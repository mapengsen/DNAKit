"""Ordered serial batch execution for DNA records."""

from __future__ import annotations

import random
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Any, Generic, Literal, TypeAlias, TypeVar, cast

from dnakit.core import DNARecord
from dnakit.core._json import to_json_compatible
from dnakit.core.provenance import Provenance
from dnakit.exceptions import BatchExecutionError, ConfigurationError, DNAKitError

T = TypeVar("T")
ErrorPolicy: TypeAlias = Literal["raise", "collect"]
ExecutionMode: TypeAlias = Literal["serial", "thread"]


@dataclass(frozen=True, slots=True)
class BatchConfig:
    """Configure bounded deterministic serial or thread-based batch execution."""

    error_policy: ErrorPolicy = "raise"
    seed: int | None = None
    jobs: int = 1
    max_records: int | None = None
    execution_mode: ExecutionMode = "serial"
    max_in_flight: int | None = None
    resume_completed_ids: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        if self.error_policy not in {"raise", "collect"}:
            raise ConfigurationError(
                "error_policy must be 'raise' or 'collect'.",
                code="INVALID_BATCH_ERROR_POLICY",
            )
        if self.seed is not None and (
            isinstance(self.seed, bool) or not isinstance(self.seed, int)
        ):
            raise ConfigurationError(
                "seed must be an integer or None.",
                code="INVALID_BATCH_SEED",
            )
        if isinstance(self.jobs, bool) or not isinstance(self.jobs, int) or self.jobs <= 0:
            raise ConfigurationError(
                "jobs must be a positive integer.",
                code="INVALID_BATCH_JOBS",
            )
        if self.execution_mode not in ("serial", "thread"):
            raise ConfigurationError(
                "execution_mode must be 'serial' or 'thread'.",
                code="INVALID_BATCH_EXECUTION_MODE",
            )
        if self.execution_mode == "serial" and self.jobs != 1:
            raise ConfigurationError(
                "Serial execution requires jobs=1.",
                code="INVALID_BATCH_JOBS",
                hint="Use execution_mode='thread' for concurrent I/O or external work.",
            )
        if self.max_records is not None and (
            isinstance(self.max_records, bool)
            or not isinstance(self.max_records, int)
            or self.max_records <= 0
        ):
            raise ConfigurationError(
                "max_records must be a positive integer or None.",
                code="INVALID_BATCH_LIMIT",
            )
        if self.max_in_flight is not None and (
            isinstance(self.max_in_flight, bool)
            or not isinstance(self.max_in_flight, int)
            or self.max_in_flight <= 0
        ):
            raise ConfigurationError(
                "max_in_flight must be a positive integer or None.",
                code="INVALID_BATCH_IN_FLIGHT_LIMIT",
            )
        if (
            self.execution_mode == "thread"
            and self.max_in_flight is not None
            and self.max_in_flight < self.jobs
        ):
            raise ConfigurationError(
                "max_in_flight must be at least jobs.",
                code="INVALID_BATCH_IN_FLIGHT_LIMIT",
            )
        if not isinstance(self.resume_completed_ids, frozenset) or any(
            not isinstance(item, str) or not item.strip() for item in self.resume_completed_ids
        ):
            raise ConfigurationError(
                "resume_completed_ids must be a frozenset of non-empty record IDs.",
                code="INVALID_BATCH_RESUME_IDS",
            )


@dataclass(frozen=True, slots=True)
class BatchContext:
    """Stable per-record execution context supplied to a batch operation."""

    input_index: int
    record_id: str
    seed: int | None

    def __post_init__(self) -> None:
        if (
            isinstance(self.input_index, bool)
            or not isinstance(self.input_index, int)
            or self.input_index < 0
        ):
            raise ConfigurationError("BatchContext input_index must be non-negative.")
        if not isinstance(self.record_id, str) or not self.record_id.strip():
            raise ConfigurationError("BatchContext record_id must be non-empty.")
        if self.seed is not None and (
            isinstance(self.seed, bool) or not isinstance(self.seed, int)
        ):
            raise ConfigurationError("BatchContext seed must be an integer or None.")


@dataclass(frozen=True, slots=True)
class BatchFailure:
    """Serializable details for one collected operation error."""

    error_type: str
    code: str | None
    message: str
    hint: str | None

    def __post_init__(self) -> None:
        for name in ("error_type", "message"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"BatchFailure {name} must be non-empty.")


@dataclass(frozen=True, slots=True)
class BatchItemResult(Generic[T]):
    """One ordered success or collected failure."""

    context: BatchContext
    value: T | None = None
    failure: BatchFailure | None = None

    def __post_init__(self) -> None:
        if self.value is not None and self.failure is not None:
            raise ConfigurationError(
                "A failed batch item cannot also contain a value.",
                code="INVALID_BATCH_ITEM_RESULT",
            )

    @property
    def succeeded(self) -> bool:
        return self.failure is None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible payload when the operation value supports it."""

        value: object = self.value
        value_to_dict = getattr(value, "to_dict", None)
        if callable(value_to_dict):
            value = value_to_dict()
        return cast(
            dict[str, Any],
            to_json_compatible({"context": self.context, "value": value, "failure": self.failure}),
        )


@dataclass(frozen=True, slots=True)
class BatchProgress:
    """Progress event emitted only through an explicit caller callback."""

    processed_count: int
    success_count: int
    failure_count: int
    last_record_id: str


@dataclass(frozen=True, slots=True)
class _ThreadOutcome(Generic[T]):
    item: BatchItemResult[T]
    error: Exception | None = None


@dataclass(frozen=True, slots=True)
class _PreparedRecord:
    """Validated input retaining its original position and derived seed."""

    record: DNARecord
    context: BatchContext


@dataclass(frozen=True, slots=True)
class BatchResult(Generic[T]):
    """Materialized ordered batch result with resolved execution settings."""

    operation: str
    items: tuple[BatchItemResult[T], ...]
    config: BatchConfig
    provenance: Provenance = field(default_factory=Provenance)

    def __post_init__(self) -> None:
        if not isinstance(self.operation, str) or not self.operation.strip():
            raise ConfigurationError(
                "operation must be a non-empty string.",
                code="INVALID_BATCH_OPERATION",
            )

    @property
    def success_count(self) -> int:
        return sum(item.succeeded for item in self.items)

    @property
    def failure_count(self) -> int:
        return len(self.items) - self.success_count

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible audit payload."""

        return {
            "operation": self.operation,
            "items": [item.to_dict() for item in self.items],
            "config": cast(dict[str, Any], to_json_compatible(self.config)),
            "success_count": self.success_count,
            "failure_count": self.failure_count,
            "provenance": self.provenance.to_dict(),
        }


BatchOperation: TypeAlias = Callable[[DNARecord, BatchContext], T]
ProgressCallback: TypeAlias = Callable[[BatchProgress], None]


def _failure(error: Exception) -> BatchFailure:
    return BatchFailure(
        error_type=type(error).__name__,
        code=error.code if isinstance(error, DNAKitError) else None,
        message=str(error),
        hint=error.hint if isinstance(error, DNAKitError) else None,
    )


def _prepare_records(
    iterator: Iterator[object],
    config: BatchConfig,
) -> Iterator[_PreparedRecord]:
    """Validate and seed every raw input before applying the resume filter.

    Deriving seeds from the unfiltered input stream guarantees that resuming a
    run produces the same context for every remaining record as a fresh run.
    ``max_records`` also bounds the raw stream, so an infinite stream containing
    only completed IDs cannot bypass the configured limit.
    """

    seed_generator = random.Random(config.seed) if config.seed is not None else None
    for input_index, record in enumerate(iterator):
        if config.max_records is not None and input_index >= config.max_records:
            raise BatchExecutionError(
                "Batch input exceeds max_records.",
                code="BATCH_RECORD_LIMIT",
                context={
                    "max_records": config.max_records,
                    "next_input_index": input_index,
                },
            )
        if not isinstance(record, DNARecord):
            raise BatchExecutionError(
                "Every batch input must be a DNARecord.",
                code="INVALID_BATCH_RECORD",
                context={"input_index": input_index, "type": type(record).__name__},
            )
        context = BatchContext(
            input_index=input_index,
            record_id=record.id,
            seed=(seed_generator.getrandbits(64) if seed_generator is not None else None),
        )
        if record.id not in config.resume_completed_ids:
            yield _PreparedRecord(record, context)


def iter_batch(
    records: Iterable[DNARecord],
    operation: BatchOperation[T],
    *,
    config: BatchConfig | None = None,
    progress: ProgressCallback | None = None,
) -> Iterator[BatchItemResult[T]]:
    """Lazily execute ``operation`` once per record in stable input order."""

    resolved = BatchConfig() if config is None else config
    if not isinstance(resolved, BatchConfig):
        raise ConfigurationError("config must be BatchConfig or None.", code="INVALID_BATCH_CONFIG")
    if not callable(operation):
        raise ConfigurationError("operation must be callable.", code="INVALID_BATCH_OPERATION")
    if progress is not None and not callable(progress):
        raise ConfigurationError(
            "progress must be callable or None.", code="INVALID_BATCH_PROGRESS"
        )
    try:
        iterator = iter(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be an iterable of DNARecord objects.",
            code="INVALID_BATCH_RECORDS",
        ) from exc
    prepared = _prepare_records(cast(Iterator[object], iterator), resolved)
    if resolved.execution_mode == "thread":
        yield from _iter_batch_threaded(
            prepared,
            operation,
            config=resolved,
            progress=progress,
        )
        return
    successes = 0
    failures = 0
    for processed, prepared_record in enumerate(prepared, start=1):
        record = prepared_record.record
        context = prepared_record.context
        item: BatchItemResult[T]
        try:
            item = BatchItemResult(context=context, value=operation(record, context))
            successes += 1
        except Exception as exc:
            if resolved.error_policy == "raise":
                raise BatchExecutionError(
                    "Batch operation failed.",
                    code="BATCH_ITEM_FAILED",
                    context={"input_index": context.input_index, "record_id": record.id},
                    hint="Use error_policy='collect' to retain per-record failures.",
                ) from exc
            item = BatchItemResult(context=context, failure=_failure(exc))
            failures += 1
        if progress is not None:
            progress(BatchProgress(processed, successes, failures, record.id))
        yield item


def _thread_item(
    record: DNARecord,
    context: BatchContext,
    operation: BatchOperation[T],
) -> _ThreadOutcome[T]:
    try:
        return _ThreadOutcome(BatchItemResult(context=context, value=operation(record, context)))
    except Exception as exc:
        return _ThreadOutcome(
            BatchItemResult(context=context, failure=_failure(exc)),
            error=exc,
        )


def _iter_batch_threaded(
    iterator: Iterator[_PreparedRecord],
    operation: BatchOperation[T],
    *,
    config: BatchConfig,
    progress: ProgressCallback | None,
) -> Iterator[BatchItemResult[T]]:
    """Execute bounded concurrent work while yielding stable input order."""

    limit = config.max_in_flight or max(config.jobs, config.jobs * 2)
    next_submit = 0
    next_yield = 0
    successes = 0
    failures = 0
    exhausted = False
    pending: dict[Future[_ThreadOutcome[T]], int] = {}
    completed: dict[int, _ThreadOutcome[T]] = {}
    with ThreadPoolExecutor(max_workers=config.jobs, thread_name_prefix="dnakit") as executor:
        while not exhausted or pending or completed:
            while not exhausted and len(pending) + len(completed) < limit:
                try:
                    prepared_record = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                future = executor.submit(
                    _thread_item,
                    prepared_record.record,
                    prepared_record.context,
                    operation,
                )
                pending[future] = next_submit
                next_submit += 1
            if pending:
                done = next(as_completed(tuple(pending)))
                completed[pending.pop(done)] = done.result()
            while next_yield in completed:
                outcome = completed.pop(next_yield)
                item = outcome.item
                if item.failure is not None and config.error_policy == "raise":
                    for future in pending:
                        future.cancel()
                    batch_error = BatchExecutionError(
                        "Batch operation failed.",
                        code="BATCH_ITEM_FAILED",
                        context={
                            "input_index": item.context.input_index,
                            "record_id": item.context.record_id,
                            "error_type": item.failure.error_type,
                        },
                        hint="Use error_policy='collect' to retain per-record failures.",
                    )
                    if outcome.error is None:
                        raise batch_error
                    raise batch_error from outcome.error
                if item.failure is None:
                    successes += 1
                else:
                    failures += 1
                if progress is not None:
                    progress(
                        BatchProgress(
                            next_yield + 1,
                            successes,
                            failures,
                            item.context.record_id,
                        )
                    )
                yield item
                next_yield += 1


def run_batch(
    records: Iterable[DNARecord],
    operation: BatchOperation[T],
    *,
    name: str,
    config: BatchConfig | None = None,
    progress: ProgressCallback | None = None,
) -> BatchResult[T]:
    """Materialize an ordered serial batch execution."""

    resolved = BatchConfig() if config is None else config
    items = tuple(iter_batch(records, operation, config=resolved, progress=progress))
    return BatchResult(name, items, resolved)


__all__ = [
    "BatchConfig",
    "BatchContext",
    "BatchFailure",
    "BatchItemResult",
    "BatchOperation",
    "BatchProgress",
    "BatchResult",
    "ErrorPolicy",
    "ProgressCallback",
    "iter_batch",
    "run_batch",
]
