"""Tests for ordered, serial, reproducible batch execution."""

from __future__ import annotations

import json
from collections.abc import Iterator

import pytest

from dnakit.batch import BatchConfig, BatchContext, BatchProgress, iter_batch, run_batch
from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import BatchExecutionError, ConfigurationError


def _records() -> list[DNARecord]:
    return [
        DNARecord(DNASequence("A"), "a"),
        DNARecord(DNASequence("CC"), "b"),
        DNARecord(DNASequence("GGG"), "c"),
    ]


def test_iter_batch_is_lazy_ordered_and_reports_progress() -> None:
    consumed: list[str] = []
    events: list[BatchProgress] = []

    def inputs() -> Iterator[DNARecord]:
        for record in _records():
            consumed.append(record.id)
            yield record

    results = iter_batch(
        inputs(),
        lambda record, context: (record.sequence.symbol_length, context.input_index),
        progress=events.append,
    )

    assert consumed == []
    first = next(results)
    assert consumed == ["a"]
    assert first.value == (1, 0)
    assert [item.value for item in results] == [(2, 1), (3, 2)]
    assert [(event.processed_count, event.last_record_id) for event in events] == [
        (1, "a"),
        (2, "b"),
        (3, "c"),
    ]


def test_seed_derivation_and_materialized_result_are_reproducible() -> None:
    def operation(_record: DNARecord, context: BatchContext) -> int | None:
        return context.seed

    config = BatchConfig(seed=42)

    first = run_batch(_records(), operation, name="seed-audit", config=config)
    second = run_batch(_records(), operation, name="seed-audit", config=config)

    assert [item.value for item in first.items] == [item.value for item in second.items]
    assert first.success_count == 3
    assert first.failure_count == 0
    payload = json.loads(json.dumps(first.to_dict()))
    assert payload["config"]["seed"] == 42
    assert payload["provenance"]["dnakit_version"]


def test_collect_policy_retains_failure_without_losing_input_order() -> None:
    def operation(record: DNARecord, _context: BatchContext) -> int:
        if record.id == "b":
            raise ValueError("bad record")
        return record.sequence.symbol_length

    result = run_batch(
        _records(),
        operation,
        name="length",
        config=BatchConfig(error_policy="collect"),
    )

    assert [item.context.record_id for item in result.items] == ["a", "b", "c"]
    assert result.success_count == 2
    assert result.failure_count == 1
    assert result.items[1].failure is not None
    assert result.items[1].failure.error_type == "ValueError"
    assert result.items[1].failure.message == "bad record"


def test_raise_policy_chains_the_original_error() -> None:
    with pytest.raises(BatchExecutionError) as exc_info:
        tuple(
            iter_batch(
                _records(),
                lambda _record, _context: 1 / 0,
            )
        )

    assert exc_info.value.code == "BATCH_ITEM_FAILED"
    assert isinstance(exc_info.value.__cause__, ZeroDivisionError)


def test_batch_limits_types_and_parallelism_are_explicit() -> None:
    with pytest.raises(ConfigurationError) as jobs_error:
        BatchConfig(jobs=2)
    assert jobs_error.value.code == "INVALID_BATCH_JOBS"
    with pytest.raises(BatchExecutionError) as limit_error:
        tuple(
            iter_batch(
                _records(),
                lambda record, _context: record.id,
                config=BatchConfig(max_records=2),
            )
        )
    assert limit_error.value.code == "BATCH_RECORD_LIMIT"
    with pytest.raises(BatchExecutionError):
        tuple(iter_batch([DNASequence("A")], lambda value, _context: value))  # type: ignore[list-item]


def test_thread_batch_preserves_order_seed_and_serial_values() -> None:
    def operation(record: DNARecord, context: BatchContext) -> tuple[str, int | None]:
        return record.id, context.seed

    serial = run_batch(_records(), operation, name="serial", config=BatchConfig(seed=42))
    threaded = run_batch(
        _records(),
        operation,
        name="threaded",
        config=BatchConfig(
            seed=42,
            jobs=3,
            execution_mode="thread",
            max_in_flight=3,
        ),
    )

    assert [item.value for item in threaded.items] == [item.value for item in serial.items]
    assert [item.context.input_index for item in threaded.items] == [0, 1, 2]


def test_thread_batch_collects_failures_in_input_order() -> None:
    def operation(record: DNARecord, _context: BatchContext) -> str:
        if record.id == "b":
            raise ValueError("bad")
        return record.id

    result = run_batch(
        _records(),
        operation,
        name="threaded-collect",
        config=BatchConfig(error_policy="collect", jobs=2, execution_mode="thread"),
    )

    assert [item.context.record_id for item in result.items] == ["a", "b", "c"]
    assert result.success_count == 2
    assert result.failure_count == 1


def test_thread_batch_raise_policy_chains_original_error() -> None:
    with pytest.raises(BatchExecutionError) as exc_info:
        tuple(
            iter_batch(
                _records(),
                lambda _record, _context: 1 / 0,
                config=BatchConfig(jobs=2, execution_mode="thread"),
            )
        )

    assert isinstance(exc_info.value.__cause__, ZeroDivisionError)


def test_batch_resume_ids_skip_already_completed_records() -> None:
    seen: list[str] = []

    def operation(record: DNARecord, _context: BatchContext) -> str:
        seen.append(record.id)
        return record.id

    result = run_batch(
        _records(),
        operation,
        name="resume",
        config=BatchConfig(resume_completed_ids=frozenset({"b"})),
    )

    assert seen == ["a", "c"]
    assert [item.context.record_id for item in result.items] == ["a", "c"]


@pytest.mark.parametrize("execution_mode,jobs", [("serial", 1), ("thread", 2)])
def test_batch_resume_preserves_original_indices_and_seeds(execution_mode: str, jobs: int) -> None:
    full = run_batch(
        _records(),
        lambda _record, context: context.seed,
        name="full",
        config=BatchConfig(seed=42, execution_mode=execution_mode, jobs=jobs),  # type: ignore[arg-type]
    )
    resumed = run_batch(
        _records(),
        lambda _record, context: context.seed,
        name="resumed",
        config=BatchConfig(
            seed=42,
            execution_mode=execution_mode,  # type: ignore[arg-type]
            jobs=jobs,
            resume_completed_ids=frozenset({"b"}),
        ),
    )

    assert [item.context.input_index for item in resumed.items] == [0, 2]
    assert [item.value for item in resumed.items] == [full.items[0].value, full.items[2].value]


def test_batch_resume_filter_does_not_bypass_validation_or_raw_input_limit() -> None:
    with pytest.raises(BatchExecutionError) as invalid:
        tuple(
            iter_batch(
                [DNASequence("A")],  # type: ignore[list-item]
                lambda value, _context: value,
                config=BatchConfig(resume_completed_ids=frozenset({"a"})),
            )
        )
    assert invalid.value.code == "INVALID_BATCH_RECORD"

    with pytest.raises(BatchExecutionError) as limit:
        tuple(
            iter_batch(
                _records(),
                lambda record, _context: record.id,
                config=BatchConfig(
                    max_records=2,
                    resume_completed_ids=frozenset({"a", "b", "c"}),
                ),
            )
        )
    assert limit.value.code == "BATCH_RECORD_LIMIT"
