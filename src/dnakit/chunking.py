"""Coordinate-aware DNA sequence chunking strategies.

The module keeps sequence extraction lazy and records enough provenance to
trace every chunk back to its source record and interval.  A BED file is
optional: without one, every input FASTA record is treated as a ``train``
region by default.  BED intervals use the standard zero-based, half-open
coordinate convention already used by DNAKit.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Any, Literal, TypeAlias

from dnakit.core.coordinates import Interval
from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, CoordinateError, InputFormatError
from dnakit.io.annotations import AnnotationDocument, AnnotationSource, read_bed
from dnakit.io.api import ReadableSource, read
from dnakit.io.config import ReadConfig
from dnakit.ops.edit import subsequence

ChunkStrategy = Literal["fixed", "sliding", "random", "multiscale", "curriculum"]
ChunkProgressStatus = Literal["started", "yielded", "completed"]
ChunkProgressCallback: TypeAlias = Callable[["ChunkProgress"], None]
ChunkBedSource: TypeAlias = AnnotationDocument | AnnotationSource


def _positive_int(value: object, field_name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(
            f"{field_name} must be a positive integer.",
            code="INVALID_CHUNKING_CONFIG",
            context={"field": field_name, "value": value},
        )


def _non_negative_int(value: object, field_name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(
            f"{field_name} must be a non-negative integer.",
            code="INVALID_CHUNKING_CONFIG",
            context={"field": field_name, "value": value},
        )


def _coerce_int_tuple(value: Sequence[int], field_name: str) -> tuple[int, ...]:
    if isinstance(value, (str, bytes)):
        raise ConfigurationError(
            f"{field_name} must be a sequence of integers.",
            code="INVALID_CHUNKING_CONFIG",
            context={"field": field_name},
        )
    try:
        resolved = tuple(value)
    except TypeError as exc:
        raise ConfigurationError(
            f"{field_name} must be a sequence of integers.",
            code="INVALID_CHUNKING_CONFIG",
            context={"field": field_name},
        ) from exc
    for index, item in enumerate(resolved):
        _positive_int(item, f"{field_name}[{index}]")
    return resolved


def _validate_text(value: object, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip() or any(char.isspace() for char in value):
        raise ConfigurationError(
            f"{field_name} must be non-empty text without whitespace.",
            code="INVALID_CHUNKING_CONFIG",
            context={"field": field_name, "value": value},
        )


@dataclass(frozen=True, slots=True)
class ChunkingConfig:
    """Validate one sequence chunking strategy.

    ``fixed`` creates non-overlapping windows.  ``sliding`` uses ``step`` and
    defaults to a one-base step.  ``random`` requires a length range and a
    sample count.  ``multiscale`` creates windows for every value in
    ``lengths``.  ``curriculum`` has the same window behavior as multiscale,
    while its optional ``stage_steps`` records how long a training stage is
    intended to run.
    """

    strategy: ChunkStrategy = "fixed"
    length: int = 1024
    step: int | None = None
    min_length: int | None = None
    max_length: int | None = None
    num_samples: int | None = None
    lengths: tuple[int, ...] = ()
    steps: tuple[int, ...] | None = None
    stage_steps: tuple[int, ...] | None = None
    include_partial: bool = False
    seed: int | None = None
    split: str = "train"
    allow_gaps: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.strategy, str) or self.strategy not in {
            "fixed",
            "sliding",
            "random",
            "multiscale",
            "curriculum",
        }:
            raise ConfigurationError(
                "Unknown sequence chunking strategy.",
                code="INVALID_CHUNKING_STRATEGY",
                context={"strategy": self.strategy},
            )
        _positive_int(self.length, "length")
        if self.step is not None:
            _positive_int(self.step, "step")
        for field_name, value in (
            ("min_length", self.min_length),
            ("max_length", self.max_length),
            ("num_samples", self.num_samples),
        ):
            if value is not None:
                _positive_int(value, field_name)
        if self.seed is not None and (
            isinstance(self.seed, bool) or not isinstance(self.seed, int)
        ):
            raise ConfigurationError(
                "seed must be an integer or None.",
                code="INVALID_CHUNKING_CONFIG",
                context={"field": "seed", "value": self.seed},
            )
        if not isinstance(self.include_partial, bool):
            raise ConfigurationError(
                "include_partial must be boolean.",
                code="INVALID_CHUNKING_CONFIG",
                context={"field": "include_partial", "value": self.include_partial},
            )
        if not isinstance(self.allow_gaps, bool):
            raise ConfigurationError(
                "allow_gaps must be boolean.",
                code="INVALID_CHUNKING_CONFIG",
                context={"field": "allow_gaps", "value": self.allow_gaps},
            )
        _validate_text(self.split, "split")

        resolved_lengths = _coerce_int_tuple(self.lengths, "lengths")
        object.__setattr__(self, "lengths", resolved_lengths)
        resolved_steps = None if self.steps is None else _coerce_int_tuple(self.steps, "steps")
        object.__setattr__(self, "steps", resolved_steps)
        resolved_stage_steps = (
            None if self.stage_steps is None else _coerce_int_tuple(self.stage_steps, "stage_steps")
        )
        object.__setattr__(self, "stage_steps", resolved_stage_steps)

        if self.strategy == "fixed":
            if self.step is not None:
                raise ConfigurationError(
                    "fixed strategy is non-overlapping; omit step.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if resolved_lengths or resolved_steps or resolved_stage_steps:
                raise ConfigurationError(
                    "lengths, steps, and stage_steps are only for multiscale or curriculum.",
                    code="INVALID_CHUNKING_CONFIG",
                )
        elif self.strategy == "sliding":
            resolved_step = 1 if self.step is None else self.step
            if resolved_step > self.length:
                raise ConfigurationError(
                    "sliding step cannot exceed window length; that would create gaps.",
                    code="INVALID_CHUNKING_CONFIG",
                    context={"length": self.length, "step": resolved_step},
                )
            object.__setattr__(self, "step", resolved_step)
            if resolved_lengths or resolved_steps or resolved_stage_steps:
                raise ConfigurationError(
                    "lengths, steps, and stage_steps are only for multiscale or curriculum.",
                    code="INVALID_CHUNKING_CONFIG",
                )
        elif self.strategy == "random":
            if self.min_length is None or self.max_length is None or self.num_samples is None:
                raise ConfigurationError(
                    "random strategy requires min_length, max_length, and num_samples.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if self.min_length > self.max_length:
                raise ConfigurationError(
                    "min_length cannot exceed max_length.",
                    code="INVALID_CHUNKING_CONFIG",
                    context={
                        "min_length": self.min_length,
                        "max_length": self.max_length,
                    },
                )
            if resolved_lengths or resolved_steps or resolved_stage_steps:
                raise ConfigurationError(
                    "lengths, steps, and stage_steps are only for multiscale or curriculum.",
                    code="INVALID_CHUNKING_CONFIG",
                )
        else:
            if not resolved_lengths:
                raise ConfigurationError(
                    f"{self.strategy} strategy requires at least one length.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if any(left >= right for left, right in pairwise(resolved_lengths)):
                raise ConfigurationError(
                    "multiscale and curriculum lengths must be strictly increasing.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if resolved_steps is not None and len(resolved_steps) != len(resolved_lengths):
                raise ConfigurationError(
                    "steps must have the same number of values as lengths.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if resolved_steps is not None and any(
                window_step > window_length
                for window_step, window_length in zip(resolved_steps, resolved_lengths, strict=True)
            ):
                raise ConfigurationError(
                    "multiscale or curriculum steps cannot exceed their window lengths.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if self.strategy == "multiscale" and resolved_stage_steps is not None:
                raise ConfigurationError(
                    "stage_steps are only valid for curriculum.",
                    code="INVALID_CHUNKING_CONFIG",
                )
            if (
                self.strategy == "curriculum"
                and resolved_stage_steps is not None
                and len(resolved_stage_steps) != len(resolved_lengths)
            ):
                raise ConfigurationError(
                    "stage_steps must have the same number of values as lengths.",
                    code="INVALID_CHUNKING_CONFIG",
                )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-friendly configuration summary."""

        return {
            "strategy": self.strategy,
            "length": self.length,
            "step": self.step,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "num_samples": self.num_samples,
            "lengths": list(self.lengths),
            "steps": None if self.steps is None else list(self.steps),
            "stage_steps": None if self.stage_steps is None else list(self.stage_steps),
            "include_partial": self.include_partial,
            "seed": self.seed,
            "split": self.split,
            "allow_gaps": self.allow_gaps,
        }


@dataclass(frozen=True, slots=True)
class LengthCurriculumStage:
    """One stage in a short-to-long training schedule."""

    index: int
    length: int
    training_steps: int | None = None

    def __post_init__(self) -> None:
        _non_negative_int(self.index, "index")
        _positive_int(self.length, "length")
        if self.training_steps is not None:
            _positive_int(self.training_steps, "training_steps")


@dataclass(frozen=True, slots=True)
class LengthCurriculum:
    """A validated short-to-long sequence length schedule.

    The schedule itself does not start training.  Use :meth:`to_config` to
    obtain a ``curriculum`` chunking configuration, then pass that config to
    :func:`iter_fasta_chunks` or :func:`iter_sequence_chunks`.
    """

    lengths: tuple[int, ...]
    stage_steps: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        resolved_lengths = _coerce_int_tuple(self.lengths, "lengths")
        if not resolved_lengths:
            raise ConfigurationError(
                "A length curriculum requires at least one length.",
                code="INVALID_LENGTH_CURRICULUM",
            )
        if any(left >= right for left, right in pairwise(resolved_lengths)):
            raise ConfigurationError(
                "Curriculum lengths must be strictly increasing.",
                code="INVALID_LENGTH_CURRICULUM",
            )
        object.__setattr__(self, "lengths", resolved_lengths)
        resolved_stage_steps = (
            None if self.stage_steps is None else _coerce_int_tuple(self.stage_steps, "stage_steps")
        )
        if resolved_stage_steps is not None and len(resolved_stage_steps) != len(resolved_lengths):
            raise ConfigurationError(
                "stage_steps must have the same number of values as lengths.",
                code="INVALID_LENGTH_CURRICULUM",
            )
        object.__setattr__(self, "stage_steps", resolved_stage_steps)

    @property
    def stages(self) -> tuple[LengthCurriculumStage, ...]:
        """Return stages in the order in which they should be trained."""

        return tuple(
            LengthCurriculumStage(
                index=index,
                length=length,
                training_steps=None if self.stage_steps is None else self.stage_steps[index],
            )
            for index, length in enumerate(self.lengths)
        )

    def to_config(
        self,
        *,
        window_step: int | None = None,
        include_partial: bool = False,
        seed: int | None = None,
        split: str = "train",
        allow_gaps: bool = False,
    ) -> ChunkingConfig:
        """Convert the schedule into a chunking configuration."""

        resolved_steps = None if window_step is None else tuple(window_step for _ in self.lengths)
        return ChunkingConfig(
            strategy="curriculum",
            lengths=self.lengths,
            steps=resolved_steps,
            stage_steps=self.stage_steps,
            include_partial=include_partial,
            seed=seed,
            split=split,
            allow_gaps=allow_gaps,
        )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-friendly curriculum summary."""

        return {
            "lengths": list(self.lengths),
            "stage_steps": None if self.stage_steps is None else list(self.stage_steps),
            "stages": [
                {
                    "index": stage.index,
                    "length": stage.length,
                    "training_steps": stage.training_steps,
                }
                for stage in self.stages
            ],
        }


def make_length_curriculum(
    lengths: Sequence[int], *, stage_steps: Sequence[int] | None = None
) -> LengthCurriculum:
    """Create a validated short-to-long curriculum."""

    return LengthCurriculum(tuple(lengths), None if stage_steps is None else tuple(stage_steps))


@dataclass(frozen=True, slots=True)
class ChunkProgress:
    """Advisory progress event emitted while chunks are generated."""

    source_id: str
    split: str
    region_index: int
    processed: int
    total: int
    status: ChunkProgressStatus


@dataclass(frozen=True, slots=True)
class SequenceChunk:
    """One extracted sequence and its source-coordinate provenance."""

    sequence: DNASequence
    id: str
    source_id: str
    split: str
    source_start: int
    source_end: int
    requested_length: int
    window_index: int
    strategy: ChunkStrategy
    region_index: int = 0
    level_index: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.sequence, DNASequence):
            raise ConfigurationError(
                "SequenceChunk.sequence must be a DNASequence.",
                code="INVALID_SEQUENCE_CHUNK",
            )
        for field_name, value in (
            ("id", self.id),
            ("source_id", self.source_id),
            ("split", self.split),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(
                    f"SequenceChunk {field_name} must be non-empty text.",
                    code="INVALID_SEQUENCE_CHUNK",
                )
        for integer_field_name, integer_value in (
            ("source_start", self.source_start),
            ("source_end", self.source_end),
            ("requested_length", self.requested_length),
        ):
            if integer_field_name == "requested_length":
                _positive_int(integer_value, integer_field_name)
            else:
                _non_negative_int(integer_value, integer_field_name)
        if self.source_end <= self.source_start:
            raise ConfigurationError(
                "SequenceChunk source_end must be greater than source_start.",
                code="INVALID_SEQUENCE_CHUNK",
            )
        if self.source_end - self.source_start > self.requested_length:
            raise ConfigurationError(
                "Extracted sequence cannot exceed requested_length.",
                code="INVALID_SEQUENCE_CHUNK",
            )
        sequence_span = self.sequence.coordinate_span
        if sequence_span is not None and sequence_span != self.source_end - self.source_start:
            raise ConfigurationError(
                "SequenceChunk coordinates do not match the extracted sequence span.",
                code="INVALID_SEQUENCE_CHUNK",
                context={
                    "sequence_span": sequence_span,
                    "source_start": self.source_start,
                    "source_end": self.source_end,
                },
            )
        _non_negative_int(self.window_index, "window_index")
        _non_negative_int(self.region_index, "region_index")
        if self.level_index is not None:
            _non_negative_int(self.level_index, "level_index")
        if not isinstance(self.strategy, str) or self.strategy not in {
            "fixed",
            "sliding",
            "random",
            "multiscale",
            "curriculum",
        }:
            raise ConfigurationError(
                "SequenceChunk strategy is not supported.",
                code="INVALID_SEQUENCE_CHUNK",
                context={"strategy": self.strategy},
            )

    @property
    def length(self) -> int:
        """Actual sequence length in coordinates."""

        return self.source_end - self.source_start

    def to_record(self) -> DNARecord:
        """Convert the chunk to a DNARecord while retaining provenance."""

        return DNARecord(
            self.sequence,
            self.id,
            description=(
                f"{self.source_id}:{self.source_start}-{self.source_end}"
                f" split={self.split} strategy={self.strategy}"
            ),
            metadata={
                "chunking": {
                    "strategy": self.strategy,
                    "split": self.split,
                    "source_id": self.source_id,
                    "source_start": self.source_start,
                    "source_end": self.source_end,
                    "requested_length": self.requested_length,
                    "window_index": self.window_index,
                    "region_index": self.region_index,
                    "level_index": self.level_index,
                }
            },
        )

    def to_dict(self) -> dict[str, Any]:
        """Return chunk metadata and symbols in a serialization-friendly form."""

        return {
            "id": self.id,
            "source_id": self.source_id,
            "split": self.split,
            "source_start": self.source_start,
            "source_end": self.source_end,
            "length": self.length,
            "requested_length": self.requested_length,
            "window_index": self.window_index,
            "region_index": self.region_index,
            "level_index": self.level_index,
            "strategy": self.strategy,
            "sequence": self.sequence.symbols,
        }


@dataclass(frozen=True, slots=True)
class _PlannedWindow:
    start: int
    end: int
    requested_length: int
    level_index: int | None = None


def _regular_windows(
    start: int,
    end: int,
    *,
    length: int,
    step: int,
    include_partial: bool,
    level_index: int | None = None,
) -> Iterator[_PlannedWindow]:
    if end <= start:
        return
    if include_partial:
        for window_start in range(start, end, step):
            window_end = min(window_start + length, end)
            yield _PlannedWindow(window_start, window_end, length, level_index)
        return
    last_start = end - length
    if last_start < start:
        return
    for window_start in range(start, last_start + 1, step):
        yield _PlannedWindow(window_start, window_start + length, length, level_index)


def _random_windows(
    start: int,
    end: int,
    *,
    config: ChunkingConfig,
    rng: random.Random,
) -> Iterator[_PlannedWindow]:
    assert config.min_length is not None
    assert config.max_length is not None
    assert config.num_samples is not None
    available = end - start
    if available < config.min_length:
        raise ConfigurationError(
            "The requested random minimum length exceeds the input region.",
            code="RANDOM_CHUNK_TOO_LONG",
            context={"region_length": available, "min_length": config.min_length},
        )
    upper_length = min(config.max_length, available)
    for _ in range(config.num_samples):
        window_length = rng.randint(config.min_length, upper_length)
        window_start = rng.randint(start, end - window_length)
        yield _PlannedWindow(window_start, window_start + window_length, window_length)


def _planned_windows(
    start: int,
    end: int,
    *,
    config: ChunkingConfig,
    rng: random.Random,
) -> Iterator[_PlannedWindow]:
    if config.strategy == "fixed":
        yield from _regular_windows(
            start,
            end,
            length=config.length,
            step=config.length,
            include_partial=config.include_partial,
        )
    elif config.strategy == "sliding":
        assert config.step is not None
        yield from _regular_windows(
            start,
            end,
            length=config.length,
            step=config.step,
            include_partial=config.include_partial,
        )
    elif config.strategy == "random":
        yield from _random_windows(start, end, config=config, rng=rng)
    else:
        for level_index, window_length in enumerate(config.lengths):
            window_step = window_length if config.steps is None else config.steps[level_index]
            yield from _regular_windows(
                start,
                end,
                length=window_length,
                step=window_step,
                include_partial=config.include_partial,
                level_index=level_index,
            )


def _count_regular_windows(
    region_length: int, *, length: int, step: int, include_partial: bool
) -> int:
    if region_length <= 0:
        return 0
    if include_partial:
        return (region_length + step - 1) // step
    if region_length < length:
        return 0
    return (region_length - length) // step + 1


def _planned_window_count(start: int, end: int, *, config: ChunkingConfig) -> int:
    region_length = end - start
    if config.strategy == "fixed":
        return _count_regular_windows(
            region_length,
            length=config.length,
            step=config.length,
            include_partial=config.include_partial,
        )
    if config.strategy == "sliding":
        assert config.step is not None
        return _count_regular_windows(
            region_length,
            length=config.length,
            step=config.step,
            include_partial=config.include_partial,
        )
    if config.strategy == "random":
        assert config.min_length is not None
        assert config.num_samples is not None
        if region_length < config.min_length:
            raise ConfigurationError(
                "The requested random minimum length exceeds the input region.",
                code="RANDOM_CHUNK_TOO_LONG",
                context={"region_length": region_length, "min_length": config.min_length},
            )
        return config.num_samples
    total = 0
    for level_index, window_length in enumerate(config.lengths):
        window_step = window_length if config.steps is None else config.steps[level_index]
        total += _count_regular_windows(
            region_length,
            length=window_length,
            step=window_step,
            include_partial=config.include_partial,
        )
    return total


def _emit_progress(callback: ChunkProgressCallback | None, event: ChunkProgress) -> None:
    if callback is None:
        return
    try:
        callback(event)
    except Exception:
        # Progress is advisory and must not change chunk generation semantics.
        return


def _validate_region(sequence: DNASequence, start: int, end: int | None) -> int:
    span = len(sequence)
    if isinstance(start, bool) or not isinstance(start, int):
        raise CoordinateError(
            "Chunk start must be an integer.",
            code="INVALID_CHUNK_COORDINATE",
            context={"start": start},
        )
    resolved_end = span if end is None else end
    if isinstance(resolved_end, bool) or not isinstance(resolved_end, int):
        raise CoordinateError(
            "Chunk end must be an integer or None.",
            code="INVALID_CHUNK_COORDINATE",
            context={"end": resolved_end},
        )
    if start < 0 or resolved_end < start or resolved_end > span:
        raise CoordinateError(
            "Chunk coordinates must be a 0-based half-open interval inside the sequence.",
            code="INVALID_CHUNK_COORDINATE",
            context={"start": start, "end": resolved_end, "sequence_length": span},
        )
    return resolved_end


def _record_identity(
    value: DNA | DNASequence | DNARecord, *, source_id: str | None
) -> tuple[DNASequence, str]:
    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        default_id = "sequence_1" if record is None else record.id
        return sequence, default_id if source_id is None else source_id
    raise TypeError("value must be DNA, DNASequence, or DNARecord.")


def _chunk_id(
    source_id: str,
    split: str,
    region_index: int,
    window_index: int,
    level_index: int | None,
) -> str:
    level = "" if level_index is None else f"|level_{level_index}"
    return f"{source_id}|{split}|region_{region_index}{level}|chunk_{window_index}"


def _iter_region_chunks(
    sequence: DNASequence,
    *,
    source_id: str,
    split: str,
    start: int,
    end: int,
    region_index: int,
    config: ChunkingConfig,
    rng: random.Random,
    progress: ChunkProgressCallback | None,
) -> Iterator[SequenceChunk]:
    total = _planned_window_count(start, end, config=config)
    _emit_progress(
        progress,
        ChunkProgress(source_id, split, region_index, 0, total, "started"),
    )
    processed = 0
    for plan in _planned_windows(start, end, config=config, rng=rng):
        selected = subsequence(sequence, plan.start, plan.end, allow_gaps=config.allow_gaps)
        chunk = SequenceChunk(
            sequence=selected,
            id=_chunk_id(source_id, split, region_index, processed, plan.level_index),
            source_id=source_id,
            split=split,
            source_start=plan.start,
            source_end=plan.end,
            requested_length=plan.requested_length,
            window_index=processed,
            strategy=config.strategy,
            region_index=region_index,
            level_index=plan.level_index,
        )
        processed += 1
        _emit_progress(
            progress,
            ChunkProgress(source_id, split, region_index, processed, total, "yielded"),
        )
        yield chunk
    _emit_progress(
        progress,
        ChunkProgress(source_id, split, region_index, processed, total, "completed"),
    )


def iter_sequence_chunks(
    value: DNA | DNASequence | DNARecord,
    *,
    config: ChunkingConfig | None = None,
    start: int = 0,
    end: int | None = None,
    source_id: str | None = None,
    split: str | None = None,
    region_index: int = 0,
    progress: ChunkProgressCallback | None = None,
) -> Iterator[SequenceChunk]:
    """Lazily split one sequence or record according to ``config``.

    The default is 1024-base, non-overlapping windows assigned to ``train``.
    Coordinates are zero-based and half-open.  Partial windows are discarded
    unless ``ChunkingConfig.include_partial`` is enabled.
    """

    resolved_config = ChunkingConfig() if config is None else config
    if not isinstance(resolved_config, ChunkingConfig):
        raise TypeError("config must be ChunkingConfig or None.")
    _non_negative_int(region_index, "region_index")
    sequence, resolved_source_id = _record_identity(value, source_id=source_id)
    resolved_split = resolved_config.split if split is None else split
    _validate_text(resolved_source_id, "source_id")
    _validate_text(resolved_split, "split")
    resolved_end = _validate_region(sequence, start, end)
    yield from _iter_region_chunks(
        sequence,
        source_id=resolved_source_id,
        split=resolved_split,
        start=start,
        end=resolved_end,
        region_index=region_index,
        config=resolved_config,
        rng=random.Random(resolved_config.seed),
        progress=progress,
    )


def _bed_regions(
    document: AnnotationDocument, *, default_split: str
) -> dict[str, tuple[tuple[int, int, str], ...]]:
    if document.format != "bed":
        raise InputFormatError(
            "Sequence chunking accepts a BED annotation document.",
            code="CHUNKING_BED_REQUIRED",
            context={"format": document.format},
        )
    grouped: dict[str, list[tuple[int, int, str]]] = {}
    for entry in document.entries:
        location = entry.feature.location
        if not isinstance(location, Interval):
            raise InputFormatError(
                "BED chunking requires simple intervals.",
                code="INVALID_BED_INTERVAL",
                context={"sequence_id": entry.sequence_id},
            )
        split = entry.feature.id or entry.feature.label or default_split
        _validate_text(split, "BED split label")
        grouped.setdefault(entry.sequence_id, []).append((location.start, location.end, split))
    return {sequence_id: tuple(regions) for sequence_id, regions in grouped.items()}


def _resolve_bed(bed: ChunkBedSource) -> AnnotationDocument:
    if isinstance(bed, AnnotationDocument):
        return bed
    return read_bed(bed)


def iter_fasta_chunks(
    source: ReadableSource,
    *,
    config: ChunkingConfig | None = None,
    bed: ChunkBedSource | None = None,
    read_config: ReadConfig | None = None,
    progress: ChunkProgressCallback | None = None,
) -> Iterator[SequenceChunk]:
    """Lazily split FASTA records, optionally restricted and labelled by BED.

    Without ``bed``, every FASTA record is processed from coordinate 0 to its
    end and uses ``config.split`` (``train`` by default).  With BED, each BED
    row selects one interval; its fourth column becomes the output split label
    (for example ``train``, ``valid``, or ``test``).  FASTA IDs and BED column
    1 values must match exactly.

    This function uses DNAKit's record reader, so one FASTA record is held in
    memory at a time.  For chromosome-scale records, pass a ``ReadConfig`` with
    sufficiently high input and sequence limits and plan memory accordingly.
    """

    resolved_config = ChunkingConfig() if config is None else config
    if not isinstance(resolved_config, ChunkingConfig):
        raise TypeError("config must be ChunkingConfig or None.")
    if read_config is not None and not isinstance(read_config, ReadConfig):
        raise TypeError("read_config must be ReadConfig or None.")
    grouped_bed = (
        None
        if bed is None
        else _bed_regions(_resolve_bed(bed), default_split=resolved_config.split)
    )
    seen_ids: set[str] = set()
    rng = random.Random(resolved_config.seed)
    with read(source, format="fasta", config=read_config) as records:
        for record in records:
            regions = None if grouped_bed is None else grouped_bed.get(record.id, ())
            if grouped_bed is not None and regions:
                seen_ids.add(record.id)
            if grouped_bed is None:
                resolved_end = _validate_region(record.sequence, 0, None)
                yield from _iter_region_chunks(
                    record.sequence,
                    source_id=record.id,
                    split=resolved_config.split,
                    start=0,
                    end=resolved_end,
                    region_index=0,
                    config=resolved_config,
                    rng=rng,
                    progress=progress,
                )
                continue
            assert regions is not None
            for region_index, (start, end, split) in enumerate(regions):
                resolved_end = _validate_region(record.sequence, start, end)
                yield from _iter_region_chunks(
                    record.sequence,
                    source_id=record.id,
                    split=split,
                    start=start,
                    end=resolved_end,
                    region_index=region_index,
                    config=resolved_config,
                    rng=rng,
                    progress=progress,
                )
    if grouped_bed is not None:
        missing = tuple(sequence_id for sequence_id in grouped_bed if sequence_id not in seen_ids)
        if missing:
            raise InputFormatError(
                "BED contains sequence IDs that were not found in FASTA.",
                code="BED_SEQUENCE_NOT_FOUND",
                context={"missing_sequence_ids": missing},
            )


__all__ = [
    "ChunkBedSource",
    "ChunkProgress",
    "ChunkProgressCallback",
    "ChunkStrategy",
    "ChunkingConfig",
    "LengthCurriculum",
    "LengthCurriculumStage",
    "SequenceChunk",
    "iter_fasta_chunks",
    "iter_sequence_chunks",
    "make_length_curriculum",
]
