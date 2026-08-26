"""Shared input and policy helpers for native similarity functions."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import TypeAlias

from dnakit.core import DNA, DNARecord, DNASequence, DNASet, Gap
from dnakit.core.facade import resolve_single_dna
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

SequenceInput: TypeAlias = DNA | DNASequence | DNARecord
NumericVector: TypeAlias = Sequence[int | float]
NamedVector: TypeAlias = Mapping[str, int | float]


def sequence_and_id(value: SequenceInput, *, role: str) -> tuple[DNASequence, str | None]:
    """Resolve a sequence input and reject untyped values consistently."""

    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        sequence_id = None if record is None else record.id
    else:
        raise ConfigurationError(
            f"{role} must be DNASequence or DNARecord.",
            code="INVALID_SIMILARITY_SEQUENCE",
            context={"role": role, "input_type": type(value).__name__},
        )
    reject_gaps(sequence, role=role)
    return sequence, sequence_id


def reject_gaps(sequence: DNASequence, *, role: str) -> None:
    """Reject every explicit Gap instead of silently omitting it."""

    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    if gaps:
        raise UnsupportedGapOperationError(
            "Native MVP similarity does not accept explicit Gap objects.",
            code="SIMILARITY_GAP_NOT_ALLOWED",
            context={
                "role": role,
                "gap_count": len(gaps),
                "unknown_gap_count": sum(gap.length is None for gap in gaps),
            },
            hint="Compare gap-free fragments explicitly; gaps are never silently omitted.",
        )


def materialize_targets(
    targets: object,
    *,
    max_targets: int,
) -> tuple[SequenceInput, ...]:
    """Materialize at most ``max_targets`` targets without exhausting iterators."""

    resolved: tuple[SequenceInput, ...]
    if isinstance(targets, DNA):
        resolved = targets.records[: max_targets + 1]
    elif isinstance(targets, (DNASequence, DNARecord)):
        resolved = (targets,)
    elif isinstance(targets, DNASet):
        resolved = targets.records[: max_targets + 1]
    elif isinstance(targets, (str, bytes)):
        raise ConfigurationError(
            "Search targets must be DNASequence/DNARecord objects, not raw text.",
            code="INVALID_SEARCH_TARGETS",
        )
    elif not isinstance(targets, Iterable):
        raise ConfigurationError(
            "Search targets must be one sequence, DNASet, or iterable of sequences.",
            code="INVALID_SEARCH_TARGETS",
        )
    else:
        collected: list[SequenceInput] = []
        for index, item in enumerate(targets):
            if not isinstance(item, (DNA, DNASequence, DNARecord)):
                raise ConfigurationError(
                    "Every search target must be DNASequence or DNARecord.",
                    code="INVALID_SEARCH_TARGET",
                    context={"target_index": index, "type": type(item).__name__},
                )
            if index >= max_targets:
                raise ConfigurationError(
                    "Search target input exceeds max_targets.",
                    code="SEARCH_TARGET_LIMIT_EXCEEDED",
                    context={
                        "target_count": max_targets + 1,
                        "target_count_is_lower_bound": True,
                        "max_targets": max_targets,
                    },
                    hint="Reduce the target input or increase max_targets explicitly.",
                )
            if isinstance(item, DNA):
                if not item.is_single:
                    raise ConfigurationError(
                        "A nested DNA search target must contain exactly one record.",
                        code="INVALID_SEARCH_TARGET",
                        context={"target_index": index, "record_count": len(item)},
                    )
                collected.append(item.record)
            else:
                collected.append(item)
        resolved = tuple(collected)
    if len(resolved) > max_targets:
        raise ConfigurationError(
            "Search target input exceeds max_targets.",
            code="SEARCH_TARGET_LIMIT_EXCEEDED",
            context={
                "target_count": max_targets + 1,
                "target_count_is_lower_bound": True,
                "max_targets": max_targets,
            },
            hint="Reduce the target input or increase max_targets explicitly.",
        )
    if any(not isinstance(item, (DNASequence, DNARecord)) for item in resolved):
        raise ConfigurationError(
            "Every search target must be DNASequence or DNARecord.",
            code="INVALID_SEARCH_TARGET",
        )
    return resolved


def validate_bool(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be boolean.", context={name: value})


def validate_positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            context={name: value},
        )


def validate_non_negative_number(value: int | float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
    ):
        raise ConfigurationError(
            f"{name} must be a finite non-negative number.",
            context={name: value},
        )
    return float(value)


__all__ = ["NamedVector", "NumericVector", "SequenceInput"]
