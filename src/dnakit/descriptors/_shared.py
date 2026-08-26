"""Shared validation and sequence traversal for native descriptors."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias

from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError

SequenceInput: TypeAlias = DNA | DNASequence | DNARecord


class DescriptorAmbiguityPolicy(str, Enum):
    """Treatment of IUPAC symbols that are not one of A, C, G, and T."""

    ERROR = "error"
    IGNORE = "ignore"


def sequence_and_id(value: SequenceInput) -> tuple[DNASequence, str | None]:
    """Resolve the sequence and optional record identifier from public input."""

    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return sequence, None if record is None else record.id
    raise ConfigurationError(
        "Descriptor input must be a DNASequence or DNARecord.",
        context={"input_type": type(value).__name__},
    )


def coerce_ambiguity_policy(
    value: DescriptorAmbiguityPolicy | str,
) -> DescriptorAmbiguityPolicy:
    """Coerce a public ambiguity policy and report stable configuration errors."""

    try:
        return (
            value
            if isinstance(value, DescriptorAmbiguityPolicy)
            else DescriptorAmbiguityPolicy(value)
        )
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Unknown descriptor ambiguity policy.",
            context={"ambiguity_policy": value},
            hint="Choose 'error' or 'ignore'.",
        ) from exc


def validate_bool(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.", context={name: value})


def validate_positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            context={name: value},
        )


def reject_ambiguity(sequence: DNASequence, policy: DescriptorAmbiguityPolicy) -> None:
    """Reject ambiguity before any partial descriptor calculation."""

    if policy is DescriptorAmbiguityPolicy.ERROR and sequence.ambiguity_count:
        raise InvalidAlphabetError(
            "This descriptor does not accept ambiguous IUPAC symbols under the error policy.",
            code="DESCRIPTOR_AMBIGUITY_NOT_ALLOWED",
            context={"ambiguity_count": sequence.ambiguity_count},
            hint="Use ambiguity_policy='ignore' to omit affected symbols or windows.",
        )


def fragments(sequence: DNASequence, *, cross_gaps: bool) -> tuple[str, ...]:
    """Return symbol fragments without crossing a gap that explicitly forbids it."""

    validate_bool(cross_gaps, "cross_gaps")
    if not cross_gaps:
        return tuple(part for part in sequence.parts if isinstance(part, str))
    resolved: list[str] = []
    current: list[str] = []
    for part in sequence.parts:
        if isinstance(part, str):
            current.append(part)
        elif part.crossable is False and current:
            resolved.append("".join(current))
            current = []
    if current:
        resolved.append("".join(current))
    return tuple(resolved)


def canonical_runs(text: str) -> Iterator[str]:
    """Yield uninterrupted A/C/G/T runs, splitting at ignored IUPAC symbols."""

    start = 0
    for index, symbol in enumerate(text):
        if symbol not in "ACGT":
            if start < index:
                yield text[start:index]
            start = index + 1
    if start < len(text):
        yield text[start:]


def iter_kmers(
    sequence: DNASequence,
    *,
    k: int,
    overlapping: bool,
    cross_gaps: bool,
) -> Iterator[str]:
    """Yield valid canonical k-mers without joining ignored ambiguity."""

    validate_positive_int(k, "k")
    validate_bool(overlapping, "overlapping")
    step = 1 if overlapping else k
    for fragment in fragments(sequence, cross_gaps=cross_gaps):
        for run in canonical_runs(fragment):
            for start in range(0, len(run) - k + 1, step):
                yield run[start : start + k]


@dataclass(frozen=True)
class SymbolView:
    """Flattened symbols with coordinate and gap-segment information."""

    text: str
    coordinate_positions: tuple[int | None, ...]
    segment_ids: tuple[int, ...]
    crossable_segment_ids: tuple[int, ...]
    unknown_gap_boundaries: frozenset[int]


def symbol_view(sequence: DNASequence) -> SymbolView:
    """Build a symbol-indexed view without silently erasing gap boundaries."""

    symbols: list[str] = []
    coordinate_positions: list[int | None] = []
    segment_ids: list[int] = []
    crossable_segment_ids: list[int] = []
    unknown_gap_boundaries: set[int] = set()
    coordinate: int | None = 0
    segment = 0
    crossable_segment = 0
    for part in sequence.parts:
        if isinstance(part, str):
            for symbol in part:
                symbols.append(symbol)
                coordinate_positions.append(coordinate)
                segment_ids.append(segment)
                crossable_segment_ids.append(crossable_segment)
                if coordinate is not None:
                    coordinate += 1
        else:
            segment += 1
            if part.crossable is False:
                crossable_segment += 1
            if part.length is None:
                unknown_gap_boundaries.add(segment)
                coordinate = None
            elif coordinate is not None:
                coordinate += part.length
    return SymbolView(
        text="".join(symbols),
        coordinate_positions=tuple(coordinate_positions),
        segment_ids=tuple(segment_ids),
        crossable_segment_ids=tuple(crossable_segment_ids),
        unknown_gap_boundaries=frozenset(unknown_gap_boundaries),
    )


__all__ = ["DescriptorAmbiguityPolicy", "SequenceInput"]
