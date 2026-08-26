"""Shared bounded traversal and coordinate helpers for pattern analyses."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TypeAlias

from dnakit.core import (
    DNA,
    CompoundLocation,
    DNARecord,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    Interval,
    Issue,
    IssueSeverity,
    Location,
    OriginClass,
    Provenance,
    ReferenceInfo,
    Strand,
    UnresolvedLocation,
)
from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.facade import resolve_single_dna
from dnakit.core.gap import Gap
from dnakit.exceptions import ConfigurationError, SequenceError
from dnakit.patterns.results import PatternResult

SequenceInput: TypeAlias = DNA | DNASequence | DNARecord

IUPAC_BASES: Mapping[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


@dataclass(frozen=True)
class Segment:
    text: str
    symbol_start: int
    coordinate_start: int | None


def resolve_sequence(value: SequenceInput) -> tuple[DNASequence, str | None]:
    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return sequence, None if record is None else record.id
    raise ConfigurationError(
        "Pattern input must be a DNASequence or DNARecord.",
        context={"input_type": type(value).__name__},
    )


def validate_positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive integer.", context={name: value})


def validate_non_negative_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer.", context={name: value})


def validate_bool(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.", context={name: value})


def segments(sequence: DNASequence) -> tuple[Segment, ...]:
    """Split at every explicit Gap and retain all resolvable coordinates."""

    resolved: list[Segment] = []
    symbol_start = 0
    coordinate: int | None = 0
    for part in sequence.parts:
        if isinstance(part, str):
            if part:
                resolved.append(Segment(part, symbol_start, coordinate))
                symbol_start += len(part)
                if coordinate is not None:
                    coordinate += len(part)
        else:
            if part.length is None:
                coordinate = None
            elif coordinate is not None:
                coordinate += part.length
    return tuple(resolved)


def _split_location(start: int, length: int, total: int) -> Location:
    if length < 0 or length > total:
        raise ConfigurationError(
            "A circular match length must be between zero and sequence length."
        )
    end = start + length
    if end <= total:
        return Interval(start, end)
    return CompoundLocation((Interval(start, total), Interval(0, end - total)))


def segment_location(
    segment: Segment,
    local_start: int,
    local_end: int,
) -> tuple[Location, Location]:
    """Map one non-wrapping segment-local interval to symbol and genomic coordinates."""

    symbol_location = Interval(
        segment.symbol_start + local_start,
        segment.symbol_start + local_end,
    )
    if segment.coordinate_start is None:
        coordinate_location: Location = UnresolvedLocation(
            "unknown-length gap precedes this pattern hit",
            (symbol_location,),
        )
    else:
        coordinate_location = Interval(
            segment.coordinate_start + local_start,
            segment.coordinate_start + local_end,
        )
    return symbol_location, coordinate_location


def circular_location(start: int, length: int, total: int) -> tuple[Location, Location, bool]:
    start %= total
    location = _split_location(start, length, total)
    return location, location, isinstance(location, CompoundLocation)


def reverse_complement_text(text: str) -> str:
    return DNASequence(text, alphabet="iupac").reverse_complement().symbols


def iupac_compatible(left: str, right: str) -> bool:
    return bool(IUPAC_BASES[left] & IUPAC_BASES[right])


def validate_iupac_text(value: str, name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    resolved = value.upper()
    if not resolved and not allow_empty:
        raise ConfigurationError(f"{name} must be non-empty.")
    invalid = sorted(set(resolved) - set(IUPAC_BASES))
    if invalid:
        raise ConfigurationError(
            f"{name} contains symbols outside the DNA IUPAC alphabet.",
            context={"invalid_symbols": invalid},
        )
    return resolved


def require_ungapped_circular(sequence: DNASequence, operation: str) -> None:
    if sequence.topology.value == "circular" and sequence.is_gapped:
        raise SequenceError(
            f"{operation} cannot resolve a circular origin across explicit gaps.",
            code="GAPPED_CIRCULAR_PATTERN_UNSUPPORTED",
        )


def require_linear(sequence: DNASequence, operation: str) -> None:
    if sequence.topology.value == "circular":
        raise SequenceError(
            f"{operation} currently requires linear topology.",
            code="CIRCULAR_PATTERN_UNSUPPORTED",
            hint="Linearize at a documented origin before running this analysis.",
        )


def pattern_provenance(
    *,
    reimplementation: bool,
    reference_name: str | None = None,
    reference_version: str | None = None,
) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=(
                ImplementationLabel.REIMPLEMENTATION
                if reimplementation
                else ImplementationLabel.NATIVE
            ),
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=(
                OriginClass.PUBLISHED_ALGORITHM if reimplementation else OriginClass.DNAKIT
            ),
        ),
        reference=(
            None
            if reference_name is None
            else ReferenceInfo(reference_name, version=reference_version)
        ),
    )


def build_result(
    sequence: DNASequence,
    sequence_id: str | None,
    *,
    name: str,
    method: str,
    algorithm_version: str,
    parameters: Mapping[str, object],
    hits: Iterable[object],
    max_matches: int,
    truncated: bool,
    provenance: Provenance,
    issues: Iterable[Issue] = (),
) -> PatternResult[object]:
    gap_count = sum(isinstance(part, Gap) for part in sequence.parts)
    unknown_gap_count = sum(
        isinstance(part, Gap) and part.length is None for part in sequence.parts
    )
    resolved_issues = list(issues)
    if truncated:
        resolved_issues.append(
            Issue(
                "PATTERN_MATCH_LIMIT_REACHED",
                IssueSeverity.WARNING,
                "Pattern hits were truncated at max_matches.",
                details={"max_matches": max_matches},
            )
        )
    return PatternResult(
        name=name,
        method=method,
        algorithm_version=algorithm_version,
        sequence_id=sequence_id,
        parameters=freeze_mapping(parameters),
        hits=tuple(hits),
        inspected_symbol_count=sequence.symbol_length,
        gap_count=gap_count,
        unknown_gap_count=unknown_gap_count,
        max_matches=max_matches,
        truncated=truncated,
        coordinate_system="0-based-half-open",
        gap_policy="split-no-crossing",
        topology=sequence.topology.value,
        provenance=provenance,
        issues=tuple(resolved_issues),
    )


def frozen(values: Mapping[str, object]) -> FrozenDict:
    return freeze_mapping(values)


def coerce_strands(value: Strand | str) -> tuple[Strand, ...]:
    try:
        strand = value if isinstance(value, Strand) else Strand(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Unknown strand selection.",
            context={"strand": value},
            hint="Choose 'forward', 'reverse', or 'both'.",
        ) from exc
    if strand is Strand.BOTH:
        return (Strand.FORWARD, Strand.REVERSE)
    if strand not in (Strand.FORWARD, Strand.REVERSE):
        raise ConfigurationError("Pattern strand must be forward, reverse, or both.")
    return (strand,)


__all__ = ["IUPAC_BASES", "SequenceInput"]
