"""Internal helpers shared by deterministic sequence operations."""

from __future__ import annotations

from collections.abc import Iterable

from dnakit.core import DNAAlphabet, DNASequence, Gap, Interval, Strandedness, Topology
from dnakit.exceptions import (
    ConfigurationError,
    CoordinateError,
    SequenceError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)


def require_sequence(sequence: DNASequence) -> DNASequence:
    """Return a validated sequence argument with a domain-specific error."""

    if not isinstance(sequence, DNASequence):
        raise SequenceError(
            "A DNASequence object is required.",
            code="INVALID_SEQUENCE_ARGUMENT",
            context={"type": type(sequence).__name__},
        )
    return sequence


def resolved_span(sequence: DNASequence, *, operation: str) -> int:
    """Return a coordinate span, refusing to guess an unknown gap length."""

    span = sequence.coordinate_span
    if span is None:
        raise UnknownLengthError(
            f"{operation} requires every Gap to have a known length.",
            code="UNRESOLVED_OPERATION_COORDINATES",
            context={"operation": operation},
            hint="Resolve unknown Gap.length values before using coordinate-based operations.",
        )
    return span


def validate_position(position: int, *, span: int, allow_end: bool = True) -> None:
    """Validate one zero-based position or insertion boundary."""

    if isinstance(position, bool) or not isinstance(position, int):
        raise CoordinateError(
            "Position must be a non-negative integer.",
            context={"position": position, "sequence_span": span},
        )
    upper_ok = position <= span if allow_end else position < span
    if position < 0 or not upper_ok:
        relation = "at most" if allow_end else "smaller than"
        raise CoordinateError(
            f"Position must be a non-negative integer {relation} the sequence span.",
            context={"position": position, "sequence_span": span},
        )


def validate_interval(start: int, end: int, *, span: int) -> Interval:
    """Build and bounds-check an internal zero-based, half-open interval."""

    interval = Interval(start, end)
    if interval.end > span:
        raise CoordinateError(
            "Interval exceeds the sequence coordinate span.",
            context={"start": start, "end": end, "sequence_span": span},
        )
    return interval


def copy_gap(gap: Gap, *, length: int) -> Gap:
    """Copy a gap while replacing only its known coordinate length."""

    return Gap(
        length,
        kind=gap.kind,
        crossable=gap.crossable,
        evidence=gap.evidence,
        metadata=gap.metadata,
    )


def interval_overlaps_gap(sequence: DNASequence, interval: Interval) -> bool:
    """Return whether a resolved interval contains any gap coordinate."""

    cursor = 0
    for part in sequence.parts:
        part_length = len(part) if isinstance(part, str) else part.length
        if part_length is None:
            raise UnknownLengthError(
                "Gap overlap cannot be resolved for an unknown-length Gap.",
                code="UNRESOLVED_GAP_COORDINATES",
            )
        part_end = cursor + part_length
        if isinstance(part, Gap) and interval.start < part_end and cursor < interval.end:
            return True
        cursor = part_end
    return False


def position_inside_gap(sequence: DNASequence, position: int) -> bool:
    """Return whether an insertion boundary lies strictly inside a gap."""

    cursor = 0
    for part in sequence.parts:
        part_length = len(part) if isinstance(part, str) else part.length
        if part_length is None:
            raise UnknownLengthError(
                "Insertion coordinates cannot be resolved for an unknown-length Gap.",
                code="UNRESOLVED_GAP_COORDINATES",
            )
        part_end = cursor + part_length
        if isinstance(part, Gap) and cursor < position < part_end:
            return True
        cursor = part_end
    return False


def reject_gap_overlap(sequence: DNASequence, interval: Interval, *, operation: str) -> None:
    """Reject edits that would silently consume assembly-gap coordinates."""

    if interval_overlaps_gap(sequence, interval):
        raise UnsupportedGapOperationError(
            f"{operation} cannot overlap a Gap.",
            code="EDIT_OVERLAPS_GAP",
            context={"operation": operation, "start": interval.start, "end": interval.end},
            hint="Edit nucleotide fragments separately or resolve the assembly gap first.",
        )


def infer_alphabet(symbols: str, *, fallback: DNAAlphabet) -> DNAAlphabet:
    """Infer the narrowest DNA alphabet for already-normalized symbols."""

    if not symbols:
        return fallback
    return DNAAlphabet.STRICT if set(symbols) <= set("ACGT") else DNAAlphabet.IUPAC


def promote_alphabet(*alphabets: DNAAlphabet) -> DNAAlphabet:
    """Return IUPAC if any contributing sequence requires it."""

    return (
        DNAAlphabet.IUPAC
        if any(alphabet is DNAAlphabet.IUPAC for alphabet in alphabets)
        else DNAAlphabet.STRICT
    )


def coerce_fragment(
    fragment: DNASequence | str,
    *,
    fallback_alphabet: DNAAlphabet,
    strandedness: Strandedness,
) -> DNASequence:
    """Convert a normalized string or linear sequence to an embeddable fragment."""

    if isinstance(fragment, str):
        alphabet = infer_alphabet(fragment, fallback=fallback_alphabet)
        return DNASequence(fragment, alphabet=alphabet, strandedness=strandedness)
    if not isinstance(fragment, DNASequence):
        raise SequenceError(
            "A sequence fragment must be DNASequence or normalized DNA text.",
            code="INVALID_SEQUENCE_FRAGMENT",
            context={"type": type(fragment).__name__},
        )
    if fragment.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "A circular DNASequence cannot be embedded as a linear fragment.",
            code="CIRCULAR_FRAGMENT_NOT_SUPPORTED",
        )
    if fragment.strandedness is not strandedness:
        raise ConfigurationError(
            "All sequence fragments must have matching strandedness.",
            code="STRANDEDNESS_MISMATCH",
            context={
                "expected": strandedness.value,
                "observed": fragment.strandedness.value,
            },
        )
    return fragment


def combine_parts(
    parts: Iterable[str | Gap],
    *,
    alphabet: DNAAlphabet,
    topology: Topology,
    strandedness: Strandedness,
) -> DNASequence:
    """Construct a fresh immutable sequence from operation output parts."""

    return DNASequence(
        parts,
        alphabet=alphabet,
        topology=topology,
        strandedness=strandedness,
    )


__all__: list[str] = []
