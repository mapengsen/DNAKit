"""Coordinate-safe edits on immutable DNA sequences."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import pairwise
from typing import Literal, TypeAlias

from dnakit.core import DNASequence, Gap, Interval, Topology
from dnakit.exceptions import ConfigurationError, CoordinateError, UnsupportedGapOperationError
from dnakit.ops._common import (
    coerce_fragment,
    combine_parts,
    copy_gap,
    infer_alphabet,
    position_inside_gap,
    promote_alphabet,
    reject_gap_overlap,
    require_sequence,
    resolved_span,
    validate_interval,
    validate_position,
)

EditKind: TypeAlias = Literal["insert", "delete", "substitute", "mask"]


@dataclass(frozen=True)
class Edit:
    """One applied edit in internal zero-based, half-open coordinates."""

    kind: EditKind
    start: int
    end: int
    replacement_parts: tuple[str | Gap, ...] = ()
    removed_parts: tuple[str | Gap, ...] = ()

    @property
    def replacement_symbols(self) -> str:
        """Replacement nucleotide symbols; inspect ``replacement_parts`` for gaps."""

        return "".join(part for part in self.replacement_parts if isinstance(part, str))

    @property
    def removed_symbols(self) -> str:
        """Removed nucleotide symbols; inspect ``removed_parts`` for gaps."""

        return "".join(part for part in self.removed_parts if isinstance(part, str))


@dataclass(frozen=True)
class EditResult:
    """A new sequence and ordered atomic edits in original coordinates."""

    sequence: DNASequence
    edits: tuple[Edit, ...]


def _split_at(sequence: DNASequence, position: int) -> tuple[list[str | Gap], list[str | Gap]]:
    """Split resolved parts at a nucleotide or gap boundary."""

    left: list[str | Gap] = []
    right: list[str | Gap] = []
    cursor = 0
    on_right = False
    for part in sequence.parts:
        part_length = len(part) if isinstance(part, str) else part.length
        assert part_length is not None
        part_end = cursor + part_length
        if on_right:
            right.append(part)
        elif position == cursor:
            right.append(part)
            on_right = True
        elif position == part_end:
            left.append(part)
        elif cursor < position < part_end:
            offset = position - cursor
            if isinstance(part, Gap):
                raise UnsupportedGapOperationError(
                    "A sequence cannot be split inside a Gap.",
                    code="EDIT_INSIDE_GAP",
                    context={"position": position},
                )
            left.append(part[:offset])
            right.append(part[offset:])
            on_right = True
        else:
            left.append(part)
        cursor = part_end
    return left, right


def _remove_interval(
    sequence: DNASequence,
    start: int,
    end: int,
) -> tuple[list[str | Gap], tuple[str | Gap, ...]]:
    """Remove a gap-free coordinate interval and return retained/removed parts."""

    left, after_start = _split_at(sequence, start)
    tail_sequence = DNASequence(
        after_start,
        alphabet=sequence.alphabet,
        topology=Topology.LINEAR,
        strandedness=sequence.strandedness,
    )
    removed_length = end - start
    removed_parts, right = _split_at(tail_sequence, removed_length)
    return [*left, *right], tuple(removed_parts)


def _require_ungapped_fragment(fragment: DNASequence, *, operation: str) -> DNASequence:
    if fragment.is_gapped:
        raise UnsupportedGapOperationError(
            f"{operation} does not accept a gapped replacement fragment.",
            code="GAPPED_EDIT_FRAGMENT",
            context={"operation": operation},
            hint="Use concat(..., gap=...) to assemble fragments around an explicit Gap.",
        )
    return fragment


def insert(
    sequence: DNASequence,
    position: int,
    fragment: DNASequence | str,
) -> EditResult:
    """Insert a normalized fragment at a zero-based coordinate boundary."""

    source = require_sequence(sequence)
    span = resolved_span(source, operation="insert")
    validate_position(position, span=span)
    if position_inside_gap(source, position):
        raise UnsupportedGapOperationError(
            "Insertion inside a Gap is not supported.",
            code="EDIT_INSIDE_GAP",
            context={"position": position},
        )
    replacement = _require_ungapped_fragment(
        coerce_fragment(
            fragment,
            fallback_alphabet=source.alphabet,
            strandedness=source.strandedness,
        ),
        operation="insert",
    )
    left, right = _split_at(source, position)
    alphabet = promote_alphabet(source.alphabet, replacement.alphabet)
    result = combine_parts(
        (*left, *replacement.parts, *right),
        alphabet=alphabet,
        topology=source.topology,
        strandedness=source.strandedness,
    )
    return EditResult(
        result,
        (Edit("insert", position, position, replacement.parts),),
    )


def delete(sequence: DNASequence, start: int, end: int) -> EditResult:
    """Delete a gap-free zero-based, half-open interval."""

    source = require_sequence(sequence)
    span = resolved_span(source, operation="delete")
    interval = validate_interval(start, end, span=span)
    reject_gap_overlap(source, interval, operation="delete")
    new_parts, removed_parts = _remove_interval(source, start, end)
    result = combine_parts(
        new_parts,
        alphabet=source.alphabet,
        topology=source.topology,
        strandedness=source.strandedness,
    )
    return EditResult(result, (Edit("delete", start, end, removed_parts=removed_parts),))


def substitute(
    sequence: DNASequence,
    start: int,
    end: int,
    fragment: DNASequence | str,
) -> EditResult:
    """Replace a gap-free zero-based, half-open interval with a fragment."""

    source = require_sequence(sequence)
    span = resolved_span(source, operation="substitute")
    interval = validate_interval(start, end, span=span)
    reject_gap_overlap(source, interval, operation="substitute")
    replacement = _require_ungapped_fragment(
        coerce_fragment(
            fragment,
            fallback_alphabet=source.alphabet,
            strandedness=source.strandedness,
        ),
        operation="substitute",
    )
    retained_parts, removed_parts = _remove_interval(source, start, end)
    interim = combine_parts(
        retained_parts,
        alphabet=source.alphabet,
        topology=Topology.LINEAR,
        strandedness=source.strandedness,
    )
    inserted = insert(interim, start, replacement)
    alphabet = promote_alphabet(source.alphabet, replacement.alphabet)
    result = combine_parts(
        inserted.sequence.parts,
        alphabet=alphabet,
        topology=source.topology,
        strandedness=source.strandedness,
    )
    return EditResult(
        result,
        (Edit("substitute", start, end, replacement.parts, removed_parts),),
    )


def trim(sequence: DNASequence, *, left: int = 0, right: int = 0) -> EditResult:
    """Remove exact coordinate lengths from both ends.

    Known gaps may be shortened by trimming.  Unknown gaps are rejected because
    the requested boundaries cannot be resolved.
    """

    source = require_sequence(sequence)
    if source.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Trimming a circular sequence is not supported in the MVP.",
            code="CIRCULAR_TRIM_NOT_SUPPORTED",
            hint="Extract an explicit linear subsequence before trimming.",
        )
    span = resolved_span(source, operation="trim")
    for name, value in (("left", left), ("right", right)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ConfigurationError(
                f"{name} trim length must be a non-negative integer.",
                code="INVALID_TRIM_LENGTH",
                context={name: value},
            )
    if left + right > span:
        raise ConfigurationError(
            "Combined trim lengths exceed the sequence span.",
            code="TRIM_EXCEEDS_SEQUENCE",
            context={"left": left, "right": right, "sequence_span": span},
        )
    result = subsequence(source, left, span - right, allow_gaps=True)
    edits: list[Edit] = []
    if left:
        edits.append(Edit("delete", 0, left, removed_parts=tuple(_slice_parts(source, 0, left))))
    if right:
        edits.append(
            Edit(
                "delete",
                span - right,
                span,
                removed_parts=tuple(_slice_parts(source, span - right, span)),
            )
        )
    return EditResult(result, tuple(edits))


def _slice_parts(sequence: DNASequence, start: int, end: int) -> list[str | Gap]:
    selected: list[str | Gap] = []
    cursor = 0
    for part in sequence.parts:
        part_length = len(part) if isinstance(part, str) else part.length
        assert part_length is not None
        part_end = cursor + part_length
        overlap_start = max(start, cursor)
        overlap_end = min(end, part_end)
        if overlap_start < overlap_end:
            local_start = overlap_start - cursor
            local_end = overlap_end - cursor
            if isinstance(part, str):
                selected.append(part[local_start:local_end])
            else:
                selected.append(copy_gap(part, length=local_end - local_start))
        cursor = part_end
    return selected


def subsequence(
    sequence: DNASequence,
    start: int,
    end: int,
    *,
    allow_gaps: bool = False,
) -> DNASequence:
    """Extract a zero-based, half-open interval, optionally across an origin.

    On a circular sequence, ``start > end`` means ``[start, span) + [0, end)``.
    ``start == end`` remains an empty interval; request ``0..span`` for a full
    circle. Linear sequences reject reversed bounds. Unknown-length gaps are
    always rejected. Known gaps are included only when ``allow_gaps=True``;
    partial overlap creates a shorter Gap with retained metadata.
    """

    source = require_sequence(sequence)
    span = resolved_span(source, operation="subsequence")
    validate_position(start, span=span)
    validate_position(end, span=span)
    intervals: tuple[Interval, ...]
    if start <= end:
        intervals = (Interval(start, end),)
    elif source.topology is Topology.CIRCULAR:
        intervals = tuple(
            interval for interval in (Interval(start, span), Interval(0, end)) if len(interval) > 0
        )
    else:
        raise CoordinateError(
            "A linear subsequence cannot have end smaller than start.",
            code="LINEAR_SUBSEQUENCE_REVERSED",
            context={"start": start, "end": end, "sequence_span": span},
            hint="Use a circular DNASequence for origin-spanning extraction.",
        )
    if not allow_gaps:
        for interval in intervals:
            reject_gap_overlap(source, interval, operation="subsequence")
    parts = [
        part
        for interval in intervals
        for part in _slice_parts(source, interval.start, interval.end)
    ]
    return combine_parts(
        parts,
        alphabet=source.alphabet,
        topology=Topology.LINEAR,
        strandedness=source.strandedness,
    )


def circular_subsequence(
    sequence: DNASequence,
    start: int,
    end: int,
    *,
    allow_gaps: bool = False,
) -> DNASequence:
    """Extract from a circular sequence, including an optional origin wrap."""

    source = require_sequence(sequence)
    if source.topology is not Topology.CIRCULAR:
        raise ConfigurationError(
            "circular_subsequence() requires topology='circular'.",
            code="CIRCULAR_TOPOLOGY_REQUIRED",
        )
    return subsequence(source, start, end, allow_gaps=allow_gaps)


def mask(
    sequence: DNASequence,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
) -> EditResult:
    """Replace one or more gap-free intervals with a one-character IUPAC mask."""

    source = require_sequence(sequence)
    if not isinstance(symbol, str) or len(symbol) != 1 or symbol not in "ACGTRYSWKMBDHVN":
        raise ConfigurationError(
            "Mask symbol must be one uppercase DNA IUPAC character.",
            code="INVALID_MASK_SYMBOL",
            context={"symbol": symbol},
        )
    span = resolved_span(source, operation="mask")
    try:
        raw_intervals = tuple(intervals)
    except TypeError as exc:
        raise ConfigurationError(
            "intervals must be an iterable of (start, end) pairs.",
            code="INVALID_MASK_INTERVALS",
        ) from exc
    resolved_items: list[Interval] = []
    for index, item in enumerate(raw_intervals):
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ConfigurationError(
                "Each mask interval must be a (start, end) pair.",
                code="INVALID_MASK_INTERVAL",
                context={"index": index},
            )
        resolved_items.append(validate_interval(item[0], item[1], span=span))
    resolved = tuple(interval for interval in resolved_items if len(interval) > 0)
    ordered = tuple(sorted(resolved, key=lambda item: (item.start, item.end)))
    if not ordered:
        result = combine_parts(
            source.parts,
            alphabet=source.alphabet,
            topology=source.topology,
            strandedness=source.strandedness,
        )
        return EditResult(result, ())
    if any(current.start < previous.end for previous, current in pairwise(ordered)):
        raise ConfigurationError(
            "Mask intervals must not overlap.",
            code="OVERLAPPING_MASK_INTERVALS",
        )
    for interval in ordered:
        reject_gap_overlap(source, interval, operation="mask")

    current = source
    edits: list[Edit] = []
    for interval in reversed(ordered):
        substituted = substitute(
            current,
            interval.start,
            interval.end,
            symbol * len(interval),
        )
        atomic = substituted.edits[0]
        edits.append(
            Edit(
                "mask",
                atomic.start,
                atomic.end,
                atomic.replacement_parts,
                atomic.removed_parts,
            )
        )
        current = substituted.sequence
    alphabet = promote_alphabet(source.alphabet, infer_alphabet(symbol, fallback=source.alphabet))
    current = combine_parts(
        current.parts,
        alphabet=alphabet,
        topology=source.topology,
        strandedness=source.strandedness,
    )
    return EditResult(current, tuple(reversed(edits)))


__all__ = [
    "Edit",
    "EditResult",
    "circular_subsequence",
    "delete",
    "insert",
    "mask",
    "subsequence",
    "substitute",
    "trim",
]
