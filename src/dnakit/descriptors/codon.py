"""Native standard-code codon statistics."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.core.sequence import DNASequence
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    coerce_ambiguity_policy,
    reject_ambiguity,
    sequence_and_id,
    symbol_view,
    validate_bool,
)
from dnakit.descriptors.results import CodonResult
from dnakit.exceptions import ConfigurationError

_STANDARD_STARTS = frozenset({"ATG"})
_STANDARD_STOPS = frozenset({"TAA", "TAG", "TGA"})


@dataclass(frozen=True)
class _CodonAudit:
    counts: Counter[str]
    ignored_ambiguity: int
    gap_interrupted: int
    incomplete_bases: int
    phase_unresolved: bool
    unresolved_downstream_bases: int


def _coordinate_codon_audit(sequence: DNASequence, frame: int) -> _CodonAudit:
    slots: dict[int, dict[int, str]] = {}
    known_gaps: list[tuple[int, int]] = []
    coordinate: int | None = 0
    unknown_gap_start: int | None = None
    unresolved_downstream = 0
    for part in sequence.parts:
        if isinstance(part, str):
            if coordinate is None:
                unresolved_downstream += len(part)
                continue
            for symbol in part:
                if coordinate >= frame:
                    relative = coordinate - frame
                    slots.setdefault(relative // 3, {})[relative % 3] = symbol
                coordinate += 1
        elif coordinate is not None:
            if part.length is None:
                unknown_gap_start = coordinate
                coordinate = None
            else:
                known_gaps.append((coordinate, coordinate + part.length))
                coordinate += part.length

    counts: Counter[str] = Counter()
    ignored_ambiguity = 0
    interrupted = 0
    incomplete = 0
    for slot, symbols in slots.items():
        codon_start = frame + slot * 3
        codon_end = codon_start + 3
        if len(symbols) == 3:
            codon = "".join(symbols[position] for position in range(3))
            if any(base not in "ACGT" for base in codon):
                ignored_ambiguity += 1
            else:
                counts[codon] += 1
            continue
        overlaps_known_gap = any(
            gap_start < codon_end and gap_end > codon_start for gap_start, gap_end in known_gaps
        )
        overlaps_unknown_gap = (
            unknown_gap_start is not None and codon_start <= unknown_gap_start < codon_end
        )
        if overlaps_known_gap or overlaps_unknown_gap:
            interrupted += 1
        else:
            incomplete += len(symbols)
    return _CodonAudit(
        counts,
        ignored_ambiguity,
        interrupted,
        incomplete,
        unknown_gap_start is not None,
        unresolved_downstream,
    )


def _concatenated_codon_audit(sequence: DNASequence, frame: int) -> _CodonAudit:
    view = symbol_view(sequence)
    codon_segments = view.crossable_segment_ids
    counts: Counter[str] = Counter()
    ignored_ambiguity = 0
    gap_interrupted = 0
    complete_candidate_count = 0
    for start in range(frame, len(view.text) - 2, 3):
        complete_candidate_count += 1
        segments = codon_segments[start : start + 3]
        if len(set(segments)) != 1:
            gap_interrupted += 1
            continue
        codon = view.text[start : start + 3]
        if any(base not in "ACGT" for base in codon):
            ignored_ambiguity += 1
        else:
            counts[codon] += 1
    eligible_after_frame = max(0, len(view.text) - frame)
    return _CodonAudit(
        counts,
        ignored_ambiguity,
        gap_interrupted,
        eligible_after_frame - complete_candidate_count * 3,
        False,
        0,
    )


def codon_statistics(
    value: SequenceInput,
    *,
    frame: int = 0,
    genetic_code: int = 1,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> CodonResult:
    """Count in-frame codons using NCBI standard genetic code table 1.

    By default, ``frame`` is a zero-based offset in the full sequence coordinate
    span: known Gap lengths advance phase, codons overlapping a Gap are omitted,
    and an unknown-length Gap makes downstream phase unresolved. With
    ``cross_gaps=True``, crossable Gaps are explicitly omitted and phase follows
    concatenated symbol coordinates; a ``crossable=False`` Gap remains a hard
    boundary. An ignored ambiguous codon is omitted without changing phase.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    if isinstance(frame, bool) or not isinstance(frame, int) or frame not in {0, 1, 2}:
        raise ConfigurationError(
            "frame must be the zero-based offset 0, 1, or 2.",
            context={"frame": frame},
        )
    if genetic_code != 1 or isinstance(genetic_code, bool):
        raise ConfigurationError(
            "MVP codon statistics currently support only genetic_code=1.",
            context={"genetic_code": genetic_code},
        )
    validate_bool(cross_gaps, "cross_gaps")

    audit = (
        _concatenated_codon_audit(sequence, frame)
        if cross_gaps
        else _coordinate_codon_audit(sequence, frame)
    )

    ordered_counts = dict(sorted(audit.counts.items()))
    codon_count = sum(ordered_counts.values())
    frequencies = {codon: count / codon_count for codon, count in ordered_counts.items()}
    start_count = sum(audit.counts[codon] for codon in _STANDARD_STARTS)
    stop_count = sum(audit.counts[codon] for codon in _STANDARD_STOPS)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return CodonResult(
        name="codon",
        method="in_frame_standard_genetic_code",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        frame=frame,
        genetic_code=genetic_code,
        counts=FrozenDict(ordered_counts),
        frequencies=FrozenDict(frequencies),
        codon_count=codon_count,
        start_count=start_count,
        stop_count=stop_count,
        start_density=start_count / codon_count if codon_count else None,
        stop_density=stop_count / codon_count if codon_count else None,
        incomplete_base_count=audit.incomplete_bases,
        ignored_ambiguity_codon_count=audit.ignored_ambiguity,
        gap_interrupted_codon_count=audit.gap_interrupted,
        phase_coordinate_system=("concatenated-symbol" if cross_gaps else "sequence-coordinate"),
        phase_unresolved_after_gap=audit.phase_unresolved,
        unresolved_downstream_base_count=audit.unresolved_downstream_bases,
    )


__all__ = ["codon_statistics"]
