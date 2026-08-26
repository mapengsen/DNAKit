"""Native exact homopolymer descriptors."""

from __future__ import annotations

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    coerce_ambiguity_policy,
    reject_ambiguity,
    sequence_and_id,
    symbol_view,
    validate_bool,
    validate_positive_int,
)
from dnakit.descriptors.results import HomopolymerResult, HomopolymerRun


def homopolymer_runs(
    value: SequenceInput,
    *,
    min_run_length: int = 2,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> HomopolymerResult:
    """Find canonical-base runs in zero-based, half-open symbol coordinates.

    Ignored IUPAC symbols always terminate a run. Gaps also terminate a run
    unless ``cross_gaps=True``; any such crossing is recorded on the run.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    validate_positive_int(min_run_length, "min_run_length")
    validate_bool(cross_gaps, "cross_gaps")
    view = symbol_view(sequence)
    run_segments = view.crossable_segment_ids if cross_gaps else view.segment_ids
    found: list[HomopolymerRun] = []
    longest = {base: 0 for base in "ACGT"}
    run_start: int | None = None

    def close_run(end: int) -> None:
        nonlocal run_start
        if run_start is None:
            return
        base = view.text[run_start]
        length = end - run_start
        longest[base] = max(longest[base], length)
        if length >= min_run_length:
            first_segment = view.segment_ids[run_start]
            last_segment = view.segment_ids[end - 1]
            crossed_gap_count = last_segment - first_segment
            crossed_unknown = any(
                boundary in view.unknown_gap_boundaries
                for boundary in range(first_segment + 1, last_segment + 1)
            )
            coordinate_start = view.coordinate_positions[run_start]
            last_coordinate = view.coordinate_positions[end - 1]
            found.append(
                HomopolymerRun(
                    base=base,
                    length=length,
                    symbol_start=run_start,
                    symbol_end=end,
                    coordinate_start=coordinate_start,
                    coordinate_end=None if last_coordinate is None else last_coordinate + 1,
                    crossed_gap_count=crossed_gap_count,
                    crossed_unknown_gap=crossed_unknown,
                )
            )
        run_start = None

    for index, symbol in enumerate(view.text):
        same_segment = run_start is None or run_segments[index] == run_segments[index - 1]
        if symbol not in "ACGT":
            close_run(index)
        elif run_start is None:
            run_start = index
        elif symbol != view.text[run_start] or not same_segment:
            close_run(index)
            run_start = index
    close_run(len(view.text))

    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return HomopolymerResult(
        name="homopolymer",
        method="exact_canonical_run",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        min_run_length=min_run_length,
        longest_length=max(longest.values(), default=0),
        longest_by_base=FrozenDict(longest),
        runs=tuple(found),
        ignored_ambiguity_count=(
            sequence.ambiguity_count if policy is DescriptorAmbiguityPolicy.IGNORE else 0
        ),
    )


__all__ = ["homopolymer_runs"]
