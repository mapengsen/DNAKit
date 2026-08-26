"""Bounded, explicit sequence-complexity and exact-repeat descriptors."""

from __future__ import annotations

import math
from collections import Counter

from dnakit.core import Gap, Topology
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

from ._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    canonical_runs,
    coerce_ambiguity_policy,
    fragments,
    reject_ambiguity,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from .results import ComplexityResult, ExactRepeatResult, RepeatRun


def linguistic_complexity(
    value: SequenceInput,
    *,
    max_word_size: int = 6,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
    max_observations: int = 10_000_000,
) -> ComplexityResult:
    """Compute observed/possible vocabulary ratios for k=1..max_word_size.

    For each k, the denominator is ``min(4**k, number of k-mer positions)``.
    The aggregate score is the product of per-k ratios, a common explicit
    linguistic-complexity definition; no DUST implementation is claimed.
    """

    sequence, sequence_id = sequence_and_id(value)
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Linguistic complexity currently requires a linear sequence.",
            code="COMPLEXITY_CIRCULAR_UNSUPPORTED",
        )
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    validate_positive_int(max_word_size, "max_word_size")
    validate_positive_int(max_observations, "max_observations")
    validate_bool(cross_gaps, "cross_gaps")
    if max_word_size > 16:
        raise ConfigurationError(
            "max_word_size cannot exceed 16.", code="COMPLEXITY_WORD_SIZE_LIMIT"
        )
    by_k: dict[str, float] = {}
    observed_by_k: dict[str, int] = {}
    possible_by_k: dict[str, int] = {}
    total_observations = 0
    for k in range(1, max_word_size + 1):
        words: set[str] = set()
        positions = 0
        for fragment in fragments(sequence, cross_gaps=cross_gaps):
            for run in canonical_runs(fragment):
                run_positions = max(0, len(run) - k + 1)
                positions += run_positions
                total_observations += run_positions
                if total_observations > max_observations:
                    raise ConfigurationError(
                        "Complexity calculation exceeds max_observations.",
                        code="COMPLEXITY_OBSERVATION_LIMIT",
                    )
                words.update(run[offset : offset + k] for offset in range(run_positions))
        possible = min(4**k, positions)
        ratio = len(words) / possible if possible else 0.0
        by_k[str(k)] = ratio
        observed_by_k[str(k)] = len(words)
        possible_by_k[str(k)] = possible
    nonempty = tuple(value for key, value in by_k.items() if possible_by_k[key] > 0)
    score = math.prod(nonempty) if nonempty else 0.0
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return ComplexityResult(
        name="linguistic_complexity",
        method="vocabulary-observed-over-possible-product",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        score=score,
        max_word_size=max_word_size,
        by_k=FrozenDict(by_k),
        observed_by_k=FrozenDict(observed_by_k),
        possible_by_k=FrozenDict(possible_by_k),
        observation_count=total_observations,
        max_observations=max_observations,
        formula="product_k(unique_kmers/min(4**k,valid_kmer_positions))",
    )


def exact_repeat_fraction(
    value: SequenceInput,
    *,
    min_unit_length: int = 1,
    max_unit_length: int = 20,
    min_repeats: int = 2,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
    max_comparisons: int = 5_000_000,
) -> ExactRepeatResult:
    """Measure union coverage of maximal exact tandem-repeat runs."""

    sequence, sequence_id = sequence_and_id(value)
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Exact tandem-repeat coverage currently requires a linear sequence.",
            code="REPEAT_CIRCULAR_UNSUPPORTED",
        )
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    for name, item in (
        ("min_unit_length", min_unit_length),
        ("max_unit_length", max_unit_length),
        ("min_repeats", min_repeats),
        ("max_comparisons", max_comparisons),
    ):
        validate_positive_int(item, name)
    validate_bool(cross_gaps, "cross_gaps")
    if min_unit_length > max_unit_length or max_unit_length > 100 or min_repeats < 2:
        raise ConfigurationError(
            "Repeat bounds require 1 <= min_unit_length <= max_unit_length <= 100 and "
            "min_repeats >= 2.",
            code="INVALID_REPEAT_DESCRIPTOR_CONFIG",
        )
    runs: list[RepeatRun] = []
    covered_positions: set[int] = set()
    comparisons = 0
    symbol_offset = 0
    for fragment in fragments(sequence, cross_gaps=cross_gaps):
        start = 0
        while start < len(fragment):
            found: tuple[str, int, int] | None = None
            for unit_length in range(min_unit_length, max_unit_length + 1):
                if start + unit_length * min_repeats > len(fragment):
                    continue
                unit = fragment[start : start + unit_length]
                if any(symbol not in "ACGT" for symbol in unit):
                    continue
                count = 1
                cursor = start + unit_length
                while cursor + unit_length <= len(fragment):
                    comparisons += 1
                    if comparisons > max_comparisons:
                        raise ConfigurationError(
                            "Repeat calculation exceeds max_comparisons.",
                            code="REPEAT_COMPARISON_LIMIT",
                        )
                    if fragment[cursor : cursor + unit_length] != unit:
                        break
                    count += 1
                    cursor += unit_length
                if count >= min_repeats:
                    found = unit, count, cursor
                    break
            if found is None:
                start += 1
                continue
            unit, count, end = found
            absolute_start = symbol_offset + start
            absolute_end = symbol_offset + end
            runs.append(RepeatRun(unit, len(unit), count, absolute_start, absolute_end))
            covered_positions.update(range(absolute_start, absolute_end))
            start = end
        symbol_offset += len(fragment)
    denominator = sequence.canonical_base_count
    repeat_count_by_unit = Counter(run.unit_length for run in runs)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return ExactRepeatResult(
        name="exact_repeat_fraction",
        method="maximal-exact-tandem-repeat-union",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        repeat_fraction=len(covered_positions) / denominator if denominator else 0.0,
        repeated_base_count=len(covered_positions),
        denominator=denominator,
        min_unit_length=min_unit_length,
        max_unit_length=max_unit_length,
        min_repeats=min_repeats,
        runs=tuple(runs),
        repeat_count_by_unit=FrozenDict(
            {str(key): value for key, value in sorted(repeat_count_by_unit.items())}
        ),
        comparisons=comparisons,
        max_comparisons=max_comparisons,
    )


__all__ = ["exact_repeat_fraction", "linguistic_complexity"]
