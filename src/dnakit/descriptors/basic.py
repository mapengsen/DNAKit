"""Native length, composition, skew, and CpG descriptors."""

from __future__ import annotations

from collections import Counter
from typing import cast

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    canonical_runs,
    coerce_ambiguity_policy,
    fragments,
    reject_ambiguity,
    sequence_and_id,
)
from dnakit.descriptors.results import (
    CompositionResult,
    ContentResult,
    CpGResult,
    LengthResult,
    SkewResult,
)


def length_features(value: SequenceInput) -> LengthResult:
    """Return lengths without guessing the span of an unknown-length gap."""

    sequence, sequence_id = sequence_and_id(value)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return LengthResult(
        name="length",
        method="explicit_sequence_parts",
        sequence_id=sequence_id,
        ambiguity_policy=None,
        cross_gaps=False,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        symbol_length=sequence.symbol_length,
        coordinate_span=sequence.coordinate_span,
        canonical_base_count=sequence.canonical_base_count,
        ambiguity_length=sequence.ambiguity_count,
        known_gap_length=sum(gap.length or 0 for gap in gaps),
    )


def base_composition(
    value: SequenceInput,
    *,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
) -> CompositionResult:
    """Count A/C/G/T; ignored ambiguity is excluded from the fraction denominator."""

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    counts = Counter(symbol for symbol in sequence.symbols if symbol in "ACGT")
    denominator = sum(counts.values())
    ordered_counts = {base: counts[base] for base in "ACGT"}
    fractions = {
        base: (ordered_counts[base] / denominator if denominator else 0.0) for base in "ACGT"
    }
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return CompositionResult(
        name="base_composition",
        method="canonical_base_count",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=False,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        counts=FrozenDict(ordered_counts),
        fractions=FrozenDict(fractions),
        denominator=denominator,
        ignored_ambiguity_count=(
            sequence.ambiguity_count if policy is DescriptorAmbiguityPolicy.IGNORE else 0
        ),
    )


def gc_at_content(
    value: SequenceInput,
    *,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
) -> ContentResult:
    """Calculate GC and AT fractions over the A/C/G/T denominator."""

    composition = base_composition(value, ambiguity_policy=ambiguity_policy)
    gc_count = cast(int, composition.counts["G"]) + cast(int, composition.counts["C"])
    at_count = cast(int, composition.counts["A"]) + cast(int, composition.counts["T"])
    denominator = composition.denominator
    return ContentResult(
        name="gc_at_content",
        method="canonical_base_fraction",
        sequence_id=composition.sequence_id,
        ambiguity_policy=composition.ambiguity_policy,
        cross_gaps=False,
        gap_count=composition.gap_count,
        unknown_gap_count=composition.unknown_gap_count,
        gc_count=gc_count,
        at_count=at_count,
        gc_fraction=gc_count / denominator if denominator else None,
        at_fraction=at_count / denominator if denominator else None,
        denominator=denominator,
        ignored_ambiguity_count=composition.ignored_ambiguity_count,
    )


def base_skew(
    value: SequenceInput,
    *,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
) -> SkewResult:
    """Calculate GC and AT skew; a zero formula denominator produces None."""

    composition = base_composition(value, ambiguity_policy=ambiguity_policy)
    a = cast(int, composition.counts["A"])
    c = cast(int, composition.counts["C"])
    g = cast(int, composition.counts["G"])
    t = cast(int, composition.counts["T"])
    gc_denominator = g + c
    at_denominator = a + t
    return SkewResult(
        name="base_skew",
        method="gc=(g-c)/(g+c);at=(a-t)/(a+t)",
        sequence_id=composition.sequence_id,
        ambiguity_policy=composition.ambiguity_policy,
        cross_gaps=False,
        gap_count=composition.gap_count,
        unknown_gap_count=composition.unknown_gap_count,
        gc_skew=(g - c) / gc_denominator if gc_denominator else None,
        at_skew=(a - t) / at_denominator if at_denominator else None,
        gc_denominator=gc_denominator,
        at_denominator=at_denominator,
        ignored_ambiguity_count=composition.ignored_ambiguity_count,
    )


def cpg_features(
    value: SequenceInput,
    *,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> CpGResult:
    """Calculate CpG statistics without crossing gaps or ambiguity by default.

    Density is ``count(CG) / eligible adjacent canonical pairs``. Observed/expected
    is ``count(CG) * canonical_base_count / (count(C) * count(G))``. Undefined
    ratios are returned as ``None`` rather than NaN or infinity.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    text_runs = tuple(
        run
        for fragment in fragments(sequence, cross_gaps=cross_gaps)
        for run in canonical_runs(fragment)
    )
    cpg_count = sum(run.count("CG") for run in text_runs)
    pair_denominator = sum(max(0, len(run) - 1) for run in text_runs)
    canonical_length = sum(len(run) for run in text_runs)
    c_count = sum(run.count("C") for run in text_runs)
    g_count = sum(run.count("G") for run in text_runs)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return CpGResult(
        name="cpg",
        method="adjacent_cg_and_length_normalized_oe",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        cpg_count=cpg_count,
        density=cpg_count / pair_denominator if pair_denominator else None,
        observed_expected=(
            cpg_count * canonical_length / (c_count * g_count) if c_count and g_count else None
        ),
        adjacent_pair_denominator=pair_denominator,
        expected_length_denominator=canonical_length,
        c_count=c_count,
        g_count=g_count,
        ignored_ambiguity_count=(
            sequence.ambiguity_count if policy is DescriptorAmbiguityPolicy.IGNORE else 0
        ),
        density_formula="count(CG)/eligible_adjacent_canonical_pairs",
        observed_expected_formula="count(CG)*canonical_base_count/(count(C)*count(G))",
    )


__all__ = ["base_composition", "base_skew", "cpg_features", "gc_at_content", "length_features"]
