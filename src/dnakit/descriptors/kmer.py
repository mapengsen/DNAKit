"""Native exact k-mer statistics."""

from __future__ import annotations

from collections import Counter

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    coerce_ambiguity_policy,
    iter_kmers,
    reject_ambiguity,
    sequence_and_id,
    validate_bool,
)
from dnakit.descriptors.results import KmerResult
from dnakit.exceptions import ConfigurationError

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical_kmer(value: str) -> str:
    """Return the lexicographically smaller strand representation."""

    if not isinstance(value, str) or not value or any(base not in "ACGT" for base in value):
        raise ConfigurationError(
            "A canonical k-mer must be a non-empty A/C/G/T string.",
            context={"kmer": value},
        )
    reverse_complement = value.translate(_COMPLEMENT)[::-1]
    return min(value, reverse_complement)


def kmer_statistics(
    value: SequenceInput,
    k: int,
    *,
    overlapping: bool = True,
    canonical: bool = False,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> KmerResult:
    """Return sparse count, frequency, and presence data for valid k-mers.

    Under ``ignore``, any window containing an IUPAC ambiguity is omitted and
    k-mers never bridge that symbol. If ``k`` exceeds every eligible run, the
    three sparse outputs are empty and the denominator is zero.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    validate_bool(canonical, "canonical")
    observed = iter_kmers(
        sequence,
        k=k,
        overlapping=overlapping,
        cross_gaps=cross_gaps,
    )
    counts = Counter(canonical_kmer(item) if canonical else item for item in observed)
    ordered_counts = dict(sorted(counts.items()))
    denominator = sum(ordered_counts.values())
    frequencies = {item: count / denominator for item, count in ordered_counts.items()}
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return KmerResult(
        name="kmer",
        method="exact_sliding_kmer",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        k=k,
        overlapping=overlapping,
        canonical=canonical,
        counts=FrozenDict(ordered_counts),
        frequencies=FrozenDict(frequencies),
        presence=tuple(ordered_counts),
        denominator=denominator,
        ignored_ambiguity_count=(
            sequence.ambiguity_count if policy is DescriptorAmbiguityPolicy.IGNORE else 0
        ),
    )


__all__ = ["canonical_kmer", "kmer_statistics"]
