"""Native Shannon entropy descriptors."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable

from dnakit.core.gap import Gap
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    coerce_ambiguity_policy,
    iter_kmers,
    reject_ambiguity,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from dnakit.descriptors.results import EntropyResult
from dnakit.exceptions import ConfigurationError


def shannon_entropy(
    value: SequenceInput,
    *,
    unit: str = "base",
    k: int = 1,
    log_base: float = 2.0,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> EntropyResult:
    """Calculate Shannon entropy for canonical bases or valid k-mers.

    An empty observation distribution has entropy ``0.0``. The log base must be
    finite, positive, and different from one.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    if unit not in {"base", "kmer"}:
        raise ConfigurationError(
            "Entropy unit must be 'base' or 'kmer'.",
            context={"unit": unit},
        )
    validate_positive_int(k, "k")
    validate_bool(cross_gaps, "cross_gaps")
    if (
        isinstance(log_base, bool)
        or not isinstance(log_base, (int, float))
        or not math.isfinite(log_base)
        or log_base <= 0
        or log_base == 1
    ):
        raise ConfigurationError(
            "log_base must be finite, positive, and different from one.",
            context={"log_base": log_base},
        )

    observations: Iterable[str]
    if unit == "base":
        observations = (symbol for symbol in sequence.symbols if symbol in "ACGT")
        resolved_k = 1
    else:
        observations = iter_kmers(
            sequence,
            k=k,
            overlapping=True,
            cross_gaps=cross_gaps,
        )
        resolved_k = k
    counts = Counter(observations)
    total = sum(counts.values())
    entropy = (
        -sum((count / total) * math.log(count / total, log_base) for count in counts.values())
        if total
        else 0.0
    )
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return EntropyResult(
        name="shannon_entropy",
        method="shannon",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        entropy=entropy,
        unit=unit,
        k=resolved_k,
        log_base=float(log_base),
        observation_count=total,
        category_count=len(counts),
        ignored_ambiguity_count=(
            sequence.ambiguity_count if policy is DescriptorAmbiguityPolicy.IGNORE else 0
        ),
    )


__all__ = ["shannon_entropy"]
