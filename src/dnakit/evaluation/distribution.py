"""Transparent descriptive distribution comparison between two DNA collections."""

from __future__ import annotations

import math
from collections import Counter

from dnakit.core import DNARecord
from dnakit.descriptors import exact_repeat_fraction, gc_at_content

from ._shared import EvaluationInput, materialize_input, mean, record_for, report, require_nonempty
from .config import DistributionEvaluationConfig
from .results import EvaluationReport


def _ecdf_distance(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    """Return the exact two-sample Kolmogorov-Smirnov D statistic."""

    if not left or not right:
        return 0.0
    support = sorted(set(left) | set(right))
    left_sorted = sorted(left)
    right_sorted = sorted(right)
    left_index = 0
    right_index = 0
    distance = 0.0
    for value in support:
        while left_index < len(left_sorted) and left_sorted[left_index] <= value:
            left_index += 1
        while right_index < len(right_sorted) and right_sorted[right_index] <= value:
            right_index += 1
        distance = max(
            distance,
            abs(left_index / len(left_sorted) - right_index / len(right_sorted)),
        )
    return distance


def _total_variation(left: Counter[str], right: Counter[str]) -> float:
    left_total = sum(left.values())
    right_total = sum(right.values())
    if left_total == 0 and right_total == 0:
        return 0.0
    support = set(left) | set(right)
    return 0.5 * math.fsum(
        abs(
            (left[key] / left_total if left_total else 0.0)
            - (right[key] / right_total if right_total else 0.0)
        )
        for key in support
    )


def _kmer_counts(records: tuple[DNARecord, ...], k: int, canonical: bool) -> Counter[str]:
    complement = str.maketrans("ACGT", "TGCA")
    counts: Counter[str] = Counter()
    for record in records:
        for part in record.sequence.parts:
            if not isinstance(part, str):
                continue
            for start in range(max(0, len(part) - k + 1)):
                word = part[start : start + k]
                if set(word) - set("ACGT"):
                    continue
                if canonical:
                    word = min(word, word.translate(complement)[::-1])
                counts[word] += 1
    return counts


def _kmer_observation_count(
    records: tuple[DNARecord, ...],
    k: int,
    *,
    initial_count: int,
    maximum: int,
) -> int:
    """Count eligible k-mers without allocating a frequency table.

    The scan stops at the first observation beyond ``maximum`` so an
    over-budget input cannot build either collection's ``Counter`` first.
    """

    observation_count = initial_count
    canonical_symbols = frozenset("ACGT")
    for record in records:
        for part in record.sequence.parts:
            if not isinstance(part, str):
                continue
            for start in range(max(0, len(part) - k + 1)):
                if set(part[start : start + k]) - canonical_symbols:
                    continue
                observation_count += 1
                if observation_count > maximum:
                    from dnakit.exceptions import ConfigurationError

                    raise ConfigurationError(
                        "Distribution comparison exceeds max_kmer_observations.",
                        code="DISTRIBUTION_KMER_LIMIT",
                        context={
                            "observation_count": observation_count,
                            "max_kmer_observations": maximum,
                        },
                    )
    return observation_count


def _motif_rates(records: tuple[DNARecord, ...], motifs: tuple[str, ...]) -> dict[str, float]:
    denominator = sum(record.sequence.symbol_length for record in records)
    return {
        motif: (
            sum(
                part.count(motif)
                for record in records
                for part in record.sequence.parts
                if isinstance(part, str)
            )
            / denominator
            if denominator
            else 0.0
        )
        for motif in motifs
    }


def _feature_distances(
    left: tuple[DNARecord, ...],
    right: tuple[DNARecord, ...],
    config: DistributionEvaluationConfig,
) -> tuple[dict[str, float], dict[str, object]]:
    distances: dict[str, float] = {}
    details: dict[str, object] = {}
    if "length" in config.features:
        left_values = tuple(float(record.sequence.symbol_length) for record in left)
        right_values = tuple(float(record.sequence.symbol_length) for record in right)
        distances["length"] = _ecdf_distance(left_values, right_values)
        details["length"] = {
            "distance": distances["length"],
            "method": "two-sample empirical-CDF supremum distance (KS D; no p-value)",
            "left_mean": mean(left_values),
            "right_mean": mean(right_values),
        }
    if "gc" in config.features:
        left_values = tuple(
            result.gc_fraction
            for record in left
            if (result := gc_at_content(record, ambiguity_policy="ignore")).gc_fraction is not None
        )
        right_values = tuple(
            result.gc_fraction
            for record in right
            if (result := gc_at_content(record, ambiguity_policy="ignore")).gc_fraction is not None
        )
        distances["gc"] = _ecdf_distance(left_values, right_values)
        details["gc"] = {
            "distance": distances["gc"],
            "method": "two-sample empirical-CDF supremum distance (KS D; no p-value)",
            "left_mean": mean(left_values),
            "right_mean": mean(right_values),
            "undefined_gc_omitted": (len(left) - len(left_values))
            + (len(right) - len(right_values)),
        }
    if "kmer" in config.features:
        left_observations = _kmer_observation_count(
            left,
            config.k,
            initial_count=0,
            maximum=config.max_kmer_observations,
        )
        observation_count = _kmer_observation_count(
            right,
            config.k,
            initial_count=left_observations,
            maximum=config.max_kmer_observations,
        )
        left_counts = _kmer_counts(left, config.k, config.canonical)
        right_counts = _kmer_counts(right, config.k, config.canonical)
        distances["kmer"] = _total_variation(left_counts, right_counts)
        details["kmer"] = {
            "distance": distances["kmer"],
            "method": "total variation over pooled exact k-mer frequency vectors",
            "left_observations": left_observations,
            "right_observations": observation_count - left_observations,
        }
    if "motif" in config.features:
        left_rates = _motif_rates(left, config.motifs)
        right_rates = _motif_rates(right, config.motifs)
        per_motif = {motif: abs(left_rates[motif] - right_rates[motif]) for motif in config.motifs}
        distances["motif"] = math.fsum(per_motif.values()) / len(per_motif)
        details["motif"] = {
            "distance": distances["motif"],
            "method": "mean absolute difference of non-overlapping motif counts per symbol",
            "left_rates": left_rates,
            "right_rates": right_rates,
            "absolute_differences": per_motif,
        }
    if "repeat" in config.features:
        left_values = tuple(
            exact_repeat_fraction(
                record,
                ambiguity_policy="ignore",
                max_comparisons=config.max_repeat_comparisons_per_sequence,
            ).repeat_fraction
            for record in left
        )
        right_values = tuple(
            exact_repeat_fraction(
                record,
                ambiguity_policy="ignore",
                max_comparisons=config.max_repeat_comparisons_per_sequence,
            ).repeat_fraction
            for record in right
        )
        distances["repeat"] = _ecdf_distance(left_values, right_values)
        details["repeat"] = {
            "distance": distances["repeat"],
            "method": "two-sample empirical-CDF supremum distance of exact-repeat union fraction",
            "left_mean": mean(left_values),
            "right_mean": mean(right_values),
        }
    return distances, details


def evaluate_distribution_similarity(
    left: EvaluationInput,
    right: EvaluationInput,
    *,
    config: DistributionEvaluationConfig | None = None,
) -> EvaluationReport:
    """Compare length, GC, k-mer, motif, and repeat distributions without p-values."""

    resolved = DistributionEvaluationConfig() if config is None else config
    if not isinstance(resolved, DistributionEvaluationConfig):
        raise TypeError("config must be DistributionEvaluationConfig or None.")
    left_items = materialize_input(left, limits=resolved.limits)
    right_items = materialize_input(right, limits=resolved.limits)
    require_nonempty(left_items, "distribution similarity")
    require_nonempty(right_items, "distribution similarity")
    if (
        sum(item.sequence.symbol_length for item in left_items + right_items)
        > resolved.limits.max_total_symbols
    ):
        from dnakit.exceptions import ConfigurationError

        raise ConfigurationError(
            "Combined distribution inputs exceed max_total_symbols.",
            code="EVALUATION_SYMBOL_LIMIT",
        )
    left_records = tuple(record_for(item) for item in left_items)
    right_records = tuple(record_for(item) for item in right_items)
    distances, details = _feature_distances(left_records, right_records, resolved)
    aggregate_distance = math.fsum(distances.values()) / len(distances)
    return report(
        name="distribution_similarity",
        method="mean-of-transparent-feature-distances",
        version="eval-distribution-v1",
        parameters={
            "features": resolved.features,
            "k": resolved.k,
            "canonical": resolved.canonical,
            "motifs": resolved.motifs,
            "aggregation": "1 - arithmetic mean of feature distances",
            "inference": "descriptive distances only; no hypothesis-test p-values",
            "motif_counting": "str.count non-overlapping",
            "gap_policy": "analyze text parts separately; no k-mer or motif crosses a Gap",
            "topology_policy": "linearized at stored origin; no circular-origin wrap",
            "ambiguity_policy": (
                "length retains IUPAC symbols; GC, k-mer, motif, and repeat calculations "
                "exclude affected noncanonical observations"
            ),
            "max_kmer_observations": resolved.max_kmer_observations,
            "max_repeat_comparisons_per_sequence": resolved.max_repeat_comparisons_per_sequence,
            "limits": resolved.limits,
        },
        metrics={
            "score": 1.0 - aggregate_distance,
            "distribution_similarity": 1.0 - aggregate_distance,
            "aggregate_distance": aggregate_distance,
            "feature_distances": distances,
            "feature_details": details,
            "left_count": len(left_records),
            "right_count": len(right_records),
        },
    )


__all__ = ["evaluate_distribution_similarity"]
