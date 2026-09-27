"""Tests for uniqueness, diversity, and redundancy reports."""

from __future__ import annotations

from typing import cast

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.evaluation import (
    DiversityEvaluationConfig,
    EvaluationLimits,
    UniquenessEvaluationConfig,
    evaluate_diversity,
    evaluate_kmer_diversity,
    evaluate_redundancy,
    evaluate_uniqueness,
)
from dnakit.exceptions import ConfigurationError


def _set(*pairs: tuple[str, str]) -> DNASet:
    return DNASet(
        DNARecord(DNASequence(symbols, alphabet="iupac"), record_id) for record_id, symbols in pairs
    )


def test_exact_and_reverse_complement_uniqueness_include_groups() -> None:
    records = _set(("a", "AAGC"), ("b", "AAGC"), ("rc", "GCTT"))
    exact = evaluate_uniqueness(records)
    reverse = evaluate_uniqueness(
        records,
        config=UniquenessEvaluationConfig(equivalence="reverse_complement"),
    )

    assert exact.metrics["uniqueness_score"] == pytest.approx(2 / 3)
    assert exact.metrics["duplicate_groups"] == (("a", "b"),)
    assert reverse.metrics["uniqueness_score"] == pytest.approx(1 / 3)


def test_iupac_and_approximate_uniqueness_are_explicit() -> None:
    iupac = evaluate_uniqueness(
        _set(("ambiguous", "AN"), ("strict", "AT")),
        config=UniquenessEvaluationConfig(equivalence="iupac"),
    )
    near = evaluate_uniqueness(
        _set(("a", "AAAA"), ("b", "AAAT"), ("c", "CCCC")),
        config=UniquenessEvaluationConfig(
            equivalence="approximate",
            approximate_method="identity",
            threshold=0.7,
        ),
    )

    assert iupac.metrics["unique_count"] == 1
    assert near.metrics["unique_count"] == 2
    assert "threshold-graph" in cast(str, near.parameters["grouping_strategy"])


def test_diversity_and_redundancy_publish_definitions_and_singleton_semantics() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "CCCC"))
    config = DiversityEvaluationConfig(method="identity", cluster_threshold=0.7)
    diversity = evaluate_diversity(records, config=config)
    redundancy = evaluate_redundancy(records, config=config)
    singleton = evaluate_diversity(_set(("only", "AAAA")), config=config)

    assert diversity.metrics["mean_pairwise_levenshtein_distance"] == pytest.approx(3.0)
    assert "cluster_count" not in diversity.metrics
    assert redundancy.metrics["near_pair_count"] == 1
    assert cast(str, redundancy.parameters["score_formula"]).startswith("mean(")
    assert singleton.metrics["score"] is None
    assert singleton.metrics["mean_pair_distance"] is None


def test_levenshtein_diversity_uses_the_published_mean_pairwise_formula() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "AATT"))
    config = DiversityEvaluationConfig()
    assert config.workers == 16

    diversity = evaluate_diversity(records, config=config)
    singleton = evaluate_diversity(_set(("only", "AAAA")), config=config)

    assert diversity.method == "mean-pairwise-levenshtein-distance"
    assert diversity.algorithm_version == "eval-diversity-levenshtein-v1"
    assert diversity.metrics["score"] == pytest.approx(4 / 3)
    assert diversity.metrics["mean_pairwise_levenshtein_distance"] == pytest.approx(4 / 3)
    assert "minimum_pairwise_levenshtein_distance" not in diversity.metrics
    assert "maximum_pairwise_levenshtein_distance" not in diversity.metrics
    assert diversity.metrics["pairwise_comparison_count"] == 3
    assert diversity.parameters["normalization"] == "none"
    assert singleton.metrics["score"] is None


def test_parallel_diversity_matches_serial_and_respects_limits() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "AATT"), ("d", "TT"))
    serial = evaluate_diversity(records, config=DiversityEvaluationConfig(workers=1))
    parallel = evaluate_diversity(records, config=DiversityEvaluationConfig(workers=2))

    assert parallel.metrics == serial.metrics
    assert parallel.parameters["workers"] == 2
    assert (
        evaluate_diversity(
            _set(("only", "AAAA")), config=DiversityEvaluationConfig(workers=2)
        ).metrics["score"]
        is None
    )

    with pytest.raises(ConfigurationError) as error:
        evaluate_diversity(
            records,
            config=DiversityEvaluationConfig(
                workers=2, limits=EvaluationLimits(max_alignment_cells=20)
            ),
        )
    assert error.value.code == "EDIT_DISTANCE_CELL_LIMIT"


@pytest.mark.parametrize("workers", [0, -1, True, 1.5])
def test_diversity_rejects_invalid_worker_counts(workers: object) -> None:
    with pytest.raises(ConfigurationError):
        DiversityEvaluationConfig(workers=workers)  # type: ignore[arg-type]


def test_kmer_diversity_reports_pairwise_and_nearest_jaccard_distances() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "CCCC"))

    result = evaluate_kmer_diversity(records, k=2)
    singleton = evaluate_kmer_diversity(_set(("only", "AAAA")), k=2)

    assert result.metrics["record_count"] == 3
    assert result.metrics["pairwise_comparison_count"] == 3
    assert result.metrics["mean_pair_distance"] == pytest.approx(5 / 6)
    assert result.metrics["mean_nearest_neighbor_distance"] == pytest.approx(2 / 3)
    assert singleton.metrics["mean_nearest_neighbor_distance"] is None
    reverse_pair = _set(("forward", "AAGC"), ("reverse", "GCTT"))
    assert evaluate_kmer_diversity(reverse_pair, k=2).metrics["mean_pair_distance"] == 0
    assert evaluate_kmer_diversity(reverse_pair, k=2, canonical=False).metrics[
        "mean_pair_distance"
    ] == pytest.approx(0.8)
    with pytest.raises(ConfigurationError):
        evaluate_kmer_diversity(records, k=2, limits=EvaluationLimits(max_pairwise_comparisons=2))


def test_diversity_uses_raw_edit_distance_for_different_lengths() -> None:
    records = _set(("a", "AAAA"), ("b", "AAA"), ("c", "T"))

    diversity = evaluate_diversity(records)

    assert diversity.metrics["score"] == pytest.approx(8 / 3)
    assert diversity.metrics["pairwise_comparison_count"] == 3


def test_collection_pairwise_limit_is_checked() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "AATT"))
    config = DiversityEvaluationConfig(
        limits=EvaluationLimits(max_pairwise_comparisons=2),
    )
    with pytest.raises(ConfigurationError) as error:
        evaluate_diversity(records, config=config)
    assert error.value.code == "EVALUATION_PAIRWISE_LIMIT"

    with pytest.raises(ConfigurationError):
        DiversityEvaluationConfig(calculation="similarity")


def test_empty_collection_has_declared_vacuous_uniqueness_and_diversity() -> None:
    uniqueness = evaluate_uniqueness(DNASet([]))
    diversity = evaluate_diversity(DNASet([]))

    assert uniqueness.metrics["score"] == 1.0
    assert diversity.metrics["score"] is None
    assert diversity.metrics["pairwise_comparison_count"] == 0
