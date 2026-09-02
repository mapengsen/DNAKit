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

    assert diversity.metrics["mean_pair_distance"] == pytest.approx(0.75)
    assert diversity.metrics["cluster_count"] == 2
    assert redundancy.metrics["near_pair_count"] == 1
    assert cast(str, redundancy.parameters["score_formula"]).startswith("mean(")
    assert singleton.metrics["score"] == 1.0
    assert singleton.metrics["mean_pair_distance"] is None


def test_levenshtein_diversity_uses_the_published_mean_pairwise_formula() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "AATT"))
    config = DiversityEvaluationConfig(calculation="levenshtein")

    diversity = evaluate_diversity(records, config=config)
    singleton = evaluate_diversity(_set(("only", "AAAA")), config=config)

    assert diversity.method == "mean-pairwise-levenshtein-distance"
    assert diversity.algorithm_version == "eval-diversity-levenshtein-v1"
    assert diversity.metrics["score"] == pytest.approx(4 / 3)
    assert diversity.metrics["mean_pairwise_levenshtein_distance"] == pytest.approx(4 / 3)
    assert diversity.metrics["minimum_pairwise_levenshtein_distance"] == 1.0
    assert diversity.metrics["maximum_pairwise_levenshtein_distance"] == 2.0
    assert diversity.metrics["pairwise_comparison_count"] == 3
    assert diversity.parameters["normalization"] == "none"
    assert singleton.metrics["score"] is None


def test_collection_pairwise_limit_is_checked() -> None:
    records = _set(("a", "AAAA"), ("b", "AAAT"), ("c", "AATT"))
    config = DiversityEvaluationConfig(
        limits=EvaluationLimits(max_pairwise_comparisons=2),
    )
    with pytest.raises(ConfigurationError) as error:
        evaluate_diversity(records, config=config)
    assert error.value.code == "EVALUATION_PAIRWISE_LIMIT"

    with pytest.raises(ConfigurationError) as levenshtein_error:
        evaluate_diversity(
            records,
            config=DiversityEvaluationConfig(
                calculation="levenshtein",
                limits=EvaluationLimits(max_pairwise_comparisons=2),
            ),
        )
    assert levenshtein_error.value.code == "EVALUATION_PAIRWISE_LIMIT"


def test_empty_collection_has_declared_vacuous_uniqueness_and_diversity() -> None:
    uniqueness = evaluate_uniqueness(DNASet([]))
    diversity = evaluate_diversity(DNASet([]))

    assert uniqueness.metrics["score"] == 1.0
    assert diversity.metrics["score"] == 1.0
    assert diversity.metrics["cluster_coverage"] == 1.0
