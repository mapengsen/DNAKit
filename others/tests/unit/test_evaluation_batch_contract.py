"""Shared DNA batch input contract across public sequence evaluation APIs."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import pytest

from dnakit import DNA
from dnakit.evaluation import (
    FragmentSimilarityConfig,
    SNNConfig,
    create_reference_library,
    evaluate_ambiguity,
    evaluate_complexity,
    evaluate_distribution_similarity,
    evaluate_diversity,
    evaluate_fragment_similarity,
    evaluate_memorization,
    evaluate_novelty,
    evaluate_quality,
    evaluate_redundancy,
    evaluate_reference_similarity,
    evaluate_snn,
    evaluate_synthesis_risk,
    evaluate_uniqueness,
    evaluate_validity,
    nearest_reference,
)
from dnakit.evaluation.results import EvaluationReport
from dnakit.exceptions import ConfigurationError


def test_per_sequence_metrics_share_dna_batch_ids_and_summary() -> None:
    queries = DNA(["ACGTACGT", "ACGTACGA", "TTTTTTTT"])
    expected_ids = queries.ids
    evaluators: tuple[Callable[[DNA], EvaluationReport], ...] = (
        evaluate_validity,
        evaluate_ambiguity,
        evaluate_complexity,
        evaluate_quality,
        evaluate_synthesis_risk,
    )

    for evaluate in evaluators:
        result = evaluate(queries)
        assert result.metrics["record_count"] == 3
        assert tuple(entry.subject_id for entry in result.entries) == expected_ids
        assert tuple(entry.input_index for entry in result.entries) == (0, 1, 2)

    validity = evaluate_validity(queries, sequence_length=8)
    assert validity.metrics["valid_fraction"] == 1.0
    complexity = evaluate_complexity(queries)
    assert complexity.metrics["score"] == pytest.approx(
        sum(cast(float, entry.metrics["score"]) for entry in complexity.entries) / 3
    )


def test_collection_metrics_accept_dna_batch_and_keep_collection_semantics() -> None:
    queries = DNA(["ACGTACGT", "ACGTACGT", "TTTTTTTT"])
    uniqueness = evaluate_uniqueness(queries)
    diversity = evaluate_diversity(queries)
    redundancy = evaluate_redundancy(queries)

    assert uniqueness.metrics["record_count"] == 3
    assert uniqueness.metrics["duplicate_groups"] == (("sequence_1", "sequence_2"),)
    assert diversity.metrics["record_count"] == 3
    assert diversity.metrics["pairwise_comparison_count"] == 3
    assert redundancy.metrics["record_count"] == 3
    assert all(not result.entries for result in (uniqueness, diversity, redundancy))


def test_reference_and_two_collection_metrics_accept_dna_batches() -> None:
    queries = DNA(["ACGTACGT", "ACGTACGA", "TTTTTTTT"])
    reference_dna = DNA(["ACGTACGT", "CCCCCCCC"])
    reference = create_reference_library(
        reference_dna, name="test", version="1", source="local:test"
    )

    for evaluate in (
        nearest_reference,
        evaluate_novelty,
        evaluate_memorization,
        evaluate_reference_similarity,
    ):
        result = evaluate(queries, reference)
        assert tuple(entry.subject_id for entry in result.entries) == queries.ids

    fragment = evaluate_fragment_similarity(
        queries, reference_dna, config=FragmentSimilarityConfig(show_progress=False)
    )
    snn = evaluate_snn(queries, reference_dna, config=SNNConfig(show_progress=False))
    distribution = evaluate_distribution_similarity(queries, reference_dna)

    assert fragment.metrics["generated_count"] == 3
    assert fragment.metrics["reference_count"] == 2
    assert not fragment.entries
    assert tuple(entry.subject_id for entry in snn.entries) == queries.ids
    assert distribution.metrics["left_count"] == 3
    assert distribution.metrics["right_count"] == 2
    assert not distribution.entries


def test_bare_string_list_error_points_to_dna_batch_constructor() -> None:
    raw_sequences: object = ["ACGT", "ACGA"]
    with pytest.raises(ConfigurationError, match=r"DNA\(\[\.\.\.\]\)") as error:
        evaluate_validity(raw_sequences)  # type: ignore[arg-type]
    assert error.value.code == "INVALID_EVALUATION_RECORD"
