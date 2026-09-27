"""Tests for versioned reference, nearest-hit, novelty, and memorization reports."""

from __future__ import annotations

import json
from importlib import import_module
from typing import cast

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.evaluation import (
    EvaluationLimits,
    ReferenceSearchConfig,
    create_reference_library,
    evaluate_memorization,
    evaluate_novelty,
    evaluate_reference_similarity,
    nearest_reference,
)
from dnakit.evaluation.results import ReferenceLibrary
from dnakit.exceptions import ConfigurationError


def _set(*pairs: tuple[str, str]) -> DNASet:
    return DNASet(DNARecord(DNASequence(symbols), record_id) for record_id, symbols in pairs)


def _reference() -> ReferenceLibrary:
    return create_reference_library(
        _set(("train-a", "AAAA"), ("train-c", "CCCC")),
        name="training",
        version="2026.1",
        source="local:test-fixture",
        date="2026-08-13",
        filters={"split": "train"},
        index_parameters={"kind": "exhaustive"},
    )


def test_reference_digest_is_deterministic_content_bound_and_json_safe() -> None:
    first = _reference()
    second = _reference()
    assert first.digest == second.digest
    assert first.digest_scope.startswith("ordered IDs")
    payload = json.loads(json.dumps(first.to_dict()))
    assert payload["record_count"] == 2
    assert "records" not in payload


def test_reference_library_requires_unique_ids() -> None:
    with pytest.raises(ConfigurationError):
        create_reference_library(
            DNASet(
                [
                    DNARecord(DNASequence("AAAA"), "same"),
                    DNARecord(DNASequence("CCCC"), "same"),
                ]
            ),
            name="bad",
            version="1",
            source="local:test",
        )


def test_nearest_reference_stable_tie_break_and_provenance() -> None:
    reference = create_reference_library(
        _set(("first", "AAAA"), ("second", "AAAA")),
        name="reference",
        version="1",
        source="local:test",
    )
    report = nearest_reference(
        _set(("query", "AAAA")),
        reference,
        config=ReferenceSearchConfig(method="exact", top_k=2),
    )

    entries = cast(list[dict[str, object]], report.to_dict()["entries"])
    metrics = cast(dict[str, object], entries[0]["metrics"])
    hits = cast(list[dict[str, object]], metrics["hits"])
    assert [hit["reference_id"] for hit in hits] == [
        "first",
        "second",
    ]
    assert report.provenance.reference is not None
    assert report.provenance.reference.checksum == reference.digest


def test_nearest_reference_top_k_heap_is_bounded_and_stable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_module = import_module("dnakit.evaluation.reference")
    real_heappush = reference_module.heapq.heappush
    largest_heap = 0

    def tracking_heappush(heap: list[object], item: object) -> None:
        nonlocal largest_heap
        real_heappush(heap, item)
        largest_heap = max(largest_heap, len(heap))

    monkeypatch.setattr(reference_module.heapq, "heappush", tracking_heappush)
    reference_count = 2_000
    reference = create_reference_library(
        DNASet(
            DNARecord(DNASequence("AAAA"), f"reference-{index:04d}")
            for index in range(reference_count)
        ),
        name="large-reference",
        version="1",
        source="local:test",
    )
    report = nearest_reference(
        _set(("query", "AAAA")),
        reference,
        config=ReferenceSearchConfig(
            method="exact",
            top_k=3,
            limits=EvaluationLimits(
                max_records=reference_count,
                max_pairwise_comparisons=reference_count,
            ),
        ),
    )

    hits = cast(tuple[object, ...], report.entries[0].metrics["hits"])
    assert [cast(dict[str, object], hit)["reference_id"] for hit in hits] == [
        "reference-0000",
        "reference-0001",
        "reference-0002",
    ]
    assert largest_heap == 3


def test_novelty_memorization_and_reference_similarity_share_reference_contract() -> None:
    reference = _reference()
    queries = _set(("copy", "AAAA"), ("novel", "GGGG"))
    config = ReferenceSearchConfig(method="identity", copy_threshold=0.9)
    novelty = evaluate_novelty(queries, reference, config=config)
    memorization = evaluate_memorization(queries, reference, config=config)
    similarity = evaluate_reference_similarity(queries, reference, config=config)

    assert novelty.entries[0].metrics["novelty"] == 0.0
    assert novelty.entries[1].metrics["novelty"] == 1.0
    assert novelty.metrics["score"] == 0.5
    assert novelty.metrics["exact_match_count"] == 1
    assert novelty.metrics["exact_match_rate"] == 0.5
    assert "is_novel" not in novelty.entries[0].metrics
    assert memorization.metrics["exact_copy_count"] == 1
    assert memorization.parameters["task_model_used"] is False
    assert similarity.metrics["mean_nearest_similarity"] == 0.5


def test_levenshtein_novelty_uses_the_published_nearest_reference_formula() -> None:
    reference = _reference()
    queries = _set(("copy", "AAAA"), ("near", "AAAT"), ("tie", "GGGG"))

    novelty = evaluate_novelty(
        queries,
        reference,
        config=ReferenceSearchConfig(novelty_calculation="levenshtein"),
    )

    assert novelty.method == "mean-nearest-reference-levenshtein-distance"
    assert novelty.algorithm_version == "eval-novelty-levenshtein-v2"
    assert novelty.metrics["score"] == pytest.approx(5 / 3)
    assert novelty.metrics["mean_nearest_levenshtein_distance"] == pytest.approx(5 / 3)
    assert novelty.metrics["exact_match_count"] == 1
    assert novelty.metrics["exact_match_rate"] == pytest.approx(1 / 3)
    assert novelty.metrics["pairwise_comparison_count"] == 6
    assert novelty.entries[0].metrics["nearest_levenshtein_distance"] == 0.0
    assert novelty.entries[1].metrics["nearest_levenshtein_distance"] == 1.0
    assert novelty.entries[2].metrics["nearest_reference_id"] == "train-a"
    assert novelty.provenance.reference is not None
    assert novelty.provenance.reference.checksum == reference.digest


def test_levenshtein_novelty_requires_a_nonempty_reference() -> None:
    empty_reference = create_reference_library(
        DNASet([]),
        name="empty",
        version="1",
        source="local:test",
    )

    with pytest.raises(ConfigurationError) as error:
        evaluate_novelty(
            _set(("query", "AAAA")),
            empty_reference,
            config=ReferenceSearchConfig(novelty_calculation="levenshtein"),
        )
    assert error.value.code == "EMPTY_EVALUATION_DATASET"


def test_reference_search_limit_precedes_pairwise_work() -> None:
    reference = _reference()
    queries = _set(("a", "AAAA"), ("b", "CCCC"))
    config = ReferenceSearchConfig(
        limits=EvaluationLimits(max_pairwise_comparisons=3),
    )
    with pytest.raises(ConfigurationError) as error:
        nearest_reference(queries, reference, config=config)
    assert error.value.code == "EVALUATION_PAIRWISE_LIMIT"


def test_reference_coverage_filter_and_copy_threshold_conflict_are_explicit() -> None:
    reference = create_reference_library(
        _set(("long", "AAAAAAAA")),
        name="reference",
        version="1",
        source="local:test",
    )
    filtered = nearest_reference(
        _set(("short", "AAAA")),
        reference,
        config=ReferenceSearchConfig(method="identity", min_coverage=0.75),
    )
    assert filtered.entries[0].metrics["nearest_reference_id"] is None
    novelty = evaluate_novelty(
        _set(("query", "AAAA")),
        reference,
        config=ReferenceSearchConfig(min_similarity=0.95, copy_threshold=0.9),
    )
    assert novelty.metrics["exact_match_rate"] == 0.0
    with pytest.raises(ConfigurationError) as error:
        evaluate_memorization(
            _set(("query", "AAAA")),
            reference,
            config=ReferenceSearchConfig(min_similarity=0.95, copy_threshold=0.9),
        )
    assert error.value.code == "REFERENCE_COPY_THRESHOLD_CONFLICT"


def test_novelty_exact_match_rate_is_independent_of_nearest_kmer_ties() -> None:
    reference = create_reference_library(
        _set(("reverse", "GGGGTTTT"), ("literal", "AAAACCCC")),
        name="reference",
        version="1",
        source="local:test",
    )
    report = evaluate_novelty(_set(("query", "AAAACCCC")), reference)

    assert report.entries[0].metrics["nearest_reference_id"] == "reverse"
    assert report.entries[0].metrics["novelty"] == 0.0
    assert report.entries[0].metrics["exact_match"] is True
    assert report.metrics["exact_match_rate"] == 1.0


def test_novelty_score_and_exact_fraction_ignore_copy_threshold() -> None:
    reference = _reference()
    queries = _set(("copy", "AAAA"), ("near", "AAAT"))
    low = evaluate_novelty(
        queries, reference, config=ReferenceSearchConfig(method="identity", copy_threshold=0.0)
    )
    high = evaluate_novelty(
        queries, reference, config=ReferenceSearchConfig(method="identity", copy_threshold=1.0)
    )

    assert low.metrics["score"] == high.metrics["score"] == pytest.approx(0.125)
    assert low.metrics["exact_match_rate"] == high.metrics["exact_match_rate"] == 0.5
    assert low.entries[1].metrics["novelty"] == pytest.approx(0.25)
    assert low.entries[1].metrics["exact_match"] is False
