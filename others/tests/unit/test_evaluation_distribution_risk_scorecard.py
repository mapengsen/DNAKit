"""Tests for distribution comparison, synthesis risk, and scorecards."""

from __future__ import annotations

import hashlib
import json
from importlib import import_module
from typing import cast

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet, Gap, MetricResult
from dnakit.evaluation import (
    DistributionEvaluationConfig,
    ScorecardConfig,
    ScoreRule,
    SynthesisRiskConfig,
    evaluate_distribution_similarity,
    evaluate_scorecard,
    evaluate_synthesis_risk,
)
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError


def _set(*pairs: tuple[str, str]) -> DNASet:
    return DNASet(DNARecord(DNASequence(symbols), record_id) for record_id, symbols in pairs)


def test_identical_distributions_have_unit_similarity_with_transparent_methods() -> None:
    records = _set(("a", "ACGTACGT"), ("b", "GGGGCCCC"))
    report = evaluate_distribution_similarity(
        records,
        records,
        config=DistributionEvaluationConfig(k=2, motifs=("CG", "GG")),
    )

    assert report.metrics["score"] == 1.0
    feature_distances = cast(dict[str, object], report.to_dict()["metrics"]["feature_distances"])
    assert set(feature_distances) == {
        "length",
        "gc",
        "kmer",
        "motif",
        "repeat",
    }
    assert cast(str, report.parameters["inference"]).startswith("descriptive distances only")


def test_distribution_rejects_empty_input_and_observation_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(ConfigurationError) as empty:
        evaluate_distribution_similarity(DNASet([]), _set(("x", "AAAA")))
    assert empty.value.code == "EMPTY_EVALUATION_DATASET"
    distribution_module = import_module("dnakit.evaluation.distribution")
    counter_calls = 0
    original_kmer_counts = distribution_module._kmer_counts

    def counted_kmer_counts(*args: object, **kwargs: object) -> object:
        nonlocal counter_calls
        counter_calls += 1
        return original_kmer_counts(*args, **kwargs)

    monkeypatch.setattr(distribution_module, "_kmer_counts", counted_kmer_counts)
    with pytest.raises(ConfigurationError) as bounded:
        evaluate_distribution_similarity(
            _set(("a", "A" * 100_000)),
            _set(("b", "AAAA")),
            config=DistributionEvaluationConfig(
                features=("kmer",),
                k=1,
                max_kmer_observations=2,
            ),
        )
    assert bounded.value.code == "DISTRIBUTION_KMER_LIMIT"
    assert bounded.value.context["observation_count"] == 3
    assert counter_calls == 0


def test_distribution_kmers_and_motifs_never_cross_gap_boundaries() -> None:
    gapped = DNASet([DNARecord(DNASequence(["A", Gap(3), "A"]), "gapped")])
    report = evaluate_distribution_similarity(
        gapped,
        _set(("separated", "AT")),
        config=DistributionEvaluationConfig(
            features=("kmer", "motif"),
            k=2,
            motifs=("AA",),
        ),
    )
    details = cast(dict[str, object], report.to_dict()["metrics"]["feature_details"])
    kmer = cast(dict[str, object], details["kmer"])
    motif = cast(dict[str, object], details["motif"])
    assert kmer["left_observations"] == 0
    assert cast(dict[str, float], motif["left_rates"])["AA"] == 0.0
    assert "no k-mer or motif crosses" in cast(str, report.parameters["gap_policy"])


def test_synthesis_risk_reports_rules_locations_and_non_experimental_disclaimer() -> None:
    risky = DNASequence("G" * 30 + "AT" * 10 + "C" * 30)
    report = evaluate_synthesis_risk(
        risky,
        config=SynthesisRiskConfig(
            window_size=10,
            window_step=5,
            homopolymer_threshold=8,
            tandem_min_unit=2,
            tandem_max_unit=2,
            tandem_min_repeats=4,
            inverted_min_arm=4,
            inverted_max_arm=4,
            inverted_max_loop=4,
        ),
    )
    entry = report.entries[0]

    assert cast(float, entry.metrics["risk_score"]) > 0
    assert entry.metrics["canonical_sequence_length"] == risky.symbol_length
    assert (
        entry.metrics["canonical_sequence_sha256"]
        == hashlib.sha256(risky.symbols.encode("ascii")).hexdigest()
    )
    assert cast(tuple[object, ...], entry.metrics["risky_gc_windows"])
    assert "not a vendor acceptance" in cast(str, report.metrics["disclaimer"])
    assert cast(str, report.parameters["structure_method"]).endswith("no folding backend")
    assert json.loads(json.dumps(report.to_dict()))["entries"][0]["issues"]


def test_synthesis_risk_rejects_gap_and_circular_topology() -> None:
    with pytest.raises(UnsupportedGapOperationError):
        evaluate_synthesis_risk(DNASequence(["AAAA", Gap(2), "TTTT"]))
    with pytest.raises(ConfigurationError) as circular:
        evaluate_synthesis_risk(DNASequence("AAAA", topology="circular"))
    assert circular.value.code == "SYNTHESIS_RISK_LINEAR_REQUIRED"


def test_synthesis_window_limit_precedes_descriptor_allocation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    synthesis_module = import_module("dnakit.evaluation.synthesis")
    descriptor_called = False

    def unexpected_window_descriptors(*args: object, **kwargs: object) -> object:
        nonlocal descriptor_called
        descriptor_called = True
        raise AssertionError("window_descriptors must not run for an over-budget input")

    monkeypatch.setattr(
        synthesis_module,
        "window_descriptors",
        unexpected_window_descriptors,
    )
    with pytest.raises(ConfigurationError) as bounded:
        evaluate_synthesis_risk(
            DNASequence("A" * 100_000),
            config=SynthesisRiskConfig(
                window_size=50,
                window_step=1,
                max_windows_per_sequence=2,
            ),
        )

    assert bounded.value.code == "SYNTHESIS_RISK_WINDOW_LIMIT"
    assert bounded.value.context["window_count"] == 100_000
    assert descriptor_called is False


def test_scorecard_supports_reports_metrics_direction_clipping_and_missing_policies() -> None:
    source = MetricResult(
        "risk",
        0.25,
        method="fixture",
        algorithm_version="1",
    )
    config = ScorecardConfig(
        rules={
            "benefit": ScoreRule(weight=2, minimum=0, maximum=10),
            "risk": ScoreRule(
                direction="lower_is_better",
                weight=1,
                minimum=0,
                maximum=1,
            ),
            "missing": ScoreRule(weight=1),
        },
        missing_policy="zero",
    )
    report = evaluate_scorecard(
        {"benefit": 15.0, "risk": source},
        config=config,
    )

    assert report.metrics["score"] == pytest.approx((2 + 0.75) / 4)
    assert report.metrics["missing_components"] == ("missing",)
    contributions = cast(list[dict[str, object]], report.to_dict()["metrics"]["contributions"])
    first = contributions[0]
    assert first["clipped"] is True
    assert cast(str, report.parameters["formula"]).startswith("sum(weight_i")

    omitted = evaluate_scorecard(
        {"benefit": 5.0},
        config=ScorecardConfig(
            rules={"benefit": ScoreRule(), "missing": ScoreRule()},
            missing_policy="omit",
        ),
    )
    assert omitted.metrics["included_weight_sum"] == 1.0
    with pytest.raises(ConfigurationError) as missing:
        evaluate_scorecard(
            {},
            config=ScorecardConfig(rules={"needed": ScoreRule()}),
        )
    assert missing.value.code == "MISSING_SCORECARD_COMPONENT"


def test_scorecard_can_select_named_evaluation_report_metric() -> None:
    distribution = evaluate_distribution_similarity(
        _set(("a", "AAAA")),
        _set(("b", "AAAA")),
        config=DistributionEvaluationConfig(features=("length",)),
    )
    scorecard = evaluate_scorecard(
        {"distribution": distribution},
        config=ScorecardConfig(rules={"distribution": ScoreRule(metric="distribution_similarity")}),
    )

    assert scorecard.metrics["score"] == 1.0


def test_scorecard_config_defensively_freezes_rule_mapping() -> None:
    rules = {"value": ScoreRule()}
    config = ScorecardConfig(rules=rules)
    rules["added"] = ScoreRule()

    assert tuple(config.rules) == ("value",)
    with pytest.raises(TypeError):
        config.rules["added"] = ScoreRule()  # type: ignore[index]
