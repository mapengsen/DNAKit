"""Tests for temporal/joint splits, leakage, and split-quality evaluation."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.datasets import (
    JointSplitConfig,
    LeakageConfig,
    SplitAssignment,
    TemporalSplitConfig,
    detect_leakage,
    evaluate_split_quality,
    joint_split,
    temporal_split,
)
from dnakit.exceptions import ConfigurationError


def _record(record_id: str, symbols: str = "AAAA", **metadata: object) -> DNARecord:
    return DNARecord(DNASequence(symbols), record_id, metadata=metadata)


def test_temporal_ratio_split_is_chronological_and_preserves_original_subset_order() -> None:
    records = [
        _record("late", date="2024-01-01"),
        _record("early", date="2020-01-01"),
        _record("middle", date="2022-01-01"),
        _record("latest", date="2025-01-01"),
    ]
    result = temporal_split(
        records,
        config=TemporalSplitConfig(ratios={"train": 0.5, "test": 0.5}),
    )

    assigned = {item.record_id: item.split for item in result.assignments}
    assert assigned == {"late": "test", "early": "train", "middle": "train", "latest": "test"}
    assert result.get("train").ids == ("early", "middle")
    assert result.get("test").ids == ("late", "latest")
    assert result.strategy == "chronological-largest-remainder"
    assert "UTC" in result.timezone_policy
    assert result.max_records == 1_000_000


def test_temporal_explicit_cutoff_boundary_is_in_earlier_split() -> None:
    records = [_record("on", date="2023-01-01"), _record("after", date="2023-01-02")]
    result = temporal_split(
        records,
        config=TemporalSplitConfig(
            ratios={"train": 0.5, "test": 0.5},
            cutoffs=("2023-01-01",),
        ),
    )
    assert {item.record_id: item.split for item in result.assignments} == {
        "on": "train",
        "after": "test",
    }


def test_temporal_timezones_are_compared_as_absolute_utc_instants() -> None:
    records = [
        _record("utc", date="2023-01-01T00:00:00Z"),
        _record("earlier-offset", date="2023-01-01T00:30:00+01:00"),
    ]
    result = temporal_split(
        records,
        config=TemporalSplitConfig(ratios={"train": 0.5, "test": 0.5}),
    )
    assert {item.record_id: item.split for item in result.assignments} == {
        "utc": "test",
        "earlier-offset": "train",
    }


def test_temporal_invalid_or_missing_dates_are_rejected() -> None:
    with pytest.raises(ConfigurationError) as missing:
        temporal_split([_record("missing")])
    assert missing.value.code == "TEMPORAL_METADATA_MISSING"
    with pytest.raises(ConfigurationError) as invalid:
        temporal_split([_record("bad", date="not-a-date")])
    assert invalid.value.code == "INVALID_TEMPORAL_METADATA"


def test_joint_split_keeps_multiple_group_constraints_atomic_and_records_heuristic() -> None:
    records = [
        _record("a1", donor="a", species="x", label=0),
        _record("a2", donor="a", species="y", label=1),
        _record("b1", donor="b", species="y", label=0),
        _record("c1", donor="c", species="z", label=1),
    ]
    result = joint_split(
        records,
        config=JointSplitConfig(
            ratios={"train": 0.75, "test": 0.25},
            group_keys=("donor", "species"),
            label_key="label",
            infeasible_policy="relax",
            ratio_tolerance=0.5,
            seed=7,
        ),
    )
    assigned = {item.record_id: item.split for item in result.assignments}
    assert assigned["a1"] == assigned["a2"] == assigned["b1"]
    assert "greedy" in result.strategy
    assert result.seed == 7
    assert result.similarity_canonical
    assert json.loads(json.dumps(result.to_dict()))["group_keys"] == ["donor", "species"]


def test_joint_split_reports_infeasible_or_explicit_relaxation() -> None:
    records = [_record(str(index), donor="same") for index in range(4)]
    config = JointSplitConfig(
        ratios={"train": 0.5, "test": 0.5},
        group_keys=("donor",),
        ratio_tolerance=0.1,
    )
    with pytest.raises(ConfigurationError) as error:
        joint_split(records, config=config)
    assert error.value.code == "JOINT_SPLIT_INFEASIBLE"

    relaxed = joint_split(
        records,
        config=JointSplitConfig(
            ratios={"train": 0.5, "test": 0.5},
            group_keys=("donor",),
            ratio_tolerance=0.1,
            infeasible_policy="relax",
        ),
    )
    assert not relaxed.feasible
    assert relaxed.relaxed
    assert relaxed.relaxations == ("ratio_tolerance",)


def test_leakage_detects_exact_and_high_similarity_cross_split_pairs() -> None:
    report = detect_leakage(
        {
            "train": DNASet([_record("a", "AAAA")]),
            "test": DNASet([_record("exact", "AAAA"), _record("near", "AAAT")]),
        },
        config=LeakageConfig(method="identity", threshold=0.7),
    )

    assert report.has_leakage
    assert report.exact_event_count == 1
    assert report.high_similarity_event_count == 1
    assert report.cross_pair_count == 2
    assert [(event.left_index, event.right_index) for event in report.events] == [(0, 0), (0, 1)]
    assert "exact dynamic-programming" in report.strictness
    assert [event.right_id for event in report.events] == ["exact", "near"]


def test_kmer_leakage_discloses_feature_based_false_negative_risk_and_limits() -> None:
    report = detect_leakage(
        {"train": DNASet([_record("a")]), "test": DNASet([_record("b", "AAAT")])},
        config=LeakageConfig(method="kmer", threshold=0.1, k=2),
    )
    assert "may miss biologically similar pairs" in report.strictness
    with pytest.raises(ConfigurationError) as pair_error:
        detect_leakage(
            {
                "a": DNASet([_record("a")]),
                "b": DNASet([_record("b"), _record("c")]),
            },
            config=LeakageConfig(max_cross_pairs=1, max_records=2, method="identity"),
        )
    assert pair_error.value.code == "LEAKAGE_RECORD_LIMIT"

    with pytest.raises(ConfigurationError) as cross_pair_error:
        detect_leakage(
            {
                "a": DNASet([_record("a")]),
                "b": DNASet([_record("b"), _record("c")]),
            },
            config=LeakageConfig(max_cross_pairs=1, max_records=3, method="identity"),
        )
    assert cross_pair_error.value.code == "LEAKAGE_PAIR_LIMIT"


def test_split_quality_definitions_cover_ratio_label_group_and_leakage() -> None:
    records = DNASet(
        [
            _record("a", label=0, donor="x"),
            _record("b", label=1, donor="x"),
            _record("c", label=0, donor="y"),
            _record("d", label=1, donor="z"),
        ]
    )
    assignments = tuple(
        SplitAssignment(index, record.id, "train" if index < 2 else "test")
        for index, record in enumerate(records)
    )
    quality = evaluate_split_quality(
        records,
        assignments,
        target_ratios={"train": 0.5, "test": 0.5},
        label_key="label",
        group_keys=("donor",),
    )

    assert quality.max_ratio_deviation == 0.0
    assert quality.label_total_variation_by_split == {"train": 0.0, "test": 0.0}
    assert quality.group_leak_count == 0
    assert quality.quality_score == 1.0

    with pytest.raises(ConfigurationError):
        evaluate_split_quality(
            records,
            assignments,
            target_ratios={"train": True, "test": 0.0},
        )

    with pytest.raises(ConfigurationError):
        evaluate_split_quality(
            records,
            assignments,
            target_ratios={"train": 0.5, "test": 0.5},
            group_keys=("",),
        )
