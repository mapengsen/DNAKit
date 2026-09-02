"""Tests for validity, ambiguity, complexity, and quality evaluation."""

from __future__ import annotations

import json
from typing import cast

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet, Gap
from dnakit.evaluation import (
    AmbiguityEvaluationConfig,
    ComplexityEvaluationConfig,
    EvaluationLimits,
    QualityEvaluationConfig,
    evaluate_ambiguity,
    evaluate_complexity,
    evaluate_quality,
    evaluate_validity,
)
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.standardize import ValidationConfig


def _record(record_id: str, symbols: str) -> DNARecord:
    return DNARecord(DNASequence(symbols, alphabet="iupac"), record_id)


def test_validity_wraps_validation_and_serializes_empty_set() -> None:
    report = evaluate_validity(
        DNASet([_record("ok", "ACGT"), _record("empty", "")]),
        config=ValidationConfig(allow_empty=False),
    )
    empty = evaluate_validity(DNASet([]))

    assert report.metrics["valid_count"] == 1
    assert report.entries[1].issues[0].code == "STD_EMPTY_SEQUENCE"
    assert empty.metrics["record_count"] == 0
    assert empty.metrics["valid_fraction"] is None
    assert json.loads(json.dumps(report.to_dict()))["algorithm_version"] == "eval-validity-v1"


def test_ambiguity_reports_positions_weights_denominator_and_gap_policy() -> None:
    sequence = DNASequence(["ARN", Gap(2), "V"], alphabet="iupac")
    report = evaluate_ambiguity(
        sequence,
        config=AmbiguityEvaluationConfig(
            max_fraction=0.5,
            gap_denominator_policy="include_known",
        ),
    )
    entry = report.entries[0]

    assert entry.metrics["positions"] == (1, 2, 3)
    assert entry.metrics["denominator"] == 6
    assert entry.metrics["weighted_ambiguity_count"] == 2.25
    assert entry.metrics["ambiguity_fraction"] == 0.5
    with pytest.raises(UnsupportedGapOperationError) as error:
        evaluate_ambiguity(
            sequence,
            config=AmbiguityEvaluationConfig(gap_denominator_policy="error"),
        )
    assert error.value.code == "AMBIGUITY_GAP_REJECTED"


def test_unknown_gap_makes_included_denominator_explicitly_undefined() -> None:
    report = evaluate_ambiguity(
        DNASequence(["AN", Gap(None), "T"], alphabet="iupac"),
        config=AmbiguityEvaluationConfig(gap_denominator_policy="include_known"),
    )

    assert report.entries[0].metrics["ambiguity_fraction"] is None
    assert report.entries[0].issues[0].code == "EVAL_AMBIGUITY_DENOMINATOR_UNKNOWN"


def test_complexity_distinguishes_low_and_balanced_sequences() -> None:
    report = evaluate_complexity(
        DNASet([_record("low", "A" * 24), _record("balanced", "ACGT" * 6)])
    )

    low_score = cast(float, report.entries[0].metrics["score"])
    balanced_score = cast(float, report.entries[1].metrics["score"])
    assert low_score < balanced_score
    assert report.entries[0].metrics["repeat_fraction"] == 1.0
    assert report.parameters["entropy_normalization"] == "min(1,max(0,H_base2/2))"


def test_quality_combines_components_and_audits_missing_phred_policy() -> None:
    good = DNARecord(
        DNASequence("ACGT" * 5),
        "good",
        letter_annotations={"phred_quality": [30] * 20},
    )
    poor = _record("poor", "N" * 20)
    report = evaluate_quality(
        DNASet([good, poor]),
        config=QualityEvaluationConfig(min_length=4),
    )

    good_score = cast(float, report.entries[0].metrics["score"])
    poor_score = cast(float, report.entries[1].metrics["score"])
    assert good_score > poor_score
    assert report.parameters["missing_phred_policy"] == "neutral completeness=1"
    assert json.loads(json.dumps(report.to_dict()))["entries"][0]["subject_id"] == "good"


def test_sequence_evaluation_enforces_materialization_limits_before_work() -> None:
    records = (_record(str(index), "ACGT") for index in range(3))
    with pytest.raises(ConfigurationError) as error:
        evaluate_validity(records, limits=EvaluationLimits(max_records=2))
    assert error.value.code == "EVALUATION_RECORD_LIMIT"


def test_quality_rejects_conflicting_nested_resource_limits() -> None:
    with pytest.raises(ConfigurationError) as error:
        QualityEvaluationConfig(
            ambiguity=AmbiguityEvaluationConfig(limits=EvaluationLimits(max_records=1)),
            complexity=ComplexityEvaluationConfig(limits=EvaluationLimits(max_records=2)),
        )
    assert error.value.code == "QUALITY_LIMIT_MISMATCH"


def test_quality_does_not_treat_undefined_ambiguity_as_perfect() -> None:
    limits = EvaluationLimits()
    report = evaluate_quality(
        DNASequence(["AN", Gap(None), "T"], alphabet="iupac"),
        config=QualityEvaluationConfig(
            ambiguity=AmbiguityEvaluationConfig(
                gap_denominator_policy="include_known",
                limits=limits,
            ),
            complexity=ComplexityEvaluationConfig(limits=limits),
        ),
    )

    entries = cast(list[dict[str, object]], report.to_dict()["entries"])
    metrics = cast(dict[str, object], entries[0]["metrics"])
    components = cast(dict[str, float], metrics["components"])
    assert components["ambiguity"] == 0.0
    assert report.parameters["undefined_ambiguity_policy"] == "ambiguity component is zero"
