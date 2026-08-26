"""Unit tests for the unified validation entry point and QC rules."""

from __future__ import annotations

import json
import math
from collections.abc import Iterator

import pytest

from dnakit.core.collection import DNASet
from dnakit.core.enums import DNAAlphabet
from dnakit.core.gap import Gap
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.standardize import (
    DatasetValidationConfig,
    ValidationConfig,
    validate,
)


def codes(report: object) -> set[str]:
    return {issue.code for issue in report.issues}  # type: ignore[attr-defined]


def test_empty_sequence_is_invalid_by_default_and_configurable() -> None:
    sequence = DNASequence("")

    report = validate(sequence)
    allowed = validate(
        sequence,
        config=ValidationConfig(allow_empty=True, sequence_length=0),
    )

    assert not report.is_valid
    assert codes(report) == {"STD_EMPTY_SEQUENCE"}
    assert allowed.is_valid


def test_validation_can_require_an_exact_symbol_length() -> None:
    sequence = DNASequence(["AC", Gap(4), "GT"])

    matching = validate(sequence, config=ValidationConfig(sequence_length=4))
    mismatching = validate(sequence, config=ValidationConfig(sequence_length=5))

    assert matching.is_valid
    assert matching.symbol_length == 4
    assert matching.coordinate_span == 8
    assert not mismatching.is_valid
    assert codes(mismatching) == {"STD_SEQUENCE_LENGTH_MISMATCH"}
    issue = mismatching.issues[0]
    assert issue.details == {"length": 4, "expected_length": 5}


def test_full_iupac_sequence_ambiguity_threshold() -> None:
    sequence = DNASequence("RYSWKMBDHVN", alphabet=DNAAlphabet.IUPAC)
    report = validate(sequence, config=ValidationConfig(max_ambiguity_fraction=0.5))

    assert not report.is_valid
    assert report.ambiguity.total_count == 11
    assert report.ambiguity.fraction == 1.0
    assert "STD_AMBIGUITY_FRACTION_HIGH" in codes(report)


def test_validate_can_apply_stricter_alphabet_than_sequence_declaration() -> None:
    sequence = DNASequence("ACN", alphabet=DNAAlphabet.IUPAC)
    report = validate(sequence, config=ValidationConfig(alphabet=DNAAlphabet.STRICT))

    assert not report.is_valid
    assert "STD_INVALID_SYMBOL" in codes(report)


def test_unknown_gap_length_and_ambiguity_denominator_are_explicit() -> None:
    sequence = DNASequence(["AN", Gap(None), "GT"], alphabet=DNAAlphabet.IUPAC)
    report = validate(
        sequence,
        config=ValidationConfig(
            allow_unknown_gap_length=False,
            ambiguity_denominator_includes_gap=True,
        ),
    )

    assert not report.is_valid
    assert report.coordinate_span is None
    assert report.ambiguity.denominator is None
    assert report.ambiguity.fraction is None
    assert {"STD_UNKNOWN_GAP_LENGTH", "STD_AMBIGUITY_DENOMINATOR_UNKNOWN"} <= codes(report)


def test_record_metadata_and_phred_quality_qc() -> None:
    record = DNARecord(
        DNASequence("ACGT"),
        "record-1",
        metadata={"sample": "S1"},
        letter_annotations={"phred_quality": (40, 30, 20, 10)},
    )
    report = validate(
        record,
        config=ValidationConfig(
            required_metadata_fields=("sample", "organism"),
            required_letter_annotations=("phred_quality",),
            minimum_mean_phred=30.0,
        ),
    )

    assert not report.is_valid
    assert report.record_id == "record-1"
    assert report.quality is not None
    assert report.quality.mean == 25.0
    assert {"STD_REQUIRED_METADATA_MISSING", "STD_QUALITY_MEAN_LOW"} <= codes(report)


def test_phred_range_check_reports_all_out_of_range_positions() -> None:
    record = DNARecord(
        DNASequence("ACGT"),
        "record-1",
        letter_annotations={"phred_quality": (40, 94, 30, -1)},
    )
    report = validate(record)

    assert not report.is_valid
    assert "STD_QUALITY_OUT_OF_RANGE" in codes(report)
    issue = next(item for item in report.issues if item.code == "STD_QUALITY_OUT_OF_RANGE")
    assert issue.details["indices"] == (1, 3)


def test_validate_accepts_collections_and_detects_duplicate_ids_with_stable_indices() -> None:
    records = [
        DNARecord(DNASequence("AC"), "same"),
        DNARecord(DNASequence("GT"), "unique"),
        DNARecord(DNASequence("TG"), "same"),
    ]

    report = validate(record for record in records)

    assert report.record_count == 3
    assert not report.ids_unique
    assert not report.is_valid
    assert report.duplicate_ids[0].id == "same"
    assert report.duplicate_ids[0].indices == (0, 2)
    assert "STD_DUPLICATE_RECORD_ID" in codes(report)
    assert report.record_reports is not None
    assert len(report.record_reports) == 3
    json.dumps(report.to_dict())


def test_validate_accepts_dataset_config_and_can_warn_on_duplicates() -> None:
    records = DNASet(
        [
            DNARecord(DNASequence("AC"), "same"),
            DNARecord(DNASequence("GT"), "same"),
        ]
    )
    report = validate(
        records,
        config=DatasetValidationConfig(
            require_unique_ids=False,
            collect_record_reports=False,
        ),
    )

    assert report.ids_unique is False
    assert report.is_valid
    assert report.record_reports is None
    duplicate_issue = next(
        issue for issue in report.issues if issue.code == "STD_DUPLICATE_RECORD_ID"
    )
    assert duplicate_issue.severity.value == "warning"


def test_validate_collection_is_invalid_when_any_record_is_invalid() -> None:
    records = [
        DNARecord(DNASequence("ACGT"), "valid"),
        DNARecord(DNASequence("AC"), "short"),
    ]

    report = validate(records, config=ValidationConfig(min_length=3))

    assert report.record_count == 2
    assert not report.is_valid
    assert report.record_reports is not None
    assert report.record_reports[0].is_valid
    assert not report.record_reports[1].is_valid
    assert "STD_INVALID_RECORDS" in codes(report)


def test_empty_dataset_is_invalid() -> None:
    report = validate([])

    assert report.record_count == 0
    assert not report.is_valid
    assert codes(report) == {"STD_EMPTY_DATASET"}


@pytest.mark.parametrize("value", [True, 1.5, "2"])
def test_validation_length_thresholds_reject_non_integer_values(value: object) -> None:
    with pytest.raises(ConfigurationError):
        ValidationConfig(min_length=value)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        ValidationConfig(sequence_length=value)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf, True, "10"])
def test_validation_quality_thresholds_require_finite_numbers(value: object) -> None:
    with pytest.raises(ConfigurationError):
        ValidationConfig(minimum_phred=value)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        ValidationConfig(minimum_mean_phred=value)  # type: ignore[arg-type]


def test_validation_config_rejects_invalid_flags_and_required_fields() -> None:
    with pytest.raises(ConfigurationError):
        ValidationConfig(allow_empty=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        ValidationConfig(required_metadata_fields="sample")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        DatasetValidationConfig(require_unique_ids=1)  # type: ignore[arg-type]


def test_validation_required_field_iterable_is_bounded(monkeypatch: pytest.MonkeyPatch) -> None:
    import dnakit.standardize.config as config_module

    monkeypatch.setattr(config_module, "_MAX_REQUIRED_FIELDS", 2)
    consumed = 0

    def fields() -> Iterator[str]:
        nonlocal consumed
        while True:
            consumed += 1
            yield f"field-{consumed}"

    with pytest.raises(ConfigurationError) as error:
        ValidationConfig(required_metadata_fields=fields())  # type: ignore[arg-type]
    assert error.value.code == "REQUIRED_FIELD_LIMIT"
    assert consumed == 3


def test_validation_entry_points_reject_falsy_non_config_objects() -> None:
    invalid_configs: tuple[object, ...] = ({}, [], 0, False)
    for value in invalid_configs:
        with pytest.raises(ConfigurationError):
            validate(DNASequence("AC"), config=value)  # type: ignore[call-overload]
        with pytest.raises(ConfigurationError):
            validate([], config=value)  # type: ignore[call-overload]
