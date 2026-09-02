"""Tests for exact and reverse-complement dataset deduplication."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence
from dnakit.datasets import DeduplicationConfig, deduplicate
from dnakit.exceptions import ConfigurationError


def _record(record_id: str, symbols: str, **metadata: object) -> DNARecord:
    alphabet = DNAAlphabet.STRICT if set(symbols) <= set("ACGT") else DNAAlphabet.IUPAC
    return DNARecord(DNASequence(symbols, alphabet=alphabet), record_id, metadata=metadata)


def test_exact_dedup_preserves_first_group_order_and_audit() -> None:
    records = [_record("a", "ACGT"), _record("b", "TTAA"), _record("c", "ACGT")]

    result = deduplicate(records)

    assert result.records.ids == ("a", "b")
    assert result.input_count == 3
    assert result.output_count == 2
    assert result.duplicate_count == 1
    assert result.removed_count == 1
    assert result.groups[0].member_ids == ("a", "c")
    assert result.groups[0].orientations == ("forward", "forward")
    assert json.loads(json.dumps(result.to_dict()))["output_ids"] == ["a", "b"]


def test_reverse_complement_dedup_reports_member_direction() -> None:
    result = deduplicate(
        [_record("forward", "AAGC"), _record("reverse", "GCTT"), _record("other", "AAAA")],
        equivalence="reverse_complement",
    )

    assert result.records.ids == ("forward", "other")
    assert result.groups[0].member_ids == ("forward", "reverse")
    assert result.groups[0].orientations == ("forward", "reverse_complement")


def test_palindrome_and_declared_alphabet_do_not_create_false_groups() -> None:
    strict = DNARecord(DNASequence("ATAT"), "strict")
    iupac = DNARecord(DNASequence("ATAT", alphabet=DNAAlphabet.IUPAC), "iupac")

    result = deduplicate([strict, iupac], equivalence="reverse_complement")

    assert result.output_count == 1
    assert result.groups[0].orientations == ("forward", "forward")


def test_last_and_best_quality_representative_policies() -> None:
    low = DNARecord(DNASequence("AC"), "low", letter_annotations={"phred_quality": (5, 5)})
    high = DNARecord(DNASequence("AC"), "high", letter_annotations={"phred_quality": (30, 30)})

    last = deduplicate([low, high], config=DeduplicationConfig(representative_policy="last"))
    best = deduplicate(
        [low, high], config=DeduplicationConfig(representative_policy="best_quality")
    )

    assert last.records.ids == ("high",)
    assert best.records.ids == ("high",)


def test_reverse_complement_orientation_respects_nonfirst_representative() -> None:
    result = deduplicate(
        [_record("forward", "AAGC"), _record("reverse", "GCTT")],
        equivalence="reverse_complement",
        config=DeduplicationConfig(representative_policy="last"),
    )

    assert result.groups[0].representative_id == "reverse"
    assert result.groups[0].orientations == ("reverse_complement", "forward")


def test_metadata_conflicts_require_an_explicit_policy() -> None:
    records = [_record("positive", "AC", label=1), _record("negative", "AC", label=0)]
    with pytest.raises(ConfigurationError) as exc_info:
        deduplicate(records, config=DeduplicationConfig(conflict_field="label"))
    assert exc_info.value.code == "DEDUPLICATION_METADATA_CONFLICT"

    dropped = deduplicate(
        records,
        config=DeduplicationConfig(conflict_field="label", conflict_policy="drop_group"),
    )
    kept = deduplicate(
        records,
        config=DeduplicationConfig(conflict_field="label", conflict_policy="keep_all"),
    )
    representative = deduplicate(
        records,
        config=DeduplicationConfig(conflict_field="label", conflict_policy="keep_representative"),
    )

    assert dropped.output_count == 0
    assert dropped.groups[0].action == "dropped"
    assert kept.records.ids == ("positive", "negative")
    assert kept.groups[0].action == "kept_all"
    assert representative.records.ids == ("positive",)
    assert representative.conflicted_group_count == 1


def test_missing_conflict_value_differs_from_an_explicit_null() -> None:
    missing = _record("missing", "AC")
    explicit_null = _record("null", "AC", label=None)

    result = deduplicate(
        [missing, explicit_null],
        config=DeduplicationConfig(conflict_field="label", conflict_policy="keep_representative"),
    )

    assert result.groups[0].conflict
    assert result.groups[0].conflict_values == (None,)
    assert result.groups[0].missing_conflict_value_count == 1


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (True, 1),
        (True, 1.0),
        (1, 1.0),
        (
            {"outer": [True, {"value": 1.0}]},
            {"outer": [1, {"value": 1}]},
        ),
    ],
)
def test_metadata_conflicts_preserve_json_types_recursively(
    left: object,
    right: object,
) -> None:
    records = [
        _record("left", "AC", label=left),
        _record("right", "AC", label=right),
    ]

    result = deduplicate(
        records,
        config=DeduplicationConfig(
            conflict_field="label",
            conflict_policy="keep_representative",
        ),
    )

    assert result.groups[0].conflict
    assert len(result.groups[0].conflict_values) == 2


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (True, 1),
        (True, 1.0),
        (1, 1.0),
        (
            {"outer": [True, {"value": 1.0}]},
            {"outer": [1, {"value": 1}]},
        ),
    ],
)
def test_metadata_merge_does_not_merge_type_distinct_values(
    left: object,
    right: object,
) -> None:
    records = [
        _record("representative", "AC"),
        _record("left", "AC", label=left),
        _record("right", "AC", label=right),
    ]

    result = deduplicate(
        records,
        config=DeduplicationConfig(merge_metadata=True),
    )

    assert "label" not in result.records[0].metadata


def test_invalid_equivalence_and_non_records_are_rejected() -> None:
    with pytest.raises(ConfigurationError):
        deduplicate([], equivalence="rotation")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        deduplicate([DNASequence("AC")])  # type: ignore[list-item]
    invalid_configs: tuple[object, ...] = ({}, [], 0, False)
    for value in invalid_configs:
        with pytest.raises(ConfigurationError):
            deduplicate([], config=value)  # type: ignore[arg-type]


def test_metadata_merge_is_stably_ordered_and_audits_configuration() -> None:
    first = _record("first", "AC", z=1, a=2)
    second = _record("second", "AC", q=3, m=4)

    result = deduplicate(
        [first, second],
        config=DeduplicationConfig(merge_metadata=True),
    )

    assert tuple(result.records[0].metadata) == ("z", "a", "m", "q")
    assert result.merge_metadata
    assert result.to_dict()["merge_metadata"] is True
