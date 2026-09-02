"""Tests for features, records, and materialized collections."""

from collections.abc import Iterator
from dataclasses import FrozenInstanceError

import pytest

from dnakit.core import (
    CompoundLocation,
    DNAAlphabet,
    DNAFeature,
    DNARecord,
    DNASequence,
    DNASet,
    Gap,
    Interval,
    Strand,
    UnresolvedLocation,
)
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError, DuplicateIDError, FeatureError, SequenceError


def test_feature_validates_fields_and_defensively_freezes_qualifiers() -> None:
    raw_qualifiers = {"note": ["first"]}
    feature = DNAFeature(
        "motif",
        Interval(1, 3),
        id="f1",
        strand=Strand.FORWARD,
        phase=0,
        qualifiers=raw_qualifiers,
    )
    raw_qualifiers["note"].append("changed")

    assert feature.qualifiers["note"] == ("first",)
    assert feature.strand is Strand.FORWARD
    with pytest.raises(FrozenInstanceError):
        feature.label = "new"  # type: ignore[misc]
    with pytest.raises(TypeError):
        hash(feature)

    for invalid_phase in (True, False, -1, 3, 1.0):
        with pytest.raises(FeatureError):
            DNAFeature("CDS", Interval(0, 1), phase=invalid_phase)  # type: ignore[arg-type]


def test_record_quality_length_counts_every_iupac_symbol() -> None:
    sequence = DNASequence("ACGN", alphabet=DNAAlphabet.IUPAC)
    record = DNARecord(
        sequence,
        "read-1",
        letter_annotations={"phred_quality": [40, 39, 38, 10]},
    )

    assert record.letter_annotations["phred_quality"] == (40, 39, 38, 10)
    with pytest.raises(SequenceError) as error:
        DNARecord(
            sequence,
            "read-2",
            letter_annotations={"phred_quality": [40, 39, 38]},
        )
    assert error.value.code == "LETTER_ANNOTATION_LENGTH_MISMATCH"


def test_record_letter_annotations_exclude_gap_span_but_not_symbols() -> None:
    sequence = DNASequence(["AC", Gap(100), "N"], alphabet=DNAAlphabet.IUPAC)
    record = DNARecord(sequence, "gapped", letter_annotations={"confidence": [1.0, 1.0, 0.5]})

    assert len(record.letter_annotations["confidence"]) == sequence.symbol_length == 3


def test_record_bounds_letter_annotation_iterables() -> None:
    consumed = 0

    def values() -> Iterator[int]:
        nonlocal consumed
        while True:
            consumed += 1
            yield 1

    with pytest.raises(SequenceError) as error:
        DNARecord(DNASequence("A"), "bounded", letter_annotations={"quality": values()})
    assert error.value.code == "LETTER_ANNOTATION_LENGTH_MISMATCH"
    assert consumed == 2


def test_core_json_rejects_cycles_and_deep_values_but_allows_shared_dag() -> None:
    cycle: list[object] = []
    cycle.append(cycle)
    with pytest.raises(ConfigurationError) as recursive:
        DNARecord(DNASequence("A"), "cycle", metadata={"value": cycle})
    assert recursive.value.code == "JSON_RECURSIVE_REFERENCE"

    deep: object = "leaf"
    for _ in range(70):
        deep = [deep]
    with pytest.raises(ConfigurationError) as depth:
        DNARecord(DNASequence("A"), "deep", metadata={"value": deep})
    assert depth.value.code == "JSON_STRUCTURE_LIMIT"

    shared = [1, 2]
    record = DNARecord(DNASequence("A"), "dag", metadata={"value": [shared, shared]})
    assert record.metadata["value"] == ((1, 2), (1, 2))


def test_record_rejects_resolved_features_outside_known_sequence_span() -> None:
    with pytest.raises(SequenceError) as error:
        DNARecord(
            DNASequence("AC"),
            "out-of-bounds",
            features=[DNAFeature("site", Interval(1, 3))],
        )
    assert error.value.code == "FEATURE_OUT_OF_BOUNDS"

    with pytest.raises(SequenceError) as compound_error:
        DNARecord(
            DNASequence("AC"),
            "compound-out-of-bounds",
            features=[DNAFeature("site", CompoundLocation([Interval(0, 1), Interval(2, 3)]))],
        )
    assert compound_error.value.code == "FEATURE_OUT_OF_BOUNDS"


def test_record_requires_unresolved_location_after_unknown_length_gap() -> None:
    sequence = DNASequence(["AC", Gap(None), "GT"])
    before_gap = DNAFeature("known", Interval(0, 2))
    unresolved = DNAFeature("unknown", UnresolvedLocation("after unknown gap"))

    record = DNARecord(sequence, "partially-resolved", features=[before_gap, unresolved])

    assert record.features == (before_gap, unresolved)
    with pytest.raises(SequenceError) as error:
        DNARecord(
            sequence,
            "incorrectly-resolved",
            features=[DNAFeature("unknown", Interval(2, 3))],
        )
    assert error.value.code == "FEATURE_LOCATION_UNRESOLVED"


def test_record_rejects_unresolved_location_anchors_outside_resolvable_span() -> None:
    with pytest.raises(SequenceError) as known_error:
        DNARecord(
            DNASequence("AC"),
            "bad-anchor",
            features=[
                DNAFeature(
                    "unknown",
                    UnresolvedLocation("uncertain boundary", [Interval(999, 1000)]),
                )
            ],
        )
    assert known_error.value.code == "FEATURE_OUT_OF_BOUNDS"

    with pytest.raises(SequenceError) as unknown_error:
        DNARecord(
            DNASequence(["AC", Gap(None), "GT"]),
            "bad-anchor-unknown-span",
            features=[
                DNAFeature(
                    "unknown",
                    UnresolvedLocation("after unknown gap", [Interval(2, 3)]),
                )
            ],
        )
    assert unknown_error.value.code == "FEATURE_LOCATION_UNRESOLVED"


def test_record_and_set_copy_input_containers() -> None:
    feature_list = [DNAFeature("site", Interval(0, 1))]
    raw_metadata = {"sample": {"batch": 2}}
    record = DNARecord(
        DNASequence("ACGT"),
        "r1",
        features=feature_list,
        metadata=raw_metadata,
    )
    source_records = [record]
    dataset = DNASet(source_records, name="example", metadata={"labels": ["a"]})
    feature_list.clear()
    source_records.clear()
    raw_metadata["sample"]["batch"] = 9

    assert len(record.features) == 1
    frozen_sample = record.metadata["sample"]
    assert isinstance(frozen_sample, FrozenDict)
    assert frozen_sample["batch"] == 2
    assert len(dataset) == 1
    assert dataset[0] is record
    assert dataset[:1].records == (record,)
    with pytest.raises(TypeError):
        hash(record)
    with pytest.raises(TypeError):
        hash(dataset)


def test_from_sequences_generates_stable_ids_and_supports_custom_factory() -> None:
    sequences = [DNASequence("A"), DNASequence("C")]

    default = DNASet.from_sequences(iter(sequences))
    custom = DNASet.from_sequences(sequences, id_factory=lambda index, sequence: f"x{index}")

    assert default.ids == ("sequence_1", "sequence_2")
    assert custom.ids == ("x0", "x1")


def test_duplicate_ids_are_reportable_but_keyed_lookup_rejects_ambiguity() -> None:
    sequence = DNASequence("AC")
    dataset = DNASet([DNARecord(sequence, "same"), DNARecord(sequence, "same")])

    assert dataset.ids == ("same", "same")
    with pytest.raises(DuplicateIDError):
        dataset.get("same")
    with pytest.raises(KeyError):
        dataset.get("missing")
