"""Tests for record operations with feature and annotation synchronization."""

import pytest

from dnakit.core import (
    CompoundLocation,
    DNAFeature,
    DNARecord,
    DNASequence,
    Gap,
    Interval,
    Strand,
    Topology,
    UnresolvedLocation,
)
from dnakit.exceptions import ConfigurationError, UnknownLengthError
from dnakit.ops import (
    canonical_origin_record,
    delete_record,
    insert_record,
    mask_record,
    reverse_complement_record,
    rotate_record,
    substitute_record,
    trim_record,
)


def _record(*, circular: bool = False) -> DNARecord:
    topology = Topology.CIRCULAR if circular else Topology.LINEAR
    return DNARecord(
        DNASequence("AACCGGTT", topology=topology),
        "record-1",
        description="fixture",
        features=[
            DNAFeature(
                "CDS",
                Interval(1, 7),
                id="cds-1",
                strand=Strand.FORWARD,
                label="coding",
                score=4.5,
                phase=1,
                qualifiers={"note": ["keep me"], "rank": 2},
                source="manual",
            ),
            DNAFeature("site", Interval(7, 8), id="tail", strand=Strand.BOTH),
        ],
        metadata={"sample": "x"},
        letter_annotations={"quality": [10, 11, 12, 13, 14, 15, 16, 17]},
    )


@pytest.mark.parametrize(
    ("policy", "expected_action", "expected_location"),
    [
        ("preserve", "resized", Interval(1, 5)),
        ("truncate", "truncated", Interval(1, 3)),
        ("split", "split", CompoundLocation((Interval(1, 3), Interval(3, 5)))),
        (
            "unresolved",
            "unresolved",
            UnresolvedLocation(
                "feature overlaps delete; exact biological extent is unresolved",
                [Interval(1, 3), Interval(3, 5)],
            ),
        ),
    ],
)
def test_delete_record_supports_overlap_policies_and_preserves_fields(
    policy: str,
    expected_action: str,
    expected_location: object,
) -> None:
    source = _record()

    result = delete_record(source, 3, 5, feature_policy=policy)  # type: ignore[arg-type]
    feature = result.record.features[0]

    assert result.sequence.symbols == "AACGTT"
    assert result.record.letter_annotations["quality"] == (10, 11, 12, 15, 16, 17)
    assert feature.location == expected_location
    assert feature.strand is Strand.FORWARD
    assert feature.qualifiers == source.features[0].qualifiers
    assert feature.id == source.features[0].id
    assert feature.type == source.features[0].type
    assert feature.label == source.features[0].label
    assert feature.score == source.features[0].score
    assert feature.phase == 1
    assert feature.source == source.features[0].source
    assert result.feature_changes[0].action == expected_action
    assert result.feature_changes[0].affected_edit_indices == (0,)
    assert result.record.metadata == source.metadata
    assert result.record.description == source.description
    assert source.sequence.symbols == "AACCGGTT"


def test_delete_record_can_remove_overlapping_feature_and_shift_later_feature() -> None:
    result = delete_record(_record(), 3, 5, feature_policy="delete")

    assert [feature.id for feature in result.record.features] == ["tail"]
    assert result.record.features[0].location == Interval(5, 6)
    assert result.feature_changes[0].action == "deleted"
    assert result.feature_changes[1].action == "shifted"
    assert result.deleted_feature_count == 1
    assert result.to_dict()["deleted_feature_count"] == 1


def test_insert_record_requires_or_drops_replacement_annotations_explicitly() -> None:
    source = _record()

    with pytest.raises(ConfigurationError) as missing:
        insert_record(source, 4, "TT")
    assert missing.value.code == "REPLACEMENT_ANNOTATIONS_REQUIRED"

    inserted = insert_record(
        source,
        4,
        "TT",
        replacement_annotations={"quality": [90, 91]},
    )
    assert inserted.sequence.symbols == "AACCTTGGTT"
    assert inserted.record.letter_annotations["quality"] == (
        10,
        11,
        12,
        13,
        90,
        91,
        14,
        15,
        16,
        17,
    )
    assert inserted.record.features[0].location == Interval(1, 9)

    dropped = insert_record(source, 4, "TT", letter_annotation_policy="drop")
    assert not dropped.record.letter_annotations
    assert dropped.letter_annotation_action == "dropped"


def test_insert_record_shifts_zero_length_feature_at_insertion_boundary() -> None:
    source = DNARecord(
        DNASequence("ACGT"),
        "point",
        features=[DNAFeature("cut", Interval(2, 2), id="cut-1")],
    )

    result = insert_record(source, 2, "AA")

    assert result.record.features[0].location == Interval(4, 4)
    assert result.feature_changes[0].action == "shifted"


def test_insert_record_excludes_boundary_insertions_from_half_open_features() -> None:
    source = DNARecord(
        DNASequence("ACGT"),
        "boundaries",
        features=[
            DNAFeature("left", Interval(0, 2), id="left"),
            DNAFeature("right", Interval(2, 4), id="right"),
        ],
    )

    result = insert_record(source, 2, "AA")

    assert result.record.features[0].location == Interval(0, 2)
    assert result.record.features[1].location == Interval(4, 6)
    assert [change.action for change in result.feature_changes] == ["preserved", "shifted"]


@pytest.mark.parametrize("operation", ["insert", "delete", "substitute"])
def test_empty_edits_do_not_remove_features_or_change_annotations(operation: str) -> None:
    source = _record()

    if operation == "insert":
        result = insert_record(source, 4, "", feature_policy="delete")
    elif operation == "delete":
        result = delete_record(source, 4, 4, feature_policy="delete")
    else:
        result = substitute_record(source, 4, 4, "", feature_policy="delete")

    assert result.sequence == source.sequence
    assert result.record.features == source.features
    assert result.record.letter_annotations == source.letter_annotations
    assert all(change.action == "preserved" for change in result.feature_changes)
    assert all(not change.affected_edit_indices for change in result.feature_changes)


def test_substitute_record_validates_annotation_keys_and_lengths() -> None:
    source = _record()

    with pytest.raises(ConfigurationError) as keys:
        substitute_record(source, 2, 4, "T", replacement_annotations={"other": [1]})
    assert keys.value.code == "REPLACEMENT_ANNOTATION_KEYS_MISMATCH"
    with pytest.raises(ConfigurationError) as length:
        substitute_record(source, 2, 4, "TT", replacement_annotations={"quality": [1]})
    assert length.value.code == "REPLACEMENT_ANNOTATION_LENGTH_MISMATCH"

    result = substitute_record(
        source,
        2,
        5,
        "T",
        feature_policy="unresolved",
        replacement_annotations={"quality": [99]},
    )
    assert result.sequence.symbols == "AATGTT"
    assert result.record.letter_annotations["quality"] == (10, 11, 99, 15, 16, 17)
    assert isinstance(result.record.features[0].location, UnresolvedLocation)


def test_mask_record_keeps_letter_alignment_and_audits_feature_policy() -> None:
    source = _record()

    result = mask_record(source, [(2, 4)], feature_policy="preserve")

    assert result.sequence.symbols == "AANNGGTT"
    assert result.record.letter_annotations == source.letter_annotations
    assert result.letter_annotation_action == "preserved"
    assert result.feature_changes[0].action == "preserved"
    assert result.parameters["feature_phase"] == "preserved_not_recomputed"


def test_trim_record_applies_two_original_coordinate_edits_in_safe_order() -> None:
    source = _record()

    result = trim_record(source, left=2, right=1, feature_policy="split")

    assert result.sequence.symbols == "CCGGT"
    assert result.record.letter_annotations["quality"] == (12, 13, 14, 15, 16)
    assert result.record.features[0].location == Interval(0, 5)
    assert result.feature_changes[0].affected_edit_indices == (0,)
    assert result.record.features[0].qualifiers == source.features[0].qualifiers
    assert result.feature_changes[1].action == "deleted"
    assert result.parameters["edit_coordinate_space"] == "original_record"


def test_zero_trim_is_audited_as_preserving_letter_annotations() -> None:
    source = _record()

    result = trim_record(source)

    assert result.record == source
    assert result.letter_annotation_action == "preserved"
    assert all(change.action == "preserved" for change in result.feature_changes)


def test_record_annotation_boundaries_at_gap_edges_map_to_symbol_offsets() -> None:
    source = DNARecord(
        DNASequence(["AA", Gap(3), "CC"]),
        "gap-edges",
        letter_annotations={"quality": [1, 2, 3, 4]},
    )

    at_left_edge = insert_record(
        source,
        2,
        "T",
        replacement_annotations={"quality": [9]},
    )
    at_right_edge = insert_record(
        source,
        5,
        "T",
        replacement_annotations={"quality": [9]},
    )

    assert at_left_edge.record.letter_annotations["quality"] == (1, 2, 9, 3, 4)
    assert at_right_edge.record.letter_annotations["quality"] == (1, 2, 9, 3, 4)


def test_reverse_complement_record_maps_location_strand_and_annotations() -> None:
    source = _record()

    result = reverse_complement_record(source)
    cds, tail = result.record.features

    assert result.sequence.symbols == "AACCGGTT"
    assert cds.location == Interval(1, 7)
    assert cds.strand is Strand.REVERSE
    assert cds.qualifiers == source.features[0].qualifiers
    assert cds.phase == source.features[0].phase
    assert tail.location == Interval(0, 1)
    assert tail.strand is Strand.BOTH
    assert result.record.letter_annotations["quality"] == tuple(range(17, 9, -1))
    assert result.feature_changes[0].action == "reverse_complemented"


def test_reverse_complement_unknown_gap_requires_explicit_feature_policy() -> None:
    source = DNARecord(
        DNASequence(["AA", Gap(None), "CC"]),
        "gapped",
        features=[DNAFeature("left", Interval(0, 2), strand=Strand.FORWARD)],
        letter_annotations={"quality": [1, 2, 3, 4]},
    )

    with pytest.raises(UnknownLengthError) as error:
        reverse_complement_record(source)
    assert error.value.code == "UNRESOLVED_REVERSE_FEATURE_COORDINATES"

    unresolved = reverse_complement_record(source, feature_policy="unresolved")
    assert isinstance(unresolved.record.features[0].location, UnresolvedLocation)
    assert unresolved.record.features[0].strand is Strand.REVERSE
    assert unresolved.record.letter_annotations["quality"] == (4, 3, 2, 1)

    deleted = reverse_complement_record(source, feature_policy="delete")
    assert not deleted.record.features


@pytest.mark.parametrize(
    ("policy", "expected_action", "location_type"),
    [
        ("split", "split", CompoundLocation),
        ("preserve", "preserved", CompoundLocation),
        ("truncate", "truncated", Interval),
        ("unresolved", "unresolved", UnresolvedLocation),
    ],
)
def test_rotate_record_maps_origin_crossing_features_and_letter_annotations(
    policy: str,
    expected_action: str,
    location_type: type[object],
) -> None:
    source = _record(circular=True)

    result = rotate_record(source, 3, feature_policy=policy)  # type: ignore[arg-type]

    assert result.sequence.symbols == "CGGTTAAC"
    assert result.record.letter_annotations["quality"] == (13, 14, 15, 16, 17, 10, 11, 12)
    assert isinstance(result.record.features[0].location, location_type)
    assert result.feature_changes[0].action == expected_action
    assert result.record.features[0].qualifiers == source.features[0].qualifiers
    assert result.parameters["effective_offset"] == 3


def test_rotate_record_delete_policy_removes_only_origin_crossing_features() -> None:
    result = rotate_record(_record(circular=True), 3, feature_policy="delete")

    assert [feature.id for feature in result.record.features] == ["tail"]
    assert result.record.features[0].location == Interval(4, 5)
    assert result.feature_changes[0].action == "deleted"


def test_rotate_record_does_not_treat_preexisting_compound_as_new_origin_split() -> None:
    source = DNARecord(
        DNASequence("AACCGGTT", topology=Topology.CIRCULAR),
        "compound",
        features=[
            DNAFeature(
                "compound",
                CompoundLocation((Interval(0, 1), Interval(4, 5))),
                id="compound-1",
            )
        ],
    )

    result = rotate_record(source, 2, feature_policy="delete")

    assert len(result.record.features) == 1
    assert result.feature_changes[0].action == "rotated"
    assert result.record.features[0].location == CompoundLocation((Interval(6, 7), Interval(2, 3)))


def test_rotate_record_maps_known_gap_boundaries_features_and_annotations() -> None:
    gap = Gap(3, metadata={"source": "assembly"})
    source = DNARecord(
        DNASequence(["AA", gap, "CC"], topology=Topology.CIRCULAR),
        "gapped-circle",
        features=[DNAFeature("right", Interval(5, 7), strand=Strand.REVERSE)],
        letter_annotations={"quality": [1, 2, 3, 4]},
    )

    result = rotate_record(source, 2)

    assert result.sequence.parts == (gap, "CCAA")
    assert result.record.features[0].location == Interval(3, 5)
    assert result.record.features[0].strand is Strand.REVERSE
    assert result.record.letter_annotations["quality"] == (3, 4, 1, 2)
    assert result.feature_changes[0].action == "rotated"


def test_canonical_origin_record_uses_sequence_rule_and_preserves_audit() -> None:
    source = DNARecord(
        DNASequence("TACACA", topology=Topology.CIRCULAR),
        "circle",
        features=[DNAFeature("site", Interval(0, 2), id="f")],
        letter_annotations={"quality": [0, 1, 2, 3, 4, 5]},
    )

    result = canonical_origin_record(source)

    assert result.sequence.symbols == "ACACAT"
    assert result.parameters["effective_offset"] == 1
    assert result.parameters["requested_offset"] is None
    assert result.parameters["algorithm"] == "booth"
    assert result.parameters["origin_rule"] == (
        "lexicographically_minimal_forward_rotation_then_smallest_offset"
    )
    assert result.record.letter_annotations["quality"] == (1, 2, 3, 4, 5, 0)
    assert result.operation == "canonical_origin"


def test_record_operations_reject_invalid_policies() -> None:
    with pytest.raises(ConfigurationError) as feature:
        delete_record(_record(), 1, 2, feature_policy="bad")  # type: ignore[arg-type]
    assert feature.value.code == "INVALID_FEATURE_OVERLAP_POLICY"
    with pytest.raises(ConfigurationError) as annotation:
        delete_record(_record(), 1, 2, letter_annotation_policy="bad")  # type: ignore[arg-type]
    assert annotation.value.code == "INVALID_LETTER_ANNOTATION_POLICY"
