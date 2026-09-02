"""Boundary tests for subsequence, editing, trimming, and masking."""

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap, GapKind, Topology
from dnakit.exceptions import (
    ConfigurationError,
    CoordinateError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)
from dnakit.ops import delete, insert, mask, subsequence, substitute, trim


def test_subsequence_uses_zero_based_half_open_coordinates() -> None:
    sequence = DNASequence("AACCGG")

    assert subsequence(sequence, 1, 5).symbols == "ACCG"
    assert subsequence(sequence, 3, 3) == DNASequence("")
    with pytest.raises(CoordinateError):
        subsequence(sequence, -1, 2)
    with pytest.raises(CoordinateError):
        subsequence(sequence, 0, 7)


def test_subsequence_gap_policy_is_explicit_and_partial_gap_is_preserved() -> None:
    gap = Gap(4, GapKind.SCAFFOLD, crossable=False, evidence=("assembly",))
    sequence = DNASequence(["AA", gap, "CC"])

    with pytest.raises(UnsupportedGapOperationError) as error:
        subsequence(sequence, 1, 7)
    assert error.value.code == "EDIT_OVERLAPS_GAP"

    result = subsequence(sequence, 1, 7, allow_gaps=True)
    assert result.parts == ("A", Gap(4, GapKind.SCAFFOLD, False, ("assembly",)), "C")
    partial = subsequence(sequence, 3, 5, allow_gaps=True)
    assert partial.parts == (Gap(2, GapKind.SCAFFOLD, False, ("assembly",)),)

    with pytest.raises(UnknownLengthError):
        subsequence(DNASequence(["AA", Gap(None), "CC"]), 0, 2)


def test_insert_delete_and_substitute_return_audited_new_objects() -> None:
    source = DNASequence("AACCGG")

    inserted = insert(source, 2, "TT")
    assert inserted.sequence.symbols == "AATTCCGG"
    assert inserted.edits[0].kind == "insert"
    assert (inserted.edits[0].start, inserted.edits[0].end) == (2, 2)
    assert inserted.edits[0].replacement_symbols == "TT"
    assert inserted.edits[0].removed_symbols == ""
    deleted = delete(source, 2, 4)
    assert deleted.sequence.symbols == "AAGG"
    assert deleted.edits[0].removed_parts == ("CC",)
    assert deleted.edits[0].removed_symbols == "CC"
    replaced = substitute(source, 2, 4, "TN")
    assert replaced.sequence.symbols == "AATNGG"
    assert replaced.sequence.alphabet is DNAAlphabet.IUPAC
    assert replaced.edits[0].removed_symbols == "CC"
    assert source.symbols == "AACCGG"


def test_empty_edits_are_well_defined() -> None:
    source = DNASequence("ACGT")

    assert delete(source, 2, 2).sequence == source
    assert substitute(source, 2, 2, "N").sequence.symbols == "ACNGT"
    empty_mask = mask(source, [])
    assert empty_mask.sequence == source
    assert empty_mask.sequence is not source
    assert empty_mask.edits == ()


def test_edit_refuses_unknown_or_interior_gap_coordinates() -> None:
    known = DNASequence(["AA", Gap(3), "CC"])
    unknown = DNASequence(["AA", Gap(None), "CC"])

    assert insert(known, 2, "T").sequence.parts == ("AAT", Gap(3), "CC")
    with pytest.raises(UnsupportedGapOperationError):
        insert(known, 3, "T")
    with pytest.raises(UnsupportedGapOperationError):
        delete(known, 1, 6)
    with pytest.raises(UnknownLengthError):
        insert(unknown, 1, "T")


def test_trim_handles_known_gaps_but_not_circular_rotation_or_unknown_gaps() -> None:
    gap = Gap(4, GapKind.SCAFFOLD)
    sequence = DNASequence(["AA", gap, "CC"])

    trimmed = trim(sequence, left=3, right=1)
    assert trimmed.sequence.parts == (
        Gap(3, GapKind.SCAFFOLD),
        "C",
    )
    assert [(edit.start, edit.end) for edit in trimmed.edits] == [(0, 3), (7, 8)]
    assert trimmed.edits[0].removed_parts == ("AA", Gap(1, GapKind.SCAFFOLD))
    assert trimmed.edits[1].removed_parts == ("C",)
    assert trim(DNASequence("ACGT"), left=4).sequence == DNASequence("")
    with pytest.raises(ConfigurationError):
        trim(sequence, left=5, right=4)
    with pytest.raises(UnknownLengthError):
        trim(DNASequence(["AA", Gap(None), "CC"]), left=1)
    with pytest.raises(ConfigurationError):
        insert(
            DNASequence("AC", topology=Topology.CIRCULAR),
            1,
            DNASequence("GT", topology=Topology.CIRCULAR),
        )
    with pytest.raises(ConfigurationError) as circular_trim:
        trim(DNASequence("AC", topology=Topology.CIRCULAR), left=1)
    assert circular_trim.value.code == "CIRCULAR_TRIM_NOT_SUPPORTED"


def test_edit_argument_type_errors_are_dnakit_errors() -> None:
    with pytest.raises(CoordinateError):
        insert(DNASequence("AC"), "1", "T")  # type: ignore[call-overload]
    with pytest.raises(CoordinateError):
        delete(DNASequence("AC"), True, 1)


def test_mask_validates_intervals_symbol_and_gap_overlap() -> None:
    source = DNASequence("AACCGGTT")

    result = mask(source, [(1, 3), (5, 7)])
    assert result.sequence.symbols == "ANNCGNNT"
    assert result.sequence.alphabet is DNAAlphabet.IUPAC
    assert [(edit.start, edit.end) for edit in result.edits] == [(1, 3), (5, 7)]
    assert [edit.removed_symbols for edit in result.edits] == ["AC", "GT"]
    with pytest.raises(ConfigurationError):
        mask(source, [(1, 4), (3, 5)])
    with pytest.raises(ConfigurationError):
        mask(source, [(1, 2)], symbol="x")
    with pytest.raises(UnsupportedGapOperationError):
        mask(DNASequence(["AA", Gap(2), "CC"]), [(1, 3)])
    with pytest.raises(ConfigurationError):
        mask(source, None)  # type: ignore[call-overload]
    with pytest.raises(ConfigurationError):
        mask(source, [(1, 2, 3)])  # type: ignore[list-item]


def test_gapped_edit_replacement_is_rejected_without_losing_audit_information() -> None:
    replacement = DNASequence(["A", Gap(2), "T"])

    with pytest.raises(UnsupportedGapOperationError) as error:
        insert(DNASequence("CC"), 1, replacement)
    assert error.value.code == "GAPPED_EDIT_FRAGMENT"
