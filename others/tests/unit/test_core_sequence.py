"""Boundary tests for immutable DNA sequence and gap values."""

from dataclasses import FrozenInstanceError

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap, GapKind, Topology
from dnakit.exceptions import (
    InvalidAlphabetError,
    SequenceError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)


def test_strict_and_iupac_alphabets_are_explicit() -> None:
    sequence = DNASequence("ACGTRYN", alphabet=DNAAlphabet.IUPAC)

    assert sequence.symbol_length == 7
    assert sequence.canonical_base_count == 4
    assert sequence.ambiguity_count == 3
    assert sequence.reverse_complement().symbols == "NRYACGT"

    with pytest.raises(InvalidAlphabetError, match="not valid"):
        DNASequence("ACGN")
    with pytest.raises(InvalidAlphabetError, match="normalize"):
        DNASequence("acgt")
    with pytest.raises(InvalidAlphabetError, match="ASCII"):
        DNASequence("核酸".encode())


def test_empty_linear_is_representable_but_empty_circular_is_rejected() -> None:
    empty = DNASequence("")

    assert empty.parts == ()
    assert empty.length == 0
    assert len(empty) == 0

    with pytest.raises(SequenceError) as error:
        DNASequence("", topology=Topology.CIRCULAR)
    assert error.value.code == "EMPTY_CIRCULAR_SEQUENCE"


def test_known_and_unknown_gap_lengths_are_not_guessed() -> None:
    known = DNASequence(["AC", Gap(5, GapKind.SCAFFOLD), "GT"])
    unknown = DNASequence(["AC", Gap(None), "GT"])

    assert known.parts == ("AC", Gap(5, GapKind.SCAFFOLD), "GT")
    assert known.symbol_length == 4
    assert known.coordinate_span == 9
    assert len(known) == 9
    assert unknown.coordinate_span is None
    assert unknown.length is None
    assert unknown.has_unknown_length
    with pytest.raises(UnknownLengthError):
        len(unknown)
    with pytest.raises(UnsupportedGapOperationError):
        known.to_string()


def test_sequence_parts_are_merged_and_transform_with_gap_order() -> None:
    gap = Gap(3, crossable=False, evidence=["paired-ends"])
    sequence = DNASequence(["A", "C", gap, "G", "A"])

    assert sequence.parts == ("AC", gap, "GA")
    assert sequence.reverse().parts == ("AG", gap, "CA")
    assert sequence.complement().parts == ("TG", gap, "CT")
    assert sequence.reverse_complement().parts == ("TC", gap, "GT")


def test_from_fragments_checks_gap_context_without_copying_flanks() -> None:
    gap = Gap(None, metadata={"source": "AGP"})
    sequence = DNASequence.from_fragments([b"AC", "GT"], [gap])

    assert sequence.parts == ("AC", gap, "GT")
    assert not hasattr(gap, "upstream")
    assert not hasattr(gap, "downstream")

    with pytest.raises(SequenceError) as error:
        DNASequence.from_fragments(["AC"], [gap])
    assert error.value.code == "FRAGMENT_GAP_COUNT_MISMATCH"


def test_core_sequence_and_nested_metadata_are_immutable() -> None:
    raw_metadata = {"evidence": ["map", {"score": 2}]}
    gap = Gap(2, metadata=raw_metadata)
    sequence = DNASequence(["AC", gap, "GT"])
    raw_metadata["evidence"].append("changed")

    frozen_evidence = gap.metadata["evidence"]
    assert isinstance(frozen_evidence, tuple)
    assert frozen_evidence[0] == "map"
    assert len(frozen_evidence) == 2
    assert hash(gap)
    assert hash(sequence)
    with pytest.raises(FrozenInstanceError):
        sequence.topology = Topology.CIRCULAR  # type: ignore[misc]
    with pytest.raises(TypeError):
        gap.metadata["new"] = "value"  # type: ignore[index]


@pytest.mark.parametrize("length", [0, -1, True, 1.5])
def test_gap_requires_positive_integer_or_none(length: object) -> None:
    with pytest.raises(SequenceError) as error:
        Gap(length)  # type: ignore[arg-type]
    assert error.value.code == "INVALID_GAP_LENGTH"


def test_circular_sequences_do_not_gain_implicit_rotation_equality() -> None:
    first = DNASequence("ACGT", topology=Topology.CIRCULAR)
    rotated = DNASequence("GTAC", topology=Topology.CIRCULAR)

    assert first != rotated
    assert first.topology is Topology.CIRCULAR
