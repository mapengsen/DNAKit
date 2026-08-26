"""Tests for deterministic circular-origin sequence operations."""

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap, GapKind, Topology
from dnakit.exceptions import (
    ConfigurationError,
    CoordinateError,
    SequenceError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)
from dnakit.ops import canonical_origin, circular_subsequence, rotate, subsequence


def _circular(parts: str | list[str | Gap]) -> DNASequence:
    return DNASequence(parts, topology=Topology.CIRCULAR)


def test_rotate_uses_left_offset_modulo_span_without_mutating_input() -> None:
    source = _circular("GATTACA")

    result = rotate(source, 9)
    negative = rotate(source, -1)

    assert result.sequence.symbols == "TTACAGA"
    assert result.requested_offset == 9
    assert result.effective_offset == 2
    assert result.sequence_span == 7
    assert result.rule == "left_rotation_offset_modulo_coordinate_span"
    assert result.parameters["direction"] == "left"
    assert negative.sequence.symbols == "AGATTAC"
    assert negative.effective_offset == 6
    assert source.symbols == "GATTACA"
    assert result.to_dict()["effective_offset"] == 2


def test_rotate_preserves_known_gap_metadata_and_rejects_gap_interior() -> None:
    gap = Gap(
        3,
        GapKind.SCAFFOLD,
        crossable=False,
        evidence=("assembly",),
        metadata={"source": "fixture"},
    )
    source = _circular(["AA", gap, "CC"])

    at_boundary = rotate(source, 2)

    assert at_boundary.sequence.parts == (gap, "CCAA")
    assert at_boundary.sequence.parts[0] is gap
    with pytest.raises(UnsupportedGapOperationError) as error:
        rotate(source, 3)
    assert error.value.code == "ROTATION_ORIGIN_INSIDE_GAP"


def test_rotate_rejects_linear_unknown_length_empty_and_bad_offsets() -> None:
    with pytest.raises(ConfigurationError) as linear:
        rotate(DNASequence("AC"), 1)
    assert linear.value.code == "CIRCULAR_TOPOLOGY_REQUIRED"
    with pytest.raises(UnknownLengthError):
        rotate(_circular(["A", Gap(None), "C"]), 0)
    with pytest.raises(ConfigurationError) as bad_offset:
        rotate(_circular("AC"), True)
    assert bad_offset.value.code == "INVALID_ROTATION_OFFSET"
    with pytest.raises(SequenceError) as empty:
        _circular("")
    assert empty.value.code == "EMPTY_CIRCULAR_SEQUENCE"


def test_canonical_origin_is_forward_lexicographic_and_tie_stable() -> None:
    source = _circular("TACACA")

    result = canonical_origin(source)
    periodic = canonical_origin(_circular("ATAT"))
    iupac = canonical_origin(
        DNASequence("NARA", alphabet=DNAAlphabet.IUPAC, topology=Topology.CIRCULAR)
    )

    rotations = [source.symbols[index:] + source.symbols[:index] for index in range(len(source))]
    assert result.sequence.symbols == min(rotations) == "ACACAT"
    assert result.effective_offset == 1
    assert periodic.sequence.symbols == "ATAT"
    assert periodic.effective_offset == 0
    assert iupac.sequence.symbols == "ANAR"
    assert result.parameters["algorithm"] == "booth"
    assert result.parameters["orientation"] == "forward_only"


def test_canonical_origin_explicitly_rejects_gaps() -> None:
    with pytest.raises(UnsupportedGapOperationError) as error:
        canonical_origin(_circular(["AC", Gap(1), "GT"]))
    assert error.value.code == "CANONICAL_ORIGIN_GAPPED_SEQUENCE"


def test_subsequence_wraps_origin_only_for_circular_topology() -> None:
    source = _circular("AACCGG")

    assert subsequence(source, 4, 2).symbols == "GGAA"
    assert circular_subsequence(source, 4, 2).symbols == "GGAA"
    assert subsequence(source, 2, 2) == DNASequence("")
    full = subsequence(source, 0, 6)
    assert full.symbols == source.symbols
    assert full.topology is Topology.LINEAR

    with pytest.raises(CoordinateError) as linear:
        subsequence(DNASequence("AACCGG"), 4, 2)
    assert linear.value.code == "LINEAR_SUBSEQUENCE_REVERSED"
    with pytest.raises(ConfigurationError) as circular_required:
        circular_subsequence(DNASequence("AC"), 1, 0)
    assert circular_required.value.code == "CIRCULAR_TOPOLOGY_REQUIRED"


def test_wrapped_subsequence_obeys_gap_and_unknown_length_policies() -> None:
    source = _circular(["AA", Gap(2, metadata={"x": 1}), "CC"])

    with pytest.raises(UnsupportedGapOperationError):
        subsequence(source, 5, 3)
    allowed = subsequence(source, 5, 3, allow_gaps=True)
    assert allowed.parts == ("CAA", Gap(1, metadata={"x": 1}))

    with pytest.raises(UnknownLengthError):
        subsequence(_circular(["AA", Gap(None), "CC"]), 1, 0, allow_gaps=True)
