"""Boundary and traceback tests for native sequence distances."""

import math

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.similarity import DistanceResult, edit_distance, hamming_distance


def _numeric_cost(value: object) -> float:
    assert isinstance(value, (int, float)) and not isinstance(value, bool)
    return float(value)


def _edit_path_cost(result: DistanceResult) -> float:
    assert result.edit_path is not None
    operation_costs = {
        "match": 0.0,
        "substitute": _numeric_cost(result.costs["substitution"]),
        "insert": _numeric_cost(result.costs["insertion"]),
        "delete": _numeric_cost(result.costs["deletion"]),
    }
    return math.fsum(operation_costs[step.operation] for step in result.edit_path)


def test_hamming_distance_reports_literal_iupac_mismatch_positions() -> None:
    left = DNARecord(DNASequence("ARYN", alphabet=DNAAlphabet.IUPAC), "left")
    right = DNARecord(DNASequence("ARYA", alphabet=DNAAlphabet.IUPAC), "right")

    result = hamming_distance(left, right, max_distance=0)

    assert result.distance == 1
    assert result.left_id == "left"
    assert result.right_id == "right"
    assert result.exceeded_max_distance
    assert [(item.position, item.left_symbol, item.right_symbol) for item in result.mismatches] == [
        (3, "N", "A")
    ]
    assert hamming_distance(DNASequence(""), DNASequence("")).distance == 0


def test_hamming_requires_equal_lengths_and_gap_free_sequences() -> None:
    with pytest.raises(ConfigurationError) as length_error:
        hamming_distance(DNASequence("A"), DNASequence("AA"))
    assert length_error.value.code == "HAMMING_LENGTH_MISMATCH"
    with pytest.raises(UnsupportedGapOperationError):
        hamming_distance(DNASequence(["A", Gap(1)]), DNASequence("AA"))


def test_levenshtein_handles_insert_delete_substitute_and_empty_sequence() -> None:
    assert edit_distance(DNASequence("ACGT"), DNASequence("AGT")).distance == 1
    assert edit_distance(DNASequence("AGT"), DNASequence("ACGT")).distance == 1
    assert edit_distance(DNASequence("ACGT"), DNASequence("ATGT")).distance == 1
    assert edit_distance(DNASequence(""), DNASequence("ACG")).distance == 3
    assert edit_distance(DNASequence("ACG"), DNASequence("")).distance == 3


def test_weighted_levenshtein_returns_exact_threshold_status() -> None:
    result = edit_distance(
        DNASequence("A"),
        DNASequence("G"),
        substitution_cost=3,
        insertion_cost=1,
        deletion_cost=1,
        max_distance=1,
    )

    assert result.distance == 2
    assert result.exceeded_max_distance
    assert result.costs == {"substitution": 3.0, "insertion": 1.0, "deletion": 1.0}


def test_levenshtein_traceback_is_ordered_and_coordinate_explicit() -> None:
    result = edit_distance(DNASequence("ACG"), DNASequence("ATG"), return_path=True)

    assert result.edit_path is not None
    assert [step.operation for step in result.edit_path] == ["match", "substitute", "match"]
    substitute = result.edit_path[1]
    assert (substitute.left_start, substitute.left_end) == (1, 2)
    assert (substitute.right_start, substitute.right_end) == (1, 2)
    assert (substitute.left_symbol, substitute.right_symbol) == ("C", "T")
    assert _edit_path_cost(result) == result.distance


def test_levenshtein_traceback_never_selects_a_nearby_non_optimal_edge() -> None:
    result = edit_distance(
        DNASequence("A"),
        DNASequence("G"),
        substitution_cost=1.0000000001,
        insertion_cost=0.5,
        deletion_cost=0.5,
        return_path=True,
    )

    assert result.distance == 1.0
    assert result.edit_path is not None
    assert [step.operation for step in result.edit_path] == ["insert", "delete"]
    assert _edit_path_cost(result) == result.distance


@pytest.mark.parametrize("return_path", [False, True])
def test_edit_distance_rejects_dp_cell_limit_before_running(
    return_path: bool,
) -> None:
    with pytest.raises(ConfigurationError) as limit_error:
        edit_distance(
            DNASequence("ACGT"),
            DNASequence("ACGT"),
            return_path=return_path,
            max_cells=24,
        )

    assert limit_error.value.code == "EDIT_DISTANCE_CELL_LIMIT"
    assert limit_error.value.context == {
        "left_length": 4,
        "right_length": 4,
        "dp_cells": 25,
        "max_cells": 24,
        "return_path": return_path,
    }


def test_edit_distance_records_dp_cell_budget_in_result() -> None:
    result = edit_distance(DNASequence("ACG"), DNASequence("AT"), max_cells=20)

    assert result.dp_cells == 12
    assert result.max_cells == 20
    assert result.to_dict()["dp_cells"] == 12


@pytest.mark.parametrize("max_cells", [0, -1, True, 1.5])
def test_edit_distance_requires_positive_integer_max_cells(max_cells: object) -> None:
    with pytest.raises(ConfigurationError):
        edit_distance(
            DNASequence("A"),
            DNASequence("A"),
            max_cells=max_cells,  # type: ignore[arg-type]
        )


def test_edit_distance_rejects_gap_and_invalid_costs_or_runtime_flags() -> None:
    with pytest.raises(UnsupportedGapOperationError):
        edit_distance(DNASequence(["A", Gap(None), "C"]), DNASequence("AC"))
    with pytest.raises(ConfigurationError):
        edit_distance(DNASequence("A"), DNASequence("G"), substitution_cost=-1)
    with pytest.raises(ConfigurationError):
        edit_distance(DNASequence("A"), DNASequence("G"), return_path=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        hamming_distance(DNASequence("A"), DNASequence("A"), max_distance=True)
