"""Literal Hamming and weighted Levenshtein distances."""

from __future__ import annotations

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.similarity._shared import (
    SequenceInput,
    sequence_and_id,
    validate_bool,
    validate_non_negative_number,
    validate_positive_int,
)
from dnakit.similarity.results import DistanceResult, EditStep, Mismatch

DEFAULT_MAX_EDIT_CELLS = 5_000_000

_DIAGONAL = 0
_DELETE = 1
_INSERT = 2


def _validate_max_distance(value: int | float | None) -> float | None:
    if value is None:
        return None
    return validate_non_negative_number(value, "max_distance")


def hamming_distance(
    left: SequenceInput,
    right: SequenceInput,
    *,
    max_distance: int | float | None = None,
) -> DistanceResult:
    """Count literal symbol differences between two equal-length sequences.

    Every IUPAC symbol compares as itself: for example, ``N == N`` and
    ``N != A``.  Explicit gaps are rejected rather than ignored.
    """

    left_sequence, left_id = sequence_and_id(left, role="left")
    right_sequence, right_id = sequence_and_id(right, role="right")
    if left_sequence.symbol_length != right_sequence.symbol_length:
        raise ConfigurationError(
            "Hamming distance requires equal symbol lengths.",
            code="HAMMING_LENGTH_MISMATCH",
            context={
                "left_length": left_sequence.symbol_length,
                "right_length": right_sequence.symbol_length,
            },
        )
    maximum = _validate_max_distance(max_distance)
    mismatches = tuple(
        Mismatch(position, left_symbol, right_symbol)
        for position, (left_symbol, right_symbol) in enumerate(
            zip(left_sequence.symbols, right_sequence.symbols, strict=True)
        )
        if left_symbol != right_symbol
    )
    distance = float(len(mismatches))
    return DistanceResult(
        name="hamming_distance",
        method="hamming",
        left_id=left_id,
        right_id=right_id,
        left_length=left_sequence.symbol_length,
        right_length=right_sequence.symbol_length,
        distance=distance,
        mismatches=mismatches,
        edit_path=None,
        costs=FrozenDict({"substitution": 1.0}),
        max_distance=maximum,
        exceeded_max_distance=maximum is not None and distance > maximum,
    )


def _levenshtein_rows(
    left: str,
    right: str,
    *,
    substitution_cost: float,
    insertion_cost: float,
    deletion_cost: float,
) -> float:
    previous = [column * insertion_cost for column in range(len(right) + 1)]
    for row, left_symbol in enumerate(left, start=1):
        current = [row * deletion_cost]
        for column, right_symbol in enumerate(right, start=1):
            diagonal_cost = 0.0 if left_symbol == right_symbol else substitution_cost
            current.append(
                min(
                    previous[column - 1] + diagonal_cost,
                    previous[column] + deletion_cost,
                    current[column - 1] + insertion_cost,
                )
            )
        previous = current
    return previous[-1]


def _levenshtein_with_path(
    left: str,
    right: str,
    *,
    substitution_cost: float,
    insertion_cost: float,
    deletion_cost: float,
) -> tuple[float, tuple[EditStep, ...]]:
    rows = len(left) + 1
    columns = len(right) + 1
    table = [[0.0] * columns for _ in range(rows)]
    predecessors = [bytearray(columns) for _ in range(rows)]
    for row in range(1, rows):
        table[row][0] = row * deletion_cost
        predecessors[row][0] = _DELETE
    for column in range(1, columns):
        table[0][column] = column * insertion_cost
        predecessors[0][column] = _INSERT
    for row in range(1, rows):
        for column in range(1, columns):
            diagonal_cost = 0.0 if left[row - 1] == right[column - 1] else substitution_cost
            best_cost = table[row - 1][column - 1] + diagonal_cost
            predecessor = _DIAGONAL
            deletion_candidate = table[row - 1][column] + deletion_cost
            if deletion_candidate < best_cost:
                best_cost = deletion_candidate
                predecessor = _DELETE
            insertion_candidate = table[row][column - 1] + insertion_cost
            if insertion_candidate < best_cost:
                best_cost = insertion_candidate
                predecessor = _INSERT
            table[row][column] = best_cost
            predecessors[row][column] = predecessor

    steps: list[EditStep] = []
    row, column = len(left), len(right)
    while row or column:
        predecessor = predecessors[row][column]
        if predecessor == _DIAGONAL and row and column:
            is_match = left[row - 1] == right[column - 1]
            steps.append(
                EditStep(
                    "match" if is_match else "substitute",
                    row - 1,
                    row,
                    column - 1,
                    column,
                    left[row - 1],
                    right[column - 1],
                )
            )
            row -= 1
            column -= 1
            continue
        if predecessor == _DELETE and row:
            steps.append(
                EditStep(
                    "delete",
                    row - 1,
                    row,
                    column,
                    column,
                    left[row - 1],
                    None,
                )
            )
            row -= 1
            continue
        if predecessor == _INSERT and column:
            steps.append(
                EditStep(
                    "insert",
                    row,
                    row,
                    column - 1,
                    column,
                    None,
                    right[column - 1],
                )
            )
            column -= 1
            continue
        raise AssertionError("Levenshtein traceback could not resolve a predecessor.")
    steps.reverse()
    return table[-1][-1], tuple(steps)


def edit_distance(
    left: SequenceInput,
    right: SequenceInput,
    *,
    substitution_cost: int | float = 1,
    insertion_cost: int | float = 1,
    deletion_cost: int | float = 1,
    max_distance: int | float | None = None,
    return_path: bool = False,
    max_cells: int = DEFAULT_MAX_EDIT_CELLS,
) -> DistanceResult:
    """Compute weighted literal Levenshtein distance by dynamic programming.

    ``max_distance`` is an auditable threshold, not an approximation: DNAKit
    still returns the exact distance and records whether it exceeds the limit.
    Without a requested path, memory is linear in the right sequence length.
    Both modes refuse more than ``max_cells`` dynamic-programming cells before
    allocating a table or starting the recurrence.
    """

    left_sequence, left_id = sequence_and_id(left, role="left")
    right_sequence, right_id = sequence_and_id(right, role="right")
    substitution = validate_non_negative_number(substitution_cost, "substitution_cost")
    insertion = validate_non_negative_number(insertion_cost, "insertion_cost")
    deletion = validate_non_negative_number(deletion_cost, "deletion_cost")
    maximum = _validate_max_distance(max_distance)
    validate_bool(return_path, "return_path")
    validate_positive_int(max_cells, "max_cells")
    dp_cells = (left_sequence.symbol_length + 1) * (right_sequence.symbol_length + 1)
    if dp_cells > max_cells:
        raise ConfigurationError(
            "Edit-distance dynamic-programming input exceeds max_cells.",
            code="EDIT_DISTANCE_CELL_LIMIT",
            context={
                "left_length": left_sequence.symbol_length,
                "right_length": right_sequence.symbol_length,
                "dp_cells": dp_cells,
                "max_cells": max_cells,
                "return_path": return_path,
            },
            hint="Raise max_cells deliberately or compare shorter sequences.",
        )
    if return_path:
        distance, path = _levenshtein_with_path(
            left_sequence.symbols,
            right_sequence.symbols,
            substitution_cost=substitution,
            insertion_cost=insertion,
            deletion_cost=deletion,
        )
    else:
        distance = _levenshtein_rows(
            left_sequence.symbols,
            right_sequence.symbols,
            substitution_cost=substitution,
            insertion_cost=insertion,
            deletion_cost=deletion,
        )
        path = None
    return DistanceResult(
        name="edit_distance",
        method="levenshtein",
        left_id=left_id,
        right_id=right_id,
        left_length=left_sequence.symbol_length,
        right_length=right_sequence.symbol_length,
        distance=distance,
        mismatches=(),
        edit_path=path,
        costs=FrozenDict(
            {
                "substitution": substitution,
                "insertion": insertion,
                "deletion": deletion,
            }
        ),
        max_distance=maximum,
        exceeded_max_distance=maximum is not None and distance > maximum,
        max_cells=max_cells,
        dp_cells=dp_cells,
    )


levenshtein_distance = edit_distance


__all__ = [
    "DEFAULT_MAX_EDIT_CELLS",
    "edit_distance",
    "hamming_distance",
    "levenshtein_distance",
]
