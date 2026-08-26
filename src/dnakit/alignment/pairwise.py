"""Needleman-Wunsch and Smith-Waterman alignment with bounded memory."""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import chain
from typing import Literal, TypeAlias, cast

from dnakit.alignment.results import AlignmentColumn, AlignmentResult
from dnakit.core import (
    Citation,
    DNARecord,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    Topology,
)
from dnakit.core._json import FrozenDict
from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

SequenceInput: TypeAlias = DNA | DNASequence | DNARecord
AlignmentMode: TypeAlias = Literal["global", "local", "semi_global"]

_STOP = 0
_DIAGONAL = 1
_DELETE = 2
_INSERT = 3
_FROM_MATCH = 1
_FROM_DELETE = 2
_FROM_INSERT = 3


@dataclass(frozen=True, slots=True)
class AlignmentConfig:
    """Scoring and resource limits for native linear-gap alignment."""

    mode: AlignmentMode = "global"
    match_score: float = 1.0
    mismatch_score: float = -1.0
    gap_score: float = -1.0
    gap_open_score: float | None = None
    gap_extend_score: float | None = None
    max_cells: int = 5_000_000

    def __post_init__(self) -> None:
        if self.mode not in ("global", "local", "semi_global"):
            raise ConfigurationError(
                "mode must be global, local, or semi_global.", code="INVALID_ALIGNMENT_MODE"
            )
        for name in ("match_score", "mismatch_score", "gap_score"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ConfigurationError(
                    f"{name} must be a finite number.", code="INVALID_ALIGNMENT_SCORE"
                )
            object.__setattr__(self, name, float(value))
        for name in ("gap_open_score", "gap_extend_score"):
            value = getattr(self, name)
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ConfigurationError(
                    f"{name} must be a finite number or None.", code="INVALID_ALIGNMENT_SCORE"
                )
            if value is not None:
                object.__setattr__(self, name, float(value))
        if (self.gap_open_score is None) != (self.gap_extend_score is None):
            raise ConfigurationError(
                "gap_open_score and gap_extend_score must be provided together.",
                code="INVALID_ALIGNMENT_SCORE",
            )
        if (
            isinstance(self.max_cells, bool)
            or not isinstance(self.max_cells, int)
            or self.max_cells <= 0
        ):
            raise ConfigurationError(
                "max_cells must be a positive integer.", code="INVALID_ALIGNMENT_CELL_LIMIT"
            )


def _sequence_and_id(value: SequenceInput, role: str) -> tuple[DNASequence, str | None]:
    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        sequence_id = None if record is None else record.id
    else:
        raise ConfigurationError(
            f"{role} must be DNASequence or DNARecord.",
            code="INVALID_ALIGNMENT_INPUT",
            context={"role": role, "type": type(value).__name__},
        )
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Pairwise alignment requires gap-free input; '-' is reserved for its output.",
            code="ALIGNMENT_INPUT_GAP_NOT_SUPPORTED",
            context={"role": role},
        )
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Pairwise alignment requires an explicitly linearized sequence.",
            code="ALIGNMENT_CIRCULAR_UNSUPPORTED",
            context={"role": role},
            hint="Rotate the sequence to the desired origin and convert it to linear topology.",
        )
    return sequence, sequence_id


def _choose(candidates: tuple[tuple[float, int], ...]) -> tuple[float, int]:
    """Choose the largest score and retain deterministic candidate order on ties."""

    best_score, best_pointer = candidates[0]
    for score, pointer in candidates[1:]:
        if score > best_score:
            best_score, best_pointer = score, pointer
    return best_score, best_pointer


def _alignment_provenance(config: AlignmentConfig) -> Provenance:
    primary = (
        Citation(
            "smith-waterman-1981",
            title="Identification of common molecular subsequences",
            doi="10.1016/0022-2836(81)90087-5",
        )
        if config.mode == "local"
        else Citation(
            "needleman-wunsch-1970",
            title="A general method applicable to the search for similarities",
            doi="10.1016/0022-2836(70)90057-4",
        )
    )
    citations: tuple[Citation, ...] = (primary,)
    if config.gap_open_score is not None:
        citations += (
            Citation(
                "gotoh-1982",
                title="An improved algorithm for matching biological sequences",
                doi="10.1016/0022-2836(82)90398-9",
            ),
        )
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.REIMPLEMENTATION,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.PUBLISHED_ALGORITHM,
            citations=citations,
        )
    )


def _fill(
    query: str,
    target: str,
    config: AlignmentConfig,
) -> tuple[list[list[float]], list[bytearray], tuple[int, int]]:
    rows = len(query) + 1
    columns = len(target) + 1
    cells = rows * columns
    if cells > config.max_cells:
        raise ConfigurationError(
            "Alignment dynamic-programming matrix exceeds max_cells.",
            code="ALIGNMENT_CELL_LIMIT",
            context={"dp_cells": cells, "max_cells": config.max_cells},
        )
    scores = [[0.0] * columns for _ in range(rows)]
    pointers = [bytearray(columns) for _ in range(rows)]
    if config.mode == "global":
        for row in range(1, rows):
            scores[row][0] = row * config.gap_score
            pointers[row][0] = _DELETE
        for column in range(1, columns):
            scores[0][column] = column * config.gap_score
            pointers[0][column] = _INSERT
    best = (0, 0)
    best_score = 0.0
    for row in range(1, rows):
        for column in range(1, columns):
            diagonal = scores[row - 1][column - 1] + (
                config.match_score
                if query[row - 1] == target[column - 1]
                else config.mismatch_score
            )
            deleted = scores[row - 1][column] + config.gap_score
            inserted = scores[row][column - 1] + config.gap_score
            candidates = (
                ((0.0, _STOP), (diagonal, _DIAGONAL), (deleted, _DELETE), (inserted, _INSERT))
                if config.mode == "local"
                else ((diagonal, _DIAGONAL), (deleted, _DELETE), (inserted, _INSERT))
            )
            score, pointer = _choose(candidates)
            scores[row][column] = score
            pointers[row][column] = pointer
            if config.mode == "local" and score > best_score:
                best_score = score
                best = (row, column)
    if config.mode == "semi_global":
        _, end_row, end_column = max(
            chain(
                ((scores[len(query)][column], len(query), column) for column in range(columns)),
                ((scores[row][len(target)], row, len(target)) for row in range(rows - 1)),
            ),
            key=lambda item: item[0],
        )
        return scores, pointers, (end_row, end_column)
    return scores, pointers, (best if config.mode == "local" else (len(query), len(target)))


def _fill_affine(
    query: str,
    target: str,
    config: AlignmentConfig,
) -> tuple[
    list[list[float]],
    tuple[list[bytearray], list[bytearray], list[bytearray]],
    tuple[int, int, int],
]:
    rows = len(query) + 1
    columns = len(target) + 1
    cells = rows * columns * 3
    if cells > config.max_cells:
        raise ConfigurationError(
            "Alignment affine dynamic-programming matrices exceed max_cells.",
            code="ALIGNMENT_CELL_LIMIT",
            context={"dp_cells": cells, "max_cells": config.max_cells},
        )
    gap_open = cast(float, config.gap_open_score)
    gap_extend = cast(float, config.gap_extend_score)
    negative = float("-inf")
    match = [[negative] * columns for _ in range(rows)]
    deleted = [[negative] * columns for _ in range(rows)]
    inserted = [[negative] * columns for _ in range(rows)]
    scores = [[negative] * columns for _ in range(rows)]
    match_pointers = [bytearray(columns) for _ in range(rows)]
    delete_pointers = [bytearray(columns) for _ in range(rows)]
    insert_pointers = [bytearray(columns) for _ in range(rows)]
    match[0][0] = scores[0][0] = 0.0
    if config.mode == "global":
        for row in range(1, rows):
            deleted[row][0] = scores[row][0] = gap_open + (row - 1) * gap_extend
            delete_pointers[row][0] = _FROM_MATCH if row == 1 else _FROM_DELETE
        for column in range(1, columns):
            inserted[0][column] = scores[0][column] = gap_open + (column - 1) * gap_extend
            insert_pointers[0][column] = _FROM_MATCH if column == 1 else _FROM_INSERT
    elif config.mode in ("local", "semi_global"):
        for row in range(rows):
            match[row][0] = scores[row][0] = 0.0
        for column in range(columns):
            match[0][column] = scores[0][column] = 0.0
    for row in range(1, rows):
        for column in range(1, columns):
            substitution = (
                config.match_score
                if query[row - 1] == target[column - 1]
                else config.mismatch_score
            )
            match_candidates: tuple[tuple[float, int], ...] = (
                (match[row - 1][column - 1] + substitution, _FROM_MATCH),
                (deleted[row - 1][column - 1] + substitution, _FROM_DELETE),
                (inserted[row - 1][column - 1] + substitution, _FROM_INSERT),
            )
            if config.mode == "local":
                match_candidates = ((0.0, _STOP), *match_candidates)
            match_previous = _choose(match_candidates)
            match[row][column] = match_previous[0]
            match_pointers[row][column] = match_previous[1]
            deleted[row][column], delete_pointers[row][column] = _choose(
                (
                    (match[row - 1][column] + gap_open, _FROM_MATCH),
                    (deleted[row - 1][column] + gap_extend, _FROM_DELETE),
                    (inserted[row - 1][column] + gap_open, _FROM_INSERT),
                )
            )
            inserted[row][column], insert_pointers[row][column] = _choose(
                (
                    (match[row][column - 1] + gap_open, _FROM_MATCH),
                    (deleted[row][column - 1] + gap_open, _FROM_DELETE),
                    (inserted[row][column - 1] + gap_extend, _FROM_INSERT),
                )
            )
            scores[row][column] = max(
                match[row][column], deleted[row][column], inserted[row][column]
            )
            if config.mode == "local":
                scores[row][column] = max(0.0, scores[row][column])
    pointers = (match_pointers, delete_pointers, insert_pointers)
    states = (
        (_FROM_MATCH, match),
        (_FROM_DELETE, deleted),
        (_FROM_INSERT, inserted),
    )
    if config.mode == "local":
        _, row, column, state = max(
            (
                (matrix[row][column], row, column, state)
                for row in range(rows)
                for column in range(columns)
                for state, matrix in states
            ),
            key=lambda item: item[0],
        )
        return scores, pointers, (row, column, state)
    if config.mode == "semi_global":
        _, row, column, state = max(
            chain(
                (
                    (matrix[-1][column], len(query), column, state)
                    for column in range(columns)
                    for state, matrix in states
                ),
                (
                    (matrix[row][-1], row, len(target), state)
                    for row in range(rows - 1)
                    for state, matrix in states
                ),
            ),
            key=lambda item: item[0],
        )
        return scores, pointers, (row, column, state)
    state = _choose(tuple((matrix[-1][-1], state) for state, matrix in states))[1]
    return scores, pointers, (len(query), len(target), state)


def _trace_affine(
    query: str,
    target: str,
    pointers: tuple[list[bytearray], list[bytearray], list[bytearray]],
    endpoint: tuple[int, int, int],
    mode: AlignmentMode,
) -> tuple[list[AlignmentColumn], list[str], list[str], int, int]:
    row, column, state = endpoint
    match_pointers, delete_pointers, insert_pointers = pointers
    reverse_columns: list[AlignmentColumn] = []
    reverse_query: list[str] = []
    reverse_target: list[str] = []
    while row > 0 or column > 0:
        if mode == "semi_global" and (row == 0 or column == 0):
            break
        if state == _FROM_MATCH:
            previous = match_pointers[row][column]
            if previous == _STOP:
                break
            query_symbol, target_symbol = query[row - 1], target[column - 1]
            operation: Literal["match", "mismatch", "insertion", "deletion"] = (
                "match" if query_symbol == target_symbol else "mismatch"
            )
            reverse_columns.append(
                AlignmentColumn(query_symbol, target_symbol, row - 1, column - 1, operation)
            )
            reverse_query.append(query_symbol)
            reverse_target.append(target_symbol)
            row -= 1
            column -= 1
        elif state == _FROM_DELETE:
            previous = delete_pointers[row][column]
            query_symbol = query[row - 1]
            reverse_columns.append(AlignmentColumn(query_symbol, "-", row - 1, None, "deletion"))
            reverse_query.append(query_symbol)
            reverse_target.append("-")
            row -= 1
        elif state == _FROM_INSERT:
            previous = insert_pointers[row][column]
            target_symbol = target[column - 1]
            reverse_columns.append(
                AlignmentColumn("-", target_symbol, None, column - 1, "insertion")
            )
            reverse_query.append("-")
            reverse_target.append(target_symbol)
            column -= 1
        else:
            raise RuntimeError("Invalid affine-alignment traceback state.")
        state = previous
    return reverse_columns, reverse_query, reverse_target, row, column


def align_pairwise(
    query: SequenceInput,
    target: SequenceInput,
    *,
    config: AlignmentConfig | None = None,
) -> AlignmentResult:
    """Align two sequences with literal-IUPAC equality and deterministic traceback."""

    resolved = AlignmentConfig() if config is None else config
    if not isinstance(resolved, AlignmentConfig):
        raise ConfigurationError(
            "config must be AlignmentConfig or None.", code="INVALID_ALIGNMENT_CONFIG"
        )
    query_sequence, query_id = _sequence_and_id(query, "query")
    target_sequence, target_id = _sequence_and_id(target, "target")
    query_text = query_sequence.symbols
    target_text = target_sequence.symbols
    if resolved.gap_open_score is None:
        scores, linear_pointers, (row, column) = _fill(query_text, target_text, resolved)
        query_end, target_end = row, column
        reverse_columns: list[AlignmentColumn] = []
        reverse_query: list[str] = []
        reverse_target: list[str] = []
        while row > 0 or column > 0:
            if resolved.mode == "semi_global" and (row == 0 or column == 0):
                break
            pointer = linear_pointers[row][column]
            if pointer == _STOP:
                break
            if pointer == _DIAGONAL:
                query_symbol, target_symbol = query_text[row - 1], target_text[column - 1]
                operation: Literal["match", "mismatch", "insertion", "deletion"] = (
                    "match" if query_symbol == target_symbol else "mismatch"
                )
                reverse_columns.append(
                    AlignmentColumn(query_symbol, target_symbol, row - 1, column - 1, operation)
                )
                reverse_query.append(query_symbol)
                reverse_target.append(target_symbol)
                row -= 1
                column -= 1
            elif pointer == _DELETE:
                query_symbol = query_text[row - 1]
                reverse_columns.append(
                    AlignmentColumn(query_symbol, "-", row - 1, None, "deletion")
                )
                reverse_query.append(query_symbol)
                reverse_target.append("-")
                row -= 1
            elif pointer == _INSERT:
                target_symbol = target_text[column - 1]
                reverse_columns.append(
                    AlignmentColumn("-", target_symbol, None, column - 1, "insertion")
                )
                reverse_query.append("-")
                reverse_target.append(target_symbol)
                column -= 1
            else:
                raise RuntimeError("Invalid alignment traceback pointer.")
        score = scores[query_end][target_end]
    else:
        scores, affine_pointers, endpoint = _fill_affine(query_text, target_text, resolved)
        query_end, target_end, _ = endpoint
        reverse_columns, reverse_query, reverse_target, row, column = _trace_affine(
            query_text,
            target_text,
            affine_pointers,
            endpoint,
            resolved.mode,
        )
        score = scores[query_end][target_end]
    columns = tuple(reversed(reverse_columns))
    matches = sum(item.operation == "match" for item in columns)
    mismatches = sum(item.operation == "mismatch" for item in columns)
    insertions = sum(item.operation == "insertion" for item in columns)
    deletions = sum(item.operation == "deletion" for item in columns)
    comparable = matches + mismatches + insertions + deletions
    query_consumed = matches + mismatches + deletions
    target_consumed = matches + mismatches + insertions
    return AlignmentResult(
        name="pairwise_alignment",
        method=resolved.mode,
        algorithm_version=(
            "affine-gap-dp-v1" if resolved.gap_open_score is not None else "linear-gap-dp-v1"
        ),
        score=score,
        aligned_query="".join(reversed(reverse_query)),
        aligned_target="".join(reversed(reverse_target)),
        query_id=query_id,
        target_id=target_id,
        query_start=row,
        query_end=query_end,
        target_start=column,
        target_end=target_end,
        matches=matches,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        identity=None if comparable == 0 else matches / comparable,
        query_coverage=1.0 if not query_text else query_consumed / len(query_text),
        target_coverage=1.0 if not target_text else target_consumed / len(target_text),
        columns=columns,
        parameters=FrozenDict(
            {
                "match_score": resolved.match_score,
                "mismatch_score": resolved.mismatch_score,
                "gap_score": resolved.gap_score,
                "gap_open_score": resolved.gap_open_score,
                "gap_extend_score": resolved.gap_extend_score,
                "max_cells": resolved.max_cells,
                "dp_cells": (len(query_text) + 1)
                * (len(target_text) + 1)
                * (3 if resolved.gap_open_score is not None else 1),
                "iupac_matching": "literal",
                "tie_break": "diagonal-delete-insert",
                "end_gap_policy": {
                    "global": "penalized",
                    "local": "excluded-local-flanks",
                    "semi_global": "free-query-and-target-ends",
                }[resolved.mode],
            }
        ),
        provenance=_alignment_provenance(resolved),
    )


__all__ = ["AlignmentConfig", "align_pairwise"]
