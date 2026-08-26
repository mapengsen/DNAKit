"""Bounded approximate substring matching by dynamic programming."""

from __future__ import annotations

import math
from collections.abc import Iterable

from dnakit.core import DNA, DNARecord, DNASequence, Strand, Topology
from dnakit.core._json import FrozenDict
from dnakit.core.facade import resolve_single_dna
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

from ._shared import SequenceInput, materialize_targets, validate_bool
from .results import ApproximateMatch, ApproximateSearchResult

DEFAULT_MAX_APPROXIMATE_CELLS = 5_000_000
DEFAULT_MAX_APPROXIMATE_MATCHES = 100_000
DEFAULT_MAX_APPROXIMATE_TARGETS = 100_000


def _finite_non_negative(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
    ):
        raise ConfigurationError(
            f"{name} must be a finite non-negative number.",
            code="INVALID_APPROXIMATE_MATCH_COST",
            context={"field": name, "value": value},
        )
    return float(value)


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code="INVALID_APPROXIMATE_MATCH_LIMIT",
        )
    return value


def _sequence_and_id_without_gap_policy(
    value: SequenceInput, *, role: str
) -> tuple[DNASequence, str | None]:
    """Resolve approximate-search inputs before applying its domain-specific Gap error."""

    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return sequence, None if record is None else record.id
    raise ConfigurationError(
        f"{role} must be DNASequence or DNARecord.",
        code="INVALID_SIMILARITY_SEQUENCE",
        context={"role": role, "input_type": type(value).__name__},
    )


def _scan_one(
    query: str,
    target: str,
    *,
    max_distance: float,
    substitution_cost: float,
    insertion_cost: float,
    deletion_cost: float,
    target_index: int,
    target_id: str | None,
    strand: Strand,
    remaining_matches: int,
) -> tuple[ApproximateMatch, ...]:
    columns = len(target) + 1
    previous = [0.0] * columns
    starts_previous = list(range(columns))
    for query_index, query_symbol in enumerate(query, start=1):
        current = [query_index * deletion_cost] + [0.0] * len(target)
        starts_current = [0] * columns
        for target_index_1, target_symbol in enumerate(target, start=1):
            candidates = (
                (
                    previous[target_index_1 - 1]
                    + (0.0 if query_symbol == target_symbol else substitution_cost),
                    starts_previous[target_index_1 - 1],
                    0,
                ),
                (
                    previous[target_index_1] + deletion_cost,
                    starts_previous[target_index_1],
                    1,
                ),
                (
                    current[target_index_1 - 1] + insertion_cost,
                    starts_current[target_index_1 - 1],
                    2,
                ),
            )
            score, start, _ = min(candidates, key=lambda item: (item[0], item[2], item[1]))
            current[target_index_1] = score
            starts_current[target_index_1] = start
        previous, current = current, previous
        starts_previous, starts_current = starts_current, starts_previous
    matches: list[ApproximateMatch] = []
    seen: set[tuple[int, int, float]] = set()
    for end, distance in enumerate(previous[1:], start=1):
        start = starts_previous[end]
        key = (start, end, distance)
        if distance <= max_distance and key not in seen:
            if len(matches) >= remaining_matches:
                raise ConfigurationError(
                    "Approximate search exceeds max_matches.",
                    code="APPROXIMATE_MATCH_LIMIT_EXCEEDED",
                    context={"match_count_lower_bound": remaining_matches + 1},
                )
            seen.add(key)
            matches.append(
                ApproximateMatch(
                    target_index=target_index,
                    target_id=target_id,
                    start=start,
                    end=end,
                    strand=strand,
                    distance=distance,
                )
            )
    return tuple(matches)


def approximate_search(
    query: SequenceInput,
    targets: SequenceInput | Iterable[SequenceInput],
    *,
    max_distance: int | float,
    substitution_cost: int | float = 1.0,
    insertion_cost: int | float = 1.0,
    deletion_cost: int | float = 1.0,
    reverse_complement: bool = False,
    max_targets: int = DEFAULT_MAX_APPROXIMATE_TARGETS,
    max_matches: int = DEFAULT_MAX_APPROXIMATE_MATCHES,
    max_cells: int = DEFAULT_MAX_APPROXIMATE_CELLS,
) -> ApproximateSearchResult:
    """Find bounded-cost query alignments ending at each target coordinate."""

    maximum = _finite_non_negative(max_distance, "max_distance")
    substitution = _finite_non_negative(substitution_cost, "substitution_cost")
    insertion = _finite_non_negative(insertion_cost, "insertion_cost")
    deletion = _finite_non_negative(deletion_cost, "deletion_cost")
    if substitution == insertion == deletion == 0.0:
        raise ConfigurationError(
            "At least one edit cost must be positive.", code="DEGENERATE_APPROXIMATE_COSTS"
        )
    validate_bool(reverse_complement, "reverse_complement")
    for name, value in (
        ("max_targets", max_targets),
        ("max_matches", max_matches),
        ("max_cells", max_cells),
    ):
        _positive_int(value, name)
    query_sequence, query_id = _sequence_and_id_without_gap_policy(query, role="query")
    if query_sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Approximate query matching requires a linear query.",
            code="APPROXIMATE_CIRCULAR_UNSUPPORTED",
        )
    if query_sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Approximate matching cannot silently omit explicit query Gap objects.",
            code="APPROXIMATE_GAP_NOT_ALLOWED",
            context={"role": "query"},
        )
    if query_sequence.symbol_length == 0:
        raise ConfigurationError(
            "Approximate matching does not accept an empty query.",
            code="EMPTY_APPROXIMATE_QUERY",
        )
    materialized = materialize_targets(targets, max_targets=max_targets)
    resolved_targets: list[tuple[DNASequence, str | None]] = []
    total_cells = 0
    for target in materialized:
        target_sequence, target_id = _sequence_and_id_without_gap_policy(target, role="target")
        if target_sequence.topology is Topology.CIRCULAR:
            raise ConfigurationError(
                "Approximate matching currently requires linear targets.",
                code="APPROXIMATE_CIRCULAR_UNSUPPORTED",
            )
        if target_sequence.is_gapped:
            raise UnsupportedGapOperationError(
                "Approximate matching cannot silently omit explicit target Gap objects.",
                code="APPROXIMATE_GAP_NOT_ALLOWED",
                context={"role": "target"},
            )
        cells = (query_sequence.symbol_length + 1) * (target_sequence.symbol_length + 1)
        total_cells += cells * (2 if reverse_complement else 1)
        if total_cells > max_cells:
            raise ConfigurationError(
                "Approximate search exceeds max_cells.",
                code="APPROXIMATE_MATCH_CELL_LIMIT",
                context={"required_cells_lower_bound": total_cells, "max_cells": max_cells},
            )
        resolved_targets.append((target_sequence, target_id))
    strands = [(query_sequence.symbols, Strand.FORWARD)]
    if reverse_complement:
        reverse = query_sequence.reverse_complement().symbols
        if reverse != query_sequence.symbols:
            strands.append((reverse, Strand.REVERSE))
    matches: list[ApproximateMatch] = []
    for target_index, (target_sequence, target_id) in enumerate(resolved_targets):
        for query_symbols, strand in strands:
            for match in _scan_one(
                query_symbols,
                target_sequence.symbols,
                max_distance=maximum,
                substitution_cost=substitution,
                insertion_cost=insertion,
                deletion_cost=deletion,
                target_index=target_index,
                target_id=target_id,
                strand=strand,
                remaining_matches=max_matches - len(matches),
            ):
                matches.append(match)
    ordered = tuple(
        sorted(
            matches, key=lambda item: (item.target_index, item.start, item.end, item.strand.value)
        )
    )
    return ApproximateSearchResult(
        name="approximate_search",
        method="semi-global-weighted-edit-distance",
        algorithm_version="dnakit-approximate-search-v1",
        query_id=query_id,
        query_length=query_sequence.symbol_length,
        target_count=len(resolved_targets),
        matches=ordered,
        max_distance=maximum,
        reverse_complement=reverse_complement,
        max_targets=max_targets,
        max_matches=max_matches,
        max_cells=max_cells,
        dp_cells=total_cells,
        parameters=FrozenDict(
            {
                "substitution_cost": substitution,
                "insertion_cost": insertion,
                "deletion_cost": deletion,
                "iupac_matching": "literal",
                "coordinate_system": "0-based-half-open",
            }
        ),
    )


__all__ = [
    "DEFAULT_MAX_APPROXIMATE_CELLS",
    "DEFAULT_MAX_APPROXIMATE_MATCHES",
    "DEFAULT_MAX_APPROXIMATE_TARGETS",
    "approximate_search",
]
