"""IUPAC-compatible and approximate deduplication."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

from dnakit._parallel import ordered_map
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core import DNARecord, DNASet, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

from ._advanced_shared import ensure_pair_limit, materialize_limited, quality_score
from .clustering import cluster_sequences
from .config import ClusterConfig, IUPACDeduplicationConfig
from .results import (
    ClusteringResult,
    IUPACDeduplicationResult,
    IUPACDuplicateGroup,
    IUPACPairRelation,
)

_IUPAC = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def _compatible(left: str, right: str) -> bool:
    if left == right:
        return True
    return len(left) == len(right) and all(
        bool(_IUPAC[left_symbol] & _IUPAC[right_symbol])
        for left_symbol, right_symbol in zip(left, right, strict=True)
    )


def _representative(
    members: tuple[int, ...],
    records: tuple[DNARecord, ...],
    policy: str,
) -> int:
    if policy == "first":
        return members[0]
    if policy == "shortest":
        return min(members, key=lambda index: (records[index].sequence.symbol_length, index))
    if policy == "longest":
        return max(members, key=lambda index: (records[index].sequence.symbol_length, -index))
    return max(members, key=lambda index: (quality_score(records[index]), -index))


def _compatibility_job(payload: tuple[int, int, str, str]) -> tuple[int, int, bool, bool]:
    left, right, first, second = payload
    return left, right, first == second, _compatible(first, second)


def deduplicate_iupac(
    records: Iterable[DNARecord],
    *,
    config: IUPACDeduplicationConfig | None = None,
    verbose: bool = True,
) -> IUPACDeduplicationResult:
    """Greedily build stable complete-link groups of pairwise-compatible IUPAC strings."""

    resolved = IUPACDeduplicationConfig() if config is None else config
    if not isinstance(resolved, IUPACDeduplicationConfig):
        raise ConfigurationError(
            "config must be IUPACDeduplicationConfig or None.",
            code="INVALID_IUPAC_DEDUPLICATION_CONFIG",
        )
    materialized = materialize_limited(records, max_records=resolved.max_records)
    comparisons = ensure_pair_limit(len(materialized), resolved.max_pairwise_comparisons)
    for index, record in enumerate(materialized):
        if record.sequence.topology is not Topology.LINEAR:
            raise ConfigurationError(
                "IUPAC-aware deduplication currently requires linear sequences.",
                code="IUPAC_DEDUPLICATION_LINEAR_REQUIRED",
                context={"input_index": index},
            )
        if any(isinstance(part, Gap) for part in record.sequence.parts):
            raise UnsupportedGapOperationError(
                "IUPAC-aware deduplication does not bridge Gap objects.",
                code="IUPAC_DEDUPLICATION_GAP_UNSUPPORTED",
            )

    compatibility = [bytearray(len(materialized) - left - 1) for left in range(len(materialized))]
    identical_pairs = 0
    compatible_pairs = 0
    conflict_pairs = 0
    pair_relations: list[IUPACPairRelation] = []
    jobs = (
        (left, right, materialized[left].sequence.symbols, materialized[right].sequence.symbols)
        for left in range(len(materialized))
        for right in range(left + 1, len(materialized))
    )
    short = max((record.sequence.symbol_length for record in materialized), default=0) <= 4096
    for left, right, identical, is_compatible in ordered_map(
        _compatibility_job, jobs, workers=1 if short else None
    ):
        compatibility[left][right - left - 1] = is_compatible
        relation: Literal["identical", "compatible", "conflict"]
        if identical:
            identical_pairs += 1
            relation = "identical"
        elif is_compatible:
            compatible_pairs += 1
            relation = "compatible"
        else:
            conflict_pairs += 1
            relation = "conflict"
        pair_relations.append(
            IUPACPairRelation(
                left,
                right,
                materialized[left].id,
                materialized[right].id,
                relation,
            )
        )

    def pair_is_compatible(left: int, right: int) -> bool:
        if left > right:
            left, right = right, left
        return bool(compatibility[left][right - left - 1])

    groups: list[list[int]] = []
    for index in range(len(materialized)):
        for group in groups:
            if all(pair_is_compatible(index, member) for member in group):
                group.append(index)
                break
        else:
            groups.append([index])

    audit: list[IUPACDuplicateGroup] = []
    selected: list[DNARecord] = []
    for group_index, raw_members in enumerate(groups):
        members = tuple(raw_members)
        representative = _representative(members, materialized, resolved.representative_policy)
        selected.append(materialized[representative])
        symbols = {materialized[index].sequence.symbols for index in members}
        audit.append(
            IUPACDuplicateGroup(
                group_index,
                materialized[representative].id,
                tuple(materialized[index].id for index in members),
                (
                    "singleton"
                    if len(members) == 1
                    else "identical"
                    if len(symbols) == 1
                    else "compatible"
                ),
            )
        )
    return IUPACDeduplicationResult(
        DNASet(selected),
        tuple(audit),
        tuple(pair_relations),
        resolved.representative_policy,
        "position-compatible iff IUPAC base-set intersection is non-empty",
        "stable-greedy-complete-link",
        len(materialized),
        len(selected),
        comparisons,
        identical_pairs,
        compatible_pairs,
        conflict_pairs,
        resolved.max_records,
        resolved.max_pairwise_comparisons,
    )


def deduplicate_approximate(
    records: Iterable[DNARecord], *, config: ClusterConfig, verbose: bool = True
) -> ClusteringResult:
    """Return threshold-graph components and representatives as near-duplicate groups."""

    return cluster_sequences(records, config=config)


__all__ = ["deduplicate_approximate", "deduplicate_iupac"]

_register_public_verbose(globals())
