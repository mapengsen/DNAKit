"""Reproducible random, stable-hash, stratified, group, and similarity splits."""

from __future__ import annotations

import hashlib
import math
import random
from collections import deque
from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import cast

from dnakit.core import DNARecord, DNASet, Topology
from dnakit.core._json import FrozenDict
from dnakit.descriptors import kmer_statistics
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

from ._metadata import metadata_value_key
from .config import SplitConfig
from .results import SplitAssignment, SplitResult, SplitSubset


def _materialize(records: Iterable[DNARecord]) -> tuple[DNARecord, ...]:
    try:
        materialized = tuple(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be an iterable of DNARecord objects.",
            code="INVALID_DATASET_RECORDS",
        ) from exc
    for index, record in enumerate(materialized):
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "Every dataset item must be a DNARecord.",
                code="INVALID_DATASET_RECORD",
                context={"input_index": index, "type": type(record).__name__},
            )
    return materialized


def _materialize_limited(records: Iterable[DNARecord], *, limit: int) -> tuple[DNARecord, ...]:
    """Materialize no more than ``limit + 1`` items from a pairwise workload."""

    try:
        iterator = iter(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be an iterable of DNARecord objects.",
            code="INVALID_DATASET_RECORDS",
        ) from exc
    materialized: list[DNARecord] = []
    for index, record in enumerate(iterator):
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "Every dataset item must be a DNARecord.",
                code="INVALID_DATASET_RECORD",
                context={"input_index": index, "type": type(record).__name__},
            )
        materialized.append(record)
        if len(materialized) > limit:
            raise ConfigurationError(
                "Similarity split input exceeds max_pairwise_records.",
                code="SIMILARITY_SPLIT_SIZE_LIMIT",
                context={
                    "consumed_record_count": len(materialized),
                    "limit": limit,
                    "record_count_is_lower_bound": True,
                },
                hint="Pre-cluster records externally or explicitly raise the safety limit.",
            )
    return tuple(materialized)


def _allocate_counts(size: int, ratios: Mapping[str, float]) -> dict[str, int]:
    exact = {name: size * ratio for name, ratio in ratios.items()}
    allocated = {name: math.floor(value) for name, value in exact.items()}
    remainder = size - sum(allocated.values())
    order = sorted(ratios, key=lambda name: (-(exact[name] - allocated[name]), name))
    for name in order[:remainder]:
        allocated[name] += 1
    return allocated


def _assign_individuals(
    indices: Sequence[int],
    ratios: Mapping[str, float],
    *,
    generator: random.Random,
    shuffle: bool,
) -> dict[int, str]:
    ordered = list(indices)
    if shuffle:
        generator.shuffle(ordered)
    return _assign_ordered_individuals(ordered, ratios)


def _assign_ordered_individuals(
    ordered: Sequence[int], ratios: Mapping[str, float]
) -> dict[int, str]:
    counts = _allocate_counts(len(ordered), ratios)
    assigned: dict[int, str] = {}
    cursor = 0
    for name in ratios:
        end = cursor + counts[name]
        for index in ordered[cursor:end]:
            assigned[index] = name
        cursor = end
    return assigned


def _stable_hash_digest(seed: int, record_id: str) -> bytes:
    """Return the versioned, cross-process hash key for one record ID."""

    payload = (
        b"dnakit.split.hash.v1\0" + str(seed).encode("ascii") + b"\0" + record_id.encode("utf-8")
    )
    return hashlib.sha256(payload).digest()


def _assign_hashed_individuals(
    records: Sequence[DNARecord],
    ratios: Mapping[str, float],
    *,
    seed: int,
    shuffle: bool,
) -> dict[int, str]:
    """Assign records by stable ID hash, independently of input order."""

    seen: set[str] = set()
    duplicate_ids: set[str] = set()
    for record in records:
        if record.id in seen:
            duplicate_ids.add(record.id)
        seen.add(record.id)
    if duplicate_ids:
        raise ConfigurationError(
            "Hash split requires unique record IDs.",
            code="HASH_SPLIT_DUPLICATE_RECORD_ID",
            context={"duplicate_record_ids": tuple(sorted(duplicate_ids))},
            hint="Assign a stable unique ID to every record before using method='hash'.",
        )

    if shuffle:
        ordered = sorted(
            range(len(records)),
            key=lambda index: (
                _stable_hash_digest(seed, records[index].id),
                records[index].id.encode("utf-8"),
            ),
        )
    else:
        ordered = sorted(range(len(records)), key=lambda index: records[index].id.encode("utf-8"))
    return _assign_ordered_individuals(ordered, ratios)


def _assign_strata(
    strata: Sequence[tuple[int, ...]],
    ratios: Mapping[str, float],
    *,
    total_size: int,
    generator: random.Random,
    shuffle: bool,
) -> dict[int, str]:
    """Allocate strata against global quotas instead of rounding each stratum alone."""

    members = [list(stratum) for stratum in strata]
    if shuffle:
        generator.shuffle(members)
        for stratum in members:
            generator.shuffle(stratum)
    targets = _allocate_counts(total_size, ratios)
    remaining = dict(targets)
    split_order = {name: index for index, name in enumerate(ratios)}
    assigned_per_stratum = [{name: 0 for name in ratios} for _ in members]
    offsets = [0 for _ in members]
    active = deque(index for index, stratum in enumerate(members) if stratum)
    assignment: dict[int, str] = {}

    while active:
        stratum_index = active.popleft()
        stratum = members[stratum_index]
        local_counts = assigned_per_stratum[stratum_index]
        candidates = [name for name in ratios if remaining[name] > 0]
        if not candidates:
            raise AssertionError("Global split quotas were exhausted before assignment completed.")
        stratum_size = len(stratum)

        def priority(
            name: str,
            stratum_size: int = stratum_size,
            local_counts: dict[str, int] = local_counts,
        ) -> tuple[float, float, int]:
            local_deficit = stratum_size * ratios[name] - local_counts[name]
            target = targets[name]
            global_pressure = remaining[name] / target if target else 0.0
            return local_deficit, global_pressure, -split_order[name]

        chosen = max(candidates, key=priority)
        record_index = stratum[offsets[stratum_index]]
        assignment[record_index] = chosen
        local_counts[chosen] += 1
        remaining[chosen] -= 1
        offsets[stratum_index] += 1
        if offsets[stratum_index] < len(stratum):
            active.append(stratum_index)

    if any(remaining.values()):
        raise AssertionError("Global split quotas were not filled by stratified assignment.")
    return assignment


class _MissingMetadataGroup:
    """Identity-only key that cannot be represented by user JSON metadata."""

    __slots__ = ()


def _metadata_units(
    records: tuple[DNARecord, ...],
    key: str,
    *,
    missing_policy: str,
) -> tuple[tuple[int, ...], ...]:
    grouped: dict[Hashable, list[int]] = {}
    for index, record in enumerate(records):
        if key not in record.metadata:
            if missing_policy == "error":
                raise ConfigurationError(
                    "A record is missing the required split metadata.",
                    code="SPLIT_METADATA_MISSING",
                    context={"input_index": index, "record_id": record.id, "metadata_key": key},
                )
            group: Hashable = _MissingMetadataGroup()
        else:
            group = metadata_value_key(record.metadata[key])
        grouped.setdefault(group, []).append(index)
    return tuple(tuple(indices) for indices in grouped.values())


def _assign_units(
    units: Sequence[tuple[int, ...]],
    ratios: Mapping[str, float],
    *,
    total_size: int,
    generator: random.Random,
    shuffle: bool,
) -> dict[int, str]:
    ordered = list(units)
    if shuffle:
        generator.shuffle(ordered)
    ordered.sort(key=len, reverse=True)
    targets = {name: total_size * ratio for name, ratio in ratios.items()}
    counts = {name: 0 for name in ratios}
    assignment: dict[int, str] = {}
    for unit in ordered:
        unit_size = len(unit)

        def objective(candidate: str, unit_size: int = unit_size) -> tuple[float, int]:
            score = 0.0
            for name in ratios:
                projected = counts[name] + (unit_size if name == candidate else 0)
                denominator = max(targets[name], 1.0)
                score += ((projected - targets[name]) ** 2) / denominator
            return score, tuple(ratios).index(candidate)

        chosen = min(ratios, key=objective)
        counts[chosen] += unit_size
        assignment.update((index, chosen) for index in unit)
    return assignment


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _similarity_units(
    records: tuple[DNARecord, ...], config: SplitConfig
) -> tuple[tuple[tuple[int, ...], ...], int]:
    if len(records) > config.max_pairwise_records:
        raise ConfigurationError(
            "Similarity split input exceeds max_pairwise_records.",
            code="SIMILARITY_SPLIT_SIZE_LIMIT",
            context={"record_count": len(records), "limit": config.max_pairwise_records},
            hint="Pre-cluster records externally or explicitly raise the safety limit.",
        )
    fingerprints: list[frozenset[str]] = []
    for index, record in enumerate(records):
        if record.sequence.topology is Topology.CIRCULAR:
            raise ConfigurationError(
                "Basic similarity split does not model circular origin wrapping.",
                code="CIRCULAR_SIMILARITY_SPLIT_UNSUPPORTED",
                context={"input_index": index, "record_id": record.id},
            )
        if record.sequence.is_gapped and config.similarity_gap_policy == "error":
            raise UnsupportedGapOperationError(
                "Basic similarity split does not silently omit or cross sequence Gaps.",
                code="SIMILARITY_SPLIT_GAP_NOT_ALLOWED",
                context={"input_index": index, "record_id": record.id},
                hint="Choose similarity_gap_policy='split' to treat fragments independently.",
            )
        statistics = kmer_statistics(
            record.sequence,
            config.similarity_k,
            ambiguity_policy=config.similarity_ambiguity_policy,
            cross_gaps=False,
        )
        fingerprints.append(frozenset(statistics.counts))
    components = _UnionFind(len(records))
    comparisons = 0
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            comparisons += 1
            if _jaccard(fingerprints[left], fingerprints[right]) >= config.similarity_threshold:
                components.union(left, right)
    grouped: dict[int, list[int]] = {}
    for index in range(len(records)):
        grouped.setdefault(components.find(index), []).append(index)
    return tuple(tuple(indices) for indices in grouped.values()), comparisons


def split(records: Iterable[DNARecord], *, config: SplitConfig) -> SplitResult:
    """Split records reproducibly while keeping declared atomic groups intact."""

    if not isinstance(config, SplitConfig):
        raise ConfigurationError(
            "config must be a SplitConfig.",
            code="INVALID_SPLIT_CONFIG",
        )
    materialized = (
        _materialize_limited(records, limit=config.max_pairwise_records)
        if config.method == "similarity"
        else _materialize(records)
    )
    ratios = config.ratios
    generator = random.Random(config.seed)
    component_count: int | None = None
    comparisons = 0
    if config.method == "random":
        assignment = _assign_individuals(
            tuple(range(len(materialized))),
            ratios,
            generator=generator,
            shuffle=config.shuffle,
        )
    elif config.method == "hash":
        assignment = _assign_hashed_individuals(
            materialized,
            ratios,
            seed=config.seed,
            shuffle=config.shuffle,
        )
    elif config.method == "stratified":
        assert config.metadata_key is not None
        strata = _metadata_units(
            materialized,
            config.metadata_key,
            missing_policy=config.missing_metadata_policy,
        )
        assignment = _assign_strata(
            strata,
            ratios,
            total_size=len(materialized),
            generator=generator,
            shuffle=config.shuffle,
        )
    else:
        if config.method == "group":
            assert config.metadata_key is not None
            units = _metadata_units(
                materialized,
                config.metadata_key,
                missing_policy=config.missing_metadata_policy,
            )
        else:
            units, comparisons = _similarity_units(materialized, config)
            component_count = len(units)
        assignment = _assign_units(
            units,
            ratios,
            total_size=len(materialized),
            generator=generator,
            shuffle=config.shuffle,
        )

    indices_by_split: dict[str, list[int]] = {name: [] for name in ratios}
    for index in range(len(materialized)):
        indices_by_split[assignment[index]].append(index)
    if not config.preserve_order:
        rank = {index: rank for rank, index in enumerate(assignment)}
        for indices in indices_by_split.values():
            indices.sort(key=rank.__getitem__)
    subsets = tuple(
        SplitSubset(name, DNASet(materialized[index] for index in indices_by_split[name]))
        for name in ratios
    )
    assignments = tuple(
        SplitAssignment(index, record.id, assignment[index])
        for index, record in enumerate(materialized)
    )
    counts = FrozenDict({name: len(indices_by_split[name]) for name in ratios})
    return SplitResult(
        subsets=subsets,
        assignments=assignments,
        method=config.method,
        ratios=cast(FrozenDict, config.ratios),
        counts=counts,
        seed=config.seed,
        preserve_order=config.preserve_order,
        metadata_key=config.metadata_key,
        component_count=component_count,
        pairwise_comparison_count=comparisons,
        similarity_method="kmer_jaccard" if config.method == "similarity" else None,
        similarity_threshold=(
            float(config.similarity_threshold) if config.method == "similarity" else None
        ),
        shuffle=config.shuffle,
        similarity_k=config.similarity_k,
        max_pairwise_records=config.max_pairwise_records,
        assignment_strategy={
            "random": "global_largest_remainder_individual",
            "hash": "sha256_record_id_rank_v1_largest_remainder",
            "stratified": "global_quota_round_robin_stratified",
            "group": "atomic_unit_greedy_target_error",
            "similarity": "connected_component_atomic_unit_greedy_target_error",
        }[config.method],
    )


__all__ = ["split"]
