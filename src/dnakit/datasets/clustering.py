"""Bounded threshold and agglomerative clustering."""

from __future__ import annotations

import heapq
import math
from collections.abc import Iterable, Sequence

from dnakit.core import DNARecord, DNASet
from dnakit.exceptions import ConfigurationError

from ._advanced_shared import (
    UnionFind,
    ensure_pair_limit,
    materialize_limited,
    quality_score,
    similarity_matrix,
    validate_pairwise_sequences,
)
from .config import (
    AdvancedRepresentativePolicy,
    ClusterConfig,
    ClusterMethod,
    HierarchicalClusteringConfig,
)
from .results import (
    Cluster,
    ClusteringResult,
    HierarchicalClusteringResult,
    LinkageStep,
    RepresentativeSelectionResult,
)


def _choose_representative(
    members: tuple[int, ...],
    records: tuple[DNARecord, ...],
    *,
    policy: AdvancedRepresentativePolicy,
    similarities: tuple[tuple[float, ...], ...] | None,
) -> int:
    if policy == "first":
        return members[0]
    if policy == "shortest":
        return min(members, key=lambda index: (records[index].sequence.symbol_length, index))
    if policy == "longest":
        return max(members, key=lambda index: (records[index].sequence.symbol_length, -index))
    if policy == "best_quality":
        return max(members, key=lambda index: (quality_score(records[index]), -index))
    if similarities is None:
        raise ConfigurationError("Medoid selection requires a similarity matrix.")
    return min(
        members,
        key=lambda candidate: (
            math.fsum(1.0 - similarities[candidate][member] for member in members),
            candidate,
        ),
    )


def cluster_sequences(
    records: Iterable[DNARecord],
    *,
    config: ClusterConfig | None = None,
) -> ClusteringResult:
    """Cluster by connected components of an exhaustive similarity-threshold graph."""

    resolved = ClusterConfig() if config is None else config
    if not isinstance(resolved, ClusterConfig):
        raise ConfigurationError("config must be ClusterConfig.", code="INVALID_CLUSTER_CONFIG")
    materialized = materialize_limited(records, max_records=resolved.max_records)
    comparisons = ensure_pair_limit(len(materialized), resolved.max_pairwise_comparisons)
    validate_pairwise_sequences(materialized, operation="cluster_sequences")
    similarities = similarity_matrix(
        materialized,
        method=resolved.method,
        k=resolved.k,
        canonical=resolved.canonical,
        max_alignment_cells=resolved.max_alignment_cells,
    )
    components = UnionFind(len(materialized))
    for left in range(len(materialized)):
        for right in range(left + 1, len(materialized)):
            if similarities[left][right] >= resolved.threshold:
                components.union(left, right)
    groups = components.groups()
    labels = [0] * len(materialized)
    clusters: list[Cluster] = []
    representatives: list[DNARecord] = []
    for cluster_index, members in enumerate(groups):
        representative = _choose_representative(
            members,
            materialized,
            policy=resolved.representative_policy,
            similarities=similarities,
        )
        for member in members:
            labels[member] = cluster_index
        representatives.append(materialized[representative])
        clusters.append(
            Cluster(
                cluster_index,
                members,
                tuple(materialized[index].id for index in members),
                representative,
                materialized[representative].id,
            )
        )
    return ClusteringResult(
        tuple(clusters),
        tuple(labels),
        DNASet(representatives),
        resolved.method,
        resolved.threshold,
        resolved.k,
        resolved.canonical,
        resolved.representative_policy,
        "exhaustive-threshold-graph-connected-components",
        resolved.seed,
        False,
        comparisons,
        resolved.max_records,
        resolved.max_pairwise_comparisons,
        resolved.max_alignment_cells,
    )


def hierarchical_cluster(
    records: Iterable[DNARecord],
    *,
    config: HierarchicalClusteringConfig | None = None,
) -> HierarchicalClusteringResult:
    """Agglomerate the closest clusters with deterministic linkage tie-breaking."""

    resolved = HierarchicalClusteringConfig() if config is None else config
    if not isinstance(resolved, HierarchicalClusteringConfig):
        raise ConfigurationError("config must be HierarchicalClusteringConfig.")
    materialized = materialize_limited(records, max_records=resolved.max_records)
    comparisons = ensure_pair_limit(len(materialized), resolved.max_pairwise_comparisons)
    validate_pairwise_sequences(materialized, operation="hierarchical_cluster")
    similarities = similarity_matrix(
        materialized,
        method=resolved.method,
        k=resolved.k,
        canonical=resolved.canonical,
        max_alignment_cells=resolved.max_alignment_cells,
    )
    active = set(range(len(materialized)))
    sizes = {index: 1 for index in active}
    distances: dict[tuple[int, int], float] = {
        (left, right): 1.0 - similarities[left][right]
        for left in range(len(materialized))
        for right in range(left + 1, len(materialized))
    }
    candidates = [(distance, left, right) for (left, right), distance in distances.items()]
    heapq.heapify(candidates)
    next_node = len(materialized)
    linkage: list[LinkageStep] = []
    distance_updates = 0

    def key(left: int, right: int) -> tuple[int, int]:
        return (left, right) if left < right else (right, left)

    while len(active) > 1:
        while candidates:
            distance, left, right = heapq.heappop(candidates)
            if left in active and right in active and distances.get((left, right)) == distance:
                break
        else:
            raise AssertionError("Hierarchical clustering candidate heap was exhausted.")
        left_size, right_size = sizes[left], sizes[right]
        others = tuple(sorted(active - {left, right}))
        active.remove(left)
        active.remove(right)
        for other in others:
            left_distance = distances[key(left, other)]
            right_distance = distances[key(right, other)]
            if resolved.linkage == "single":
                new_distance = min(left_distance, right_distance)
            elif resolved.linkage == "complete":
                new_distance = max(left_distance, right_distance)
            else:
                new_distance = (left_size * left_distance + right_size * right_distance) / (
                    left_size + right_size
                )
            new_key = key(next_node, other)
            distances[new_key] = new_distance
            heapq.heappush(candidates, (new_distance, *new_key))
            distance_updates += 1
        for other in others:
            distances.pop(key(left, other))
            distances.pop(key(right, other))
        distances.pop(key(left, right))
        sizes[next_node] = left_size + right_size
        active.add(next_node)
        linkage.append(
            LinkageStep(
                len(linkage),
                left,
                right,
                next_node,
                distance,
                sizes[next_node],
            )
        )
        next_node += 1
    return HierarchicalClusteringResult(
        tuple(record.id for record in materialized),
        tuple(linkage),
        resolved.method,
        resolved.linkage,
        "1 - configured pair similarity",
        resolved.k,
        resolved.canonical,
        comparisons,
        distance_updates,
        resolved.max_records,
        resolved.max_pairwise_comparisons,
        resolved.max_alignment_cells,
    )


def select_representatives(
    records: Iterable[DNARecord],
    labels: Sequence[int],
    *,
    policy: AdvancedRepresentativePolicy = "first",
    medoid_method: ClusterMethod = "identity",
    k: int = 3,
    canonical: bool = False,
    max_records: int = 1_000,
    max_pairwise_comparisons: int = 500_000,
    max_alignment_cells: int = 5_000_000,
) -> RepresentativeSelectionResult:
    """Select one stable representative per user-supplied cluster label."""

    validation = ClusterConfig(
        method=medoid_method,
        k=k,
        canonical=canonical,
        representative_policy=policy,
        max_records=max_records,
        max_pairwise_comparisons=max_pairwise_comparisons,
        max_alignment_cells=max_alignment_cells,
    )
    materialized = materialize_limited(records, max_records=validation.max_records)
    label_tuple = tuple(labels)
    if len(label_tuple) != len(materialized) or any(
        isinstance(label, bool) or not isinstance(label, int) or label < 0 for label in label_tuple
    ):
        raise ConfigurationError("labels must be non-negative integers aligned with records.")
    grouped: dict[int, list[int]] = {}
    for index, label in enumerate(label_tuple):
        grouped.setdefault(label, []).append(index)
    comparisons = 0
    similarities: tuple[tuple[float, ...], ...] | None = None
    if policy == "medoid":
        comparisons = ensure_pair_limit(len(materialized), validation.max_pairwise_comparisons)
        validate_pairwise_sequences(materialized, operation="select_representatives")
        similarities = similarity_matrix(
            materialized,
            method=validation.method,
            k=validation.k,
            canonical=validation.canonical,
            max_alignment_cells=validation.max_alignment_cells,
        )
    chosen = tuple(
        _choose_representative(
            tuple(members),
            materialized,
            policy=policy,
            similarities=similarities,
        )
        for _, members in sorted(grouped.items())
    )
    return RepresentativeSelectionResult(
        DNASet(materialized[index] for index in chosen),
        tuple(materialized[index].id for index in chosen),
        policy,
        len(grouped),
        validation.method if policy == "medoid" else None,
        comparisons,
    )


__all__ = ["cluster_sequences", "hierarchical_cluster", "select_representatives"]
