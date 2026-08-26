"""Immutable, JSON-compatible summaries of dataset operations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal, cast

from dnakit.core import DNASet
from dnakit.core._json import FrozenDict, JSONValue, to_json_compatible
from dnakit.exceptions import ConfigurationError

DedupAction = Literal["deduplicated", "kept_all", "dropped"]
Orientation = Literal["forward", "reverse_complement"]


@dataclass(frozen=True, slots=True)
class DedupGroup:
    """One equivalence class and its chosen handling."""

    group_index: int
    representative_id: str
    member_ids: tuple[str, ...]
    orientations: tuple[Orientation, ...]
    conflict: bool = False
    conflict_values: tuple[JSONValue, ...] = ()
    missing_conflict_value_count: int = 0
    action: DedupAction = "deduplicated"
    rotation_offsets: tuple[int | None, ...] = ()
    rotation_offset_definition: str | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.group_index, bool)
            or not isinstance(self.group_index, int)
            or self.group_index < 0
        ):
            raise ConfigurationError("DedupGroup group_index must be non-negative.")
        if not self.representative_id or not self.member_ids:
            raise ConfigurationError("DedupGroup IDs must be non-empty.")
        if len(self.member_ids) != len(self.orientations):
            raise ConfigurationError("DedupGroup orientations must align with member_ids.")
        if self.representative_id not in self.member_ids:
            raise ConfigurationError("DedupGroup representative must be a member.")
        if self.action not in {"deduplicated", "kept_all", "dropped"}:
            raise ConfigurationError("DedupGroup action is invalid.")
        if self.rotation_offsets and len(self.rotation_offsets) != len(self.member_ids):
            raise ConfigurationError("DedupGroup rotation offsets must align with member_ids.")
        if any(
            value is not None
            and (isinstance(value, bool) or not isinstance(value, int) or value < 0)
            for value in self.rotation_offsets
        ):
            raise ConfigurationError("DedupGroup rotation offsets must be non-negative or None.")
        if (
            isinstance(self.missing_conflict_value_count, bool)
            or not isinstance(self.missing_conflict_value_count, int)
            or self.missing_conflict_value_count < 0
        ):
            raise ConfigurationError("Missing conflict value count must be non-negative.")


@dataclass(frozen=True, slots=True)
class DeduplicationResult:
    """Non-redundant records, equivalence groups, and conflict audit."""

    records: DNASet
    groups: tuple[DedupGroup, ...]
    equivalence: str
    representative_policy: str
    conflict_field: str | None
    conflict_policy: str
    input_count: int
    output_count: int
    duplicate_count: int
    conflicted_group_count: int
    removed_count: int
    merge_metadata: bool = False

    def __post_init__(self) -> None:
        if self.output_count != len(self.records):
            raise ConfigurationError("Deduplication output_count does not match records.")
        if self.input_count < 0 or self.output_count < 0 or self.output_count > self.input_count:
            raise ConfigurationError("Deduplication counts are inconsistent.")
        if self.removed_count != self.input_count - self.output_count:
            raise ConfigurationError("Deduplication removed_count is inconsistent.")
        if self.duplicate_count != self.input_count - len(self.groups):
            raise ConfigurationError("Deduplication duplicate_count is inconsistent.")
        if not isinstance(self.merge_metadata, bool):
            raise ConfigurationError("Deduplication merge_metadata must be a boolean.")

    def to_dict(self) -> dict[str, Any]:
        """Return an audit payload without duplicating full sequences."""

        payload = {
            "equivalence": self.equivalence,
            "representative_policy": self.representative_policy,
            "conflict_field": self.conflict_field,
            "conflict_policy": self.conflict_policy,
            "merge_metadata": self.merge_metadata,
            "input_count": self.input_count,
            "output_count": self.output_count,
            "duplicate_count": self.duplicate_count,
            "conflicted_group_count": self.conflicted_group_count,
            "removed_count": self.removed_count,
            "output_ids": self.records.ids,
            "groups": self.groups,
        }
        return cast(dict[str, Any], to_json_compatible(payload))


@dataclass(frozen=True, slots=True)
class SplitAssignment:
    """Assignment of one original record position to one named split."""

    input_index: int
    record_id: str
    split: str

    def __post_init__(self) -> None:
        if (
            isinstance(self.input_index, bool)
            or not isinstance(self.input_index, int)
            or self.input_index < 0
        ):
            raise ConfigurationError("SplitAssignment input_index must be non-negative.")
        if not self.record_id or not self.split:
            raise ConfigurationError("SplitAssignment names must be non-empty.")


@dataclass(frozen=True, slots=True)
class SplitSubset:
    """One named materialized dataset subset."""

    name: str
    records: DNASet

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigurationError("SplitSubset name must be non-empty.")


@dataclass(frozen=True, slots=True)
class SplitResult:
    """Materialized subsets plus deterministic assignment metadata."""

    subsets: tuple[SplitSubset, ...]
    assignments: tuple[SplitAssignment, ...]
    method: str
    ratios: FrozenDict
    counts: FrozenDict
    seed: int
    preserve_order: bool
    metadata_key: str | None = None
    component_count: int | None = None
    pairwise_comparison_count: int = 0
    similarity_method: str | None = None
    similarity_threshold: float | None = None
    shuffle: bool = True
    similarity_k: int | None = None
    max_pairwise_records: int | None = None
    assignment_strategy: str | None = None

    def __post_init__(self) -> None:
        names = tuple(subset.name for subset in self.subsets)
        if len(names) != len(set(names)) or set(names) != set(self.ratios):
            raise ConfigurationError("SplitResult subset names must match ratio names.")
        if set(self.counts) != set(names):
            raise ConfigurationError("SplitResult count names must match subset names.")
        if any(self.counts[name] != len(self.get(name)) for name in names):
            raise ConfigurationError("SplitResult counts do not match materialized subsets.")
        if sum(cast(int, self.counts[name]) for name in names) != len(self.assignments):
            raise ConfigurationError("SplitResult assignments do not match subset counts.")
        if not isinstance(self.shuffle, bool):
            raise ConfigurationError("SplitResult shuffle must be a boolean.")
        for field_name in ("similarity_k", "max_pairwise_records"):
            value = getattr(self, field_name)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int) or value <= 0
            ):
                raise ConfigurationError(f"SplitResult {field_name} must be positive or None.")
        if self.assignment_strategy is not None and (
            not isinstance(self.assignment_strategy, str) or not self.assignment_strategy.strip()
        ):
            raise ConfigurationError("SplitResult assignment_strategy must be non-empty or None.")

    def get(self, name: str) -> DNASet:
        """Return one subset by name."""

        for subset in self.subsets:
            if subset.name == name:
                return subset.records
        raise KeyError(name)

    def to_dict(self) -> dict[str, Any]:
        """Return an audit payload containing IDs rather than full sequences."""

        payload = {
            "method": self.method,
            "ratios": self.ratios,
            "counts": self.counts,
            "seed": self.seed,
            "shuffle": self.shuffle,
            "preserve_order": self.preserve_order,
            "metadata_key": self.metadata_key,
            "component_count": self.component_count,
            "pairwise_comparison_count": self.pairwise_comparison_count,
            "similarity_method": self.similarity_method,
            "similarity_threshold": self.similarity_threshold,
            "similarity_k": self.similarity_k,
            "max_pairwise_records": self.max_pairwise_records,
            "assignment_strategy": self.assignment_strategy,
            "subsets": {subset.name: subset.records.ids for subset in self.subsets},
            "assignments": self.assignments,
        }
        return cast(dict[str, Any], to_json_compatible(payload))


@dataclass(frozen=True, slots=True)
class IUPACDuplicateGroup:
    """Complete-link group whose member symbols are pairwise compatible."""

    group_index: int
    representative_id: str
    member_ids: tuple[str, ...]
    relation: Literal["singleton", "identical", "compatible"]


@dataclass(frozen=True, slots=True)
class IUPACPairRelation:
    """Auditable relation for one exhaustively compared input pair."""

    left_index: int
    right_index: int
    left_id: str
    right_id: str
    relation: Literal["identical", "compatible", "conflict"]


@dataclass(frozen=True, slots=True)
class IUPACDeduplicationResult:
    records: DNASet
    groups: tuple[IUPACDuplicateGroup, ...]
    pair_relations: tuple[IUPACPairRelation, ...]
    representative_policy: str
    compatibility_rule: str
    grouping_strategy: str
    input_count: int
    output_count: int
    pairwise_comparison_count: int
    identical_pair_count: int
    compatible_pair_count: int
    conflict_pair_count: int
    max_records: int
    max_pairwise_comparisons: int

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class Cluster:
    cluster_index: int
    member_indices: tuple[int, ...]
    member_ids: tuple[str, ...]
    representative_index: int
    representative_id: str


@dataclass(frozen=True, slots=True)
class ClusteringResult:
    clusters: tuple[Cluster, ...]
    labels: tuple[int, ...]
    representatives: DNASet
    method: str
    threshold: float
    k: int
    canonical: bool
    representative_policy: str
    clustering_strategy: str
    seed: int
    seed_used: bool
    pairwise_comparison_count: int
    max_records: int
    max_pairwise_comparisons: int
    max_alignment_cells: int

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class NeuralClusteringResult:
    """K-means groups obtained from sequence-level DNA model representations."""

    clusters: tuple[Cluster, ...]
    labels: tuple[int, ...]
    representatives: DNASet
    centers: tuple[tuple[float, ...], ...]
    model_name: str
    checkpoint_path: str | None
    pooling: str
    embedding_dimension: int
    clustering_dimension: int
    n_clusters: int
    inertia: float
    silhouette_score: float | None
    normalize: bool
    pca_components: int | None
    pca_explained_variance_ratio: tuple[float, ...]
    seed: int
    n_init: int
    max_iter: int
    iteration_count: int
    input_count: int
    clustering_strategy: str = "l2-normalize-optional-pca-kmeans++-lloyd"

    def __post_init__(self) -> None:
        if len(self.labels) != self.input_count:
            raise ConfigurationError("Neural clustering labels do not match input_count.")
        if len(self.clusters) != self.n_clusters or len(self.centers) != self.n_clusters:
            raise ConfigurationError("Neural clustering cluster counts are inconsistent.")
        if len(self.representatives) != self.n_clusters:
            raise ConfigurationError("Neural clustering representatives are inconsistent.")
        if any(len(center) != self.clustering_dimension for center in self.centers):
            raise ConfigurationError("Neural clustering center dimensions are inconsistent.")
        numeric = [self.inertia, *(value for center in self.centers for value in center)]
        if self.silhouette_score is not None:
            numeric.append(self.silhouette_score)
        numeric.extend(self.pca_explained_variance_ratio)
        if any(not math.isfinite(float(value)) for value in numeric):
            raise ConfigurationError("Neural clustering metrics must be finite.")

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible labels, centers, scores, and provenance."""

        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class LinkageStep:
    step_index: int
    left_node: int
    right_node: int
    new_node: int
    distance: float
    member_count: int


@dataclass(frozen=True, slots=True)
class HierarchicalClusteringResult:
    record_ids: tuple[str, ...]
    linkage: tuple[LinkageStep, ...]
    method: str
    linkage_method: str
    distance_definition: str
    k: int
    canonical: bool
    pairwise_comparison_count: int
    linkage_distance_update_count: int
    max_records: int
    max_pairwise_comparisons: int
    max_alignment_cells: int

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class RepresentativeSelectionResult:
    representatives: DNASet
    representative_ids: tuple[str, ...]
    policy: str
    cluster_count: int
    medoid_metric: str | None
    pairwise_comparison_count: int

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class TemporalSplitResult:
    subsets: tuple[SplitSubset, ...]
    assignments: tuple[SplitAssignment, ...]
    counts: FrozenDict
    ratios: FrozenDict
    metadata_key: str
    cutoffs: tuple[str, ...]
    strategy: str
    preserve_order: bool
    timezone_policy: str
    max_records: int

    def get(self, name: str) -> DNASet:
        for subset in self.subsets:
            if subset.name == name:
                return subset.records
        raise KeyError(name)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class JointSplitResult:
    subsets: tuple[SplitSubset, ...]
    assignments: tuple[SplitAssignment, ...]
    counts: FrozenDict
    target_counts: FrozenDict
    achieved_ratios: FrozenDict
    group_keys: tuple[str, ...]
    label_key: str | None
    similarity_threshold: float | None
    similarity_k: int
    similarity_canonical: bool
    seed: int
    strategy: str
    feasible: bool
    relaxed: bool
    relaxations: tuple[str, ...]
    max_ratio_deviation: float
    ratio_tolerance: float
    atomic_unit_count: int
    pairwise_comparison_count: int
    max_records: int
    max_pairwise_comparisons: int
    infeasible_policy: str

    def get(self, name: str) -> DNASet:
        for subset in self.subsets:
            if subset.name == name:
                return subset.records
        raise KeyError(name)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class LeakageEvent:
    left_split: str
    right_split: str
    left_id: str
    right_id: str
    left_index: int
    right_index: int
    similarity: float
    exact: bool


@dataclass(frozen=True, slots=True)
class LeakageReport:
    events: tuple[LeakageEvent, ...]
    method: str
    threshold: float
    strictness: str
    exact_event_count: int
    high_similarity_event_count: int
    cross_pair_count: int
    max_records: int
    max_cross_pairs: int
    max_events: int
    truncated: bool
    k: int
    canonical: bool
    max_alignment_cells: int

    @property
    def has_leakage(self) -> bool:
        return bool(self.events)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class SplitQualityResult:
    total_count: int
    counts: FrozenDict
    target_ratios: FrozenDict
    achieved_ratios: FrozenDict
    ratio_deviations: FrozenDict
    max_ratio_deviation: float
    label_key: str | None
    label_total_variation_by_split: FrozenDict
    group_keys: tuple[str, ...]
    group_leak_count: int
    leakage_event_count: int | None
    quality_score: float
    score_definition: str
    parameters: FrozenDict

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = [
    "Cluster",
    "ClusteringResult",
    "DedupAction",
    "DedupGroup",
    "DeduplicationResult",
    "HierarchicalClusteringResult",
    "IUPACDeduplicationResult",
    "IUPACDuplicateGroup",
    "IUPACPairRelation",
    "JointSplitResult",
    "LeakageEvent",
    "LeakageReport",
    "LinkageStep",
    "NeuralClusteringResult",
    "Orientation",
    "RepresentativeSelectionResult",
    "SplitAssignment",
    "SplitQualityResult",
    "SplitResult",
    "SplitSubset",
    "TemporalSplitResult",
]
