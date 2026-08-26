"""Validated configuration for deterministic dataset operations."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.representations.models import RepresentationConfig

DedupEquivalence: TypeAlias = Literal[
    "exact",
    "reverse_complement",
    "circular",
    "circular_reverse_complement",
]
RepresentativePolicy: TypeAlias = Literal["first", "last", "best_quality"]
ConflictPolicy: TypeAlias = Literal["error", "drop_group", "keep_representative", "keep_all"]
SplitMethod: TypeAlias = Literal["random", "hash", "stratified", "group", "similarity"]
MissingMetadataPolicy: TypeAlias = Literal["error", "separate"]
SimilarityAmbiguityPolicy: TypeAlias = Literal["error", "ignore"]
SimilarityGapPolicy: TypeAlias = Literal["error", "split"]
ClusterMethod: TypeAlias = Literal["identity", "edit", "kmer", "fingerprint"]
LinkageMethod: TypeAlias = Literal["single", "complete", "average"]
AdvancedRepresentativePolicy: TypeAlias = Literal[
    "first", "shortest", "longest", "best_quality", "medoid"
]
InfeasiblePolicy: TypeAlias = Literal["error", "relax"]


@dataclass(frozen=True, slots=True)
class DeduplicationConfig:
    """Control representative selection and duplicate-label conflicts."""

    representative_policy: RepresentativePolicy = "first"
    conflict_field: str | None = None
    conflict_policy: ConflictPolicy = "error"
    merge_metadata: bool = False

    def __post_init__(self) -> None:
        if self.representative_policy not in {"first", "last", "best_quality"}:
            raise ConfigurationError(
                "Unknown representative_policy.",
                code="INVALID_REPRESENTATIVE_POLICY",
                context={"representative_policy": self.representative_policy},
            )
        if self.conflict_policy not in {
            "error",
            "drop_group",
            "keep_representative",
            "keep_all",
        }:
            raise ConfigurationError(
                "Unknown conflict_policy.",
                code="INVALID_CONFLICT_POLICY",
                context={"conflict_policy": self.conflict_policy},
            )
        if self.conflict_field is not None and (
            not isinstance(self.conflict_field, str) or not self.conflict_field.strip()
        ):
            raise ConfigurationError(
                "conflict_field must be None or a non-empty metadata key.",
                code="INVALID_CONFLICT_FIELD",
            )
        if not isinstance(self.merge_metadata, bool):
            raise ConfigurationError(
                "merge_metadata must be a boolean.",
                code="INVALID_DEDUPLICATION_CONFIG",
            )


def _default_ratios() -> dict[str, float]:
    return {"train": 0.8, "valid": 0.1, "test": 0.1}


@dataclass(frozen=True, slots=True)
class SplitConfig:
    """Control random, stable-hash, stratified, group-aware, or k-mer similarity splitting."""

    method: SplitMethod = "random"
    ratios: Mapping[str, float] = field(default_factory=_default_ratios)
    seed: int = 0
    shuffle: bool = True
    preserve_order: bool = True
    metadata_key: str | None = None
    missing_metadata_policy: MissingMetadataPolicy = "error"
    similarity_k: int = 3
    similarity_threshold: float = 0.8
    similarity_ambiguity_policy: SimilarityAmbiguityPolicy = "error"
    similarity_gap_policy: SimilarityGapPolicy = "error"
    max_pairwise_records: int = 5_000

    def __post_init__(self) -> None:
        if self.method not in {"random", "hash", "stratified", "group", "similarity"}:
            raise ConfigurationError(
                "Unknown dataset split method.",
                code="INVALID_SPLIT_METHOD",
                context={"method": self.method},
            )
        if not isinstance(self.ratios, Mapping):
            raise ConfigurationError(
                "ratios must be a mapping of split names to fractions.",
                code="INVALID_SPLIT_RATIOS",
            )
        ratio_items = tuple(self.ratios.items())
        if len(ratio_items) < 2:
            raise ConfigurationError(
                "At least two split ratios are required.",
                code="INVALID_SPLIT_RATIOS",
            )
        normalized_ratios: dict[str, float] = {}
        for name, value in ratio_items:
            if not isinstance(name, str) or not name.strip():
                raise ConfigurationError(
                    "Split names must be non-empty strings.",
                    code="INVALID_SPLIT_RATIOS",
                )
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ConfigurationError(
                    "Every split ratio must be a finite positive number.",
                    code="INVALID_SPLIT_RATIOS",
                    context={"split": name, "ratio": value},
                )
            normalized_ratios[name] = float(value)
        ratio_sum = sum(normalized_ratios.values())
        if not math.isclose(ratio_sum, 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ConfigurationError(
                "Split ratios must sum to 1.0.",
                code="INVALID_SPLIT_RATIOS",
                context={"sum": ratio_sum},
            )
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise ConfigurationError(
                "seed must be an integer.",
                code="INVALID_SPLIT_SEED",
            )
        for field_name in ("shuffle", "preserve_order"):
            if not isinstance(getattr(self, field_name), bool):
                raise ConfigurationError(
                    f"{field_name} must be a boolean.",
                    code="INVALID_SPLIT_CONFIG",
                )
        if self.method in {"stratified", "group"} and (
            not isinstance(self.metadata_key, str) or not self.metadata_key.strip()
        ):
            raise ConfigurationError(
                f"metadata_key is required for method={self.method!r}.",
                code="SPLIT_METADATA_KEY_REQUIRED",
            )
        if self.metadata_key is not None and (
            not isinstance(self.metadata_key, str) or not self.metadata_key.strip()
        ):
            raise ConfigurationError(
                "metadata_key must be None or a non-empty string.",
                code="INVALID_SPLIT_METADATA_KEY",
            )
        if self.missing_metadata_policy not in {"error", "separate"}:
            raise ConfigurationError(
                "missing_metadata_policy must be 'error' or 'separate'.",
                code="INVALID_MISSING_METADATA_POLICY",
            )
        if (
            isinstance(self.similarity_k, bool)
            or not isinstance(self.similarity_k, int)
            or self.similarity_k <= 0
        ):
            raise ConfigurationError(
                "similarity_k must be a positive integer.",
                code="INVALID_SIMILARITY_SPLIT_CONFIG",
            )
        if (
            isinstance(self.similarity_threshold, bool)
            or not isinstance(self.similarity_threshold, (int, float))
            or not math.isfinite(self.similarity_threshold)
            or not 0 <= self.similarity_threshold <= 1
        ):
            raise ConfigurationError(
                "similarity_threshold must be finite and within [0, 1].",
                code="INVALID_SIMILARITY_SPLIT_CONFIG",
            )
        if self.similarity_ambiguity_policy not in {"error", "ignore"}:
            raise ConfigurationError(
                "similarity_ambiguity_policy must be 'error' or 'ignore'.",
                code="INVALID_SIMILARITY_SPLIT_CONFIG",
            )
        if self.similarity_gap_policy not in {"error", "split"}:
            raise ConfigurationError(
                "similarity_gap_policy must be 'error' or 'split'.",
                code="INVALID_SIMILARITY_SPLIT_CONFIG",
            )
        if (
            isinstance(self.max_pairwise_records, bool)
            or not isinstance(self.max_pairwise_records, int)
            or self.max_pairwise_records <= 0
        ):
            raise ConfigurationError(
                "max_pairwise_records must be a positive integer.",
                code="INVALID_SIMILARITY_SPLIT_CONFIG",
            )
        object.__setattr__(self, "ratios", FrozenDict(normalized_ratios))


def _validate_probability(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not 0.0 <= value <= 1.0
    ):
        raise ConfigurationError(
            f"{name} must be finite and within [0, 1].",
            code="INVALID_ADVANCED_DATASET_CONFIG",
        )
    return float(value)


def _validate_bounded_positive_int(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0 or value > maximum:
        raise ConfigurationError(
            f"{name} must be an integer in [1, {maximum}].",
            code="INVALID_ADVANCED_DATASET_LIMIT",
        )
    return value


@dataclass(frozen=True, slots=True)
class IUPACDeduplicationConfig:
    """Bound pairwise work for all-pairs IUPAC compatibility grouping."""

    representative_policy: AdvancedRepresentativePolicy = "first"
    max_records: int = 2_000
    max_pairwise_comparisons: int = 1_000_000

    def __post_init__(self) -> None:
        if self.representative_policy not in {
            "first",
            "shortest",
            "longest",
            "best_quality",
            "medoid",
        }:
            raise ConfigurationError("Unknown advanced representative policy.")
        if self.representative_policy == "medoid":
            raise ConfigurationError(
                "IUPAC compatibility groups do not define a numeric medoid distance.",
                code="IUPAC_MEDOID_UNDEFINED",
            )
        _validate_bounded_positive_int(self.max_records, "max_records", 10_000)
        _validate_bounded_positive_int(
            self.max_pairwise_comparisons,
            "max_pairwise_comparisons",
            50_000_000,
        )


@dataclass(frozen=True, slots=True)
class ClusterConfig:
    """Configuration for exhaustive bounded threshold-graph clustering."""

    method: ClusterMethod = "identity"
    threshold: float = 0.9
    k: int = 3
    canonical: bool = False
    representative_policy: AdvancedRepresentativePolicy = "first"
    seed: int = 0
    max_records: int = 1_000
    max_pairwise_comparisons: int = 500_000
    max_alignment_cells: int = 5_000_000

    def __post_init__(self) -> None:
        if self.method not in {"identity", "edit", "kmer", "fingerprint"}:
            raise ConfigurationError("Unknown clustering method.", code="INVALID_CLUSTER_METHOD")
        object.__setattr__(self, "threshold", _validate_probability(self.threshold, "threshold"))
        _validate_bounded_positive_int(self.k, "k", 12)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        if self.representative_policy not in {
            "first",
            "shortest",
            "longest",
            "best_quality",
            "medoid",
        }:
            raise ConfigurationError("Unknown advanced representative policy.")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise ConfigurationError("seed must be an integer.")
        _validate_bounded_positive_int(self.max_records, "max_records", 5_000)
        _validate_bounded_positive_int(
            self.max_pairwise_comparisons,
            "max_pairwise_comparisons",
            10_000_000,
        )
        _validate_bounded_positive_int(
            self.max_alignment_cells,
            "max_alignment_cells",
            100_000_000,
        )


@dataclass(frozen=True, slots=True)
class HierarchicalClusteringConfig:
    """Configuration for exact agglomerative clustering over a bounded matrix."""

    method: ClusterMethod = "identity"
    linkage: LinkageMethod = "average"
    k: int = 3
    canonical: bool = False
    max_records: int = 500
    max_pairwise_comparisons: int = 124_750
    max_alignment_cells: int = 5_000_000

    def __post_init__(self) -> None:
        if self.method not in {"identity", "edit", "kmer", "fingerprint"}:
            raise ConfigurationError("Unknown clustering method.", code="INVALID_CLUSTER_METHOD")
        if self.linkage not in {"single", "complete", "average"}:
            raise ConfigurationError("Unknown linkage method.", code="INVALID_LINKAGE_METHOD")
        _validate_bounded_positive_int(self.k, "k", 12)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        _validate_bounded_positive_int(self.max_records, "max_records", 1_000)
        _validate_bounded_positive_int(
            self.max_pairwise_comparisons,
            "max_pairwise_comparisons",
            499_500,
        )
        _validate_bounded_positive_int(
            self.max_alignment_cells,
            "max_alignment_cells",
            100_000_000,
        )


@dataclass(frozen=True, slots=True)
class NeuralClusteringConfig:
    """Configure DNA model representations followed by seeded k-means."""

    representation: RepresentationConfig = field(default_factory=RepresentationConfig)
    n_clusters: int = 2
    n_init: int = 10
    max_iter: int = 300
    tolerance: float = 1e-4
    normalize: bool = True
    pca_components: int | None = None
    silhouette_sample_size: int = 5_000
    seed: int = 0
    max_records: int = 10_000

    def __post_init__(self) -> None:
        if not isinstance(self.representation, RepresentationConfig):
            raise ConfigurationError(
                "representation must be RepresentationConfig.",
                code="INVALID_NEURAL_CLUSTER_CONFIG",
            )
        _validate_bounded_positive_int(self.n_clusters, "n_clusters", 100_000)
        _validate_bounded_positive_int(self.n_init, "n_init", 1_000)
        _validate_bounded_positive_int(self.max_iter, "max_iter", 100_000)
        if (
            isinstance(self.tolerance, bool)
            or not isinstance(self.tolerance, (int, float))
            or not math.isfinite(float(self.tolerance))
            or not 0 < float(self.tolerance) <= 1.0
        ):
            raise ConfigurationError(
                "tolerance must be finite and in (0, 1].",
                code="INVALID_NEURAL_CLUSTER_CONFIG",
            )
        object.__setattr__(self, "tolerance", float(self.tolerance))
        if not isinstance(self.normalize, bool):
            raise ConfigurationError(
                "normalize must be boolean.",
                code="INVALID_NEURAL_CLUSTER_CONFIG",
            )
        if self.pca_components is not None:
            _validate_bounded_positive_int(self.pca_components, "pca_components", 100_000)
        _validate_bounded_positive_int(
            self.silhouette_sample_size,
            "silhouette_sample_size",
            1_000_000,
        )
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise ConfigurationError(
                "seed must be an integer.",
                code="INVALID_NEURAL_CLUSTER_CONFIG",
            )
        _validate_bounded_positive_int(self.max_records, "max_records", 1_000_000)


@dataclass(frozen=True, slots=True)
class TemporalSplitConfig:
    """Chronological split configuration using ISO-8601 metadata strings."""

    metadata_key: str = "date"
    ratios: Mapping[str, float] = field(default_factory=lambda: {"train": 0.8, "test": 0.2})
    cutoffs: tuple[str, ...] | None = None
    preserve_order: bool = True
    max_records: int = 1_000_000

    def __post_init__(self) -> None:
        if not isinstance(self.metadata_key, str) or not self.metadata_key.strip():
            raise ConfigurationError("metadata_key must be a non-empty string.")
        validated = SplitConfig(ratios=self.ratios, preserve_order=self.preserve_order)
        object.__setattr__(self, "ratios", validated.ratios)
        _validate_bounded_positive_int(self.max_records, "max_records", 10_000_000)
        if self.cutoffs is not None:
            cutoffs = tuple(self.cutoffs)
            if len(cutoffs) != len(self.ratios) - 1 or any(
                not isinstance(value, str) or not value.strip() for value in cutoffs
            ):
                raise ConfigurationError(
                    "cutoffs must contain one non-empty ISO timestamp per split boundary.",
                    code="INVALID_TEMPORAL_CUTOFFS",
                )
            object.__setattr__(self, "cutoffs", cutoffs)


@dataclass(frozen=True, slots=True)
class JointSplitConfig:
    """Greedy multi-constraint split with explicit infeasibility handling."""

    ratios: Mapping[str, float] = field(default_factory=_default_ratios)
    group_keys: tuple[str, ...] = ()
    label_key: str | None = None
    similarity_threshold: float | None = None
    similarity_k: int = 3
    similarity_canonical: bool = True
    seed: int = 0
    ratio_tolerance: float = 0.05
    infeasible_policy: InfeasiblePolicy = "error"
    max_records: int = 2_000
    max_pairwise_comparisons: int = 1_000_000

    def __post_init__(self) -> None:
        validated = SplitConfig(ratios=self.ratios, seed=self.seed)
        object.__setattr__(self, "ratios", validated.ratios)
        keys = tuple(self.group_keys)
        if any(not isinstance(key, str) or not key.strip() for key in keys):
            raise ConfigurationError("group_keys must contain non-empty strings.")
        if len(keys) != len(set(keys)):
            raise ConfigurationError("group_keys must be unique.")
        object.__setattr__(self, "group_keys", keys)
        if self.label_key is not None and (
            not isinstance(self.label_key, str) or not self.label_key.strip()
        ):
            raise ConfigurationError("label_key must be non-empty or None.")
        if self.similarity_threshold is not None:
            object.__setattr__(
                self,
                "similarity_threshold",
                _validate_probability(self.similarity_threshold, "similarity_threshold"),
            )
        if not keys and self.similarity_threshold is None:
            raise ConfigurationError(
                "At least one group key or similarity constraint is required.",
                code="JOINT_SPLIT_CONSTRAINT_REQUIRED",
            )
        _validate_bounded_positive_int(self.similarity_k, "similarity_k", 12)
        if not isinstance(self.similarity_canonical, bool):
            raise ConfigurationError("similarity_canonical must be boolean.")
        object.__setattr__(
            self, "ratio_tolerance", _validate_probability(self.ratio_tolerance, "ratio_tolerance")
        )
        if self.infeasible_policy not in {"error", "relax"}:
            raise ConfigurationError("infeasible_policy must be 'error' or 'relax'.")
        _validate_bounded_positive_int(self.max_records, "max_records", 5_000)
        _validate_bounded_positive_int(
            self.max_pairwise_comparisons,
            "max_pairwise_comparisons",
            10_000_000,
        )


@dataclass(frozen=True, slots=True)
class LeakageConfig:
    """Cross-split exhaustive comparison controls."""

    method: ClusterMethod = "identity"
    threshold: float = 0.9
    k: int = 3
    canonical: bool = False
    max_records: int = 5_000
    max_cross_pairs: int = 1_000_000
    max_events: int = 100_000
    max_alignment_cells: int = 5_000_000

    def __post_init__(self) -> None:
        if self.method not in {"identity", "edit", "kmer", "fingerprint"}:
            raise ConfigurationError("Unknown leakage method.")
        object.__setattr__(self, "threshold", _validate_probability(self.threshold, "threshold"))
        _validate_bounded_positive_int(self.k, "k", 12)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        _validate_bounded_positive_int(self.max_records, "max_records", 20_000)
        _validate_bounded_positive_int(self.max_cross_pairs, "max_cross_pairs", 50_000_000)
        _validate_bounded_positive_int(self.max_events, "max_events", 5_000_000)
        _validate_bounded_positive_int(
            self.max_alignment_cells,
            "max_alignment_cells",
            100_000_000,
        )


__all__ = [
    "AdvancedRepresentativePolicy",
    "ClusterConfig",
    "ClusterMethod",
    "ConflictPolicy",
    "DedupEquivalence",
    "DeduplicationConfig",
    "HierarchicalClusteringConfig",
    "IUPACDeduplicationConfig",
    "InfeasiblePolicy",
    "JointSplitConfig",
    "LeakageConfig",
    "LinkageMethod",
    "MissingMetadataPolicy",
    "NeuralClusteringConfig",
    "RepresentativePolicy",
    "SimilarityAmbiguityPolicy",
    "SimilarityGapPolicy",
    "SplitConfig",
    "SplitMethod",
    "TemporalSplitConfig",
]
