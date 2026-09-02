"""Validated configuration for transparent DNA evaluation methods."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, TypeAlias

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.representations import RepresentationConfig

GapDenominatorPolicy: TypeAlias = Literal["exclude", "include_known", "error"]
PairSimilarityMethod: TypeAlias = Literal["exact", "identity", "edit", "kmer"]
DiversityCalculation: TypeAlias = Literal["similarity", "levenshtein"]
NoveltyCalculation: TypeAlias = Literal["similarity", "levenshtein"]
UniquenessEquivalence: TypeAlias = Literal[
    "exact",
    "reverse_complement",
    "circular",
    "circular_reverse_complement",
    "iupac",
    "approximate",
]
ScoreDirection: TypeAlias = Literal["higher_is_better", "lower_is_better"]
MissingScorePolicy: TypeAlias = Literal["error", "omit", "zero"]
KmerAmbiguityPolicy: TypeAlias = Literal["error", "ignore"]


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _probability(value: object, name: str) -> float:
    resolved = _finite(value, name)
    if not 0.0 <= resolved <= 1.0:
        raise ConfigurationError(f"{name} must be within [0, 1].")
    return resolved


def _positive_int(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 < value <= maximum:
        raise ConfigurationError(f"{name} must be an integer within [1, {maximum}].")
    return value


def _non_negative_int(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer within [0, {maximum}].")
    return value


def _freeze_weights(
    values: Mapping[str, float],
    *,
    name: str,
    allowed: frozenset[str] | None = None,
) -> FrozenDict:
    if not isinstance(values, Mapping) or not values:
        raise ConfigurationError(f"{name} must be a non-empty mapping.")
    resolved: dict[str, float] = {}
    for key, value in values.items():
        if not isinstance(key, str) or not key.strip():
            raise ConfigurationError(f"{name} keys must be non-empty strings.")
        if allowed is not None and key not in allowed:
            raise ConfigurationError(f"{name} contains an unknown component: {key!r}.")
        weight = _finite(value, f"{name}[{key!r}]")
        if weight < 0:
            raise ConfigurationError(f"{name} values must be non-negative.")
        resolved[key] = weight
    if math.fsum(resolved.values()) <= 0:
        raise ConfigurationError(f"{name} must contain at least one positive weight.")
    return FrozenDict(resolved)


@dataclass(frozen=True, slots=True)
class EvaluationLimits:
    """Shared caps for materialization and exhaustive calculations."""

    max_records: int = 1_000
    max_total_symbols: int = 10_000_000
    max_pairwise_comparisons: int = 500_000
    max_alignment_cells: int = 5_000_000

    def __post_init__(self) -> None:
        _positive_int(self.max_records, "max_records", 1_000_000)
        _positive_int(self.max_total_symbols, "max_total_symbols", 1_000_000_000)
        _positive_int(self.max_pairwise_comparisons, "max_pairwise_comparisons", 50_000_000)
        _positive_int(self.max_alignment_cells, "max_alignment_cells", 100_000_000)


def _ambiguity_weights() -> dict[str, float]:
    return {
        "R": 0.5,
        "Y": 0.5,
        "S": 0.5,
        "W": 0.5,
        "K": 0.5,
        "M": 0.5,
        "B": 0.75,
        "D": 0.75,
        "H": 0.75,
        "V": 0.75,
        "N": 1.0,
    }


@dataclass(frozen=True, slots=True)
class AmbiguityEvaluationConfig:
    max_fraction: float = 0.05
    symbol_weights: Mapping[str, float] = field(default_factory=_ambiguity_weights)
    gap_denominator_policy: GapDenominatorPolicy = "exclude"
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_fraction", _probability(self.max_fraction, "max_fraction"))
        if self.gap_denominator_policy not in {"exclude", "include_known", "error"}:
            raise ConfigurationError("Unknown gap_denominator_policy.")
        allowed = frozenset("RYSWKMBDHVN")
        weights = _freeze_weights(self.symbol_weights, name="symbol_weights", allowed=allowed)
        if set(weights) != allowed:
            raise ConfigurationError("symbol_weights must define every ambiguous IUPAC symbol.")
        if any(
            isinstance(value, (int, float)) and not isinstance(value, bool) and value > 1.0
            for value in weights.values()
        ):
            raise ConfigurationError("Ambiguity symbol weights must be within [0, 1].")
        object.__setattr__(self, "symbol_weights", weights)
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


def _complexity_weights() -> dict[str, float]:
    return {
        "entropy": 1.0,
        "linguistic": 1.0,
        "repeat_cleanliness": 1.0,
        "homopolymer_cleanliness": 1.0,
    }


@dataclass(frozen=True, slots=True)
class ComplexityEvaluationConfig:
    max_word_size: int = 6
    max_repeat_unit: int = 20
    min_repeat_count: int = 2
    acceptable_homopolymer_length: int = 6
    weights: Mapping[str, float] = field(default_factory=_complexity_weights)
    max_observations_per_sequence: int = 10_000_000
    max_comparisons_per_sequence: int = 5_000_000
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        _positive_int(self.max_word_size, "max_word_size", 16)
        _positive_int(self.max_repeat_unit, "max_repeat_unit", 100)
        if self.min_repeat_count < 2:
            raise ConfigurationError("min_repeat_count must be at least 2.")
        _positive_int(self.min_repeat_count, "min_repeat_count", 1_000_000)
        _positive_int(
            self.acceptable_homopolymer_length, "acceptable_homopolymer_length", 1_000_000
        )
        _positive_int(
            self.max_observations_per_sequence, "max_observations_per_sequence", 100_000_000
        )
        _positive_int(
            self.max_comparisons_per_sequence, "max_comparisons_per_sequence", 100_000_000
        )
        object.__setattr__(
            self,
            "weights",
            _freeze_weights(
                self.weights,
                name="weights",
                allowed=frozenset(_complexity_weights()),
            ),
        )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


def _quality_weights() -> dict[str, float]:
    return {
        "validity": 2.0,
        "ambiguity": 1.0,
        "complexity": 1.0,
        "length": 1.0,
        "completeness": 1.0,
    }


@dataclass(frozen=True, slots=True)
class QualityEvaluationConfig:
    min_length: int = 1
    max_length: int | None = None
    ambiguity: AmbiguityEvaluationConfig = field(default_factory=AmbiguityEvaluationConfig)
    complexity: ComplexityEvaluationConfig = field(default_factory=ComplexityEvaluationConfig)
    weights: Mapping[str, float] = field(default_factory=_quality_weights)
    warning_score: float = 0.75
    pass_score: float = 0.9

    def __post_init__(self) -> None:
        _non_negative_int(self.min_length, "min_length", 1_000_000_000)
        if self.max_length is not None:
            _positive_int(self.max_length, "max_length", 1_000_000_000)
            if self.max_length < self.min_length:
                raise ConfigurationError("max_length cannot be smaller than min_length.")
        if not isinstance(self.ambiguity, AmbiguityEvaluationConfig):
            raise ConfigurationError("ambiguity must be AmbiguityEvaluationConfig.")
        if not isinstance(self.complexity, ComplexityEvaluationConfig):
            raise ConfigurationError("complexity must be ComplexityEvaluationConfig.")
        if self.ambiguity.limits != self.complexity.limits:
            raise ConfigurationError(
                "Quality ambiguity and complexity configs must use identical limits.",
                code="QUALITY_LIMIT_MISMATCH",
            )
        object.__setattr__(
            self,
            "weights",
            _freeze_weights(self.weights, name="weights", allowed=frozenset(_quality_weights())),
        )
        warning = _probability(self.warning_score, "warning_score")
        passing = _probability(self.pass_score, "pass_score")
        if warning > passing:
            raise ConfigurationError("warning_score cannot exceed pass_score.")
        object.__setattr__(self, "warning_score", warning)
        object.__setattr__(self, "pass_score", passing)


@dataclass(frozen=True, slots=True)
class UniquenessEvaluationConfig:
    equivalence: UniquenessEquivalence = "exact"
    approximate_method: PairSimilarityMethod = "kmer"
    threshold: float = 0.9
    k: int = 5
    canonical: bool = True
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        if self.equivalence not in {
            "exact",
            "reverse_complement",
            "circular",
            "circular_reverse_complement",
            "iupac",
            "approximate",
        }:
            raise ConfigurationError("Unknown uniqueness equivalence.")
        if self.approximate_method not in {"exact", "identity", "edit", "kmer"}:
            raise ConfigurationError("Unknown approximate_method.")
        object.__setattr__(self, "threshold", _probability(self.threshold, "threshold"))
        _positive_int(self.k, "k", 32)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


@dataclass(frozen=True, slots=True)
class DiversityEvaluationConfig:
    method: PairSimilarityMethod = "kmer"
    k: int = 5
    canonical: bool = True
    cluster_threshold: float = 0.9
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)
    calculation: DiversityCalculation = "similarity"
    show_progress: bool = False

    def __post_init__(self) -> None:
        if self.method not in {"exact", "identity", "edit", "kmer"}:
            raise ConfigurationError("Unknown diversity similarity method.")
        _positive_int(self.k, "k", 32)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        object.__setattr__(
            self, "cluster_threshold", _probability(self.cluster_threshold, "cluster_threshold")
        )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")
        if self.calculation not in {"similarity", "levenshtein"}:
            raise ConfigurationError("Unknown diversity calculation.")
        if not isinstance(self.show_progress, bool):
            raise ConfigurationError("show_progress must be boolean.")


@dataclass(frozen=True, slots=True)
class ReferenceSearchConfig:
    method: PairSimilarityMethod = "kmer"
    k: int = 7
    canonical: bool = True
    top_k: int = 1
    min_similarity: float = 0.0
    min_coverage: float = 0.0
    copy_threshold: float = 0.9
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)
    novelty_calculation: NoveltyCalculation = "similarity"
    show_progress: bool = False

    def __post_init__(self) -> None:
        if self.method not in {"exact", "identity", "edit", "kmer"}:
            raise ConfigurationError("Unknown reference similarity method.")
        _positive_int(self.k, "k", 32)
        _positive_int(self.top_k, "top_k", 100_000)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        object.__setattr__(
            self, "min_similarity", _probability(self.min_similarity, "min_similarity")
        )
        object.__setattr__(self, "min_coverage", _probability(self.min_coverage, "min_coverage"))
        object.__setattr__(
            self, "copy_threshold", _probability(self.copy_threshold, "copy_threshold")
        )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")
        if self.novelty_calculation not in {"similarity", "levenshtein"}:
            raise ConfigurationError("Unknown novelty calculation.")
        if not isinstance(self.show_progress, bool):
            raise ConfigurationError("show_progress must be boolean.")


@dataclass(frozen=True, slots=True)
class DistributionEvaluationConfig:
    features: tuple[str, ...] = ("length", "gc", "kmer", "motif", "repeat")
    k: int = 3
    canonical: bool = True
    motifs: tuple[str, ...] = ("CG",)
    max_kmer_observations: int = 10_000_000
    max_repeat_comparisons_per_sequence: int = 1_000_000
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        allowed = {"length", "gc", "kmer", "motif", "repeat"}
        features = tuple(self.features)
        if not features or len(set(features)) != len(features) or set(features) - allowed:
            raise ConfigurationError("features must be unique supported distribution features.")
        object.__setattr__(self, "features", features)
        _positive_int(self.k, "k", 12)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("canonical must be boolean.")
        motifs = tuple(self.motifs)
        if len(set(motifs)) != len(motifs) or any(
            not isinstance(motif, str) or not motif or set(motif) - set("ACGT") for motif in motifs
        ):
            raise ConfigurationError("motifs must contain unique, non-empty A/C/G/T strings.")
        if "motif" in features and not motifs:
            raise ConfigurationError("At least one motif is required for motif distribution.")
        object.__setattr__(self, "motifs", motifs)
        _positive_int(self.max_kmer_observations, "max_kmer_observations", 100_000_000)
        _positive_int(
            self.max_repeat_comparisons_per_sequence,
            "max_repeat_comparisons_per_sequence",
            100_000_000,
        )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


@dataclass(frozen=True, slots=True)
class FrechetDistanceConfig:
    """Configure a Fréchet distance over DNA foundation-model representations."""

    representation: RepresentationConfig = field(default_factory=RepresentationConfig)
    normalize: bool = True
    max_cross_gram_elements: int = 1_000_000
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        if not isinstance(self.representation, RepresentationConfig):
            raise ConfigurationError(
                "representation must be RepresentationConfig.",
                code="INVALID_FRECHET_DISTANCE_CONFIG",
            )
        if not isinstance(self.normalize, bool):
            raise ConfigurationError(
                "normalize must be boolean.",
                code="INVALID_FRECHET_DISTANCE_CONFIG",
            )
        _positive_int(
            self.max_cross_gram_elements,
            "max_cross_gram_elements",
            100_000_000,
        )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError(
                "limits must be EvaluationLimits.",
                code="INVALID_FRECHET_DISTANCE_CONFIG",
            )


@dataclass(frozen=True, slots=True)
class FragmentSimilarityConfig:
    """Configure the DNA adaptation of the MOSES Frag metric."""

    k: int = 3
    canonical: bool = True
    ambiguity_policy: KmerAmbiguityPolicy = "ignore"
    max_kmer_observations: int = 10_000_000
    show_progress: bool = True
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        _positive_int(self.k, "k", 32)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError(
                "canonical must be boolean.",
                code="INVALID_FRAGMENT_SIMILARITY_CONFIG",
            )
        if self.ambiguity_policy not in {"error", "ignore"}:
            raise ConfigurationError(
                "ambiguity_policy must be 'error' or 'ignore'.",
                code="INVALID_FRAGMENT_SIMILARITY_CONFIG",
            )
        _positive_int(
            self.max_kmer_observations,
            "max_kmer_observations",
            100_000_000,
        )
        if not isinstance(self.show_progress, bool):
            raise ConfigurationError(
                "show_progress must be boolean.",
                code="INVALID_FRAGMENT_SIMILARITY_CONFIG",
            )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError(
                "limits must be EvaluationLimits.",
                code="INVALID_FRAGMENT_SIMILARITY_CONFIG",
            )


@dataclass(frozen=True, slots=True)
class SNNConfig:
    """Configure similarity to nearest neighbor over DNA bit fingerprints."""

    k: int = 7
    n_bits: int = 1_024
    canonical: bool = True
    seed: int = 0
    ambiguity_policy: KmerAmbiguityPolicy = "ignore"
    max_fingerprint_elements: int = 10_000_000
    show_progress: bool = True
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        _positive_int(self.k, "k", 32)
        _positive_int(self.n_bits, "n_bits", 1_000_000)
        if not isinstance(self.canonical, bool):
            raise ConfigurationError(
                "canonical must be boolean.",
                code="INVALID_SNN_CONFIG",
            )
        if (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, int)
            or not 0 <= self.seed < 2**64
        ):
            raise ConfigurationError(
                "seed must be an integer in [0, 2**64).",
                code="INVALID_SNN_CONFIG",
            )
        if self.ambiguity_policy not in {"error", "ignore"}:
            raise ConfigurationError(
                "ambiguity_policy must be 'error' or 'ignore'.",
                code="INVALID_SNN_CONFIG",
            )
        _positive_int(
            self.max_fingerprint_elements,
            "max_fingerprint_elements",
            100_000_000,
        )
        if not isinstance(self.show_progress, bool):
            raise ConfigurationError(
                "show_progress must be boolean.",
                code="INVALID_SNN_CONFIG",
            )
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError(
                "limits must be EvaluationLimits.",
                code="INVALID_SNN_CONFIG",
            )


@dataclass(frozen=True, slots=True)
class SynthesisRiskConfig:
    global_gc_min: float = 0.30
    global_gc_max: float = 0.70
    local_gc_min: float = 0.20
    local_gc_max: float = 0.80
    window_size: int = 50
    window_step: int = 10
    homopolymer_threshold: int = 8
    tandem_min_unit: int = 2
    tandem_max_unit: int = 20
    tandem_min_repeats: int = 4
    inverted_min_arm: int = 10
    inverted_max_arm: int = 30
    inverted_max_loop: int = 80
    max_windows_per_sequence: int = 1_000_000
    max_pattern_comparisons_per_sequence: int = 2_000_000
    max_matches_per_sequence: int = 100_000
    limits: EvaluationLimits = field(default_factory=EvaluationLimits)

    def __post_init__(self) -> None:
        for name in ("global_gc_min", "global_gc_max", "local_gc_min", "local_gc_max"):
            object.__setattr__(self, name, _probability(getattr(self, name), name))
        if self.global_gc_min > self.global_gc_max or self.local_gc_min > self.local_gc_max:
            raise ConfigurationError("Minimum GC thresholds cannot exceed maximum thresholds.")
        for name, maximum in (
            ("window_size", 1_000_000),
            ("window_step", 1_000_000),
            ("homopolymer_threshold", 1_000_000),
            ("tandem_min_unit", 100),
            ("tandem_max_unit", 100),
            ("tandem_min_repeats", 1_000_000),
            ("inverted_min_arm", 1_000),
            ("inverted_max_arm", 1_000),
            ("max_windows_per_sequence", 10_000_000),
            ("max_pattern_comparisons_per_sequence", 100_000_000),
            ("max_matches_per_sequence", 1_000_000),
        ):
            _positive_int(getattr(self, name), name, maximum)
        _non_negative_int(self.inverted_max_loop, "inverted_max_loop", 1_000_000)
        if self.tandem_min_unit > self.tandem_max_unit:
            raise ConfigurationError("tandem_min_unit cannot exceed tandem_max_unit.")
        if self.tandem_min_repeats < 2:
            raise ConfigurationError("tandem_min_repeats must be at least 2.")
        if self.inverted_min_arm > self.inverted_max_arm:
            raise ConfigurationError("inverted_min_arm cannot exceed inverted_max_arm.")
        if not isinstance(self.limits, EvaluationLimits):
            raise ConfigurationError("limits must be EvaluationLimits.")


@dataclass(frozen=True, slots=True)
class ScoreRule:
    """One explicit min-max normalization rule for a scorecard component."""

    direction: ScoreDirection = "higher_is_better"
    weight: float = 1.0
    minimum: float = 0.0
    maximum: float = 1.0
    metric: str | None = None

    def __post_init__(self) -> None:
        if self.direction not in {"higher_is_better", "lower_is_better"}:
            raise ConfigurationError("Unknown score direction.")
        weight = _finite(self.weight, "weight")
        minimum = _finite(self.minimum, "minimum")
        maximum = _finite(self.maximum, "maximum")
        if weight <= 0 or minimum >= maximum:
            raise ConfigurationError("ScoreRule requires weight > 0 and minimum < maximum.")
        if self.metric is not None and (
            not isinstance(self.metric, str) or not self.metric.strip()
        ):
            raise ConfigurationError("ScoreRule metric must be non-empty or None.")
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "minimum", minimum)
        object.__setattr__(self, "maximum", maximum)


@dataclass(frozen=True, slots=True)
class ScorecardConfig:
    rules: Mapping[str, ScoreRule]
    missing_policy: MissingScorePolicy = "error"
    warning_score: float = 0.5
    pass_score: float = 0.8

    def __post_init__(self) -> None:
        if not isinstance(self.rules, Mapping) or not self.rules:
            raise ConfigurationError("rules must be a non-empty mapping.")
        rules: dict[str, ScoreRule] = {}
        for name, rule in self.rules.items():
            if not isinstance(name, str) or not name.strip() or not isinstance(rule, ScoreRule):
                raise ConfigurationError(
                    "Scorecard rules require non-empty names and ScoreRule values."
                )
            rules[name] = rule
        object.__setattr__(self, "rules", MappingProxyType(rules))
        if self.missing_policy not in {"error", "omit", "zero"}:
            raise ConfigurationError("Unknown missing_policy.")
        warning = _probability(self.warning_score, "warning_score")
        passing = _probability(self.pass_score, "pass_score")
        if warning > passing:
            raise ConfigurationError("warning_score cannot exceed pass_score.")
        object.__setattr__(self, "warning_score", warning)
        object.__setattr__(self, "pass_score", passing)


__all__ = [
    "AmbiguityEvaluationConfig",
    "ComplexityEvaluationConfig",
    "DistributionEvaluationConfig",
    "DiversityEvaluationConfig",
    "EvaluationLimits",
    "FragmentSimilarityConfig",
    "FrechetDistanceConfig",
    "QualityEvaluationConfig",
    "ReferenceSearchConfig",
    "SNNConfig",
    "ScoreRule",
    "ScorecardConfig",
    "SynthesisRiskConfig",
    "UniquenessEvaluationConfig",
]
