"""Transparent, deterministic DNA sequence and dataset evaluation APIs."""

from dnakit.evaluation.collection import (
    evaluate_diversity,
    evaluate_redundancy,
    evaluate_uniqueness,
)
from dnakit.evaluation.config import (
    AmbiguityEvaluationConfig,
    ComplexityEvaluationConfig,
    DistributionEvaluationConfig,
    DiversityEvaluationConfig,
    EvaluationLimits,
    FragmentSimilarityConfig,
    FrechetDistanceConfig,
    QualityEvaluationConfig,
    ReferenceSearchConfig,
    ScorecardConfig,
    ScoreRule,
    SNNConfig,
    SynthesisRiskConfig,
    UniquenessEvaluationConfig,
)
from dnakit.evaluation.distribution import evaluate_distribution_similarity
from dnakit.evaluation.frechet import evaluate_frechet_distance
from dnakit.evaluation.generative import evaluate_fragment_similarity, evaluate_snn
from dnakit.evaluation.reference import (
    create_reference_library,
    evaluate_memorization,
    evaluate_novelty,
    evaluate_reference_similarity,
    nearest_reference,
)
from dnakit.evaluation.results import EvaluationEntry, EvaluationReport, ReferenceLibrary
from dnakit.evaluation.scorecard import ScoreInput, evaluate_scorecard
from dnakit.evaluation.sequence import (
    evaluate_ambiguity,
    evaluate_complexity,
    evaluate_quality,
    evaluate_validity,
)
from dnakit.evaluation.synthesis import evaluate_synthesis_risk

__all__ = [
    "AmbiguityEvaluationConfig",
    "ComplexityEvaluationConfig",
    "DistributionEvaluationConfig",
    "DiversityEvaluationConfig",
    "EvaluationEntry",
    "EvaluationLimits",
    "EvaluationReport",
    "FragmentSimilarityConfig",
    "FrechetDistanceConfig",
    "QualityEvaluationConfig",
    "ReferenceLibrary",
    "ReferenceSearchConfig",
    "SNNConfig",
    "ScoreInput",
    "ScoreRule",
    "ScorecardConfig",
    "SynthesisRiskConfig",
    "UniquenessEvaluationConfig",
    "create_reference_library",
    "evaluate_ambiguity",
    "evaluate_complexity",
    "evaluate_distribution_similarity",
    "evaluate_diversity",
    "evaluate_fragment_similarity",
    "evaluate_frechet_distance",
    "evaluate_memorization",
    "evaluate_novelty",
    "evaluate_quality",
    "evaluate_redundancy",
    "evaluate_reference_similarity",
    "evaluate_scorecard",
    "evaluate_snn",
    "evaluate_synthesis_risk",
    "evaluate_uniqueness",
    "evaluate_validity",
    "nearest_reference",
]
