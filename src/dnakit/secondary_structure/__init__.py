"""Secondary-structure parsing and explicit NUPACK integration."""

from .dotbracket import (
    analyze_dot_bracket,
    ensemble_defect_from_probabilities,
    pair_probability_metrics,
    target_structure_probability,
)
from .nupack import NupackAdapter, probe_nupack
from .results import (
    AccessibilityWindow,
    ComplexConcentration,
    NupackComplexResult,
    NupackTubeResult,
    PairProbabilityResult,
    PredictedSecondaryStructure,
    SecondaryBasePair,
    SecondaryStem,
    SecondaryStructureSummary,
)

__all__ = [
    "AccessibilityWindow",
    "ComplexConcentration",
    "NupackAdapter",
    "NupackComplexResult",
    "NupackTubeResult",
    "PairProbabilityResult",
    "PredictedSecondaryStructure",
    "SecondaryBasePair",
    "SecondaryStem",
    "SecondaryStructureSummary",
    "analyze_dot_bracket",
    "ensemble_defect_from_probabilities",
    "pair_probability_metrics",
    "probe_nupack",
    "target_structure_probability",
]
