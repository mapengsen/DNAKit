"""Public DNA normalization and validation API."""

from dnakit._processing.dispatch import register as _register_inputs
from dnakit._verbosity import register_public as _register_public_verbose

from .config import (
    AmbiguityPolicy,
    DatasetValidationConfig,
    NormalizationConfig,
    UPolicy,
    ValidationConfig,
)
from .gaps import (
    AGPAssemblyConfig,
    AGPAssemblyResult,
    AGPAssemblySegment,
    AGPComponentValue,
    CircularBoundaryPolicy,
    GapNormalizable,
    GapNormalizationChange,
    GapNormalizationConfig,
    GapNormalizationResult,
    normalize_gaps,
    sequence_from_agp,
)
from .normalize import RawSequenceInput, normalize
from .results import (
    AmbiguityReport,
    DatasetValidationReport,
    DuplicateID,
    InputPosition,
    InvalidSymbol,
    NormalizationChange,
    NormalizationResult,
    NormalizationStep,
    ProbabilityResolution,
    QualitySummary,
    RawInputSnapshot,
    SymbolCount,
    SymbolOccurrence,
    ValidationReport,
)
from .validate import validate, validate_set

_register_inputs(globals(), ("normalize",), mode="normalize")
_register_inputs(globals(), ("normalize_gaps",), mode="sequence")
_register_inputs(globals(), ("validate",), mode="report")
_register_inputs(globals(), ("validate_set",), mode="collection")

__all__ = [
    "AGPAssemblyConfig",
    "AGPAssemblyResult",
    "AGPAssemblySegment",
    "AGPComponentValue",
    "AmbiguityPolicy",
    "AmbiguityReport",
    "CircularBoundaryPolicy",
    "DatasetValidationConfig",
    "DatasetValidationReport",
    "DuplicateID",
    "GapNormalizable",
    "GapNormalizationChange",
    "GapNormalizationConfig",
    "GapNormalizationResult",
    "InputPosition",
    "InvalidSymbol",
    "NormalizationChange",
    "NormalizationConfig",
    "NormalizationResult",
    "NormalizationStep",
    "ProbabilityResolution",
    "QualitySummary",
    "RawInputSnapshot",
    "RawSequenceInput",
    "SymbolCount",
    "SymbolOccurrence",
    "UPolicy",
    "ValidationConfig",
    "ValidationReport",
    "normalize",
    "normalize_gaps",
    "sequence_from_agp",
    "validate",
    "validate_set",
]

_register_public_verbose(globals())
