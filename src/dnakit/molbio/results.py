"""Immutable JSON-serializable molecular-biology simulation results."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal, cast

from dnakit.core import DNASequence, Issue, Provenance, Strand
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError
from dnakit.patterns.results import GuideCandidate

EndPolarity = Literal["blunt", "5prime", "3prime"]
EndSide = Literal["left", "right"]


class _Serializable:
    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, kw_only=True)
class MolBioResult(_Serializable):
    """Common audit envelope carried directly by every public result."""

    method: str
    algorithm_version: str
    parameters: FrozenDict
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __post_init__(self) -> None:
        for name in ("method", "algorithm_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"{name} must be a non-empty string.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("parameters must be a FrozenDict.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("provenance must be Provenance.")
        if any(not isinstance(issue, Issue) for issue in self.issues):
            raise ConfigurationError("issues must contain Issue objects.")


@dataclass(frozen=True)
class EndDescriptor(_Serializable):
    """One abstract double-stranded fragment end in top-strand coordinates."""

    polarity: EndPolarity
    overhang_sequence_5to3: str
    cohesive_key: str
    side: EndSide
    five_prime_phosphorylated: bool
    top_cut: int | None
    bottom_cut: int | None
    source: str

    def __post_init__(self) -> None:
        if self.polarity not in ("blunt", "5prime", "3prime"):
            raise ConfigurationError("Unknown end polarity.")
        if self.side not in ("left", "right"):
            raise ConfigurationError("End side must be left or right.")
        if set(self.overhang_sequence_5to3) - set("ACGTRYSWKMBDHVN"):
            raise ConfigurationError("End overhang must contain uppercase DNA IUPAC symbols.")
        if self.polarity == "blunt" and self.overhang_sequence_5to3:
            raise ConfigurationError("A blunt end cannot carry an overhang sequence.")
        if self.polarity != "blunt" and not self.overhang_sequence_5to3:
            raise ConfigurationError("A sticky end requires an overhang sequence.")
        if set(self.cohesive_key) - set("ACGTRYSWKMBDHVN"):
            raise ConfigurationError("End cohesive_key must contain uppercase DNA IUPAC symbols.")
        if self.polarity == "blunt" and self.cohesive_key:
            raise ConfigurationError("A blunt end must use an empty cohesive_key.")
        if self.polarity != "blunt" and not self.cohesive_key:
            raise ConfigurationError("A sticky end requires a cohesive_key.")
        if not isinstance(self.five_prime_phosphorylated, bool):
            raise ConfigurationError("five_prime_phosphorylated must be boolean.")
        for name in ("top_cut", "bottom_cut"):
            value = getattr(self, name)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int) or value < 0
            ):
                raise ConfigurationError(f"{name} must be a non-negative integer or None.")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ConfigurationError("End source must be non-empty.")


@dataclass(frozen=True)
class DigestCut(_Serializable):
    enzymes: tuple[str, ...]
    top_cut: int
    bottom_cut: int
    polarity: EndPolarity
    overhang_sequence_5to3: str


@dataclass(frozen=True)
class DigestFragment(_Serializable):
    id: str
    sequence: DNASequence
    source_start: int
    source_end: int
    wraps_origin: bool
    left_end: EndDescriptor
    right_end: EndDescriptor


@dataclass(frozen=True)
class RestrictionDigestResult(MolBioResult):
    source_length: int
    source_topology: str
    cuts: tuple[DigestCut, ...]
    fragments: tuple[DigestFragment, ...]
    coordinate_system: str = "0-based-half-open"


@dataclass(frozen=True)
class EndTypeResult(MolBioResult):
    end: EndDescriptor


@dataclass(frozen=True)
class LigationCompatibilityResult(MolBioResult):
    compatible: bool
    reason: str
    left: EndDescriptor
    right: EndDescriptor


@dataclass(frozen=True)
class LigationResult(MolBioResult):
    product: DNASequence
    fragment_ids: tuple[str, ...]
    junction_count: int
    circularized: bool


@dataclass(frozen=True)
class PrimerBindingHit(_Serializable):
    primer_name: str
    strand: Strand
    start: int
    end: int
    matched_template: str
    mismatch_positions_5to3: tuple[int, ...]
    three_prime_mismatch_count: int
    wraps_origin: bool

    @property
    def mismatch_count(self) -> int:
        return len(self.mismatch_positions_5to3)


@dataclass(frozen=True)
class PrimerMatchResult(MolBioResult):
    primer_sequence: str
    template_length: int
    template_topology: str
    hits: tuple[PrimerBindingHit, ...]
    truncated: bool
    coordinate_system: str = "0-based-half-open"


@dataclass(frozen=True)
class Amplicon(_Serializable):
    sequence: DNASequence
    template_start: int
    template_end: int
    wraps_origin: bool
    forward_binding: PrimerBindingHit
    reverse_binding: PrimerBindingHit


@dataclass(frozen=True)
class PCRResult(MolBioResult):
    template_length: int
    template_topology: str
    amplicons: tuple[Amplicon, ...]
    truncated: bool
    coordinate_system: str = "0-based-half-open"


@dataclass(frozen=True)
class ConditionalAnalysis(_Serializable):
    capability: str
    status: str
    available: bool
    backend_required: str
    automatic_probe: bool
    automatic_execution: bool
    reason: str
    execution_performed: bool = False
    structure_found: bool | None = None
    tm_celsius: float | None = None
    delta_g_kcal_per_mol: float | None = None
    delta_h_kcal_per_mol: float | None = None
    delta_s_cal_per_k_mol: float | None = None
    ascii_structure: str | None = None
    backend_name: str | None = None
    backend_version: str | None = None

    def __post_init__(self) -> None:
        for name in ("capability", "status", "backend_required", "reason"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Conditional analysis {name} must be non-empty.")
        for name in (
            "available",
            "automatic_probe",
            "automatic_execution",
            "execution_performed",
        ):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"Conditional analysis {name} must be boolean.")
        if self.available != self.execution_performed:
            raise ConfigurationError(
                "Conditional analysis availability must match explicit execution status."
            )
        if self.structure_found is not None and not isinstance(self.structure_found, bool):
            raise ConfigurationError("structure_found must be boolean or None.")
        values = (
            self.tm_celsius,
            self.delta_g_kcal_per_mol,
            self.delta_h_kcal_per_mol,
            self.delta_s_cal_per_k_mol,
        )
        if any(value is not None and not math.isfinite(value) for value in values):
            raise ConfigurationError("Structure thermodynamic values must be finite or None.")
        if self.available and any(value is None for value in values):
            raise ConfigurationError("An available structure analysis requires all numeric values.")
        if not self.available and any(value is not None for value in values):
            raise ConfigurationError("An unavailable structure analysis cannot contain values.")
        if self.available and self.structure_found is None:
            raise ConfigurationError("An available structure analysis requires structure_found.")
        if not self.available and (
            self.structure_found is not None
            or self.ascii_structure is not None
            or self.backend_name is not None
            or self.backend_version is not None
        ):
            raise ConfigurationError("An unavailable structure analysis cannot claim backend data.")
        if self.available and (
            not isinstance(self.backend_name, str) or not self.backend_name.strip()
        ):
            raise ConfigurationError("An available structure analysis requires backend_name.")
        for name in ("backend_name", "backend_version"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(f"{name} must be a non-empty string or None.")
        if self.ascii_structure is not None and (
            not isinstance(self.ascii_structure, str) or len(self.ascii_structure) > 100_000
        ):
            raise ConfigurationError("ascii_structure must be text within 100,000 characters.")


@dataclass(frozen=True)
class PrimerPropertiesResult(MolBioResult):
    primer_sequence: str
    gc_fraction: float
    tm_celsius: float
    tm_method: str
    hairpin: ConditionalAnalysis
    self_dimer: ConditionalAnalysis
    heterodimer: ConditionalAnalysis | None
    paired_primer_sequence: str | None
    paired_gc_fraction: float | None
    paired_tm_celsius: float | None
    paired_hairpin: ConditionalAnalysis | None
    paired_self_dimer: ConditionalAnalysis | None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not math.isfinite(self.gc_fraction) or not 0.0 <= self.gc_fraction <= 1.0:
            raise ConfigurationError("gc_fraction must be finite and between 0 and 1.")
        if not math.isfinite(self.tm_celsius):
            raise ConfigurationError("tm_celsius must be finite.")
        missing_pair = self.paired_primer_sequence is None
        if missing_pair != (self.paired_gc_fraction is None) or missing_pair != (
            self.paired_tm_celsius is None
        ):
            raise ConfigurationError("Paired primer sequence, GC, and Tm must be set together.")
        if self.paired_gc_fraction is not None and (
            not math.isfinite(self.paired_gc_fraction) or not 0.0 <= self.paired_gc_fraction <= 1.0
        ):
            raise ConfigurationError("paired_gc_fraction must be finite and between 0 and 1.")
        if self.paired_tm_celsius is not None and not math.isfinite(self.paired_tm_celsius):
            raise ConfigurationError("paired_tm_celsius must be finite.")
        if missing_pair != (self.paired_hairpin is None) or missing_pair != (
            self.paired_self_dimer is None
        ):
            raise ConfigurationError(
                "Paired hairpin and self-dimer analyses must follow paired primer presence."
            )
        if missing_pair != (self.heterodimer is None):
            raise ConfigurationError("Heterodimer analysis must follow paired primer presence.")
        analyses = (
            self.hairpin,
            self.self_dimer,
            self.heterodimer,
            self.paired_hairpin,
            self.paired_self_dimer,
        )
        if any(
            analysis is not None and not isinstance(analysis, ConditionalAnalysis)
            for analysis in analyses
        ):
            raise ConfigurationError("Primer structure fields must be ConditionalAnalysis values.")
        expected_capabilities = (
            (self.hairpin, "hairpin"),
            (self.self_dimer, "self_dimer"),
            (self.heterodimer, "heterodimer"),
            (self.paired_hairpin, "hairpin"),
            (self.paired_self_dimer, "self_dimer"),
        )
        if any(
            analysis is not None and analysis.capability != capability
            for analysis, capability in expected_capabilities
        ):
            raise ConfigurationError("Primer structure capability labels are inconsistent.")


@dataclass(frozen=True)
class PrimerDesignRequest(_Serializable):
    template: DNASequence
    target_start: int
    target_end: int
    primer_length_range: tuple[int, int]
    tm_range_celsius: tuple[float, float]
    gc_range: tuple[float, float]
    product_length_range: tuple[int, int]
    excluded_regions: tuple[tuple[int, int], ...]
    candidate_count: int
    thermodynamic_conditions: FrozenDict

    def __post_init__(self) -> None:
        if not isinstance(self.template, DNASequence):
            raise ConfigurationError("Primer design template must be DNASequence.")
        length = self.template.symbol_length
        if not 0 <= self.target_start < self.target_end <= length:
            raise ConfigurationError("Primer design target interval is invalid.")
        for name, integer_bounds in (
            ("primer_length_range", self.primer_length_range),
            ("product_length_range", self.product_length_range),
        ):
            if (
                len(integer_bounds) != 2
                or any(
                    isinstance(item, bool) or not isinstance(item, int) for item in integer_bounds
                )
                or integer_bounds[0] <= 0
                or integer_bounds[0] > integer_bounds[1]
            ):
                raise ConfigurationError(f"{name} is invalid.")
        for name, float_bounds in (
            ("tm_range_celsius", self.tm_range_celsius),
            ("gc_range", self.gc_range),
        ):
            if (
                len(float_bounds) != 2
                or any(not math.isfinite(item) for item in float_bounds)
                or float_bounds[0] > float_bounds[1]
            ):
                raise ConfigurationError(f"{name} is invalid.")
        if not 0.0 <= self.gc_range[0] <= self.gc_range[1] <= 1.0:
            raise ConfigurationError("gc_range must lie between zero and one.")
        for start, end in self.excluded_regions:
            if not 0 <= start < end <= length:
                raise ConfigurationError("An excluded primer-design interval is invalid.")
        if (
            isinstance(self.candidate_count, bool)
            or not isinstance(self.candidate_count, int)
            or self.candidate_count <= 0
        ):
            raise ConfigurationError("candidate_count must be positive.")
        if not isinstance(self.thermodynamic_conditions, FrozenDict):
            raise ConfigurationError("thermodynamic_conditions must be FrozenDict.")


@dataclass(frozen=True)
class PrimerDesignCandidate(_Serializable):
    """One bounded Primer3 primer-pair candidate in template coordinates."""

    rank: int
    left_primer_sequence: str
    right_primer_sequence: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    product_size: int
    left_tm_celsius: float
    right_tm_celsius: float
    left_gc_fraction: float
    right_gc_fraction: float
    pair_penalty: float | None

    def __post_init__(self) -> None:
        if isinstance(self.rank, bool) or not isinstance(self.rank, int) or self.rank < 0:
            raise ConfigurationError("Primer candidate rank must be a non-negative integer.")
        for name in ("left_primer_sequence", "right_primer_sequence"):
            value = getattr(self, name)
            if (
                not isinstance(value, str)
                or not value
                or len(value) > 1_000
                or set(value) - set("ACGT")
            ):
                raise ConfigurationError(f"{name} must be bounded canonical DNA.")
        coordinates = (self.left_start, self.left_end, self.right_start, self.right_end)
        if any(isinstance(value, bool) or not isinstance(value, int) for value in coordinates):
            raise ConfigurationError("Primer candidate coordinates must be integers.")
        if not 0 <= self.left_start < self.left_end <= self.right_start < self.right_end:
            raise ConfigurationError("Primer candidate intervals are invalid or overlap.")
        if (
            isinstance(self.product_size, bool)
            or not isinstance(self.product_size, int)
            or self.product_size <= 0
            or self.product_size != self.right_end - self.left_start
        ):
            raise ConfigurationError("Primer candidate product_size is inconsistent.")
        for name in ("left_tm_celsius", "right_tm_celsius"):
            if not math.isfinite(getattr(self, name)):
                raise ConfigurationError(f"{name} must be finite.")
        for name in ("left_gc_fraction", "right_gc_fraction"):
            value = getattr(self, name)
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ConfigurationError(f"{name} must be finite and between zero and one.")
        if self.pair_penalty is not None and not math.isfinite(self.pair_penalty):
            raise ConfigurationError("pair_penalty must be finite or None.")


@dataclass(frozen=True)
class PrimerDesignInterfaceResult(MolBioResult):
    request: PrimerDesignRequest
    backend_name: str
    status: str
    execution_performed: bool
    candidates: tuple[PrimerDesignCandidate, ...]
    reason: str

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.request, PrimerDesignRequest):
            raise ConfigurationError("Primer design result requires a PrimerDesignRequest.")
        for name in ("backend_name", "status", "reason"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Primer design {name} must be non-empty.")
        if not isinstance(self.execution_performed, bool):
            raise ConfigurationError("execution_performed must be boolean.")
        if any(not isinstance(candidate, PrimerDesignCandidate) for candidate in self.candidates):
            raise ConfigurationError("candidates must contain PrimerDesignCandidate values.")
        if len(self.candidates) > self.request.candidate_count:
            raise ConfigurationError("Candidate output exceeds the requested candidate count.")
        if self.execution_performed != self.status.startswith("execution-complete"):
            raise ConfigurationError("Primer design status disagrees with execution_performed.")


@dataclass(frozen=True)
class AssemblyStep(_Serializable):
    index: int
    left_fragment_id: str
    right_fragment_id: str
    operation: str
    junction_sequence: str


@dataclass(frozen=True)
class AssemblyResult(MolBioResult):
    assembly_method: str
    product: DNASequence
    steps: tuple[AssemblyStep, ...]
    circularized: bool
    complete: bool


@dataclass(frozen=True)
class OffTargetHit(_Serializable):
    guide_sequence: str
    reference_id: str
    strand: Strand
    start: int
    end: int
    mismatch_positions_5to3: tuple[int, ...]
    wraps_origin: bool


@dataclass(frozen=True)
class CrisprScanResult(MolBioResult):
    candidates: tuple[GuideCandidate, ...]
    off_targets: tuple[OffTargetHit, ...]
    candidate_truncated: bool
    off_target_truncated: bool
    efficiency_prediction_performed: bool


@dataclass(frozen=True)
class SequenceChange(_Serializable):
    position: int
    before: str
    after: str
    reason: str


@dataclass(frozen=True)
class RuleOptimizationResult(MolBioResult):
    original: DNASequence
    optimized: DNASequence
    changes: tuple[SequenceChange, ...]
    initial_score: float
    final_score: float
    iterations: int
    constraints_satisfied: bool


@dataclass(frozen=True)
class CodonOptimizationResult(MolBioResult):
    original: DNASequence
    optimized: DNASequence
    original_translation: str
    optimized_translation: str
    changes: tuple[SequenceChange, ...]
    original_cai: float
    optimized_cai: float
    original_gc_fraction: float
    optimized_gc_fraction: float
    usage_table_name: str
    usage_table_version: str
    usage_table_checksum: str


@dataclass(frozen=True)
class Mutation(_Serializable):
    position: int
    reference: str
    alternate: str


@dataclass(frozen=True)
class SequenceVariant(_Serializable):
    id: str
    sequence: DNASequence
    mutations: tuple[Mutation, ...]


@dataclass(frozen=True)
class MutationLibraryResult(MolBioResult):
    source: DNASequence
    variants: tuple[SequenceVariant, ...]
    total_possible_variants: int
    sampled: bool
    seed: int


__all__ = [
    "Amplicon",
    "AssemblyResult",
    "AssemblyStep",
    "CodonOptimizationResult",
    "ConditionalAnalysis",
    "CrisprScanResult",
    "DigestCut",
    "DigestFragment",
    "EndDescriptor",
    "EndTypeResult",
    "LigationCompatibilityResult",
    "LigationResult",
    "MolBioResult",
    "Mutation",
    "MutationLibraryResult",
    "OffTargetHit",
    "PCRResult",
    "PrimerBindingHit",
    "PrimerDesignCandidate",
    "PrimerDesignInterfaceResult",
    "PrimerDesignRequest",
    "PrimerMatchResult",
    "PrimerPropertiesResult",
    "RestrictionDigestResult",
    "RuleOptimizationResult",
    "SequenceChange",
    "SequenceVariant",
]
