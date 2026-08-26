"""Immutable secondary-structure and ensemble result objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import BackendInfo, Provenance
from dnakit.core._json import to_json_compatible


class _SerializableResult:
    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class SecondaryBasePair(_SerializableResult):
    first_global_index: int
    second_global_index: int
    first_strand_index: int
    first_local_index: int
    second_strand_index: int
    second_local_index: int
    bracket_type: str
    inter_strand: bool


@dataclass(frozen=True)
class SecondaryStem(_SerializableResult):
    base_pairs: tuple[SecondaryBasePair, ...]
    length: int
    inter_strand: bool


@dataclass(frozen=True)
class SecondaryStructureSummary(_SerializableResult):
    strands_5to3: tuple[str, ...]
    dot_bracket: str
    base_pairs: tuple[SecondaryBasePair, ...]
    stems: tuple[SecondaryStem, ...]
    structure_type: str
    base_pair_count: int
    paired_base_fraction: float
    stem_lengths: tuple[int, ...]
    hairpin_count: int
    hairpin_loop_lengths: tuple[int, ...]
    max_contiguous_pair_count: int
    three_prime_window: int
    three_prime_dimer: bool
    three_prime_dimer_max_contiguous_pairs: int
    method: str


@dataclass(frozen=True)
class AccessibilityWindow(_SerializableResult):
    start: int
    end: int
    mean_unpaired_probability: float


@dataclass(frozen=True)
class PairProbabilityResult(_SerializableResult):
    strands_5to3: tuple[str, ...]
    pair_probabilities: tuple[tuple[float, ...], ...]
    pairing_probabilities_by_base: tuple[float, ...]
    unpaired_probabilities_by_base: tuple[float, ...]
    accessibility_window_size: int
    accessibility_windows: tuple[AccessibilityWindow, ...]
    most_accessible_window_start: int | None
    method: str
    applicability: str


@dataclass(frozen=True)
class PredictedSecondaryStructure(_SerializableResult):
    summary: SecondaryStructureSummary
    free_energy_kcal_per_mol: float
    stack_free_energy_kcal_per_mol: float | None


@dataclass(frozen=True)
class NupackComplexResult(_SerializableResult):
    strands_5to3: tuple[str, ...]
    material: str
    temperature_celsius: float
    monovalent_molar: float
    magnesium_molar: float
    ensemble: str
    partition_function_log: float
    ensemble_free_energy_kcal_per_mol: float
    mfe_structures: tuple[PredictedSecondaryStructure, ...]
    pair_probabilities: PairProbabilityResult
    suboptimal_structures: tuple[PredictedSecondaryStructure, ...]
    boltzmann_samples: tuple[str, ...]
    ensemble_size: int
    target_structure: SecondaryStructureSummary | None
    target_structure_probability: float | None
    target_ensemble_defect: float | None
    method: str
    backend: BackendInfo
    provenance: Provenance


@dataclass(frozen=True)
class ComplexConcentration(_SerializableResult):
    name: str
    strand_names: tuple[str, ...]
    concentration_molar: float
    is_target: bool


@dataclass(frozen=True)
class NupackTubeResult(_SerializableResult):
    strand_names: tuple[str, ...]
    sequences_5to3: tuple[str, ...]
    input_concentrations_molar: tuple[float, ...]
    complex_concentrations: tuple[ComplexConcentration, ...]
    target_complex_name: str
    target_strand_names: tuple[str, ...]
    target_complex_concentration_molar: float
    complex_fraction_denominator_molar: float
    target_complex_fraction: float
    non_target_complex_fraction: float
    fraction_bases_unpaired: float
    max_complex_size: int
    method: str
    backend: BackendInfo
    provenance: Provenance


__all__ = [
    "AccessibilityWindow",
    "ComplexConcentration",
    "NupackComplexResult",
    "NupackTubeResult",
    "PairProbabilityResult",
    "PredictedSecondaryStructure",
    "SecondaryBasePair",
    "SecondaryStem",
    "SecondaryStructureSummary",
]
