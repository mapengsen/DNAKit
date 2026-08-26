"""Immutable, serializable thermodynamic results with explicit units."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import BackendInfo, ImplementationLabel, Issue, Provenance
from dnakit.core._json import FrozenDict, to_json_compatible

from .config import ThermodynamicConditions


class _SerializableResult:
    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class MolecularWeightResult(_SerializableResult):
    value_dalton: float
    value_kilodalton: float
    strand_count: int
    sequence_length: int
    five_prime_phosphorylated: bool
    method: str
    algorithm_version: str
    applicability: str
    parameters: FrozenDict
    provenance: Provenance


@dataclass(frozen=True)
class ExtinctionCoefficientResult(_SerializableResult):
    value_m_inverse_cm_inverse: float
    wavelength_nm: int
    sequence_length: int
    method: str
    algorithm_version: str
    applicability: str
    parameters: FrozenDict
    provenance: Provenance


@dataclass(frozen=True)
class SaltCorrectionResult(_SerializableResult):
    delta_s_cal_per_k_mol: float
    sequence_length: int
    model: str
    model_version: str
    conditions: ThermodynamicConditions
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class StackingStep(_SerializableResult):
    index: int
    top_5to3: str
    bottom_3to5: str
    delta_h_kcal_per_mol: float
    delta_s_cal_per_k_mol: float
    delta_g_kcal_per_mol: float


@dataclass(frozen=True)
class StackingResult(_SerializableResult):
    sequence: str
    complement_3to5: str
    temperature_celsius: float
    parameter_set: str
    steps: tuple[StackingStep, ...]
    total_delta_h_kcal_per_mol: float
    total_delta_s_cal_per_k_mol: float
    total_delta_g_kcal_per_mol: float
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class NearestNeighborResult(_SerializableResult):
    sequence: str
    complement_5to3: str
    sequence_length: int
    parameter_set: str
    conditions: ThermodynamicConditions
    stacking_steps: tuple[StackingStep, ...]
    stacking_delta_h_kcal_per_mol: float
    stacking_delta_s_cal_per_k_mol: float
    initiation_delta_h_kcal_per_mol: float
    initiation_delta_s_cal_per_k_mol: float
    symmetry_delta_h_kcal_per_mol: float
    symmetry_delta_s_cal_per_k_mol: float
    salt_delta_s_cal_per_k_mol: float
    delta_h_kcal_per_mol: float
    delta_s_cal_per_k_mol: float
    delta_g_kcal_per_mol: float
    tm_celsius: float
    self_complementary: bool
    concentration_divisor: int
    gas_constant_cal_per_k_mol: float
    reference_sodium_molar: float
    tm_equation: str
    method: str
    algorithm_version: str
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class MeltingTemperatureResult(_SerializableResult):
    tm_celsius: float
    sequence_length: int
    method: str
    algorithm_version: str
    conditions: ThermodynamicConditions
    parameter_set: str | None
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class DuplexStabilityResult(_SerializableResult):
    sequence_a_5to3: str
    sequence_b_5to3: str
    fully_complementary: bool
    stable_at_temperature: bool
    stability_criterion: str
    delta_g_kcal_per_mol: float
    tm_celsius: float
    conditions: ThermodynamicConditions
    model: str
    applicability: str
    thermodynamics: NearestNeighborResult | Primer3ThermodynamicResult
    provenance: Provenance


@dataclass(frozen=True)
class WindowTmPoint(_SerializableResult):
    start: int
    end: int
    sequence: str
    tm_celsius: float


@dataclass(frozen=True)
class WindowTmResult(_SerializableResult):
    sequence_length: int
    window_size: int
    step: int
    method: str
    conditions: ThermodynamicConditions
    windows: tuple[WindowTmPoint, ...]
    min_tm_celsius: float | None
    max_tm_celsius: float | None
    max_windows: int
    coordinate_system: str
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class ConditionalCapability(_SerializableResult):
    requirement_id: str
    capability: str
    status: str
    execution_supported: bool
    automatic_install: bool
    automatic_download: bool
    compatible_backends: tuple[str, ...]
    backend_info: BackendInfo | None
    reason: str


@dataclass(frozen=True)
class Primer3ThermodynamicResult(_SerializableResult):
    """One explicitly requested calculation from user-installed Primer3 CLI tools."""

    capability: str
    sequences_5to3: tuple[str, ...]
    structure_found: bool | None
    tm_celsius: float
    delta_g_kcal_per_mol: float | None
    delta_h_kcal_per_mol: float | None
    delta_s_cal_per_k_mol: float | None
    ascii_structure: str | None
    conditions: ThermodynamicConditions
    max_loop: int | None
    method: str
    algorithm_version: str
    backend: BackendInfo
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __post_init__(self) -> None:
        import math

        if self.capability not in {"tm", "hairpin", "self_dimer", "heterodimer"}:
            raise ValueError("Unknown Primer3 thermodynamic capability.")
        expected_sequences = 2 if self.capability == "heterodimer" else 1
        if len(self.sequences_5to3) != expected_sequences or any(
            not sequence or len(sequence) > 60 or set(sequence) - set("ACGT")
            for sequence in self.sequences_5to3
        ):
            raise ValueError("Primer3 sequences must contain 1-60 nt of canonical DNA.")
        if self.structure_found is not None and not isinstance(self.structure_found, bool):
            raise ValueError("structure_found must be boolean or None.")
        if not math.isfinite(self.tm_celsius):
            raise ValueError("tm_celsius must be finite.")
        thermodynamic_values = (
            self.delta_g_kcal_per_mol,
            self.delta_h_kcal_per_mol,
            self.delta_s_cal_per_k_mol,
        )
        if any(value is not None and not math.isfinite(value) for value in thermodynamic_values):
            raise ValueError("Primer3 thermodynamic values must be finite or None.")
        if self.capability == "tm" and any(value is not None for value in thermodynamic_values):
            raise ValueError("Tm-only results cannot claim structure thermodynamics.")
        if self.capability != "tm" and any(value is None for value in thermodynamic_values):
            raise ValueError("Structure results require dG, dH, and dS values.")
        if self.max_loop is not None and (
            isinstance(self.max_loop, bool)
            or not isinstance(self.max_loop, int)
            or not 1 <= self.max_loop <= 30
        ):
            raise ValueError("max_loop must be None or an integer in [1, 30].")
        if (self.capability == "tm") != (self.max_loop is None):
            raise ValueError(
                "Tm results require max_loop=None; structure results require max_loop."
            )
        if self.ascii_structure is not None and (
            not isinstance(self.ascii_structure, str) or len(self.ascii_structure) > 100_000
        ):
            raise ValueError("ascii_structure must be text within 100,000 characters.")
        for name in ("method", "algorithm_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string.")
        if not isinstance(self.conditions, ThermodynamicConditions):
            raise ValueError("conditions must be ThermodynamicConditions.")
        if not isinstance(self.backend, BackendInfo) or not self.backend.available:
            raise ValueError("Primer3 result requires an available BackendInfo.")
        if self.capability not in self.backend.capabilities:
            raise ValueError("Primer3 result capability is absent from BackendInfo.")
        if not isinstance(self.provenance, Provenance):
            raise ValueError("provenance must be Provenance.")
        if self.provenance.backend != self.backend:
            raise ValueError("Primer3 result backend and provenance backend must match.")
        if self.provenance.implementation.label is not ImplementationLabel.ADAPTER:
            raise ValueError("Primer3 results must identify an adapter implementation.")
        if any(not isinstance(issue, Issue) for issue in self.issues):
            raise ValueError("issues must contain Issue objects.")


@dataclass(frozen=True)
class OpticalModification(_SerializableResult):
    """One explicitly counted label or chemical modification correction."""

    name: str
    count: int = 1
    extinction_coefficient_260_delta_m_inverse_cm_inverse: float = 0.0
    molecular_weight_delta_dalton: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Optical modification name must be non-empty text.")
        if isinstance(self.count, bool) or not isinstance(self.count, int) or self.count <= 0:
            raise ValueError("Optical modification count must be a positive integer.")
        for field_name in (
            "extinction_coefficient_260_delta_m_inverse_cm_inverse",
            "molecular_weight_delta_dalton",
        ):
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError(f"{field_name} must be finite numeric data.")


@dataclass(frozen=True)
class LabelAbsorbanceCorrection(_SerializableResult):
    """Measured dye-channel absorbance contribution to subtract from A260."""

    name: str
    absorbance_at_label_max: float
    a260_correction_factor: float

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Label correction name must be non-empty text.")
        for field_name in ("absorbance_at_label_max", "a260_correction_factor"):
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0.0
            ):
                raise ValueError(f"{field_name} must be finite and non-negative.")


@dataclass(frozen=True)
class OpticalPropertiesResult(_SerializableResult):
    strand_type: str
    sequences_5to3: tuple[str, ...]
    sequence_lengths: tuple[int, ...]
    native_extinction_coefficient_260_m_inverse_cm_inverse: float
    modification_extinction_coefficient_260_m_inverse_cm_inverse: float
    extinction_coefficient_260_m_inverse_cm_inverse: float
    molecular_weight_dalton: float
    one_od260_nmol: float
    one_od260_microgram: float
    method: str
    modifications: tuple[OpticalModification, ...]
    assumptions: tuple[str, ...]
    provenance: Provenance


@dataclass(frozen=True)
class A260ConcentrationResult(_SerializableResult):
    measured_a260: float
    label_a260_subtracted: float
    corrected_a260: float
    path_length_cm: float
    dilution_factor: float
    extinction_coefficient_260_m_inverse_cm_inverse: float
    molar_concentration_molar: float
    molar_concentration_micromolar: float
    molar_concentration_nanomolar: float
    mass_concentration_g_per_l: float
    mass_concentration_ng_per_microliter: float
    volume_liter: float | None
    amount_mol: float | None
    mass_microgram: float | None
    label_corrections: tuple[LabelAbsorbanceCorrection, ...]
    method: str
    provenance: Provenance


@dataclass(frozen=True)
class OligoQuantityResult(_SerializableResult):
    molecular_weight_dalton: float
    volume_liter: float | None
    molar_concentration_molar: float | None
    mass_concentration_g_per_l: float | None
    amount_mol: float
    amount_nmol: float
    mass_g: float
    mass_microgram: float
    input_kind: str
    method: str


@dataclass(frozen=True)
class BindingEquilibriumResult(_SerializableResult):
    temperature_celsius: float
    delta_g_kcal_per_mol: float
    association_constant_m_inverse: float
    dissociation_constant_molar: float
    total_strand_concentration_molar: float
    free_strand_concentration_molar: float
    duplex_concentration_molar: float
    duplex_fraction: float
    self_complementary: bool
    model: str
    thermodynamics: NearestNeighborResult
    provenance: Provenance


@dataclass(frozen=True)
class MeltingCurvePoint(_SerializableResult):
    temperature_celsius: float
    duplex_fraction: float
    delta_g_kcal_per_mol: float
    association_constant_m_inverse: float


@dataclass(frozen=True)
class TheoreticalMeltingCurveResult(_SerializableResult):
    sequence: str
    temperatures_celsius: tuple[float, ...]
    points: tuple[MeltingCurvePoint, ...]
    midpoint_temperature_celsius: float | None
    conditions: ThermodynamicConditions
    model: str
    applicability: str
    provenance: Provenance


@dataclass(frozen=True)
class TerminalStabilityResult(_SerializableResult):
    sequence: str
    window_size: int
    five_prime_sequence: str
    three_prime_sequence: str
    five_prime_delta_g_kcal_per_mol: float
    three_prime_delta_g_kcal_per_mol: float
    less_stable_end: str
    conditions: ThermodynamicConditions
    five_prime_thermodynamics: NearestNeighborResult
    three_prime_thermodynamics: NearestNeighborResult
    provenance: Provenance


@dataclass(frozen=True)
class CosolventCorrectionResult(_SerializableResult):
    uncorrected_tm_celsius: float
    dmso_percent: float
    dmso_factor_celsius_per_percent: float
    dmso_delta_tm_celsius: float
    formamide_molar: float
    gc_fraction: float
    formamide_delta_tm_celsius: float
    corrected_tm_celsius: float
    method: str
    applicability: str
    provenance: Provenance


__all__ = [
    "A260ConcentrationResult",
    "BindingEquilibriumResult",
    "ConditionalCapability",
    "CosolventCorrectionResult",
    "DuplexStabilityResult",
    "LabelAbsorbanceCorrection",
    "MeltingCurvePoint",
    "MeltingTemperatureResult",
    "MolecularWeightResult",
    "NearestNeighborResult",
    "OligoQuantityResult",
    "OpticalModification",
    "OpticalPropertiesResult",
    "Primer3ThermodynamicResult",
    "SaltCorrectionResult",
    "StackingResult",
    "StackingStep",
    "TerminalStabilityResult",
    "TheoreticalMeltingCurveResult",
    "WindowTmPoint",
    "WindowTmResult",
]
