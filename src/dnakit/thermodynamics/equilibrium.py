"""Duplex equilibrium, theoretical melting-curve and terminal-stability helpers."""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from contextlib import suppress
from itertools import islice, pairwise
from typing import TypeAlias

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError

from ._shared import canonical_linear_symbols, native_provenance
from .calculations import nearest_neighbor
from .config import NearestNeighborConfig, ThermodynamicConditions
from .results import (
    BindingEquilibriumResult,
    CosolventCorrectionResult,
    MeltingCurvePoint,
    TerminalStabilityResult,
    TheoreticalMeltingCurveResult,
)

MeltingCurveProgress: TypeAlias = Callable[[int, int], None]
_GAS_CONSTANT_KCAL_PER_K_MOL = 0.00198720425864083
_MAX_CURVE_POINTS = 100_001


def _conditions(value: ThermodynamicConditions | None) -> ThermodynamicConditions:
    resolved = ThermodynamicConditions() if value is None else value
    if not isinstance(resolved, ThermodynamicConditions):
        raise ConfigurationError(
            "conditions must be ThermodynamicConditions or None.",
            code="INVALID_THERMODYNAMIC_CONDITIONS",
        )
    return resolved


def _association_constant(delta_g_kcal_per_mol: float, temperature_celsius: float) -> float:
    temperature_kelvin = temperature_celsius + 273.15
    exponent = -delta_g_kcal_per_mol / (_GAS_CONSTANT_KCAL_PER_K_MOL * temperature_kelvin)
    return math.exp(min(700.0, max(-700.0, exponent)))


def binding_equilibrium(
    sequence: DNASequence,
    *,
    complement: DNASequence | None = None,
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
) -> BindingEquilibriumResult:
    """Calculate Ka, Kd and exact two-state duplex fraction from mass balance."""

    resolved = _conditions(conditions)
    thermodynamics = nearest_neighbor(
        sequence,
        complement=complement,
        conditions=resolved,
        config=config,
    )
    association = _association_constant(
        thermodynamics.delta_g_kcal_per_mol,
        resolved.temperature_celsius,
    )
    total = resolved.strand_concentration_molar
    if thermodynamics.self_complementary:
        free = 2.0 * total / (1.0 + math.sqrt(1.0 + 8.0 * association * total))
        duplex = max(0.0, (total - free) / 2.0)
        fraction = 0.0 if total == 0.0 else 2.0 * duplex / total
    else:
        each_strand_total = total / 2.0
        free = (
            2.0 * each_strand_total / (1.0 + math.sqrt(1.0 + 4.0 * association * each_strand_total))
        )
        duplex = max(0.0, each_strand_total - free)
        fraction = 0.0 if each_strand_total == 0.0 else duplex / each_strand_total
    return BindingEquilibriumResult(
        temperature_celsius=resolved.temperature_celsius,
        delta_g_kcal_per_mol=thermodynamics.delta_g_kcal_per_mol,
        association_constant_m_inverse=association,
        dissociation_constant_molar=1.0 / association,
        total_strand_concentration_molar=total,
        free_strand_concentration_molar=free,
        duplex_concentration_molar=duplex,
        duplex_fraction=min(1.0, max(0.0, fraction)),
        self_complementary=thermodynamics.self_complementary,
        model="ideal-two-state-perfect-duplex-mass-balance-v1",
        thermodynamics=thermodynamics,
        provenance=thermodynamics.provenance,
    )


def _curve_conditions(
    base: ThermodynamicConditions, temperature_celsius: float
) -> ThermodynamicConditions:
    return ThermodynamicConditions(
        temperature_celsius=temperature_celsius,
        sodium_molar=base.sodium_molar,
        potassium_molar=base.potassium_molar,
        magnesium_molar=base.magnesium_molar,
        dntp_molar=base.dntp_molar,
        strand_concentration_molar=base.strand_concentration_molar,
        dmso_percent=base.dmso_percent,
        dmso_factor_celsius_per_percent=base.dmso_factor_celsius_per_percent,
        formamide_molar=base.formamide_molar,
        salt_model=base.salt_model,
    )


def theoretical_melting_curve(
    sequence: DNASequence,
    temperatures_celsius: Iterable[float],
    *,
    complement: DNASequence | None = None,
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
    progress: MeltingCurveProgress | None = None,
) -> TheoreticalMeltingCurveResult:
    """Evaluate the ideal two-state duplex fraction at bounded temperatures."""

    if isinstance(temperatures_celsius, (str, bytes)) or not isinstance(
        temperatures_celsius, Iterable
    ):
        raise ConfigurationError(
            "temperatures_celsius must be an iterable.", code="INVALID_MELTING_CURVE"
        )
    raw_temperatures = tuple(islice(iter(temperatures_celsius), _MAX_CURVE_POINTS + 1))
    if not 2 <= len(raw_temperatures) <= _MAX_CURVE_POINTS:
        raise ConfigurationError(
            f"A melting curve requires 2-{_MAX_CURVE_POINTS} temperature points.",
            code="INVALID_MELTING_CURVE",
        )
    temperatures: list[float] = []
    for value in raw_temperatures:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0.0 <= value <= 100.0
        ):
            raise ConfigurationError(
                "Every melting-curve temperature must be finite and in [0, 100] C.",
                code="INVALID_MELTING_CURVE_TEMPERATURE",
            )
        temperatures.append(float(value))
    if any(right <= left for left, right in pairwise(temperatures)):
        raise ConfigurationError(
            "Melting-curve temperatures must be strictly increasing.",
            code="UNSORTED_MELTING_CURVE_TEMPERATURES",
        )
    if progress is not None and not callable(progress):
        raise ConfigurationError(
            "progress must be callable or None.", code="INVALID_MELTING_CURVE_PROGRESS"
        )
    resolved = _conditions(conditions)
    points: list[MeltingCurvePoint] = []
    curve_provenance = None
    for index, temperature in enumerate(temperatures, start=1):
        equilibrium = binding_equilibrium(
            sequence,
            complement=complement,
            conditions=_curve_conditions(resolved, temperature),
            config=config,
        )
        points.append(
            MeltingCurvePoint(
                temperature_celsius=temperature,
                duplex_fraction=equilibrium.duplex_fraction,
                delta_g_kcal_per_mol=equilibrium.delta_g_kcal_per_mol,
                association_constant_m_inverse=equilibrium.association_constant_m_inverse,
            )
        )
        curve_provenance = equilibrium.provenance
        if progress is not None:
            with suppress(Exception):
                progress(index, len(temperatures))
    midpoint: float | None = None
    for left, right in pairwise(points):
        if left.duplex_fraction == 0.5:
            midpoint = left.temperature_celsius
            break
        if left.duplex_fraction >= 0.5 >= right.duplex_fraction:
            span = left.duplex_fraction - right.duplex_fraction
            fraction = 0.0 if span == 0.0 else (left.duplex_fraction - 0.5) / span
            midpoint = left.temperature_celsius + fraction * (
                right.temperature_celsius - left.temperature_celsius
            )
            break
    if midpoint is None and points[-1].duplex_fraction == 0.5:
        midpoint = points[-1].temperature_celsius
    sequence_symbols = canonical_linear_symbols(
        sequence,
        operation="theoretical_melting_curve",
        min_length=2,
        max_length=60,
    )
    return TheoreticalMeltingCurveResult(
        sequence=sequence_symbols,
        temperatures_celsius=tuple(temperatures),
        points=tuple(points),
        midpoint_temperature_celsius=midpoint,
        conditions=resolved,
        model="ideal-two-state-perfect-duplex-mass-balance-v1",
        applicability=(
            "Equilibrium two-state curve for the native fully complementary nearest-neighbor "
            "model; it is not an instrument response, kinetic trace, or heat-capacity model."
        ),
        provenance=native_provenance() if curve_provenance is None else curve_provenance,
    )


def terminal_stability(
    sequence: DNASequence,
    *,
    window_size: int = 5,
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
) -> TerminalStabilityResult:
    """Compare nearest-neighbor stability of equal 5-prime and 3-prime windows."""

    symbols = canonical_linear_symbols(
        sequence,
        operation="terminal_stability",
        min_length=2,
        max_length=60,
    )
    if (
        isinstance(window_size, bool)
        or not isinstance(window_size, int)
        or not 2 <= window_size <= len(symbols)
    ):
        raise ConfigurationError(
            "window_size must be an integer in [2, sequence length].",
            code="INVALID_TERMINAL_STABILITY_WINDOW",
        )
    resolved = _conditions(conditions)
    five_sequence = symbols[:window_size]
    three_sequence = symbols[-window_size:]
    five = nearest_neighbor(DNASequence(five_sequence), conditions=resolved, config=config)
    three = nearest_neighbor(DNASequence(three_sequence), conditions=resolved, config=config)
    if math.isclose(five.delta_g_kcal_per_mol, three.delta_g_kcal_per_mol, abs_tol=1e-12):
        less_stable = "equal"
    elif five.delta_g_kcal_per_mol > three.delta_g_kcal_per_mol:
        less_stable = "5-prime"
    else:
        less_stable = "3-prime"
    return TerminalStabilityResult(
        sequence=symbols,
        window_size=window_size,
        five_prime_sequence=five_sequence,
        three_prime_sequence=three_sequence,
        five_prime_delta_g_kcal_per_mol=five.delta_g_kcal_per_mol,
        three_prime_delta_g_kcal_per_mol=three.delta_g_kcal_per_mol,
        less_stable_end=less_stable,
        conditions=resolved,
        five_prime_thermodynamics=five,
        three_prime_thermodynamics=three,
        provenance=five.provenance,
    )


def cosolvent_tm_correction(
    sequence: DNASequence,
    uncorrected_tm_celsius: float,
    *,
    dmso_percent: float = 0.0,
    dmso_factor_celsius_per_percent: float = 0.6,
    formamide_molar: float = 0.0,
) -> CosolventCorrectionResult:
    """Apply explicit Primer3 empirical DMSO and formamide Tm corrections."""

    symbols = canonical_linear_symbols(
        sequence,
        operation="cosolvent_tm_correction",
        min_length=1,
        max_length=1_000_000,
    )
    values = {
        "uncorrected_tm_celsius": uncorrected_tm_celsius,
        "dmso_percent": dmso_percent,
        "dmso_factor_celsius_per_percent": dmso_factor_celsius_per_percent,
        "formamide_molar": formamide_molar,
    }
    resolved_values: dict[str, float] = {}
    for name, value in values.items():
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ConfigurationError(
                f"{name} must be finite numeric data.",
                code="INVALID_COSOLVENT_CORRECTION",
            )
        resolved_values[name] = float(value)
    if not 0.0 <= resolved_values["dmso_percent"] <= 100.0:
        raise ConfigurationError(
            "dmso_percent must be in [0, 100].", code="INVALID_COSOLVENT_CORRECTION"
        )
    if not 0.0 <= resolved_values["dmso_factor_celsius_per_percent"] <= 2.0:
        raise ConfigurationError(
            "dmso_factor_celsius_per_percent must be in [0, 2].",
            code="INVALID_COSOLVENT_CORRECTION",
        )
    if not 0.0 <= resolved_values["formamide_molar"] <= 30.0:
        raise ConfigurationError(
            "formamide_molar must be in [0, 30].", code="INVALID_COSOLVENT_CORRECTION"
        )
    gc_fraction = sum(base in "GC" for base in symbols) / len(symbols)
    dmso_delta = -(
        resolved_values["dmso_factor_celsius_per_percent"] * resolved_values["dmso_percent"]
    )
    formamide_delta = (0.453 * gc_fraction - 2.88) * resolved_values["formamide_molar"]
    corrected = resolved_values["uncorrected_tm_celsius"] + dmso_delta + formamide_delta
    return CosolventCorrectionResult(
        uncorrected_tm_celsius=resolved_values["uncorrected_tm_celsius"],
        dmso_percent=resolved_values["dmso_percent"],
        dmso_factor_celsius_per_percent=resolved_values["dmso_factor_celsius_per_percent"],
        dmso_delta_tm_celsius=dmso_delta,
        formamide_molar=resolved_values["formamide_molar"],
        gc_fraction=gc_fraction,
        formamide_delta_tm_celsius=formamide_delta,
        corrected_tm_celsius=corrected,
        method="primer3-manual-empirical-cosolvent-correction-v1",
        applicability=(
            "Empirical additive correction only; DMSO and formamide effects are not a "
            "mechanistic nearest-neighbor free-energy model."
        ),
        provenance=native_provenance(citation=False),
    )


__all__ = [
    "MeltingCurveProgress",
    "binding_equilibrium",
    "cosolvent_tm_correction",
    "terminal_stability",
    "theoretical_melting_curve",
]
