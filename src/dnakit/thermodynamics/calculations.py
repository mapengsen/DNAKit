"""Strict native physicochemical and thermodynamic calculations for canonical DNA."""

from __future__ import annotations

import math
from typing import Literal

from dnakit.core import (
    Citation,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    Strandedness,
)
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

from ._shared import canonical_linear_symbols, native_provenance
from .backends import Primer3CLIAdapter, validate_primer3_result
from .config import NearestNeighborConfig, ThermodynamicConditions, TmMethod
from .parameters import PARAMETER_SETS, NearestNeighborParameterSet
from .results import (
    DuplexStabilityResult,
    ExtinctionCoefficientResult,
    MeltingTemperatureResult,
    MolecularWeightResult,
    NearestNeighborResult,
    SaltCorrectionResult,
    StackingResult,
    StackingStep,
    WindowTmPoint,
    WindowTmResult,
)

_BASE_MASSES = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.20}
_BASE_MASSES_DALTON = FrozenDict(_BASE_MASSES)
_UNPHOSPHORYLATED_TERMINAL_CORRECTION_DALTON = -61.96
_FIVE_PRIME_PHOSPHATE_ADDITION_DALTON = 79.0
_INDIVIDUAL_BASE_EXTINCTION_260_VALUES = {
    "A": 15_400.0,
    "C": 7_400.0,
    "G": 11_500.0,
    "T": 8_700.0,
}
_INDIVIDUAL_BASE_EXTINCTION_260 = FrozenDict(_INDIVIDUAL_BASE_EXTINCTION_260_VALUES)
_NEAREST_NEIGHBOR_EXTINCTION_260_VALUES = {
    "AA": 27_400.0,
    "AC": 21_200.0,
    "AG": 25_000.0,
    "AT": 22_800.0,
    "CA": 21_200.0,
    "CC": 14_600.0,
    "CG": 18_000.0,
    "CT": 15_200.0,
    "GA": 25_200.0,
    "GC": 17_600.0,
    "GG": 21_600.0,
    "GT": 20_000.0,
    "TA": 23_400.0,
    "TC": 16_200.0,
    "TG": 19_000.0,
    "TT": 16_800.0,
}
_NEAREST_NEIGHBOR_EXTINCTION_260 = FrozenDict(_NEAREST_NEIGHBOR_EXTINCTION_260_VALUES)
_GAS_CONSTANT_CAL_PER_K_MOL = 1.9872
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_NATIVE_APPLICABILITY = (
    "Linear, ungapped, canonical DNA; fully complementary duplex; 2-60 nt; "
    "SantaLucia 1998 unified DNA/DNA parameters; total monovalent Na+ plus K+; "
    "no mismatches, dangling ends, Mg2+, dNTP, cosolvent, or modification model."
)


def _conditions(value: ThermodynamicConditions | None) -> ThermodynamicConditions:
    resolved = ThermodynamicConditions() if value is None else value
    if not isinstance(resolved, ThermodynamicConditions):
        raise ConfigurationError(
            "conditions must be ThermodynamicConditions or None.",
            code="INVALID_THERMODYNAMIC_CONDITIONS",
        )
    return resolved


def _nn_config(value: NearestNeighborConfig | None) -> NearestNeighborConfig:
    resolved = NearestNeighborConfig() if value is None else value
    if not isinstance(resolved, NearestNeighborConfig):
        raise ConfigurationError(
            "config must be NearestNeighborConfig or None.",
            code="INVALID_NEAREST_NEIGHBOR_CONFIG",
        )
    return resolved


def _require_native_ions(conditions: ThermodynamicConditions) -> None:
    if (
        conditions.magnesium_molar != 0.0
        or conditions.dntp_molar != 0.0
        or conditions.dmso_percent != 0.0
        or conditions.formamide_molar != 0.0
    ):
        raise ConfigurationError(
            "The native model does not combine Mg2+, dNTP, DMSO, or formamide with its "
            "monovalent correction.",
            code="UNSUPPORTED_NATIVE_ION_MODEL",
            context={
                "magnesium_molar": conditions.magnesium_molar,
                "dntp_molar": conditions.dntp_molar,
                "dmso_percent": conditions.dmso_percent,
                "formamide_molar": conditions.formamide_molar,
            },
            hint=(
                "Use zero Mg2+/dNTP/cosolvents, an explicit empirical correction, "
                "or a validated external backend."
            ),
        )


def _reverse_complement(symbols: str) -> str:
    return symbols.translate(_COMPLEMENT)[::-1]


def _complement_3to5(symbols: str) -> str:
    return symbols.translate(_COMPLEMENT)


def molecular_weight(
    sequence: DNASequence,
    *,
    strand: Literal["single", "double"] = "single",
    five_prime_phosphorylated: bool = False,
    max_sequence_length: int = 1_000_000,
) -> MolecularWeightResult:
    """Calculate anhydrous unmodified DNA oligonucleotide molecular weight."""

    if strand not in {"single", "double"}:
        raise ConfigurationError(
            "strand must be 'single' or 'double'.",
            code="INVALID_MOLECULAR_WEIGHT_STRAND",
        )
    if not isinstance(five_prime_phosphorylated, bool):
        raise ConfigurationError(
            "five_prime_phosphorylated must be a boolean.",
            code="INVALID_TERMINAL_MODIFICATION",
        )
    if (
        isinstance(max_sequence_length, bool)
        or not isinstance(max_sequence_length, int)
        or not 1 <= max_sequence_length <= 10_000_000
    ):
        raise ConfigurationError(
            "max_sequence_length must be an integer in [1, 10000000].",
            code="INVALID_THERMODYNAMIC_LIMIT",
        )
    symbols = canonical_linear_symbols(
        sequence,
        operation="molecular_weight",
        min_length=1,
        max_length=max_sequence_length,
    )

    def strand_mass(strand_symbols: str) -> float:
        value = math.fsum(_BASE_MASSES[base] for base in strand_symbols)
        value += _UNPHOSPHORYLATED_TERMINAL_CORRECTION_DALTON
        if five_prime_phosphorylated:
            value += _FIVE_PRIME_PHOSPHATE_ADDITION_DALTON
        return value

    value = strand_mass(symbols)
    strand_count = 1 if strand == "single" else 2
    if strand == "double":
        value += strand_mass(_reverse_complement(symbols))
    return MolecularWeightResult(
        value_dalton=value,
        value_kilodalton=value / 1000.0,
        strand_count=strand_count,
        sequence_length=len(symbols),
        five_prime_phosphorylated=five_prime_phosphorylated,
        method="anhydrous-deoxynucleotide-residue-sum",
        algorithm_version="dnakit-molecular-weight-v1",
        applicability=(
            "Linear, ungapped, canonical, unmodified DNA with hydroxyl termini; "
            "optional 5-prime monophosphate. Double strand adds the complete "
            "Watson-Crick reverse-complement strand."
        ),
        parameters=FrozenDict(
            {
                "base_masses_dalton": _BASE_MASSES_DALTON,
                "unphosphorylated_terminal_correction_dalton": (
                    _UNPHOSPHORYLATED_TERMINAL_CORRECTION_DALTON
                ),
                "five_prime_phosphate_addition_dalton": (_FIVE_PRIME_PHOSPHATE_ADDITION_DALTON),
                "mass_table_version": "anhydrous-dna-oligo-residues-v1",
                "mass_table_units": "Da per residue",
                "mass_table_source": (
                    "Standard anhydrous DNA oligonucleotide residue formula: "
                    "sum(dA,dC,dG,dT) - 61.96 Da"
                ),
            }
        ),
        provenance=native_provenance(citation=False),
    )


def extinction_coefficient_260nm(
    sequence: DNASequence,
    *,
    max_sequence_length: int = 1_000_000,
) -> ExtinctionCoefficientResult:
    """Calculate the theoretical molar extinction coefficient of an ssDNA oligo at 260 nm."""

    if (
        isinstance(max_sequence_length, bool)
        or not isinstance(max_sequence_length, int)
        or not 1 <= max_sequence_length <= 10_000_000
    ):
        raise ConfigurationError(
            "max_sequence_length must be an integer in [1, 10000000].",
            code="INVALID_THERMODYNAMIC_LIMIT",
        )
    symbols = canonical_linear_symbols(
        sequence,
        operation="extinction_coefficient_260nm",
        min_length=1,
        max_length=max_sequence_length,
    )
    if sequence.strandedness is not Strandedness.SINGLE:
        raise ConfigurationError(
            "extinction_coefficient_260nm requires a single-stranded DNA sequence.",
            code="EXTINCTION_COEFFICIENT_SINGLE_STRAND_ONLY",
            hint=(
                "Calculate each unhybridized oligonucleotide strand separately; "
                "duplex absorbance requires a duplex-specific hypochromicity model."
            ),
        )

    if len(symbols) == 1:
        value = _INDIVIDUAL_BASE_EXTINCTION_260_VALUES[symbols]
    else:
        dinucleotide_sum = math.fsum(
            _NEAREST_NEIGHBOR_EXTINCTION_260_VALUES[symbols[index : index + 2]]
            for index in range(len(symbols) - 1)
        )
        internal_base_sum = math.fsum(
            _INDIVIDUAL_BASE_EXTINCTION_260_VALUES[base] for base in symbols[1:-1]
        )
        value = dinucleotide_sum - internal_base_sum

    provenance = Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.REIMPLEMENTATION,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.PUBLISHED_ALGORITHM,
            citations=(
                Citation(
                    "warshaw-tinoco1966",
                    title="Optical properties of sixteen dinucleoside phosphates",
                    doi="10.1016/0022-2836(66)90115-X",
                ),
                Citation(
                    "cantor-warshaw-shapiro1970",
                    title=(
                        "Oligonucleotide interactions. III. Circular dichroism studies "
                        "of the conformation of deoxyoligonucleotides"
                    ),
                    doi="10.1002/bip.1970.360090909",
                ),
            ),
        )
    )
    return ExtinctionCoefficientResult(
        value_m_inverse_cm_inverse=value,
        wavelength_nm=260,
        sequence_length=len(symbols),
        method="nearest-neighbor-hypochromicity",
        algorithm_version="dnakit-extinction-coefficient-260nm-v1",
        applicability=(
            "Linear, ungapped, canonical, single-stranded, unmodified DNA at 260 nm. "
            "The result is theoretical and does not replace an experimental A260 measurement; "
            "fluorophores and other modifications require separate correction."
        ),
        parameters=FrozenDict(
            {
                "parameter_set": "warshaw-tinoco-cantor-dna-260nm-v1",
                "formula": "sum(pair_i, i=1..N-1) - sum(base_i, i=2..N-1)",
                "reference_temperature_celsius": 25.0,
                "reference_ph": 7.0,
                "individual_base_m_inverse_cm_inverse": (_INDIVIDUAL_BASE_EXTINCTION_260),
                "nearest_neighbor_m_inverse_cm_inverse": (_NEAREST_NEIGHBOR_EXTINCTION_260),
                "single_nucleotide_rule": "use individual-base coefficient",
            }
        ),
        provenance=provenance,
    )


def salt_correction(
    sequence_length: int,
    *,
    conditions: ThermodynamicConditions | None = None,
) -> SaltCorrectionResult:
    """Return the SantaLucia 1998 monovalent entropy correction."""

    resolved = _conditions(conditions)
    _require_native_ions(resolved)
    if isinstance(sequence_length, bool) or not isinstance(sequence_length, int):
        raise ConfigurationError(
            "sequence_length must be an integer.",
            code="INVALID_THERMODYNAMIC_SEQUENCE_LENGTH",
        )
    if not 2 <= sequence_length <= 60:
        raise ConfigurationError(
            "Native salt correction is limited to sequence lengths from 2 to 60 nt.",
            code="THERMODYNAMIC_SEQUENCE_LENGTH_OUT_OF_DOMAIN",
            context={"sequence_length": sequence_length},
        )
    correction = 0.368 * (sequence_length - 1) * math.log(resolved.monovalent_molar)
    return SaltCorrectionResult(
        delta_s_cal_per_k_mol=correction,
        sequence_length=sequence_length,
        model=resolved.salt_model,
        model_version="santalucia1998-equation-v1",
        conditions=resolved,
        applicability=(
            "Total monovalent Na+ plus K+ correction; length 2-60 nt; "
            "Mg2+/dNTP/cosolvents excluded."
        ),
        provenance=native_provenance(),
    )


def _stacking_steps(
    symbols: str,
    *,
    temperature_celsius: float,
    parameter_set: NearestNeighborParameterSet,
) -> tuple[StackingStep, ...]:
    temperature_kelvin = temperature_celsius + 273.15
    steps: list[StackingStep] = []
    complement = _complement_3to5(symbols)
    for index in range(len(symbols) - 1):
        top = symbols[index : index + 2]
        parameter = parameter_set.stacking[top]
        delta_g = parameter.delta_h_kcal_per_mol - (
            temperature_kelvin * parameter.delta_s_cal_per_k_mol / 1000.0
        )
        steps.append(
            StackingStep(
                index=index,
                top_5to3=top,
                bottom_3to5=complement[index : index + 2],
                delta_h_kcal_per_mol=parameter.delta_h_kcal_per_mol,
                delta_s_cal_per_k_mol=parameter.delta_s_cal_per_k_mol,
                delta_g_kcal_per_mol=delta_g,
            )
        )
    return tuple(steps)


def stacking_interactions(
    sequence: DNASequence,
    *,
    temperature_celsius: float = 37.0,
    config: NearestNeighborConfig | None = None,
) -> StackingResult:
    """Return each canonical nearest-neighbor stacking contribution."""

    resolved_config = _nn_config(config)
    temperature = ThermodynamicConditions(
        temperature_celsius=temperature_celsius,
        sodium_molar=1.0,
    ).temperature_celsius
    symbols = canonical_linear_symbols(
        sequence,
        operation="stacking_interactions",
        min_length=2,
        max_length=resolved_config.max_sequence_length,
    )
    parameter_set = PARAMETER_SETS[resolved_config.parameter_set]
    steps = _stacking_steps(
        symbols,
        temperature_celsius=temperature,
        parameter_set=parameter_set,
    )
    delta_h = math.fsum(step.delta_h_kcal_per_mol for step in steps)
    delta_s = math.fsum(step.delta_s_cal_per_k_mol for step in steps)
    return StackingResult(
        sequence=symbols,
        complement_3to5=_complement_3to5(symbols),
        temperature_celsius=temperature,
        parameter_set=parameter_set.version,
        steps=steps,
        total_delta_h_kcal_per_mol=delta_h,
        total_delta_s_cal_per_k_mol=delta_s,
        total_delta_g_kcal_per_mol=(delta_h - (temperature + 273.15) * delta_s / 1000.0),
        applicability=(
            "Canonical DNA/DNA stacking steps only; excludes initiation, symmetry, and salt."
        ),
        provenance=native_provenance(),
    )


def _nearest_neighbor_symbols(
    symbols: str,
    *,
    conditions: ThermodynamicConditions,
    parameter_set: NearestNeighborParameterSet,
) -> NearestNeighborResult:
    steps = _stacking_steps(
        symbols,
        temperature_celsius=conditions.temperature_celsius,
        parameter_set=parameter_set,
    )
    stacking_h = math.fsum(step.delta_h_kcal_per_mol for step in steps)
    stacking_s = math.fsum(step.delta_s_cal_per_k_mol for step in steps)
    terminal_parameters = (
        parameter_set.terminal_at if base in "AT" else parameter_set.terminal_gc
        for base in (symbols[0], symbols[-1])
    )
    terminal_tuple = tuple(terminal_parameters)
    initiation_h = math.fsum(item.delta_h_kcal_per_mol for item in terminal_tuple)
    initiation_s = math.fsum(item.delta_s_cal_per_k_mol for item in terminal_tuple)
    self_complementary = symbols == _reverse_complement(symbols)
    symmetry_h = parameter_set.symmetry.delta_h_kcal_per_mol if self_complementary else 0.0
    symmetry_s = parameter_set.symmetry.delta_s_cal_per_k_mol if self_complementary else 0.0
    salt_s = 0.368 * (len(symbols) - 1) * math.log(conditions.monovalent_molar)
    delta_h = math.fsum((stacking_h, initiation_h, symmetry_h))
    delta_s = math.fsum((stacking_s, initiation_s, symmetry_s, salt_s))
    temperature_kelvin = conditions.temperature_celsius + 273.15
    delta_g = delta_h - temperature_kelvin * delta_s / 1000.0
    divisor = 1 if self_complementary else 4
    concentration_term = _GAS_CONSTANT_CAL_PER_K_MOL * math.log(
        conditions.strand_concentration_molar / divisor
    )
    denominator = delta_s + concentration_term
    if denominator == 0.0:
        raise ConfigurationError(
            "The selected conditions produce a singular Tm denominator.",
            code="SINGULAR_TM_CONDITIONS",
        )
    tm_celsius = (1000.0 * delta_h / denominator) - 273.15
    return NearestNeighborResult(
        sequence=symbols,
        complement_5to3=_reverse_complement(symbols),
        sequence_length=len(symbols),
        parameter_set=parameter_set.version,
        conditions=conditions,
        stacking_steps=steps,
        stacking_delta_h_kcal_per_mol=stacking_h,
        stacking_delta_s_cal_per_k_mol=stacking_s,
        initiation_delta_h_kcal_per_mol=initiation_h,
        initiation_delta_s_cal_per_k_mol=initiation_s,
        symmetry_delta_h_kcal_per_mol=symmetry_h,
        symmetry_delta_s_cal_per_k_mol=symmetry_s,
        salt_delta_s_cal_per_k_mol=salt_s,
        delta_h_kcal_per_mol=delta_h,
        delta_s_cal_per_k_mol=delta_s,
        delta_g_kcal_per_mol=delta_g,
        tm_celsius=tm_celsius,
        self_complementary=self_complementary,
        concentration_divisor=divisor,
        gas_constant_cal_per_k_mol=_GAS_CONSTANT_CAL_PER_K_MOL,
        reference_sodium_molar=parameter_set.reference_sodium_molar,
        tm_equation="1000*dH/(dS+R*ln(Ct/divisor))-273.15",
        method="nearest-neighbor-perfect-duplex",
        algorithm_version="dnakit-santalucia1998-v1",
        applicability=_NATIVE_APPLICABILITY,
        provenance=native_provenance(),
    )


def nearest_neighbor(
    sequence: DNASequence,
    *,
    complement: DNASequence | None = None,
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
) -> NearestNeighborResult:
    """Calculate complete Watson-Crick duplex ΔH, ΔS, ΔG, and Tm."""

    resolved_conditions = _conditions(conditions)
    _require_native_ions(resolved_conditions)
    resolved_config = _nn_config(config)
    symbols = canonical_linear_symbols(
        sequence,
        operation="nearest_neighbor",
        min_length=2,
        max_length=resolved_config.max_sequence_length,
    )
    if complement is not None:
        complement_symbols = canonical_linear_symbols(
            complement,
            operation="nearest_neighbor complement",
            min_length=2,
            max_length=resolved_config.max_sequence_length,
        )
        expected = _reverse_complement(symbols)
        if complement_symbols != expected:
            raise ConfigurationError(
                "The native nearest-neighbor model supports only a full reverse complement.",
                code="MISMATCHED_DUPLEX_UNSUPPORTED",
                context={"expected_complement_5to3": expected},
                hint="Use a validated mismatch-capable external backend.",
            )
    return _nearest_neighbor_symbols(
        symbols,
        conditions=resolved_conditions,
        parameter_set=PARAMETER_SETS[resolved_config.parameter_set],
    )


def melting_temperature(
    sequence: DNASequence,
    *,
    method: TmMethod = "nearest_neighbor",
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
) -> MeltingTemperatureResult:
    """Calculate Tm using either the bounded Wallace rule or native NN model."""

    resolved_conditions = _conditions(conditions)
    if method == "wallace":
        if config is not None:
            raise ConfigurationError(
                "NearestNeighborConfig is not applicable to the Wallace rule.",
                code="TM_CONFIG_NOT_APPLICABLE",
            )
        symbols = canonical_linear_symbols(
            sequence,
            operation="Wallace melting_temperature",
            min_length=2,
            max_length=13,
        )
        tm = 2.0 * sum(symbol in "AT" for symbol in symbols) + 4.0 * sum(
            symbol in "GC" for symbol in symbols
        )
        return MeltingTemperatureResult(
            tm_celsius=tm,
            sequence_length=len(symbols),
            method="wallace-2at-4gc",
            algorithm_version="wallace-short-oligo-v1",
            conditions=resolved_conditions,
            parameter_set=None,
            applicability=(
                "Empirical 2(A+T)+4(G+C) rule restricted here to canonical 2-13 nt DNA; "
                "recorded ionic and concentration conditions are not modeled."
            ),
            provenance=native_provenance(citation=False),
        )
    if method != "nearest_neighbor":
        raise ConfigurationError(
            "method must be 'wallace' or 'nearest_neighbor'.",
            code="UNKNOWN_TM_METHOD",
            context={"method": method},
        )
    result = nearest_neighbor(
        sequence,
        conditions=resolved_conditions,
        config=config,
    )
    return MeltingTemperatureResult(
        tm_celsius=result.tm_celsius,
        sequence_length=result.sequence_length,
        method=result.method,
        algorithm_version=result.algorithm_version,
        conditions=result.conditions,
        parameter_set=result.parameter_set,
        applicability=result.applicability,
        provenance=result.provenance,
    )


def duplex_stability(
    sequence_a: DNASequence,
    sequence_b: DNASequence,
    *,
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
    backend: Literal["native", "primer3-cli"] = "native",
    adapter: Primer3CLIAdapter | None = None,
    max_loop: int = 30,
    output_structure: bool = False,
) -> DuplexStabilityResult:
    """Evaluate duplex stability with a native or explicit Primer3 backend.

    The default native path is restricted to a complete Watson-Crick reverse
    complement. ``backend="primer3-cli"`` is an explicit opt-in path for
    mismatches and dangling ends; it requires a caller-configured CLI adapter,
    never searches, installs, or downloads Primer3, and records the resolved
    executable in result provenance.
    """

    if not isinstance(backend, str) or backend not in {"native", "primer3-cli"}:
        raise ConfigurationError(
            "backend must be 'native' or 'primer3-cli'.",
            code="UNKNOWN_DUPLEX_STABILITY_BACKEND",
            context={"backend": backend},
        )
    if adapter is not None and not isinstance(adapter, Primer3CLIAdapter):
        raise ConfigurationError(
            "adapter must be Primer3CLIAdapter or None.",
            code="INVALID_DUPLEX_STABILITY_ADAPTER",
        )
    if isinstance(max_loop, bool) or not isinstance(max_loop, int) or not 1 <= max_loop <= 30:
        raise ConfigurationError(
            "max_loop must be an integer in [1, 30].",
            code="INVALID_DUPLEX_STABILITY_MAX_LOOP",
        )
    if not isinstance(output_structure, bool):
        raise ConfigurationError(
            "output_structure must be boolean.",
            code="INVALID_DUPLEX_STABILITY_STRUCTURE_OPTION",
        )
    if backend == "native" and adapter is not None:
        raise ConfigurationError(
            "adapter is only valid when backend='primer3-cli'.",
            code="UNUSED_DUPLEX_STABILITY_ADAPTER",
        )
    if backend == "native" and (max_loop != 30 or output_structure):
        raise ConfigurationError(
            "max_loop and output_structure are only valid with backend='primer3-cli'.",
            code="UNUSED_DUPLEX_STABILITY_STRUCTURE_OPTION",
        )
    if backend == "primer3-cli" and config is not None:
        raise ConfigurationError(
            "NearestNeighborConfig is only valid with backend='native'.",
            code="UNUSED_DUPLEX_STABILITY_NATIVE_CONFIG",
        )
    if backend == "primer3-cli":
        resolved_adapter = Primer3CLIAdapter() if adapter is None else adapter
        resolved_conditions = _conditions(conditions)
        requested_a = canonical_linear_symbols(
            sequence_a,
            operation="Primer3 duplex stability sequence A",
            min_length=1,
            max_length=60,
        )
        requested_b = canonical_linear_symbols(
            sequence_b,
            operation="Primer3 duplex stability sequence B",
            min_length=1,
            max_length=60,
        )
        raw_primer3_result = resolved_adapter.heterodimer(
            sequence_a,
            sequence_b,
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        primer3_result = validate_primer3_result(
            raw_primer3_result,
            capability="heterodimer",
            sequences_5to3=(requested_a, requested_b),
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
            error_code="MISMATCHED_PRIMER3_DUPLEX_RESULT",
        )
        delta_g = primer3_result.delta_g_kcal_per_mol
        if delta_g is None:  # Defensive guard for third-party adapter contract drift.
            raise ConfigurationError(
                "Primer3 heterodimer result omitted delta_g_kcal_per_mol.",
                code="INVALID_PRIMER3_DUPLEX_RESULT",
            )
        symbols_a, symbols_b = primer3_result.sequences_5to3
        fully_complementary = symbols_b == _reverse_complement(symbols_a)
        return DuplexStabilityResult(
            sequence_a_5to3=symbols_a,
            sequence_b_5to3=symbols_b,
            fully_complementary=fully_complementary,
            stable_at_temperature=(
                bool(primer3_result.structure_found)
                and primer3_result.tm_celsius > primer3_result.conditions.temperature_celsius
            ),
            stability_criterion=(
                "primer3_structure_found and "
                "heterodimer_tm_celsius > configured_temperature_celsius"
            ),
            delta_g_kcal_per_mol=delta_g,
            tm_celsius=primer3_result.tm_celsius,
            conditions=primer3_result.conditions,
            model=primer3_result.algorithm_version,
            applicability=(
                "Canonical linear DNA of 1-60 nt per strand; Primer3 heterodimer "
                "model permits mismatch and dangling-end structures selected by the "
                "backend. User-supplied alignments and chemical modifications are excluded."
            ),
            thermodynamics=primer3_result,
            provenance=primer3_result.provenance,
        )

    result = nearest_neighbor(
        sequence_a,
        complement=sequence_b,
        conditions=conditions,
        config=config,
    )
    return DuplexStabilityResult(
        sequence_a_5to3=result.sequence,
        sequence_b_5to3=result.complement_5to3,
        fully_complementary=True,
        stable_at_temperature=result.tm_celsius > result.conditions.temperature_celsius,
        stability_criterion="tm_celsius > configured_temperature_celsius",
        delta_g_kcal_per_mol=result.delta_g_kcal_per_mol,
        tm_celsius=result.tm_celsius,
        conditions=result.conditions,
        model=result.algorithm_version,
        applicability=_NATIVE_APPLICABILITY,
        thermodynamics=result,
        provenance=result.provenance,
    )


def window_tm(
    sequence: DNASequence,
    window_size: int,
    *,
    step: int = 1,
    method: TmMethod = "nearest_neighbor",
    conditions: ThermodynamicConditions | None = None,
    config: NearestNeighborConfig | None = None,
    max_windows: int = 100_000,
) -> WindowTmResult:
    """Calculate a bounded local Tm profile in 0-based half-open coordinates."""

    for name, value in (
        ("window_size", window_size),
        ("step", step),
        ("max_windows", max_windows),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ConfigurationError(
                f"{name} must be a positive integer.",
                code="INVALID_WINDOW_TM_CONFIG",
                context={"field": name, "value": value},
            )
    if method not in {"wallace", "nearest_neighbor"}:
        raise ConfigurationError(
            "method must be 'wallace' or 'nearest_neighbor'.",
            code="UNKNOWN_TM_METHOD",
        )
    if method == "wallace" and config is not None:
        raise ConfigurationError(
            "NearestNeighborConfig is not applicable to a Wallace-rule Tm profile.",
            code="TM_CONFIG_NOT_APPLICABLE",
        )
    if max_windows > 1_000_000:
        raise ConfigurationError(
            "max_windows cannot exceed the hard safety ceiling of 1000000.",
            code="INVALID_WINDOW_TM_CONFIG",
        )
    resolved_config: NearestNeighborConfig | None = None
    if method == "wallace":
        method_max = 13
    else:
        resolved_config = _nn_config(config)
        method_max = resolved_config.max_sequence_length
    if not 2 <= window_size <= method_max:
        raise ConfigurationError(
            f"window_size must be in [2, {method_max}] for method={method!r}.",
            code="WINDOW_TM_SIZE_OUT_OF_DOMAIN",
        )
    symbols = canonical_linear_symbols(
        sequence,
        operation="window_tm",
        min_length=window_size,
        max_length=10_000_000,
    )
    count = ((len(symbols) - window_size) // step) + 1
    if count > max_windows:
        raise ConfigurationError(
            "Requested local Tm profile exceeds max_windows.",
            code="WINDOW_TM_LIMIT_EXCEEDED",
            context={"window_count": count, "max_windows": max_windows},
            hint="Increase step, decrease input length, or explicitly raise max_windows.",
        )
    resolved_conditions = _conditions(conditions)
    if method == "nearest_neighbor":
        _require_native_ions(resolved_conditions)
        assert resolved_config is not None
        parameter_set = PARAMETER_SETS[resolved_config.parameter_set]
    else:
        parameter_set = PARAMETER_SETS["santalucia1998-v1"]
    points: list[WindowTmPoint] = []
    for start in range(0, len(symbols) - window_size + 1, step):
        subsequence = symbols[start : start + window_size]
        if method == "wallace":
            tm = float(
                2 * sum(symbol in "AT" for symbol in subsequence)
                + 4 * sum(symbol in "GC" for symbol in subsequence)
            )
        else:
            tm = _nearest_neighbor_symbols(
                subsequence,
                conditions=resolved_conditions,
                parameter_set=parameter_set,
            ).tm_celsius
        points.append(WindowTmPoint(start, start + window_size, subsequence, tm))
    tm_values = tuple(point.tm_celsius for point in points)
    return WindowTmResult(
        sequence_length=len(symbols),
        window_size=window_size,
        step=step,
        method=method,
        conditions=resolved_conditions,
        windows=tuple(points),
        min_tm_celsius=min(tm_values, default=None),
        max_tm_celsius=max(tm_values, default=None),
        max_windows=max_windows,
        coordinate_system="0-based-half-open-symbol",
        applicability=(
            "Local fixed-window computational Tm profile; not an experimental melting curve. "
            + (
                "Wallace rule conditions are recorded but not modeled."
                if method == "wallace"
                else _NATIVE_APPLICABILITY
            )
        ),
        provenance=(
            native_provenance(citation=False) if method == "wallace" else native_provenance()
        ),
    )


__all__ = [
    "duplex_stability",
    "extinction_coefficient_260nm",
    "melting_temperature",
    "molecular_weight",
    "nearest_neighbor",
    "salt_correction",
    "stacking_interactions",
    "window_tm",
]
