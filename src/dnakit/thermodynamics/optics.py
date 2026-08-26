"""Optical properties and unit-safe oligonucleotide concentration conversions."""

from __future__ import annotations

import math
from collections.abc import Iterable
from itertools import islice
from typing import Literal

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError

from ._shared import canonical_linear_symbols
from .calculations import extinction_coefficient_260nm, molecular_weight
from .results import (
    A260ConcentrationResult,
    LabelAbsorbanceCorrection,
    OligoQuantityResult,
    OpticalModification,
    OpticalPropertiesResult,
)

_AVERAGE_DUPLEX_EXTINCTION_PER_BASE_PAIR = 13_200.0
_MAX_CORRECTIONS = 1_000


def _finite(
    value: object,
    name: str,
    *,
    minimum: float | None = None,
    strictly_positive: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"{name} must be finite numeric data.",
            code="INVALID_OPTICAL_VALUE",
            context={"field": name, "value": value},
        )
    result = float(value)
    if strictly_positive and result <= 0.0:
        raise ConfigurationError(
            f"{name} must be greater than zero.",
            code="INVALID_OPTICAL_VALUE",
            context={"field": name, "value": value},
        )
    if minimum is not None and result < minimum:
        raise ConfigurationError(
            f"{name} must be at least {minimum}.",
            code="INVALID_OPTICAL_VALUE",
            context={"field": name, "value": value},
        )
    return result


def _bounded_tuple(
    values: Iterable[object], expected_type: type[object], name: str
) -> tuple[object, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Iterable):
        raise ConfigurationError(f"{name} must be an iterable.", code="INVALID_OPTICAL_VALUE")
    resolved = tuple(islice(iter(values), _MAX_CORRECTIONS + 1))
    if len(resolved) > _MAX_CORRECTIONS:
        raise ConfigurationError(
            f"{name} exceeds the {_MAX_CORRECTIONS}-item safety limit.",
            code="OPTICAL_CORRECTION_LIMIT_EXCEEDED",
        )
    if any(not isinstance(item, expected_type) for item in resolved):
        raise ConfigurationError(
            f"Every {name} item must be {expected_type.__name__}.",
            code="INVALID_OPTICAL_VALUE",
        )
    return resolved


def optical_properties(
    sequence: DNASequence,
    *,
    strand_type: Literal["single", "double"] = "single",
    complement: DNASequence | None = None,
    duplex_method: Literal["average-base-pair", "strand-sum-hypochromicity"] = (
        "average-base-pair"
    ),
    hypochromicity_fraction: float | None = None,
    modifications: Iterable[OpticalModification] = (),
) -> OpticalPropertiesResult:
    """Calculate epsilon260, molecular weight, nmol/OD and mass/OD.

    ``average-base-pair`` uses 13,200 M^-1 cm^-1 per base pair and is the
    conventional 1 OD260 = 50 microgram/mL duplex approximation.
    ``strand-sum-hypochromicity`` requires the caller to provide the fractional
    absorbance reduction on duplex formation instead of silently assuming one.
    """

    if strand_type not in {"single", "double"}:
        raise ConfigurationError(
            "strand_type must be 'single' or 'double'.",
            code="INVALID_OPTICAL_STRAND_TYPE",
        )
    symbols = canonical_linear_symbols(
        sequence,
        operation="optical_properties",
        min_length=1,
        max_length=1_000_000,
    )
    sequence_single = DNASequence(symbols)
    resolved_modifications = _bounded_tuple(modifications, OpticalModification, "modifications")
    typed_modifications = tuple(
        item for item in resolved_modifications if isinstance(item, OpticalModification)
    )

    sequences: tuple[str, ...]
    assumptions: tuple[str, ...]
    if strand_type == "single":
        if complement is not None:
            raise ConfigurationError(
                "complement is only valid for double-stranded optical properties.",
                code="UNUSED_OPTICAL_COMPLEMENT",
            )
        if hypochromicity_fraction is not None:
            raise ConfigurationError(
                "hypochromicity_fraction is only valid for a duplex strand-sum model.",
                code="UNUSED_HYPOCHROMICITY",
            )
        extinction = extinction_coefficient_260nm(sequence_single)
        native_epsilon = extinction.value_m_inverse_cm_inverse
        mass = molecular_weight(sequence_single).value_dalton
        sequences = (symbols,)
        method = extinction.method
        assumptions = (
            "Single-stranded, linear, canonical and unmodified DNA before explicit corrections.",
            "One OD260 unit means A260 multiplied by sample volume in mL equals one, "
            "using a 1 cm reference path.",
        )
        provenance = extinction.provenance
    else:
        expected_complement = sequence_single.reverse_complement().symbols
        if complement is None:
            complement_symbols = expected_complement
        else:
            complement_symbols = canonical_linear_symbols(
                complement,
                operation="optical_properties complement",
                min_length=1,
                max_length=1_000_000,
            )
            if complement_symbols != expected_complement:
                raise ConfigurationError(
                    "The native duplex optical model requires a full reverse complement.",
                    code="MISMATCHED_OPTICAL_DUPLEX",
                    context={"expected_complement_5to3": expected_complement},
                )
        if duplex_method == "average-base-pair":
            if hypochromicity_fraction is not None:
                raise ConfigurationError(
                    "average-base-pair does not accept hypochromicity_fraction.",
                    code="UNUSED_HYPOCHROMICITY",
                )
            native_epsilon = _AVERAGE_DUPLEX_EXTINCTION_PER_BASE_PAIR * len(symbols)
            method = "average-dsdna-base-pair-extinction"
            assumptions = (
                "Fully complementary duplex with average epsilon260 of 13,200 M^-1 cm^-1 "
                "per base pair.",
                "This average conversion is not a sequence-specific hypochromicity prediction.",
                "One OD260 unit means A260 multiplied by sample volume in mL equals one, "
                "using a 1 cm reference path.",
            )
        elif duplex_method == "strand-sum-hypochromicity":
            if hypochromicity_fraction is None:
                raise ConfigurationError(
                    "strand-sum-hypochromicity requires hypochromicity_fraction.",
                    code="MISSING_HYPOCHROMICITY",
                )
            fraction = _finite(
                hypochromicity_fraction,
                "hypochromicity_fraction",
                minimum=0.0,
            )
            if fraction >= 1.0:
                raise ConfigurationError(
                    "hypochromicity_fraction must be in [0, 1).",
                    code="INVALID_HYPOCHROMICITY",
                )
            primary = extinction_coefficient_260nm(sequence_single)
            secondary = extinction_coefficient_260nm(DNASequence(complement_symbols))
            native_epsilon = (
                primary.value_m_inverse_cm_inverse + secondary.value_m_inverse_cm_inverse
            ) * (1.0 - fraction)
            method = "sequence-specific-strand-sum-with-explicit-hypochromicity"
            assumptions = (
                "Fully complementary duplex.",
                "The hypochromicity fraction is caller-supplied experimental or literature input.",
                "One OD260 unit means A260 multiplied by sample volume in mL equals one, "
                "using a 1 cm reference path.",
            )
        else:
            raise ConfigurationError(
                "Unknown duplex_method.",
                code="UNKNOWN_DUPLEX_OPTICAL_METHOD",
                context={"duplex_method": duplex_method},
            )
        mass_result = molecular_weight(sequence_single, strand="double")
        mass = mass_result.value_dalton
        provenance = mass_result.provenance
        sequences = (symbols, complement_symbols)

    modification_epsilon = math.fsum(
        item.count * item.extinction_coefficient_260_delta_m_inverse_cm_inverse
        for item in typed_modifications
    )
    modification_mass = math.fsum(
        item.count * item.molecular_weight_delta_dalton for item in typed_modifications
    )
    corrected_epsilon = native_epsilon + modification_epsilon
    corrected_mass = mass + modification_mass
    if corrected_epsilon <= 0.0 or corrected_mass <= 0.0:
        raise ConfigurationError(
            "Modification corrections must leave positive epsilon260 and molecular weight.",
            code="INVALID_OPTICAL_CORRECTION_TOTAL",
        )
    return OpticalPropertiesResult(
        strand_type=strand_type,
        sequences_5to3=sequences,
        sequence_lengths=tuple(len(item) for item in sequences),
        native_extinction_coefficient_260_m_inverse_cm_inverse=native_epsilon,
        modification_extinction_coefficient_260_m_inverse_cm_inverse=modification_epsilon,
        extinction_coefficient_260_m_inverse_cm_inverse=corrected_epsilon,
        molecular_weight_dalton=corrected_mass,
        one_od260_nmol=1_000_000.0 / corrected_epsilon,
        one_od260_microgram=1_000.0 * corrected_mass / corrected_epsilon,
        method=method,
        modifications=typed_modifications,
        assumptions=assumptions,
        provenance=provenance,
    )


def concentration_from_a260(
    measured_a260: float,
    properties: OpticalPropertiesResult,
    *,
    path_length_cm: float = 1.0,
    dilution_factor: float = 1.0,
    label_corrections: Iterable[LabelAbsorbanceCorrection] = (),
    volume_liter: float | None = None,
) -> A260ConcentrationResult:
    """Apply dye correction and Beer-Lambert law to an A260 measurement."""

    absorbance = _finite(measured_a260, "measured_a260", minimum=0.0)
    path = _finite(path_length_cm, "path_length_cm", strictly_positive=True)
    dilution = _finite(dilution_factor, "dilution_factor", strictly_positive=True)
    if not isinstance(properties, OpticalPropertiesResult):
        raise ConfigurationError(
            "properties must be OpticalPropertiesResult.", code="INVALID_OPTICAL_PROPERTIES"
        )
    resolved_corrections = _bounded_tuple(
        label_corrections, LabelAbsorbanceCorrection, "label_corrections"
    )
    typed_corrections = tuple(
        item for item in resolved_corrections if isinstance(item, LabelAbsorbanceCorrection)
    )
    subtracted = math.fsum(
        item.absorbance_at_label_max * item.a260_correction_factor for item in typed_corrections
    )
    corrected = absorbance - subtracted
    if corrected < 0.0:
        raise ConfigurationError(
            "Label correction exceeds the measured A260.",
            code="NEGATIVE_CORRECTED_A260",
            context={"measured_a260": absorbance, "subtracted_a260": subtracted},
        )
    epsilon = properties.extinction_coefficient_260_m_inverse_cm_inverse
    molar = corrected * dilution / (epsilon * path)
    mass_g_per_l = molar * properties.molecular_weight_dalton
    resolved_volume = (
        None
        if volume_liter is None
        else _finite(volume_liter, "volume_liter", strictly_positive=True)
    )
    amount = None if resolved_volume is None else molar * resolved_volume
    mass_microgram = (
        None if resolved_volume is None else mass_g_per_l * resolved_volume * 1_000_000.0
    )
    return A260ConcentrationResult(
        measured_a260=absorbance,
        label_a260_subtracted=subtracted,
        corrected_a260=corrected,
        path_length_cm=path,
        dilution_factor=dilution,
        extinction_coefficient_260_m_inverse_cm_inverse=epsilon,
        molar_concentration_molar=molar,
        molar_concentration_micromolar=molar * 1_000_000.0,
        molar_concentration_nanomolar=molar * 1_000_000_000.0,
        mass_concentration_g_per_l=mass_g_per_l,
        mass_concentration_ng_per_microliter=mass_g_per_l * 1_000.0,
        volume_liter=resolved_volume,
        amount_mol=amount,
        mass_microgram=mass_microgram,
        label_corrections=typed_corrections,
        method="beer-lambert-a260-with-explicit-label-correction",
        provenance=properties.provenance,
    )


def convert_oligo_quantity(
    molecular_weight_dalton: float,
    *,
    volume_liter: float | None = None,
    molar_concentration_molar: float | None = None,
    mass_concentration_g_per_l: float | None = None,
    amount_mol: float | None = None,
    mass_g: float | None = None,
) -> OligoQuantityResult:
    """Convert exactly one concentration or amount input into all four forms."""

    molecular_weight = _finite(
        molecular_weight_dalton, "molecular_weight_dalton", strictly_positive=True
    )
    provided = {
        "molar_concentration_molar": molar_concentration_molar,
        "mass_concentration_g_per_l": mass_concentration_g_per_l,
        "amount_mol": amount_mol,
        "mass_g": mass_g,
    }
    selected = tuple(name for name, value in provided.items() if value is not None)
    if len(selected) != 1:
        raise ConfigurationError(
            "Provide exactly one concentration or amount input.",
            code="OLIGO_QUANTITY_INPUT_COUNT",
            context={"provided_fields": selected},
        )
    input_kind = selected[0]
    raw_value = _finite(provided[input_kind], input_kind, minimum=0.0)
    resolved_volume = (
        None
        if volume_liter is None
        else _finite(volume_liter, "volume_liter", strictly_positive=True)
    )
    molar: float | None
    mass_concentration: float | None
    if (input_kind.endswith("per_l") or input_kind.endswith("molar")) and resolved_volume is None:
        raise ConfigurationError(
            "volume_liter is required when converting a concentration to an amount.",
            code="OLIGO_QUANTITY_VOLUME_REQUIRED",
        )
    if input_kind == "molar_concentration_molar":
        molar = raw_value
        mass_concentration = molar * molecular_weight
        assert resolved_volume is not None
        amount = molar * resolved_volume
        mass = amount * molecular_weight
    elif input_kind == "mass_concentration_g_per_l":
        mass_concentration = raw_value
        molar = mass_concentration / molecular_weight
        assert resolved_volume is not None
        mass = mass_concentration * resolved_volume
        amount = mass / molecular_weight
    elif input_kind == "amount_mol":
        amount = raw_value
        mass = amount * molecular_weight
        molar = None if resolved_volume is None else amount / resolved_volume
        mass_concentration = None if resolved_volume is None else mass / resolved_volume
    else:
        mass = raw_value
        amount = mass / molecular_weight
        molar = None if resolved_volume is None else amount / resolved_volume
        mass_concentration = None if resolved_volume is None else mass / resolved_volume
    return OligoQuantityResult(
        molecular_weight_dalton=molecular_weight,
        volume_liter=resolved_volume,
        molar_concentration_molar=molar,
        mass_concentration_g_per_l=mass_concentration,
        amount_mol=amount,
        amount_nmol=amount * 1_000_000_000.0,
        mass_g=mass,
        mass_microgram=mass * 1_000_000.0,
        input_kind=input_kind,
        method="molar-mass-volume-dimensional-conversion",
    )


__all__ = [
    "concentration_from_a260",
    "convert_oligo_quantity",
    "optical_properties",
]
