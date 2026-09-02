"""Tests for optical conversions and ideal duplex equilibrium extensions."""

from __future__ import annotations

from itertools import pairwise

import pytest

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.thermodynamics import (
    LabelAbsorbanceCorrection,
    OpticalModification,
    ThermodynamicConditions,
    binding_equilibrium,
    concentration_from_a260,
    convert_oligo_quantity,
    cosolvent_tm_correction,
    nearest_neighbor,
    optical_properties,
    terminal_stability,
    theoretical_melting_curve,
)


def test_single_and_double_strand_optical_properties_have_explicit_units() -> None:
    single = optical_properties(DNASequence("ACGT"))
    double = optical_properties(DNASequence("ACGT"), strand_type="double")

    assert single.extinction_coefficient_260_m_inverse_cm_inverse == pytest.approx(40_300.0)
    assert single.molecular_weight_dalton == pytest.approx(1_173.84)
    assert single.one_od260_nmol == pytest.approx(1_000_000.0 / 40_300.0)
    assert single.one_od260_microgram == pytest.approx(1_000.0 * 1_173.84 / 40_300.0)
    assert double.sequences_5to3 == ("ACGT", "ACGT")
    assert double.extinction_coefficient_260_m_inverse_cm_inverse == 52_800.0
    assert double.molecular_weight_dalton == pytest.approx(2_347.68)


def test_optical_modification_and_label_absorbance_corrections_are_explicit() -> None:
    properties = optical_properties(
        DNASequence("ACGT"),
        modifications=(
            OpticalModification(
                "fluorophore",
                count=2,
                extinction_coefficient_260_delta_m_inverse_cm_inverse=1_000.0,
                molecular_weight_delta_dalton=100.0,
            ),
        ),
    )
    concentration = concentration_from_a260(
        0.5,
        properties,
        label_corrections=(
            LabelAbsorbanceCorrection(
                "fluorophore",
                absorbance_at_label_max=0.2,
                a260_correction_factor=0.1,
            ),
        ),
    )

    assert properties.extinction_coefficient_260_m_inverse_cm_inverse == 42_300.0
    assert properties.molecular_weight_dalton == pytest.approx(1_373.84)
    assert concentration.label_a260_subtracted == pytest.approx(0.02)
    assert concentration.corrected_a260 == pytest.approx(0.48)
    assert concentration.molar_concentration_molar == pytest.approx(0.48 / 42_300.0)


def test_a260_and_quantity_conversions_round_trip() -> None:
    properties = optical_properties(DNASequence("ACGT"))
    measured = concentration_from_a260(
        0.403,
        properties,
        volume_liter=0.001,
    )
    converted = convert_oligo_quantity(
        properties.molecular_weight_dalton,
        volume_liter=0.001,
        molar_concentration_molar=measured.molar_concentration_molar,
    )
    from_mass = convert_oligo_quantity(
        properties.molecular_weight_dalton,
        volume_liter=0.001,
        mass_g=converted.mass_g,
    )

    assert measured.molar_concentration_micromolar == pytest.approx(10.0)
    assert measured.mass_concentration_ng_per_microliter == pytest.approx(11.7384)
    assert measured.amount_mol == pytest.approx(1e-8)
    assert converted.amount_nmol == pytest.approx(10.0)
    assert converted.mass_microgram == pytest.approx(11.7384)
    assert from_mass.molar_concentration_molar == pytest.approx(1e-5)


def test_quantity_conversion_requires_one_source_and_volume_for_concentration() -> None:
    with pytest.raises(ConfigurationError) as missing_volume:
        convert_oligo_quantity(1_000.0, molar_concentration_molar=1e-6)
    assert missing_volume.value.code == "OLIGO_QUANTITY_VOLUME_REQUIRED"

    with pytest.raises(ConfigurationError) as too_many:
        convert_oligo_quantity(1_000.0, amount_mol=1e-9, mass_g=1e-6)
    assert too_many.value.code == "OLIGO_QUANTITY_INPUT_COUNT"


def test_binding_constants_and_two_state_melting_curve_are_consistent() -> None:
    conditions = ThermodynamicConditions(
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    )
    equilibrium = binding_equilibrium(DNASequence("GTGCAT"), conditions=conditions)
    progress: list[tuple[int, int]] = []
    curve = theoretical_melting_curve(
        DNASequence("GTGCAT"),
        range(0, 51, 5),
        conditions=conditions,
        progress=lambda completed, total: progress.append((completed, total)),
    )
    nearest = nearest_neighbor(DNASequence("GTGCAT"), conditions=conditions)

    assert equilibrium.association_constant_m_inverse * equilibrium.dissociation_constant_molar == (
        pytest.approx(1.0)
    )
    assert 0.0 <= equilibrium.duplex_fraction <= 1.0
    assert all(
        left.duplex_fraction > right.duplex_fraction for left, right in pairwise(curve.points)
    )
    assert curve.midpoint_temperature_celsius == pytest.approx(nearest.tm_celsius, abs=0.1)
    assert progress == [(index, 11) for index in range(1, 12)]


def test_terminal_cosolvent_and_total_monovalent_corrections() -> None:
    sodium = nearest_neighbor(
        DNASequence("GTGCAT"),
        conditions=ThermodynamicConditions(sodium_molar=0.05, potassium_molar=0.0),
    )
    potassium = nearest_neighbor(
        DNASequence("GTGCAT"),
        conditions=ThermodynamicConditions(sodium_molar=0.0, potassium_molar=0.05),
    )
    terminal = terminal_stability(
        DNASequence("AACCGGTT"),
        conditions=ThermodynamicConditions(sodium_molar=1.0),
    )
    corrected = cosolvent_tm_correction(
        DNASequence("ACGT"),
        60.0,
        dmso_percent=5.0,
        formamide_molar=1.0,
    )

    assert potassium.tm_celsius == pytest.approx(sodium.tm_celsius)
    assert terminal.less_stable_end == "equal"
    assert terminal.five_prime_sequence == "AACCG"
    assert terminal.three_prime_sequence == "CGGTT"
    assert corrected.dmso_delta_tm_celsius == -3.0
    assert corrected.formamide_delta_tm_celsius == pytest.approx(-2.6535)
    assert corrected.corrected_tm_celsius == pytest.approx(54.3465)
