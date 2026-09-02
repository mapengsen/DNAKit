"""Tests for bounded native thermodynamic calculations."""

from __future__ import annotations

import json
import math
import stat
import sys
from pathlib import Path

import pytest

from dnakit.core import (
    BackendInfo,
    DNAAlphabet,
    DNASequence,
    ExecutionMode,
    Gap,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Topology,
)
from dnakit.core._json import FrozenDict
from dnakit.core.provenance import Provenance
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    InvalidAlphabetError,
    UnsupportedGapOperationError,
)
from dnakit.thermodynamics import (
    NearestNeighborConfig,
    Primer3CLIAdapter,
    Primer3ThermodynamicResult,
    ThermodynamicConditions,
    duplex_stability,
    extinction_coefficient_260nm,
    melting_temperature,
    molecular_weight,
    nearest_neighbor,
    salt_correction,
    stacking_interactions,
    window_tm,
)


def test_molecular_weight_uses_explicit_anhydrous_residue_table() -> None:
    single = molecular_weight(DNASequence("ACGT"))
    double = molecular_weight(DNASequence("A"), strand="double")
    phosphorylated = molecular_weight(DNASequence("A"), five_prime_phosphorylated=True)

    assert single.value_dalton == pytest.approx(1173.84)
    assert single.value_kilodalton == pytest.approx(1.17384)
    assert double.value_dalton == pytest.approx((313.21 - 61.96) + (304.20 - 61.96))
    assert phosphorylated.value_dalton == pytest.approx(313.21 - 61.96 + 79.0)
    assert single.parameters["base_masses_dalton"] == {
        "A": 313.21,
        "C": 289.18,
        "G": 329.21,
        "T": 304.2,
    }


def test_extinction_coefficient_260nm_matches_published_acgt_example() -> None:
    result = extinction_coefficient_260nm(DNASequence("ACGT"))

    assert result.value_m_inverse_cm_inverse == pytest.approx(40_300.0)
    assert result.wavelength_nm == 260
    assert result.sequence_length == 4
    pair_parameters = result.parameters["nearest_neighbor_m_inverse_cm_inverse"]
    base_parameters = result.parameters["individual_base_m_inverse_cm_inverse"]
    assert isinstance(pair_parameters, FrozenDict)
    assert isinstance(base_parameters, FrozenDict)
    assert pair_parameters["AC"] == 21_200.0
    assert base_parameters["C"] == 7_400.0
    assert tuple(citation.key for citation in result.provenance.implementation.citations) == (
        "warshaw-tinoco1966",
        "cantor-warshaw-shapiro1970",
    )
    assert json.loads(json.dumps(result.to_dict()))["value_m_inverse_cm_inverse"] == 40_300.0


def test_extinction_coefficient_handles_sequence_order_and_short_oligos() -> None:
    assert extinction_coefficient_260nm(
        DNASequence("TGCA")
    ).value_m_inverse_cm_inverse == pytest.approx(38_900.0)
    assert extinction_coefficient_260nm(
        DNASequence("AA")
    ).value_m_inverse_cm_inverse == pytest.approx(27_400.0)
    assert extinction_coefficient_260nm(
        DNASequence("A")
    ).value_m_inverse_cm_inverse == pytest.approx(15_400.0)


def test_extinction_coefficient_rejects_double_stranded_and_invalid_limits() -> None:
    with pytest.raises(ConfigurationError) as strand_error:
        extinction_coefficient_260nm(DNASequence("ACGT", strandedness="double"))
    assert strand_error.value.code == "EXTINCTION_COEFFICIENT_SINGLE_STRAND_ONLY"

    for invalid_limit in (True, 0, 10_000_001):
        with pytest.raises(ConfigurationError) as limit_error:
            extinction_coefficient_260nm(
                DNASequence("A"),
                max_sequence_length=invalid_limit,
            )
        assert limit_error.value.code == "INVALID_THERMODYNAMIC_LIMIT"


def test_wallace_rule_is_bounded_and_records_unmodeled_conditions() -> None:
    result = melting_temperature(DNASequence("AACG"), method="wallace")

    assert result.tm_celsius == 12.0
    assert result.sequence_length == 4
    assert "conditions are not modeled" in result.applicability
    with pytest.raises(ConfigurationError) as error:
        melting_temperature(DNASequence("A" * 14), method="wallace")
    assert error.value.code == "THERMODYNAMIC_SEQUENCE_TOO_LONG"
    with pytest.raises(ConfigurationError) as config_error:
        melting_temperature(
            DNASequence("ACGT"),
            method="wallace",
            config=NearestNeighborConfig(),
        )
    assert config_error.value.code == "TM_CONFIG_NOT_APPLICABLE"


def test_santalucia_example_has_published_stacking_and_total_parameters() -> None:
    conditions = ThermodynamicConditions(
        temperature_celsius=37.0,
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    )
    result = nearest_neighbor(DNASequence("GTGCAT"), conditions=conditions)

    assert result.stacking_delta_h_kcal_per_mol == pytest.approx(-42.4)
    assert result.stacking_delta_s_cal_per_k_mol == pytest.approx(-112.6)
    assert result.initiation_delta_h_kcal_per_mol == pytest.approx(2.4)
    assert result.initiation_delta_s_cal_per_k_mol == pytest.approx(1.3)
    assert result.delta_h_kcal_per_mol == pytest.approx(-40.0)
    assert result.delta_s_cal_per_k_mol == pytest.approx(-111.3)
    assert result.salt_delta_s_cal_per_k_mol == 0.0
    assert [step.top_5to3 for step in result.stacking_steps] == [
        "GT",
        "TG",
        "GC",
        "CA",
        "AT",
    ]


def test_nearest_neighbor_self_complement_symmetry_and_units() -> None:
    result = nearest_neighbor(
        DNASequence("ACGT"),
        conditions=ThermodynamicConditions(sodium_molar=1.0),
    )

    assert result.self_complementary
    assert result.concentration_divisor == 1
    assert result.symmetry_delta_h_kcal_per_mol == 0.0
    assert result.symmetry_delta_s_cal_per_k_mol == -1.4
    assert result.gas_constant_cal_per_k_mol == 1.9872
    assert result.tm_equation == "1000*dH/(dS+R*ln(Ct/divisor))-273.15"
    expected_g = result.delta_h_kcal_per_mol - (310.15 * result.delta_s_cal_per_k_mol / 1000.0)
    assert result.delta_g_kcal_per_mol == pytest.approx(expected_g)
    assert json.loads(json.dumps(result.to_dict()))["parameter_set"] == "santalucia1998-v1"


def test_salt_correction_matches_published_monovalent_entropy_equation() -> None:
    conditions = ThermodynamicConditions(sodium_molar=0.05)
    result = salt_correction(10, conditions=conditions)

    assert result.delta_s_cal_per_k_mol == pytest.approx(0.368 * 9 * math.log(0.05))
    assert result.conditions.sodium_molar == 0.05
    assert result.model_version == "santalucia1998-equation-v1"


def test_stacking_result_excludes_terminal_symmetry_and_salt_terms() -> None:
    result = stacking_interactions(DNASequence("AA"), temperature_celsius=37.0)

    assert len(result.steps) == 1
    assert result.steps[0].bottom_3to5 == "TT"
    assert result.total_delta_h_kcal_per_mol == -7.9
    assert result.total_delta_s_cal_per_k_mol == -22.2
    assert result.total_delta_g_kcal_per_mol == pytest.approx(-1.01467)
    assert "excludes initiation" in result.applicability


def test_duplex_stability_requires_exact_reverse_complement() -> None:
    stable = duplex_stability(DNASequence("GTGCAT"), DNASequence("ATGCAC"))
    assert stable.fully_complementary
    assert stable.delta_g_kcal_per_mol == stable.thermodynamics.delta_g_kcal_per_mol
    assert stable.stable_at_temperature == (
        stable.tm_celsius > stable.conditions.temperature_celsius
    )
    assert stable.stability_criterion == "tm_celsius > configured_temperature_celsius"

    with pytest.raises(ConfigurationError) as mismatch:
        duplex_stability(DNASequence("GTGCAT"), DNASequence("ATGCAA"))
    assert mismatch.value.code == "MISMATCHED_DUPLEX_UNSUPPORTED"


def test_duplex_stability_explicit_primer3_backend_supports_mismatch(
    tmp_path: Path,
) -> None:
    executable = tmp_path / "ntthal"
    executable.write_text(
        (
            f"#!{sys.executable}\n"
            "print('Calculated thermodynamical parameters for dimer:\\t'"
            "      'dS = -100\\tdH = -30000\\tdG = -2500\\tt = 55')\n"
            "print('SEQ\\tACGT')\n"
        ),
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    adapter = Primer3CLIAdapter(ntthal_path=executable)

    result = duplex_stability(
        DNASequence("GTGCAT"),
        DNASequence("ATGCAA"),
        backend="primer3-cli",
        adapter=adapter,
    )

    assert not result.fully_complementary
    assert result.model == "dnakit-primer3-cli-adapter-v1"
    assert isinstance(result.thermodynamics, Primer3ThermodynamicResult)
    assert result.thermodynamics.capability == "heterodimer"
    assert result.provenance.implementation.label == "adapter"

    dangling = duplex_stability(
        DNASequence("GTGCAT"),
        DNASequence("TGCAC"),
        backend="primer3-cli",
        adapter=adapter,
        output_structure=True,
    )
    assert not dangling.fully_complementary
    assert isinstance(dangling.thermodynamics, Primer3ThermodynamicResult)
    assert dangling.thermodynamics.ascii_structure is not None


def test_duplex_stability_backend_arguments_are_strict() -> None:
    with pytest.raises(ConfigurationError, match="backend must"):
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAC"),
            backend="unknown",  # type: ignore[arg-type]
        )


def test_duplex_stability_rejects_mismatched_adapter_result() -> None:
    class MismatchedAdapter(Primer3CLIAdapter):
        def heterodimer(
            self,
            sequence_a: DNASequence,
            sequence_b: DNASequence,
            **kwargs: object,
        ) -> Primer3ThermodynamicResult:
            del sequence_a, sequence_b, kwargs
            info = BackendInfo(
                "primer3-cli",
                version="test",
                license_expression="GPL-2.0-or-later",
                capabilities=("heterodimer",),
            )
            return Primer3ThermodynamicResult(
                capability="heterodimer",
                sequences_5to3=("AAAA", "TTTT"),
                structure_found=True,
                tm_celsius=1.0,
                delta_g_kcal_per_mol=-1.0,
                delta_h_kcal_per_mol=-1.0,
                delta_s_cal_per_k_mol=-1.0,
                ascii_structure=None,
                conditions=ThermodynamicConditions(),
                max_loop=30,
                method="fixture",
                algorithm_version="fixture-v1",
                backend=info,
                provenance=Provenance(
                    implementation=ImplementationInfo(
                        label=ImplementationLabel.ADAPTER,
                        execution_mode=ExecutionMode.EXTERNAL,
                        origin_class=OriginClass.INTEGRATION,
                    ),
                    backend=info,
                ),
                issues=(),
            )

    with pytest.raises(BackendExecutionError) as error:
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAA"),
            backend="primer3-cli",
            adapter=MismatchedAdapter(),
        )
    assert error.value.code == "MISMATCHED_PRIMER3_DUPLEX_RESULT"
    with pytest.raises(ConfigurationError, match="only valid"):
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAC"),
            adapter=Primer3CLIAdapter(),
        )
    with pytest.raises(ConfigurationError, match="max_loop"):
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAC"),
            max_loop=0,
        )
    with pytest.raises(ConfigurationError, match="only valid"):
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAC"),
            output_structure=True,
        )
    with pytest.raises(ConfigurationError, match="only valid"):
        duplex_stability(
            DNASequence("GTGCAT"),
            DNASequence("ATGCAC"),
            backend="primer3-cli",
            config=NearestNeighborConfig(),
        )


def test_window_tm_coordinates_counts_bounds_and_resource_limit() -> None:
    wallace = window_tm(DNASequence("AACGTT"), 4, step=2, method="wallace")

    assert [(item.start, item.end, item.sequence) for item in wallace.windows] == [
        (0, 4, "AACG"),
        (2, 6, "CGTT"),
    ]
    assert [item.tm_celsius for item in wallace.windows] == [12.0, 12.0]
    assert wallace.coordinate_system == "0-based-half-open-symbol"

    with pytest.raises(ConfigurationError) as resource_error:
        window_tm(DNASequence("A" * 100), 2, max_windows=10)
    assert resource_error.value.code == "WINDOW_TM_LIMIT_EXCEEDED"
    with pytest.raises(ConfigurationError):
        window_tm(DNASequence("AAAA"), 2, max_windows=1_000_001)


@pytest.mark.parametrize(
    "sequence",
    [
        DNASequence("AN", alphabet=DNAAlphabet.IUPAC),
        DNASequence(["A", Gap(1), "T"]),
        DNASequence("AT", topology=Topology.CIRCULAR),
    ],
)
def test_native_calculations_reject_out_of_domain_sequence_states(
    sequence: DNASequence,
) -> None:
    error_type = (
        InvalidAlphabetError
        if sequence.ambiguity_count
        else UnsupportedGapOperationError
        if sequence.is_gapped
        else ConfigurationError
    )
    with pytest.raises(error_type):
        nearest_neighbor(sequence)


def test_empty_too_long_ions_and_invalid_runtime_configs_are_rejected() -> None:
    with pytest.raises(ConfigurationError):
        molecular_weight(DNASequence(""))
    with pytest.raises(ConfigurationError) as long_error:
        nearest_neighbor(DNASequence("A" * 61))
    assert long_error.value.code == "THERMODYNAMIC_SEQUENCE_TOO_LONG"
    with pytest.raises(ConfigurationError) as ion_error:
        nearest_neighbor(
            DNASequence("ACGT"),
            conditions=ThermodynamicConditions(magnesium_molar=0.001),
        )
    assert ion_error.value.code == "UNSUPPORTED_NATIVE_ION_MODEL"

    invalid_configs: tuple[object, ...] = ({}, [], 0, False)
    for value in invalid_configs:
        with pytest.raises(ConfigurationError):
            nearest_neighbor(DNASequence("AC"), config=value)  # type: ignore[arg-type]
        with pytest.raises(ConfigurationError):
            nearest_neighbor(DNASequence("AC"), conditions=value)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "factory",
    [
        lambda: ThermodynamicConditions(temperature_celsius=float("nan")),
        lambda: ThermodynamicConditions(sodium_molar=float("inf")),
        lambda: ThermodynamicConditions(strand_concentration_molar=0.0),
        lambda: ThermodynamicConditions(temperature_celsius=-273.15),
        lambda: NearestNeighborConfig(max_sequence_length=True),
        lambda: NearestNeighborConfig(max_sequence_length=1),
        lambda: NearestNeighborConfig(max_sequence_length=61),
    ],
)
def test_non_finite_and_out_of_domain_conditions_are_rejected(factory: object) -> None:
    with pytest.raises(ConfigurationError):
        factory()  # type: ignore[operator]
