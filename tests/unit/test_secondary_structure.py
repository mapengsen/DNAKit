"""Tests for backend-neutral secondary structure and the conditional NUPACK adapter."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from dnakit.core import BackendInfo, DNASequence
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.secondary_structure import (
    NupackAdapter,
    analyze_dot_bracket,
    ensemble_defect_from_probabilities,
    pair_probability_metrics,
    probe_nupack,
    target_structure_probability,
)


def test_dot_bracket_reports_hairpin_stem_loop_and_maximum_pairing() -> None:
    summary = analyze_dot_bracket(
        (DNASequence("ATCCTAGTTATAGGAT"),),
        "((((((....))))))",
    )

    assert summary.structure_type == "hairpin"
    assert summary.base_pair_count == 6
    assert summary.stem_lengths == (6,)
    assert summary.hairpin_count == 1
    assert summary.hairpin_loop_lengths == (4,)
    assert summary.max_contiguous_pair_count == 6
    assert not summary.three_prime_dimer


def test_dot_bracket_reports_interstrand_and_three_prime_dimer() -> None:
    summary = analyze_dot_bracket(
        (DNASequence("CCC"), DNASequence("GGG")),
        "(((+)))",
    )

    assert summary.structure_type == "heterodimer"
    assert summary.stems[0].inter_strand
    assert summary.three_prime_dimer
    assert summary.three_prime_dimer_max_contiguous_pairs == 3

    unbound = analyze_dot_bracket(
        (DNASequence("CCC"), DNASequence("GGG")),
        "...+...",
    )
    assert unbound.structure_type == "unbound-strands"


def test_pair_probabilities_accessibility_defect_and_target_probability() -> None:
    matrix = (
        (0.7, 0.0, 0.0, 0.3),
        (0.0, 0.8, 0.2, 0.0),
        (0.0, 0.2, 0.8, 0.0),
        (0.3, 0.0, 0.0, 0.7),
    )
    probabilities = pair_probability_metrics(
        (DNASequence("ACGT"),),
        matrix,
        accessibility_window_size=2,
    )
    target = analyze_dot_bracket((DNASequence("ACGT"),), "(())")

    assert probabilities.pairing_probabilities_by_base == pytest.approx((0.3, 0.2, 0.2, 0.3))
    assert probabilities.unpaired_probabilities_by_base == pytest.approx((0.7, 0.8, 0.8, 0.7))
    assert probabilities.most_accessible_window_start == 1
    assert ensemble_defect_from_probabilities(target, probabilities) == pytest.approx(0.75)
    assert target_structure_probability(-2.0, -2.0) == 1.0
    assert 0.0 < target_structure_probability(-1.0, -2.0) < 1.0

    mismatched_target = analyze_dot_bracket((DNASequence("TGCA"),), "(())")
    with pytest.raises(ConfigurationError) as mismatch:
        ensemble_defect_from_probabilities(mismatched_target, probabilities)
    assert mismatch.value.code == "ENSEMBLE_DEFECT_STRAND_MISMATCH"


def test_dot_bracket_and_probability_input_validation_is_explicit() -> None:
    with pytest.raises(ConfigurationError) as structure_error:
        analyze_dot_bracket((DNASequence("ACGT"),), "((..)")
    assert structure_error.value.code == "DOT_BRACKET_LENGTH_MISMATCH"

    with pytest.raises(ConfigurationError) as probability_error:
        pair_probability_metrics(
            (DNASequence("AC"),),
            ((0.5, 0.4), (0.4, 0.5)),
        )
    assert probability_error.value.code == "INVALID_PAIR_PROBABILITY_ROW_SUM"


def test_nupack_probe_never_installs_or_downloads() -> None:
    info = probe_nupack()

    assert info.name == "nupack"
    assert info.metadata["automatic_install"] is False
    assert info.metadata["automatic_download"] is False
    assert info.metadata["import_executed"] is False
    if not info.available:
        with pytest.raises(BackendUnavailableError) as error:
            NupackAdapter().analyze_complex((DNASequence("ACGT"),))
        assert error.value.code == "NUPACK_UNAVAILABLE"


@dataclass
class _FakeStructure:
    structure: str
    energy: float
    stack_energy: float | None


class _FakeArray:
    def tolist(self) -> list[list[float]]:
        return [
            [0.8, 0.0, 0.0, 0.2],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.2, 0.0, 0.0, 0.8],
        ]


class _FakePairMatrix:
    def to_array(self) -> _FakeArray:
        return _FakeArray()


@dataclass(frozen=True)
class _FakeStrand:
    sequence: str
    name: str


class _FakeComplex:
    def __init__(self, strands: list[_FakeStrand], *, name: str) -> None:
        self.strands = tuple(strands)
        self.name = name


@dataclass(frozen=True)
class _FakeSetSpec:
    max_size: int
    include: list[_FakeComplex]


@dataclass(eq=False)
class _FakeTube:
    strands: dict[_FakeStrand, float]
    complexes: _FakeSetSpec
    name: str


@dataclass(frozen=True)
class _FakeTubeResult:
    complex_concentrations: dict[_FakeComplex, float]
    fraction_bases_unpaired: float


class _FakePartition:
    def log(self) -> float:
        return 2.0


class _FakeNupack:
    Strand = _FakeStrand
    Complex = _FakeComplex
    SetSpec = _FakeSetSpec
    Tube = _FakeTube

    @staticmethod
    def Model(**parameters: object) -> object:
        return parameters

    @staticmethod
    def pfunc(**parameters: object) -> tuple[_FakePartition, float]:
        del parameters
        return (_FakePartition(), -2.0)

    @staticmethod
    def mfe(**parameters: object) -> list[_FakeStructure]:
        del parameters
        return [_FakeStructure("....", -1.5, None)]

    @staticmethod
    def pairs(**parameters: object) -> _FakePairMatrix:
        del parameters
        return _FakePairMatrix()

    @staticmethod
    def subopt(**parameters: object) -> list[_FakeStructure]:
        del parameters
        return [_FakeStructure("....", -1.5, None)]

    @staticmethod
    def sample(**parameters: object) -> list[str]:
        count = parameters["num_sample"]
        assert isinstance(count, int)
        return ["...."] * count

    @staticmethod
    def ensemble_size(**parameters: object) -> int:
        del parameters
        return 7

    @staticmethod
    def structure_probability(**parameters: object) -> float:
        del parameters
        return 0.25

    @staticmethod
    def defect(**parameters: object) -> float:
        del parameters
        return 0.1

    @staticmethod
    def tube_analysis(**parameters: object) -> dict[_FakeTube, _FakeTubeResult]:
        tubes = parameters["tubes"]
        assert isinstance(tubes, list) and len(tubes) == 1
        tube = tubes[0]
        assert isinstance(tube, _FakeTube)
        target = tube.complexes.include[0]
        first_strand = next(iter(tube.strands))
        monomer = _FakeComplex([first_strand], name=f"({first_strand.name})")
        return {
            tube: _FakeTubeResult(
                complex_concentrations={target: 4e-7, monomer: 1e-7},
                fraction_bases_unpaired=0.2,
            )
        }


def test_nupack_adapter_maps_bounded_outputs_without_real_installation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adapter = NupackAdapter()
    adapter._info = BackendInfo(
        "nupack",
        version="4-test",
        license_expression="LicenseRef-NUPACK",
        capabilities=probe_nupack().capabilities,
        available=True,
    )
    monkeypatch.setattr(adapter, "_module", lambda: _FakeNupack())

    result = adapter.analyze_complex(
        (DNASequence("ACGT"),),
        target_structure="....",
        num_samples=2,
        accessibility_window_size=2,
    )

    assert result.partition_function_log == 2.0
    assert result.material == "dna"
    assert result.ensemble_free_energy_kcal_per_mol == -2.0
    assert result.mfe_structures[0].summary.dot_bracket == "...."
    assert result.boltzmann_samples == ("....", "....")
    assert result.ensemble_size == 7
    assert result.target_structure_probability == 0.25
    assert result.target_ensemble_defect == 0.1
    assert result.backend.version == "4-test"

    tube = adapter.analyze_tube(
        {"a": DNASequence("CCC"), "b": DNASequence("GGG")},
        {"a": 1e-6, "b": 1e-6},
        target_strand_names=("a", "b"),
    )
    assert tube.sequences_5to3 == ("CCC", "GGG")
    assert tube.target_strand_names == ("a", "b")
    assert tube.target_complex_concentration_molar == pytest.approx(4e-7)
    assert tube.complex_fraction_denominator_molar == pytest.approx(5e-7)
    assert tube.target_complex_fraction == pytest.approx(0.8)
    assert tube.non_target_complex_fraction == pytest.approx(0.2)
    assert tube.fraction_bases_unpaired == 0.2


def test_nupack_adapter_rejects_invalid_user_numbers_as_configuration_errors() -> None:
    adapter = NupackAdapter()
    with pytest.raises(ConfigurationError) as error:
        adapter.analyze_complex(
            (DNASequence("ACGT"),),
            suboptimal_energy_gap_kcal_per_mol=float("nan"),
        )
    assert error.value.code == "INVALID_NUPACK_CONFIGURATION"
