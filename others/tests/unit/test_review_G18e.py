"""Contract tests for the dadi and tskit review wrappers in dnakit.population (group G18e).

Each wrapped upstream function is checked against a hand calculation or a direct call to the
upstream dadi function with the same inputs. Inputs are written inside tmp_path, and no network is
used. The tskit module has no wrapped function; its test only checks that it exports nothing.
"""

from __future__ import annotations

import inspect
import itertools
import random
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.population import dadi_review_functions as dadi_review
from dnakit.population import tskit_review_functions as tskit_review

np = pytest.importorskip("numpy")
dadi = pytest.importorskip("dadi")

SAMPLE_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
chr1\t100\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|0\t0|1\t1|1\t0|0
chr1\t200\t.\tC\tT\t.\tPASS\tAA=C\tGT\t0|1\t1|1\t0|0\t0|0
chr1\t300\t.\tG\tA\t.\tPASS\tAA=G\tGT\t0|0\t0|0\t0|1\t1|1
chr1\t400\t.\tT\tC\t.\tPASS\tAA=T\tGT\t1|0\t0|1\t0|0\t1|1
"""
POPINFO = "s1 popA\ns2 popA\ns3 popB\ns4 popB\n"
SLIM_POP1 = """#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 2 m1 250 0.2
2 5 m2 400 0.0
Genomes:
p1:0 A 0 1
p1:1 A 1
p1:2 A
"""
SLIM_POP2 = """#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 3 m1 300 0.2
Genomes:
p2:0 A 0
p2:1 A 1
p2:2 A 0 1
"""
JC_RATES = [[0.0 if i == j else 1 / 3 for j in range(4)] for i in range(4)]
TRI_KEYS = ["".join(p) for p in itertools.product("ACGT", repeat=3)]


def _expected_counts(pts: list[int], ns: tuple[int, ...], params: list[float]) -> Any:
    """Expected 1D spectrum of a two-epoch model, scaled to counts."""

    model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
    return np.ma.filled(model(tuple(params), ns, pts) * 200, 0.0)


def _noisy_counts(rng: Any, mean: Any) -> list[float]:
    return [float(value) for value in rng.poisson(mean).astype(float)]


def _data_and_boots(num_boot: int = 4) -> tuple[list[float], list[list[float]]]:
    rng = np.random.default_rng(7)
    mean = _expected_counts([20, 30, 40], (10,), [1.5, 0.5])
    data = _noisy_counts(rng, mean)
    boots = [_noisy_counts(rng, mean) for _ in range(num_boot)]
    return data, boots


def _three_epoch_data() -> tuple[list[float], list[list[float]]]:
    """Counts from a three-epoch model, for the nested tests."""

    rng = np.random.default_rng(3)
    model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.three_epoch)
    mean = np.ma.filled(model((2.0, 0.5, 0.2, 0.3), (10,), [20, 30, 40]) * 200, 0.0)
    return _noisy_counts(rng, mean), [_noisy_counts(rng, mean) for _ in range(4)]


def test_public_functions_take_keyword_only_verbose() -> None:
    for name in dadi_review.__all__:
        signature = inspect.signature(getattr(dadi_review, name))
        verbose = signature.parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True
    assert tskit_review.__all__ == []


def test_total_substitution_rate_hand_calculation() -> None:
    jc = dadi_review.total_substitution_rate(JC_RATES, [0.25] * 4)
    assert jc["status"] == "completed"
    assert jc["total_rate"] == pytest.approx(1.0)
    uniform = [[0.0 if i == j else 1 / 63 for j in range(64)] for i in range(64)]
    result = dadi_review.total_substitution_rate(uniform, [1 / 64] * 64)
    assert result["num_states"] == 64
    assert result["total_rate"] == pytest.approx(1.0)


def test_total_substitution_rate_rejects_bad_input() -> None:
    with pytest.raises(ConfigurationError):
        dadi_review.total_substitution_rate([[0.0, 1.0], [1.0, 0.0]], [0.5, 0.5])
    with pytest.raises(ConfigurationError):
        dadi_review.total_substitution_rate(JC_RATES, [0.5, 0.5, 0.5, 0.5])


def test_ms_command_line_hand_strings() -> None:
    single = dadi_review.ms_command_line(2.0, [5], "-eN 0.5 2.0", 10)
    assert single["command"] == "ms 5 10 -t 2.000000  -eN 0.5 2.0"
    multi = dadi_review.ms_command_line(
        2.0, [3, 4], "-eN 0.5 2.0", 10, recombination=1.0, seeds=(1, 2, 3)
    )
    assert multi["command"] == (
        "ms 7 10 -t 2.000000 -I 2 3 4 -eN 0.5 2.0 -r 1.000000 20 -seeds 1 2 3"
    )
    with pytest.raises(ConfigurationError):
        dadi_review.ms_command_line(2.0, [5], "-eN 0.5\n2.0", 10)


def test_ms_core_command_2d_matches_upstream() -> None:
    split = dadi_review.ms_core_command_2d("split_mig", [1.0, 2.0, 0.5, 0.25])
    assert split["core"] == (
        "-n 1 1.000000 -n 2 2.000000 -ma x 0.500000 0.500000 x -ej 0.250000 2 1 -en 0.250000 1 1"
    )
    im = dadi_review.ms_core_command_2d("IM", [0.4, 1.0, 2.0, 0.5, 0.1, 0.2])
    assert im["core"] == dadi.Demographics2D.IM_mscore((0.4, 1.0, 2.0, 0.5, 0.1, 0.2))
    pre = dadi_review.ms_core_command_2d("IM_pre", [1.5, 0.3, 0.4, 1.0, 2.0, 0.5, 0.1, 0.2])
    assert pre["core"] == dadi.Demographics2D.IM_pre_mscore(
        (1.5, 0.3, 0.4, 1.0, 2.0, 0.5, 0.1, 0.2)
    )
    with pytest.raises(ConfigurationError):
        dadi_review.ms_core_command_2d("IM", [0.0, 1.0, 2.0, 0.5, 0.1, 0.2])


def test_expected_neutral_sfs_2d_matches_upstream() -> None:
    result = dadi_review.expected_neutral_sfs_2d((3, 3), 20)
    direct = dadi.Demographics2D.snm(None, (3, 3), 20)
    assert len(result["spectrum"]) == 4
    assert result["spectrum"][0][0] is None
    assert result["spectrum"][3][3] is None
    assert result["spectrum"][0][1] == pytest.approx(float(direct[0, 1]))
    assert result["spectrum"][1][2] == pytest.approx(float(direct[1, 2]))


def test_trinucleotide_misid_table_matches_upstream(tmp_path: Path) -> None:
    rates = [[0.0 if i == j else 1 / 63 for j in range(64)] for i in range(64)]
    frequencies = {key: 1 / 64 for key in TRI_KEYS}
    target = tmp_path / "misid.txt"
    result = dadi_review.trinucleotide_misid_table(
        target,
        expected_substitutions=0.1,
        rate_matrix=rates,
        trinucleotide_frequencies=frequencies,
    )
    assert result["rows"] == 256
    direct = tmp_path / "direct.txt"
    dadi.Misc.make_fux_table(str(direct), 0.1, np.array(rates, dtype=float), frequencies)
    assert target.read_text() == direct.read_text()
    values = [float(line.split()[-1]) for line in target.read_text().splitlines()]
    assert all(0.0 <= value <= 1.0 for value in values)
    with pytest.raises(ConfigurationError):
        dadi_review.trinucleotide_misid_table(
            target,
            expected_substitutions=0.1,
            rate_matrix=rates,
            trinucleotide_frequencies=frequencies,
        )


def test_slim_output_data_dict_matches_upstream(tmp_path: Path) -> None:
    pop1 = tmp_path / "p1.slimout"
    pop2 = tmp_path / "p2.slimout"
    pop1.write_text(SLIM_POP1)
    pop2.write_text(SLIM_POP2)
    result = dadi_review.slim_output_data_dict([pop1, pop2])
    assert result["populations"] == ["0", "1"]
    assert result["sample_sizes"] == [3, 3]
    assert result["num_sites"] == 4
    first = result["sites"][0]
    assert first["key"] == "SLiM_100.0"
    assert first["calls"] == {"0": [2, 1], "1": [1, 2]}
    with open(pop1) as one, open(pop2) as two:
        direct, _ = dadi.Misc.dd_from_SLiM_files([one, two])
    assert {site["key"] for site in result["sites"]} == set(direct)
    filtered = dadi_review.slim_output_data_dict([pop1, pop2], mutation_types=["m2"])
    assert filtered["num_sites"] == 1
    assert filtered["sites"][0]["key"] == "SLiM_400.5"


def _write_vcf_inputs(tmp_path: Path) -> tuple[Path, Path]:
    vcf = tmp_path / "variants.vcf"
    popinfo = tmp_path / "popinfo.txt"
    vcf.write_text(SAMPLE_VCF)
    popinfo.write_text(POPINFO)
    return vcf, popinfo


def test_vcf_subsample_bootstraps_are_seeded_and_match_upstream(tmp_path: Path) -> None:
    vcf, popinfo = _write_vcf_inputs(tmp_path)
    kwargs = {
        "subsample": {"popA": 2, "popB": 2},
        "pop_ids": ["popA", "popB"],
        "num_boot": 2,
        "chunk_size": 200,
        "seed": 3,
    }
    python_state = random.getstate()
    numpy_state = np.random.get_state()
    first = dadi_review.vcf_subsample_bootstrap_spectra(vcf, popinfo, **kwargs)
    assert random.getstate() == python_state
    assert np.array_equal(np.random.get_state()[1], numpy_state[1])
    second = dadi_review.vcf_subsample_bootstrap_spectra(vcf, popinfo, **kwargs)
    assert first["spectra"] == second["spectra"]
    assert len(first["spectra"]) == 2
    assert len(first["spectra"][0]) == 5
    random.seed(3)
    np.random.seed(3)
    direct = dadi.Misc.bootstraps_subsample_vcf(
        str(vcf), str(popinfo), {"popA": 2, "popB": 2}, 2, 200, ["popA", "popB"]
    )
    expected = [[None if np.ma.is_masked(v) else float(v) for v in row] for row in direct[0]]
    assert first["spectra"][0] == expected


def test_godambe_and_gim_match_upstream() -> None:
    data, boots = _data_and_boots()
    model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
    grids = [20, 30, 40]
    observed = dadi.Spectrum(np.array(data))
    spectra = [dadi.Spectrum(np.array(item)) for item in boots]
    gim = dadi_review.parameter_uncertainty_gim(
        data, boots, model="two_epoch", p0=[1.5, 0.5], grid_points=grids
    )
    direct = dadi.Godambe.GIM_uncert(model, grids, spectra, [1.5, 0.5], observed)
    assert gim["uncertainties"] == pytest.approx([float(value) for value in direct])
    assert len(gim["uncertainties"]) == 3
    godambe = dadi_review.godambe_matrices(
        data, boots, model="two_epoch", p0=[1.5, 0.5], grid_points=grids
    )
    matrix, hessian, variability, _ = dadi.Godambe.get_godambe(
        model, grids, spectra, [1.5, 0.5], observed, 0.01
    )
    assert np.allclose(godambe["godambe_information"], matrix, rtol=1e-9)
    assert np.allclose(godambe["hessian"], hessian, rtol=1e-9)
    assert np.allclose(godambe["variability"], variability, rtol=1e-9)


def test_nested_model_statistics_match_upstream() -> None:
    data, boots = _three_epoch_data()
    model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.three_epoch)
    grids = [20, 30, 40]
    observed = dadi.Spectrum(np.array(data))
    spectra = [dadi.Spectrum(np.array(item)) for item in boots]
    p0 = [1.0, 0.5, 0.2, 0.3]
    lrt = dadi_review.nested_model_statistic(
        data, boots, test="lrt", model="three_epoch", p0=p0, nested_indices=[0], grid_points=grids
    )
    direct_lrt = dadi.Godambe.LRT_adjust(model, grids, spectra, p0, observed, [0])
    assert lrt["adjustment"] == pytest.approx(float(direct_lrt))
    score = dadi_review.nested_model_statistic(
        data,
        boots,
        test="score",
        model="three_epoch",
        p0=p0,
        nested_indices=[0],
        grid_points=grids,
        adjusted_and_original=True,
    )
    adjusted, original = dadi.Godambe.score_stat(
        model, grids, spectra, p0, observed, [0], adj_and_org=True
    )
    assert score["statistic"] == pytest.approx(float(adjusted))
    assert score["unadjusted"] == pytest.approx(float(original))
    wald = dadi_review.nested_model_statistic(
        data,
        boots,
        test="wald",
        model="three_epoch",
        p0=p0,
        nested_indices=[0],
        grid_points=grids,
        full_params=[2.0, 0.5, 0.2, 0.3],
    )
    direct_wald = dadi.Godambe.Wald_stat(
        model, grids, spectra, p0, observed, [0], [2.0, 0.5, 0.2, 0.3]
    )
    assert wald["statistic"] == pytest.approx(float(direct_wald))


def test_nested_model_statistic_rejects_bad_options() -> None:
    data, boots = _three_epoch_data()
    p0 = [1.0, 0.5, 0.2, 0.3]
    with pytest.raises(ConfigurationError):
        dadi_review.nested_model_statistic(
            data, boots, test="lrt", model="three_epoch", p0=p0, nested_indices=[0, 1, 2, 3]
        )
    with pytest.raises(ConfigurationError):
        dadi_review.nested_model_statistic(
            data, boots, test="wald", model="three_epoch", p0=p0, nested_indices=[0]
        )
    with pytest.raises(ConfigurationError):
        dadi_review.nested_model_statistic(
            data,
            boots,
            test="lrt",
            model="three_epoch",
            p0=p0,
            nested_indices=[0],
            adjusted_and_original=True,
        )


def test_fit_to_residual_matches_upstream() -> None:
    data, _ = _data_and_boots(num_boot=1)
    target = [value * 0.5 for value in data]
    result = dadi_review.fit_demographic_model_to_residual(
        data, target, model="two_epoch", p0=[0.8, 0.3], grid_points=(20, 30, 40), max_iterations=2
    )
    model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
    xopt, fopt, _, _, _, _, _ = dadi.Inference.optimize_log_resid(
        [0.8, 0.3],
        dadi.Spectrum(np.array(data)),
        model,
        dadi.Spectrum(np.array(target)),
        [20, 30, 40],
        verbose=0,
        maxiter=2,
        full_output=True,
    )
    assert result["params"] == pytest.approx([float(value) for value in xopt])
    assert result["objective"] == pytest.approx(float(fopt))
    assert len(result["params"]) == 2


def test_bootstrap_and_data_shape_errors() -> None:
    data, boots = _data_and_boots(num_boot=2)
    with pytest.raises(ConfigurationError):
        dadi_review.parameter_uncertainty_gim(
            data, [boots[0][:-1]], model="two_epoch", p0=[1.5, 0.5]
        )
    with pytest.raises(ConfigurationError):
        dadi_review.parameter_uncertainty_gim(data, boots, model="not_a_model", p0=[1.5, 0.5])
