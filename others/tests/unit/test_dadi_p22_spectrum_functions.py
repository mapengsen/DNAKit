"""Tests for dadi_p22_spectrum_functions: spectrum operations, statistics, low coverage, diffusion.

Each result is compared with the dadi call on the same input. Random operations are checked for
reproducibility with a seed and for leaving NumPy's global random state untouched. A socket guard
blocks the network.
"""

from __future__ import annotations

import importlib.util
import inspect
import socket
import sys

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import dadi_p22_spectrum_functions as module

requires_pkg = pytest.mark.skipif(
    importlib.util.find_spec("dadi") is None, reason="dadi is not installed"
)
COVERAGE = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _snm(n=10):
    import dadi

    return dadi.Demographics1D.snm((), (n,), 40)


def _nested1d(n=10):
    values = _snm(n).tolist()
    values[0] = values[-1] = None
    return values


def _joint():
    import dadi
    import numpy as np

    fs = dadi.Demographics2D.split_mig((1, 1, 0.2, 1), (6, 6), 20)
    return fs, np.ma.getdata(fs).tolist()


@requires_pkg
@pytest.mark.parametrize(
    "statistic", ["S", "Watterson_theta", "theta_L", "Zengs_E", "pi", "Tajima_D"]
)
def test_statistics_match_dadi(statistic):
    import dadi

    got = module.calculate_sfs_statistic(_nested1d(), statistic, verbose=False)["value"]
    assert got == pytest.approx(float(getattr(dadi.Spectrum(_snm()), statistic)()))


@requires_pkg
def test_fst_matches_dadi_and_one_population_statistics_reject_joint_spectra():
    import dadi
    import numpy as np

    fs, nested = _joint()
    got = module.calculate_sfs_statistic(nested, "Fst", verbose=False)["value"]
    assert got == pytest.approx(float(dadi.Spectrum(np.ma.getdata(fs)).Fst()))
    with pytest.raises(ConfigurationError):
        module.calculate_sfs_statistic(nested, "Tajima_D", verbose=False)
    with pytest.raises(ConfigurationError):
        module.calculate_sfs_statistic(nested, "kurtosis", verbose=False)


@requires_pkg
def test_mask_project_subset_reorder_log_match_dadi():
    import dadi
    import numpy as np

    fs, nested = _joint()
    plain = dadi.Spectrum(np.ma.getdata(fs))

    def values(out):
        return np.array([[np.nan if v is None else v for v in row] for row in out["spectrum"]])

    corners = module.mask_spectrum_corners(nested, verbose=False)
    assert corners["spectrum"][0][0] is None and corners["spectrum"][6][6] is None
    projected = module.project_spectrum(nested, [4, 3], verbose=False)
    assert projected["sample_sizes"] == [4, 3]
    assert np.allclose(
        values(projected), np.ma.filled(plain.project([4, 3]), np.nan), equal_nan=True
    )
    kept = module.subset_spectrum_populations(nested, keep=[2], verbose=False)
    summed = module.subset_spectrum_populations(nested, marginalize_over=[1], verbose=False)
    expected = plain.filter_pops([2])
    assert kept["spectrum"][1:-1] == pytest.approx(list(np.ma.getdata(expected)[1:-1]))
    assert summed["spectrum"] == kept["spectrum"]
    reordered = module.reorder_spectrum_populations(nested, [2, 1], verbose=False)
    assert np.allclose(
        values(reordered), np.ma.filled(plain.reorder_pops([2, 1]), np.nan), equal_nan=True
    )
    logged = module.log_spectrum(_nested1d(), verbose=False)
    assert logged["spectrum"][2] == pytest.approx(float(np.log(_snm()[2])))
    with pytest.raises(ConfigurationError):
        module.subset_spectrum_populations(nested, verbose=False)
    with pytest.raises(ConfigurationError):
        module.subset_spectrum_populations(nested, keep=[1], marginalize_over=[2], verbose=False)
    with pytest.raises(ConfigurationError):
        module.subset_spectrum_populations(nested, keep=[3], verbose=False)
    with pytest.raises(ConfigurationError):
        module.project_spectrum(nested, [9, 3], verbose=False)
    with pytest.raises(ConfigurationError):
        module.reorder_spectrum_populations(nested, [1, 1], verbose=False)


@requires_pkg
def test_unfold_matches_dadi_and_needs_a_folded_spectrum():
    import dadi
    import numpy as np

    folded = dadi.Spectrum(_snm()).fold()
    nested = [
        None if m else float(v)
        for v, m in zip(np.ma.getdata(folded), np.ma.getmaskarray(folded), strict=True)
    ]
    out = module.unfold_spectrum(nested, verbose=False)
    expected = folded.unfold()
    assert out["folded"] is False
    assert [x for x in out["spectrum"] if x is not None] == pytest.approx(
        [float(x) for x in np.ma.compressed(expected)]
    )
    assert out["spectrum"][0] is None and out["spectrum"][-1] is None


@requires_pkg
def test_resampling_is_seeded_and_restores_global_state():
    import numpy as np

    np.random.seed(11)
    before = np.random.get_state()[1].copy()
    first = module.resample_spectrum(_nested1d(), seed=5, verbose=False)
    after = np.random.get_state()[1]
    assert (before == after).all()
    again = module.resample_spectrum(_nested1d(), seed=5, verbose=False)
    other = module.resample_spectrum(_nested1d(), seed=6, verbose=False)
    assert first["spectrum"] == again["spectrum"] and first["spectrum"] != other["spectrum"]
    assert all(float(v).is_integer() for v in first["spectrum"] if v is not None)
    fixed = module.resample_spectrum(
        _nested1d(), "fixed_size", num_sites=500, seed=2, verbose=False
    )
    assert sum(v for v in fixed["spectrum"] if v is not None) <= 500
    assert fixed["num_sites"] == 500
    nonmasked = module.resample_spectrum(
        _nested1d(), "fixed_size", num_sites=500, only_nonmasked=True, seed=2, verbose=False
    )
    assert sum(v for v in nonmasked["spectrum"] if v is not None) == 500
    with pytest.raises(ConfigurationError):
        module.resample_spectrum(_nested1d(), "fixed_size", verbose=False)
    with pytest.raises(ConfigurationError):
        module.resample_spectrum(_nested1d(), "bootstrap", verbose=False)


@requires_pkg
def test_ms_file_parsing_and_errors(tmp_path):
    import dadi

    path = tmp_path / "toy.ms"
    path.write_text(
        "ms 4 2 -t 2.0 -I 2 2 2\n123 456 789\n\n//\nsegsites: 2\npositions: 0.1 0.5\n"
        "01\n10\n01\n10\n\n"
        "//\nsegsites: 1\npositions: 0.3\n1\n0\n0\n1\n"
    )
    out = module.spectrum_from_simulation_file(str(path), verbose=False)
    expected = dadi.Spectrum.from_ms_file(str(path))
    assert out["sample_sizes"] == [2, 2] and out["command"].startswith("ms 4 2")
    assert out["spectrum"][1][1] == pytest.approx(float(expected[1, 1]))
    summed = module.spectrum_from_simulation_file(str(path), average=False, verbose=False)
    assert summed["spectrum"][1][1] == pytest.approx(2 * out["spectrum"][1][1])
    split = module.spectrum_from_simulation_file(str(path), pop_assignments=[4], verbose=False)
    assert split["sample_sizes"] == [4]
    with pytest.raises(ConfigurationError):
        module.spectrum_from_simulation_file(str(tmp_path / "none.ms"), verbose=False)
    with pytest.raises(ConfigurationError):
        module.spectrum_from_simulation_file(str(path), file_format="vcf", verbose=False)


@requires_pkg
def test_fragment_snp_data(tmp_path):
    path = tmp_path / "snps.txt"
    path.write_text(
        "Ingroup\tOutgroup\tAllele1\tA\tB\tAllele2\tA\tB\tGene\tPosition\n"
        "TCT\tTCT\tC\t4\t3\tT\t2\t2\tg1\tchr1_100\n"
        "ACA\tACA\tC\t1\t6\tA\t3\t1\tg1\tchr1_200\n"
        "GCG\tGCG\tC\t1\t6\tG\t3\t1\tg2\tchr2_50\n"
    )
    out = module.fragment_snp_data(str(path), 100000, include_snp_ids=True, verbose=False)
    assert out["num_snps"] == 3 and sum(c["num_snps"] for c in out["chunks"]) == 3
    assert {i for c in out["chunks"] for i in c["snp_ids"]} == {
        "g1_chr1_100",
        "g1_chr1_200",
        "g2_chr2_50",
    }
    with pytest.raises(ConfigurationError):
        module.fragment_snp_data(str(path), 0, verbose=False)
    with pytest.raises(ConfigurationError):
        module.fragment_snp_data(str(tmp_path / "x"), 10, verbose=False)


@requires_pkg
def test_extrapolated_spectrum_matches_dadi():
    import dadi
    import numpy as np

    for log in (False, True):
        out = module.expected_sfs_extrapolated(
            "two_epoch",
            [2.0, 0.1],
            [8],
            grid_points=[30, 40, 50],
            log_extrapolation=log,
            verbose=False,
        )
        wrap = dadi.Numerics.make_extrap_log_func if log else dadi.Numerics.make_extrap_func
        expected = wrap(dadi.Demographics1D.two_epoch)((2.0, 0.1), (8,), [30, 40, 50])
        got = np.array([np.nan if v is None else v for v in out["spectrum"]])
        assert np.allclose(got, np.ma.filled(expected, np.nan), equal_nan=True)
    with pytest.raises(ConfigurationError):
        module.expected_sfs_extrapolated(
            "two_epoch", [2.0, 0.1], [8], grid_points=[30, 40], verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.expected_sfs_extrapolated("split_mig", [1, 1, 0.1, 1], [8], verbose=False)


@requires_pkg
def test_log_likelihood_derivatives_match_godambe():
    import dadi
    import numpy as np

    data = dadi.Demographics1D.two_epoch((2.0, 0.3), (10,), 40) * 500
    nested = [float(v) for v in np.ma.getdata(data)]
    grids = [20, 30, 40]
    function = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
    observed = dadi.Spectrum(np.array(nested))

    def objective(p, spec):
        return dadi.Inference.ll_multinom(function(tuple(p), (10,), grids), spec)

    hess = module.log_likelihood_derivatives(nested, "two_epoch", [1.5, 0.2], verbose=False)
    grad = module.log_likelihood_derivatives(
        nested, "two_epoch", [1.5, 0.2], derivative="gradient", verbose=False
    )
    assert np.allclose(
        hess["value"], dadi.Godambe.get_hess(objective, [1.5, 0.2], 0.01, args=(observed,))
    )
    assert np.allclose(
        grad["value"],
        np.ravel(dadi.Godambe.get_grad(objective, [1.5, 0.2], 0.01, args=(observed,))),
    )
    poisson = module.log_likelihood_derivatives(
        nested, "two_epoch", [1.5, 0.2], multinom=False, verbose=False
    )
    assert poisson["value"] != hess["value"]
    with pytest.raises(ConfigurationError):
        module.log_likelihood_derivatives(nested, "split_mig", [1, 1], verbose=False)
    with pytest.raises(ConfigurationError):
        module.log_likelihood_derivatives(
            nested, "two_epoch", [1.5, 0.2], derivative="jacobian", verbose=False
        )


@requires_pkg
def test_integrate_diffusion_phi_one_to_five_populations():
    import dadi
    import numpy as np

    xx = dadi.Numerics.default_grid(16)
    phi = dadi.PhiManip.phi_1D(xx)
    one = module.integrate_diffusion_phi(
        1, 0.1, nu=[2.0], grid_points=16, return_phi=True, sample_sizes=[6], verbose=False
    )
    expected = dadi.Integration.one_pop(phi, xx, 0.1, nu=2.0)
    assert np.allclose(one["phi"], expected)
    assert one["spectrum"][1] == pytest.approx(
        float(dadi.Spectrum.from_phi(expected, [6], (xx,))[1])
    )
    two = module.integrate_diffusion_phi(
        2,
        0.1,
        nu=[2.0, 1.0],
        migration=[[0, 1.0], [0.5, 0]],
        gamma=[0.0, -1.0],
        grid_points=16,
        return_phi=True,
        verbose=False,
    )
    expected2 = dadi.Integration.two_pops(
        dadi.PhiManip.phi_1D_to_2D(xx, phi),
        xx,
        0.1,
        nu1=2.0,
        nu2=1.0,
        m12=1.0,
        m21=0.5,
        gamma1=0.0,
        gamma2=-1.0,
    )
    assert np.allclose(two["phi"], expected2)
    for pops, grid in ((3, 12), (4, 10), (5, 8)):
        out = module.integrate_diffusion_phi(
            pops, 0.02, grid_points=grid, sample_sizes=[4] * pops, verbose=False
        )
        assert out["phi_shape"] == [grid] * pops and len(out["spectrum"]) == 5
    x_linked = module.integrate_diffusion_phi(
        1, 0.1, grid_points=16, x_linked=True, return_phi=True, verbose=False
    )
    x_expected = dadi.Integration.one_pop_X(dadi.PhiManip.phi_1D_X(xx), xx, 0.1)
    assert np.allclose(x_linked["phi"], x_expected)
    with pytest.raises(ConfigurationError):
        module.integrate_diffusion_phi(2, 0.1, x_linked=True, verbose=False)
    with pytest.raises(ConfigurationError):
        module.integrate_diffusion_phi(3, 0.01, grid_points=40, return_phi=True, verbose=False)
    with pytest.raises(ConfigurationError):
        module.integrate_diffusion_phi(6, 0.1, verbose=False)
    with pytest.raises(ConfigurationError):
        module.integrate_diffusion_phi(2, 0.1, migration=[[0, 1]], verbose=False)


@requires_pkg
def test_low_coverage_functions():
    import numpy as np
    from dadi.LowPass import LowPass

    cov = {"A": np.array([np.arange(8), COVERAGE["A"]])}
    no_call = module.gatk_no_call_probability(COVERAGE["A"], 6, verbose=False)["probabilities"]
    assert no_call == pytest.approx(
        list(LowPass.probability_of_no_call_1D_GATK_multisample(cov["A"], 6, 0))
    )
    assert no_call[0] == 1.0
    enough = module.probability_enough_individuals_covered(COVERAGE["A"], 6, 4, verbose=False)[
        "probability"
    ]
    assert enough == pytest.approx(
        float(LowPass.probability_enough_individuals_covered(cov["A"], 6, 4))
    )

    state = np.random.get_state()[1].copy()
    reads = module.simulate_low_coverage_reads(
        COVERAGE, [0, 1, 2, 1], [8], 6, seed=3, verbose=False
    )
    assert (np.random.get_state()[1] == state).all()
    again = module.simulate_low_coverage_reads(
        COVERAGE, [0, 1, 2, 1], [8], 6, seed=3, verbose=False
    )
    assert reads == again
    assert len(reads["reference_reads"]) == 6 and len(reads["reference_reads"][0]) == 4
    assert all(r == 0 for r in [row[0] for row in reads["alternate_reads"]])  # individual 1 is 0/0
    assert all(row[2] == 0 for row in reads["reference_reads"])  # individual 3 is 1/1

    calls = [[0, 1, 2, 99], [1, 1, 0, 2], [99, 99, 99, 1]]
    sub = module.subsample_genotypes(calls, 4, seed=1, verbose=False)
    assert sub["num_loci"] == 2 and all(len(row) == 2 for row in sub["genotypes"])
    assert sub == module.subsample_genotypes(calls, 4, seed=1, verbose=False)

    spectrum = module.simulate_gatk_calling_spectrum(
        COVERAGE, [3], [6], [4], 300, seed=1, verbose=False
    )
    assert (
        sum(spectrum["probabilities"]) == pytest.approx(1.0) and len(spectrum["probabilities"]) == 5
    )
    assert spectrum == module.simulate_gatk_calling_spectrum(
        COVERAGE, [3], [6], [4], 300, seed=1, verbose=False
    )

    model = module.expected_sfs_low_pass(
        "two_epoch",
        [2.0, 0.1],
        COVERAGE,
        [6],
        [4],
        grid_points=[20, 30, 40],
        simulations_per_frequency=20,
        seed=1,
        verbose=False,
    )
    assert model["sample_sizes"] == [4] and len(model["spectrum"]) == 5
    with pytest.raises(ConfigurationError):
        module.gatk_no_call_probability(COVERAGE["A"], 5, verbose=False)
    with pytest.raises(ConfigurationError):
        module.gatk_no_call_probability([0.5, 0.1], 6, verbose=False)
    with pytest.raises(ConfigurationError):
        module.simulate_low_coverage_reads(COVERAGE, [0, 1], [8], 6, verbose=False)
    with pytest.raises(ConfigurationError):
        module.simulate_gatk_calling_spectrum(COVERAGE, [3], [6], [8], 10, verbose=False)
    with pytest.raises(ConfigurationError):
        module.subsample_genotypes([[0, 1, 5]], 2, verbose=False)
    with pytest.raises(ConfigurationError):
        module.expected_sfs_low_pass("two_epoch", [2, 0.1], COVERAGE, [6], [4, 4], verbose=False)


def test_missing_dadi_raises_backend_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError):
        module.calculate_sfs_statistic([0, 1, 1, 0], "S", verbose=False)


def test_public_functions_have_verbose_flag():
    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
