"""Tests for dadi_gap_functions: spectrum, residual and DFE plots, DFE mixture spectra, Demes
inference and the depth-of-coverage distribution.

Plots are drawn with dadi on the Agg backend into tmp_path and checked for the PNG signature and
the returned data. The DFE mixture tests build dadi's selection caches on tiny grids (sample sizes
4 and 4), so the linearity identities hold exactly and the runs stay in seconds. Upstream defects
found while testing are recorded as strict expected failures and in the disposition skip reasons.
A socket guard blocks any connection, and dadi's caches run in one process (see the fixture).
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import socket
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import dadi_gap_functions as module

requires_pkg = pytest.mark.skipif(
    importlib.util.find_spec("dadi") is None,
    reason="dadi is not installed",
)


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)
    # dadi's caches fork worker processes when several CPUs are visible; one CPU keeps them in-process.
    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 1)


def _is_png(path: Path) -> bool:
    return path.is_file() and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


SFS_1D = [None, 4.0, 2.0, 1.25, 1.0, None]
SFS_2D = [[None, 3.0, 1.0], [2.0, 1.0, 0.5], [None, 0.5, None]]
SFS_3D = [
    [[None, 2.0, 1.0], [1.5, 1.0, 0.5], [0.8, 0.4, 0.2]],
    [[1.0, 0.5, 0.3], [0.6, 0.3, 0.2], [0.4, 0.2, 0.1]],
    [[0.5, 0.3, 0.1], [0.3, 0.2, 0.1], [0.2, 0.1, None]],
]
RESID_2D = [[None, 0.5, -0.2], [-0.1, 0.0, 0.3], [None, 0.2, None]]
TRI = [[None, 1.0, 2.0], [1.0, 0.5, None], [2.0, None, 0.0]]
TWO_LOCUS = [[[None, 1.0], [1.0, 0.5]], [[1.0, 0.5], [0.5, 0.2]]]


# ---------------------------------------------------------------------------
# Spectrum plots
# ---------------------------------------------------------------------------


@requires_pkg
def test_plot_sfs_1d_writes_png_and_returns_spectrum(tmp_path):
    target = tmp_path / "one.png"
    result = module.plot_sfs_1d(SFS_1D, str(target), verbose=False)
    assert _is_png(target)
    assert result["image_path"] == str(target.resolve())
    assert result["image_format"] == "png"
    assert result["sample_sizes"] == [5]
    assert result["spectrum"] == SFS_1D


@requires_pkg
def test_plot_sfs_1d_comparison_multinom_and_poisson(tmp_path):
    model = [None, 4.0, 2.0, 1.25, 1.0, None]
    data = [None, 3.5, 2.5, 1.0, 1.0, None]
    multinom = module.plot_sfs_1d_comparison(model, data, str(tmp_path / "m.png"), verbose=False)
    poisson = module.plot_sfs_1d_comparison(
        model, data, str(tmp_path / "p.png"), likelihood="poisson", residual="linear", verbose=False
    )
    assert _is_png(tmp_path / "m.png") and _is_png(tmp_path / "p.png")
    assert multinom["likelihood"] == "multinom" and multinom["residual"] == "Anscombe"
    assert poisson["likelihood"] == "poisson" and poisson["residual"] == "linear"
    assert multinom["model"] == model and multinom["data"] == data


@requires_pkg
def test_plot_sfs_2d_and_residual(tmp_path):
    result = module.plot_sfs_2d(
        SFS_2D, str(tmp_path / "2d.png"), vmin=0.25, vmax=4.0, pop_ids=["A", "B"], verbose=False
    )
    assert _is_png(tmp_path / "2d.png")
    assert result["sample_sizes"] == [2, 2]
    assert result["spectrum"] == SFS_2D
    residual = module.plot_sfs_2d_residual(
        RESID_2D, str(tmp_path / "r.png"), resid_range=0.5, pop_ids=["A", "B"], verbose=False
    )
    assert _is_png(tmp_path / "r.png")
    assert residual["residual"] == RESID_2D and residual["resid_range"] == 0.5


@requires_pkg
def test_plot_sfs_2d_comparison_and_meta_residual(tmp_path):
    comparison = module.plot_sfs_2d_comparison(
        SFS_2D, SFS_2D, str(tmp_path / "cmp.png"), likelihood="poisson", verbose=False
    )
    assert _is_png(tmp_path / "cmp.png")
    assert comparison["likelihood"] == "poisson"
    assert comparison["model"] == SFS_2D and comparison["data"] == SFS_2D
    meta = module.plot_sfs_2d_meta_residual(
        RESID_2D,
        [[None, 0.1, -0.2], [-0.1, 0.0, 0.3], [None, 0.2, None]],
        str(tmp_path / "meta.png"),
        verbose=False,
    )
    assert _is_png(tmp_path / "meta.png")
    assert meta["synonymous_residual"] == RESID_2D


@requires_pkg
def test_plot_sfs_3d_family(tmp_path):
    single = module.plot_sfs_3d(SFS_3D, str(tmp_path / "3d.png"), verbose=False)
    assert _is_png(tmp_path / "3d.png")
    assert single["sample_sizes"] == [2, 2, 2] and single["spectrum"] == SFS_3D
    comparison = module.plot_sfs_3d_comparison(
        SFS_3D, SFS_3D, str(tmp_path / "3dc.png"), residual="linear", verbose=False
    )
    assert _is_png(tmp_path / "3dc.png") and comparison["residual"] == "linear"
    pairwise = module.plot_sfs_3d_pairwise(
        SFS_3D, str(tmp_path / "pair.png"), pop_ids=["A", "B", "C"], verbose=False
    )
    assert _is_png(tmp_path / "pair.png") and pairwise["spectrum"] == SFS_3D


@requires_pkg
def test_plot_trispectrum_argument_checks_run_before_drawing(tmp_path):
    with pytest.raises(ConfigurationError, match="folded"):
        module.plot_trispectrum(TRI, str(tmp_path / "t.png"), folded="yes", verbose=False)
    with pytest.raises(ConfigurationError, match="titles must be strings"):
        module.plot_trispectrum_comparison(
            TRI, TRI, str(tmp_path / "t2.png"), title1=3, verbose=False
        )
    assert not (tmp_path / "t.png").exists()


@pytest.mark.xfail(
    strict=True,
    raises=ConfigurationError,
    reason=(
        "upstream dadi.Triallele.plotting.plot_single_trispectrum passes norm=LogNorm together with "
        "vmin/vmax; matplotlib 3.11.2 raises ValueError: Passing a Normalize instance simultaneously "
        "with vmin/vmax is not supported"
    ),
)
@requires_pkg
def test_plot_trispectrum_draws_png(tmp_path):
    result = module.plot_trispectrum(TRI, str(tmp_path / "t.png"), verbose=False)
    assert _is_png(tmp_path / "t.png") and result["folded"] is False


@pytest.mark.xfail(
    strict=True,
    raises=ConfigurationError,
    reason=(
        "upstream dadi.Triallele.plotting.plot_trispectrum_comp passes norm=LogNorm together with "
        "vmin/vmax; matplotlib 3.11.2 raises ValueError: Passing a Normalize instance simultaneously "
        "with vmin/vmax is not supported"
    ),
)
@requires_pkg
def test_plot_trispectrum_comparison_draws_png(tmp_path):
    result = module.plot_trispectrum_comparison(TRI, TRI, str(tmp_path / "tc.png"), verbose=False)
    assert _is_png(tmp_path / "tc.png") and result["folded"] is False


@requires_pkg
def test_plot_two_locus_comparison_draws_png(tmp_path):
    result = module.plot_two_locus_comparison(
        TWO_LOCUS, TWO_LOCUS, str(tmp_path / "tl.png"), resid_range=1.0, verbose=False
    )
    assert _is_png(tmp_path / "tl.png")
    assert result["resid_range"] == 1.0
    assert len(result["model"]) == len(result["data"]) == 2
    assert all(len(plane) == 2 and all(len(row) == 2 for row in plane) for plane in result["model"])


@pytest.mark.xfail(
    strict=True,
    raises=ConfigurationError,
    reason=(
        "upstream dadi.TwoLocus.plotting.plot_3d calls ax.colorbar() unconditionally; the call fails "
        "with AttributeError: 'Axes3D' object has no attribute 'colorbar'"
    ),
)
@requires_pkg
def test_plot_two_locus_spectrum_draws_png(tmp_path):
    result = module.plot_two_locus_spectrum(TWO_LOCUS, str(tmp_path / "tl1.png"), verbose=False)
    assert _is_png(tmp_path / "tl1.png") and result["spectrum"] == TWO_LOCUS


@requires_pkg
def test_two_locus_slice_arguments_are_checked_before_drawing(tmp_path):
    with pytest.raises(ConfigurationError, match="max_sum"):
        module.plot_two_locus_spectrum(
            TWO_LOCUS, str(tmp_path / "x.png"), max_sum=-1, verbose=False
        )
    assert not (tmp_path / "x.png").exists()


# ---------------------------------------------------------------------------
# DFE bivariate plots
# ---------------------------------------------------------------------------


@requires_pkg
def test_plot_dfe_bivariate_gamma_and_lognormal(tmp_path):
    gamma = module.plot_dfe_bivariate(
        "biv_ind_gamma", [0.2, 10.0], str(tmp_path / "g.png"), grid_size=8, verbose=False
    )
    assert _is_png(tmp_path / "g.png")
    assert gamma["distribution"] == "biv_ind_gamma" and gamma["params"] == [0.2, 10.0]
    grid = gamma["gamma_grid"]
    assert len(grid) == 8 and all(value < 0 for value in grid)
    # Negative gammas on a log grid, from -1e-3 (weakest) down to -100 (strongest), as dadi expects.
    assert grid[0] == pytest.approx(-1e-3) and grid[-1] == pytest.approx(-100.0)
    lognormal = module.plot_dfe_bivariate(
        "biv_lognormal", [1.0, 0.5, 0.2], str(tmp_path / "l.png"), grid_size=8, verbose=False
    )
    assert _is_png(tmp_path / "l.png")
    assert lognormal["distribution"] == "biv_lognormal" and lognormal["params"] == [1.0, 0.5, 0.2]


@requires_pkg
def test_plot_dfe_bivariate_point_pos(tmp_path):
    result = module.plot_dfe_bivariate_point_pos(
        "biv_lognormal",
        [1.0, 0.5, 0.0, 0.1, 5.0, 0.1, 3.0],
        str(tmp_path / "pp.png"),
        rho=0.2,
        grid_size=8,
        verbose=False,
    )
    assert _is_png(tmp_path / "pp.png")
    assert result["rho"] == 0.2 and result["params"] == [1.0, 0.5, 0.0, 0.1, 5.0, 0.1, 3.0]


@requires_pkg
def test_plot_dfe_bivariate_rejects_unknown_distribution(tmp_path):
    with pytest.raises(ConfigurationError, match="distribution"):
        module.plot_dfe_bivariate("biv_normal", [1.0], str(tmp_path / "x.png"), verbose=False)
    assert not (tmp_path / "x.png").exists()


# ---------------------------------------------------------------------------
# DFE mixture spectra (dadi selection caches)
# ---------------------------------------------------------------------------

DEMO = [2.0, 2.0, 0.2, 0.0]
MIX_OPTS = dict(
    sample_sizes=[4, 4], grid_points=(12, 14, 16), gamma_points=5, gamma_max=200.0, verbose=False
)


def _mixture(variant, dfe_params, theta=1.0):
    return module.dfe_mixture_spectrum(
        variant, DEMO, dfe_params=dfe_params, theta=theta, **MIX_OPTS
    )


@requires_pkg
def test_mixture_weight_interpolates_between_components():
    # The mixture is (1 - w) * lognormal spectrum + w * bivariate spectrum. The two end points are
    # computed with weight 0 and 1, so the half-weight spectrum must be their average.
    only_1d = _mixture("mixture", [1.0, 0.5, 0.2, 0.0])["spectrum"]
    only_2d = _mixture("mixture", [1.0, 0.5, 0.2, 1.0])["spectrum"]
    half = _mixture("mixture", [1.0, 0.5, 0.2, 0.5])
    assert half["variant"] == "mixture" and half["sample_sizes"] == [4, 4]
    for row_half, row_1d, row_2d in zip(half["spectrum"], only_1d, only_2d, strict=True):
        for value_half, value_1d, value_2d in zip(row_half, row_1d, row_2d, strict=True):
            if value_half is None:
                assert value_1d is None and value_2d is None
            else:
                assert value_half == pytest.approx(
                    0.5 * value_1d + 0.5 * value_2d, rel=1e-6, abs=1e-9
                )


@requires_pkg
def test_mixture_spectrum_scales_linearly_with_theta():
    one = _mixture("mixture", [1.0, 0.5, 0.2, 0.3], theta=1.0)["spectrum"]
    two = _mixture("mixture", [1.0, 0.5, 0.2, 0.3], theta=2.0)["spectrum"]
    for row_one, row_two in zip(one, two, strict=True):
        for value_one, value_two in zip(row_one, row_two, strict=True):
            if value_one is None:
                assert value_two is None
            else:
                assert value_two == pytest.approx(2.0 * value_one, rel=1e-9)


@requires_pkg
def test_symmetric_point_pos_and_vourlaki_spectra_have_expected_shape():
    symmetric = _mixture("mixture_symmetric_point_pos", [1.0, 0.5, 0.2, 0.1, 5.0, 0.3])
    assert symmetric["variant"] == "mixture_symmetric_point_pos"
    assert len(symmetric["spectrum"]) == 5 and all(len(row) == 5 for row in symmetric["spectrum"])
    assert symmetric["gamma_points"] == 5 and symmetric["grid_points"] == [12, 14, 16]
    vourlaki = _mixture("Vourlaki_mixture", [0.2, 10.0, 0.0, 10.0, 0.0, 0.0])
    assert vourlaki["variant"] == "Vourlaki_mixture"
    assert len(vourlaki["spectrum"]) == 5 and all(len(row) == 5 for row in vourlaki["spectrum"])
    values = [v for row in vourlaki["spectrum"] for v in row if v is not None]
    assert values and max(values) > 0


@pytest.mark.xfail(
    strict=True,
    raises=ConfigurationError,
    reason=(
        "upstream dadi.DFE.Cache2D_mod.mixture_point_pos calls integrate_point_pos(params2, None, "
        "sel_dist2, theta, None), which passes None as rho; the call fails with TypeError: "
        "unsupported operand type(s) for *: 'NoneType' and 'float'"
    ),
)
@requires_pkg
def test_mixture_point_pos_variant_draws_spectrum():
    result = _mixture("mixture_point_pos", [1.0, 0.5, 0.2, 0.1, 5.0, 0.1, 3.0, 0.3])
    assert result["variant"] == "mixture_point_pos"


@pytest.mark.parametrize(
    ("variant", "dfe_params", "message"),
    [
        ("mixture", [1.0, 0.5, 0.2], "needs 4"),
        ("mixture_symmetric_point_pos", [1.0, 0.5, 0.2], "needs 6"),
        ("unknown", [1.0, 0.5, 0.2, 0.3], "variant"),
    ],
)
@requires_pkg
def test_mixture_checks_variant_and_parameter_count(variant, dfe_params, message):
    with pytest.raises(ConfigurationError, match=message):
        _mixture(variant, dfe_params)


@requires_pkg
def test_mixture_checks_demography_length():
    with pytest.raises(ConfigurationError, match="four values"):
        module.dfe_mixture_spectrum(
            "mixture", [2.0, 2.0, 0.2], [4, 4], [1.0, 0.5, 0.2, 0.3], verbose=False
        )


def test_mixture_needs_dadi(monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError) as info:
        _mixture("mixture", [1.0, 0.5, 0.2, 0.3])
    assert "bio-dadi" in (info.value.hint or "")


# ---------------------------------------------------------------------------
# Demes inference
# ---------------------------------------------------------------------------


@requires_pkg
def test_demes_fit_rejects_missing_inputs_and_bad_method(tmp_path):
    graph = tmp_path / "graph.yaml"
    options = tmp_path / "options.yaml"
    graph.write_text("description: x\ntime_units: generations\ndemes: []\n", encoding="utf-8")
    options.write_text("parameters: []\n", encoding="utf-8")
    sfs = [[None, 1.0], [1.0, None]]
    with pytest.raises(ConfigurationError, match="does not exist"):
        module.fit_demes_graph_parameters(
            tmp_path / "missing.yaml", options, sfs, ["A", "B"], [8, 10, 12], verbose=False
        )
    with pytest.raises(ConfigurationError, match="method"):
        module.fit_demes_graph_parameters(
            graph, options, sfs, ["A", "B"], [8, 10, 12], method="newton", verbose=False
        )


@pytest.mark.xfail(
    strict=True,
    raises=ConfigurationError,
    reason=(
        "upstream dadi.Demes.Inference.optimize reads the options file with ruamel.yaml.load, which "
        'ruamel.yaml 0.19.1 (installed) removed; the call fails with AttributeError: "load()" has '
        "been removed"
    ),
)
@requires_pkg
def test_demes_fit_returns_parameters_in_bounds(tmp_path):
    graph = tmp_path / "model.yaml"
    options = tmp_path / "options.yaml"
    graph.write_text(
        "description: split\ntime_units: generations\ndemes:\n"
        "  - name: anc\n    epochs:\n      - start_size: 1000\n        end_time: 500\n"
        "  - name: A\n    ancestors: [anc]\n    epochs:\n      - start_size: 1000\n"
        "  - name: B\n    ancestors: [anc]\n    epochs:\n      - start_size: 2000\n",
        encoding="utf-8",
    )
    options.write_text(
        "parameters:\n  - name: nuB\n    values:\n      - demes:\n          B:\n"
        "            epochs:\n              0: start_size\n    lower_bound: 1\n    upper_bound: 10000\n",
        encoding="utf-8",
    )
    sfs = [
        [None, 0.5, 0.2, 0.1, None],
        [0.6, 0.9, 0.4, 0.2, 0.1],
        [0.2, 0.3, 0.5, 0.2, 0.1],
        [0.1, 0.2, 0.1, 0.3, 0.1],
        [None, 0.1, 0.1, 0.1, None],
    ]
    result = module.fit_demes_graph_parameters(
        graph, options, sfs, ["A", "B"], [8, 10, 12], maxiter=3, verbose=False
    )
    assert result["parameter_names"] == ["nuB"]
    assert 1.0 <= result["params"][0] <= 10000.0


# ---------------------------------------------------------------------------
# Depth-of-coverage distribution
# ---------------------------------------------------------------------------


@requires_pkg
def test_coverage_depth_distribution_known_answer():
    # Population A depths over both sites: 1, 2, 2, 3 -> frequencies over depths 0..3.
    # Population B depths over both sites: 0, 1 -> frequencies 0.5 at depths 0 and 1.
    result = module.coverage_depth_distribution(
        {"A": [[1, 2], [2, 3]], "B": [[0], [1]]}, verbose=False
    )
    assert result["sites"] == 2
    assert result["distribution"]["A"]["depth"] == [0, 1, 2, 3]
    assert result["distribution"]["A"]["frequency"] == pytest.approx([0.0, 0.25, 0.5, 0.25])
    assert result["distribution"]["B"]["depth"] == [0, 1]
    assert result["distribution"]["B"]["frequency"] == pytest.approx([0.5, 0.5])


@pytest.mark.parametrize(
    ("site_depths", "message"),
    [
        ({"A": [[1], [2]], "B": [[1]]}, "same, non-zero number of sites"),
        ({"A": [[-1]]}, "depth"),
        ({}, "at least one population"),
        ({"A": [[]]}, "at least one depth"),
    ],
)
@requires_pkg
def test_coverage_depth_distribution_rejects_bad_input(site_depths, message):
    with pytest.raises(ConfigurationError, match=message):
        module.coverage_depth_distribution(site_depths, verbose=False)


# ---------------------------------------------------------------------------
# Output and dependency rules
# ---------------------------------------------------------------------------


@requires_pkg
def test_plot_refuses_existing_output_and_non_png(tmp_path):
    target = tmp_path / "exists.png"
    target.write_bytes(b"x")
    with pytest.raises(ConfigurationError, match="already exists"):
        module.plot_sfs_1d(SFS_1D, str(target), verbose=False)
    with pytest.raises(ConfigurationError, match=r"\.png"):
        module.plot_sfs_1d(SFS_1D, str(tmp_path / "plot.pdf"), verbose=False)
    assert target.read_bytes() == b"x"


@requires_pkg
def test_plot_checks_shape_and_pop_ids(tmp_path):
    with pytest.raises(ConfigurationError, match="axes"):
        module.plot_sfs_2d(SFS_1D, str(tmp_path / "a.png"), verbose=False)
    with pytest.raises(ConfigurationError, match="pop_ids"):
        module.plot_sfs_2d(SFS_2D, str(tmp_path / "b.png"), pop_ids=["A"], verbose=False)


def test_plot_needs_dadi(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError) as info:
        module.plot_sfs_1d(SFS_1D, str(tmp_path / "x.png"), verbose=False)
    assert "bio-dadi" in (info.value.hint or "")
    assert not (tmp_path / "x.png").exists()
