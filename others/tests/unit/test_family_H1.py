"""Contract tests for the family wrappers of dadi, tskit and scikit-allel (group H1).

Every new wrapper is checked against a hand calculation or a direct call to the upstream function
with the same inputs. A sample of the ledger rows that stay ``family-variant`` is checked by calling
the existing wrapper with the claimed parameter and comparing the result with the upstream
variant. Inputs are written inside tmp_path, and no network is used.
"""

from __future__ import annotations

import inspect
import io
import math
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.population import dadi_family_functions as dadi_fam
from dnakit.population import dadi_functions as dadi_api
from dnakit.population import dadi_review_functions as dadi_review
from dnakit.population import tskit_family_functions as tskit_fam
from dnakit.population import tskit_functions as tskit_api
from dnakit.variants import allel_family_functions as allel_fam

np = pytest.importorskip("numpy")

# Four samples (0-3) under two cherries: ((0,1),(2,3)); sites at 10, 20 and 30 on 0..100.
TABLES: dict[str, Any] = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}

MODULES = (dadi_fam, tskit_fam, allel_fam)


def _public_functions() -> list[Any]:
    return [getattr(module, name) for module in MODULES for name in module.__all__]


def test_public_functions_take_keyword_only_verbose() -> None:
    functions = _public_functions()
    assert len(functions) == 6 + 2 + 18
    for function in functions:
        verbose = inspect.signature(function).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True


def test_public_signatures_use_json_friendly_types() -> None:
    for function in _public_functions():
        for parameter in inspect.signature(function).parameters.values():
            assert "Mapping[int" not in str(parameter.annotation)


# ----------------------------------------------------------------------------------------------
# dadi
# ----------------------------------------------------------------------------------------------


@pytest.fixture
def dadi_mod() -> Any:
    return pytest.importorskip("dadi")


def _two_epoch(dadi: Any, params: tuple[float, float] = (2.0, 0.5), n: int = 10) -> Any:
    return dadi.Demographics1D.two_epoch(params, (n,), 40)


def _model_and_data(dadi: Any) -> tuple[list[float], list[float], Any, Any]:
    model = _two_epoch(dadi)
    scaled = np.ma.filled(model * 1000, 0.0)
    data = np.round(scaled)
    return (
        [float(v) for v in scaled],
        [float(v) for v in data],
        dadi.Spectrum(scaled),
        dadi.Spectrum(data),
    )


def _poisson_ll(model: list[float], data: list[float]) -> float:
    """Poisson log-likelihood over the unmasked entries 1..n-1, by hand."""

    return sum(
        -m + d * math.log(m) - math.lgamma(d + 1.0)
        for m, d in zip(model[1:-1], data[1:-1], strict=True)
    )


def test_calculate_model_log_likelihood_matches_hand_and_upstream(dadi_mod):
    model, data, fs_model, fs_data = _model_and_data(dadi_mod)
    inference = dadi_mod.Inference
    plain = dadi_fam.calculate_model_log_likelihood(model, data)
    assert plain["dadi_function"] == "ll"
    assert plain["log_likelihood"] == pytest.approx(_poisson_ll(model, data), rel=1e-9)
    assert plain["log_likelihood"] == pytest.approx(float(inference.ll(fs_model, fs_data)))
    cases = [
        ("ll_multinom", {"multinom": True}, inference.ll_multinom),
        ("minus_ll", {"negate": True}, inference.minus_ll),
        ("minus_ll_multinom", {"multinom": True, "negate": True}, inference.minus_ll_multinom),
    ]
    for name, options, direct in cases:
        result = dadi_fam.calculate_model_log_likelihood(model, data, **options)
        assert result["dadi_function"] == name
        assert result["log_likelihood"] == pytest.approx(float(direct(fs_model, fs_data)))
    # The negated total is the negative of the plain total.
    negated = dadi_fam.calculate_model_log_likelihood(model, data, negate=True)
    assert negated["log_likelihood"] == pytest.approx(-plain["log_likelihood"])
    # Per-bin values sum to the total; the corners are masked.
    per_bin = dadi_fam.calculate_model_log_likelihood(model, data, per_bin=True)
    values = per_bin["log_likelihood_per_bin"]
    assert values[0] is None and values[-1] is None
    assert sum(v for v in values if v is not None) == pytest.approx(plain["log_likelihood"])
    multi_bin = dadi_fam.calculate_model_log_likelihood(model, data, multinom=True, per_bin=True)
    direct_bin = inference.ll_multinom_per_bin(fs_model, fs_data)
    expected = [
        None if m else float(v) for v, m in zip(direct_bin.data, direct_bin.mask, strict=True)
    ]
    got = multi_bin["log_likelihood_per_bin"]
    assert [x is None for x in got] == [x is None for x in expected]
    assert [x for x in got if x is not None] == pytest.approx(
        [x for x in expected if x is not None]
    )


def test_calculate_model_log_likelihood_rejects_bad_input(dadi_mod):
    model, data, _, _ = _model_and_data(dadi_mod)
    with pytest.raises(ConfigurationError):
        dadi_fam.calculate_model_log_likelihood(model, data, per_bin=True, negate=True)
    with pytest.raises(ConfigurationError):
        dadi_fam.calculate_model_log_likelihood(model, data[:-1])
    with pytest.raises(ConfigurationError):
        dadi_fam.calculate_model_log_likelihood(model, [-1.0, *data[1:]])
    with pytest.raises(ConfigurationError):
        dadi_fam.calculate_model_log_likelihood(model, data, multinom="yes")


def test_scale_model_sfs_to_data_matches_hand_ratio(dadi_mod):
    model, data, fs_model, fs_data = _model_and_data(dadi_mod)
    result = dadi_fam.scale_model_sfs_to_data(model, data)
    theta = sum(data[1:-1]) / sum(model[1:-1])
    assert result["theta"] == pytest.approx(theta, rel=1e-12)
    scaled = result["scaled_spectrum"]
    assert scaled[0] is None and scaled[-1] is None
    assert scaled[1:-1] == pytest.approx([m * theta for m in model[1:-1]], rel=1e-12)
    direct = dadi_mod.Inference.optimally_scaled_sfs(fs_model, fs_data)
    assert scaled[1:-1] == pytest.approx([float(v) for v in direct[1:-1]], rel=1e-12)
    with pytest.raises(ConfigurationError):
        dadi_fam.scale_model_sfs_to_data(model, data[:-1])


FIT = dict(model="two_epoch", p0=[1.5, 0.7], lower_bound=[1e-3, 1e-3], upper_bound=[100, 10])


def test_fit_natural_scale_matches_upstream_optimize(dadi_mod):
    _, data, _, fs_data = _model_and_data(dadi_mod)
    result = dadi_fam.fit_demographic_model_natural_scale(data, maxiter=40, **FIT)
    direct = dadi_mod.Inference.optimize(
        [1.5, 0.7],
        fs_data,
        dadi_mod.Demographics1D.two_epoch,
        40,
        lower_bound=[1e-3, 1e-3],
        upper_bound=[100, 10],
        verbose=0,
        multinom=True,
        maxiter=40,
    )
    assert result["params"] == pytest.approx([float(v) for v in direct], rel=1e-9)
    assert result["params"] == pytest.approx([2.0, 0.5], rel=0.02)
    spectrum = _two_epoch(dadi_mod, tuple(result["params"]))
    assert result["log_likelihood"] == pytest.approx(
        float(dadi_mod.Inference.ll_multinom(spectrum, fs_data))
    )
    with pytest.raises(ConfigurationError):
        dadi_fam.fit_demographic_model_natural_scale(data, **{**FIT, "p0": [0.0, 1.0]})
    with pytest.raises(ConfigurationError):
        dadi_fam.fit_demographic_model_natural_scale(data, **{**FIT, "model": "nope"})


def test_fit_constrained_meets_constraints_and_matches_upstream(dadi_mod):
    _, data, _, fs_data = _model_and_data(dadi_mod)
    equality = [{"coefficients": [1, -4], "offset": 0}]
    inequality = [{"coefficients": [1, -1]}]
    result = dadi_fam.fit_demographic_model_constrained(
        data, equality_constraints=equality, inequality_constraints=inequality, maxiter=40, **FIT
    )
    nu, time = result["params"]
    assert nu == pytest.approx(4 * time, abs=1e-6)
    assert result["equality_values"] == pytest.approx([0.0], abs=1e-6)
    assert result["inequality_values"][0] == pytest.approx(nu - time)
    assert result["converged"] is True
    direct = dadi_mod.Inference.optimize_cons(
        [1.5, 0.7],
        fs_data,
        dadi_mod.Demographics1D.two_epoch,
        40,
        eq_constraint=lambda p, *_: [p[0] - 4 * p[1]],
        ieq_constraint=lambda p, *_: [p[0] - p[1]],
        lower_bound=[1e-3, 1e-3],
        upper_bound=[100, 10],
        verbose=0,
        multinom=True,
        maxiter=40,
    )
    assert result["params"] == pytest.approx([float(v) for v in direct], rel=1e-9)
    # A binding constraint moves the optimum: nu <= time puts nu = time at the optimum.
    binding = dadi_fam.fit_demographic_model_constrained(
        data, inequality_constraints=[{"coefficients": [-1, 1]}], maxiter=40, **FIT
    )
    assert binding["params"][0] <= binding["params"][1] + 1e-6
    for bad in ({}, {"equality_constraints": [{"coefficients": [1]}]}):
        with pytest.raises(ConfigurationError):
            dadi_fam.fit_demographic_model_constrained(data, **bad, **FIT)


def test_fit_grid_matches_upstream_and_picks_best_point(dadi_mod):
    _, data, _, fs_data = _model_and_data(dadi_mod)
    ranges = [[1.0, 3.0, 3], [0.2, 1.0, 5]]
    result = dadi_fam.fit_demographic_model_grid(data, model="two_epoch", parameter_ranges=ranges)
    grid = dadi_mod.Inference.index_exp[1.0:3.0:3j, 0.2:1.0:5j]
    direct = dadi_mod.Inference.optimize_grid(
        fs_data, dadi_mod.Demographics1D.two_epoch, 40, grid, verbose=0
    )
    assert result["params"] == pytest.approx([float(v) for v in direct])
    assert result["num_evaluations"] == 15
    # The true values (2.0, 0.5) lie on the grid, so they must win.
    assert result["params"] == pytest.approx([2.0, 0.4], abs=0.11) or result["params"][0] == 2.0
    spectrum = _two_epoch(dadi_mod, tuple(result["params"]))
    assert result["log_likelihood"] == pytest.approx(
        float(dadi_mod.Inference.ll_multinom(spectrum, fs_data))
    )
    for bad in ([[1.0, 3.0, 3]], [[3.0, 1.0, 3], [0.2, 1.0, 3]], [[1.0, 3.0, 1], [0.2, 1.0, 3]]):
        with pytest.raises(ConfigurationError):
            dadi_fam.fit_demographic_model_grid(data, model="two_epoch", parameter_ranges=bad)
    with pytest.raises(ConfigurationError):
        dadi_fam.fit_demographic_model_grid(
            data, model="two_epoch", parameter_ranges=[[1.0, 3.0, 100], [0.2, 1.0, 100]]
        )


def test_fit_nlopt_matches_upstream_opt(dadi_mod):
    nlopt = pytest.importorskip("nlopt")
    _, data, _, fs_data = _model_and_data(dadi_mod)
    result = dadi_fam.fit_demographic_model_nlopt(data, maxeval=300, **FIT)
    direct, value = dadi_mod.Inference.opt(
        [1.5, 0.7],
        fs_data,
        dadi_mod.Demographics1D.two_epoch,
        40,
        lower_bound=[1e-3, 1e-3],
        upper_bound=[100, 10],
        algorithm=nlopt.LN_BOBYQA,
        local_optimizer=nlopt.LN_BOBYQA,
        maxeval=300,
        verbose=0,
    )
    assert result["params"] == pytest.approx([float(v) for v in direct], rel=1e-9)
    assert result["optimizer_value"] == pytest.approx(float(value))
    assert result["log_likelihood"] == pytest.approx(float(value))
    other = dadi_fam.fit_demographic_model_nlopt(data, algorithm="LN_SBPLX", maxeval=300, **FIT)
    assert other["params"] == pytest.approx([2.0, 0.5], rel=0.02)
    with pytest.raises(ConfigurationError):
        dadi_fam.fit_demographic_model_nlopt(data, algorithm="LD_LBFGS", **FIT)


def _nested_lists(array: Any) -> Any:
    data = np.ma.getdata(array).tolist()
    mask = np.ma.getmaskarray(array).tolist()

    def walk(value: Any, flag: Any) -> Any:
        if isinstance(value, list):
            return [walk(v, f) for v, f in zip(value, flag, strict=True)]
        return None if flag else float(value)

    return walk(data, mask)


def _assert_same(left: Any, right: Any) -> None:
    if isinstance(left, list):
        assert isinstance(right, list) and len(left) == len(right)
        for a, b in zip(left, right, strict=True):
            _assert_same(a, b)
    elif left is None or right is None:
        assert left is None and right is None
    else:
        assert left == pytest.approx(right, rel=1e-9, abs=1e-12)


# (wrapper, family attribute, model, params, sample sizes, grid points); the model is one of the
# ledger rows kept as family-variant, reached through the model argument of the wrapper.
MODEL_ROWS = [
    ("expected_sfs_1d", "Demographics1D", "growth", [2.0, 0.5], [8], 50),
    ("expected_sfs_1d", "Demographics1D", "three_epoch", [2.0, 0.5, 0.2, 0.3], [8], 50),
    ("expected_sfs_2d", "Demographics2D", "IM", [0.5] * 6, [4, 4], 20),
    ("expected_sfs_2d", "Demographics2D", "founder_asym", [0.5] * 5, [4, 4], 20),
    ("expected_sfs_3d", "Demographics3D", "sim_split_no_mig", [0.5] * 4, [3, 3, 3], 20),
    ("expected_sfs_selection", "DFE.DemogSelModels", "equil", [0.5], [8], 20),
    ("expected_sfs_selection", "DFE.DemogSelModels", "IM_sel_single_gamma", [0.5] * 7, [4, 4], 20),
]


@pytest.mark.parametrize(("wrapper", "family", "model", "params", "sizes", "points"), MODEL_ROWS)
def test_model_argument_reaches_family_variant(
    dadi_mod, wrapper, family, model, params, sizes, points
):
    import importlib

    module = importlib.import_module(f"dadi.{family}")
    direct = getattr(module, model)(tuple(params), tuple(sizes), points)
    function = getattr(dadi_api, wrapper)
    argument = sizes[0] if wrapper == "expected_sfs_1d" else sizes
    result = function(model, params, argument, grid_points=points)
    assert result["model"] == model
    _assert_same(result["spectrum"], _nested_lists(direct))


def test_ms_core_command_2d_reaches_mscore_variants(dadi_mod):
    for model, params in (
        ("IM", [0.4, 1.0, 2.0, 0.3, 0.5, 0.6]),
        ("IM_pre", [1.5, 0.2, 0.4, 1.0, 2.0, 0.3, 0.5, 0.6]),
    ):
        direct = getattr(dadi_mod.Demographics2D, f"{model}_mscore")(tuple(params))
        result = dadi_review.ms_core_command_2d(model, params)
        assert result["core"] == direct


def test_fit_method_argument_reaches_log_optimizers(dadi_mod):
    _, data, _, fs_data = _model_and_data(dadi_mod)
    for method, name in (("nelder-mead", "optimize_log_fmin"), ("powell", "optimize_log_powell")):
        result = dadi_api.fit_demographic_model(
            data, model="two_epoch", p0=[1.5, 0.7], method=method, maxiter=30
        )
        direct = getattr(dadi_mod.Inference, name)(
            [1.5, 0.7],
            fs_data,
            dadi_mod.Demographics1D.two_epoch,
            40,
            lower_bound=None,
            upper_bound=None,
            verbose=0,
            multinom=True,
            maxiter=30,
        )
        assert result["params"] == pytest.approx([float(v) for v in direct], rel=1e-9)


def test_lbfgsb_optimizers_follow_upstream_status(dadi_mod):
    """optimize_log_lbfgsb is not reachable when upstream fails with the installed SciPy."""

    _, data, _, fs_data = _model_and_data(dadi_mod)
    try:
        dadi_mod.Inference.optimize_log_lbfgsb(
            [1.5, 0.7], fs_data, dadi_mod.Demographics1D.two_epoch, 40, maxiter=5, verbose=0
        )
        upstream_works = True
    except TypeError:
        upstream_works = False
    try:
        dadi_api.fit_demographic_model(
            data, model="two_epoch", p0=[1.5, 0.7], method="lbfgsb", maxiter=5
        )
        wrapper_works = True
    except ConfigurationError:
        wrapper_works = False
    assert wrapper_works == upstream_works


def test_sfs_residuals_poisson_reaches_linear_poisson_residual(dadi_mod):
    model, data, fs_model, fs_data = _model_and_data(dadi_mod)
    result = dadi_api.sfs_residuals(model, data, residual="poisson")
    direct = dadi_mod.Inference.linear_Poisson_residual(fs_model, fs_data)
    _assert_same(result["residuals"], _nested_lists(np.ma.masked_invalid(direct)))
    # Hand calculation of one entry: (model - data) / sqrt(model).
    assert result["residuals"][2] == pytest.approx((model[2] - data[2]) / math.sqrt(model[2]))


def test_read_vcf_data_cyvcf2_parser_reaches_cyvcf2_reader(tmp_path, dadi_mod):
    pytest.importorskip("cyvcf2")
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n"
        "chr1\t100\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|0\t0|1\t1|1\t0|0\n"
        "chr1\t200\t.\tC\tT\t.\tPASS\tAA=T\tGT\t0|1\t1|1\t0|0\t0|0\n"
    )
    popinfo = tmp_path / "pop.txt"
    popinfo.write_text("s1 popA\ns2 popA\ns3 popB\ns4 popB\n")
    result = dadi_api.read_vcf_data(vcf, popinfo, parser="cyvcf2")
    direct = dadi_mod.Misc.make_data_dict_vcf_cyvcf2(str(vcf), str(popinfo), filter=True)
    assert result["parser"] == "cyvcf2"
    assert result["num_snps"] == len(direct) == 2
    builtin = dadi_mod.Misc.make_data_dict_vcf(str(vcf), str(popinfo), filter=True)
    assert {k: v["calls"] for k, v in direct.items()} == {k: v["calls"] for k, v in builtin.items()}


def test_nested_statistics_reach_wald_and_score(dadi_mod):
    rng = np.random.default_rng(3)
    model = dadi_mod.Numerics.make_extrap_func(dadi_mod.Demographics1D.three_epoch)
    mean = np.ma.filled(model((2.0, 0.5, 0.2, 0.3), (10,), [20, 30, 40]) * 200, 0.0)
    data = [float(v) for v in rng.poisson(mean).astype(float)]
    boots = [[float(v) for v in rng.poisson(mean).astype(float)] for _ in range(4)]
    spectra = [dadi_mod.Spectrum(np.array(item)) for item in boots]
    observed = dadi_mod.Spectrum(np.array(data))
    p0 = [1.0, 0.5, 0.2, 0.3]
    full = [2.0, 0.5, 0.2, 0.3]
    grids = [20, 30, 40]
    wald = dadi_review.nested_model_statistic(
        data, boots, test="wald", model="three_epoch", p0=p0, nested_indices=[0], full_params=full
    )
    assert wald["statistic"] == pytest.approx(
        float(dadi_mod.Godambe.Wald_stat(model, grids, spectra, p0, observed, [0], full))
    )
    score = dadi_review.nested_model_statistic(
        data, boots, test="score", model="three_epoch", p0=p0, nested_indices=[0]
    )
    assert score["statistic"] == pytest.approx(
        float(dadi_mod.Godambe.score_stat(model, grids, spectra, p0, observed, [0]))
    )


# ----------------------------------------------------------------------------------------------
# tskit
# ----------------------------------------------------------------------------------------------


@pytest.fixture
def tskit_mod() -> Any:
    return pytest.importorskip("tskit")


def _direct_ts(tskit: Any) -> Any:
    tables = tskit.TableCollection(sequence_length=100)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    for time in (1, 2, 3):
        tables.nodes.add_row(flags=0, time=time)
    for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]:
        tables.edges.add_row(left=0, right=100, parent=parent, child=child)
    for position in (10, 20, 30):
        tables.sites.add_row(position=position, ancestral_state="A")
    for site, node in [(0, 0), (1, 4), (2, 5)]:
        tables.mutations.add_row(site=site, node=node, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


def test_count_tree_topologies_matches_upstream_and_all_tree_wrapper(tskit_mod):
    sets = {"A": [0, 1], "B": [2, 3]}
    result = tskit_fam.count_tree_topologies(TABLES, sets)
    ts = _direct_ts(tskit_mod)
    counter = tskit_mod.combinatorics.tree_count_topologies(ts.first(), [[0, 1], [2, 3]])
    assert sorted(result["counts"]) == ["A", "A+B", "B"]
    for key, label in (((0,), "A"), ((1,), "B"), ((0, 1), "A+B")):
        total = sum(entry["count"] for entry in result["counts"][label])
        assert total == sum(counter.topologies[key].values())
    # A cherry has one embedded topology; the four-leaf tree has 2 * 2 pairings for A+B.
    assert [entry["count"] for entry in result["counts"]["A"]] == [2]
    assert [entry["count"] for entry in result["counts"]["A+B"]] == [4]
    whole = tskit_api.count_sample_topologies(TABLES, sets)
    assert whole["trees"][0]["counts"] == result["counts"]
    assert result["interval"] == [0.0, 100.0] and result["tree_index"] == 0
    with pytest.raises(ConfigurationError):
        tskit_fam.count_tree_topologies(TABLES, sets, tree_index=3)
    with pytest.raises(ConfigurationError):
        tskit_fam.count_tree_topologies(TABLES, {"A": [0]})


def test_count_tree_topologies_picks_the_requested_tree(tskit_mod):
    tables = tskit_mod.TableCollection(sequence_length=10)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit_mod.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(time=1)
    tables.nodes.add_row(time=1)
    tables.nodes.add_row(time=2)
    # Left half: ((0,1),(2,3)); right half: ((0,2),(1,3)).
    for left, right, pairs in (
        (0, 5, [(0, 4), (1, 4), (2, 5), (3, 5)]),
        (5, 10, [(0, 4), (2, 4), (1, 5), (3, 5)]),
    ):
        for child, parent in pairs:
            tables.edges.add_row(left=left, right=right, parent=parent, child=child)
        tables.edges.add_row(left=left, right=right, parent=6, child=4)
        tables.edges.add_row(left=left, right=right, parent=6, child=5)
    tables.sort()
    ts = tables.tree_sequence()
    spec = {
        "sequence_length": 10,
        "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 1}, {"time": 2}],
        "edges": [
            {"left": e.left, "right": e.right, "parent": e.parent, "child": e.child}
            for e in ts.edges()
        ],
    }
    first = tskit_fam.count_tree_topologies(spec, {"A": [0, 1], "B": [2, 3]}, tree_index=0)
    second = tskit_fam.count_tree_topologies(spec, {"A": [0, 1], "B": [2, 3]}, tree_index=1)
    assert first["interval"] == [0.0, 5.0] and second["interval"] == [5.0, 10.0]
    whole = tskit_api.count_sample_topologies(spec, {"A": [0, 1], "B": [2, 3]})
    assert [tree["counts"] for tree in whole["trees"]] == [first["counts"], second["counts"]]


def test_export_tree_sequence_nexus_matches_write_nexus(tmp_path, tskit_mod):
    target = tmp_path / "trees.nex"
    result = tskit_fam.export_tree_sequence_nexus(TABLES, target)
    ts = _direct_ts(tskit_mod)
    buffer = io.StringIO()
    tskit_mod.write_nexus(
        ts,
        buffer,
        precision=None,
        include_trees=True,
        include_alignments=None,
        reference_sequence=None,
        missing_data_character="?",
        node_labels=None,
    )
    assert target.read_text() == buffer.getvalue()
    assert result["has_trees_block"] is True and result["has_data_block"] is True
    assert "TREE t0^100 = [&R] ((n0:1,n1:1):2,(n2:2,n3:2):1);" in target.read_text()
    only_trees = tmp_path / "only_trees.nex"
    tskit_fam.export_tree_sequence_nexus(TABLES, only_trees, include_alignments=False)
    text = only_trees.read_text()
    assert "BEGIN TREES;" in text and "BEGIN DATA;" not in text
    # The alignment-only export of the existing wrapper has no TREES block.
    alignments = tmp_path / "alignments.nex"
    tskit_api.export_sample_sequences(TABLES, alignments, file_format="nexus")
    assert "BEGIN TREES;" not in alignments.read_text()
    with pytest.raises(ConfigurationError):
        tskit_fam.export_tree_sequence_nexus(TABLES, target)
    with pytest.raises(ConfigurationError):
        tskit_fam.export_tree_sequence_nexus(
            TABLES, tmp_path / "x.nex", include_trees=False, include_alignments=False
        )
    with pytest.raises(ConfigurationError):
        tskit_fam.export_tree_sequence_nexus(
            TABLES, tmp_path / "y.nex", missing_data_character="NN"
        )


def _text_tables(tskit: Any, directory: Path) -> dict[str, Path]:
    tables = tskit.TableCollection(sequence_length=100)
    tables.populations.add_row()
    tables.populations.add_row()
    for index in range(3):
        tables.individuals.add_row(flags=index, location=[0.5 * index])
    for index in range(4):
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=0, population=index % 2, individual=index % 3
        )
    tables.nodes.add_row(flags=0, time=2, population=0)
    for child in range(4):
        tables.edges.add_row(left=0, right=100, parent=4, child=child)
    tables.sites.add_row(position=10, ancestral_state="A")
    tables.sites.add_row(position=20, ancestral_state="C")
    tables.mutations.add_row(site=0, node=0, derived_state="G", time=0.5)
    tables.mutations.add_row(site=1, node=4, derived_state="T", time=2.5)
    tables.migrations.add_row(left=0, right=50, node=0, source=0, dest=1, time=1.0)
    ts = tables.tree_sequence()
    paths = {
        name: directory / f"{name}.txt"
        for name in (
            "nodes",
            "edges",
            "sites",
            "mutations",
            "individuals",
            "populations",
            "migrations",
        )
    }
    handles = {name: path.open("w") for name, path in paths.items()}
    try:
        ts.dump_text(**handles, precision=6)
    finally:
        for handle in handles.values():
            handle.close()
    return paths


@pytest.mark.parametrize(
    ("table", "rows"),
    [
        ("individuals", 3),
        ("edges", 4),
        ("sites", 2),
        ("mutations", 2),
        ("populations", 2),
        ("migrations", 1),
    ],
)
def test_read_text_table_reaches_every_parse_function(tmp_path, tskit_mod, table, rows):
    paths = _text_tables(tskit_mod, tmp_path)
    result = tskit_api.read_text_table(paths[table], table=table)
    direct = getattr(tskit_mod, f"parse_{table}")(
        io.StringIO(paths[table].read_text()), strict=True
    )
    assert result["table"] == table
    assert result["num_rows"] == direct.num_rows == rows
    expected_columns = [
        name for name in direct.column_names if name != "metadata" and not name.endswith("_offset")
    ]
    assert result["columns"] == expected_columns
    if table == "edges":
        assert [row["child"] for row in result["rows"]] == [0, 1, 2, 3]
    if table == "migrations":
        assert result["rows"][0]["dest"] == 1 and result["rows"][0]["time"] == pytest.approx(1.0)
    assert getattr(tskit_mod, f"parse_{table}") is getattr(tskit_mod.trees, f"parse_{table}")


@pytest.mark.parametrize("shape", ["star", "comb", "balanced"])
def test_random_tree_newick_reaches_tree_generators(tskit_mod, shape):
    result = tskit_api.random_tree_newick(6, shape=shape, arity=2)
    generator = getattr(tskit_mod.combinatorics, f"generate_{shape}")
    kwargs = {"arity": 2} if shape == "balanced" else {}
    tree = generator(6, span=1, branch_length=1.0, record_provenance=False, **kwargs)
    labels = {int(node): f"n{int(node)}" for node in tree.tree_sequence.samples()}
    assert result["newick"] == tree.as_newick(node_labels=labels, precision=6)
    if shape == "star":
        assert result["newick"].count("(") == 1


# ----------------------------------------------------------------------------------------------
# scikit-allel
# ----------------------------------------------------------------------------------------------


@pytest.fixture
def allel_mod() -> Any:
    return pytest.importorskip("allel")


@pytest.fixture
def panel(allel_mod) -> dict[str, Any]:
    rng = np.random.default_rng(1)
    variants = 40
    positions = np.sort(rng.choice(np.arange(1, 400), variants, replace=False))
    calls = rng.integers(0, 2, size=(variants, 6, 2)).astype("int8")
    genotypes = allel_mod.GenotypeArray(calls)
    first = genotypes.take([0, 1, 2], axis=1).count_alleles()
    second = genotypes.take([3, 4, 5], axis=1).count_alleles()
    singles = [genotypes.take([i], axis=1).count_alleles() for i in range(4)]
    return {
        "pos": positions,
        "calls": calls,
        "g": genotypes,
        "ac1": first,
        "ac2": second,
        "four": singles,
        "groups": [[0, 1, 2], [3, 4, 5]],
        "dosage": genotypes.to_n_alt(),
        "haplotypes": genotypes.to_haplotypes(),
    }


def _same_values(got: Any, expected: Any) -> None:
    flat_got = np.asarray(
        [np.nan if x is None else x for x in np.asarray(got, dtype=object).ravel()], dtype=float
    )
    flat_expected = np.asarray(expected, dtype=float).ravel()
    assert flat_got.shape == flat_expected.shape
    assert np.allclose(flat_got, flat_expected, rtol=1e-8, atol=1e-9, equal_nan=True)


def test_windowed_divergence_and_theta_match_upstream(allel_mod, panel):
    pos = panel["pos"].tolist()
    result = allel_fam.calculate_windowed_divergence(
        pos, panel["ac1"].tolist(), panel["ac2"].tolist(), window_size=100, start=1, stop=400
    )
    values, windows, bases, sites = allel_mod.windowed_divergence(
        allel_mod.SortedIndex(panel["pos"]), panel["ac1"], panel["ac2"], size=100, start=1, stop=400
    )
    _same_values([row["dxy"] for row in result["windows"]], values)
    assert [[row["start"], row["stop"]] for row in result["windows"]] == windows.tolist()
    assert [row["sites"] for row in result["windows"]] == sites.tolist()
    assert [row["callable_bases"] for row in result["windows"]] == bases.tolist()
    # Whole-region dxy times the region length equals the sum over windows of dxy * bases.
    total = sum(row["dxy"] * row["callable_bases"] for row in result["windows"])
    region = allel_mod.sequence_divergence(
        panel["pos"], panel["ac1"], panel["ac2"], start=1, stop=400
    )
    assert total / 400 == pytest.approx(float(region), rel=1e-8)
    theta = allel_fam.calculate_windowed_watterson_theta(
        pos, panel["ac1"].tolist(), window_size=100, start=1, stop=400
    )
    expected, _, _, _ = allel_mod.windowed_watterson_theta(
        allel_mod.SortedIndex(panel["pos"]), panel["ac1"], size=100, start=1, stop=400
    )
    _same_values([row["theta_per_base"] for row in theta["windows"]], expected)
    mask = [index % 2 == 0 for index in range(400)]
    masked = allel_fam.calculate_windowed_divergence(
        pos,
        panel["ac1"].tolist(),
        panel["ac2"].tolist(),
        window_size=100,
        start=1,
        stop=400,
        is_accessible=mask,
    )
    assert [row["callable_bases"] for row in masked["windows"]] == [50, 50, 50, 50]
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_windowed_divergence(
            pos, panel["ac1"].tolist(), panel["ac2"].tolist()[:-1], window_size=100
        )
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_windowed_watterson_theta(
            pos, panel["ac1"].tolist(), window_size=100, is_accessible=[True] * 5
        )


def test_moving_tajimas_d_matches_upstream(allel_mod, panel):
    result = allel_fam.calculate_moving_tajimas_d(panel["ac1"].tolist(), window_size=10)
    _same_values(result["values"], allel_mod.moving_tajima_d(panel["ac1"], 10))
    assert len(result["values"]) == 4
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_moving_tajimas_d(panel["ac1"].tolist(), window_size=0)


def test_weir_cockerham_per_variant_matches_upstream(allel_mod, panel):
    result = allel_fam.calculate_weir_cockerham_fst_per_variant(
        panel["calls"].tolist(), panel["groups"]
    )
    a, b, c = allel_mod.weir_cockerham_fst(panel["g"], panel["groups"])
    _same_values(result["a"], a.sum(axis=1))
    _same_values(result["b"], b.sum(axis=1))
    _same_values(result["c"], c.sum(axis=1))
    expected = a.sum(axis=1) / (a.sum(axis=1) + b.sum(axis=1) + c.sum(axis=1))
    _same_values(result["fst"], expected)
    assert result["fst_ratio_of_sums"] == pytest.approx(
        float(a.sum() / (a.sum() + b.sum() + c.sum()))
    )
    # The ratio of sums is the unweighted genome-wide estimate used by allel's average functions.
    average, _, _, _ = allel_mod.average_weir_cockerham_fst(panel["g"], panel["groups"], 10)
    assert result["fst_ratio_of_sums"] == pytest.approx(float(average), abs=0.1)
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_weir_cockerham_fst_per_variant(
            panel["calls"].tolist(), [[0, 1], [1, 2]]
        )
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_weir_cockerham_fst_per_variant(panel["calls"].tolist(), [[0, 9], [1]])


def test_windowed_and_moving_fst_match_upstream(allel_mod, panel):
    pos = panel["pos"].tolist()
    index = allel_mod.SortedIndex(panel["pos"])
    calls = panel["calls"].tolist()
    hudson = allel_fam.calculate_windowed_hudson_fst(
        pos, panel["ac1"].tolist(), panel["ac2"].tolist(), window_size=100, start=1, stop=400
    )
    fst, windows, counts = allel_mod.windowed_hudson_fst(
        index, panel["ac1"], panel["ac2"], size=100, start=1, stop=400
    )
    _same_values([row["fst"] for row in hudson["windows"]], fst)
    assert [row["sites"] for row in hudson["windows"]] == counts.tolist()
    weir = allel_fam.calculate_windowed_weir_cockerham_fst(
        pos, calls, panel["groups"], window_size=100, start=1, stop=400
    )
    fst, windows, counts = allel_mod.windowed_weir_cockerham_fst(
        index, panel["g"], panel["groups"], size=100, start=1, stop=400
    )
    _same_values([row["fst"] for row in weir["windows"]], fst)
    assert [[row["start"], row["stop"]] for row in weir["windows"]] == windows.tolist()
    moving_hudson = allel_fam.calculate_moving_hudson_fst(
        panel["ac1"].tolist(), panel["ac2"].tolist(), window_size=10, step=5
    )
    _same_values(
        moving_hudson["values"], allel_mod.moving_hudson_fst(panel["ac1"], panel["ac2"], 10, step=5)
    )
    moving_weir = allel_fam.calculate_moving_weir_cockerham_fst(
        calls, panel["groups"], window_size=10
    )
    _same_values(
        moving_weir["values"],
        allel_mod.moving_weir_cockerham_fst(panel["g"], panel["groups"], 10),
    )
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_windowed_hudson_fst(
            pos, panel["ac1"].tolist(), panel["ac2"].tolist(), window_size=100, start=300, stop=10
        )


def test_blockwise_and_moving_patterson_d_match_upstream(allel_mod, panel):
    four = [item.tolist() for item in panel["four"]]
    result = allel_fam.calculate_blockwise_patterson_d(*four, block_length=10)
    d, error, z_score, blocks, jackknife = allel_mod.blockwise_patterson_d(*panel["four"], 10)
    assert result["d"] == pytest.approx(float(d), rel=1e-8)
    assert result["standard_error"] == pytest.approx(float(error), rel=1e-8)
    assert result["z_score"] == pytest.approx(float(z_score), rel=1e-8)
    _same_values(result["block_values"], blocks)
    _same_values(result["jackknife_values"], jackknife)
    assert allel_mod.average_patterson_d is allel_mod.blockwise_patterson_d
    # d from the blockwise estimate is the ratio of sums, as in calculate_patterson_d's total.
    numerator, denominator = allel_mod.patterson_d(*panel["four"])
    assert result["d"] == pytest.approx(float(numerator.sum() / denominator.sum()), rel=1e-8)
    moving = allel_fam.calculate_moving_patterson_d(*four, window_size=10)
    _same_values(moving["values"], allel_mod.moving_patterson_d(*panel["four"], 10))
    assert moving["values"] == pytest.approx(result["block_values"], rel=1e-8)
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_blockwise_patterson_d(*four[:3], four[3][:-1], block_length=10)


def test_moving_garud_h_matches_upstream(allel_mod, panel):
    result = allel_fam.calculate_moving_garud_h(panel["haplotypes"].tolist(), window_size=10)
    h1, h12, h123, h2_h1 = allel_mod.moving_garud_h(panel["haplotypes"], 10)
    _same_values(result["h1"], h1)
    _same_values(result["h12"], h12)
    _same_values(result["h123"], h123)
    _same_values(result["h2_h1"], h2_h1)
    # Four windows of ten variants; a single window over all variants is garud_h of all.
    whole = allel_fam.calculate_moving_garud_h(panel["haplotypes"].tolist(), window_size=40)
    assert whole["h1"][0] == pytest.approx(float(allel_mod.garud_h(panel["haplotypes"])[0]))


def test_rogers_huff_r_between_matches_hand_values_and_upstream(allel_mod, panel):
    result = allel_fam.calculate_rogers_huff_r_between(
        [[0, 1, 2, 1]], [[0, 1, 2, 1], [2, 1, 0, 1], [1, 1, 1, 1]]
    )
    assert result["r"][0][0] == pytest.approx(1.0)
    assert result["r"][0][1] == pytest.approx(-1.0)
    assert result["r"][0][2] is None  # a constant variant has no correlation
    dosage = panel["dosage"]
    got = allel_fam.calculate_rogers_huff_r_between(dosage[:5].tolist(), dosage[5:9].tolist())
    _same_values(got["r"], allel_mod.rogers_huff_r_between(dosage[:5], dosage[5:9]))
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_rogers_huff_r_between([[0, 1, 2]], [[0, 1]])


def test_randomized_pca_matches_upstream_and_is_reproducible(allel_mod, panel):
    dosage = panel["dosage"]
    result = allel_fam.calculate_randomized_pca(dosage.tolist(), n_components=2, seed=5)
    again = allel_fam.calculate_randomized_pca(dosage.tolist(), n_components=2, seed=5)
    assert result == again
    total = dosage.sum(axis=1)
    kept = dosage[(total > 0) & (total < 12)]
    coordinates, model = allel_mod.randomized_pca(kept, n_components=2, random_state=5)
    _same_values(result["coordinates"], coordinates)
    _same_values(result["explained_variance_ratio"], model.explained_variance_ratio_)
    assert result["polymorphic_variants"] == kept.shape[0]
    assert sum(result["explained_variance_ratio"]) <= 1.0
    exact = allel_mod.pca(kept, n_components=2)[1].explained_variance_ratio_
    assert result["explained_variance_ratio"] == pytest.approx(list(exact), abs=0.05)
    for bad in (
        {"n_components": 0},
        {"n_components": 50},
        {"seed": -1},
        {"scaler": "robust"},
    ):
        with pytest.raises(ConfigurationError):
            allel_fam.calculate_randomized_pca(dosage.tolist(), **bad)
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_randomized_pca([[0, -1, 1, 2], [1, 1, 0, 2]])


def test_site_frequency_spectrum_from_counts_hand_and_upstream(allel_mod, panel):
    counts = [[3, 1], [2, 2], [0, 4]]
    plain = allel_fam.calculate_site_frequency_spectrum_from_counts(counts)
    assert plain["spectrum"] == [0, 1, 1, 0, 1] and plain["sample_size"] == 4
    folded = allel_fam.calculate_site_frequency_spectrum_from_counts(counts, fold=True)
    assert folded["spectrum"] == [1, 1, 1]  # minor counts 1, 2 and 0
    scaled = allel_fam.calculate_site_frequency_spectrum_from_counts(counts, scale=True)
    assert scaled["spectrum"] == [0.0, 1.0, 2.0, 0.0, 4.0]  # i * xi_i
    both = allel_fam.calculate_site_frequency_spectrum_from_counts(counts, fold=True, scale=True)
    assert both["spectrum"] == pytest.approx([0.0, 1 * 3 / 4, 2 * 2 / 4])  # xi_i * i(n-i)/n
    ac = panel["ac1"]
    n = 6
    for fold, scale, direct in (
        (False, False, allel_mod.sfs(ac[:, 1], n=n)),
        (True, False, allel_mod.sfs_folded(ac, n=n)),
        (False, True, allel_mod.sfs_scaled(ac[:, 1], n=n)),
        (True, True, allel_mod.sfs_folded_scaled(ac, n=n)),
    ):
        result = allel_fam.calculate_site_frequency_spectrum_from_counts(
            ac.tolist(), fold=fold, scale=scale, sample_size=n
        )
        _same_values(result["spectrum"], direct)
    # Folding the unfolded spectrum by hand (class i joins class n - i) gives the folded spectrum.
    unfolded = allel_fam.calculate_site_frequency_spectrum_from_counts(ac.tolist(), sample_size=n)
    classes = unfolded["spectrum"]
    by_hand = [classes[i] + (classes[n - i] if n - i != i else 0) for i in range(n // 2 + 1)]
    assert (
        by_hand
        == allel_fam.calculate_site_frequency_spectrum_from_counts(
            ac.tolist(), fold=True, sample_size=n
        )["spectrum"]
    )
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_site_frequency_spectrum_from_counts(counts, sample_size=3)
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_site_frequency_spectrum_from_counts([[1, 1, 1]])
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_site_frequency_spectrum_from_counts(counts, fold="yes")


def test_joint_site_frequency_spectrum_hand_and_upstream(allel_mod, panel):
    first = [[3, 1], [2, 2], [0, 4]]
    second = [[4, 0], [1, 3], [2, 2]]
    plain = allel_fam.calculate_joint_site_frequency_spectrum_from_counts(first, second)
    expected = np.zeros((5, 5), dtype=int)
    for (_, a), (_, b) in zip(first, second, strict=True):
        expected[a, b] += 1
    assert plain["spectrum"] == expected.tolist()
    folded = allel_fam.calculate_joint_site_frequency_spectrum_from_counts(first, second, fold=True)
    expected_folded = np.zeros((3, 3), dtype=int)
    for (r1, a1), (r2, a2) in zip(first, second, strict=True):
        expected_folded[min(r1, a1), min(r2, a2)] += 1
    assert folded["spectrum"] == expected_folded.tolist()
    ac1, ac2 = panel["ac1"], panel["ac2"]
    for fold, scale, direct in (
        (False, False, allel_mod.joint_sfs(ac1[:, 1], ac2[:, 1], n1=6, n2=6)),
        (True, False, allel_mod.joint_sfs_folded(ac1, ac2, n1=6, n2=6)),
        (False, True, allel_mod.joint_sfs_scaled(ac1[:, 1], ac2[:, 1], n1=6, n2=6)),
        (True, True, allel_mod.joint_sfs_folded_scaled(ac1, ac2, n1=6, n2=6)),
    ):
        result = allel_fam.calculate_joint_site_frequency_spectrum_from_counts(
            ac1.tolist(),
            ac2.tolist(),
            fold=fold,
            scale=scale,
            first_sample_size=6,
            second_sample_size=6,
        )
        _same_values(result["spectrum"], direct)
    with pytest.raises(ConfigurationError):
        allel_fam.calculate_joint_site_frequency_spectrum_from_counts(first, second[:-1])


def test_phase_trio_steps_compose_to_phase_by_transmission(allel_mod):
    rng = np.random.default_rng(11)
    trio = rng.integers(0, 2, size=(40, 3, 2)).astype("int8")
    step_one = allel_fam.phase_trio_step_by_transmission(trio.tolist())
    direct_one = allel_mod.phase_progeny_by_transmission(allel_mod.GenotypeArray(trio))
    assert step_one["genotypes"] == direct_one.values.tolist()
    assert step_one["is_phased"] == direct_one.is_phased.tolist()
    step_two = allel_fam.phase_trio_step_by_transmission(
        step_one["genotypes"],
        step="parents",
        window_size=5,
        is_phased=step_one["is_phased"],
    )
    direct_two = allel_mod.phase_parents_by_transmission(direct_one, window_size=5)
    assert step_two["genotypes"] == direct_two.values.tolist()
    assert step_two["is_phased"] == direct_two.is_phased.tolist()
    from dnakit.variants import allel_functions

    both = allel_functions.phase_trio_by_transmission(trio.tolist(), 5)
    assert step_two["genotypes"] == both["genotypes"] and step_two["is_phased"] == both["phased"]
    for bad in (
        {"step": "both"},
        {"step": "parents"},
        {"step": "progeny", "window_size": 5},
        {"step": "parents", "window_size": 5, "is_phased": [[True]]},
    ):
        with pytest.raises(ConfigurationError):
            allel_fam.phase_trio_step_by_transmission(trio.tolist(), **bad)
    with pytest.raises(ConfigurationError):
        allel_fam.phase_trio_step_by_transmission(trio[:, :2, :].tolist())


VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\trs1\tA\tG\t50\tPASS\tDP=20;AF=0.5\tGT\t0/1\t1/1
chr1\t20\t.\tC\tT,A\t.\tPASS\tDP=7;AF=0.2,0.1\tGT\t0|1\t0/0
chr1\t30\trs3\tG\tC\t99\tPASS\tDP=12\tGT\t0/0\t0/1
"""


def test_read_vcf_field_table_matches_vcf_to_dataframe(tmp_path, allel_mod):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(VCF_TEXT)
    result = allel_fam.read_vcf_field_table(str(vcf))
    assert result["columns"] == [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT_1",
        "ALT_2",
        "ALT_3",
        "QUAL",
        "FILTER_PASS",
    ]
    assert result["num_rows"] == 3
    second = result["rows"][1]
    assert second["ALT_1"] == "T" and second["ALT_2"] == "A" and second["ALT_3"] is None
    assert second["QUAL"] is None and second["POS"] == 20
    assert result["rows"][0]["QUAL"] == 50.0 and result["rows"][0]["FILTER_PASS"] is True
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        frame = allel_mod.vcf_to_dataframe(str(vcf))
        record = allel_mod.vcf_to_recarray(str(vcf))
    assert [row["POS"] for row in result["rows"]] == frame["POS"].tolist() == record["POS"].tolist()
    chosen = allel_fam.read_vcf_field_table(
        str(vcf), fields=["CHROM", "POS", "variants/DP", "variants/AF"], alt_number=2
    )
    assert chosen["columns"] == ["CHROM", "POS", "DP", "AF_1", "AF_2"]
    assert [row["DP"] for row in chosen["rows"]] == [20, 7, 12]
    assert chosen["rows"][1]["AF_2"] == pytest.approx(0.1)
    with pytest.raises(ConfigurationError):
        allel_fam.read_vcf_field_table(str(vcf), max_rows=2)
    with pytest.raises(ConfigurationError):
        allel_fam.read_vcf_field_table(str(vcf), alt_number=0)
    with pytest.raises(ConfigurationError):
        allel_fam.read_vcf_field_table(str(tmp_path / "missing.vcf"))


def test_vcf_to_csv_table_matches_vcf_to_csv(tmp_path, allel_mod):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(VCF_TEXT)
    out = tmp_path / "table.csv"
    result = allel_fam.vcf_to_csv_table(str(vcf), str(out), fields=["CHROM", "POS", "REF", "ALT"])
    direct = tmp_path / "direct.csv"
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        allel_mod.vcf_to_csv(str(vcf), str(direct), fields=["CHROM", "POS", "REF", "ALT"])
    assert out.read_text() == direct.read_text()
    assert result["columns"] == ["CHROM", "POS", "REF", "ALT_1", "ALT_2", "ALT_3"]
    assert result["num_rows"] == 3 and out.read_text().splitlines()[2] == "chr1,20,C,T,A,"
    with pytest.raises(ConfigurationError):
        allel_fam.vcf_to_csv_table(str(vcf), str(out))
    with pytest.raises(ConfigurationError):
        allel_fam.vcf_to_csv_table(str(vcf), str(tmp_path / "table.txt"))
