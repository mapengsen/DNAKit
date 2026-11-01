"""Tests for dadi_p22_twolocus_functions: two-locus spectrum operations and spectrum files.

Each result is compared with the dadi call on the same spectrum. A socket guard blocks the network.
"""

from __future__ import annotations

import importlib.util
import inspect
import socket
import sys

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import dadi_p22_twolocus_functions as module

requires_pkg = pytest.mark.skipif(
    importlib.util.find_spec("dadi") is None, reason="dadi is not installed"
)
N = 6


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _cube(seed=1, n=N):
    import numpy as np

    return np.random.default_rng(seed).random((n + 1,) * 3)


def _masked_nested(seed=1, n=N):
    return module.two_locus_spectrum_transform(
        _cube(seed, n).tolist(), "mask_infeasible", verbose=False
    )["spectrum"]


def _array(nested):
    import numpy as np

    flat = np.array(nested, dtype=float)
    return flat


@requires_pkg
def test_transform_matches_dadi_numerics():
    import numpy as np
    from dadi.TwoLocus import TLSpectrum_mod, numerics

    nested = _masked_nested()

    def fresh():  # dadi's fold edits its input's mask in place; use a copy per reference
        return TLSpectrum_mod.TLSpectrum(_cube())

    references = {
        "fold": lambda: fresh().fold(),
        "fold_ancestral": lambda: numerics.fold_ancestral(fresh()),
        "fold_lr": lambda: numerics.fold_lr(fresh()),
        "misidentification": lambda: numerics.misidentification(fresh(), 0.1),
        "project": lambda: numerics.project(fresh(), 4),
    }
    for operation, reference in references.items():
        expected = reference()
        out = module.two_locus_spectrum_transform(
            nested, operation, sample_size=4, p_misid=0.1, verbose=False
        )
        got = _array(out["spectrum"])
        assert np.allclose(got, np.ma.filled(expected, np.nan), equal_nan=True), operation
    assert module.two_locus_spectrum_transform(nested, "fold", verbose=False)["folded"] is True
    assert (
        module.two_locus_spectrum_transform(nested, "project", sample_size=4, verbose=False)[
            "sample_size"
        ]
        == 4
    )


@requires_pkg
def test_transform_masks_unfold_and_array_to_spectrum():
    import numpy as np

    raw = _cube().tolist()
    masked = module.two_locus_spectrum_transform(raw, "mask_infeasible", verbose=False)["spectrum"]
    array = _array(masked)
    assert np.isnan(array[0, 0, 0]) and np.isnan(array[6, 6, 6]) and not np.isnan(array[1, 1, 1])
    folded = module.two_locus_spectrum_transform(masked, "fold", verbose=False)
    unfolded = module.two_locus_spectrum_transform(
        folded["spectrum"], "unfold", folded=True, verbose=False
    )
    assert unfolded["folded"] is False
    as_spectrum = module.two_locus_spectrum_transform(raw, "array_to_spectrum", verbose=False)
    assert (
        _array(as_spectrum["spectrum"])[0, 0, 0] != _array(as_spectrum["spectrum"])[0, 0, 0]
    )  # NaN: masked
    with pytest.raises(ConfigurationError):
        module.two_locus_spectrum_transform(masked, "unfold", verbose=False)
    with pytest.raises(ConfigurationError):
        module.two_locus_spectrum_transform(masked, "project", verbose=False)
    with pytest.raises(ConfigurationError):
        module.two_locus_spectrum_transform(masked, "misidentification", verbose=False)
    with pytest.raises(ConfigurationError):
        module.two_locus_spectrum_transform(masked, "project", sample_size=9, verbose=False)
    with pytest.raises(ConfigurationError):
        module.two_locus_spectrum_transform([[1, 2], [3, 4]], "fold", verbose=False)


@requires_pkg
def test_marginals_and_mean_r2_match_dadi():
    import numpy as np
    from dadi.TwoLocus import TLSpectrum_mod, numerics

    nested = _masked_nested()
    fs = TLSpectrum_mod.TLSpectrum(_cube())
    marginal = module.two_locus_marginal_spectra(nested, verbose=False)
    assert np.allclose(_array(marginal["locus_a"])[1:-1], np.ma.getdata(fs.marginalA())[1:-1])
    assert np.allclose(_array(marginal["locus_b"])[1:-1], np.ma.getdata(fs.marginalB())[1:-1])
    single = module.two_locus_marginal_spectra(nested, method="single_locus", verbose=False)
    a, _b = numerics.to_single_locus(fs)
    assert np.allclose(_array(single["locus_a"])[1:-1], np.ma.getdata(a)[1:-1])
    for method in ("numerics", "spectrum_method"):
        value = module.two_locus_mean_r2(nested, method=method, verbose=False)["mean_r2"]
        assert value == pytest.approx(float(numerics.mean_r2(fs)))
    with pytest.raises(ConfigurationError):
        module.two_locus_mean_r2(nested, method="other", verbose=False)


@requires_pkg
def test_ld_per_bin_values():
    import numpy as np

    out = module.two_locus_ld_per_bin(4, verbose=False)
    d, r2 = _array(out["D"]), _array(out["r2"])
    assert d.shape == (5, 5, 5)
    # AB = 2, Ab = 0, aB = 0 of n = 4: pA = pB = 0.5, pAB = 0.5, D = 0.25, r2 = 1
    assert d[2, 0, 0] == pytest.approx(0.25) and r2[2, 0, 0] == pytest.approx(1.0)
    assert np.isnan(r2[0, 0, 0]) or r2[0, 0, 0] != r2[0, 0, 0]
    with pytest.raises(ConfigurationError):
        module.two_locus_ld_per_bin(1, verbose=False)


@requires_pkg
def test_log_likelihoods_match_dadi_inference():
    import numpy as np
    from dadi.TwoLocus import TLSpectrum_mod, inference

    model_nested = _masked_nested(1)
    data_nested = _masked_nested(2)
    model = TLSpectrum_mod.TLSpectrum(_cube(1))
    data = TLSpectrum_mod.TLSpectrum(_cube(2))
    plain = module.two_locus_log_likelihood(model_nested, data_nested, verbose=False)
    assert plain["log_likelihood"] == pytest.approx(float(inference.ll(model, data)))
    multi = module.two_locus_log_likelihood(model_nested, data_nested, multinom=True, verbose=False)
    assert multi["log_likelihood"] == pytest.approx(float(inference.ll_multinom(model, data)))
    per_bin = module.two_locus_log_likelihood(
        model_nested, data_nested, per_bin=True, verbose=False
    )
    expected = inference.ll_per_bin(model, data)
    assert np.allclose(
        _array(per_bin["log_likelihood_per_bin"]), np.ma.filled(expected, np.nan), equal_nan=True
    )
    both = module.two_locus_log_likelihood_over_rho_bins(
        [model_nested, model_nested], [data_nested, data_nested], verbose=False
    )
    assert both["log_likelihood"] == pytest.approx(2 * plain["log_likelihood"])
    both_multi = module.two_locus_log_likelihood_over_rho_bins(
        [model_nested, model_nested], [data_nested, data_nested], multinom=True, verbose=False
    )
    assert both_multi["log_likelihood"] == pytest.approx(2 * multi["log_likelihood"])
    scaling = module.scale_two_locus_spectrum(model_nested, data_nested, verbose=False)
    assert scaling["scaling"] == pytest.approx(float(inference.optimal_sfs_scaling(model, data)))
    scaled = _array(scaling["scaled_model"])
    assert np.nansum(scaled) == pytest.approx(
        float(np.ma.sum(inference.optimally_scaled_sfs(model, data)))
    )
    with pytest.raises(ConfigurationError):
        module.two_locus_log_likelihood(model_nested, _masked_nested(2, 5), verbose=False)
    with pytest.raises(ConfigurationError):
        module.two_locus_log_likelihood_over_rho_bins([model_nested], [], verbose=False)


@requires_pkg
def test_condition_and_sample_density():
    import numpy as np
    from dadi.TwoLocus import TLSpectrum_mod, numerics

    nested = _masked_nested()
    fs = TLSpectrum_mod.TLSpectrum(_cube())
    out = module.condition_two_locus_spectrum(nested, 2, verbose=False)
    assert np.allclose(_array(out["spectrum"]), numerics.condition(fs, 2))
    assert _array(out["spectrum"]).shape == (3, N - 2 + 1)
    with pytest.raises(ConfigurationError):
        module.condition_two_locus_spectrum(nested, 0, verbose=False)
    phi = np.random.default_rng(3).random((6, 6, 6))
    sampled = module.sample_two_locus_density(phi.tolist(), 4, verbose=False)
    expected = numerics.sample_cached(phi, 4, numerics.grid(5))
    assert np.allclose(_array(sampled["spectrum"]), np.ma.filled(expected, np.nan), equal_nan=True)
    with pytest.raises(ConfigurationError):
        module.sample_two_locus_density(np.ones((60, 60, 60)).tolist(), 4, verbose=False)


@requires_pkg
def test_spectrum_file_round_trip(tmp_path):
    import numpy as np

    nested = _masked_nested()
    target = tmp_path / "tl.fs"
    written = module.write_multilocus_spectrum_file(
        nested,
        str(target),
        kind="two_locus",
        extrap_x=0.01,
        extrap_t=0.005,
        comment_lines=["note"],
        verbose=False,
    )
    assert written["bytes"] == target.stat().st_size
    read = module.read_multilocus_spectrum_file(str(target), kind="two_locus", verbose=False)
    assert read["comments"] == ["note"] and read["extrap_x"] == 0.01 and read["extrap_t"] == 0.005
    assert np.allclose(_array(read["spectrum"]), _array(nested), equal_nan=True)
    with pytest.raises(ConfigurationError):
        module.write_multilocus_spectrum_file(nested, str(target), kind="two_locus", verbose=False)
    tri = [
        [None, None, None, None],
        [None, 1.5, None, None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    tri_target = tmp_path / "tri.fs"
    module.write_multilocus_spectrum_file(tri, str(tri_target), folded_major=True, verbose=False)
    tri_read = module.read_multilocus_spectrum_file(
        str(tri_target), mask_infeasible=False, verbose=False
    )
    assert tri_read["folded_major"] is True and tri_read["folded_ancestral"] is False
    assert tri_read["spectrum"][1][1] == 1.5 and tri_read["spectrum"][0][0] is None
    with pytest.raises(ConfigurationError):
        module.read_multilocus_spectrum_file(str(tmp_path / "missing.fs"), verbose=False)
    with pytest.raises(ConfigurationError):
        module.read_multilocus_spectrum_file(str(tri_target), kind="two_locus", verbose=False)


def test_missing_dadi_raises_backend_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError):
        module.two_locus_ld_per_bin(4, verbose=False)


def test_public_functions_have_verbose_flag():
    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
