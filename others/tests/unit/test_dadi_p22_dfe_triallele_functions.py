"""Tests for dadi_p22_dfe_triallele_functions: DFE-integrated spectra and the triallelic functions.

Results are compared with direct calls of the dadi classes on the same inputs. Tiny caches keep the
tests fast; a socket guard blocks the network.
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import socket
import sys

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import dadi_p22_dfe_triallele_functions as module

requires_pkg = pytest.mark.skipif(
    importlib.util.find_spec("dadi") is None, reason="dadi is not installed"
)

CACHE = {"grid_points": [16, 20, 24], "gamma_points": 5, "gamma_max": 50.0, "verbose": False}


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)
    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 1)


def _filled(rows):
    import numpy as np

    return np.array([[np.nan if v is None else v for v in row] for row in rows], dtype=float)


@requires_pkg
def test_dfe_negative_1d_matches_cache1d():
    import dadi
    import numpy as np
    from dadi.DFE import Cache1D_mod, DemogSelModels, PDFs

    out = module.dfe_integrated_spectrum(
        "negative_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0], theta=2.0, **CACHE
    )
    cache = Cache1D_mod.Cache1D(
        (2.0, 0.1),
        [6],
        DemogSelModels.two_epoch_sel,
        [16, 20, 24],
        gamma_bounds=(1e-4, 50.0),
        gamma_pts=5,
        cpus=1,
    )
    expected = cache.integrate([0.2, 100.0], None, PDFs.gamma, 2.0, None)
    assert out["kind"] == "negative_1d" and out["distribution"] == "gamma"
    got = np.array([np.nan if v is None else v for v in out["spectrum"]])
    assert np.allclose(got[1:-1], np.ma.getdata(expected)[1:-1])
    assert out["provenance"]["package"].startswith("dadi ")
    assert dadi.__name__ == "dadi"


@requires_pkg
def test_dfe_point_pos_1d_positive_mass_changes_spectrum():
    base = module.dfe_integrated_spectrum(
        "point_pos_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0, 0.0, 5.0], **CACHE
    )
    pos = module.dfe_integrated_spectrum(
        "point_pos_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0, 0.2, 5.0], **CACHE
    )
    neg = module.dfe_integrated_spectrum(
        "negative_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0], **CACHE
    )
    assert base["spectrum"][1:-1] == pytest.approx(neg["spectrum"][1:-1])
    assert pos["spectrum"][3] != pytest.approx(base["spectrum"][3])


@requires_pkg
def test_dfe_two_dimensional_kinds_agree_for_symmetric_parameters():
    args = ("split_mig_sel", [1, 1, 0.1, 1], [4, 4])
    point = module.dfe_integrated_spectrum(
        "bivariate_point_pos_2d", *args, [-1, 1, 0.5, 0.1, 3, 0.1, 3], **CACHE
    )
    sym = module.dfe_integrated_spectrum(
        "bivariate_symmetric_point_pos_2d", *args, [-1, 1, 0.5, 0.1, 3], **CACHE
    )
    plain = module.dfe_integrated_spectrum(
        "bivariate_2d", *args, [-1, 1, 0.5], exterior_int=False, **CACHE
    )
    assert _filled(point["spectrum"])[1:-1, 1:-1] == pytest.approx(
        _filled(sym["spectrum"])[1:-1, 1:-1]
    )
    assert len(plain["spectrum"]) == 5 and len(plain["spectrum"][0]) == 5
    assert plain["distribution"] == "biv_lognormal"


@requires_pkg
def test_dfe_rejects_bad_input():
    with pytest.raises(ConfigurationError):
        module.dfe_integrated_spectrum("nope", "equil", [], [4], [1.0, 1.0], verbose=False)
    with pytest.raises(ConfigurationError):
        module.dfe_integrated_spectrum(
            "negative_1d", "split_mig_sel", [1], [4], [1.0], verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.dfe_integrated_spectrum("point_pos_1d", "equil", [], [4], [0.5], verbose=False)
    with pytest.raises(ConfigurationError):
        module.dfe_integrated_spectrum(
            "negative_1d", "equil", [], [4], [1.0, 1.0], distribution="biv_lognormal", verbose=False
        )


@requires_pkg
def test_triallelic_neutral_equilibrium_matches_dadi_with_numpy_float_patched(monkeypatch):
    import numpy as np
    from dadi.Triallele import demographics

    out = module.triallelic_neutral_equilibrium(30, 8, verbose=False)
    monkeypatch.setattr(np, "float", float, raising=False)
    upstream = demographics.equilibrium([], 8, 30)
    got = _filled(out["spectrum"])
    assert np.allclose(got, np.ma.filled(upstream, np.nan), equal_nan=True)
    assert out["sample_size"] == 8 and out["folded_major"] is False
    folded = module.triallelic_neutral_equilibrium(30, 8, folded=True, verbose=False)
    assert folded["folded_major"] is True
    with_density = module.triallelic_neutral_equilibrium(12, 5, return_density=True, verbose=False)
    assert len(with_density["phi"]) == 13
    again = module.triallelic_neutral_equilibrium(12, 5, phi=with_density["phi"], verbose=False)
    assert again["spectrum"] == with_density["spectrum"]
    with pytest.raises(ConfigurationError):
        module.triallelic_neutral_equilibrium(12, None, folded=True, verbose=False)
    with pytest.raises(ConfigurationError):
        module.triallelic_neutral_equilibrium(12, None, phi=with_density["phi"], verbose=False)


@requires_pkg
def test_alternate_mutation_spectrum_known_values():
    out = module.alternate_mutation_triallelic_spectrum(6, verbose=False)
    spectrum = out["spectrum"]
    # fs[i, j] = 2n / (n - 2) / ((n - na - 1)(n - na)(n - na + 1)), na = n - i - j (Jenkins 2014)
    assert spectrum[1][1] == pytest.approx(2 * 6 / 4 / (1 * 2 * 3))  # na = 4
    assert spectrum[0][1] is None and spectrum[5][5] is None
    assert sum(v for row in spectrum for v in row if v is not None) == pytest.approx(1.0)
    with pytest.raises(ConfigurationError):
        module.alternate_mutation_triallelic_spectrum(2, verbose=False)


@requires_pkg
def test_triallelic_transform_and_segregating_sites():
    import numpy as np

    spectrum = module.triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
    folded = module.triallelic_spectrum_transform(spectrum, "fold_major", verbose=False)
    assert folded["folded_major"] is True
    unfolded = module.triallelic_spectrum_transform(
        folded["spectrum"], "unfold", folded_major=True, verbose=False
    )
    assert unfolded["folded_major"] is False
    ancestral = module.triallelic_spectrum_transform(spectrum, "fold_ancestral", verbose=False)
    assert ancestral["folded_ancestral"] is True
    logged = module.triallelic_spectrum_transform(spectrum, "log", verbose=False)
    first = next(v for v in spectrum[1] if v is not None)
    assert next(v for v in logged["spectrum"][1] if v is not None) == pytest.approx(
        float(np.log(first))
    )
    masked = module.triallelic_spectrum_transform(
        [[1.0] * 5 for _ in range(5)], "mask_infeasible", verbose=False
    )
    assert masked["spectrum"][0][0] is None and masked["spectrum"][4][4] is None
    assert masked["spectrum"][1][1] == 1.0
    with pytest.raises(ConfigurationError):
        module.triallelic_spectrum_transform(spectrum, "unfold", verbose=False)
    with pytest.raises(ConfigurationError):
        module.triallelic_spectrum_transform(spectrum, "square", verbose=False)
    total = module.triallelic_segregating_sites(spectrum, verbose=False)["segregating_sites"]
    assert total == pytest.approx(sum(v for row in spectrum for v in row if v is not None))
    infeasible = module.triallelic_segregating_sites(
        [[1.0] * 4 for _ in range(4)], apply_infeasible_mask=True, verbose=False
    )
    assert infeasible["segregating_sites"] == 1.0  # only (1, 1) has i + j < n = 3


def test_missing_dadi_raises_backend_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError):
        module.triallelic_neutral_equilibrium(12, 5, verbose=False)
    with pytest.raises(BackendUnavailableError):
        module.dfe_integrated_spectrum("negative_1d", "equil", [], [4], [1.0, 1.0], verbose=False)


def test_public_functions_have_verbose_flag():
    import inspect

    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
