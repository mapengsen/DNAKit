"""Tests for dadi_verify_functions: triallelic plots, the two-locus 3D plot, the DFE point-mass mixture
and Demes inference, each run where the upstream call fails.

Each implementation is compared with dadi's own arithmetic on the same input: the fold with dadi's
``_fold``, the 3D voxel faces with the faces dadi's ``plot_3d`` builds (its failing ``colorbar`` call
is stubbed only inside the test), the DFE mixture with dadi's integration methods called with the
rho argument they need, and the Demes inference with dadi's ``optimize`` whose YAML loader is replaced
by the current ruamel.yaml API. The upstream defects themselves are asserted as well. The masks are
taken from dadi's own spectrum objects, which normalise them. Plots are PNG files in tmp_path. Socket
connections are blocked; dadi's caches run in one process.
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import socket

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.population import dadi_verify_functions as module

requires_dadi = pytest.mark.skipif(
    importlib.util.find_spec("dadi") is None or importlib.util.find_spec("demes") is None,
    reason="dadi and demes are not installed",
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)
    # dadi's caches fork worker processes when several CPUs are visible; one CPU keeps them in-process.
    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 1)


def _triallelic(size=5):
    """An unfolded triallelic spectrum: entries with i + j <= size - 1 are valid, the rest are masked."""
    return [[float(3 + 2 * i + j) if i + j <= size - 1 else None for j in range(size)] for i in range(size)]


def _nested(array):
    """A masked array as the nested list the public functions take (None for masked entries)."""
    import numpy as np

    mask = np.ma.getmaskarray(array)
    data = np.ma.getdata(array)
    if data.ndim == 1:
        return [None if mask[i] else float(data[i]) for i in range(data.shape[0])]
    return [_nested(array_slice) for array_slice in [np.ma.masked_array(data[i], mask=mask[i]) for i in range(data.shape[0])]]


def _png(path):
    return path.read_bytes()[:8] == PNG_SIGNATURE


@requires_dadi
def test_triallelic_plot_writes_png_and_drawn_values(tmp_path):
    import numpy as np

    spectrum = _triallelic(5)
    target = tmp_path / "tri.png"
    result = module.plot_triallelic_spectrum(spectrum, str(target), colorbar=True, verbose=False)
    assert _png(target)
    assert result["image_path"] == str(target)
    positive = [v for row in spectrum for v in row if v is not None and v > 0]
    assert result["vmin"] == min(positive)
    assert result["vmax"] == max(positive)
    drawn = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]], dtype=float)
    expected = np.array([[np.nan if v is None else v for v in row] for row in spectrum], dtype=float)
    assert np.allclose(drawn, np.where(expected >= result["vmin"], expected, np.nan), equal_nan=True)


@requires_dadi
def test_fold_matches_dadi_fold():
    import dadi
    import numpy as np
    from dadi.Triallele import plotting as upstream

    spectrum = dadi.Spectrum(np.ma.masked_invalid(np.array([[np.nan if v is None else v for v in row] for row in _triallelic(6)])))
    assert not np.ma.getmaskarray(spectrum)[1, 2]
    upstream_folded = upstream._fold(spectrum.copy())
    same_input = np.ma.array(np.ma.getdata(spectrum).copy(), mask=np.ma.getmaskarray(spectrum).copy())
    ours = module._fold(np, same_input)
    assert np.array_equal(np.ma.getmaskarray(ours), np.ma.getmaskarray(upstream_folded))
    keep = ~np.ma.getmaskarray(ours)
    assert np.allclose(np.ma.getdata(ours)[keep], np.ma.getdata(upstream_folded)[keep])


@requires_dadi
def test_folded_plot_rejects_an_already_folded_spectrum(tmp_path):
    spectrum = _triallelic(5)
    spectrum[1][2] = None
    with pytest.raises(ConfigurationError, match="already folded"):
        module.plot_triallelic_spectrum(spectrum, str(tmp_path / "x.png"), folded=True, verbose=False)


@requires_dadi
def test_upstream_single_trispectrum_fails_in_matplotlib():
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    from dadi.Triallele import plotting as upstream
    from dadi.Triallele.TriSpectrum_mod import TriSpectrum

    raw = np.ma.masked_invalid(np.array([[np.nan if v is None else v for v in row] for row in _triallelic(5)]))
    spectrum = TriSpectrum(raw)
    with pytest.raises(ValueError, match="Normalize"):
        upstream.plot_single_trispectrum(spectrum, folded=False)


@requires_dadi
def test_comparison_residual_matches_the_formula(tmp_path):
    import numpy as np

    first = _triallelic(5)
    second = [[None if v is None else v * 1.25 for v in row] for row in first]
    target = tmp_path / "comp.png"
    result = module.plot_triallelic_spectrum_comparison(first, second, str(target), verbose=False)
    assert _png(target)
    residual = np.array([[np.nan if v is None else v for v in row] for row in result["residual"]], dtype=float)
    a = np.array([[np.nan if v is None else v for v in row] for row in first], dtype=float)
    b = np.array([[np.nan if v is None else v for v in row] for row in second], dtype=float)
    expected = (a - b) / np.sqrt(a)
    assert np.allclose(residual, expected, equal_nan=True)


@requires_dadi
def test_two_locus_faces_match_upstream_plot_3d(tmp_path, monkeypatch):
    import numpy as np
    from mpl_toolkits.mplot3d import art3d, axes3d
    from dadi.TwoLocus import plotting as upstream
    from dadi.TwoLocus.TLSpectrum_mod import TLSpectrum

    raw = np.ma.masked_invalid(
        np.array(
            [[[np.nan if i + j + k >= 6 else float(1 + i + 2 * j + 3 * k) for k in range(4)] for j in range(4)] for i in range(4)]
        )
    )
    fs = TLSpectrum(raw.copy())
    values = _nested(fs)
    recorded = {"ours": [], "upstream": []}
    original = art3d.Poly3DCollection
    current = {"name": "ours"}

    class Recorder(original):
        def __init__(self, polys, *args, **kwargs):
            recorded[current["name"]].append(np.array(polys, dtype=float))
            super().__init__(polys, *args, **kwargs)

    monkeypatch.setattr(art3d, "Poly3DCollection", Recorder)
    monkeypatch.setattr(axes3d.Axes3D, "colorbar", lambda self, *a, **k: None, raising=False)

    result = module.plot_two_locus_spectrum_3d(values, str(tmp_path / "tl.png"), verbose=False)
    assert _png(tmp_path / "tl.png")
    current["name"] = "upstream"
    upstream.plot_3d(fs, vmin=result["vmin"], vmax=result["vmax"], show=False)
    assert len(recorded["ours"]) == len(recorded["upstream"]) == 1
    assert np.array_equal(recorded["ours"][0], recorded["upstream"][0])
    assert result["drawn_entries"] > 0


@requires_dadi
def test_upstream_plot_3d_fails_on_axes3d_colorbar():
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    from dadi.TwoLocus import plotting as upstream
    from dadi.TwoLocus.TLSpectrum_mod import TLSpectrum

    raw = np.ma.masked_invalid(np.array([[[float(1 + i + j + k) for k in range(4)] for j in range(4)] for i in range(4)]))
    fs = TLSpectrum(raw)
    with pytest.raises(AttributeError, match="colorbar"):
        upstream.plot_3d(fs, show=False)


@requires_dadi
def test_dfe_mixture_point_pos_matches_dadi_components():
    import numpy as np
    from dadi.DFE import Cache1D_mod, Cache2D_mod, DemogSelModels, PDFs

    demo = (1.2, 0.8, 0.3, 0.5)
    values = [-1.1, 0.6, 0.4, 0.2, 0.5, 0.3, 0.2, 0.6]
    pdf_params, (rho, ppos1, gamma1, ppos2, gamma2, p2d) = values[:-6], values[-6:]
    ns, pts, theta = [4, 4], [10, 12, 14], 1.0
    kwargs = dict(gamma_bounds=(1e-4, 2000.0), gamma_pts=5, additional_gammas=[gamma1, gamma2], verbose=False)
    s1 = Cache1D_mod.Cache1D(demo, ns, DemogSelModels.split_mig_sel_single_gamma, pts, **kwargs)
    s2 = Cache2D_mod.Cache2D(demo, ns, DemogSelModels.split_mig_sel, pts, **kwargs)
    fs1 = s1.integrate_point_pos(list(pdf_params) + [ppos1, gamma1], None, PDFs.lognormal, theta, Npos=1, pts=None)
    fs2 = s2.integrate_point_pos(
        list(pdf_params) + [rho, ppos1, gamma1, ppos2, gamma2], None, PDFs.biv_lognormal, theta, rho, None
    )
    expected = (1 - p2d) * fs1 + p2d * fs2
    result = module.dfe_mixture_point_pos_spectrum(
        list(demo), ns, values, theta, grid_points=pts, gamma_points=5, gamma_max=2000.0, verbose=False
    )
    got = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]], dtype=float)
    expected_mask = np.ma.getmaskarray(expected)
    assert np.array_equal(np.isnan(got), expected_mask)
    keep = ~expected_mask
    assert np.allclose(got[keep], np.ma.getdata(expected)[keep], rtol=1e-8, atol=1e-12)


@requires_dadi
def test_upstream_mixture_point_pos_fails_on_none_rho():
    from dadi.DFE import Cache1D_mod, Cache2D_mod, DemogSelModels, PDFs
    from dadi.DFE import Cache2D_mod as mixtures

    demo = (1.2, 0.8, 0.3, 0.5)
    values = [-1.1, 0.6, 0.4, 0.2, 0.5, 0.3, 0.2, 0.6]
    ns, pts = [4, 4], [10, 12, 14]
    kwargs = dict(gamma_bounds=(1e-4, 2000.0), gamma_pts=5, additional_gammas=[0.5, 0.2], verbose=False)
    s1 = Cache1D_mod.Cache1D(demo, ns, DemogSelModels.split_mig_sel_single_gamma, pts, **kwargs)
    s2 = Cache2D_mod.Cache2D(demo, ns, DemogSelModels.split_mig_sel, pts, **kwargs)
    with pytest.raises(TypeError):
        mixtures.mixture_point_pos(values, None, s1, s2, PDFs.lognormal, PDFs.biv_lognormal, 1.0, None)


DEMES_GRAPH = """description: split with two demes
time_units: generations
demes:
  - name: anc
    epochs:
      - start_size: 1000
        end_time: 500
  - name: A
    ancestors: [anc]
    epochs:
      - start_size: 800
        end_time: 0
  - name: B
    ancestors: [anc]
    epochs:
      - start_size: 600
        end_time: 0
"""

DEMES_OPTIONS = """parameters:
  - name: nuA
    values:
      - demes:
          A:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
  - name: nuB
    values:
      - demes:
          B:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
"""


def _repaired_options_loader(fname):
    import ruamel.yaml

    with open(fname, "r") as handle:
        options = ruamel.yaml.YAML(typ="safe").load(handle)
    if "parameters" not in options.keys():
        raise ValueError("parameters to fit must be specified")
    return options


def _sfs_nested(spectrum):
    import numpy as np

    mask = np.ma.getmaskarray(spectrum)
    return [[None if mask[i, j] else float(spectrum[i, j]) for j in range(spectrum.shape[1])] for i in range(spectrum.shape[0])]


def _demes_inputs(tmp_path):
    import dadi

    graph = tmp_path / "graph.yaml"
    graph.write_text(DEMES_GRAPH)
    options = tmp_path / "options.yaml"
    options.write_text(DEMES_OPTIONS)
    data = dadi.Demographics2D.split_mig([1.1, 0.7, 0.4, 0.6], [4, 4], 10)
    data.pop_ids = ["A", "B"]
    return graph, options, data


@requires_dadi
def test_demes_optimize_matches_dadi_with_repaired_loader(tmp_path, monkeypatch):
    import numpy as np
    from dadi.Demes import Inference

    graph, options, data = _demes_inputs(tmp_path)
    result = module.optimize_demes_graph_parameters(
        graph, options, _sfs_nested(data), ["A", "B"], [10, 12, 14], maxiter=3, verbose=False
    )
    monkeypatch.setattr(Inference, "_get_params_dict", _repaired_options_loader)
    with open("/dev/null", "w") as sink:
        names, xopt, fopt = Inference.optimize(
            str(graph), str(options), data, [10, 12, 14], maxiter=3, verbose=0, output_stream=sink
        )
    assert result["parameter_names"] == names
    assert np.allclose(result["parameters"], xopt, rtol=1e-9)
    assert abs(result["objective"] - fopt) < 1e-9


@requires_dadi
def test_upstream_demes_optimize_fails_on_removed_yaml_load(tmp_path):
    from dadi.Demes import Inference

    graph, options, data = _demes_inputs(tmp_path)
    with open("/dev/null", "w") as sink:
        with pytest.raises(AttributeError, match="load"):
            Inference.optimize(str(graph), str(options), data, [10, 12, 14], maxiter=1, verbose=0, output_stream=sink)


@requires_dadi
def test_demes_output_is_a_new_file(tmp_path):
    import demes

    graph, options, data = _demes_inputs(tmp_path)
    out = tmp_path / "fitted.yaml"
    result = module.optimize_demes_graph_parameters(
        graph, options, _sfs_nested(data), ["A", "B"], [10, 12, 14], maxiter=2, output_path=out, verbose=False
    )
    assert result["output_path"] == str(out)
    assert demes.load(str(out)).demes[1].name == "A"
    with pytest.raises(ConfigurationError):
        module.optimize_demes_graph_parameters(
            graph, options, _sfs_nested(data), ["A", "B"], [10, 12, 14], maxiter=2, output_path=out, verbose=False
        )
