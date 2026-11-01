"""Tests for dadi_gap2_functions: the Vourlaki (2022) validation figure and the MayaVi 3D spectrum plot.

The Vourlaki figure is drawn with dadi on small caches (sample sizes 4 and 4, three small grids)
so the test runs in seconds; the figure, the spectra it returns and the parameter table are checked.
MayaVi is not installed in the test environment. Its guard path runs everywhere; the plot logic
runs against a fake ``mayavi.mlab`` module that records the calls; the real render is skipped.
A socket guard blocks any connection.
"""

from __future__ import annotations

import importlib.util
import math
import multiprocessing
import socket
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import dadi_gap2_functions as module


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


class _FakeMlab:
    """Records the mlab calls made by plot_3d_spectrum_mayavi and writes a PNG on savefig."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple, dict]] = []
        self.barchart_args: tuple = ()
        self.axes_kwargs: dict = {}
        self.axes_obj = SimpleNamespace(axes=SimpleNamespace(label_format="unset"))

    def figure(self, *args, **kwargs):
        self.calls.append(("figure", args, kwargs))
        return "figure-1"

    def clf(self, figure):
        self.calls.append(("clf", (figure,), {}))

    def barchart(self, xs, ys, zs, values, **kwargs):
        self.calls.append(("barchart", (), kwargs))
        self.barchart_args = (list(xs), list(ys), list(zs), [float(v) for v in values])

    def axes(self, **kwargs):
        self.axes_kwargs = kwargs
        return self.axes_obj

    def savefig(self, path, figure=None):
        self.calls.append(("savefig", (path,), {"figure": figure}))
        Path(path).write_bytes(b"\x89PNG\r\n\x1a\n" + b"fake")

    def close(self, figure):
        self.calls.append(("close", (figure,), {}))


# ---------------------------------------------------------------------------
# Vourlaki (2022) validation figure
# ---------------------------------------------------------------------------


@requires_pkg
def test_vourlaki_validation_writes_png_and_spectra(tmp_path):
    target = tmp_path / "nested" / "validation.png"
    result = module.vourlaki_visual_validation_png(
        str(target), sample_sizes=(4, 4), grid_points=(12, 14, 16), gamma_points=5, verbose=False
    )
    assert _is_png(target)
    assert result["path"] == str(target.resolve())
    data = result["data"]
    assert data["sample_sizes"] == [4, 4]
    assert data["grid_points"] == [12, 14, 16]
    assert data["demographic_params"] == [2, 2, 0.2, 0]
    assert (data["alpha"], data["beta"], data["gamma_pos"]) == (0.2, 10.0, 10.0)
    panels = data["panels"]
    assert [panel["panel"] for panel in panels] == ["A", "B", "C", "D"]
    # The four scenarios of dadi's visual_validation: (ppos_wild, pchange, pchange_pos).
    assert [(p["ppos_wild"], p["pchange"], p["pchange_pos"]) for p in panels] == [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, 0.9),
    ]
    for panel in panels:
        spectrum = panel["spectrum"]
        assert len(spectrum) == 5 and all(len(row) == 5 for row in spectrum)
        values = [v for row in spectrum for v in row if v is not None]
        # The 12-16 point grids of this test leave small integration round-off, so entries may dip
        # slightly below zero; they must stay finite and small against the largest entry.
        assert values and all(math.isfinite(v) for v in values)
        assert min(values) >= -0.01 * max(values) and max(values) > 0
        # The marginal spectra have n + 1 entries for each population.
        assert len(panel["wild_marginal"]) == 5
        assert len(panel["domesticate_marginal"]) == 5


@requires_pkg
def test_vourlaki_validation_refuses_existing_output(tmp_path):
    target = tmp_path / "validation.png"
    target.write_bytes(b"existing")
    with pytest.raises(ConfigurationError, match="already exists"):
        module.vourlaki_visual_validation_png(str(target), sample_sizes=(4, 4), verbose=False)
    assert target.read_bytes() == b"existing"


@pytest.mark.parametrize(
    ("output", "kwargs", "message"),
    [
        ("figure.pdf", {}, "must end with .png"),
        ("", {}, "file path"),
        ("ok.png", {"sample_sizes": (1, 4)}, "sample"),
        ("ok.png", {"sample_sizes": (4,)}, "sample_sizes"),
        ("ok.png", {"grid_points": (12, 14)}, "grid_points"),
        ("ok.png", {"gamma_points": 2}, "gamma_points"),
    ],
)
@requires_pkg
def test_vourlaki_validation_rejects_bad_arguments(tmp_path, output, kwargs, message):
    path = str(tmp_path / output) if output else output
    with pytest.raises(ConfigurationError, match=message):
        module.vourlaki_visual_validation_png(path, verbose=False, **kwargs)


def test_vourlaki_validation_needs_dadi(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError) as info:
        module.vourlaki_visual_validation_png(str(tmp_path / "v.png"), verbose=False)
    assert "bio-dadi" in (info.value.hint or "")


# ---------------------------------------------------------------------------
# MayaVi three-population spectrum
# ---------------------------------------------------------------------------


def test_mayavi_missing_raises_backend_error_with_hint(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "mayavi.mlab", None)
    spectrum = [[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]]
    with pytest.raises(BackendUnavailableError) as info:
        module.plot_3d_spectrum_mayavi(spectrum, str(tmp_path / "m.png"), verbose=False)
    assert "mayavi" in (info.value.hint or "").lower()
    assert not (tmp_path / "m.png").exists()


def test_mayavi_plot_with_fake_mlab_draws_positive_entries(tmp_path, monkeypatch):
    fake = _FakeMlab()
    monkeypatch.setitem(sys.modules, "mayavi.mlab", fake)
    # Masked corner (None) and the value 1.0 below vmin=2.0 are not drawn.
    spectrum = [
        [[None, 1.0], [2.0, 4.0]],
        [[8.0, 2.0], [4.0, 2.0]],
    ]
    target = tmp_path / "three.png"
    result = module.plot_3d_spectrum_mayavi(
        spectrum, str(target), vmin=2.0, vmax=8.0, pop_ids=["AA", "BB", "CC"], verbose=False
    )
    assert _is_png(target)
    data = result["data"]
    assert data["pop_ids"] == ["AA", "BB", "CC"]
    assert data["sample_sizes"] == [1, 1, 1]
    assert (data["vmin"], data["vmax"]) == (2.0, 8.0)
    # Entries with value >= 2.0 are drawn: 2.0, 4.0, 8.0, 2.0, 4.0, 2.0 (six bars).
    assert data["plotted_entries"] == 6
    xs, ys, zs, values = fake.barchart_args
    assert len(xs) == len(ys) == len(zs) == len(values) == 6
    # Colour value is log(v / vmin) / log(vmax / vmin): 2 -> 0, 4 -> 0.5, 8 -> 1.
    assert sorted(values) == pytest.approx([0.0, 0.0, 0.0, 0.5, 0.5, 1.0], abs=1e-12)
    assert fake.axes_kwargs["xlabel"] == "AA"
    assert fake.axes_kwargs["ylabel"] == "BB"
    assert fake.axes_kwargs["zlabel"] == "CC"
    assert fake.axes_obj.axes.label_format == ""
    names = [call[0] for call in fake.calls]
    assert names[-2:] == ["savefig", "close"]


def test_mayavi_colour_range_defaults_to_positive_extremes(tmp_path, monkeypatch):
    fake = _FakeMlab()
    monkeypatch.setitem(sys.modules, "mayavi.mlab", fake)
    spectrum = [[[0.0, 1.0], [2.0, 8.0]]]
    result = module.plot_3d_spectrum_mayavi(spectrum, str(tmp_path / "d.png"), verbose=False)
    assert (result["data"]["vmin"], result["data"]["vmax"]) == (1.0, 8.0)
    assert result["data"]["pop_ids"] == ["pop0", "pop1", "pop2"]


@pytest.mark.parametrize(
    ("spectrum", "kwargs", "message"),
    [
        ([[[0.0, 0.0], [0.0, 0.0]]], {}, "positive"),
        ([[[1.0, 2.0], [3.0, 4.0]]], {"vmin": 0.0}, "vmin"),
        ([[[1.0, 2.0], [3.0, 4.0]]], {"vmin": 4.0, "vmax": 2.0}, "below vmax"),
        ([[1.0, 2.0], [3.0, 4.0]], {}, "axes"),
        ([[[1.0, 2.0], [3.0, 4.0]]], {"pop_ids": ["A", "B"]}, "pop_ids"),
    ],
)
def test_mayavi_rejects_bad_input_before_drawing(tmp_path, monkeypatch, spectrum, kwargs, message):
    fake = _FakeMlab()
    monkeypatch.setitem(sys.modules, "mayavi.mlab", fake)
    with pytest.raises(ConfigurationError, match=message):
        module.plot_3d_spectrum_mayavi(spectrum, str(tmp_path / "bad.png"), verbose=False, **kwargs)
    assert fake.calls == []


@pytest.mark.skipif(
    importlib.util.find_spec("mayavi") is None,
    reason="MayaVi is not installed here; the real 3D render needs mayavi and a display",
)
def test_mayavi_real_render(tmp_path):
    spectrum = [[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]]
    result = module.plot_3d_spectrum_mayavi(spectrum, str(tmp_path / "real.png"), verbose=False)
    assert _is_png(tmp_path / "real.png")
    assert result["data"]["plotted_entries"] == 8
