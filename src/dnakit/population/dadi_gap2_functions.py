"""dadi wrappers for the rows that ``dadi_gap_functions`` does not cover.

``vourlaki_visual_validation_png`` redraws dadi's ``DFE.Vourlaki2022.visual_validation`` figure
(four Vourlaki 2022 DFE scenarios with their two-population spectra and marginal spectra). Upstream
writes ``validation.pdf`` into the working directory and calls ``plt.show()``; this version takes an
explicit PNG path, keeps the same model and parameter sets, and returns the spectra it drew.
``plot_3d_spectrum_mayavi`` draws a three-population spectrum as a MayaVi 3D bar chart and saves it
as a PNG. MayaVi is an optional package and is not part of the dadi extra; without it the function
raises ``BackendUnavailableError`` with an install hint.

dadi, NumPy, matplotlib and MayaVi are imported only when a function runs. No network is used.
Output files are created new: an existing path is refused and its parent directory is created.
"""

from __future__ import annotations

import importlib
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _integer, _number, _provenance
from dnakit.population.dadi_functions import _nested
from dnakit.population.dadi_gap_functions import _array, _pop_ids, _pyplot

_PNG_DPI = 150
_MAX_SAMPLE = 40
_MAX_GRID_POINT = 200
_MAX_GAMMA_POINTS = 60
_MAYAVI_HINT = (
    "Install MayaVi with: python -m pip install mayavi. The 3D plot needs a display or an "
    "offscreen VTK build."
)
# Upstream scenario order: (ppos_wild, pchange, pchange_pos) for panels A to D.
_VALIDATION_PANELS = (
    ("A", 0.0, 0.0, 0.0),
    ("B", 1.0, 0.0, 0.0),
    ("C", 0.0, 1.0, 0.0),
    ("D", 0.0, 1.0, 0.9),
)
_VOURLAKI_DEMOGRAPHY = (2, 2, 0.2, 0)  # nu1, nu2, T, m (as in dadi's visual_validation)
_VOURLAKI_ALPHA = 0.2
_VOURLAKI_BETA = 10.0
_VOURLAKI_GAMMA_POS = 10.0
_VOURLAKI_GAMMA_BOUNDS = (1e-4, 200.0)


def _output_png(output_path: object) -> Path:
    """A new PNG path; its parent directory is created. Existing files are never overwritten."""

    if not isinstance(output_path, (str, Path)) or not str(output_path).strip():
        raise ConfigurationError("output_path must be a file path.")
    target = Path(output_path).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _sample_sizes(value: object) -> tuple[int, int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError("sample_sizes must list two integer sample sizes.")
    return (
        _integer(value[0], "sample size", minimum=2, maximum=_MAX_SAMPLE),
        _integer(value[1], "sample size", minimum=2, maximum=_MAX_SAMPLE),
    )


def _grid(value: object) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigurationError("grid_points must list three integer grid sizes.")
    return [_integer(item, "grid point", minimum=8, maximum=_MAX_GRID_POINT) for item in value]


def vourlaki_visual_validation_png(
    output_path: str,
    *,
    sample_sizes: Sequence[int] = (10, 10),
    grid_points: Sequence[int] = (30, 35, 40),
    gamma_points: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw dadi's Vourlaki (2022) validation figure to a PNG and return the spectra it shows.

    The demography is dadi's visual-validation model (nu1 = nu2 = 2, T = 0.2, m = 0) with
    alpha = 0.2 and beta = 10 for the negative gamma DFE, and a positive point mass at gamma = 10.
    Four scenarios (panels A to D) vary ppos_wild, pchange and pchange_pos in ``Vourlaki_mixture``.
    Each panel is a two-population spectrum; the lower panels show the wild and domesticated
    marginal spectra on a log scale. ``sample_sizes`` are the two sample sizes, ``grid_points``
    the three integration grids and ``gamma_points`` the number of cached gamma values. Defaults
    match upstream. The returned data holds each panel's spectrum and its parameters.
    """

    dadi = _dadi_module("dadi")
    np = require("numpy", _DADI_EXTRA)
    plotting_vourlaki = _vourlaki_module()
    target = _output_png(output_path)
    ns = _sample_sizes(sample_sizes)
    pts = _grid(grid_points)
    gammas = _integer(gamma_points, "gamma_points", minimum=5, maximum=_MAX_GAMMA_POINTS)
    demo = list(_VOURLAKI_DEMOGRAPHY)
    cache1d = _dadi_cache("Cache1D")
    cache2d = _dadi_cache("Cache2D")
    selection = _dadi_module("dadi.DFE.DemogSelModels")
    spectra_by_panel: list[dict[str, Any]] = []
    plt = _pyplot()
    figure = None
    try:
        s1 = cache1d(
            tuple(demo),
            list(ns),
            selection.split_mig_sel_single_gamma,
            pts,
            gamma_bounds=_VOURLAKI_GAMMA_BOUNDS,
            gamma_pts=gammas,
            additional_gammas=[_VOURLAKI_GAMMA_POS],
            verbose=False,
        )
        s2 = cache2d(
            tuple(demo),
            list(ns),
            selection.split_mig_sel,
            pts,
            gamma_bounds=_VOURLAKI_GAMMA_BOUNDS,
            gamma_pts=gammas,
            additional_gammas=[_VOURLAKI_GAMMA_POS],
            verbose=False,
        )
        figure = plt.figure(figsize=(10, 8))
        figure.clear()
        axis_wild = figure.add_subplot(3, 2, 5)
        axis_domesticate = figure.add_subplot(3, 2, 6)
        for index, (label, wild, change, change_pos) in enumerate(_VALIDATION_PANELS):
            params = [_VOURLAKI_ALPHA, _VOURLAKI_BETA, wild, _VOURLAKI_GAMMA_POS, change, change_pos]
            fs = plotting_vourlaki(params, None, s1, s2, 1.0, None)
            axis = figure.add_subplot(3, 2, index + 1)
            dadi.Plotting.plot_single_2d_sfs(fs, ax=axis, show=False)
            axis.set_title(
                f"ppos_wild={wild:.1f}, pchange={change:.1f}, pchange_pos={change_pos:.1f}"
            )
            axis.text(0.02, 0.98, label, fontsize="x-large", va="top", transform=axis.transAxes)
            wild_marginal = fs.filter_pops(tokeep=[1])
            domesticate_marginal = fs.filter_pops(tokeep=[2])
            axis_wild.semilogy(wild_marginal, "-o", label=label)
            axis_domesticate.semilogy(domesticate_marginal, "-o", label=label)
            spectra_by_panel.append(
                {
                    "panel": label,
                    "ppos_wild": wild,
                    "pchange": change,
                    "pchange_pos": change_pos,
                    "spectrum": _nested(np, fs),
                    "wild_marginal": _nested(np, wild_marginal),
                    "domesticate_marginal": _nested(np, domesticate_marginal),
                }
            )
        axis_wild.set_title("wild")
        axis_domesticate.set_title("domesticate")
        axis_wild.legend()
        axis_domesticate.legend()
        figure.tight_layout()
        figure.savefig(str(target), format="png", dpi=_PNG_DPI)
    except (ValueError, TypeError, IndexError, KeyError, AttributeError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi could not draw the Vourlaki validation figure: {exc}") from exc
    finally:
        if figure is not None:
            plt.close(figure)
    return {
        "path": str(target),
        "data": {
            "sample_sizes": list(ns),
            "grid_points": pts,
            "gamma_points": gammas,
            "demographic_params": demo,
            "alpha": _VOURLAKI_ALPHA,
            "beta": _VOURLAKI_BETA,
            "gamma_pos": _VOURLAKI_GAMMA_POS,
            "panels": spectra_by_panel,
            "provenance": _provenance(
                "dadi",
                "Vourlaki 2022 DFE validation scenarios, redrawn from dadi's visual_validation",
            ),
        },
    }


def _dadi_module(dotted: str) -> Any:
    return require(dotted, _DADI_EXTRA)


def _vourlaki_module() -> Any:
    """The Vourlaki mixture function (dadi.DFE.Vourlaki2022.Vourlaki_mixture)."""

    return _dadi_module("dadi.DFE.Vourlaki2022").Vourlaki_mixture


def _dadi_cache(name: str) -> Any:
    module = "dadi.DFE.Cache1D_mod" if name == "Cache1D" else "dadi.DFE.Cache2D_mod"
    return getattr(_dadi_module(module), name)


def _mayavi_mlab() -> Any:
    try:
        return importlib.import_module("mayavi.mlab")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"mayavi.mlab could not be imported: {exc}", hint=_MAYAVI_HINT
        ) from exc


def plot_3d_spectrum_mayavi(
    spectrum: Sequence[Any],
    output_path: str,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    pop_ids: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a three-population spectrum as a MayaVi 3D bar chart and save it as a PNG.

    ``spectrum`` is a 3D nested list of counts (masked entries as None). Each entry from ``vmin``
    up is a bar, coloured by its logarithm scaled between ``vmin`` and ``vmax``. Both default to
    the smallest and largest positive counts; ``vmin`` must be above zero (the colour scale is
    logarithmic). This follows dadi's ``plot_3d_spectrum_mayavi`` and uses MayaVi's ``mlab`` API
    directly. ``pop_ids`` names the three axes. The returned data holds the plotted range and the
    spectrum, not the figure. MayaVi is required (``mayavi.mlab``); without it a
    BackendUnavailableError is raised.
    """

    np = require("numpy", _DADI_EXTRA)
    array = _array(np, spectrum, "spectrum", (3,), signed=False)
    names = _pop_ids(pop_ids, 3) if pop_ids is not None else ["pop0", "pop1", "pop2"]
    low, high = _colour_range(np, array, vmin, vmax)
    target = _output_png(output_path)
    mlab = _mayavi_mlab()
    toplot = (~np.ma.getmaskarray(array)) & (np.ma.getdata(array) >= low)
    values = np.ma.getdata(array)[toplot].astype(float)
    normalized = np.minimum((np.log(values) - math.log(low)) / (math.log(high) - math.log(low)), 1.0)
    xs, ys, zs = np.indices(array.shape)
    figure = mlab.figure(bgcolor=(1, 1, 1), size=(800, 800))
    mlab.clf(figure)
    mlab.barchart(
        xs[toplot],
        ys[toplot],
        zs[toplot],
        normalized,
        colormap="hsv",
        scale_mode="none",
        lateral_scale=1,
        figure=figure,
    )
    axes = mlab.axes(
        xlabel=names[0], ylabel=names[1], zlabel=names[2], figure=figure, color=(0, 0, 0)
    )
    axes.axes.label_format = ""
    mlab.savefig(str(target), figure=figure)
    mlab.close(figure)
    return {
        "path": str(target),
        "data": {
            "pop_ids": names,
            "sample_sizes": [int(size - 1) for size in array.shape],
            "vmin": low,
            "vmax": high,
            "plotted_entries": int(toplot.sum()),
            "spectrum": _nested(np, array),
            "provenance": _provenance("dadi", "three-population spectrum drawn with MayaVi mlab"),
        },
    }


def _colour_range(np: Any, array: Any, vmin: object, vmax: object) -> tuple[float, float]:
    present = np.ma.getdata(array)[~np.ma.getmaskarray(array)]
    positive = present[present > 0]
    if positive.size == 0:
        raise ConfigurationError("spectrum must contain at least one positive, unmasked count.")
    low = float(positive.min()) if vmin is None else _number(vmin, "vmin", positive=True)
    high = float(positive.max()) if vmax is None else _number(vmax, "vmax", positive=True)
    if low >= high:
        raise ConfigurationError("vmin must be below vmax.")
    return low, high


__all__ = [
    "plot_3d_spectrum_mayavi",
    "vourlaki_visual_validation_png",
]

_register_public_verbose(globals(), tuple(__all__))
