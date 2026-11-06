"""dadi operations that DNAKit runs without the calls that fail in the installed packages.

Four dadi functions cannot run as written with the installed matplotlib 3.11, dadi 2.4.4, demes 0.2.3
and ruamel.yaml 0.19.1. Each one is implemented here with the same arithmetic and the same drawing:

- ``dadi.Triallele.plotting.plot_single_trispectrum`` calls ``pcolor`` with ``norm=LogNorm`` and
  ``vmin``/``vmax`` together, which matplotlib 3.11 rejects. The plot is drawn with the norm alone.
- ``dadi.Triallele.plotting.plot_trispectrum_comp`` draws its two panels with that same call; it is
  implemented on the same single-panel drawing.
- ``dadi.TwoLocus.plotting.plot_3d`` calls ``ax.colorbar()``, which an ``Axes3D`` does not have. The
  voxels are drawn as upstream does and the colour bar is added with ``fig.colorbar``.
- ``dadi.DFE.Cache2D_mod.mixture_point_pos`` passes ``None`` as the ``rho`` argument of
  ``Cache2D.integrate_point_pos``, which multiplies it and fails. ``rho`` is passed instead.
- ``dadi.Demes.Inference.optimize`` reads its options file with ``ruamel.yaml.load``, which
  ruamel.yaml 0.19 removed. The options are read with ``YAML(typ="safe")``; every other step is
  dadi's own helper.

Plots are written to new PNG files and return the path together with the spectra that were drawn.
dadi, NumPy, SciPy and matplotlib are imported only when a function runs. Output files are never
overwritten.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _integer, _number, _provenance
from dnakit.population.dadi_functions import _nested, _new_file, _numbers, _sizes

_PLOT_DPI = 100
_MAX_GRID_SIZE = 200
_MAX_CACHE_GAMMA_POINTS = 60
_MAX_CACHE_SAMPLE = 40
_MAX_DEMES_ITERATIONS = 10_000
_DEMES_METHODS = ("fmin", "powell", "lbfgsb")
_OPERATION_ERRORS = (ValueError, TypeError, IndexError, KeyError, AttributeError, ZeroDivisionError)
_TRISPECTRUM_CMAP = "cubehelix_r"
_RESIDUAL_CMAP = "RdBu"
_SCRUNCH = (0.2, 0.8)


def _dadi_module(dotted: str) -> Any:
    return require(dotted, _DADI_EXTRA)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _DADI_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _DADI_EXTRA)


def _cell(value: object, name: str) -> Any:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_cell(item, name) for item in value]
    if value is None:
        return math.nan
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must hold finite numbers, or None for masked entries.")
    return float(value)


def _masked(np: Any, value: object, name: str, ndims: tuple[int, ...]) -> Any:
    """A nested list as a masked array; None marks an entry that dadi masks."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(
            f"{name} must be a non-empty nested list; masked entries are None."
        )
    try:
        raw = np.array(_cell(value, name), dtype=float)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be a rectangular nested list.") from exc
    if raw.ndim not in ndims:
        raise ConfigurationError(f"{name} must have {list(ndims)} axes, found {raw.ndim}.")
    data = np.ma.masked_invalid(raw)
    present = np.ma.getdata(data)[~np.ma.getmaskarray(data)]
    if (present < 0).any():
        raise ConfigurationError(f"{name} must contain only non-negative values or None.")
    return data


def _reals(value: object, name: str) -> list[float]:
    """Finite numbers of any sign (the lognormal mu of a DFE is usually negative)."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of numbers.")
    return [_number_any(item, f"{name} entry") for item in value]


def _number_any(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _colour_limits(
    np: Any, arrays: Sequence[Any], vmin: object, vmax: object
) -> tuple[float, float]:
    """Log-scale limits: vmin defaults to the smallest positive count, vmax to the largest."""
    present = [np.ma.getdata(item)[~np.ma.getmaskarray(item)] for item in arrays]
    values = np.concatenate(present) if present else np.array([])
    if values.size == 0:
        raise ConfigurationError("the spectrum has no unmasked entries to draw.")
    positive = values[values > 0]
    if vmin is None:
        if positive.size == 0:
            raise ConfigurationError("the spectrum has no positive entries for a log colour scale.")
        low = float(positive.min())
    else:
        low = _number(vmin, "vmin", positive=True)
    high = float(values.max()) if vmax is None else _number(vmax, "vmax", positive=True)
    if not high > low:
        raise ConfigurationError("vmax must be greater than vmin.")
    return low, high


def _colormap(matplotlib: Any, name: object) -> Any:
    if not isinstance(name, str) or not name:
        raise ConfigurationError("cmap must be the name of a matplotlib colormap.")
    try:
        return matplotlib.colormaps[name]
    except KeyError as exc:
        raise ConfigurationError(f"cmap {name!r} is not a matplotlib colormap.") from exc


def _fold(np: Any, sfs: Any) -> Any:
    """Fold a triallelic spectrum, as dadi's Triallele plotting ``_fold`` does (same masks)."""
    if sfs.ndim != 2 or sfs.shape[0] != sfs.shape[1] or sfs.shape[0] < 3:
        raise ConfigurationError("a triallelic spectrum must be a square array of at least 3 by 3.")
    if np.ma.getmaskarray(sfs)[1, 2]:
        raise ConfigurationError("the spectrum is already folded (entry [1, 2] is masked).")
    n = sfs.shape[0]
    summed = sfs + np.ma.transpose(sfs)
    data = np.array(np.ma.getdata(summed), dtype=float)
    for ii in range(n):
        data[ii, ii] = data[ii, ii] / 2
    mask = np.array(np.ma.getmaskarray(summed), dtype=bool)
    mask[0, :] = True
    mask[:, 0] = True
    for ii in range(n):
        mask[ii, ii + 1 :] = True
        mask[ii, n - 1 - ii :] = True
    return np.ma.array(data, mask=mask)


def _shown(np: Any, sfs: Any, low: float) -> Any:
    """The entries that are drawn: unmasked and at least vmin."""
    below = np.ma.getdata(sfs) < low
    return np.ma.masked_where(below | np.ma.getmaskarray(sfs), np.ma.getdata(sfs))


def _draw_single(
    matplotlib: Any, np: Any, ax: Any, sfs: Any, *, folded: bool, low: float, high: float, cmap: Any
) -> tuple[Any, Any]:
    shown = _shown(np, sfs, low)
    grid = np.transpose(shown) if folded else shown
    norm = matplotlib.colors.LogNorm(vmin=low, vmax=high)
    mappable = ax.pcolor(grid, norm=norm, cmap=cmap)
    ax.set_xlim(0, sfs.shape[1])
    if folded:
        ax.set_ylim(0, sfs.shape[0] / 2)
        ax.set_xlabel("Major derived allele frequency")
        ax.set_ylabel("Minor derived allele frequency")
    else:
        ax.set_ylim(0, sfs.shape[0])
        ax.set_xlabel("Frequency")
        ax.set_ylabel("Frequency")
    return mappable, shown


def _png_target(image_path: object) -> Path:
    target = _new_file(image_path, "image_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    return target


def plot_triallelic_spectrum(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    folded: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = False,
    cmap: str = _TRISPECTRUM_CMAP,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a triallelic spectrum as a log-scale heat map and write it to a new PNG file.

    ``spectrum`` is a square nested list of counts; ``None`` marks an entry that is masked. With
    ``folded=True`` the spectrum is folded first, as dadi's ``plot_single_trispectrum`` does. Entries
    below ``vmin`` are not drawn. The colour map is the upstream one (``cubehelix_r``); upstream
    ignores its ``cmap`` argument, so here it is honoured. When ``vmin`` is omitted the smallest
    positive count is used, because a log scale cannot start at zero.
    """
    plt = _pyplot()
    matplotlib = require("matplotlib", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    target = _png_target(image_path)
    sfs = _masked(np, spectrum, "spectrum", (2,))
    is_folded = bool(folded)
    if is_folded:
        sfs = _fold(np, sfs)
    low, high = _colour_limits(np, [sfs], vmin, vmax)
    colours = _colormap(matplotlib, cmap)
    figure, ax = plt.subplots()
    try:
        mappable, shown = _draw_single(
            matplotlib, np, ax, sfs, folded=is_folded, low=low, high=high, cmap=colours
        )
        if colorbar:
            bar = figure.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, shrink=0.85)
            bar.set_label("Count", fontsize=10)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"matplotlib could not draw the spectrum: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_triallelic_spectrum",
        "image_path": str(target),
        "image_format": "png",
        "folded": is_folded,
        "vmin": low,
        "vmax": high,
        "cmap": cmap,
        "spectrum": _nested(np, shown),
        "provenance": _provenance(
            "dadi",
            "triallelic spectrum drawn as in dadi plot_single_trispectrum with a LogNorm alone; matplotlib",
        ),
    }


def plot_triallelic_spectrum_comparison(
    first: Sequence[Any],
    second: Sequence[Any],
    image_path: str | Path,
    *,
    folded: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    resid_range: float | None = None,
    colorbar: bool = False,
    titles: Sequence[str] = ("sfs1", "sfs2", "Residual"),
    cmap: str = _TRISPECTRUM_CMAP,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw two triallelic spectra and their residual in three stacked panels, as a new PNG file.

    The residual is ``(first - second) / sqrt(first)``, as in dadi's ``plot_trispectrum_comp``. The
    two spectra share one colour scale; the residual uses a diverging colour map centred on zero,
    spanning ``-resid_range`` to ``resid_range`` or, when that is omitted, the largest absolute value.
    """
    plt = _pyplot()
    matplotlib = require("matplotlib", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    target = _png_target(image_path)
    first_sfs = _masked(np, first, "first", (2,))
    second_sfs = _masked(np, second, "second", (2,))
    if first_sfs.shape != second_sfs.shape:
        raise ConfigurationError("first and second must have the same shape.")
    if not isinstance(titles, Sequence) or isinstance(titles, str) or len(titles) != 3:
        raise ConfigurationError("titles must list three panel titles.")
    is_folded = bool(folded)
    if is_folded:
        first_sfs = _fold(np, first_sfs)
        second_sfs = _fold(np, second_sfs)
    low, high = _colour_limits(np, [first_sfs, second_sfs], vmin, vmax)
    colours = _colormap(matplotlib, cmap)
    with np.errstate(divide="ignore", invalid="ignore"):
        residual = (first_sfs - second_sfs) / np.sqrt(first_sfs)
    residual = np.ma.masked_invalid(residual)
    if is_folded:
        residual = _fold(np, residual)
    if resid_range is None:
        span = float(np.max(np.abs(np.ma.getdata(residual)[~np.ma.getmaskarray(residual)])))
    else:
        span = _number(resid_range, "resid_range", positive=True)
    figure, axes = plt.subplots(3, 1, figsize=(6, 14))
    try:
        mappables = []
        for ax, sfs, title in zip(axes[:2], (first_sfs, second_sfs), titles[:2], strict=True):
            mappable, _ = _draw_single(
                matplotlib, np, ax, sfs, folded=is_folded, low=low, high=high, cmap=colours
            )
            ax.set_title(title)
            mappables.append(mappable)
        ax3 = axes[2]
        ax3.set_title(titles[2])
        grid = np.transpose(residual) if is_folded else residual
        residual_map = ax3.pcolor(grid, vmin=-span, vmax=span, cmap=_RESIDUAL_CMAP)
        ax3.set_xlim(0, first_sfs.shape[1])
        if is_folded:
            ax3.set_ylim(0, first_sfs.shape[0] / 2)
        else:
            ax3.set_ylim(0, first_sfs.shape[0])
        if colorbar:
            for ax, mappable in zip(axes[:2], mappables, strict=True):
                bar = figure.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, shrink=0.85)
                bar.set_label("Count", fontsize=10)
            bar = figure.colorbar(residual_map, ax=ax3, fraction=0.046, pad=0.04, shrink=0.85)
            bar.set_label("Residual", fontsize=10)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"matplotlib could not draw the comparison: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_triallelic_spectrum_comparison",
        "image_path": str(target),
        "image_format": "png",
        "folded": is_folded,
        "vmin": low,
        "vmax": high,
        "resid_range": span,
        "first": _nested(np, _shown(np, first_sfs, low)),
        "second": _nested(np, _shown(np, second_sfs, low)),
        "residual": _nested(np, residual),
        "provenance": _provenance(
            "dadi",
            "triallelic comparison as in dadi plot_trispectrum_comp; the panels use the single-panel draw",
        ),
    }


def _voxel_faces(ii: int, jj: int, kk: int) -> list[list[list[float]]]:
    """The six faces of voxel (ii, jj, kk), in the order dadi's TwoLocus plot_3d builds them."""
    return [
        [
            [ii - 0.5, jj + 0.5, kk + 0.5],
            [ii + 0.5, jj + 0.5, kk + 0.5],
            [ii + 0.5, jj - 0.5, kk + 0.5],
            [ii - 0.5, jj - 0.5, kk + 0.5],
        ],
        [
            [ii - 0.5, jj + 0.5, kk - 0.5],
            [ii + 0.5, jj + 0.5, kk - 0.5],
            [ii + 0.5, jj - 0.5, kk - 0.5],
            [ii - 0.5, jj - 0.5, kk - 0.5],
        ],
        [
            [ii - 0.5, jj + 0.5, kk + 0.5],
            [ii + 0.5, jj + 0.5, kk + 0.5],
            [ii + 0.5, jj + 0.5, kk - 0.5],
            [ii - 0.5, jj + 0.5, kk - 0.5],
        ],
        [
            [ii - 0.5, jj - 0.5, kk + 0.5],
            [ii + 0.5, jj - 0.5, kk + 0.5],
            [ii + 0.5, jj - 0.5, kk - 0.5],
            [ii - 0.5, jj - 0.5, kk - 0.5],
        ],
        [
            [ii + 0.5, jj - 0.5, kk + 0.5],
            [ii + 0.5, jj + 0.5, kk + 0.5],
            [ii + 0.5, jj + 0.5, kk - 0.5],
            [ii + 0.5, jj - 0.5, kk - 0.5],
        ],
        [
            [ii - 0.5, jj - 0.5, kk + 0.5],
            [ii - 0.5, jj + 0.5, kk + 0.5],
            [ii - 0.5, jj + 0.5, kk - 0.5],
            [ii - 0.5, jj - 0.5, kk - 0.5],
        ],
    ]


def _visible_faces(toplot: Any, ii: int, jj: int, kk: int) -> list[int]:
    """Which of the six faces of a drawn voxel touch no other drawn voxel (dadi's rule)."""
    nx, ny, nz = toplot.shape
    checks = [
        (kk < nz - 1 and toplot[ii, jj, kk + 1]),
        (kk > 0 and toplot[ii, jj, kk - 1]),
        (jj < ny - 1 and toplot[ii, jj + 1, kk]),
        (jj > 0 and toplot[ii, jj - 1, kk]),
        (ii < nx - 1 and toplot[ii + 1, jj, kk]),
        (ii > 0 and toplot[ii - 1, jj, kk]),
    ]
    return [face for face, hidden in enumerate(checks) if not hidden]


def plot_two_locus_spectrum_3d(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    max_sum: int | None = None,
    max_index: int | None = None,
    colorbar: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a three-dimensional two-locus spectrum as coloured voxels in a new PNG file.

    ``spectrum`` is a cubic nested list of counts with ``None`` for masked entries. Voxels below
    ``vmin`` are not drawn. ``max_sum`` hides voxels whose index sum exceeds it, and ``max_index``
    hides voxels whose three indices all exceed it, as in dadi's ``plot_3d``. Colours follow dadi's
    scrunched ``cubehelix_r`` scale on a log axis; the colour bar is drawn on the same scale.
    """
    plt = _pyplot()
    matplotlib = require("matplotlib", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    target = _png_target(image_path)
    fs = _masked(np, spectrum, "spectrum", (3,))
    low, high = _colour_limits(np, [fs], vmin, vmax)
    if max_sum is not None:
        max_sum = _integer(max_sum, "max_sum", minimum=0)
    if max_index is not None:
        max_index = _integer(max_index, "max_index", minimum=0)
    shown = _shown(np, fs, low)
    present = ~np.ma.getmaskarray(shown)
    # dadi swaps the first and third axes before drawing; the index sums are unchanged by that.
    drawn = np.transpose(present, (2, 1, 0)).copy()
    data = np.transpose(np.ma.getdata(fs), (2, 1, 0)).copy()
    if max_sum is not None:
        ii, jj, kk = np.indices(drawn.shape)
        drawn &= (ii + jj + kk) <= max_sum
    if max_index is not None:
        ii, jj, kk = np.indices(drawn.shape)
        drawn &= ~((ii > max_index) & (jj > max_index) & (kk > max_index))
    normalized = np.zeros(drawn.shape, dtype=float)
    cells = np.argwhere(drawn)
    if cells.size:
        values = data[drawn]
        normalized[drawn] = np.minimum(
            (np.log(values) - np.log(low)) / (np.log(high) - np.log(low)), 1
        )
    cube = matplotlib.colormaps["cubehelix_r"]
    scrunched_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "cubehelix_r_scrunched", cube(np.linspace(*_SCRUNCH, 256))
    )
    colours = cube(0.2 + 0.6 * normalized)
    polys, polycolors = [], []
    for ii, jj, kk in cells:
        ii, jj, kk = int(ii), int(jj), int(kk)
        for face in _visible_faces(drawn, ii, jj, kk):
            polys.append(_voxel_faces(ii, jj, kk)[face])
            polycolors.append(colours[ii, jj, kk])
    figure = plt.figure()
    ax = figure.add_subplot(projection="3d")
    try:
        poly3d = require("mpl_toolkits.mplot3d.art3d", _DADI_EXTRA)
        polycoll = poly3d.Poly3DCollection(
            polys, facecolor=polycolors, edgecolor="k", linewidths=0.5
        )
        ax.add_collection3d(polycoll)
        ax.set_xlim3d(-0.5, drawn.shape[0] - 0.5)
        ax.set_ylim3d(-0.5, drawn.shape[1] - 0.5)
        ax.set_zlim3d(-0.5, drawn.shape[2] - 0.5)
        ax.set_xlabel("aB", horizontalalignment="left")
        ax.set_ylabel("Ab", verticalalignment="bottom")
        ax.set_zlabel("AB", verticalalignment="bottom")
        ax.view_init(elev=30.0, azim=45)
        if colorbar:
            mappable = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.LogNorm(vmin=low, vmax=high), cmap=scrunched_cmap
            )
            mappable.set_array([])
            bar = figure.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
            bar.set_label("Count", fontsize=10)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"matplotlib could not draw the 3D spectrum: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_two_locus_spectrum_3d",
        "image_path": str(target),
        "image_format": "png",
        "vmin": low,
        "vmax": high,
        "max_sum": max_sum,
        "max_index": max_index,
        "drawn_entries": int(drawn.sum()),
        "spectrum": _nested(np, _shown(np, fs, low)),
        "provenance": _provenance(
            "dadi",
            "two-locus spectrum drawn as voxels as in dadi TwoLocus plot_3d, colour bar added with fig.colorbar",
        ),
    }


def dfe_mixture_point_pos_spectrum(
    demographic_params: Sequence[float],
    sample_sizes: Sequence[int],
    dfe_params: Sequence[float],
    theta: float = 1.0,
    *,
    grid_points: Sequence[int] = (20, 25, 30),
    gamma_points: int = 20,
    gamma_max: float = 2000.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected two-population spectrum for dadi's point-mass DFE mixture (``mixture_point_pos``).

    ``demographic_params`` are nu1, nu2, T and m for dadi's ``split_mig_sel`` demography.
    ``dfe_params`` are eight values: the lognormal mu and sigma, then rho, the positive proportion
    and gamma for population 1, the positive proportion and gamma for population 2, and the weight of
    the two-dimensional part. The result is ``(1 - weight) * spectrum_1D + weight * spectrum_2D``, where
    the two-dimensional part receives ``rho`` as its ``rho`` argument (upstream passes ``None``).
    """
    cache_mod = _dadi_module("dadi.DFE.Cache1D_mod")
    np = require("numpy", _DADI_EXTRA)
    demo = _numbers(demographic_params, "demographic_params")
    if len(demo) != 4:
        raise ConfigurationError("demographic_params must be nu1, nu2, T, m (four values).")
    ns = _sizes(sample_sizes, 2, _MAX_CACHE_SAMPLE, "sample_sizes")
    pts = _grid_points(grid_points)
    gammas = _integer(gamma_points, "gamma_points", minimum=5, maximum=_MAX_CACHE_GAMMA_POINTS)
    ceiling = _number(gamma_max, "gamma_max", positive=True)
    values = _reals(dfe_params, "dfe_params")
    if len(values) != 8:
        raise ConfigurationError(f"mixture_point_pos needs 8 dfe_params, got {len(values)}.")
    scale = _number(theta, "theta", positive=True)
    pdf_params = list(values[:-6])
    rho, ppos1, gamma_pos1, ppos2, gamma_pos2, p2d = values[-6:]
    if not 0 <= p2d <= 1:
        raise ConfigurationError("the weight of the two-dimensional part must be between 0 and 1.")
    selection = _dadi_module("dadi.DFE.DemogSelModels")
    cache1d = cache_mod.Cache1D
    cache2d = _dadi_module("dadi.DFE.Cache2D_mod").Cache2D
    pdfs = _dadi_module("dadi.DFE.PDFs")
    point_gammas = [gamma_pos1, gamma_pos2]
    try:
        s1 = cache1d(
            tuple(demo),
            list(ns),
            selection.split_mig_sel_single_gamma,
            pts,
            gamma_bounds=(1e-4, ceiling),
            gamma_pts=gammas,
            additional_gammas=point_gammas,
            verbose=False,
        )
        s2 = cache2d(
            tuple(demo),
            list(ns),
            selection.split_mig_sel,
            pts,
            gamma_bounds=(1e-4, ceiling),
            gamma_pts=gammas,
            additional_gammas=point_gammas,
            verbose=False,
        )
        fs1 = s1.integrate_point_pos(
            pdf_params + [ppos1, gamma_pos1], None, pdfs.lognormal, scale, Npos=1, pts=None
        )
        fs2 = s2.integrate_point_pos(
            pdf_params + [rho, ppos1, gamma_pos1, ppos2, gamma_pos2],
            None,
            pdfs.biv_lognormal,
            scale,
            rho,
            None,
        )
        spectrum = (1 - p2d) * fs1 + p2d * fs2
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"dadi could not compute mixture_point_pos: {exc}") from exc
    return {
        "status": "completed",
        "operation": "dfe_mixture_point_pos_spectrum",
        "variant": "mixture_point_pos",
        "demographic_params": demo,
        "sample_sizes": list(ns),
        "dfe_params": values,
        "theta": scale,
        "grid_points": pts,
        "gamma_points": gammas,
        "gamma_max": ceiling,
        "spectrum": _nested(np, spectrum),
        "provenance": _provenance(
            "dadi",
            "point-mass DFE mixture, split_mig_sel caches; rho passed to the 2D integration",
        ),
    }


def _grid_points(value: object) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigurationError("grid_points must list three integer grid sizes.")
    return [_integer(item, "grid point", minimum=8, maximum=_MAX_GRID_SIZE) for item in value]


def _read_inference_options(path: Path) -> dict[str, Any]:
    """Read dadi's inference options with the current ruamel.yaml API (dadi calls the removed load)."""
    yaml = require("ruamel.yaml", _DADI_EXTRA)
    with open(path, encoding="utf-8") as handle:
        options = yaml.YAML(typ="safe").load(handle)
    if not isinstance(options, dict) or "parameters" not in options:
        raise ConfigurationError("parameters to fit must be specified in the inference options.")
    return options


def optimize_demes_graph_parameters(
    graph_path: str | Path,
    inference_options_path: str | Path,
    sfs: Sequence[Any],
    pop_ids: Sequence[str],
    grid_points: Sequence[int],
    *,
    maxiter: int = 1000,
    perturb: float = 0,
    log: bool = True,
    method: str = "fmin",
    uL: float | None = None,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit the parameters of a demes graph to a two- or three-population spectrum with dadi.

    This is dadi's ``Demes.Inference.optimize`` with the same optimisation steps (parameter set-up,
    bounds, constraints, the chosen optimiser and the log parameterisation). ``graph_path`` is a
    demes YAML file and ``inference_options_path`` lists the parameters to fit. ``sfs`` is the observed
    spectrum as a nested list with ``pop_ids`` naming its populations; ``grid_points`` gives three
    grid sizes. ``output_path`` (optional, new file) receives the fitted demes model as YAML.
    """
    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    scipy_optimize = require("scipy.optimize", _DADI_EXTRA)
    inference = _dadi_module("dadi.Demes.Inference")
    graph = Path(graph_path).expanduser()
    options_file = Path(inference_options_path).expanduser()
    if not graph.is_file() or not options_file.is_file():
        raise ConfigurationError("graph_path and inference_options_path must be existing files.")
    if isinstance(pop_ids, (str, bytes)) or not isinstance(pop_ids, Sequence):
        raise ConfigurationError("pop_ids must be a list of population names.")
    ids = [str(item) for item in pop_ids]
    values = _masked(np, sfs, "sfs", (1, 2, 3))
    data = dadi.Spectrum(values)
    if len(ids) != data.ndim or len(set(ids)) != data.ndim or any(not item for item in ids):
        raise ConfigurationError(f"pop_ids must name each of the {data.ndim} populations once.")
    data.pop_ids = ids
    pts = _grid_points(grid_points)
    if method not in _DEMES_METHODS:
        raise ConfigurationError(f"method must be one of {sorted(_DEMES_METHODS)}.")
    iterations = _integer(maxiter, "maxiter", minimum=1, maximum=_MAX_DEMES_ITERATIONS)
    if isinstance(perturb, bool) or not isinstance(perturb, Real) or perturb < 0:
        raise ConfigurationError("perturb must be a non-negative number.")
    if not isinstance(log, bool):
        raise ConfigurationError("log must be true or false.")
    if uL is not None:
        uL = _number(uL, "uL", positive=True)
    target = _new_file(output_path, "output_path") if output_path is not None else None
    try:
        builder = inference._get_demes_dict(str(graph))
        options = _read_inference_options(options_file)
        param_names, p0, lower_bound, upper_bound = inference._set_up_params_and_bounds(
            options, builder
        )
        constraints = inference._set_up_constraints(options, param_names)
        if perturb > 0:
            p0 = inference._perturb_params_constrained(
                p0, perturb, lower_bound, upper_bound, constraints
            )
        if log:
            p0 = np.log(p0) + 1
            obj_fun = inference._object_func_log
        else:
            obj_fun = inference._object_func
        args = (
            data,
            pts,
            builder,
            options,
            lower_bound,
            upper_bound,
            constraints,
            0,
            uL,
            False,
            None,
        )
        if method == "fmin":
            xopt, fopt, _, _, _ = scipy_optimize.fmin(
                obj_fun,
                p0,
                args=args,
                disp=False,
                maxiter=iterations,
                maxfun=iterations,
                full_output=True,
            )
        elif method == "powell":
            xopt, fopt, _, _, _, _ = scipy_optimize.fmin_powell(
                obj_fun,
                p0,
                args=args,
                disp=False,
                maxiter=iterations,
                maxfun=iterations,
                full_output=True,
            )
        else:
            bounds = list(zip(np.log(lower_bound) + 1, np.log(upper_bound) + 1, strict=False))
            xopt, fopt, _ = scipy_optimize.fmin_l_bfgs_b(
                obj_fun,
                p0,
                bounds=bounds,
                epsilon=1e-2,
                args=args,
                iprint=-1,
                pgtol=1e-5,
                maxiter=iterations,
                maxfun=iterations,
                approx_grad=True,
            )
        if log:
            xopt = np.exp(xopt - 1)
        fitted = None
        if target is not None:
            fitted_builder = inference._update_builder(builder, options, xopt)
            demes = require("demes", _DADI_EXTRA)
            demes.dump(demes.Graph.fromdict(fitted_builder), str(target))
            fitted = str(target)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"dadi could not optimise the demes graph: {exc}") from exc
    return {
        "status": "completed",
        "operation": "optimize_demes_graph_parameters",
        "method": method,
        "log": bool(log),
        "maxiter": iterations,
        "parameter_names": list(param_names),
        "parameters": [float(item) for item in np.asarray(xopt).ravel()],
        "objective": float(fopt),
        "output_path": fitted,
        "provenance": _provenance(
            "dadi",
            "Demes.Inference.optimize steps; inference options read with ruamel.yaml YAML(typ='safe')",
        ),
    }


__all__ = [
    "dfe_mixture_point_pos_spectrum",
    "optimize_demes_graph_parameters",
    "plot_triallelic_spectrum",
    "plot_triallelic_spectrum_comparison",
    "plot_two_locus_spectrum_3d",
]

_register_public_verbose(globals())
