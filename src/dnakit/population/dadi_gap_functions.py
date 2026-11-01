"""dadi wrappers for plots, DFE mixture models, Demes inference and depth-of-coverage distributions.

Plots are drawn by dadi.Plotting (and the TwoLocus, Triallele and DFE plotting modules) with the
matplotlib Agg backend. Each plot is written to a new PNG file, and the function returns that path
together with the spectra or parameters that were drawn, never a plotting object. The DFE mixture
models build dadi's selection caches from a documented demographic model on the fly. Demes
inference reads a demes graph and a dadi inference-options YAML file given by path. The
depth-of-coverage distribution takes per-site depths in memory.

dadi, NumPy and matplotlib are imported only when a function runs. Output files are never overwritten.
"""

from __future__ import annotations

import importlib
import io
import math
from collections.abc import Callable, Mapping, Sequence
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _integer, _number, _provenance
from dnakit.population.dadi_functions import _choice, _nested, _new_file, _numbers, _sizes

_LIKELIHOODS = ("multinom", "poisson")
_RESIDUALS = ("Anscombe", "linear")
_DFE_VARIANTS = ("mixture", "mixture_symmetric_point_pos", "mixture_point_pos", "Vourlaki_mixture")
_DFE_DISTRIBUTIONS = ("biv_lognormal", "biv_ind_gamma")
_DEMES_METHODS = ("fmin", "powell", "lbfgsb")
_MAX_GRID_SIZE = 200
_MAX_CACHE_GAMMA_POINTS = 60
_MAX_CACHE_SAMPLE = 40
_MAX_DEMES_ITERATIONS = 10_000
_MAX_COVERAGE_SITES = 1_000_000
_PLOT_DPI = 100
_OPERATION_ERRORS = (ValueError, TypeError, IndexError, KeyError, AttributeError, ZeroDivisionError)


def _dadi_and_numpy() -> tuple[Any, Any]:
    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    return dadi, np


def _dadi_module(dotted: str) -> Any:
    try:
        return importlib.import_module(dotted)
    except ImportError as exc:
        raise BackendUnavailableError(
            f"{dotted} could not be imported: {exc}",
            hint='Install the optional dependency with: python -m pip install "dnakit[bio-dadi]"',
        ) from exc


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _DADI_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _DADI_EXTRA)


def _fill(value: object, name: str) -> Any:
    if value is None:
        return math.nan
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_fill(item, name) for item in value]
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(
            f"{name} must be a nested list of finite numbers; masked entries are None."
        )
    return float(value)


def _array(np: Any, value: object, name: str, dims: tuple[int, ...], *, signed: bool) -> Any:
    """A nested list as a masked array; None marks an entry that dadi masks."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty nested list; masked entries are None.")
    try:
        raw = np.array(_fill(value, name), dtype=float)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be a rectangular nested list.") from exc
    if raw.ndim not in dims:
        raise ConfigurationError(f"{name} must have {list(dims)} axes, found {raw.ndim}.")
    present = raw[~np.isnan(raw)]
    if not signed and (present < 0).any():
        raise ConfigurationError(f"{name} must contain only non-negative values or None.")
    return np.ma.masked_invalid(raw)


def _pop_ids(value: object, ndim: int) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("pop_ids must be a list of population names.")
    ids = [str(item) for item in value]
    if len(ids) != ndim or len(set(ids)) != ndim or any(not item for item in ids):
        raise ConfigurationError(f"pop_ids must name each of the {ndim} populations once.")
    return ids


def _spectrum(
    dadi: Any,
    np: Any,
    value: object,
    name: str,
    dims: tuple[int, ...] = (1, 2, 3),
    *,
    signed: bool = False,
    pop_ids: object = None,
) -> Any:
    spectrum = dadi.Spectrum(_array(np, value, name, dims, signed=signed))
    if pop_ids is not None:
        spectrum.pop_ids = _pop_ids(pop_ids, spectrum.ndim)
    return spectrum


def _grid_points(value: object) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigurationError("grid_points must list three integer grid sizes.")
    return [_integer(item, "grid point", minimum=8, maximum=_MAX_GRID_SIZE) for item in value]


def _render(image_path: object, operation: str, draw: Callable[[Any], None]) -> Path:
    """Run ``draw`` and save the figure it creates as a PNG at a new path."""

    target = _new_file(image_path, "image_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    plt = _pyplot()
    before = set(plt.get_fignums())
    try:
        try:
            draw(plt)
        except _OPERATION_ERRORS as exc:
            raise ConfigurationError(f"dadi could not draw {operation}: {exc}") from exc
        created = [number for number in plt.get_fignums() if number not in before]
        figure = plt.figure(created[-1]) if created else plt.gcf()
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    finally:
        for number in plt.get_fignums():
            if number not in before:
                plt.close(number)
    return target


def _plot_result(
    operation: str, target: Path, data: dict[str, Any], convention: str
) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "image_path": str(target),
        "image_format": "png",
        **data,
        "provenance": _provenance("dadi", convention + "; drawn with matplotlib (Agg backend)"),
    }


def _likelihood_function(plotting: Any, family: str, likelihood: str) -> Any:
    kind = _choice(likelihood, _LIKELIHOODS, "likelihood")
    name = {"multinom": f"{family}_multinom", "poisson": f"{family}_Poisson"}[kind]
    return getattr(plotting, name)


def _negative_gamma_grid(np: Any, gamma_min: float, gamma_max: float, size: int) -> Any:
    low = _number(gamma_min, "gamma_min", positive=True)
    high = _number(gamma_max, "gamma_max", positive=True)
    if low >= high:
        raise ConfigurationError("gamma_min must be below gamma_max.")
    points = _integer(size, "grid_size", minimum=5, maximum=_MAX_GRID_SIZE)
    return -np.logspace(np.log10(low), np.log10(high), points)


def plot_sfs_1d(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a one-population site frequency spectrum to a PNG, using dadi's ``plot_1d_fs``.

    ``spectrum`` is a 1D list of counts (masked corners as None). The result holds the PNG path
    and the spectrum that was drawn.
    """

    dadi, np = _dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, "spectrum", (1,))
    target = _render(
        image_path,
        "plot_sfs_1d",
        lambda plt: dadi.Plotting.plot_1d_fs(fs, show=False),
    )
    return _plot_result(
        "plot_sfs_1d",
        target,
        {"sample_sizes": [int(n) for n in fs.sample_sizes], "spectrum": _nested(np, fs)},
        "one-population spectrum",
    )


def plot_sfs_1d_comparison(
    model: Sequence[Any],
    data: Sequence[Any],
    image_path: str | Path,
    *,
    likelihood: str = "multinom",
    residual: str = "Anscombe",
    plot_masked: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a model spectrum against a data spectrum with residuals to a PNG, using dadi.

    ``likelihood`` is ``multinom`` (dadi's ``plot_1d_comp_multinom``) or ``poisson``
    (``plot_1d_comp_Poisson``). ``residual`` is ``Anscombe`` (default) or ``linear``.
    """

    dadi, np = _dadi_and_numpy()
    kind = _choice(likelihood, _LIKELIHOODS, "likelihood")
    resid = _choice(residual, _RESIDUALS, "residual")
    model_fs = _spectrum(dadi, np, model, "model", (1,))
    data_fs = _spectrum(dadi, np, data, "data", (1,))
    function = _likelihood_function(dadi.Plotting, "plot_1d_comp", kind)
    target = _render(
        image_path,
        "plot_sfs_1d_comparison",
        lambda plt: function(
            model_fs, data_fs, residual=resid, plot_masked=_flag(plot_masked, "plot_masked"), show=False
        ),
    )
    return _plot_result(
        "plot_sfs_1d_comparison",
        target,
        {
            "likelihood": kind,
            "residual": resid,
            "model": _nested(np, model_fs),
            "data": _nested(np, data_fs),
        },
        f"model and data spectra, {kind} likelihood residuals",
    )


def plot_sfs_2d(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    pop_ids: Sequence[str] | None = None,
    colorbar: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a two-population joint site frequency spectrum as a heat map, using dadi.

    ``spectrum`` is a 2D list of counts (masked corners as None). ``vmin`` and ``vmax`` set the
    colour range; ``pop_ids`` names the two populations on the axes.
    """

    dadi, np = _dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, "spectrum", (2,), pop_ids=pop_ids)
    low, high = _optional_range(vmin, vmax)
    target = _render(
        image_path,
        "plot_sfs_2d",
        lambda plt: dadi.Plotting.plot_single_2d_sfs(
            fs, vmin=low, vmax=high, pop_ids=pop_ids, colorbar=_flag_colorbar(colorbar), show=False
        ),
    )
    return _plot_result(
        "plot_sfs_2d",
        target,
        {"sample_sizes": [int(n) for n in fs.sample_sizes], "spectrum": _nested(np, fs)},
        "two-population spectrum",
    )


def _optional_range(vmin: object, vmax: object) -> tuple[float | None, float | None]:
    low = None if vmin is None else _number(vmin, "vmin")
    high = None if vmax is None else _number(vmax, "vmax")
    if low is not None and high is not None and low >= high:
        raise ConfigurationError("vmin must be below vmax.")
    return low, high


def _flag_colorbar(value: object) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError("colorbar must be true or false.")
    return value


def plot_sfs_2d_residual(
    residual: Sequence[Any],
    image_path: str | Path,
    *,
    resid_range: float | None = None,
    pop_ids: Sequence[str] | None = None,
    colorbar: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a two-population residual spectrum (signed values) as a heat map, using dadi.

    ``residual`` is a 2D nested list of finite numbers, masked entries as None. ``resid_range``
    saturates the colour scale at plus and minus that value.
    """

    dadi, np = _dadi_and_numpy()
    array = _spectrum(dadi, np, residual, "residual", (2,), signed=True, pop_ids=pop_ids)
    bound = None if resid_range is None else _number(resid_range, "resid_range", positive=True)
    target = _render(
        image_path,
        "plot_sfs_2d_residual",
        lambda plt: dadi.Plotting.plot_2d_resid(
            array, resid_range=bound, pop_ids=pop_ids, colorbar=_flag_colorbar(colorbar)
        ),
    )
    return _plot_result(
        "plot_sfs_2d_residual",
        target,
        {"residual": _nested(np, array), "resid_range": bound},
        "two-population residuals",
    )


def plot_sfs_2d_comparison(
    model: Sequence[Any],
    data: Sequence[Any],
    image_path: str | Path,
    *,
    likelihood: str = "multinom",
    vmin: float | None = None,
    vmax: float | None = None,
    resid_range: float | None = None,
    pop_ids: Sequence[str] | None = None,
    residual: str = "Anscombe",
    adjust: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw model, data and residual two-population spectra side by side, using dadi.

    ``likelihood`` picks dadi's ``plot_2d_comp_multinom`` (default) or ``plot_2d_comp_Poisson``.
    ``residual`` is ``Anscombe`` (default) or ``linear``.
    """

    dadi, np = _dadi_and_numpy()
    kind = _choice(likelihood, _LIKELIHOODS, "likelihood")
    resid = _choice(residual, _RESIDUALS, "residual")
    model_fs = _spectrum(dadi, np, model, "model", (2,), pop_ids=pop_ids)
    data_fs = _spectrum(dadi, np, data, "data", (2,), pop_ids=pop_ids)
    low, high = _optional_range(vmin, vmax)
    bound = None if resid_range is None else _number(resid_range, "resid_range", positive=True)
    function = _likelihood_function(dadi.Plotting, "plot_2d_comp", kind)
    target = _render(
        image_path,
        "plot_sfs_2d_comparison",
        lambda plt: function(
            model_fs,
            data_fs,
            vmin=low,
            vmax=high,
            resid_range=bound,
            pop_ids=pop_ids,
            residual=resid,
            adjust=_flag_colorbar(adjust),
            show=False,
        ),
    )
    return _plot_result(
        "plot_sfs_2d_comparison",
        target,
        {
            "likelihood": kind,
            "residual": resid,
            "model": _nested(np, model_fs),
            "data": _nested(np, data_fs),
        },
        f"two-population model and data spectra, {kind} residuals",
    )


def plot_sfs_2d_meta_residual(
    synonymous_residual: Sequence[Any],
    nonsynonymous_residual: Sequence[Any],
    image_path: str | Path,
    *,
    resid_range: float | None = None,
    pop_ids: Sequence[str] | None = None,
    adjust: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the meta-residual comparison of two two-population residual spectra, using dadi.

    The first residual is the synonymous one and the second the nonsynonymous one, as in dadi's
    ``plot_2d_meta_resid``. Both must be 2D residual lists with masked entries as None.
    """

    dadi, np = _dadi_and_numpy()
    syn = _spectrum(dadi, np, synonymous_residual, "synonymous_residual", (2,), signed=True)
    nonsyn = _spectrum(
        dadi, np, nonsynonymous_residual, "nonsynonymous_residual", (2,), signed=True
    )
    if pop_ids is not None:
        syn.pop_ids = _pop_ids(pop_ids, 2)
        nonsyn.pop_ids = _pop_ids(pop_ids, 2)
    bound = None if resid_range is None else _number(resid_range, "resid_range", positive=True)
    target = _render(
        image_path,
        "plot_sfs_2d_meta_residual",
        lambda plt: dadi.Plotting.plot_2d_meta_resid(
            syn,
            nonsyn,
            resid_range=bound,
            pop_ids=pop_ids,
            adjust=_flag_colorbar(adjust),
            show=False,
        ),
    )
    return _plot_result(
        "plot_sfs_2d_meta_residual",
        target,
        {
            "synonymous_residual": _nested(np, syn),
            "nonsynonymous_residual": _nested(np, nonsyn),
            "resid_range": bound,
        },
        "synonymous and nonsynonymous two-population residuals",
    )


def plot_sfs_3d(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    pop_ids: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a three-population spectrum as a logarithmic heat map of slices, using dadi.

    ``spectrum`` is a 3D list of counts (masked entries as None). The colour range is
    ``vmin`` to ``vmax`` when given.
    """

    dadi, np = _dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, "spectrum", (3,), pop_ids=pop_ids)
    low, high = _optional_range(vmin, vmax)
    target = _render(
        image_path,
        "plot_sfs_3d",
        lambda plt: dadi.Plotting.plot_3d_spectrum(
            fs, vmin=low, vmax=high, pop_ids=pop_ids, show=False
        ),
    )
    return _plot_result(
        "plot_sfs_3d",
        target,
        {"sample_sizes": [int(n) for n in fs.sample_sizes], "spectrum": _nested(np, fs)},
        "three-population spectrum",
    )


def plot_sfs_3d_comparison(
    model: Sequence[Any],
    data: Sequence[Any],
    image_path: str | Path,
    *,
    likelihood: str = "multinom",
    vmin: float | None = None,
    vmax: float | None = None,
    resid_range: float | None = None,
    pop_ids: Sequence[str] | None = None,
    residual: str = "Anscombe",
    adjust: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw model, data and residual three-population spectra as slice panels, using dadi.

    ``likelihood`` picks dadi's ``plot_3d_comp_multinom`` (default) or ``plot_3d_comp_Poisson``.
    ``residual`` is ``Anscombe`` (default) or ``linear``.
    """

    dadi, np = _dadi_and_numpy()
    kind = _choice(likelihood, _LIKELIHOODS, "likelihood")
    resid = _choice(residual, _RESIDUALS, "residual")
    model_fs = _spectrum(dadi, np, model, "model", (3,), pop_ids=pop_ids)
    data_fs = _spectrum(dadi, np, data, "data", (3,), pop_ids=pop_ids)
    low, high = _optional_range(vmin, vmax)
    bound = None if resid_range is None else _number(resid_range, "resid_range", positive=True)
    function = _likelihood_function(dadi.Plotting, "plot_3d_comp", kind)
    target = _render(
        image_path,
        "plot_sfs_3d_comparison",
        lambda plt: function(
            model_fs,
            data_fs,
            vmin=low,
            vmax=high,
            resid_range=bound,
            pop_ids=pop_ids,
            residual=resid,
            adjust=_flag_colorbar(adjust),
            show=False,
        ),
    )
    return _plot_result(
        "plot_sfs_3d_comparison",
        target,
        {
            "likelihood": kind,
            "residual": resid,
            "model": _nested(np, model_fs),
            "data": _nested(np, data_fs),
        },
        f"three-population model and data spectra, {kind} residuals",
    )


def plot_sfs_3d_pairwise(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    pop_ids: Sequence[str] | None = None,
    adjust: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the three pairwise two-population marginals of a three-population spectrum, using dadi.

    ``spectrum`` is a 3D list of counts (masked entries as None). Each pair of populations is
    shown as a marginal two-population spectrum.
    """

    dadi, np = _dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, "spectrum", (3,), pop_ids=pop_ids)
    low, high = _optional_range(vmin, vmax)
    target = _render(
        image_path,
        "plot_sfs_3d_pairwise",
        lambda plt: dadi.Plotting.plot_3d_pairwise(
            fs, vmin=low, vmax=high, pop_ids=pop_ids, adjust=_flag_colorbar(adjust), show=False
        ),
    )
    return _plot_result(
        "plot_sfs_3d_pairwise",
        target,
        {"sample_sizes": [int(n) for n in fs.sample_sizes], "spectrum": _nested(np, fs)},
        "three-population spectrum, pairwise marginals",
    )


def plot_trispectrum(
    spectrum: Sequence[Any],
    image_path: str | Path,
    *,
    folded: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a triallelic (three-allele) spectrum as a heat map, using dadi's Triallele plotting.

    ``spectrum`` is a 2D list of counts with masked entries as None. ``folded=True`` folds the
    spectrum before drawing, as dadi's ``plot_single_trispectrum`` does.
    """

    dadi, np = _dadi_and_numpy()
    tri = _dadi_module("dadi.Triallele.TriSpectrum_mod").TriSpectrum
    plotting = _dadi_module("dadi.Triallele.plotting")
    array = _array(np, spectrum, "spectrum", (2,), signed=False)
    sfs = tri(array)
    low, high = _optional_range(vmin, vmax)
    fold = _flag(folded, "folded")
    target = _render(
        image_path,
        "plot_trispectrum",
        lambda plt: plotting.plot_single_trispectrum(
            sfs, folded=fold, vmin=low, vmax=high, colorbar=_flag_colorbar(colorbar)
        ),
    )
    return _plot_result(
        "plot_trispectrum",
        target,
        {"folded": fold, "spectrum": _nested(np, sfs)},
        "triallelic spectrum",
    )


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def plot_trispectrum_comparison(
    first: Sequence[Any],
    second: Sequence[Any],
    image_path: str | Path,
    *,
    folded: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    resid_range: float | None = None,
    title1: str = "sfs1",
    title2: str = "sfs2",
    title3: str = "Residual",
    colorbar: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw two triallelic spectra and their residual side by side, using dadi's Triallele plotting.

    ``first`` and ``second`` are 2D lists of counts with masked entries as None. Panel titles are
    ``title1``, ``title2`` and ``title3``.
    """

    dadi, np = _dadi_and_numpy()
    tri = _dadi_module("dadi.Triallele.TriSpectrum_mod").TriSpectrum
    plotting = _dadi_module("dadi.Triallele.plotting")
    sfs1 = tri(_array(np, first, "first", (2,), signed=False))
    sfs2 = tri(_array(np, second, "second", (2,), signed=False))
    low, high = _optional_range(vmin, vmax)
    bound = None if resid_range is None else _number(resid_range, "resid_range", positive=True)
    fold = _flag(folded, "folded")
    target = _render(
        image_path,
        "plot_trispectrum_comparison",
        lambda plt: plotting.plot_trispectrum_comp(
            sfs1,
            sfs2,
            folded=fold,
            vmin=low,
            vmax=high,
            resid_range=bound,
            colorbar=_flag_colorbar(colorbar),
            title1=_title(title1),
            title2=_title(title2),
            title3=_title(title3),
        ),
    )
    return _plot_result(
        "plot_trispectrum_comparison",
        target,
        {"folded": fold, "first": _nested(np, sfs1), "second": _nested(np, sfs2)},
        "two triallelic spectra and their residual",
    )


def _title(value: object) -> str:
    if not isinstance(value, str):
        raise ConfigurationError("panel titles must be strings.")
    return value


def plot_two_locus_spectrum(
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
    """Draw a two-locus (three-dimensional) spectrum as a 3D scatter plot, using dadi's TwoLocus.

    ``spectrum`` is a 3D list of counts with masked entries as None. ``max_sum`` slices the
    spectrum to points closer to the origin than that diagonal plane; ``max_index`` keeps only
    entries that many steps from an axis plane.
    """

    dadi, np = _dadi_and_numpy()
    tl = _dadi_module("dadi.TwoLocus.TLSpectrum_mod").TLSpectrum
    plotting = _dadi_module("dadi.TwoLocus.plotting")
    fs = tl(_array(np, spectrum, "spectrum", (3,), signed=False))
    low, high = _optional_range(vmin, vmax)
    cut = _optional_int(max_sum, "max_sum")
    index = _optional_int(max_index, "max_index")
    target = _render(
        image_path,
        "plot_two_locus_spectrum",
        lambda plt: plotting.plot_3d(
            fs,
            vmin=low,
            vmax=high,
            max_sum=cut,
            max_index=index,
            show=False,
            colorbar=_flag_colorbar(colorbar),
        ),
    )
    return _plot_result(
        "plot_two_locus_spectrum",
        target,
        {"spectrum": _nested(np, fs)},
        "two-locus spectrum",
    )


def _optional_int(value: object, name: str) -> int | None:
    return None if value is None else _integer(value, name, minimum=0, maximum=_MAX_GRID_SIZE)


def plot_two_locus_comparison(
    model: Sequence[Any],
    data: Sequence[Any],
    image_path: str | Path,
    *,
    resid_range: float = 1.0,
    max_sum: int | None = None,
    max_index: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a two-locus model spectrum against data as 3D panels with residuals, using dadi.

    ``resid_range`` saturates the residual colour scale. ``max_sum`` and ``max_index`` slice the
    spectra as in ``plot_two_locus_spectrum``.
    """

    dadi, np = _dadi_and_numpy()
    tl = _dadi_module("dadi.TwoLocus.TLSpectrum_mod").TLSpectrum
    plotting = _dadi_module("dadi.TwoLocus.plotting")
    model_fs = tl(_array(np, model, "model", (3,), signed=False))
    data_fs = tl(_array(np, data, "data", (3,), signed=False))
    bound = _number(resid_range, "resid_range", positive=True)
    cut = _optional_int(max_sum, "max_sum")
    index = _optional_int(max_index, "max_index")
    target = _render(
        image_path,
        "plot_two_locus_comparison",
        lambda plt: plotting.plot_3d_comp(
            model_fs, data_fs, resid_range=bound, max_sum=cut, max_index=index, show=False
        ),
    )
    return _plot_result(
        "plot_two_locus_comparison",
        target,
        {"resid_range": bound, "model": _nested(np, model_fs), "data": _nested(np, data_fs)},
        "two-locus model and data spectra",
    )


def plot_dfe_bivariate(
    distribution: str,
    params: Sequence[float],
    image_path: str | Path,
    *,
    gamma_min: float = 1e-3,
    gamma_max: float = 100.0,
    grid_size: int = 50,
    colorbar: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a bivariate distribution of fitness effects (negative gammas) as a heat map, using dadi.

    ``distribution`` is ``biv_lognormal`` (``params`` = mu, sigma, rho in log scale) or
    ``biv_ind_gamma`` (``params`` = alpha, beta, rho ignored). The gamma grid is logarithmic
    from ``gamma_min`` to ``gamma_max`` on both axes, drawn as negative values as dadi expects.
    """

    dadi, np = _dadi_and_numpy()
    name = _choice(distribution, _DFE_DISTRIBUTIONS, "distribution")
    values = _numbers(params, "params")
    pdfs = _dadi_module("dadi.DFE.PDFs")
    grid = _negative_gamma_grid(np, gamma_min, gamma_max, grid_size)
    plotting = _dadi_module("dadi.DFE.Plotting")
    target = _render(
        image_path,
        "plot_dfe_bivariate",
        lambda plt: plotting.plot_biv_dfe(
            grid, grid, getattr(pdfs, name), values, colorbar=_flag_colorbar(colorbar)
        ),
    )
    return _plot_result(
        "plot_dfe_bivariate",
        target,
        {
            "distribution": name,
            "params": values,
            "gamma_grid": [float(item) for item in grid],
        },
        f"bivariate DFE, {name}, negative gamma grid",
    )


def plot_dfe_bivariate_point_pos(
    distribution: str,
    params: Sequence[float],
    image_path: str | Path,
    *,
    rho: float = 0.0,
    gamma_min: float = 1e-3,
    gamma_max: float = 100.0,
    grid_size: int = 50,
    colorbar: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a bivariate DFE with a positive point mass as a heat map, using dadi.

    ``params`` holds the distribution parameters followed by the last four values
    (proportion positive in population 1, positive gamma in population 1, proportion positive in
    population 2, positive gamma in population 2), as dadi's ``plot_biv_point_pos_dfe`` expects.
    ``rho`` is the correlation passed to the distribution. ``distribution`` is ``biv_lognormal``
    or ``biv_ind_gamma``.
    """

    dadi, np = _dadi_and_numpy()
    name = _choice(distribution, _DFE_DISTRIBUTIONS, "distribution")
    values = _numbers(params, "params", minimum=5)
    correlation = _number(rho, "rho", minimum=-1.0)
    pdfs = _dadi_module("dadi.DFE.PDFs")
    grid = _negative_gamma_grid(np, gamma_min, gamma_max, grid_size)
    plotting = _dadi_module("dadi.DFE.Plotting")
    target = _render(
        image_path,
        "plot_dfe_bivariate_point_pos",
        lambda plt: plotting.plot_biv_point_pos_dfe(
            grid,
            grid,
            getattr(pdfs, name),
            values,
            rho=correlation,
            colorbar=_flag_colorbar(colorbar),
        ),
    )
    return _plot_result(
        "plot_dfe_bivariate_point_pos",
        target,
        {
            "distribution": name,
            "params": values,
            "rho": correlation,
            "gamma_grid": [float(item) for item in grid],
        },
        f"bivariate DFE with positive point mass, {name}, negative gamma grid",
    )


def dfe_mixture_spectrum(
    variant: str,
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
    """Expected two-population spectrum for a DFE mixture model, using dadi's selection caches.

    The demography is dadi's ``split_mig_sel_single_gamma`` (one-dimensional cache) and
    ``split_mig_sel`` (two-dimensional cache) with ``demographic_params`` = nu1, nu2, T, m. The
    caches are built from ``grid_points`` (three grid sizes), ``gamma_points`` log-spaced selection
    coefficients up to ``gamma_max``, and the point-mass gammas of the chosen variant.

    ``variant`` picks the mixture: ``mixture`` (lognormal 1D plus bivariate lognormal 2D,
    ``dfe_params`` = mu, sigma, rho, weight_2d), ``mixture_symmetric_point_pos`` (adds
    rho-, proportion-, gamma-, weight; 6 values), ``mixture_point_pos`` (8 values, separate
    positive proportion and gamma per population), or ``Vourlaki_mixture`` (gamma and independent
    bivariate gamma: alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos).
    """

    dadi, np = _dadi_and_numpy()
    name = _choice(variant, _DFE_VARIANTS, "variant")
    demo = _numbers(demographic_params, "demographic_params")
    if len(demo) != 4:
        raise ConfigurationError("demographic_params must be nu1, nu2, T, m (four values).")
    ns = _sizes(sample_sizes, 2, _MAX_CACHE_SAMPLE, "sample_sizes")
    pts = _grid_points(grid_points)
    gammas = _integer(gamma_points, "gamma_points", minimum=5, maximum=_MAX_CACHE_GAMMA_POINTS)
    ceiling = _number(gamma_max, "gamma_max", positive=True)
    values = _numbers(dfe_params, "dfe_params")
    scale = _number(theta, "theta", positive=True)
    expected, point_gammas = _dfe_layout(name, values)
    if len(values) != expected:
        raise ConfigurationError(f"{name} needs {expected} dfe_params, got {len(values)}.")
    selection = _dadi_module("dadi.DFE.DemogSelModels")
    mixtures = _dadi_module("dadi.DFE.Cache2D_mod")
    cache1d = _dadi_module("dadi.DFE.Cache1D_mod").Cache1D
    cache2d = _dadi_module("dadi.DFE.Cache2D_mod").Cache2D
    pdfs = _dadi_module("dadi.DFE.PDFs")
    if name == "Vourlaki_mixture":
        dist1, dist2 = pdfs.gamma, pdfs.biv_ind_gamma
    else:
        dist1, dist2 = pdfs.lognormal, pdfs.biv_lognormal
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
        if name == "Vourlaki_mixture":
            spectrum = _dadi_module("dadi.DFE.Vourlaki2022").Vourlaki_mixture(
                values, None, s1, s2, scale, None
            )
        else:
            function = getattr(mixtures, name)
            spectrum = function(values, None, s1, s2, dist1, dist2, scale, None)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"dadi could not compute {name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "dfe_mixture_spectrum",
        "variant": name,
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
            "expected spectrum of a DFE mixture; split_mig_sel demography with selection caches",
        ),
    }


def _dfe_layout(variant: str, values: list[float]) -> tuple[int, list[float]]:
    """Number of parameters the variant takes, and the positive point-mass gammas it needs."""

    if variant == "mixture":
        return 4, []
    if variant == "mixture_symmetric_point_pos":
        return 6, [values[-2]] if len(values) >= 2 else []
    if variant == "mixture_point_pos":
        return 8, [values[-4], values[-2]] if len(values) >= 4 else []
    return 6, [values[3]] if len(values) >= 4 else []


def fit_demes_graph_parameters(
    graph: str | Path,
    inference_options: str | Path,
    sfs: Sequence[Any],
    pop_ids: Sequence[str],
    grid_points: Sequence[int],
    *,
    method: str = "fmin",
    maxiter: int = 1000,
    log: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit the free parameters of a demes graph to a frequency spectrum, using dadi's Demes inference.

    ``graph`` is a demes YAML file and ``inference_options`` a dadi inference-options YAML file
    (a ``parameters`` list naming each free value, its location in the graph, and optional
    bounds). ``sfs`` is the observed spectrum as a nested list of counts (masked entries as None),
    and ``pop_ids`` names its populations, which must match demes in the graph. ``grid_points`` is
    three integration grid sizes. The optimisation is dadi's ``Demes.Inference.optimize`` with a
    Poisson-free multinomial likelihood (no ``uL``), run with the given ``method``
    (``fmin``, ``powell`` or ``lbfgsb``) and ``maxiter``.
    """

    dadi, np = _dadi_and_numpy()
    inference = _dadi_module("dadi.Demes.Inference")
    graph_path = input_path(graph)
    options_path = input_path(inference_options)
    spectrum = _spectrum(dadi, np, sfs, "sfs", (1, 2, 3), pop_ids=pop_ids)
    pts = _grid_points(grid_points)
    chosen = _choice(method, _DEMES_METHODS, "method")
    iterations = _integer(maxiter, "maxiter", minimum=1, maximum=_MAX_DEMES_ITERATIONS)
    use_log = _flag(log, "log")
    try:
        names, xopt, fopt = inference.optimize(
            str(graph_path),
            str(options_path),
            spectrum,
            pts,
            maxiter=iterations,
            verbose=0,
            log=use_log,
            method=chosen,
            output_stream=io.StringIO(),
            output=None,
        )
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"dadi Demes inference failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "fit_demes_graph_parameters",
        "graph": str(graph_path),
        "inference_options": str(options_path),
        "pop_ids": list(pop_ids),
        "sample_sizes": [int(n) for n in spectrum.sample_sizes],
        "grid_points": pts,
        "method": chosen,
        "maxiter": iterations,
        "log": use_log,
        "parameter_names": [str(name) for name in names],
        "params": [float(value) for value in np.atleast_1d(xopt)],
        "fopt": _finite_or_none(fopt),
        "provenance": _provenance(
            "dadi", "demes graph fitted by dadi Demes.Inference.optimize; fopt as dadi returns it"
        ),
    }


def _finite_or_none(value: object) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None


def coverage_depth_distribution(
    site_depths: Mapping[str, Sequence[Sequence[int]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Depth-of-coverage distribution for each population, using dadi's LowPass ``compute_cov_dist``.

    ``site_depths`` maps each population ID to a list over sites, and each site gives the
    sequencing depths of the individuals in that population. Every population must list the same
    number of sites. The result gives, per population, the depths 0..max and their normalised
    frequencies, as dadi returns them.
    """

    np = require("numpy", _DADI_EXTRA)
    compute = _dadi_module("dadi.LowPass.LowPass").compute_cov_dist
    if not isinstance(site_depths, Mapping) or not site_depths:
        raise ConfigurationError("site_depths must map at least one population ID to its sites.")
    pops = [str(pop) for pop in site_depths]
    if len(set(pops)) != len(pops) or any(not pop for pop in pops):
        raise ConfigurationError("population IDs must be unique, non-empty strings.")
    per_pop = [site_depths[pop] for pop in site_depths]
    counts = {len(sites) for sites in per_pop if not isinstance(sites, (str, bytes))}
    if len(counts) != 1 or not per_pop[0]:
        raise ConfigurationError("every population must list the same, non-zero number of sites.")
    sites_n = counts.pop()
    if sites_n > _MAX_COVERAGE_SITES:
        raise ConfigurationError(f"at most {_MAX_COVERAGE_SITES} sites are accepted.")
    data: dict[str, dict[str, Any]] = {}
    for index in range(sites_n):
        coverage = {}
        for pop, sites in zip(pops, per_pop, strict=True):
            depths = sites[index]
            if isinstance(depths, (str, bytes)) or not isinstance(depths, Sequence) or not depths:
                raise ConfigurationError("each site must list at least one depth per population.")
            coverage[pop] = np.array(
                [_integer(depth, "depth", minimum=0) for depth in depths], dtype=int
            )
        data[str(index)] = {"coverage": coverage}
    try:
        distribution = compute(data, pops)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"dadi could not compute the coverage distribution: {exc}") from exc
    result = {}
    for pop in pops:
        depth_row, frequency_row = distribution[pop][0], distribution[pop][1]
        result[pop] = {
            "depth": [int(value) for value in depth_row],
            "frequency": [float(value) for value in frequency_row],
        }
    return {
        "status": "completed",
        "operation": "coverage_depth_distribution",
        "sites": sites_n,
        "distribution": result,
        "provenance": _provenance("dadi", "depth-of-coverage distribution, LowPass compute_cov_dist"),
    }


__all__ = [
    "coverage_depth_distribution",
    "dfe_mixture_spectrum",
    "fit_demes_graph_parameters",
    "plot_dfe_bivariate",
    "plot_dfe_bivariate_point_pos",
    "plot_sfs_1d",
    "plot_sfs_1d_comparison",
    "plot_sfs_2d",
    "plot_sfs_2d_comparison",
    "plot_sfs_2d_meta_residual",
    "plot_sfs_2d_residual",
    "plot_sfs_3d",
    "plot_sfs_3d_comparison",
    "plot_sfs_3d_pairwise",
    "plot_trispectrum",
    "plot_trispectrum_comparison",
    "plot_two_locus_comparison",
    "plot_two_locus_spectrum",
]

_register_public_verbose(globals(), tuple(__all__))
