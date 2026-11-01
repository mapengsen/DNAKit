"""scikit-allel wrappers for haplotype, site-frequency and distance plots, and for ROH calls.

The plot functions call scikit-allel's own drawing functions with matplotlib's Agg backend. Each
one writes a new PNG file (150 dpi) at ``output_path``, never overwrites an existing file, and
returns that path together with the numbers that were drawn, never a plotting object. The ROH
wrapper calls scikit-allel's Poisson HMM caller, which needs the optional ``protopunica`` package.
The VCF-to-Zarr conversion and the HTML repr helpers of scikit-allel are not wrapped (see the
manifest notes).

scikit-allel, NumPy, matplotlib and seaborn are imported only when a function runs. No function
uses the network.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_ALLEL_EXTRA = "bio-scikit-allel"
_PLOT_EXTRA = "bio-plots"
_PLOT_DPI = 150
_OPERATION_ERRORS = (ValueError, TypeError, IndexError, KeyError, NotImplementedError)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _PLOT_EXTRA)


def _allel_and_numpy() -> tuple[Any, Any]:
    allel = require("allel", _ALLEL_EXTRA)
    np = require("numpy", _ALLEL_EXTRA)
    return allel, np


def _numbers(
    np: Any,
    value: object,
    name: str,
    dims: tuple[int, ...],
    *,
    integer: bool = False,
    minimum: float | None = None,
) -> Any:
    """A list (or nested list) of finite numbers as a NumPy array with one of the allowed ndims."""

    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular list of numbers.") from exc
    if array.ndim not in dims or array.size == 0:
        axes = " or ".join(str(dim) for dim in dims)
        raise ConfigurationError(f"{name} must be a non-empty list with {axes} axes.")
    if not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must contain only finite numbers.")
    if minimum is not None and (array < minimum).any():
        raise ConfigurationError(f"{name} must not contain values smaller than {minimum}.")
    if integer:
        if (array != np.floor(array)).any():
            raise ConfigurationError(f"{name} must contain whole numbers.")
        return array.astype(np.int64)
    return array


def _haplotype_matrix(np: Any, value: object, name: str) -> Any:
    """Allele indices, variants x haplotypes (0 and 1 for biallelic sites)."""

    return _numbers(np, value, name, (2,), integer=True, minimum=0)


def _count(value: object, name: str, *, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    return value


def _index(value: object, size: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < size:
        raise ConfigurationError(f"{name} must be an integer from 0 to {size - 1}.")
    return value


def _new_png(output_path: object) -> Path:
    """A new PNG path; the parent directory is created and existing files are never overwritten."""

    if not isinstance(output_path, (str, Path)) or not str(output_path).strip():
        raise ConfigurationError("output_path must be a non-empty file path.")
    target = Path(output_path).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    if target.exists():
        raise ConfigurationError(f"output_path already exists and is not overwritten: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _render(output_path: object, operation: str, draw: Callable[[], Any]) -> tuple[Path, Any]:
    """Run ``draw`` (scikit-allel returns an Axes or a Figure), save it, close new figures."""

    target = _new_png(output_path)
    plt = _pyplot()
    before = set(plt.get_fignums())
    try:
        try:
            drawn = draw()
        except _OPERATION_ERRORS as exc:
            raise ConfigurationError(f"scikit-allel could not draw {operation}: {exc}") from exc
        figure = drawn if hasattr(drawn, "savefig") else drawn.figure
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    finally:
        for number in plt.get_fignums():
            if number not in before:
                plt.close(number)
    return target, drawn


def _plain(value: Any) -> Any:
    """Convert NumPy values to JSON-serialisable Python values; non-finite floats become None."""

    if hasattr(value, "tolist") and not isinstance(value, (list, tuple, dict)):
        value = value.tolist()
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _result(target: Path, data: dict[str, Any]) -> dict[str, Any]:
    return {"path": str(target), "data": _plain(data)}


def roh_poissonhmm(
    genotypes: list[list[int]],
    positions: list[int],
    *,
    phet_roh: float = 0.001,
    phet_nonroh_low: float = 0.0025,
    phet_nonroh_high: float = 0.01,
    transition: float = 0.001,
    window_size: int = 1000,
    min_roh: int = 0,
    is_accessible: list[bool] | None = None,
    contig_size: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call runs of homozygosity (ROH) in one individual with scikit-allel's Poisson HMM.

    ``genotypes`` is a variants x ploidy list of allele indices (for a diploid, two columns) and
    ``positions`` the matching 1-based variant positions. The chromosome is cut into accessible
    windows of ``window_size`` bases, the heterozygous calls per window feed a Poisson hidden Markov
    model with the heterozygote rate ``phet_roh`` inside a ROH and the two rates
    ``phet_nonroh_low`` and ``phet_nonroh_high`` outside it, and the ROH blocks are returned with
    the fraction of the contig they cover (FROH). This needs the optional ``protopunica`` package;
    without it a BackendUnavailableError with an install hint is raised.
    """

    allel, np = _allel_and_numpy()
    genotype_array = _numbers(np, genotypes, "genotypes", (2,), integer=True, minimum=0)
    position_array = _numbers(np, positions, "positions", (1,), integer=True, minimum=0)
    if position_array.shape[0] != genotype_array.shape[0]:
        raise ConfigurationError("positions must have one entry per genotype row.")
    accessible = None
    if is_accessible is not None:
        # One flag per base of the contig; scikit-allel then takes the contig size from its length.
        accessible = _numbers(np, is_accessible, "is_accessible", (1,)).astype(bool)
    if contig_size is not None:
        contig_size = _count(contig_size, "contig_size")
    if accessible is None and contig_size is None:
        raise ConfigurationError("Provide is_accessible (one flag per base) or contig_size.")
    window = _count(window_size, "window_size")
    _count(min_roh, "min_roh", minimum=0)
    require("protopunica", _ALLEL_EXTRA)
    blocks, froh = allel.roh_poissonhmm(
        genotype_array,
        position_array,
        phet_roh=phet_roh,
        phet_nonroh=(phet_nonroh_low, phet_nonroh_high),
        transition=transition,
        window_size=window,
        min_roh=min_roh,
        is_accessible=accessible,
        contig_size=contig_size,
    )
    return {
        "status": "completed",
        "roh_blocks": _plain(blocks.to_dict(orient="records")),
        "block_count": len(blocks),
        "froh": float(froh),
        "window_size": window,
        "phet_roh": phet_roh,
        "phet_nonroh": [phet_nonroh_low, phet_nonroh_high],
        "transition": transition,
        "min_roh": min_roh,
    }


def _matrix_size(values: Any) -> int:
    """Samples (or variants) in a square matrix, or in its condensed upper-triangle form."""

    if values.ndim == 2:
        if values.shape[0] != values.shape[1]:
            raise ConfigurationError(
                "A square matrix must have the same number of rows and columns."
            )
        return values.shape[0]
    length = int(values.shape[0])
    size = round((1 + math.sqrt(1 + 8 * length)) / 2)
    if size * (size - 1) // 2 != length:
        raise ConfigurationError("A condensed matrix must hold n*(n-1)/2 values for some n.")
    return size


def _sfs_plot(
    upstream: str,
    spectrum: object,
    output_path: object,
    *,
    yscale: str,
    bins: int | list[float] | None,
    n: int | None,
    clip_endpoints: bool,
    label: str | None,
    plot_kwargs: dict[str, Any] | None,
) -> dict[str, Any]:
    allel, np = _allel_and_numpy()
    values = _numbers(np, spectrum, "spectrum", (1,), minimum=0)
    if not isinstance(yscale, str) or not yscale:
        raise ConfigurationError(
            "yscale must be a matplotlib scale name such as 'log' or 'linear'."
        )
    bin_value: Any = None
    if isinstance(bins, list):
        bin_value = _numbers(np, bins, "bins", (1,))
    elif bins is not None:
        bin_value = _count(bins, "bins")
    if n is not None:
        n = _count(n, "n")
    kwargs = dict(plot_kwargs or {})
    target, ax = _render(
        output_path,
        upstream,
        lambda: getattr(allel, upstream)(
            values,
            yscale=yscale,
            bins=bin_value,
            n=n,
            clip_endpoints=clip_endpoints,
            label=label,
            plot_kwargs=kwargs,
        ),
    )
    line = ax.lines[0]
    return _result(
        target,
        {
            "spectrum": values,
            "x": line.get_xdata(),
            "y": line.get_ydata(),
            "xlabel": ax.get_xlabel(),
            "ylabel": ax.get_ylabel(),
            "yscale": ax.get_yscale(),
            "n": n,
            "bins": bin_value,
            "clip_endpoints": clip_endpoints,
            "label": label,
        },
    )


def _joint_plot(
    upstream: str,
    spectrum: object,
    output_path: object,
    *,
    imshow_kwargs: dict[str, Any] | None,
) -> dict[str, Any]:
    allel, np = _allel_and_numpy()
    values = _numbers(np, spectrum, "spectrum", (2,), minimum=0)
    kwargs = dict(imshow_kwargs or {})
    target, ax = _render(
        output_path, upstream, lambda: getattr(allel, upstream)(values, imshow_kwargs=kwargs)
    )
    return _result(
        target,
        {
            "spectrum": values,
            "shape": list(values.shape),
            "drawn_matrix": ax.images[0].get_array(),
            "xlabel": ax.get_xlabel(),
            "ylabel": ax.get_ylabel(),
        },
    )


def plot_sfs(
    spectrum: list[float],
    output_path: str,
    *,
    yscale: str = "log",
    bins: int | list[float] | None = None,
    n: int | None = None,
    clip_endpoints: bool = True,
    label: str | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a one-dimensional site frequency spectrum with scikit-allel's plot_sfs.

    ``spectrum`` holds the site counts by derived allele count (index 0 to the number of
    chromosomes). ``bins`` groups the allele counts (an int or bin edges), ``n`` (chromosomes
    sampled) switches the x axis to derived allele frequency, and ``clip_endpoints`` drops the
    first and last entries as scikit-allel does by default. Returns the PNG path and the points
    drawn (x, y) with the axis labels.
    """

    return _sfs_plot(
        "plot_sfs",
        spectrum,
        output_path,
        yscale=yscale,
        bins=bins,
        n=n,
        clip_endpoints=clip_endpoints,
        label=label,
        plot_kwargs=plot_kwargs,
    )


def plot_sfs_folded(
    spectrum: list[float],
    output_path: str,
    *,
    yscale: str = "log",
    bins: int | list[float] | None = None,
    n: int | None = None,
    clip_endpoints: bool = True,
    label: str | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a folded site frequency spectrum (minor allele counts) with plot_sfs_folded.

    The arguments and the returned data are those of plot_sfs; scikit-allel labels the x axis
    'minor allele count'.
    """

    return _sfs_plot(
        "plot_sfs_folded",
        spectrum,
        output_path,
        yscale=yscale,
        bins=bins,
        n=n,
        clip_endpoints=clip_endpoints,
        label=label,
        plot_kwargs=plot_kwargs,
    )


def plot_sfs_scaled(
    spectrum: list[float],
    output_path: str,
    *,
    yscale: str = "log",
    bins: int | list[float] | None = None,
    n: int | None = None,
    clip_endpoints: bool = True,
    label: str | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a scaled site frequency spectrum with scikit-allel's plot_sfs_scaled.

    scikit-allel draws the values it is given and only changes the y-axis label to 'scaled site
    frequency'. To plot scaled values, pass the output of scikit-allel's sfs_scaled. The result
    holds the values as drawn.
    """

    return _sfs_plot(
        "plot_sfs_scaled",
        spectrum,
        output_path,
        yscale=yscale,
        bins=bins,
        n=n,
        clip_endpoints=clip_endpoints,
        label=label,
        plot_kwargs=plot_kwargs,
    )


def plot_sfs_folded_scaled(
    spectrum: list[float],
    output_path: str,
    *,
    yscale: str = "log",
    bins: int | list[float] | None = None,
    n: int | None = None,
    clip_endpoints: bool = True,
    label: str | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a folded, scaled site frequency spectrum with scikit-allel's plot_sfs_folded_scaled.

    As for plot_sfs_scaled, the values are drawn as given; only the labels differ.
    """

    return _sfs_plot(
        "plot_sfs_folded_scaled",
        spectrum,
        output_path,
        yscale=yscale,
        bins=bins,
        n=n,
        clip_endpoints=clip_endpoints,
        label=label,
        plot_kwargs=plot_kwargs,
    )


def plot_joint_sfs(
    spectrum: list[list[float]],
    output_path: str,
    *,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a joint site frequency spectrum of two populations as a log-scaled image.

    ``spectrum`` is the matrix of site counts indexed by the derived allele count in population 1
    (rows) and population 2 (columns). scikit-allel draws it with population 1 on the x axis. The
    result holds the PNG path, the input, and the matrix passed to the image (transposed).
    """

    return _joint_plot("plot_joint_sfs", spectrum, output_path, imshow_kwargs=imshow_kwargs)


def plot_joint_sfs_folded(
    spectrum: list[list[float]],
    output_path: str,
    *,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a folded joint site frequency spectrum (minor allele counts) with scikit-allel."""

    return _joint_plot("plot_joint_sfs_folded", spectrum, output_path, imshow_kwargs=imshow_kwargs)


def plot_joint_sfs_scaled(
    spectrum: list[list[float]],
    output_path: str,
    *,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a scaled joint site frequency spectrum with scikit-allel's plot_joint_sfs_scaled.

    As for the 1D scaled plot, scikit-allel draws the values as given (the labels are those of
    the unscaled plot); pass the output of scikit-allel's joint_sfs_scaled.
    """

    return _joint_plot("plot_joint_sfs_scaled", spectrum, output_path, imshow_kwargs=imshow_kwargs)


def plot_joint_sfs_folded_scaled(
    spectrum: list[list[float]],
    output_path: str,
    *,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a folded, scaled joint site frequency spectrum with scikit-allel."""

    return _joint_plot(
        "plot_joint_sfs_folded_scaled", spectrum, output_path, imshow_kwargs=imshow_kwargs
    )


def fig_voight_painting(
    haplotypes: list[list[int]],
    output_path: str,
    *,
    index: int | None = None,
    palette: str = "colorblind",
    height_factor: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw scikit-allel's figure of shared haplotype prefixes on both flanks of one variant.

    ``haplotypes`` is a biallelic 0/1 matrix (variants x haplotypes). The flanks are cut at
    ``index`` (default: the middle variant). The figure shows the two Voight paintings above the
    EHH decay curves of each flank. The result holds the PNG path, the index used, and the
    painting of each flank from scikit-allel's voight_painting.
    """

    allel, np = _allel_and_numpy()
    require("seaborn", _ALLEL_EXTRA)
    matrix = _haplotype_matrix(np, haplotypes, "haplotypes")
    centre = matrix.shape[0] // 2 if index is None else _index(index, matrix.shape[0], "index")
    target, _ = _render(
        output_path,
        "fig_voight_painting",
        lambda: allel.fig_voight_painting(
            matrix, index=centre, palette=palette, height_factor=height_factor
        ),
    )
    left, _ = allel.voight_painting(matrix[: centre + 1][::-1])
    right, _ = allel.voight_painting(matrix[centre:])
    return _result(
        target,
        {
            "index": centre,
            "n_variants": int(matrix.shape[0]),
            "n_haplotypes": int(matrix.shape[1]),
            "palette": palette,
            "left_painting": left,
            "right_painting": right,
        },
    )


def plot_voight_painting(
    painting: list[list[int]],
    output_path: str,
    *,
    palette: str = "colorblind",
    flank: str = "right",
    height_factor: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw one Voight painting, in which each shared haplotype prefix is a coloured block.

    ``painting`` is the painting array of scikit-allel's voight_painting (variants x haplotypes;
    0 marks singleton haplotypes). ``flank`` 'left' reverses the variant axis. The result holds
    the PNG path, the painting, the flank and the number of colours (the largest group id).
    """

    if flank not in ("right", "left"):
        raise ConfigurationError("flank must be 'right' or 'left'.")
    allel, np = _allel_and_numpy()
    require("seaborn", _ALLEL_EXTRA)
    values = _numbers(np, painting, "painting", (2,), integer=True, minimum=0)
    target, _ = _render(
        output_path,
        "plot_voight_painting",
        lambda: allel.plot_voight_painting(
            values, palette=palette, flank=flank, height_factor=height_factor
        ),
    )
    return _result(
        target,
        {
            "painting": values,
            "flank": flank,
            "n_colors": int(values.max()),
            "palette": palette,
            "height_factor": height_factor,
        },
    )


def plot_haplotype_frequencies(
    haplotypes: list[list[int]],
    output_path: str,
    *,
    palette: str = "Paired",
    singleton_color: str = "w",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the frequencies of the distinct haplotypes of a region as one bar of coloured spans.

    ``haplotypes`` is a variants x haplotypes allele-index matrix. Each distinct haplotype that
    occurs more than once gets a palette colour and each singleton gets ``singleton_color``. The
    result holds the PNG path, the count of each distinct haplotype (scikit-allel's distinct_counts
    order) and the number of singletons.
    """

    allel, np = _allel_and_numpy()
    require("seaborn", _ALLEL_EXTRA)
    matrix = _haplotype_matrix(np, haplotypes, "haplotypes")
    array = allel.HaplotypeArray(matrix)
    counts = array.distinct_counts()
    target, _ = _render(
        output_path,
        "plot_haplotype_frequencies",
        lambda: allel.plot_haplotype_frequencies(
            array, palette=palette, singleton_color=singleton_color
        ),
    )
    return _result(
        target,
        {
            "n_variants": int(matrix.shape[0]),
            "n_haplotypes": int(matrix.shape[1]),
            "haplotype_counts": counts,
            "singleton_count": int(np.count_nonzero(counts == 1)),
            "palette": palette,
            "singleton_color": singleton_color,
        },
    )


def plot_moving_haplotype_frequencies(
    positions: list[int],
    haplotypes: list[list[int]],
    window_size: int,
    output_path: str,
    *,
    start: int = 0,
    stop: int | None = None,
    n: int | None = None,
    palette: str = "Paired",
    singleton_color: str = "w",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw haplotype frequencies in consecutive windows of ``window_size`` variants.

    ``positions`` are the 1-based positions of the rows of ``haplotypes`` (variants x haplotypes).
    Each column of the plot is a window and each row a frequency rank of the haplotypes in that
    window (rank 1 is the most frequent non-singleton); ``n`` keeps only the n most frequent ranks.
    ``start`` and ``stop`` are variant indices. The result holds the PNG path, the first and last
    position of each window, and the rank matrix as drawn (windows x haplotypes).
    """

    allel, np = _allel_and_numpy()
    require("seaborn", _ALLEL_EXTRA)
    position = _numbers(np, positions, "positions", (1,), integer=True, minimum=0)
    matrix = _haplotype_matrix(np, haplotypes, "haplotypes")
    if position.shape[0] != matrix.shape[0]:
        raise ConfigurationError("positions must have one entry per variant row of haplotypes.")
    size = _count(window_size, "window_size")
    if n is not None:
        n = _count(n, "n")
    array = allel.HaplotypeArray(matrix)
    target, ax = _render(
        output_path,
        "plot_moving_haplotype_frequencies",
        lambda: allel.plot_moving_haplotype_frequencies(
            position,
            array,
            size,
            start=start,
            stop=stop,
            n=n,
            palette=palette,
            singleton_color=singleton_color,
        ),
    )
    windows = allel.moving_statistic(
        position,
        statistic=lambda values: (values[0], values[-1]),
        size=size,
        start=start,
        stop=stop,
    )
    rank_matrix = np.asarray(ax.collections[0].get_array()).T
    return _result(
        target,
        {
            "window_size": size,
            "start": start,
            "stop": stop,
            "n": n,
            "window_bounds": windows,
            "rank_matrix": rank_matrix,
            "palette": palette,
            "singleton_color": singleton_color,
        },
    )


def plot_variant_locator(
    positions: list[int],
    output_path: str,
    *,
    step: int | None = None,
    start: int | None = None,
    stop: int | None = None,
    flip: bool = False,
    line_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw lines from the genome position of each variant to its place in the variant index.

    ``positions`` is one sorted contig of 1-based positions. ``step`` draws every step-th variant;
    by default scikit-allel uses one line per hundredth of the variants, here never less than
    every variant (scikit-allel fails with a zero step for fewer than 100 variants). ``start`` and
    ``stop`` bound the genome region. The result holds the PNG path, the step, the x limits and
    the number of lines drawn.
    """

    allel, np = _allel_and_numpy()
    position = _numbers(np, positions, "positions", (1,), integer=True, minimum=0)
    if position.size > 1 and (np.diff(position) < 0).any():
        raise ConfigurationError("positions must be sorted in ascending order.")
    step_value = max(1, position.size // 100) if step is None else _count(step, "step")
    kwargs = dict(line_kwargs or {})
    target, ax = _render(
        output_path,
        "plot_variant_locator",
        lambda: allel.plot_variant_locator(
            position, step=step_value, start=start, stop=stop, flip=flip, line_kwargs=kwargs
        ),
    )
    x_low, x_high = ax.get_xlim()
    return _result(
        target,
        {
            "n_variants": int(position.size),
            "step": step_value,
            "start": x_low,
            "stop": x_high,
            "flip": flip,
            "lines_drawn": len(ax.lines),
        },
    )


def plot_pairwise_distance(
    distances: list[float],
    output_path: str,
    *,
    labels: list[str] | None = None,
    colorbar: bool = True,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a pairwise distance matrix between samples as an image.

    ``distances`` is the condensed form (the upper triangle, row by row, n*(n-1)/2 values) or a
    square n x n matrix. ``labels`` gives one name per sample. The result holds the PNG path, the
    number of samples, the labels, and the square matrix as drawn.
    """

    allel, np = _allel_and_numpy()
    values = _numbers(np, distances, "distances", (1, 2))
    size = _matrix_size(values)
    if labels is not None and (not isinstance(labels, list) or len(labels) != size):
        raise ConfigurationError(f"labels must hold {size} names, one per sample.")
    kwargs = dict(imshow_kwargs or {})
    target, ax = _render(
        output_path,
        "plot_pairwise_distance",
        lambda: allel.plot_pairwise_distance(
            values, labels=labels, colorbar=colorbar, imshow_kwargs=kwargs
        ),
    )
    return _result(
        target,
        {
            "n_samples": size,
            "labels": labels,
            "colorbar": colorbar,
            "square_matrix": ax.images[0].get_array(),
        },
    )


def plot_pairwise_ld(
    ld_values: list[float],
    output_path: str,
    *,
    colorbar: bool = True,
    imshow_kwargs: dict[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw linkage disequilibrium (r squared) between all pairs of variants as a triangle image.

    ``ld_values`` is the condensed form (the upper triangle, row by row, of n*(n-1)/2 values for n
    variants) or a square matrix. Like scikit-allel, the image shows the lower triangle with the
    variant order flipped. The result holds the PNG path, the number of variants, the input and
    the matrix as drawn.
    """

    allel, np = _allel_and_numpy()
    values = _numbers(np, ld_values, "ld_values", (1, 2))
    size = _matrix_size(values)
    kwargs = dict(imshow_kwargs or {})
    target, ax = _render(
        output_path,
        "plot_pairwise_ld",
        lambda: allel.plot_pairwise_ld(values, colorbar=colorbar, imshow_kwargs=kwargs),
    )
    return _result(
        target,
        {
            "n_variants": size,
            "ld_values": values,
            "colorbar": colorbar,
            "drawn_matrix": ax.images[0].get_array(),
        },
    )


__all__ = [
    "fig_voight_painting",
    "plot_haplotype_frequencies",
    "plot_joint_sfs",
    "plot_joint_sfs_folded",
    "plot_joint_sfs_folded_scaled",
    "plot_joint_sfs_scaled",
    "plot_moving_haplotype_frequencies",
    "plot_pairwise_distance",
    "plot_pairwise_ld",
    "plot_sfs",
    "plot_sfs_folded",
    "plot_sfs_folded_scaled",
    "plot_sfs_scaled",
    "plot_variant_locator",
    "plot_voight_painting",
    "roh_poissonhmm",
]

_register_public_verbose(globals())
