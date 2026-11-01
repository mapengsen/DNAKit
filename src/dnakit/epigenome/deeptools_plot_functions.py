"""deepTools plot wrappers and plot-data helpers: fragment-length density, GC-bias plot, matrix
heatmap,
enrichment bars, and blacklist-aware template lengths.

Each wrapper calls the deepTools function that does the work: ``bamPEFragmentSize.getDensity``,
``computeGCBias.plotGCbias``, ``heatmapper.read_matrix_file`` with ``plotHeatmap.plotMatrix``,
``plotEnrichment.plotEnrichment``, ``utilities.getTLen`` and ``mapReduce.blSubtract``. Plot
functions
write a new PNG file (creating its parent directory), and return the absolute path with the numbers
that
were drawn; they never return a figure. Helpers that only lay out or draw part of a figure are not
wrapped.

deepTools, matplotlib, pysam and deeptoolsintervals are imported only when a function runs. No
network is used.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import statistics
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.epigenome.deeptools_functions import _integer, _number, _plain, _provenance, _text
from dnakit.exceptions import ConfigurationError

_DEEPTOOLS_EXTRA = "bio-deeptools"
_MAX_LENGTH = 2_000_000_000
_MAX_FRAGMENTS = 10_000_000
_MAX_REGION = 10_000_000
_MAX_COUNT = 1_000_000_000_000
_MAX_READS = 100_000_000
_MAX_FEATURES = 10_000
_MAX_SAMPLES = 1_000
_MAX_GC_ROWS = 10_000_000
_GC_WINDOW = (0.2, 0.7)
_AVERAGE_TYPES = ("mean", "median", "sum", "min", "max", "std")
_PLOT_TYPES = ("lines", "se", "std", "overlapped_lines", "fill")
_WHAT_TO_SHOW = ("plot, heatmap and colorbar", "plot and heatmap", "heatmap and colorbar")
_DEFAULT_ENRICHMENT_COLOR = "#9E4A06"
_OPERATION_ERRORS = (
    AssertionError,
    AttributeError,
    IndexError,
    KeyError,
    OSError,
    RuntimeError,
    TypeError,
    ValueError,
    ZeroDivisionError,
)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _DEEPTOOLS_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _DEEPTOOLS_EXTRA)


def _plotting_module(dotted: str) -> Any:
    """Import a deepTools module; plotHeatmap changes numpy error handling and interactive mode on
    import.

    Those two global settings are restored after the import so that importing this package does not
    silence floating-point warnings for the rest of the process.
    """

    numpy = require("numpy", _DEEPTOOLS_EXTRA)
    plt = _pyplot()
    errors = numpy.geterr()
    interactive = plt.isinteractive()
    try:
        return require(dotted, _DEEPTOOLS_EXTRA)
    finally:
        numpy.seterr(**errors)
        if interactive:
            plt.ion()


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choice(value: object, options: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of: {', '.join(options)}.")
    return value


def _short_text(value: object, name: str, limit: int = 500) -> str:
    if (
        not isinstance(value, str)
        or len(value) > limit
        or any(char in value for char in "\x00\r\n")
    ):
        raise ConfigurationError(
            f"{name} must be a single-line string of at most {limit} characters."
        )
    return value


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _output_png(output_path: object) -> Path:
    """A new PNG path; the parent directory is created. Existing files are never overwritten."""

    if not isinstance(output_path, (str, Path)):
        raise ConfigurationError("output_path must be a file path.")
    target = Path(_text(str(output_path), "output_path")).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _sequence(value: object, name: str, maximum: int) -> Sequence[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if not 1 <= len(value) <= maximum:
        raise ConfigurationError(f"{name} must contain 1 to {maximum} entries.")
    return value


def _numbers(value: object, name: str, minimum: float, maximum: float, limit: int) -> list[float]:
    return [_number(item, name, minimum, maximum) for item in _sequence(value, name, limit)]


def _labels(value: object, name: str, limit: int) -> list[str]:
    return [_short_text(item, name, limit=255) for item in _sequence(value, name, limit)]


def _draw(plt: Any, target: Path, operation: str, call: Callable[[], Any]) -> None:
    """Run one deepTools drawing call, close the figures it leaves open, and report its failures."""

    before = set(plt.get_fignums())
    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            call()
    except SystemExit as exc:
        lines = [line for line in captured.getvalue().splitlines() if line.strip()]
        detail = lines[-1] if lines else str(exc.code)
        raise ConfigurationError(f"deepTools rejected the {operation}: {detail}") from exc
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"deepTools could not draw the {operation}: {exc}") from exc
    finally:
        for number in plt.get_fignums():
            if number not in before:
                plt.close(number)
    if not target.is_file():
        raise ConfigurationError(f"deepTools did not write {target}.")


def _frequency_table(value: object) -> list[list[float]]:
    """Rows of (fragments, reads, observed/expected ratio), one per fragment length as computeGCBias
    writes them."""

    rows = _sequence(value, "frequencies", _MAX_GC_ROWS)
    if len(rows) < 3:
        raise ConfigurationError("frequencies must have at least 3 rows.")
    table: list[list[float]] = []
    for row in rows:
        if isinstance(row, (str, bytes, Mapping)) or not isinstance(row, Sequence) or len(row) != 3:
            raise ConfigurationError(
                "each frequencies row must hold three numbers: fragments, reads and "
                "observed/expected ratio."
            )
        fragments = _number(row[0], "frequencies", 0, _MAX_COUNT)
        reads = _number(row[1], "frequencies", 0, _MAX_COUNT)
        ratio = _number(row[2], "frequencies", 0, _MAX_COUNT)
        if ratio <= 0:
            raise ConfigurationError("the observed/expected ratio in frequencies must be positive.")
        table.append([fragments, reads, ratio])
    return table


def _reads_table(value: object) -> list[list[float]]:
    """Rows of (reads, GC fraction) for regions of one size, as computeGCBias tabulates them."""

    rows = _sequence(value, "reads_per_gc", _MAX_GC_ROWS)
    table: list[list[float]] = []
    for row in rows:
        if isinstance(row, (str, bytes, Mapping)) or not isinstance(row, Sequence) or len(row) != 2:
            raise ConfigurationError(
                "each reads_per_gc row must hold two numbers: reads and GC fraction."
            )
        reads = _number(row[0], "reads_per_gc", 0, _MAX_COUNT)
        gc = _number(row[1], "reads_per_gc", 0, 1)
        table.append([reads, gc])
    return table


def fragment_length_density(
    lengths: list[float],
    min_value: float,
    max_value: float,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Density of fragment lengths in 100 bins over [min_value, max_value], as deepTools' getDensity
    computes it.

    Lengths outside the range are not counted. The histogram is drawn in a scratch figure that is
    closed
    again, so the caller's current matplotlib figure is not cleared.
    """

    low = _number(min_value, "min_value", 0, _MAX_LENGTH)
    high = _number(max_value, "max_value", 0, _MAX_LENGTH)
    if high <= low:
        raise ConfigurationError("max_value must be greater than min_value.")
    values = _numbers(lengths, "lengths", 0, _MAX_LENGTH, _MAX_FRAGMENTS)
    density_module = _plotting_module("deeptools.bamPEFragmentSize")
    plt = _pyplot()
    scratch = plt.figure()
    try:
        density, edges = density_module.getDensity(values, low, high)
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(
            f"deepTools could not compute the fragment-length density: {exc}"
        ) from exc
    finally:
        plt.close(scratch)
    return {
        "status": "completed",
        "operation": "fragment_length_density",
        "fragments": len(values),
        "range": [low, high],
        "bins": len(density),
        "bin_edges": _plain(edges),
        "density": _plain(density),
        "provenance": _provenance(
            "deeptools",
            "100-bin histogram over the range with density=True (bamPEFragmentSize.getDensity)",
            note="fragments outside the range are not counted",
        ),
    }


def plot_gc_bias(
    frequencies: list[list[float]],
    reads_per_gc: list[list[float]],
    region_size: int,
    output_path: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw deepTools' GC-bias plot: reads per GC bin (box plots) and the observed/expected log2
    ratio.

    ``frequencies`` has one row per fragment length: (fragments, reads, observed/expected ratio), as
    computeGCBias writes its frequency table. ``reads_per_gc`` has one row per region: (reads, GC
    fraction).
    Only GC bins from 0.2 to 0.7 are drawn. Saves a PNG at a new path and returns its path and the
    plotted values.
    """

    table = _frequency_table(frequencies)
    reads = _reads_table(reads_per_gc)
    size = _integer(region_size, "region_size", 1, _MAX_REGION)
    gc_bias = _plotting_module("deeptools.computeGCBias")
    numpy = require("numpy", _DEEPTOOLS_EXTRA)
    frequency_array = numpy.array(table, dtype=float)
    reads_array = numpy.array(reads, dtype=float)
    target = _output_png(output_path)
    plt = _pyplot()
    _draw(
        plt,
        target,
        "GC-bias plot",
        lambda: gc_bias.plotGCbias(
            str(target), frequency_array, reads_array, size, image_format="png"
        ),
    )
    binned, left_edges = gc_bias.bin_by(reads_array[:, 0], reads_array[:, 1], nbins=100)
    bins = []
    for edge, values in zip(left_edges, binned, strict=True):
        if _GC_WINDOW[0] <= edge <= _GC_WINDOW[1]:
            bins.append(
                {
                    "gc": float(edge),
                    "regions": len(values),
                    "median_reads": float(numpy.median(values)) if len(values) else None,
                }
            )
    rows = len(table)
    return {
        "path": str(target),
        "data": {
            "operation": "plot_gc_bias",
            "region_size": size,
            "gc_bins": bins,
            "gc_fraction": _plain(numpy.linspace(0, 1, rows)),
            "log2_observed_expected": _plain(numpy.log2(frequency_array[:, 2])),
            "window": list(_GC_WINDOW),
            "provenance": _provenance(
                "deeptools", "computeGCBias.plotGCbias, bins from computeGCBias.bin_by"
            ),
        },
    }


def plot_matrix_heatmap(
    matrix_path: str,
    output_path: str,
    *,
    plot_title: str = "",
    colormap: str = "binary",
    missing_data_color: str = "black",
    alpha: float = 1.0,
    z_min: float | None = None,
    z_max: float | None = None,
    y_min: float | None = None,
    y_max: float | None = None,
    average_type: str = "median",
    start_label: str = "TSS",
    end_label: str = "TES",
    heatmap_height_cm: float = 25,
    heatmap_width_cm: float = 7.5,
    per_group: bool = False,
    what_to_show: str = "plot, heatmap and colorbar",
    plot_type: str = "lines",
    interpolation_method: str = "auto",
    dpi: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Heatmap (and summary profile) of a deepTools computeMatrix output, drawn with deepTools'
    plotMatrix.

    ``matrix_path`` is the gzipped matrix that computeMatrix writes. It is read with deepTools'
    heatmapper
    class, and the heatmap is saved as a new PNG. When z_min or z_max is None, deepTools picks the
    1st and
    98th percentiles of the matrix.
    """

    source = _existing_file(matrix_path, "matrix_path")
    title = _short_text(plot_title, "plot_title")
    cmap = _short_text(colormap, "colormap", limit=100)
    missing = _short_text(missing_data_color, "missing_data_color", limit=100)
    transparency = _number(alpha, "alpha", 0, 1)
    low = None if z_min is None else _number(z_min, "z_min", -_MAX_COUNT, _MAX_COUNT)
    high = None if z_max is None else _number(z_max, "z_max", -_MAX_COUNT, _MAX_COUNT)
    if low is not None and high is not None and high <= low:
        raise ConfigurationError("z_max must be greater than z_min.")
    lower_y = None if y_min is None else _number(y_min, "y_min", -_MAX_COUNT, _MAX_COUNT)
    upper_y = None if y_max is None else _number(y_max, "y_max", -_MAX_COUNT, _MAX_COUNT)
    statistic = _choice(average_type, _AVERAGE_TYPES, "average_type")
    start = _short_text(start_label, "start_label", limit=100)
    end = _short_text(end_label, "end_label", limit=100)
    height = _number(heatmap_height_cm, "heatmap_height_cm", 0.5, 100)
    width = _number(heatmap_width_cm, "heatmap_width_cm", 0.5, 100)
    grouped = _flag(per_group, "per_group")
    show = _choice(what_to_show, _WHAT_TO_SHOW, "what_to_show")
    line_type = _choice(plot_type, _PLOT_TYPES, "plot_type")
    interpolation = _text(interpolation_method, "interpolation_method", limit=100)
    resolution = _integer(dpi, "dpi", 50, 1200)
    plot_module = _plotting_module("deeptools.plotHeatmap")
    heatmapper = require("deeptools.heatmapper", _DEEPTOOLS_EXTRA)
    hm = heatmapper.heatmapper()
    try:
        hm.read_matrix_file(str(source))
    except _OPERATION_ERRORS as exc:
        raise ConfigurationError(f"{source.name} is not a computeMatrix output: {exc}") from exc
    target = _output_png(output_path)
    color_map = {
        "colorMap": [cmap],
        "missingDataColor": missing,
        "alpha": transparency,
        "colorList": None,
        "colorNumber": 255,
    }
    plt = _pyplot()
    _draw(
        plt,
        target,
        "matrix heatmap",
        lambda: plot_module.plotMatrix(
            hm,
            str(target),
            colorMapDict=color_map,
            plotTitle=title,
            zMin=None if low is None else [low],
            zMax=None if high is None else [high],
            yMin=lower_y,
            yMax=upper_y,
            averageType=statistic,
            startLabel=start,
            endLabel=end,
            heatmapHeight=height,
            heatmapWidth=width,
            perGroup=grouped,
            whatToShow=show,
            plot_type=line_type,
            image_format="png",
            dpi=resolution,
            interpolation_method=interpolation,
        ),
    )
    boundaries = [int(value) for value in hm.matrix.group_boundaries]
    sample_bounds = [int(value) for value in hm.matrix.sample_boundaries]
    return {
        "path": str(target),
        "data": {
            "operation": "plot_matrix_heatmap",
            "matrix": str(source),
            "samples": [str(label) for label in hm.matrix.sample_labels],
            "groups": [str(label) for label in hm.matrix.group_labels],
            "regions": boundaries[-1],
            "group_boundaries": boundaries,
            "group_sizes": [boundaries[i + 1] - boundaries[i] for i in range(len(boundaries) - 1)],
            "sample_boundaries": sample_bounds,
            "columns_per_sample": [
                sample_bounds[i + 1] - sample_bounds[i] for i in range(len(sample_bounds) - 1)
            ],
            "parameters": _plain(hm.parameters),
            "colormap": cmap,
            "z_min": low,
            "z_max": high,
            "average_type": statistic,
            "plot_type": line_type,
            "what_to_show": show,
        },
    }


def plot_enrichment_counts(
    features: list[str],
    sample_labels: list[str],
    feature_counts: list[dict[str, int]],
    total_counts: list[int],
    output_path: str,
    *,
    per_sample: bool = False,
    variable_scales: bool = False,
    plot_title: str = "",
    colors: list[str] | None = None,
    alpha: float = 0.9,
    plot_width_cm: float = 20,
    plot_height_cm: float = 20,
    panels_per_row: int = 4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enrichment bar plot from per-sample feature counts, drawn with deepTools' plotEnrichment.

    ``feature_counts`` has one dict per sample that maps each feature to its read count, and
    ``total_counts`` holds each sample's total read count. Bars show each feature's percentage of
    the total.
    Counting reads over a BAM and BED file is the job of the plot_enrichment facade, not of this
    function.
    """

    names = _labels(features, "features", _MAX_FEATURES)
    labels = _labels(sample_labels, "sample_labels", _MAX_SAMPLES)
    rows = _sequence(feature_counts, "feature_counts", _MAX_SAMPLES)
    if len(rows) != len(labels):
        raise ConfigurationError("feature_counts needs one dict per sample label.")
    counts_rows: list[dict[str, int]] = []
    for row in rows:
        if not isinstance(row, Mapping):
            raise ConfigurationError(
                "each feature_counts entry must be a dict of feature to count."
            )
        counts: dict[str, int] = {}
        for name in names:
            if name not in row:
                raise ConfigurationError(f"feature_counts is missing the feature {name}.")
            counts[name] = _integer(row[name], "feature count", 0, _MAX_READS)
        counts_rows.append(counts)
    totals_list = _sequence(total_counts, "total_counts", _MAX_SAMPLES)
    if len(totals_list) != len(labels):
        raise ConfigurationError("total_counts needs one value per sample label.")
    totals = [_integer(item, "total_counts", 1, _MAX_READS) for item in totals_list]
    bars = len(names) if per_sample else len(labels)
    if colors is None:
        bar_colors = [_DEFAULT_ENRICHMENT_COLOR] * bars
    else:
        bar_colors = [
            _short_text(item, "colors", limit=100) for item in _sequence(colors, "colors", bars)
        ]
        if len(bar_colors) < bars:
            raise ConfigurationError(f"{bars} colors are needed for this layout.")
    flag = _flag(variable_scales, "variable_scales")
    grouped = _flag(per_sample, "per_sample")
    transparency = _number(alpha, "alpha", 0, 1)
    width = _number(plot_width_cm, "plot_width_cm", 1, 200)
    height = _number(plot_height_cm, "plot_height_cm", 1, 200)
    panels = _integer(panels_per_row, "panels_per_row", 1, 100)
    title = _short_text(plot_title, "plot_title")
    enrichment = _plotting_module("deeptools.plotEnrichment")
    target = _output_png(output_path)
    args = argparse.Namespace(
        bamfiles=list(labels),
        perSample=grouped,
        labels=list(labels),
        colors=bar_colors,
        alpha=transparency,
        variableScales=flag,
        plotHeight=height,
        plotWidth=width,
        plotTitle=title,
        numPlotsPerRow=panels,
        plotFile=str(target),
        plotFileFormat="png",
    )
    plt = _pyplot()
    _draw(
        plt,
        target,
        "enrichment plot",
        lambda: enrichment.plotEnrichment(args, counts_rows, totals, names),
    )
    samples = []
    for label, counts, total in zip(labels, counts_rows, totals, strict=True):
        samples.append(
            {
                "sample": label,
                "total_counts": total,
                "features": {
                    name: {"count": counts[name], "percent": 100.0 * counts[name] / total}
                    for name in names
                },
            }
        )
    return {
        "path": str(target),
        "data": {
            "operation": "plot_enrichment_counts",
            "features": names,
            "per_sample": grouped,
            "bars_per_panel": bars,
            "samples": samples,
            "provenance": _provenance(
                "deeptools", "percent of total reads per feature (plotEnrichment)"
            ),
        },
    }


def template_lengths_from_bam(
    bam_path: str,
    *,
    not_abs: bool = False,
    max_reads: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Template length of each alignment in a BAM or SAM file, as deepTools' getTLen computes it.

    For a read with a TLEN field this is the absolute TLEN, or the signed TLEN when ``not_abs`` is
    true.
    For a read without TLEN it is the sum of the CIGAR operations M, D, = and X (codes 0, 2, 7 and
    8),
    which is what deepTools does. Reading stops after ``max_reads`` alignments.
    """

    source = _existing_file(bam_path, "bam_path")
    flag = _flag(not_abs, "not_abs")
    limit = _integer(max_reads, "max_reads", 1, _MAX_READS)
    utilities = require("deeptools.utilities", _DEEPTOOLS_EXTRA)
    pysam = require("pysam", _DEEPTOOLS_EXTRA)
    lengths: list[int] = []
    truncated = False
    try:
        with pysam.AlignmentFile(str(source), "r") as handle:
            for read in handle:
                if len(lengths) == limit:
                    truncated = True
                    break
                lengths.append(int(utilities.getTLen(read, notAbs=flag)))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{source.name} could not be read as alignments: {exc}") from exc
    summary = None
    if lengths:
        summary = {"min": min(lengths), "median": statistics.median(lengths), "max": max(lengths)}
    return {
        "status": "completed",
        "operation": "template_lengths_from_bam",
        "bam": str(source),
        "not_abs": flag,
        "alignments": len(lengths),
        "truncated": truncated,
        "template_lengths": lengths,
        "summary": summary,
        "provenance": _provenance(
            "deeptools",
            "getTLen: absolute TLEN when present, otherwise the sum of CIGAR operations "
            "M, D, = and X "
            "(codes 0, 2, 7, 8)",
        ),
    }


def subtract_blacklist_intervals(
    blacklist_path: str,
    chrom: str,
    start: int,
    end: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """The parts of [start, end) on one chromosome that no BED blacklist interval overlaps.

    Uses deepTools' blSubtract on a deeptoolsintervals tree of the blacklist file. Intervals are
    0-based,
    half-open. The overlapping blacklist intervals are returned as they are stored in the file.
    """

    source = _existing_file(blacklist_path, "blacklist_path")
    name = _text(chrom, "chrom", limit=255)
    left = _integer(start, "start", 0, _MAX_LENGTH)
    right = _integer(end, "end", 1, _MAX_LENGTH)
    if right <= left:
        raise ConfigurationError("end must be greater than start.")
    intervals = require("deeptoolsintervals", _DEEPTOOLS_EXTRA)
    mapreduce = require("deeptools.mapReduce", _DEEPTOOLS_EXTRA)
    try:
        tree = intervals.GTF(str(source))
    except Exception as exc:  # deeptoolsintervals raises plain exceptions for malformed BED files
        raise ConfigurationError(f"{source.name} could not be read as a BED file: {exc}") from exc
    overlaps = tree.findOverlaps(name, left, right) or []
    blacklisted = [[int(item[0]), int(item[1])] for item in overlaps]
    remaining = [[int(a), int(b)] for a, b in mapreduce.blSubtract(tree, name, [left, right])]
    return {
        "status": "completed",
        "operation": "subtract_blacklist_intervals",
        "blacklist": str(source),
        "chrom": name,
        "start": left,
        "end": right,
        "blacklisted_overlaps": blacklisted,
        "remaining": remaining,
        "remaining_bases": sum(b - a for a, b in remaining),
        "provenance": _provenance(
            "deeptools",
            "blSubtract on a deeptoolsintervals tree of the BED file",
            coordinate_system="0-based half-open",
        ),
    }


__all__ = [
    "fragment_length_density",
    "plot_enrichment_counts",
    "plot_gc_bias",
    "plot_matrix_heatmap",
    "subtract_blacklist_intervals",
    "template_lengths_from_bam",
]

_register_public_verbose(globals())
