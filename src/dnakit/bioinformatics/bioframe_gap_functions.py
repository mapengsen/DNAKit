"""bioframe operations as DNAKit wrappers: UCSC assembly tables, pairix reads, bigBed/bigWig writers
and interval plots.

Network wrappers call bioframe's fetch functions and record the service, the request, a UTC
retrieval time and the data. ``read_pairix`` needs the pypairix Python bindings. ``to_bigbed`` and
``to_bigwig`` are facades over the UCSC programs ``bedToBigBed`` and ``bedGraphToBigWig``, which the
user installs; they take explicit file paths. Plot functions write one PNG to a new path and return
the drawn data. They draw on a matplotlib Figure object, so no pyplot global state is used.

Missing optional packages raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Ordinary tests never use the network.
"""

from __future__ import annotations

import datetime as dt
import math
import re
import shutil
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_BIOFRAME_EXTRA = "bio-bioframe"
_PAIRIX_EXTRA = "bio-pairix"
_PLOT_EXTRA = "bio-plots"
_PROVIDERS = ("ucsc", "local")
_SERVICES = {
    "ucsc": "UCSC Genome Browser database tables (hgdownload.soe.ucsc.edu)",
    "local": "bioframe bundled assembly tables (no network access)",
}
_ASSEMBLY = re.compile(r"[A-Za-z0-9_.-]{1,64}")
_BED_FORMAT = re.compile(r"bed[1-9][0-9]*(\+[0-9]*)?")
_MAX_PLOT_INTERVALS = 100_000
_MAX_TAIL = 4000
_DEFAULT_FACECOLOR = "skyblue"
_DEFAULT_EDGECOLOR = "dimgray"
_HINT_BEDTOBIGBED = (
    "Install the UCSC bedToBigBed program (conda install -c bioconda ucsc-bedtobigbed, "
    "or download it from "
    "https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_BEDGRAPHTOBIGWIG = (
    "Install the UCSC bedGraphToBigWig program "
    "(conda install -c bioconda ucsc-bedgraphtobigwig, or download it "
    "from https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _bioframe_version() -> str:
    try:
        return version("bioframe")
    except (
        PackageNotFoundError
    ):  # pragma: no cover - bioframe is installed whenever its wrappers run.
        return "unknown"


def _plain(value: Any) -> Any:
    """Convert numpy and pandas scalars to JSON-friendly values; NaN and infinity become None."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item") and not isinstance(value, (str, Mapping)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _records(frame: Any) -> list[dict[str, Any]]:
    return [
        {str(key): _plain(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _assembly(value: object) -> str:
    if not isinstance(value, str) or not _ASSEMBLY.fullmatch(value):
        raise ConfigurationError(
            "db must be an assembly name such as hg38 (letters, digits, '_', '.' or '-' only)."
        )
    return value


def _provider(value: object) -> str:
    if value not in _PROVIDERS:
        raise ConfigurationError(f"provider must be one of {', '.join(_PROVIDERS)}.")
    return str(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new_path(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _optional_output(value: object, name: str) -> Path | None:
    return None if value is None else _new_path(value, name)


def _integers(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError(f"{name} must map chromosome names to lengths.")
    sizes: dict[str, int] = {}
    for chrom, length in value.items():
        _text(chrom, "chromosome name")
        if isinstance(length, bool) or not isinstance(length, int) or length < 1:
            raise ConfigurationError(f"The length of {chrom} in {name} must be a positive integer.")
        sizes[chrom] = length
    return sizes


def _wrap_fetch_error(exc: Exception, what: str) -> BackendExecutionError:
    return BackendExecutionError(f"bioframe could not {what}: {exc}")


def fetch_centromeres(
    db: str,
    *,
    provider: str = "ucsc",
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch the centromere position of every chromosome of an assembly (bioframe
    fetch_centromeres).
    """

    assembly = _assembly(db)
    source = _provider(provider)
    target = _optional_output(output_path, "output_path")
    bioframe = require("bioframe", _BIOFRAME_EXTRA)
    try:
        frame = bioframe.fetch_centromeres(assembly, provider=source)
    except (OSError, ValueError, KeyError) as exc:
        raise _wrap_fetch_error(
            exc, f"fetch centromeres for {assembly} from provider {source!r}"
        ) from exc
    if target is not None:
        frame.to_csv(target, sep="\t", index=False)
    return {
        "status": "completed",
        "service": _SERVICES[source],
        "request": {"db": assembly, "provider": source},
        "retrieved_at": _utc_now(),
        "bioframe_version": _bioframe_version(),
        "row_count": len(frame),
        "data": _records(frame),
        "output_path": None if target is None else str(target),
    }


def fetch_chromsizes(
    db: str,
    *,
    provider: str = "ucsc",
    as_bed: bool = False,
    filter_chroms: bool = True,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch the chromosome lengths of an assembly as a name-to-length map or a BED3 table (bioframe
    fetch_chromsizes).
    """

    assembly = _assembly(db)
    source = _provider(provider)
    as_table = _flag(as_bed, "as_bed")
    filtered = _flag(filter_chroms, "filter_chroms")
    target = _optional_output(output_path, "output_path")
    bioframe = require("bioframe", _BIOFRAME_EXTRA)
    try:
        sizes = bioframe.fetch_chromsizes(
            assembly, provider=source, as_bed=as_table, filter_chroms=filtered
        )
    except (OSError, ValueError, KeyError) as exc:
        raise _wrap_fetch_error(
            exc, f"fetch chromosome sizes for {assembly} from provider {source!r}"
        ) from exc
    if as_table:
        data: Any = _records(sizes)
        if target is not None:
            sizes.to_csv(target, sep="\t", index=False, header=False)
    else:
        data = {str(name): int(_plain(length)) for name, length in sizes.items()}
        if target is not None:
            sizes.to_csv(target, sep="\t", header=False)
    return {
        "status": "completed",
        "service": _SERVICES[source],
        "request": {
            "db": assembly,
            "provider": source,
            "as_bed": as_table,
            "filter_chroms": filtered,
        },
        "retrieved_at": _utc_now(),
        "bioframe_version": _bioframe_version(),
        "row_count": len(data),
        "data": data,
        "output_path": None if target is None else str(target),
    }


def read_pairix(
    path: str,
    region1: str,
    *,
    region2: str | None = None,
    chromsizes: Mapping[str, int] | None = None,
    usecols: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the pairs of a pairix-indexed file that fall in one region, or between two regions
    (bioframe read_pairix).
    """

    source = _existing(path, "path")
    first = _text(region1, "region1")
    second = None if region2 is None else _text(region2, "region2")
    sizes = None if chromsizes is None else _integers(chromsizes, "chromsizes")
    columns = None
    if usecols is not None:
        if isinstance(usecols, str) or not isinstance(usecols, Sequence) or not usecols:
            raise ConfigurationError("usecols must be a non-empty sequence of column names.")
        columns = [_text(name, "usecols item") for name in usecols]
    pypairix = require("pypairix", _PAIRIX_EXTRA)
    require("cytoolz", _PAIRIX_EXTRA)
    bioframe = require("bioframe", _BIOFRAME_EXTRA)
    pandas = require("pandas", _BIOFRAME_EXTRA)
    series = None if sizes is None else pandas.Series(sizes, dtype="int64")
    failures = (
        OSError,
        ValueError,
        KeyError,
        RuntimeError,
        getattr(pypairix, "PairixError", RuntimeError),
    )
    try:
        frame = bioframe.read_pairix(
            str(source), first, region2=second, chromsizes=series, usecols=columns
        )
    except failures as exc:
        raise BackendExecutionError(f"pypairix could not read {source} for {first}: {exc}") from exc
    return {
        "status": "completed",
        "backend": "pypairix",
        "bioframe_version": _bioframe_version(),
        "request": {"path": str(source), "region1": first, "region2": second, "usecols": columns},
        "row_count": len(frame),
        "columns": [str(name) for name in frame.columns],
        "data": _records(frame),
    }


def _program(name: str, explicit: str | None, hint: str) -> str:
    if explicit is not None:
        candidate = Path(_text(explicit, "program path")).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name} path does not exist: {candidate}", hint=hint)
        return str(candidate)
    found = shutil.which(name)
    if found is None:
        raise BackendUnavailableError(f"{name} was not found on PATH.", hint=hint)
    return found


def _run_tool(
    name: str,
    args: list[str],
    output: Path,
    explicit: str | None,
    hint: str,
) -> dict[str, Any]:
    executable = _program(name, explicit, hint)
    result = execute_bounded_command(
        executable,
        args,
        backend_id=f"ucsc-{name.lower()}",
        cwd=output.parent,
        verbose=False,
    )
    produced = output.is_file()
    return {
        "status": "completed" if result.return_code == 0 and produced else "failed",
        "backend": name,
        "executable": executable,
        "command": [name, *args],
        "return_code": result.return_code,
        "outputs": [str(output)] if produced else [],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


def to_bigbed(
    bed: str,
    chromsizes: str,
    output_bigbed: str,
    *,
    bed_format: str | None = None,
    bedtobigbed_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a sorted BED file as an indexed bigBed file with UCSC bedToBigBed (explicit chromsizes
    file).
    """

    bed_file = _existing(bed, "bed")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigbed, "output_bigbed")
    args: list[str] = []
    if bed_format is not None:
        if not isinstance(bed_format, str) or not _BED_FORMAT.fullmatch(bed_format):
            raise ConfigurationError(
                "bed_format must look like bed3, bed6 or bed6+4 (bedToBigBed -type)."
            )
        args.append(f"-type={bed_format}")
    args += [str(bed_file), str(sizes_file), str(out)]
    return _run_tool("bedToBigBed", args, out, bedtobigbed_path, _HINT_BEDTOBIGBED)


def to_bigwig(
    bedgraph: str,
    chromsizes: str,
    output_bigwig: str,
    *,
    bedgraphtobigwig_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a sorted bedGraph file as an indexed bigWig file with UCSC bedGraphToBigWig (explicit
    chromsizes file).
    """

    graph_file = _existing(bedgraph, "bedgraph")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigwig, "output_bigwig")
    args = [str(graph_file), str(sizes_file), str(out)]
    return _run_tool("bedGraphToBigWig", args, out, bedgraphtobigwig_path, _HINT_BEDGRAPHTOBIGWIG)


def _figsize(value: object) -> tuple[float, float]:
    if isinstance(value, str) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError("figsize must be two positive numbers (width, height per track).")
    width, height = value
    for item in (width, height):
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not 0 < item <= 100:
            raise ConfigurationError("figsize values must be numbers in (0, 100].")
    return float(width), float(height)


def _numbers(values: object, name: str) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty sequence of numbers.")
    if len(values) > _MAX_PLOT_INTERVALS:
        raise ConfigurationError(f"{name} can hold at most {_MAX_PLOT_INTERVALS} values.")
    result: list[float] = []
    for value in values:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ConfigurationError(f"{name} must contain finite numbers only.")
        result.append(value)
    return result


def _figure_classes() -> tuple[Any, Any, Any]:
    figure = require("matplotlib.figure", _PLOT_EXTRA)
    canvas = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    patches = require("matplotlib.patches", _PLOT_EXTRA)
    return figure, canvas, patches


def _draw_intervals(
    axis: Any,
    patches: Any,
    starts: Sequence[float],
    ends: Sequence[float],
    levels: Sequence[float],
    labels: Sequence[str | None],
    colors: Sequence[str],
) -> None:
    """Draw one 0.6-high rectangle per interval at its level, with an optional centred label."""

    for start, end, level, label, color in zip(starts, ends, levels, labels, colors, strict=True):
        axis.add_patch(
            patches.Rectangle(
                (start, level - 0.3),
                end - start,
                0.6,
                facecolor=color,
                edgecolor=_DEFAULT_EDGECOLOR,
            )
        )
        if label is not None:
            axis.text(
                (start + end) / 2,
                level,
                label,
                horizontalalignment="center",
                verticalalignment="center",
            )


def _finish_axis(
    axis: Any, top: float, xlim: tuple[float, float], show_coords: bool, title: str | None
) -> None:
    axis.set_ylim(-0.5, top + 0.5)
    axis.set_xlim(xlim[0], xlim[1])
    axis.set_frame_on(False)
    axis.set_yticks([])
    if not show_coords:
        axis.set_xticks([])
    if title is not None:
        axis.set_title(title)


def _save_png(figure: Any, canvas: Any, path: Path) -> None:
    canvas(figure)
    figure.savefig(str(path), format="png", dpi=120)


def _read_intervals(pandas: Any, source: Path, label_column: int | None) -> Any:
    try:
        frame = pandas.read_csv(
            source, sep="\t", header=None, dtype={0: str}, keep_default_na=False
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(
            f"intervals must be a tab-separated file without a header: {exc}"
        ) from exc
    if frame.shape[0] == 0 or frame.shape[0] > _MAX_PLOT_INTERVALS:
        raise ConfigurationError(f"intervals must hold 1 to {_MAX_PLOT_INTERVALS} rows.")
    if frame.shape[1] < 3:
        raise ConfigurationError(
            "intervals needs at least three tab-separated columns: chrom, start, end."
        )
    if label_column is not None and (
        isinstance(label_column, bool)
        or not isinstance(label_column, int)
        or not 0 <= label_column < frame.shape[1]
    ):
        raise ConfigurationError(
            f"label_column must be a column index from 0 to {frame.shape[1] - 1}."
        )
    try:
        frame[1] = frame[1].astype("int64")
        frame[2] = frame[2].astype("int64")
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(
            "The start and end columns of intervals must be integers."
        ) from exc
    if (frame[1] < 0).any() or (frame[1] > frame[2]).any():
        raise ConfigurationError("Every interval needs 0 <= start <= end.")
    return frame


def plot_intervals(
    intervals: str,
    output_png: str,
    *,
    label_column: int | None = None,
    color: str = _DEFAULT_FACECOLOR,
    show_coords: bool = False,
    figsize: Sequence[float] = (10.0, 2.0),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a BED-like interval file as one stacked track per chromosome in a PNG (bioframe
    plot_intervals).
    """

    source = _existing(intervals, "intervals")
    out = _new_path(output_png, "output_png")
    width, track_height = _figsize(figsize)
    fill = _text(color, "color")
    coords = _flag(show_coords, "show_coords")
    pandas = require("pandas", _BIOFRAME_EXTRA)
    numpy = require("numpy", _BIOFRAME_EXTRA)
    bioframe = require("bioframe", _BIOFRAME_EXTRA)
    frame = _read_intervals(pandas, source, label_column)
    figure_module, canvas_module, patches = _figure_classes()

    chroms = sorted(set(frame[0]))
    figure = figure_module.Figure(figsize=(width, track_height * len(chroms)))
    tracks: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for index, chrom in enumerate(chroms, start=1):
        part = frame[frame[0] == chrom]
        starts = [int(value) for value in part[1]]
        ends = [int(value) for value in part[2]]
        labels: list[str | None] = (
            [str(value) for value in part[label_column]]
            if label_column is not None
            else [None] * len(starts)
        )
        levels = [
            int(v)
            for v in bioframe.core.arrops.stack_intervals(
                numpy.asarray(starts, dtype=numpy.int64), numpy.asarray(ends, dtype=numpy.int64)
            )
        ]
        axis = figure.add_subplot(len(chroms), 1, index)
        _draw_intervals(axis, patches, starts, ends, levels, labels, [fill] * len(starts))
        _finish_axis(axis, max(levels), (-0.5, max(ends) + 0.5), coords, chrom)
        tracks.append({"chrom": chrom, "interval_count": len(starts), "max_level": max(levels)})
        rows.extend(
            {"chrom": chrom, "start": s, "end": e, "level": lv, "label": lb}
            for s, e, lv, lb in zip(starts, ends, levels, labels, strict=True)
        )
    _save_png(figure, canvas_module.FigureCanvasAgg, out)
    return {
        "status": "completed",
        "backend": "bioframe + matplotlib",
        "image_path": str(out),
        "data": {
            "source": str(source),
            "interval_count": len(rows),
            "chromosomes": chroms,
            "tracks": tracks,
            "intervals": rows,
        },
    }


def plot_intervals_arr(
    starts: Sequence[float],
    ends: Sequence[float],
    output_png: str,
    *,
    levels: Sequence[float] | None = None,
    labels: Sequence[str | None] | None = None,
    colors: str | Sequence[str] | None = None,
    xlim: Sequence[float] | None = None,
    show_coords: bool = False,
    figsize: Sequence[float] = (10.0, 2.0),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot one set of intervals given as start and end arrays in a PNG (bioframe
    plot_intervals_arr).
    """

    start_values = _numbers(starts, "starts")
    end_values = _numbers(ends, "ends")
    if len(start_values) != len(end_values):
        raise ConfigurationError("starts and ends must have the same length.")
    if any(s > e for s, e in zip(start_values, end_values, strict=True)):
        raise ConfigurationError("Every interval needs start <= end.")
    count = len(start_values)
    out = _new_path(output_png, "output_png")
    width, track_height = _figsize(figsize)
    coords = _flag(show_coords, "show_coords")
    bioframe = require("bioframe", _BIOFRAME_EXTRA)
    numpy = require("numpy", _BIOFRAME_EXTRA)
    if levels is None:
        level_values: list[float] = [
            int(v)
            for v in bioframe.core.arrops.stack_intervals(
                numpy.asarray(start_values, dtype=numpy.float64),
                numpy.asarray(end_values, dtype=numpy.float64),
            )
        ]
    else:
        level_values = _numbers(levels, "levels")
        if len(level_values) != count:
            raise ConfigurationError("levels must have one value per interval.")
    if labels is None:
        label_values: list[str | None] = [None] * count
    else:
        if isinstance(labels, str) or not isinstance(labels, Sequence) or len(labels) != count:
            raise ConfigurationError("labels must be a sequence with one entry per interval.")
        label_values = [None if item is None else str(item) for item in labels]
    if colors is None:
        color_values = [_DEFAULT_FACECOLOR] * count
    elif isinstance(colors, str):
        color_values = [_text(colors, "colors")] * count
    else:
        if not isinstance(colors, Sequence) or len(colors) != count:
            raise ConfigurationError("colors must be one colour name or one per interval.")
        color_values = [_text(item, "colors item") for item in colors]
    if xlim is None:
        x_low, x_high = -0.5, max(end_values) + 0.5
    else:
        bounds = _numbers(xlim, "xlim")
        if len(bounds) != 2 or bounds[0] >= bounds[1]:
            raise ConfigurationError("xlim must be two numbers with lower < upper.")
        x_low, x_high = bounds
    figure_module, canvas_module, patches = _figure_classes()
    figure = figure_module.Figure(figsize=(width, track_height))
    axis = figure.add_subplot(1, 1, 1)
    _draw_intervals(
        axis, patches, start_values, end_values, level_values, label_values, color_values
    )
    top = max(level_values)
    _finish_axis(axis, top, (x_low, x_high), coords, None)
    _save_png(figure, canvas_module.FigureCanvasAgg, out)
    return {
        "status": "completed",
        "backend": "bioframe + matplotlib",
        "image_path": str(out),
        "data": {
            "interval_count": count,
            "starts": start_values,
            "ends": end_values,
            "levels": level_values,
            "labels": label_values,
            "colors": color_values,
            "xlim": [x_low, x_high],
            "ylim": [-0.5, top + 0.5],
        },
    }


__all__ = [
    "fetch_centromeres",
    "fetch_chromsizes",
    "plot_intervals",
    "plot_intervals_arr",
    "read_pairix",
    "to_bigbed",
    "to_bigwig",
]
_register_public_verbose(globals())
