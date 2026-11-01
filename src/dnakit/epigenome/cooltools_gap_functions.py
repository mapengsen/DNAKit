"""cooltools operations not covered by the existing DNAKit cooltools helpers (Fig. 2 coverage gap).

Download and listing functions wrap cooltools' example-dataset helpers, which fetch a dataset table and its files
over HTTP(S) and cache them in a data folder that the caller names. The listing can also read a local table, so it
runs offline. The saddle plot wraps cooltools.api.saddle.saddleplot, which draws through pyplot's current figure;
the wrapper lets it create the figure, saves that figure to a new PNG and closes it. The command-line helpers of
cooltools.cli are not wrapped (see the disposition file).
"""

from __future__ import annotations

import datetime as _dt
import math
import re
import warnings
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-cooltools"
_PLOT_EXTRA = "bio-plots"
_SAVE_DPI = 150
_URL = re.compile(r"https?://\S+")
_SERVICE_DATASETS = (
    "cooltools example datasets (dataset table and files from the cooltools repository)"
)
_SERVICE_FILE = "remote file over HTTP(S), fetched with cooltools' download_file"
_SERVICE_PLOT = "cooltools.api.saddle.saddleplot rendered to PNG"


def _module(name: str, extra: str = _EXTRA) -> Any:
    return require(name, extra)


def _retrieved() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _version() -> str | None:
    return getattr(_module("cooltools"), "__version__", None)


def _text(value: object, name: str, *, limit: int = 512) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > limit:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return value.strip()


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _source(value: str, name: str) -> str:
    """An http(s) URL, or the path of an existing local dataset table."""
    if _URL.fullmatch(value) is not None:
        return value
    local = Path(value).expanduser().resolve()
    if local.is_file():
        return str(local)
    raise ConfigurationError(f"{name} must be an http(s) URL or an existing local file: {value}")


def _new_file(value: object, name: str, *, suffix: str | None = None) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new path: {path}")
    if suffix is not None and path.suffix.lower() != suffix:
        raise ConfigurationError(f"{name} must end with {suffix}: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _upstream(name: str, call: Callable[[], Any]) -> Any:
    """Run one cooltools call; configuration errors pass through, other failures become BackendExecutionError."""
    try:
        return call()
    except ConfigurationError:
        raise
    except (
        Exception
    ) as exc:  # cooltools reports unknown datasets and checksum mismatches as plain exceptions
        raise BackendExecutionError(f"cooltools {name} failed: {exc}") from exc


def _paths(result: Any) -> list[str]:
    """cooltools returns '' for no file, one path, or a list of paths."""
    if result in (None, ""):
        return []
    if isinstance(result, (list, tuple)):
        return [str(item) for item in result]
    return [str(result)]


def _datasets(rows: Any) -> list[dict[str, Any]]:
    return [
        {
            str(key): (value if isinstance(value, (int, float)) else str(value))
            for key, value in row.items()
        }
        for row in rows
    ]


def _document(
    operation: str, service: str, request: dict[str, Any], data: dict[str, Any]
) -> dict[str, Any]:
    return {
        "status": "complete",
        "operation": operation,
        "service": service,
        "backend": "cooltools",
        "version": _version(),
        "request": request,
        "retrieved_utc": _retrieved(),
        "data": data,
    }


def download_cooltools_example_data(
    name: str = "all",
    *,
    data_dir: str,
    cache: bool = True,
    ignore_checksum: bool = False,
    info_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download one cooltools example dataset, or all of them, into a named folder (cooltools download_data).

    Needs network access unless the files are already cached in data_dir (cache true). Each file is checked
    against the md5 checksum in the dataset table unless ignore_checksum is true. info_url may be an http(s)
    URL or a local dataset table; by default cooltools' published table is used.
    """
    dataset = _text(name, "name", limit=128)
    folder = Path(_text(data_dir, "data_dir", limit=4096)).expanduser().resolve()
    if folder.exists() and not folder.is_dir():
        raise ConfigurationError(f"data_dir is not a folder: {folder}")
    folder.mkdir(parents=True, exist_ok=True)
    cache_flag = _flag(cache, "cache")
    checksum_flag = _flag(ignore_checksum, "ignore_checksum")
    io_module = _module("cooltools.lib.io")
    table = (
        io_module.URL_DATA if info_url is None else _source(_text(info_url, "info_url"), "info_url")
    )
    result = _upstream(
        "download_data",
        lambda: io_module.download_data(
            dataset,
            cache=cache_flag,
            data_dir=str(folder),
            ignore_checksum=checksum_flag,
            _url_info=table,
        ),
    )
    files = _paths(result)
    return _document(
        "download_cooltools_example_data",
        _SERVICE_DATASETS,
        {
            "name": dataset,
            "data_dir": str(folder),
            "cache": cache_flag,
            "ignore_checksum": checksum_flag,
            "info_url": table,
        },
        {"files": files, "file_count": len(files), "data_dir": str(folder)},
    )


def list_cooltools_example_datasets(
    *,
    url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the cooltools example-dataset table that print_available_datasets prints (cooltools print_available_datasets).

    Each row holds the dataset key, comment, source link, local file name and md5 checksum. The table is read
    from cooltools' published URL unless url names an http(s) URL or a local table; a local table needs no network.
    """
    io_module = _module("cooltools.lib.io")
    table = io_module.URL_DATA if url is None else _source(_text(url, "url"), "url")
    rows = _upstream("print_available_datasets", lambda: io_module._get_datasets_info(table))
    datasets = _datasets(rows)
    return _document(
        "list_cooltools_example_datasets",
        _SERVICE_DATASETS,
        {"url": table},
        {"datasets": datasets, "count": len(datasets)},
    )


def download_cooltools_file(
    url: str,
    output_path: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download one file from an http(s) URL to a new local path (cooltools download_file).

    The output path must not exist yet; the result records its size and md5 checksum. Needs network access.
    """
    source = _text(url, "url", limit=2048)
    if _URL.fullmatch(source) is None:
        raise ConfigurationError("url must be an http or https URL.")
    target = _new_file(output_path, "output_path")
    io_module = _module("cooltools.lib.io")
    _upstream("download_file", lambda: io_module.download_file(source, str(target)))
    if not target.is_file() or target.stat().st_size == 0:
        raise BackendExecutionError(f"cooltools download_file wrote no data to {target}")
    return _document(
        "download_cooltools_file",
        _SERVICE_FILE,
        {"url": source, "output_path": str(target)},
        {
            "output_path": str(target),
            "bytes": target.stat().st_size,
            "md5": _upstream("get_md5sum", lambda: io_module.get_md5sum(str(target))),
        },
    )


def _pair(value: Sequence[float], name: str) -> Sequence[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError(f"{name} must hold exactly two numbers.")
    return value


def _track_records(rows: Sequence[Sequence[object]]) -> list[list[object]]:
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence) or not rows:
        raise ConfigurationError(
            "track_rows must be a non-empty list of [chrom, start, end, value] rows."
        )
    records: list[list[object]] = []
    for index, row in enumerate(rows):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 4:
            raise ConfigurationError(
                f"track_rows[{index}] must hold four values: chrom, start, end, value."
            )
        chrom, start, end, value = row
        if not isinstance(chrom, str) or not chrom.strip():
            raise ConfigurationError(f"track_rows[{index}] chrom must be a non-empty string.")
        for label, number in (("start", start), ("end", end)):
            if isinstance(number, bool) or not isinstance(number, int) or number < 0:
                raise ConfigurationError(
                    f"track_rows[{index}] {label} must be a non-negative integer."
                )
        if end <= start:  # both are ints here
            raise ConfigurationError(f"track_rows[{index}] must have end greater than start.")
        numeric = math.nan if value is None else _number(value, f"track_rows[{index}] value")
        records.append([chrom.strip(), start, end, numeric])
    return records


def _matrix_values(rows: Sequence[Sequence[object]], n_bins: int) -> list[list[float]]:
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence) or not rows:
        raise ConfigurationError("saddle_matrix must be a non-empty square list of rows.")
    side = len(rows)
    if side != n_bins + 2:
        # cooltools.api.saddle returns the flanked matrix; without the flanks saddleplot fails with a shape mismatch.
        raise ConfigurationError(
            f"saddle_matrix must be {n_bins + 2} x {n_bins + 2}, with the two outlier flanks that "
            f"cooltools.api.saddle returns; got {side} x {side} for {n_bins} bins."
        )
    values: list[list[float]] = []
    for index, row in enumerate(rows):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != side:
            raise ConfigurationError("saddle_matrix must be square.")
        values.append(
            [math.nan if item is None else _number(item, f"saddle_matrix[{index}]") for item in row]
        )
    return values


def plot_cooltools_saddle(
    track_rows: Sequence[Sequence[object]],
    saddle_matrix: Sequence[Sequence[object]],
    n_bins: int,
    output_png: str,
    *,
    vrange: Sequence[float] | None = None,
    qrange: Sequence[float] | None = (0.0, 1.0),
    scale: str = "log",
    vmin: float = 0.5,
    vmax: float = 2.0,
    cmap: str = "coolwarm",
    title: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a cooltools saddle plot of a track against a saddle matrix into a new PNG (cooltools saddleplot).

    track_rows are [chrom, start, end, value] rows; saddle_matrix is (n_bins + 2) square and includes the two
    outlier flanks that cooltools.api.saddle returns. Exactly one of vrange (value limits) or qrange (quantile limits) must be given:
    pass qrange=None when using vrange. Returns the PNG path and the digitised track counts and heatmap values.
    """
    records = _track_records(track_rows)
    if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 1:
        raise ConfigurationError("n_bins must be a positive integer.")
    matrix = _matrix_values(saddle_matrix, n_bins)
    if (vrange is None) == (qrange is None):
        raise ConfigurationError(
            "Supply exactly one of vrange or qrange (pass qrange=None with vrange)."
        )
    value_range = (
        None if vrange is None else tuple(_number(v, "vrange") for v in _pair(vrange, "vrange"))
    )
    quantile_range = (
        None if qrange is None else tuple(_number(v, "qrange") for v in _pair(qrange, "qrange"))
    )
    if quantile_range is not None and not (0.0 <= quantile_range[0] <= quantile_range[1] <= 1.0):
        raise ConfigurationError("qrange must satisfy 0 <= low <= high <= 1.")
    if value_range is not None and value_range[0] >= value_range[1]:
        raise ConfigurationError("vrange must satisfy low < high.")
    if scale not in ("log", "linear"):
        raise ConfigurationError("scale must be 'log' or 'linear'.")
    limit_low = _number(vmin, "vmin")
    limit_high = _number(vmax, "vmax")
    if limit_low >= limit_high or (scale == "log" and limit_low <= 0):
        raise ConfigurationError("vmin must be below vmax, and above zero for a log scale.")
    colour = _text(cmap, "cmap", limit=64)
    heading = None if title is None else _text(title, "title", limit=200)
    target = _new_file(output_png, "output_png", suffix=".png")

    pandas = _module("pandas")
    numpy = _module("numpy")
    saddle_api = _module("cooltools.api.saddle")
    plt = _module("matplotlib.pyplot", _PLOT_EXTRA)
    track = pandas.DataFrame(records, columns=["chrom", "start", "end", "value"])
    saddle = numpy.array(matrix, dtype=float)
    digitized, bin_edges = saddle_api.digitize(
        track, n_bins, vrange=value_range, qrange=quantile_range
    )
    codes = digitized.iloc[:, 3].astype("int64").to_numpy()

    figure = None
    open_before = set(plt.get_fignums())
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            grid = saddle_api.saddleplot(
                track,
                saddle,
                n_bins,
                vrange=value_range,
                qrange=quantile_range,
                cmap=colour,
                scale=scale,
                vmin=limit_low,
                vmax=limit_high,
                title=heading,
            )
            figure = grid["ax_heatmap"].figure
            figure.savefig(str(target), format="png", dpi=_SAVE_DPI)
    except Exception as exc:  # upstream draws with pyplot; report any failure as a backend error
        raise BackendExecutionError(f"cooltools saddleplot failed: {exc}") from exc
    finally:
        # saddleplot leaves its figure open when it fails part-way, so close every figure this call opened.
        if figure is not None:
            plt.close(figure)
        for number in set(plt.get_fignums()) - open_before:
            plt.close(number)

    heat = saddle[1:-1, 1:-1]  # the plotted cells, without the outlier flanks
    return {
        "status": "completed",
        "operation": "plot_cooltools_saddle",
        "service": _SERVICE_PLOT,
        "backend": "cooltools",
        "version": _version(),
        "image_path": str(target),
        "data": {
            "n_bins": n_bins,
            "vrange": None if value_range is None else list(value_range),
            "qrange": None if quantile_range is None else list(quantile_range),
            "scale": scale,
            "vmin": limit_low,
            "vmax": limit_high,
            "cmap": colour,
            "track_rows": len(records),
            "bin_edges": [float(edge) for edge in bin_edges],
            "bin_counts": [int((codes == code).sum()) for code in range(1, n_bins + 1)],
            "left_outlier_count": int((codes == 0).sum()),
            "right_outlier_count": int((codes == n_bins + 1).sum()),
            "missing_count": int((codes == -1).sum()),
            "heatmap_shape": [int(heat.shape[0]), int(heat.shape[1])],
            "heatmap_values": [
                [None if not math.isfinite(item) else float(item) for item in row]
                for row in heat.tolist()
            ],
        },
        "upstream_warnings": [str(item.message) for item in caught],
    }


__all__ = [
    "download_cooltools_example_data",
    "download_cooltools_file",
    "list_cooltools_example_datasets",
    "plot_cooltools_saddle",
]

_register_public_verbose(globals())
