"""Cooler file description, chromosome sizes, region plots, pixel dumps, triplet loading and pair sorting.

Each public function imports cooler (or pysam, pandas, matplotlib) only when it is called.
Missing dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Inputs are checked before upstream calls; upstream ``ValueError``, ``KeyError`` and
``OSError`` are reported as ``ConfigurationError``.

Functions that write a file create a new path and never overwrite one. Results are plain
dictionaries; NaN and infinite values are returned as ``None``.
"""

from __future__ import annotations

import math
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_PYSAM_EXTRA = "bio-reads"
_PLOT_EXTRA = "bio-plots"
_MAX_CHROMOSOMES = 10_000
_MAX_CHROM_LENGTH = 2_000_000_000
_MAX_BIN_SIZE = 10_000_000
_MAX_BINS = 1_000_000
_MAX_REGION_BINS = 300
_MAX_DUMP_ROWS = 100_000
_MAX_ROWS = 10_000_000
_SCALES = ("linear", "log2", "log10")
_FORMATS = ("coo", "bg2")
_INDEXES = ("none", "tabix")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _path(value: object, name: str) -> Path:
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


def _chromsizes(value: object) -> dict[str, int]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CHROMOSOMES:
        raise ConfigurationError(
            f"chromsizes must map 1 to {_MAX_CHROMOSOMES} chromosome names to lengths."
        )
    sizes: dict[str, int] = {}
    for name, length in value.items():
        _text(name, "chromosome name")
        if name in sizes:
            raise ConfigurationError(f"chromsizes repeats the chromosome {name!r}.")
        sizes[name] = _integer(length, f"length of {name}", 1, _MAX_CHROM_LENGTH)
    return sizes


def _plain(value: Any) -> Any:
    """Convert numpy and bytes values to JSON-friendly Python values; NaN and inf become None."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not hasattr(value, "item"):
        return _plain(value.tolist())
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _records(frame: Any) -> list[dict[str, Any]]:
    return [
        {str(key): _plain(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _open_cooler(cooler: Any, path: Path) -> Any:
    try:
        return cooler.Cooler(str(path))
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {path} ({exc})") from exc


def _weighted(clr: Any, balanced: bool) -> None:
    if balanced and "weight" not in clr.bins()[:1].columns:
        raise ConfigurationError(
            "This cooler has no balancing weights; run balance_contact_matrix() first."
        )


def _extent_check(clr: Any, text: str, limit: int) -> None:
    try:
        lo, hi = clr.extent(text)
    except (KeyError, ValueError) as exc:
        raise ConfigurationError(f"Region is not valid for this cooler: {text} ({exc})") from exc
    if hi - lo > limit:
        raise ConfigurationError(
            f"Region spans {hi - lo} bins, above max_region_bins={limit}; choose a smaller region."
        )


# ---------------------------------------------------------------- describe


def describe_cooler(cool_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Return the info block, binsize and chromosome count of a cooler, using cooler (cooler info)."""

    source = _path(cool_path, "cool_path")
    cooler = require("cooler", _COOLER_EXTRA)
    util = require("cooler.util", _COOLER_EXTRA)
    clr = _open_cooler(cooler, source)
    try:
        info = dict(clr.info)
        binsize = int(clr.binsize)
        chromosomes = list(clr.chromnames)
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {source} ({exc})") from exc
    metadata = info.pop("metadata", None)
    return {
        "status": "completed",
        "operation": "describe_cooler",
        "cool_path": str(source),
        "binsize": binsize,
        "bin_count": int(info["nbins"]),
        "chromosome_count": len(chromosomes),
        "chromosomes": [str(name) for name in chromosomes],
        "info": _plain(util.attrs_to_jsonable(info)),
        "metadata": None if metadata is None else _plain(util.attrs_to_jsonable(dict(metadata))),
        "provenance": _provenance(
            "cooler",
            "info block of the cooler file (cooler.Cooler.info) without the metadata block",
            note="metadata, when present, is returned separately",
        ),
    }


# ---------------------------------------------------------------- chromosome sizes


def read_cooler_chromsizes(
    cool_path: str | Path,
    *,
    from_bins: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the chromosome lengths of a cooler from its chroms table, or inferred from its bins."""

    source = _path(cool_path, "cool_path")
    flag = _flag(from_bins, "from_bins")
    cooler = require("cooler", _COOLER_EXTRA)
    util = require("cooler.util", _COOLER_EXTRA)
    clr = _open_cooler(cooler, source)
    try:
        stored = [(str(name), int(length)) for name, length in clr.chromsizes.items()]
        binsize = int(clr.binsize)
        bins = clr.bins()[:] if flag else None
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {source} ({exc})") from exc
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "read_cooler_chromsizes",
        "cool_path": str(source),
        "source": "bins" if flag else "chroms",
        "binsize": binsize,
        "chromosome_count": len(stored),
        "chromosomes": [{"name": name, "length": length} for name, length in stored],
    }
    if bins is not None:
        inferred = util.get_chromsizes(bins)
        inferred_rows = [(str(name), int(length)) for name, length in inferred.items()]
        inferred_binsize = util.get_binsize(bins)
        result["chromosomes"] = [{"name": name, "length": length} for name, length in inferred_rows]
        result["chromosome_count"] = len(inferred_rows)
        result["bin_binsize"] = None if inferred_binsize is None else int(inferred_binsize)
        result["uniform_bins"] = inferred_binsize is not None
        result["matches_chroms_table"] = inferred_rows == stored
    result["provenance"] = _provenance(
        "cooler",
        "chroms table of the cooler (cooler.Cooler.chromsizes)"
        if bins is None
        else "last bin end per chromosome of the bins table (cooler.util.get_chromsizes)",
        note="bin_binsize is cooler.util.get_binsize; it is null when the bins are non-uniform",
    )
    return result


# ---------------------------------------------------------------- plot


def plot_contact_matrix_region(
    cool_path: str | Path,
    region: str,
    output_path: str | Path,
    *,
    region2: str | None = None,
    balanced: bool = False,
    field: str = "count",
    scale: str = "log10",
    cmap: str = "YlOrRd",
    zmin: float | None = None,
    zmax: float | None = None,
    dpi: int = 150,
    max_region_bins: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Render a contact-matrix region of a cooler to a PNG (cooler show --out) and return its values."""

    source = _path(cool_path, "cool_path")
    row_text = _text(region, "region")
    col_text = row_text if region2 is None else _text(region2, "region2")
    target = _new_path(output_path, "output_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    use_weights = _flag(balanced, "balanced")
    if scale not in _SCALES:
        raise ConfigurationError(f"scale must be one of {_SCALES}.")
    if not isinstance(field, str) or not field.strip():
        raise ConfigurationError("field must be a non-empty column name.")
    if cmap not in _colormaps():
        raise ConfigurationError(f"cmap is not a matplotlib colormap name: {cmap}")
    if zmin is not None and (isinstance(zmin, bool) or not isinstance(zmin, (int, float))):
        raise ConfigurationError("zmin must be a number or None.")
    if zmax is not None and (isinstance(zmax, bool) or not isinstance(zmax, (int, float))):
        raise ConfigurationError("zmax must be a number or None.")
    pixel_dpi = _integer(dpi, "dpi", 50, 600)
    limit = _integer(max_region_bins, "max_region_bins", 1, _MAX_REGION_BINS)
    cooler = require("cooler", _COOLER_EXTRA)
    np = require("numpy", "bio-analysis")
    util = require("cooler.util", _COOLER_EXTRA)
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    figure_module = require("matplotlib.figure", _PLOT_EXTRA)
    canvas_module = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    clr = _open_cooler(cooler, source)
    _weighted(clr, use_weights)
    _extent_check(clr, row_text, limit)
    _extent_check(clr, col_text, limit)
    try:
        row_chrom, row_lo, row_hi = util.parse_region(row_text, clr.chromsizes)
        col_chrom, col_lo, col_hi = util.parse_region(col_text, clr.chromsizes)
        matrix = np.asarray(clr.matrix(balance=use_weights, field=field).fetch(row_text, col_text))
    except (KeyError, ValueError) as exc:
        raise ConfigurationError(f"Region is not valid for this cooler: {exc}") from exc
    values = matrix.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        if scale == "log2":
            values = np.log2(values)
        elif scale == "log10":
            values = np.log10(values)
    figure = figure_module.Figure(figsize=(11, 10))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    image = axes.imshow(
        values,
        interpolation="none",
        extent=[col_lo, col_hi, row_hi, row_lo],
        vmin=zmin,
        vmax=zmax,
        cmap=cmap,
    )
    axes.set_ylabel(f"{row_chrom} coordinate")
    axes.set_xlabel(f"{col_chrom} coordinate")
    colorbar = figure.colorbar(image, ax=axes)
    colorbar.set_label(
        {
            "linear": "relative contact frequency",
            "log2": "log 2 ( relative contact frequency )",
            "log10": "log 10 ( relative contact frequency )",
        }[scale]
    )
    figure.savefig(str(target), dpi=pixel_dpi, format="png")
    del matplotlib
    return {
        "status": "completed",
        "operation": "plot_contact_matrix_region",
        "cool_path": str(source),
        "region": row_text,
        "region2": col_text,
        "output_path": str(target),
        "output_bytes": target.stat().st_size,
        "binsize": int(clr.binsize),
        "balanced": use_weights,
        "field": field,
        "scale": scale,
        "shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "matrix": [[_plain(float(value)) for value in row] for row in matrix],
        "provenance": _provenance(
            "cooler",
            "matrix fetched with cooler and drawn with matplotlib imshow (cooler show --out)",
            scale="values are log-transformed for the image only; matrix is the raw fetch",
            coordinate_system="0-based half-open bins; axis extent in genomic coordinates",
        ),
    }


def _colormaps() -> set[str]:
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    return set(matplotlib.colormaps)


# ---------------------------------------------------------------- pixel dump


def dump_cooler_pixels(
    cool_path: str | Path,
    *,
    region: str | None = None,
    region2: str | None = None,
    balanced: bool = False,
    join: bool = False,
    limit: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dump stored pixels, optionally restricted to a region, with balanced values and bin annotation."""

    source = _path(cool_path, "cool_path")
    use_weights = _flag(balanced, "balanced")
    use_join = _flag(join, "join")
    rows_limit = _integer(limit, "limit", 1, _MAX_DUMP_ROWS)
    row_text = None if region is None else _text(region, "region")
    col_text = None if region2 is None else _text(region2, "region2")
    if col_text is not None and row_text is None:
        raise ConfigurationError("region2 needs region.")
    cooler = require("cooler", _COOLER_EXTRA)
    clr = _open_cooler(cooler, source)
    _weighted(clr, use_weights)
    if row_text is not None:
        _extent_check(clr, row_text, _MAX_REGION_BINS * 1_000)
        if col_text is not None:
            _extent_check(clr, col_text, _MAX_REGION_BINS * 1_000)
    try:
        selector = clr.matrix(balance=use_weights, as_pixels=True, join=use_join)
        if row_text is None:
            frame = selector[:]
        else:
            frame = selector.fetch(row_text) if col_text is None else selector.fetch(row_text, col_text)
        total = len(frame)
    except (KeyError, ValueError) as exc:
        raise ConfigurationError(f"Pixels could not be read for this region: {exc}") from exc
    shown = frame.iloc[:rows_limit]
    return {
        "status": "completed",
        "operation": "dump_cooler_pixels",
        "cool_path": str(source),
        "region": row_text,
        "region2": col_text,
        "balanced": use_weights,
        "join": use_join,
        "pixel_count": int(total),
        "returned": len(shown),
        "columns": [str(column) for column in shown.columns],
        "pixels": _records(shown),
        "provenance": _provenance(
            "cooler",
            "stored pixels in bin order via cooler.Cooler.matrix(as_pixels=True)",
            balanced="balanced = weight1 * weight2 * count when balanced is true",
            order="rows follow the stored bin1_id, bin2_id order; limit keeps the first rows",
        ),
    }


# ---------------------------------------------------------------- triplet loading


def load_contact_triplets(
    records_path: str | Path,
    cool_path: str | Path,
    *,
    chromsizes: Mapping[str, int],
    binsize: int,
    file_format: str = "coo",
    one_based: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a bin-pair triplet (coo) or 2D bedGraph (bg2) contact file into a cooler, using cooler load."""

    source = _path(records_path, "records_path")
    output = _new_path(cool_path, "cool_path")
    sizes = _chromsizes(chromsizes)
    size = _integer(binsize, "binsize", 1, _MAX_BIN_SIZE)
    if file_format not in _FORMATS:
        raise ConfigurationError(f"file_format must be one of {_FORMATS}.")
    shifted = _flag(one_based, "one_based")
    layout = sum(-(-length // size) for length in sizes.values())
    if layout > _MAX_BINS:
        raise ConfigurationError(f"The bin layout has {layout} bins, above the limit {_MAX_BINS}.")
    cooler = require("cooler", _COOLER_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    bins = cooler.binnify(pandas.Series(sizes, dtype="int64"), size)
    if file_format == "coo":
        columns = ["bin1_id", "bin2_id", "count"]
        dtypes = {"bin1_id": "int64", "bin2_id": "int64", "count": "int64"}
        pipeline = cooler.create.sanitize_pixels(
            bins, is_one_based=shifted, tril_action="reflect", sort=True
        )
    else:
        columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]
        dtypes = {
            "chrom1": str,
            "start1": "int64",
            "end1": "int64",
            "chrom2": str,
            "start2": "int64",
            "end2": "int64",
            "count": "int64",
        }
        pipeline = cooler.create.sanitize_records(
            bins, schema="bg2", is_one_based=shifted, tril_action="reflect", sort=True
        )
    try:
        with pandas.read_csv(
            source,
            sep="\t",
            comment="#",
            header=None,
            names=columns,
            dtype=dtypes,
            iterator=True,
            chunksize=1_000_000,
        ) as reader:
            cooler.create.create_from_unordered(
                str(output),
                bins,
                map(pipeline, reader),
                columns=["count"],
                dtypes={},
                ensure_sorted=False,
                symmetric_upper=True,
                triucheck=True,
                temp_dir=str(output.parent),
                delete_temp=True,
            )
        info = cooler.Cooler(str(output)).info
    except (ValueError, IndexError, KeyError, OSError) as exc:
        output.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler rejected the {file_format} input: {exc}") from exc
    return {
        "status": "completed",
        "operation": "load_contact_triplets",
        "records": str(source),
        "cool_path": str(output),
        "file_format": file_format,
        "one_based": shifted,
        "binsize": size,
        "bins": layout,
        "nonzero_pixels": int(info["nnz"]),
        "total_contacts": int(info["sum"]),
        "provenance": _provenance(
            "cooler",
            "records merged by bin pair; lower-triangle records reflected; duplicates summed "
            "(cooler.create.create_from_unordered with sanitize_pixels or sanitize_records)",
            coordinate_system="coo: bin ids; bg2: genomic start1/start2 are assigned to bins",
        ),
    }


# ---------------------------------------------------------------- pair sorting


def sort_pairs_file(
    pairs_path: str | Path,
    output_path: str | Path,
    *,
    chromsizes: Mapping[str, int],
    index: str = "none",
    verbose: bool = True,
) -> dict[str, Any]:
    """Orient each pair as upper triangular, sort by position and optionally index (cooler csort, Python)."""

    source = _path(pairs_path, "pairs_path")
    sizes = _chromsizes(chromsizes)
    mode = _text(index, "index")
    if mode not in _INDEXES:
        raise ConfigurationError(
            "index must be 'none' or 'tabix'; cooler's default pairix index needs the pairix "
            "program, which DNAKit does not run."
        )
    target = _new_path(output_path, "output_path")
    if mode == "tabix" and not target.name.endswith(".gz"):
        raise ConfigurationError("output_path must end with .gz when index is 'tabix'.")
    rank = {name: position for position, name in enumerate(sizes, 1)}
    kept: list[tuple[str, int, str, int, str]] = []
    seen = dropped = flipped = 0
    try:
        with source.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                seen += 1
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    raise ConfigurationError(f"A pair line has fewer than 5 columns: {line.strip()!r}")
                chrom1, chrom2 = parts[1], parts[3]
                if chrom1 not in rank or chrom2 not in rank:
                    dropped += 1
                    continue
                pos1, pos2 = int(parts[2]), int(parts[4])
                if rank[chrom1] > rank[chrom2] or (rank[chrom1] == rank[chrom2] and pos1 > pos2):
                    parts[1], parts[3] = parts[3], parts[1]
                    parts[2], parts[4] = parts[4], parts[2]
                    chrom1, chrom2, pos1, pos2 = parts[1], parts[3], int(parts[2]), int(parts[4])
                    flipped += 1
                kept.append((chrom1, pos1, chrom2, pos2, "\t".join(parts)))
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable pairs file: {source} ({exc})") from exc
    kept.sort(key=lambda row: (row[0].encode(), row[1], row[2].encode(), row[3], row[4].encode()))
    index_file: str | None = None
    if mode == "none":
        with target.open("w", encoding="utf-8") as handle:
            for row in kept:
                handle.write(row[4] + "\n")
    else:
        pysam = require("pysam", _PYSAM_EXTRA)
        with tempfile.TemporaryDirectory(dir=str(target.parent)) as scratch:
            plain = Path(scratch) / "sorted.pairs"
            with plain.open("w", encoding="utf-8") as handle:
                for row in kept:
                    handle.write(row[4] + "\n")
            pysam.tabix_compress(str(plain), str(target), force=False)
        try:
            pysam.tabix_index(
                str(target),
                seq_col=1,
                start_col=2,
                end_col=2,
                meta_char="#",
                zerobased=False,
                force=False,
            )
        except (OSError, ValueError) as exc:
            target.unlink(missing_ok=True)
            raise ConfigurationError(f"tabix could not index the sorted pairs: {exc}") from exc
        index_file = str(target) + ".tbi"
    return {
        "status": "completed",
        "operation": "sort_pairs_file",
        "pairs": str(source),
        "output": str(target),
        "index": mode,
        "index_file": index_file,
        "chromosome_count": len(sizes),
        "records_read": seen,
        "records_dropped": dropped,
        "records_flipped": flipped,
        "records_written": len(kept),
        "provenance": _provenance(
            "cooler",
            "mates flipped to upper triangle by the chromosome order, then sorted by chrom1, pos1, "
            "chrom2, pos2 in byte order (cooler csort, tabix sort preset)",
            index="index=tabix uses pysam.tabix_index (cols 2 and 3); cooler's pairix preset is not run",
            header="comment lines starting with '#' are not written, as in cooler csort",
        ),
    }


__all__ = [
    "describe_cooler",
    "dump_cooler_pixels",
    "load_contact_triplets",
    "plot_contact_matrix_region",
    "read_cooler_chromsizes",
    "sort_pairs_file",
]

_register_public_verbose(globals(), tuple(__all__))
