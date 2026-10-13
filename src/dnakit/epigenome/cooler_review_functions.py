"""Single-cell pixel tables, group listings and bin-size inference from cooler.

Each public function imports cooler only when it is called. Missing dependencies raise
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint. Inputs are checked
before upstream calls; upstream ``ValueError``, ``KeyError`` and ``OSError`` are reported as
``ConfigurationError``.

Output files are created as new paths and are never overwritten.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_MAX_CHROMOSOMES = 10_000
_MAX_CHROM_LENGTH = 2_000_000_000
_MAX_BINSIZE = 10_000_000
_MAX_BINS = 200_000
_MAX_CELLS = 256
_MAX_PAIRS_PER_CELL = 5_000_000
_CELL_NAME = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new(value: object, name: str) -> Path:
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
    for key, length in value.items():
        name = _text(key, "chromsizes chromosome")
        if name in sizes:
            raise ConfigurationError(f"chromsizes repeats the chromosome {name!r}.")
        sizes[name] = _integer(length, f"length of {name}", 1, _MAX_CHROM_LENGTH)
    return sizes


def _cells(value: object) -> dict[str, Path]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CELLS:
        raise ConfigurationError(f"pairs_by_cell must map 1 to {_MAX_CELLS} cell names to files.")
    cells: dict[str, Path] = {}
    for key, path in value.items():
        if not isinstance(key, str) or not _CELL_NAME.match(key):
            raise ConfigurationError(
                "Cell names must be 1 to 64 letters, digits, '_', '.' or '-' characters."
            )
        cells[key] = _existing(path, f"pairs_by_cell[{key!r}]")
    return cells


def _layout(sizes: Mapping[str, int], binsize: int) -> tuple[dict[str, int], int]:
    """Bin offset of each chromosome and the total bin count (last bin may be shorter)."""

    offsets: dict[str, int] = {}
    total = 0
    for name, length in sizes.items():
        offsets[name] = total
        total += -(-length // binsize)
    if total > _MAX_BINS:
        raise ConfigurationError(
            f"The layout has {total} bins, above the limit of {_MAX_BINS}; use a larger binsize."
        )
    return offsets, total


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _read_pairs(path: Path) -> Any:
    pandas = require("pandas", _COOLER_EXTRA)
    try:
        frame = pandas.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            usecols=[1, 2, 3, 4],
            names=["chrom1", "pos1", "chrom2", "pos2"],
            dtype={"chrom1": str, "chrom2": str, "pos1": "int64", "pos2": "int64"},
            nrows=_MAX_PAIRS_PER_CELL + 1,
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(
            f"Input is not a readable .pairs file with at least 5 columns: {path} ({exc})"
        ) from exc
    if len(frame) > _MAX_PAIRS_PER_CELL:
        raise ConfigurationError(f"A cell has more than {_MAX_PAIRS_PER_CELL} pairs.")
    if frame.empty:
        raise ConfigurationError(f"The .pairs file has no pairs: {path}")
    return frame


def _bin_ids(
    frame: Any, chrom: str, pos: str, sizes: Mapping[str, int], offsets: Any, size: int
) -> Any:
    """Global bin index of each (chrom, pos) row; positions must lie inside the chromosome."""

    unknown = ~frame[chrom].isin(list(sizes))
    if unknown.any():
        raise ConfigurationError(
            f"A pair uses a chromosome missing from chromsizes: {frame.loc[unknown, chrom].iloc[0]}"
        )
    lengths = frame[chrom].map(sizes)
    if ((frame[pos] < 0) | (frame[pos] >= lengths)).any():
        raise ConfigurationError("A pair position lies outside its chromosome.")
    return frame[chrom].map(offsets) + frame[pos] // size


def create_single_cell_cooler(
    output_path: str | Path,
    pairs_by_cell: Mapping[str, str | Path],
    chromsizes: Mapping[str, int],
    binsize: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a single-cell (.scool) cooler from one .pairs file per cell (cooler create_scool)."""

    target = _new(output_path, "output_path")
    sizes = _chromsizes(chromsizes)
    size = _integer(binsize, "binsize", 1, _MAX_BINSIZE)
    cells = _cells(pairs_by_cell)
    offsets, bin_count = _layout(sizes, size)
    pandas = require("pandas", _COOLER_EXTRA)
    cooler = require("cooler", _COOLER_EXTRA)
    bins = cooler.binnify(pandas.Series(sizes, dtype="int64"), size)
    bin_offsets = pandas.Series(offsets, dtype="int64")
    pixel_tables: dict[str, Any] = {}
    summary: list[dict[str, Any]] = []
    for name in sorted(cells):
        pairs = _read_pairs(cells[name])
        first = _bin_ids(pairs, "chrom1", "pos1", sizes, bin_offsets, size)
        second = _bin_ids(pairs, "chrom2", "pos2", sizes, bin_offsets, size)
        lo = first.where(first <= second, second).astype("int64")
        hi = first.where(first > second, second).astype("int64")
        table = (
            pandas.DataFrame({"bin1_id": lo, "bin2_id": hi})
            .groupby(["bin1_id", "bin2_id"])
            .size()
            .reset_index(name="count")
            .astype({"count": "int64"})
        )
        pixel_tables[name] = table
        summary.append({"name": name, "pair_count": len(pairs), "nonzero_pixels": len(table)})
    try:
        cooler.create_scool(str(target), bins, pixel_tables)
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"cooler could not write the single-cell file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "create_single_cell_cooler",
        "output_path": str(target),
        "binsize": size,
        "bin_count": bin_count,
        "cell_count": len(summary),
        "cells": summary,
        "provenance": _provenance(
            "cooler",
            "per-cell pixel counts from .pairs files, binned and stored with cooler.create_scool",
            pixels="upper triangle; counts are summed per bin pair",
        ),
    }


def list_cooler_groups(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """List the groups that hold a cooler data collection in an HDF5 file, using cooler."""

    source = _existing(path, "path")
    api = require("cooler.api", _COOLER_EXTRA)
    try:
        groups = sorted(str(group) for group in api.list_coolers(str(source)))
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable HDF5 file: {source} ({exc})") from exc
    return {
        "status": "completed",
        "operation": "list_cooler_groups",
        "path": str(source),
        "group_count": len(groups),
        "groups": groups,
        "provenance": _provenance(
            "cooler",
            "groups that hold a cooler data collection (cooler.api.list_coolers)",
            order="sorted by group path",
        ),
    }


def list_scool_cells(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """List the cells of a single-cell (.scool) cooler file, using cooler."""

    source = _existing(path, "path")
    fileops = require("cooler.fileops", _COOLER_EXTRA)
    try:
        is_scool = fileops.is_scool_file(str(source))
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable HDF5 file: {source} ({exc})") from exc
    if not is_scool:
        raise ConfigurationError(f"Input is not a single-cell cooler file: {source}")
    try:
        groups = sorted(str(group) for group in fileops.list_scool_cells(str(source)))
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Could not list the cells of {source} ({exc})") from exc
    cells = [{"name": group.rsplit("/", 1)[-1], "group": group} for group in groups]
    return {
        "status": "completed",
        "operation": "list_scool_cells",
        "path": str(source),
        "cell_count": len(cells),
        "cells": cells,
        "provenance": _provenance(
            "cooler",
            "cell groups of a single-cell file (cooler.fileops.list_scool_cells)",
            order="sorted by group path",
        ),
    }


def infer_cooler_binsize(cool_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Infer the bin size of a cooler from its bin table, using cooler.util.get_binsize."""

    source = _existing(cool_path, "cool_path")
    cooler = require("cooler", _COOLER_EXTRA)
    util = require("cooler.util", _COOLER_EXTRA)
    try:
        bins = cooler.Cooler(str(source)).bins()[:]
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {source} ({exc})") from exc
    size = util.get_binsize(bins)
    return {
        "status": "completed",
        "operation": "infer_cooler_binsize",
        "cool_path": str(source),
        "bin_count": len(bins),
        "uniform": size is not None,
        "binsize": None if size is None else int(size),
        "provenance": _provenance(
            "cooler",
            "bin size shared by all bins except the last bin of each chromosome "
            "(cooler.util.get_binsize)",
            non_uniform="binsize is None when the bins have different sizes",
        ),
    }


__all__ = [
    "create_single_cell_cooler",
    "infer_cooler_binsize",
    "list_cooler_groups",
    "list_scool_cells",
]

_register_public_verbose(globals(), tuple(__all__))
