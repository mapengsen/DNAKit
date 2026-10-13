"""Genome tiling, overlap counting, table-text parsing, and Jaccard similarity through pyranges1.

Functions read local BED files or in-memory mappings and return JSON-compatible records. Coordinates
are 0-based half-open (BED convention). pyranges1 is imported only when a function runs; a missing
package raises BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _integer,
    _lengths,
    _package,
    _pandas_records,
    _rows_limit,
    _scan_bed,
    _text,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_STRAND_BEHAVIOR = ("auto", "same", "opposite", "ignore")


def _ranges_records(
    ranges: Any, limit: int, what: str, *, value_key: str | None = None
) -> list[dict[str, Any]]:
    records = []
    for row in _pandas_records(ranges, limit, what):
        record: dict[str, Any] = {
            "chrom": str(row["Chromosome"]),
            "start": int(row["Start"]),
            "end": int(row["End"]),
        }
        if "Strand" in row and row["Strand"] is not None:
            record["strand"] = str(row["Strand"])
        if value_key is not None:
            record[value_key] = row[value_key]
        records.append(record)
    return sorted(records, key=lambda r: (r["chrom"], r["start"], r["end"]))


def tile_pyranges1_genome(
    chromsizes: Mapping[str, int],
    *,
    tile_size: int,
    full_last_tile: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tile each chromosome of a size mapping into fixed-width intervals with pyranges1."""

    lengths = _lengths(chromsizes, "chromsizes")
    width = _integer(tile_size, "tile_size", 1, 1_000_000_000)
    if not isinstance(full_last_tile, bool):
        raise ConfigurationError("full_last_tile must be a boolean.")
    limit = _rows_limit(max_rows)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    tiles = pyranges1.tile_genome(lengths, width, full_last_tile=full_last_tile)
    records = _ranges_records(tiles, limit, "genome tiles")
    return {
        "status": "completed",
        "operation": "tile_pyranges1_genome",
        "tile_size": width,
        "tiles": records,
        "tile_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "tiles are 0-based half-open and start at 0 on each chromosome",
        },
    }


def count_pyranges1_overlaps(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the intervals of another BED file that overlap each source interval with pyranges1."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    left = pyranges1.read_bed(str(path))
    right = pyranges1.read_bed(str(other_path))
    counted = pyranges1.count_overlaps({"overlaps": right}, left)
    records = _ranges_records(counted, limit, "overlap counts", value_key="overlaps")
    return {
        "status": "completed",
        "operation": "count_pyranges1_overlaps",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "overlaps counts intervals of the other file for each source interval",
        },
    }


def parse_pyranges1_string(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a headed multi-line table (Chromosome, Start, End, optional Strand) with pyranges1."""

    body = _text(text, "text", maximum=65_536)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    try:
        ranges = pyranges1.from_string(body)
    except Exception as exc:  # pyranges1 raises ValueError for missing columns or bad numbers
        raise ConfigurationError(f"pyranges1 could not parse the text: {exc}") from exc
    records = _ranges_records(ranges, _DEFAULT_ROWS, "parsed intervals")
    return {
        "status": "completed",
        "operation": "parse_pyranges1_string",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "the first line is the header; Start and End are taken as written",
        },
    }


def jaccard_pyranges1_intervals(
    source: str | Path,
    other: str | Path,
    *,
    strand_behavior: str = "auto",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the base-pair Jaccard index of two BED files with pyranges1.stats.jaccard."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    if strand_behavior not in _STRAND_BEHAVIOR:
        raise ConfigurationError(f"strand_behavior must be one of {', '.join(_STRAND_BEHAVIOR)}.")
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    stats = require("pyranges1.stats", "bio-pyranges1")
    left = pyranges1.read_bed(str(path))
    right = pyranges1.read_bed(str(other_path))
    value = stats.jaccard(left, right, strand_behavior=strand_behavior)
    return {
        "status": "completed",
        "operation": "jaccard_pyranges1_intervals",
        "jaccard": float(value),
        "strand_behavior": strand_behavior,
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "shared bases divided by bases covered by either file",
        },
    }


__all__ = [
    "count_pyranges1_overlaps",
    "jaccard_pyranges1_intervals",
    "parse_pyranges1_string",
    "tile_pyranges1_genome",
]

_register_public_verbose(globals(), tuple(__all__))
