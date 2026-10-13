"""Interval operations on BED files and region strings through bioframe.

Each public function reads local BED3 to BED6 files or an in-memory chromosome-size mapping,
calls one bioframe operation, and returns JSON-compatible records. Coordinates are 0-based
half-open unless a docstring says otherwise. bioframe is imported only when a function runs;
a missing package raises BackendUnavailableError with the install hint.
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

_HOW = ("inner", "left", "right", "outer")


def _read_bed(bioframe: Any, path: Path, limit: int) -> Any:
    return bioframe.read_table(str(path), schema=f"bed{_scan_bed(path, limit)}")


def overlap_bioframe_intervals(
    source: str | Path,
    other: str | Path,
    *,
    how: str = "inner",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pair every overlapping interval of two BED files with bioframe.overlap."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    if how not in _HOW:
        raise ConfigurationError(f"how must be one of {', '.join(_HOW)}.")
    bioframe = require("bioframe", "bio-bioframe")
    pairs = bioframe.overlap(
        _read_bed(bioframe, path, limit), _read_bed(bioframe, other_path, limit), how=how
    )
    records = _pandas_records(pairs, limit, "Overlap pairs", ordered=True)
    return {
        "status": "completed",
        "operation": "overlap_bioframe_intervals",
        "how": how,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "0-based half-open; suffix _ marks columns from the other file",
        },
    }


def subtract_bioframe_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the bases covered by another BED file from each interval with bioframe.subtract."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    bioframe = require("bioframe", "bio-bioframe")
    remainder = bioframe.subtract(
        _read_bed(bioframe, path, limit), _read_bed(bioframe, other_path, limit)
    )
    intervals = _pandas_records(remainder, limit, "Remaining intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "subtract_bioframe_intervals",
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "0-based half-open; only the pieces left after removing other remain",
        },
    }


def expand_bioframe_intervals(
    source: str | Path,
    *,
    pad: int,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Grow each BED interval by pad bases on both sides with bioframe.expand."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    amount = _integer(pad, "pad", 0, 1_000_000_000)
    bioframe = require("bioframe", "bio-bioframe")
    expanded = bioframe.expand(_read_bed(bioframe, path, limit), pad=amount)
    intervals = _pandas_records(expanded, limit, "Expanded intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "expand_bioframe_intervals",
        "pad": amount,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "start moves left by pad and end moves right by pad; no clipping at 0",
        },
    }


def cluster_bioframe_intervals(
    source: str | Path,
    *,
    min_dist: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign each BED interval to a cluster of overlapping or nearby intervals with bioframe."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    distance = _integer(min_dist, "min_dist", 0, 1_000_000_000)
    bioframe = require("bioframe", "bio-bioframe")
    clustered = bioframe.cluster(_read_bed(bioframe, path, limit), min_dist=distance)
    records = _pandas_records(clustered, limit, "Clustered intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "cluster_bioframe_intervals",
        "min_dist": distance,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "intervals on one chromosome join a cluster when within min_dist",
        },
    }


def binnify_bioframe_chromosomes(
    chromsizes: Mapping[str, int],
    *,
    binsize: int,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Divide each chromosome of a size mapping into fixed-width bins with bioframe.binnify."""

    lengths = _lengths(chromsizes, "chromsizes")
    width = _integer(binsize, "binsize", 1, 1_000_000_000)
    limit = _rows_limit(max_rows)
    pandas = require("pandas", "bio-bioframe")
    bioframe = require("bioframe", "bio-bioframe")
    bins = bioframe.binnify(pandas.Series(lengths, dtype="int64"), width)
    records = _pandas_records(bins, limit, "Genome bins", ordered=True)
    return {
        "status": "completed",
        "operation": "binnify_bioframe_chromosomes",
        "binsize": width,
        "bins": records,
        "bin_count": len(records),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "bins are 0-based half-open; the last bin of a chromosome may be shorter",
        },
    }


def parse_bioframe_region(region: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a UCSC-style region string such as chr1:101-200 with bioframe.parse_region_string."""

    text = _text(region, "region", maximum=200)
    bioframe = require("bioframe", "bio-bioframe")
    try:
        chrom, start, end = bioframe.parse_region_string(text)
        start, end = int(start), int(end)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(
            f"region is not a UCSC-style chrom:start-end string: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "parse_bioframe_region",
        "region": text,
        "chrom": chrom,
        "start": start,
        "end": end,
        "provenance": {
            "package": _package("bioframe"),
            "convention": "numbers are returned as written, commas removed; no coordinate shift",
        },
    }


__all__ = [
    "binnify_bioframe_chromosomes",
    "cluster_bioframe_intervals",
    "expand_bioframe_intervals",
    "overlap_bioframe_intervals",
    "parse_bioframe_region",
    "subtract_bioframe_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
