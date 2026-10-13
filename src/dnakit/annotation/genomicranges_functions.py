"""Genomicranges wrappers for GTF gene ranges, promoter and terminator windows, and interval tiling.

Each public function runs one upstream genomicranges operation on local input and returns a
JSON-serializable dictionary with ``status``, ``operation``, result fields, and ``provenance``.

Coordinates in the returned records are 0-based half-open (BED). Upstream genomicranges stores
ranges 1-based and inclusive, and its GTF reader keeps the GTF end minus one; these wrappers
restore the GTF end so that a feature covers the same bases as in the GTF file.

Functions are importable as ``dnakit.annotation.genomicranges_functions.<name>``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_COMMENT_PREFIXES = ("#", "track", "browser")
_SITES = ("TSS", "TES")


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _choice(value: object, name: str, allowed: Iterable[str]) -> str:
    options = tuple(allowed)
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {', '.join(options)}.")
    return value


def _plain(value: Any) -> Any:
    """Convert numpy, biocframe, and missing values to JSON-compatible Python objects."""

    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return _plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _record_order(record: Mapping[str, Any]) -> tuple[tuple[int, Any], ...]:
    """Type-safe sort key so rows come back in one order regardless of hash seed."""

    keys: list[tuple[int, Any]] = []
    for value in record.values():
        if value is None:
            keys.append((3, ""))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            keys.append((0, float(value)))
        else:
            keys.append((2, str(value)))
    return tuple(keys)


def _ordered(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(records, key=_record_order)


def _checked(rows: list[Any], limit: int, what: str) -> list[Any]:
    if len(rows) > limit:
        raise ConfigurationError(f"{what} has {len(rows)} rows; max_rows={limit} is the limit.")
    return rows


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("genomicranges"), "convention": convention}


_BED_CONVENTION = "0-based half-open coordinates (BED)"


def _bed_rows(path: Path, limit: int) -> list[tuple[str, int, int, str | None]]:
    """Read chrom, start, end, and optional strand from a BED file."""

    rows: list[tuple[str, int, int, str | None]] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            if not text or text.startswith(_COMMENT_PREFIXES):
                continue
            fields = text.split("\t")
            if len(fields) < 3:
                raise ConfigurationError(f"{path}:{number} needs at least three BED columns.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(f"{path}:{number} has non-integer coordinates.") from exc
            if start < 0 or end <= start:
                raise ConfigurationError(f"{path}:{number} must satisfy 0 <= start < end.")
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else None
            rows.append((fields[0], start, end, strand))
    if not rows:
        raise ConfigurationError(f"{path} contains no intervals.")
    _checked(rows, limit, "Interval file")
    return rows


# GTF gene ranges.


def read_gtf_ranges(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GTF file with genomicranges and return one record per feature with its attributes.

    Each record has ``chrom``, 0-based ``start``, exclusive ``end``, ``strand``, and every GTF
    attribute column (for example ``feature``, ``gene_id``, and ``gene_name``) as a string.
    """

    source = input_path(path)
    limit = _rows_limit(max_rows)
    genomicranges = require("genomicranges", "bio-genomicranges")
    ranges = genomicranges.read_gtf(str(source))
    if len(ranges) > limit:
        raise ConfigurationError(f"GTF file has {len(ranges)} rows; max_rows={limit} is the limit.")
    seqnames = list(ranges.get_seqnames())
    starts = list(ranges.get_start())
    ends = list(ranges.get_end())
    strands = list(ranges.get_strand())
    mcols = ranges.get_mcols()
    columns = {name: list(mcols.get_column(name)) for name in mcols.get_column_names()}
    features = []
    for index in range(len(ranges)):
        strand_code = int(strands[index])
        record: dict[str, Any] = {
            "chrom": str(seqnames[index]),
            # Upstream starts are 1-based GTF starts; upstream ends are GTF ends minus one.
            "start": int(starts[index]) - 1,
            "end": int(ends[index]) + 1,
            "strand": "+" if strand_code == 1 else "-" if strand_code == -1 else None,
        }
        for name, values in columns.items():
            record[str(name)] = _plain(values[index])
        if record["start"] < 0:
            raise ConfigurationError("A GTF feature starts before position 1.")
        features.append(record)
    return {
        "status": "completed",
        "operation": "read_gtf_ranges",
        "feature_count": len(features),
        "features": _ordered(features),
        "provenance": _provenance(
            _BED_CONVENTION + "; upstream GTF end minus one is corrected to the GTF end"
        ),
    }


# Promoter and terminator windows.


def upstream_downstream_regions(
    source: str | Path,
    *,
    upstream: int = 1000,
    downstream: int = 0,
    site: str = "TSS",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return promoter (``TSS``) or terminator (``TES``) windows around each BED interval.

    Strand ``+`` and missing strands use the start side for ``TSS``; strand ``-`` uses the end
    side. ``upstream`` bases lie before the site and ``downstream`` bases after it, both in
    the transcription direction. Windows that would start before position 0 are clipped.
    """

    path = input_path(source)
    before = _integer(upstream, "upstream", 0, 10**9)
    after = _integer(downstream, "downstream", 0, 10**9)
    mode = _choice(site, "site", _SITES)
    limit = _rows_limit(max_rows)
    intervals = _bed_rows(path, limit)
    numpy = require("numpy", "bio-genomicranges")
    utils = require("genomicranges.utils", "bio-genomicranges")
    # Upstream works on 1-based inclusive ranges with integer strand codes (+1, -1).
    starts = numpy.array([start + 1 for _, start, _, _ in intervals], dtype=numpy.int64)
    ends = numpy.array([end for _, _, end, _ in intervals], dtype=numpy.int64)
    strands = numpy.array(
        [-1 if strand == "-" else 1 for _, _, _, strand in intervals], dtype=numpy.int64
    )
    new_starts, widths = utils.compute_up_down(starts, ends, strands, before, after, mode)
    windows = []
    for (chrom, start, end, strand), new_start, width in zip(
        intervals, new_starts.tolist(), widths.tolist(), strict=True
    ):
        window_start = max(0, int(new_start) - 1)
        window_end = int(new_start) + int(width) - 1
        windows.append(
            {
                "chrom": chrom,
                "start": window_start,
                "end": window_end,
                "strand": strand,
                "source_start": start,
                "source_end": end,
            }
        )
    return {
        "status": "completed",
        "operation": "upstream_downstream_regions",
        "site": mode,
        "upstream": before,
        "downstream": after,
        "window_count": len(windows),
        "windows": _ordered(windows),
        "provenance": _provenance(
            _BED_CONVENTION + "; upstream returns widths, converted to window ends here"
        ),
    }


# Interval tiling.


def split_interval_bins(
    start: int,
    end: int,
    bin_size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tile the half-open range ``[start, end)`` into consecutive bins of ``bin_size`` bases.

    The last bin may be shorter. The bins cover the range exactly once.
    """

    first = _integer(start, "start", 0, 10**12)
    last = _integer(end, "end", first + 1, 10**12)
    size = _integer(bin_size, "bin_size", 1, 10**12)
    utils = require("genomicranges.utils", "bio-genomicranges")
    raw = utils.split_intervals(first, last - 1, size)
    bins = []
    for index, (bin_start, offset) in enumerate(raw):
        bin_end = int(bin_start) + int(offset) + 1
        bins.append({"index": index, "start": int(bin_start), "end": bin_end})
    if sum(item["end"] - item["start"] for item in bins) != last - first:
        raise ConfigurationError("genomicranges returned bins that do not cover the range.")
    return {
        "status": "completed",
        "operation": "split_interval_bins",
        "bin_count": len(bins),
        "bins": bins,
        "provenance": _provenance(_BED_CONVENTION + "; the last bin may be shorter"),
    }


def slide_interval_windows(
    start: int,
    end: int,
    width: int,
    step: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return sliding windows of ``width`` bases moved by ``step`` across ``[start, end)``.

    Windows lie fully inside the range. If the range is shorter than ``width`` it is returned
    as one window.
    """

    first = _integer(start, "start", 0, 10**12)
    last = _integer(end, "end", first + 1, 10**12)
    window = _integer(width, "width", 1, 10**12)
    move = _integer(step, "step", 1, 10**12)
    utils = require("genomicranges.utils", "bio-genomicranges")
    raw = utils.slide_intervals(first, last - 1, window, move)
    windows = []
    for index, (window_start, offset) in enumerate(raw):
        window_end = int(window_start) + int(offset) + 1
        windows.append({"index": index, "start": int(window_start), "end": window_end})
    return {
        "status": "completed",
        "operation": "slide_interval_windows",
        "window_count": len(windows),
        "windows": windows,
        "provenance": _provenance(_BED_CONVENTION + f"; width={window}; step={move}"),
    }


__all__ = [
    "read_gtf_ranges",
    "slide_interval_windows",
    "split_interval_bins",
    "upstream_downstream_regions",
]

_register_public_verbose(globals(), tuple(__all__))
