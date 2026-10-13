"""bigWig header reading, bedGraph-to-bigWig writing, and binned statistics through pyBigWig.

All functions use pyBigWig.open on local files. Intervals are 0-based half-open. The write function
creates a new bigWig file and refuses to overwrite one. pyBigWig is imported only when a
function runs; a missing package raises BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _BIGWIG_STATISTICS,
    _data_lines,
    _integer,
    _lengths,
    _package,
    _plain,
    _text,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def _open(pybigwig: Any, path: Path, mode: str) -> Any:
    handle = pybigwig.open(str(path), mode)
    if handle is None:
        raise ConfigurationError(f"pyBigWig could not open {path.name} in mode {mode!r}.")
    return handle


def read_bigwig_header(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read the chromosome lengths and header summary of a bigWig file with pyBigWig."""

    source = input_path(path)
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open(pybigwig, source, "r")
    try:
        if not handle.isBigWig():
            raise ConfigurationError(f"{source.name} is not a bigWig file.")
        lengths = handle.chroms()
        header = handle.header()
    finally:
        handle.close()
    chromosomes = [{"chrom": name, "length": int(lengths[name])} for name in sorted(lengths)]
    return {
        "status": "completed",
        "operation": "read_bigwig_header",
        "chromosomes": chromosomes,
        "chromosome_count": len(chromosomes),
        "header": {str(key): _plain(value) for key, value in sorted(header.items())},
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "chromosome lengths are as stored in the bigWig header",
        },
    }


def write_bigwig_intervals(
    bedgraph: str | Path,
    output: str | Path,
    *,
    chromsizes: Mapping[str, int],
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a bigWig file from bedGraph rows (chrom, start, end, value) with pyBigWig."""

    source = input_path(bedgraph)
    lengths = _lengths(chromsizes, "chromsizes")
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {target.parent}")
    rows: list[tuple[str, int, int, float]] = []
    for number, fields in _data_lines(source):
        if len(fields) < 4:
            raise ConfigurationError(f"bedGraph line {number} needs chrom, start, end, value.")
        chrom = fields[0]
        try:
            start, end, value = int(fields[1]), int(fields[2]), float(fields[3])
        except ValueError as exc:
            raise ConfigurationError(f"bedGraph line {number} has a non-numeric field.") from exc
        if chrom not in lengths:
            raise ConfigurationError(f"chromosome {chrom} is missing from chromsizes.")
        if not 0 <= start < end <= lengths[chrom]:
            raise ConfigurationError(f"bedGraph line {number} lies outside {chrom}.")
        rows.append((chrom, start, end, value))
    if not rows:
        raise ConfigurationError("bedGraph contains no intervals.")
    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open(pybigwig, target, "w")
    try:
        handle.addHeader([(name, lengths[name]) for name in sorted(lengths)])
        handle.addEntries(
            [row[0] for row in rows],
            [row[1] for row in rows],
            ends=[row[2] for row in rows],
            values=[row[3] for row in rows],
        )
    except RuntimeError as exc:  # pyBigWig rejects overlapping or unsorted entries
        raise ConfigurationError(f"pyBigWig could not write the intervals: {exc}") from exc
    finally:
        handle.close()
    return {
        "status": "completed",
        "operation": "write_bigwig_intervals",
        "path": str(target),
        "interval_count": len(rows),
        "chromosomes": sorted(lengths),
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "intervals must not overlap; values are written as float",
        },
    }


def summarize_bigwig_bins(
    path: str | Path,
    *,
    chrom: str,
    start: int,
    end: int,
    bins: int,
    statistic: str = "mean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a bigWig region into equal bins with pyBigWig stats(nBins=...)."""

    source = input_path(path)
    name = _text(chrom, "chrom", maximum=256)
    begin = _integer(start, "start", 0, 10**12)
    finish = _integer(end, "end", 1, 10**12)
    count = _integer(bins, "bins", 1, 100_000)
    if begin >= finish:
        raise ConfigurationError("start must be less than end.")
    if statistic not in _BIGWIG_STATISTICS:
        raise ConfigurationError(f"statistic must be one of {sorted(_BIGWIG_STATISTICS)}.")
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open(pybigwig, source, "r")
    try:
        lengths = handle.chroms()
        if name not in lengths:
            raise ConfigurationError(f"chromosome {name} is not in {source.name}.")
        if finish > lengths[name]:
            raise ConfigurationError(f"end {finish} is beyond {name} length {lengths[name]}.")
        values = handle.stats(name, begin, finish, type=statistic, nBins=count)
    finally:
        handle.close()
    return {
        "status": "completed",
        "operation": "summarize_bigwig_bins",
        "chrom": name,
        "start": begin,
        "end": finish,
        "bins": count,
        "statistic": statistic,
        "values": [_plain(value) for value in values],
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "equal-width bins over [start, end); None marks an empty bin",
        },
    }


__all__ = [
    "read_bigwig_header",
    "summarize_bigwig_bins",
    "write_bigwig_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
