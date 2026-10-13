"""Fixed-column GTF parsing, attribute expansion, and GTF writing through gtfparse.

Functions read local GTF files and write one new GTF file. Coordinates are kept as stored (1-based,
inclusive). gtfparse is imported only when a function runs; a missing package raises
BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _checked,
    _gtf_limit,
    _gtf_rows,
    _package,
    _pandas_records,
    _rows_limit,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def parse_gtfparse_fixed_columns(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the nine GTF columns with the raw attribute text using gtfparse.parse_gtf."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    try:
        frame = gtfparse.parse_gtf(str(path), split_attributes=False)
    except Exception as exc:  # gtfparse raises for malformed GTF lines
        raise ConfigurationError(f"gtfparse could not read {path.name}: {exc}") from exc
    records = _pandas_records(frame, limit, "GTF rows")
    return {
        "status": "completed",
        "operation": "parse_gtfparse_fixed_columns",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "start and end are 1-based inclusive; attribute is the raw column text",
        },
    }


def expand_gtfparse_attributes(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GTF file and expand its attribute column into one column per key with gtfparse."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    try:
        frame = gtfparse.parse_gtf_and_expand_attributes(str(path))
    except Exception as exc:  # gtfparse raises ParsingError and pandas errors
        raise ConfigurationError(f"gtfparse could not read {path.name}: {exc}") from exc
    records = _gtf_rows(frame, limit, "GTF rows")
    columns = sorted(str(column) for column in frame.columns)
    return {
        "status": "completed",
        "operation": "expand_gtfparse_attributes",
        "columns": columns,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "absent attribute keys are None; coordinates are 1-based inclusive",
        },
    }


def write_gtfparse_table(
    source: str | Path,
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a GTF file and write it back out as a new GTF file with gtfparse.write_gtf."""

    path = input_path(source)
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {target.parent}")
    gtfparse = require("gtfparse", "bio-gtfparse")
    try:
        frame = gtfparse.parse_gtf(str(path), split_attributes=False)
        gtfparse.write_gtf(frame, str(target))
    except Exception as exc:  # gtfparse raises for malformed GTF lines or unwritable paths
        raise ConfigurationError(f"gtfparse could not write {target.name}: {exc}") from exc
    lines = target.read_text(encoding="utf-8").splitlines()
    _checked(lines, 10_000_000, "written GTF")
    return {
        "status": "completed",
        "operation": "write_gtfparse_table",
        "path": str(target),
        "line_count": len(lines),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "frame is written as parsed; '.' frames come back as 0",
        },
    }


__all__ = [
    "expand_gtfparse_attributes",
    "parse_gtfparse_fixed_columns",
    "write_gtfparse_table",
]

_register_public_verbose(globals(), tuple(__all__))
