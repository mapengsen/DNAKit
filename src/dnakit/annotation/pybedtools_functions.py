"""Feature coordinates, BEDPE end tagging, and chromosome-size files through pybedtools.

pybedtools drives the BEDTools binary (installed at /usr/bin/bedtools). Functions read local files
only. Coordinates are 0-based half-open. pybedtools is imported only when a function runs; a missing
package raises BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _bedtools_call,
    _checked,
    _data_lines,
    _integer,
    _lengths,
    _package,
    _rows_limit,
    _scan_bed,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def _columns_must_be(path: Path, expected: int, what: str) -> None:
    seen = 0
    for number, fields in _data_lines(path):
        if len(fields) != expected:
            raise ConfigurationError(
                f"{what} line {number} must have {expected} columns; found {len(fields)}."
            )
        seen += 1
    if not seen:
        raise ConfigurationError(f"{path.name} contains no data lines.")


def _feature_records(bedtool: Any, limit: int, what: str) -> list[dict[str, Any]]:
    records = []
    for line in bedtool:
        fields = str(line).rstrip("\n").split("\t")
        records.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "name": fields[3] if len(fields) > 3 else None,
                "strand": fields[5] if len(fields) > 5 else None,
            }
        )
    _checked(records, limit, what)
    return sorted(records, key=lambda r: (r["chrom"], r["start"], r["end"], r["name"] or ""))


def five_prime_pybedtools_coordinates(
    source: str | Path,
    *,
    upstream: int = 500,
    downstream: int = 500,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the window around each feature's 5' end using pybedtools featurefuncs.five_prime."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    before = _integer(upstream, "upstream", 0, 1_000_000_000)
    after = _integer(downstream, "downstream", 0, 1_000_000_000)
    _scan_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    windows = _bedtools_call(
        "five_prime",
        lambda: pybedtools.BedTool(str(path)).each(
            featurefuncs.five_prime, upstream=before, downstream=after
        ),
    )
    records = _feature_records(windows, limit, "5' windows")
    return {
        "status": "completed",
        "operation": "five_prime_pybedtools_coordinates",
        "upstream": before,
        "downstream": after,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "5' end is start on + or end on -; window is clipped by upstream code",
        },
    }


def midpoint_pybedtools_coordinates(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the 1-bp midpoint of each BED feature with pybedtools featurefuncs.midpoint."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    midpoints = _bedtools_call(
        "midpoint", lambda: pybedtools.BedTool(str(path)).each(featurefuncs.midpoint)
    )
    records = _feature_records(midpoints, limit, "midpoints")
    return {
        "status": "completed",
        "operation": "midpoint_pybedtools_coordinates",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "each midpoint is a 1-bp interval at floor((start + end) / 2)",
        },
    }


def _tags(rows: Any) -> list[dict[str, Any]]:
    tags = []
    for row in rows:
        fields = list(row)
        overlap = int(fields[-1])
        if overlap > 0:
            tags.append({"query": fields[9], "overlap_bp": overlap})
    return sorted(tags, key=lambda tag: (tag["query"], tag["overlap_bp"]))


def tag_pybedtools_bedpe(
    bedpe: str | Path,
    queries: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag each BEDPE end with the named query intervals it overlaps (pybedtools contrib)."""

    path = input_path(bedpe)
    query_path = input_path(queries)
    limit = _rows_limit(max_rows)
    _columns_must_be(path, 10, "BEDPE")
    _columns_must_be(query_path, 4, "query BED4")
    pybedtools = require("pybedtools", "bio-pybedtools")
    long_range = require("pybedtools.contrib.long_range_interaction", "bio-pybedtools")

    def run() -> list[dict[str, Any]]:
        # Each end is a groupby iterator that is valid only until the next outer item.
        tagged, _, _ = long_range.tag_bedpe(
            pybedtools.BedTool(str(path)),
            {"query": pybedtools.BedTool(str(query_path))},
            verbose=False,
        )
        return [
            {"name": name, "left": _tags(left), "right": _tags(right)}
            for name, left, right in tagged
        ]

    pairs = _bedtools_call("tag_bedpe", run)
    _checked(pairs, limit, "BEDPE pairs")
    pairs.sort(key=lambda pair: str(pair["name"]))
    return {
        "status": "completed",
        "operation": "tag_pybedtools_bedpe",
        "pairs": pairs,
        "pair_count": len(pairs),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "query names come from column 4 of the query BED4; overlap is in bp",
        },
    }


def write_pybedtools_chromsizes(
    chromsizes: dict[str, int],
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a chromosome-size file with pybedtools; the output path must not exist yet."""

    lengths = _lengths(chromsizes, "chromsizes")
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {target.parent}")
    pybedtools = require("pybedtools", "bio-pybedtools")
    _bedtools_call(
        "chromsizes_to_file",
        lambda: pybedtools.chromsizes_to_file(
            {chrom: (0, length) for chrom, length in sorted(lengths.items())}, fn=str(target)
        ),
    )
    chromosomes = sorted(
        (
            {"chrom": fields[0], "length": int(fields[1])}
            for fields in (
                line.split("\t") for line in target.read_text(encoding="utf-8").splitlines()
            )
            if len(fields) == 2
        ),
        key=lambda row: row["chrom"],
    )
    return {
        "status": "completed",
        "operation": "write_pybedtools_chromsizes",
        "path": str(target),
        "chromosomes": chromosomes,
        "chromosome_count": len(chromosomes),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "two-column tab file: chromosome name and length",
        },
    }


__all__ = [
    "five_prime_pybedtools_coordinates",
    "midpoint_pybedtools_coordinates",
    "tag_pybedtools_bedpe",
    "write_pybedtools_chromsizes",
]

_register_public_verbose(globals(), tuple(__all__))
