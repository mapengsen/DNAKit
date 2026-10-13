"""Overlap lengths, windows, score scaling, BEDPE cis/trans calls, GFF/BED conversion, and shuffled
Jaccard through pybedtools (review group G18g).

pybedtools drives the BEDTools binary (installed at /usr/bin/bedtools). Functions read local files
only. Coordinates are 0-based half-open unless a function returns GFF text, which is 1-based
inclusive. The shuffled Jaccard takes an explicit seed, and no function changes global state.
pybedtools is imported only when a function runs; a missing package raises
BackendUnavailableError with the install hint.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _bed_rows,
    _bedtools_call,
    _checked,
    _data_lines,
    _integer,
    _package,
    _plain,
    _rows_limit,
    _text,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MAX_COORDINATE = 2_000_000_000  # pybedtools.overlap takes C ints
_MAX_WIDTH = 1_000_000_000
_LABEL_PATTERN = re.compile(r"[A-Za-z0-9_.-]{1,64}")


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _check_bed(
    path: Path, limit: int, *, width: int | None = None, numeric: int | None = None
) -> None:
    """Validate BED3 to BED6 rows; optionally require a column count and a numeric column."""

    _bed_rows(path, limit)
    if width is None and numeric is None:
        return
    for number, fields in _data_lines(path):
        if width is not None and len(fields) != width:
            raise ConfigurationError(
                f"{path.name} line {number} must have {width} columns; found {len(fields)}."
            )
        if numeric is not None:
            if len(fields) <= numeric:
                raise ConfigurationError(f"{path.name} line {number} needs column {numeric + 1}.")
            try:
                value = float(fields[numeric])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number} column {numeric + 1} must be numeric."
                ) from exc
            if not math.isfinite(value):
                raise ConfigurationError(
                    f"{path.name} line {number} column {numeric + 1} must be finite."
                )


def _check_gff(path: Path, limit: int) -> None:
    count = 0
    for number, fields in _data_lines(path):
        count += 1
        if count > limit:
            raise ConfigurationError(f"{path.name} has more than max_rows={limit} data lines.")
        if len(fields) != 9:
            raise ConfigurationError(f"GFF line {number} must have 9 columns; found {len(fields)}.")
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError as exc:
            raise ConfigurationError(f"GFF line {number} needs integer start and end.") from exc
        if not 1 <= start <= end:
            raise ConfigurationError(f"GFF line {number} must satisfy 1 <= start <= end.")
    if not count:
        raise ConfigurationError(f"{path.name} contains no GFF features.")


def _check_columns(path: Path, expected: int, what: str) -> None:
    seen = 0
    for number, fields in _data_lines(path):
        if len(fields) != expected:
            raise ConfigurationError(
                f"{what} line {number} must have {expected} columns; found {len(fields)}."
            )
        seen += 1
    if not seen:
        raise ConfigurationError(f"{path.name} contains no data lines.")


def _check_named_bed(path: Path, limit: int, what: str) -> None:
    """Validate BED3 to BED6 rows and require a name column (BED4 or wider)."""

    _check_bed(path, limit)
    for number, fields in _data_lines(path):
        if len(fields) < 4:
            raise ConfigurationError(f"{what} line {number} needs a name in column 4.")


def _check_genome(path: Path) -> None:
    seen = 0
    for number, fields in _data_lines(path):
        if len(fields) != 2:
            raise ConfigurationError(f"genome line {number} must be chrom and length.")
        try:
            length = int(fields[1])
        except ValueError as exc:
            raise ConfigurationError(f"genome line {number} length must be an integer.") from exc
        if length < 1:
            raise ConfigurationError(f"genome line {number} length must be positive.")
        seen += 1
    if not seen:
        raise ConfigurationError(f"{path.name} contains no chromosome sizes.")


def _bed_records(bedtool: Any, limit: int, what: str) -> list[dict[str, Any]]:
    records = []
    for line in bedtool:
        fields = str(line).rstrip("\n").split("\t")
        records.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "name": fields[3] if len(fields) > 3 else None,
                "score": fields[4] if len(fields) > 4 else None,
                "strand": fields[5] if len(fields) > 5 else None,
            }
        )
    _checked(records, limit, what)
    return sorted(records, key=lambda r: (r["chrom"], r["start"], r["end"], r["name"] or ""))


def _gff_records(bedtool: Any, limit: int) -> list[dict[str, Any]]:
    records = []
    for line in bedtool:
        fields = str(line).rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ConfigurationError("bed2gff returned a line without nine GFF columns.")
        records.append(
            {
                "seqid": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "score": fields[5],
                "strand": fields[6],
                "frame": fields[7],
                "attributes": fields[8],
            }
        )
    _checked(records, limit, "GFF records")
    return sorted(records, key=lambda r: (r["seqid"], r["start"], r["end"], r["attributes"]))


def _window(
    source: str | Path,
    *,
    function: str,
    operation: str,
    max_rows: int,
    upstream: int | None = None,
    downstream: int | None = None,
    width: int | None = None,
    convention: str,
) -> dict[str, Any]:
    path = input_path(source)
    limit = _rows_limit(max_rows)
    _check_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    options: dict[str, int] = {}
    if upstream is not None and downstream is not None:
        options = {"upstream": upstream, "downstream": downstream}
    elif width is not None:
        options = {"width": width}
    windows = _bedtools_call(
        function,
        lambda: pybedtools.BedTool(str(path)).each(getattr(featurefuncs, function), **options),
    )
    records = _bed_records(windows, limit, f"{function} windows")
    result: dict[str, Any] = {
        "status": "completed",
        "operation": operation,
        "records": records,
        "record_count": len(records),
        "provenance": {"package": _package("pybedtools"), "convention": convention},
    }
    result.update(options)
    return result


def overlap_pybedtools_bp(
    start1: int,
    end1: int,
    start2: int,
    end2: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the signed overlap in bp of two half-open intervals with pybedtools.overlap."""

    first = _integer(start1, "start1", 0, _MAX_COORDINATE)
    first_end = _integer(end1, "end1", 0, _MAX_COORDINATE)
    second = _integer(start2, "start2", 0, _MAX_COORDINATE)
    second_end = _integer(end2, "end2", 0, _MAX_COORDINATE)
    if first >= first_end or second >= second_end:
        raise ConfigurationError("each interval must satisfy start < end.")
    pybedtools = require("pybedtools", "bio-pybedtools")
    # pybedtools.overlap takes the two starts before the two ends.
    signed = int(pybedtools.overlap(first, second, first_end, second_end))
    return {
        "status": "completed",
        "operation": "overlap_pybedtools_bp",
        "overlap_bp": signed,
        "overlaps": signed > 0,
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "half-open 0-based intervals; touching intervals give 0; "
            "a negative overlap_bp is the gap between them",
        },
    }


def normalize_pybedtools_scores_to_length(
    source: str | Path,
    *,
    index: int = 4,
    scalar: float = 0.001,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Divide the value at a BED column by the interval length and multiply by scalar."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    column = _integer(index, "index", 0, 5)
    factor = _finite(scalar, "scalar")
    _check_bed(path, limit, numeric=column)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    scaled = _bedtools_call(
        "normalized_to_length",
        lambda: pybedtools.BedTool(str(path)).each(
            featurefuncs.normalized_to_length, idx=column, scalar=factor
        ),
    )
    records = _bed_records(scaled, limit, "normalized intervals")
    return {
        "status": "completed",
        "operation": "normalize_pybedtools_scores_to_length",
        "index": column,
        "scalar": factor,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "value * scalar / (end - start); the default 0.001 gives a per-kb rate",
        },
    }


def scale_pybedtools_bedgraph_values(
    source: str | Path,
    *,
    scalar: float,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multiply the value column of a four-column bedGraph by scalar with pybedtools."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    factor = _finite(scalar, "scalar")
    _check_bed(path, limit, width=4, numeric=3)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    scaled = _bedtools_call(
        "bedgraph_scale",
        lambda: pybedtools.BedTool(str(path)).each(featurefuncs.bedgraph_scale, scalar=factor),
    )
    records = _bed_records(scaled, limit, "scaled bedGraph intervals")
    return {
        "status": "completed",
        "operation": "scale_pybedtools_bedgraph_values",
        "scalar": factor,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "bedGraph value is written to the fourth column as a float string",
        },
    }


def center_pybedtools_coordinates(
    source: str | Path,
    *,
    width: int = 100,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the window of width bp centred on each interval (pybedtools featurefuncs.center)."""

    window = _integer(width, "width", 1, _MAX_WIDTH)
    return _window(
        source,
        function="center",
        operation="center_pybedtools_coordinates",
        max_rows=max_rows,
        width=window,
        convention=(
            "window is 2 * (width // 2) bp around the midpoint, and width 1 gives the base at "
            "the midpoint; intervals shorter than width are unchanged; a window that would start "
            "before 0 is clamped to start 1 by upstream"
        ),
    )


def three_prime_pybedtools_coordinates(
    source: str | Path,
    *,
    upstream: int = 500,
    downstream: int = 500,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the window around each feature's 3' end using pybedtools featurefuncs.three_prime."""

    before = _integer(upstream, "upstream", 0, 1_000_000_000)
    after = _integer(downstream, "downstream", 0, 1_000_000_000)
    return _window(
        source,
        function="three_prime",
        operation="three_prime_pybedtools_coordinates",
        max_rows=max_rows,
        upstream=before,
        downstream=after,
        convention=(
            "3' end is stop on + or start on -; window is upstream bp before and downstream bp "
            "after it, in the feature's strand direction; unstranded features act as +"
        ),
    )


def filter_pybedtools_intervals_by_length(
    source: str | Path,
    *,
    greater_than: int | None = None,
    less_than: int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep BED intervals whose length is above and/or below strict thresholds (pybedtools)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    if greater_than is None and less_than is None:
        raise ConfigurationError("set greater_than, less_than, or both.")
    rules: dict[str, int | None] = {"greater_than": None, "less_than": None}
    if greater_than is not None:
        rules["greater_than"] = _integer(greater_than, "greater_than", 0, _MAX_COORDINATE)
    if less_than is not None:
        rules["less_than"] = _integer(less_than, "less_than", 0, _MAX_COORDINATE)
    _check_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")

    def run() -> Any:
        bedtool = pybedtools.BedTool(str(path))
        if rules["greater_than"] is not None:
            bedtool = bedtool.filter(featurefuncs.greater_than, size=rules["greater_than"])
        if rules["less_than"] is not None:
            bedtool = bedtool.filter(featurefuncs.less_than, size=rules["less_than"])
        return bedtool

    kept = _bed_records(_bedtools_call("filter", run), limit, "filtered intervals")
    return {
        "status": "completed",
        "operation": "filter_pybedtools_intervals_by_length",
        "greater_than": rules["greater_than"],
        "less_than": rules["less_than"],
        "records": kept,
        "record_count": len(kept),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "length is end - start; both thresholds are strict (> and <)",
        },
    }


def gff_to_bed_pybedtools(
    source: str | Path,
    *,
    name_field: str | int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GFF features to BED6 with pybedtools featurefuncs.gff2bed."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _check_gff(path, limit)
    options: dict[str, Any] = {}
    if name_field is not None:
        if isinstance(name_field, int):
            options["name_field"] = _integer(name_field, "name_field", 0, 8)
        else:
            options["name_field"] = _text(name_field, "name_field", maximum=64)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    converted = _bedtools_call(
        "gff2bed",
        lambda: pybedtools.BedTool(str(path)).each(featurefuncs.gff2bed, **options),
    )
    records = _bed_records(converted, limit, "BED features")
    return {
        "status": "completed",
        "operation": "gff_to_bed_pybedtools",
        "name_field": options.get("name_field"),
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "BED6 start is 0-based; the name is the GFF ID unless name_field picks "
            "an attribute or column, and a missing name is '.'",
        },
    }


def bed_to_gff_pybedtools(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert BED3 to BED6 intervals to GFF records with pybedtools featurefuncs.bed2gff."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _check_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    featurefuncs = require("pybedtools.featurefuncs", "bio-pybedtools")
    converted = _bedtools_call(
        "bed2gff",
        lambda: pybedtools.BedTool(str(path)).each(featurefuncs.bed2gff),
    )
    records = _gff_records(converted, limit)
    return {
        "status": "completed",
        "operation": "bed_to_gff_pybedtools",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "GFF start is the BED start plus 1 (1-based); source and feature are "
            "'.'; name goes in the Name attribute",
        },
    }


def cis_trans_pybedtools_interactions(
    bedpe: str | Path,
    queries: Mapping[str, str | Path],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag BEDPE ends with named BED4 query sets and list cis and distal partners (pybedtools)."""

    path = input_path(bedpe)
    limit = _rows_limit(max_rows)
    if not isinstance(queries, Mapping) or not 2 <= len(queries) <= 8:
        # bedtools intersect drops the -names column when only one query file is given.
        raise ConfigurationError("queries must map two to eight label names to BED4 files.")
    _check_columns(path, 10, "BEDPE")
    query_paths: dict[str, Path] = {}
    for label, source in sorted(queries.items()):
        name = _text(label, "query label", maximum=64)
        if not _LABEL_PATTERN.fullmatch(name):
            raise ConfigurationError(f"query label {name!r} must use letters, digits, _ . or -.")
        query_paths[name] = input_path(source)
        _check_named_bed(query_paths[name], limit, f"query {name}")
    pybedtools = require("pybedtools", "bio-pybedtools")
    long_range = require("pybedtools.contrib.long_range_interaction", "bio-pybedtools")

    def run() -> list[dict[str, Any]]:
        # tag_bedpe yields groupby iterators, so the table is built before the next pair is read.
        iterator, count, extra = long_range.tag_bedpe(
            pybedtools.BedTool(str(path)),
            {name: pybedtools.BedTool(str(p)) for name, p in query_paths.items()},
            verbose=False,
        )
        frame = long_range.cis_trans_interactions(iterator, count, extra, verbose=False)
        return frame.to_dict("records")

    rows = _bedtools_call("cis_trans_interactions", run)
    records = [{str(key): _plain(value) for key, value in row.items()} for row in rows]
    _checked(records, limit, "BEDPE cis/trans rows")
    records.sort(key=lambda row: tuple(str(row[key]) for key in sorted(row)))
    return {
        "status": "completed",
        "operation": "cis_trans_pybedtools_interactions",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "cis_* names a query on the same BEDPE end as the target; distal_* names "
            "a query on the other end; '.' marks a missing partner",
        },
    }


def _jaccard(values: dict[str, Any]) -> dict[str, Any]:
    return {
        "intersection_bp": int(values["intersection"]),
        "union_bp": int(values["union-intersection"]),
        "jaccard": float(values["jaccard"]),
        "n_intersections": int(values["n_intersections"]),
    }


def random_jaccard_pybedtools_intervals(
    a: str | Path,
    b: str | Path,
    genome: str | Path,
    *,
    seed: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the Jaccard index of two BED files with one seeded random shuffle of the first."""

    path_a = input_path(a)
    path_b = input_path(b)
    genome_path = input_path(genome)
    limit = _rows_limit(max_rows)
    seed_value = _integer(seed, "seed", 0, 2**31 - 1)
    _check_bed(path_a, limit)
    _check_bed(path_b, limit)
    _check_genome(genome_path)
    pybedtools = require("pybedtools", "bio-pybedtools")
    stats = require("pybedtools.stats", "bio-pybedtools")

    def observed() -> dict[str, Any]:
        return (
            pybedtools.BedTool(str(path_a)).sort().jaccard(pybedtools.BedTool(str(path_b)).sort())
        )

    def shuffled() -> dict[str, Any]:
        return stats.random_jaccard(
            pybedtools.BedTool(str(path_a)),
            pybedtools.BedTool(str(path_b)),
            str(genome_path),
            {"seed": seed_value},
            {},
        )

    return {
        "status": "completed",
        "operation": "random_jaccard_pybedtools_intervals",
        "seed": seed_value,
        "observed": _jaccard(_bedtools_call("jaccard", observed)),
        "shuffled": _jaccard(_bedtools_call("random_jaccard", shuffled)),
        "provenance": {
            "package": _package("pybedtools"),
            "convention": "A is shuffled once with bedtools shuffle -seed on the genome file; "
            "one draw from the null, not a p-value",
        },
    }


__all__ = [
    "bed_to_gff_pybedtools",
    "center_pybedtools_coordinates",
    "cis_trans_pybedtools_interactions",
    "filter_pybedtools_intervals_by_length",
    "gff_to_bed_pybedtools",
    "normalize_pybedtools_scores_to_length",
    "overlap_pybedtools_bp",
    "random_jaccard_pybedtools_intervals",
    "scale_pybedtools_bedgraph_values",
    "three_prime_pybedtools_coordinates",
]

_register_public_verbose(globals(), tuple(__all__))
