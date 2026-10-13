"""Review wrappers for bioframe interval, region, sequence, file, and assembly utilities.

Each public function reads local BED, FASTA, SAM or BAM, tabix, bigBed, chromsize, or
in-memory interval inputs and calls one bioframe operation. Coordinates are 0-based half-open
unless a docstring says otherwise. Region strings are read as written (no coordinate shift),
which matches the existing ``parse_bioframe_region`` wrapper. bioframe is imported only when a
function runs; a missing package raises BackendUnavailableError with the install hint.
Network-backed bioframe paths (UCSC chromsize and centromere downloads, UCSC mRNA fetches) are
not exposed here; every input must be a local file or an in-memory mapping.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _data_lines,
    _flag,
    _integer,
    _lengths,
    _package,
    _pandas_records,
    _plain,
    _rows_limit,
    _scan_bed,
    _text,
)
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_BED_SCHEMAS = ("bed3", "bed4", "bed5", "bed6", "bed9", "bed12")
_BED_WIDTHS = frozenset({3, 4, 5, 6, 9, 12})
_SANITIZE_ACTIONS = (None, "drop", "flip")
_RUN_COLUMNS = ("name", "score", "strand")
_RELATIVE_TO = ("midpoints", "endpoints")
_MRNA_COLUMNS = ("tName", "tStart", "tEnd")
_COORD_MAX = 10**12
_REGION_MAX = 200
_SEQUENCE_MAX = 10_000_000
_DIGEST_CUT_OFFSET = 2


# Helpers. Upstream failures on bad coordinates, names, or file contents become ConfigurationError.


def _bioframe() -> Any:
    return require("bioframe", "bio-bioframe")


def _run(label: str, call: Any) -> Any:
    try:
        return call()
    except (IndexError, KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{label} rejected the input: {exc}") from exc


def _choice(value: object, name: str, allowed: tuple[Any, ...]) -> Any:
    if value not in allowed:
        options = ", ".join(repr(item) for item in allowed)
        raise ConfigurationError(f"{name} must be one of {options}.")
    return value


def _selector(value: object, name: str) -> list[str] | str | None:
    if value is None or value == "all":
        return value
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be None, 'all', or a non-empty list of strings.")
    return [_text(item, name, maximum=64) for item in value]


def _optional_int(value: object, name: str) -> int | None:
    return None if value is None else _integer(value, name, 0, _COORD_MAX)


def _region_bounds(
    chrom: object, start: object, end: object
) -> tuple[str | None, int | None, int | None]:
    """Validate an optional chrom, start, end triple used by indexed readers."""

    name = None if chrom is None else _text(chrom, "chrom", maximum=256)
    begin = _optional_int(start, "start")
    finish = _optional_int(end, "end")
    if name is None and (begin is not None or finish is not None):
        raise ConfigurationError("start and end need a chrom.")
    if begin is not None and finish is not None and begin >= finish:
        raise ConfigurationError("start must be less than end.")
    return name, begin, finish


def _bed(bioframe: Any, path: Path, limit: int) -> Any:
    """Read a BED3 to BED6 file; the column count selects the bioframe schema."""

    schema = f"bed{_scan_bed(path, limit)}"
    return _run("read_table", lambda: bioframe.read_table(str(path), schema=schema, comment="#"))


def _view(bioframe: Any, path: Path, limit: int) -> Any:
    frame = _bed(bioframe, path, limit)
    if "name" not in frame.columns:
        raise ConfigurationError(f"{path.name} is a view and needs names: use BED4 to BED6.")
    return frame


def _records(
    frame: Any,
    limit: int,
    what: str,
    *,
    ordered: bool = False,
    fields: tuple[str, ...] = ("start", "end"),
) -> list[dict[str, Any]]:
    """Records as JSON-compatible dicts; coordinates that pandas upcast to float become int."""

    rows = _pandas_records(frame, limit, what, ordered=ordered)
    for row in rows:
        for field in fields:
            value = row.get(field)
            if isinstance(value, float) and math.isfinite(value) and value.is_integer():
                row[field] = int(value)
    return rows


def _fasta_records(bioframe: Any, source: str | Path) -> Any:
    path = input_path(source)
    # pysam writes <fasta>.fai beside the file on first use when the index is missing.
    return _run("load_fasta", lambda: bioframe.load_fasta(str(path)))


def _check(operation: str, check: str, passed: bool, count: int, convention: str) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "check": check,
        "passed": passed,
        "record_count": count,
        "provenance": {"package": _package("bioframe"), "convention": convention},
    }


# Validation checks on BED files, views, and intervals.


def validate_bioframe_bed(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a BED3 to BED6 file against the bioframe bedframe rules with is_bedframe."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    passed = _run("is_bedframe", lambda: bool(bioframe.is_bedframe(frame, raise_errors=False)))
    return _check(
        "validate_bioframe_bed",
        "bedframe",
        passed,
        len(frame),
        "passed is false when a column type or start/end rule fails; no row is changed",
    )


def check_bioframe_cataloged(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that every interval name in a BED4+ file names a region of a view (is_cataloged)."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    passed = _run(
        "is_cataloged",
        lambda: bool(bioframe.is_cataloged(frame, views, raise_errors=False, df_view_col="name")),
    )
    return _check(
        "check_bioframe_cataloged",
        "cataloged",
        passed,
        len(frame),
        "the name column of source must match region names of view; coordinates are not used",
    )


def check_bioframe_contained(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that every interval of a BED file lies inside some region of a view (is_contained)."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    passed = _run(
        "is_contained",
        lambda: bool(bioframe.is_contained(frame, views, raise_errors=False)),
    )
    return _check(
        "check_bioframe_contained",
        "contained",
        passed,
        len(frame),
        "an interval passes when it lies entirely within one view region (end-exclusive)",
    )


def check_bioframe_covering(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the intervals of a BED file cover every base of a view (is_covering)."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    passed = _run("is_covering", lambda: bool(bioframe.is_covering(frame, views)))
    return _check(
        "check_bioframe_covering",
        "covering",
        passed,
        len(frame),
        "passes when the union of source intervals spans each view region without gaps",
    )


def check_bioframe_overlapping(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report whether any two intervals of a BED file overlap (is_overlapping)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    passed = _run("is_overlapping", lambda: bool(bioframe.is_overlapping(frame)))
    return _check(
        "check_bioframe_overlapping",
        "overlapping",
        passed,
        len(frame),
        "passed is true when at least one pair of intervals shares a base",
    )


def check_bioframe_sorted(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report whether a BED file is already in bioframe sort order (is_sorted)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    passed = _run("is_sorted", lambda: bool(bioframe.is_sorted(frame)))
    return _check(
        "check_bioframe_sorted",
        "sorted",
        passed,
        len(frame),
        "passed is true when sort_bioframe_intervals would leave the file unchanged",
    )


def check_bioframe_tiling(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the intervals of a BED file tile every region of a view (is_tiling)."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    passed = _run(
        "is_tiling",
        lambda: bool(bioframe.is_tiling(frame, views, df_view_col=None)),
    )
    return _check(
        "check_bioframe_tiling",
        "tiling",
        passed,
        len(frame),
        "passes when the intervals abut from start to end of each view region, without overlap",
    )


def check_bioframe_viewframe(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that a BED file is a valid view frame: unique names, valid chrom, start, end."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    passed = _run("is_viewframe", lambda: bool(bioframe.is_viewframe(frame, raise_errors=False)))
    return _check(
        "check_bioframe_viewframe",
        "viewframe",
        passed,
        len(frame),
        "a view frame needs a name column (BED4 or wider) with unique names and valid intervals",
    )


# Region strings and coordinates.


def resolve_bioframe_region(
    region: str,
    *,
    chromsizes: Mapping[str, int] | None = None,
    check_bounds: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a region string into chrom, start, end, filling a chromosome-only name from sizes."""

    text = _text(region, "region", maximum=_REGION_MAX)
    sizes = None if chromsizes is None else _lengths(chromsizes, "chromsizes")
    bounds = _flag(check_bounds, "check_bounds")
    bioframe = _bioframe()
    chrom, start, end = _run(
        "parse_region",
        lambda: bioframe.parse_region(text, sizes, check_bounds=bounds),
    )
    if end is None:
        raise ConfigurationError("region has no end; pass chromsizes to fill the chromosome end.")
    begin, finish = int(start), int(end)
    return {
        "status": "completed",
        "operation": "resolve_bioframe_region",
        "region": text,
        "chrom": chrom,
        "start": begin,
        "end": finish,
        "length": finish - begin,
        "provenance": {
            "package": _package("bioframe"),
            "convention": "numbers are read as written, commas removed; a chromosome name alone "
            "spans 0 to chromsizes[chrom]",
        },
    }


def format_bioframe_region(
    chrom: str,
    start: int,
    end: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a chrom, start, end triple as a UCSC region string with to_ucsc_string."""

    name = _text(chrom, "chrom", maximum=256)
    begin = _integer(start, "start", 0, _COORD_MAX)
    finish = _integer(end, "end", 1, _COORD_MAX)
    if begin >= finish:
        raise ConfigurationError("start must be less than end.")
    bioframe = _bioframe()
    region = _run("to_ucsc_string", lambda: bioframe.to_ucsc_string((name, begin, finish)))
    return {
        "status": "completed",
        "operation": "format_bioframe_region",
        "region": region,
        "chrom": name,
        "start": begin,
        "end": finish,
        "provenance": {
            "package": _package("bioframe"),
            "convention": "start and end are written as given; no coordinate shift",
        },
    }


# Interval cleaning, ordering, and set operations.


def sanitize_bioframe_intervals(
    source: str | Path,
    *,
    start_exceed_end_action: str | None = None,
    drop_null: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clean BED3 to BED6 intervals with sanitize_bedframe, flipping or dropping start>end rows."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    action = _choice(start_exceed_end_action, "start_exceed_end_action", _SANITIZE_ACTIONS)
    drop = _flag(drop_null, "drop_null")
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    cleaned = _run(
        "sanitize_bedframe",
        lambda: bioframe.sanitize_bedframe(frame, drop_null=drop, start_exceed_end_action=action),
    )
    rows = _records(cleaned, limit, "Sanitized intervals")
    return {
        "status": "completed",
        "operation": "sanitize_bioframe_intervals",
        "start_exceed_end_action": action,
        "drop_null": drop,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "None raises on start>end; flip swaps start and end; drop removes "
            "those rows; null rows are removed only when drop_null is true",
        },
    }


def sort_bioframe_intervals(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort BED3 to BED6 intervals by chromosome name text, then start, with sort_bedframe."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    ordered = _run("sort_bedframe", lambda: bioframe.sort_bedframe(frame))
    rows = _records(ordered, limit, "Sorted intervals")
    return {
        "status": "completed",
        "operation": "sort_bioframe_intervals",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "chromosome names compare as text (chr10 before chr2), then start",
        },
    }


def setdiff_bioframe_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop each interval of source that overlaps any interval of other (setdiff)."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    reference = _bed(bioframe, other_path, limit)
    kept = _run("setdiff", lambda: bioframe.setdiff(frame, reference))
    rows = _records(kept, limit, "Remaining intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "setdiff_bioframe_intervals",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "whole intervals are removed; partly overlapping intervals are not "
            "clipped (subtract_bioframe_intervals clips)",
        },
    }


def trim_bioframe_intervals(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip each interval of source to the regions of a view with trim; drop what lies outside."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    trimmed = _run("trim", lambda: bioframe.trim(frame, view_df=views))
    rows = [row for row in _records(trimmed, limit, "Trimmed intervals") if row.get("chrom")]
    return {
        "status": "completed",
        "operation": "trim_bioframe_intervals",
        "intervals": rows,
        "interval_count": len(rows),
        "dropped_count": len(trimmed) - len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "each kept interval is clipped to the view region it overlaps; "
            "intervals with no overlap are dropped and counted in dropped_count",
        },
    }


def complement_bioframe_intervals(
    source: str | Path,
    view: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the parts of each view region that no interval of source covers (complement)."""

    path = input_path(source)
    view_path = input_path(view)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    gaps = _run("complement", lambda: bioframe.complement(frame, view_df=views))
    rows = _records(gaps, limit, "Gap intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "complement_bioframe_intervals",
        "gaps": rows,
        "gap_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "gaps are 0-based half-open; view_region names the view region each gap "
            "belongs to",
        },
    }


def count_bioframe_overlaps(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the intervals of other that overlap each interval of source (count_overlaps)."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    reference = _bed(bioframe, other_path, limit)
    counted = _run("count_overlaps", lambda: bioframe.count_overlaps(frame, reference))
    rows = _records(counted, limit, "Overlap counts")
    return {
        "status": "completed",
        "operation": "count_bioframe_overlaps",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "count is the number of other intervals overlapping each source "
            "interval; source order is kept",
        },
    }


def mark_bioframe_runs(
    source: str | Path,
    column: str,
    *,
    allow_overlaps: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number runs of neighbouring intervals that share one column value (mark_runs)."""

    path = input_path(source)
    field = _choice(column, "column", _RUN_COLUMNS)
    overlaps = _flag(allow_overlaps, "allow_overlaps")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    if field not in frame.columns:
        raise ConfigurationError(f"column {field} needs a BED file that has that column.")
    marked = _run(
        "mark_runs",
        lambda: bioframe.mark_runs(frame, field, allow_overlaps=overlaps),
    )
    rows = _records(marked, limit, "Run-marked intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "mark_bioframe_runs",
        "column": field,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "run numbers start at 0 per chromosome and count up while neighbouring "
            "intervals share the column value; overlaps are rejected unless allowed",
        },
    }


def merge_bioframe_runs(
    source: str | Path,
    column: str,
    *,
    allow_overlaps: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge neighbouring intervals that share one column value into single runs (merge_runs)."""

    path = input_path(source)
    field = _choice(column, "column", _RUN_COLUMNS)
    overlaps = _flag(allow_overlaps, "allow_overlaps")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    if field not in frame.columns:
        raise ConfigurationError(f"column {field} needs a BED file that has that column.")
    merged = _run(
        "merge_runs",
        lambda: bioframe.merge_runs(frame, field, allow_overlaps=overlaps),
    )
    rows = _records(merged, limit, "Merged runs", ordered=True)
    return {
        "status": "completed",
        "operation": "merge_bioframe_runs",
        "column": field,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "a run spans from the first start to the last end of neighbouring "
            "intervals with the same column value",
        },
    }


def pair_bioframe_intervals_by_distance(
    source: str | Path,
    *,
    min_sep: int,
    max_sep: int,
    min_intervening: int | None = None,
    max_intervening: int | None = None,
    relative_to: str = "midpoints",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List unique interval pairs whose separation falls in a window (pair_by_distance)."""

    path = input_path(source)
    low = _integer(min_sep, "min_sep", 0, _COORD_MAX)
    high = _integer(max_sep, "max_sep", 1, _COORD_MAX)
    if low >= high:
        raise ConfigurationError("min_sep must be less than max_sep.")
    near = _optional_int(min_intervening, "min_intervening")
    far = _optional_int(max_intervening, "max_intervening")
    anchor = _choice(relative_to, "relative_to", _RELATIVE_TO)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    pairs = _run(
        "pair_by_distance",
        lambda: bioframe.pair_by_distance(
            frame,
            min_sep=low,
            max_sep=high,
            min_intervening=near,
            max_intervening=far,
            relative_to=anchor,
        ),
    )
    rows = _records(
        pairs,
        limit,
        "Interval pairs",
        ordered=True,
        fields=("start_1", "end_1", "start_2", "end_2"),
    )
    return {
        "status": "completed",
        "operation": "pair_bioframe_intervals_by_distance",
        "relative_to": anchor,
        "pairs": rows,
        "pair_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "separation is measured between midpoints by default and is "
            "min_sep <= distance < max_sep; suffixes _1 and _2 mark the two intervals",
        },
    }


def select_bioframe_region(
    source: str | Path,
    region: str,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the intervals of a BED file that overlap a region string (select)."""

    path = input_path(source)
    text = _text(region, "region", maximum=_REGION_MAX)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    selected = _run("select", lambda: bioframe.select(frame, text))
    rows = _records(selected, limit, "Selected intervals")
    return {
        "status": "completed",
        "operation": "select_bioframe_region",
        "region": text,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "an interval is selected when it overlaps the region; a bare "
            "chromosome name selects every interval on it",
        },
    }


def assign_bioframe_views(
    source: str | Path,
    view: str | Path,
    *,
    drop_unassigned: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Label each interval of a BED file with the view region that contains it (assign_view)."""

    path = input_path(source)
    view_path = input_path(view)
    drop = _flag(drop_unassigned, "drop_unassigned")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    views = _view(bioframe, view_path, limit)
    assigned = _run(
        "assign_view",
        lambda: bioframe.assign_view(frame, views, drop_unassigned=drop),
    )
    rows = _records(assigned, limit, "View assignments")
    return {
        "status": "completed",
        "operation": "assign_bioframe_views",
        "drop_unassigned": drop,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "view_region holds the view name of each interval; null when no "
            "view region contains it",
        },
    }


def read_bioframe_table(
    source: str | Path,
    *,
    schema: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a tab-delimited BED3 to BED12 table with bioframe.read_table and a checked schema."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    widths = set()
    count = 0
    for _, fields in _data_lines(path):
        count += 1
        if count > limit:
            raise ConfigurationError(f"{path.name} has more than max_rows={limit} data lines.")
        widths.add(len(fields))
    if not count:
        raise ConfigurationError(f"{path.name} contains no data lines.")
    if len(widths) != 1:
        raise ConfigurationError(f"{path.name} rows need one column count; found {sorted(widths)}.")
    found = widths.pop()
    if schema is None:
        if found not in _BED_WIDTHS:
            raise ConfigurationError(
                f"{path.name} has {found} columns; BED3, BED4, BED5, BED6, BED9 and BED12 "
                "are supported."
            )
        name = f"bed{found}"
    else:
        name = _choice(schema, "schema", _BED_SCHEMAS)
        if found != int(name[3:]):
            raise ConfigurationError(f"schema {name} needs {name[3:]} columns; found {found}.")
    bioframe = _bioframe()
    frame = _run("read_table", lambda: bioframe.read_table(str(path), schema=name, comment="#"))
    rows = _records(frame, limit, "Table rows")
    return {
        "status": "completed",
        "operation": "read_bioframe_table",
        "schema": name,
        "rows": rows,
        "row_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "lines starting with # are skipped; track and browser lines are not",
        },
    }


def write_bioframe_bed(
    source: str | Path,
    output: str | Path,
    *,
    schema: str = "infer",
    require_sorted: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Re-write a BED file through bioframe.to_bed to a new path after validating its fields."""

    path = input_path(source)
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    dest_schema = _choice(schema, "schema", ("infer", *_BED_SCHEMAS))
    sorted_only = _flag(require_sorted, "require_sorted")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    try:
        _run(
            "to_bed",
            lambda: bioframe.to_bed(
                frame, path=str(target), schema=dest_schema, require_sorted=sorted_only
            ),
        )
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    return {
        "status": "completed",
        "operation": "write_bioframe_bed",
        "output": str(target),
        "schema": dest_schema,
        "record_count": len(frame),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "tab-separated with no header; require_sorted rejects unsorted rows",
        },
    }


# Sequence statistics and restriction digests.


def digest_bioframe_fasta(
    source: str | Path,
    enzyme: str,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut every FASTA sequence into restriction fragments with bioframe.digest (Biopython)."""

    name = _text(enzyme, "enzyme", maximum=64)
    if not name.isidentifier():
        raise ConfigurationError("enzyme must be a Biopython Restriction name, such as EcoRI.")
    limit = _rows_limit(max_rows)
    require("Bio.Restriction", "bio-formats")
    bioframe = _bioframe()
    records = _fasta_records(bioframe, source)
    fragments = _run("digest", lambda: bioframe.digest(records, name))
    # bioframe.digest adds 1 to Biopython's search position, which already points one base past
    # the cut (1-based). The 0-based boundary is therefore two less than bioframe reports.
    ends = fragments.groupby("chrom")["end"].transform("max")
    fragments = fragments.copy()
    fragments.loc[fragments["start"] > 0, "start"] -= _DIGEST_CUT_OFFSET
    fragments.loc[fragments["end"] < ends, "end"] -= _DIGEST_CUT_OFFSET
    rows = _records(fragments, limit, "Restriction fragments", ordered=True)
    return {
        "status": "completed",
        "operation": "digest_bioframe_fasta",
        "enzyme": name,
        "fragments": rows,
        "fragment_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "fragments tile each sequence from 0 to its length; interior cuts are "
            "0-based boundaries after the cut base (EcoRI cuts G^AATTC after the G)",
        },
    }


def sequence_gc_fraction(
    sequence: str,
    *,
    mapped_only: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the fraction of G and C in one DNA string with bioframe.seq_gc."""

    text = sequence_text(sequence)
    if len(text) > _SEQUENCE_MAX:
        raise ConfigurationError(f"sequence has more than {_SEQUENCE_MAX} bases.")
    mapped = _flag(mapped_only, "mapped_only")
    bioframe = _bioframe()
    value = _run("seq_gc", lambda: bioframe.seq_gc(text, mapped_only=mapped))
    fraction = None if math.isnan(float(value)) else float(value)
    return {
        "status": "completed",
        "operation": "sequence_gc_fraction",
        "length": len(text),
        "gc_fraction": fraction,
        "mapped_only": mapped,
        "provenance": {
            "package": _package("bioframe"),
            "convention": "G and C over all bases; with mapped_only, N is excluded from the "
            "denominator; null when no base is counted",
        },
    }


def frac_gc_bioframe_intervals(
    source: str | Path,
    fasta: str | Path,
    *,
    mapped_only: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add the GC fraction of each BED interval of a FASTA sequence with bioframe.frac_gc."""

    path = input_path(source)
    mapped = _flag(mapped_only, "mapped_only")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    records = _fasta_records(bioframe, fasta)
    scored = _run("frac_gc", lambda: bioframe.frac_gc(frame, records, mapped_only=mapped))
    scored = scored.rename(columns={"GC": "frac_gc"})
    rows = _records(scored, limit, "GC intervals")
    return {
        "status": "completed",
        "operation": "frac_gc_bioframe_intervals",
        "mapped_only": mapped,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "frac_gc is G plus C over the interval's bases; with mapped_only, N is "
            "excluded; null when the interval has no counted base",
        },
    }


def frac_mapped_bioframe_intervals(
    source: str | Path,
    fasta: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add the fraction of non-N bases of each BED interval of a FASTA with bioframe.frac_mapped."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    records = _fasta_records(bioframe, fasta)
    scored = _run("frac_mapped", lambda: bioframe.frac_mapped(frame, records))
    rows = _records(scored, limit, "Mapped intervals")
    return {
        "status": "completed",
        "operation": "frac_mapped_bioframe_intervals",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "frac_mapped is the non-N share of the sequence bases inside the "
            "interval; 0 when the interval lies outside the sequence",
        },
    }


def frac_gene_coverage_bioframe(
    source: str | Path,
    mrna: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Measure how much of each interval mRNA alignments cover (frac_gene_coverage, local table)."""

    path = input_path(source)
    mrna_path = input_path(mrna)
    limit = _rows_limit(max_rows)
    pandas = require("pandas", "bio-bioframe")
    bioframe = _bioframe()
    frame = _bed(bioframe, path, limit)
    table = _run("read mRNA table", lambda: pandas.read_table(mrna_path))
    missing = [column for column in _MRNA_COLUMNS if column not in table.columns]
    if missing:
        raise ConfigurationError(
            f"mRNA table needs columns {', '.join(_MRNA_COLUMNS)}; missing {', '.join(missing)}."
        )
    scored = _run("frac_gene_coverage", lambda: bioframe.frac_gene_coverage(frame, table))
    rows = _records(scored, limit, "Gene coverage intervals")
    return {
        "status": "completed",
        "operation": "frac_gene_coverage_bioframe",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "coverage counts bases covered by mRNA alignments; count is the number "
            "of overlapping mRNA alignments; the mRNA table is a local tab-separated file "
            "with a header",
        },
    }


def summarize_bioframe_fasta(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the sequence names and lengths of an indexed FASTA file with bioframe.load_fasta."""

    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    records = _fasta_records(bioframe, source)
    if len(records) > limit:
        raise ConfigurationError(
            f"FASTA has {len(records)} sequences; max_rows={limit} is the limit."
        )
    sequences = [
        {
            "name": str(name),
            # Lengths come from the pysam index; sequences are not read into memory.
            "length": int(record.ff.get_reference_length(record.ref)),
        }
        for name, record in records.items()
    ]
    return {
        "status": "completed",
        "operation": "summarize_bioframe_fasta",
        "sequences": sequences,
        "sequence_count": len(sequences),
        "total_length": sum(item["length"] for item in sequences),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "sequences keep FASTA order; the .fai index is created beside the file "
            "when missing",
        },
    }


# Genome assemblies, chromosome sizes, and arms.


def list_bioframe_assemblies(*, verbose: bool = True) -> dict[str, Any]:
    """List the genome assemblies whose metadata ships with bioframe (assemblies_available)."""

    bioframe = _bioframe()
    table = _run("assemblies_available", bioframe.assemblies_available)
    rows = [
        {
            "organism": _plain(item["organism"]),
            "provider": _plain(item["provider"]),
            "provider_build": _plain(item["provider_build"]),
            "release_year": _plain(item["release_year"]),
            "default_roles": [str(role) for role in item["default_roles"]],
            "default_units": [str(unit) for unit in item["default_units"]],
            "url": _plain(item["url"]),
        }
        for item in table.to_dict("records")
    ]
    rows.sort(key=lambda row: (str(row["provider"]), str(row["provider_build"])))
    return {
        "status": "completed",
        "operation": "list_bioframe_assemblies",
        "assemblies": rows,
        "assembly_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "metadata is bundled with bioframe; no network access is used",
        },
    }


def describe_bioframe_assembly(
    name: str,
    *,
    roles: Sequence[str] | str | None = None,
    units: Sequence[str] | str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the sequences of a bundled genome assembly, such as ucsc.hg38 (assembly_info)."""

    text = _text(name, "name", maximum=64)
    role_filter = _selector(roles, "roles")
    unit_filter = _selector(units, "units")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    info = _run(
        "assembly_info",
        lambda: bioframe.assembly_info(text, roles=role_filter, units=unit_filter),
    )
    sequences = _pandas_records(info.seqinfo, limit, "Assembly sequences")
    return {
        "status": "completed",
        "operation": "describe_bioframe_assembly",
        "name": text,
        "organism": info.organism,
        "provider": info.provider,
        "provider_build": info.provider_build,
        "release_year": _plain(info.release_year),
        "url": info.url,
        "sequence_count": len(sequences),
        "sequences": sequences,
        "provenance": {
            "package": _package("bioframe"),
            "convention": "default roles and units apply unless roles or units are given; "
            "sequences keep metadata file order",
        },
    }


def make_bioframe_chromarms(
    chromsizes: Mapping[str, int],
    midpoints: Mapping[str, int],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split each chromosome at its centromere midpoint into p and q arms (make_chromarms)."""

    sizes = _lengths(chromsizes, "chromsizes")
    if not isinstance(midpoints, Mapping) or not midpoints:
        raise ConfigurationError("midpoints must be a non-empty mapping from chrom to midpoint.")
    mids = {
        _text(chrom, "midpoints chromosome"): _integer(mid, "midpoint", 1, _COORD_MAX)
        for chrom, mid in midpoints.items()
    }
    if set(mids) != set(sizes):
        raise ConfigurationError("midpoints must name every chromosome in chromsizes exactly.")
    for chrom in sizes:
        if mids[chrom] >= sizes[chrom]:
            raise ConfigurationError(f"midpoint for {chrom} must lie inside its length.")
    limit = _rows_limit(max_rows)
    pandas = require("pandas", "bio-bioframe")
    bioframe = _bioframe()
    names = sorted(sizes)
    chroms = pandas.DataFrame({"chrom": names, "length": [sizes[name] for name in names]})
    centers = pandas.DataFrame({"chrom": names, "mid": [mids[name] for name in names]})
    arms = _run("make_chromarms", lambda: bioframe.make_chromarms(chroms, centers))
    rows = _records(arms, limit, "Chromosome arms", ordered=True)
    return {
        "status": "completed",
        "operation": "make_bioframe_chromarms",
        "arms": rows,
        "arm_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "each chromosome becomes <chrom>_p from 0 to its midpoint and <chrom>_q "
            "from the midpoint to its length",
        },
    }


def read_bioframe_chromsizes(
    source: str | Path,
    *,
    filter_chroms: bool = True,
    natsort: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local chrom.sizes or chromInfo file with bioframe.read_chromsizes."""

    path = input_path(source)
    keep = _flag(filter_chroms, "filter_chroms")
    natural = _flag(natsort, "natsort")
    limit = _rows_limit(max_rows)
    bioframe = _bioframe()
    sizes = _run(
        "read_chromsizes",
        lambda: bioframe.read_chromsizes(str(path), filter_chroms=keep, natsort=natural),
    )
    if len(sizes) > limit:
        raise ConfigurationError(
            f"chromsizes has {len(sizes)} rows; max_rows={limit} is the limit."
        )
    chromosomes = [{"chrom": str(name), "length": int(length)} for name, length in sizes.items()]
    return {
        "status": "completed",
        "operation": "read_bioframe_chromsizes",
        "chromsizes": chromosomes,
        "chromosome_count": len(chromosomes),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "filter_chroms keeps chr1 to chr22, chrX, chrY and chrM; natsort orders "
            "names naturally",
        },
    }


# Indexed and alignment readers.


def read_bioframe_tabix(
    source: str | Path,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a bgzipped, tabix-indexed BED file, whole or by region (read_tabix via pysam)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    if not any(Path(f"{path}{suffix}").is_file() for suffix in (".tbi", ".csi")):
        raise ConfigurationError(
            f"{path.name} has no .tbi or .csi index; build one with pysam.tabix_index."
        )
    name, begin, finish = _region_bounds(chrom, start, end)
    require("pysam", "bio-reads")
    bioframe = _bioframe()
    frame = _run("read_tabix", lambda: bioframe.read_tabix(str(path), name, begin, finish))
    width = len(frame.columns)
    if list(frame.columns) == list(range(width)) and width in _BED_WIDTHS:
        frame.columns = list(bioframe.SCHEMAS[f"bed{width}"])
    rows = _records(frame, limit, "Tabix intervals")
    return {
        "status": "completed",
        "operation": "read_bioframe_tabix",
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "start and end are 0-based half-open; a region query returns records "
            "overlapping chrom:start-end",
        },
    }


def read_bioframe_alignments(
    source: str | Path,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read SAM or BAM alignment records with bioframe.read_alignments (pysam)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    name, begin, finish = _region_bounds(chrom, start, end)
    require("pysam", "bio-reads")
    bioframe = _bioframe()
    frame = _run(
        "read_alignments",
        lambda: bioframe.read_alignments(str(path), name, begin, finish),
    )
    rows = _pandas_records(frame, limit, "Alignment records")
    return {
        "status": "completed",
        "operation": "read_bioframe_alignments",
        "alignments": rows,
        "alignment_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "POS is 0-based (the SAM column minus one); region queries need an "
            "indexed BAM",
        },
    }


def read_bioframe_bigbed(
    source: str | Path,
    chrom: str,
    *,
    start: int | None = None,
    end: int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the intervals of one chromosome from a local bigBed file (read_bigbed via pyBigWig)."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    name, begin, finish = _region_bounds(chrom, start, end)
    require("pyBigWig", "bio-pybigwig")
    bioframe = _bioframe()
    frame = _run(
        "read_bigbed",
        lambda: bioframe.read_bigbed(str(path), name, start=begin, end=finish, engine="pybigwig"),
    )
    rows = _records(frame, limit, "bigBed intervals")
    return {
        "status": "completed",
        "operation": "read_bioframe_bigbed",
        "chrom": name,
        "intervals": rows,
        "interval_count": len(rows),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "start and end are 0-based half-open; rest holds the remaining BED "
            "fields as tab-separated text",
        },
    }


__all__ = [
    "assign_bioframe_views",
    "check_bioframe_cataloged",
    "check_bioframe_contained",
    "check_bioframe_covering",
    "check_bioframe_overlapping",
    "check_bioframe_sorted",
    "check_bioframe_tiling",
    "check_bioframe_viewframe",
    "complement_bioframe_intervals",
    "count_bioframe_overlaps",
    "describe_bioframe_assembly",
    "digest_bioframe_fasta",
    "format_bioframe_region",
    "frac_gc_bioframe_intervals",
    "frac_gene_coverage_bioframe",
    "frac_mapped_bioframe_intervals",
    "list_bioframe_assemblies",
    "make_bioframe_chromarms",
    "mark_bioframe_runs",
    "merge_bioframe_runs",
    "pair_bioframe_intervals_by_distance",
    "read_bioframe_alignments",
    "read_bioframe_bigbed",
    "read_bioframe_chromsizes",
    "read_bioframe_tabix",
    "read_bioframe_table",
    "resolve_bioframe_region",
    "sanitize_bioframe_intervals",
    "select_bioframe_region",
    "sequence_gc_fraction",
    "setdiff_bioframe_intervals",
    "sort_bioframe_intervals",
    "summarize_bioframe_fasta",
    "trim_bioframe_intervals",
    "validate_bioframe_bed",
    "write_bioframe_bed",
]

_register_public_verbose(globals(), tuple(__all__))
