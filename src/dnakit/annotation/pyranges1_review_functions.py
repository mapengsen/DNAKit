"""Review wrappers for pyranges1 readers, interval concatenation, ORFs and overlap statistics.

- Readers: ``read_pyranges1_gff()``, ``read_pyranges1_gff3()``, ``read_pyranges1_gtf_full()`` and
  ``read_pyranges1_bigwig()`` return one JSON record per genomic feature or signal interval.
- Intervals: ``concatenate_pyranges1_intervals()`` stacks several BED files into one sorted list.
- ORFs: ``calculate_pyranges1_cds_frame()`` assigns the reading frame of CDS features grouped by an
  attribute; ``extend_pyranges1_orfs()`` extends BED6 CDS intervals to open reading frames using a
  FASTA file.
- Statistics: ``calculate_pyranges1_forbes()``, ``calculate_pyranges1_mcc()`` and
  ``calculate_pyranges1_relative_distance()`` compare interval sets.

Coordinates are 0-based half-open, as pyranges1 reports them. pyranges1 is imported only when a
function runs; a missing package raises ``BackendUnavailableError`` with the install hint.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _integer,
    _lengths,
    _package,
    _plain,
    _rows_limit,
    _scan_bed,
    _text,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-pyranges1"
_STRAND_BEHAVIOR = ("auto", "same", "opposite", "ignore")
_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,49}\Z")
_MAX_BIGWIG_BYTES = 200_000_000
_RENAME = {
    "Chromosome": "chrom",
    "Start": "start",
    "End": "end",
    "Strand": "strand",
    "Value": "value",
}


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _identifier(value: object, name: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise ConfigurationError(f"{name} must be an identifier such as 'transcript_id'.")
    return value


def _call(operation: str, function: Callable[[], Any]) -> Any:
    """Run a pyranges1 call and report its failures as input errors."""

    try:
        return function()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"pyranges1 could not {operation}: {exc}") from exc


def _sort_key(record: Mapping[str, object]) -> tuple[str, int, int, str]:
    return (
        str(record.get("chrom", "")),
        int(record.get("start") or 0),
        int(record.get("end") or 0),
        json.dumps(record, sort_keys=True, default=str),
    )


def _frame_records(frame: Any, limit: int, what: str) -> list[dict[str, Any]]:
    """Rows of a PyRanges table as JSON records with chrom, start, end and strand names."""

    if len(frame) > limit:
        raise ConfigurationError(f"{what} has {len(frame)} rows; max_rows={limit} is the limit.")
    records = []
    for row in frame.to_dict("records"):
        records.append(
            {_RENAME.get(str(key), str(key)): _plain(value) for key, value in row.items()}
        )
    return sorted(records, key=_sort_key)


def _result(operation: str, fields: Mapping[str, object], convention: str) -> dict[str, object]:
    return {
        "status": "completed",
        "operation": operation,
        "provenance": {"package": _package("pyranges1"), "convention": convention},
        **fields,
    }


def _record_fields(records: list[dict[str, Any]]) -> dict[str, object]:
    return {"records": records, "record_count": len(records)}


_GTF_CONVENTION = (
    "Coordinates are 0-based half-open; attribute columns keep the names written by the reader. "
    "Start is shifted by one from the 1-based file, as pyranges1 does."
)


def read_pyranges1_gff(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
    verbose: bool = True,
) -> dict[str, object]:
    """Read a GTF or GFF file with pyranges1 into one record per feature with its attributes."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    duplicates = _flag(duplicate_attr, "duplicate_attr")
    ignore = _flag(ignore_bad, "ignore_bad")
    pyranges1 = require("pyranges1", _EXTRA)
    frame = _call(
        "read the GTF/GFF file",
        lambda: pyranges1.read_gff(
            source, nrows=limit + 1, duplicate_attr=duplicates, ignore_bad=ignore
        ),
    )
    records = _frame_records(frame, limit, "GTF/GFF features")
    return _result(
        "read_pyranges1_gff",
        _record_fields(records),
        _GTF_CONVENTION,
    )


def read_pyranges1_gff3(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Read a GFF3 file with pyranges1; attribute keys become record fields."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pyranges1 = require("pyranges1", _EXTRA)
    frame = _call(
        "read the GFF3 file",
        lambda: pyranges1.read_gff3(source, nrows=limit + 1),
    )
    records = _frame_records(frame, limit, "GFF3 features")
    return _result("read_pyranges1_gff3", _record_fields(records), _GTF_CONVENTION)


def read_pyranges1_gtf_full(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
    verbose: bool = True,
) -> dict[str, object]:
    """Read a GTF file with pyranges1's full attribute parser, one record per feature."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    duplicates = _flag(duplicate_attr, "duplicate_attr")
    ignore = _flag(ignore_bad, "ignore_bad")
    readers = require("pyranges1.readers", _EXTRA)
    frame = _call(
        "read the GTF file",
        lambda: readers.read_gtf_full(
            source, nrows=limit + 1, duplicate_attr=duplicates, ignore_bad=ignore
        ),
    )
    records = _frame_records(frame, limit, "GTF features")
    return _result("read_pyranges1_gtf_full", _record_fields(records), _GTF_CONVENTION)


def read_pyranges1_bigwig(
    path: str | Path,
    *,
    max_file_bytes: int = 50_000_000,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Read every interval of a bigWig track as chrom, start, end and value records."""

    source = input_path(path)
    size = source.stat().st_size
    byte_limit = _integer(max_file_bytes, "max_file_bytes", 1, _MAX_BIGWIG_BYTES)
    if size > byte_limit:
        raise ConfigurationError(f"bigWig file has {size} bytes; max_file_bytes is {byte_limit}.")
    limit = _rows_limit(max_rows)
    require("pyBigWig", "bio-pybigwig")
    pyranges1 = require("pyranges1", _EXTRA)
    frame = _call("read the bigWig file", lambda: pyranges1.read_bigwig(source))
    records = _frame_records(frame, limit, "bigWig intervals")
    return _result(
        "read_pyranges1_bigwig",
        _record_fields(records),
        "Intervals are 0-based half-open with the track value. The whole file is read before the "
        "max_rows check, so max_file_bytes bounds the input size.",
    )


def concatenate_pyranges1_intervals(
    paths: Sequence[str | Path],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Concatenate several BED files into one sorted list of intervals with pyranges1."""

    if isinstance(paths, (str, bytes)) or not isinstance(paths, Sequence) or not paths:
        raise ConfigurationError("paths must be a non-empty list of BED file paths.")
    if len(paths) > 50:
        raise ConfigurationError("paths may contain at most 50 files.")
    limit = _rows_limit(max_rows)
    sources = [input_path(item) for item in paths]
    for source in sources:
        _scan_bed(source, limit)
    pyranges1 = require("pyranges1", _EXTRA)
    frames = [pyranges1.read_bed(str(source)) for source in sources]
    combined = _call("concatenate the intervals", lambda: pyranges1.concat(frames))
    records = _frame_records(combined, limit, "concatenated intervals")
    return _result(
        "concatenate_pyranges1_intervals",
        {**_record_fields(records), "source_count": len(sources)},
        "BED files are read with pyranges1.read_bed and stacked; output is sorted by chromosome, "
        "start, end and the remaining columns.",
    )


def calculate_pyranges1_cds_frame(
    path: str | Path,
    *,
    feature: str = "CDS",
    group_by: str | Sequence[str] = ("transcript_id",),
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Assign the reading frame of each CDS feature from the CDS lengths in its group."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    if not isinstance(feature, str) or not feature.strip() or len(feature) > 50:
        raise ConfigurationError("feature must be a non-empty string of at most 50 characters.")
    groups = [group_by] if isinstance(group_by, str) else list(group_by)
    if not 1 <= len(groups) <= 5:
        raise ConfigurationError("group_by must name between 1 and 5 attribute columns.")
    names = [_identifier(item, "group_by") for item in groups]
    pyranges1 = require("pyranges1", _EXTRA)
    orfs = require("pyranges1.orfs", _EXTRA)
    features = _call("read the GTF/GFF file", lambda: pyranges1.read_gff(source, nrows=limit + 1))
    selected = features[features["Feature"] == feature]
    if len(selected) == 0:
        raise ConfigurationError(f"No feature named {feature!r} in {source.name}.")
    missing = [name for name in names if name not in selected.columns]
    if missing:
        raise ConfigurationError(f"group_by columns not found: {missing}.")
    framed = _call(
        "calculate CDS frames",
        lambda: orfs.calculate_frame(selected, group_by=names, frame_col="cds_frame"),
    )
    rows = []
    for row in framed.to_dict("records"):
        record: dict[str, Any] = {
            "chrom": _plain(row["Chromosome"]),
            "start": _plain(row["Start"]),
            "end": _plain(row["End"]),
            "strand": _plain(row["Strand"]),
        }
        for name in names:
            record[name] = _plain(row[name])
        record["frame"] = _plain(row["cds_frame"])
        rows.append(record)
    rows.sort(key=_sort_key)
    return _result(
        "calculate_pyranges1_cds_frame",
        _record_fields(rows),
        "frame is the sum of the earlier CDS lengths in the same group modulo 3, as pyranges1 "
        "computes it. Intervals are sorted within each group by start.",
    )


def extend_pyranges1_orfs(
    path: str | Path,
    fasta_path: str | Path,
    *,
    group_by: str = "Name",
    direction: str | None = None,
    keep_off_bounds: bool = False,
    chunk_size: int = 900,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Extend BED6 CDS intervals to complete open reading frames using a FASTA sequence."""

    source = input_path(path)
    fasta = input_path(fasta_path)
    limit = _rows_limit(max_rows)
    width = _scan_bed(source, limit)
    if width != 6:
        raise ConfigurationError("CDS input must be BED6 (with a Name and a Strand column).")
    name = _identifier(group_by, "group_by")
    if direction not in (None, "up", "down"):
        raise ConfigurationError("direction must be None (both ways), 'up' or 'down'.")
    keep = _flag(keep_off_bounds, "keep_off_bounds")
    step = _integer(chunk_size, "chunk_size", 3, 3_000)
    if step % 3:
        raise ConfigurationError("chunk_size must be a multiple of three.")
    require("pyfaidx", "bio-pyfaidx")
    pyranges1 = require("pyranges1", _EXTRA)
    orfs = require("pyranges1.orfs", _EXTRA)
    cds = pyranges1.read_bed(str(source))
    extended = _call(
        "extend the ORFs",
        lambda: orfs.extend_orfs(
            cds,
            str(fasta),
            group_by=name,
            direction=direction,
            keep_off_bounds=keep,
            chunk_size=step,
        ),
    )
    records = _frame_records(extended, limit, "ORF intervals")
    return _result(
        "extend_pyranges1_orfs",
        _record_fields(records),
        "Intervals are extended to the nearest in-frame ATG upstream and stop codon (TAG, TGA, "
        "TAA) downstream. pyfaidx writes a .fai index next to the FASTA file on first use.",
    )


def calculate_pyranges1_forbes(
    source: str | Path,
    other: str | Path,
    *,
    chromsizes: Mapping[str, int],
    strand_behavior: str = "auto",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Forbes coefficient of two BED interval sets on a genome of given chromosome sizes."""

    left_path = input_path(source)
    right_path = input_path(other)
    limit = _rows_limit(max_rows)
    if strand_behavior not in _STRAND_BEHAVIOR:
        raise ConfigurationError(f"strand_behavior must be one of {', '.join(_STRAND_BEHAVIOR)}.")
    sizes = _lengths(chromsizes, "chromsizes")
    _scan_bed(left_path, limit)
    _scan_bed(right_path, limit)
    pyranges1 = require("pyranges1", _EXTRA)
    stats = require("pyranges1.stats", _EXTRA)
    left = pyranges1.read_bed(str(left_path))
    right = pyranges1.read_bed(str(right_path))
    value = _call(
        "compute the Forbes coefficient",
        lambda: stats.forbes(left, right, chromsizes=sizes, strand_behavior=strand_behavior),
    )
    return _result(
        "calculate_pyranges1_forbes",
        {"forbes": _plain(float(value)), "strand_behavior": strand_behavior},
        "Forbes coefficient from pyranges1.stats.forbes; chromsizes give the genome length.",
    )


def calculate_pyranges1_mcc(
    paths: Sequence[str | Path],
    labels: Sequence[str],
    *,
    use_strand: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Matthews correlation coefficient of base-pair overlaps between labelled BED files."""

    if (
        isinstance(paths, (str, bytes))
        or not isinstance(paths, Sequence)
        or not 2 <= len(paths) <= 10
    ):
        raise ConfigurationError("paths must be a list of 2 to 10 BED file paths.")
    if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence):
        raise ConfigurationError("labels must be a list of strings, one per BED file.")
    if len(labels) != len(paths):
        raise ConfigurationError("labels and paths must have the same length.")
    names = [_text(item, "label", maximum=100) for item in labels]
    if len(set(names)) != len(names):
        raise ConfigurationError("labels must be unique.")
    flag = _flag(use_strand, "use_strand")
    limit = _rows_limit(max_rows)
    sources = [input_path(item) for item in paths]
    for source in sources:
        _scan_bed(source, limit)
    pyranges1 = require("pyranges1", _EXTRA)
    stats = require("pyranges1.stats", _EXTRA)
    frames = [pyranges1.read_bed(str(source)) for source in sources]
    table = _call(
        "compute the MCC table",
        lambda: stats.mcc(frames, labels=names, use_strand=flag),
    )
    rows = [
        {
            "truth_label": str(row["T"]),
            "prediction_label": str(row["F"]),
            "tp": _plain(row["TP"]),
            "fp": _plain(row["FP"]),
            "tn": _plain(row["TN"]),
            "fn": _plain(row["FN"]),
            "mcc": _plain(row["MCC"]),
        }
        for row in table.to_dict("records")
    ]
    rows.sort(key=lambda row: (row["truth_label"], row["prediction_label"]))
    return _result(
        "calculate_pyranges1_mcc",
        {"records": rows, "record_count": len(rows), "use_strand": flag},
        "Each row treats truth_label as the reference set and prediction_label as the predicted "
        "set; counts are base pairs from pyranges1.stats.mcc.",
    )


def calculate_pyranges1_relative_distance(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, object]:
    """Spatial correlation between two BED interval sets as relative-distance bins."""

    left_path = input_path(source)
    right_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(left_path, limit)
    _scan_bed(right_path, limit)
    pyranges1 = require("pyranges1", _EXTRA)
    stats = require("pyranges1.stats", _EXTRA)
    left = pyranges1.read_bed(str(left_path))
    right = pyranges1.read_bed(str(right_path))
    table = _call(
        "compute relative distances",
        lambda: stats.relative_distance(left, right),
    )
    rows = [
        {
            "reldist": _plain(row["reldist"]),
            "count": _plain(row["count"]),
            "total": _plain(row["total"]),
            "fraction": _plain(row["fraction"]),
        }
        for row in table.to_dict("records")
    ]
    return _result(
        "calculate_pyranges1_relative_distance",
        {"records": rows, "record_count": len(rows)},
        "Rows come from pyranges1.stats.relative_distance: reldist, count, total and fraction.",
    )


__all__ = [
    "calculate_pyranges1_cds_frame",
    "calculate_pyranges1_forbes",
    "calculate_pyranges1_mcc",
    "calculate_pyranges1_relative_distance",
    "concatenate_pyranges1_intervals",
    "extend_pyranges1_orfs",
    "read_pyranges1_bigwig",
    "read_pyranges1_gff",
    "read_pyranges1_gff3",
    "read_pyranges1_gtf_full",
]

_register_public_verbose(globals(), tuple(__all__))
