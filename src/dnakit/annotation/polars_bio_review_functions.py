"""Review wrappers for polars-bio header, schema, SQL, alignment, FASTQ, bigBed, and A2M/A3M files.

Each public function reads or writes local files with ``polars_bio``, runs one upstream operation,
and returns a JSON-serializable dictionary with ``status``, ``operation``, result fields, and
``provenance``. Coordinates follow the BED convention (0-based, half-open) unless a docstring
says otherwise. Rows are sorted by their values so repeated runs return identical lists, except
for ``query_genomic_tables`` when the SQL text has an ``ORDER BY``; then the query order is kept.

Functions are importable as ``dnakit.annotation.polars_bio_review_functions.<name>``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_MAX_QUERY = 10_000
_MAX_SQL_TABLES = 16
_BED_CONVENTION = "0-based half-open coordinates (BED)"
_TABLE_NAME = re.compile(r"^[a-z_][a-z0-9_]{0,62}$")
_SQL_START = re.compile(r"^(select|with)\b", re.IGNORECASE)
_ORDER_BY = re.compile(r"\border\s+by\b", re.IGNORECASE)
_ALIGNMENT_FORMATS = ("sam", "bam", "cram")
_HEADER_READERS = {"vcf": "scan_vcf", "bcf": "scan_bcf", "sam": "scan_sam", "bam": "scan_bam"}
_DESCRIBERS = {
    "vcf": "describe_vcf",
    "bcf": "describe_bcf",
    "sam": "describe_sam",
    "bam": "describe_bam",
    "cool": "describe_cool",
    "sto": "describe_sto",
}
_ALIGNMENT_SEQUENCE_READERS = {"a2m": "read_a2m", "a3m": "read_a3m"}
_SQL_REGISTRARS = {
    "bed": "register_bed",
    "gtf": "register_gtf",
    "gff": "register_gff",
    "vcf": "register_vcf",
    "bcf": "register_bcf",
    "sam": "register_sam",
    "bam": "register_bam",
    "fasta": "register_fasta",
    "fastq": "register_fastq",
    "a2m": "register_a2m",
    "a3m": "register_a3m",
    "sto": "register_sto",
    "pairs": "register_pairs",
    "cool": "register_cool",
    "bigwig": "register_bigwig",
    "bigbed": "register_bigbed",
}


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int = 512) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _names(value: Sequence[str] | None, name: str) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise ConfigurationError(f"{name} must be a list of field names.")
    items = [_text(item, name, maximum=128) for item in value]
    if not items:
        raise ConfigurationError(f"{name} must not be empty when given.")
    return items


def _new_output(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output file already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _plain(value: Any) -> Any:
    """Convert numpy, polars, and missing values to JSON-compatible Python objects."""

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
        elif isinstance(value, list):
            keys.append((2, str(value)))
        else:
            keys.append((2, str(value)))
    return tuple(keys)


def _ordered(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(records, key=_record_order)


def _records(
    frame: Any, limit: int, what: str, *, keep_order: bool = False
) -> list[dict[str, Any]]:
    """Materialize a polars frame as JSON-compatible records within the row limit."""

    if hasattr(frame, "collect"):
        frame = frame.collect()
    if frame.height > limit:
        raise ConfigurationError(f"{what} has {frame.height} rows; max_rows={limit} is the limit.")
    rows = [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]
    return rows if keep_order else _ordered(rows)


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("polars-bio"), "convention": convention}


# Headers and schemas.


def read_file_header(
    path: str | Path,
    *,
    file_format: str = "vcf",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the header of a VCF, BCF, SAM, or BAM file as polars-bio parses it."""

    source = input_path(path)
    fmt = _choice(file_format, "file_format", tuple(_HEADER_READERS))
    pb = require("polars_bio", "bio-polarsbio")
    frame = getattr(pb, _HEADER_READERS[fmt])(str(source))
    metadata = pb.get_metadata(frame)
    header = _plain(metadata.get("header") or {})
    header.pop("_datafusion_table_name", None)
    return {
        "status": "completed",
        "operation": "read_file_header",
        "file_format": fmt,
        "coordinate_system_zero_based": _plain(metadata.get("coordinate_system_zero_based")),
        "header": header,
        "provenance": _provenance("header values as parsed by polars-bio; the lazy scan is used"),
    }


def describe_genomic_file(
    path: str | Path,
    *,
    file_format: str = "vcf",
    reference_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe the schema or header entries of a VCF, BCF, SAM, BAM, CRAM, cool, or STO file."""

    source = input_path(path)
    fmt = _choice(file_format, "file_format", ("vcf", "bcf", "sam", "bam", "cram", "cool", "sto"))
    if fmt == "cram":
        if reference_path is None:
            raise ConfigurationError("reference_path is required to describe a CRAM file.")
        reference = input_path(reference_path)
    elif reference_path is not None:
        raise ConfigurationError("reference_path applies only to CRAM files.")
    pb = require("polars_bio", "bio-polarsbio")
    if fmt == "cram":
        frame = pb.describe_cram(str(source), reference_path=str(reference))
    else:
        frame = getattr(pb, _DESCRIBERS[fmt])(str(source))
    entries = _records(frame, _DEFAULT_ROWS, "Description")
    return {
        "status": "completed",
        "operation": "describe_genomic_file",
        "file_format": fmt,
        "entry_count": len(entries),
        "entries": entries,
        "provenance": _provenance(
            "schema or header entries; SAM, BAM, and CRAM tag types come from a sample of records"
        ),
    }


# Alignments and interval records.


def read_cram_records(
    path: str | Path,
    reference_path: str | Path,
    *,
    tag_fields: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read CRAM alignments, decoded against a reference FASTA, into one record per read."""

    source = input_path(path)
    reference = input_path(reference_path)
    tags = _names(tag_fields, "tag_fields")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_cram(
        str(source),
        reference_path=str(reference),
        tag_fields=tags,
        use_zero_based=True,
    )
    alignments = _records(frame, limit, "CRAM file")
    return {
        "status": "completed",
        "operation": "read_cram_records",
        "record_count": len(alignments),
        "alignments": alignments,
        "provenance": _provenance(_BED_CONVENTION + "; CRAM decoded with the reference FASTA"),
    }


def read_bigbed_intervals(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every interval of a bigBed file with its tab-separated fields after the coordinates."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_bigbed(str(source), schema="rest", use_zero_based=True)
    rows = _records(frame, limit, "bigBed file")
    intervals = []
    for row in rows:
        rest = row.pop("rest", None)
        row["fields"] = rest.split("\t") if rest else []
        intervals.append(row)
    return {
        "status": "completed",
        "operation": "read_bigbed_intervals",
        "interval_count": len(intervals),
        "intervals": _ordered(intervals),
        "provenance": _provenance(_BED_CONVENTION + "; fields are the stored columns after end"),
    }


def read_alignment_sequences(
    path: str | Path,
    *,
    file_format: str = "a2m",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an A2M or A3M multiple-sequence alignment, one record per aligned sequence."""

    source = input_path(path)
    fmt = _choice(file_format, "file_format", tuple(_ALIGNMENT_SEQUENCE_READERS))
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = getattr(pb, _ALIGNMENT_SEQUENCE_READERS[fmt])(str(source))
    rows = _records(frame, limit, "Alignment file")
    lengths = {len(row.get("sequence") or "") for row in rows}
    return {
        "status": "completed",
        "operation": "read_alignment_sequences",
        "file_format": fmt,
        "sequence_count": len(rows),
        "alignment_columns": lengths.pop() if len(lengths) == 1 else None,
        "sequences": rows,
        "provenance": _provenance(
            "gaps '-' and '.' and A3M lowercase insertions are kept as stored"
        ),
    }


# Writers.


def write_fastq_reads(
    reads: Sequence[Mapping[str, Any]],
    destination: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write reads with name, sequence, and Phred+33 quality_scores to a new FASTQ file."""

    target = _new_output(destination)
    if isinstance(reads, (str, bytes, Mapping)) or not isinstance(reads, Iterable):
        raise ConfigurationError("reads must be a list of records with name and sequence.")
    records = list(reads)
    if not records or len(records) > _MAX_ROWS:
        raise ConfigurationError(f"reads must hold between 1 and {_MAX_ROWS} records.")
    names: list[str] = []
    descriptions: list[str | None] = []
    sequences: list[str] = []
    qualities: list[str] = []
    for index, read in enumerate(records, 1):
        if not isinstance(read, Mapping):
            raise ConfigurationError(f"read {index} must be a mapping.")
        name = _text(read.get("name"), "read name", maximum=256)
        if any(character.isspace() for character in name):
            raise ConfigurationError(f"read name {name!r} cannot contain whitespace.")
        description = read.get("description")
        if description is not None:
            description = _text(description, "description", maximum=1024)
            if "\n" in description or "\r" in description:
                raise ConfigurationError("description cannot contain line breaks.")
        sequence = sequence_text(read.get("sequence"))
        quality = read.get("quality_scores")
        if (
            not isinstance(quality, str)
            or len(quality) != len(sequence)
            or any(not "!" <= character <= "~" for character in quality)
        ):
            raise ConfigurationError(
                f"read {name!r} needs quality_scores of Phred+33 characters, one per base."
            )
        names.append(name)
        descriptions.append(description)
        sequences.append(sequence)
        qualities.append(quality)
    pb = require("polars_bio", "bio-polarsbio")
    pl = require("polars", "bio-polarsbio")
    frame = pl.DataFrame(
        {
            "name": names,
            "description": pl.Series("description", descriptions, dtype=pl.String),
            "sequence": sequences,
            "quality_scores": qualities,
        }
    )
    written = int(pb.write_fastq(frame, str(target)))
    return {
        "status": "completed",
        "operation": "write_fastq_reads",
        "records_written": written,
        "destination": str(target),
        "provenance": _provenance(
            "sequences are written uppercase; qualities are written as given"
        ),
    }


def write_alignment_file(
    source: str | Path,
    destination: str | Path,
    *,
    source_format: str = "sam",
    output_format: str = "bam",
    reference_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read SAM, BAM, or CRAM alignments and write them as SAM, BAM, or CRAM in record order."""

    source_path = input_path(source)
    src_fmt = _choice(source_format, "source_format", _ALIGNMENT_FORMATS)
    out_fmt = _choice(output_format, "output_format", _ALIGNMENT_FORMATS)
    reference: Path | None = None
    if "cram" in (src_fmt, out_fmt):
        if reference_path is None:
            raise ConfigurationError("reference_path is required for CRAM input or output.")
        reference = input_path(reference_path)
    target = _new_output(destination)
    pb = require("polars_bio", "bio-polarsbio")
    if src_fmt == "cram":
        frame = pb.read_cram(str(source_path), reference_path=str(reference), use_zero_based=True)
    elif src_fmt == "sam":
        frame = pb.read_sam(str(source_path), use_zero_based=True)
    else:
        frame = pb.read_bam(str(source_path), use_zero_based=True)
    if out_fmt == "cram":
        written = pb.write_cram(frame, str(target), reference_path=str(reference))
    elif out_fmt == "sam":
        written = pb.write_sam(frame, str(target))
    else:
        written = pb.write_bam(frame, str(target))
    return {
        "status": "completed",
        "operation": "write_alignment_file",
        "source_format": src_fmt,
        "output_format": out_fmt,
        "records_written": int(written),
        "destination": str(target),
        "provenance": _provenance(
            _BED_CONVENTION + "; header lines come from the source; @HD SO is set by the writer"
        ),
    }


# SQL over local genomic files.


def _table_name(value: object) -> str:
    name = _text(value, "table name", maximum=63)
    if not _TABLE_NAME.match(name):
        raise ConfigurationError(
            f"table name {name!r} must be lowercase letters, digits, or underscores."
        )
    return name


def _select_statement(value: object) -> str:
    text = _text(value, "query", maximum=_MAX_QUERY).strip()
    body = text[:-1].rstrip() if text.endswith(";") else text
    if ";" in body or not _SQL_START.match(body):
        raise ConfigurationError(
            "query must be one SELECT or WITH statement over the registered tables."
        )
    return body


def _session_has_table(pb: Any, name: str) -> bool:
    try:
        pb.ctx.table(name)
    except KeyError:
        return False
    return True


def query_genomic_tables(
    tables: Mapping[str, tuple[str | Path, str]],
    query: str,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one read-only SQL SELECT over local genomic files, each registered under a table name.

    ``tables`` maps a lowercase table name to a ``(path, file_format)`` pair. The tables are
    registered in polars-bio's session for this call only and removed afterwards, even on error.
    """

    if not isinstance(tables, Mapping) or not 1 <= len(tables) <= _MAX_SQL_TABLES:
        raise ConfigurationError(f"tables must map 1 to {_MAX_SQL_TABLES} table names to files.")
    sources: list[tuple[str, Path, str]] = []
    for alias, spec in tables.items():
        name = _table_name(alias)
        if not isinstance(spec, (tuple, list)) or len(spec) != 2:
            raise ConfigurationError(f"table {name!r} must be a (path, file_format) pair.")
        path = input_path(spec[0])
        fmt = _choice(spec[1], f"file format of {name}", tuple(_SQL_REGISTRARS))
        sources.append((name, path, fmt))
    statement = _select_statement(query)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    for name, _, _ in sources:
        if _session_has_table(pb, name):
            raise ConfigurationError(f"A table named {name!r} is already registered.")
    registered: list[str] = []
    try:
        for name, path, fmt in sources:
            getattr(pb, _SQL_REGISTRARS[fmt])(str(path), name=name)
            registered.append(name)
        materialized = pb.sql(statement).collect()
        columns = list(materialized.columns)
        rows = _records(
            materialized,
            limit,
            "SQL result",
            keep_order=bool(_ORDER_BY.search(statement)),
        )
    finally:
        for name in registered:
            pb.ctx.deregister_table(name)
    return {
        "status": "completed",
        "operation": "query_genomic_tables",
        "table_names": [name for name, _, _ in sources],
        "columns": columns,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance(
            "SQL runs in polars-bio's session; registered tables are removed after the call"
        ),
    }


__all__ = [
    "describe_genomic_file",
    "query_genomic_tables",
    "read_alignment_sequences",
    "read_bigbed_intervals",
    "read_cram_records",
    "read_file_header",
    "write_alignment_file",
    "write_fastq_reads",
]

_register_public_verbose(globals(), tuple(__all__))
