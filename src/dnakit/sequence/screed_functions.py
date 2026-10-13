"""Offline screed helpers for exporting databases, parsing FASTX files and formatting records.

screed databases are SQLite files written by ``build_screed_database``. The
functions here read them back out, parse plain FASTA/FASTQ files with screed's
own iterators, or format explicit records with screed's record writer. Every
path is local; nothing is downloaded.
"""

from __future__ import annotations

import io
import sqlite3
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.sequence.backends import (
    _close,
    _flag,
    _integer,
    _member,
    _new_output,
    _provenance,
    _screed_target,
    _source,
    _text,
)

_SCREED_EXTRA = "bio-screed"
_FASTX_FORMATS = ("fasta", "fastq")
_DEFAULT_MAX_RECORDS = 1_000
_MAX_RECORDS_LIMIT = 100_000
_MAX_NAME = 512
_MAX_SEQUENCE = 1_000_000
_MAX_DESCRIPTION = 4_096


def _record_count(screed: Any, target: Path) -> int:
    try:
        handle = screed.ScreedDB(str(target))
    except (OSError, ValueError, sqlite3.Error) as exc:
        raise ConfigurationError(f"screed could not open {target.name}: {exc}") from exc
    try:
        return len(list(handle.keys()))
    except (OSError, TypeError, ValueError, sqlite3.Error) as exc:
        raise ConfigurationError(f"screed could not read {target.name}: {exc}") from exc
    finally:
        _close(handle)


def export_screed_database(
    database: str | Path,
    output: str | Path,
    *,
    output_format: str = "fasta",
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write every record of a screed database to a new FASTA or FASTQ file through screed."""

    target = _screed_target(database)
    if not target.is_file():
        raise ConfigurationError(f"screed database does not exist: {target}")
    fmt = _member(output_format, _FASTX_FORMATS, "output_format")
    destination = _new_output(output, "output")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    screed = require("screed", _SCREED_EXTRA)
    count = _record_count(screed, target)
    if count == 0:
        raise ConfigurationError(f"screed database {target.name} holds no records to export.")
    if count > limit:
        raise ConfigurationError(
            f"{target.name} holds {count} records, above max_records={limit}; "
            "increase the explicit limit."
        )
    exporter = screed.ToFasta if fmt == "fasta" else screed.ToFastq
    try:
        written = exporter(str(target), str(destination))
    except (OSError, TypeError, ValueError, sqlite3.Error) as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(f"screed could not export {target.name}: {exc}") from exc
    if written != count:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(
            f"screed wrote {written} records but the database holds {count}; output removed."
        )
    return {
        "status": "completed",
        "operation": "export_screed_database",
        "database": str(target),
        "output": str(destination),
        "output_format": fmt,
        "record_count": written,
        "provenance": _provenance(
            "screed",
            writer="screed.ToFasta" if fmt == "fasta" else "screed.ToFastq",
            line_wrap="sequence lines wrapped at 80 columns by screed",
            missing_quality="FASTQ export fills absent qualities with screed's placeholder",
        ),
    }


def _screed_entry(raw: Any, file_format: str) -> dict[str, Any]:
    sequence = str(raw.sequence)
    description = getattr(raw, "description", None)
    if description is None:
        # screed's FASTQ parser stores the header text under 'annotations'.
        description = getattr(raw, "annotations", "")
    entry: dict[str, Any] = {
        "name": str(raw.name),
        "description": str(description or ""),
        "length": len(sequence),
        "sequence": sequence,
    }
    if file_format == "fastq":
        entry["quality"] = str(raw.quality)
    return entry


def read_screed_fastx_records(
    source: str | Path,
    *,
    file_format: str = "fasta",
    parse_description: bool = True,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse an uncompressed FASTA or FASTQ file with screed's record iterator."""

    path = _source(source)
    fmt = _member(file_format, _FASTX_FORMATS, "file_format")
    describe = _flag(parse_description, "parse_description")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    screed = require("screed", _SCREED_EXTRA)
    iterator = screed.fasta.fasta_iter if fmt == "fasta" else screed.fastq.fastq_iter
    records: list[dict[str, Any]] = []
    try:
        with path.open("rb") as handle:
            for raw in iterator(handle, parse_description=describe):
                if len(records) >= limit:
                    raise ConfigurationError(
                        f"{path.name} has more than max_records={limit} records; "
                        "increase the explicit limit."
                    )
                records.append(_screed_entry(raw, fmt))
    except ConfigurationError:
        raise
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"screed could not parse {path.name} as {fmt}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_screed_fastx_records",
        "source": str(path),
        "format": fmt,
        "record_count": len(records),
        "records": records,
        "provenance": _provenance(
            "screed",
            parser="screed.fasta.fasta_iter" if fmt == "fasta" else "screed.fastq.fastq_iter",
            compression="uncompressed files only",
            description_rule="FASTQ header text (screed annotations) is reported as description",
        ),
    }


def _token(value: object, name: str, *, maximum: int) -> str:
    text = _text(value, name, maximum=maximum)
    if any(character.isspace() for character in text):
        raise ConfigurationError(f"{name} cannot contain whitespace.")
    return text


def _record_for_write(screed: Any, item: object, file_format: str, index: int) -> Any:
    if not isinstance(item, Mapping):
        raise ConfigurationError(f"records[{index}] must be a mapping with name and sequence.")
    name = _token(item.get("name"), f"records[{index}].name", maximum=_MAX_NAME)
    sequence = _token(item.get("sequence"), f"records[{index}].sequence", maximum=_MAX_SEQUENCE)
    if not all(
        character.isascii() and (character.isalpha() or character in "-*") for character in sequence
    ):
        raise ConfigurationError(
            f"records[{index}].sequence may contain only ASCII letters, '-' and '*'."
        )
    fields: dict[str, Any] = {"name": name, "sequence": sequence}
    description = item.get("description", "")
    if not isinstance(description, str) or "\n" in description or "\r" in description:
        raise ConfigurationError(f"records[{index}].description must be one line of text.")
    if len(description) > _MAX_DESCRIPTION:
        raise ConfigurationError(f"records[{index}].description is too long.")
    if description:
        fields["description"] = description
    quality = item.get("quality")
    if file_format == "fastq":
        if (
            not isinstance(quality, str)
            or len(quality) != len(sequence)
            or not all(33 <= ord(character) <= 126 for character in quality)
        ):
            raise ConfigurationError(
                f"records[{index}].quality must be printable ASCII with one symbol per base."
            )
        fields["quality"] = quality
    elif quality is not None:
        raise ConfigurationError("quality is written only when output_format is 'fastq'.")
    return screed.Record(**fields)


def format_fastx_records(
    records: Sequence[Mapping[str, Any]],
    *,
    output_format: str = "fasta",
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format explicit name, description, sequence and quality records as FASTA or FASTQ text."""

    fmt = _member(output_format, _FASTX_FORMATS, "output_format")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    if isinstance(records, (str, bytes)) or not isinstance(records, Sequence) or not records:
        raise ConfigurationError("records must be a non-empty list of name/sequence mappings.")
    if len(records) > limit:
        raise ConfigurationError(f"records has more than max_records={limit} entries.")
    screed = require("screed", _SCREED_EXTRA)
    buffer = io.BytesIO()
    for index, item in enumerate(records):
        record = _record_for_write(screed, item, fmt, index)
        screed.screedRecord.write_fastx(record, buffer)
    return {
        "status": "completed",
        "operation": "format_fastx_records",
        "output_format": fmt,
        "record_count": len(records),
        "text": buffer.getvalue().decode("utf-8"),
        "provenance": _provenance(
            "screed",
            writer="screed.screedRecord.write_fastx",
            layout="each sequence is written on one line; no line wrapping is applied",
        ),
    }


__all__ = [
    "export_screed_database",
    "format_fastx_records",
    "read_screed_fastx_records",
]

_register_public_verbose(globals(), tuple(__all__))
