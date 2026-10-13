"""Offline dnaio helpers for FASTA/FASTQ records, paired reads and read-name checks.

The functions read explicit local files or name lists through dnaio. Record
and pair counts are bounded, and malformed input is reported as a
ConfigurationError rather than a partial result.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.sequence.backends import _integer, _provenance, _source, _text

_DNAIO_EXTRA = "bio-dnaio"
_DEFAULT_MAX_RECORDS = 1_000
_MAX_RECORDS_LIMIT = 100_000
_MAX_NAMES = 10_000
_MAX_NAME = 1_024


def _read_dict(record: Any) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "name": str(record.name),
        "id": str(record.id),
        "comment": str(record.comment or ""),
        "sequence": str(record.sequence),
        "length": len(record.sequence),
    }
    if record.qualities is not None:
        entry["qualities"] = str(record.qualities)
    return entry


def read_fastx_records(
    source: str | Path,
    *,
    mate_source: str | Path | None = None,
    interleaved: bool = False,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read FASTA or FASTQ records, or paired reads, through dnaio.open with a record limit."""

    path = _source(source)
    if not isinstance(interleaved, bool):
        raise ConfigurationError("interleaved must be True or False.")
    if interleaved and mate_source is not None:
        raise ConfigurationError("interleaved=True cannot be combined with mate_source.")
    mate = None if mate_source is None else _source(mate_source)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    dnaio = require("dnaio", _DNAIO_EXTRA)

    if mate is not None:
        layout = "paired-end"
    elif interleaved:
        layout = "interleaved"
    else:
        layout = "single-end"
    try:
        if mate is not None:
            reader_context = dnaio.open(str(path), str(mate))
        elif interleaved:
            reader_context = dnaio.open(str(path), interleaved=True)
        else:
            reader_context = dnaio.open(str(path))
    except (OSError, ValueError, dnaio.FileFormatError) as exc:
        raise ConfigurationError(f"dnaio could not open {path.name}: {exc}") from exc

    entries: list[dict[str, Any]] = []
    try:
        with reader_context as reader:
            for item in reader:
                if len(entries) >= limit:
                    raise ConfigurationError(
                        f"{path.name} has more than max_records={limit} entries; "
                        "increase the explicit limit."
                    )
                if layout == "single-end":
                    entries.append(_read_dict(item))
                else:
                    first, second = item
                    entries.append({"read1": _read_dict(first), "read2": _read_dict(second)})
    except ConfigurationError:
        raise
    except (OSError, TypeError, ValueError, dnaio.FileFormatError) as exc:
        raise ConfigurationError(f"dnaio could not read {path.name}: {exc}") from exc

    fmt: str | None = None
    if entries:
        first_entry = entries[0] if layout == "single-end" else entries[0]["read1"]
        fmt = "fastq" if "qualities" in first_entry else "fasta"
    return {
        "status": "completed",
        "operation": "read_fastx_records",
        "source": str(path),
        "mate_source": None if mate is None else str(mate),
        "layout": layout,
        "format": fmt,
        "entry_unit": "record" if layout == "single-end" else "pair",
        "entry_count": len(entries),
        "records": entries,
        "provenance": _provenance(
            "dnaio",
            parser="dnaio.open",
            paired_check="dnaio raises when mate names disagree while reading pairs",
        ),
    }


def check_read_mates(names: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Check whether read names are mates using dnaio.records_are_mates on name-only records."""

    if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or len(names) < 2:
        raise ConfigurationError("names must be a list of at least two read names.")
    if len(names) > _MAX_NAMES:
        raise ConfigurationError(f"names has more than {_MAX_NAMES} entries.")
    values: list[str] = []
    for name in names:
        text = _text(name, "read name", maximum=_MAX_NAME)
        if "\n" in text or "\r" in text:
            raise ConfigurationError("read names cannot contain line breaks.")
        values.append(text)
    dnaio = require("dnaio", _DNAIO_EXTRA)
    try:
        # Only names are compared, so the sequence placeholder is empty and qualities absent.
        records = [dnaio.SequenceRecord(name, "", None) for name in values]
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"dnaio rejected a read name: {exc}") from exc
    return {
        "status": "completed",
        "operation": "check_read_mates",
        "read_count": len(values),
        "mates": bool(dnaio.records_are_mates(*records)),
        "names": values,
        "provenance": _provenance(
            "dnaio",
            rule="dnaio.records_are_mates compares names only",
        ),
    }


__all__ = [
    "check_read_mates",
    "read_fastx_records",
]

_register_public_verbose(globals(), tuple(__all__))
