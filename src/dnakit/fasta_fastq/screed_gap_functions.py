"""screed operations: build an indexed screed database from FASTA or FASTQ and dump it back to text.

screed stores each record in SQLite so that sequences can be retrieved by name. ``create_screed_database``
writes the database to an explicit path (screed's own command line writes next to the input file), and the
two dump functions write FASTA or FASTQ text to new files. Quality values are required for FASTQ output.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-screed"


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _database_record_count(database: Path) -> int:
    dict_table = require("screed.DBConstants", _EXTRA)._DICT_TABLE
    connection = sqlite3.connect(str(database))
    try:
        return int(connection.execute(f"SELECT COUNT(*) FROM {dict_table}").fetchone()[0])
    finally:
        connection.close()


def _first_record_fields(database: Path) -> set[str]:
    openscreed = require("screed.openscreed", _EXTRA)
    db = openscreed.ScreedDB(str(database))
    try:
        for value in db.itervalues():
            return set(value.keys())
        return set()
    finally:
        db.close()


def create_screed_database(
    sequence_file: str,
    output_database: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index a FASTA or FASTQ file into a screed database at an explicit path (screed createscreed.make_db).

    The format is detected from the first record. The database name gets screed's ``_screed`` suffix when it
    lacks it. The output must not exist. Returns the database path and the number of indexed records.
    """
    createscreed = require("screed.createscreed", _EXTRA)
    openscreed = require("screed.openscreed", _EXTRA)
    fasta = require("screed.fasta", _EXTRA)
    fastq = require("screed.fastq", _EXTRA)
    dbconstants = require("screed.DBConstants", _EXTRA)
    source = _existing(sequence_file, "sequence_file")
    target = _new(output_database, "output_database")
    if not target.name.endswith(dbconstants.fileExtension):
        target = target.with_name(target.name + dbconstants.fileExtension)
        if target.exists():
            raise ConfigurationError(f"output_database must not exist yet: {target}")
    iterator = openscreed.Open(str(source), parse_description=True)
    field_types = {
        fastq.fastq_iter.__name__: fastq.FieldTypes,
        fasta.fasta_iter.__name__: fasta.FieldTypes,
    }.get(iterator.iter_fn.__name__)
    if field_types is None:
        raise ConfigurationError(f"Unsupported sequence format for screed: {source}")
    createscreed.create_db(str(target), field_types, iterator)
    count = _database_record_count(target)
    return {
        "status": "completed",
        "backend": "screed",
        "input_format": "FASTQ" if field_types is fastq.FieldTypes else "FASTA",
        "database": str(target),
        "record_count": count,
    }


def _dump(database: str, output: str, *, fields_required: tuple[str, ...], writer: str) -> dict[str, Any]:
    conversion = require("screed.conversion", _EXTRA)
    source = _existing(database, "database")
    target = _new(output, "output")
    if _database_record_count(source) == 0:
        raise ConfigurationError(f"The screed database has no records: {source}")
    fields = _first_record_fields(source)
    missing = [field for field in fields_required if field not in fields]
    if missing:
        raise ConfigurationError(
            f"The screed database has no {', '.join(missing)} field, so it cannot be written as {writer}."
        )
    function = getattr(conversion, "ToFastq" if writer == "FASTQ" else "ToFasta")
    count = int(function(str(source), str(target)))
    return {
        "status": "completed",
        "backend": "screed",
        "database": str(source),
        "format": writer,
        "output": str(target),
        "record_count": count,
    }


def dump_screed_fasta(database: str, output_fasta: str, *, verbose: bool = True) -> dict[str, Any]:
    """Write every record of a screed database as FASTA text to a new file (screed conversion.ToFasta)."""
    return _dump(database, output_fasta, fields_required=("name", "sequence"), writer="FASTA")


def dump_screed_fastq(database: str, output_fastq: str, *, verbose: bool = True) -> dict[str, Any]:
    """Write every record of a screed database as FASTQ text to a new file (screed conversion.ToFastq).

    Needs a database built from FASTQ input, because the quality values are part of the record.
    """
    return _dump(database, output_fastq, fields_required=("name", "sequence", "quality"), writer="FASTQ")


__all__ = [
    "create_screed_database",
    "dump_screed_fasta",
    "dump_screed_fastq",
]
_register_public_verbose(globals())
