"""Biopython ``BioSQL`` operations added by integration pack P19.

BioSQL stores sequence records in a relational database (SQLite, MySQL or PostgreSQL) created from the
BioSQL schema (https://github.com/biosql/biosql). Records live in named sub-databases (``namespace``,
the ``biodatabase`` table). The functions connect with ``driver`` ``sqlite3`` (a database file, the
default), ``mysql.connector`` or ``psycopg2`` (``host``, ``user``, ``password``, ``port``); the driver
package for MySQL or PostgreSQL is installed by the user. A new SQLite database can be created from the
BioSQL schema file (``biosqldb-sqlite.sql``) with ``schema_sql`` in ``biosql_load_records``. The schema
file is not shipped with Biopython and is not redistributed by DNAKit. No taxonomy look-up on the NCBI
server is done while loading.
"""

from __future__ import annotations

import importlib.metadata
import re
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_DRIVERS = ("sqlite3", "mysql.connector", "psycopg2")
_ID_TYPES = ("display_id", "accession", "version", "primary_id", "name", "gi")
_FORMAT = re.compile(r"[a-z][a-z0-9\-]{1,30}")
_NAMESPACE = re.compile(r"[A-Za-z0-9_.\-]{1,100}")
_MAX_RECORDS = 1_000_000
_MAX_SEQUENCE = 100_000_000
_MAX_LISTED = 1_000_000


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _namespace(value: object) -> str:
    if not isinstance(value, str) or not _NAMESPACE.fullmatch(value):
        raise ConfigurationError("namespace must use letters, digits, '_', '.' or '-'.")
    return value


@contextmanager
def _server(
    driver: str,
    database: str | Path,
    host: str | None,
    user: str | None,
    password: str | None,
    port: int | None,
    *,
    schema_sql: str | Path | None = None,
) -> Iterator[Any]:
    """Open a BioSQL server; create the SQLite database from ``schema_sql`` when it does not exist."""

    if driver not in _DRIVERS:
        raise ConfigurationError(f"driver must be one of {list(_DRIVERS)}.")
    biosql = require("BioSQL.BioSeqDatabase", _EXTRA)
    if driver == "sqlite3":
        path = Path(database).expanduser().resolve()
        creating = schema_sql is not None
        if creating and path.exists():
            raise ConfigurationError(f"database already exists, schema_sql would not be applied: {path}")
        if not creating and not path.is_file():
            raise ConfigurationError(f"BioSQL SQLite database not found: {path}")
        if creating and not path.parent.is_dir():
            raise ConfigurationError(f"The folder for the database does not exist: {path.parent}")
        server = biosql.open_database(driver="sqlite3", db=str(path))
        if creating:
            try:
                server.load_database_sql(str(cli.existing_file(schema_sql, "schema_sql")))
                server.commit()
            except Exception:
                server.close()
                path.unlink(missing_ok=True)
                raise
    else:
        if schema_sql is not None:
            raise ConfigurationError("schema_sql is only supported for SQLite databases.")
        if not isinstance(database, str) or not database:
            raise ConfigurationError("database must be the name of the server database.")
        settings: dict[str, Any] = {"db": database}
        for key, value in (("host", host), ("user", user), ("password", password), ("port", port)):
            if value is not None:
                settings[key] = value
        try:
            server = biosql.open_database(driver=driver, **settings)
        except ImportError as exc:
            raise ConfigurationError(f"The database driver {driver} is not installed: {exc}") from exc
    try:
        yield server
    finally:
        server.close()


def biosql_load_records(
    database: str | Path,
    namespace: str,
    records_path: str | Path,
    *,
    record_format: str = "genbank",
    driver: str = "sqlite3",
    schema_sql: str | Path | None = None,
    host: str | None = None,
    user: str | None = None,
    password: str | None = None,
    port: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load sequence records from a file into a BioSQL sub-database.

    Reads ``records_path`` with ``Bio.SeqIO`` (``record_format``: genbank, embl, fasta, swiss, ...), creates the
    sub-database ``namespace`` when it is new (``BioSeqDatabase.new_database``) and loads the records with
    ``BioSeqDatabase.load`` (``BioSQL.Loader.DatabaseLoader.load_seqrecord``), then commits. For SQLite,
    pass ``schema_sql`` (path of ``biosqldb-sqlite.sql``) to create a new database file; an existing database
    is never recreated. Returns the number of records loaded.
    """

    space = _namespace(namespace)
    source = input_path(records_path)
    if not isinstance(record_format, str) or not _FORMAT.fullmatch(record_format):
        raise ConfigurationError("record_format must be a Bio.SeqIO format name.")
    seqio = require("Bio.SeqIO", _EXTRA)
    with _server(driver, database, host, user, password, port, schema_sql=schema_sql) as server:
        existed = space in server
        db = server[space] if existed else server.new_database(space)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with source.open(encoding="utf-8") as handle:
                    count = db.load(seqio.parse(handle, record_format), fetch_NCBI_taxonomy=False)
            server.commit()
        except Exception as exc:
            server.rollback()
            raise BackendExecutionError(f"BioSQL could not load the records: {exc}") from exc
        total = len(db)
    return {
        "status": "completed",
        "operation": "biosql_load_records",
        "database": str(database),
        "namespace": space,
        "created_namespace": not existed,
        "records_loaded": int(count),
        "records_in_namespace": int(total),
        "provenance": _provenance(parser="BioSQL.BioSeqDatabase.BioSeqDatabase.load; BioSQL.Loader.DatabaseLoader"),
    }


def _summary(record: Any, include_sequence: bool, max_sequence_length: int) -> dict[str, Any]:
    annotations = {}
    for key, value in record.annotations.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            annotations[key] = value
        elif isinstance(value, (list, tuple)) and all(isinstance(item, (str, int, float)) for item in value):
            annotations[key] = list(value)
    text = None
    length = len(record.seq)
    if include_sequence:
        if length > max_sequence_length:
            raise ConfigurationError(
                f"the sequence has {length} bases; raise max_sequence_length or set include_sequence=False."
            )
        text = str(record.seq)
    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "length": length,
        "sequence": text,
        "feature_count": len(record.features),
        "feature_types": sorted({feature.type for feature in record.features}),
        "annotations": annotations,
        "dbxrefs": list(record.dbxrefs),
    }


def biosql_get_record(
    database: str | Path,
    namespace: str,
    identifier: str,
    *,
    id_type: str = "display_id",
    include_sequence: bool = True,
    max_sequence_length: int = 1_000_000,
    output_path: str | Path | None = None,
    output_format: str = "genbank",
    driver: str = "sqlite3",
    host: str | None = None,
    user: str | None = None,
    password: str | None = None,
    port: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve one record from a BioSQL sub-database by identifier.

    ``id_type`` selects the identifier kind: ``display_id`` (``get_Seq_by_id``), ``accession``
    (``get_Seq_by_acc``), ``version`` (``get_Seq_by_ver``), ``primary_id``, ``name`` or ``gi``
    (``BioSeqDatabase.lookup``). Returns the record id, name, description, length, sequence, feature
    summary, annotations and cross-references; with ``output_path`` (a new file) the record is also
    written in ``output_format``.
    """

    space = _namespace(namespace)
    if id_type not in _ID_TYPES:
        raise ConfigurationError(f"id_type must be one of {list(_ID_TYPES)}.")
    if not isinstance(identifier, str) or not identifier:
        raise ConfigurationError("identifier must be a non-empty string.")
    if not isinstance(include_sequence, bool):
        raise ConfigurationError("include_sequence must be True or False.")
    limit = max_sequence_length
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= _MAX_SEQUENCE:
        raise ConfigurationError(f"max_sequence_length must be an integer in [1, {_MAX_SEQUENCE}].")
    target = cli.new_file(output_path, "output_path") if output_path is not None else None
    if target is not None and (not isinstance(output_format, str) or not _FORMAT.fullmatch(output_format)):
        raise ConfigurationError("output_format must be a Bio.SeqIO format name.")
    with _server(driver, database, host, user, password, port) as server:
        if space not in server:
            raise ConfigurationError(f"namespace {space!r} is not in the database.")
        db = server[space]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                record = {
                    "display_id": lambda: db.get_Seq_by_id(identifier),
                    "accession": lambda: db.get_Seq_by_acc(identifier),
                    "version": lambda: db.get_Seq_by_ver(identifier),
                }.get(id_type, lambda: db.lookup(**{id_type: identifier}))()
        except (KeyError, IndexError) as exc:
            raise ConfigurationError(f"no record with {id_type} {identifier!r} in {space!r}.") from exc
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary = _summary(record, include_sequence, limit)
            if target is not None:
                seqio = require("Bio.SeqIO", _EXTRA)
                try:
                    seqio.write(record, str(target), output_format)
                except ValueError as exc:
                    target.unlink(missing_ok=True)
                    raise ConfigurationError(f"the record cannot be written as {output_format}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "biosql_get_record",
        "namespace": space,
        "id_type": id_type,
        "record": summary,
        "output": str(target) if target is not None else None,
        "provenance": _provenance(
            parser="BioSQL.BioSeqDatabase.BioSeqDatabase.get_Seq_by_id / get_Seq_by_acc / get_Seq_by_ver / lookup"
        ),
    }


def biosql_get_subsequence(
    database: str | Path,
    namespace: str,
    identifier: str,
    start: int,
    end: int,
    *,
    id_type: str = "display_id",
    driver: str = "sqlite3",
    host: str | None = None,
    user: str | None = None,
    password: str | None = None,
    port: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a slice of a stored sequence without loading the whole record.

    The slice ``[start, end)`` (0-based, half-open) is cut inside the database with SQL ``SUBSTR``
    (``BioSQL.BioSeqDatabase.Adaptor.get_subseq_as_string``). ``identifier`` and ``id_type`` select the
    record as in ``biosql_get_record``.
    """

    space = _namespace(namespace)
    if id_type not in _ID_TYPES:
        raise ConfigurationError(f"id_type must be one of {list(_ID_TYPES)}.")
    if not isinstance(identifier, str) or not identifier:
        raise ConfigurationError("identifier must be a non-empty string.")
    for value, name in ((start, "start"), (end, "end")):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ConfigurationError(f"{name} must be a non-negative integer.")
    if end <= start:
        raise ConfigurationError("end must be greater than start.")
    if end - start > _MAX_SEQUENCE:
        raise ConfigurationError(f"the slice may be at most {_MAX_SEQUENCE} bases.")
    with _server(driver, database, host, user, password, port) as server:
        if space not in server:
            raise ConfigurationError(f"namespace {space!r} is not in the database.")
        db = server[space]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                record = db.lookup(**{id_type: identifier})
        except (KeyError, IndexError) as exc:
            raise ConfigurationError(f"no record with {id_type} {identifier!r} in {space!r}.") from exc
        seqid = record._primary_id  # noqa: SLF001 - internal bioentry id of the DBSeqRecord
        length = len(record.seq)
        if end > length:
            raise ConfigurationError(f"end {end} is beyond the sequence length {length}.")
        text = server.adaptor.get_subseq_as_string(seqid, start, end)
    return {
        "status": "completed",
        "operation": "biosql_get_subsequence",
        "namespace": space,
        "identifier": identifier,
        "start": start,
        "end": end,
        "sequence_length": length,
        "sequence": text,
        "provenance": _provenance(parser="BioSQL.BioSeqDatabase.Adaptor.get_subseq_as_string"),
    }


def biosql_list_entries(
    database: str | Path,
    *,
    namespace: str | None = None,
    max_entries: int = 10_000,
    driver: str = "sqlite3",
    host: str | None = None,
    user: str | None = None,
    password: str | None = None,
    port: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the sub-databases of a BioSQL database, or the entries of one sub-database.

    Without ``namespace`` returns the sub-database names with their entry counts
    (``Adaptor.list_biodatabase_names``). With ``namespace`` returns the internal entry ids and the display
    ids (``Adaptor.list_bioentry_ids``, ``list_bioentry_display_ids``), at most ``max_entries``.
    """

    if isinstance(max_entries, bool) or not isinstance(max_entries, int) or not 1 <= max_entries <= _MAX_LISTED:
        raise ConfigurationError(f"max_entries must be an integer in [1, {_MAX_LISTED}].")
    with _server(driver, database, host, user, password, port) as server:
        adaptor = server.adaptor
        if namespace is None:
            names = list(adaptor.list_biodatabase_names())
            return {
                "status": "completed",
                "operation": "biosql_list_entries",
                "namespaces": [
                    {"name": name, "entry_count": len(adaptor.list_bioentry_ids(adaptor.fetch_dbid_by_dbname(name)))}
                    for name in names
                ],
                "provenance": _provenance(parser="BioSQL.BioSeqDatabase.Adaptor.list_biodatabase_names"),
            }
        space = _namespace(namespace)
        if space not in server:
            raise ConfigurationError(f"namespace {space!r} is not in the database.")
        dbid = adaptor.fetch_dbid_by_dbname(space)
        ids = list(adaptor.list_bioentry_ids(dbid))
        display = list(adaptor.list_bioentry_display_ids(dbid))
    return {
        "status": "completed",
        "operation": "biosql_list_entries",
        "namespace": space,
        "entry_count": len(ids),
        "truncated": len(ids) > max_entries,
        "entry_ids": ids[:max_entries],
        "display_ids": display[:max_entries],
        "provenance": _provenance(
            parser="BioSQL.BioSeqDatabase.Adaptor.list_bioentry_ids / list_bioentry_display_ids"
        ),
    }


__all__ = [
    "biosql_get_record",
    "biosql_get_subsequence",
    "biosql_list_entries",
    "biosql_load_records",
]

_register_public_verbose(globals(), tuple(__all__))
