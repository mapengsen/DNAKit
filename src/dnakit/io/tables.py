"""Strict, bounded, atomic export for result objects and two-dimensional tables."""

from __future__ import annotations

import csv
import io
import json
import math
import os
import platform
import re
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Any, BinaryIO, Literal, Protocol, TextIO, cast, runtime_checkable

from dnakit.core._json import FrozenDict, JSONScalar, freeze_mapping, to_json_compatible
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    DNAKitError,
    InputFormatError,
)

from ._advanced_common import write_text_path
from .parquet import ParquetBackendStatus, require_parquet_backend

TableFormat = Literal["csv", "tsv", "json", "parquet"]
ParquetCompression = Literal["none", "snappy", "gzip", "brotli", "lz4", "zstd"]
TableScalar = JSONScalar
TableColumnType = Literal["any", "string", "integer", "number", "boolean"]
_MAX_TABLE_SCHEMA_COLUMNS = 10_000


@runtime_checkable
class DictResult(Protocol):
    """Structural protocol for serializable DNAKit result objects."""

    def to_dict(self) -> dict[str, Any]: ...


def _validate_columns(columns: Iterable[str]) -> tuple[str, ...]:
    if isinstance(columns, (str, bytes)):
        raise ConfigurationError(
            "Table columns must be an iterable of column names.",
            code="INVALID_TABLE_SCHEMA",
        )
    try:
        resolved = tuple(islice(iter(columns), _MAX_TABLE_SCHEMA_COLUMNS + 1))
    except TypeError as exc:
        raise ConfigurationError(
            "Table columns must be an iterable of column names.",
            code="INVALID_TABLE_SCHEMA",
        ) from exc
    if len(resolved) > _MAX_TABLE_SCHEMA_COLUMNS:
        raise ConfigurationError(
            "Table schema exceeds the hard column limit.",
            code="INVALID_TABLE_SCHEMA",
            context={"max_columns": _MAX_TABLE_SCHEMA_COLUMNS},
        )
    if not resolved or any(
        not isinstance(column, str)
        or not column
        or "\n" in column
        or "\r" in column
        or "\t" in column
        for column in resolved
    ):
        raise ConfigurationError(
            "Table columns must be non-empty single-line strings without tabs.",
            code="INVALID_TABLE_SCHEMA",
        )
    if len(set(resolved)) != len(resolved):
        raise ConfigurationError("Table columns must be unique.", code="INVALID_TABLE_SCHEMA")
    return resolved


def _validate_scalar(value: object, *, row_index: int, column: str) -> TableScalar:
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float) and math.isfinite(value):
        return value
    raise InputFormatError(
        "Table cells must be JSON scalar values.",
        code="INVALID_TABLE_CELL",
        context={
            "row_index": row_index,
            "column": column,
            "value_type": type(value).__name__,
        },
    )


def _validate_schema_scalar(
    value: object, *, schema: TableSchema, row_index: int, column: str
) -> TableScalar:
    scalar = _validate_scalar(value, row_index=row_index, column=column)
    if scalar is None:
        if column not in cast(tuple[str, ...], schema.nullable):
            raise InputFormatError(
                "A non-nullable table column contains a missing value.",
                code="TABLE_NULL_NOT_ALLOWED",
                context={"row_index": row_index, "column": column},
            )
        return None
    expected = schema.column_type(column)
    valid = (
        expected == "any"
        or (expected == "string" and isinstance(scalar, str))
        or (expected == "integer" and isinstance(scalar, int) and not isinstance(scalar, bool))
        or (
            expected == "number"
            and isinstance(scalar, (int, float))
            and not isinstance(scalar, bool)
        )
        or (expected == "boolean" and isinstance(scalar, bool))
    )
    if not valid:
        raise InputFormatError(
            "Table cell does not match its declared column type.",
            code="TABLE_COLUMN_TYPE_MISMATCH",
            context={
                "row_index": row_index,
                "column": column,
                "expected": expected,
                "actual": type(scalar).__name__,
            },
        )
    return scalar


@dataclass(frozen=True, slots=True)
class TableSchema:
    """Explicit stable column order and optional semantic schema version."""

    columns: tuple[str, ...]
    schema_version: str = "dnakit.table.v1"
    column_types: Mapping[str, TableColumnType] = field(
        default_factory=lambda: cast(Mapping[str, TableColumnType], FrozenDict())
    )
    nullable: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "columns", _validate_columns(self.columns))
        if not isinstance(self.schema_version, str) or not self.schema_version.strip():
            raise ConfigurationError(
                "Table schema_version must be non-empty.", code="INVALID_TABLE_SCHEMA"
            )
        if not isinstance(self.column_types, Mapping):
            raise ConfigurationError(
                "Table column_types must be a mapping.", code="INVALID_TABLE_SCHEMA"
            )
        raw_types = dict(self.column_types)
        if any(column not in self.columns for column in raw_types) or any(
            value not in ("any", "string", "integer", "number", "boolean")
            for value in raw_types.values()
        ):
            raise ConfigurationError(
                "Table column_types keys and values must match the declared schema.",
                code="INVALID_TABLE_SCHEMA",
            )
        object.__setattr__(self, "column_types", freeze_mapping(raw_types))
        nullable = self.columns if self.nullable is None else tuple(self.nullable)
        if (
            any(not isinstance(column, str) for column in nullable)
            or len(set(nullable)) != len(nullable)
            or any(column not in self.columns for column in nullable)
        ):
            raise ConfigurationError(
                "Table nullable columns must be unique declared columns.",
                code="INVALID_TABLE_SCHEMA",
            )
        object.__setattr__(self, "nullable", nullable)

    def column_type(self, column: str) -> TableColumnType:
        """Return one column's strict scalar type, defaulting to ``any``."""

        return self.column_types.get(column, "any")


@dataclass(frozen=True, slots=True)
class TableExportResult:
    format: TableFormat
    row_count: int
    columns: tuple[str, ...]
    schema_version: str
    target_path: str
    byte_count: int
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if self.format not in {"csv", "tsv", "json", "parquet"}:
            raise ConfigurationError(
                "Table export result has an invalid format.", code="INVALID_TABLE_RESULT"
            )
        _validate_columns(self.columns)
        if not isinstance(self.schema_version, str) or not self.schema_version:
            raise ConfigurationError(
                "Table export schema_version must be non-empty.",
                code="INVALID_TABLE_RESULT",
            )
        if not isinstance(self.target_path, str) or not self.target_path:
            raise ConfigurationError(
                "Table export target_path must be non-empty.", code="INVALID_TABLE_RESULT"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in (self.row_count, self.byte_count)
        ):
            raise ConfigurationError(
                "Table export counts must be non-negative integers.",
                code="INVALID_TABLE_RESULT",
            )
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError(
                "Table export parameters must be an immutable mapping.",
                code="INVALID_TABLE_RESULT",
            )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class TableReadResult:
    """Immutable rows and audit metadata returned by :func:`read_table`."""

    format: TableFormat
    rows: tuple[FrozenDict, ...]
    schema: TableSchema
    source_path: str | None
    byte_count: int
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if self.format not in {"csv", "tsv", "json", "parquet"}:
            raise ConfigurationError(
                "Table read result has an invalid format.", code="INVALID_TABLE_RESULT"
            )
        if not isinstance(self.rows, tuple) or any(
            not isinstance(row, FrozenDict) for row in self.rows
        ):
            raise ConfigurationError(
                "Table read rows must be an immutable tuple of FrozenDict values.",
                code="INVALID_TABLE_RESULT",
            )
        if not isinstance(self.schema, TableSchema):
            raise ConfigurationError(
                "Table read schema must be a TableSchema.", code="INVALID_TABLE_RESULT"
            )
        if any(tuple(row) != self.schema.columns for row in self.rows):
            raise ConfigurationError(
                "Table read row keys and order must match its schema.",
                code="INVALID_TABLE_RESULT",
            )
        if self.source_path is not None and (
            not isinstance(self.source_path, str) or not self.source_path
        ):
            raise ConfigurationError(
                "Table read source_path must be non-empty text or None.",
                code="INVALID_TABLE_RESULT",
            )
        if (
            isinstance(self.byte_count, bool)
            or not isinstance(self.byte_count, int)
            or self.byte_count < 0
        ):
            raise ConfigurationError(
                "Table read byte_count must be a non-negative integer.",
                code="INVALID_TABLE_RESULT",
            )
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError(
                "Table read parameters must be an immutable mapping.",
                code="INVALID_TABLE_RESULT",
            )

    @property
    def row_count(self) -> int:
        return len(self.rows)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def _schema_payload(schema: TableSchema) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": schema.schema_version,
        "columns": list(schema.columns),
    }
    if schema.column_types:
        payload["column_types"] = dict(schema.column_types)
    if schema.nullable != schema.columns:
        payload["nullable"] = list(cast(tuple[str, ...], schema.nullable))
    return payload


def _native_table_backend() -> dict[str, object]:
    return {
        "name": "python-stdlib",
        "version": platform.python_version(),
        "implementation": "native",
        "license": "PSF-2.0",
    }


def _parquet_backend(status: ParquetBackendStatus, *, compression: object) -> dict[str, object]:
    return {
        "name": status.distribution,
        "version": status.version,
        "implementation": status.implementation,
        "license": "Apache-2.0",
        "compression": compression,
    }


class _BoundedTableTextWriter:
    """Count UTF-8 bytes before forwarding table text writes."""

    def __init__(self, raw: TextIO, *, limit: int) -> None:
        self._raw = raw
        self._limit = limit
        self._count = 0

    def write(self, value: str) -> int:
        size = len(value.encode("utf-8"))
        if self._count + size > self._limit:
            raise InputFormatError(
                "Table export exceeds max_output_bytes.",
                code="TABLE_OUTPUT_LIMIT_EXCEEDED",
                context={"max_output_bytes": self._limit},
            )
        written = self._raw.write(value)
        self._count += size
        return written

    def flush(self) -> None:
        self._raw.flush()

    def __getattr__(self, name: str) -> object:
        return getattr(self._raw, name)


class _BoundedTableBinaryWriter:
    """Count physical bytes before forwarding Parquet writes."""

    def __init__(self, raw: BinaryIO, *, limit: int) -> None:
        self._raw = raw
        self._limit = limit
        self._count = 0

    def write(self, value: bytes | bytearray | memoryview) -> int:
        size = len(value)
        if self._count + size > self._limit:
            raise InputFormatError(
                "Table export exceeds max_output_bytes.",
                code="TABLE_OUTPUT_LIMIT_EXCEEDED",
                context={"max_output_bytes": self._limit},
            )
        written = self._raw.write(value)
        self._count += written
        return written

    def flush(self) -> None:
        self._raw.flush()

    def __getattr__(self, name: str) -> object:
        return getattr(self._raw, name)


def _iter_rows(
    rows: Iterable[Mapping[str, object] | Sequence[object] | DictResult],
    *,
    schema: TableSchema,
    max_rows: int,
    max_cell_characters: int,
) -> Iterable[dict[str, TableScalar]]:
    for row_index, row in enumerate(rows):
        if row_index >= max_rows:
            raise InputFormatError(
                "Table export exceeds max_rows.",
                code="TABLE_ROW_LIMIT_EXCEEDED",
                context={"max_rows": max_rows},
            )
        if isinstance(row, DictResult):
            raw: object = row.to_dict()
        else:
            raw = row
        width = len(schema.columns)
        if isinstance(raw, Mapping):
            keys = tuple(islice(iter(raw), width + 1))
            if (
                len(keys) != width
                or any(not isinstance(key, str) for key in keys)
                or set(keys) != set(schema.columns)
            ):
                raise InputFormatError(
                    "Mapping row keys must match the table schema exactly.",
                    code="TABLE_SCHEMA_MISMATCH",
                    context={"row_index": row_index},
                )
            values = tuple(raw[column] for column in schema.columns)
        elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
            values = tuple(islice(iter(raw), width + 1))
            if len(values) != width:
                raise InputFormatError(
                    "Sequence row length must equal the table column count.",
                    code="TABLE_SCHEMA_MISMATCH",
                    context={"row_index": row_index},
                )
        else:
            raise TypeError(
                "Each table row must be a mapping, non-text sequence, or object with to_dict()."
            )
        resolved: dict[str, TableScalar] = {}
        for column, value in zip(schema.columns, values, strict=True):
            scalar = _validate_schema_scalar(
                value, schema=schema, row_index=row_index, column=column
            )
            if isinstance(scalar, str) and len(scalar) > max_cell_characters:
                raise InputFormatError(
                    "Table cell exceeds max_cell_characters.",
                    code="TABLE_CELL_LIMIT_EXCEEDED",
                    context={
                        "row_index": row_index,
                        "column": column,
                        "max_cell_characters": max_cell_characters,
                    },
                )
            resolved[column] = scalar
        yield resolved


def _write_delimited(
    handle: TextIO,
    rows: Iterable[dict[str, TableScalar]],
    *,
    columns: tuple[str, ...],
    delimiter: str,
    null_value: str,
) -> int:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    count = 0
    for row in rows:
        values: dict[str, TableScalar] = {}
        for column, value in row.items():
            if isinstance(value, str) and value == null_value:
                raise InputFormatError(
                    "A string table cell conflicts with null_value.",
                    code="TABLE_NULL_SENTINEL_COLLISION",
                    context={"row_index": count, "column": column, "null_value": null_value},
                )
            values[column] = null_value if value is None else value
        writer.writerow(values)
        count += 1
    return count


def _write_json(
    handle: TextIO,
    rows: Iterable[dict[str, TableScalar]],
    *,
    schema: TableSchema,
) -> int:
    handle.write(
        json.dumps(
            _schema_payload(schema),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )[:-1]
        + ',"rows":['
    )
    count = 0
    for row in rows:
        if count:
            handle.write(",")
        handle.write(json.dumps(row, ensure_ascii=False, sort_keys=False, separators=(",", ":")))
        count += 1
    handle.write("]}\n")
    return count


def _write_parquet(
    rows: Iterable[dict[str, TableScalar]],
    path: Path,
    *,
    schema: TableSchema,
    overwrite: bool,
    max_output_bytes: int,
    compression: ParquetCompression,
) -> int:
    require_parquet_backend()
    try:
        import pyarrow as pa  # type: ignore[import-untyped]
        import pyarrow.parquet as pq  # type: ignore[import-untyped]
    except (ImportError, AttributeError) as exc:
        raise BackendExecutionError(
            "PyArrow was detected but could not be imported.",
            code="PARQUET_BACKEND_IMPORT_FAILED",
        ) from exc
    materialized = list(rows)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    if not path.parent.exists():
        raise FileNotFoundError(f"Output parent directory does not exist: {path.parent}")
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(descriptor)
    temporary = Path(name)
    try:
        columns = {column: [row[column] for row in materialized] for column in schema.columns}
        table = pa.table(columns)
        metadata = dict(table.schema.metadata or {})
        metadata[b"dnakit_schema_version"] = schema.schema_version.encode("utf-8")
        metadata[b"dnakit_table_schema"] = json.dumps(
            _schema_payload(schema),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        table = table.replace_schema_metadata(metadata)
        with temporary.open("wb") as raw:
            bounded = _BoundedTableBinaryWriter(raw, limit=max_output_bytes)
            pq.write_table(
                table,
                bounded,
                compression=None if compression == "none" else compression,
            )
            bounded.flush()
            os.fsync(raw.fileno())
    except InputFormatError:
        temporary.unlink(missing_ok=True)
        raise
    except Exception as exc:
        temporary.unlink(missing_ok=True)
        raise BackendExecutionError(
            "PyArrow failed to serialize the bounded table.",
            code="PARQUET_BACKEND_EXECUTION_FAILED",
            context={"reason": str(exc)},
        ) from exc
    try:
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
            temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return len(materialized)


def export_table(
    rows: Iterable[Mapping[str, object] | Sequence[object] | DictResult],
    target: str | os.PathLike[str],
    *,
    format: TableFormat,
    schema: TableSchema,
    overwrite: bool = False,
    max_rows: int = 1_000_000,
    max_columns: int = 10_000,
    max_cell_characters: int = 10_000_000,
    max_output_bytes: int = 1_000_000_000,
    null_value: str = r"\N",
    parquet_compression: ParquetCompression = "snappy",
) -> TableExportResult:
    """Export one bounded table without inferring or silently widening its schema."""

    if format not in {"csv", "tsv", "json", "parquet"}:
        raise ConfigurationError("Unknown table export format.", code="INVALID_TABLE_FORMAT")
    if not isinstance(schema, TableSchema):
        raise TypeError("schema must be a TableSchema.")
    for name, value in (
        ("max_rows", max_rows),
        ("max_columns", max_columns),
        ("max_cell_characters", max_cell_characters),
        ("max_output_bytes", max_output_bytes),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ConfigurationError(
                f"{name} must be positive.",
                code="INVALID_TABLE_LIMIT",
                context={"field": name},
            )
    if not isinstance(null_value, str) or not null_value:
        raise ConfigurationError(
            "null_value must be non-empty text.", code="INVALID_TABLE_NULL_VALUE"
        )
    if len(null_value) > max_cell_characters:
        raise ConfigurationError(
            "null_value exceeds max_cell_characters.", code="INVALID_TABLE_NULL_VALUE"
        )
    if not isinstance(parquet_compression, str) or parquet_compression not in (
        "none",
        "snappy",
        "gzip",
        "brotli",
        "lz4",
        "zstd",
    ):
        raise ConfigurationError(
            "parquet_compression is unsupported.", code="INVALID_PARQUET_COMPRESSION"
        )
    if len(schema.columns) > max_columns:
        raise InputFormatError(
            "Table schema exceeds max_columns.",
            code="TABLE_COLUMN_LIMIT_EXCEEDED",
            context={"max_columns": max_columns},
        )
    path = Path(target)
    prepared = _iter_rows(
        rows,
        schema=schema,
        max_rows=max_rows,
        max_cell_characters=max_cell_characters,
    )
    if format == "parquet":
        parquet_status = require_parquet_backend()
        count = _write_parquet(
            prepared,
            path,
            schema=schema,
            overwrite=overwrite,
            max_output_bytes=max_output_bytes,
            compression=parquet_compression,
        )
        backend = _parquet_backend(parquet_status, compression=parquet_compression)
    else:

        def writer(handle: TextIO) -> int:
            bounded = cast(TextIO, _BoundedTableTextWriter(handle, limit=max_output_bytes))
            if format == "json":
                return _write_json(bounded, prepared, schema=schema)
            return _write_delimited(
                bounded,
                prepared,
                columns=schema.columns,
                delimiter="," if format == "csv" else "\t",
                null_value=null_value,
            )

        count = write_text_path(path, writer, overwrite=overwrite, create_parents=False)
        backend = _native_table_backend()
    return TableExportResult(
        format,
        count,
        schema.columns,
        schema.schema_version,
        str(path),
        path.stat().st_size,
        freeze_mapping(
            {
                "max_rows": max_rows,
                "max_columns": max_columns,
                "max_cell_characters": max_cell_characters,
                "max_output_bytes": max_output_bytes,
                "null_value": null_value,
                "parquet_compression": parquet_compression,
                "backend": backend,
            }
        ),
    )


def _validate_table_limits(
    *,
    max_rows: int,
    max_columns: int,
    max_cell_characters: int,
    max_file_bytes: int,
    max_decoded_bytes: int,
) -> None:
    for name, value in (
        ("max_rows", max_rows),
        ("max_columns", max_columns),
        ("max_cell_characters", max_cell_characters),
        ("max_file_bytes", max_file_bytes),
        ("max_decoded_bytes", max_decoded_bytes),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ConfigurationError(
                f"{name} must be positive.",
                code="INVALID_TABLE_LIMIT",
                context={"field": name},
            )


@contextmanager
def _open_table_text(
    source: str | os.PathLike[str] | TextIO, *, max_file_bytes: int
) -> Iterator[tuple[TextIO, str | None]]:
    if isinstance(source, (str, os.PathLike)):
        path = Path(source)
        try:
            stat = path.stat()
            if stat.st_size > max_file_bytes:
                raise InputFormatError(
                    "Table input exceeds max_file_bytes.",
                    code="TABLE_FILE_LIMIT_EXCEEDED",
                    context={"max_file_bytes": max_file_bytes, "source": str(path)},
                )
            with path.open("rb") as probe:
                if probe.read(2) == b"\x1f\x8b":
                    raise InputFormatError(
                        "Generic table reading does not accept gzip-compressed input.",
                        code="COMPRESSED_TABLE_UNSUPPORTED",
                    )
            with path.open("r", encoding="utf-8", newline="") as handle:
                yield handle, str(path)
        except InputFormatError:
            raise
        except (OSError, UnicodeError) as exc:
            raise InputFormatError(
                "Could not open or decode table input.",
                code="TABLE_READ_FAILED",
                context={"source": str(path), "reason": str(exc)},
            ) from exc
        return
    if not hasattr(source, "read"):
        raise TypeError("source must be a path or readable text stream.")
    stream = source
    stream_probe = stream.read(0)
    if not isinstance(stream_probe, str):
        raise TypeError("generic CSV, TSV, and JSON readers require a decoded text stream.")
    yield stream, None


def _read_bounded_text(
    handle: TextIO, *, max_file_bytes: int, limit_error_code: str = "TABLE_FILE_LIMIT_EXCEEDED"
) -> tuple[str, int]:
    chunks: list[str] = []
    byte_count = 0
    while True:
        chunk = handle.read(min(65_536, max_file_bytes - byte_count + 1))
        if not chunk:
            break
        byte_count += len(chunk.encode("utf-8"))
        if byte_count > max_file_bytes:
            raise InputFormatError(
                "Table input exceeds max_file_bytes.",
                code=limit_error_code,
                context={"max_file_bytes": max_file_bytes},
            )
        chunks.append(chunk)
    return "".join(chunks), byte_count


def _parse_delimited_scalar(
    raw: str,
    *,
    schema: TableSchema,
    column: str,
    row_index: int,
    missing_values: tuple[str, ...],
) -> TableScalar:
    if raw in missing_values:
        return _validate_schema_scalar(None, schema=schema, row_index=row_index, column=column)
    expected = schema.column_type(column)
    try:
        if expected in {"any", "string"}:
            value: object = raw
        elif expected == "integer":
            if re.fullmatch(r"[+-]?\d+", raw) is None:
                raise ValueError
            value = int(raw)
        elif expected == "number":
            value = float(raw)
            if not math.isfinite(value):
                raise ValueError
        else:
            lowered = raw.lower()
            if lowered not in {"true", "false"}:
                raise ValueError
            value = lowered == "true"
    except (OverflowError, ValueError) as exc:
        raise InputFormatError(
            "Delimited table cell cannot be decoded as its declared type.",
            code="TABLE_COLUMN_TYPE_MISMATCH",
            context={"row_index": row_index, "column": column, "expected": expected},
        ) from exc
    return _validate_schema_scalar(value, schema=schema, row_index=row_index, column=column)


def _validated_read_row(
    row: Mapping[str, object],
    *,
    schema: TableSchema,
    row_index: int,
    max_cell_characters: int,
) -> FrozenDict:
    if set(row) != set(schema.columns) or any(not isinstance(key, str) for key in row):
        raise InputFormatError(
            "Table row keys must match the declared schema exactly.",
            code="TABLE_SCHEMA_MISMATCH",
            context={"row_index": row_index},
        )
    resolved: dict[str, TableScalar] = {}
    for column in schema.columns:
        value = _validate_schema_scalar(
            row[column], schema=schema, row_index=row_index, column=column
        )
        if isinstance(value, str) and len(value) > max_cell_characters:
            raise InputFormatError(
                "Table cell exceeds max_cell_characters.",
                code="TABLE_CELL_LIMIT_EXCEEDED",
                context={"row_index": row_index, "column": column},
            )
        resolved[column] = value
    return freeze_mapping(resolved)


def _read_delimited_table(
    text: str,
    *,
    delimiter: str,
    schema: TableSchema,
    missing_values: tuple[str, ...],
    max_rows: int,
    max_cell_characters: int,
) -> tuple[FrozenDict, ...]:
    previous_limit = csv.field_size_limit()
    rows: list[FrozenDict] = []
    try:
        csv.field_size_limit(max_cell_characters)
        reader = csv.reader(io.StringIO(text), delimiter=delimiter, strict=True)
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise InputFormatError("Table input is empty.", code="INVALID_TABLE_HEADER") from exc
        if header != schema.columns:
            raise InputFormatError(
                "Delimited table header must match the declared schema exactly.",
                code="TABLE_SCHEMA_MISMATCH",
                context={"actual_columns": header, "expected_columns": schema.columns},
            )
        for row_index, values in enumerate(reader):
            if row_index >= max_rows:
                raise InputFormatError(
                    "Table input exceeds max_rows.",
                    code="TABLE_ROW_LIMIT_EXCEEDED",
                    context={"max_rows": max_rows},
                )
            if len(values) != len(schema.columns):
                raise InputFormatError(
                    "Delimited table row field count does not match its header.",
                    code="TABLE_SCHEMA_MISMATCH",
                    context={"row_index": row_index},
                )
            resolved = {
                column: _parse_delimited_scalar(
                    raw,
                    schema=schema,
                    column=column,
                    row_index=row_index,
                    missing_values=missing_values,
                )
                for column, raw in zip(schema.columns, values, strict=True)
            }
            rows.append(
                _validated_read_row(
                    resolved,
                    schema=schema,
                    row_index=row_index,
                    max_cell_characters=max_cell_characters,
                )
            )
    except csv.Error as exc:
        raise InputFormatError(
            "Delimited table is malformed or exceeds max_cell_characters.",
            code="INVALID_TABLE_DELIMITED",
            context={"reason": str(exc)},
        ) from exc
    finally:
        csv.field_size_limit(previous_limit)
    return tuple(rows)


def _read_json_table(
    text: str,
    *,
    schema: TableSchema,
    max_rows: int,
    max_cell_characters: int,
) -> tuple[FrozenDict, ...]:
    def strict_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON object key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> object:
        raise ValueError(f"non-standard JSON constant: {value}")

    try:
        payload = json.loads(
            text,
            object_pairs_hook=strict_object,
            parse_constant=reject_constant,
        )
    except (json.JSONDecodeError, RecursionError, ValueError) as exc:
        context: dict[str, object] = {}
        if isinstance(exc, json.JSONDecodeError):
            context = {"line_number": exc.lineno, "column_number": exc.colno}
        raise InputFormatError(
            "JSON table input is invalid.",
            code="INVALID_TABLE_JSON",
            context=context,
        ) from exc
    if not isinstance(payload, Mapping) or set(payload) - {
        "schema_version",
        "columns",
        "column_types",
        "nullable",
        "rows",
    }:
        raise InputFormatError("JSON table root/schema is invalid.", code="INVALID_TABLE_JSON")
    metadata = {key: value for key, value in payload.items() if key != "rows"}
    if metadata != _schema_payload(schema):
        raise InputFormatError(
            "JSON table metadata does not match the declared schema.",
            code="TABLE_SCHEMA_MISMATCH",
        )
    raw_rows = payload.get("rows")
    if not isinstance(raw_rows, list):
        raise InputFormatError("JSON table rows must be an array.", code="INVALID_TABLE_JSON")
    if len(raw_rows) > max_rows:
        raise InputFormatError(
            "Table input exceeds max_rows.",
            code="TABLE_ROW_LIMIT_EXCEEDED",
            context={"max_rows": max_rows},
        )
    rows: list[FrozenDict] = []
    for row_index, raw in enumerate(raw_rows):
        if not isinstance(raw, Mapping):
            raise InputFormatError(
                "Each table row must be an object.",
                code="INVALID_TABLE_ROW",
                context={"row_index": row_index},
            )
        rows.append(
            _validated_read_row(
                cast(Mapping[str, object], raw),
                schema=schema,
                row_index=row_index,
                max_cell_characters=max_cell_characters,
            )
        )
    return tuple(rows)


def _read_parquet_table(
    source: str | os.PathLike[str],
    *,
    schema: TableSchema,
    max_rows: int,
    max_cell_characters: int,
    max_file_bytes: int,
    max_decoded_bytes: int,
) -> tuple[tuple[FrozenDict, ...], int, str, dict[str, object]]:
    status = require_parquet_backend()
    path = Path(source)
    try:
        stat = path.stat()
        if stat.st_size > max_file_bytes:
            raise InputFormatError(
                "Table input exceeds max_file_bytes.",
                code="TABLE_FILE_LIMIT_EXCEEDED",
                context={"max_file_bytes": max_file_bytes},
            )
        with path.open("rb") as probe:
            if probe.read(2) == b"\x1f\x8b":
                raise InputFormatError(
                    "Generic Parquet reading does not accept gzip wrappers.",
                    code="COMPRESSED_TABLE_UNSUPPORTED",
                )
        import pyarrow.parquet as pq

        parquet = pq.ParquetFile(path)
        if parquet.metadata.num_rows > max_rows:
            raise InputFormatError(
                "Table input exceeds max_rows.",
                code="TABLE_ROW_LIMIT_EXCEEDED",
                context={"max_rows": max_rows},
            )
        decoded_size = sum(
            parquet.metadata.row_group(row_group).column(column).total_uncompressed_size
            for row_group in range(parquet.metadata.num_row_groups)
            for column in range(parquet.metadata.num_columns)
        )
        if decoded_size > max_decoded_bytes:
            raise InputFormatError(
                "Parquet table exceeds max_decoded_bytes before materialization.",
                code="TABLE_DECODED_LIMIT_EXCEEDED",
                context={
                    "decoded_bytes": decoded_size,
                    "max_decoded_bytes": max_decoded_bytes,
                },
            )
        if tuple(parquet.schema_arrow.names) != schema.columns:
            raise InputFormatError(
                "Parquet columns must match the declared schema exactly.",
                code="TABLE_SCHEMA_MISMATCH",
            )
        metadata = parquet.schema_arrow.metadata or {}
        stored_version = metadata.get(b"dnakit_schema_version")
        if stored_version is not None and stored_version.decode("utf-8") != schema.schema_version:
            raise InputFormatError(
                "Parquet schema version does not match the declared schema.",
                code="TABLE_SCHEMA_MISMATCH",
            )
        stored_schema = metadata.get(b"dnakit_table_schema")
        if stored_schema is not None:
            try:
                embedded_schema = json.loads(stored_schema.decode("utf-8"))
            except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
                raise InputFormatError(
                    "Parquet DNAKit schema metadata is invalid.",
                    code="INVALID_TABLE_SCHEMA_METADATA",
                ) from exc
            if embedded_schema != _schema_payload(schema):
                raise InputFormatError(
                    "Parquet metadata does not match the declared schema.",
                    code="TABLE_SCHEMA_MISMATCH",
                )
        rows: list[FrozenDict] = []
        for batch in parquet.iter_batches(batch_size=min(max_rows, 4_096)):
            for raw in batch.to_pylist():
                row_index = len(rows)
                rows.append(
                    _validated_read_row(
                        cast(Mapping[str, object], raw),
                        schema=schema,
                        row_index=row_index,
                        max_cell_characters=max_cell_characters,
                    )
                )
        compressions = sorted(
            {
                parquet.metadata.row_group(row_group).column(column).compression.lower()
                for row_group in range(parquet.metadata.num_row_groups)
                for column in range(parquet.metadata.num_columns)
            }
        )
        backend = _parquet_backend(status, compression=compressions)
        return tuple(rows), stat.st_size, str(path), backend
    except InputFormatError:
        raise
    except (OSError, UnicodeError) as exc:
        raise InputFormatError(
            "Could not open Parquet table input.",
            code="TABLE_READ_FAILED",
            context={"source": str(path), "reason": str(exc)},
        ) from exc
    except DNAKitError:
        raise
    except Exception as exc:
        raise BackendExecutionError(
            "PyArrow failed to read the bounded table.",
            code="PARQUET_BACKEND_EXECUTION_FAILED",
            context={"reason": str(exc)},
        ) from exc


def read_table(
    source: str | os.PathLike[str] | TextIO,
    *,
    format: TableFormat,
    schema: TableSchema,
    missing_values: Sequence[str] = (r"\N",),
    max_rows: int = 1_000_000,
    max_columns: int = 10_000,
    max_cell_characters: int = 10_000_000,
    max_file_bytes: int = 1_000_000_000,
    max_decoded_bytes: int = 1_000_000_000,
) -> TableReadResult:
    """Strictly read one bounded two-dimensional table with an explicit schema."""

    if format not in {"csv", "tsv", "json", "parquet"}:
        raise ConfigurationError("Unknown table read format.", code="INVALID_TABLE_FORMAT")
    if not isinstance(schema, TableSchema):
        raise TypeError("schema must be a TableSchema.")
    _validate_table_limits(
        max_rows=max_rows,
        max_columns=max_columns,
        max_cell_characters=max_cell_characters,
        max_file_bytes=max_file_bytes,
        max_decoded_bytes=max_decoded_bytes,
    )
    if len(schema.columns) > max_columns:
        raise InputFormatError(
            "Table schema exceeds max_columns.",
            code="TABLE_COLUMN_LIMIT_EXCEEDED",
            context={"max_columns": max_columns},
        )
    resolved_missing = tuple(missing_values)
    if any(not isinstance(value, str) for value in resolved_missing) or len(
        set(resolved_missing)
    ) != len(resolved_missing):
        raise ConfigurationError(
            "missing_values must contain unique strings.", code="INVALID_TABLE_MISSING_VALUES"
        )
    source_path: str | None
    if format == "parquet":
        if not isinstance(source, (str, os.PathLike)):
            raise TypeError("Parquet table reading requires a local filesystem path.")
        rows, byte_count, source_path, backend = _read_parquet_table(
            source,
            schema=schema,
            max_rows=max_rows,
            max_cell_characters=max_cell_characters,
            max_file_bytes=max_file_bytes,
            max_decoded_bytes=max_decoded_bytes,
        )
    else:
        with _open_table_text(source, max_file_bytes=max_file_bytes) as (handle, source_path):
            text_value, byte_count = _read_bounded_text(
                handle,
                max_file_bytes=min(max_file_bytes, max_decoded_bytes),
                limit_error_code=(
                    "TABLE_DECODED_LIMIT_EXCEEDED"
                    if max_decoded_bytes < max_file_bytes
                    else "TABLE_FILE_LIMIT_EXCEEDED"
                ),
            )
        if format == "json":
            rows = _read_json_table(
                text_value,
                schema=schema,
                max_rows=max_rows,
                max_cell_characters=max_cell_characters,
            )
        else:
            rows = _read_delimited_table(
                text_value,
                delimiter="," if format == "csv" else "\t",
                schema=schema,
                missing_values=resolved_missing,
                max_rows=max_rows,
                max_cell_characters=max_cell_characters,
            )
        backend = _native_table_backend()
    return TableReadResult(
        format,
        rows,
        schema,
        source_path,
        byte_count,
        freeze_mapping(
            {
                "missing_values": resolved_missing,
                "max_rows": max_rows,
                "max_columns": max_columns,
                "max_cell_characters": max_cell_characters,
                "max_file_bytes": max_file_bytes,
                "max_decoded_bytes": max_decoded_bytes,
                "backend": backend,
            }
        ),
    )


def export_result(
    result: DictResult,
    target: str | os.PathLike[str],
    *,
    format: TableFormat = "json",
    overwrite: bool = False,
) -> TableExportResult:
    """Export one flat result object's scalar fields as a one-row table."""

    if not isinstance(result, DictResult):
        raise TypeError("result must provide to_dict().")
    raw_row = result.to_dict()
    if not isinstance(raw_row, Mapping) or not raw_row:
        raise InputFormatError(
            "Result to_dict() must return one non-empty flat mapping.",
            code="INVALID_RESULT_EXPORT",
        )
    raw_items = tuple(islice(raw_row.items(), _MAX_TABLE_SCHEMA_COLUMNS + 1))
    if len(raw_items) > _MAX_TABLE_SCHEMA_COLUMNS:
        raise InputFormatError(
            "Result to_dict() exceeds the hard column limit.",
            code="INVALID_RESULT_EXPORT",
            context={"max_columns": _MAX_TABLE_SCHEMA_COLUMNS},
        )
    row: dict[str, object] = {}
    for key, value in raw_items:
        if not isinstance(key, str):
            raise InputFormatError(
                "Result to_dict() keys must be strings.", code="INVALID_RESULT_EXPORT"
            )
        if value is None or isinstance(value, (bool, int, float, str)):
            row[key] = value
        else:
            try:
                row[key] = json.dumps(
                    to_json_compatible(value),
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
            except (DNAKitError, RecursionError, TypeError, ValueError) as exc:
                raise InputFormatError(
                    "Result fields must be JSON-compatible.",
                    code="INVALID_RESULT_EXPORT",
                    context={"field": key},
                ) from exc
    columns = tuple(row)
    return export_table(
        (row,),
        target,
        format=format,
        schema=TableSchema(columns, schema_version="dnakit.result-table.v1"),
        overwrite=overwrite,
    )


__all__ = [
    "DictResult",
    "ParquetCompression",
    "TableColumnType",
    "TableExportResult",
    "TableFormat",
    "TableReadResult",
    "TableScalar",
    "TableSchema",
    "export_result",
    "export_table",
    "read_table",
]
