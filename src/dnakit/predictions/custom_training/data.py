"""Bounded, format-independent data ingestion for custom property training."""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import mimetypes
import re
import sqlite3
import xml.etree.ElementTree as ET
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Any, Protocol, cast, runtime_checkable

from dnakit.core._json import to_json_compatible
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    InputFormatError,
)
from dnakit.io import ReadConfig, read

from .llm import (
    OpenAICompatibleSchemaClient,
    SchemaLLMClient,
    mapping_from_llm_proposal,
)
from .schema import DataInspection, DataMapping, LLMSchemaConfig

_TEXT_PROBE_BYTES = 65_536
_IUPAC = frozenset("ACGTRYSWKMBDHVNU")
_SEQUENCE_NAMES = (
    "sequence",
    "seq",
    "dna",
    "dna_sequence",
    "nucleotide_sequence",
    "nucleotides",
)
_ID_NAMES = ("id", "index", "sample_id", "sequence_id", "seq_id", "identifier", "name")
_SPLIT_NAMES = ("split", "subset", "partition", "set")
_LABEL_NAMES = ("label", "target", "class", "score", "property", "value", "y")
_FORMAT_SUFFIXES = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".tab": "tsv",
    ".json": "json",
    ".jsonl": "jsonl",
    ".ndjson": "jsonl",
    ".yaml": "yaml",
    ".yml": "yaml",
    ".xml": "xml",
    ".xlsx": "excel",
    ".xls": "excel",
    ".parquet": "parquet",
    ".pq": "parquet",
    ".sqlite": "sqlite",
    ".sqlite3": "sqlite",
    ".db": "sqlite",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".ffn": "fasta",
    ".fq": "fastq",
    ".fastq": "fastq",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbff": "genbank",
    ".genbank": "genbank",
}


@runtime_checkable
class DataReader(Protocol):
    """Extension point for formats outside DNAKit's built-in reader set."""

    def read(
        self,
        source: Path,
        *,
        mapping: DataMapping,
        max_records: int,
        max_file_bytes: int,
    ) -> Iterable[Mapping[str, object]]:
        """Yield bounded record mappings from ``source``."""


@dataclass(frozen=True, slots=True)
class PropertyData:
    """Validated rows ready for target encoding and deterministic splitting."""

    record_ids: tuple[str, ...]
    sequences: tuple[str, ...]
    labels: tuple[tuple[object, ...], ...]
    label_names: tuple[str, ...]
    supplied_splits: tuple[str | None, ...]
    inspection: DataInspection
    dropped_records: int = 0
    features: tuple[tuple[float, ...], ...] = ()

    @property
    def record_count(self) -> int:
        return len(self.record_ids)


def _source_path(source: str | Path) -> Path:
    try:
        path = Path(source).expanduser().resolve()
    except (TypeError, ValueError, OSError) as exc:
        raise ConfigurationError(
            "source must be a valid local filesystem path.",
            code="INVALID_CUSTOM_DATA_SOURCE",
        ) from exc
    if not path.is_file():
        raise InputFormatError(
            "The custom-property data source is not a file.",
            code="CUSTOM_DATA_SOURCE_NOT_FOUND",
            context={"source": str(path)},
        )
    return path


def _check_source_size(path: Path, maximum: int) -> int:
    try:
        size = path.stat().st_size
    except OSError as exc:
        raise InputFormatError(
            "Could not inspect the custom-property data source.",
            code="CUSTOM_DATA_READ_FAILED",
            context={"source": str(path)},
        ) from exc
    if size > maximum:
        raise InputFormatError(
            "The custom-property data source exceeds max_file_bytes.",
            code="CUSTOM_DATA_FILE_LIMIT",
            context={"source_size_bytes": size, "max_file_bytes": maximum},
        )
    return size


def _inner_suffix(path: Path) -> str:
    suffixes = [value.lower() for value in path.suffixes]
    if suffixes and suffixes[-1] == ".gz" and len(suffixes) > 1:
        return suffixes[-2]
    return suffixes[-1] if suffixes else ""


def _probe_bytes(path: Path) -> bytes:
    try:
        if path.suffix.lower() == ".gz":
            with gzip.open(path, "rb") as handle:
                return handle.read(_TEXT_PROBE_BYTES)
        with path.open("rb") as handle:
            return handle.read(_TEXT_PROBE_BYTES)
    except (OSError, gzip.BadGzipFile) as exc:
        raise InputFormatError(
            "Could not read the custom-property data source.",
            code="CUSTOM_DATA_READ_FAILED",
            context={"source": str(path)},
        ) from exc


def _looks_textual(raw: bytes) -> bool:
    if not raw:
        return True
    if b"\x00" in raw:
        return False
    try:
        raw.decode("utf-8")
    except UnicodeDecodeError:
        return False
    control = sum(value < 9 or 13 < value < 32 for value in raw)
    return control / len(raw) < 0.01


def detect_data_format(path: str | Path) -> str:
    """Detect a built-in format from extension, magic bytes, and a bounded probe."""

    path = _source_path(path)
    suffix = _inner_suffix(path)
    if suffix in _FORMAT_SUFFIXES:
        return _FORMAT_SUFFIXES[suffix]
    raw = _probe_bytes(path)
    if raw.startswith(b"PAR1"):
        return "parquet"
    if raw.startswith(b"SQLite format 3\x00"):
        return "sqlite"
    if raw.startswith(b"PK\x03\x04"):
        return "excel"
    if not _looks_textual(raw):
        return "unknown"
    text = raw.decode("utf-8-sig", errors="strict").lstrip()
    if text.startswith(">"):
        return "fasta"
    if text.startswith("LOCUS"):
        return "genbank"
    if text.startswith("<?xml") or text.startswith("<"):
        return "xml"
    if text.startswith("[") or text.startswith("{"):
        try:
            json.loads(text)
        except json.JSONDecodeError:
            return "jsonl" if text.startswith("{") else "json"
        return "json"
    lines = text.splitlines()
    if len(lines) >= 4 and lines[0].startswith("@") and lines[2].startswith("+"):
        return "fastq"
    try:
        dialect = csv.Sniffer().sniff(text[:32_768], delimiters=",\t;|")
    except csv.Error:
        return "unknown"
    return (
        "tsv" if dialect.delimiter == "\t" else "csv" if dialect.delimiter == "," else "delimited"
    )


@contextmanager
def _open_text(path: Path, *, encoding: str, max_file_bytes: int) -> Iterator[io.TextIOBase]:
    raw: Any
    try:
        raw = (
            gzip.open(path, "rb")  # noqa: SIM115 - TextIOWrapper owns and closes it.
            if path.suffix.lower() == ".gz"
            else path.open("rb")
        )
        bounded = _BoundedBinaryReader(raw, max_file_bytes)
        text = io.TextIOWrapper(cast(Any, bounded), encoding=encoding, newline="")
    except (LookupError, OSError) as exc:
        raise InputFormatError(
            "Could not open the custom-property text source.",
            code="CUSTOM_DATA_READ_FAILED",
            context={"source": str(path), "encoding": encoding},
        ) from exc
    try:
        yield text
    except UnicodeError as exc:
        raise InputFormatError(
            "The custom-property source could not be decoded.",
            code="CUSTOM_DATA_DECODE_FAILED",
            context={"encoding": encoding},
        ) from exc
    finally:
        text.close()


class _BoundedBinaryReader(io.RawIOBase):
    def __init__(self, raw: Any, limit: int) -> None:
        self.raw = raw
        self.limit = limit
        self.count = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        remaining = self.limit - self.count
        if remaining < 1:
            probe = self.raw.read(1)
            if probe:
                raise InputFormatError(
                    "Decoded custom-property input exceeds max_file_bytes.",
                    code="CUSTOM_DATA_FILE_LIMIT",
                    context={"max_file_bytes": self.limit},
                )
            return 0
        chunk = self.raw.read(min(len(buffer), remaining + 1))
        if len(chunk) > remaining:
            raise InputFormatError(
                "Decoded custom-property input exceeds max_file_bytes.",
                code="CUSTOM_DATA_FILE_LIMIT",
                context={"max_file_bytes": self.limit},
            )
        size = len(chunk)
        buffer[:size] = chunk
        self.count += size
        return size


def _bounded_records(
    records: Iterable[Mapping[str, object]],
    *,
    max_records: int,
) -> tuple[dict[str, object], ...]:
    materialized: list[dict[str, object]] = []
    try:
        iterator = iter(records)
    except TypeError as exc:
        raise InputFormatError(
            "The data reader did not return an iterable.", code="INVALID_CUSTOM_DATA_READER"
        ) from exc
    for index, record in enumerate(iterator):
        if index >= max_records:
            raise InputFormatError(
                "The custom-property source exceeds max_records.",
                code="CUSTOM_DATA_RECORD_LIMIT",
                context={"max_records": max_records},
            )
        if not isinstance(record, Mapping):
            raise InputFormatError(
                "Every custom-property source record must be a mapping.",
                code="INVALID_CUSTOM_DATA_RECORD",
                context={"record_index": index, "type": type(record).__name__},
            )
        materialized.append({str(key): value for key, value in record.items()})
    if not materialized:
        raise InputFormatError(
            "The custom-property source contains no records.", code="EMPTY_CUSTOM_DATA"
        )
    return tuple(materialized)


def _read_delimited(
    path: Path,
    mapping: DataMapping,
    *,
    maximum: int,
) -> Iterable[Mapping[str, object]]:
    delimiter = mapping.delimiter
    if delimiter is None:
        delimiter = "\t" if mapping.format == "tsv" else ","
        if mapping.format == "delimited":
            probe = _probe_bytes(path).decode(mapping.encoding, errors="strict")
            try:
                delimiter = csv.Sniffer().sniff(probe, delimiters=",\t;|").delimiter
            except csv.Error as exc:
                raise InputFormatError(
                    "Could not infer a delimiter for the custom-property source.",
                    code="CUSTOM_DATA_DELIMITER_UNKNOWN",
                    hint="Set DataMapping(delimiter=...).",
                ) from exc
    with _open_text(path, encoding=mapping.encoding, max_file_bytes=maximum) as handle:
        for _ in range(mapping.header_row):
            if handle.readline() == "":
                raise InputFormatError(
                    "header_row lies beyond the end of the source.",
                    code="INVALID_CUSTOM_DATA_HEADER",
                )
        previous_limit = csv.field_size_limit()
        try:
            csv.field_size_limit(min(maximum, 10_000_000))
            reader = csv.DictReader(handle, delimiter=delimiter)
            if reader.fieldnames is None or any(
                name is None or not name for name in reader.fieldnames
            ):
                raise InputFormatError(
                    "The delimited source has an empty or missing header.",
                    code="INVALID_CUSTOM_DATA_HEADER",
                )
            if len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise InputFormatError(
                    "The delimited source header contains duplicate fields.",
                    code="INVALID_CUSTOM_DATA_HEADER",
                )
            yield from reader
        finally:
            csv.field_size_limit(previous_limit)


def _resolve_record_path(value: object, path: str | None) -> object:
    current = value
    if path is None:
        return current
    for part in path.split("."):
        if not isinstance(current, Mapping) or part not in current:
            raise InputFormatError(
                "record_path does not exist in the structured source.",
                code="CUSTOM_DATA_RECORD_PATH_MISSING",
                context={"record_path": path, "missing_part": part},
            )
        current = current[part]
    return current


def _structured_records(value: object, record_path: str | None) -> Iterable[Mapping[str, object]]:
    selected = _resolve_record_path(value, record_path)
    if isinstance(selected, Mapping):
        list_values = [item for item in selected.values() if isinstance(item, list)]
        selected = list_values[0] if record_path is None and len(list_values) == 1 else [selected]
    if not isinstance(selected, list):
        raise InputFormatError(
            "The structured source record container must be an array.",
            code="INVALID_CUSTOM_DATA_RECORDS",
            hint="Set record_path to the array containing supervised records.",
        )
    for index, item in enumerate(selected):
        if not isinstance(item, Mapping):
            raise InputFormatError(
                "Each structured source record must be an object.",
                code="INVALID_CUSTOM_DATA_RECORD",
                context={"record_index": index},
            )
        yield cast(Mapping[str, object], item)


def _read_json(path: Path, mapping: DataMapping, *, maximum: int) -> Iterable[Mapping[str, object]]:
    with _open_text(path, encoding=mapping.encoding, max_file_bytes=maximum) as handle:
        try:
            value = json.load(handle)
        except json.JSONDecodeError as exc:
            raise InputFormatError(
                "The custom-property JSON source is invalid.",
                code="INVALID_CUSTOM_DATA_JSON",
                context={"line": exc.lineno, "column": exc.colno},
            ) from exc
    return _structured_records(value, mapping.record_path)


def _read_jsonl(
    path: Path, mapping: DataMapping, *, maximum: int
) -> Iterable[Mapping[str, object]]:
    with _open_text(path, encoding=mapping.encoding, max_file_bytes=maximum) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise InputFormatError(
                    "The custom-property JSONL source is invalid.",
                    code="INVALID_CUSTOM_DATA_JSON",
                    context={"line": line_number, "column": exc.colno},
                ) from exc
            if not isinstance(value, Mapping):
                raise InputFormatError(
                    "Each JSONL line must be an object.",
                    code="INVALID_CUSTOM_DATA_RECORD",
                    context={"line": line_number},
                )
            yield cast(Mapping[str, object], value)


def _read_yaml(path: Path, mapping: DataMapping, *, maximum: int) -> Iterable[Mapping[str, object]]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - PyYAML is a base dependency.
        raise BackendUnavailableError(
            "YAML input requires PyYAML.", code="MISSING_CUSTOM_DATA_DEPENDENCY"
        ) from exc
    with _open_text(path, encoding=mapping.encoding, max_file_bytes=maximum) as handle:
        try:
            value = yaml.safe_load(handle)
        except yaml.YAMLError as exc:
            raise InputFormatError(
                "The custom-property YAML source is invalid.",
                code="INVALID_CUSTOM_DATA_YAML",
            ) from exc
    return _structured_records(value, mapping.record_path)


def _xml_record(element: ET.Element) -> dict[str, object]:
    row: dict[str, object] = {f"@{key}": value for key, value in element.attrib.items()}
    if element.text and element.text.strip():
        row["text"] = element.text.strip()
    for child in element:
        value: object = _xml_record(child) if child else (child.text or "").strip()
        if child.tag in row:
            existing = row[child.tag]
            row[child.tag] = [*existing, value] if isinstance(existing, list) else [existing, value]
        else:
            row[child.tag] = value
    return row


class _DoctypeRejectingReader:
    """Proxy a text stream while rejecting XML document type declarations."""

    _MARKER = "<!doctype"

    def __init__(self, handle: io.TextIOBase) -> None:
        self._handle = handle
        self._tail = ""

    def read(self, size: int = -1) -> str:
        chunk = self._handle.read(size)
        combined = (self._tail + chunk).lower()
        if self._MARKER in combined:
            raise InputFormatError(
                "XML document type declarations are not accepted.",
                code="UNSAFE_CUSTOM_DATA_XML",
            )
        self._tail = combined[-(len(self._MARKER) - 1) :]
        return chunk


def _read_xml(path: Path, mapping: DataMapping, *, maximum: int) -> Iterable[Mapping[str, object]]:
    with _open_text(path, encoding=mapping.encoding, max_file_bytes=maximum) as handle:
        try:
            guarded_handle = _DoctypeRejectingReader(handle)
            root = ET.parse(guarded_handle).getroot()
        except ET.ParseError as exc:
            raise InputFormatError(
                "The custom-property XML source is invalid.",
                code="INVALID_CUSTOM_DATA_XML",
            ) from exc
    elements: Sequence[ET.Element]
    if mapping.record_path is None:
        elements = list(root)
    else:
        parts = mapping.record_path.split(".")
        elements = [root]
        if parts and parts[0] == root.tag:
            parts = parts[1:]
        for part in parts:
            elements = [child for parent in elements for child in parent if child.tag == part]
    if not elements:
        raise InputFormatError(
            "No XML records matched record_path.", code="CUSTOM_DATA_RECORD_PATH_MISSING"
        )
    return (_xml_record(element) for element in elements)


def _read_excel(path: Path, mapping: DataMapping) -> Iterable[Mapping[str, object]]:
    try:
        import pandas as pd  # type: ignore[import-untyped]
    except ImportError as exc:
        raise BackendUnavailableError(
            "Excel input requires pandas and openpyxl.",
            code="MISSING_CUSTOM_DATA_DEPENDENCY",
            hint='Install the training extra with: python -m pip install "dnakit[training]"',
        ) from exc
    try:
        frame = pd.read_excel(
            path,
            sheet_name=0 if mapping.sheet_name is None else mapping.sheet_name,
            header=mapping.header_row,
        )
    except Exception as exc:
        raise BackendExecutionError(
            "Could not read the custom-property Excel source.",
            code="CUSTOM_DATA_READ_FAILED",
            context={"error_type": type(exc).__name__},
        ) from exc
    if not hasattr(frame, "to_dict"):
        raise InputFormatError(
            "Excel sheet selection did not resolve one table.",
            code="INVALID_CUSTOM_DATA_RECORDS",
        )
    return cast(Iterable[Mapping[str, object]], frame.to_dict(orient="records"))


def _read_parquet(path: Path) -> Iterable[Mapping[str, object]]:
    try:
        import pyarrow.parquet as pq  # type: ignore[import-untyped]
    except ImportError as exc:
        raise BackendUnavailableError(
            "Parquet input requires PyArrow.",
            code="MISSING_CUSTOM_DATA_DEPENDENCY",
            hint='Install the I/O extra with: python -m pip install "dnakit[io]"',
        ) from exc
    try:
        return cast(Iterable[Mapping[str, object]], pq.read_table(path).to_pylist())
    except Exception as exc:
        raise BackendExecutionError(
            "Could not read the custom-property Parquet source.",
            code="CUSTOM_DATA_READ_FAILED",
            context={"error_type": type(exc).__name__},
        ) from exc


def _sqlite_tables(connection: sqlite3.Connection) -> tuple[str, ...]:
    cursor = connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table' "
        "AND name NOT LIKE 'sqlite_%' ORDER BY name"
    )
    return tuple(str(row[0]) for row in cursor.fetchall())


def _read_sqlite(
    path: Path,
    mapping: DataMapping,
    *,
    max_records: int,
) -> Iterable[Mapping[str, object]]:
    try:
        connection = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        raise InputFormatError(
            "Could not open the SQLite source in read-only mode.",
            code="CUSTOM_DATA_READ_FAILED",
        ) from exc
    try:
        tables = _sqlite_tables(connection)
        table = mapping.table_name
        if table is None:
            if len(tables) != 1:
                raise InputFormatError(
                    "SQLite input requires table_name when the database does not have "
                    "exactly one table.",
                    code="CUSTOM_DATA_TABLE_REQUIRED",
                    context={"tables": tables},
                )
            table = tables[0]
        if table not in tables:
            raise InputFormatError(
                "The selected SQLite table does not exist.",
                code="CUSTOM_DATA_TABLE_MISSING",
                context={"table_name": table, "tables": tables},
            )
        escaped = table.replace('"', '""')
        cursor = connection.execute(f'SELECT * FROM "{escaped}" LIMIT ?', (max_records + 1,))
        names = tuple(str(item[0]) for item in cast(Sequence[Sequence[object]], cursor.description))
        rows = cursor.fetchall()
        return tuple(dict(zip(names, row, strict=True)) for row in rows)
    except sqlite3.Error as exc:
        raise BackendExecutionError(
            "Could not read the selected SQLite table.",
            code="CUSTOM_DATA_READ_FAILED",
        ) from exc
    finally:
        connection.close()


def _header_values(description: str, separator: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for token in re.split(rf"(?:{re.escape(separator)}|;|\s)+", description.strip()):
        key, marker, value = token.partition("=")
        if marker and key and value and key not in fields:
            fields[key] = value
    return fields


def _read_sequence_file(
    path: Path,
    mapping: DataMapping,
    *,
    data_format: str,
    max_records: int,
    max_file_bytes: int,
) -> Iterable[Mapping[str, object]]:
    config = ReadConfig(
        encoding=mapping.encoding,
        max_records=max_records,
        max_input_bytes=max_file_bytes,
    )
    source = read(path, format=data_format, config=config, mode="stream")
    with source:
        return tuple(
            {
                "id": record.id,
                "sequence": record.sequence.symbols,
                "description": record.description,
                "header": _header_values(record.description, mapping.fasta_header_separator),
                "metadata": to_json_compatible(record.metadata),
            }
            for record in source
        )


def _read_builtin(
    path: Path,
    mapping: DataMapping,
    *,
    detected_format: str,
    max_records: int,
    max_file_bytes: int,
) -> tuple[dict[str, object], ...]:
    data_format = detected_format if mapping.format == "auto" else mapping.format
    if data_format == "unknown" or data_format == "auto":
        raise InputFormatError(
            "DNAKit could not determine how to parse this source.",
            code="CUSTOM_DATA_FORMAT_UNKNOWN",
            hint="Provide DataMapping(format=...) or a custom DataReader.",
        )
    if data_format in {"csv", "tsv", "delimited"}:
        records = _read_delimited(path, mapping, maximum=max_file_bytes)
    elif data_format == "json":
        records = _read_json(path, mapping, maximum=max_file_bytes)
    elif data_format == "jsonl":
        records = _read_jsonl(path, mapping, maximum=max_file_bytes)
    elif data_format == "yaml":
        records = _read_yaml(path, mapping, maximum=max_file_bytes)
    elif data_format == "xml":
        records = _read_xml(path, mapping, maximum=max_file_bytes)
    elif data_format == "excel":
        records = _read_excel(path, mapping)
    elif data_format == "parquet":
        records = _read_parquet(path)
    elif data_format == "sqlite":
        records = _read_sqlite(path, mapping, max_records=max_records)
    elif data_format in {"fasta", "fastq", "genbank"}:
        records = _read_sequence_file(
            path,
            mapping,
            data_format=data_format,
            max_records=max_records,
            max_file_bytes=max_file_bytes,
        )
    else:  # pragma: no cover - DataMapping validates all choices.
        raise AssertionError(f"Unsupported custom data format: {data_format}")
    return _bounded_records(records, max_records=max_records)


def _flatten_keys(value: Mapping[str, object], prefix: str = "") -> tuple[str, ...]:
    keys: list[str] = []
    for raw_key, child in value.items():
        key = f"{prefix}.{raw_key}" if prefix else str(raw_key)
        keys.append(key)
        if isinstance(child, Mapping):
            keys.extend(_flatten_keys(cast(Mapping[str, object], child), key))
    return tuple(keys)


def _field(record: Mapping[str, object], path: str) -> object:
    if path in record:
        return record[path]
    current: object = record
    for part in path.split("."):
        if not isinstance(current, Mapping) or part not in current:
            raise KeyError(path)
        current = current[part]
    return current


def _normalized_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _named_column(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    normalized = {_normalized_name(column): column for column in columns}
    return next((normalized[name] for name in candidates if name in normalized), None)


def _dna_score(records: Sequence[Mapping[str, object]], field: str) -> float:
    scores: list[float] = []
    for record in records:
        try:
            value = _field(record, field)
        except KeyError:
            continue
        if not isinstance(value, str):
            continue
        sequence = "".join(value.split()).upper()
        if len(sequence) < 4:
            continue
        scores.append(sum(symbol in _IUPAC for symbol in sequence) / len(sequence))
    return sum(scores) / len(scores) if scores else 0.0


def _local_mapping(
    records: Sequence[Mapping[str, object]],
    *,
    detected_format: str,
) -> DataMapping:
    columns = tuple(dict.fromkeys(key for record in records for key in _flatten_keys(record)))
    sequence_field = _named_column(columns, _SEQUENCE_NAMES)
    if sequence_field is None and columns:
        scored = sorted(
            ((_dna_score(records, column), column) for column in columns),
            reverse=True,
        )
        if scored[0][0] >= 0.9:
            sequence_field = scored[0][1]
    sequence_field = sequence_field or "sequence"
    id_field = _named_column(columns, _ID_NAMES)
    split_field = _named_column(columns, _SPLIT_NAMES)
    excluded = {sequence_field, id_field, split_field, "description", "metadata", "header"}
    labels = [
        column
        for candidate in _LABEL_NAMES
        if (column := _named_column(columns, (candidate,))) is not None and column not in excluded
    ]
    if not labels:
        labels = [
            column
            for column in columns
            if column not in excluded
            and not column.startswith(("header.", "metadata."))
            and all(not isinstance(record.get(column), (dict, list, tuple)) for record in records)
        ]
    resolved_format = "auto" if detected_format == "unknown" else detected_format
    return DataMapping(
        format=cast(Any, resolved_format),
        sequence_field=sequence_field,
        label_fields=tuple(dict.fromkeys(labels)),
        id_field=id_field,
        split_field=split_field,
    )


def _safe_sample(value: object, *, include_values: bool) -> object:
    if not include_values:
        return {"type": type(value).__name__}
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, str):
        return value[:512]
    if isinstance(value, Mapping):
        return {
            str(key)[:128]: _safe_sample(child, include_values=True)
            for key, child in list(value.items())[:50]
        }
    if isinstance(value, (list, tuple)):
        return [_safe_sample(child, include_values=True) for child in value[:20]]
    return {"type": type(value).__name__}


def _llm_profile(
    path: Path,
    *,
    detected_format: str,
    source_size: int,
    records: Sequence[Mapping[str, object]],
    estimated_records: int | None,
    config: LLMSchemaConfig,
) -> tuple[dict[str, object], int, int]:
    columns = tuple(dict.fromkeys(key for record in records for key in _flatten_keys(record)))
    field_types: dict[str, list[str]] = {}
    for column in columns[:1_000]:
        observed: set[str] = set()
        for record in records:
            try:
                observed.add(type(_field(record, column)).__name__)
            except KeyError:
                continue
        field_types[column] = sorted(observed)
    shared_row_limit = min(config.sample_rows, len(records))
    if config.include_sample_values:
        if estimated_records is None:
            shared_row_limit = min(shared_row_limit, max(0, len(records) - 1))
        elif estimated_records <= shared_row_limit:
            shared_row_limit = max(0, estimated_records - 1)
    samples = [
        {
            str(key)[:128]: _safe_sample(value, include_values=config.include_sample_values)
            for key, value in list(record.items())[:100]
        }
        for record in records[:shared_row_limit]
    ]
    serialized = json.dumps(samples, ensure_ascii=False)
    if len(serialized) > config.max_sample_characters:
        serialized = serialized[: config.max_sample_characters]
        samples_payload: object = {
            "truncated_json_prefix": serialized,
            "note": "The sample was truncated at the configured character limit.",
        }
    else:
        samples_payload = samples
    raw = _probe_bytes(path)
    textual = _looks_textual(raw)
    text_excerpt: object = None
    if not records and textual:
        decoded = raw.decode("utf-8-sig", errors="strict")
        if config.include_sample_values:
            safe_excerpt_length = min(config.max_sample_characters, len(decoded) // 2)
            text_excerpt = decoded[:safe_excerpt_length]
        else:
            lines = decoded.splitlines()[: config.sample_rows]
            text_excerpt = {
                "line_count_in_probe": len(decoded.splitlines()),
                "sample_line_lengths": [len(line) for line in lines],
                "separator_counts": {
                    repr(separator): sum(line.count(separator) for line in lines)
                    for separator in (",", "\t", ";", "|", " ")
                },
            }
    profile = {
        "file_name": path.name,
        "extension": "".join(path.suffixes[-2:]),
        "media_type": mimetypes.guess_type(path.name)[0],
        "size_bytes": source_size,
        "detected_format": detected_format,
        "probe_size_bytes": len(raw),
        "textual_probe": textual,
        "columns_or_paths": list(columns[:1_000]),
        "field_types": field_types,
        "sample_records": samples_payload,
        "unparsed_text_excerpt": text_excerpt,
        "sample_policy": {
            "whole_file_shared": False,
            "maximum_records": config.sample_rows,
            "maximum_characters": config.max_sample_characters,
            "values_included": config.include_sample_values,
            "complete_source_values_shared": False,
        },
    }
    shared_characters = min(len(serialized), config.max_sample_characters)
    if isinstance(text_excerpt, str):
        shared_characters += len(text_excerpt)
    shared_characters = min(shared_characters, config.max_sample_characters)
    return profile, len(samples), shared_characters


def inspect_property_data(
    source: str | Path,
    *,
    mapping: DataMapping | None = None,
    llm: LLMSchemaConfig | None = None,
    llm_client: SchemaLLMClient | None = None,
    reader: DataReader | None = None,
    max_records: int = 1_000_000,
    max_file_bytes: int = 1_000_000_000,
) -> DataInspection:
    """Inspect a local source and select an explicit, auditable data mapping.

    When ``llm`` is provided, only the bounded profile described by its sampling
    options is sent to the configured API.  The complete source is never uploaded.
    """

    path = _source_path(source)
    size = _check_source_size(path, max_file_bytes)
    detected = "custom" if reader is not None else detect_data_format(path)
    probe_mapping = mapping or DataMapping(
        format=cast(Any, "auto" if detected in {"custom", "unknown"} else detected)
    )
    samples: tuple[dict[str, object], ...] = ()
    estimated_records: int | None = None
    sample_warning: list[str] = []
    try:
        if reader is not None:
            preview_limit = min(max_records, (llm.sample_rows + 1 if llm else 8))
            samples = _bounded_records(
                islice(
                    reader.read(
                        path,
                        mapping=probe_mapping,
                        max_records=max_records,
                        max_file_bytes=max_file_bytes,
                    ),
                    preview_limit,
                ),
                max_records=preview_limit,
            )
            if llm is not None and len(samples) <= llm.sample_rows:
                estimated_records = len(samples)
        elif detected != "unknown":
            all_probe_records = _read_builtin(
                path,
                probe_mapping,
                detected_format=detected,
                max_records=max_records,
                max_file_bytes=max_file_bytes,
            )
            samples = all_probe_records[: (llm.sample_rows + 1 if llm else 8)]
            estimated_records = len(all_probe_records)
    except (BackendExecutionError, BackendUnavailableError, InputFormatError) as exc:
        if mapping is not None or llm is None:
            raise
        sample_warning.append(f"Local preview was unavailable: {exc.code}")
    local = mapping or _local_mapping(samples, detected_format=detected)
    recommendations: tuple[str, ...] = ()
    llm_warnings: tuple[str, ...] = ()
    origin: str = "explicit" if mapping is not None else "local"
    llm_model: str | None = None
    shared_rows = 0
    shared_characters = 0
    if llm is not None and mapping is None:
        profile, shared_rows, shared_characters = _llm_profile(
            path,
            detected_format=detected,
            source_size=size,
            records=samples,
            estimated_records=estimated_records,
            config=llm,
        )
        client = OpenAICompatibleSchemaClient() if llm_client is None else llm_client
        if not isinstance(client, SchemaLLMClient):
            raise ConfigurationError(
                "llm_client must implement SchemaLLMClient.",
                code="INVALID_LLM_SCHEMA_CLIENT",
            )
        proposal = client.analyze(profile, config=llm)
        local, recommendations, llm_warnings = mapping_from_llm_proposal(
            proposal,
            fallback=local,
        )
        origin = "llm"
        llm_model = llm.model
    columns = tuple(dict.fromkeys(key for record in samples for key in _flatten_keys(record)))
    warnings = [*sample_warning, *llm_warnings]
    if not local.label_fields:
        warnings.append(
            "No target fields were identified; set DataMapping.label_fields explicitly."
        )
    return DataInspection(
        source_path=str(path),
        detected_format=detected,
        source_size_bytes=size,
        columns=columns,
        estimated_records=estimated_records,
        mapping=local,
        mapping_origin=cast(Any, origin),
        recommendations=recommendations,
        warnings=tuple(warnings),
        llm_model=llm_model,
        shared_sample_records=shared_rows,
        shared_sample_characters=shared_characters,
    )


def _transform_sequence(value: object, mapping: DataMapping, *, record_index: int) -> str:
    if not isinstance(value, str):
        raise InputFormatError(
            "The mapped sequence field must contain text.",
            code="INVALID_CUSTOM_SEQUENCE",
            context={"record_index": record_index, "field": mapping.sequence_field},
        )
    sequence = value
    for transform in mapping.sequence_transforms:
        if transform == "strip":
            sequence = sequence.strip()
        elif transform == "uppercase":
            sequence = sequence.upper()
        elif transform == "remove_whitespace":
            sequence = "".join(sequence.split())
        elif transform == "u_to_t":
            sequence = sequence.replace("U", "T").replace("u", "t")
        elif transform == "non_iupac_to_n":
            sequence = "".join(symbol if symbol.upper() in _IUPAC else "N" for symbol in sequence)
    if not sequence:
        raise InputFormatError(
            "A mapped DNA sequence is empty after cleaning.",
            code="INVALID_CUSTOM_SEQUENCE",
            context={"record_index": record_index},
        )
    invalid = sorted(set(sequence.upper()) - _IUPAC)
    if invalid:
        raise InputFormatError(
            "A mapped sequence contains symbols outside the DNA/RNA IUPAC alphabet.",
            code="INVALID_CUSTOM_SEQUENCE",
            context={"record_index": record_index, "symbols": tuple(invalid[:50])},
            hint="Add the 'non_iupac_to_n' transform only if masking is scientifically valid.",
        )
    return sequence.upper().replace("U", "T")


def _missing(value: object, missing_values: Sequence[str]) -> bool:
    if value is None or (isinstance(value, str) and value.strip() in missing_values):
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    try:
        return bool(value != value)
    except (TypeError, ValueError):
        return False


def load_property_data(
    source: str | Path,
    *,
    inspection: DataInspection | None = None,
    mapping: DataMapping | None = None,
    llm: LLMSchemaConfig | None = None,
    llm_client: SchemaLLMClient | None = None,
    reader: DataReader | None = None,
    max_records: int = 1_000_000,
    max_file_bytes: int = 1_000_000_000,
    require_labels: bool = True,
    missing_label_policy: str = "error",
    require_sequences: bool = True,
) -> PropertyData:
    """Load mapped records and validate IDs, sequences, labels, and optional splits."""

    resolved_inspection = inspection or inspect_property_data(
        source,
        mapping=mapping,
        llm=llm,
        llm_client=llm_client,
        reader=reader,
        max_records=max_records,
        max_file_bytes=max_file_bytes,
    )
    path = _source_path(source)
    selected = resolved_inspection.mapping
    if require_sequences and selected.sequence_field is None:
        raise InputFormatError(
            "This model requires a mapped sequence field.", code="CUSTOM_DATA_FIELD_MISSING"
        )
    if reader is None:
        records = _read_builtin(
            path,
            selected,
            detected_format=resolved_inspection.detected_format,
            max_records=max_records,
            max_file_bytes=max_file_bytes,
        )
    else:
        records = _bounded_records(
            reader.read(
                path,
                mapping=selected,
                max_records=max_records,
                max_file_bytes=max_file_bytes,
            ),
            max_records=max_records,
        )
    labels = selected.label_fields
    if require_labels and not labels:
        raise InputFormatError(
            "No label fields are configured for custom-property training.",
            code="CUSTOM_DATA_LABELS_MISSING",
            hint="Set DataMapping(label_fields=(...)) or enable LLM schema analysis.",
        )
    record_ids: list[str] = []
    sequences: list[str] = []
    target_rows: list[tuple[object, ...]] = []
    splits: list[str | None] = []
    feature_rows: list[tuple[float, ...]] = []
    dropped = 0
    for index, record in enumerate(records):
        try:
            sequence_value = (
                _field(record, selected.sequence_field)
                if selected.sequence_field is not None
                else None
            )
        except KeyError as exc:
            raise InputFormatError(
                "The mapped sequence field does not exist.",
                code="CUSTOM_DATA_FIELD_MISSING",
                context={"record_index": index, "field": selected.sequence_field},
            ) from exc
        row_labels: list[object] = []
        missing_fields: list[str] = []
        for field in labels:
            try:
                value = _field(record, field)
            except KeyError:
                missing_fields.append(field)
                continue
            if _missing(value, selected.missing_values):
                missing_fields.append(field)
            row_labels.append(value)
        if missing_fields:
            if missing_label_policy == "drop":
                dropped += 1
                continue
            raise InputFormatError(
                "A custom-property record has missing target values.",
                code="CUSTOM_DATA_LABEL_MISSING",
                context={"record_index": index, "fields": tuple(missing_fields)},
                hint="Clean the source or use missing_label_policy='drop'.",
            )
        if selected.id_field is None:
            record_id = f"record_{index + 1}"
        else:
            try:
                raw_id = _field(record, selected.id_field)
            except KeyError as exc:
                raise InputFormatError(
                    "The mapped ID field does not exist.",
                    code="CUSTOM_DATA_FIELD_MISSING",
                    context={"record_index": index, "field": selected.id_field},
                ) from exc
            record_id = str(raw_id).strip()
            if not record_id:
                raise InputFormatError(
                    "A mapped record ID is empty.",
                    code="INVALID_CUSTOM_RECORD_ID",
                    context={"record_index": index},
                )
        split_value: str | None = None
        if selected.split_field is not None:
            try:
                raw_split = _field(record, selected.split_field)
            except KeyError as exc:
                raise InputFormatError(
                    "The mapped split field does not exist.",
                    code="CUSTOM_DATA_FIELD_MISSING",
                    context={"record_index": index, "field": selected.split_field},
                ) from exc
            split_value = str(raw_split).strip().lower()
        record_ids.append(record_id)
        sequences.append(
            _transform_sequence(sequence_value, selected, record_index=index)
            if selected.sequence_field is not None
            else ""
        )
        numeric: list[float] = []
        for feature in selected.feature_fields:
            try:
                raw_feature = _field(record, feature)
                if _missing(raw_feature, selected.missing_values):
                    raise ValueError("missing feature")
                number = float(cast(Any, raw_feature))
                if not math.isfinite(number):
                    raise ValueError("non-finite feature")
            except (KeyError, TypeError, ValueError, OverflowError) as exc:
                raise InputFormatError(
                    "Numeric feature fields must exist and contain finite numbers.",
                    code="INVALID_CUSTOM_FEATURE",
                    context={"record_index": index, "field": feature},
                ) from exc
            numeric.append(number)
        feature_rows.append(tuple(numeric))
        target_rows.append(tuple(row_labels))
        splits.append(split_value)
    if not record_ids:
        raise InputFormatError(
            "No usable records remain after data validation.",
            code="CUSTOM_DATA_TOO_SMALL",
            context={"usable_records": len(record_ids), "dropped_records": dropped},
        )
    duplicates = sorted(value for value, count in Counter(record_ids).items() if count > 1)
    if duplicates:
        raise InputFormatError(
            "Mapped record IDs must be unique.",
            code="DUPLICATE_CUSTOM_RECORD_ID",
            context={"examples": tuple(duplicates[:20])},
        )
    return PropertyData(
        tuple(record_ids),
        tuple(sequences),
        tuple(target_rows),
        tuple(labels),
        tuple(splits),
        resolved_inspection,
        dropped,
        tuple(feature_rows),
    )


__all__ = [
    "DataReader",
    "PropertyData",
    "detect_data_format",
    "inspect_property_data",
    "load_property_data",
]
