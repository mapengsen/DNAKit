"""Bounded chunk iteration and persistent random access for plain FASTA."""

from __future__ import annotations

import hashlib
import io
import json
import os
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal, TypeVar, cast

from dnakit.core._json import to_json_compatible
from dnakit.core.enums import Topology
from dnakit.core.record import DNARecord
from dnakit.core.sequence import IUPAC_SYMBOLS, DNASequence
from dnakit.exceptions import ConfigurationError, DuplicateIDError, InputFormatError

from ._advanced_common import write_text_path
from .api import read_one

T = TypeVar("T")
_INDEX_SCHEMA = "dnakit.fasta-index.v1"
_FASTQ_INDEX_SCHEMA = "dnakit.fastq-index.v1"


def iter_chunks(values: Iterable[T], *, chunk_size: int) -> Iterator[tuple[T, ...]]:
    """Yield stable input-order chunks without materializing the input."""

    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int):
        raise TypeError("chunk_size must be an integer.")
    if chunk_size < 1:
        raise ConfigurationError(
            "chunk_size must be positive.",
            code="INVALID_CHUNK_SIZE",
            context={"chunk_size": chunk_size},
        )
    chunk: list[T] = []
    for value in values:
        chunk.append(value)
        if len(chunk) == chunk_size:
            yield tuple(chunk)
            chunk.clear()
    if chunk:
        yield tuple(chunk)


@dataclass(frozen=True, slots=True)
class FastaIndexEntry:
    """Byte range and decoded symbol count for one FASTA record."""

    record_id: str
    byte_start: int
    byte_end: int
    sequence_length: int

    def __post_init__(self) -> None:
        if (
            not isinstance(self.record_id, str)
            or not self.record_id
            or any(symbol.isspace() for symbol in self.record_id)
        ):
            raise ConfigurationError(
                "FASTA index record_id must be non-empty.", code="INVALID_FASTA_INDEX"
            )
        values = (self.byte_start, self.byte_end, self.sequence_length)
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values
        ):
            raise ConfigurationError(
                "FASTA index offsets and lengths must be non-negative integers.",
                code="INVALID_FASTA_INDEX",
            )
        if self.byte_end <= self.byte_start:
            raise ConfigurationError(
                "FASTA index byte range must be non-empty.", code="INVALID_FASTA_INDEX"
            )


@dataclass(frozen=True, slots=True)
class FastaIndex:
    """Verified immutable index for an uncompressed FASTA path."""

    source_path: str
    source_sha256: str
    source_size: int
    source_mtime_ns: int
    entries: tuple[FastaIndexEntry, ...]
    schema_version: str = _INDEX_SCHEMA
    _by_id: Mapping[str, FastaIndexEntry] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.schema_version != _INDEX_SCHEMA:
            raise ConfigurationError(
                "Unsupported FASTA index schema.", code="UNSUPPORTED_FASTA_INDEX_SCHEMA"
            )
        if not isinstance(self.source_path, str) or not self.source_path:
            raise ConfigurationError(
                "FASTA index source_path must be non-empty.", code="INVALID_FASTA_INDEX"
            )
        if not re_full_sha256(self.source_sha256):
            raise ConfigurationError(
                "FASTA index checksum must be lowercase SHA-256.", code="INVALID_FASTA_INDEX"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in (self.source_size, self.source_mtime_ns)
        ):
            raise ConfigurationError(
                "FASTA index source size and mtime must be non-negative integers.",
                code="INVALID_FASTA_INDEX",
            )
        if not isinstance(self.entries, tuple) or any(
            not isinstance(entry, FastaIndexEntry) for entry in self.entries
        ):
            raise ConfigurationError(
                "FASTA index entries have invalid types.", code="INVALID_FASTA_INDEX"
            )
        lookup: dict[str, FastaIndexEntry] = {}
        previous_end = 0
        for entry in self.entries:
            if entry.record_id in lookup:
                raise DuplicateIDError(
                    "FASTA index contains duplicate record IDs.",
                    context={"record_id": entry.record_id},
                )
            if entry.byte_start < previous_end or entry.byte_end > self.source_size:
                raise ConfigurationError(
                    "FASTA index byte ranges are invalid or overlap.", code="INVALID_FASTA_INDEX"
                )
            lookup[entry.record_id] = entry
            previous_end = entry.byte_end
        object.__setattr__(self, "_by_id", MappingProxyType(lookup))

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(entry.record_id for entry in self.entries)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_path": self.source_path,
            "source_sha256": self.source_sha256,
            "source_size": self.source_size,
            "source_mtime_ns": self.source_mtime_ns,
            "entries": [cast(dict[str, Any], to_json_compatible(entry)) for entry in self.entries],
        }

    def fetch(
        self,
        record_id: str,
        *,
        start: int | None = None,
        end: int | None = None,
        strand: Literal["+", "-", "forward", "reverse"] = "+",
        max_record_bytes: int = 100_000_000,
    ) -> DNARecord:
        """Fetch one record or one 0-based half-open subsequence by ID."""

        if not isinstance(record_id, str) or not record_id:
            raise TypeError("record_id must be non-empty text.")
        try:
            entry = self._by_id[record_id]
        except KeyError as exc:
            raise KeyError(record_id) from exc
        if (
            isinstance(max_record_bytes, bool)
            or not isinstance(max_record_bytes, int)
            or max_record_bytes < 1
        ):
            raise ConfigurationError(
                "max_record_bytes must be positive.", code="INVALID_FASTA_INDEX_LIMIT"
            )
        byte_count = entry.byte_end - entry.byte_start
        if byte_count > max_record_bytes:
            raise InputFormatError(
                "Indexed FASTA record exceeds max_record_bytes.",
                code="FASTA_RECORD_LIMIT_EXCEEDED",
                context={
                    "record_id": record_id,
                    "byte_count": byte_count,
                    "max_record_bytes": max_record_bytes,
                },
            )
        path = Path(self.source_path)
        try:
            stat = path.stat()
            if stat.st_size != self.source_size or stat.st_mtime_ns != self.source_mtime_ns:
                raise InputFormatError(
                    "Indexed FASTA changed after index verification.",
                    code="STALE_FASTA_INDEX",
                    context={"source": str(path)},
                )
            with path.open("rb") as handle:
                handle.seek(entry.byte_start)
                payload = handle.read(byte_count)
            final_stat = path.stat()
            if final_stat.st_size != stat.st_size or final_stat.st_mtime_ns != stat.st_mtime_ns:
                raise InputFormatError(
                    "Indexed FASTA changed while a record was being fetched.",
                    code="STALE_FASTA_INDEX",
                    context={"source": str(path)},
                )
        except OSError as exc:
            raise InputFormatError(
                "Could not read indexed FASTA source.",
                code="INPUT_READ_FAILED",
                context={"source": str(path), "reason": str(exc)},
            ) from exc
        if len(payload) != byte_count:
            raise InputFormatError(
                "Indexed FASTA byte range is truncated.",
                code="STALE_FASTA_INDEX",
                context={"record_id": record_id},
            )
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise InputFormatError(
                "Indexed FASTA record is not UTF-8.",
                code="INPUT_DECODE_FAILED",
                context={"record_id": record_id},
            ) from exc
        record = read_one(io.StringIO(text), format="fasta")
        if record.id != record_id or record.sequence.symbol_length != entry.sequence_length:
            raise InputFormatError(
                "Indexed FASTA content no longer matches the index.",
                code="STALE_FASTA_INDEX",
                context={"record_id": record_id},
            )
        if start is None and end is None and strand in {"+", "forward"}:
            return record
        resolved_start = 0 if start is None else start
        resolved_end = entry.sequence_length if end is None else end
        for name, value in (("start", resolved_start), ("end", resolved_end)):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{name} must be an integer or None.")
        if (
            resolved_start < 0
            or resolved_end < resolved_start
            or resolved_end > entry.sequence_length
        ):
            raise InputFormatError(
                "Indexed FASTA query coordinates are outside the record.",
                code="FASTA_QUERY_OUT_OF_BOUNDS",
                context={
                    "record_id": record_id,
                    "start": resolved_start,
                    "end": resolved_end,
                    "sequence_length": entry.sequence_length,
                },
            )
        if strand not in {"+", "-", "forward", "reverse"}:
            raise ConfigurationError(
                "strand must be '+', '-', 'forward', or 'reverse'.",
                code="INVALID_FASTA_QUERY_STRAND",
            )
        selected = record.sequence.to_string()[resolved_start:resolved_end]
        sequence = DNASequence(
            selected,
            alphabet=record.sequence.alphabet,
            topology=Topology.LINEAR,
            strandedness=record.sequence.strandedness,
        )
        if strand in {"-", "reverse"}:
            sequence = sequence.reverse_complement()
        return DNARecord(
            sequence,
            record.id,
            description=record.description,
            metadata={
                "fasta_index_query": {
                    "start": resolved_start,
                    "end": resolved_end,
                    "strand": "reverse" if strand in {"-", "reverse"} else "forward",
                    "coordinate_system": "0-based-half-open",
                }
            },
        )


def re_full_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(symbol in "0123456789abcdef" for symbol in value)
    )


def _digest_path(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_non_negative_integer(value: object, *, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field} is not a non-negative integer")
    return value


def _json_string(value: object, *, field: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field} is not a string")
    return value


def _json_object_without_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _read_bounded_index_text(path: Path, *, max_index_bytes: int, limit_error_code: str) -> str:
    try:
        with path.open("rb") as handle:
            payload = handle.read(max_index_bytes + 1)
    except OSError as exc:
        raise InputFormatError(
            "Could not read index input.",
            code="INPUT_READ_FAILED",
            context={"source": str(path), "reason": str(exc)},
        ) from exc
    if len(payload) > max_index_bytes:
        raise InputFormatError(
            "Index input exceeds max_index_bytes.",
            code=limit_error_code,
            context={"max_index_bytes": max_index_bytes},
        )
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "Index input is not UTF-8.",
            code="INPUT_DECODE_FAILED",
            context={"source": str(path)},
        ) from exc


def build_fasta_index(
    source_path: str | os.PathLike[str],
    index_path: str | os.PathLike[str] | None = None,
    *,
    overwrite: bool = False,
    max_records: int = 10_000_000,
    max_line_length: int = 10_000_000,
) -> FastaIndex:
    """Scan and persist an ID index for an uncompressed UTF-8 FASTA file."""

    for name, value in (("max_records", max_records), ("max_line_length", max_line_length)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ConfigurationError(f"{name} must be positive.", code="INVALID_FASTA_INDEX_LIMIT")
    source = Path(source_path).resolve()
    if not source.is_file():
        raise InputFormatError(
            "FASTA index source must be an existing file.",
            code="INPUT_OPEN_FAILED",
            context={"source": str(source)},
        )
    initial_stat = source.stat()
    with source.open("rb") as probe:
        if probe.read(2) == b"\x1f\x8b":
            raise InputFormatError(
                "Persistent FASTA indexing supports only uncompressed files.",
                code="COMPRESSED_FASTA_INDEX_UNSUPPORTED",
            )
    entries: list[FastaIndexEntry] = []
    seen: set[str] = set()
    current_id: str | None = None
    current_start = 0
    current_length = 0
    offset = 0
    with source.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if len(raw.rstrip(b"\r\n")) > max_line_length:
                raise InputFormatError(
                    "FASTA line exceeds max_line_length.",
                    code="FASTA_LINE_TOO_LONG",
                    context={"line_number": line_number, "max_line_length": max_line_length},
                )
            if raw.startswith(b">"):
                if current_id is not None:
                    entries.append(
                        FastaIndexEntry(current_id, current_start, offset, current_length)
                    )
                if len(entries) >= max_records:
                    raise InputFormatError(
                        "FASTA input exceeds max_records.",
                        code="FASTA_INDEX_RECORD_LIMIT_EXCEEDED",
                        context={"max_records": max_records},
                    )
                try:
                    header = raw[1:].decode("utf-8").strip()
                except UnicodeDecodeError as exc:
                    raise InputFormatError(
                        "FASTA header is not UTF-8.",
                        code="INPUT_DECODE_FAILED",
                        context={"line_number": line_number},
                    ) from exc
                if not header:
                    raise InputFormatError(
                        "FASTA header must contain an ID.",
                        code="EMPTY_FASTA_HEADER",
                        context={"line_number": line_number},
                    )
                current_id = header.split(maxsplit=1)[0]
                if current_id in seen:
                    raise DuplicateIDError(
                        "FASTA indexing requires unique record IDs.",
                        context={"record_id": current_id, "line_number": line_number},
                    )
                seen.add(current_id)
                current_start = offset
                current_length = 0
            else:
                body = raw.rstrip(b"\r\n")
                if not body:
                    offset += len(raw)
                    continue
                if current_id is None:
                    raise InputFormatError(
                        "FASTA sequence data appears before the first header.",
                        code="FASTA_MISSING_HEADER",
                        context={"line_number": line_number},
                    )
                try:
                    symbols = body.decode("ascii").upper()
                except UnicodeDecodeError as exc:
                    raise InputFormatError(
                        "FASTA sequence is not ASCII DNA text.",
                        code="INVALID_SEQUENCE_CONTENT",
                        context={"line_number": line_number},
                    ) from exc
                if any(symbol not in IUPAC_SYMBOLS for symbol in symbols):
                    raise InputFormatError(
                        "FASTA sequence contains invalid DNA symbols.",
                        code="INVALID_SEQUENCE_CONTENT",
                        context={"line_number": line_number},
                    )
                current_length += len(symbols)
            offset += len(raw)
    if current_id is not None:
        entries.append(FastaIndexEntry(current_id, current_start, offset, current_length))
    digest = _digest_path(source)
    final_stat = source.stat()
    if (
        initial_stat.st_size != final_stat.st_size
        or initial_stat.st_mtime_ns != final_stat.st_mtime_ns
    ):
        raise InputFormatError(
            "FASTA source changed while its index was being built.",
            code="FASTA_CHANGED_DURING_INDEX",
            context={"source": str(source)},
        )
    index = FastaIndex(
        str(source),
        digest,
        final_stat.st_size,
        final_stat.st_mtime_ns,
        tuple(entries),
    )
    destination = (
        Path(index_path) if index_path is not None else Path(str(source) + ".dnakit.fai.json")
    )

    def writer(handle: Any) -> int:
        json.dump(
            index.to_dict(), handle, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        handle.write("\n")
        return len(entries)

    write_text_path(destination, writer, overwrite=overwrite, create_parents=False)
    return index


def load_fasta_index(
    index_path: str | os.PathLike[str],
    *,
    source_path: str | os.PathLike[str] | None = None,
    verify_checksum: bool = True,
    max_index_bytes: int = 64 * 1024 * 1024,
) -> FastaIndex:
    """Load and optionally checksum-verify a persistent DNAKit FASTA index."""

    if not isinstance(verify_checksum, bool):
        raise TypeError("verify_checksum must be a boolean.")
    path = Path(index_path)
    if (
        isinstance(max_index_bytes, bool)
        or not isinstance(max_index_bytes, int)
        or max_index_bytes < 1
    ):
        raise ConfigurationError(
            "max_index_bytes must be positive.", code="INVALID_FASTA_INDEX_LIMIT"
        )
    try:
        payload = json.loads(
            _read_bounded_index_text(
                path,
                max_index_bytes=max_index_bytes,
                limit_error_code="FASTA_INDEX_FILE_LIMIT_EXCEEDED",
            ),
            object_pairs_hook=_json_object_without_duplicates,
        )
    except InputFormatError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError, ValueError) as exc:
        raise InputFormatError(
            "Could not decode FASTA index JSON.",
            code="INVALID_FASTA_INDEX",
            context={"source": str(path)},
        ) from exc
    if not isinstance(payload, Mapping) or not isinstance(payload.get("entries"), list):
        raise InputFormatError("FASTA index root/schema is invalid.", code="INVALID_FASTA_INDEX")
    try:
        entries = tuple(
            FastaIndexEntry(
                _json_string(item["record_id"], field="record_id"),
                _json_non_negative_integer(item["byte_start"], field="byte_start"),
                _json_non_negative_integer(item["byte_end"], field="byte_end"),
                _json_non_negative_integer(item["sequence_length"], field="sequence_length"),
            )
            for item in payload["entries"]
            if isinstance(item, Mapping)
        )
        if len(entries) != len(payload["entries"]):
            raise ValueError("non-object entry")
        stored_source = _json_string(payload["source_path"], field="source_path")
        resolved_source = (
            str(Path(source_path).resolve()) if source_path is not None else stored_source
        )
        index = FastaIndex(
            resolved_source,
            _json_string(payload["source_sha256"], field="source_sha256"),
            _json_non_negative_integer(payload["source_size"], field="source_size"),
            _json_non_negative_integer(payload["source_mtime_ns"], field="source_mtime_ns"),
            entries,
            _json_string(payload["schema_version"], field="schema_version"),
        )
    except (KeyError, TypeError, ValueError, ConfigurationError, DuplicateIDError) as exc:
        raise InputFormatError(
            "FASTA index fields are invalid.", code="INVALID_FASTA_INDEX"
        ) from exc
    source = Path(index.source_path)
    try:
        stat = source.stat()
        if stat.st_size != index.source_size:
            raise InputFormatError(
                "FASTA source size no longer matches its index.",
                code="STALE_FASTA_INDEX",
                context={"source": str(source)},
            )
        if not verify_checksum and stat.st_mtime_ns != index.source_mtime_ns:
            raise InputFormatError(
                "FASTA source mtime no longer matches its index.",
                code="STALE_FASTA_INDEX",
                context={"source": str(source)},
            )
        if verify_checksum:
            if _digest_path(source) != index.source_sha256:
                raise InputFormatError(
                    "FASTA source checksum no longer matches its index.",
                    code="STALE_FASTA_INDEX",
                    context={"source": str(source)},
                )
            if stat.st_mtime_ns != index.source_mtime_ns:
                index = FastaIndex(
                    index.source_path,
                    index.source_sha256,
                    index.source_size,
                    stat.st_mtime_ns,
                    index.entries,
                    index.schema_version,
                )
    except InputFormatError:
        raise
    except OSError as exc:
        raise InputFormatError(
            "Indexed FASTA source is unavailable.",
            code="INPUT_OPEN_FAILED",
            context={"source": str(source)},
        ) from exc
    return index


@dataclass(frozen=True, slots=True)
class FastqIndexEntry:
    """Byte range and symbol count for one strict four-line FASTQ record."""

    record_id: str
    byte_start: int
    byte_end: int
    sequence_length: int

    def __post_init__(self) -> None:
        if (
            not isinstance(self.record_id, str)
            or not self.record_id
            or any(symbol.isspace() for symbol in self.record_id)
        ):
            raise ConfigurationError(
                "FASTQ index record_id must be non-empty without whitespace.",
                code="INVALID_FASTQ_INDEX",
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in (self.byte_start, self.byte_end, self.sequence_length)
        ):
            raise ConfigurationError(
                "FASTQ index offsets and lengths must be non-negative integers.",
                code="INVALID_FASTQ_INDEX",
            )
        if self.byte_end <= self.byte_start:
            raise ConfigurationError(
                "FASTQ index byte range must be non-empty.", code="INVALID_FASTQ_INDEX"
            )


@dataclass(frozen=True, slots=True)
class FastqIndex:
    """Verified immutable index for an uncompressed four-line FASTQ path."""

    source_path: str
    source_sha256: str
    source_size: int
    source_mtime_ns: int
    entries: tuple[FastqIndexEntry, ...]
    phred_offset: int = 33
    schema_version: str = _FASTQ_INDEX_SCHEMA
    _by_id: Mapping[str, FastqIndexEntry] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.schema_version != _FASTQ_INDEX_SCHEMA:
            raise ConfigurationError(
                "Unsupported FASTQ index schema.", code="UNSUPPORTED_FASTQ_INDEX_SCHEMA"
            )
        if not isinstance(self.source_path, str) or not self.source_path:
            raise ConfigurationError(
                "FASTQ index source_path must be non-empty.", code="INVALID_FASTQ_INDEX"
            )
        if not re_full_sha256(self.source_sha256):
            raise ConfigurationError(
                "FASTQ index checksum must be lowercase SHA-256.", code="INVALID_FASTQ_INDEX"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in (self.source_size, self.source_mtime_ns)
        ):
            raise ConfigurationError(
                "FASTQ index source size and mtime must be non-negative integers.",
                code="INVALID_FASTQ_INDEX",
            )
        if (
            isinstance(self.phred_offset, bool)
            or not isinstance(self.phred_offset, int)
            or not 33 <= self.phred_offset <= 64
        ):
            raise ConfigurationError(
                "FASTQ index phred_offset must be between 33 and 64.",
                code="INVALID_FASTQ_INDEX",
            )
        if not isinstance(self.entries, tuple) or any(
            not isinstance(entry, FastqIndexEntry) for entry in self.entries
        ):
            raise ConfigurationError(
                "FASTQ index entries have invalid types.", code="INVALID_FASTQ_INDEX"
            )
        lookup: dict[str, FastqIndexEntry] = {}
        previous_end = 0
        for entry in self.entries:
            if entry.record_id in lookup:
                raise DuplicateIDError(
                    "FASTQ index contains duplicate record IDs.",
                    context={"record_id": entry.record_id},
                )
            if entry.byte_start != previous_end or entry.byte_end > self.source_size:
                raise ConfigurationError(
                    "FASTQ index byte ranges must be valid, contiguous, and non-overlapping.",
                    code="INVALID_FASTQ_INDEX",
                )
            lookup[entry.record_id] = entry
            previous_end = entry.byte_end
        if previous_end != self.source_size:
            raise ConfigurationError(
                "FASTQ index entries must cover the complete source.",
                code="INVALID_FASTQ_INDEX",
            )
        object.__setattr__(self, "_by_id", MappingProxyType(lookup))

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(entry.record_id for entry in self.entries)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_path": self.source_path,
            "source_sha256": self.source_sha256,
            "source_size": self.source_size,
            "source_mtime_ns": self.source_mtime_ns,
            "phred_offset": self.phred_offset,
            "entries": [cast(dict[str, Any], to_json_compatible(entry)) for entry in self.entries],
        }

    def fetch(
        self,
        record_id: str,
        *,
        start: int | None = None,
        end: int | None = None,
        strand: Literal["+", "-", "forward", "reverse"] = "+",
        max_record_bytes: int = 100_000_000,
    ) -> DNARecord:
        """Fetch a FASTQ record or quality-synchronized 0-based subsequence."""

        if not isinstance(record_id, str) or not record_id:
            raise TypeError("record_id must be non-empty text.")
        try:
            entry = self._by_id[record_id]
        except KeyError as exc:
            raise KeyError(record_id) from exc
        if (
            isinstance(max_record_bytes, bool)
            or not isinstance(max_record_bytes, int)
            or max_record_bytes < 1
        ):
            raise ConfigurationError(
                "max_record_bytes must be positive.", code="INVALID_FASTQ_INDEX_LIMIT"
            )
        byte_count = entry.byte_end - entry.byte_start
        if byte_count > max_record_bytes:
            raise InputFormatError(
                "Indexed FASTQ record exceeds max_record_bytes.",
                code="FASTQ_RECORD_LIMIT_EXCEEDED",
                context={"record_id": record_id, "max_record_bytes": max_record_bytes},
            )
        path = Path(self.source_path)
        try:
            stat = path.stat()
            if stat.st_size != self.source_size or stat.st_mtime_ns != self.source_mtime_ns:
                raise InputFormatError(
                    "Indexed FASTQ changed after index verification.",
                    code="STALE_FASTQ_INDEX",
                    context={"source": str(path)},
                )
            with path.open("rb") as handle:
                handle.seek(entry.byte_start)
                payload = handle.read(byte_count)
            final_stat = path.stat()
            if final_stat.st_size != stat.st_size or final_stat.st_mtime_ns != stat.st_mtime_ns:
                raise InputFormatError(
                    "Indexed FASTQ changed while a record was being fetched.",
                    code="STALE_FASTQ_INDEX",
                    context={"source": str(path)},
                )
        except InputFormatError:
            raise
        except OSError as exc:
            raise InputFormatError(
                "Could not read indexed FASTQ source.",
                code="INPUT_READ_FAILED",
                context={"source": str(path), "reason": str(exc)},
            ) from exc
        if len(payload) != byte_count:
            raise InputFormatError(
                "Indexed FASTQ byte range is truncated.",
                code="STALE_FASTQ_INDEX",
                context={"record_id": record_id},
            )
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise InputFormatError(
                "Indexed FASTQ record is not UTF-8.",
                code="INPUT_DECODE_FAILED",
                context={"record_id": record_id},
            ) from exc
        record = read_one(
            io.StringIO(text),
            format="fastq",
            config=_fastq_read_config(self.phred_offset, entry.sequence_length, byte_count),
        )
        if record.id != record_id or record.sequence.symbol_length != entry.sequence_length:
            raise InputFormatError(
                "Indexed FASTQ content no longer matches the index.",
                code="STALE_FASTQ_INDEX",
                context={"record_id": record_id},
            )
        if start is None and end is None and strand in {"+", "forward"}:
            return record
        resolved_start, resolved_end = _resolve_fastq_query(
            record_id=record_id,
            sequence_length=entry.sequence_length,
            start=start,
            end=end,
            strand=strand,
        )
        selected = record.sequence.to_string()[resolved_start:resolved_end]
        sequence = DNASequence(
            selected,
            alphabet=record.sequence.alphabet,
            topology=Topology.LINEAR,
            strandedness=record.sequence.strandedness,
        )
        qualities = tuple(record.letter_annotations["phred_quality"])[resolved_start:resolved_end]
        if strand in {"-", "reverse"}:
            sequence = sequence.reverse_complement()
            qualities = tuple(reversed(qualities))
        return DNARecord(
            sequence,
            record.id,
            description=record.description,
            metadata={
                "fastq_index_query": {
                    "start": resolved_start,
                    "end": resolved_end,
                    "strand": "reverse" if strand in {"-", "reverse"} else "forward",
                    "coordinate_system": "0-based-half-open",
                }
            },
            letter_annotations={"phred_quality": qualities},
        )


def _fastq_read_config(phred_offset: int, sequence_length: int, byte_count: int) -> Any:
    from .config import ReadConfig

    return ReadConfig(
        phred_offset=phred_offset,
        max_sequence_symbols=max(1, sequence_length),
        max_input_bytes=max(1, byte_count),
    )


def _resolve_fastq_query(
    *,
    record_id: str,
    sequence_length: int,
    start: int | None,
    end: int | None,
    strand: str,
) -> tuple[int, int]:
    resolved_start = 0 if start is None else start
    resolved_end = sequence_length if end is None else end
    for name, value in (("start", resolved_start), ("end", resolved_end)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer or None.")
    if resolved_start < 0 or resolved_end < resolved_start or resolved_end > sequence_length:
        raise InputFormatError(
            "Indexed FASTQ query coordinates are outside the record.",
            code="FASTQ_QUERY_OUT_OF_BOUNDS",
            context={
                "record_id": record_id,
                "start": resolved_start,
                "end": resolved_end,
                "sequence_length": sequence_length,
            },
        )
    if strand not in {"+", "-", "forward", "reverse"}:
        raise ConfigurationError(
            "strand must be '+', '-', 'forward', or 'reverse'.",
            code="INVALID_FASTQ_QUERY_STRAND",
        )
    return resolved_start, resolved_end


def build_fastq_index(
    source_path: str | os.PathLike[str],
    index_path: str | os.PathLike[str] | None = None,
    *,
    overwrite: bool = False,
    phred_offset: int = 33,
    max_records: int = 10_000_000,
    max_line_length: int = 10_000_000,
    max_record_bytes: int = 100_000_000,
    max_source_bytes: int = 1_000_000_000,
) -> FastqIndex:
    """Scan and persist an ID index for strict uncompressed four-line FASTQ."""

    for name, value in (
        ("max_records", max_records),
        ("max_line_length", max_line_length),
        ("max_record_bytes", max_record_bytes),
        ("max_source_bytes", max_source_bytes),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ConfigurationError(f"{name} must be positive.", code="INVALID_FASTQ_INDEX_LIMIT")
    if (
        isinstance(phred_offset, bool)
        or not isinstance(phred_offset, int)
        or not 33 <= phred_offset <= 64
    ):
        raise ConfigurationError(
            "phred_offset must be between 33 and 64.", code="INVALID_FASTQ_INDEX_LIMIT"
        )
    source = Path(source_path).resolve()
    if not source.is_file():
        raise InputFormatError(
            "FASTQ index source must be an existing file.",
            code="INPUT_OPEN_FAILED",
            context={"source": str(source)},
        )
    initial_stat = source.stat()
    if initial_stat.st_size > max_source_bytes:
        raise InputFormatError(
            "FASTQ source exceeds max_source_bytes.",
            code="FASTQ_SOURCE_LIMIT_EXCEEDED",
            context={"max_source_bytes": max_source_bytes},
        )
    entries: list[FastqIndexEntry] = []
    seen: set[str] = set()
    with source.open("rb") as handle:
        if handle.read(2) == b"\x1f\x8b":
            raise InputFormatError(
                "Persistent FASTQ indexing supports only uncompressed files.",
                code="COMPRESSED_FASTQ_INDEX_UNSUPPORTED",
            )
        handle.seek(0)
        while True:
            start = handle.tell()
            header = handle.readline()
            if not header:
                break
            if len(entries) >= max_records:
                raise InputFormatError(
                    "FASTQ input exceeds max_records.",
                    code="FASTQ_INDEX_RECORD_LIMIT_EXCEEDED",
                    context={"max_records": max_records},
                )
            lines = (header, handle.readline(), handle.readline(), handle.readline())
            if any(not line for line in lines):
                raise InputFormatError(
                    "FASTQ record ended before all four lines were present.",
                    code="TRUNCATED_FASTQ_RECORD",
                    context={"byte_start": start},
                )
            for offset, line in enumerate(lines):
                if len(line.rstrip(b"\r\n")) > max_line_length:
                    raise InputFormatError(
                        "FASTQ line exceeds max_line_length.",
                        code="FASTQ_LINE_TOO_LONG",
                        context={"byte_start": start, "record_line": offset + 1},
                    )
            end = handle.tell()
            if end - start > max_record_bytes:
                raise InputFormatError(
                    "FASTQ record exceeds max_record_bytes.",
                    code="FASTQ_RECORD_LIMIT_EXCEEDED",
                    context={"byte_start": start, "max_record_bytes": max_record_bytes},
                )
            payload = b"".join(lines)
            try:
                text = payload.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise InputFormatError(
                    "FASTQ index source is not UTF-8.", code="INPUT_DECODE_FAILED"
                ) from exc
            record = read_one(
                io.StringIO(text),
                format="fastq",
                config=_fastq_read_config(phred_offset, max_line_length, len(payload)),
            )
            if record.id in seen:
                raise DuplicateIDError(
                    "FASTQ indexing requires unique record IDs.",
                    context={"record_id": record.id},
                )
            seen.add(record.id)
            entries.append(FastqIndexEntry(record.id, start, end, record.sequence.symbol_length))
    digest = _digest_path(source)
    final_stat = source.stat()
    if (
        initial_stat.st_size != final_stat.st_size
        or initial_stat.st_mtime_ns != final_stat.st_mtime_ns
    ):
        raise InputFormatError(
            "FASTQ source changed while its index was being built.",
            code="FASTQ_CHANGED_DURING_INDEX",
            context={"source": str(source)},
        )
    index = FastqIndex(
        str(source),
        digest,
        final_stat.st_size,
        final_stat.st_mtime_ns,
        tuple(entries),
        phred_offset,
    )
    destination = (
        Path(index_path) if index_path is not None else Path(str(source) + ".dnakit.fqi.json")
    )

    def writer(handle: Any) -> int:
        json.dump(
            index.to_dict(), handle, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        handle.write("\n")
        return len(entries)

    write_text_path(destination, writer, overwrite=overwrite, create_parents=False)
    return index


def load_fastq_index(
    index_path: str | os.PathLike[str],
    *,
    source_path: str | os.PathLike[str] | None = None,
    verify_checksum: bool = True,
    max_index_bytes: int = 64 * 1024 * 1024,
    max_entries: int = 10_000_000,
) -> FastqIndex:
    """Load and checksum-verify a persistent DNAKit FASTQ index."""

    if not isinstance(verify_checksum, bool):
        raise TypeError("verify_checksum must be a boolean.")
    for name, value in (("max_index_bytes", max_index_bytes), ("max_entries", max_entries)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ConfigurationError(f"{name} must be positive.", code="INVALID_FASTQ_INDEX_LIMIT")
    path = Path(index_path)
    try:
        payload = json.loads(
            _read_bounded_index_text(
                path,
                max_index_bytes=max_index_bytes,
                limit_error_code="FASTQ_INDEX_FILE_LIMIT_EXCEEDED",
            ),
            object_pairs_hook=_json_object_without_duplicates,
        )
    except InputFormatError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError, ValueError) as exc:
        raise InputFormatError(
            "Could not decode FASTQ index JSON.",
            code="INVALID_FASTQ_INDEX",
            context={"source": str(path)},
        ) from exc
    required = {
        "schema_version",
        "source_path",
        "source_sha256",
        "source_size",
        "source_mtime_ns",
        "phred_offset",
        "entries",
    }
    if (
        not isinstance(payload, Mapping)
        or set(payload) != required
        or not isinstance(payload.get("entries"), list)
    ):
        raise InputFormatError("FASTQ index root/schema is invalid.", code="INVALID_FASTQ_INDEX")
    raw_entries = cast(list[object], payload["entries"])
    if len(raw_entries) > max_entries:
        raise InputFormatError(
            "FASTQ index exceeds max_entries.",
            code="FASTQ_INDEX_ENTRY_LIMIT_EXCEEDED",
            context={"max_entries": max_entries},
        )
    try:
        entries: list[FastqIndexEntry] = []
        for item in raw_entries:
            if not isinstance(item, Mapping) or set(item) != {
                "record_id",
                "byte_start",
                "byte_end",
                "sequence_length",
            }:
                raise ValueError("invalid FASTQ index entry schema")
            entries.append(
                FastqIndexEntry(
                    _json_string(item["record_id"], field="record_id"),
                    _json_non_negative_integer(item["byte_start"], field="byte_start"),
                    _json_non_negative_integer(item["byte_end"], field="byte_end"),
                    _json_non_negative_integer(item["sequence_length"], field="sequence_length"),
                )
            )
        stored_source = _json_string(payload["source_path"], field="source_path")
        resolved_source = (
            str(Path(source_path).resolve()) if source_path is not None else stored_source
        )
        index = FastqIndex(
            resolved_source,
            _json_string(payload["source_sha256"], field="source_sha256"),
            _json_non_negative_integer(payload["source_size"], field="source_size"),
            _json_non_negative_integer(payload["source_mtime_ns"], field="source_mtime_ns"),
            tuple(entries),
            _json_non_negative_integer(payload["phred_offset"], field="phred_offset"),
            _json_string(payload["schema_version"], field="schema_version"),
        )
    except (KeyError, TypeError, ValueError, ConfigurationError, DuplicateIDError) as exc:
        raise InputFormatError(
            "FASTQ index fields are invalid.", code="INVALID_FASTQ_INDEX"
        ) from exc
    source = Path(index.source_path)
    try:
        stat = source.stat()
        if stat.st_size != index.source_size:
            raise InputFormatError(
                "FASTQ source size no longer matches its index.", code="STALE_FASTQ_INDEX"
            )
        if stat.st_mtime_ns != index.source_mtime_ns:
            raise InputFormatError(
                "FASTQ source mtime no longer matches its index.", code="STALE_FASTQ_INDEX"
            )
        if verify_checksum:
            if _digest_path(source) != index.source_sha256:
                raise InputFormatError(
                    "FASTQ source checksum no longer matches its index.",
                    code="STALE_FASTQ_INDEX",
                )
            final_stat = source.stat()
            if final_stat.st_size != stat.st_size or final_stat.st_mtime_ns != stat.st_mtime_ns:
                raise InputFormatError(
                    "FASTQ source changed during checksum verification.",
                    code="STALE_FASTQ_INDEX",
                )
    except InputFormatError:
        raise
    except OSError as exc:
        raise InputFormatError(
            "Indexed FASTQ source is unavailable.",
            code="INPUT_OPEN_FAILED",
            context={"source": str(source)},
        ) from exc
    return index


__all__ = [
    "FastaIndex",
    "FastaIndexEntry",
    "FastqIndex",
    "FastqIndexEntry",
    "build_fasta_index",
    "build_fastq_index",
    "iter_chunks",
    "load_fasta_index",
    "load_fastq_index",
]
