"""Public file and stream I/O entry points."""

from __future__ import annotations

import codecs
import gzip
import hashlib
import io
import os
import tempfile
from collections.abc import Callable, Iterable, Iterator
from contextlib import suppress
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Literal, TextIO, TypeAlias, cast, overload

from dnakit.core._json import to_json_compatible
from dnakit.core.collection import DNASet
from dnakit.core.enums import ExecutionMode, ImplementationLabel, OriginClass
from dnakit.core.facade import DNA
from dnakit.core.provenance import ArtifactRef, ImplementationInfo, Provenance
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, InputFormatError

from ._formats import (
    iter_delimited,
    iter_fasta,
    iter_fastq,
    iter_json_array,
    iter_json_lines,
    write_delimited,
    write_fasta,
    write_fastq,
    write_json_array,
    write_json_lines,
)
from .config import ReadConfig, WriteConfig
from .genbank import iter_genbank, write_genbank
from .results import GeneratedID, WriteResult
from .source import RecordSource

PathSource: TypeAlias = str | os.PathLike[str]
ReadableSource: TypeAlias = PathSource | TextIO | BinaryIO
WritableTarget: TypeAlias = PathSource | TextIO | BinaryIO
WritableRecords: TypeAlias = DNA | DNASequence | DNARecord | Iterable[DNASequence | DNARecord]
ReadMode: TypeAlias = Literal["stream", "dna"]

_FORMAT_ALIASES = {
    "fa": "fasta",
    "fasta": "fasta",
    "fna": "fasta",
    "ffn": "fasta",
    "faa": "fasta",
    "fq": "fastq",
    "fastq": "fastq",
    "csv": "csv",
    "tsv": "tsv",
    "tab": "tsv",
    "json": "json",
    "jsonl": "jsonl",
    "ndjson": "jsonl",
    "gb": "genbank",
    "gbk": "genbank",
    "gbff": "genbank",
    "genbank": "genbank",
}
_UNSUPPORTED_FORMATS = {
    "parquet": "Parquet support requires the optional advanced table backend.",
    "gff": "GFF support belongs to the advanced annotation codec.",
    "gff3": "GFF3 support belongs to the advanced annotation codec.",
    "bed": "BED support belongs to the advanced annotation codec.",
    "agp": "AGP support belongs to the advanced assembly codec.",
    "fai": "Indexed FASTA access belongs to the advanced indexing adapter.",
}
_MEDIA_TYPES = {
    "fasta": "text/x-fasta",
    "fastq": "text/x-fastq",
    "csv": "text/csv",
    "tsv": "text/tab-separated-values",
    "json": "application/json",
    "jsonl": "application/x-ndjson",
    "genbank": "text/x-genbank",
}


class _BoundedTextReader:
    """Proxy a decoded stream while bounding bytes returned to parsers."""

    def __init__(self, raw: TextIO, *, encoding: str, limit: int, source: str | None) -> None:
        self._raw = raw
        self._encoding = encoding
        self._encoder = codecs.getincrementalencoder(encoding)()
        self._limit = limit
        self._source = source
        self._count = 0

    def _checked(self, value: str) -> str:
        self._count += len(self._encoder.encode(value))
        if self._count > self._limit:
            raise InputFormatError(
                "Decoded input exceeds max_input_bytes.",
                code="INPUT_BYTE_LIMIT_EXCEEDED",
                context={"source": self._source, "max_input_bytes": self._limit},
            )
        return value

    def _bounded_size(self, size: int | None) -> int:
        remaining = self._limit - self._count
        if size is None or size < 0:
            return remaining + 1
        return min(size, remaining + 1)

    def read(self, size: int = -1) -> str:
        return self._checked(self._raw.read(self._bounded_size(size)))

    def readline(self, size: int = -1) -> str:
        return self._checked(self._raw.readline(self._bounded_size(size)))

    def __iter__(self) -> _BoundedTextReader:
        return self

    def __next__(self) -> str:
        value = self.readline()
        if not value:
            raise StopIteration
        return value

    def __getattr__(self, name: str) -> object:
        return getattr(self._raw, name)


class _BoundedBinaryWriter:
    """Count physical output bytes before forwarding each binary write."""

    def __init__(self, raw: BinaryIO, *, limit: int) -> None:
        self._raw = raw
        self._limit = limit
        self._count = 0

    def write(self, value: bytes | bytearray | memoryview) -> int:
        size = len(value)
        if self._count + size > self._limit:
            raise InputFormatError(
                "Serialized output exceeds max_output_bytes.",
                code="OUTPUT_BYTE_LIMIT_EXCEEDED",
                context={"max_output_bytes": self._limit},
            )
        written = self._raw.write(value)
        self._count += written
        return written

    def flush(self) -> None:
        self._raw.flush()

    def __getattr__(self, name: str) -> object:
        return getattr(self._raw, name)


class _BoundedTextWriter:
    """Count encoded output bytes before forwarding writes to a text stream."""

    def __init__(self, raw: TextIO, *, encoding: str, limit: int) -> None:
        self._raw = raw
        self._encoding = encoding
        self._limit = limit
        self._count = 0

    def write(self, value: str) -> int:
        size = len(value.encode(self._encoding))
        if self._count + size > self._limit:
            raise InputFormatError(
                "Serialized output exceeds max_output_bytes.",
                code="OUTPUT_BYTE_LIMIT_EXCEEDED",
                context={"max_output_bytes": self._limit},
            )
        written = self._raw.write(value)
        self._count += size
        return written

    def flush(self) -> None:
        self._raw.flush()

    def __getattr__(self, name: str) -> object:
        return getattr(self._raw, name)


def _path(value: object) -> Path | None:
    if isinstance(value, (str, os.PathLike)):
        return Path(value)
    return None


def _stream_name(stream: object) -> str | None:
    value = getattr(stream, "name", None)
    if not isinstance(value, (str, os.PathLike)):
        return None
    name = os.fspath(value)
    return name if isinstance(name, str) and name else None


def _normalize_format(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("format must be a string or None.")
    normalized = value.strip().lower().lstrip(".")
    if normalized.endswith(".gz"):
        normalized = normalized[:-3]
    if normalized in _FORMAT_ALIASES:
        return _FORMAT_ALIASES[normalized]
    if normalized in _UNSUPPORTED_FORMATS:
        raise InputFormatError(
            f"Format {value!r} is not supported by the MVP.",
            code="UNSUPPORTED_FORMAT",
            context={"format": normalized},
            hint=_UNSUPPORTED_FORMATS[normalized],
        )
    raise InputFormatError(
        f"Unknown DNAKit I/O format {value!r}.",
        code="UNSUPPORTED_FORMAT",
        context={"format": normalized},
        hint="Choose FASTA, FASTQ, GenBank, CSV, TSV, JSON, or JSONL.",
    )


def _infer_format(name: str | None, explicit: str | None) -> str:
    if explicit is not None:
        return _normalize_format(explicit)
    if name is None:
        raise InputFormatError(
            "format is required for a stream without a usable filename.",
            code="FORMAT_REQUIRED",
            hint=("Pass format='fasta', 'fastq', 'genbank', 'csv', 'tsv', 'json', or 'jsonl'."),
        )
    filename = Path(name).name.lower()
    if filename.endswith(".gz"):
        filename = filename[:-3]
    suffix = Path(filename).suffix.lstrip(".")
    if not suffix:
        raise InputFormatError(
            "Could not infer a format from the filename.",
            code="FORMAT_REQUIRED",
            context={"source": name},
        )
    return _normalize_format(suffix)


def _use_gzip(name: str | None, compression: str) -> bool:
    if compression == "gzip":
        return True
    if compression == "none":
        return False
    return bool(name and name.lower().endswith(".gz"))


def _open_read_text(
    source: ReadableSource, config: ReadConfig
) -> tuple[TextIO, Callable[[], None], str | None]:
    path = _path(source)
    name = str(path) if path is not None else _stream_name(source)
    compressed = _use_gzip(name, config.compression)
    if path is not None:
        try:
            if compressed:
                handle: TextIO = gzip.open(  # noqa: SIM115 - lazy owner is RecordSource.
                    path, "rt", encoding=config.encoding, newline=""
                )
            else:
                handle = path.open("r", encoding=config.encoding, newline="")
        except OSError as exc:
            raise InputFormatError(
                "Could not open input path.",
                code="INPUT_OPEN_FAILED",
                context={"source": str(path), "reason": str(exc)},
            ) from exc
        bounded = _BoundedTextReader(
            handle, encoding=config.encoding, limit=config.max_input_bytes, source=name
        )
        return cast(TextIO, bounded), handle.close, name

    if not hasattr(source, "read"):
        raise TypeError("source must be a path or a readable text/binary stream.")
    close_source = config.close_source is True
    probe = source.read(0)
    if isinstance(probe, str):
        text_handle = cast(TextIO, source)
        if compressed:
            if close_source:
                text_handle.close()
            raise ConfigurationError(
                "A text stream cannot be decompressed by DNAKit.",
                code="TEXT_STREAM_COMPRESSION_CONFLICT",
                hint="Pass an already decompressed text stream with compression='none'.",
            )
        bounded = _BoundedTextReader(
            text_handle, encoding=config.encoding, limit=config.max_input_bytes, source=name
        )
        return (
            cast(TextIO, bounded),
            text_handle.close if close_source else (lambda: None),
            name,
        )
    if not isinstance(probe, bytes):
        raise TypeError("source.read(0) must return str or bytes.")

    binary_handle = cast(BinaryIO, source)
    if compressed:
        gzip_handle = gzip.GzipFile(fileobj=binary_handle, mode="rb")
        gzip_wrapper = io.TextIOWrapper(gzip_handle, encoding=config.encoding, newline="")

        def close_binary_gzip() -> None:
            try:
                gzip_wrapper.close()
            finally:
                if close_source:
                    binary_handle.close()

        bounded = _BoundedTextReader(
            gzip_wrapper, encoding=config.encoding, limit=config.max_input_bytes, source=name
        )
        return cast(TextIO, bounded), close_binary_gzip, name

    binary_wrapper = io.TextIOWrapper(binary_handle, encoding=config.encoding, newline="")

    def release_binary() -> None:
        if binary_wrapper.closed:
            return
        if close_source:
            binary_wrapper.close()
        else:
            binary_wrapper.detach()

    bounded = _BoundedTextReader(
        binary_wrapper, encoding=config.encoding, limit=config.max_input_bytes, source=name
    )
    return cast(TextIO, bounded), release_binary, name


def _reader(handle: TextIO, format: str, config: ReadConfig) -> Iterator[DNARecord]:
    if format == "fasta":
        return iter_fasta(handle, config)
    if format == "fastq":
        return iter_fastq(handle, config)
    if format in {"csv", "tsv"}:
        delimiter = config.delimiter or ("," if format == "csv" else "\t")
        return iter_delimited(handle, config, delimiter=delimiter)
    if format == "json":
        return iter_json_array(handle, config)
    if format == "jsonl":
        return iter_json_lines(handle, config)
    if format == "genbank":
        return iter_genbank(handle, config)
    raise AssertionError(f"Unreachable normalized format: {format}")


def _guard_reader(
    iterator: Iterator[DNARecord], *, source_name: str | None, max_records: int
) -> Iterator[DNARecord]:
    try:
        for index, record in enumerate(iterator, start=1):
            if index > max_records:
                raise InputFormatError(
                    "Input exceeds the configured record limit.",
                    code="INPUT_RECORD_LIMIT_EXCEEDED",
                    context={"source": source_name, "max_records": max_records},
                )
            yield record
    except (InputFormatError, ConfigurationError):
        raise
    except UnicodeError as exc:
        raise InputFormatError(
            "Input could not be decoded with the configured encoding.",
            code="INPUT_DECODE_FAILED",
            context={"source": source_name, "reason": str(exc)},
        ) from exc
    except OSError as exc:
        raise InputFormatError(
            "Input stream could not be read or decompressed.",
            code="INPUT_READ_FAILED",
            context={"source": source_name, "reason": str(exc)},
        ) from exc


@overload
def read(
    source: ReadableSource,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
    mode: Literal["stream"] = "stream",
) -> RecordSource: ...


@overload
def read(
    source: ReadableSource,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
    mode: Literal["dna"],
) -> DNA: ...


def read(
    source: ReadableSource,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
    mode: ReadMode = "stream",
) -> RecordSource | DNA:
    """Read DNA with one name, choosing a lazy stream or unified DNA object."""

    resolved = ReadConfig() if config is None else config
    if not isinstance(resolved, ReadConfig):
        raise TypeError("config must be ReadConfig or None.")
    if mode not in {"stream", "dna"}:
        raise ConfigurationError(
            "mode must be 'stream' or 'dna'.",
            code="INVALID_READ_MODE",
            context={"mode": mode},
        )
    path = _path(source)
    name = str(path) if path is not None else _stream_name(source)
    resolved_format = _infer_format(name, format)
    handle, close_callback, source_name = _open_read_text(source, resolved)
    try:
        iterator = _guard_reader(
            _reader(handle, resolved_format, resolved),
            source_name=source_name,
            max_records=resolved.max_records,
        )
    except BaseException:
        close_callback()
        raise
    records = RecordSource(
        iterator,
        close_callback=close_callback,
        source_name=source_name,
        format=resolved_format,
    )
    if mode == "dna":
        return DNA(records.collect())
    return records


def read_one(
    source: ReadableSource,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
) -> DNARecord:
    """Read exactly one record, rejecting empty and multi-record inputs."""

    with read(source, format=format, config=config, mode="stream") as records:
        try:
            first = next(records)
        except StopIteration as exc:
            raise InputFormatError(
                "Expected exactly one record, but the input was empty.",
                code="EXPECTED_ONE_RECORD",
                context={"record_count": 0},
            ) from exc
        try:
            next(records)
        except StopIteration:
            return first
        raise InputFormatError(
            "Expected exactly one record, but the input contained multiple records.",
            code="EXPECTED_ONE_RECORD",
            context={"record_count": "at_least_2"},
        )


def read_set(
    source: ReadableSource,
    *,
    format: str | None = None,
    config: ReadConfig | None = None,
) -> DNASet:
    """Explicitly materialize every input record into a repeatable DNASet."""

    return read(source, format=format, config=config, mode="stream").collect()


def _records_with_ids(
    records: WritableRecords,
    config: WriteConfig,
    generated_ids: list[GeneratedID],
) -> Iterator[DNARecord]:
    if isinstance(records, (DNASequence, DNARecord)):
        values: Iterable[DNASequence | DNARecord] = (records,)
    elif isinstance(records, Iterable):
        values = records
    else:
        raise TypeError("records must be DNASequence, DNARecord, or an iterable of them.")
    for index, value in enumerate(values):
        if isinstance(value, DNARecord):
            yield value
        elif isinstance(value, DNASequence):
            if config.anonymous_id_policy == "error":
                raise InputFormatError(
                    "Anonymous DNASequence output is disabled by configuration.",
                    code="ANONYMOUS_SEQUENCE_ID_REQUIRED",
                    context={"input_index": index},
                )
            generated_id = f"{config.anonymous_id_prefix}{index + 1}"
            generated_ids.append(GeneratedID(index, generated_id))
            yield DNARecord(value, generated_id)
        else:
            raise TypeError(
                f"records[{index}] must be DNASequence or DNARecord, not {type(value).__name__}."
            )


def _writer(
    handle: TextIO,
    records: Iterable[DNARecord],
    format: str,
    config: WriteConfig,
) -> int:
    if format == "fasta":
        return write_fasta(handle, records, config)
    if format == "fastq":
        return write_fastq(handle, records, config)
    if format in {"csv", "tsv"}:
        delimiter = config.delimiter or ("," if format == "csv" else "\t")
        return write_delimited(handle, records, config, delimiter=delimiter)
    if format == "json":
        return write_json_array(handle, records, config)
    if format == "jsonl":
        return write_json_lines(handle, records, config)
    if format == "genbank":
        return write_genbank(handle, records, config)
    raise AssertionError(f"Unreachable normalized format: {format}")


def _write_to_stream(
    target: TextIO | BinaryIO,
    records: Iterable[DNARecord],
    format: str,
    config: WriteConfig,
    *,
    compressed: bool,
) -> int:
    close_target = config.close_target is True
    try:
        cast(Any, target).write("")
        is_text_target = True
    except TypeError:
        is_text_target = False
    if is_text_target:
        text_target = cast(TextIO, target)
        if compressed:
            if close_target:
                text_target.close()
            raise ConfigurationError(
                "A text target cannot be gzip-compressed by DNAKit.",
                code="TEXT_STREAM_COMPRESSION_CONFLICT",
                hint="Pass a binary stream or an output path ending in .gz.",
            )
        try:
            destination = (
                cast(
                    TextIO,
                    _BoundedTextWriter(
                        text_target,
                        encoding=config.encoding,
                        limit=config.max_output_bytes,
                    ),
                )
                if config.max_output_bytes is not None
                else text_target
            )
            count = _writer(destination, records, format, config)
            text_target.flush()
            return count
        finally:
            if close_target:
                text_target.close()

    binary_target = cast(BinaryIO, target)
    bounded_target = (
        cast(
            BinaryIO,
            _BoundedBinaryWriter(binary_target, limit=config.max_output_bytes),
        )
        if config.max_output_bytes is not None
        else binary_target
    )
    raw: BinaryIO | None = None
    try:
        raw = (
            cast(
                BinaryIO,
                gzip.GzipFile(
                    filename="",
                    fileobj=bounded_target,
                    mode="wb",
                    compresslevel=config.compression_level,
                    mtime=0,
                ),
            )
            if compressed
            else bounded_target
        )
        wrapper = io.TextIOWrapper(raw, encoding=config.encoding, newline="")
    except BaseException:
        if compressed and raw is not None:
            with suppress(BaseException):
                raw.close()
        if close_target:
            binary_target.close()
        raise
    try:
        count = _writer(wrapper, records, format, config)
        wrapper.flush()
        return count
    finally:
        if compressed:
            try:
                wrapper.close()
            finally:
                if close_target:
                    binary_target.close()
        elif close_target:
            wrapper.close()
        elif not wrapper.closed:
            wrapper.detach()


def _artifact(path: Path, format: str) -> ArtifactRef:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = path.stat()
    return ArtifactRef(
        relative_path=os.path.relpath(path.resolve(), Path.cwd()),
        media_type=_MEDIA_TYPES[format],
        schema_version="dnakit-io-v1",
        sha256=digest.hexdigest(),
        byte_size=stat.st_size,
        created_at=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
    )


def _write_path(
    path: Path,
    records: Iterable[DNARecord],
    format: str,
    config: WriteConfig,
) -> tuple[int, ArtifactRef]:
    if path.exists() and not config.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    parent = path.parent
    if not parent.exists():
        if config.create_parents:
            parent.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Output parent directory does not exist: {parent}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=parent
    )
    temporary_path = Path(temporary_name)
    os.close(descriptor)
    compressed = _use_gzip(str(path), config.compression)
    try:
        with temporary_path.open("wb") as binary_target:
            stream_config = WriteConfig(
                encoding=config.encoding,
                overwrite=config.overwrite,
                create_parents=config.create_parents,
                compression=config.compression,
                compression_level=config.compression_level,
                close_target=False,
                line_width=config.line_width,
                phred_offset=config.phred_offset,
                delimiter=config.delimiter,
                anonymous_id_policy=config.anonymous_id_policy,
                anonymous_id_prefix=config.anonymous_id_prefix,
                feature_policy=config.feature_policy,
                json_indent=config.json_indent,
                max_output_bytes=config.max_output_bytes,
            )
            count = _write_to_stream(
                binary_target, records, format, stream_config, compressed=compressed
            )
            binary_target.flush()
            os.fsync(binary_target.fileno())
        if config.overwrite:
            os.replace(temporary_path, path)
        else:
            os.link(temporary_path, path)
            temporary_path.unlink()
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return count, _artifact(path, format)


def write(
    records: WritableRecords,
    target: WritableTarget,
    *,
    format: str | None = None,
    config: WriteConfig | None = None,
) -> WriteResult:
    """Serialize records with stable anonymous IDs and safe path replacement."""

    resolved = WriteConfig() if config is None else config
    if not isinstance(resolved, WriteConfig):
        raise TypeError("config must be WriteConfig or None.")
    path = _path(target)
    name = str(path) if path is not None else _stream_name(target)
    resolved_format = _infer_format(name, format)
    generated_ids: list[GeneratedID] = []
    prepared = _records_with_ids(records, resolved, generated_ids)
    if path is not None:
        count, artifact = _write_path(path, prepared, resolved_format, resolved)
        byte_count: int | None = artifact.byte_size
    else:
        if not hasattr(target, "write"):
            raise TypeError("target must be a path or writable text/binary stream.")
        count = _write_to_stream(
            cast(TextIO | BinaryIO, target),
            prepared,
            resolved_format,
            resolved,
            compressed=_use_gzip(name, resolved.compression),
        )
        artifact = None
        byte_count = None
    return WriteResult(
        resolved_format,
        count,
        byte_count=byte_count,
        generated_ids=generated_ids,
        target_artifact=artifact,
        parameters=cast(dict[str, object], to_json_compatible(resolved)),
        provenance=Provenance(
            implementation=ImplementationInfo(
                label=ImplementationLabel.REIMPLEMENTATION,
                execution_mode=ExecutionMode.INTERNAL,
                origin_class=OriginClass.STANDARD,
            )
        ),
    )


__all__ = [
    "PathSource",
    "ReadMode",
    "ReadableSource",
    "WritableRecords",
    "WritableTarget",
    "read",
    "read_one",
    "read_set",
    "write",
]
