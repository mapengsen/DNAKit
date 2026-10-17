"""Shared plumbing of the io_formats file-format engines.

Every format (``map``, ``ped``, ``paf``, ``newick``, ``nexus``, ``phylip``, ``twobit``, ``zarr``) is a
small :class:`Format` subclass that knows how to parse a source into a *document*, expose its
*records*, and serialise a document again.  The 96 public ``<operation>_<format>`` functions and the
generic ``convert``/``filter``/... functions are bound to these engines by :func:`run`
(``_impl_dispatch.py``), so each operation has one implementation per format, not one per function.

Input handling is bounded: files and gzip streams are read up to ``max_bytes``, nothing is written
unless an explicit ``output`` / ``output_dir`` option is given, and no network access happens.
"""

# ruff: noqa: E501, RUF012

from __future__ import annotations

import base64
import hashlib
import math
import os
import re
import zlib
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

DEFAULT_MAX_BYTES = 256 * 1024 * 1024
DEFAULT_MAX_RECORDS = 200_000
DEFAULT_PAGE = 1000

KNOWN_SUFFIXES = frozenset(
    {
        ".map",
        ".ped",
        ".paf",
        ".nwk",
        ".newick",
        ".tre",
        ".tree",
        ".nex",
        ".nexus",
        ".nxs",
        ".phy",
        ".phylip",
        ".2bit",
        ".zarr",
        ".gz",
    }
)


def bad_input(
    message: str, *, code: str = "IO_FORMATS_BAD_INPUT", hint: str | None = None
) -> InputFormatError:
    return InputFormatError(message, code=code, hint=hint)


def bad_option(message: str, *, hint: str | None = None) -> ConfigurationError:
    return ConfigurationError(message, code="IO_FORMATS_BAD_OPTION", hint=hint)


# ---------------------------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------------------------


class Opts:
    """Typed, checked access to the ``options`` mapping; unknown keys are rejected by :meth:`finish`."""

    def __init__(
        self, options: Mapping[str, object] | None, allowed: Sequence[str], where: str
    ) -> None:
        if options is not None and not isinstance(options, Mapping):
            raise bad_option(
                f"{where}: options must be a mapping.", hint="Pass options={'key': value}."
            )
        self.raw: dict[str, object] = dict(options or {})
        self.allowed = set(allowed) | {"backend"}
        self.where = where
        self.used: set[str] = set()

    def _take(self, key: str) -> object:
        self.used.add(key)
        return self.raw.get(key)

    def has(self, key: str) -> bool:
        return self.raw.get(key) is not None

    def get(self, key: str, default: object = None) -> object:
        value = self._take(key)
        return default if value is None else value

    def text(
        self, key: str, default: str | None = None, choices: Sequence[str] | None = None
    ) -> str | None:
        value = self._take(key)
        if value is None:
            return default
        if not isinstance(value, (str, os.PathLike)):
            raise bad_option(f"{self.where}: option '{key}' must be a string.")
        value = str(value)
        if choices is not None:
            lowered = value.lower()
            match = next((c for c in choices if c.lower() == lowered), None)
            if match is None:
                raise bad_option(
                    f"{self.where}: option '{key}' must be one of {list(choices)}, got {value!r}."
                )
            return match
        return value

    def integer(
        self, key: str, default: int | None, minimum: int | None = 0, maximum: int | None = None
    ) -> int | None:
        value = self._take(key)
        if value is None:
            return default
        if isinstance(value, bool) or not isinstance(value, int):
            if isinstance(value, float) and value.is_integer():
                value = int(value)
            else:
                raise bad_option(f"{self.where}: option '{key}' must be an integer.")
        if (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
            raise bad_option(
                f"{self.where}: option '{key}'={value} is outside [{minimum}, {maximum}]."
            )
        return value

    def number(
        self,
        key: str,
        default: float | None = None,
        minimum: float | None = None,
        maximum: float | None = None,
    ) -> float | None:
        value = self._take(key)
        if value is None:
            return default
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise bad_option(f"{self.where}: option '{key}' must be a finite number.")
        if (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
            raise bad_option(
                f"{self.where}: option '{key}'={value} is outside [{minimum}, {maximum}]."
            )
        return float(value)

    def boolean(self, key: str, default: bool = False) -> bool:
        value = self._take(key)
        if value is None:
            return default
        if not isinstance(value, bool):
            raise bad_option(f"{self.where}: option '{key}' must be true or false.")
        return value

    def names(self, key: str) -> list[str] | None:
        """A list of strings, or a comma-separated string; ``None`` when absent."""

        value = self._take(key)
        if value is None:
            return None
        if isinstance(value, str):
            return [part.strip() for part in value.split(",") if part.strip()]
        if isinstance(value, Sequence) and all(
            isinstance(item, (str, int)) and not isinstance(item, bool) for item in value
        ):
            return [str(item) for item in value]
        raise bad_option(
            f"{self.where}: option '{key}' must be a list of strings or a comma-separated string."
        )

    def name_set(self, key: str) -> frozenset[str] | None:
        """Like :meth:`names` but as a set (fast membership tests for large lists)."""

        values = self.names(key)
        return None if values is None else frozenset(values)

    def mapping(self, key: str) -> dict[str, object] | None:
        value = self._take(key)
        if value is None:
            return None
        if not isinstance(value, Mapping):
            raise bad_option(f"{self.where}: option '{key}' must be a mapping.")
        return {str(k): v for k, v in value.items()}

    def finish(self) -> None:
        unknown = sorted(set(self.raw) - self.allowed)
        if unknown:
            raise bad_option(
                f"{self.where}: unknown option(s) {unknown}.",
                hint=f"Accepted options: {sorted(self.allowed - {'backend'})}.",
            )


# ---------------------------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------------------------


def _bounded_gunzip(data: bytes, limit: int) -> bytes:
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        out = decoder.decompress(data, limit + 1)
        if decoder.unconsumed_tail or len(out) > limit:
            raise bad_input(
                f"Decompressed input is larger than max_bytes={limit}.", hint="Raise max_bytes."
            )
        if not decoder.eof:
            raise bad_input("Truncated gzip stream.")
    except zlib.error as exc:
        raise bad_input(f"Invalid gzip stream: {exc}") from exc
    return out


def looks_like_path(value: str) -> bool:
    return bool(value) and "\n" not in value and "\x00" not in value and len(value) <= 4096


def read_source(source: object, max_bytes: int = DEFAULT_MAX_BYTES) -> tuple[bytes, Path | None]:
    """Return the bytes of a path, text or bytes input (gzip is unpacked) and the path, if any."""

    path: Path | None = None
    if isinstance(source, bytes | bytearray | memoryview):
        data = bytes(source)
    elif isinstance(source, os.PathLike) or (
        isinstance(source, str) and looks_like_path(source) and os.path.exists(source)
    ):
        path = Path(source).expanduser()
        if not path.is_file():
            raise bad_input(
                f"Input path is not a readable file: {path}", code="IO_FORMATS_FILE_NOT_FOUND"
            )
        size = path.stat().st_size
        if size > max_bytes:
            raise ConfigurationError(
                f"{path.name} has {size} bytes, above max_bytes={max_bytes}.",
                code="IO_FORMATS_TOO_LARGE",
                hint="Raise the max_bytes option or use a streaming/indexing operation.",
            )
        data = path.read_bytes()
    elif isinstance(source, str):
        if (
            looks_like_path(source)
            and Path(source).suffix.lower() in KNOWN_SUFFIXES
            and not source.lstrip().startswith(("(", "#", ">"))
        ):
            raise bad_input(
                f"Input file does not exist: {source}",
                code="IO_FORMATS_FILE_NOT_FOUND",
                hint="Pass an existing path or the file content itself.",
            )
        data = source.encode("utf-8")
    else:
        raise bad_input(
            f"Unsupported input of type {type(source).__name__}; expected a path, text or bytes.",
            hint="Pass a file path, the file content as str/bytes, or (for write) a list of records.",
        )
    if len(data) > max_bytes:
        raise ConfigurationError(
            f"Input has {len(data)} bytes, above max_bytes={max_bytes}.",
            code="IO_FORMATS_TOO_LARGE",
            hint="Raise the max_bytes option.",
        )
    if data[:2] == b"\x1f\x8b":
        data = _bounded_gunzip(data, max_bytes)
    return data, path


def decode_text(data: bytes) -> str:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise bad_input(f"Input is not valid UTF-8 text: {exc}") from exc
    if "\x00" in text:
        raise bad_input("Input contains NUL bytes; it is not a text file of this format.")
    return text


def line_table(data: bytes) -> Iterator[tuple[int, int, str]]:
    """Yield ``(byte_offset, byte_length_incl_newline, text_without_newline)`` for each line."""

    offset = 0
    for raw in data.splitlines(keepends=True):
        try:
            text = raw.rstrip(b"\r\n").decode("utf-8-sig" if offset == 0 else "utf-8")
        except UnicodeDecodeError as exc:
            raise bad_input(f"Input is not valid UTF-8 text near byte {offset}: {exc}") from exc
        if "\x00" in text:
            raise bad_input("Input contains NUL bytes; it is not a text file of this format.")
        yield offset, len(raw), text
        offset += len(raw)


_INT = re.compile(r"[+-]?\d+\Z")


def to_int(token: str) -> int:
    if not _INT.match(token):
        raise ValueError(f"{token!r} is not an integer")
    return int(token)


def to_float(token: str) -> float:
    if "_" in token:
        raise ValueError(f"{token!r} is not a number")
    value = float(token)
    if not math.isfinite(value):
        raise ValueError(f"{token!r} is not finite")
    return value


def fmt_number(value: float) -> str:
    return format(value, ".15g")


def chrom_key(name: str) -> tuple[int, int, str]:
    """Natural chromosome order: 1..22, X, Y, XY, M/MT, then anything else alphabetically."""

    short = name[3:] if name.lower().startswith("chr") else name
    if short.isdigit():
        return (0, int(short), "")
    special = {"X": 23, "Y": 24, "XY": 25, "M": 26, "MT": 26}
    if short.upper() in special:
        return (1, special[short.upper()], "")
    return (2, 0, name)


def parse_region(text: str) -> tuple[str, int, int]:
    """``name:start-end`` -> ``(name, start, end)`` (the caller defines the coordinate convention)."""

    match = re.fullmatch(r"(.+):([\d,]+)-([\d,]+)", text.strip())
    if not match:
        raise bad_option(f"Region {text!r} must look like 'name:start-end'.")
    start, end = int(match.group(2).replace(",", "")), int(match.group(3).replace(",", ""))
    if end < start:
        raise bad_option(f"Region {text!r} has end < start.")
    return match.group(1), start, end


def safe_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._") or "part"


def sha256_hex(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


# ---------------------------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------------------------


@dataclass
class Payload:
    """A serialised result: text, bytes, or a directory store written on disk (zarr)."""

    fmt: str
    text: str | None = None
    data: bytes | None = None
    path: str | None = None
    count: int = 0
    info: dict[str, object] = field(default_factory=dict)

    def blob(self) -> bytes:
        return self.data if self.data is not None else (self.text or "").encode("utf-8")


def _write_file(path: Path, blob: bytes, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise ConfigurationError(
            f"Refusing to overwrite existing file: {path}",
            code="IO_FORMATS_EXISTS",
            hint="Set overwrite=True or choose another output path.",
        )
    if path.is_dir():
        raise ConfigurationError(
            f"Output path is a directory: {path}", code="IO_FORMATS_BAD_OUTPUT"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)


def emit(payload: Payload, opts: Opts) -> dict[str, object]:
    """Turn a payload into a JSON-friendly result, writing ``output`` when requested."""

    output = opts.text("output")
    overwrite = opts.boolean("overwrite", False)
    result: dict[str, object] = {"format": payload.fmt, "record_count": payload.count}
    result.update(payload.info)
    if payload.path is not None:
        result["output_path"] = payload.path
        return result
    blob = payload.blob()
    result["bytes"] = len(blob)
    result["sha256"] = sha256_hex(blob)
    if output:
        target = Path(output).expanduser()
        _write_file(target, blob, overwrite)
        result["output_path"] = str(target)
    elif payload.data is not None:
        result["data_base64"] = base64.b64encode(blob).decode("ascii")
    else:
        result["text"] = payload.text
    return result


def emit_parts(
    fmt: str, parts: list[tuple[str, Payload]], suffix: str, opts: Opts
) -> dict[str, object]:
    """Result of a split: one entry per part; parts are written under ``output_dir`` when given."""

    directory = opts.text("output_dir")
    overwrite = opts.boolean("overwrite", False)
    entries: list[dict[str, object]] = []
    for number, (label, part) in enumerate(parts, 1):
        entry: dict[str, object] = {"name": label, "record_count": part.count}
        entry.update(part.info)
        if part.path is not None:
            entry["output_path"] = part.path
        else:
            blob = part.blob()
            entry["bytes"] = len(blob)
            entry["sha256"] = sha256_hex(blob)
            if directory:
                target = (
                    Path(directory).expanduser() / f"part_{number:04d}_{safe_label(label)}{suffix}"
                )
                _write_file(target, blob, overwrite)
                entry["output_path"] = str(target)
            elif part.data is not None:
                entry["data_base64"] = base64.b64encode(blob).decode("ascii")
            else:
                entry["text"] = part.text
        entries.append(entry)
    return {"format": fmt, "part_count": len(entries), "parts": entries}


def page(records: Sequence[object], opts: Opts, total: int | None = None) -> dict[str, object]:
    """Cursor paging used by the stream operations."""

    start = opts.integer("start", 0, 0)
    limit = opts.integer("limit", DEFAULT_PAGE, 1, 100_000)
    assert start is not None and limit is not None
    total = len(records) if total is None else total
    chunk = list(records[start : start + limit])
    end = start + len(chunk)
    return {
        "start": start,
        "limit": limit,
        "returned": len(chunk),
        "total": total,
        "next_start": end if end < total else None,
        "done": end >= total,
        "records": chunk,
    }


# ---------------------------------------------------------------------------------------------
# engine base class
# ---------------------------------------------------------------------------------------------

COMMON_OPTIONS: dict[str, dict[str, str]] = {
    "read": {
        "max_bytes": "int, default 268435456: refuse inputs larger than this",
        "max_records": "int, default 200000: raise instead of returning more records (use stream_* to page)",
    },
    "write": {
        "max_bytes": "int, default 268435456: refuse text/file inputs larger than this",
        "output": "path to write (optional; the content is returned inline otherwise)",
        "overwrite": "bool, default False: allow replacing an existing output file",
    },
    "validate": {"max_bytes": "int, default 268435456: refuse inputs larger than this"},
    "filter": {
        "max_bytes": "int, default 268435456",
        "output": "path to write (optional)",
        "overwrite": "bool, default False",
    },
    "sort": {
        "max_bytes": "int, default 268435456",
        "output": "path to write (optional)",
        "overwrite": "bool, default False",
    },
    "merge": {
        "max_bytes": "int, default 268435456",
        "output": "path to write (optional)",
        "overwrite": "bool, default False",
    },
    "split": {
        "max_bytes": "int, default 268435456",
        "output_dir": "directory to write the parts to as part_0001_<label><ext> (optional)",
        "overwrite": "bool, default False",
    },
    "stream": {
        "max_bytes": "int, default 268435456",
        "start": "int, default 0: index of the first record of the page",
        "limit": "int, default 1000: records per page (max 100000)",
    },
    "summarize": {"max_bytes": "int, default 268435456"},
    "index": {"max_bytes": "int, default 268435456"},
    "convert": {
        "max_bytes": "int, default 268435456",
        "to": "target format name (required)",
        "output": "path to write (optional)",
        "overwrite": "bool, default False",
    },
}


class Format:
    """Base class: subclasses implement parse/dump/records and the operation hooks they support."""

    name = ""
    title = ""
    noun = "records"
    suffix = ""
    kind = "text"  # "text" | "binary" | "store"
    #: format-specific options, ``{operation: {key: description}}``
    options: dict[str, dict[str, str]] = {}
    #: one-paragraph semantics per operation (shown in the docstrings)
    semantics: dict[str, str] = {}
    #: option keys that are recognised only to be rejected with a specific message, ``{operation: (keys...)}``
    rejected_options: dict[str, tuple[str, ...]] = {}
    #: conversion targets besides the format itself
    targets: tuple[str, ...] = ("json",)
    #: operations implemented as ``<op>_op(source, opts)`` instead of through parse/dump
    native_ops: frozenset[str] = frozenset()

    # -- parsing / serialising ------------------------------------------------------------
    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> Any:
        raise NotImplementedError

    def from_records(self, records: object, opts: Opts) -> Any:
        raise bad_input(f"Cannot build {self.title} from {type(records).__name__}.")

    def dump(self, doc: Any, opts: Opts) -> Payload:
        raise NotImplementedError

    def records(self, doc: Any) -> list[dict[str, object]]:
        raise NotImplementedError

    def describe(self, doc: Any) -> dict[str, object]:
        return {}

    def count(self, doc: Any) -> int:
        return len(self.records(doc))

    # -- loading --------------------------------------------------------------------------
    def coerce(self, source: object, opts: Opts, errors: list[str] | None = None) -> Any:
        if isinstance(source, (Mapping, list, tuple)):
            return self.from_records(source, opts)
        data, _ = read_source(
            source, opts.integer("max_bytes", DEFAULT_MAX_BYTES, 1) or DEFAULT_MAX_BYTES
        )
        return self.parse(data, errors, opts)

    def coerce_for_write(self, source: object, opts: Opts) -> Any:
        return self.coerce(source, opts)

    def bad(self, errors: list[str] | None, where: str, message: str) -> None:
        text = f"{where}: {message}" if where else message
        if errors is None:
            raise bad_input(
                f"Invalid {self.title}: {text}", code=f"IO_FORMATS_BAD_{self.name.upper()}"
            )
        errors.append(text)

    # -- operation hooks (return doc / payload) -------------------------------------------
    def validate_doc(
        self, doc: Any, opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        return {}

    def filter(self, doc: Any, opts: Opts) -> Any:
        raise NotImplementedError

    def sort(self, doc: Any, opts: Opts) -> Any:
        raise NotImplementedError

    def merge(self, docs: list[Any], opts: Opts) -> Any:
        raise NotImplementedError

    def split(self, doc: Any, opts: Opts) -> list[tuple[str, Any]]:
        raise NotImplementedError

    def summarize(self, doc: Any, opts: Opts) -> dict[str, object]:
        raise NotImplementedError

    def index(self, doc: Any, opts: Opts) -> dict[str, object]:
        raise NotImplementedError

    def convert(self, doc: Any, target: str, opts: Opts) -> Payload:
        raise NotImplementedError


def json_payload(fmt: str, value: object, count: int) -> Payload:
    import json

    return Payload(
        fmt="json",
        text=json.dumps(value, ensure_ascii=False, indent=2),
        count=count,
        info={"source_format": fmt},
    )


def sources_list(input: object) -> list[object]:
    """The inputs of a merge: a sequence of sources, or ``{"inputs": [...]}``."""

    if isinstance(input, Mapping):
        input = input.get("inputs", input.get("items"))
    if isinstance(input, (str, bytes, os.PathLike)) or not isinstance(input, Sequence) or not input:
        raise bad_input(
            "merge needs a non-empty list of inputs (paths, texts or bytes).",
            hint="Pass input=[source1, source2, ...].",
        )
    return list(input)
