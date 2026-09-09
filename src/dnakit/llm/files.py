"""Bounded, explicit local text-file input for reusable LLM processing steps."""

from __future__ import annotations

import codecs
import os
import stat
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

_BINARY_SIGNATURES = (
    b"%PDF-",
    b"PK\x03\x04",
    b"PK\x05\x06",
    b"PK\x07\x08",
    b"\x1f\x8b",
    b"\x89PNG\r\n\x1a\n",
    b"\xff\xd8\xff",
    b"GIF87a",
    b"GIF89a",
    b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1",
    b"\x7fELF",
)


def _file_path(file_path: str | os.PathLike[str]) -> Path:
    try:
        value = os.fspath(file_path)
        if not isinstance(value, str) or not value.strip() or "\x00" in value or "://" in value:
            raise ValueError
        return Path(value).expanduser()
    except (TypeError, ValueError, RuntimeError):
        raise InputFormatError(
            "file_path must be a non-empty local path, not a URL or file handle.",
            code="INVALID_LLM_FILE_PATH",
        ) from None


def _nonblocking_open(path: str, flags: int) -> int:
    # If a checked regular file is replaced by a FIFO before open(), do not wait
    # for a writer. The opened descriptor is checked again before any read.
    return os.open(path, flags | getattr(os, "O_NONBLOCK", 0))


def read_llm_file(
    file_path: str | os.PathLike[str],
    *,
    encoding: str = "utf-8-sig",
    max_bytes: int = 64_000,
) -> str:
    """Read one complete bounded local text file, without any API call.

    Supports text irrespective of extension, including LICENSE, code, CSV, JSON,
    Markdown, and plain-text sequence files. UTF-8 (with optional BOM) is the
    default; other encodings must be selected explicitly. No directory traversal,
    PDF/Office extraction, archive expansion, or binary attachments are performed.
    Oversized/empty/binary files fail rather than being truncated or split. Only
    regular files are accepted, including explicitly selected symlinks to them.
    """

    if (
        isinstance(max_bytes, bool)
        or not isinstance(max_bytes, int)
        or not 1 <= max_bytes <= 100_000_000
    ):
        raise ConfigurationError(
            "max_bytes must be an integer in [1, 100000000].", code="INVALID_LLM_FILE_CONFIG"
        )
    try:
        if not isinstance(encoding, str) or not encoding.strip():
            raise LookupError
        # Empty bytes.decode() can skip codec lookup entirely. Resolve the
        # decoder explicitly and reject byte-to-byte transformation codecs.
        decoder = codecs.getincrementaldecoder(encoding)()
        if not isinstance(decoder.decode(b"", final=False), str):
            raise LookupError
    except (LookupError, TypeError, UnicodeError):
        raise ConfigurationError(
            "encoding must name a supported text encoding.", code="INVALID_LLM_FILE_CONFIG"
        ) from None
    path = _file_path(file_path)

    def validate_file(info: os.stat_result) -> None:
        if not stat.S_ISREG(info.st_mode):
            raise InputFormatError(
                "file_path must refer to one regular file, not a directory or special device.",
                code="INVALID_LLM_FILE_TYPE",
            )
        if info.st_size > max_bytes:
            raise InputFormatError(
                "The file exceeds max_file_bytes; no content was sent.",
                code="LLM_FILE_LIMIT",
                context={"max_file_bytes": max_bytes},
                hint="Select a smaller file or explicitly increase the file and request limits.",
            )

    try:
        validate_file(path.stat())
        with open(path, "rb", opener=_nonblocking_open) as handle:
            validate_file(os.fstat(handle.fileno()))
            raw = handle.read(max_bytes + 1)
    except FileNotFoundError:
        raise InputFormatError(
            "The selected LLM input file does not exist.", code="LLM_FILE_NOT_FOUND"
        ) from None
    except OSError:
        raise InputFormatError(
            "The selected LLM input file could not be read.", code="LLM_FILE_READ_ERROR"
        ) from None
    if len(raw) > max_bytes:
        raise InputFormatError(
            "The file grew beyond max_file_bytes while reading; no content was sent.",
            code="LLM_FILE_LIMIT",
        )
    if raw.startswith(_BINARY_SIGNATURES):
        raise InputFormatError(
            "Binary documents are not supported by the text-file input reader.",
            code="UNSUPPORTED_LLM_FILE_FORMAT",
            hint="Extract text from PDF/Office files explicitly before using input_data.",
        )
    try:
        content = raw.decode(encoding)
    except UnicodeError:
        raise InputFormatError(
            "The file cannot be decoded using file_encoding.",
            code="LLM_FILE_ENCODING_ERROR",
            hint="Specify the correct text encoding; binary files need separate text extraction.",
        ) from None
    if any(
        (ord(character) < 32 and character not in "\t\n\r\f") or 0x7F <= ord(character) <= 0x9F
        for character in content
    ):
        raise InputFormatError(
            "The file contains binary control characters rather than ordinary text.",
            code="UNSUPPORTED_LLM_FILE_FORMAT",
        )
    if not content.strip():
        raise InputFormatError("The selected file contains no text.", code="EMPTY_LLM_FILE")
    return content
