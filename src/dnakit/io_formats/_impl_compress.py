"""Deterministic gzip codec behind ``compress_<format>`` / ``decompress_<format>``.

The sixteen public functions share one implementation: the bytes of the input (text, bytes, a file, or - for
``zarr`` - a directory archived as tar) are gzip-compressed with a fixed header (mtime 0) so that the same input
always gives the same output, and decompression returns the exact bytes again.  Nothing is written unless
``options['output']`` is given.
"""

# ruff: noqa: E501

from __future__ import annotations

import base64
import binascii
import gzip
import io
import os
import stat
import tarfile
import zlib
from collections.abc import Iterator, Mapping
from pathlib import Path, PurePosixPath

from dnakit._candidate_runtime import candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io_formats._impl_common import (
    DEFAULT_MAX_BYTES,
    KNOWN_SUFFIXES,
    Opts,
    _write_file,
    bad_input,
    looks_like_path,
    sha256_hex,
)

MODULE = "io_formats"
FORMATS = ("map", "newick", "nexus", "paf", "ped", "phylip", "twobit", "zarr")
_GZIP_MAGIC = b"\x1f\x8b"
_ALLOWED = {
    "compress": ("max_bytes", "level", "output", "overwrite"),
    "decompress": ("max_bytes", "output", "overwrite"),
}


def _too_large(what: str, size: int, limit: int) -> ConfigurationError:
    return ConfigurationError(
        f"{what} has {size} bytes, above max_bytes={limit}.",
        code="IO_FORMATS_TOO_LARGE",
        hint="Raise the max_bytes option.",
    )


def _not_found(source: object) -> InputFormatError:
    return bad_input(
        f"Input file does not exist: {source}",
        code="IO_FORMATS_FILE_NOT_FOUND",
        hint="Pass an existing path or the content itself.",
    )


def _is_existing_path(source: object) -> bool:
    if isinstance(source, os.PathLike):
        return True
    return isinstance(source, str) and looks_like_path(source) and os.path.exists(source)


def _looks_like_missing_file(source: str) -> bool:
    """A one-line string that names a file of a known suffix (a typo'd path, not content)."""

    return (
        looks_like_path(source)
        and Path(source).suffix.lower() in KNOWN_SUFFIXES | {".tar", ".fa", ".fasta", ".fna"}
        and not source.lstrip().startswith(("(", "#", ">"))
    )


def _read_file(path: Path, limit: int) -> bytes:
    if not path.is_file():
        raise _not_found(path)
    size = path.stat().st_size
    if size > limit:
        raise _too_large(path.name, size, limit)
    return path.read_bytes()


# ---------------------------------------------------------------------------------------------
# tar archives of directories (zarr stores)
# ---------------------------------------------------------------------------------------------


def _walk(directory: str, prefix: str) -> Iterator[tuple[str, str, bool]]:
    for name in sorted(os.listdir(directory)):
        full = os.path.join(directory, name)
        arcname = prefix + name
        mode = os.lstat(full).st_mode
        if stat.S_ISLNK(mode):
            raise bad_input(
                f"Symbolic links cannot be archived: {full}",
                code="IO_FORMATS_UNSAFE_ARCHIVE",
                hint="Replace the link by the file or directory it points to.",
            )
        if stat.S_ISDIR(mode):
            yield arcname, full, True
            yield from _walk(full, arcname + "/")
        elif stat.S_ISREG(mode):
            yield arcname, full, False
        else:
            raise bad_input(
                f"Only regular files and directories can be archived: {full}",
                code="IO_FORMATS_UNSAFE_ARCHIVE",
            )


def _tar_directory(root: Path, limit: int) -> bytes:
    """A reproducible tar archive of ``root`` (sorted entries, no times/owners, paths relative to ``root``)."""

    buffer = io.BytesIO()
    total = 0
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.GNU_FORMAT, encoding="utf-8") as archive:
        for arcname, full, is_dir in _walk(str(root), ""):
            info = tarfile.TarInfo(arcname)
            info.mtime, info.uid, info.gid, info.uname, info.gname = 0, 0, 0, "", ""
            if is_dir:
                info.type, info.mode = tarfile.DIRTYPE, 0o755
                archive.addfile(info)
                continue
            info.size, info.mode = os.path.getsize(full), 0o644
            total += info.size
            if total > limit:
                raise _too_large(f"The directory {root}", total, limit)
            with open(full, "rb") as handle:
                archive.addfile(info, handle)
    return buffer.getvalue()


def _safe_parts(name: str) -> tuple[str, ...]:
    """The path components of a member name, rejecting absolute paths and ``..``."""

    posix = name.replace("\\", "/")
    pure = PurePosixPath(posix)
    if posix.startswith("/") or pure.is_absolute() or (len(posix) > 1 and posix[1] == ":"):
        raise bad_input(
            f"Unsafe archive member (absolute path): {name!r}", code="IO_FORMATS_UNSAFE_ARCHIVE"
        )
    parts = tuple(part for part in pure.parts if part != ".")
    if ".." in parts:
        raise bad_input(
            f"Unsafe archive member (path escapes the target with '..'): {name!r}",
            code="IO_FORMATS_UNSAFE_ARCHIVE",
        )
    if "\x00" in name:
        raise bad_input("Unsafe archive member (NUL in name).", code="IO_FORMATS_UNSAFE_ARCHIVE")
    return parts


def _open_tar(data: bytes) -> tarfile.TarFile | None:
    """The tar archive held by ``data``, or ``None`` when the bytes are not a (readable) tar archive."""

    try:
        archive = tarfile.open(fileobj=io.BytesIO(data), mode="r:", encoding="utf-8")  # noqa: SIM115
        archive.getmembers()
    except (tarfile.TarError, EOFError, OSError, UnicodeError):
        return None
    return archive


def _extract_tar(archive: tarfile.TarFile, target: Path, overwrite: bool, limit: int) -> int:
    """Extract regular files and directories of a tar archive into ``target``; return the file count.

    Everything is validated before anything is written.  Links, devices and other special members are rejected
    (a link could point outside the target), as are absolute paths and ``..`` components.
    """

    members = archive.getmembers()
    kinds: dict[tuple[str, ...], str] = {}
    plan: list[tuple[tarfile.TarInfo, tuple[str, ...]]] = []
    total = 0
    for member in members:
        parts = _safe_parts(member.name)
        if member.isdir():
            kind = "dir"
        elif member.isreg():
            kind = "file"
            total += member.size
            if total > limit:
                raise _too_large("The archive content", total, limit)
        else:
            raise bad_input(
                f"Unsafe archive member {member.name!r}: links and special files are not extracted.",
                code="IO_FORMATS_UNSAFE_ARCHIVE",
            )
        if not parts:
            continue  # the "." entry of the root
        for depth in range(1, len(parts)):
            if kinds.get(parts[:depth]) == "file":
                raise bad_input(
                    f"Archive member {member.name!r} lies below a file.", code="IO_FORMATS_BAD_ARCHIVE"
                )
        previous = kinds.get(parts)
        if previous is not None and (previous != kind or kind == "file"):
            raise bad_input(
                f"Archive member {member.name!r} occurs twice.", code="IO_FORMATS_BAD_ARCHIVE"
            )
        kinds[parts] = kind
        plan.append((member, parts))

    if target.exists():
        if not target.is_dir():
            raise ConfigurationError(
                f"Output path is not a directory: {target}", code="IO_FORMATS_BAD_OUTPUT"
            )
        if any(target.iterdir()) and not overwrite:
            raise ConfigurationError(
                f"Refusing to extract into the non-empty directory: {target}",
                code="IO_FORMATS_EXISTS",
                hint="Set overwrite=True or choose another output directory.",
            )
    target.mkdir(parents=True, exist_ok=True)
    base = os.path.realpath(target)
    files = 0
    for member, parts in plan:
        destination = os.path.join(base, *parts)
        if os.path.commonpath([base, os.path.realpath(destination)]) != base:
            raise bad_input(
                f"Unsafe archive member {member.name!r}: it resolves outside the output directory.",
                code="IO_FORMATS_UNSAFE_ARCHIVE",
            )
        try:
            if member.isdir():
                os.makedirs(destination, exist_ok=True)
                continue
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            source = archive.extractfile(member)
            assert source is not None
            with open(destination, "wb") as handle:
                handle.write(source.read())
            files += 1
        except OSError as exc:
            raise ConfigurationError(
                f"Cannot extract {member.name!r} into {target}: {exc}", code="IO_FORMATS_BAD_OUTPUT"
            ) from exc
    return files


# ---------------------------------------------------------------------------------------------
# compress
# ---------------------------------------------------------------------------------------------


def _compress_source(fmt: str, source: object, limit: int) -> tuple[bytes, bool]:
    """The bytes to compress and whether they are a tar archive of a directory."""

    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
        if len(data) > limit:
            raise _too_large("The input", len(data), limit)
        return data, False
    if _is_existing_path(source):
        path = Path(source).expanduser()  # type: ignore[arg-type]
        if path.is_dir():
            if fmt != "zarr":
                raise bad_input(
                    f"Input path is a directory: {path}",
                    code="IO_FORMATS_BAD_INPUT",
                    hint="Only the zarr functions accept a directory (it is archived as tar before compression).",
                )
            return _tar_directory(path, limit), True
        return _read_file(path, limit), False
    if isinstance(source, str):
        if _looks_like_missing_file(source):
            raise _not_found(source)
        try:
            data = source.encode("utf-8")
        except UnicodeEncodeError as exc:
            raise bad_input(f"The text cannot be encoded as UTF-8: {exc}") from exc
        if len(data) > limit:
            raise _too_large("The input", len(data), limit)
        return data, False
    raise bad_input(
        f"Unsupported input of type {type(source).__name__}; expected text, bytes or a path.",
        hint="Pass the file content as str/bytes or the path of a file"
        + (" or directory (zarr)." if fmt == "zarr" else "."),
    )


def _compress(fmt: str, function: str, source: object, opts: Opts) -> dict[str, object]:
    limit = opts.integer("max_bytes", DEFAULT_MAX_BYTES, 1) or DEFAULT_MAX_BYTES
    level = opts.integer("level", 6, 1, 9) or 6
    output = opts.text("output")
    overwrite = opts.boolean("overwrite", False)
    data, archived = _compress_source(fmt, source, limit)
    packed = gzip.compress(data, compresslevel=level, mtime=0)
    result: dict[str, object] = {
        "format": fmt,
        "encoding": "gzip+base64",
        "data": base64.b64encode(packed).decode("ascii"),
        "original_bytes": len(data),
        "compressed_bytes": len(packed),
        "sha256_original": sha256_hex(data),
    }
    if archived:
        result["archive"] = "tar"
    if output:
        target = Path(output).expanduser()
        _write_file(target, packed, overwrite)
        result["output_path"] = str(target)
    return result


# ---------------------------------------------------------------------------------------------
# decompress
# ---------------------------------------------------------------------------------------------


def _gzip_bytes(data: bytes) -> bytes:
    """A gzip stream given either directly or as base64 text."""

    if data[:2] == _GZIP_MAGIC:
        return data
    try:
        raw = base64.b64decode(b"".join(data.split()), validate=True)
    except (binascii.Error, ValueError) as exc:
        raise bad_input(
            "The input is neither a gzip stream nor valid base64 text.",
            code="IO_FORMATS_BAD_BASE64",
            hint="Pass the 'data' of a compress_* result, the bytes of a .gz file, or the path of a .gz file.",
        ) from exc
    if raw[:2] != _GZIP_MAGIC:
        raise bad_input(
            "The base64 payload is not a gzip stream.",
            code="IO_FORMATS_BAD_GZIP",
            hint="Pass the 'data' of a compress_* result.",
        )
    return raw


def _gunzip(data: bytes, limit: int) -> bytes:
    pieces: list[bytes] = []
    total = 0
    rest = data
    while rest:
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            chunk = decoder.decompress(rest, limit + 1 - total)
        except zlib.error as exc:
            raise bad_input(
                f"Invalid gzip stream: {exc}",
                code="IO_FORMATS_BAD_GZIP",
                hint="The payload must come from gzip/compress_* and be complete.",
            ) from exc
        total += len(chunk)
        if total > limit or decoder.unconsumed_tail:
            raise _too_large("The decompressed data", total, limit)
        if not decoder.eof:
            raise bad_input(
                "Truncated gzip stream.",
                code="IO_FORMATS_BAD_GZIP",
                hint="The payload must be complete.",
            )
        pieces.append(chunk)
        rest = decoder.unused_data
        if rest and rest[:2] != _GZIP_MAGIC:
            raise bad_input(
                "Unexpected data after the end of the gzip stream.", code="IO_FORMATS_BAD_GZIP"
            )
    return b"".join(pieces)


def _decompress_source(source: object, limit: int) -> bytes:
    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
    elif _is_existing_path(source):
        path = Path(source).expanduser()  # type: ignore[arg-type]
        if path.is_dir():
            raise bad_input(f"Input path is a directory, not a .gz file: {path}")
        data = _read_file(path, limit)
    elif isinstance(source, str):
        if _looks_like_missing_file(source):
            raise _not_found(source)
        data = source.encode("utf-8")
    else:
        raise bad_input(
            f"Unsupported input of type {type(source).__name__}; expected base64 text, gzip bytes or a .gz path.",
            hint="Pass the 'data' of a compress_* result.",
        )
    if len(data) > limit:
        raise _too_large("The input", len(data), limit)
    return _gzip_bytes(data)


def _decompress(fmt: str, function: str, source: object, opts: Opts) -> dict[str, object]:
    limit = opts.integer("max_bytes", DEFAULT_MAX_BYTES, 1) or DEFAULT_MAX_BYTES
    output = opts.text("output")
    overwrite = opts.boolean("overwrite", False)
    packed = _decompress_source(source, limit)
    data = _gunzip(packed, limit)
    result: dict[str, object] = {
        "format": fmt,
        "original_bytes": len(data),
        "compressed_bytes": len(packed),
        "sha256": sha256_hex(data),
    }
    archive = _open_tar(data) if fmt == "zarr" else None
    if archive is not None:  # the archive of a directory made by compress_zarr
        result["archive"] = "tar"
        if output:
            target = Path(output).expanduser()
            result["file_count"] = _extract_tar(archive, target, overwrite, limit)
            result["output_path"] = str(target)
        else:
            result["data_base64"] = base64.b64encode(data).decode("ascii")
        return result
    if output:
        target = Path(output).expanduser()
        _write_file(target, data, overwrite)
        result["output_path"] = str(target)
        return result
    try:
        result["text"] = data.decode("utf-8")
    except UnicodeDecodeError:
        result["data_base64"] = base64.b64encode(data).decode("ascii")
    return result


# ---------------------------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------------------------


def run(function_name: str, input: object, options: Mapping[str, object] | None) -> dict[str, object]:
    """Execute ``compress_<format>`` or ``decompress_<format>`` and wrap the result in the envelope."""

    operation, _, fmt = function_name.partition("_")
    if operation not in _ALLOWED or fmt not in FORMATS:
        raise ConfigurationError(f"Unknown io_formats operation {function_name!r}.")
    opts = Opts(options, _ALLOWED[operation], f"io_formats.{function_name}")
    opts.finish()
    work = _compress if operation == "compress" else _decompress
    result = work(fmt, function_name, input, opts)
    return candidate_ok(MODULE, function_name, input, result, options=options)  # type: ignore[return-value]
