"""Shared bounded text helpers for advanced codecs."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

from dnakit.exceptions import InputFormatError


@contextmanager
def open_text_source(source: str | os.PathLike[str] | TextIO) -> Iterator[TextIO]:
    """Borrow streams and own path handles."""

    if isinstance(source, (str, os.PathLike)):
        path = Path(source)
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                yield handle
        except OSError as exc:
            raise InputFormatError(
                "Could not open the advanced text input.",
                code="INPUT_OPEN_FAILED",
                context={"source": str(path), "reason": str(exc)},
            ) from exc
        return
    if not hasattr(source, "read"):
        raise TypeError("source must be a path or readable text stream.")
    probe = source.read(0)
    if not isinstance(probe, str):
        raise TypeError("advanced annotation codecs require a decoded text stream.")
    yield source


def checked_lines(handle: TextIO, *, max_line_length: int) -> Iterator[tuple[int, str]]:
    """Yield newline-free lines while enforcing a hard character limit."""

    if isinstance(max_line_length, bool) or not isinstance(max_line_length, int):
        raise TypeError("max_line_length must be an integer.")
    if max_line_length < 1:
        raise ValueError("max_line_length must be positive.")
    for line_number, raw in enumerate(handle, start=1):
        line = raw.rstrip("\r\n")
        if len(line) > max_line_length:
            raise InputFormatError(
                "An input line exceeds the configured resource limit.",
                code="LINE_TOO_LONG",
                context={
                    "line_number": line_number,
                    "line_length": len(line),
                    "max_line_length": max_line_length,
                },
            )
        yield line_number, line


def write_text_path(
    target: str | os.PathLike[str] | TextIO,
    writer: Callable[[TextIO], int],
    *,
    overwrite: bool,
    create_parents: bool,
) -> int:
    """Write text atomically for paths and directly to borrowed streams."""

    if not isinstance(overwrite, bool) or not isinstance(create_parents, bool):
        raise TypeError("overwrite and create_parents must be booleans.")
    if not isinstance(target, (str, os.PathLike)):
        if not hasattr(target, "write"):
            raise TypeError("target must be a path or writable text stream.")
        count = writer(target)
        target.flush()
        return count

    path = Path(target)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    if not path.parent.exists():
        if not create_parents:
            raise FileNotFoundError(f"Output parent directory does not exist: {path.parent}")
        path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(name)
    os.close(descriptor)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            count = writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
            temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return count


__all__ = ["checked_lines", "open_text_source", "write_text_path"]
