"""Reopen FASTA records by offsets for random splitting without a sequence spool."""

from __future__ import annotations

import io
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import DNAKitError, InputFormatError
from dnakit.io._formats import iter_fasta

from .inputs import input_config, is_path


def absolute_error(error: InputFormatError, path: Path, offset: int) -> InputFormatError:
    """Keep file-relative diagnostics when a parser consumes only one region."""
    line = error.context.get("line_number")
    if not offset or not isinstance(line, int):
        return InputFormatError(
            error.message, code=error.code, context=error.context, hint=error.hint
        )
    lines, remaining, previous_cr = 0, offset, False
    with path.open("rb") as handle:
        while remaining:
            block = handle.read(min(1024**2, remaining))
            if not block:
                break
            lines += block.count(b"\n") + block.count(b"\r") - block.count(b"\r\n")
            if previous_cr and block.startswith(b"\n"):
                lines -= 1
            previous_cr = block.endswith(b"\r")
            remaining -= len(block)
    return InputFormatError(
        error.message,
        code=error.code,
        context={**error.context, "line_number": line + lines},
        hint=error.hint,
    )


@dataclass(frozen=True)
class IndexedRecords:
    locations: tuple[tuple[Path, int, int], ...]
    stamps: tuple[tuple[Path, int, int, int], ...]
    layouts: tuple[tuple[int, bytes], ...] = ()

    def validate_sources(self) -> None:
        for path, inode, size, modified in self.stamps:
            stat = path.stat()
            if (stat.st_ino, stat.st_size, stat.st_mtime_ns) != (inode, size, modified):
                raise InputFormatError("Input changed during splitting.", code="INPUT_CHANGED")

    def __iter__(self) -> Iterator[DNARecord]:
        self.validate_sources()
        for path, start, end in self.locations:
            with path.open("rb") as source:
                source.seek(start)
                text = source.read(end - start).decode("utf-8")
            try:
                yield from iter_fasta(io.StringIO(text, newline=""), input_config(path, None))
            except InputFormatError as exc:
                raise absolute_error(exc, path, start) from exc
        self.validate_sources()

    def select(self, indices: list[int]) -> IndexedRecords:
        return IndexedRecords(
            tuple(self.locations[index] for index in indices),
            self.stamps,
            tuple(self.layouts[index] for index in indices) if self.layouts else (),
        )


def try_index(
    value: Any, sources: tuple[Path, ...]
) -> tuple[list[DNARecord], IndexedRecords] | None:
    def only_paths(item: Any) -> bool:
        return is_path(item) or (
            isinstance(item, (list, tuple)) and all(only_paths(part) for part in item)
        )

    if not only_paths(value) or any(
        path.suffix.lower() not in {".fa", ".fasta", ".fna", ".ffn", ".faa"}
        or path.stat().st_size == 0
        for path in sources
    ):
        return None
    from dnakit.ops._fasta_direction import _Fallback, _index

    headers, locations, stamps, layouts = [], [], [], []
    try:
        for source in sources:
            stat = source.stat()
            stamps.append((source, stat.st_ino, stat.st_size, stat.st_mtime_ns))
            records, _ = _index(source, 80)
            start = 0
            for record in records:
                locations.append((source, start, record.end))
                headers.append(DNARecord(DNASequence(""), record.identifier))
                layouts.append((record.bases, record.header))
                start = record.end
    except (_Fallback, DNAKitError, UnicodeError):
        return None
    return headers, IndexedRecords(tuple(locations), tuple(stamps), tuple(layouts))
