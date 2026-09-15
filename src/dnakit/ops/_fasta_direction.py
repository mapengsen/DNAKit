"""Bounded file-region processing for large, regularly wrapped genomic FASTA.

Workers read and write disjoint regions themselves. Only offsets cross process
boundaries; the parent never materializes or pickles chromosome-sized records.
Unusual layouts and explicit reader configurations retain the general codec.
"""

from __future__ import annotations

import mmap
import multiprocessing
import os
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import cast

from rich.console import Console

from dnakit._parallel import available_memory, resolve_workers
from dnakit._processing.inputs import input_config, output_path, protect_inputs
from dnakit._progress import Progress
from dnakit._verbosity import progress_enabled
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import DNAKitError, InputFormatError
from dnakit.io import ReadConfig, WriteConfig, WriteResult
from dnakit.io._formats import _header, _sequence
from dnakit.io.api import _artifact

_MIN_BYTES = 16 * 1024**2
_BLOCK_BYTES = 4 * 1024**2
_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn", ".faa"}
_COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")


class _Fallback(Exception):
    """Retry through the general reader, including its detailed diagnostics."""


@dataclass(frozen=True)
class _Record:
    start: int
    end: int
    bases: int
    line_bases: int
    newline_bytes: int
    header: bytes
    identifier: str
    output_start: int

    def byte_offset(self, base: int) -> int:
        lines, column = divmod(base, self.line_bases)
        return min(self.end, self.start + lines * (self.line_bases + self.newline_bytes) + column)


def _index(source: Path, line_width: int) -> tuple[list[_Record], int]:
    records: list[_Record] = []
    output_size = 0
    with (
        source.open("rb") as handle,
        mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data,
    ):
        if data[:1] != b">":
            raise _Fallback
        marks = [0]
        max_records = ReadConfig().max_records
        from dnakit._io_progress import read_observer

        observer = read_observer.get()
        cursor = 0
        while cursor < len(data):
            limit = min(len(data), cursor + 8 * 1024**2)
            position = data.find(b"\n>", cursor, min(len(data), limit + 1))
            while position != -1 and position < limit:
                marks.append(position + 1)
                if len(marks) > max_records:
                    raise _Fallback
                position = data.find(b"\n>", position + 2, min(len(data), limit + 1))
            cursor = limit
            if observer is not None:
                observer(str(source), cursor)
        marks.append(len(data))
        for mark, end in pairwise(marks):
            header_end = data.find(b"\n", mark, min(end, mark + 1024**2))
            if header_end == -1:
                if end - mark > 1024**2:
                    raise _Fallback
                header_end = end
            if data.find(b"\r", mark, max(mark, header_end - 1)) != -1:
                raise _Fallback
            identifier, description = _header(
                data[mark:header_end].decode("utf-8").rstrip("\r"), line_number=1, marker="FASTA"
            )
            header = (">" + identifier + (" " + description if description else "") + "\n").encode()
            start = min(end, header_end + 1)
            newline = data.find(b"\n", start, end)
            crlf = newline != -1 and data[newline - 1 : newline] == b"\r"
            newline_bytes = 2 if crlf else 1
            if crlf and data[end - 1 : end] == b"\n" and data[end - 2 : end] != b"\r\n":
                raise _Fallback
            width = end - start if newline == -1 else newline - start - int(crlf)
            if width == 0 and start != end:
                raise _Fallback
            if start == end:
                bases, width = 0, 1
            elif newline == -1:
                bases = width
            else:
                full, rest = divmod(end - start, width + newline_bytes)
                trailing = newline_bytes if rest and data[end - 1 : end] == b"\n" else 0
                if rest < trailing or rest - trailing > width:
                    raise _Fallback
                bases = full * width + rest - trailing
            output_size += len(header)
            records.append(
                _Record(start, end, bases, width, newline_bytes, header, identifier, output_size)
            )
            output_size += bases + (bases + line_width - 1) // line_width
    return records, output_size


@dataclass(frozen=True)
class _Job:
    source: Path
    target: Path
    record: _Record
    start: int
    end: int
    line_width: int
    reader: ReadConfig
    reverse: bool = True
    complement: bool = True


def _process(job: _Job) -> int:
    record = job.record
    first = record.byte_offset(record.bases - job.end if job.reverse else job.start)
    last = record.byte_offset(record.bases - job.start if job.reverse else job.end)
    with job.source.open("rb") as source:
        source.seek(first)
        raw = source.read(last - first).translate(None, b"\r\n")
    if len(raw) != job.end - job.start:
        return -1  # Irregular wrapping: retry using the authoritative codec.
    try:
        sequence = _sequence(
            raw.decode("utf-8"), job.reader, record_id=record.identifier, line_number=1
        )
    except (DNAKitError, UnicodeError):
        return -1
    # Input is validated once. Complementing valid IUPAC symbols preserves validity.
    bases = sequence.symbols.encode("ascii")
    if job.reverse:
        bases = bases[::-1]
    if job.complement:
        bases = bases.translate(_COMPLEMENT)
    output = (
        b"\n".join(bases[i : i + job.line_width] for i in range(0, len(bases), job.line_width))
        + b"\n"
    )
    with job.target.open("r+b") as target:
        target.seek(record.output_start + job.start + job.start // job.line_width)
        target.write(output)
    return len(bases)


def try_reverse_fasta(
    source: Path,
    target: Path | None,
    *,
    read_config: ReadConfig | None,
    write_config: WriteConfig | None,
    sequence_column: str | None,
    workers: int,
    show_progress: bool,
    reverse: bool = True,
    complement: bool = True,
) -> WriteResult | None:
    """Return None when the general file adapter should handle this request."""
    if (
        read_config is not None
        or sequence_column is not None
        or source.suffix.lower() not in _SUFFIXES
    ):
        return None
    config = cast(object, write_config or WriteConfig())
    if not isinstance(config, WriteConfig):
        return None
    source = source.expanduser().resolve()
    operation = (
        "reverse_complement" if reverse and complement else "reverse" if reverse else "complement"
    )
    target = (output_path(source, operation) if target is None else target).expanduser().resolve()
    if (
        target.suffix.lower() not in _SUFFIXES
        or config.compression == "gzip"
        or config.max_output_bytes is not None
        or config.encoding.lower().replace("_", "-") not in {"utf-8", "utf8"}
        or source.stat().st_size < _MIN_BYTES
    ):
        return None
    protect_inputs((source,), target)
    if target.exists() and not config.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {target}")
    original_stat = source.stat()
    try:
        records, total_bytes = _index(source, config.line_width)
    except (_Fallback, UnicodeError, DNAKitError):
        return None
    if config.create_parents:
        target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    temporary = Path(name)
    os.close(descriptor)
    block = max(config.line_width, _BLOCK_BYTES // config.line_width * config.line_width)
    reader = input_config(source, None)

    def jobs() -> Iterator[_Job]:
        for record in records:
            for start in range(0, record.bases, block):
                yield _Job(
                    source,
                    temporary,
                    record,
                    start,
                    min(record.bases, start + block),
                    config.line_width,
                    reader,
                    reverse,
                    complement,
                )

    task_count = sum((record.bases + block - 1) // block for record in records)
    count = min(
        resolve_workers(workers),
        8 if original_stat.st_size <= 512 * 1024**2 else resolve_workers(workers),
        max(1, task_count),
        max(1, available_memory() // max(128 * 1024**2, block * 12)),
    )
    try:
        with temporary.open("wb") as handle:
            handle.truncate(total_bytes)
            for record in records:
                handle.seek(record.output_start - len(record.header))
                handle.write(record.header)
        with Progress(
            console=Console(stderr=True), disable=not progress_enabled(show_progress)
        ) as display:
            task = display.add_task(f"{operation} (bases)", total=sum(r.bases for r in records))

            def collect(results: Iterator[int]) -> None:
                for completed in results:
                    if completed < 0:
                        raise _Fallback
                    display.advance(task, completed)

            if count == 1:
                collect(map(_process, jobs()))
            else:
                with multiprocessing.get_context("spawn").Pool(count) as pool:
                    collect(pool.imap(_process, jobs(), chunksize=1))
                    pool.close()
                    pool.join()
        current_stat = source.stat()
        if (original_stat.st_ino, original_stat.st_size, original_stat.st_mtime_ns) != (
            current_stat.st_ino,
            current_stat.st_size,
            current_stat.st_mtime_ns,
        ):
            raise InputFormatError("Input changed during processing.", code="INPUT_CHANGED")
        with temporary.open("r+b") as sync_handle:
            os.fsync(sync_handle.fileno())
        if config.overwrite:
            os.replace(temporary, target)
        else:
            os.link(temporary, target)
        artifact = _artifact(target, "fasta")
        return WriteResult(
            "fasta",
            len(records),
            byte_count=artifact.byte_size,
            target_artifact=artifact,
            parameters=cast(dict[str, object], to_json_compatible(config)),
        )
    except _Fallback:
        return None
    finally:
        temporary.unlink(missing_ok=True)
