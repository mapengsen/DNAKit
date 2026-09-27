"""Streaming row I/O for ordinary files without a DNAKit table envelope.

CSV/TSV values remain strings on read. JSON and Parquet preserve their own
types. Use ``read_table``/``export_table`` when an explicit typed schema and
resource limits are required instead.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import sys
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from types import TracebackType
from typing import Any, cast

from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.io._formats import write_fasta
from dnakit.io.artifacts import _staged_path
from dnakit.io.config import WriteConfig
from dnakit.io.parquet import require_parquet_backend

_ALIASES = {"fa": "fasta", "fna": "fasta", "ndjson": "jsonl", "pq": "parquet"}
_FORMATS = {"csv", "tsv", "json", "jsonl", "parquet", "fasta"}


def _format(path: Path) -> str:
    suffixes = path.suffixes
    if suffixes == [".gz"]:
        raise ValueError(f"Missing format before gzip suffix: {path}")
    suffix = suffixes[-2] if suffixes[-1:] == [".gz"] else path.suffix
    name = suffix.lstrip(".").lower()
    name = _ALIASES.get(name, name)
    if name not in _FORMATS or (name == "parquet" and path.suffix == ".gz"):
        raise ValueError(f"Unsupported row file format: {path}")
    return name


def iter_rows(
    source: str | os.PathLike[str],
    *,
    batch_size: int = 4096,
    verbose: bool = True,
) -> Iterator[dict[str, Any]]:
    """Iterate CSV, TSV, JSONL or Parquet rows; JSON arrays are loaded in memory.

    Text files may be gzip compressed. Parquet requires the optional pyarrow
    backend. A partially consumed iterator can be explicitly closed.
    """
    path = Path(source)
    format = _format(path)
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")
    if format == "parquet":
        require_parquet_backend()
        import pyarrow.parquet as pq

        with cast(Any, pq).ParquetFile(path) as parquet:
            for batch in parquet.iter_batches(batch_size=batch_size):
                yield from batch.to_pylist()
        return
    if format == "fasta":
        raise ValueError("Read FASTA records with dnakit.read().")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        if format in {"csv", "tsv"}:
            yield from csv.DictReader(stream, delimiter="," if format == "csv" else "\t")
        else:
            values = (
                json.load(stream)
                if format == "json"
                else (json.loads(line) for line in stream if line.strip())
            )
            if isinstance(values, dict):
                raise ValueError("Expected an array of rows; use load_json() for objects.")
            for row in values:
                if not isinstance(row, dict):
                    raise ValueError("Each row must be an object.")
                yield row


class RowWriter:
    """Write rows to one or several formats in a single pass.

    Targets may be CSV, TSV, JSON array, JSONL, Parquet or FASTA paths.
    Column order is inferred from the first row or supplied explicitly.
    Empty input requires explicit columns. Each row must contain exactly
    those columns. CSV/TSV serialize None as an empty cell. FASTA exports only
    ``id`` and ``sequence``; use a second tabular target to retain metadata.

    All outputs are staged until the context exits successfully. Each file is
    replaced atomically, but multiple final replacements are not a filesystem
    transaction. Parquet buffers at most ``batch_size`` rows and fixes its
    schema from the first batch; use export_table for an explicit typed schema.
    """

    def __init__(
        self,
        target: str | os.PathLike[str] | Sequence[str | os.PathLike[str]],
        *,
        columns: Sequence[str] | None = None,
        overwrite: bool = False,
        batch_size: int = 4096,
    ) -> None:
        self.paths = tuple(
            Path(p) for p in ([target] if isinstance(target, (str, os.PathLike)) else target)
        )
        if not self.paths or len({p.resolve() for p in self.paths}) != len(self.paths):
            raise ValueError("Output paths must be non-empty and unique.")
        self.formats = tuple(_format(p) for p in self.paths)
        self.columns = None if columns is None else tuple(columns)
        if self.columns is not None and (
            not self.columns
            or len(set(self.columns)) != len(self.columns)
            or any(not isinstance(c, str) or not c for c in self.columns)
        ):
            raise ValueError("columns must contain unique non-empty names.")
        if batch_size < 1:
            raise ValueError("batch_size must be positive.")
        self.overwrite = overwrite
        self.batch_size = batch_size
        self.count = 0
        self._stack: ExitStack | None = None
        self._sinks: list[dict[str, Any]] = []
        self._started = False

    def __enter__(self) -> RowWriter:
        if self._stack is not None or self._sinks:
            raise RuntimeError("RowWriter cannot be entered twice.")
        stack = ExitStack()
        try:
            for path, format in zip(self.paths, self.formats, strict=True):
                if format == "parquet":
                    require_parquet_backend()
                temporary = stack.enter_context(_staged_path(path, overwrite=self.overwrite))
                sink: dict[str, Any] = {"format": format, "path": temporary, "buffer": []}
                if format != "parquet":
                    opener = gzip.open if path.suffix == ".gz" else open
                    sink["stream"] = stack.enter_context(
                        opener(temporary, "wt", encoding="utf-8", newline="")
                    )
                self._sinks.append(sink)
        except BaseException:
            stack.__exit__(*sys.exc_info())
            raise
        self._stack = stack
        return self

    def _start(self) -> None:
        if self._started:
            return
        if not self.columns:
            raise ValueError("Supply columns when writing an empty collection.")
        if any(not isinstance(c, str) or not c for c in self.columns):
            raise ValueError("columns must contain non-empty string names.")
        for sink in self._sinks:
            format = sink["format"]
            if format in {"csv", "tsv"}:
                sink["writer"] = csv.DictWriter(
                    sink["stream"],
                    fieldnames=self.columns,
                    delimiter="," if format == "csv" else "\t",
                )
                sink["writer"].writeheader()
            elif format == "json":
                sink["stream"].write("[")
            elif format == "fasta" and not {"id", "sequence"}.issubset(self.columns):
                raise ValueError("FASTA export requires id and sequence columns.")
        self._started = True

    def _flush_parquet(self, sink: dict[str, Any]) -> None:
        import pyarrow as pa
        import pyarrow.parquet as pq

        writer = sink.get("writer")
        if not sink["buffer"] and writer is not None:
            return
        table = (
            pa.Table.from_pylist(sink["buffer"], schema=writer.schema if writer else None)
            if sink["buffer"]
            else pa.table({name: pa.array([], type=pa.string()) for name in self.columns or ()})
        )
        if writer is None:
            writer = cast(Any, pq).ParquetWriter(sink["path"], table.schema)
            sink["writer"] = writer
            assert self._stack is not None
            self._stack.callback(writer.close)
        writer.write_table(table)
        sink["buffer"].clear()

    def write(self, row: Mapping[str, Any]) -> None:
        """Append one row to every target without retaining preceding rows."""
        if self._stack is None:
            raise RuntimeError("Use RowWriter inside a with block.")
        if self.columns is None:
            self.columns = tuple(row)
        self._start()
        if set(row) != set(self.columns or ()):
            raise ValueError("Row keys do not match the output columns.")
        value = {key: row[key] for key in self.columns or ()}
        for sink in self._sinks:
            format = sink["format"]
            if format in {"csv", "tsv"}:
                sink["writer"].writerow(value)
            elif format in {"json", "jsonl"}:
                if format == "json" and self.count:
                    sink["stream"].write(",")
                json.dump(value, sink["stream"], ensure_ascii=False, allow_nan=False)
                sink["stream"].write("\n")
            elif format == "fasta":
                record = DNARecord(DNASequence(value["sequence"], alphabet="iupac"), value["id"])
                write_fasta(
                    sink["stream"], [record], WriteConfig(line_width=max(1, len(record.sequence)))
                )
            else:
                sink["buffer"].append(value)
                if len(sink["buffer"]) >= self.batch_size:
                    self._flush_parquet(sink)
        self.count += 1

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool | None:
        assert self._stack is not None
        stack = self._stack
        try:
            if exc_type is None:
                self._start()
                for sink in self._sinks:
                    if sink["format"] == "json":
                        sink["stream"].write("]\n")
                    elif sink["format"] == "parquet":
                        self._flush_parquet(sink)
        except BaseException:
            stack.__exit__(*sys.exc_info())
            raise
        else:
            return stack.__exit__(exc_type, exc, traceback)
        finally:
            self._stack = None


def write_rows(
    rows: Iterable[Mapping[str, Any]],
    target: str | os.PathLike[str] | Sequence[str | os.PathLike[str]],
    *,
    columns: Sequence[str] | None = None,
    overwrite: bool = False,
    batch_size: int = 4096,
    verbose: bool = True,
) -> int:
    """Write a row iterator once to one or several paths; return its row count."""
    with RowWriter(target, columns=columns, overwrite=overwrite, batch_size=batch_size) as writer:
        for row in rows:
            writer.write(row)
    return writer.count
