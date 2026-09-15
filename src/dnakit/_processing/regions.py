"""Process complete FASTA record regions without sending sequences through IPC."""

from __future__ import annotations

import inspect
import multiprocessing
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, cast

from rich.console import Console

from dnakit._parallel import (
    _function,
    available_memory,
    current_workers,
    pack,
    unpack,
    worker_scope,
)
from dnakit._progress import Progress
from dnakit._verbosity import _enabled, progress_enabled
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import DNAKitError, InputFormatError
from dnakit.io import WriteConfig, WriteResult
from dnakit.io.api import _artifact

_MIN_BYTES = 16 * 1024**2
_SHARD_BYTES = 4 * 1024**2


def _region_job(payload: bytes) -> bytes:
    try:
        return pack((True, _run_region_job(unpack(payload))))
    except Exception as exc:
        return pack((False, exc))


def _run_region_job(job: Any) -> tuple[str, int]:
    module, name, args, kwargs, input_name, mode, options, source, first, last = job
    from .sequence import run_sequence

    with source.open("rb") as handle:
        handle.seek(first)
        text = handle.read(last - first).decode("utf-8")
    function = _function(module, name)
    bound = inspect.signature(function).bind(*args, **kwargs)
    token = _enabled.set(False)
    try:
        with worker_scope(1):
            result = run_sequence(function, bound, input_name, mode, options, _region_text=text)
    except InputFormatError as exc:
        from .indexed_records import absolute_error

        raise absolute_error(exc, source, first) from exc
    finally:
        _enabled.reset(token)
    return str(options["output_path"]), result.record_count


def try_regions(
    function: Any,
    bound: Any,
    input_name: str,
    mode: str,
    options: dict[str, Any],
    sources: tuple[Path, ...],
    target: Path,
) -> WriteResult | None:
    # Keep richer/custom codecs and resource limits in their authoritative paths.
    if (
        current_workers() == 1
        or bound.arguments.get("rng") is not None
        or bound.arguments.get("progress") is not None
        or len(sources) != 1
        or not isinstance(bound.arguments[input_name], (str, os.PathLike))
        or sources[0].suffix.lower() not in {".fa", ".fasta", ".fna", ".ffn", ".faa"}
        or target.suffix.lower() not in {".fa", ".fasta", ".fna", ".jsonl"}
        or options.get("read_config") is not None
        or options.get("sequence_column") is not None
        or options.get("topology") is not None
        or mode not in {"sequence", "record", "text", "normalize", "report", "chunks"}
        or sources[0].stat().st_size < _MIN_BYTES
    ):
        return None
    config = options.get("write_config") or WriteConfig()
    if (
        not isinstance(config, WriteConfig)
        or config.compression == "gzip"
        or config.max_output_bytes is not None
        or config.encoding.lower().replace("_", "-") not in {"utf-8", "utf8"}
    ):
        return None
    show = options.get("show_progress", True)
    if not isinstance(show, bool):
        return None  # Normal adapter raises the public configuration error.
    from dnakit.ops._fasta_direction import _Fallback, _index

    original = sources[0].stat()
    try:
        records, _ = _index(sources[0], config.line_width)
    except (_Fallback, DNAKitError, UnicodeError):
        return None
    regions = []
    first = 0
    for record in records:
        if record.end - first >= _SHARD_BYTES:
            regions.append((first, record.end))
            first = record.end
    if first < original.st_size:
        regions.append((first, original.st_size))
    if len(regions) < 2:
        return None
    if target.exists() and not config.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {target}")
    if config.create_parents:
        target.parent.mkdir(parents=True, exist_ok=True)
    largest = max(end - start for start, end in regions)
    # Reports can contain a per-position object for every ambiguous base.
    memory_per_worker = max(512 * 1024**2, largest * (256 if mode == "report" else 12))
    count = min(current_workers(), len(regions), 8, max(1, available_memory() // memory_per_worker))
    with tempfile.TemporaryDirectory(
        prefix=f".{target.name}.regions-", dir=target.parent
    ) as folder:
        stage = Path(folder)
        jobs = [
            (
                function.__module__,
                function.__name__,
                bound.args,
                bound.kwargs,
                input_name,
                mode,
                {
                    **options,
                    "workers": 1,
                    "show_progress": False,
                    "output_path": stage / f"{i}{target.suffix}",
                },
                sources[0],
                first,
                last,
            )
            for i, (first, last) in enumerate(regions)
        ]
        with Progress(console=Console(stderr=True), disable=not progress_enabled(show)) as display:
            task = display.add_task(function.__name__, total=len(regions))

            def collect(results: Any) -> list[tuple[str, int]]:
                completed = []
                for payload in results:
                    success, result = unpack(payload)
                    if not success:
                        raise result
                    completed.append(result)
                    display.advance(task)
                return completed

            if count == 1:
                completed = collect(map(_region_job, map(pack, jobs)))
            else:
                with multiprocessing.get_context("spawn").Pool(count) as pool:
                    completed = collect(pool.imap(_region_job, map(pack, jobs), chunksize=1))
                    pool.close()
                    pool.join()
        current = sources[0].stat()
        if (original.st_ino, original.st_size, original.st_mtime_ns) != (
            current.st_ino,
            current.st_size,
            current.st_mtime_ns,
        ):
            raise InputFormatError("Input changed during processing.", code="INPUT_CHANGED")
        combined = stage / "combined"
        with combined.open("wb") as output:
            for path, _ in completed:
                with open(path, "rb") as fragment:
                    shutil.copyfileobj(fragment, output, length=1024**2)
                Path(path).unlink()
            output.flush()
            os.fsync(output.fileno())
        if config.overwrite:
            os.replace(combined, target)
        else:
            os.link(combined, target)
        fmt = "jsonl" if target.suffix.lower() == ".jsonl" else "fasta"
        artifact = _artifact(target, fmt)
        return WriteResult(
            fmt,
            sum(count for _, count in completed),
            byte_count=artifact.byte_size,
            target_artifact=artifact,
            parameters=(
                None
                if mode in {"text", "report"}
                else cast(dict[str, object], to_json_compatible(config))
            ),
        )
