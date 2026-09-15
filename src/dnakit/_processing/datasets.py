"""Collection operations and transactional exports for inferred file inputs."""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from collections.abc import Callable, Mapping
from contextlib import ExitStack, closing
from inspect import BoundArguments
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import BarColumn, SpinnerColumn, TextColumn, TimeElapsedColumn

from dnakit._progress import Progress
from dnakit._verbosity import progress_enabled, status_enabled
from dnakit.core import DNA, DNASet
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError
from dnakit.io import WriteConfig, write

from .inputs import iter_inputs, output_path, paths_in, protect_inputs, table_config


def _write_subset(job: Any) -> None:
    records, path, config, *counters = job
    if counters:
        from .progress import tracked_records

        records = tracked_records(records, counters[0])
    write(records, path, config=config, verbose=False)
    if counters:
        completed = counters[0].read_text().split()[0]
        counters[0].write_text(f"{completed} done")


def run_collection(
    function: Callable[..., Any],
    bound: BoundArguments,
    input_name: str,
    mode: str,
    options: dict[str, Any],
) -> Any:
    value = bound.arguments[input_name]
    sources = paths_in(value)
    config = options.get("read_config")
    column = options.get("sequence_column")
    show = options.get("show_progress", True)
    if not isinstance(show, bool):
        raise ConfigurationError("show_progress must be boolean.")
    show = progress_enabled(show)
    destination = options.get("output_path")
    if sources and destination is None:
        destination = output_path(sources[0], function.__name__, "")
    target = Path(destination).expanduser().resolve() if destination is not None else None
    writer_config = options.get("write_config") or WriteConfig()
    if not isinstance(writer_config, WriteConfig):
        raise TypeError("write_config must be WriteConfig or None.")
    if target is not None:
        protect_inputs(sources, target)
        if any(target in source.parents for source in sources):
            raise ConfigurationError("An output directory cannot contain its input files.")
        if target.exists() and not writer_config.overwrite:
            raise FileExistsError(f"Refusing to overwrite existing output directory: {target}")
    console = Console(stderr=True)
    if status_enabled(show) and sources:
        console.print(f"Starting {function.__name__}: {sources[0]}", markup=False, soft_wrap=True)
    if (
        sources
        and target is not None
        and function.__name__ == "split"
        and getattr(bound.arguments.get("config"), "method", None) in {"random", "kfold"}
    ):
        from .file_split import random_file_split

        result = random_file_split(
            function, bound, input_name, target, sources, writer_config, config, column, show
        )
        if status_enabled(show):
            console.print(f"split complete: {target}", markup=False, soft_wrap=True)
        return result
    with ExitStack() as stack:
        display = stack.enter_context(
            Progress(
                SpinnerColumn(),
                TextColumn("{task.description}: {task.completed:.0f} records"),
                BarColumn(),
                TimeElapsedColumn(),
                console=console,
                disable=not (show and sources),
            )
        )
        task = display.add_task("Reading inputs", total=None)

        def tracked(stream: Any) -> Any:
            for record in stream:
                display.advance(task)
                yield record
            display.update(task, description=function.__name__)

        if mode == "mapping":
            if not isinstance(value, Mapping):
                raise TypeError("This operation requires a mapping of split names to inputs.")
            converted = {}
            for name, item in value.items():
                stream = stack.enter_context(
                    closing(iter_inputs(item, read_config=config, sequence_column=column))
                )
                converted[name] = DNASet(tracked(stream))
            bound.arguments[input_name] = converted
        else:
            stream = stack.enter_context(
                closing(iter_inputs(value, read_config=config, sequence_column=column))
            )
            bound.arguments[input_name] = (
                DNASet(tracked(stream)) if mode == "set" else tracked(stream)
            )
        result = function(*bound.args, **bound.kwargs)
    if target is not None:
        _export(result, target, function.__name__, sources, writer_config, config, column)
    if status_enabled(show) and sources:
        console.print(f"{function.__name__} complete: {target}", markup=False, soft_wrap=True)
    return result


def _export(
    result: Any,
    target: Path,
    operation: str,
    sources: tuple[Path, ...],
    config: WriteConfig,
    read_config: Any,
    sequence_column: str | None,
    *,
    progress: Any = None,
) -> None:
    """Install a complete output folder only after all subsets and reports succeed."""
    if config.create_parents:
        target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{target.name}.", dir=target.parent) as directory:
        stage = Path(directory) / "output"
        stage.mkdir()
        outputs: dict[str, Any] = {}
        if hasattr(result, "subsets"):
            for subset in result.subsets:
                if not re.fullmatch(r"[A-Za-z0-9_-]+", subset.name):
                    raise ConfigurationError(
                        "Split output names must use letters, digits, '_' or '-'."
                    )
                outputs[subset.name] = subset.records
        elif isinstance(result, (DNASet, DNA)):
            outputs[operation] = result.records if isinstance(result, DNA) else result
        elif hasattr(result, "records"):
            outputs[operation] = result.records
        if hasattr(result, "representatives"):
            outputs["representatives"] = result.representatives
        jobs: list[Any] = []
        for name, records in outputs.items():
            # JSONL retains all annotations; CSV inputs retain their original table columns.
            csv_input = any(
                suffix.lower() in {".csv", ".tsv"} for path in sources for suffix in path.suffixes
            )
            extension = (".csv" if csv_input else ".jsonl") + (
                ".gz" if config.compression == "gzip" else ""
            )
            path = stage / f"{name}{extension}"
            selected = table_config(sources, path, config, read_config, sequence_column)
            jobs.append((records, path, selected))
        from dnakit._parallel import ordered_map

        from .indexed_records import IndexedRecords
        from .split_exports import write_indexed_outputs

        with ExitStack() as stack:
            if progress is not None and progress.enabled:
                jobs = [(*job, Path(directory) / f"{job[1].name}.progress") for job in jobs]
                stack.enter_context(progress.writing(jobs))
            # These iterables contain offsets only; never pickle materialized datasets.
            if write_indexed_outputs(jobs):
                pass
            elif jobs and all(isinstance(job[0], IndexedRecords) for job in jobs):
                for _ in ordered_map(_write_subset, jobs, batch_size=1):
                    pass
            else:
                for job in jobs:
                    _write_subset(job)
        if hasattr(result, "validate_sources"):
            result.validate_sources()
        payload = result.to_dict() if hasattr(result, "to_dict") else to_json_compatible(result)
        (stage / "result.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        backup = Path(directory) / "previous"
        existed = target.exists()
        if existed:
            if not config.overwrite:
                raise FileExistsError(f"Output already exists: {target}")
            os.rename(target, backup)
        try:
            os.rename(stage, target)
        except BaseException:
            if existed:
                os.rename(backup, target)
            raise
        if backup.is_dir():
            shutil.rmtree(backup)
