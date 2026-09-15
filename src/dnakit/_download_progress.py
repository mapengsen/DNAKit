"""Shared console feedback for public download APIs."""

from __future__ import annotations

from collections.abc import Callable
from functools import wraps
from inspect import signature
from pathlib import Path
from typing import Any, ParamSpec, TypeVar

from rich.console import Console
from rich.progress import (
    BarColumn,
    DownloadColumn,
    MofNCompleteColumn,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

from dnakit._progress import Progress
from dnakit._verbosity import progress_enabled, status_enabled
from dnakit.exceptions import ConfigurationError

P = ParamSpec("P")
R = TypeVar("R")


def with_download_progress(
    target_parameter: str, *, callback_parameter: str = "progress"
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Show default feedback while leaving explicit callbacks in control.

    Output goes to stderr so stdout remains usable for JSON and pipelines.
    Dataset callbacks also silence nested file displays, producing one session.
    """

    def decorate(function: Callable[P, R]) -> Callable[P, R]:
        parameters = signature(function)

        @wraps(function)
        def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
            bound = parameters.bind(*args, **kwargs)
            bound.apply_defaults()
            enabled = bound.arguments.get(
                "show_progress", getattr(bound.arguments.get("config"), "show_progress", True)
            )
            if not isinstance(enabled, bool):
                raise ConfigurationError("show_progress must be boolean.")
            enabled = progress_enabled(enabled)
            if not enabled or bound.arguments[callback_parameter] is not None:
                return function(*args, **kwargs)

            target = str(Path(bound.arguments[target_parameter]).expanduser().resolve())
            console = Console(stderr=True)
            if status_enabled(enabled):
                console.print(
                    f"[DNAKit] Starting download: {target}",
                    markup=False,
                    highlight=False,
                    soft_wrap=True,
                )
            is_query = callback_parameter == "query_progress"
            with Progress(
                SpinnerColumn(),
                TextColumn("{task.description}", markup=False),
                BarColumn(),
                TaskProgressColumn(),
                MofNCompleteColumn() if is_query else DownloadColumn(),
                TimeElapsedColumn(),
                console=console,
            ) as display:
                task = display.add_task("Connecting / resolving", total=None)
                previous_key: tuple[object, ...] | None = None

                def report(event: Any) -> None:
                    nonlocal previous_key
                    phase = getattr(event, "phase", "download")
                    index = getattr(event, "file_index", 1)
                    count = getattr(event, "file_count", 1)
                    key = (phase, index, getattr(event, "current_file", None))
                    description = {
                        "resolve": "Resolving assembly",
                        "download": f"Downloading {index}/{count}",
                        "extract": "Extracting FASTA",
                    }.get(phase, phase)
                    if is_query:
                        description = "Fetching sequences"
                    completed = event.completed if is_query else event.bytes_completed
                    total = event.total if is_query else event.total_bytes
                    if key != previous_key:
                        display.reset(task, description=description, total=total)
                        previous_key = key
                    display.update(task, completed=completed, total=total)

                bound.arguments[callback_parameter] = report
                try:
                    result = function(*bound.args, **bound.kwargs)
                except BaseException:
                    if status_enabled(enabled):
                        console.print(
                            f"[DNAKit] Download failed or interrupted: {target}",
                            markup=False,
                            highlight=False,
                            soft_wrap=True,
                        )
                    raise
                current = display.tasks[0]
                if current.total is None:
                    display.update(task, total=current.completed)
            destination = getattr(result, "fasta_path", target)
            if status_enabled(enabled):
                console.print(
                    f"[DNAKit] Download complete: {destination}",
                    markup=False,
                    highlight=False,
                    soft_wrap=True,
                )
            return result

        return wrapped

    return decorate
