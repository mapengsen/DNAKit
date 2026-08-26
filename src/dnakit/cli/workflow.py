"""Standalone local CLI for strict JSON/YAML DNAKit workflows."""

from __future__ import annotations

import argparse
import json
import os
from collections.abc import Sequence
from pathlib import Path

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn

from dnakit.exceptions import DNAKitError
from dnakit.workflows import WorkflowProgress, WorkflowRunResult, run_workflow


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m dnakit.cli.workflow",
        description="Run a strict, bounded, local DNAKit JSON/YAML workflow.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="validate and execute one workflow")
    run.add_argument("config", type=Path, help="workflow .json/.yaml path")
    run.add_argument("--dry-run", action="store_true", help="validate without writing")
    run.add_argument("--resume", action="store_true", help="verify hashes and resume outputs")
    run.add_argument("--no-progress", action="store_true", help="disable the Rich progress bar")
    return parser


def _run_with_progress(config: Path, *, dry_run: bool, resume: bool) -> WorkflowRunResult:
    progress = Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        transient=True,
    )
    task = progress.add_task("DNAKit workflow", total=1)

    def update(event: WorkflowProgress) -> None:
        total = max(1, event.total)
        completed = event.index if event.status != "started" else max(0, event.index - 1)
        progress.update(
            task,
            total=total,
            completed=min(completed, total),
            description=f"{event.step_id}: {event.status}",
        )

    with progress:
        return run_workflow(config, dry_run=dry_run, resume=resume, progress=update)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the standalone workflow command and return a process exit status."""

    arguments = _parser().parse_args(None if argv is None else list(argv))
    console = Console(stderr=True)
    try:
        result = (
            run_workflow(
                arguments.config,
                dry_run=arguments.dry_run,
                resume=arguments.resume,
            )
            if arguments.no_progress or not os.isatty(2)
            else _run_with_progress(
                arguments.config,
                dry_run=arguments.dry_run,
                resume=arguments.resume,
            )
        )
    except (DNAKitError, OSError) as error:
        console.print(f"[red]DNAKit workflow failed:[/red] {error}")
        return 2
    print(
        json.dumps(
            result.to_dict(),
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
    )
    return 0 if result.status != "failed" else 1


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["main"]
