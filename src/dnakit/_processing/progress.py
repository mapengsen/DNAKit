"""Byte-weighted split progress with stage-specific time estimates."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from threading import Event, Thread
from time import monotonic
from typing import Any

from rich.console import Console
from rich.progress import (
    BarColumn,
    DownloadColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

from dnakit._io_progress import read_observer
from dnakit._progress import Progress


class SplitProgress:
    def __init__(self, sources: tuple[Path, ...], enabled: bool):
        self.enabled = enabled
        self.sizes = {str(path): path.stat().st_size for path in sources}
        self.positions: dict[str, int] = {}
        self.display = Progress(
            TextColumn("{task.description}", markup=False),
            BarColumn(),
            TaskProgressColumn(),
            DownloadColumn(),
            TransferSpeedColumn(),
            TextColumn("elapsed"),
            TimeElapsedColumn(),
            TextColumn("ETA"),
            TimeRemainingColumn(),
            console=Console(stderr=True),
            disable=not enabled,
        )
        self.read_task = self.display.add_task(
            "Reading split inputs", total=sum(self.sizes.values())
        )

    def __enter__(self) -> SplitProgress:
        self.display.start()
        return self

    def __exit__(self, *args: Any) -> None:
        self.display.stop()

    def read(self, source: str, completed: int) -> None:
        if source in self.sizes:
            self.positions[source] = min(self.sizes[source], completed)
            self.display.update(self.read_task, completed=sum(self.positions.values()))

    @contextmanager
    def reading(self) -> Any:
        self.positions.clear()
        self.display.reset(self.read_task, total=sum(self.sizes.values()))
        token = read_observer.set(self.read if self.enabled else None)
        try:
            yield
        finally:
            read_observer.reset(token)

    def read_done(self) -> None:
        self.display.update(self.read_task, completed=sum(self.sizes.values()))

    @contextmanager
    def writing(self, jobs: list[Any]) -> Any:
        """Poll tiny per-worker counters; workers never render competing progress bars."""
        tasks = []
        for records, target, _config, counter in jobs:
            weights = record_weights(records)
            task = self.display.add_task(f"Writing {target.name}", total=sum(weights))
            tasks.append((task, counter, target.name, sum(weights)))
        stop = Event()
        states: dict[int, str] = {}
        starts = {task: self.display.tasks[task].start_time for task, *_ in tasks}

        def refresh() -> None:
            for task, counter, name, total in tasks:
                try:
                    value, state = counter.read_text().split()
                    if states.get(task) != state:
                        if state == "checking":
                            self.display.reset(task, total=None, description=f"Checking {name}")
                            self.display.tasks[task].start_time = starts[task]
                        elif state == "done":
                            self.display.update(task, total=total, description=f"Saved {name}")
                        states[task] = state
                    if state != "checking":
                        self.display.update(task, completed=int(value))
                except (OSError, ValueError):
                    continue

        def monitor() -> None:
            while not stop.wait(0.2):
                refresh()

        thread = Thread(target=monitor, name="dnakit-split-progress", daemon=True)
        if self.enabled:
            thread.start()
        try:
            yield
        finally:
            stop.set()
            if self.enabled:
                thread.join()
            refresh()


def record_weights(records: Any) -> tuple[int, ...]:
    # Indexed records: (path, start, end). Spool records: (offset, serialized length).
    return tuple(
        location[2] - location[1] if len(location) == 3 else location[1]
        for location in records.locations
    )


def tracked_records(records: Any, counter: Path) -> Any:
    completed = 0
    reported, reported_at = 0, monotonic()
    counter.write_text("0 writing")
    for record, weight in zip(records, record_weights(records), strict=True):
        yield record
        completed += weight
        now = monotonic()
        if completed - reported >= 1024**2 or now - reported_at >= 0.2:
            counter.write_text(f"{completed} writing")
            reported, reported_at = completed, now
    counter.write_text(f"{completed} checking")
