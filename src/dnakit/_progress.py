"""Consistent elapsed time and stage-local ETA for existing console progress."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sized
from typing import Any, TypeVar

from rich.console import Console
from rich.progress import Progress as RichProgress
from rich.progress import TextColumn, TimeElapsedColumn, TimeRemainingColumn

from dnakit._verbosity import progress_enabled

_T = TypeVar("_T")


class Progress(RichProgress):
    """Keep caller-selected units and add missing time columns exactly once."""

    def __init__(self, *columns: Any, **kwargs: Any) -> None:
        selected = list(columns or self.get_default_columns())
        if not any(isinstance(column, TimeElapsedColumn) for column in selected):
            index = next(
                (i for i, column in enumerate(selected) if isinstance(column, TimeRemainingColumn)),
                len(selected),
            )
            selected.insert(index, TimeElapsedColumn())
        if not any(isinstance(column, TimeRemainingColumn) for column in selected):
            selected.append(TimeRemainingColumn())
        labeled: list[Any] = []
        for column in selected:
            label = (
                "elapsed"
                if isinstance(column, TimeElapsedColumn)
                else "ETA"
                if isinstance(column, TimeRemainingColumn)
                else None
            )
            if label is not None and (
                not labeled or getattr(labeled[-1], "text_format", None) != label
            ):
                labeled.append(TextColumn(label))
            labeled.append(column)
        kwargs["disable"] = not progress_enabled(not kwargs.get("disable", False))
        if kwargs.get("console") is None:
            kwargs["console"] = Console(stderr=True)
        super().__init__(*labeled, **kwargs)


def track(
    sequence: Iterable[_T],
    description: str = "Working...",
    total: float | None = None,
    *,
    disable: bool = False,
    console: Console | None = None,
) -> Iterator[_T]:
    """Track sized inputs without pre-consuming streams to guess their total."""
    if total is None and isinstance(sequence, Sized):
        total = len(sequence)
    with Progress(console=console, disable=disable) as progress:
        task = progress.add_task(description, total=total)
        for item in sequence:
            yield item
            progress.advance(task)
