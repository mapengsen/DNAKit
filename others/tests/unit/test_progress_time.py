"""Time display, lazy tracking, and global quiet-mode contracts."""

from io import StringIO

from rich.console import Console
from rich.progress import SpinnerColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn

from dnakit._progress import Progress, track
from dnakit._verbosity import _enabled


def test_known_total_has_eta_and_phase_reset_discards_old_speed() -> None:
    output = StringIO()
    now = [0.0]
    with Progress(
        console=Console(file=output, width=160), auto_refresh=False, get_time=lambda: now[0]
    ) as progress:
        task_id = progress.add_task("Downloading", total=100)
        now[0] = 1
        progress.update(task_id, completed=10)
        now[0] = 2
        progress.update(task_id, completed=20)
        task = progress.tasks[task_id]
        assert task.elapsed == 2
        assert task.time_remaining is not None and task.time_remaining > 0
        progress.refresh()
        progress.reset(task_id, description="Extracting", total=200)
        assert task.completed == 0
        assert task.time_remaining is None
    assert "elapsed" in output.getvalue() and "ETA" in output.getvalue()


def test_unknown_total_shows_elapsed_without_guessing_remaining() -> None:
    now = [0.0]
    with Progress(
        SpinnerColumn(),
        TextColumn("Waiting for backend"),
        auto_refresh=False,
        console=Console(file=StringIO()),
        get_time=lambda: now[0],
    ) as progress:
        task_id = progress.add_task("backend", total=None)
        now[0] = 12
        progress.advance(task_id)
        task = progress.tasks[task_id]
        assert task.elapsed == 12
        assert task.total is None and task.time_remaining is None
        assert any(isinstance(column, TimeElapsedColumn) for column in progress.columns)


def test_existing_time_columns_are_not_duplicated() -> None:
    progress = Progress(
        TextColumn("elapsed"),
        TimeElapsedColumn(),
        TextColumn("ETA"),
        TimeRemainingColumn(),
        disable=True,
    )
    assert len(progress.columns) == 4


def test_tracking_does_not_materialize_unknown_stream() -> None:
    consumed = []

    def source():
        for value in range(10):
            consumed.append(value)
            yield value

    stream = track(source(), console=Console(file=StringIO()), disable=True)
    assert consumed == []
    assert next(stream) == 0 and consumed == [0]
    assert list(stream) == list(range(1, 10))


def test_global_quiet_mode_suppresses_tracker_output() -> None:
    output = StringIO()
    token = _enabled.set(False)
    try:
        assert list(track([1, 2, 3], console=Console(file=output))) == [1, 2, 3]
    finally:
        _enabled.reset(token)
    assert output.getvalue() == ""
