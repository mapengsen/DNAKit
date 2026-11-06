"""Shared pytest configuration."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Generator
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError


def _missing_optional_dependency(error: BaseException) -> str | None:
    """Return a skip reason when ``error`` means an optional package or program is absent.

    The core test job installs only ``.[dev]``, so tests that exercise an optional
    backend must skip there instead of failing. DNAKit reports a missing backend as
    ``BackendUnavailableError``; a direct import reports ``ModuleNotFoundError``
    (``importlib.metadata.PackageNotFoundError`` is a subclass). A missing ``dnakit``
    module is a real failure and is never skipped.
    """

    if isinstance(error, BackendUnavailableError):
        return f"optional backend unavailable: {error}"
    if isinstance(error, ModuleNotFoundError) and error.name:
        if error.name.split(".")[0] == "dnakit":
            return None
        return f"optional dependency not installed: {error.name}"
    return None


@pytest.hookimpl(wrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[None]
) -> Generator[None, pytest.TestReport, pytest.TestReport]:
    """Report a missing optional dependency in setup or the test body as a skip."""

    report = yield
    if report.failed and call.when in ("setup", "call") and call.excinfo is not None:
        reason = _missing_optional_dependency(call.excinfo.value)
        if reason is not None:
            path, line, _ = item.location
            report.outcome = "skipped"
            report.longrepr = (path, (line or 0) + 1, f"Skipped: {reason}")
    return report


# paper/write_paper/ is a local, git-ignored folder. Some tests reuse its synthetic Fig. 2
# data; a fresh checkout (such as CI) does not have it.
_PAPER_FOLDER = Path(__file__).resolve().parents[2] / "paper" / "write_paper"


def _uses_paper_folder(module: object) -> bool:
    """Return whether a test module holds a module-level path inside the paper folder."""

    for value in vars(module).values():
        if isinstance(value, (str, Path)):
            path = Path(value)
            if path.is_absolute() and path.is_relative_to(_PAPER_FOLDER):
                return True
    return False


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip tests that read the local paper folder when it is absent."""

    if _PAPER_FOLDER.is_dir():
        return
    skip = pytest.mark.skip(reason="needs the git-ignored local paper/write_paper data")
    checked: dict[int, bool] = {}
    for item in items:
        module = getattr(item, "module", None)
        if module is None:
            continue
        if id(module) not in checked:
            checked[id(module)] = _uses_paper_folder(module)
        if checked[id(module)]:
            item.add_marker(skip)


def _snapshot() -> tuple[type[logging.Logger], list[logging.Handler], int]:
    root = logging.getLogger()
    return logging.getLoggerClass(), list(root.handlers), root.level


def _restore(state: tuple[type[logging.Logger], list[logging.Handler], int]) -> None:
    logger_class, handlers, level = state
    root = logging.getLogger()
    if logging.getLoggerClass() is not logger_class:
        logging.setLoggerClass(logger_class)
    for handler in list(root.handlers):
        if handler not in handlers:
            root.removeHandler(handler)
    root.setLevel(level)


@pytest.hookimpl(hookwrapper=True)
def pytest_make_collect_report(collector):
    """Undo global logging changes made while a test module imports optional packages.

    Some packages (for example MACS3) install a logger class and a root handler at import time.
    Importing them during collection would otherwise change every logger created afterwards,
    such as the ones PyTorch Lightning creates when it is first imported.
    """

    state = _snapshot()
    yield
    _restore(state)


@pytest.fixture(autouse=True)
def _preserve_global_logging():
    """Undo the same global logging changes when a test imports such a package lazily."""

    state = _snapshot()
    yield
    _restore(state)


@pytest.fixture(scope="session", autouse=True)
def _remove_stray_artifacts():
    """Remove files that optional packages write into the working directory during tests.

    Only artifacts that did not exist when the session started are removed.
    """

    names = ("lightning_logs", ".pybiomart.sqlite")
    root = Path.cwd()
    existed = {name for name in names if (root / name).exists()}
    yield
    for name in names:
        path = root / name
        if name in existed or not path.exists():
            continue
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        else:
            path.unlink(missing_ok=True)
