"""Shared pytest configuration."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pytest


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
