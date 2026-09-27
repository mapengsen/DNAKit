"""Real contribution checks plus explicit checker/timeout failure behavior."""

from __future__ import annotations

import ast
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path
from unittest.mock import Mock

import pytest

import dnakit.tools.development_quality as quality
from dnakit.exceptions import ConfigurationError

ROOT = Path(__file__).resolve().parents[3]


def _require_checkers() -> None:
    try:
        quality._command("ruff", ROOT)
        quality._command("mypy", ROOT)
    except ConfigurationError:
        pytest.skip("Actual contribution checks require Ruff and Mypy")


def test_real_ruff_formats_candidate_and_strict_mypy_accepts_types(tmp_path: Path) -> None:
    _require_checkers()
    source = quality.prepare_generated_source(
        "def metric(value:int,*,offset:int=0)->int:\n return value+offset\n",
        "metric",
        "Return the supplied value plus an optional offset.",
    )

    formatted, report = quality.check_python_quality(
        source,
        root=ROOT,
        directory=tmp_path,
        name="metric",
        verbose=False,
    )

    assert "def metric(value: int, *, offset: int = 0) -> int:" in formatted
    assert "    return value + offset\n" in formatted
    assert (
        ast.get_docstring(ast.parse(formatted))
        == "Return the supplied value plus an optional offset."
    )
    assert (tmp_path / "metric.py").read_text(encoding="utf-8") == formatted
    assert report["formatted"] and set(report["checks"]) == {"ruff", "mypy"}
    assert all(item["passed"] for item in report["checks"].values())
    assert all(Path(item["log"]).is_file() for item in report["checks"].values())


def test_real_strict_mypy_rejects_an_incorrect_return_type(tmp_path: Path) -> None:
    _require_checkers()
    source = quality.prepare_generated_source(
        'def metric(value: int) -> int:\n    return "wrong type"\n',
        "metric",
        "Count a value.",
    )

    with pytest.raises(ConfigurationError, match="TOOL_DEVELOPMENT_QUALITY_FAILED"):
        quality.check_python_quality(
            source, root=ROOT, directory=tmp_path, name="metric", verbose=False
        )

    assert "All checks passed" in (tmp_path / "ruff.log").read_text(encoding="utf-8")
    assert "return-value" in (tmp_path / "mypy.log").read_text(encoding="utf-8")


def test_missing_checkers_fail_without_installing_or_running_a_subprocess(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    monkeypatch.setattr(os, "access", lambda path, mode: False)
    run = Mock()
    monkeypatch.setattr(subprocess, "run", run)

    with pytest.raises(ConfigurationError, match="TOOL_DEVELOPMENT_CHECKER_MISSING"):
        quality.check_python_quality(
            "def metric(): pass",
            root=tmp_path,
            directory=tmp_path / "audit",
            name="metric",
            verbose=False,
        )

    run.assert_not_called()


def test_checker_timeout_is_reported_as_a_quality_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(quality, "_command", lambda tool, root: [tool])
    monkeypatch.setattr(subprocess, "run", Mock(side_effect=subprocess.TimeoutExpired("ruff", 60)))
    with pytest.raises(ConfigurationError, match="TOOL_DEVELOPMENT_QUALITY_FAILED"):
        quality.check_python_quality(
            "def metric(): pass", root=tmp_path, directory=tmp_path, name="metric", verbose=False
        )


@pytest.mark.parametrize(
    "declaration",
    [
        "__all__ = ['other']",
        "__all__ = ['metric']\n__all__ = ['metric']",
        "__all__ = dynamically_computed()",
    ],
)
def test_explicit_exports_must_identify_only_the_generated_function(declaration: str) -> None:
    with pytest.raises(ConfigurationError, match="GENERATED_TOOL_INVALID"):
        quality.prepare_generated_source(
            f"def metric():\n    pass\n{declaration}\n", "metric", "Metric."
        )


def test_existing_documentation_future_imports_and_exports_are_preserved() -> None:
    source = (
        '"""Existing module documentation."""\nfrom __future__ import annotations\n\n'
        'def metric(value: int) -> int:\n    return value\n\n__all__ = ["metric"]\n'
    )
    assert quality.prepare_generated_source(source, "metric", "Replacement description") == source


def test_candidate_filename_cannot_escape_its_audit_directory(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="GENERATED_TOOL_INVALID"):
        quality.check_python_quality(
            "def metric(): pass",
            root=tmp_path,
            directory=tmp_path / "audit",
            name="../escape",
            verbose=False,
        )
    assert not (tmp_path / "escape.py").exists()
