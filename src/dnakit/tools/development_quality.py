"""Format and check reviewable Python candidates using DNAKit's development tools.

This version adds explicit exports, Ruff formatting/lint and strict Mypy checks.
No dependency installer is invoked; missing checkers leave candidates unactivated.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from dnakit._progress import track
from dnakit.exceptions import ConfigurationError


def _command(tool: str, root: Path) -> list[str]:
    if importlib.util.find_spec(tool) is not None:
        return [sys.executable, "-m", tool]
    executable = shutil.which(tool)
    if executable:
        return [executable]
    # A source checkout may keep its development CLI in a local virtualenv.
    for directory in (root, Path(__file__).resolve().parents[3]):
        for relative in (f".venv/bin/{tool}", f".venv/Scripts/{tool}.exe"):
            candidate = directory / relative
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return [str(candidate)]
    raise ConfigurationError(
        f"The {tool} checker is required before activating contribution-ready code.",
        code="TOOL_DEVELOPMENT_CHECKER_MISSING",
        hint="Install the DNAKit dev extra in the development environment explicitly.",
    )


def _run(
    command: list[str], *, root: Path, directory: Path, label: str, source: str | None = None
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    source_paths = [root / "src"]
    package_source = Path(__file__).resolve().parents[2]
    # Installed typed packages are discovered by Mypy itself. Adding the whole
    # site-packages directory to MYPYPATH is forbidden by Mypy.
    if package_source.name == "src" and package_source not in source_paths:
        source_paths.append(package_source)
    environment["MYPYPATH"] = os.pathsep.join(str(path) for path in source_paths)
    try:
        result = subprocess.run(
            command,
            input=source,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=directory,
            env=environment,
            timeout=60,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ConfigurationError(
            f"Candidate {label} check could not finish: {type(exc).__name__}.",
            code="TOOL_DEVELOPMENT_QUALITY_FAILED",
        ) from None
    (directory / f"{label}.log").write_text(result.stdout[-100_000:], encoding="utf-8")
    return result


def prepare_generated_source(source: str, name: str, description: str) -> str:
    """Add package-level documentation and exactly one explicit public export."""
    try:
        module = ast.parse(source)
    except SyntaxError as exc:
        raise ConfigurationError(
            "Candidate source is not valid Python.", code="GENERATED_TOOL_INVALID"
        ) from exc
    declarations = [
        node
        for node in module.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        and any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in (node.targets if isinstance(node, ast.Assign) else [node.target])
        )
    ]
    if declarations:
        if len(declarations) != 1:
            raise ConfigurationError(
                "Only one __all__ declaration is allowed.", code="GENERATED_TOOL_INVALID"
            )
        try:
            exports = ast.literal_eval(declarations[0].value)  # type: ignore[arg-type]
        except (ValueError, TypeError):
            exports = None
        if exports != [name] and exports != (name,):
            raise ConfigurationError(
                "__all__ must export only the generated function.", code="GENERATED_TOOL_INVALID"
            )
    else:
        source = source.rstrip() + f"\n\n__all__ = [{name!r}]\n"
    if not ast.get_docstring(module):
        source = json.dumps(description, ensure_ascii=False) + "\n\n" + source
    return source


def check_python_quality(
    source: str, *, root: Path, directory: Path, name: str, verbose: bool = True
) -> tuple[str, dict[str, Any]]:
    """Return formatted source only when syntax, style and strict typing pass."""
    if not isinstance(name, str) or not name.isidentifier():
        raise ConfigurationError(
            "Candidate name must be a Python identifier.", code="GENERATED_TOOL_INVALID"
        )
    directory.mkdir(parents=True, exist_ok=True)
    candidate = directory / f"{name}.py"
    ruff, mypy = _command("ruff", root), _command("mypy", root)
    formatted = _run(
        [
            *ruff,
            "format",
            "--isolated",
            "--line-length",
            "100",
            "--target-version",
            "py310",
            "--stdin-filename",
            str(candidate),
            "-",
        ],
        root=root,
        directory=directory,
        label="format",
        source=source,
    )
    if formatted.returncode:
        raise ConfigurationError(
            "Candidate formatting failed.", code="TOOL_DEVELOPMENT_QUALITY_FAILED"
        )
    candidate.write_text(formatted.stdout, encoding="utf-8")
    checks = [
        (
            "ruff",
            [
                *ruff,
                "check",
                "--isolated",
                "--select",
                "B,E,F,I,RUF,SIM,UP",
                "--line-length",
                "100",
                "--target-version",
                "py310",
                str(candidate),
            ],
        ),
        (
            "mypy",
            [
                *mypy,
                "--strict",
                "--follow-imports=silent",
                "--python-version",
                "3.10",
                "--cache-dir",
                str(directory / ".mypy_cache"),
                str(candidate),
            ],
        ),
    ]
    report: dict[str, Any] = {"formatted": True, "checks": {}}
    for label, command in track(
        checks, description="Check DNAKit contribution", disable=not verbose
    ):
        result = _run(command, root=root, directory=directory, label=label)
        report["checks"][label] = {
            "passed": result.returncode == 0,
            "log": str(directory / f"{label}.log"),
        }
        if result.returncode:
            raise ConfigurationError(
                f"Candidate failed {label}: {result.stdout[-12_000:]}",
                code="TOOL_DEVELOPMENT_QUALITY_FAILED",
            )
    return candidate.read_text(encoding="utf-8"), report
