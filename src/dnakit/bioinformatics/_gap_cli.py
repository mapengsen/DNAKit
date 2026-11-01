"""Private helpers shared by the command-line facades of the DNA-function gap batches.

A facade runs one fixed program with an explicit argument list (no shell) through
``dnakit.backends.external.execute_bounded_command``. The program is resolved from an explicit path, then from
PATH, then from the directory of the running interpreter (where pip installs console scripts). A missing program
raises ``BackendUnavailableError`` carrying the facade's install hint. Facade results are JSON-serialisable
dictionaries: status, backend, version, command, return code, existing outputs and a bounded log tail.
"""

from __future__ import annotations

import re
import shutil
import sys
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError

LOG_TAIL_CHARS = 4000
DEFAULT_TIMEOUT_SECONDS = 3600.0
MAX_OUTPUT_BYTES = 5_000_000
MAX_ARTIFACT_BYTES = 2_000_000_000
_VERSION = re.compile(r"\d+(?:\.\d+)+")


def program(name: str, *, explicit: str | None, hint: str) -> str:
    """Resolve an executable: explicit path, then PATH, then the interpreter's own bin directory."""
    if explicit is not None:
        candidate = Path(explicit).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name} executable does not exist: {candidate}", hint=hint)
        return str(candidate.resolve())
    found = shutil.which(name)
    if found is None:
        sibling = Path(sys.executable).parent / name
        if sibling.is_file():
            return str(sibling.resolve())
        raise BackendUnavailableError(f"{name} was not found on PATH.", hint=hint)
    return str(Path(found).resolve())


@lru_cache(maxsize=64)
def _version_of(executable: str) -> str | None:
    try:
        result = execute_bounded_command(
            executable,
            ["--version"],
            backend_id=Path(executable).name,
            cwd=Path(executable).parent,
            timeout_seconds=60.0,
            max_output_bytes=100_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    for line in (result.output or "").splitlines():
        match = _VERSION.search(line)
        if match is not None:
            return match.group(0)
    return None


def run(
    name: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    explicit: str | None,
    hint: str,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Run one fixed program; ``status`` is completed only when the exit code is 0 and every expected file exists."""
    executable = program(name, explicit=explicit, hint=hint)
    result = execute_bounded_command(
        executable,
        list(arguments),
        backend_id=name,
        cwd=cwd,
        timeout_seconds=timeout_seconds,
        max_output_bytes=MAX_OUTPUT_BYTES,
        monitored_output_paths=list(expected),
        max_monitored_output_bytes=MAX_ARTIFACT_BYTES,
        environment_overrides=dict(environment) if environment is not None else None,
        verbose=False,
    )
    produced = [str(path) for path in expected if path.is_file()]
    complete = result.return_code == 0 and len(produced) == len(expected)
    return {
        "status": "completed" if complete else "failed",
        "backend": name,
        "version": _version_of(executable),
        "command": [name, *arguments],
        "return_code": result.return_code,
        "outputs": produced,
        "log_tail": (result.output or "")[-LOG_TAIL_CHARS:],
    }


def existing_file(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def new_file(value: str | Path, name: str) -> Path:
    """A path that does not exist yet and whose folder exists; facades never overwrite outputs."""
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def new_directory(value: str | Path, name: str) -> Path:
    """Create a fresh output directory; an existing directory is refused."""
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    path.mkdir(parents=True)
    return path


def text(value: object, name: str) -> str:
    """Non-empty text that cannot be read as a command-line flag."""
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if value.startswith("-") or "\x00" in value:
        raise ConfigurationError(f"{name} must not start with '-' or contain NUL characters.")
    return value


def integer(value: object, name: str, *, minimum: int | None = None, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def number(value: object, name: str, *, minimum: float | None = None, maximum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return float(value)


def choice(value: object, options: Sequence[str], name: str) -> str:
    if value not in options:
        raise ConfigurationError(f"{name} must be one of {', '.join(options)}.")
    return str(value)


def option(flag: str, value: object | None) -> list[str]:
    """``[flag, str(value)]`` when a value is given, otherwise nothing (optional fixed flags)."""
    return [] if value is None else [flag, str(value)]
