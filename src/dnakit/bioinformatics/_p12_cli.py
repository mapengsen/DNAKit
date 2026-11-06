"""Private runner shared by the P12 command-line facades (Centrifuge, DeepVariant, HUMAnN, NanoPlot, seqtk).

It extends ``dnakit.bioinformatics._gap_cli`` (program resolution, argument validators) with what the P12 programs
need on top: a version probe that is configurable per program, capture of standard output into a new file with a
bounded size (several programs print their result to stdout), a hard check that the exit code is 0 and that every
declared output exists, and small parsers for tab-separated tables. Nothing is executed through a shell and no
program is ever started with an argument that the calling facade did not assemble from its own whitelist.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import re
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._gap_cli import (
    LOG_TAIL_CHARS,
    choice,
    existing_file,
    integer,
    new_directory,
    new_file,
    number,
    option,
    program,
    text,
)
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
    DNAKitError,
)

TIMEOUT_SECONDS = 86_400.0
MAX_LOG_BYTES = 5_000_000
MAX_STDOUT_BYTES = 2_000_000_000
_VERSION = re.compile(r"\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*\d)?")


@dataclass(frozen=True)
class Tool:
    """One external program: how to find it, how to read its version and what to tell a user who lacks it."""

    name: str
    executable: str
    hint: str
    tested_version: str
    citation: str
    version_arguments: tuple[str, ...] = ("--version",)
    version_pattern: str | None = None


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def optional_text(value: object | None, name: str) -> str | None:
    return None if value is None else text(value, name)


def free_text(value: object, name: str) -> str:
    """A non-empty value that may begin with '-' (a regular expression, a negative number)."""
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def existing_path(value: str | Path, name: str) -> Path:
    """An existing file or directory."""
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise ConfigurationError(f"{name} must exist: {path}")
    return path


def existing_directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def existing_files(values: object, name: str, *, minimum: int = 1) -> list[Path]:
    if isinstance(values, (str, Path)) or not isinstance(values, Sequence) or len(values) < minimum:
        raise ConfigurationError(f"{name} must be a list of at least {minimum} file path(s).")
    return [existing_file(item, name) for item in values]


def threads_value(value: object, name: str = "threads") -> int:
    return integer(value, name, minimum=1, maximum=256)


def timeout_value(value: object) -> float:
    return number(value, "timeout_seconds", minimum=1.0, maximum=TIMEOUT_SECONDS)


@lru_cache(maxsize=64)
def _version_of(executable: str, arguments: tuple[str, ...], pattern: str | None) -> str | None:
    try:
        result = execute_bounded_command(
            executable,
            list(arguments),
            backend_id=Path(executable).name,
            cwd=Path(executable).parent,
            timeout_seconds=60.0,
            max_output_bytes=200_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    expression = re.compile(pattern) if pattern else _VERSION
    for line in (result.output or "").splitlines():
        match = expression.search(line)
        if match is not None:
            return match.group(1) if expression.groups else match.group(0)
    return None


def _tail(value: str) -> str:
    return value[-LOG_TAIL_CHARS:]


def _stream_stdout(
    command: Sequence[str],
    *,
    cwd: Path,
    target: Path,
    environment: Mapping[str, str] | None,
    timeout_seconds: float,
) -> tuple[int, str]:
    """Run ``command`` with standard output written straight to a new file and standard error kept as the log."""
    env = os.environ.copy()
    env.update({"LANG": "C", "LC_ALL": "C"})
    if environment:
        env.update(environment)
    with tempfile.TemporaryFile() as errors, target.open("xb") as out:
        process = subprocess.Popen(
            list(command),
            cwd=str(cwd),
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=errors,
            shell=False,
            close_fds=True,
            start_new_session=os.name == "posix",
        )
        deadline = time.monotonic() + timeout_seconds
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{Path(command[0]).name} exceeded its time limit.")
                if (
                    target.stat().st_size > MAX_STDOUT_BYTES
                    or os.fstat(errors.fileno()).st_size > MAX_LOG_BYTES
                ):
                    raise BackendExecutionError(
                        f"{Path(command[0]).name} exceeded its output size limit."
                    )
                time.sleep(0.01)
        finally:
            if process.poll() is None:
                with contextlib.suppress(ProcessLookupError, PermissionError):
                    os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        errors.seek(max(0, os.fstat(errors.fileno()).st_size - LOG_TAIL_CHARS * 2))
        log = errors.read().decode("utf-8", errors="replace")
    return int(process.returncode), log


def invoke(
    tool: Tool,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path] = (),
    stdout_to: Path | None = None,
    explicit: str | None = None,
    timeout_seconds: float = TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
    accept_codes: Sequence[int] = (0,),
    executable_name: str | None = None,
) -> dict[str, Any]:
    """Run one program with a fixed argument vector and return a JSON-serialisable record.

    ``expected`` files must exist afterwards (``stdout_to`` is expected automatically). A failing exit code or a
    missing expected file raises ``BackendExecutionError`` carrying the log tail; a missing program raises
    ``BackendUnavailableError`` with the install hint of ``tool``.
    """
    name = executable_name or tool.executable
    executable = program(name, explicit=explicit, hint=tool.hint)
    timeout = timeout_value(timeout_seconds)
    arguments = [str(item) for item in arguments]
    wanted = list(expected)
    if stdout_to is not None:
        wanted.append(stdout_to)
        code, log = _stream_stdout(
            [executable, *arguments],
            cwd=cwd,
            target=stdout_to,
            environment=environment,
            timeout_seconds=timeout,
        )
    else:
        result = execute_bounded_command(
            executable,
            arguments,
            backend_id=name,
            cwd=cwd,
            timeout_seconds=timeout,
            max_output_bytes=MAX_LOG_BYTES,
            monitored_output_paths=[
                path for path in wanted if path.is_relative_to(cwd) and path != cwd
            ],
            max_monitored_output_bytes=MAX_STDOUT_BYTES,
            environment_overrides=dict(environment) if environment else None,
            verbose=False,
        )
        code, log = result.return_code, result.output or ""
    if code not in accept_codes:
        raise BackendExecutionError(
            f"{tool.name} exited with code {code}.",
            code=f"{tool.name.upper().replace('-', '_')}_EXECUTION_FAILED",
            context={"command": [name, *arguments], "log_tail": _tail(log)},
        )
    missing = [str(path) for path in wanted if not path.exists()]
    if missing:
        raise BackendExecutionError(
            f"{tool.name} did not write the expected output(s): {missing}",
            code=f"{tool.name.upper().replace('-', '_')}_OUTPUT_MISSING",
            context={"command": [name, *arguments], "log_tail": _tail(log)},
        )
    return {
        "status": "completed",
        "backend": tool.name,
        "version": _version_of(executable, tool.version_arguments, tool.version_pattern),
        "command": [name, *arguments],
        "return_code": code,
        "outputs": [str(path) for path in wanted],
        "log_tail": _tail(log),
        "provenance": {
            "tool": tool.name,
            "tested_version": tool.tested_version,
            "executable": executable,
            "citation": tool.citation,
        },
    }


def finish(record: dict[str, Any], operation: str, **extra: Any) -> dict[str, Any]:
    """Add the DNAKit operation name and facade-specific fields to an ``invoke`` record."""
    record["operation"] = operation
    record.update(extra)
    return record


# ---------------------------------------------------------------------------
# Small parsers
# ---------------------------------------------------------------------------


def convert_number(value: str) -> Any:
    try:
        return int(value)
    except ValueError:
        pass
    try:
        parsed = float(value)
    except ValueError:
        return value
    return parsed if parsed == parsed and parsed not in (float("inf"), float("-inf")) else value


def open_text(path: Path) -> Any:
    return (
        gzip.open(path, "rt", encoding="utf-8", errors="replace")
        if path.suffix == ".gz"
        else path.open("r", encoding="utf-8", errors="replace")
    )


def read_table(
    path: Path,
    *,
    delimiter: str = "\t",
    header: bool = True,
    names: Sequence[str] | None = None,
    comment: str | None = None,
    limit: int = 10_000,
    text_columns: frozenset[str] = frozenset(),
) -> list[dict[str, Any]]:
    """Rows of a delimited text file as dictionaries (at most ``limit`` rows; numbers converted)."""
    rows: list[dict[str, Any]] = []
    columns = list(names) if names is not None else None
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or (
                comment and line.startswith(comment) and not (header and columns is None)
            ):
                continue
            fields = line.split(delimiter)
            if columns is None:
                columns = [item.lstrip(comment or "") if comment else item for item in fields]
                if header:
                    continue
                columns = [f"column_{i + 1}" for i in range(len(fields))]
            if len(rows) >= limit:
                break
            rows.append(
                {
                    (columns[i] if i < len(columns) else f"column_{i + 1}"): (
                        value
                        if (i < len(columns) and columns[i] in text_columns)
                        else convert_number(value)
                    )
                    for i, value in enumerate(fields)
                }
            )
    return rows


def count_lines(path: Path) -> int:
    with open_text(path) as handle:
        return sum(1 for line in handle if line.strip())


__all__ = [
    "LOG_TAIL_CHARS",
    "Tool",
    "choice",
    "convert_number",
    "count_lines",
    "existing_directory",
    "existing_file",
    "existing_files",
    "existing_path",
    "finish",
    "flag",
    "free_text",
    "integer",
    "invoke",
    "new_directory",
    "new_file",
    "number",
    "open_text",
    "option",
    "optional_text",
    "program",
    "read_table",
    "text",
    "threads_value",
    "timeout_value",
]
