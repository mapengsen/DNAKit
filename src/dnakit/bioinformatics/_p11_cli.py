"""Private helpers shared by the command-line facades of integration pack P11 (Bakta, IQ-TREE, RepeatMasker, minigraph).

Programs run without a shell through an explicit argument list. Standard output is streamed into a file (or, when the
program writes its own result files, into the log) inside a fresh output directory that is never overwritten; standard
error goes to the log file. The executable is resolved by ``dnakit.bioinformatics._gap_cli.program`` (explicit path, PATH,
the interpreter's bin directory) and a missing program raises ``BackendUnavailableError`` with an install hint.
"""

from __future__ import annotations

import contextlib
import math
import os
import signal
import subprocess
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics import _gap_cli as cli
from dnakit.exceptions import BackendExecutionError, BackendTimeoutError, ConfigurationError

LOG_NAME = "tool.log"
MAX_ARTIFACT_BYTES = 2_000_000_000
MAX_LOG_BYTES = 20_000_000
ROW_LIMIT = 1000


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def switch(arguments: list[str], option: str, value: object, name: str) -> None:
    """Append a valueless option when the Python flag is True."""
    if flag(value, name):
        arguments.append(option)


def opt_int(
    arguments: list[str],
    option: str,
    value: int | None,
    name: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> None:
    if value is not None:
        arguments.extend((option, str(cli.integer(value, name, minimum=minimum, maximum=maximum))))


def opt_num(
    arguments: list[str],
    option: str,
    value: float | None,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> None:
    if value is not None:
        arguments.extend(
            (option, format(cli.number(value, name, minimum=minimum, maximum=maximum), ".12g"))
        )


def opt_text(arguments: list[str], option: str, value: str | None, name: str) -> None:
    if value is not None:
        arguments.extend((option, cli.text(value, name)))


def opt_file(arguments: list[str], option: str, value: str | Path | None, name: str) -> None:
    if value is not None:
        arguments.extend((option, str(cli.existing_file(value, name))))


def opt_dir(arguments: list[str], option: str, value: str | Path | None, name: str) -> None:
    if value is not None:
        path = Path(value).expanduser().resolve()
        if not path.is_dir():
            raise ConfigurationError(f"{name} must be an existing directory: {path}")
        arguments.extend((option, str(path)))


def file_path(value: str | Path, name: str) -> str:
    return str(cli.existing_file(value, name))


def file_paths(values: Sequence[str | Path], name: str) -> list[str]:
    if isinstance(values, (str, Path)) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of files.")
    return [file_path(item, name) for item in values]


def token(value: str, name: str) -> str:
    """A name usable inside a file name or a model string: no whitespace-only values, no leading dash, no NUL."""
    return cli.text(value, name)


def file_stem(value: str, name: str) -> str:
    """A file-name prefix: letters, digits, underscore, dot and dash; no path separators."""
    text = cli.text(value, name)
    if not all(character.isalnum() or character in "_.-" for character in text):
        raise ConfigurationError(
            f"{name} may contain only letters, digits, '_', '.' and '-': {text!r}"
        )
    return text


# ---------------------------------------------------------------------------
# Running programs
# ---------------------------------------------------------------------------


def stream(
    executable: str,
    arguments: Sequence[str],
    directory: Path,
    *,
    stdout_name: str | None = None,
    log_name: str = LOG_NAME,
    timeout_seconds: float = cli.DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
) -> int:
    """Run one program; standard output goes to ``stdout_name`` (or to the log), standard error to the log.

    Returns the exit code. A timeout kills the process group and raises ``BackendTimeoutError``.
    """
    timeout = cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    log_path = directory / log_name
    out_path = None if stdout_name is None else directory / stdout_name
    env = {**os.environ, "LANG": "C", "LC_ALL": "C", **(dict(environment) if environment else {})}
    with contextlib.ExitStack() as stack:
        log = stack.enter_context(log_path.open("xb"))
        stdout = stack.enter_context(out_path.open("xb")) if out_path is not None else log
        try:
            process = subprocess.Popen(
                [executable, *arguments],
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=log,
                cwd=directory,
                shell=False,
                env=env,
                start_new_session=os.name == "posix",
            )
        except OSError as exc:
            raise BackendExecutionError(
                "External command could not be started.",
                code="P11_START_FAILED",
                context={"error": str(exc)},
            ) from exc
        deadline = time.monotonic() + timeout
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{Path(executable).name} exceeded its time limit.")
                if (out_path is not None and out_path.stat().st_size > MAX_ARTIFACT_BYTES) or (
                    log_path.stat().st_size > MAX_LOG_BYTES
                ):
                    raise BackendExecutionError(
                        f"{Path(executable).name} exceeded its output size limit."
                    )
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
    return int(process.returncode)


def tail(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[-cli.LOG_TAIL_CHARS :]


def run_checked(
    operation: str,
    executable: str,
    arguments: Sequence[str],
    directory: Path,
    *,
    stdout_name: str | None = None,
    timeout_seconds: float = cli.DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
) -> None:
    """Run and raise ``BackendExecutionError`` (with the log tail) when the exit code is not 0."""
    code = stream(
        executable,
        arguments,
        directory,
        stdout_name=stdout_name,
        timeout_seconds=timeout_seconds,
        environment=environment,
    )
    if code != 0:
        raise BackendExecutionError(
            f"{operation} failed (exit code {code}).",
            code="P11_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": tail(directory / LOG_NAME)},
        )


def version_of(executable: str) -> str | None:
    """The ``--version`` number of an executable (cached), or None when it cannot be read."""
    return cli._version_of(executable)


def require_outputs(operation: str, directory: Path, outputs: Mapping[str, Path]) -> None:
    missing = [name for name, path in outputs.items() if not path.is_file()]
    if missing:
        raise BackendExecutionError(
            f"{operation} finished but did not write: {', '.join(missing)}.",
            code="P11_OUTPUT_MISSING",
            context={"operation": operation, "log_tail": tail(directory / LOG_NAME)},
        )


def result(
    operation: str,
    backend: str,
    version: str | None,
    command: Sequence[str],
    directory: Path,
    outputs: Mapping[str, Path],
    *,
    tool: str,
    tested_version: str,
    citation: str,
    **extra: object,
) -> dict[str, Any]:
    """The JSON-serialisable result dictionary of one completed run."""
    payload: dict[str, Any] = {
        "status": "completed",
        "backend": backend,
        "operation": operation,
        "version": version,
        "command": list(command),
        "output_directory": str(directory),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "log_tail": tail(directory / LOG_NAME),
        "provenance": {"tool": tool, "tested_version": tested_version, "citation": citation},
    }
    payload.update(extra)
    return payload


def count_lines(path: Path, prefix: str | None = None) -> int:
    """Number of non-empty lines, or of lines starting with ``prefix``."""
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if prefix is None:
                total += 1 if line.strip() else 0
            elif line.startswith(prefix):
                total += 1
    return total


def table(path: Path, limit: int = ROW_LIMIT) -> tuple[list[list[str]], int]:
    """Tab-separated data lines (comment lines excluded): at most ``limit`` rows and the total count."""
    rows: list[list[str]] = []
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            total += 1
            if len(rows) < limit:
                rows.append(line.rstrip("\n").split("\t"))
    return rows, total


def number(value: str) -> Any:
    """``int`` or finite ``float`` when the text is a number, otherwise the text itself (``nan`` stays text)."""
    for kind in (int, float):
        try:
            parsed = kind(value)
        except ValueError:
            continue
        if isinstance(parsed, float) and not math.isfinite(parsed):
            return value
        return parsed
    return value
