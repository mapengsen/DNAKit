"""Private runner shared by the pack-P14 command-line facades (BLAST+, BamTools, HISAT2, Kraken 2, Primer3, SPAdes,
WhatsHap, fastp).

It extends ``dnakit.bioinformatics._gap_cli`` with the two things those programs need: a program-specific version
flag (BLAST+ uses ``-version``, SPAdes ``--version`` on ``spades.py`` ...) and a ``stdout_file`` for programs that
print their result to standard output (``hisat2-inspect``, ``bamtools count``, ``kraken2-inspect``, ``ntdpal``).
Standard output goes to a new file, standard error is kept as the bounded log. No shell is used; the argument list is
built by the public function from typed, validated parameters.
"""

from __future__ import annotations

import os
import re
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from contextlib import suppress
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
    DNAKitError,
)

LOG_TAIL_CHARS = _gap_cli.LOG_TAIL_CHARS
DEFAULT_TIMEOUT_SECONDS = _gap_cli.DEFAULT_TIMEOUT_SECONDS
_VERSION = re.compile(r"\d+(?:\.\d+)+")


@lru_cache(maxsize=64)
def version_of(executable: str, version_args: tuple[str, ...] = ("--version",)) -> str | None:
    """First dotted number printed by ``executable <version_args>``, or None when the program does not print one."""
    try:
        result = execute_bounded_command(
            executable,
            list(version_args),
            backend_id=Path(executable).name,
            cwd=Path(executable).parent,
            timeout_seconds=60.0,
            max_output_bytes=100_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    match = _VERSION.search(result.output or "")
    return match.group(0) if match is not None else None


def run(
    name: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path] = (),
    explicit: str | None,
    hint: str,
    version_args: Sequence[str] = ("--version",),
    stdout_file: Path | None = None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
    program_name: str | None = None,
) -> dict[str, Any]:
    """Run one fixed program; ``status`` is ``completed`` only when the exit code is 0 and every expected file exists.

    ``name`` is the executable that is looked up (``program_name`` overrides it for the result's backend label).
    With ``stdout_file`` the standard output is streamed into that new file (listed in ``outputs``) and the log
    tail holds standard error only.
    """
    executable = _gap_cli.program(name, explicit=explicit, hint=hint)
    version = version_of(executable, tuple(version_args))
    label = program_name or name
    if stdout_file is None:
        result = execute_bounded_command(
            executable,
            list(arguments),
            backend_id=label,
            cwd=cwd,
            timeout_seconds=timeout_seconds,
            max_output_bytes=_gap_cli.MAX_OUTPUT_BYTES,
            monitored_output_paths=list(expected),
            max_monitored_output_bytes=_gap_cli.MAX_ARTIFACT_BYTES,
            environment_overrides=dict(environment) if environment is not None else None,
            verbose=False,
        )
        return_code = result.return_code
        log = result.output or ""
        files = list(expected)
    else:
        return_code, log = _run_to_stdout_file(
            executable,
            arguments,
            cwd=cwd,
            stdout_file=stdout_file,
            timeout_seconds=timeout_seconds,
            environment=environment,
            label=label,
        )
        files = [stdout_file, *expected]
    produced = [str(path) for path in files if path.is_file()]
    complete = return_code == 0 and len(produced) == len(files)
    return {
        "status": "completed" if complete else "failed",
        "backend": label,
        "version": version,
        "command": [label, *arguments],
        "return_code": return_code,
        "outputs": produced,
        "log_tail": log[-LOG_TAIL_CHARS:],
    }


def _run_to_stdout_file(
    executable: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    stdout_file: Path,
    timeout_seconds: float,
    environment: Mapping[str, str] | None,
    label: str,
) -> tuple[int, str]:
    env = os.environ.copy()
    env.update({"LANG": "C", "LC_ALL": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    if environment:
        env.update(environment)
    with stdout_file.open("xb") as out, tempfile.TemporaryFile() as err:
        process = subprocess.Popen(
            [executable, *arguments],
            cwd=str(cwd),
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=err,
            shell=False,
            start_new_session=os.name == "posix",
        )
        deadline = time.monotonic() + timeout_seconds
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{label} exceeded its time limit.")
                if stdout_file.stat().st_size > _gap_cli.MAX_ARTIFACT_BYTES:
                    raise BackendExecutionError(f"{label} exceeded its output size limit.")
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                with suppress(ProcessLookupError, PermissionError):
                    os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        err.seek(0)
        log = err.read(_gap_cli.MAX_OUTPUT_BYTES).decode("utf-8", errors="replace")
    return process.returncode, log


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def nucleotides(value: object, name: str, *, alphabet: str = "ACGTNacgtn") -> str:
    """Non-empty text made only of nucleotide letters."""
    if (
        not isinstance(value, str)
        or not value
        or any(character not in alphabet for character in value)
    ):
        raise ConfigurationError(f"{name} must be a non-empty string of the letters {alphabet}.")
    return value


def existing_directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def line_text(value: object, name: str) -> str:
    """Single-line text without NUL; may start with '-' (use it only where a flag value is always supplied)."""
    if not isinstance(value, str) or not value or "\x00" in value or "\n" in value or "\r" in value:
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def tsv_rows(path: Path, *, limit: int = 100_000) -> list[list[str]]:
    """Tab-separated rows of a small output file (comment lines starting with '#' are skipped)."""
    rows: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
            if len(rows) >= limit:
                break
    return rows
