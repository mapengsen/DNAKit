"""Private helpers of the pack-P13 command-line facades (BEDTools, BWA, MAFFT, SURVIVOR, hifiasm).

``dnakit.bioinformatics._gap_cli.run`` captures standard output in memory and merges it with standard error. Several of
these programs write their result to standard output (BEDTools, ``bwa aln``/``samse``/``sampe``/``fastmap``/``pemerge``,
MAFFT), so :func:`run_to_file` streams standard output to a new file and keeps standard error separate. Both helpers run
one fixed program with an explicit argument list (no shell) under a time limit and an output-size limit, and return the
same JSON-serialisable result dictionary as the shared layer: ``status``, ``backend``, ``version``, ``command``,
``return_code``, ``outputs`` and a bounded ``log_tail``.
"""

from __future__ import annotations

import contextlib
import os
import re
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
    DNAKitError,
)

LOG_TAIL_CHARS = _cli.LOG_TAIL_CHARS
_VERSION = re.compile(r"\d+(?:\.\d+)+")


def version_of(executable: str, flags: Sequence[str] = ("--version",)) -> str | None:
    """First dotted version number printed (to either stream, whatever the exit code) by ``executable <flags>``."""
    try:
        result = execute_bounded_command(
            executable,
            list(flags),
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
    expected: Sequence[Path],
    explicit: str | None,
    hint: str,
    version_flags: Sequence[str] | None = ("--version",),
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
    success_codes: Sequence[int] = (0,),
) -> dict[str, Any]:
    """Like ``_gap_cli.run`` but with a configurable version flag and accepted exit codes."""
    executable = _cli.program(name, explicit=explicit, hint=hint)
    result = execute_bounded_command(
        executable,
        list(arguments),
        backend_id=name,
        cwd=cwd,
        timeout_seconds=timeout_seconds,
        max_output_bytes=_cli.MAX_OUTPUT_BYTES,
        monitored_output_paths=list(expected),
        max_monitored_output_bytes=_cli.MAX_ARTIFACT_BYTES,
        environment_overrides=dict(environment) if environment is not None else None,
        verbose=False,
    )
    produced = [str(path) for path in expected if path.is_file()]
    complete = result.return_code in success_codes and len(produced) == len(expected)
    return {
        "status": "completed" if complete else "failed",
        "backend": name,
        "version": None if version_flags is None else version_of(executable, version_flags),
        "command": [name, *arguments],
        "return_code": result.return_code,
        "outputs": produced,
        "log_tail": (result.output or "")[-LOG_TAIL_CHARS:],
    }


def run_to_file(
    name: str,
    arguments: Sequence[str],
    *,
    output: Path | None,
    explicit: str | None,
    hint: str,
    cwd: Path | None = None,
    stdin_path: Path | None = None,
    extra_outputs: Sequence[Path] = (),
    version_flags: Sequence[str] | None = ("--version",),
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
    success_codes: Sequence[int] = (0,),
) -> dict[str, Any]:
    """Run one program with optional standard input from ``stdin_path``.

    With ``output`` the standard output is streamed to that new file and standard error becomes ``log_tail``. With
    ``output=None`` both streams go to ``log_tail`` (the program writes its results to ``extra_outputs``).
    """
    executable = _cli.program(name, explicit=explicit, hint=hint)
    if timeout_seconds <= 0:
        raise ConfigurationError("timeout_seconds must be positive.")
    env = os.environ.copy()
    env.update({"LANG": "C", "LC_ALL": "C"})
    if environment:
        env.update(environment)
    command = [executable, *arguments]
    if cwd is not None:
        work = cwd
    elif output is not None:
        work = output.parent
    else:
        raise ConfigurationError("cwd is required when output is None.")
    deadline = time.monotonic() + float(timeout_seconds)
    stdout_handle: Any = None
    stdin_handle: Any = subprocess.DEVNULL
    opened: list[Any] = []
    try:
        if stdin_path is not None:
            stdin_handle = stdin_path.open("rb")
            opened.append(stdin_handle)
        log_file = tempfile.TemporaryFile()  # noqa: SIM115 - closed explicitly below
        opened.append(log_file)
        if output is not None:
            stdout_handle = output.open("xb")
            opened.append(stdout_handle)
            stderr_handle: Any = log_file
        else:
            stdout_handle = log_file
            stderr_handle = subprocess.STDOUT
        process = subprocess.Popen(
            command,
            stdin=stdin_handle,
            stdout=stdout_handle,
            stderr=stderr_handle,
            cwd=str(work),
            env=env,
            shell=False,
            start_new_session=os.name == "posix",
        )
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(
                        f"{name} exceeded its time limit of {timeout_seconds} s."
                    )
                if output is not None and output.stat().st_size > _cli.MAX_ARTIFACT_BYTES:
                    raise BackendExecutionError(f"{name} exceeded the output size limit.")
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
        log_file.seek(0)
        data = log_file.read()
        log = data[-LOG_TAIL_CHARS:].decode("utf-8", errors="replace")
    finally:
        for handle in opened:
            with contextlib.suppress(OSError):
                handle.close()
    wanted = [
        path for path in ((output,) if output is not None else ()) if path is not None
    ] + list(extra_outputs)
    produced = [str(path) for path in wanted if path.is_file()]
    complete = process.returncode in success_codes and len(produced) == len(wanted)
    shown = [name, *arguments]
    if stdin_path is not None:
        shown += ["<", str(stdin_path)]
    if output is not None:
        shown += [">", str(output)]
    return {
        "status": "completed" if complete else "failed",
        "backend": name,
        "version": None if version_flags is None else version_of(executable, version_flags),
        "command": shown,
        "return_code": process.returncode,
        "outputs": produced,
        "log_tail": log,
    }


def count_lines(path: Path, *, skip_prefixes: Sequence[str] = ()) -> int:
    """Number of non-empty lines of a text output that do not start with one of ``skip_prefixes``."""
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not any(line.startswith(prefix) for prefix in skip_prefixes):
                count += 1
    return count
