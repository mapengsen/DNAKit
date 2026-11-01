"""Private helpers of the GEMINI, Merqury and modkit command-line facades (integration pack P10).

``dnakit.bioinformatics._gap_cli`` runs a program and merges standard output into a bounded log. The programs of this pack
often write their table to standard output, so this module adds a runner that streams standard output into a new file with a
size limit (no shell, process group isolated), plus small readers for the tables the facades report back.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import signal
import subprocess
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import BackendExecutionError, BackendTimeoutError, ConfigurationError

PREVIEW_ROWS = 20


def stdout_to_file(
    name: str,
    arguments: Sequence[str],
    *,
    output: Path,
    cwd: Path,
    explicit: str | None,
    hint: str,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    environment: Mapping[str, str] | None = None,
    extra_expected: Sequence[Path] = (),
    program_path: str | None = None,
) -> dict[str, Any]:
    """Run one program, streaming its standard output into the new file ``output`` and keeping standard error as the log.

    ``status`` is completed when the exit code is 0 and ``output`` and every ``extra_expected`` file exist. The result has the
    same keys as ``_gap_cli.run``.
    """
    executable = program_path or _cli.program(name, explicit=explicit, hint=hint)
    if output.exists():
        raise ConfigurationError(f"{output} must not exist yet.")
    env = os.environ.copy()
    env.update({"LANG": "C", "LC_ALL": "C"})
    if environment:
        env.update(environment)
    diagnostic = output.with_name(output.name + ".stderr.tmp")
    with output.open("xb") as stdout, diagnostic.open("xb") as stderr:
        process = subprocess.Popen(
            [executable, *arguments],
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            cwd=str(cwd),
            env=env,
            shell=False,
            start_new_session=os.name == "posix",
        )
        deadline = time.monotonic() + timeout_seconds
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{name} exceeded its time limit.")
                if (
                    output.stat().st_size > _cli.MAX_ARTIFACT_BYTES
                    or diagnostic.stat().st_size > _cli.MAX_OUTPUT_BYTES
                ):
                    raise BackendExecutionError(f"{name} exceeded its output size limit.")
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
    log = diagnostic.read_text(encoding="utf-8", errors="replace")
    diagnostic.unlink()
    expected = [output, *extra_expected]
    produced = [str(path) for path in expected if path.is_file()]
    complete = process.returncode == 0 and len(produced) == len(expected)
    return {
        "status": "completed" if complete else "failed",
        "backend": name,
        "version": _cli._version_of(executable),
        "command": [name, *arguments],
        "return_code": process.returncode,
        "outputs": produced,
        "log_tail": log[-_cli.LOG_TAIL_CHARS :],
    }


def open_text(path: Path):
    """Open a text or gzip/BGZF file for reading."""
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def count_records(path: Path, *, comment: str | None = "#") -> int:
    """Number of non-empty lines that do not start with the comment prefix."""
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.strip() and not (comment and line.startswith(comment)):
                count += 1
    return count


def table_preview(
    path: Path,
    *,
    delimiter: str = "\t",
    header: bool = False,
    comment: str | None = "#",
    limit: int = PREVIEW_ROWS,
) -> dict[str, Any]:
    """Row count plus the first ``limit`` rows (lists of strings) of a delimited table; ``header`` uses the first line as names."""
    names: list[str] | None = None
    rows: list[list[str]] = []
    total = 0
    with open_text(path) as handle:
        for line in handle:
            text = line.rstrip("\n").rstrip("\r")
            if not text.strip() or (comment and text.startswith(comment)):
                continue
            fields = text.split(delimiter)
            if header and names is None:
                names = fields
                continue
            total += 1
            if len(rows) < limit:
                rows.append(fields)
    return {"columns": names, "row_count": total, "rows": rows, "rows_shown": len(rows)}


def keyword_only_text(value: object, name: str) -> str:
    """Free text (SQL, motif) that cannot be read as a flag and has no NUL or line break."""
    text = _cli.text(value, name)
    if "\n" in text or "\r" in text:
        raise ConfigurationError(f"{name} must be a single line.")
    return text


def provenance(
    tool: str, tested_version: str, version: str | None, citation: str
) -> dict[str, Any]:
    return {
        "tool": tool,
        "version": version,
        "tested_version": tested_version,
        "citation": citation,
    }
