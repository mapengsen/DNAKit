"""Private runner shared by the AUGUSTUS facades of integration pack P04.

Every AUGUSTUS program and helper script is run as one fixed executable with an explicit argument list (no shell).
Standard input can come from a file and standard output can go to a new file, because many AUGUSTUS scripts are
filters (``script < in > out``). The executable is resolved from an explicit ``augustus_bin_dir`` (the folder that
holds ``augustus`` and its scripts), then from PATH, then from the folder of the running interpreter. The child
process gets that folder in front of PATH (so ``#!/usr/bin/env perl`` and the sibling programs resolve), the folder
in ``PERL5LIB`` (the scripts load ``SplicedAlignment.pm`` and ``helpMod.pm`` from it) and ``AUGUSTUS_CONFIG_PATH``.
"""

from __future__ import annotations

import contextlib
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
    DNAKitError,
)

TESTED_VERSION = "3.5.0"
CITATION = "https://doi.org/10.1093/nar/gkl200"
INSTALL_HINT = (
    "Install AUGUSTUS 3.5.0 or later (for example `conda install -c bioconda augustus`, or build it from "
    "https://github.com/Gaius-Augustus/Augustus), then put its bin folder on PATH or pass augustus_bin_dir. "
    "Most helper scripts need Perl; some need Python 3, samtools or the UCSC tools."
)
TIMEOUT_SECONDS = 24 * 3600.0
VERSION_TIMEOUT_SECONDS = 60.0
MAX_LOG_BYTES = 5_000_000
MAX_CAPTURE_BYTES = 50_000_000
MAX_ARTIFACT_BYTES = 2_000_000_000
LOG_TAIL_CHARS = 4000
PREVIEW_LINES = 10
_VERSION_LINE = re.compile(r"\((\d+(?:\.\d+)+)\)")
_VERSIONS: dict[str, str | None] = {}

NAME = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.+-]*$")


def bin_dir(value: object) -> str | None:
    """Validate ``augustus_bin_dir``: ``None`` (use PATH) or an existing folder."""
    if value is None or value == "":
        return None
    path = Path(str(value)).expanduser()
    if not path.is_dir():
        raise BackendUnavailableError(
            f"augustus_bin_dir is not a folder: {path}", hint=INSTALL_HINT
        )
    return str(path.resolve())


def executable(tool: str, directory: str | None) -> str:
    explicit = None if directory is None else str(Path(directory) / tool)
    return cli.program(tool, explicit=explicit, hint=INSTALL_HINT)


def config_directory(value: object, exe: str) -> str | None:
    """The AUGUSTUS_CONFIG_PATH for a run: explicit, then the environment, then ``<bin>/../config``."""
    if value is not None:
        path = Path(str(value)).expanduser()
        if not path.is_dir():
            raise ConfigurationError(f"config_path must be an existing folder: {path}")
        return str(path.resolve())
    if os.environ.get("AUGUSTUS_CONFIG_PATH"):
        return None
    sibling = Path(exe).resolve().parent.parent / "config"
    return str(sibling) if (sibling / "species").is_dir() else None


def version(directory: str | None, exe: str) -> str | None:
    """Version of the AUGUSTUS installation that ``exe`` belongs to (``augustus --version``)."""
    parent = str(Path(exe).resolve().parent)
    if parent not in _VERSIONS:
        _VERSIONS[parent] = None
        try:
            main = executable("augustus", parent)
            result = execute_bounded_command(
                main,
                ("--version",),
                backend_id="augustus",
                cwd=Path(main).parent,
                timeout_seconds=VERSION_TIMEOUT_SECONDS,
                max_output_bytes=100_000,
                verbose=False,
            )
            match = _VERSION_LINE.search(result.output or "")
            _VERSIONS[parent] = match.group(1) if match else None
        except DNAKitError:
            _VERSIONS[parent] = None
    return _VERSIONS[parent]


def tail(text: str) -> str:
    return text[-LOG_TAIL_CHARS:]


@dataclass(frozen=True)
class Completed:
    """Result of one bounded run."""

    tool: str
    return_code: int
    log: str
    executable: str
    version: str | None
    command: tuple[str, ...]
    stdout_path: Path | None


def run(
    tool: str,
    arguments: Sequence[str],
    *,
    directory: str | None,
    cwd: Path,
    stdin: Path | None = None,
    stdout: Path | None = None,
    capture: bool = False,
    config_path: str | None = None,
    accept: Sequence[int] = (0,),
    expected: Sequence[Path] = (),
    environment: Mapping[str, str] | None = None,
    timeout_seconds: float = TIMEOUT_SECONDS,
) -> Completed:
    """Run one fixed program.

    Standard output goes to ``stdout`` (a new file this facade owns) or, with ``capture=True``, to a temporary file
    that the caller reads through ``Completed.stdout_path`` (and must delete with :func:`captured_text`). Files in
    ``expected`` must exist afterwards. ``cwd`` is the folder the program runs in; programs that write files next
    to their working directory write them there.
    """
    exe = executable(tool, directory)
    program_version = version(directory, exe)
    config = config_directory(config_path, exe)
    parent = str(Path(exe).resolve().parent)
    env = dict(os.environ)
    env["PATH"] = parent + os.pathsep + env.get("PATH", "")
    env["PERL5LIB"] = parent + (os.pathsep + env["PERL5LIB"] if env.get("PERL5LIB") else "")
    if config is not None:
        env["AUGUSTUS_CONFIG_PATH"] = config
    if environment:
        env.update(environment)
    command = (exe, *arguments)
    created = stdout is not None
    with tempfile.TemporaryDirectory(prefix="dnakit-augustus-") as scratch:
        scratch_dir = Path(scratch)
        target = stdout if stdout is not None else scratch_dir / "stdout"
        diagnostic = scratch_dir / "stderr"
        limit = MAX_ARTIFACT_BYTES if created else MAX_CAPTURE_BYTES
        input_handle = stdin.open("rb") if stdin is not None else None
        try:
            with target.open("xb") as out_handle, diagnostic.open("xb") as err_handle:
                process = subprocess.Popen(
                    command,
                    stdin=input_handle if input_handle is not None else subprocess.DEVNULL,
                    stdout=out_handle,
                    stderr=err_handle,
                    cwd=cwd,
                    env=env,
                    shell=False,
                    start_new_session=os.name == "posix",
                )
                deadline = time.monotonic() + timeout_seconds
                try:
                    while process.poll() is None:
                        if time.monotonic() > deadline:
                            raise BackendTimeoutError(f"{tool} exceeded its time limit.")
                        if (
                            target.stat().st_size > limit
                            or diagnostic.stat().st_size > MAX_LOG_BYTES
                        ):
                            raise BackendExecutionError(f"{tool} exceeded its output size limit.")
                        time.sleep(0.02)
                finally:
                    if process.poll() is None:
                        if os.name == "posix":
                            with contextlib.suppress(ProcessLookupError):
                                os.killpg(process.pid, signal.SIGKILL)
                        else:
                            process.kill()
                    process.wait()
        except BaseException:
            if created and stdout is not None:
                with contextlib.suppress(OSError):
                    stdout.unlink()
            raise
        finally:
            if input_handle is not None:
                input_handle.close()
        log = diagnostic.read_text(encoding="utf-8", errors="replace")
        failed = process.returncode not in accept
        missing = [str(path) for path in expected if not path.exists()]
        if failed or missing:
            if created and stdout is not None:
                with contextlib.suppress(OSError):
                    stdout.unlink()
            raise BackendExecutionError(
                f"AUGUSTUS {tool} exited with code {process.returncode}"
                + (f" and did not write {missing}." if missing and not failed else "."),
                code="AUGUSTUS_EXECUTION_FAILED",
                context={"tool": tool, "log_tail": tail(log)},
            )
        kept: Path | None = stdout
        if capture:
            descriptor, copy_name = tempfile.mkstemp(prefix="dnakit-augustus-out-")
            os.close(descriptor)
            kept = Path(copy_name)
            shutil.copyfile(target, kept)
    return Completed(
        tool, int(process.returncode), log, exe, program_version, (tool, *arguments), kept
    )


def captured_text(completed: Completed) -> str:
    """Read (and delete) the standard output captured by ``run(..., capture=True)``."""
    if completed.stdout_path is None:
        return ""
    try:
        return completed.stdout_path.read_text(encoding="utf-8", errors="replace")
    finally:
        with contextlib.suppress(OSError):
            completed.stdout_path.unlink()


def result(
    operation: str,
    completed: Completed,
    *,
    outputs: Mapping[str, Path | str] | None = None,
    **extra: object,
) -> dict[str, Any]:
    """The JSON-serialisable result dictionary of a facade."""
    payload: dict[str, Any] = {
        "status": "completed",
        "backend": "augustus",
        "operation": operation,
        "program": completed.tool,
        "version": completed.version,
        "command": list(completed.command),
        "return_code": completed.return_code,
        "outputs": {name: str(path) for name, path in (outputs or {}).items() if path is not None},
        "log_tail": tail(completed.log),
        "provenance": {
            "tool": "AUGUSTUS",
            "version": completed.version,
            "tested_version": TESTED_VERSION,
            "executable": completed.executable,
            "citation": CITATION,
        },
    }
    payload.update(extra)
    return payload


def line_summary(path: Path, *, skip_comments: bool = True) -> dict[str, Any]:
    """Line count and the first lines of a text output (comment lines starting with ``#`` are not counted)."""
    count = 0
    preview: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if skip_comments and line.startswith("#"):
                continue
            if not line.strip():
                continue
            count += 1
            if len(preview) < PREVIEW_LINES:
                preview.append(line.rstrip("\n")[:300])
    return {"lines": count, "preview": preview}


def fasta_count(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.startswith(">"))


# ---------------------------------------------------------------------------
# Argument validation shared by the facades
# ---------------------------------------------------------------------------


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def name_token(value: object, name: str) -> str:
    """A species or file-name token: letters, digits, '_', '.', '+', '-'; never a flag."""
    text = cli.text(value, name)
    if NAME.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} may contain only letters, digits, '_', '.', '+' and '-': {text!r}"
        )
    return text


def optional_file(value: object, name: str) -> Path | None:
    return None if value is None else cli.existing_file(str(value), name)


def optional_int(
    value: object, name: str, *, minimum: int | None = None, maximum: int | None = None
) -> int | None:
    return None if value is None else cli.integer(value, name, minimum=minimum, maximum=maximum)


def optional_number(
    value: object, name: str, *, minimum: float | None = None, maximum: float | None = None
) -> float | None:
    return None if value is None else cli.number(value, name, minimum=minimum, maximum=maximum)


def valued(flag_name: str, value: object | None) -> list[str]:
    """``[flag, str(value)]`` for a fixed ``--flag value`` option, nothing when the value is None."""
    return [] if value is None else [flag_name, str(value)]


def eq(flag_name: str, value: object | None) -> list[str]:
    """``["--flag=value"]`` for a fixed ``--flag=value`` option, nothing when the value is None."""
    return [] if value is None else [f"{flag_name}={value}"]


def switch(flag_name: str, on: bool) -> list[str]:
    return [flag_name] if on else []


def boolean_word(
    flag_name: str, value: bool | None, words: tuple[str, str] = ("true", "false")
) -> list[str]:
    """``--flag=true``/``--flag=false`` (AUGUSTUS style) for an optional bool."""
    if value is None:
        return []
    return [f"{flag_name}={words[0] if flag(value, flag_name) else words[1]}"]


def simple(
    operation: str,
    tool: str,
    arguments: Sequence[str],
    *,
    directory: str | None,
    output: Path,
    stdin: Path | None = None,
    to_stdout: bool = True,
    kind: str = "text",
    accept: Sequence[int] = (0,),
    config_path: str | None = None,
    **extra: object,
) -> dict[str, Any]:
    """Run one filter-style program whose single output is ``output`` (standard output or a file it writes)."""
    if to_stdout:
        completed = run(
            tool,
            arguments,
            directory=directory,
            cwd=output.parent,
            stdin=stdin,
            stdout=output,
            accept=accept,
            config_path=config_path,
        )
    else:
        completed = run(
            tool,
            arguments,
            directory=directory,
            cwd=output.parent,
            stdin=stdin,
            capture=True,
            expected=[output],
            accept=accept,
            config_path=config_path,
        )
        captured_text(completed)
    summary: dict[str, Any] = (
        {"records": fasta_count(output)} if kind == "fasta" else line_summary(output)
    )
    summary.update(extra)
    return result(operation, completed, outputs={"output": output}, **summary)
