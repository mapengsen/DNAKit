"""Shared lazy dependencies, validation, and bounded external execution."""

from __future__ import annotations

import contextlib
import importlib
import math
import os
import shutil
import signal
import subprocess
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.core import DNA, DNARecord, DNASequence
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)


def require(module: str, extra: str) -> Any:
    """Import a dependency without changing the user's installation."""
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(
            f"{module} could not be imported: {exc}",
            hint=f'Install the optional dependency with: python -m pip install "dnakit[{extra}]"',
        ) from exc


def symbol(module: ModuleType | Any, name: str) -> Any:
    candidate = getattr(module, name, None)
    if candidate is None:
        raise BackendUnavailableError(f"Installed {module.__name__} does not provide {name}.")
    return candidate


def positive(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def sequence_text(value: str | DNASequence | DNARecord | DNA, *, gaps: bool = False) -> str:
    if isinstance(value, DNA):
        value = value.record
    if isinstance(value, DNARecord):
        value = value.sequence
    if isinstance(value, DNASequence):
        if any(not isinstance(part, str) for part in value.parts):
            raise ConfigurationError("Resolve explicit Gap objects before using this API.")
        value = "".join(part for part in value.parts if isinstance(part, str))
    if not isinstance(value, str) or not value:
        raise ConfigurationError("A non-empty DNA sequence is required.")
    value = value.upper()
    alphabet = set("ACGTRYSWKMBDHVN" + ("-" if gaps else ""))
    if set(value) - alphabet:
        raise ConfigurationError("Sequence contains unsupported DNA symbols.")
    return value


def named_sequence_items(sequences: DNA | Mapping[str, Any]) -> tuple[tuple[str, Any], ...]:
    """Preserve IDs and order for a DNA batch or an existing named mapping."""

    if isinstance(sequences, DNA):
        items = tuple((record.id, record.sequence) for record in sequences)
        if len({record_id for record_id, _ in items}) != len(items):
            raise ConfigurationError("DNA record IDs must be unique for named sequence APIs.")
        return items
    if not isinstance(sequences, Mapping):
        raise ConfigurationError("sequences must be DNA([...]) or an ID-to-sequence mapping.")
    return tuple(sequences.items())


def matrix(values: Any, *, min_rows: int = 1) -> Any:
    np = require("numpy", "bio-analysis")
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("Expected a rectangular numeric matrix.") from exc
    if array.ndim != 2 or array.shape[0] < min_rows or array.shape[1] < 1:
        raise ConfigurationError(f"Expected a matrix with at least {min_rows} rows and one column.")
    if not np.isfinite(array).all():
        raise ConfigurationError("Matrix must contain only finite values.")
    return array


def input_path(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}")
    return path


@dataclass(frozen=True)
class ExecutionConfig:
    """Limits per external process; sequence output is streamed directly to disk."""

    threads: int = 1
    timeout_seconds: float = 3600.0
    max_output_bytes: int = 1_000_000
    max_artifact_bytes: int = 1_000_000_000

    def __post_init__(self) -> None:
        positive(self.threads, "threads")
        positive(self.max_output_bytes, "max_output_bytes")
        positive(self.max_artifact_bytes, "max_artifact_bytes")
        if self.max_output_bytes > 100_000_000 or self.max_artifact_bytes > 10_000_000_000:
            raise ConfigurationError("Execution limits exceed supported log/artifact maxima.")
        if not math.isfinite(self.timeout_seconds) or not 0 < self.timeout_seconds <= 86_400:
            raise ConfigurationError("timeout_seconds must be in (0, 86400].")


@dataclass(frozen=True)
class ToolResult:
    """Persistent output directory, produced files, and diagnostic command log."""

    backend: str
    directory: Path
    artifacts: Mapping[str, Path]
    log: str


def workspace(value: str | Path) -> Path:
    """Claim a fresh output directory; never overwrite previous results."""
    path = Path(value).expanduser().resolve()
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise ConfigurationError(f"Output directory already exists: {path}") from exc
    return path


def executable(name: str) -> str:
    resolved = shutil.which(name)
    if resolved is None:
        raise BackendUnavailableError(
            f"External executable {name!r} was not found.",
            hint="Install this executable separately and place it on PATH, or pass its full path.",
        )
    return str(Path(resolved).resolve())


def run(
    program: str,
    arguments: Iterable[str],
    directory: Path,
    config: ExecutionConfig,
    *,
    outputs: Iterable[Path] = (),
) -> str:
    outputs = tuple(outputs)
    result = execute_bounded_command(
        executable(program),
        tuple(arguments),
        backend_id=Path(program).name,
        cwd=directory,
        timeout_seconds=config.timeout_seconds,
        max_output_bytes=config.max_output_bytes,
        monitored_output_paths=tuple(outputs),
        max_monitored_output_bytes=config.max_artifact_bytes,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"{program} failed with exit code {result.return_code}.",
            context={"directory": str(directory), "log": result.output},
        )
    if missing := [str(path) for path in outputs if not path.is_file()]:
        raise BackendExecutionError(f"{program} did not produce expected files: {missing}")
    return result.output


def artifacts(directory: Path) -> dict[str, Path]:
    return {
        str(path.relative_to(directory)): path for path in directory.rglob("*") if path.is_file()
    }


def run_to_file(
    program: str, arguments: Iterable[str], output: Path, config: ExecutionConfig
) -> str:
    """Stream stdout to a fresh artifact, keep stderr separate, enforce time/size limits."""
    command = (executable(program), *arguments)
    diagnostic = output.with_name(output.name + ".stderr")
    with output.open("xb") as stdout, diagnostic.open("xb") as stderr:
        process = subprocess.Popen(
            command,
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            cwd=output.parent,
            shell=False,
            start_new_session=os.name == "posix",
        )
        deadline = time.monotonic() + config.timeout_seconds
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{program} exceeded its time limit.")
                if (
                    output.stat().st_size > config.max_artifact_bytes
                    or diagnostic.stat().st_size > config.max_output_bytes
                ):
                    raise BackendExecutionError(f"{program} exceeded its output size limit.")
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
    if (
        output.stat().st_size > config.max_artifact_bytes
        or diagnostic.stat().st_size > config.max_output_bytes
    ):
        raise BackendExecutionError(f"{program} exceeded its output size limit.")
    log = diagnostic.read_text(encoding="utf-8", errors="replace")
    if process.returncode != 0:
        raise BackendExecutionError(
            f"{program} exited with code {process.returncode}.", context={"log": log}
        )
    return log
