"""Passive registrations and bounded execution for user-supplied CLI tools."""

from __future__ import annotations

import os
import re
import shutil
import signal
import stat
import subprocess
import threading
import time
from collections.abc import Iterable
from contextlib import suppress
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import BinaryIO, cast

from dnakit.core import BackendInfo
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)

_VERSION = re.compile(r"(?P<version>\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?)")
_MAX_EXTERNAL_ARGUMENTS = 256
_MAX_EXTERNAL_ARGUMENT_BYTES = 65_536
_MAX_MONITORED_OUTPUTS = 16


def _terminate_process_tree(process: subprocess.Popen[bytes]) -> None:
    """Best-effort termination of a version command and its descendants."""

    if process.poll() is not None:
        return
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            return
        except OSError:
            process.kill()
    else:  # pragma: no cover - Linux is the required validation platform.
        process.kill()


@dataclass(frozen=True, slots=True)
class BoundedCommandResult:
    """Captured metadata from one explicitly requested, bounded CLI execution."""

    return_code: int
    output: str
    elapsed_seconds: float


def _validated_command_arguments(arguments: Iterable[str]) -> tuple[str, ...]:
    if isinstance(arguments, (str, bytes)):
        raise ConfigurationError("External command arguments must be an iterable of strings.")
    try:
        resolved = tuple(islice(iter(arguments), _MAX_EXTERNAL_ARGUMENTS + 1))
    except TypeError as exc:
        raise ConfigurationError(
            "External command arguments must be an iterable of strings."
        ) from exc
    if len(resolved) > _MAX_EXTERNAL_ARGUMENTS:
        raise ConfigurationError(
            "External command exceeds the argument-count limit.",
            code="EXTERNAL_ARGUMENT_LIMIT_EXCEEDED",
            context={"max_arguments": _MAX_EXTERNAL_ARGUMENTS},
        )
    if any(
        not isinstance(argument, str) or not argument or "\x00" in argument for argument in resolved
    ):
        raise ConfigurationError("External command arguments must be non-empty NUL-free strings.")
    argument_bytes = sum(len(argument.encode("utf-8")) + 1 for argument in resolved)
    if argument_bytes > _MAX_EXTERNAL_ARGUMENT_BYTES:
        raise ConfigurationError(
            "External command exceeds the encoded argument-size limit.",
            code="EXTERNAL_ARGUMENT_LIMIT_EXCEEDED",
            context={
                "max_argument_bytes": _MAX_EXTERNAL_ARGUMENT_BYTES,
                "observed_argument_bytes": argument_bytes,
            },
        )
    return resolved


def execute_bounded_command(
    executable_path: str | Path,
    arguments: Iterable[str],
    *,
    backend_id: str,
    cwd: str | Path,
    timeout_seconds: float = 300.0,
    max_output_bytes: int = 1_000_000,
    monitored_output_paths: Iterable[str | Path] = (),
    max_monitored_output_bytes: int = 100_000_000,
) -> BoundedCommandResult:
    """Run an explicit executable without a shell under time and output limits.

    Domain adapters remain responsible for constructing a fixed command
    whitelist. This executor only accepts an absolute resolved executable,
    bounds the argument vector and captured output, isolates a POSIX process
    group, and monitors declared output artifacts inside ``cwd``.
    """

    if not isinstance(backend_id, str) or not backend_id.strip():
        raise ConfigurationError("backend_id must be a non-empty string.")
    try:
        executable = Path(executable_path).expanduser().resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise BackendUnavailableError(
            "Explicit external backend executable is unavailable.",
            context={"backend_id": backend_id, "executable": str(executable_path)},
        ) from exc
    try:
        executable_stat = executable.stat()
    except OSError as exc:  # pragma: no cover - resolve already verified the path.
        raise BackendUnavailableError(
            "Explicit external backend executable is unavailable.",
            context={"backend_id": backend_id, "executable": str(executable)},
        ) from exc
    if not stat.S_ISREG(executable_stat.st_mode) or not os.access(executable, os.X_OK):
        raise BackendUnavailableError(
            "Explicit external backend path must be an executable regular file.",
            context={"backend_id": backend_id, "executable": str(executable)},
        )

    try:
        working_directory = Path(cwd).expanduser().resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError("External command cwd must be an existing directory.") from exc
    if not working_directory.is_dir():
        raise ConfigurationError("External command cwd must be an existing directory.")
    resolved_arguments = _validated_command_arguments(arguments)
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or not 0 < timeout_seconds <= 86_400
    ):
        raise ConfigurationError("timeout_seconds must be in (0, 86400].")
    for name, value, maximum in (
        ("max_output_bytes", max_output_bytes, 100_000_000),
        ("max_monitored_output_bytes", max_monitored_output_bytes, 10_000_000_000),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
            raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")

    if isinstance(monitored_output_paths, (str, bytes, Path)):
        monitored_values: Iterable[str | Path] = (monitored_output_paths,)
    else:
        monitored_values = monitored_output_paths
    try:
        monitored = tuple(
            islice(
                (Path(path).expanduser() for path in monitored_values),
                _MAX_MONITORED_OUTPUTS + 1,
            )
        )
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("monitored_output_paths must contain filesystem paths.") from exc
    if len(monitored) > _MAX_MONITORED_OUTPUTS:
        raise ConfigurationError(
            "Too many monitored external output paths.",
            code="EXTERNAL_OUTPUT_PATH_LIMIT_EXCEEDED",
            context={"max_paths": _MAX_MONITORED_OUTPUTS},
        )
    resolved_monitored: tuple[Path, ...] = tuple(
        path.resolve(strict=False) if path.is_absolute() else (working_directory / path).resolve()
        for path in monitored
    )
    if any(
        path == working_directory or not path.is_relative_to(working_directory)
        for path in resolved_monitored
    ):
        raise ConfigurationError(
            "Monitored external outputs must stay inside cwd.",
            code="UNSAFE_EXTERNAL_OUTPUT_PATH",
        )

    environment = os.environ.copy()
    environment.update({"LANG": "C", "LC_ALL": "C"})
    command = (str(executable), *resolved_arguments)
    started = time.monotonic()
    try:
        process = subprocess.Popen(
            command,
            cwd=str(working_directory),
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=False,
            close_fds=True,
            start_new_session=os.name == "posix",
        )
    except OSError as exc:
        raise BackendExecutionError(
            "External backend command could not be started.",
            code="BACKEND_COMMAND_EXECUTION_FAILED",
            context={"backend_id": backend_id, "error_type": type(exc).__name__},
        ) from exc
    if process.stdout is None:  # pragma: no cover - guaranteed by stdout=PIPE.
        process.kill()
        process.wait()
        raise BackendExecutionError("External backend output pipe could not be created.")

    stdout = cast(BinaryIO, process.stdout)
    output = bytearray()
    output_too_large = threading.Event()

    def read_bounded_output() -> None:
        try:
            while len(output) <= max_output_bytes:
                remaining = max_output_bytes + 1 - len(output)
                chunk = stdout.read(min(65_536, remaining))
                if not chunk:
                    return
                output.extend(chunk)
                if len(output) > max_output_bytes:
                    output_too_large.set()
                    return
        except (OSError, ValueError):
            return

    reader = threading.Thread(target=read_bounded_output, daemon=True)
    reader.start()
    timed_out = False
    artifact_too_large: Path | None = None
    monitored_output_bytes = 0
    deadline = started + float(timeout_seconds)
    try:
        while process.poll() is None:
            if output_too_large.is_set():
                _terminate_process_tree(process)
                break
            monitored_output_bytes = 0
            for path in resolved_monitored:
                try:
                    monitored_output_bytes += path.stat().st_size
                    if monitored_output_bytes > max_monitored_output_bytes:
                        artifact_too_large = path
                        _terminate_process_tree(process)
                        break
                except FileNotFoundError:
                    continue
                except OSError as exc:
                    _terminate_process_tree(process)
                    raise BackendExecutionError(
                        "External backend output artifact could not be inspected.",
                        code="BACKEND_ARTIFACT_INSPECTION_FAILED",
                        context={"backend_id": backend_id, "path": str(path)},
                    ) from exc
            if artifact_too_large is not None:
                break
            if time.monotonic() >= deadline:
                timed_out = True
                _terminate_process_tree(process)
                break
            time.sleep(0.005)
        process.wait()
    finally:
        if process.poll() is None:
            _terminate_process_tree(process)
            process.wait()
        reader.join(timeout=1.0)
        if reader.is_alive() and os.name == "posix":
            with suppress(OSError, ProcessLookupError):
                os.killpg(process.pid, signal.SIGKILL)
        stdout.close()
        reader.join(timeout=1.0)

    monitored_output_bytes = 0
    for path in resolved_monitored:
        try:
            monitored_output_bytes += path.stat().st_size
            if monitored_output_bytes > max_monitored_output_bytes:
                artifact_too_large = path
                break
        except FileNotFoundError:
            continue
        except OSError as exc:
            raise BackendExecutionError(
                "External backend output artifact could not be inspected.",
                code="BACKEND_ARTIFACT_INSPECTION_FAILED",
                context={"backend_id": backend_id, "path": str(path)},
            ) from exc
    if timed_out:
        raise BackendTimeoutError(
            "External backend command exceeded its timeout.",
            context={"backend_id": backend_id, "timeout_seconds": float(timeout_seconds)},
        )
    if output_too_large.is_set() or len(output) > max_output_bytes:
        raise BackendExecutionError(
            "External backend command output exceeded max_output_bytes.",
            code="BACKEND_COMMAND_OUTPUT_LIMIT",
            context={
                "backend_id": backend_id,
                "max_output_bytes": max_output_bytes,
                "observed_bytes": len(output),
            },
        )
    if artifact_too_large is not None:
        raise BackendExecutionError(
            "External backend artifact exceeded max_monitored_output_bytes.",
            code="BACKEND_ARTIFACT_OUTPUT_LIMIT",
            context={
                "backend_id": backend_id,
                "path": str(artifact_too_large),
                "max_output_bytes": max_monitored_output_bytes,
                "observed_bytes": monitored_output_bytes,
            },
        )
    return BoundedCommandResult(
        return_code=process.returncode,
        output=bytes(output).decode("utf-8", errors="replace"),
        elapsed_seconds=time.monotonic() - started,
    )


@dataclass(frozen=True, slots=True)
class ExternalCLIAdapter:
    """Explicit opt-in adapter handle that never executes during discovery or loading."""

    backend_id: str
    executable: str
    capabilities: frozenset[str]
    license_expression: str

    def __post_init__(self) -> None:
        if not all(
            isinstance(value, str) and value.strip() for value in (self.backend_id, self.executable)
        ):
            raise ConfigurationError("External backend identifiers must be non-empty strings.")
        if not self.capabilities:
            raise ConfigurationError("External backend capabilities must not be empty.")

    def probe(self) -> BackendInfo:
        return passive_external_probe(
            self.backend_id,
            self.executable,
            self.capabilities,
            self.license_expression,
        )

    def version(
        self,
        *,
        timeout_seconds: float = 5.0,
        max_output_bytes: int = 1_000_000,
    ) -> BackendInfo:
        """Explicitly execute only the backend's version command with a timeout."""

        info = self.probe()
        if not info.available or info.executable_path is None:
            raise BackendUnavailableError(
                "External backend executable is unavailable.",
                context={"backend_id": self.backend_id, "executable": self.executable},
            )
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or not 0 < timeout_seconds <= 60
        ):
            raise ConfigurationError("timeout_seconds must be in (0, 60].")
        if (
            isinstance(max_output_bytes, bool)
            or not isinstance(max_output_bytes, int)
            or not 1 <= max_output_bytes <= 10_000_000
        ):
            raise ConfigurationError("max_output_bytes must be an integer in [1, 10000000].")
        flag = "version" if self.backend_id == "mmseqs2" else "--version"
        deadline = time.monotonic() + float(timeout_seconds)
        try:
            process = subprocess.Popen(
                (info.executable_path, flag),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                start_new_session=os.name == "posix",
            )
        except OSError as exc:
            raise BackendExecutionError(
                "External backend version command could not be started.",
                code="BACKEND_VERSION_EXECUTION_FAILED",
                context={"backend_id": self.backend_id, "error_type": type(exc).__name__},
            ) from exc
        if process.stdout is None:  # pragma: no cover - guaranteed by stdout=PIPE.
            process.kill()
            process.wait()
            raise BackendExecutionError("External backend output pipe could not be created.")
        stdout = cast(BinaryIO, process.stdout)
        output = bytearray()
        output_too_large = threading.Event()

        def read_bounded_output() -> None:
            while len(output) <= max_output_bytes:
                remaining = max_output_bytes + 1 - len(output)
                chunk = stdout.read(min(65_536, remaining))
                if not chunk:
                    return
                output.extend(chunk)
                if len(output) > max_output_bytes:
                    output_too_large.set()
                    return

        reader = threading.Thread(target=read_bounded_output, daemon=True)
        reader.start()
        timed_out = False
        try:
            while process.poll() is None:
                if output_too_large.is_set():
                    _terminate_process_tree(process)
                    break
                if time.monotonic() >= deadline:
                    timed_out = True
                    _terminate_process_tree(process)
                    break
                time.sleep(0.005)
            process.wait()
        finally:
            if process.poll() is None:
                _terminate_process_tree(process)
                process.wait()
            reader.join(timeout=1.0)
            stdout.close()
        if timed_out:
            raise BackendTimeoutError(
                "External backend version command exceeded its timeout.",
                context={
                    "backend_id": self.backend_id,
                    "timeout_seconds": float(timeout_seconds),
                },
            )
        if output_too_large.is_set() or len(output) > max_output_bytes:
            raise BackendExecutionError(
                "External backend version output exceeded max_output_bytes.",
                code="BACKEND_VERSION_OUTPUT_LIMIT",
                context={
                    "backend_id": self.backend_id,
                    "max_output_bytes": max_output_bytes,
                    "observed_bytes": len(output),
                },
            )
        text = bytes(output).decode("utf-8", errors="replace").strip()
        match = _VERSION.search(text)
        return BackendInfo(
            info.name,
            version=None if match is None else match.group("version"),
            executable_path=info.executable_path,
            license_expression=info.license_expression,
            capabilities=info.capabilities,
            available=process.returncode == 0,
            metadata={
                "probe_mode": "explicit-version-command",
                "return_code": process.returncode,
                "output_excerpt": text[:1_000],
                "max_output_bytes": max_output_bytes,
            },
        )


def passive_external_probe(
    backend_id: str,
    executable: str,
    capabilities: frozenset[str],
    license_expression: str,
) -> BackendInfo:
    """Locate an executable without running it or reading sequence data."""

    path = shutil.which(executable)
    return BackendInfo(
        backend_id,
        executable_path=path,
        license_expression=license_expression,
        capabilities=capabilities,
        available=path is not None,
        metadata={
            "probe_mode": "path-only",
            "execution": "explicit opt-in only",
            "redistributed": False,
        },
    )


EXTERNAL_CLI_ADAPTERS = (
    ExternalCLIAdapter("blastn", "blastn", frozenset({"sequence-search"}), "LicenseRef-NCBI"),
    ExternalCLIAdapter(
        "mmseqs2",
        "mmseqs",
        frozenset({"sequence-search", "clustering"}),
        "GPL-3.0-or-later",
    ),
    ExternalCLIAdapter(
        "sourmash",
        "sourmash",
        frozenset({"sketch", "sequence-search", "similarity"}),
        "BSD-3-Clause",
    ),
    ExternalCLIAdapter(
        "dashing",
        "dashing",
        frozenset({"sketch", "similarity"}),
        "GPL-3.0-only",
    ),
)


__all__ = [
    "EXTERNAL_CLI_ADAPTERS",
    "BoundedCommandResult",
    "ExternalCLIAdapter",
    "execute_bounded_command",
    "passive_external_probe",
]
