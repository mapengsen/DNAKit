"""Opt-in local Python execution for the Agent's data-processing fallback.

This is a subprocess runner, not a sandbox. Code has the current user's access
to files, network and installed packages. No credentials or environment snapshot
are written by the runner. Use only for tasks authorized by the caller.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class PythonExecutionConfig:
    """Caller-owned execution limits; the model cannot change these settings.

    Each run creates a unique audit directory under ``workspace/.dnakit_agent``.
    Scripts use the current Python interpreter and workspace as their cwd.
    Timeout is per attempt, including validation code written by the model.
    Output tails sent to the LLM are bounded; full stdout/stderr stay on disk.
    ``workspace`` sets the working directory, not a filesystem access boundary.
    """

    workspace: str | Path
    timeout_seconds: int = 120
    max_attempts: int = 5
    max_output_bytes: int = 16_000

    def __post_init__(self) -> None:
        if not isinstance(self.workspace, (str, Path)) or not str(self.workspace).strip():
            raise ConfigurationError("workspace must be a path.", code="AGENT_EXECUTION_CONFIG")
        for name, upper in (
            ("timeout_seconds", 86_400),
            ("max_attempts", 100),
            ("max_output_bytes", 64_000),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= upper:
                raise ConfigurationError(
                    f"{name} must be an integer in [1, {upper}].", code="AGENT_EXECUTION_CONFIG"
                )


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _stop(process: subprocess.Popen[bytes]) -> None:
    if os.name == "posix":
        with suppress(ProcessLookupError):
            os.killpg(process.pid, signal.SIGKILL)
    elif process.poll() is None:
        process.kill()
    process.wait()


def _tail(path: Path, limit: int) -> tuple[str, bool]:
    if not path.exists():
        return "", False
    with path.open("rb") as stream:
        size = stream.seek(0, 2)
        stream.seek(max(0, size - limit))
        return stream.read(limit).decode("utf-8", errors="replace"), size > limit


class _PythonExecutor:
    def __init__(self, config: PythonExecutionConfig) -> None:
        self.config = config
        self.workspace = Path(config.workspace).expanduser().resolve()
        root = self.workspace / ".dnakit_agent"
        root.mkdir(parents=True, exist_ok=True)
        self.directory = Path(tempfile.mkdtemp(prefix="run_", dir=root))
        self.attempts = 0
        _write_json(
            self.directory / "settings.json",
            {
                "workspace": str(self.workspace),
                "python": sys.executable,
                "timeout_seconds": config.timeout_seconds,
                "max_attempts": config.max_attempts,
                "max_output_bytes": config.max_output_bytes,
            },
        )

    def record(self, entry: dict[str, object]) -> None:
        with (self.directory / "history.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(entry, ensure_ascii=False) + "\n")

    def finish(self, status: str, text: str) -> None:
        _write_json(
            self.directory / "summary.json",
            {
                "status": status,
                "text": text,
                "code_attempts": self.attempts,
            },
        )

    def execute(
        self, *, code: str, reason: str, inputs: object, expected_outputs: list[str], step: int
    ) -> dict[str, object]:
        if self.attempts >= self.config.max_attempts:
            raise ConfigurationError(
                "Python execution attempt limit reached; inspect saved errors or use other tools.",
                code="AGENT_CODE_ATTEMPT_LIMIT",
            )
        if len(code.encode("utf-8")) > 200_000:
            raise ConfigurationError("Python code exceeds 200 KB.", code="AGENT_CODE_SIZE_LIMIT")
        if len(expected_outputs) > 100:
            raise ConfigurationError(
                "At most 100 output paths are allowed.", code="AGENT_CODE_OUTPUT_LIMIT"
            )
        # Resolve before creating the attempt so invalid paths have no execution effects.
        targets = []
        for name in expected_outputs:
            path = Path(name).expanduser()
            targets.append((path if path.is_absolute() else self.workspace / path).resolve())
        self.attempts += 1
        directory = self.directory / f"attempt_{self.attempts:03d}"
        directory.mkdir()
        script = directory / "script.py"
        input_path = directory / "inputs.json"
        script.write_text(code, encoding="utf-8")
        _write_json(input_path, inputs)
        _write_json(
            directory / "request.json",
            {
                "step": step,
                "reason": reason,
                "expected_outputs": [str(p) for p in targets],
            },
        )
        stdout_path, stderr_path = directory / "stdout.log", directory / "stderr.log"
        start = time.monotonic()
        status = "succeeded"
        returncode = None
        execution_error = None
        try:
            with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                process = subprocess.Popen(
                    [sys.executable, "-u", str(script), str(input_path)],
                    cwd=self.workspace,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=os.name == "posix",
                )
                try:
                    returncode = process.wait(timeout=self.config.timeout_seconds)
                    if returncode != 0:
                        status = "failed"
                except subprocess.TimeoutExpired:
                    status = "timeout"
                except BaseException:
                    status = "interrupted"
                    raise
                finally:
                    # Also reap descendants left behind by a completed script on POSIX.
                    _stop(process)
                    returncode = process.returncode
        except OSError as exc:
            status, execution_error = "failed", str(exc)
        finally:
            outputs = []
            for target in targets:
                try:
                    size = target.stat().st_size if target.is_file() else 0
                    outputs.append(
                        {"path": str(target), "nonempty_file": size > 0, "size_bytes": size}
                    )
                except OSError as exc:
                    outputs.append({"path": str(target), "nonempty_file": False, "error": str(exc)})
            if status == "succeeded" and any(not item["nonempty_file"] for item in outputs):
                status = "missing_outputs"
            stdout, stdout_truncated = _tail(stdout_path, self.config.max_output_bytes)
            stderr, stderr_truncated = _tail(stderr_path, self.config.max_output_bytes)
            result = {
                "status": status,
                "returncode": returncode,
                "elapsed_seconds": round(time.monotonic() - start, 3),
                "attempt": self.attempts,
                "directory": str(directory),
                "script_path": str(script),
                "inputs_path": str(input_path),
                "stdout_path": str(stdout_path),
                "stderr_path": str(stderr_path),
                "stdout": stdout,
                "stderr": stderr,
                "stdout_truncated": stdout_truncated,
                "stderr_truncated": stderr_truncated,
                "outputs": outputs,
                "execution_error": execution_error,
            }
            _write_json(directory / "result.json", result)
        return result


__all__ = ["PythonExecutionConfig"]
