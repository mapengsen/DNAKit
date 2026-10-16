"""Run Agent-callable DNAKit functions in a separate process for the graphical interface.

The interface talks to ``python -m dnakit.gui._worker`` with one JSON request per line on stdin
and one JSON reply per line on stdout. A separate process lets the interface stop a call that
runs too long, or that the user cancels, without restarting the interface. Anything a function
prints goes to stderr, so stdout carries only the protocol.
"""

from __future__ import annotations

import json
import os
import queue
import subprocess
import sys
import threading
from pathlib import Path
from typing import IO, Any, cast

READY_TIMEOUT_S = 180


class ToolWorker:
    """One long-lived worker process that holds the DNAKit tool registry."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._process: subprocess.Popen[str] | None = None
        self._replies: queue.Queue[str] = queue.Queue()

    def _start(self) -> None:
        env = {**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"}
        process = subprocess.Popen(
            [sys.executable, "-m", "dnakit.gui._worker"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            bufsize=1,
            env=env,
        )
        self._process = process
        self._replies = queue.Queue()
        threading.Thread(target=_pump, args=(process, self._replies), daemon=True).start()
        if self._wait(READY_TIMEOUT_S) is None:
            self.stop()
            raise RuntimeError("The DNAKit worker did not start.")

    def _wait(self, timeout: float | None) -> dict[str, Any] | None:
        try:
            line = self._replies.get(timeout=timeout)
        except queue.Empty:
            return None
        if not line:
            raise RuntimeError("The DNAKit worker stopped; run the call again.")
        return cast(dict[str, Any], json.loads(line))

    def stop(self) -> None:
        """Kill the worker; a running call returns an error and the next call restarts it."""

        process = self._process
        if process is not None and process.poll() is None:
            process.kill()
            process.wait(10)

    def call(
        self,
        name: str,
        arguments: dict[str, Any],
        cwd: Path,
        *,
        timeout: float | None,
        allow_side_effects: bool = False,
    ) -> tuple[bool, Any]:
        """Run one registry tool; return ``(True, result)`` or ``(False, error message)``."""

        with self._lock:
            if self._process is None or self._process.poll() is not None:
                self._start()
            process = self._process
            assert process is not None and process.stdin is not None
            request = {
                "name": name,
                "arguments": arguments,
                "cwd": str(cwd),
                "allow_side_effects": allow_side_effects,
            }
            process.stdin.write(json.dumps(request, ensure_ascii=False) + "\n")
            process.stdin.flush()
            try:
                reply = self._wait(timeout)
            except RuntimeError as exc:
                self.stop()
                return False, str(exc)
            if reply is None:
                self.stop()
                return False, f"Stopped after {timeout:.0f} s (time limit of this interface)."
        if reply["status"] == "ok":
            return True, reply["result"]
        return False, reply["error"]


def _pump(process: subprocess.Popen[str], replies: queue.Queue[str]) -> None:
    stdout = cast(IO[str], process.stdout)
    for line in stdout:
        replies.put(line)
    replies.put("")


def main() -> None:
    """Serve requests from the graphical interface until stdin closes."""

    replies = os.fdopen(os.dup(1), "w", buffering=1, encoding="utf-8")
    os.dup2(2, 1)
    sys.stdout = sys.stderr

    from dnakit.tools.registry import default_tool_registry

    registry = default_tool_registry(verbose=False)
    replies.write(json.dumps({"status": "ready"}) + "\n")
    for line in sys.stdin:
        request = json.loads(line)
        try:
            os.chdir(request["cwd"])
            result = registry.execute(
                request["name"],
                request["arguments"],
                allow_side_effects=bool(request.get("allow_side_effects", False)),
            )
            reply = json.dumps({"status": "ok", "result": result}, ensure_ascii=False, default=str)
        except Exception as exc:  # every DNAKit or Python error is reported to the interface
            reply = json.dumps({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
        replies.write(reply + "\n")


if __name__ == "__main__":
    main()
