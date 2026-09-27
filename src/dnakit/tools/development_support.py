"""Durable development transactions shared by generated tools and text overlays.

Candidates live outside the installed package. Activation writes reviewed files
with precondition hashes, backups and a recoverable journal; it never publishes.
"""

from __future__ import annotations

import difflib
import hashlib
import importlib
import json
import os
import re
import sys
import tempfile
import threading
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit._progress import track
from dnakit.exceptions import ConfigurationError

_LOCK_STATE = threading.local()


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(data: object) -> bytes:
    return (
        json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _file_hash(path: Path) -> str | None:
    return _digest(path.read_bytes()) if path.is_file() else None


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _source_root() -> Path:
    """Require the running DNAKit to originate from a writable source tree."""
    root = Path(__file__).resolve().parents[3]
    project = root / "pyproject.toml"
    if (
        "site-packages" in root.parts
        or not Path(__file__).resolve().is_relative_to(root / "src/dnakit/tools")
        or not project.is_file()
        or not (root / "src/dnakit/tools/agent.py").is_file()
    ):
        raise ConfigurationError(
            "Tool development requires a writable DNAKit source checkout.",
            code="AGENT_DEVELOPMENT_SOURCE_REQUIRED",
            hint="Use an editable installation of DNAKit; site-packages is not rewritten.",
        )
    parser = importlib.import_module("tomllib" if sys.version_info >= (3, 11) else "tomli")
    if parser.loads(project.read_text(encoding="utf-8")).get("project", {}).get("name") != "dnakit":
        raise ConfigurationError(
            "The source project is not DNAKit.", code="AGENT_DEVELOPMENT_SOURCE_REQUIRED"
        )
    try:
        with tempfile.TemporaryFile(dir=root / "src/dnakit"):
            pass
    except OSError as exc:
        raise ConfigurationError(
            "The DNAKit source tree is not writable.", code="AGENT_DEVELOPMENT_SOURCE_REQUIRED"
        ) from exc
    return root


@contextmanager
def _write_lock(root: Path) -> Iterator[None]:
    """OS-owned lock is released on interruption, including process termination."""
    held: set[Path] = getattr(_LOCK_STATE, "held", set())
    if root in held:
        yield
        return
    directory = root / ".dnakit/development"
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / ".lock").open("a+b") as stream:
        if os.name == "nt":
            msvcrt = importlib.import_module("msvcrt")

            stream.seek(0)
            stream.write(b"0")
            stream.flush()
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            _LOCK_STATE.held = held | {root}
            yield
        finally:
            _LOCK_STATE.held = held
            if os.name == "nt":
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


class DevelopmentJob:
    """One candidate, at most three revisions, with immutable validation cases."""

    def __init__(
        self, root: Path, kind: str, identity: str, action: dict[str, Any], verbose: bool
    ) -> None:
        self.root = root.resolve()
        self.verbose = verbose
        identity_hash = _digest(_json_bytes([kind, identity, action.get("reason")]))[:24]
        self.change_id = f"{kind}-{identity_hash}"
        if not re.fullmatch(r"[a-z_]+-[a-f0-9]{24}", self.change_id):
            raise ConfigurationError(
                "Invalid development identity.", code="AGENT_DEVELOPMENT_INVALID"
            )
        self.directory = self.root / ".dnakit/development" / self.change_id
        self.state: dict[str, Any] = {}
        self.replay: dict[str, Any] | None = None
        self._watched: dict[Path, tuple[Path, str | None]] = {}
        self._action_hash = _digest(_json_bytes(action))
        with _write_lock(self.root):
            state_path = self.directory / "state.json"
            self._state_hash = _file_hash(state_path)
            self.state = json.loads(state_path.read_text()) if state_path.is_file() else {}
            if self.state.get("phase") == "activating":
                self._recover_activation()
            initial = self.state.get("initial_action", action)
            for field in ("tests", "scenarios", "contract", "target", "name", "base_sha256"):
                if initial.get(field) != action.get(field):
                    raise ConfigurationError(
                        f"Candidate {field} cannot change between revisions.",
                        code="AGENT_DEVELOPMENT_FIXED_SCENARIOS",
                    )
            if self.state.get("phase") == "activated":
                if self.state.get("action_hash") != self._action_hash:
                    raise ConfigurationError(
                        "This candidate is already activated; use a new version reason.",
                        code="AGENT_DEVELOPMENT_ALREADY_ACTIVE",
                    )
                if not self._matches("new_sha256"):
                    raise ConfigurationError(
                        "Activated files changed; candidate cannot be replayed.",
                        code="AGENT_DEVELOPMENT_CHANGED",
                    )
                self.replay = self.state["result"]
            elif self.state.get("phase") in {"validating", "executing", "baseline", "candidate"}:
                # An actual call may have partial external effects. Never silently rerun it.
                self.replay = self.reject(
                    "Previous validation was interrupted; inspect outputs before a new candidate.",
                    code="AGENT_DEVELOPMENT_VALIDATION_UNCERTAIN",
                )
            elif self.state.get("phase") in {"rolled_back", "rolling_back"}:
                raise ConfigurationError(
                    "This candidate was rolled back; use a new version reason.",
                    code="AGENT_DEVELOPMENT_ROLLED_BACK",
                )
            elif self.state.get("action_hash") == self._action_hash and self.state.get("result"):
                self.replay = self.state["result"]
            self.attempt = int(self.state.get("attempt", 0))
            self.resuming = self.state.get("action_hash") == self._action_hash and self.state.get(
                "phase"
            ) in {"tested", "validated", "baseline_completed"}
            if self.replay is None:
                if self.attempt >= 3 and not self.resuming:
                    raise ConfigurationError(
                        "This candidate reached three revisions.",
                        code="AGENT_DEVELOPMENT_ROUND_LIMIT",
                    )
                if not self.resuming:
                    self.attempt += 1
                self.state.update(
                    version=1,
                    change_id=self.change_id,
                    identity=identity,
                    kind=kind,
                    initial_action=initial,
                    action_hash=self._action_hash,
                    attempt=self.attempt,
                    revision=int(self.state.get("revision", 0)) + 1,
                    phase=self.state.get("phase") if self.resuming else "staged",
                )
                self.state.pop("result", None)
                self.directory.mkdir(parents=True, exist_ok=True)
                _atomic_write(
                    self.directory / f"candidate-{self.attempt}.json", _json_bytes(action)
                )
                self._save()

    def _save(self) -> None:
        with _write_lock(self.root):
            self._check_current()
            data = _json_bytes(self.state)
            _atomic_write(self.directory / "state.json", data)
            self._state_hash = _digest(data)

    def _check_current(self) -> None:
        if _file_hash(self.directory / "state.json") != self._state_hash:
            raise ConfigurationError(
                "Another process advanced this candidate; reload its state.",
                code="AGENT_DEVELOPMENT_CONFLICT",
            )

    def phase(self, name: str, **data: Any) -> None:
        self.state.update(phase=name)
        self.state.setdefault("events", []).append({"phase": name, "attempt": self.attempt, **data})
        self._save()
        if self.verbose:
            print(f"[DNAKit development] {self.change_id}: {name} ({self.attempt}/3)", flush=True)

    def watch(self, paths: list[Path]) -> None:
        for path in paths:
            resolved = self._target(path)
            self._watched.setdefault(Path(os.path.abspath(path)), (resolved, _file_hash(resolved)))

    def _target(self, path: Path) -> Path:
        resolved = path.resolve()
        if not resolved.is_relative_to(self.root) or resolved == self.root:
            raise ConfigurationError(
                "Activation path must stay inside the source tree.", code="AGENT_DEVELOPMENT_PATH"
            )
        return resolved

    def _matches(self, field: str) -> bool:
        return bool(self.state.get("writes")) and all(
            _file_hash(self._target(Path(item["path"]))) == item[field]
            for item in self.state["writes"]
        )

    def _receipt(self) -> dict[str, Any]:
        return {
            "status": "activated",
            "change_id": self.change_id,
            "audit_dir": str(self.directory),
            "artifacts": [item["path"] for item in self.state["writes"]],
        }

    def _recover_activation(self) -> None:
        if self._matches("new_sha256"):
            self.state.update(phase="activated", result=self._receipt())
            self._save()
            return
        # Undo an incomplete multi-file transaction only while every file is known.
        for item in self.state["writes"]:
            if _file_hash(self._target(Path(item["path"]))) not in (
                item["new_sha256"],
                item["old_sha256"],
            ):
                raise ConfigurationError(
                    "Interrupted activation conflicts with edits; inspect the audit directory.",
                    code="AGENT_DEVELOPMENT_CONFLICT",
                )
        self._restore()
        self.state.update(
            phase="rejected",
            result={
                "status": "rejected",
                "reason": "Interrupted activation was rolled back.",
                "change_id": self.change_id,
                "audit_dir": str(self.directory),
            },
        )
        self._save()

    def _restore(self) -> None:
        for item in self.state["writes"]:
            path = Path(item["path"])
            if self._target(path) != path or _file_hash(path) not in {
                item["old_sha256"],
                item["new_sha256"],
            }:
                raise ConfigurationError(
                    "Rollback conflicts with a later edit; transaction journal was retained.",
                    code="AGENT_DEVELOPMENT_CONFLICT",
                )
        for item in reversed(self.state["writes"]):
            path = self._target(Path(item["path"]))
            if item["old_sha256"] is None:
                path.unlink(missing_ok=True)
            else:
                backup = self.directory / item["backup"]
                data = backup.read_bytes()
                if _digest(data) != item["old_sha256"]:
                    raise ConfigurationError(
                        "Rollback backup checksum mismatch.", code="AGENT_DEVELOPMENT_CHANGED"
                    )
                _atomic_write(path, data)

    def activate(self, writes: Mapping[Path, bytes]) -> dict[str, Any]:
        if not writes:
            raise ConfigurationError("No activation files.", code="AGENT_DEVELOPMENT_INVALID")
        with _write_lock(self.root):
            self._check_current()
            for original, watched in self._watched.items():
                resolved = self._target(original)
                if (resolved, _file_hash(resolved)) != watched:
                    raise ConfigurationError(
                        "A validation input changed before activation.",
                        code="AGENT_DEVELOPMENT_CONFLICT",
                    )
            pending: list[dict[str, Any]] = []
            patches: list[str] = []
            for index, (raw_path, data) in enumerate(writes.items()):
                path = self._target(raw_path)
                old = path.read_bytes() if path.is_file() else None
                old_hash = _digest(old) if old is not None else None
                expected_state = self._watched.get(Path(os.path.abspath(raw_path)))
                if expected_state is not None and (path, old_hash) != expected_state:
                    raise ConfigurationError(
                        "Source changed while candidate was being validated.",
                        code="AGENT_DEVELOPMENT_CONFLICT",
                    )
                backup = f"before-{index}.bin"
                if old is not None:
                    _atomic_write(self.directory / backup, old)
                patches.extend(
                    difflib.unified_diff(
                        (old or b"").decode("utf-8", errors="replace").splitlines(keepends=True),
                        data.decode("utf-8", errors="replace").splitlines(keepends=True),
                        fromfile=str(path),
                        tofile=str(path),
                    )
                )
                pending.append(
                    {
                        "path": str(path),
                        "backup": backup,
                        "old_sha256": old_hash,
                        "new_sha256": _digest(data),
                    }
                )
            self.state["writes"] = pending
            _atomic_write(self.directory / "changes.diff", "".join(patches).encode())
            self.phase("activating")
            try:
                for item, data in track(
                    list(zip(pending, writes.values(), strict=True)),
                    description="Activate validated tool",
                    disable=not self.verbose,
                ):
                    path = Path(item["path"])
                    if self._target(path) != path or _file_hash(path) != item["old_sha256"]:
                        raise ConfigurationError(
                            "Activation target changed after journaling.",
                            code="AGENT_DEVELOPMENT_CONFLICT",
                        )
                    _atomic_write(path, data)
            except Exception:
                self._restore()
                self.phase("rejected")
                raise
            result = self._receipt()
            self.complete(result)
            return result

    def complete(self, result: dict[str, Any]) -> None:
        self.state.update(phase="activated", result=result)
        self._save()

    def reject(self, reason: str, **data: Any) -> dict[str, Any]:
        if self.state.get("phase") == "activated":
            # A post-commit reporting error must not erase the durable commit receipt.
            return dict(self.state["result"])
        if self.state.get("phase") in {"activating", "rolling_back"}:
            # Keep the journal recoverable when an I/O error also prevented rollback.
            return {
                "status": "uncertain",
                "reason": reason,
                "change_id": self.change_id,
                "audit_dir": str(self.directory),
                **data,
            }
        result = {
            "status": "rejected",
            "reason": reason,
            "change_id": self.change_id,
            "audit_dir": str(self.directory),
            **data,
        }
        self.state.update(phase="rejected", result=result)
        self._save()
        return result


def rollback_tool_change(
    change_id: str, *, allow_side_effects: bool = False, verbose: bool = True
) -> dict[str, object]:
    """Restore a development transaction if no later edit overlaps its files.

    Use the change_id saved in a generate_tool/optimize_tool observation. Revert
    dependent changes in reverse order; running independent MCP servers must restart.
    """
    if allow_side_effects is not True:
        raise ConfigurationError(
            "Rollback requires allow_side_effects=True.", code="AGENT_TOOL_CONFIRMATION_REQUIRED"
        )
    if not isinstance(change_id, str) or not re.fullmatch(r"[a-z_]+-[a-f0-9]{24}", change_id):
        raise ConfigurationError("Invalid change_id.", code="AGENT_DEVELOPMENT_INVALID")
    root = _source_root()
    with _write_lock(root):
        directory = root / ".dnakit/development" / change_id
        state = json.loads((directory / "state.json").read_text())
        job = object.__new__(DevelopmentJob)
        job.root, job.directory, job.state = root, directory, state
        job._state_hash = _file_hash(directory / "state.json")
        job.change_id, job.verbose = change_id, verbose
        if state.get("phase") == "rolled_back":
            _clear_development_caches()
            return {"status": "rolled_back", "change_id": change_id}
        if state.get("phase") == "rolling_back":
            intact = all(
                _file_hash(job._target(Path(item["path"])))
                in (item["old_sha256"], item["new_sha256"])
                for item in state["writes"]
            )
        else:
            intact = state.get("phase") == "activated" and job._matches("new_sha256")
        if not intact:
            raise ConfigurationError(
                "Files changed after activation; rollback would overwrite later work.",
                code="AGENT_DEVELOPMENT_CONFLICT",
            )
        # Persist an inverse journal so interrupted rollbacks can be diagnosed.
        state["phase"] = "rolling_back"
        job._save()
        job._restore()
        _refresh_rolled_back_code(state)
        _clear_development_caches()
        state["phase"] = "rolled_back"
        job._save()
    if verbose:
        print(f"[DNAKit development] Rolled back {change_id}", flush=True)
    return {"status": "rolled_back", "change_id": change_id, "audit_dir": str(directory)}


def _clear_development_caches() -> None:
    from dnakit.skills.catalog import _catalog

    from .registry import default_tool_registry

    default_tool_registry.cache_clear()
    _catalog.cache_clear()


def _refresh_rolled_back_code(state: dict[str, Any]) -> None:
    """Refresh live exports before committing the rollback completion marker."""
    from dnakit import generated

    generated._refresh()
    if state.get("kind") != "repair":
        return
    result = state.get("result", {})
    tool_name = str(result.get("tool_name", state["identity"]))
    if tool_name.startswith("dnakit.generated."):
        return
    import inspect

    from .runtime_updates import refresh_repaired_function

    module_name, _, export_name = tool_name.rpartition(".")
    function = inspect.unwrap(getattr(importlib.import_module(module_name), export_name))
    source_path = inspect.getsourcefile(function)
    if source_path:
        refresh_repaired_function(tool_name, source_path, function.__name__)
