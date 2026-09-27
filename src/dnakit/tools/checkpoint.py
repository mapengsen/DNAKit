"""File-backed checkpoints for explicitly requested Agent run recovery.

Snapshots cover identifiable explicit files, bounded directory trees, tool code,
and declared Python artifacts. They cannot discover arbitrary Python's implicit
dependencies or make external effects atomic with writing a checkpoint. Pending
actions therefore remain explicitly uncertain after interruption.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import os
import sys
import tempfile
import time
from contextlib import suppress
from copy import deepcopy
from pathlib import Path
from types import CodeType
from typing import Any

from dnakit._progress import track
from dnakit._verbosity import progress_enabled
from dnakit._version import __version__
from dnakit.exceptions import ConfigurationError

from .execution import PythonExecutionConfig
from .registry import AgentToolRegistry

_FORMAT = 1
_DEVELOPMENT_ACTIONS = frozenset({"generate_tool", "optimize_tool"})
_MAX_TREE_ENTRIES = 10_000
_MAX_CHECKPOINT_BYTES = 256_000_000
_OUTPUT_KEYS = frozenset(
    {
        "output",
        "output_path",
        "output_dir",
        "output_file",
        "output_folder",
        "save_path",
        "destination",
        "destination_path",
        "destination_dir",
        "dest",
        "out_dir",
    }
)


def _json(value: object) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def _digest(value: object) -> str:
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reject_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON number: {value}")


def _tree_error(error: OSError) -> None:
    raise error


def _code_structure(code: CodeType) -> dict[str, Any]:
    def constant(value: object) -> object:
        if isinstance(value, CodeType):
            return _code_structure(value)
        if isinstance(value, (tuple, frozenset)):
            values = [constant(item) for item in value]
            return sorted(values, key=repr) if isinstance(value, frozenset) else values
        # Constants are compiler literals; repr has stable contents, unlike the
        # reference/interning flags emitted by marshal.dumps(code).
        return {"type": type(value).__name__, "value": repr(value)}

    return {
        "bytecode": code.co_code.hex(),
        "constants": [constant(v) for v in code.co_consts],
        "names": code.co_names,
        "variables": code.co_varnames,
        "freevars": code.co_freevars,
        "cellvars": code.co_cellvars,
        "arguments": [code.co_argcount, code.co_posonlyargcount, code.co_kwonlyargcount],
        "flags": code.co_flags,
    }


class _AgentCheckpoint:
    """One run directory, with a validated prefix copied from the latest run.

    The caller keeps its own history copy and invokes ``record`` after appending
    each observation. Only successful observations populate ``calls``. Runtime
    credentials and model configuration are deliberately absent from this API.
    Activated tool-development history validates its declared artifacts, without
    granting permission to perform new development or counting as tool execution.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        request: str,
        *,
        registry: AgentToolRegistry,
        code_execution: PythonExecutionConfig | None,
        allow_side_effects: bool,
        reuse_results: bool,
        verbose: bool = True,
        context: object = None,
    ) -> None:
        if not isinstance(cache_dir, (str, Path)) or not str(cache_dir).strip():
            raise ConfigurationError(
                "cache_dir must be a non-empty filesystem path.", code="AGENT_CACHE_DIR_INVALID"
            )
        if not isinstance(reuse_results, bool):
            raise ConfigurationError(
                "reuse_results must be boolean.", code="AGENT_REUSE_RESULTS_INVALID"
            )
        self._registry = registry
        self._code_execution = code_execution
        self._allow_side_effects = allow_side_effects
        self._cwd = Path.cwd().resolve()
        self._workspace = (
            Path(code_execution.workspace).expanduser().resolve()
            if code_execution is not None
            else None
        )
        try:
            self._cache_root = Path(cache_dir).expanduser().resolve()
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(
                "cache_dir is not a usable filesystem path.", code="AGENT_CACHE_DIR_INVALID"
            ) from exc
        self._identity = {
            "request": request,
            "cwd": str(self._cwd),
            "workspace": str(self._workspace) if self._workspace is not None else None,
            "knowledge_fingerprint": _digest(context),
        }
        self.history: list[dict[str, Any]] = []
        self.calls: dict[str, dict[str, Any]] = {}
        self._proofs: list[dict[str, Any]] = []
        self.pending: dict[str, Any] | None = None
        self.uncertain_actions: list[dict[str, Any]] = []
        self.reused_steps = 0
        self.invalidated_from: int | None = None
        self.source_directory: Path | None = None
        self._status = "running"
        self._text = ""
        task_directory = self._cache_root / "agent_runs" / _digest(self._identity)
        try:
            task_directory.mkdir(parents=True, exist_ok=True)
            candidates = list(task_directory.glob("run_*/checkpoint.json")) if reuse_results else []
            if candidates:
                source = max(
                    candidates, key=lambda path: (path.stat().st_mtime_ns, path.parent.name)
                )
                self._restore(source, verbose=verbose)
            self.directory = Path(
                tempfile.mkdtemp(prefix=f"run_{time.time_ns()}_", dir=task_directory)
            )
        except OSError as exc:
            raise ConfigurationError(
                "Cannot create or read the Agent checkpoint directory.",
                code="AGENT_CHECKPOINT_IO_ERROR",
                context={"directory": str(task_directory)},
            ) from exc
        self._save()

    def update_registry(self, registry: AgentToolRegistry) -> None:
        """Use the runner's refreshed catalog for subsequent tool fingerprints.

        Existing proofs retain their original fingerprints, so a changed tool
        definition still invalidates its earlier calls on the next restoration.
        """
        if not isinstance(registry, AgentToolRegistry):
            raise ConfigurationError(
                "registry must be AgentToolRegistry.", code="AGENT_REGISTRY_INVALID"
            )
        self._registry = registry

    def _development_paths(self, result: object) -> tuple[Path, ...] | None:
        if not isinstance(result, dict) or result.get("status") != "activated":
            return None
        artifacts = result.get("artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            return None
        paths: list[Path] = []
        for name in artifacts:
            if not isinstance(name, str) or not name.strip() or "\0" in name or "://" in name:
                return None
            try:
                path = Path(name).expanduser()
                paths.append((path if path.is_absolute() else self._cwd / path).absolute())
            except (OSError, ValueError, RuntimeError):
                return None
        return tuple(paths)

    def _tool_fingerprint(self, name: str) -> str | None:
        try:
            tool = self._registry.resolve(name)
            function = inspect.unwrap(tool.function)
            source_file = inspect.getsourcefile(function)
            source = (
                _file_digest(Path(source_file))
                if source_file and Path(source_file).is_file()
                else None
            )
            code = getattr(function, "__code__", None)
            bytecode = _code_structure(code) if isinstance(code, CodeType) else None
            if source is None and bytecode is None:
                return None
            metadata: object = None
            if tool.name in {
                "dnakit.integrations.list_mcp_tools",
                "dnakit.integrations.call_mcp_tool",
            }:
                from .descriptions import read_descriptions

                metadata = read_descriptions()["mcp"]
            elif tool.name.startswith("dnakit.skills."):
                from dnakit.skills.catalog import _locations

                metadata = _file_digest(_locations()[1])
            return _digest(
                {
                    "package_version": __version__,
                    "manifest": tool.to_dict(),
                    "source": source,
                    "bytecode": bytecode,
                    **({"local_metadata": metadata} if metadata is not None else {}),
                }
            )
        except (ConfigurationError, OSError, TypeError, ValueError):
            return None

    def _ignored_tree(self, path: Path) -> bool:
        resolved = path.resolve()
        audit = (self._workspace or self._cwd) / ".dnakit_agent"
        return resolved.is_relative_to(self._cache_root) or resolved.is_relative_to(audit)

    def _snapshot(self, path: Path) -> dict[str, Any]:
        try:
            if path.is_file():
                before = path.stat()
                digest = _file_digest(path)
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    return {"kind": "unstable"}
                return {
                    "kind": "file",
                    "size": after.st_size,
                    "sha256": digest,
                    "resolved": str(path.resolve()),
                }
            if path.is_dir():
                if self._ignored_tree(path):
                    return {"kind": "unbounded_directory"}
                members: dict[str, Any] = {}
                directories: list[str] = []
                visited = 0
                for directory, names, files in os.walk(
                    path,
                    followlinks=False,
                    onerror=_tree_error,
                ):
                    base = Path(directory)
                    if self._ignored_tree(base):
                        names[:] = []
                        continue
                    visited += len(names) + len(files) + 1
                    if visited > _MAX_TREE_ENTRIES:
                        return {"kind": "unbounded_directory"}
                    names[:] = sorted(name for name in names if not self._ignored_tree(base / name))
                    if any((base / name).is_symlink() for name in names):
                        return {"kind": "unbounded_directory"}
                    directories.extend(str((base / name).relative_to(path)) for name in names)
                    for name in sorted(files):
                        child = base / name
                        if self._ignored_tree(child):
                            continue
                        if not child.is_file():
                            return {"kind": "special_file"}
                        value = self._snapshot(child)
                        if value.get("kind") != "file":
                            return value
                        members[str(child.relative_to(path))] = value
                return {
                    "kind": "directory",
                    "members": members,
                    "directories": directories,
                    "resolved": str(path.resolve()),
                }
            return {"kind": "missing"}
        except (OSError, ValueError, RuntimeError):
            return {"kind": "unreadable"}

    def _files(self, value: object, *, base: Path, inputs: bool = False) -> dict[str, Any]:
        paths: set[Path] = set()

        def visit(item: object) -> None:
            if isinstance(item, dict):
                for key, child in item.items():
                    if not (inputs and str(key).casefold() in _OUTPUT_KEYS):
                        visit(child)
            elif isinstance(item, (list, tuple)):
                for child in item:
                    visit(child)
            elif isinstance(item, (str, Path)):
                text = str(item)
                if not text or len(text) > 4096 or any(c in text for c in ("\n", "\0", "://")):
                    return
                try:
                    path = Path(text).expanduser()
                    path = path if path.is_absolute() else base / path
                    if path.exists():
                        paths.add(path.absolute())
                except (OSError, ValueError, RuntimeError):
                    pass

        visit(value)
        return {str(path): self._snapshot(path) for path in sorted(paths)}

    def _output_paths(self, value: object, *, base: Path) -> list[str]:
        paths: set[str] = set()

        def visit(item: object, output: bool = False) -> None:
            if isinstance(item, dict):
                for key, child in item.items():
                    visit(child, output or str(key).casefold() in _OUTPUT_KEYS)
            elif isinstance(item, (list, tuple)):
                for child in item:
                    visit(child, output)
            elif output and isinstance(item, str) and item and "://" not in item:
                try:
                    path = Path(item).expanduser()
                    paths.add(str(path if path.is_absolute() else base / path))
                except (OSError, ValueError, RuntimeError):
                    pass

        visit(value)
        return sorted(paths)

    def _files_valid(self, files: object, *, outputs: bool = False) -> bool:
        if not isinstance(files, dict):
            return False
        for name, previous in files.items():
            if not isinstance(name, str) or not isinstance(previous, dict):
                return False
            path = Path(name)
            if outputs and previous.get("kind") == "directory":
                # Later workflow steps may add reports under an output directory.
                # Existing recorded artifacts must still match, but additions do
                # not invalidate a step that already produced its own files.
                if not path.is_dir() or str(path.resolve()) != previous.get("resolved"):
                    return False
                members = previous.get("members")
                if not isinstance(members, dict) or any(
                    self._snapshot(path / child) != value for child, value in members.items()
                ):
                    return False
                if any(not (path / child).is_dir() for child in previous.get("directories", [])):
                    return False
            elif (
                previous.get("kind") not in {"file", "directory"}
                or self._snapshot(path) != previous
            ):
                return False
        return True

    def _valid(self, entry: dict[str, Any], proof: dict[str, Any]) -> bool:
        action, observation = entry["action"], entry["observation"]
        kind = action.get("action")
        if observation.get("ok") is not True:
            return True
        if kind in _DEVELOPMENT_ACTIONS:
            # Activation is already applied to the source tree. Reading its
            # receipt does not enable another generation/optimization action.
            paths = self._development_paths(observation.get("result"))
            outputs = proof.get("outputs")
            return (
                paths is not None
                and isinstance(outputs, dict)
                and set(outputs) == {str(path) for path in paths}
                and self._files_valid(proof.get("inputs"))
                and self._files_valid(outputs, outputs=True)
            )
        if kind == "search":
            try:
                for tool in observation.get("tools", []):
                    self._registry.resolve(tool["name"])
            except (ConfigurationError, KeyError, TypeError):
                return False
            return True
        if kind in {"describe", "call", "inspect_tool_source"}:
            name = action.get("tool_name")
            if not isinstance(name, str):
                return False
            fingerprint = self._tool_fingerprint(name)
            if fingerprint is None or fingerprint != proof.get("tool_fingerprint"):
                return False
            tool = self._registry.resolve(name)
            if kind == "call" and (
                not tool.agent_compatible
                or (tool.requires_confirmation and not self._allow_side_effects)
            ):
                return False
        if kind == "inspect_tool_source":
            result = observation.get("result")
            if not isinstance(result, dict) or result.get("name") != action["tool_name"]:
                return False
            path, digest = result.get("source_path"), result.get("base_sha256")
            outputs = proof.get("outputs")
            return (
                isinstance(path, str)
                and isinstance(digest, str)
                and isinstance(outputs, dict)
                and set(outputs) == {path}
                and isinstance(outputs[path], dict)
                and outputs[path].get("sha256") == digest
                and self._files_valid(outputs)
            )
        if kind == "execute_python" and (
            self._code_execution is None
            or not self._allow_side_effects
            or proof.get("runtime") != self._runtime()
        ):
            return False
        if kind in {"call", "execute_python"}:
            if proof.get("authorized_side_effects") and not self._allow_side_effects:
                return False
            return self._files_valid(proof.get("inputs")) and self._files_valid(
                proof.get("outputs"),
                outputs=True,
            )
        return True

    @staticmethod
    def _runtime() -> dict[str, str]:
        return {
            "python": sys.version,
            "executable": sys.executable,
            "package_version": __version__,
        }

    def _restore(self, source: Path, *, verbose: bool) -> None:
        try:
            if source.stat().st_size > _MAX_CHECKPOINT_BYTES:
                raise ValueError("Checkpoint exceeds the size limit")
            payload = json.loads(
                source.read_text(encoding="utf-8"), parse_constant=_reject_constant
            )
            if (
                not isinstance(payload, dict)
                or payload.get("format") != _FORMAT
                or payload.get("identity") != self._identity
                or not isinstance(payload.get("history"), list)
                or not isinstance(payload.get("proofs"), list)
                or len(payload["history"]) != len(payload["proofs"])
                or not isinstance(payload.get("uncertain_actions", []), list)
            ):
                raise ValueError("Invalid checkpoint structure or version")
            for number, (entry, proof) in enumerate(
                zip(payload["history"], payload["proofs"], strict=True),
                1,
            ):
                if (
                    not isinstance(entry, dict)
                    or type(entry.get("step")) is not int
                    or entry["step"] != number
                    or not isinstance(entry.get("action"), dict)
                    or not isinstance(entry.get("observation"), dict)
                    or not isinstance(proof, dict)
                ):
                    raise ValueError("Invalid checkpoint history")
                action, observation = entry["action"], entry["observation"]
                kind = action.get("action")
                if not isinstance(kind, str) or not isinstance(observation.get("ok"), bool):
                    raise ValueError("Invalid action or observation")
                if observation["ok"]:
                    if kind in {"call", "describe", "inspect_tool_source"} and not isinstance(
                        action.get("tool_name"), str
                    ):
                        raise ValueError("Invalid tool action")
                    if kind == "search":
                        tools = observation.get("tools")
                        if not isinstance(tools, list) or any(
                            not isinstance(tool, dict) or not isinstance(tool.get("name"), str)
                            for tool in tools
                        ):
                            raise ValueError("Invalid tool discovery results")
                    if kind == "inspect_tool_source":
                        result = observation.get("result")
                        if not isinstance(result, dict) or any(
                            not isinstance(result.get(key), str)
                            for key in ("name", "source_path", "base_sha256", "source")
                        ):
                            raise ValueError("Invalid inspected source evidence")
                    if (
                        kind
                        in {"call", "execute_python", "inspect_tool_source"} | _DEVELOPMENT_ACTIONS
                    ):
                        if "result" not in observation:
                            raise ValueError("Successful action has no result")
                        for files in (proof.get("inputs"), proof.get("outputs")):
                            if not isinstance(files, dict):
                                raise ValueError("Invalid file validation evidence")
                            for name, snapshot in files.items():
                                if not isinstance(name, str) or not isinstance(snapshot, dict):
                                    raise ValueError("Invalid file snapshot")
                                if snapshot.get("kind") == "directory":
                                    members = snapshot.get("members")
                                    directories = snapshot.get("directories")
                                    if (
                                        not isinstance(members, dict)
                                        or any(
                                            not isinstance(child, str)
                                            or not isinstance(value, dict)
                                            for child, value in members.items()
                                        )
                                        or not isinstance(directories, list)
                                        or any(not isinstance(child, str) for child in directories)
                                    ):
                                        raise ValueError("Invalid directory snapshot")
                if proof.get("cache_key") is not None and not isinstance(proof["cache_key"], str):
                    raise ValueError("Invalid call cache key")
            uncertain = payload.get("uncertain_actions", [])
            pending = payload.get("pending")
            if pending is not None:
                uncertain = [*uncertain, pending]
            if any(
                not isinstance(item, dict)
                or not isinstance(item.get("action"), dict)
                or (item.get("cache_key") is not None and not isinstance(item["cache_key"], str))
                for item in uncertain
            ):
                raise ValueError("Invalid pending action")
        except (OSError, ValueError, TypeError, KeyError, RecursionError) as exc:
            raise ConfigurationError(
                "Agent checkpoint is unreadable or malformed; "
                "inspect it or use reuse_results=False.",
                code="AGENT_CHECKPOINT_INVALID",
                context={"path": str(source)},
            ) from exc
        self.source_directory = source.parent
        self.uncertain_actions = deepcopy(uncertain)
        entries = zip(payload["history"], payload["proofs"], strict=True)
        for entry, proof in track(
            entries,
            total=len(payload["history"]),
            description="Validating saved Agent steps",
            disable=not progress_enabled(verbose),
        ):
            if not self._valid(entry, proof):
                self.invalidated_from = entry["step"]
                break
            self.history.append(deepcopy(entry))
            self._proofs.append(deepcopy(proof))
            observation = entry["observation"]
            if (
                entry["action"].get("action") in {"call", "execute_python"}
                and observation.get("ok") is True
            ):
                self.reused_steps += 1
                key = proof.get("cache_key")
                if isinstance(key, str):
                    self.calls[key] = deepcopy(observation)

    def begin(
        self,
        action: dict[str, Any],
        *,
        inputs: object = None,
        cache_key: str | None = None,
    ) -> None:
        """Durably record execution intent before a tool or subprocess can act."""
        if self.pending is not None:
            raise ConfigurationError(
                "A previous action has not been recorded.", code="AGENT_CHECKPOINT_PENDING"
            )
        base = self._workspace if action.get("action") == "execute_python" else self._cwd
        self.pending = {
            "action": deepcopy(action),
            "cache_key": cache_key,
            "inputs": self._files(inputs, base=base or self._cwd, inputs=True),
            "output_paths": self._output_paths(inputs, base=base or self._cwd),
            "started_at_ns": time.time_ns(),
        }
        if action.get("action") == "call":
            self.pending["tool_fingerprint"] = self._tool_fingerprint(action["tool_name"])
        self._save()

    def record(self, entry: dict[str, Any], *, cache_key: str | None = None) -> None:
        """Append one completed observation and clear this run's pending action."""
        if entry.get("step") != len(self.history) + 1:
            raise ConfigurationError(
                "Checkpoint history steps must be consecutive.",
                code="AGENT_CHECKPOINT_STEP_INVALID",
            )
        action, observation = entry["action"], entry["observation"]
        kind = action.get("action")
        pending = self.pending if self.pending and self.pending.get("action") == action else {}
        proof: dict[str, Any] = {
            "cache_key": cache_key if cache_key is not None else pending.get("cache_key"),
            "inputs": deepcopy(pending.get("inputs", {})),
            "outputs": {},
            "authorized_side_effects": self._allow_side_effects,
        }
        if observation.get("ok") is True:
            if kind in {"call", "describe", "inspect_tool_source"}:
                proof["tool_fingerprint"] = pending.get(
                    "tool_fingerprint"
                ) or self._tool_fingerprint(action["tool_name"])
            result = observation.get("result")
            if kind == "inspect_tool_source":
                if not isinstance(result, dict) or not isinstance(result.get("source_path"), str):
                    raise ConfigurationError(
                        "Source inspection requires a source path and content hash.",
                        code="AGENT_REPAIR_SOURCE_INVALID",
                    )
                path = result["source_path"]
                proof["outputs"] = {path: self._snapshot(Path(path))}
            elif kind in _DEVELOPMENT_ACTIONS:
                artifacts = self._development_paths(result)
                if artifacts is None:
                    raise ConfigurationError(
                        "Successful tool development requires activated status and artifact paths.",
                        code="AGENT_DEVELOPMENT_RESULT_INVALID",
                    )
                # Do not scan audit_dir: append-only transaction logs can change
                # independently of the activated implementation and descriptions.
                # Explicit missing artifacts stay in the proof and prevent reuse.
                proof["outputs"] = {str(path): self._snapshot(path) for path in artifacts}
            elif kind == "call":
                proof["outputs"] = self._files(result, base=self._cwd)
                for name in pending.get("output_paths", []):
                    proof["outputs"][name] = self._snapshot(Path(name))
            elif kind == "execute_python" and isinstance(result, dict):
                paths = [
                    result.get(key)
                    for key in ("script_path", "inputs_path", "stdout_path", "stderr_path")
                ]
                paths.extend(action.get("expected_outputs", []))
                paths.append(result.get("outputs"))
                proof["outputs"] = self._files(paths, base=self._workspace or self._cwd)
                proof["runtime"] = self._runtime()
                for name in action.get("expected_outputs", []):
                    target = Path(name).expanduser()
                    if not target.is_absolute():
                        target = (self._workspace or self._cwd) / target
                    proof["outputs"][str(target)] = self._snapshot(target)
            # In-place file operations are validated against their resulting state.
            for path in proof["inputs"].keys() & proof["outputs"].keys():
                proof["inputs"][path] = proof["outputs"][path]
            key = proof["cache_key"]
            if kind in {"call", "execute_python"} and isinstance(key, str):
                self.calls[key] = deepcopy(observation)
            elif (
                kind in _DEVELOPMENT_ACTIONS and isinstance(key, str) and self._valid(entry, proof)
            ):
                # The development transaction reconciled its durable receipt.
                # Only verified activation resolves a matching interrupted job;
                # failed or missing-artifact results leave uncertainty intact.
                self.uncertain_actions = [
                    item
                    for item in self.uncertain_actions
                    if not (
                        item.get("cache_key") == key
                        and item.get("action", {}).get("action") == kind
                    )
                ]
        self.history.append(deepcopy(entry))
        self._proofs.append(proof)
        self.pending = None
        self._save()

    def finish(self, status: str, text: str = "") -> None:
        """Persist final status without model configuration or credentials."""
        self._status, self._text = status, text
        self._save()

    def _save(self) -> None:
        payload = {
            "format": _FORMAT,
            "identity": self._identity,
            "package_version": __version__,
            "history": self.history,
            "proofs": self._proofs,
            "pending": self.pending,
            "uncertain_actions": self.uncertain_actions,
            "status": self._status,
            "text": self._text,
            "updated_at_ns": time.time_ns(),
            "source_directory": str(self.source_directory) if self.source_directory else None,
        }
        temporary: str | None = None
        try:
            encoded = _json(payload).encode("utf-8")
            if len(encoded) > _MAX_CHECKPOINT_BYTES:
                raise ValueError("Checkpoint exceeds the size limit")
            descriptor, temporary = tempfile.mkstemp(prefix=".checkpoint_", dir=self.directory)
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(encoded)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.directory / "checkpoint.json")
        except (OSError, ValueError, TypeError) as exc:
            raise ConfigurationError(
                "Could not save the Agent checkpoint.",
                code="AGENT_CHECKPOINT_WRITE_FAILED",
                context={"directory": str(self.directory)},
            ) from exc
        finally:
            if temporary is not None:
                with suppress(OSError):
                    Path(temporary).unlink(missing_ok=True)
