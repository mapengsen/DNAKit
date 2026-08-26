"""Create, persist, and verify reproducible run manifests."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from dnakit.core import ArtifactRef, Issue, Provenance, RunManifest
from dnakit.exceptions import ConfigurationError, InputFormatError

RunStatus = Literal["running", "succeeded", "failed", "cancelled"]
_MAX_MANIFEST_BYTES = 100_000_000
_MAX_MANIFEST_DEPTH = 64
_MAX_MANIFEST_NODES = 1_000_000


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise InputFormatError(
                "Run manifest contains a duplicate JSON key.",
                code="INVALID_RUN_MANIFEST",
                context={"key": key},
            )
        result[key] = value
    return result


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _validate_manifest_shape(value: object) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(value, 0)]
    while stack:
        item, depth = stack.pop()
        nodes += 1
        if nodes > _MAX_MANIFEST_NODES or depth > _MAX_MANIFEST_DEPTH:
            raise InputFormatError(
                "Run manifest exceeds structural limits.",
                code="RUN_MANIFEST_STRUCTURE_LIMIT",
                context={
                    "max_nodes": _MAX_MANIFEST_NODES,
                    "max_depth": _MAX_MANIFEST_DEPTH,
                },
            )
        if isinstance(item, Mapping):
            stack.extend((key, depth + 1) for key in item)
            stack.extend((child, depth + 1) for child in item.values())
        elif isinstance(item, list):
            stack.extend((child, depth + 1) for child in item)


def artifact_from_path(
    path: str | os.PathLike[str],
    *,
    media_type: str,
    schema_version: str,
    max_bytes: int | None = None,
) -> ArtifactRef:
    """Hash one regular file without copying its contents into the manifest."""

    if not isinstance(media_type, str) or not media_type.strip():
        raise ConfigurationError("media_type must be a non-empty string.")
    if not isinstance(schema_version, str) or not schema_version.strip():
        raise ConfigurationError("schema_version must be a non-empty string.")
    if max_bytes is not None and (
        isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes < 1
    ):
        raise ConfigurationError("max_bytes must be a positive integer or None.")
    try:
        resolved = Path(path).expanduser().absolute()
    except TypeError as exc:
        raise ConfigurationError("Artifact path must be a filesystem path.") from exc
    if not resolved.is_file() or resolved.is_symlink():
        raise InputFormatError(
            "Artifact path must be a regular non-symlink file.",
            code="INVALID_ARTIFACT_PATH",
            context={"path": str(resolved)},
        )
    before = resolved.stat()
    if max_bytes is not None and before.st_size > max_bytes:
        raise InputFormatError(
            "Artifact exceeds max_bytes.",
            code="ARTIFACT_SIZE_LIMIT",
            context={"path": str(resolved), "byte_size": before.st_size, "max_bytes": max_bytes},
        )
    digest = hashlib.sha256()
    byte_size = 0
    with resolved.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            byte_size += len(chunk)
            if max_bytes is not None and byte_size > max_bytes:
                raise InputFormatError(
                    "Artifact grew beyond max_bytes while being hashed.",
                    code="ARTIFACT_SIZE_LIMIT",
                    context={"path": str(resolved), "max_bytes": max_bytes},
                )
            digest.update(chunk)
    stat = resolved.stat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
    if identity_before != identity_after or byte_size != stat.st_size:
        raise InputFormatError(
            "Artifact changed while it was being hashed.",
            code="ARTIFACT_CHANGED_DURING_HASH",
            context={"path": str(resolved)},
        )
    return ArtifactRef(
        relative_path=os.path.relpath(resolved, Path.cwd()),
        media_type=media_type,
        schema_version=schema_version,
        sha256=digest.hexdigest(),
        byte_size=byte_size,
        created_at=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
    )


@dataclass(slots=True)
class RunManifestBuilder:
    """Mutable builder whose public product is an immutable RunManifest."""

    run_id: str
    command: tuple[str, ...]
    resolved_config: Mapping[str, object]
    seed: int | None = None
    seed_derivation: str | None = None
    provenance: Provenance = field(default_factory=Provenance)
    started_at: str = field(default_factory=_utc_now)
    inputs: list[ArtifactRef] = field(default_factory=list)
    outputs: list[ArtifactRef] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)

    def __post_init__(self) -> None:
        RunManifest(
            self.run_id,
            self.command,
            self.resolved_config,
            self.provenance,
            self.started_at,
            "running",
            seed=self.seed,
            seed_derivation=self.seed_derivation,
        )

    def add_input(self, artifact: ArtifactRef) -> None:
        if not isinstance(artifact, ArtifactRef):
            raise ConfigurationError("artifact must be ArtifactRef.")
        self.inputs.append(artifact)

    def add_output(self, artifact: ArtifactRef) -> None:
        if not isinstance(artifact, ArtifactRef):
            raise ConfigurationError("artifact must be ArtifactRef.")
        self.outputs.append(artifact)

    def add_issue(self, issue: Issue) -> None:
        if not isinstance(issue, Issue):
            raise ConfigurationError("issue must be Issue.")
        self.issues.append(issue)

    def build(
        self,
        *,
        status: RunStatus = "running",
        finished_at: str | None = None,
    ) -> RunManifest:
        if status not in ("running", "succeeded", "failed", "cancelled"):
            raise ConfigurationError("Unknown run status.", code="INVALID_RUN_STATUS")
        if status == "running" and finished_at is not None:
            raise ConfigurationError(
                "A running manifest cannot have finished_at.", code="INVALID_RUN_TIMESTAMPS"
            )
        if status != "running" and finished_at is None:
            finished_at = _utc_now()
        return RunManifest(
            self.run_id,
            self.command,
            self.resolved_config,
            self.provenance,
            self.started_at,
            status,
            seed=self.seed,
            seed_derivation=self.seed_derivation,
            inputs=tuple(self.inputs),
            outputs=tuple(self.outputs),
            finished_at=finished_at,
            issues=tuple(self.issues),
        )


def save_manifest(
    manifest: RunManifest,
    path: str | os.PathLike[str],
    *,
    overwrite: bool = False,
) -> ArtifactRef:
    """Persist a manifest atomically as canonical JSON."""

    if not isinstance(manifest, RunManifest):
        raise ConfigurationError("manifest must be RunManifest.")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be boolean.")
    target = Path(path)
    if target.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite manifest: {target}")
    if not target.parent.is_dir():
        raise FileNotFoundError(f"Manifest parent does not exist: {target.parent}")
    payload = manifest.to_dict()
    _validate_manifest_shape(payload)
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    if len(encoded) > _MAX_MANIFEST_BYTES:
        raise ConfigurationError(
            "Run manifest exceeds the byte limit.",
            code="RUN_MANIFEST_SIZE_LIMIT",
            context={"byte_size": len(encoded), "max_bytes": _MAX_MANIFEST_BYTES},
        )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, target)
        else:
            os.link(temporary, target)
            temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return artifact_from_path(
        target,
        media_type="application/vnd.dnakit.run-manifest+json",
        schema_version="dnakit-run-manifest-v1",
    )


def load_manifest(path: str | os.PathLike[str]) -> dict[str, object]:
    """Load and minimally validate a persisted manifest without executing it."""

    try:
        source = Path(path).expanduser().absolute()
    except TypeError as exc:
        raise InputFormatError(
            "Run manifest path must be a filesystem path.", code="INVALID_RUN_MANIFEST"
        ) from exc
    if source.is_symlink() or not source.is_file():
        raise InputFormatError(
            "Run manifest path must be a regular non-symlink file.",
            code="INVALID_RUN_MANIFEST",
            context={"path": str(source)},
        )
    byte_size = source.stat().st_size
    if byte_size > _MAX_MANIFEST_BYTES:
        raise InputFormatError(
            "Run manifest exceeds the byte limit.",
            code="RUN_MANIFEST_SIZE_LIMIT",
            context={"byte_size": byte_size, "max_bytes": _MAX_MANIFEST_BYTES},
        )
    try:
        with source.open("rb") as handle:
            raw = handle.read(_MAX_MANIFEST_BYTES + 1)
        if len(raw) > _MAX_MANIFEST_BYTES:
            raise InputFormatError(
                "Run manifest grew beyond the byte limit while being read.",
                code="RUN_MANIFEST_SIZE_LIMIT",
                context={"max_bytes": _MAX_MANIFEST_BYTES},
            )
        payload: object = json.loads(raw.decode("utf-8"), object_pairs_hook=_unique_json_object)
    except InputFormatError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise InputFormatError(
            "Run manifest cannot be decoded.",
            code="INVALID_RUN_MANIFEST",
            context={"path": str(source)},
        ) from exc
    _validate_manifest_shape(payload)
    if not isinstance(payload, dict):
        raise InputFormatError("Run manifest root must be an object.", code="INVALID_RUN_MANIFEST")
    required = {"run_id", "command", "resolved_config", "provenance", "started_at", "status"}
    missing = sorted(required - payload.keys())
    if missing:
        raise InputFormatError(
            "Run manifest is missing required fields.",
            code="INVALID_RUN_MANIFEST",
            context={"missing": missing},
        )
    if not isinstance(payload["run_id"], str) or not payload["run_id"].strip():
        raise InputFormatError("Run manifest run_id is invalid.", code="INVALID_RUN_MANIFEST")
    command = payload["command"]
    if (
        not isinstance(command, list)
        or not command
        or any(not isinstance(argument, str) or not argument for argument in command)
    ):
        raise InputFormatError("Run manifest command is invalid.", code="INVALID_RUN_MANIFEST")
    if not isinstance(payload["resolved_config"], dict) or not isinstance(
        payload["provenance"], dict
    ):
        raise InputFormatError(
            "Run manifest configuration or provenance is invalid.",
            code="INVALID_RUN_MANIFEST",
        )
    if not isinstance(payload["started_at"], str) or not payload["started_at"].strip():
        raise InputFormatError("Run manifest started_at is invalid.", code="INVALID_RUN_MANIFEST")
    status = payload["status"]
    if status not in ("running", "succeeded", "failed", "cancelled"):
        raise InputFormatError("Run manifest status is invalid.", code="INVALID_RUN_MANIFEST")
    finished_at = payload.get("finished_at")
    if finished_at is not None and (not isinstance(finished_at, str) or not finished_at.strip()):
        raise InputFormatError("Run manifest finished_at is invalid.", code="INVALID_RUN_MANIFEST")
    if (status == "running") != (finished_at is None):
        raise InputFormatError(
            "Run manifest status and finished_at are inconsistent.",
            code="INVALID_RUN_MANIFEST",
        )
    seed = payload.get("seed")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise InputFormatError("Run manifest seed is invalid.", code="INVALID_RUN_MANIFEST")
    seed_derivation = payload.get("seed_derivation")
    if seed_derivation is not None and (
        not isinstance(seed_derivation, str) or not seed_derivation.strip()
    ):
        raise InputFormatError(
            "Run manifest seed_derivation is invalid.", code="INVALID_RUN_MANIFEST"
        )
    for field_name in ("inputs", "outputs", "issues"):
        value = payload.get(field_name, [])
        if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
            raise InputFormatError(
                f"Run manifest {field_name} is invalid.", code="INVALID_RUN_MANIFEST"
            )
    return {str(key): value for key, value in payload.items()}


__all__ = ["RunManifestBuilder", "artifact_from_path", "load_manifest", "save_manifest"]
