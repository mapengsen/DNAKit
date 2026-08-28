"""Small safe cache for JSON-compatible DNAKit results."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections.abc import Mapping
from contextlib import suppress
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, cast

from dnakit.core import ArtifactRef
from dnakit.core._json import freeze_json, to_json_compatible
from dnakit.exceptions import CacheError, ConfigurationError

_NAMESPACE = re.compile(r"^[a-zA-Z0-9_.-]+$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_ROOT_MARKER = ".dnakit-cache-v1"
_NAMESPACE_MARKER = ".dnakit-cache-namespace-v1"
DEFAULT_MAX_ENTRY_BYTES = 20_000_000
_MAX_CACHE_DEPTH = 64
_MAX_CACHE_NODES = 1_000_000


def _dnakit_version() -> str:
    try:
        return version("dnakit")
    except PackageNotFoundError:
        return "unknown"


def _validate_json_shape(value: object) -> None:
    nodes = [0]
    active: set[int] = set()

    def visit(item: object, depth: int) -> None:
        nodes[0] += 1
        if depth > _MAX_CACHE_DEPTH or nodes[0] > _MAX_CACHE_NODES:
            raise ConfigurationError(
                "Cache value exceeds structural limits.", code="CACHE_STRUCTURE_LIMIT"
            )
        if isinstance(item, (Mapping, tuple, list)):
            identity = id(item)
            if identity in active:
                raise ConfigurationError(
                    "Cache value contains a recursive object.", code="CACHE_STRUCTURE_LIMIT"
                )
            active.add(identity)
            try:
                if isinstance(item, Mapping):
                    for key in item:
                        visit(key, depth + 1)
                        visit(item[key], depth + 1)
                else:
                    for child in item:
                        visit(child, depth + 1)
            finally:
                active.remove(identity)

    visit(value, 0)


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise CacheError(
                "Cache entry contains a duplicate JSON key.",
                code="CACHE_DECODE_ERROR",
                context={"key": key},
            )
        result[key] = value
    return result


def _canonical_bytes(value: object) -> bytes:
    try:
        _validate_json_shape(value)
        compatible = to_json_compatible(freeze_json(value))
        return json.dumps(
            compatible,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except RecursionError as exc:
        raise ConfigurationError(
            "Cache values exceed structural limits.", code="CACHE_STRUCTURE_LIMIT"
        ) from exc
    except (TypeError, ValueError, ConfigurationError) as exc:
        raise ConfigurationError(
            "Cache values must be finite JSON-compatible data.", code="INVALID_CACHE_VALUE"
        ) from exc


@dataclass(frozen=True, slots=True)
class CacheKey:
    """Versioned content hash for all computation-defining inputs."""

    namespace: str
    digest: str
    schema_version: str = "dnakit-cache-key-v1"

    def __post_init__(self) -> None:
        if not isinstance(self.namespace, str) or not _NAMESPACE.fullmatch(self.namespace):
            raise ConfigurationError(
                "Cache namespace must be a safe non-empty identifier.",
                code="INVALID_CACHE_NAMESPACE",
            )
        if not isinstance(self.digest, str) or not _DIGEST.fullmatch(self.digest):
            raise ConfigurationError(
                "Cache digest must contain 64 lowercase hexadecimal digits.",
                code="INVALID_CACHE_DIGEST",
            )
        if not isinstance(self.schema_version, str) or not self.schema_version.strip():
            raise ConfigurationError("Cache key schema_version must be non-empty.")

    @classmethod
    def from_components(
        cls,
        namespace: str,
        components: Mapping[str, object],
        *,
        schema_version: str = "dnakit-cache-key-v1",
    ) -> CacheKey:
        if not isinstance(namespace, str) or not _NAMESPACE.fullmatch(namespace):
            raise ConfigurationError(
                "Cache namespace must be a safe non-empty identifier.",
                code="INVALID_CACHE_NAMESPACE",
            )
        if not isinstance(components, Mapping):
            raise ConfigurationError(
                "Cache components must be a mapping.", code="INVALID_CACHE_COMPONENTS"
            )
        if not isinstance(schema_version, str) or not schema_version.strip():
            raise ConfigurationError("Cache key schema_version must be non-empty.")
        payload = _canonical_bytes(
            {
                "schema_version": schema_version,
                "dnakit_version": _dnakit_version(),
                "namespace": namespace,
                "components": components,
            }
        )
        return cls(namespace, hashlib.sha256(payload).hexdigest(), schema_version)


@dataclass(frozen=True, slots=True)
class CacheClearReport:
    """Audited count and bytes removed from one cache namespace or all namespaces."""

    namespace: str | None
    removed_count: int
    removed_bytes: int

    def __post_init__(self) -> None:
        if self.namespace is not None and (
            not isinstance(self.namespace, str) or not _NAMESPACE.fullmatch(self.namespace)
        ):
            raise ConfigurationError("Cache report namespace is invalid.")
        for name in ("removed_count", "removed_bytes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Cache report {name} must be non-negative.")


class JSONCache:
    """Filesystem cache using canonical JSON and atomic same-directory replacement."""

    def __init__(
        self,
        root: str | os.PathLike[str],
        *,
        max_entry_bytes: int = DEFAULT_MAX_ENTRY_BYTES,
    ) -> None:
        if (
            isinstance(max_entry_bytes, bool)
            or not isinstance(max_entry_bytes, int)
            or not 1 <= max_entry_bytes <= 1_000_000_000
        ):
            raise ConfigurationError(
                "max_entry_bytes must be an integer in [1, 1000000000].",
                code="INVALID_CACHE_ENTRY_LIMIT",
            )
        self.max_entry_bytes = max_entry_bytes
        try:
            self.root = Path(root).resolve()
        except TypeError as exc:
            raise ConfigurationError(
                "Cache root must be a filesystem path.", code="INVALID_CACHE_ROOT"
            ) from exc
        if self.root in {Path("/"), Path.home().resolve(), Path.cwd().resolve()}:
            raise CacheError(
                "Cache root must be a dedicated subdirectory, not a broad filesystem root.",
                code="UNSAFE_CACHE_ROOT",
                context={"root": str(self.root)},
            )
        self.root.mkdir(parents=True, exist_ok=True)
        if not self.root.is_dir():
            raise CacheError("Cache root is not a directory.", code="INVALID_CACHE_ROOT")
        marker = self.root / _ROOT_MARKER
        nonempty = next(self.root.iterdir(), None) is not None
        if nonempty and not marker.is_file():
            raise CacheError(
                "Existing cache root has no DNAKit cache marker.",
                code="UNSAFE_CACHE_ROOT",
                context={"root": str(self.root)},
                hint="Choose a new empty directory dedicated to DNAKit cache data.",
            )
        marker.touch(exist_ok=True)

    def _namespace_dir(self, namespace: str, *, create: bool) -> Path:
        directory = self.root / namespace
        if directory.is_symlink():
            raise CacheError(
                "Cache namespace must not be a symbolic link.",
                code="UNSAFE_CACHE_NAMESPACE",
                context={"namespace": namespace},
            )
        marker = directory / _NAMESPACE_MARKER
        if create:
            directory.mkdir(parents=False, exist_ok=True)
            if any(directory.iterdir()) and not marker.is_file():
                raise CacheError(
                    "Existing namespace has no DNAKit cache marker.",
                    code="UNSAFE_CACHE_NAMESPACE",
                    context={"namespace": namespace},
                )
            marker.touch(exist_ok=True)
        elif directory.exists() and not marker.is_file():
            raise CacheError(
                "Cache namespace marker is missing.",
                code="UNSAFE_CACHE_NAMESPACE",
                context={"namespace": namespace},
            )
        return directory

    def _path(self, key: CacheKey) -> Path:
        if not isinstance(key, CacheKey):
            raise ConfigurationError("key must be CacheKey.", code="INVALID_CACHE_KEY")
        return self._namespace_dir(key.namespace, create=False) / f"{key.digest}.json"

    def get(self, key: CacheKey) -> object | None:
        path = self._path(key)
        if path.is_symlink():
            raise CacheError(
                "Cache entry must not be a symbolic link.",
                code="UNSAFE_CACHE_ENTRY",
                context={"path": str(path)},
            )
        if not path.exists():
            return None
        try:
            byte_size = path.stat().st_size
            if byte_size > self.max_entry_bytes:
                raise CacheError(
                    "Cache entry exceeds max_entry_bytes.",
                    code="CACHE_ENTRY_SIZE_LIMIT",
                    context={
                        "path": str(path),
                        "byte_size": byte_size,
                        "max_entry_bytes": self.max_entry_bytes,
                    },
                )
            with path.open("rb") as handle:
                raw = handle.read(self.max_entry_bytes + 1)
            if len(raw) > self.max_entry_bytes:
                raise CacheError(
                    "Cache entry grew beyond max_entry_bytes while being read.",
                    code="CACHE_ENTRY_SIZE_LIMIT",
                    context={"path": str(path), "max_entry_bytes": self.max_entry_bytes},
                )
            envelope = json.loads(raw, object_pairs_hook=_unique_json_object)
            _validate_json_shape(envelope)
            if not isinstance(envelope, dict):
                raise ValueError("not an object")
            if (
                envelope.get("cache_schema") != "dnakit-json-cache-v1"
                or envelope.get("key_schema") != key.schema_version
                or envelope.get("namespace") != key.namespace
                or envelope.get("digest") != key.digest
            ):
                raise CacheError(
                    "Cached envelope does not match the requested key.",
                    code="CACHE_KEY_MISMATCH",
                    context={"path": str(path)},
                )
            payload = envelope["payload"]
            expected = envelope["payload_sha256"]
            actual = hashlib.sha256(_canonical_bytes(payload)).hexdigest()
            if actual != expected:
                raise CacheError(
                    "Cached payload checksum does not match.",
                    code="CACHE_INTEGRITY_ERROR",
                    context={"path": str(path), "expected": expected, "actual": actual},
                )
            return to_json_compatible(freeze_json(payload))
        except CacheError:
            raise
        except (
            OSError,
            KeyError,
            TypeError,
            ValueError,
            ConfigurationError,
            json.JSONDecodeError,
            RecursionError,
        ) as exc:
            raise CacheError(
                "Cache entry cannot be decoded.",
                code="CACHE_DECODE_ERROR",
                context={"path": str(path)},
            ) from exc

    def put(self, key: CacheKey, value: object) -> ArtifactRef:
        path = self._path(key)
        path = self._namespace_dir(key.namespace, create=True) / path.name
        _validate_json_shape(value)
        try:
            payload = cast(Any, to_json_compatible(freeze_json(value)))
        except RecursionError as exc:
            raise ConfigurationError(
                "Cache value exceeds structural limits.", code="CACHE_STRUCTURE_LIMIT"
            ) from exc
        payload_sha256 = hashlib.sha256(_canonical_bytes(payload)).hexdigest()
        envelope = {
            "cache_schema": "dnakit-json-cache-v1",
            "key_schema": key.schema_version,
            "namespace": key.namespace,
            "digest": key.digest,
            "payload_sha256": payload_sha256,
            "payload": payload,
        }
        encoded = _canonical_bytes(envelope)
        if len(encoded) > self.max_entry_bytes:
            raise CacheError(
                "Encoded cache entry exceeds max_entry_bytes.",
                code="CACHE_ENTRY_SIZE_LIMIT",
                context={
                    "byte_size": len(encoded),
                    "max_entry_bytes": self.max_entry_bytes,
                },
            )
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{key.digest}.", suffix=".tmp", dir=path.parent
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(encoded)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, path)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
        stat = path.stat()
        return ArtifactRef(
            relative_path=os.path.relpath(path, Path.cwd()),
            media_type="application/vnd.dnakit.cache+json",
            schema_version="dnakit-json-cache-v1",
            sha256=hashlib.sha256(encoded).hexdigest(),
            byte_size=stat.st_size,
            created_at=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
        )

    def invalidate(self, key: CacheKey) -> bool:
        path = self._path(key)
        if not path.exists():
            return False
        try:
            path.unlink()
        except OSError as exc:
            raise CacheError(
                "Cache entry could not be removed.",
                code="CACHE_REMOVE_ERROR",
                context={"path": str(path)},
            ) from exc
        with suppress(OSError):
            path.parent.rmdir()
        return True

    def clear(self, namespace: str | None = None) -> CacheClearReport:
        if namespace is not None and (
            not isinstance(namespace, str) or not _NAMESPACE.fullmatch(namespace)
        ):
            raise ConfigurationError(
                "Cache namespace must be a safe identifier.", code="INVALID_CACHE_NAMESPACE"
            )
        roots = (
            [self._namespace_dir(namespace, create=False)]
            if namespace is not None
            else sorted(self.root.iterdir())
        )
        count = 0
        byte_count = 0
        for directory in roots:
            if (
                not directory.is_dir()
                or directory.is_symlink()
                or directory.parent != self.root
                or not (directory / _NAMESPACE_MARKER).is_file()
            ):
                continue
            for path in sorted(directory.glob("*.json")):
                if path.parent != directory:
                    continue
                byte_count += path.stat().st_size
                path.unlink()
                count += 1
            (directory / _NAMESPACE_MARKER).unlink(missing_ok=True)
            with suppress(OSError):
                directory.rmdir()
        return CacheClearReport(namespace, count, byte_count)


__all__ = ["DEFAULT_MAX_ENTRY_BYTES", "CacheClearReport", "CacheKey", "JSONCache"]
