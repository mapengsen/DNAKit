"""Thread-safe lazy registry for optional computational backends."""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from importlib.metadata import EntryPoint, entry_points
from threading import RLock

from dnakit.core import BackendInfo
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
)

BackendProbe = Callable[[], BackendInfo]
BackendLoader = Callable[[], object]

_BACKEND_ID = re.compile(r"^[a-z][a-z0-9_.-]*$")


def _validate_backend_id(backend_id: object) -> str:
    if not isinstance(backend_id, str) or not _BACKEND_ID.fullmatch(backend_id):
        raise ConfigurationError(
            "backend_id must be a lowercase dotted, dashed, or underscored identifier.",
            code="INVALID_BACKEND_ID",
            context={"backend_id": backend_id},
        )
    return backend_id


@dataclass(frozen=True, slots=True)
class BackendRegistration:
    """Lazy probe and loader associated with one stable backend identifier."""

    backend_id: str
    probe: BackendProbe
    loader: BackendLoader | None = None
    source: str = "manual"

    def __post_init__(self) -> None:
        _validate_backend_id(self.backend_id)
        if not callable(self.probe):
            raise ConfigurationError("probe must be callable.", code="INVALID_BACKEND_PROBE")
        if self.loader is not None and not callable(self.loader):
            raise ConfigurationError(
                "loader must be callable or None.", code="INVALID_BACKEND_LOADER"
            )
        if not isinstance(self.source, str) or not self.source.strip():
            raise ConfigurationError(
                "source must be a non-empty string.", code="INVALID_BACKEND_SOURCE"
            )


def _entry_point_registration(point: EntryPoint) -> BackendRegistration:
    def load() -> object:
        candidate = point.load()
        return candidate() if isinstance(candidate, type) else candidate

    def probe() -> BackendInfo:
        candidate = load()
        probe_method = getattr(candidate, "probe", None)
        if not callable(probe_method):
            raise BackendExecutionError(
                "Backend entry point does not expose a callable probe().",
                code="BACKEND_PROBE_MISSING",
                context={"backend_id": point.name, "entry_point": point.value},
            )
        info = probe_method()
        if not isinstance(info, BackendInfo):
            raise BackendExecutionError(
                "Backend probe did not return BackendInfo.",
                code="INVALID_BACKEND_PROBE_RESULT",
                context={"backend_id": point.name},
            )
        return info

    return BackendRegistration(point.name, probe, load, source=f"entry-point:{point.value}")


class BackendRegistry:
    """Register and resolve backend adapters without importing them eagerly."""

    def __init__(self) -> None:
        self._registrations: dict[str, BackendRegistration] = {}
        self._lock = RLock()

    def register(self, registration: BackendRegistration, *, replace: bool = False) -> None:
        if not isinstance(registration, BackendRegistration):
            raise ConfigurationError(
                "registration must be BackendRegistration.", code="INVALID_BACKEND_REGISTRATION"
            )
        if not isinstance(replace, bool):
            raise ConfigurationError("replace must be boolean.", code="INVALID_BACKEND_REPLACE")
        with self._lock:
            if registration.backend_id in self._registrations and not replace:
                raise ConfigurationError(
                    "Backend identifier is already registered.",
                    code="DUPLICATE_BACKEND",
                    context={"backend_id": registration.backend_id},
                )
            self._registrations[registration.backend_id] = registration

    def unregister(self, backend_id: str) -> bool:
        resolved_id = _validate_backend_id(backend_id)
        with self._lock:
            return self._registrations.pop(resolved_id, None) is not None

    def registrations(self) -> tuple[BackendRegistration, ...]:
        with self._lock:
            return tuple(self._registrations[key] for key in sorted(self._registrations))

    def __iter__(self) -> Iterator[str]:
        return iter(tuple(item.backend_id for item in self.registrations()))

    def __len__(self) -> int:
        with self._lock:
            return len(self._registrations)

    def _get(self, backend_id: str) -> BackendRegistration:
        resolved_id = _validate_backend_id(backend_id)
        with self._lock:
            registration = self._registrations.get(resolved_id)
        if registration is None:
            raise BackendUnavailableError(
                "Requested backend is not registered.",
                code="BACKEND_NOT_REGISTERED",
                context={"backend_id": backend_id},
                hint="Inspect dnakit backends and install or register the required adapter.",
            )
        return registration

    def probe(self, backend_id: str) -> BackendInfo:
        registration = self._get(backend_id)
        try:
            result = registration.probe()
        except BackendUnavailableError:
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "Backend probe failed.",
                code="BACKEND_PROBE_FAILED",
                context={"backend_id": backend_id, "error_type": type(exc).__name__},
            ) from exc
        if not isinstance(result, BackendInfo):
            raise BackendExecutionError(
                "Backend probe did not return BackendInfo.",
                code="INVALID_BACKEND_PROBE_RESULT",
                context={"backend_id": backend_id},
            )
        return result

    def probe_all(self) -> tuple[BackendInfo, ...]:
        return tuple(self.probe(backend_id) for backend_id in self)

    def load(self, backend_id: str, *, capability: str | None = None) -> object:
        if capability is not None and (not isinstance(capability, str) or not capability.strip()):
            raise ConfigurationError(
                "capability must be a non-empty string or None.",
                code="INVALID_BACKEND_CAPABILITY",
            )
        registration = self._get(backend_id)
        info = self.probe(backend_id)
        if not info.available:
            raise BackendUnavailableError(
                "Registered backend is not available in this environment.",
                context={"backend_id": backend_id},
            )
        if capability is not None and capability not in info.capabilities:
            raise BackendUnavailableError(
                "Backend does not advertise the requested capability.",
                code="BACKEND_CAPABILITY_UNAVAILABLE",
                context={"backend_id": backend_id, "capability": capability},
            )
        if registration.loader is None:
            raise BackendUnavailableError(
                "Backend registration has no adapter loader.",
                code="BACKEND_LOADER_UNAVAILABLE",
                context={"backend_id": backend_id},
            )
        try:
            return registration.loader()
        except Exception as exc:
            raise BackendExecutionError(
                "Backend adapter could not be loaded.",
                code="BACKEND_LOAD_FAILED",
                context={"backend_id": backend_id, "error_type": type(exc).__name__},
            ) from exc

    def discover(self, *, group: str = "dnakit.backends", replace: bool = False) -> tuple[str, ...]:
        if not isinstance(group, str) or not group.strip():
            raise ConfigurationError("group must be non-empty.", code="INVALID_ENTRY_POINT_GROUP")
        if not isinstance(replace, bool):
            raise ConfigurationError("replace must be boolean.", code="INVALID_BACKEND_REPLACE")
        discovered: list[str] = []
        selected = entry_points().select(group=group)
        for point in sorted(selected, key=lambda item: item.name):
            self.register(_entry_point_registration(point), replace=replace)
            discovered.append(point.name)
        return tuple(discovered)


backend_registry = BackendRegistry()

__all__ = ["BackendRegistration", "BackendRegistry", "backend_registry"]
