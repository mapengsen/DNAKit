"""Passive metadata describing an optional computational backend."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, init=False)
class BackendInfo:
    """Backend probe result without importing a domain protocol or adapter."""

    name: str
    version: str | None
    executable_path: str | None
    package_location: str | None
    license_expression: str | None
    capabilities: frozenset[str]
    available: bool
    metadata: FrozenDict

    def __init__(
        self,
        name: str,
        *,
        version: str | None = None,
        executable_path: str | None = None,
        package_location: str | None = None,
        license_expression: str | None = None,
        capabilities: Iterable[str] = (),
        available: bool = True,
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("Backend name must be a non-empty string.")
        for field_name, value in (
            ("version", version),
            ("executable_path", executable_path),
            ("package_location", package_location),
            ("license_expression", license_expression),
        ):
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(
                    f"Backend {field_name} must be None or a non-empty string.",
                    context={field_name: value},
                )
        capability_set = frozenset(capabilities)
        if any(not isinstance(item, str) or not item.strip() for item in capability_set):
            raise ConfigurationError("Backend capabilities must be non-empty strings.")
        if not isinstance(available, bool):
            raise ConfigurationError("Backend available must be a boolean.")

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "executable_path", executable_path)
        object.__setattr__(self, "package_location", package_location)
        object.__setattr__(self, "license_expression", license_expression)
        object.__setattr__(self, "capabilities", capability_set)
        object.__setattr__(self, "available", available)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["BackendInfo"]
