"""Structured, serializable non-fatal diagnostics."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.coordinates import CompoundLocation, Interval, Location, UnresolvedLocation
from dnakit.core.enums import IssueSeverity
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, init=False)
class Issue:
    """A stable machine-readable finding that does not require raising."""

    code: str
    severity: IssueSeverity
    message: str
    location: Location | None
    details: FrozenDict

    def __init__(
        self,
        code: str,
        severity: IssueSeverity | str,
        message: str,
        location: Location | None = None,
        details: Mapping[str, object] | None = None,
    ) -> None:
        if not isinstance(code, str) or not code.strip():
            raise ConfigurationError("Issue code must be a non-empty string.")
        if not isinstance(message, str) or not message.strip():
            raise ConfigurationError("Issue message must be a non-empty string.")
        try:
            resolved_severity = (
                severity if isinstance(severity, IssueSeverity) else IssueSeverity(severity)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Unknown issue severity.",
                context={"severity": severity},
                hint=f"Choose one of: {', '.join(item.value for item in IssueSeverity)}.",
            ) from exc
        if location is not None and not isinstance(
            location, (Interval, CompoundLocation, UnresolvedLocation)
        ):
            raise ConfigurationError("Issue location must be None or an internal Location.")

        object.__setattr__(self, "code", code)
        object.__setattr__(self, "severity", resolved_severity)
        object.__setattr__(self, "message", message)
        object.__setattr__(self, "location", location)
        object.__setattr__(self, "details", freeze_mapping(details))

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["Issue"]
