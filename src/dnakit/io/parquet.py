"""Conditional availability probe for the optional PyArrow Parquet adapter."""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Any, cast

from dnakit.core._json import to_json_compatible
from dnakit.exceptions import BackendUnavailableError


@dataclass(frozen=True, slots=True)
class ParquetBackendStatus:
    available: bool
    distribution: str = "pyarrow"
    version: str | None = None
    implementation: str = "conditional-adapter"

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def parquet_backend_status() -> ParquetBackendStatus:
    """Probe without importing PyArrow or mutating the environment."""

    try:
        available = importlib.util.find_spec("pyarrow") is not None
    except (ImportError, AttributeError, ValueError):
        available = False
    resolved_version: str | None = None
    if available:
        try:
            resolved_version = version("pyarrow")
        except PackageNotFoundError:
            available = False
    return ParquetBackendStatus(available, version=resolved_version)


def require_parquet_backend() -> ParquetBackendStatus:
    """Return backend metadata or raise one stable actionable error."""

    status = parquet_backend_status()
    if not status.available:
        raise BackendUnavailableError(
            "The optional PyArrow Parquet backend is not installed.",
            code="PARQUET_BACKEND_UNAVAILABLE",
            hint=(
                "Install a compatible pyarrow build in the local environment before using Parquet."
            ),
        )
    return status


__all__ = ["ParquetBackendStatus", "parquet_backend_status", "require_parquet_backend"]
