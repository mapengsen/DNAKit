"""Small validation and conversion helpers shared by the pack-38 modules (internal, not public API)."""

from __future__ import annotations

import math
from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path
from dnakit.exceptions import ConfigurationError


def integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    result = float(value)
    if not math.isfinite(result) or not minimum <= result <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return result


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def choice(value: object, name: str, allowed: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def new_path(value: object, name: str) -> Path:
    """A path for a file to be created: the parent must exist and the file must not."""

    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; refusing to overwrite: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory of {name} does not exist: {path.parent}")
    return path


def plain(value: Any) -> Any:
    """Convert numpy/pandas scalars to Python values; NaN and infinity become None."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item") and not isinstance(value, (str, Mapping)):
        try:
            value = value.item()
        except (ValueError, TypeError):
            return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value.__class__.__name__ in {"Interval", "Timestamp"}:
        return str(value)
    return value


def records(frame: Any, limit: int | None = None) -> list[dict[str, Any]]:
    rows = frame if limit is None else frame.head(limit)
    return [{str(k): plain(v) for k, v in row.items()} for row in rows.to_dict(orient="records")]


def provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}
