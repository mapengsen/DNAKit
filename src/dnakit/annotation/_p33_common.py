"""Shared argument checks and result helpers for the P33 pybedtools and pyranges1 modules.

Private to the package: the public functions live in ``pybedtools_p33_functions`` and
``pyranges1_p33_functions``. Nothing here imports pybedtools, pyranges1 or pandas at import time.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError

DEFAULT_ROWS = 1000
MAX_ROWS = 1_000_000
_NAME = re.compile(r"[A-Za-z0-9_.:+-]{1,64}")


def package_version(distribution: str) -> str | None:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return None


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return value


def number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return float(value)


def label(value: object, name: str) -> str:
    if not isinstance(value, str) or not _NAME.fullmatch(value):
        raise ConfigurationError(f"{name} must be 1-64 letters, digits or _.:+-")
    return value


def text(value: object, name: str, maximum: int = 4096) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if "\x00" in value or "\n" in value or "\r" in value:
        raise ConfigurationError(f"{name} must not contain control characters.")
    return value


def rows_limit(value: object) -> int:
    return integer(value, "max_rows", 1, MAX_ROWS)


def new_output(value: object, name: str = "output_path") -> Path | None:
    """Output path that must not exist yet; its parent directory must exist."""

    if value is None:
        return None
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}")
    return path


def plain(value: Any) -> Any:
    """Convert numpy and pandas scalars to JSON-friendly values; NaN and infinity become None."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item") and not isinstance(value, (str, Mapping)):
        try:
            value = value.item()
        except (ValueError, TypeError):
            return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def check_genome_file(path: Path) -> dict[str, int]:
    """Read a two-column chromosome-size file (name, length) and return it as a dict."""

    sizes: dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for number_, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ConfigurationError(f"{path.name} line {number_} must be name and length.")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ConfigurationError(f"{path.name} line {number_}: length is not an integer.") from exc
            if length < 1:
                raise ConfigurationError(f"{path.name} line {number_}: length must be positive.")
            sizes[fields[0]] = length
    if not sizes:
        raise ConfigurationError(f"{path.name} holds no chromosome sizes.")
    return sizes
