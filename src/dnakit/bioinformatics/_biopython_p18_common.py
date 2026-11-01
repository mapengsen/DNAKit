"""Private helpers shared by the four Biopython pack-18 modules.

The modules ``dnakit.sequence.biopython_p18_functions``,
``dnakit.alignment.biopython_p18_alignment_functions``,
``dnakit.phylogeny.biopython_p18_tree_functions`` and
``dnakit.regulation.biopython_p18_motif_data_functions`` use these small validators so that they
raise the same ``ConfigurationError`` messages and report the same provenance block.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path
from dnakit.exceptions import ConfigurationError

EXTRA = "bio-formats"
PACKAGE = "biopython"
MAX_SEQUENCE = 10_000_000
MAX_RECORDS = 100_000
MAX_FILE_BYTES = 512 * 1024 * 1024
MAX_TEXT = 5_000_000
_LONG_LIST = 1_000


def provenance(**notes: str) -> dict[str, Any]:
    """Return the package name, the installed version and the given notes."""

    try:
        installed: str | None = importlib.metadata.version(PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": PACKAGE, "package_version": installed, **notes}


def integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def optional_integer(value: object, name: str, *, minimum: int, maximum: int) -> int | None:
    return None if value is None else integer(value, name, minimum=minimum, maximum=maximum)


def number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    result = float(value)
    if result != result or result in (float("inf"), float("-inf")):
        raise ConfigurationError(f"{name} must be a finite number.")
    return result


def choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def text(value: object, name: str, *, limit: int = MAX_TEXT, blank: bool = False) -> str:
    if not isinstance(value, str) or (not blank and not value.strip()):
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} is longer than {limit} characters.")
    return value


def identifier(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or any(c.isspace() for c in value):
        raise ConfigurationError(f"{name} must be a non-empty string without spaces.")
    return value


def source_file(path: object, name: str = "path") -> Path:
    """Return an existing input file below the size limit."""

    if not isinstance(path, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    source = input_path(path)
    if source.stat().st_size > MAX_FILE_BYTES:
        raise ConfigurationError(f"{source.name} is larger than {MAX_FILE_BYTES} bytes.")
    return source


def new_output(value: object, name: str = "output_path") -> Path:
    """Return a path that does not exist yet and whose directory exists."""

    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new file: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"{name} parent directory does not exist: {path.parent}")
    return path


def open_new(path: Path, *, binary: bool) -> Any:
    """Open a new file for writing (exclusive creation), in binary or UTF-8 text mode."""

    if binary:
        return path.open("xb")
    return path.open("x", encoding="utf-8")


def plain(value: object) -> Any:
    """Turn Biopython containers into JSON-friendly values; long lists are summarised."""

    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return (
            value if value == value and value not in (float("inf"), float("-inf")) else str(value)
        )
    if isinstance(value, bytes):
        return value.decode("latin-1")
    if isinstance(value, Mapping):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [plain(item) for item in sorted(value, key=str)]
    if isinstance(value, (list, tuple)):
        if len(value) > _LONG_LIST:
            return {"length": len(value), "first": [plain(item) for item in value[:20]]}
        return [plain(item) for item in value]
    if hasattr(value, "tolist"):
        return plain(value.tolist())
    return str(value)


def optional_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return result if result == result and result not in (float("inf"), float("-inf")) else None
