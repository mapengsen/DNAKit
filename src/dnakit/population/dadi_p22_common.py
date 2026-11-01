"""Shared input checks and result helpers for the dadi P22 modules (private, no public functions).

Spectra travel as nested lists; an entry that dadi masks is ``None``. dadi and NumPy are imported
only when a function runs.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _integer, _provenance
from dnakit.population.dadi_functions import _nested, _new_file

OPERATION_ERRORS = (ValueError, TypeError, IndexError, KeyError, AttributeError, ZeroDivisionError)
MAX_ARRAY_ENTRIES = 2_000_000


def dadi_and_numpy() -> tuple[Any, Any]:
    return require("dadi", _DADI_EXTRA), require("numpy", _DADI_EXTRA)


def number(
    value: object,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    positive: bool = False,
) -> float:
    """A finite real number of any sign unless bounded."""

    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    result = float(value)
    if positive and result <= 0:
        raise ConfigurationError(f"{name} must be greater than zero.")
    if minimum is not None and result < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum:g}.")
    if maximum is not None and result > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum:g}.")
    return result


def integer(value: object, name: str, *, minimum: int = 0, maximum: int | None = None) -> int:
    return _integer(value, name, minimum=minimum, maximum=maximum)


def reals(
    value: object, name: str, *, length: int | None = None, minimum_length: int = 1
) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of numbers.")
    if len(value) < minimum_length or (length is not None and len(value) != length):
        need = f"exactly {length}" if length is not None else f"at least {minimum_length}"
        raise ConfigurationError(f"{name} must hold {need} numbers; got {len(value)}.")
    return [number(item, f"{name} entry") for item in value]


def choice(value: object, options: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {sorted(options)}.")
    return value


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _fill(value: object, name: str) -> Any:
    if value is None:
        return math.nan
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_fill(item, name) for item in value]
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(
            f"{name} must be a nested list of finite numbers; masked entries are None."
        )
    return float(value)


def masked_array(
    np: Any,
    value: object,
    name: str,
    ndims: Sequence[int] = (1, 2, 3),
    *,
    nonnegative: bool = True,
    minimum_length: int = 3,
    square: bool = False,
) -> Any:
    """A nested list as a masked array, ``None`` marking masked entries."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(
            f"{name} must be a non-empty nested list; masked entries are None."
        )
    try:
        raw = np.array(_fill(value, name), dtype=float)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be a rectangular nested list.") from exc
    if raw.ndim not in ndims:
        raise ConfigurationError(f"{name} must have {list(ndims)} axes, found {raw.ndim}.")
    if raw.size > MAX_ARRAY_ENTRIES:
        raise ConfigurationError(f"{name} is too large ({raw.size} entries).")
    if min(raw.shape) < minimum_length:
        raise ConfigurationError(
            f"{name} needs every axis to have at least {minimum_length} entries."
        )
    if square and len(set(raw.shape)) != 1:
        raise ConfigurationError(f"{name} must have the same length on every axis.")
    present = raw[~np.isnan(raw)]
    if nonnegative and (present < 0).any():
        raise ConfigurationError(f"{name} must contain only non-negative values or None.")
    missing = np.isnan(raw)
    return np.ma.array(np.where(missing, 0.0, raw), mask=missing)


def nested(np: Any, value: Any) -> Any:
    return _nested(np, value)


def scalar(value: Any) -> float | None:
    result = float(value)
    return result if math.isfinite(result) else None


def new_file(destination: object, name: str) -> Path:
    return _new_file(destination, name)


def existing_file(path: object, name: str) -> Path:
    if not isinstance(path, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(str(path)).expanduser().resolve()
    if not target.is_file():
        raise ConfigurationError(f"{name} does not exist: {target}")
    return target


def result(operation: str, convention: str, data: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        **data,
        "provenance": _provenance("dadi", convention),
    }


def guarded(operation: str, call: Any) -> Any:
    """Run a dadi call and turn its argument errors into ``ConfigurationError``."""

    try:
        return call()
    except OPERATION_ERRORS as exc:
        raise ConfigurationError(f"dadi could not run {operation}: {exc}") from exc
