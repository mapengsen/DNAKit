"""Type-preserving keys for immutable JSON metadata values."""

from __future__ import annotations

from collections.abc import Hashable

from dnakit.core._json import FrozenDict, JSONValue


def metadata_value_key(value: JSONValue) -> Hashable:
    """Return a recursive key that does not apply Python numeric coercion rules."""

    if value is None:
        return ("null",)
    if isinstance(value, bool):
        return ("bool", value)
    if isinstance(value, int):
        return ("int", value)
    if isinstance(value, float):
        return ("float", value)
    if isinstance(value, str):
        return ("string", value)
    if isinstance(value, tuple):
        return ("array", tuple(metadata_value_key(item) for item in value))
    if isinstance(value, FrozenDict):
        return (
            "object",
            tuple((key, metadata_value_key(item)) for key, item in sorted(value.items())),
        )
    raise TypeError(f"Unsupported frozen JSON metadata type: {type(value).__name__}")


__all__ = ["metadata_value_key"]
