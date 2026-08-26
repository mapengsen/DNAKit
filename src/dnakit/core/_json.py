"""Immutable JSON-compatible containers for core value objects."""

from __future__ import annotations

import json
import math
from collections.abc import Iterator, Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
from types import MappingProxyType
from typing import TypeAlias, Union

from dnakit.exceptions import ConfigurationError

JSONScalar: TypeAlias = bool | int | float | str | None
JSONValue: TypeAlias = Union[JSONScalar, tuple["JSONValue", ...], "FrozenDict"]
_MAX_JSON_DEPTH = 64
_MAX_JSON_NODES = 1_000_000


class _TraversalBudget:
    __slots__ = ("active", "nodes")

    def __init__(self) -> None:
        self.nodes = 0
        self.active: set[int] = set()

    def count(self, depth: int) -> None:
        self.nodes += 1
        if depth > _MAX_JSON_DEPTH or self.nodes > _MAX_JSON_NODES:
            raise ConfigurationError(
                "JSON-compatible value exceeds structural limits.",
                code="JSON_STRUCTURE_LIMIT",
                context={"max_depth": _MAX_JSON_DEPTH, "max_nodes": _MAX_JSON_NODES},
            )

    def enter(self, value: object) -> None:
        identity = id(value)
        if identity in self.active:
            raise ConfigurationError(
                "JSON-compatible value contains a recursive reference.",
                code="JSON_RECURSIVE_REFERENCE",
            )
        self.active.add(identity)

    def leave(self, value: object) -> None:
        self.active.remove(id(value))


class FrozenDict(Mapping[str, JSONValue]):
    """A recursively immutable and hashable mapping with string keys."""

    __slots__ = ("_data", "_hash")

    _data: Mapping[str, JSONValue]
    _hash: int

    def __init__(self, values: Mapping[str, object] | None = None) -> None:
        source = {} if values is None else values
        frozen = _freeze_json(source, _TraversalBudget(), 0)
        if not isinstance(frozen, FrozenDict):  # pragma: no cover - source is a Mapping.
            raise ConfigurationError("FrozenDict source must be a mapping.")
        object.__setattr__(self, "_data", frozen._data)
        object.__setattr__(self, "_hash", frozen._hash)

    def __getitem__(self, key: str) -> JSONValue:
        return self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mapping):
            return NotImplemented
        return dict(self.items()) == dict(other.items())

    def __repr__(self) -> str:
        return f"FrozenDict({self._data!r})"

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("FrozenDict is immutable.")


def _freeze_json(value: object, budget: _TraversalBudget, depth: int) -> JSONValue:
    budget.count(depth)
    if isinstance(value, Enum):
        return _freeze_json(value.value, budget, depth + 1)
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ConfigurationError(
                "JSON-compatible floating-point values must be finite.",
                context={"value": repr(value)},
            )
        return value
    if isinstance(value, Mapping):
        budget.enter(value)
        try:
            data: dict[str, JSONValue] = {}
            for key in value:
                budget.count(depth + 1)
                if not isinstance(key, str):
                    raise ConfigurationError(
                        "JSON-compatible mapping keys must be strings.",
                        context={"key_type": type(key).__name__},
                    )
                data[key] = _freeze_json(value[key], budget, depth + 1)
            frozen = object.__new__(FrozenDict)
            object.__setattr__(frozen, "_data", MappingProxyType(data))
            object.__setattr__(frozen, "_hash", hash(tuple(sorted(data.items()))))
            return frozen
        finally:
            budget.leave(value)
    if isinstance(value, (tuple, list)):
        budget.enter(value)
        try:
            return tuple(_freeze_json(item, budget, depth + 1) for item in value)
        finally:
            budget.leave(value)
    raise ConfigurationError(
        "Value is not JSON-compatible.",
        context={"value_type": type(value).__name__},
        hint="Use None, booleans, finite numbers, strings, lists, or string-keyed mappings.",
    )


def freeze_json(value: object) -> JSONValue:
    """Validate and recursively freeze a bounded JSON-compatible value."""

    return _freeze_json(value, _TraversalBudget(), 0)


def freeze_mapping(values: Mapping[str, object] | None) -> FrozenDict:
    """Return an immutable defensive copy of a JSON-compatible mapping."""

    return FrozenDict(values)


def to_json_compatible(value: object) -> object:
    """Convert core values into plain objects accepted by ``json.dumps``."""

    return _to_json_compatible(value, _TraversalBudget(), 0)


def _to_json_compatible(value: object, budget: _TraversalBudget, depth: int) -> object:
    budget.count(depth)
    if isinstance(value, Enum):
        return _to_json_compatible(value.value, budget, depth + 1)
    if isinstance(value, Mapping):
        budget.enter(value)
        try:
            result: dict[str, object] = {}
            for key in value:
                budget.count(depth + 1)
                result[str(key)] = _to_json_compatible(value[key], budget, depth + 1)
            return result
        finally:
            budget.leave(value)
    if isinstance(value, frozenset):
        budget.enter(value)
        try:
            converted = [_to_json_compatible(item, budget, depth + 1) for item in value]
            return sorted(
                converted,
                key=lambda item: json.dumps(
                    item,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            )
        finally:
            budget.leave(value)
    if isinstance(value, (tuple, list)):
        budget.enter(value)
        try:
            return [_to_json_compatible(item, budget, depth + 1) for item in value]
        finally:
            budget.leave(value)
    if is_dataclass(value) and not isinstance(value, type):
        budget.enter(value)
        try:
            return {
                field.name: _to_json_compatible(getattr(value, field.name), budget, depth + 1)
                for field in fields(value)
            }
        finally:
            budget.leave(value)
    return value


__all__ = [
    "FrozenDict",
    "JSONScalar",
    "JSONValue",
    "freeze_json",
    "freeze_mapping",
    "to_json_compatible",
]
