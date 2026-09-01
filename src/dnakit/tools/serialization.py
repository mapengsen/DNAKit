"""Bounded conversion of DNAKit results into MCP-safe JSON values."""

from __future__ import annotations

import base64
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from pathlib import PurePath

from dnakit.core import DNA
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class AgentSerializationLimits:
    """Hard limits applied before a result is returned to an Agent."""

    max_depth: int = 64
    max_nodes: int = 200_000
    max_items: int = 20_000
    max_bytes: int = 5_000_000

    def __post_init__(self) -> None:
        values = {
            "max_depth": self.max_depth,
            "max_nodes": self.max_nodes,
            "max_items": self.max_items,
            "max_bytes": self.max_bytes,
        }
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in values.values()
        ):
            raise ConfigurationError(
                "Agent serialization limits must be positive integers.",
                code="AGENT_SERIALIZATION_LIMIT_INVALID",
                context=values,
            )


class _Budget:
    def __init__(self, limits: AgentSerializationLimits) -> None:
        self.limits = limits
        self.nodes = 0
        self.active: set[int] = set()

    def count(self, depth: int) -> None:
        self.nodes += 1
        if depth > self.limits.max_depth or self.nodes > self.limits.max_nodes:
            raise ConfigurationError(
                "Agent result exceeds structural serialization limits.",
                code="AGENT_RESULT_STRUCTURE_LIMIT",
                context={
                    "max_depth": self.limits.max_depth,
                    "max_nodes": self.limits.max_nodes,
                },
            )

    def enter(self, value: object) -> None:
        identity = id(value)
        if identity in self.active:
            raise ConfigurationError(
                "Agent result contains a recursive reference.",
                code="AGENT_RESULT_RECURSIVE",
            )
        self.active.add(identity)

    def leave(self, value: object) -> None:
        self.active.remove(id(value))


def _sequence_values(
    value: Iterable[object],
    *,
    budget: _Budget,
    depth: int,
) -> list[object]:
    items: list[object] = []
    iterator = iter(value)
    try:
        for index, item in enumerate(iterator):
            if index >= budget.limits.max_items:
                raise ConfigurationError(
                    "Agent result contains too many collection items.",
                    code="AGENT_RESULT_ITEM_LIMIT",
                    context={"max_items": budget.limits.max_items},
                    hint="Use a bounded DNAKit input or write the complete result to a file.",
                )
            items.append(_convert(item, budget=budget, depth=depth + 1))
    finally:
        close = getattr(value, "close", None)
        if callable(close):
            close()
    return items


def _convert(value: object, *, budget: _Budget, depth: int) -> object:
    budget.count(depth)
    if isinstance(value, Enum):
        return _convert(value.value, budget=budget, depth=depth + 1)
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ConfigurationError(
                "Agent results cannot contain non-finite floating-point values.",
                code="AGENT_RESULT_NONFINITE",
                context={"value": repr(value)},
            )
        return value
    if isinstance(value, PurePath):
        return str(value)
    if isinstance(value, bytes):
        return {
            "encoding": "base64",
            "data": base64.b64encode(value).decode("ascii"),
        }
    if isinstance(value, DNA):
        return {
            "records": _convert(value.records, budget=budget, depth=depth + 1),
            "name": value.name,
            "source": value.source,
            "version": value.version,
            "metadata": _convert(
                value.collection_metadata,
                budget=budget,
                depth=depth + 1,
            ),
        }
    if isinstance(value, Mapping):
        budget.enter(value)
        try:
            if len(value) > budget.limits.max_items:
                raise ConfigurationError(
                    "Agent result mapping contains too many entries.",
                    code="AGENT_RESULT_ITEM_LIMIT",
                    context={"max_items": budget.limits.max_items},
                )
            return {
                str(key): _convert(item, budget=budget, depth=depth + 1)
                for key, item in value.items()
            }
        finally:
            budget.leave(value)
    if isinstance(value, (list, tuple)):
        budget.enter(value)
        try:
            if len(value) > budget.limits.max_items:
                raise ConfigurationError(
                    "Agent result sequence contains too many items.",
                    code="AGENT_RESULT_ITEM_LIMIT",
                    context={"max_items": budget.limits.max_items},
                )
            return [_convert(item, budget=budget, depth=depth + 1) for item in value]
        finally:
            budget.leave(value)
    if isinstance(value, (set, frozenset)):
        budget.enter(value)
        try:
            converted = [_convert(item, budget=budget, depth=depth + 1) for item in value]
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

    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        budget.enter(value)
        try:
            return _convert(to_dict(), budget=budget, depth=depth + 1)
        finally:
            budget.leave(value)
    if is_dataclass(value) and not isinstance(value, type):
        budget.enter(value)
        try:
            return {
                field.name: _convert(
                    getattr(value, field.name),
                    budget=budget,
                    depth=depth + 1,
                )
                for field in fields(value)
            }
        finally:
            budget.leave(value)

    module_name = type(value).__module__
    if module_name.startswith(("numpy", "torch")):
        to_list = getattr(value, "tolist", None)
        if callable(to_list):
            return _convert(to_list(), budget=budget, depth=depth + 1)
        item = getattr(value, "item", None)
        if callable(item):
            return _convert(item(), budget=budget, depth=depth + 1)
    if isinstance(value, Iterable):
        budget.enter(value)
        try:
            return _sequence_values(value, budget=budget, depth=depth)
        finally:
            budget.leave(value)
    raise ConfigurationError(
        "DNAKit result cannot be represented as bounded JSON.",
        code="AGENT_RESULT_UNSERIALIZABLE",
        context={"result_type": f"{type(value).__module__}.{type(value).__qualname__}"},
        hint="Use the Python API for results containing live backend or callback objects.",
    )


def serialize_for_agent(
    value: object,
    *,
    limits: AgentSerializationLimits | None = None,
) -> object:
    """Return a bounded plain JSON value suitable for an MCP tool response."""

    resolved_limits = AgentSerializationLimits() if limits is None else limits
    if not isinstance(resolved_limits, AgentSerializationLimits):
        raise ConfigurationError(
            "limits must be AgentSerializationLimits or None.",
            code="AGENT_SERIALIZATION_LIMIT_INVALID",
        )
    converted = _convert(value, budget=_Budget(resolved_limits), depth=0)
    encoded = json.dumps(
        converted,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    if len(encoded) > resolved_limits.max_bytes:
        raise ConfigurationError(
            "Agent result exceeds the serialized byte limit.",
            code="AGENT_RESULT_BYTE_LIMIT",
            context={
                "result_bytes": len(encoded),
                "max_bytes": resolved_limits.max_bytes,
            },
            hint="Reduce the requested result or use a DNAKit file-output function.",
        )
    return converted


__all__ = ["AgentSerializationLimits", "serialize_for_agent"]
