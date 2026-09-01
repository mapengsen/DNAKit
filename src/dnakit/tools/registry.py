"""Search, describe, and execute the Agent-facing DNAKit tool catalog."""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from collections.abc import Mapping
from functools import lru_cache

from dnakit.exceptions import ConfigurationError

from .adapters import coerce_arguments
from .catalog import CATEGORY_LABELS, AgentToolSpec, build_tool_catalog
from .serialization import AgentSerializationLimits, serialize_for_agent

_TOKEN = re.compile(r"[\w]+", re.UNICODE)
_MAX_ARGUMENT_BYTES = 2_000_000


class AgentToolRegistry:
    """Immutable searchable registry over public DNAKit functions."""

    def __init__(
        self,
        tools: tuple[AgentToolSpec, ...] | None = None,
        *,
        serialization_limits: AgentSerializationLimits | None = None,
    ) -> None:
        discovered = build_tool_catalog() if tools is None else tools
        names = [tool.name for tool in discovered]
        if len(names) != len(set(names)):
            duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
            raise ConfigurationError(
                "Agent tool names must be unique.",
                code="AGENT_TOOL_DUPLICATE",
                context={"tools": duplicates},
            )
        self._tools = tuple(sorted(discovered, key=lambda tool: tool.name))
        self._by_name = {tool.name: tool for tool in self._tools}
        aliases: defaultdict[str, list[AgentToolSpec]] = defaultdict(list)
        for tool in self._tools:
            aliases[tool.name.removeprefix("dnakit.")].append(tool)
            aliases[tool.name.rsplit(".", 1)[-1]].append(tool)
        self._aliases = {name: tuple(values) for name, values in aliases.items()}
        if serialization_limits is not None and not isinstance(
            serialization_limits, AgentSerializationLimits
        ):
            raise ConfigurationError(
                "serialization_limits must be AgentSerializationLimits or None.",
                code="AGENT_SERIALIZATION_LIMIT_INVALID",
            )
        self._serialization_limits = (
            AgentSerializationLimits() if serialization_limits is None else serialization_limits
        )

    @property
    def tools(self) -> tuple[AgentToolSpec, ...]:
        return self._tools

    def stats(self) -> dict[str, object]:
        """Return aggregate catalog counts without loading optional backends."""

        compatible = sum(tool.agent_compatible for tool in self._tools)
        return {
            "total": len(self._tools),
            "agent_compatible": compatible,
            "python_only": len(self._tools) - compatible,
            "categories": len({tool.category for tool in self._tools}),
        }

    def categories(self) -> list[dict[str, object]]:
        """Return stable category names and callable counts."""

        totals = Counter(tool.category for tool in self._tools)
        compatible = Counter(tool.category for tool in self._tools if tool.agent_compatible)
        result: list[dict[str, object]] = []
        for category in sorted(totals):
            labels = CATEGORY_LABELS.get(category, (category, category))
            result.append(
                {
                    "category": category,
                    "label_zh": labels[0],
                    "label_en": labels[1],
                    "total": totals[category],
                    "agent_compatible": compatible[category],
                }
            )
        return result

    def resolve(self, name: str) -> AgentToolSpec:
        """Resolve a full name or an unambiguous category/function alias."""

        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError(
                "Agent tool name must be a non-empty string.",
                code="AGENT_TOOL_NAME_INVALID",
            )
        normalized = name.strip()
        exact = self._by_name.get(normalized)
        if exact is not None:
            return exact
        matches = self._aliases.get(normalized, ())
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ConfigurationError(
                "Agent tool alias is ambiguous.",
                code="AGENT_TOOL_AMBIGUOUS",
                context={"alias": normalized, "matches": [tool.name for tool in matches]},
                hint="Use the complete dnakit.<category>.<function> name.",
            )
        suggestions = [
            tool.name for tool in self.search(normalized, limit=5, include_python_only=True)
        ]
        raise ConfigurationError(
            "Agent tool was not found.",
            code="AGENT_TOOL_NOT_FOUND",
            context={"tool": normalized, "suggestions": suggestions},
        )

    def describe(self, name: str) -> dict[str, object]:
        """Return one full Agent tool manifest."""

        return self.resolve(name).to_dict()

    def search(
        self,
        query: str = "",
        *,
        category: str | None = None,
        limit: int = 20,
        include_python_only: bool = False,
    ) -> tuple[AgentToolSpec, ...]:
        """Search names, descriptions, and bilingual category labels."""

        if not isinstance(query, str):
            raise ConfigurationError(
                "Agent tool search query must be a string.",
                code="AGENT_TOOL_QUERY_INVALID",
            )
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
            raise ConfigurationError(
                "Agent tool search limit must be an integer in [1, 100].",
                code="AGENT_TOOL_LIMIT_INVALID",
            )
        if category is not None and (
            not isinstance(category, str) or category not in CATEGORY_LABELS
        ):
            raise ConfigurationError(
                "Unknown Agent tool category.",
                code="AGENT_TOOL_CATEGORY_INVALID",
                context={"category": category},
            )
        normalized = query.strip().casefold()
        query_tokens = set(_TOKEN.findall(normalized))
        ranked: list[tuple[int, str, AgentToolSpec]] = []
        for tool in self._tools:
            if category is not None and tool.category != category:
                continue
            if not include_python_only and not tool.agent_compatible:
                continue
            labels = CATEGORY_LABELS.get(tool.category, (tool.category, tool.category))
            short_name = tool.name.rsplit(".", 1)[-1]
            haystack = " ".join(
                (
                    tool.name,
                    short_name.replace("_", " "),
                    tool.description,
                    labels[0],
                    labels[1],
                )
            ).casefold()
            if not normalized:
                score = 1
            else:
                haystack_tokens = set(_TOKEN.findall(haystack))
                score = 0
                if normalized == tool.name.casefold():
                    score += 1_000
                if normalized == short_name.casefold():
                    score += 800
                if normalized in haystack:
                    score += 300
                score += 40 * len(query_tokens & haystack_tokens)
                score += sum(10 for token in query_tokens if token in haystack)
            if score > 0:
                ranked.append((score, tool.name, tool))
        ranked.sort(key=lambda item: (-item[0], item[1]))
        return tuple(item[2] for item in ranked[:limit])

    def execute(
        self,
        name: str,
        arguments: Mapping[str, object] | None = None,
        *,
        allow_side_effects: bool = False,
    ) -> object:
        """Execute one compatible public function and return bounded JSON."""

        tool = self.resolve(name)
        if not tool.agent_compatible:
            raise ConfigurationError(
                "This DNAKit function requires live Python objects and is not Agent-callable.",
                code="AGENT_TOOL_PYTHON_ONLY",
                context={
                    "tool": tool.name,
                    "incompatibilities": list(tool.incompatibilities),
                },
                hint="Call this function through the regular Python API.",
            )
        if tool.requires_confirmation and not allow_side_effects:
            raise ConfigurationError(
                "This DNAKit tool may write files, download models, or run an external program.",
                code="AGENT_TOOL_CONFIRMATION_REQUIRED",
                context={"tool": tool.name, "effect": tool.effect},
                hint="Review the tool and set allow_side_effects=true to authorize this call.",
            )
        resolved_arguments: Mapping[str, object] = {} if arguments is None else arguments
        if not isinstance(resolved_arguments, Mapping) or any(
            not isinstance(key, str) for key in resolved_arguments
        ):
            raise ConfigurationError(
                "Agent tool arguments must be an object with string keys.",
                code="AGENT_ARGUMENTS_INVALID",
            )
        try:
            encoded_arguments = json.dumps(
                resolved_arguments,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Agent tool arguments must contain only JSON values.",
                code="AGENT_ARGUMENTS_INVALID",
            ) from exc
        if len(encoded_arguments) > _MAX_ARGUMENT_BYTES:
            raise ConfigurationError(
                "Agent tool arguments exceed the byte limit.",
                code="AGENT_ARGUMENT_BYTE_LIMIT",
                context={
                    "argument_bytes": len(encoded_arguments),
                    "max_bytes": _MAX_ARGUMENT_BYTES,
                },
            )
        positional, keywords = coerce_arguments(tool.call_plan, resolved_arguments)
        result = tool.function(*positional, **keywords)
        return serialize_for_agent(result, limits=self._serialization_limits)


@lru_cache(maxsize=1)
def default_tool_registry() -> AgentToolRegistry:
    """Return the process-wide immutable DNAKit Agent registry."""

    return AgentToolRegistry()


__all__ = ["AgentToolRegistry", "default_tool_registry"]
