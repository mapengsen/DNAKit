"""Compact MCP server exposing the complete DNAKit Agent tool registry."""

from __future__ import annotations

import importlib
import time
from collections.abc import Callable, Mapping
from typing import Protocol, cast

from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError

from .catalog import CATEGORY_LABELS
from .registry import AgentToolRegistry, default_tool_registry
from .serialization import serialize_for_agent


class MCPServer(Protocol):
    """Small protocol used to keep FastMCP an optional dependency."""

    def tool(
        self,
        *,
        name: str,
        description: str,
        annotations: Mapping[str, object],
    ) -> Callable[[Callable[..., object]], Callable[..., object]]: ...

    def run(self, *, show_banner: bool | None = None) -> None: ...


MCPServerFactory = Callable[..., MCPServer]


def _fastmcp_factory() -> MCPServerFactory:
    try:
        module = importlib.import_module("fastmcp")
    except ModuleNotFoundError as exc:
        raise BackendUnavailableError(
            "FastMCP is required to run the DNAKit Agent server.",
            code="DNAKIT_MCP_UNAVAILABLE",
            hint='Install it with: python -m pip install "dnakit[agent]"',
        ) from exc
    factory = getattr(module, "FastMCP", None)
    if factory is None or not callable(factory):
        raise BackendUnavailableError(
            "The installed fastmcp package does not expose FastMCP.",
            code="DNAKIT_MCP_INCOMPATIBLE",
            hint="Install a DNAKit-compatible fastmcp release.",
        )
    return cast(MCPServerFactory, factory)


def _error_payload(tool_name: str | None, error: Exception) -> dict[str, object]:
    if isinstance(error, DNAKitError):
        context = serialize_for_agent(dict(error.context))
        return {
            "ok": False,
            "tool": tool_name,
            "error": {
                "type": type(error).__name__,
                "code": error.code,
                "message": error.message,
                "context": context,
                "hint": error.hint,
            },
        }
    return {
        "ok": False,
        "tool": tool_name,
        "error": {
            "type": type(error).__name__,
            "code": "AGENT_TOOL_EXECUTION_ERROR",
            "message": str(error),
            "context": {},
            "hint": None,
        },
    }


def create_server(registry: AgentToolRegistry | None = None) -> MCPServer:
    """Create a compact MCP server without importing optional model backends."""

    resolved_registry = default_tool_registry() if registry is None else registry
    if not isinstance(resolved_registry, AgentToolRegistry):
        raise ConfigurationError(
            "registry must be AgentToolRegistry or None.",
            code="AGENT_REGISTRY_INVALID",
        )
    server = _fastmcp_factory()(
        name="DNAKit",
        instructions=(
            "Use search_dnakit_tools to find a capability, describe_dnakit_tool to inspect "
            "its JSON schema, then call_dnakit_tool with the full tool name. File-writing, "
            "model-download, and external-process tools require allow_side_effects=true."
        ),
    )

    @server.tool(
        name="dnakit_catalog_stats",
        description="Return aggregate counts for the installed DNAKit Agent tool catalog.",
        annotations={
            "readOnlyHint": True,
            "destructiveHint": False,
            "idempotentHint": True,
            "openWorldHint": False,
        },
    )
    def catalog_stats() -> dict[str, object]:
        return resolved_registry.stats()

    @server.tool(
        name="list_dnakit_categories",
        description="List DNAKit tool categories with Chinese and English labels.",
        annotations={
            "readOnlyHint": True,
            "destructiveHint": False,
            "idempotentHint": True,
            "openWorldHint": False,
        },
    )
    def list_categories() -> dict[str, object]:
        return {"categories": resolved_registry.categories()}

    @server.tool(
        name="list_dnakit_tools",
        description="List a bounded page of DNAKit functions available to the Agent.",
        annotations={
            "readOnlyHint": True,
            "destructiveHint": False,
            "idempotentHint": True,
            "openWorldHint": False,
        },
    )
    def list_tools(
        category: str | None = None,
        offset: int = 0,
        limit: int = 50,
        include_python_only: bool = False,
    ) -> dict[str, object]:
        try:
            if category is not None and category not in CATEGORY_LABELS:
                raise ConfigurationError(
                    "Unknown DNAKit tool category.",
                    code="AGENT_TOOL_CATEGORY_INVALID",
                    context={"category": category},
                )
            if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
                raise ConfigurationError(
                    "offset must be a non-negative integer.",
                    code="AGENT_TOOL_OFFSET_INVALID",
                )
            if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
                raise ConfigurationError(
                    "limit must be an integer in [1, 100].",
                    code="AGENT_TOOL_LIMIT_INVALID",
                )
            matches = [
                tool
                for tool in resolved_registry.tools
                if (category is None or tool.category == category)
                and (include_python_only or tool.agent_compatible)
            ]
            page = matches[offset : offset + limit]
            return {
                "ok": True,
                "total": len(matches),
                "offset": offset,
                "limit": limit,
                "tools": [tool.to_summary() for tool in page],
            }
        except Exception as exc:
            return _error_payload(None, exc)

    @server.tool(
        name="search_dnakit_tools",
        description=(
            "Search DNAKit public functions by capability, function name, or category. "
            "Use English scientific terms when possible."
        ),
        annotations={
            "readOnlyHint": True,
            "destructiveHint": False,
            "idempotentHint": True,
            "openWorldHint": False,
        },
    )
    def search_tools(
        query: str,
        category: str | None = None,
        limit: int = 20,
        include_python_only: bool = False,
    ) -> dict[str, object]:
        try:
            matches = resolved_registry.search(
                query,
                category=category,
                limit=limit,
                include_python_only=include_python_only,
            )
            return {
                "ok": True,
                "count": len(matches),
                "tools": [tool.to_summary() for tool in matches],
            }
        except Exception as exc:
            return _error_payload(None, exc)

    @server.tool(
        name="describe_dnakit_tool",
        description="Return the complete input schema and safety metadata for one DNAKit tool.",
        annotations={
            "readOnlyHint": True,
            "destructiveHint": False,
            "idempotentHint": True,
            "openWorldHint": False,
        },
    )
    def describe_tool(tool_name: str) -> dict[str, object]:
        try:
            return {"ok": True, "tool": resolved_registry.describe(tool_name)}
        except Exception as exc:
            return _error_payload(tool_name, exc)

    @server.tool(
        name="call_dnakit_tool",
        description=(
            "Call one public DNAKit function with JSON arguments. Run describe_dnakit_tool "
            "first. Set allow_side_effects only after reviewing file, model, or process effects."
        ),
        annotations={
            "readOnlyHint": False,
            "destructiveHint": True,
            "idempotentHint": False,
            "openWorldHint": True,
        },
    )
    def call_tool(
        tool_name: str,
        arguments: dict[str, object] | None = None,
        allow_side_effects: bool = False,
    ) -> dict[str, object]:
        started = time.perf_counter()
        try:
            result = resolved_registry.execute(
                tool_name,
                arguments,
                allow_side_effects=allow_side_effects,
            )
            return {
                "ok": True,
                "tool": resolved_registry.resolve(tool_name).name,
                "elapsed_ms": round((time.perf_counter() - started) * 1_000, 3),
                "result": result,
            }
        except Exception as exc:
            return _error_payload(tool_name, exc)

    return server


def main() -> None:
    """Run the local stdio MCP server used by Agent clients."""

    create_server().run(show_banner=False)


if __name__ == "__main__":
    main()


__all__ = ["MCPServer", "create_server", "main"]
