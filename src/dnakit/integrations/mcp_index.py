"""Offline search over the DNA- and gene-related tools of catalogued MCP servers.

The index is a snapshot of tools/list from the third-party servers in
mcp_servers.json, reduced to tools whose main input or output is DNA or genes.
Searching it starts no server and reads no credentials; call a tool with
call_mcp_tool after installing its server.
"""

from __future__ import annotations

import copy
import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, cast

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.tools.retrieval import BM25Config, _bm25_scores

_INDEX = Path(__file__).with_name("mcp_tool_index.json")
_CATALOG = Path(__file__).with_name("mcp_servers.json")
_SPLIT = re.compile(r"[._\-/]+|(?<=[a-z])(?=[A-Z])")


@lru_cache(maxsize=1)
def _index() -> dict[str, Any]:
    return cast(dict[str, Any], json.loads(_INDEX.read_text(encoding="utf-8")))


@lru_cache(maxsize=1)
def _servers() -> dict[str, dict[str, Any]]:
    data = json.loads(_CATALOG.read_text(encoding="utf-8"))
    return {server["id"]: server for server in data["servers"]}


def _call(tool: dict[str, Any]) -> dict[str, Any]:
    """How to invoke an indexed tool through dnakit.integrations.call_mcp_tool."""
    via = tool.get("via")
    if via is None:
        return {
            "function": "dnakit.integrations.call_mcp_tool",
            "server_id": tool["server_id"],
            "tool_name": tool["name"],
            "arguments": "the tool's own arguments (see input_schema)",
        }
    return {
        "function": "dnakit.integrations.call_mcp_tool",
        "server_id": tool["server_id"],
        "tool_name": via["tool_name"],
        "arguments": {
            via["name_argument"]: tool["name"],
            via["arguments_argument"]: "the tool's own arguments (see input_schema)",
        },
    }


def _summary(tool: dict[str, Any]) -> dict[str, Any]:
    server = _servers()[tool["server_id"]]
    schema = tool.get("input_schema") or {}
    summary = {
        "server_id": tool["server_id"],
        "tool_name": tool["name"],
        "category": tool["category"],
        "summary": tool["summary"],
        "required_arguments": list(schema.get("required", [])),
        "listing": tool["method"],
        "server_requires": {
            "environment": list(server.get("required_environment", [])),
            "installation": server.get("installation", ""),
        },
        "call": _call(tool),
    }
    if "requires" in tool:
        summary["tool_requires"] = copy.deepcopy(tool["requires"])
    return summary


def _document(tool: dict[str, Any]) -> str:
    server = _servers()[tool["server_id"]]
    name_words = " ".join(_SPLIT.split(tool["name"]))
    return " ".join(
        [
            tool["name"],
            name_words,
            name_words,
            tool["summary"],
            " ".join(tool.get("keywords", [])),
            tool["category"].replace("_", " "),
            tool.get("description", ""),
            tool["server_id"],
            server.get("name", ""),
        ]
    )


@lru_cache(maxsize=1)
def _documents() -> tuple[tuple[str, str], ...]:
    return tuple(
        (str(position), _document(tool)) for position, tool in enumerate(_index()["tools"])
    )


def _check_filters(server_id: str | None, category: str | None) -> None:
    if server_id is not None and (not isinstance(server_id, str) or server_id not in _servers()):
        raise ConfigurationError(
            "Unknown MCP server; use list_mcp_servers().",
            code="MCP_INDEX_SERVER_INVALID",
            context={"server_id": server_id},
        )
    if category is not None and (
        not isinstance(category, str) or category not in _index()["categories"]
    ):
        raise ConfigurationError(
            "Unknown MCP tool category.",
            code="MCP_INDEX_CATEGORY_INVALID",
            context={"category": category, "valid": sorted(_index()["categories"])},
        )


def search_mcp_tools(
    query: str,
    *,
    server_id: str | None = None,
    category: str | None = None,
    limit: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rank indexed DNA/gene tools of third-party MCP servers with local BM25.

    Only tools judged DNA- or gene-related are indexed; server housekeeping,
    literature, drug and protein-only tools are left out. query is English
    text; an empty query with server_id or category lists those tools in
    index order. No server is started and no network is used.
    """
    if not isinstance(query, str) or len(query) > 1000:
        raise ConfigurationError(
            "query must be text of at most 1000 characters.", code="MCP_INDEX_QUERY_INVALID"
        )
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
        raise ConfigurationError(
            "limit must be an integer in [1, 100].", code="MCP_INDEX_LIMIT_INVALID"
        )
    _check_filters(server_id, category)
    tools = _index()["tools"]

    def wanted(tool: dict[str, Any]) -> bool:
        return (server_id is None or tool["server_id"] == server_id) and (
            category is None or tool["category"] == category
        )

    if query.strip():
        scores = _bm25_scores(query, list(_documents()), BM25Config())
        order = sorted(scores, key=lambda key: (-scores[key], int(key)))
        matches = [tools[int(key)] for key in order if wanted(tools[int(key)])]
    elif server_id is not None or category is not None:
        matches = [tool for tool in tools if wanted(tool)]
    else:
        raise ConfigurationError(
            "Give a query, a server_id or a category.", code="MCP_INDEX_QUERY_INVALID"
        )
    return {
        "query": query,
        "count": min(len(matches), limit),
        "total_matches": len(matches),
        "indexed_tools": len(tools),
        "tools": [_summary(tool) for tool in matches[:limit]],
    }


def describe_mcp_tool(server_id: str, tool_name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return one indexed MCP tool's description, input schema and server facts.

    The schema is the snapshot taken when the index was built; call_mcp_tool
    validates against the live schema of the installed server.
    """
    _check_filters(server_id, None)
    for tool in _index()["tools"]:
        if tool["server_id"] == server_id and tool["name"] == tool_name:
            server = _servers()[server_id]
            result = _summary(tool)
            result.update(
                {
                    "description": tool.get("description", ""),
                    "keywords": list(tool.get("keywords", [])),
                    "input_schema": copy.deepcopy(tool.get("input_schema") or {}),
                    "annotations": copy.deepcopy(tool.get("annotations") or {}),
                    "server": {
                        key: copy.deepcopy(server.get(key))
                        for key in (
                            "name",
                            "repository",
                            "license",
                            "source_commit",
                            "command",
                            "arguments",
                            "transports",
                            "installation",
                            "required_environment",
                            "secret_environment",
                            "effects",
                        )
                    },
                }
            )
            return result
    excluded = next(
        (
            item
            for item in _index()["excluded"]
            if item["server_id"] == server_id and item["name"] == tool_name
        ),
        None,
    )
    raise ConfigurationError(
        "This tool is not in the DNA/gene MCP tool index.",
        code="MCP_INDEX_TOOL_NOT_INDEXED",
        context={
            "server_id": server_id,
            "tool_name": tool_name,
            "excluded_as": None if excluded is None else excluded["category"],
            "reason": None if excluded is None else excluded["reason"],
        },
        hint="list_mcp_tools(server_id) lists every tool of an installed server.",
    )


def mcp_tool_index_summary(*, verbose: bool = True) -> dict[str, Any]:
    """Count listed, indexed and excluded MCP tools per server and category.

    listed_tools is what each server returned from tools/list (or its source
    at the pinned commit); indexed_tools keeps the DNA/gene-related subset.
    """
    data = _index()
    by_category: dict[str, int] = {code: 0 for code in data["categories"]}
    for tool in data["tools"]:
        by_category[tool["category"]] += 1
    excluded_by_reason: dict[str, int] = {}
    for item in data["excluded"]:
        excluded_by_reason[item["category"]] = excluded_by_reason.get(item["category"], 0) + 1
    servers = copy.deepcopy(data["servers"])
    return {
        "built_date": data["built_date"],
        "scope": data["scope"],
        "catalog_servers": len(_servers()),
        "listed_servers": sum(1 for server in servers if server["listed_tools"] is not None),
        "servers_with_indexed_tools": sum(1 for server in servers if server["indexed_tools"]),
        "listed_tools": sum(server["listed_tools"] or 0 for server in servers),
        "indexed_tools": len(data["tools"]),
        "excluded_tools": len(data["excluded"]),
        "indexed_by_category": by_category,
        "excluded_by_reason": dict(sorted(excluded_by_reason.items())),
        "categories": copy.deepcopy(data["categories"]),
        "servers": servers,
    }


__all__ = ["describe_mcp_tool", "mcp_tool_index_summary", "search_mcp_tools"]
_register_public_verbose(globals(), tuple(__all__))
