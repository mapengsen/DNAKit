"""Offline index of DNA/gene tools from the catalogued third-party MCP servers."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, cast

import pytest

import dnakit.integrations as integrations
from dnakit.exceptions import ConfigurationError
from dnakit.integrations import (
    describe_mcp_tool,
    integration_inventory,
    list_mcp_servers,
    mcp_tool_index_summary,
    search_mcp_tools,
)
from dnakit.tools import server as server_module
from dnakit.tools.registry import AgentToolRegistry

INDEX = Path(integrations.__file__).with_name("mcp_tool_index.json")
EXCLUDED_CODES = {
    "server_meta",
    "literature",
    "chemistry_drug",
    "protein_only",
    "clinical_nongenetic",
    "ontology_nongene",
    "cell_level",
    "generic_platform",
    "other",
}


@pytest.fixture(scope="module")
def data() -> dict[str, Any]:
    return cast(dict[str, Any], json.loads(INDEX.read_text(encoding="utf-8")))


def test_every_catalog_server_has_a_listing_record(data) -> None:
    catalog = {server["id"]: server for server in list_mcp_servers(verbose=False)}
    records = {server["server_id"]: server for server in data["servers"]}
    assert records.keys() == catalog.keys()
    for server_id, record in records.items():
        assert record["source_commit"] == catalog[server_id]["source_commit"], server_id
        assert record["listed_tools"] is not None, server_id
        assert record["listing"], server_id


def test_indexed_and_excluded_tools_partition_each_listing(data) -> None:
    indexed: dict[str, int] = {}
    excluded: dict[str, int] = {}
    for tool in data["tools"]:
        indexed[tool["server_id"]] = indexed.get(tool["server_id"], 0) + 1
    for item in data["excluded"]:
        excluded[item["server_id"]] = excluded.get(item["server_id"], 0) + 1
        assert item["category"] in EXCLUDED_CODES
        assert item["reason"]
    for record in data["servers"]:
        server_id = record["server_id"]
        assert record["indexed_tools"] == indexed.get(server_id, 0), server_id
        assert indexed.get(server_id, 0) + excluded.get(server_id, 0) == record["listed_tools"]
    keys = [(tool["server_id"], tool["name"]) for tool in data["tools"]]
    keys += [(item["server_id"], item["name"]) for item in data["excluded"]]
    assert len(keys) == len(set(keys))


def test_indexed_entries_are_complete(data) -> None:
    for tool in data["tools"]:
        assert tool["category"] in data["categories"], tool["name"]
        assert tool["summary"] and len(tool["summary"]) <= 200
        assert isinstance(tool["input_schema"], dict)
        assert tool["method"] in {"live-stdio", "live-http", "source", "source-catalog"}
        if tool["method"] == "source-catalog":
            assert tool["via"]["tool_name"], tool["name"]


def test_summary_counts_match_the_index(data) -> None:
    summary = mcp_tool_index_summary(verbose=False)
    assert summary["catalog_servers"] == 49
    assert summary["listed_servers"] == 49
    assert summary["indexed_tools"] == len(data["tools"])
    assert summary["excluded_tools"] == len(data["excluded"])
    assert summary["listed_tools"] == len(data["tools"]) + len(data["excluded"])
    assert sum(summary["indexed_by_category"].values()) == summary["indexed_tools"]
    inventory = integration_inventory(verbose=False)
    assert inventory["mcp_indexed_dna_gene_tools"] == summary["indexed_tools"]
    assert inventory["mcp_listed_tools"] == summary["listed_tools"]


@pytest.mark.parametrize(
    ("query", "server_id"),
    [
        ("addgene plasmid", "addgene"),
        ("juicebox hi-c contact map", "juicebox"),
        ("polygenic risk score", None),
        ("gnomad allele frequency", None),
    ],
)
def test_search_finds_tools_only_reachable_through_mcp(query, server_id) -> None:
    result = search_mcp_tools(query, limit=10, verbose=False)
    assert result["count"] > 0
    if server_id is not None:
        assert server_id in {tool["server_id"] for tool in result["tools"]}
    top = result["tools"][0]
    assert top["call"]["function"] == "dnakit.integrations.call_mcp_tool"


def test_search_filters_and_listing_without_query() -> None:
    listed = search_mcp_tools("", server_id="gnomad", limit=100, verbose=False)
    assert listed["count"] == listed["total_matches"] > 0
    assert {tool["server_id"] for tool in listed["tools"]} == {"gnomad"}
    variants = search_mcp_tools("frequency", category="variant", limit=20, verbose=False)
    assert all(tool["category"] == "variant" for tool in variants["tools"])
    with pytest.raises(ConfigurationError) as caught:
        search_mcp_tools("", verbose=False)
    assert caught.value.code == "MCP_INDEX_QUERY_INVALID"
    with pytest.raises(ConfigurationError) as caught:
        search_mcp_tools("gene", server_id="no-such-server", verbose=False)
    assert caught.value.code == "MCP_INDEX_SERVER_INVALID"
    with pytest.raises(ConfigurationError) as caught:
        search_mcp_tools("gene", category="drugs", verbose=False)
    assert caught.value.code == "MCP_INDEX_CATEGORY_INVALID"
    with pytest.raises(ConfigurationError):
        search_mcp_tools("gene", limit=0, verbose=False)


def test_dispatched_tools_describe_the_wrapper_call(data) -> None:
    wrapped = next(tool for tool in data["tools"] if tool.get("via"))
    described = describe_mcp_tool(wrapped["server_id"], wrapped["name"], verbose=False)
    call = described["call"]
    assert call["tool_name"] == wrapped["via"]["tool_name"]
    assert call["arguments"][wrapped["via"]["name_argument"]] == wrapped["name"]
    assert described["input_schema"] == wrapped["input_schema"]
    assert described["server"]["repository"]


def test_describe_explains_excluded_tools(data) -> None:
    item = data["excluded"][0]
    with pytest.raises(ConfigurationError) as caught:
        describe_mcp_tool(item["server_id"], item["name"], verbose=False)
    assert caught.value.code == "MCP_INDEX_TOOL_NOT_INDEXED"
    assert caught.value.context["excluded_as"] == item["category"]


class _FakeServer:
    def __init__(self, **settings: object) -> None:
        self.tools: dict[str, Callable[..., object]] = {}

    def tool(
        self, *, name: str, description: str, annotations: Mapping[str, object]
    ) -> Callable[[Callable[..., object]], Callable[..., object]]:
        def register(function: Callable[..., object]) -> Callable[..., object]:
            self.tools[name] = function
            return function

        return register

    def run(self, *, show_banner: bool | None = None) -> None:
        del show_banner


def test_agent_server_search_appends_mcp_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = _FakeServer()
    monkeypatch.setattr(
        "dnakit.tools.server.importlib.import_module",
        lambda name: type("FakeFastMCP", (), {"FastMCP": lambda **kwargs: fake}),
    )
    server_module.create_server(AgentToolRegistry(), verbose=False)
    search = fake.tools["search_dnakit_tools"]
    result = cast(dict[str, Any], search("addgene plasmid"))
    assert result["ok"] is True
    assert 0 < len(result["mcp_tools"]) <= 5
    assert "addgene" in {tool["server_id"] for tool in result["mcp_tools"]}
    assert "mcp_tools" not in cast(dict[str, Any], search("addgene plasmid", mcp_limit=0))
    assert "mcp_tools" not in cast(dict[str, Any], search("plasmid", category="databases"))
    bad = cast(dict[str, Any], search("plasmid", mcp_limit=99))
    assert bad["error"]["code"] == "AGENT_TOOL_LIMIT_INVALID"
