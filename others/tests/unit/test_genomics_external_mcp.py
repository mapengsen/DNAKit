"""Real stdio protocol tests and offline MCP bridge failure-path checks."""

from __future__ import annotations

import asyncio
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.integrations import mcp as bridge

_SERVER = """
import asyncio
import os
import sys
from pathlib import Path
from fastmcp import FastMCP
from fastmcp.exceptions import ToolError
server = FastMCP("DNAKit local protocol fixture")
secret = os.environ.get("SERVICE_TOKEN", "")
counter = 0
print("credential: " + secret, file=sys.stderr)
@server.tool
def add(a: int, b: int) -> dict:
    Path("called").write_text("yes")
    return {"sum": a + b, "credential": secret}
@server.tool
def remote_failure() -> str:
    raise ToolError("failed " + secret)
@server.tool
def plain_failure() -> str:
    return "HTTP error occurred: 500 - " + secret
@server.tool
def json_failure() -> dict:
    return {"error": "provider failure " + secret}
@server.tool
def increment(value: int) -> dict:
    global counter
    counter += value
    return {"counter": counter}
@server.tool
async def slow() -> str:
    await asyncio.sleep(10)
    return "finished"
server.run(transport="stdio", show_banner=False)
"""


@pytest.fixture
def local_server(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    pytest.importorskip("fastmcp")
    (tmp_path / "fixture_server.py").write_text(_SERVER)
    definition = {
        "id": "fixture",
        "name": "Fixture",
        "command": sys.executable,
        "arguments": ["fixture_server.py"],
        "required_files": ["fixture_server.py"],
        "required_environment": [],
        "optional_environment": ["SERVICE_TOKEN"],
        "secret_environment": ["SERVICE_TOKEN"],
        "installation": "local test fixture",
        "application_error_prefixes": ["HTTP error occurred:"],
        "application_error_keys": ["error", "errors"],
    }
    monkeypatch.setattr(bridge, "_definitions", lambda: [definition])
    monkeypatch.setenv("DNAKIT_TEST_TOKEN", "a-unique-secret-for-mcp-tests")
    monkeypatch.setenv(
        "DNAKIT_TEST_PYTHONPATH",
        os.pathsep.join(
            str(Path(path).resolve()) for path in sys.path if path and Path(path).is_dir()
        ),
    )
    return {
        "server_id": "fixture",
        "workdir": str(tmp_path),
        "environment_variables": {
            "SERVICE_TOKEN": "DNAKIT_TEST_TOKEN",
            "PYTHONPATH": "DNAKIT_TEST_PYTHONPATH",
        },
        "timeout_seconds": 30,
        "verbose": False,
    }


def test_catalog_is_offline_and_verified(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bridge, "_sdk", lambda: pytest.fail("SDK must not be imported"))
    servers = bridge.list_mcp_servers(verbose=False)
    assert {server["id"] for server in servers} >= {
        "gigwa-mcp",
        "chorus",
        "mavedb-link",
        "biomcp",
        "juicebox",
        "genechat",
        "genepattern",
    }
    assert len({server["id"] for server in servers}) == len(servers)
    assert all(server["command"] not in {"uvx", "npx", "pip"} for server in servers)
    assert all(server["sources"] for server in servers)
    mavedb = next(s for s in servers if s["id"] == "mavedb-link")
    assert mavedb["http_url"] == "https://mavedb-link.genefoundry.org/mcp"


def test_expanded_catalog_has_pinned_entrypoints_and_no_install_runners() -> None:
    catalog = {server["id"]: server for server in bridge.list_mcp_servers(verbose=False)}
    assert catalog.keys() >= {
        "brapi",
        "ensembl",
        "biothings",
        "gget",
        "ucsc-genome",
        "opentargets",
        "plant-genomics",
        "biocontext-kb",
        "grebi",
        "ebi-tools",
        "gwas-mcp",
        "ncbi-database",
        "pubmed",
        "encode-toolkit",
        "tooluniverse",
        "kiro-genomics",
        "pubmed-anthropic",
        "opentargets-platform",
        "10x-genomics",
    }
    for server in catalog.values():
        assert len(server["source_commit"]) == 40
        assert all(char in "0123456789abcdef" for char in server["source_commit"])
        assert server["sources"] and server["license"] and server["attribution"]
        if "stdio" in server["transports"]:
            assert server["command"] not in {"uv", "uvx", "npx", "npm", "bunx", "pip"}
            assert isinstance(server["command"], str)
        else:
            assert server["command"] is None and server["http_url"].startswith("https://")
        assert not set(server.get("fixed_environment", {})) & set(server["secret_environment"])
    assert catalog["ucsc-genome"]["arguments"] == ["ucsc-genome-mcp.py"]
    assert catalog["encode-toolkit"]["arguments"] == ["-m", "encode_connector.server.main"]
    assert catalog["gget"]["arguments"] == ["stdio"]
    assert {"kiro-microbiology", "kiro-molbio"}.isdisjoint(catalog)
    for identifier in ["pubmed-anthropic", "opentargets-platform", "10x-genomics"]:
        assert catalog[identifier]["project_id"] == "anthropic-life-sciences"


def test_fixed_transport_environment_cannot_inherit_or_accept_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = bridge._definition("biocontext-kb")
    monkeypatch.setenv("MCP_ENVIRONMENT", "PRODUCTION")
    env, secrets = bridge._environment(server, None)
    assert env["MCP_ENVIRONMENT"] == "DEVELOPMENT" and secrets == []
    with pytest.raises(ConfigurationError, match="mapping"):
        bridge._environment(server, {"MCP_ENVIRONMENT": "MY_TRANSPORT"})


def test_http_only_recipe_rejects_stdio_before_startup() -> None:
    with pytest.raises(ConfigurationError, match="transport"):
        bridge._connection(bridge._definition("grebi"), "stdio", None, None, None)
    config, secrets = bridge._connection(bridge._definition("grebi"), "http", None, None, None)
    assert config == {"url": "https://wwwdev.ebi.ac.uk/kg/api/v1/mcp"}
    assert secrets == []


def test_missing_socks_support_is_actionable_without_exposing_proxy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("HTTPS_PROXY", "socks5://proxy-secret@example.test:1080")
    monkeypatch.setattr(bridge, "proxy_bypass", lambda _: False)
    monkeypatch.setattr(
        bridge, "getproxies", lambda: {"https": "socks5://proxy-secret@example.test:1080"}
    )

    def missing(name: str) -> None:
        assert name == "socksio"
        raise ImportError("missing proxy dependency")

    monkeypatch.setattr(bridge.importlib, "import_module", missing)
    with pytest.raises(BackendUnavailableError, match="SOCKS") as error:
        bridge._http_proxy_support("https://mcp.example.test/mcp")
    assert "httpx[socks]" in str(error.value)
    assert "proxy-secret" not in str(error.value)


def test_https_proxy_precedence_and_no_proxy_do_not_require_socks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        bridge, "getproxies", lambda: {"https": "http://proxy.test", "all": "socks5://other.test"}
    )
    monkeypatch.setattr(bridge, "proxy_bypass", lambda _: False)
    monkeypatch.setattr(
        bridge.importlib, "import_module", lambda _: pytest.fail("SOCKS dependency not needed")
    )
    bridge._http_proxy_support("https://mcp.example.test/mcp")
    monkeypatch.setattr(bridge, "getproxies", lambda: {"all": "socks5://other.test"})
    monkeypatch.setattr(bridge, "proxy_bypass", lambda _: True)
    bridge._http_proxy_support("https://mcp.example.test/mcp")


def test_real_stdio_tool_listing_and_call(
    local_server: dict[str, Any], capsys: pytest.CaptureFixture
) -> None:
    listing = bridge.list_mcp_tools(**local_server)
    assert {tool["name"] for tool in listing["tools"]} == {
        "add",
        "remote_failure",
        "slow",
        "increment",
        "plain_failure",
        "json_failure",
    }
    schema = next(tool["inputSchema"] for tool in listing["tools"] if tool["name"] == "add")
    assert schema["properties"]["a"]["type"] == "integer"
    result = bridge.call_mcp_tool(tool_name="add", arguments={"a": 2, "b": 3}, **local_server)
    assert result["status"] == "completed" and result["isError"] is False
    assert result["structuredContent"]["sum"] == 5
    assert result["structuredContent"]["credential"] == "[REDACTED]"
    assert "a-unique-secret-for-mcp-tests" not in json.dumps(result)
    assert (Path(local_server["workdir"]) / "called").read_text() == "yes"
    captured = capsys.readouterr()
    assert "a-unique-secret-for-mcp-tests" not in captured.out + captured.err


def test_real_remote_errors_are_not_success(local_server: dict[str, Any]) -> None:
    result = bridge.call_mcp_tool(tool_name="remote_failure", arguments={}, **local_server)
    assert result["isError"] is True and result["status"] == "remote_error"
    assert result["content"][0]["text"] == "failed [REDACTED]"


def test_real_source_verified_text_errors_are_not_success(local_server: dict[str, Any]) -> None:
    result = bridge.call_mcp_tool(tool_name="plain_failure", arguments={}, **local_server)
    assert result["isError"] is True and result["status"] == "remote_error"
    assert result["upstream_isError"] is False
    assert result["error_detection"] == "verified_upstream_text_prefix"
    assert "a-unique-secret-for-mcp-tests" not in json.dumps(result)


def test_real_source_verified_json_errors_stop_workflow(local_server: dict[str, Any]) -> None:
    result = bridge.call_mcp_workflow(
        calls=[
            {"tool_name": "json_failure", "arguments": {}},
            {"tool_name": "increment", "arguments": {"value": 1}},
        ],
        **local_server,
    )
    assert result["status"] == "remote_error" and result["executed_steps"] == 1
    assert result["steps"][0]["upstream_isError"] is False
    assert result["steps"][0]["error_detection"] == "verified_upstream_json_error"
    assert "a-unique-secret-for-mcp-tests" not in json.dumps(result)


def test_real_workflow_keeps_state_and_stops_at_remote_error(local_server: dict[str, Any]) -> None:
    result = bridge.call_mcp_workflow(
        calls=[
            {"tool_name": "increment", "arguments": {"value": 2}},
            {"tool_name": "increment", "arguments": {"value": 3}},
            {"tool_name": "remote_failure", "arguments": {}},
            {"tool_name": "add", "arguments": {"a": 4, "b": 5}},
        ],
        **local_server,
    )
    assert result["status"] == "remote_error" and result["isError"] is True
    assert result["requested_steps"] == 4 and result["executed_steps"] == 3
    assert result["steps"][0]["structuredContent"]["counter"] == 2
    assert result["steps"][1]["structuredContent"]["counter"] == 5
    assert result["steps"][2]["isError"] is True
    assert "a-unique-secret-for-mcp-tests" not in json.dumps(result)
    assert not (Path(local_server["workdir"]) / "called").exists()


@pytest.mark.parametrize(
    "calls",
    [[], [{}], [{"tool_name": "add", "arguments": []}], [{"tool_name": "", "arguments": {}}]],
)
def test_invalid_workflow_precedes_connection(calls: list[dict[str, Any]]) -> None:
    with pytest.raises(ConfigurationError):
        bridge.call_mcp_workflow("biomcp", calls, verbose=False)


@pytest.mark.parametrize(
    "name,arguments", [("missing", {}), ("add", {"a": "2", "b": 3}), ("add", {"a": 2})]
)
def test_real_schema_validation_precedes_invocation(
    local_server: dict[str, Any],
    name: str,
    arguments: dict[str, Any],
) -> None:
    with pytest.raises(ConfigurationError):
        bridge.call_mcp_tool(tool_name=name, arguments=arguments, **local_server)
    assert not (Path(local_server["workdir"]) / "called").exists()


def test_missing_sdk_has_actionable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    def missing(name: str) -> None:
        raise ImportError("not installed")

    monkeypatch.setattr(bridge.importlib, "import_module", missing)
    with pytest.raises(BackendUnavailableError, match="FastMCP"):
        bridge._sdk()


def test_unknown_missing_and_source_workdir(monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(ConfigurationError, match="Unknown MCP"):
        bridge.list_mcp_tools("no-such-server", verbose=False)
    with pytest.raises(BackendUnavailableError):
        bridge.list_mcp_tools("biomcp", executable_path="/no/such/program", verbose=False)
    with pytest.raises(ConfigurationError, match="workdir"):
        bridge.list_mcp_tools("mavedb-link", executable_path=sys.executable, verbose=False)
    monkeypatch.setattr(bridge.shutil, "which", lambda _: "/usr/bin/uvx")
    with pytest.raises(ConfigurationError, match="installers"):
        bridge.list_mcp_tools("biomcp", executable_path="uvx", verbose=False)


def test_named_credentials_and_pairs(monkeypatch: pytest.MonkeyPatch) -> None:
    server = bridge._definition("gigwa-mcp")
    monkeypatch.setenv("GIGWA_URL", "https://example.test/gigwa")
    monkeypatch.setenv("GIGWA_USER", "user")
    monkeypatch.delenv("GIGWA_PASS", raising=False)
    with pytest.raises(ConfigurationError, match="both"):
        bridge._environment(server, None)
    with pytest.raises(ConfigurationError, match="mapping"):
        bridge._environment(server, {"GIGWA_PASS": "literal password with spaces"})
    with pytest.raises(ConfigurationError):
        bridge._environment(server, {"LD_PRELOAD": "MY_LIBRARY"})


def test_event_loop_rejection_does_not_create_coroutine() -> None:
    called = []

    async def operation() -> dict[str, Any]:
        called.append(True)
        return {}

    async def scenario() -> None:
        with pytest.raises(ConfigurationError, match="event loop"):
            bridge._sync(operation, 1)

    asyncio.run(scenario())
    assert called == []


def test_timeout_and_sdk_exception_are_sanitized() -> None:
    async def slow() -> dict[str, Any]:
        await asyncio.sleep(1)
        return {}

    with pytest.raises(BackendTimeoutError):
        bridge._sync(slow, 0.01)

    async def broken() -> dict[str, Any]:
        raise RuntimeError("secret-credential-from-sdk")

    with pytest.raises(BackendExecutionError) as error:
        bridge._sync(broken, 1)
    assert "secret-credential-from-sdk" not in str(error.value)
    assert error.value.__suppress_context__ is True


def test_sdk_diagnostics_are_suppressed_only_inside_bridge(
    caplog: pytest.LogCaptureFixture,
) -> None:
    logger = logging.getLogger("mcp.client.stdio")
    with bridge._quiet_sdk_logs():
        logger.error("secret-credential-in-protocol-parser")
    logger.warning("unrelated-client-diagnostic")
    assert "secret-credential" not in caplog.text
    assert "unrelated-client-diagnostic" in caplog.text


def test_schema_references_do_not_fetch_urls() -> None:
    with pytest.raises(BackendExecutionError, match="schema"):
        bridge._validate_arguments({}, {"$ref": "https://example.test/schema.json"})
    bridge._validate_arguments(
        {"n": 3},
        {
            "type": "object",
            "properties": {"n": {"$ref": "#/$defs/n"}},
            "$defs": {"n": {"type": "integer"}},
        },
    )


@pytest.mark.parametrize(
    "arguments", [{"value": float("nan")}, {"value": object()}, {1: "bad-key"}]
)
def test_non_json_arguments_precede_connection(arguments: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        bridge.call_mcp_tool("biomcp", "unused", arguments, verbose=False)


def test_limits_and_unverified_http_rejected() -> None:
    with pytest.raises(ConfigurationError):
        bridge.list_mcp_tools("biomcp", timeout_seconds=float("inf"), verbose=False)
    with pytest.raises(ConfigurationError, match="HTTP"):
        bridge.list_mcp_tools("biomcp", transport="http", verbose=False)
    with pytest.raises(BackendExecutionError, match="limit"):
        bridge._json_size({"value": "large"}, 1)
