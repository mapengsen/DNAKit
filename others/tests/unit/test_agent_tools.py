"""Tests for compact Agent discovery and MCP-safe DNAKit execution."""

from __future__ import annotations

import importlib
import inspect
from collections.abc import Callable, Mapping
from typing import cast

import pytest

import dnakit.tools.server as server_module
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import PropertyPredictionConfig
from dnakit.tools import AgentToolRegistry, default_tool_registry
from dnakit.tools.adapters import coerce_value
from dnakit.tools.catalog import PUBLIC_TOOL_MODULES


def test_catalog_contains_every_explicitly_exported_public_function() -> None:
    registry = default_tool_registry()
    expected: set[str] = set()
    for module_name in PUBLIC_TOOL_MODULES:
        module = importlib.import_module(module_name)
        expected.update(
            f"{module_name}.{name}"
            for name in getattr(module, "__all__", ())
            if isinstance(name, str) and inspect.isfunction(getattr(module, name, None))
        )

    assert {tool.name for tool in registry.tools} == expected
    stats = registry.stats()
    total = cast(int, stats["total"])
    assert total == len(expected)
    assert total >= 322
    assert stats["agent_compatible"] == total - 3
    assert stats["python_only"] == 3
    assert stats["categories"] == 28


def test_catalog_excludes_callbacks_and_the_low_level_command_executor() -> None:
    registry = default_tool_registry()

    for name in (
        "dnakit.batch.iter_batch",
        "dnakit.batch.run_batch",
        "dnakit.backends.execute_bounded_command",
    ):
        manifest = registry.describe(name)
        assert manifest["agent_compatible"] is False
        with pytest.raises(ConfigurationError) as caught:
            registry.execute(name, {})
        assert caught.value.code == "AGENT_TOOL_PYTHON_ONLY"


def test_search_describe_and_execute_native_function() -> None:
    registry = default_tool_registry()
    matches = registry.search("molecular weight")

    assert matches[0].name == "dnakit.thermodynamics.molecular_weight"
    manifest = registry.describe("thermodynamics.molecular_weight")
    assert manifest["agent_compatible"] is True
    assert manifest["input_schema"]["required"] == ["sequence"]  # type: ignore[index]

    result = registry.execute(
        "dnakit.thermodynamics.molecular_weight",
        {"sequence": "ACGT", "strand": "single"},
    )
    assert isinstance(result, dict)
    assert result["value_dalton"] == pytest.approx(1173.84)
    assert result["strand_count"] == 1


def test_agent_adapter_constructs_nested_prediction_configuration() -> None:
    value = coerce_value(
        {
            "model": "enformer",
            "task": "human_tracks",
            "batch_size": 2,
            "show_progress": False,
        },
        PropertyPredictionConfig,
        path="config",
    )

    assert isinstance(value, PropertyPredictionConfig)
    assert value.model == "enformer"
    assert value.task == "human_tracks"
    assert value.batch_size == 2
    assert value.show_progress is False


def test_enformer_benchmark_prediction_is_registered_as_model_tool() -> None:
    manifest = default_tool_registry().describe("predictions.predict_enformer_benchmark")

    assert manifest["agent_compatible"] is True
    assert manifest["effect"] == "model"
    assert manifest["requires_confirmation"] is True


def test_agent_result_contains_plain_json_values() -> None:
    result = default_tool_registry().execute(
        "dnakit.ops.reverse_complement",
        {"value": "AACG"},
    )

    assert isinstance(result, dict)
    records = result["records"]
    assert isinstance(records, list)
    sequence = records[0]["sequence"]
    assert sequence["parts"] == ["CGTT"]
    assert sequence["alphabet"] == "strict"
    assert sequence["topology"] == "linear"


def test_serialized_results_can_feed_follow_up_agent_tools() -> None:
    registry = default_tool_registry()
    index = registry.execute(
        "dnakit.similarity.build_sketch_index",
        {
            "records": [
                {"id": "a", "sequence": "ACGTACGT"},
                {"id": "b", "sequence": "ACGTTCGT"},
            ],
            "k": 3,
            "num_hashes": 20,
        },
    )
    neighbors = registry.execute(
        "dnakit.similarity.nearest_neighbors",
        {"query": "ACGTACGT", "index": index, "top_k": 2},
    )

    assert isinstance(neighbors, dict)
    hits = cast(list[dict[str, object]], neighbors["hits"])
    assert [hit["record_id"] for hit in hits] == ["a", "b"]

    alignment = registry.execute(
        "dnakit.alignment.align_pairwise",
        {"query": "ACGT", "target": "ACGA"},
    )
    plot = registry.execute(
        "dnakit.visualization.plot_alignment",
        {"result": alignment},
    )
    assert isinstance(plot, dict)
    assert plot["kind"] == "alignment"
    assert str(plot["svg"]).startswith("<svg")


def test_side_effecting_tools_require_explicit_authorization() -> None:
    registry = default_tool_registry()

    with pytest.raises(ConfigurationError) as caught:
        registry.execute("dnakit.visualization.save_svg", {})

    assert caught.value.code == "AGENT_TOOL_CONFIRMATION_REQUIRED"


class _FakeMCPServer:
    def __init__(self, **settings: object) -> None:
        self.settings = settings
        self.tools: dict[str, Callable[..., object]] = {}
        self.ran = False
        self.show_banner: bool | None = None

    def tool(
        self,
        *,
        name: str,
        description: str,
        annotations: Mapping[str, object],
    ) -> Callable[[Callable[..., object]], Callable[..., object]]:
        del description, annotations

        def register(function: Callable[..., object]) -> Callable[..., object]:
            self.tools[name] = function
            return function

        return register

    def run(self, *, show_banner: bool | None = None) -> None:
        self.ran = True
        self.show_banner = show_banner


def test_compact_server_registers_discovery_and_call_tools(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = _FakeMCPServer()
    registry = AgentToolRegistry()

    def fake_import(name: str) -> object:
        assert name == "fastmcp"
        return type("FakeFastMCPModule", (), {"FastMCP": lambda **kwargs: fake})

    monkeypatch.setattr("dnakit.tools.server.importlib.import_module", fake_import)
    server = server_module.create_server(registry)

    assert server is fake
    assert set(fake.tools) == {
        "call_dnakit_tool",
        "describe_dnakit_tool",
        "dnakit_catalog_stats",
        "list_dnakit_categories",
        "list_dnakit_tools",
        "search_dnakit_tools",
    }
    response = fake.tools["call_dnakit_tool"](
        "dnakit.thermodynamics.molecular_weight",
        {"sequence": "ACGT"},
    )
    assert isinstance(response, dict)
    assert response["ok"] is True


def test_server_reports_missing_optional_dependency(monkeypatch: pytest.MonkeyPatch) -> None:
    registry = AgentToolRegistry()

    def missing_fastmcp(name: str) -> object:
        assert name == "fastmcp"
        raise ModuleNotFoundError(name)

    monkeypatch.setattr("dnakit.tools.server.importlib.import_module", missing_fastmcp)
    with pytest.raises(BackendUnavailableError) as caught:
        server_module.create_server(registry)

    assert caught.value.code == "DNAKIT_MCP_UNAVAILABLE"
    assert "dnakit[agent]" in str(caught.value)


def test_server_entrypoint_disables_the_networked_banner_check(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = _FakeMCPServer()
    monkeypatch.setattr(server_module, "create_server", lambda: fake)

    server_module.main()

    assert fake.ran is True
    assert fake.show_banner is False
