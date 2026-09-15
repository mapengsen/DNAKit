"""Tests for compact Agent discovery and MCP-safe DNAKit execution."""

from __future__ import annotations

import importlib
import inspect
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, cast

import pytest

import dnakit.tools.server as server_module
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import PropertyPredictionConfig
from dnakit.tools import AgentToolRegistry, EmbeddingConfig, default_tool_registry
from dnakit.tools.adapters import coerce_value
from dnakit.tools.catalog import PUBLIC_TOOL_MODULES


def test_agent_reverse_complement_preserves_file_path(tmp_path: Path) -> None:
    source = tmp_path / "input.fna"
    target = tmp_path / "output.fna"
    source.write_text(">chr1\nAACGTN\n>chr2\nCCAA\n", encoding="utf-8")
    registry = default_tool_registry()
    registry.execute(
        "dnakit.ops.reverse_complement",
        {"value": str(source), "output_path": str(target), "workers": 1},
    )
    from dnakit.io import read

    records = list(read(target))
    assert [(record.id, record.sequence.to_string()) for record in records] == [
        ("chr1", "NACGTT"),
        ("chr2", "TTGG"),
    ]


def test_agent_split_accepts_multiple_files_and_writer_config(tmp_path: Path) -> None:
    import json

    sources = [tmp_path / "original.fna", tmp_path / "reverse.fna"]
    for source, sequence in zip(sources, ["AACGT", "ACGTT"], strict=True):
        source.write_text("".join(f">seq{i}\n{sequence}\n" for i in range(5)))
    target = tmp_path / "splits"
    registry = default_tool_registry()
    registry.execute(
        "dnakit.datasets.split",
        {
            "records": [str(source) for source in sources],
            "config": {
                "method": "random",
                "ratios": {"train": 0.8, "valid": 0.1, "test": 0.1},
                "seed": 42,
            },
            "output_path": str(target),
            "write_config": {"overwrite": False},
            "workers": 1,
        },
    )
    rows = {
        name: [json.loads(line) for line in (target / f"{name}.jsonl").read_text().splitlines()]
        for name in ["train", "valid", "test"]
    }
    assert {name: len(values) for name, values in rows.items()} == {
        "train": 8,
        "valid": 1,
        "test": 1,
    }
    assert sorted(row["sequence"] for values in rows.values() for row in values) == sorted(
        ["AACGT"] * 5 + ["ACGTT"] * 5
    )


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
    assert stats["categories"] == len({tool.category for tool in registry.tools})


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


def test_input_schema_retains_live_function_instructions() -> None:
    registry = default_tool_registry()
    for name in (
        "dnakit.references.download_genome",
        "dnakit.references.resolve_genome_assembly",
        "dnakit.ops.reverse_complement",
    ):
        module_name, function_name = name.rsplit(".", 1)
        function = getattr(importlib.import_module(module_name), function_name)
        tool = registry.resolve(name)
        assert tool.function is function
        docstring = inspect.getdoc(function)
        assert docstring is not None
        assert tool.call_plan.input_schema["description"] == docstring[:8_000]

    schema = cast(dict[str, Any], registry.describe("references.download_genome")["input_schema"])
    assert 'query="hg38"' in schema["description"]
    assert "not an LLM key" in schema["description"]
    assert schema["properties"]["query"]["type"] == "string"


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


def test_wrapped_metadata_filter_accepts_json_records_and_criteria() -> None:
    result = default_tool_registry().execute(
        "dnakit.io.filter_by_metadata",
        {
            "dataset": [
                {"id": "keep", "sequence": "ACGT", "metadata": {"group": "selected"}},
                {"id": "drop", "sequence": "TGCA", "metadata": {"group": "other"}},
            ],
            "criteria": {"group": "selected"},
        },
    )
    assert isinstance(result, dict)
    records = cast(list[dict[str, object]], result["records"])
    assert [record["id"] for record in records] == ["keep"]


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
        "search_dnakit_skills",
        "read_dnakit_skill",
    }
    response = fake.tools["call_dnakit_tool"](
        "dnakit.thermodynamics.molecular_weight",
        {"sequence": "ACGT"},
    )
    assert isinstance(response, dict)
    assert response["ok"] is True

    skills = fake.tools["search_dnakit_skills"]("变异注释", limit=3)
    assert isinstance(skills, dict)
    assert skills["ok"] is True
    assert len(skills["skills"]) == 3
    workflow = fake.tools["read_dnakit_skill"]("dnakit", limit=100)
    assert isinstance(workflow, dict)
    assert workflow["ok"] is True
    assert workflow["content"].startswith("---\nname: dnakit\n")
    assert workflow["next_offset"] == 100


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


def test_mcp_search_modes_match_python_and_do_not_execute(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = _FakeMCPServer()
    monkeypatch.setattr(server_module, "_fastmcp_factory", lambda: lambda **kw: fake)
    name = "dnakit.references.download_genome"
    registry = AgentToolRegistry(
        (default_tool_registry(verbose=False).resolve(name),),
        embedding_config=EmbeddingConfig(
            model="offline-model",
            base_url="https://embedding.example/v1",
            api_key="offline-key",
            cache_dir=None,
        ),
    )
    monkeypatch.setattr(
        "dnakit.tools.embeddings._embed", lambda texts, config: [(1.0, 0.0)] * len(texts)
    )
    server_module.create_server(registry, verbose=False)
    for mode in ("keyword", "semantic", "hybrid"):
        result = cast(
            dict[str, Any], fake.tools["search_dnakit_tools"]("download genome", mode=mode)
        )
        assert result["ok"] is True
        assert result["tools"] == [
            tool.to_summary() for tool in registry.search("download genome", mode=mode)
        ]
    result = cast(dict[str, Any], fake.tools["search_dnakit_tools"]("download", mode="bad"))
    assert result["error"]["code"] == "AGENT_SEARCH_MODE_INVALID"
    server_module.create_server(AgentToolRegistry(), verbose=False)
    result = cast(dict[str, Any], fake.tools["search_dnakit_tools"]("download", mode="semantic"))
    assert result["error"]["code"] == "MISSING_EMBEDDING_CONFIG"


def test_mcp_entrypoint_loads_embedding_settings_from_environment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fake = _FakeMCPServer()
    registries: list[AgentToolRegistry] = []

    def create(registry: AgentToolRegistry) -> _FakeMCPServer:
        registries.append(registry)
        return fake

    monkeypatch.setattr(server_module, "create_server", create)
    monkeypatch.setenv("DNAKIT_EMBEDDING_MODEL", "offline-model")
    monkeypatch.setenv("DNAKIT_EMBEDDING_BASE_URL", "https://embedding.example/v1")
    monkeypatch.setenv("DNAKIT_EMBEDDING_CACHE_DIR", str(tmp_path / "vectors"))
    server_module.main(verbose=False)
    assert fake.ran is True
    embedding_config = registries[0].embedding_config
    assert embedding_config is not None
    assert embedding_config.model == "offline-model"
    assert not (tmp_path / "vectors").exists()
