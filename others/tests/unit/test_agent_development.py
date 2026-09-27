"""Offline end-to-end development, registry refresh, budget and permission tests."""

from __future__ import annotations

import json
from collections.abc import Iterator
from copy import deepcopy
from pathlib import Path
from typing import Any, cast
from unittest.mock import Mock

import pytest

import dnakit.generated as generated
import dnakit.tools.agent as agent
import dnakit.tools.descriptions as descriptions
import dnakit.tools.development_support as support
import dnakit.tools.server as server
from dnakit.core._json import FrozenDict
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.llm import LLMConfig, LLMResult
from dnakit.tools import (
    AgentSerializationLimits,
    AgentToolRegistry,
    EmbeddingConfig,
    RAGConfig,
    default_tool_registry,
    rollback_tool_change,
    run_agent,
)
from dnakit.tools.catalog import _build_spec

CONFIG = LLMConfig(provider="deepseek", model="offline-test", api_key="offline")
NAME = "dnakit.generated.count_literal_motif"
CODE = '''from __future__ import annotations

def count_literal_motif(sequence: str, motif: str, *, verbose: bool = True) -> int:
    """Count literal, non-overlapping motif occurrences; matching is case sensitive."""
    if not motif:
        raise ValueError("motif cannot be empty")
    return sequence.count(motif)
'''


def _candidate() -> dict[str, Any]:
    return {
        "action": "generate_tool",
        "name": "count_literal_motif",
        "reason": "The selected registry cannot count literal non-overlapping motifs.",
        "contract": {
            "description": "Count non-overlapping literal motif occurrences, case sensitive.",
            "parameters": {
                "type": "object",
                "properties": {"sequence": {"type": "string"}, "motif": {"type": "string"}},
                "required": ["sequence", "motif"],
            },
            "returns": {"type": "integer"},
            "dependencies": [],
            "effects": ["local"],
            "success_criteria": ["Return the exact non-overlapping occurrence count."],
        },
        "code": CODE,
        "tests": [
            {"arguments": {"sequence": "CGCGAA", "motif": "CG"}, "expected": 2},
            {"arguments": {"sequence": "", "motif": "A"}, "expected": 0},
            {"arguments": {"sequence": "AAAA", "motif": "AA"}, "expected": 2},
        ],
    }


@pytest.fixture
def checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    directory = tmp_path / "source/src/dnakit/generated"
    directory.mkdir(parents=True)
    (directory / "manifest.json").write_text('{"version":1,"tools":{}}')
    root = directory.parents[2]
    original = generated._DIRECTORY
    monkeypatch.setattr(generated, "_DIRECTORY", directory)
    monkeypatch.setattr(agent, "_source_root", lambda: root)
    monkeypatch.setattr(support, "_source_root", lambda: root)
    monkeypatch.setattr(descriptions, "_DESCRIPTION_PATH", descriptions.description_path(root))
    generated._refresh()
    default_tool_registry.cache_clear()
    try:
        yield root
    finally:
        monkeypatch.setattr(generated, "_DIRECTORY", original)
        generated._refresh()
        default_tool_registry.cache_clear()


def _reply(value: dict[str, Any]) -> LLMResult:
    return LLMResult(
        text="json",
        model="offline",
        response_id=None,
        finish_reason="stop",
        output_format="json",
        _json_data=FrozenDict(value),
    )


def _model(monkeypatch: pytest.MonkeyPatch, actions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    pending = iter(actions)
    inputs: list[dict[str, Any]] = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        inputs.append(deepcopy(kwargs["input_data"]))
        return _reply(next(pending))

    monkeypatch.setattr(agent, "call_llm", reply)
    return inputs


def _finish(status: str = "blocked") -> dict[str, str]:
    return {"action": "finish", "status": status, "text": "Observed results."}


@pytest.mark.parametrize(
    "generate,optimize,action",
    [
        (False, False, "generate_tool"),
        (False, False, "optimize_tool"),
        (True, False, "optimize_tool"),
        (False, True, "generate_tool"),
    ],
)
def test_flags_are_independent_and_model_cannot_enable_them(
    checkout: Path, monkeypatch: pytest.MonkeyPatch, generate: bool, optimize: bool, action: str
) -> None:
    _model(
        monkeypatch,
        [
            {"action": action, "generate_tools": True, "optimize_tools": True},
            _finish(),
        ],
    )
    result = run_agent(
        "A task",
        config=CONFIG,
        registry=AgentToolRegistry(()),
        allow_side_effects=True,
        generate_tools=generate,
        optimize_tools=optimize,
        verbose=False,
    )
    assert result.steps[0]["observation"]["error"]["code"] == "AGENT_DEVELOPMENT_DISABLED"  # type: ignore[index]
    assert not (checkout / ".dnakit").exists()


@pytest.mark.parametrize("flag", ["generate_tools", "optimize_tools"])
def test_development_requires_caller_permission_and_boolean_flags(
    monkeypatch: pytest.MonkeyPatch, flag: str
) -> None:
    model = Mock()
    monkeypatch.setattr(agent, "call_llm", model)
    options: dict[str, Any] = {flag: True}
    with pytest.raises(ConfigurationError, match="AGENT_TOOL_CONFIRMATION_REQUIRED"):
        run_agent("task", config=CONFIG, **options, verbose=False)
    options[flag] = "true"
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFIG"):
        run_agent("task", config=CONFIG, **options, verbose=False)
    model.assert_not_called()


def test_no_generation_without_prior_search(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _model(monkeypatch, [_candidate(), _finish()])
    result = run_agent(
        "task",
        config=CONFIG,
        registry=AgentToolRegistry(()),
        allow_side_effects=True,
        generate_tools=True,
        verbose=False,
    )
    assert result.steps[0]["observation"]["error"]["code"] == "AGENT_DEVELOPMENT_SEARCH_REQUIRED"  # type: ignore[index]
    assert not (checkout / ".dnakit").exists()


class _MCP:
    def __init__(self, **kwargs: Any) -> None:
        self.functions: dict[str, Any] = {}

    def tool(self, *, name: str, **kwargs: Any) -> Any:
        def register(function: Any) -> Any:
            self.functions[name] = function
            return function

        return register


def test_generated_python_and_mcp_match_and_agent_immediately_uses_it(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    embeddings = EmbeddingConfig(model="unused", base_url="https://unused.invalid/v1", api_key="x")
    limits = AgentSerializationLimits(max_items=12)
    registry = AgentToolRegistry((), embedding_config=embeddings, serialization_limits=limits)
    refreshed: list[AgentToolRegistry] = []
    original = AgentToolRegistry._refreshed

    def refresh(self: AgentToolRegistry, additions: tuple[str, ...] = ()) -> AgentToolRegistry:
        result = original(self, additions)
        refreshed.append(result)
        return result

    monkeypatch.setattr(AgentToolRegistry, "_refreshed", refresh)
    inputs = _model(
        monkeypatch,
        [
            {"action": "search", "query": "literal motif"},
            _candidate(),
            {"approved": True, "issues": []},
            {"action": "describe", "tool_name": NAME},
            {
                "action": "call",
                "tool_name": NAME,
                "arguments": {"sequence": "CGCGAA", "motif": "CG"},
            },
            _finish("completed"),
        ],
    )
    result = run_agent(
        "Count literal CG in CGCGAA",
        config=CONFIG,
        registry=registry,
        generate_tools=True,
        allow_side_effects=True,
        max_steps=6,
        cache_dir=checkout / "cache",
        verbose=False,
    )
    assert result.status == "completed"
    assert result.steps[3]["observation"]["result"] == 2  # type: ignore[index]
    assert len(inputs) == 6
    assert [tool.name for tool in refreshed[-1].tools] == [NAME]
    assert refreshed[-1].embedding_config is embeddings
    assert refreshed[-1]._serialization_limits is limits
    function = cast(Any, generated).count_literal_motif
    assert function("CGCGAA", "CG") == 2
    monkeypatch.setattr(server, "_fastmcp_factory", lambda: _MCP)
    mcp = server.create_server(refreshed[-1], verbose=False)
    assert mcp.functions["describe_dnakit_tool"](NAME)["ok"]  # type: ignore[attr-defined]
    assert (
        mcp.functions["call_dnakit_tool"](  # type: ignore[attr-defined]
            NAME, {"sequence": "CGCGAA", "motif": "CG"}
        )["result"]
        == 2
    )

    # A fresh run with development disabled can reuse the now-installed tool's results.
    inputs = _model(monkeypatch, [_finish("completed")])
    resumed = run_agent(
        "Count literal CG in CGCGAA",
        config=CONFIG,
        registry=refreshed[-1],
        allow_side_effects=True,
        cache_dir=checkout / "cache",
        reuse_results=True,
        verbose=False,
    )
    assert resumed.reused_steps == 1
    assert len(inputs) == 1
    receipt = result.steps[1]["observation"]["result"]  # type: ignore[index]
    rollback_tool_change(receipt["change_id"], allow_side_effects=True, verbose=False)
    assert not hasattr(generated, "count_literal_motif")
    assert not original(refreshed[-1]).tools


def test_internal_review_uses_the_same_model_budget(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _model(monkeypatch, [{"action": "search", "query": "motif"}, _candidate()])
    with pytest.raises(BackendExecutionError, match="AGENT_STEP_LIMIT"):
        run_agent(
            "task",
            config=CONFIG,
            registry=AgentToolRegistry(()),
            allow_side_effects=True,
            generate_tools=True,
            max_steps=2,
            verbose=False,
        )
    assert len(inputs) == 2
    assert generated._entries() == {}
    assert list((checkout / ".dnakit/development").glob("*/candidate-1.json"))


def _absolute(value: int, *, verbose: bool = True) -> int:
    """Transform value."""
    return abs(value)


def test_optimization_evaluation_is_blind_and_retrieval_sees_new_description(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = _build_spec("dnakit.test", "absolute", _absolute)
    registry = AgentToolRegistry((tool,))
    decision = {
        "action": "optimize_tool",
        "target": {"kind": "native", "name": tool.name},
        "reason": "Transform does not explain the absolute value operation.",
        "changes": {"description": "Return the absolute value (magnitude) of an integer."},
        "scenarios": [
            {
                "request": "Find magnitude of minus three",
                "expected_arguments": {"value": -3},
                "expected": 3,
            }
        ],
    }
    actions: Iterator[dict[str, Any]] = iter(
        [
            {"action": "search", "query": "absolute"},
            decision,
            {"action": "search", "query": "magnitude"},
            {"action": "describe", "tool_name": tool.name},
            {"action": "call", "tool_name": tool.name, "arguments": {"value": -3}},
            _finish("completed"),
        ]
    )
    evaluations = []

    def respond(*args: Any, **kwargs: Any) -> LLMResult:
        data = kwargs["input_data"]
        if data.get("stage") == "planning":
            evaluations.append(data)
            assert set(data) == {"stage", "request", "tools", "skill"}
            assert "expected" not in json.dumps(data)
            text = data["tools"][0]["description"]
            return _reply(
                {"tool_name": tool.name, "arguments": {"value": -3 if "magnitude" in text else 3}}
            )
        if data.get("candidate", {}).get("stage") == "review":
            return _reply({"approved": True, "semantic_preserved": True, "issues": []})
        return _reply(next(actions))

    monkeypatch.setattr(agent, "call_llm", respond)
    result = run_agent(
        "Find magnitude of minus three",
        config=CONFIG,
        registry=registry,
        optimize_tools=True,
        allow_side_effects=True,
        max_steps=9,
        verbose=False,
    )
    assert result.status == "completed"
    assert len(evaluations) == 2
    assert result.steps[1]["observation"]["ok"] is True  # type: ignore[index]
    assert result.steps[2]["observation"]["tools"][0]["name"] == tool.name  # type: ignore[index]
    assert result.steps[-1]["observation"]["result"] == 3  # type: ignore[index]
    assert _absolute(-3) == 3


def test_development_can_share_local_rag_and_checkpoints(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    docs = checkout / "knowledge"
    docs.mkdir()
    (docs / "motifs.md").write_text("Literal motifs use case sensitive non-overlapping matching.")
    inputs = _model(monkeypatch, [_finish("completed")])
    result = run_agent(
        "How are literal motifs matched?",
        config=CONFIG,
        registry=AgentToolRegistry(()),
        rag=RAGConfig(directory=docs),
        cache_dir=checkout / "cache",
        generate_tools=True,
        optimize_tools=True,
        allow_side_effects=True,
        verbose=False,
    )
    assert result.knowledge_context
    assert inputs[0]["knowledge_context"]
    assert not (checkout / ".dnakit").exists()


def test_internal_model_error_propagates_and_preserves_pending_candidate(
    checkout: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    replies: Iterator[LLMResult | ConfigurationError] = iter(
        [
            _reply({"action": "search", "query": "motif"}),
            _reply(_candidate()),
            ConfigurationError("secret credential in transport message", code="MODEL_DOWN"),
        ]
    )

    def respond(*args: Any, **kwargs: Any) -> LLMResult:
        value = next(replies)
        if isinstance(value, Exception):
            raise value
        return value

    monkeypatch.setattr(agent, "call_llm", respond)
    with pytest.raises(ConfigurationError, match="MODEL_DOWN"):
        run_agent(
            "task",
            config=CONFIG,
            registry=AgentToolRegistry(()),
            allow_side_effects=True,
            generate_tools=True,
            cache_dir=checkout / "cache",
            verbose=False,
        )
    checkpoint = next((checkout / "cache").glob("agent_runs/*/*/checkpoint.json"))
    saved = checkpoint.read_text()
    payload = json.loads(saved)
    assert payload["status"] == "interrupted"
    assert payload["pending"]["action"]["action"] == "generate_tool"
    assert "secret credential" not in saved
    assert not generated._entries()


def test_description_activation_keeps_unrelated_successful_call_cache() -> None:
    affected = json.dumps(["dnakit.ops.first", {"value": 1}])
    unrelated = json.dumps(["dnakit.ops.second", {"value": 2}])
    calls: dict[str, dict[str, object]] = {
        affected: {"ok": True, "result": 1},
        unrelated: {"ok": True, "result": 2},
    }
    names = agent._invalidate_developed_calls(
        calls,
        {
            "action": "optimize_tool",
            "target": {"kind": "native", "name": "dnakit.ops.first"},
        },
    )
    assert names == {"dnakit.ops.first"}
    assert set(calls) == {unrelated}


def test_missing_generated_dependency_does_not_disable_other_tools(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        generated,
        "_entries",
        lambda: {
            "optional_absolute": {
                "effect": "local",
                "description": "Optional absolute value",
                "contract": {"dependencies": ["missing_dnakit_optional_dependency"]},
            }
        },
    )
    optional = _build_spec("dnakit.generated", "optional_absolute", _absolute)
    native = _build_spec("dnakit.test", "absolute", _absolute)
    registry = AgentToolRegistry((optional, native))
    assert not optional.agent_compatible
    assert "missing_dnakit_optional_dependency" in optional.incompatibilities[0]
    assert registry.search("absolute") == (native,)
    assert registry.execute(native.name, {"value": -3}) == 3
