"""Agent repair permissions, fresh source evidence, runtime refresh and recovery."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
from collections.abc import Callable
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, cast
from unittest.mock import Mock

import pytest

import dnakit.ops as ops
import dnakit.tools.agent as agent
import dnakit.tools.registry as registries
import dnakit.tools.repair as repair
import dnakit.tools.server as server
from dnakit.core._json import FrozenDict
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.llm import LLMConfig, LLMResult
from dnakit.tools.catalog import AgentToolSpec, _build_spec, build_tool_catalog
from dnakit.tools.checkpoint import _AgentCheckpoint
from dnakit.tools.development_support import DevelopmentJob
from dnakit.tools.embeddings import EmbeddingConfig
from dnakit.tools.registry import AgentToolRegistry
from dnakit.tools.serialization import AgentSerializationLimits

_NAME = "dnakit.ops.increment_fixture"
_CONFIG = LLMConfig(provider="deepseek", model="offline", api_key="offline")
_CODE = '''from __future__ import annotations

def increment_fixture(value: int, *, verbose: bool = True) -> int:
    """Increment the supplied integer."""
    return value - 1
'''


@pytest.fixture
def checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, AgentToolRegistry]:
    root = tmp_path / "checkout"
    source = root / "src/dnakit/ops/repair_fixture.py"
    source.parent.mkdir(parents=True)
    source.write_text(_CODE)
    module_name = "dnakit.ops._agent_repair_fixture"
    spec = importlib.util.spec_from_file_location(module_name, source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, module_name, module)
    spec.loader.exec_module(module)
    monkeypatch.setattr(ops, "increment_fixture", module.increment_fixture, raising=False)
    monkeypatch.setattr(agent, "_source_root", lambda: root)

    def catalog(*, verbose: bool = True) -> tuple[AgentToolSpec, ...]:
        return (_build_spec("dnakit.ops", "increment_fixture", cast(Any, ops).increment_fixture),)

    monkeypatch.setattr(registries, "build_tool_catalog", catalog)
    return root, source, AgentToolRegistry(catalog())


def _model(
    monkeypatch: pytest.MonkeyPatch,
    actions: list[dict[str, Any] | Callable[[], dict[str, Any]]],
) -> list[dict[str, Any]]:
    pending = iter(actions)
    inputs = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        inputs.append(deepcopy(kwargs["input_data"]))
        action = next(pending)
        value = action() if callable(action) else action
        return LLMResult(
            text="json",
            model="offline",
            response_id=None,
            finish_reason="stop",
            output_format="json",
            _json_data=FrozenDict(value),
        )

    monkeypatch.setattr(agent, "call_llm", reply)
    return inputs


def _finish(status: str = "blocked") -> dict[str, str]:
    return {"action": "finish", "status": status, "text": "Observed offline results."}


def _prepare() -> list[dict[str, Any]]:
    return [
        {"action": "search", "query": _NAME},
        {"action": "describe", "tool_name": _NAME},
        {"action": "inspect_tool_source", "tool_name": _NAME},
    ]


def _repair(source: Path) -> dict[str, Any]:
    return {
        "action": "optimize_tool",
        "mode": "repair",
        "target": {"kind": "native", "name": _NAME},
        "reason": "Increment returned a decrement.",
        "base_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "code": _CODE.split("\n\n", 1)[1].replace("value - 1", "value + 1"),
        "tests": [
            {"arguments": {"value": 0}, "expected": 1},
            {"arguments": {"value": 4}, "expected": 5},
        ],
    }


def _run(registry: AgentToolRegistry, **kwargs: Any) -> agent.AgentResult:
    options: dict[str, Any] = {
        "config": _CONFIG,
        "registry": registry,
        "allow_side_effects": True,
        "optimize_tools": True,
        "optimization_mode": "repair",
        "verbose": False,
    }
    options.update(kwargs)
    return agent.run_agent("Fix increment and verify it", **options)


def _observation(result: agent.AgentResult, index: int) -> dict[str, Any]:
    return cast(dict[str, Any], result.steps[index]["observation"])


def _activate_stub(source: Path, *, ask: bool = False) -> Callable[..., dict[str, Any]]:
    def activate(action: dict[str, Any], **kwargs: Any) -> dict[str, Any]:
        if ask:
            kwargs["ask"](
                "Review the compatible repair", {"stage": "review", "code": action["code"]}
            )
        source.write_text(_CODE.replace("value - 1", "value + 1"))
        return {
            "status": "activated",
            "tool_name": _NAME,
            "source_path": str(source),
            "function_name": "increment_fixture",
            "artifacts": [str(source)],
        }

    return activate


@pytest.mark.parametrize("mode,optimize", [("description_only", True), ("repair", False)])
def test_source_inspection_cannot_enable_repair(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
    optimize: bool,
) -> None:
    _, _, registry = checkout
    _model(monkeypatch, [_prepare()[2], _finish()])
    result = _run(registry, optimization_mode=mode, optimize_tools=optimize)
    assert _observation(result, 0)["error"]["code"] == "AGENT_REPAIR_DISABLED"


@pytest.mark.parametrize(
    "prepare,expected",
    [([], "AGENT_TOOL_NOT_DISCOVERED"), (_prepare()[:1], "AGENT_TOOL_NOT_DESCRIBED")],
)
def test_inspection_requires_search_and_description(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
    prepare: list[dict[str, Any]],
    expected: str,
) -> None:
    _, _, registry = checkout
    _model(monkeypatch, [*prepare, _prepare()[2], _finish()])
    result = _run(registry)
    assert _observation(result, len(prepare))["error"]["code"] == expected


def test_repair_requires_opt_in_valid_mode_and_source_inspection(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, source, registry = checkout
    model = Mock()
    monkeypatch.setattr(agent, "call_llm", model)
    invalid_options: list[tuple[dict[str, Any], str]] = [
        ({"allow_side_effects": False}, "AGENT_TOOL_CONFIRMATION_REQUIRED"),
        ({"optimization_mode": "unsafe"}, "AGENT_DEVELOPMENT_CONFIG"),
    ]
    for options, code in invalid_options:
        with pytest.raises(ConfigurationError, match=code):
            _run(registry, **options)
    model.assert_not_called()
    _model(monkeypatch, [*_prepare()[:2], _repair(source), _finish()])
    result = _run(registry)
    assert _observation(result, 2)["error"]["code"] == "AGENT_REPAIR_INSPECTION_REQUIRED"
    _model(monkeypatch, [_prepare()[0], _repair(source), _finish()])
    result = _run(registry, optimization_mode="description_only")
    assert _observation(result, 1)["error"]["code"] == "AGENT_REPAIR_DISABLED"


@pytest.mark.parametrize("change_file", [False, True])
def test_changed_source_or_invented_hash_never_reaches_repair(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
    change_file: bool,
) -> None:
    _, source, registry = checkout
    action = _repair(source)

    def stale() -> dict[str, Any]:
        if change_file:
            source.write_text(_CODE + "\n# Concurrent source edit.\n")
        else:
            action["base_sha256"] = "0" * 64
        return action

    repair_call = Mock(side_effect=AssertionError("must not execute"))
    monkeypatch.setattr(repair, "repair_tool", repair_call)
    _model(monkeypatch, [*_prepare(), stale, _finish()])
    result = _run(registry)
    assert _observation(result, 3)["error"]["code"] == "AGENT_REPAIR_SOURCE_CHANGED"
    repair_call.assert_not_called()


def test_repaired_implementation_is_used_with_filters_effects_and_embeddings_preserved(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, source, registry = checkout
    embeddings = EmbeddingConfig(
        model="unused", base_url="https://unused.invalid", api_key="unused"
    )
    limits = AgentSerializationLimits(max_items=17)
    selected = replace(registry.resolve(_NAME), effect="filesystem", requires_confirmation=True)
    registry = AgentToolRegistry(
        (selected,), embedding_config=embeddings, serialization_limits=limits
    )
    repair_call = Mock(side_effect=_activate_stub(source))
    monkeypatch.setattr(repair, "repair_tool", repair_call)
    refreshed = []
    original_refresh = AgentToolRegistry._refreshed

    def refresh(
        self: AgentToolRegistry, additions: tuple[str, ...] = (), **kwargs: Any
    ) -> AgentToolRegistry:
        result = original_refresh(self, additions, **kwargs)
        refreshed.append(result)
        return result

    monkeypatch.setattr(AgentToolRegistry, "_refreshed", refresh)
    action = _repair(source)
    action.pop("mode")  # Configured repair is the default for native targets.
    inputs = _model(
        monkeypatch,
        [
            *_prepare(),
            action,
            {"action": "call", "tool_name": _NAME, "arguments": {"value": 5}},
            _prepare()[1],
            {"action": "call", "tool_name": _NAME, "arguments": {"value": 5}},
            _finish("completed"),
        ],
    )
    result = _run(registry)
    assert result.status == "completed"
    assert _observation(result, 4)["error"]["code"] == "AGENT_TOOL_NOT_DESCRIBED"
    assert _observation(result, 6)["result"] == 6
    assert repair_call.call_args.args[0]["mode"] == "repair"
    assert inputs[0]["development"]["optimization_mode"] == "repair"
    assert [tool.name for tool in refreshed[-1].tools] == [_NAME]
    assert refreshed[-1].embedding_config is embeddings
    assert refreshed[-1]._serialization_limits is limits
    assert refreshed[-1].resolve(_NAME).effect == "filesystem"
    assert refreshed[-1].resolve(_NAME).requires_confirmation
    assert cast(Any, ops).increment_fixture(5) == 6
    assert registry.resolve(_NAME).function(5) == 4  # Caller-owned snapshots remain immutable.


def test_custom_registry_is_not_replaced_by_same_named_catalog_function(
    checkout: tuple[Path, Path, AgentToolRegistry],
) -> None:
    _, _, registry = checkout

    def custom(value: int, *, verbose: bool = True) -> int:
        """A caller-selected implementation."""
        return value * 10

    custom_spec = replace(registry.resolve(_NAME), function=custom)
    chosen = AgentToolRegistry((custom_spec,))
    assert chosen._refreshed().resolve(_NAME).function is custom
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_INVALID"):
        chosen._refreshed(replacements=("dnakit.ops.not_selected",))


def test_missing_catalog_replacement_requires_a_repaired_matching_public_export(
    checkout: tuple[Path, Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    _, _, registry = checkout
    monkeypatch.setattr(registries, "build_tool_catalog", lambda **kwargs: ())
    with pytest.raises(ConfigurationError, match="matching repaired export"):
        registry._refreshed(replacements=(_NAME,))

    def unrelated(value: int, *, verbose: bool = True) -> int:
        return value

    monkeypatch.setattr(ops, "increment_fixture", unrelated)
    with pytest.raises(ConfigurationError, match="matching repaired export"):
        registry._refreshed(replacements=(_NAME,))


def test_real_custom_native_repair_refreshes_python_agent_and_mcp_without_expanding_registry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "complete_checkout"
    installed = Path(__file__).resolve().parents[3]
    shutil.copytree(
        installed / "src",
        root / "src",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
    )
    shutil.copyfile(installed / "pyproject.toml", root / "pyproject.toml")
    module_name = "dnakit.agent_repair_contribution_fixture"
    tool_name = f"{module_name}.magnitude"
    source = root / "src/dnakit/agent_repair_contribution_fixture.py"
    original_function = '''def magnitude(value: int, *, verbose: bool = True) -> int:
    """Return magnitude, rejecting zero according to the fixture contract."""
    if value == 0:
        raise ValueError("Zero is excluded")
    return value
'''
    source.write_text(
        '"""Independent custom-native repair fixture."""\n\n'
        "from __future__ import annotations\n\n"
        + original_function
        + '\n\ndef untouched() -> int:\n    return 17\n\n__all__ = ["magnitude", "untouched"]\n'
    )
    specification = importlib.util.spec_from_file_location(module_name, source)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    monkeypatch.setitem(sys.modules, module_name, module)
    specification.loader.exec_module(module)
    monkeypatch.setattr(agent, "_source_root", lambda: root)
    embeddings = EmbeddingConfig(
        model="unused", base_url="https://unused.invalid", api_key="unused"
    )
    limits = AgentSerializationLimits(max_items=23)
    registry = AgentToolRegistry(
        (_build_spec(module_name, "magnitude", module.magnitude),),
        embedding_config=embeddings,
        serialization_limits=limits,
    )
    assert all(tool.name != tool_name for tool in build_tool_catalog(verbose=False))
    refreshed = []
    original_refresh = AgentToolRegistry._refreshed

    def refresh(
        self: AgentToolRegistry, additions: tuple[str, ...] = (), **kwargs: Any
    ) -> AgentToolRegistry:
        result = original_refresh(self, additions, **kwargs)
        refreshed.append(result)
        return result

    monkeypatch.setattr(AgentToolRegistry, "_refreshed", refresh)
    fixed = original_function.replace(
        "verbose: bool = True", "verbose: bool = True, offset: int = 0"
    ).replace("return value", "return abs(value) + offset")
    action = {
        "action": "optimize_tool",
        "mode": "repair",
        "target": {"kind": "native", "name": tool_name},
        "reason": "Fix negative magnitudes and add an optional offset while preserving old calls.",
        "base_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "code": fixed,
        "tests": [
            {"arguments": {"value": 3}, "expected": 3},
            {"arguments": {"value": -3}, "expected": 3},
            {"arguments": {"value": 3, "offset": 2}, "expected": 5},
            {"arguments": {"value": 0}, "expected_error": {"type": "ValueError"}},
        ],
    }
    inputs = _model(
        monkeypatch,
        [
            {"action": "search", "query": tool_name},
            {"action": "describe", "tool_name": tool_name},
            {"action": "inspect_tool_source", "tool_name": tool_name},
            action,
            {"approved": True, "issues": []},
            {"action": "describe", "tool_name": tool_name},
            {"action": "call", "tool_name": tool_name, "arguments": {"value": -3, "offset": 2}},
            _finish("completed"),
        ],
    )
    result = _run(registry, max_steps=8)
    assert result.status == "completed", result.steps
    receipt = _observation(result, 3)["result"]
    assert receipt["status"] == "activated"
    assert receipt["validation"]["baseline"]["cases"][1]["observed"] == -3
    assert receipt["validation"]["candidate"]["cases"][1]["observed"] == 3
    assert receipt["validation"]["baseline"]["cases"][2]["status"] == "unsupported_optional"
    assert receipt["validation"]["candidate"]["cases"][2]["executed"]
    assert receipt["compatibility"]["added_optional_parameters"] == ["offset"]
    assert len(inputs) == 8 and "offset" in source.read_text()
    current = refreshed[-1]
    assert [tool.name for tool in current.tools] == [tool_name]
    assert current.embedding_config is embeddings and current._serialization_limits is limits
    assert current.resolve(tool_name).effect == registry.resolve(tool_name).effect
    assert (
        current.resolve(tool_name).requires_confirmation
        == registry.resolve(tool_name).requires_confirmation
    )
    description = _observation(result, 4)["tool"]
    assert description["input_schema"]["properties"]["offset"] == {"type": "integer", "default": 0}
    assert "offset" not in description["input_schema"]["required"]
    assert _observation(result, 5)["result"] == module.magnitude(-3, offset=2) == 5
    assert module.untouched() == 17

    class RecordingMCP:
        def __init__(self, **kwargs: Any) -> None:
            self.functions: dict[str, Any] = {}

        def tool(self, *, name: str, **kwargs: Any) -> Any:
            def register(function: Any) -> Any:
                self.functions[name] = function
                return function

            return register

    monkeypatch.setattr(server, "_fastmcp_factory", lambda: RecordingMCP)
    mcp = cast(Any, server.create_server(current, verbose=False))
    mcp_result = mcp.functions["call_dnakit_tool"](tool_name, {"value": -3, "offset": 2})
    assert mcp_result["ok"] and mcp_result["result"] == 5
    assert mcp.functions["describe_dnakit_tool"](tool_name)["tool"]["input_schema"]["properties"][
        "offset"
    ]


def test_repair_mode_keeps_explicit_text_optimization_and_rejects_remote_repair(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import dnakit.tools.optimization as optimization

    _, source, registry = checkout
    text_action = {
        "action": "optimize_tool",
        "mode": "description_only",
        "target": {"kind": "native", "name": _NAME},
    }
    remote_action = {**_repair(source), "target": {"kind": "mcp", "name": "remote"}}
    optimize = Mock(return_value={"status": "rejected", "reason": "No improvement."})
    monkeypatch.setattr(optimization, "optimize_tool", optimize)
    _model(monkeypatch, [_prepare()[0], text_action, remote_action, _finish()])
    result = _run(registry)
    optimize.assert_called_once()
    assert _observation(result, 2)["error"]["code"] == "AGENT_REPAIR_TARGET_UNSUPPORTED"
    assert source.read_text() == _CODE


def test_repair_reviews_use_the_shared_model_budget(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, source, registry = checkout
    monkeypatch.setattr(repair, "repair_tool", _activate_stub(source, ask=True))
    _model(monkeypatch, [*_prepare(), _repair(source)])
    with pytest.raises(BackendExecutionError, match="AGENT_STEP_LIMIT"):
        _run(registry, max_steps=4)
    assert source.read_text() == _CODE


def _checkpoint(
    root: Path, registry: AgentToolRegistry, *, reuse: bool = False
) -> _AgentCheckpoint:
    return _AgentCheckpoint(
        root / "cache",
        "Fix increment and verify it",
        registry=registry,
        code_execution=None,
        allow_side_effects=True,
        reuse_results=reuse,
        verbose=False,
        context=[],
    )


def test_source_inspection_checkpoint_rejects_stale_hash_and_malformed_evidence(
    checkout: tuple[Path, Path, AgentToolRegistry],
) -> None:
    root, source, registry = checkout
    checkpoint = _checkpoint(root, registry)
    inspection = repair.inspect_native_tool(_NAME, root=root, registry=registry)
    checkpoint.record(
        {"step": 1, "action": _prepare()[2], "observation": {"ok": True, "result": inspection}}
    )
    restored = _checkpoint(root, registry, reuse=True)
    assert restored.history == checkpoint.history
    source.write_text(_CODE + "\n# Changed since inspection.\n")
    changed = _checkpoint(root, registry, reuse=True)
    assert changed.invalidated_from == 1 and changed.history == []
    path = changed.directory / "checkpoint.json"
    payload = json.loads(path.read_text())
    payload["history"] = checkpoint.history
    payload["history"][0]["observation"]["result"] = "corrupted"
    payload["proofs"] = [{}]
    path.write_text(json.dumps(payload))
    with pytest.raises(ConfigurationError, match="AGENT_CHECKPOINT_INVALID"):
        _checkpoint(root, registry, reuse=True)


def test_interrupted_activated_repair_replays_receipt_without_old_source_or_execution(
    checkout: tuple[Path, Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root, source, registry = checkout
    action = _repair(source)
    checkpoint = _checkpoint(root, registry)
    checkpoint.record(
        {
            "step": 1,
            "action": _prepare()[0],
            "observation": {"ok": True, "tools": [registry.resolve(_NAME).to_summary()]},
        }
    )
    key = json.dumps(action, sort_keys=True, ensure_ascii=False, allow_nan=False)
    checkpoint.begin(action, cache_key=key)
    checkpoint.finish("interrupted")
    job = DevelopmentJob(root, "repair", _NAME, action, False)
    receipt = job.activate({source: _CODE.replace("value - 1", "value + 1").encode()})
    metadata = {
        "tool_name": _NAME,
        "source_path": str(source),
        "function_name": "increment_fixture",
    }
    job.state["repair_metadata"] = metadata
    job.complete({**receipt, **metadata})
    # repair_tool must use its durable receipt before checking the now-outdated base hash.
    model_inputs = _model(
        monkeypatch,
        [
            action,
            _prepare()[1],
            {"action": "call", "tool_name": _NAME, "arguments": {"value": 4}},
            _finish("completed"),
        ],
    )
    result = _run(registry, cache_dir=root / "cache", reuse_results=True)
    assert result.status == "completed" and _observation(result, 3)["result"] == 5
    assert model_inputs[0]["resume"]["uncertain_actions"]
    assert not model_inputs[-1]["resume"]["uncertain_actions"]
