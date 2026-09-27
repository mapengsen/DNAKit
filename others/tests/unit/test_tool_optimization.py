"""Offline description optimization with real callable execution and file activation."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, cast
from unittest.mock import Mock

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools.catalog import _build_spec
from dnakit.tools.descriptions import (
    apply_native_description,
    description_path,
    read_descriptions,
    schema_contract,
)
from dnakit.tools.development_support import DevelopmentJob
from dnakit.tools.optimization import optimize_tool
from dnakit.tools.registry import AgentToolRegistry

NAME = "dnakit.ops.optimization_fixture"


@pytest.fixture
def fixture_registry() -> tuple[AgentToolRegistry, list[int]]:
    observed: list[int] = []

    def increment(value: int = 1) -> int:
        """Apply an operation."""
        observed.append(value)
        return value + 1

    return AgentToolRegistry(
        (_build_spec("dnakit.ops", "optimization_fixture", increment),)
    ), observed


def _action() -> dict[str, Any]:
    return {
        "action": "optimize_tool",
        "target": {"kind": "native", "name": NAME},
        "reason": "Clarify the required integer argument without changing the operation.",
        "changes": {
            "description": "Add one to the supplied integer.",
            "parameters": {"value": "The integer to increment."},
        },
        "scenarios": [
            {"request": "Ambiguous request", "expected_arguments": {"value": 3}, "expected": 4},
            {"request": "Existing request", "expected_arguments": {"value": 7}, "expected": 8},
        ],
    }


def _planner(*, improve: bool = True, regression: bool = False) -> Mock:
    def ask(prompt: str, data: dict[str, Any]) -> dict[str, Any]:
        if data["stage"] == "review":
            return {"approved": True, "semantic_preserved": True, "issues": []}
        assert "expected" not in data and "expected_arguments" not in data
        old = data["request"] == "Existing request"
        value = 7 if old else 3
        if not old and data["variant"] == "baseline" and improve:
            value = 2
        if old and data["variant"] == "candidate" and regression:
            value = 8
        return {"tool_name": data["tools"][0]["name"], "arguments": {"value": value}}

    return Mock(side_effect=ask)


def test_native_improvement_executes_correct_plans_and_preserves_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, calls = fixture_registry
    before = deepcopy(registry.describe(NAME))
    planner = _planner()
    result = optimize_tool(_action(), root=tmp_path, registry=registry, ask=planner, verbose=False)
    assert result["status"] == "activated"
    assert calls == [7, 3, 7]
    assert result["validation"]["baseline"][0]["executed"] is False
    monkeypatch.setattr("dnakit.tools.descriptions._DESCRIPTION_PATH", description_path(tmp_path))
    changed = apply_native_description(registry.resolve(NAME))
    changed_schema = cast(dict[str, Any], changed.call_plan.input_schema)
    assert changed.description == "Add one to the supplied integer."
    assert changed_schema["properties"]["value"]["description"] == "The integer to increment."
    assert schema_contract(changed.call_plan.input_schema) == schema_contract(
        before["input_schema"]
    )
    changed_schema["properties"]["value"]["description"] = "Local mutation"
    assert registry.describe(NAME) == before
    replay = optimize_tool(
        _action(),
        root=tmp_path,
        registry=registry,
        ask=Mock(side_effect=AssertionError),
        verbose=False,
    )
    assert replay["change_id"] == result["change_id"] and calls == [7, 3, 7]


@pytest.mark.parametrize("improve,regression", [(False, False), (True, True)])
def test_self_approval_cannot_replace_improvement_or_allow_regression(
    tmp_path: Path,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
    improve: bool,
    regression: bool,
) -> None:
    registry, _ = fixture_registry
    result = optimize_tool(
        _action(),
        root=tmp_path,
        registry=registry,
        ask=_planner(improve=improve, regression=regression),
        verbose=False,
    )
    assert result["status"] == "rejected"
    assert result["regression"] is regression
    assert not description_path(tmp_path).exists()


def test_schema_change_is_rejected_before_review_or_execution(
    tmp_path: Path,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, calls = fixture_registry
    action = _action()
    action["changes"] = {"parameters": {"value": {"type": "string"}}}
    ask = Mock(side_effect=AssertionError)
    result = optimize_tool(action, root=tmp_path, registry=registry, ask=ask, verbose=False)
    assert result["status"] == "rejected" and calls == []
    assert result["error"]["code"] == "TOOL_DESCRIPTION_CONTRACT_CHANGE"
    ask.assert_not_called()


def test_version_merge_and_repeated_application_keep_examples_and_schema_text_consistent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, _ = fixture_registry
    monkeypatch.setattr("dnakit.tools.descriptions._DESCRIPTION_PATH", description_path(tmp_path))
    initial = _action()
    initial["changes"].update(
        full_description="Add one to the supplied integer.\n\nOnly integers are supported.",
        examples=[{"arguments": {"value": 3}, "result": 4}],
    )
    result = optimize_tool(initial, root=tmp_path, registry=registry, ask=_planner(), verbose=False)
    assert result["status"] == "activated"
    once = apply_native_description(registry.resolve(NAME))
    twice = apply_native_description(once)
    assert once.to_dict() == twice.to_dict()
    assert twice.full_description.count("DNAKit usage examples:") == 1
    assert twice.call_plan.input_schema["description"] == twice.full_description
    changed = _action()
    changed["reason"] = "Refine the summary while preserving the previous detailed instructions."
    changed["changes"] = {"description": "Return the supplied integer plus one."}
    revised_registry = AgentToolRegistry((once,))
    result = optimize_tool(
        changed, root=tmp_path, registry=revised_registry, ask=_planner(), verbose=False
    )
    assert result["status"] == "activated"
    current = apply_native_description(registry.resolve(NAME))
    assert current.description == "Return the supplied integer plus one."
    assert current.full_description.startswith(current.description)
    assert "Only integers are supported." in current.full_description
    assert current.full_description.count("DNAKit usage examples:") == 1
    assert len(read_descriptions(description_path(tmp_path))["native"][NAME]["versions"]) == 2


def test_parameter_only_revision_changes_search_text_and_expires_on_effect_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, _ = fixture_registry
    original = registry.resolve(NAME)
    action = _action()
    action["changes"] = {"parameters": {"value": "The integer to increment, including negatives."}}
    result = optimize_tool(action, root=tmp_path, registry=registry, ask=_planner(), verbose=False)
    assert result["status"] == "activated"
    monkeypatch.setattr("dnakit.tools.descriptions._DESCRIPTION_PATH", description_path(tmp_path))
    changed = apply_native_description(original)
    assert changed.description == original.description
    assert changed.full_description == original.full_description
    assert original.search_text is None
    assert changed.search_text is not None and "including negatives" in changed.search_text
    assert apply_native_description(changed).search_text == changed.search_text
    changed_contract = replace(original, effect="network", requires_confirmation=True)
    assert apply_native_description(changed_contract) is changed_contract


def test_mcp_overlay_changes_listing_only_and_keeps_live_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    import dnakit.integrations.mcp as mcp

    registry, _ = fixture_registry
    original: dict[str, Any] = {
        "name": "increment",
        "description": "Apply an operation.",
        "inputSchema": {
            "type": "object",
            "properties": {"value": {"type": "integer", "default": 1}},
            "required": ["value"],
            "additionalProperties": False,
        },
    }
    observed: list[dict[str, Any]] = []

    def request(
        server_id: str, *, tool_name: str | None, arguments: Any, **kwargs: Any
    ) -> dict[str, Any]:
        if tool_name is None:
            return {"server_id": server_id, "tools": [deepcopy(original)]}
        mcp._validate_arguments(arguments, original["inputSchema"])
        observed.append(deepcopy(arguments))
        return {"isError": False, "structuredContent": arguments["value"] + 1}

    monkeypatch.setattr(mcp, "_request", request)
    monkeypatch.setattr("dnakit.tools.descriptions._DESCRIPTION_PATH", description_path(tmp_path))
    action = _action()
    action["target"] = {"kind": "mcp", "server_id": "fixture", "name": "increment"}
    result = optimize_tool(action, root=tmp_path, registry=registry, ask=_planner(), verbose=False)
    assert result["status"] == "activated" and len(observed) == 3
    listing = mcp.list_mcp_tools("fixture", verbose=False)
    assert listing["tools"][0]["description"] == "Add one to the supplied integer."
    assert original["description"] == "Apply an operation."
    with pytest.raises(ConfigurationError):
        mcp.call_mcp_tool("fixture", "increment", {"value": "wrong"}, verbose=False)
    assert len(observed) == 3


def _skill(root: Path) -> tuple[dict[str, Any], str]:
    directory = root / "skills" / "fixture-skill"
    directory.mkdir(parents=True)
    original = (
        "---\nname: fixture-skill\ndescription: Apply an operation.\nlicense: MIT\n---\n"
        "Use a threshold of 0.05.\n```bash\npython sample.py --limit 5\n```\n"
        "Source: https://example.org/source\n"
    )
    (directory / "SKILL.md").write_text(original)
    entry = {
        "id": "fixture-skill",
        "description": "Apply an operation.",
        "summary_zh": "处理输入整数",
        "license": "MIT",
        "upstream_url": "https://example.org/source",
        "attribution": "Author",
        "upstream_revision": "revision",
        "files": [
            {
                "path": "SKILL.md",
                "sha256": hashlib.sha256(original.encode()).hexdigest(),
                "bytes": len(original.encode()),
            }
        ],
    }
    (root / "registry").mkdir()
    (root / "registry/catalog.json").write_text(
        json.dumps({"schema_version": 1, "skills": [entry]})
    )
    candidate = original.replace(
        "description: Apply an operation.", "description: Increment the supplied integer."
    )
    action = _action()
    action["target"] = {"kind": "skill", "name": "fixture-skill"}
    action["changes"] = {"text": candidate, "summary_zh": "给输入整数加一"}
    for scenario in action["scenarios"]:
        scenario["tool_name"] = NAME
    return action, original


def test_skill_activation_updates_hash_and_summary_without_losing_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    import dnakit.skills.catalog as catalog

    action, original = _skill(tmp_path)
    registry, calls = fixture_registry
    monkeypatch.setattr(
        catalog, "_locations", lambda: (tmp_path / "skills", tmp_path / "registry/catalog.json")
    )
    catalog._catalog.cache_clear()
    assert catalog.read_skill_file("fixture-skill")["content"] == original
    result = optimize_tool(action, root=tmp_path, registry=registry, ask=_planner(), verbose=False)
    assert result["status"] == "activated" and calls == [7, 3, 7]
    updated = catalog.read_skill_file("fixture-skill")
    assert updated["content"] == action["changes"]["text"]
    assert updated["skill"]["description"] == "Increment the supplied integer."
    assert updated["skill"]["license"] == "MIT"
    assert updated["attribution"] == "Author"
    catalog._catalog.cache_clear()


@pytest.mark.parametrize(
    "old,new", [("0.05", "0.10"), ("--limit 5", "--limit 9"), ("license: MIT", "license: BSD")]
)
def test_skill_threshold_command_or_license_edits_are_rejected(
    tmp_path: Path,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
    old: str,
    new: str,
) -> None:
    action, original = _skill(tmp_path)
    action["changes"]["text"] = action["changes"]["text"].replace(old, new)
    registry, calls = fixture_registry
    result = optimize_tool(
        action,
        root=tmp_path,
        registry=registry,
        ask=Mock(side_effect=AssertionError),
        verbose=False,
    )
    assert result["status"] == "rejected" and calls == []
    assert (tmp_path / "skills/fixture-skill/SKILL.md").read_text() == original


def test_interrupted_validation_is_not_silently_reexecuted(
    tmp_path: Path,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, calls = fixture_registry
    interrupted = Mock(
        side_effect=[
            {"approved": True, "semantic_preserved": True, "issues": []},
            KeyboardInterrupt(),
        ]
    )
    with pytest.raises(KeyboardInterrupt):
        optimize_tool(_action(), root=tmp_path, registry=registry, ask=interrupted, verbose=False)
    ask = Mock(side_effect=AssertionError)
    result = optimize_tool(_action(), root=tmp_path, registry=registry, ask=ask, verbose=False)
    assert result["status"] == "rejected"
    assert result["code"] == "AGENT_DEVELOPMENT_VALIDATION_UNCERTAIN" and calls == []
    ask.assert_not_called()


def test_validated_candidate_resumes_activation_without_repeating_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, calls = fixture_registry
    original_activate = DevelopmentJob.activate
    monkeypatch.setattr(DevelopmentJob, "activate", Mock(side_effect=KeyboardInterrupt))
    with pytest.raises(KeyboardInterrupt):
        optimize_tool(_action(), root=tmp_path, registry=registry, ask=_planner(), verbose=False)
    assert calls == [7, 3, 7]
    monkeypatch.setattr(DevelopmentJob, "activate", original_activate)
    ask = Mock(side_effect=AssertionError)
    result = optimize_tool(_action(), root=tmp_path, registry=registry, ask=ask, verbose=False)
    assert result["status"] == "activated" and calls == [7, 3, 7]
    ask.assert_not_called()


def test_concurrent_metadata_edit_is_preserved(
    tmp_path: Path,
    fixture_registry: tuple[AgentToolRegistry, list[int]],
) -> None:
    registry, _ = fixture_registry
    plan = _planner()
    path = description_path(tmp_path)
    other = {"schema_version": 1, "native": {}, "mcp": {}, "other_editor": True}

    def ask(prompt: str, data: dict[str, Any]) -> dict[str, Any]:
        if data["stage"] == "review":
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(other))
        return cast(dict[str, Any], plan(prompt, data))

    result = optimize_tool(_action(), root=tmp_path, registry=registry, ask=ask, verbose=False)
    assert result["status"] == "rejected"
    assert result["error"]["code"] == "AGENT_DEVELOPMENT_CONFLICT"
    assert read_descriptions(path) == other
