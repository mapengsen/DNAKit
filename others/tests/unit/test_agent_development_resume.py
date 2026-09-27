"""Development receipts are reusable only while activated artifacts still match."""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools import AgentToolRegistry, default_tool_registry
from dnakit.tools.checkpoint import _AgentCheckpoint

WEIGHT = "dnakit.thermodynamics.molecular_weight"


def test_local_mcp_overlay_changes_invalidate_cached_schema_descriptions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    name = "dnakit.integrations.list_mcp_tools"
    registry = AgentToolRegistry((default_tool_registry(verbose=False).resolve(name),))
    metadata = {"schema_version": 1, "native": {}, "mcp": {}}
    monkeypatch.setattr("dnakit.tools.descriptions.read_descriptions", lambda: metadata)
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    checkpoint.record(
        {
            "step": 1,
            "action": {"action": "describe", "tool_name": name},
            "observation": {"ok": True, "tool": registry.describe(name)},
        }
    )
    checkpoint.finish("completed")
    metadata["mcp"] = {"upstream": {"tool": {"active": "new-description"}}}
    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)
    assert restored.invalidated_from == 1
    assert restored.history == []


@pytest.fixture
def registry() -> AgentToolRegistry:
    return AgentToolRegistry((default_tool_registry(verbose=False).resolve(WEIGHT),))


def _checkpoint(
    directory: Path,
    registry: AgentToolRegistry,
    *,
    reuse: bool = False,
    permission: bool = True,
) -> _AgentCheckpoint:
    return _AgentCheckpoint(
        directory,
        "Develop a tool and use the verified implementation",
        registry=registry,
        code_execution=None,
        allow_side_effects=permission,
        reuse_results=reuse,
        verbose=False,
    )


def _record_call(checkpoint: _AgentCheckpoint, registry: AgentToolRegistry, sequence: str) -> None:
    arguments = {"sequence": sequence}
    action = {"action": "call", "tool_name": WEIGHT, "arguments": arguments}
    key = json.dumps([WEIGHT, arguments], sort_keys=True)
    checkpoint.begin(action, inputs=arguments, cache_key=key)
    checkpoint.record(
        {
            "step": len(checkpoint.history) + 1,
            "action": action,
            "observation": {"ok": True, "result": registry.execute(WEIGHT, arguments)},
        },
        cache_key=key,
    )


def _record_development(
    checkpoint: _AgentCheckpoint,
    kind: str,
    artifacts: list[Path],
    *,
    inputs: dict[str, Any] | None = None,
) -> None:
    action = {"action": kind, "reason": "Observed capability or documentation gap"}
    key = json.dumps(action, sort_keys=True)
    checkpoint.begin(action, inputs=inputs, cache_key=key)
    checkpoint.record(
        {
            "step": len(checkpoint.history) + 1,
            "action": action,
            "observation": {
                "ok": True,
                "result": {
                    "status": "activated",
                    "change_id": "offline-change",
                    "audit_dir": str(checkpoint.directory / "development_audit"),
                    "artifacts": [str(path) for path in artifacts],
                    "tool_name": WEIGHT,
                },
            },
        },
        cache_key=key,
    )


@pytest.mark.parametrize("kind", ["generate_tool", "optimize_tool"])
def test_activated_history_remains_readable_without_development_permission(
    tmp_path: Path,
    registry: AgentToolRegistry,
    kind: str,
) -> None:
    artifact = tmp_path / "implementation.py"
    artifact.write_text("def generated():\n    return 1\n", encoding="utf-8")
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    _record_development(checkpoint, kind, [artifact])
    checkpoint.finish("interrupted")
    audit = checkpoint.directory / "development_audit"
    audit.mkdir()
    (audit / "review.txt").write_text("Additional review", encoding="utf-8")

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True, permission=False)

    assert restored.history == checkpoint.history
    assert restored.calls == {} and restored.reused_steps == 0
    assert restored.invalidated_from is None
    payload = json.loads((restored.directory / "checkpoint.json").read_text(encoding="utf-8"))
    assert set(payload["proofs"][0]["outputs"]) == {str(artifact)}
    assert payload["proofs"][0]["outputs"][str(artifact)]["sha256"]


@pytest.mark.parametrize("kind", ["generate_tool", "optimize_tool"])
@pytest.mark.parametrize("change", ["implementation", "manifest", "rollback"])
def test_artifact_edit_deletion_and_rollback_truncate_development_and_later_calls(
    tmp_path: Path,
    registry: AgentToolRegistry,
    kind: str,
    change: str,
) -> None:
    implementation = tmp_path / "implementation.py"
    manifest = tmp_path / "manifest.json"
    descriptions = tmp_path / "descriptions.json"
    implementation.write_text("def generated():\n    return 1\n", encoding="utf-8")
    manifest.write_text('{"name": "generated"}\n', encoding="utf-8")
    descriptions.write_text('{"description": "new wording"}\n', encoding="utf-8")
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    _record_call(checkpoint, registry, "ACGT")
    _record_development(checkpoint, kind, [implementation, manifest, descriptions])
    _record_call(checkpoint, registry, "TTTT")
    checkpoint.finish("completed")

    if change == "implementation":
        previous = implementation.stat()
        implementation.write_text("def generated():\n    return 2\n", encoding="utf-8")
        os.utime(implementation, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    elif change == "manifest":
        manifest.unlink()
    else:
        descriptions.write_text('{"description": "old wording"}\n', encoding="utf-8")

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)

    assert restored.invalidated_from == 2
    assert restored.history == checkpoint.history[:1]
    assert restored.reused_steps == 1 and len(restored.calls) == 1


@pytest.mark.parametrize("kind", ["generate_tool", "optimize_tool"])
def test_interrupted_development_remains_uncertain_across_restarts(
    tmp_path: Path,
    registry: AgentToolRegistry,
    kind: str,
) -> None:
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    action = {"action": kind, "reason": "Develop after observed gap"}
    key = json.dumps(action, sort_keys=True)
    checkpoint.begin(action, cache_key=key)
    (tmp_path / "partially_activated.py").write_text("# partial state\n", encoding="utf-8")
    checkpoint.finish("interrupted", "KeyboardInterrupt")

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)
    assert restored.history == [] and restored.calls == {} and restored.reused_steps == 0
    assert restored.pending is None
    assert len(restored.uncertain_actions) == 1
    assert restored.uncertain_actions[0]["cache_key"] == key
    assert restored.uncertain_actions[0]["action"] == action
    restored.finish("blocked", "Await development transaction recovery")

    again = _checkpoint(tmp_path / "cache", registry, reuse=True, permission=False)
    assert again.uncertain_actions == restored.uncertain_actions
    assert again.history == [] and again.reused_steps == 0


def test_missing_declared_artifact_cannot_be_reused(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    _record_development(checkpoint, "generate_tool", [tmp_path / "missing.py"])
    checkpoint.finish("completed")

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)

    assert restored.invalidated_from == 1 and restored.history == []
    assert restored.calls == {} and restored.reused_steps == 0


@pytest.mark.parametrize("artifact_exists", [False, True])
def test_reconciled_activation_clears_uncertainty_only_with_valid_artifacts(
    tmp_path: Path,
    registry: AgentToolRegistry,
    artifact_exists: bool,
) -> None:
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    action = {"action": "generate_tool", "reason": "Generate missing capability"}
    checkpoint.begin(action, cache_key="interrupted-job")
    checkpoint.finish("interrupted", "KeyboardInterrupt")
    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)
    artifact = tmp_path / "implementation.py"
    if artifact_exists:
        artifact.write_text("def generated():\n    return 1\n", encoding="utf-8")
    restored.begin(action, cache_key="interrupted-job")
    restored.record(
        {
            "step": 1,
            "action": action,
            "observation": {
                "ok": True,
                "result": {"status": "activated", "artifacts": [str(artifact)]},
            },
        },
        cache_key="interrupted-job",
    )

    assert bool(restored.uncertain_actions) is not artifact_exists
    again = _checkpoint(tmp_path / "cache", registry, reuse=True)
    assert bool(again.uncertain_actions) is not artifact_exists
    assert len(again.history) == int(artifact_exists)
    assert again.calls == {} and again.reused_steps == 0


def test_failed_development_remains_diagnostic_history_and_is_not_cached(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    action = {"action": "optimize_tool", "reason": "Correct an ambiguous description"}
    checkpoint.begin(action, cache_key="failed-development")
    checkpoint.record(
        {
            "step": 1,
            "action": action,
            "observation": {"ok": False, "error": {"code": "DEVELOPMENT_REVIEW_FAILED"}},
        },
        cache_key="failed-development",
    )

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True, permission=False)

    assert restored.history == checkpoint.history
    assert restored.calls == {} and restored.reused_steps == 0
    assert restored.pending is None and restored.uncertain_actions == []


def test_invalid_activation_result_preserves_pending_receipt(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    action = {"action": "generate_tool", "reason": "Generate a missing tool"}
    checkpoint.begin(action, cache_key="development-key")
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_RESULT_INVALID"):
        checkpoint.record(
            {
                "step": 1,
                "action": action,
                "observation": {"ok": True, "result": {"status": "activated", "artifacts": []}},
            },
            cache_key="development-key",
        )

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)

    assert restored.history == []
    assert restored.uncertain_actions[0]["cache_key"] == "development-key"


def test_optimization_validates_the_activated_version_of_an_in_place_input(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    artifact = tmp_path / "descriptions.json"
    artifact.write_text('{"description": "old wording"}\n', encoding="utf-8")
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    action = {"action": "optimize_tool", "reason": "Clarify a parameter"}
    checkpoint.begin(action, inputs={"source_path": str(artifact)}, cache_key="optimize")
    artifact.write_text('{"description": "new wording"}\n', encoding="utf-8")
    checkpoint.record(
        {
            "step": 1,
            "action": action,
            "observation": {
                "ok": True,
                "result": {"status": "activated", "artifacts": [str(artifact)]},
            },
        },
        cache_key="optimize",
    )

    restored = _checkpoint(tmp_path / "cache", registry, reuse=True)

    assert restored.history == checkpoint.history and restored.invalidated_from is None


def test_refreshed_registry_fingerprints_new_calls_without_rewriting_old_proofs(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    refreshed = AgentToolRegistry(
        (replace(registry.resolve(WEIGHT), description="Clarified molecular-weight description"),)
    )
    checkpoint = _checkpoint(tmp_path / "old_and_new", registry)
    _record_call(checkpoint, registry, "ACGT")
    checkpoint.update_registry(refreshed)
    _record_call(checkpoint, refreshed, "TTTT")
    payload = json.loads((checkpoint.directory / "checkpoint.json").read_text(encoding="utf-8"))
    assert payload["proofs"][0]["tool_fingerprint"] != payload["proofs"][1]["tool_fingerprint"]
    restored = _checkpoint(tmp_path / "old_and_new", refreshed, reuse=True)
    assert restored.invalidated_from == 1 and restored.history == []

    fresh = _checkpoint(tmp_path / "only_new", registry)
    fresh.update_registry(refreshed)
    _record_call(fresh, refreshed, "TTTT")
    restored = _checkpoint(tmp_path / "only_new", refreshed, reuse=True)
    assert restored.reused_steps == 1 and restored.history == fresh.history


def test_corrupt_development_proof_is_reported_as_a_checkpoint_error(
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    artifact = tmp_path / "implementation.py"
    artifact.write_text("def generated():\n    return 1\n", encoding="utf-8")
    checkpoint = _checkpoint(tmp_path / "cache", registry)
    _record_development(checkpoint, "generate_tool", [artifact])
    path = checkpoint.directory / "checkpoint.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["proofs"][0]["outputs"] = []
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ConfigurationError, match="AGENT_CHECKPOINT_INVALID"):
        _checkpoint(tmp_path / "cache", registry, reuse=True)
