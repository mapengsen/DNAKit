"""Offline, real-process regression coverage for compatible native tool repairs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from dnakit import generated
from dnakit.exceptions import ConfigurationError
from dnakit.tools import repair as repair_module
from dnakit.tools.catalog import _build_spec
from dnakit.tools.development_support import DevelopmentJob, rollback_tool_change
from dnakit.tools.registry import AgentToolRegistry
from dnakit.tools.repair import inspect_native_tool, repair_tool

NAME = "dnakit.repair_fixture.magnitude"
FUNCTION = '''def magnitude(value: int, *, verbose: bool = True) -> int:
    """Return the magnitude; zero is rejected by this fixture's contract."""
    if value == 0:
        raise ValueError("Zero is excluded")
    return value
'''
FIXED = FUNCTION.replace("    return value", "    return abs(value)")
MODULE = (
    "from __future__ import annotations\n\nUNRELATED = 17\n\n"
    + FUNCTION
    + """
def untouched() -> int:
    return UNRELATED

__all__ = ["magnitude", "untouched"]
"""
)


@pytest.fixture
def source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, AgentToolRegistry]:
    package = tmp_path / "src/dnakit"
    package.mkdir(parents=True)
    # A small real package supplies the tool while installed DNAKit supplies its
    # normal adapter dependencies. Full source checkouts copy all their modules.
    installed_package = Path(__file__).resolve().parents[3] / "src/dnakit"
    (package / "__init__.py").write_text(f"__path__.append({str(installed_package)!r})\n")
    path = package / "repair_fixture.py"
    path.write_text(MODULE)
    specification = importlib.util.spec_from_file_location("dnakit.repair_fixture", path)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    monkeypatch.setitem(sys.modules, specification.name, module)
    specification.loader.exec_module(module)
    registry = AgentToolRegistry(
        (_build_spec("dnakit.repair_fixture", "magnitude", module.magnitude),)
    )
    return path, registry


def _action(path: Path, *, code: str = FIXED, name: str = NAME) -> dict[str, Any]:
    return {
        "action": "optimize_tool",
        "mode": "repair",
        "target": {"kind": "native", "name": name},
        "reason": (
            "Negative magnitudes incorrectly remain negative; preserve positive and zero behavior."
        ),
        "base_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "code": code,
        "tests": [
            {"arguments": {"value": 4}, "expected": 4},
            {"arguments": {"value": -4}, "expected": 4},
            {"arguments": {"value": 0}, "expected_error": {"type": "ValueError"}},
        ],
    }


def _review() -> Mock:
    return Mock(return_value={"approved": True, "issues": []})


def _repair(
    root: Path, registry: AgentToolRegistry, action: dict[str, Any], reviewer: Mock | None = None
) -> dict[str, Any]:
    return repair_tool(
        action, root=root, registry=registry, ask=reviewer or _review(), verbose=False
    )


def test_inspect_returns_registered_source_and_exact_module_fingerprint(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    result = inspect_native_tool(NAME, root=tmp_path, registry=registry)
    assert result["source_path"] == str(path)
    assert result["relative_path"] == "src/dnakit/repair_fixture.py"
    assert result["source"] == FUNCTION
    assert result["module_source"] == MODULE
    assert result["base_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert result["function_name"] == "magnitude" and result["effect"] == "local"


def test_repair_proves_old_success_error_preservation_and_real_bug_fix(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    result = _repair(tmp_path, registry, _action(path))
    assert result["status"] == "activated", result
    assert path.read_text() == MODULE.replace(FUNCTION, FIXED)
    validation = result["validation"]
    assert validation["baseline"]["cases"][1]["observed"] == -4
    assert validation["candidate"]["cases"][1]["observed"] == 4
    assert validation["baseline"]["cases"][2]["matched"] is True
    assert validation["candidate"]["cases"][2]["error"]["type"] == "ValueError"
    assert validation["preserved_successes"] == [1]
    assert validation["improved_cases"] == [2]
    assert result["kind"] == "native"
    assert result["source_path"] == str(path) and result["code"] == FIXED
    assert len(result["artifacts"]) == 6  # source plus the reviewable contribution files
    # Runtime mutation belongs to the caller; this module never reloads live classes/functions.
    assert registry.execute(NAME, {"value": -4}) == -4


def test_optional_extension_is_classified_as_unsupported_without_fake_baseline_execution(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    code = FUNCTION.replace("verbose: bool = True", "verbose: bool = True, positive: bool = False")
    code = code.replace("return value", "return abs(value) if positive else value")
    action = _action(path, code=code)
    action["tests"][1]["arguments"]["positive"] = True
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "activated", result
    baseline = result["validation"]["baseline"]["cases"][1]
    assert baseline["status"] == "unsupported_optional" and baseline["executed"] is False
    assert result["validation"]["candidate"]["cases"][1]["executed"] is True
    assert result["compatibility"]["added_optional_parameters"] == ["positive"]


def test_new_optional_parameter_needs_a_successful_demonstration(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    code = FUNCTION.replace("verbose: bool = True", "verbose: bool = True, positive: bool = False")
    code = code.replace("return value", 'raise ValueError("Optional feature is unimplemented")')
    # Preserve the demonstrated old success while the new parameter only errors.
    code = code.replace(
        '    raise ValueError("Optional feature is unimplemented")',
        '    if positive:\n        raise ValueError("Optional feature is unimplemented")\n'
        "    return value",
    )
    action = _action(path, code=code)
    action["tests"][1] = {
        "arguments": {"value": -4, "positive": True},
        "expected_error": {"type": "ValueError"},
    }
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected"
    assert "demonstrated fix" in result["reason"]


@pytest.mark.parametrize(
    "code",
    [
        FIXED.replace("value: int", "value: str"),
        FIXED.replace("verbose: bool = True", "verbose: bool = False"),
        FIXED.replace(") -> int:", ") -> float:"),
        FIXED.replace("value: int, *", "value: int, extra: int, *"),
        FIXED + "\nUNRELATED = 99\n",
        FIXED.replace("def magnitude", "def renamed"),
    ],
)
def test_breaking_or_unrelated_edits_never_reach_model_review(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], code: str
) -> None:
    path, registry = source
    reviewer = _review()
    result = _repair(tmp_path, registry, _action(path, code=code), reviewer)
    assert result["status"] == "rejected"
    assert path.read_text() == MODULE
    reviewer.assert_not_called()


def test_new_network_scope_is_rejected_before_execution(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    code = FIXED.replace(
        "    return abs(value)", "    import urllib.request\n    return abs(value)"
    )
    reviewer = _review()
    result = _repair(tmp_path, registry, _action(path, code=code), reviewer)
    assert result["status"] == "rejected" and result["code"] == "TOOL_REPAIR_EFFECT_EXPANSION"
    reviewer.assert_not_called()


def test_successful_review_cannot_override_regression_or_missing_improvement(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    result = _repair(
        tmp_path, registry, _action(path, code=FIXED.replace("abs(value)", "abs(value) + 1"))
    )
    assert result["status"] == "rejected"
    assert result["validation"]["candidate"]["cases"][0]["matched"] is False
    assert path.read_text() == MODULE


def test_identical_behavior_has_no_demonstrated_improvement(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    action = _action(path, code=FUNCTION)
    action["tests"][1]["expected"] = -4
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected"
    assert "demonstrated fix" in result["reason"]


def test_missing_old_success_is_not_enough_for_activation(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    action = _action(path)
    action["tests"].pop(0)
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected"
    assert "successful old call" in result["reason"]


def test_stale_inspection_and_external_source_are_rejected(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    action = _action(path)
    path.write_text(MODULE.replace("UNRELATED = 17", "UNRELATED = 18"))
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected" and result["code"] == "TOOL_REPAIR_STALE_SOURCE"
    with pytest.raises(ConfigurationError, match="inside this checkout"):
        inspect_native_tool(NAME, root=tmp_path / "different_checkout", registry=registry)


def test_rejected_review_never_executes_baseline_or_candidate(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    monkeypatch.setattr(
        "dnakit.tools.repair._run_variant", Mock(side_effect=AssertionError("must not execute"))
    )
    reviewer = Mock(return_value={"approved": False, "issues": ["Inadequate evidence"]})
    result = _repair(tmp_path, registry, _action(path), reviewer)
    assert result["status"] == "rejected" and path.read_text() == MODULE


def test_activation_replay_never_reexecutes_or_requires_the_old_source_fingerprint(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    action = _action(path)
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "activated", result
    monkeypatch.setattr(
        "dnakit.tools.repair._run_variant", Mock(side_effect=AssertionError("must replay"))
    )
    replay = _repair(tmp_path, registry, action, Mock(side_effect=AssertionError("must replay")))
    assert replay == result


def test_completed_validation_resumes_activation_without_repeating_calls(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    action = _action(path)
    with monkeypatch.context() as interrupted:
        interrupted.setattr(
            DevelopmentJob, "activate", Mock(side_effect=InterruptedError("before writes"))
        )
        with pytest.raises(InterruptedError):
            _repair(tmp_path, registry, action)
    monkeypatch.setattr(
        "dnakit.tools.repair._run_variant",
        Mock(side_effect=AssertionError("must reuse validation")),
    )
    result = _repair(
        tmp_path, registry, action, Mock(side_effect=AssertionError("must reuse review"))
    )
    assert result["status"] == "activated"


def test_completed_baseline_resumes_only_candidate_and_reuses_quality(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    action = _action(path)
    original_phase = DevelopmentJob.phase

    def interrupt_before_candidate(self: DevelopmentJob, name: str, **data: Any) -> None:
        if name == "candidate":
            raise InterruptedError("Baseline was persisted; candidate has not started")
        original_phase(self, name, **data)

    with monkeypatch.context() as interrupted:
        interrupted.setattr(DevelopmentJob, "phase", interrupt_before_candidate)
        with pytest.raises(InterruptedError):
            _repair(tmp_path, registry, action)
    runner = Mock(wraps=repair_module._run_variant)
    monkeypatch.setattr(repair_module, "_run_variant", runner)
    monkeypatch.setattr(
        repair_module, "_quality", Mock(side_effect=AssertionError("reuse quality"))
    )
    result = _repair(tmp_path, registry, action, Mock(side_effect=AssertionError("reuse review")))
    assert result["status"] == "activated", result
    assert [call.args[2] for call in runner.call_args_list] == ["candidate"]


@pytest.mark.parametrize("tamper", ["formatted_code", "candidate_sha256"])
def test_resume_rejects_changed_code_or_validation_evidence(
    tmp_path: Path,
    source: tuple[Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
    tamper: str,
) -> None:
    path, registry = source
    action = _action(path)
    directories: list[Path] = []

    def interrupt(self: DevelopmentJob, writes: dict[Path, bytes]) -> dict[str, object]:
        directories.append(self.directory)
        raise InterruptedError("Before activation")

    with monkeypatch.context() as interrupted:
        interrupted.setattr(DevelopmentJob, "activate", interrupt)
        with pytest.raises(InterruptedError):
            _repair(tmp_path, registry, action)
    state_path = directories[0] / "state.json"
    state = json.loads(state_path.read_text())
    if tamper == "formatted_code":
        state["repair_quality"]["formatted_code"] += "\n"
    else:
        for event in state["events"]:
            if event["phase"] == "validated":
                event["validation"]["candidate_sha256"] = "0" * 64
    state_path.write_text(json.dumps(state))
    monkeypatch.setattr(
        repair_module, "_run_variant", Mock(side_effect=AssertionError("Do not repeat calls"))
    )
    result = _repair(tmp_path, registry, action, Mock(side_effect=AssertionError("No review")))
    assert result["status"] == "rejected"
    assert path.read_text() == MODULE


def test_unrelated_existing_type_errors_do_not_block_repair(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    path.write_text(MODULE + '\ndef existing_bad() -> int:\n    return "existing issue"\n')
    result = _repair(tmp_path, registry, _action(path))
    assert result["status"] == "activated", result
    assert 'return "existing issue"' in path.read_text()


def test_new_type_error_is_rejected_even_when_fixed_cases_would_pass(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    code = FIXED.replace(
        "    return abs(value)",
        '    if value == 7:\n        return "wrong type"\n    return abs(value)',
    )
    reviewer = _review()
    result = _repair(tmp_path, registry, _action(path, code=code), reviewer)
    assert result["status"] == "rejected" and result["code"] == "TOOL_REPAIR_QUALITY_FAILED"
    assert path.read_text() == MODULE
    reviewer.assert_not_called()


def test_missing_dependency_preserves_source_and_staged_candidate(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    # Existing implementation already attempts this optional dependency. A repair
    # may not claim correctness merely by removing an unavailable environment.
    broken = MODULE.replace(
        "    return value", "    import repair_missing_dependency\n    return value"
    )
    path.write_text(broken)
    action = _action(path)
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected"
    assert result["validation"]["baseline"]["cases"][0]["environment_failure"] is True
    assert path.read_text() == broken
    assert Path(result["audit_dir"], "candidate-01.py").exists()


@pytest.mark.parametrize(
    "failure",
    [
        'raise TimeoutError("Service unavailable")',
        "from dnakit.exceptions import ConfigurationError\n"
        '        raise ConfigurationError("Missing credential", code="MISSING_LLM_API_KEY")',
        'raise ValueError("Wrapped service error") from ConnectionError("Service unavailable")',
    ],
)
def test_environment_failures_cannot_be_removed_to_claim_a_successful_fix(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], failure: str
) -> None:
    path, registry = source
    broken = MODULE.replace(
        "    return value", f"    if value < 0:\n        {failure}\n    return value"
    )
    path.write_text(broken)
    result = _repair(tmp_path, registry, _action(path))
    assert result["status"] == "rejected", result
    assert result["validation"]["baseline"]["cases"][1]["environment_failure"] is True
    assert "candidate" not in result["validation"]
    assert path.read_text() == broken


def test_explicit_environment_error_contract_is_preserved_as_normal_regression(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry]
) -> None:
    path, registry = source
    path.write_text(MODULE.replace("raise ValueError", "raise TimeoutError"))
    action = _action(path, code=FIXED.replace("raise ValueError", "raise TimeoutError"))
    action["tests"][2]["expected_error"] = {"type": "TimeoutError"}
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "activated", result
    validation = result["validation"]
    assert validation["baseline"]["cases"][2]["environment_failure"] is False
    assert validation["candidate"]["cases"][2]["matched"] is True
    assert validation["improved_cases"] == [2]


def test_candidate_timeout_cannot_activate(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    monkeypatch.setattr("dnakit.tools.repair._MAX_SECONDS", 2.0)
    action = _action(
        path, code=FIXED.replace("    return abs(value)", "    while True:\n        pass")
    )
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "rejected"
    assert result["validation"]["candidate"]["status"] == "timeout"
    assert path.read_text() == MODULE


def test_original_source_and_contribution_files_rollback_together(
    tmp_path: Path, source: tuple[Path, AgentToolRegistry], monkeypatch: pytest.MonkeyPatch
) -> None:
    path, registry = source
    result = _repair(tmp_path, registry, _action(path))
    assert result["status"] == "activated", result
    monkeypatch.setattr("dnakit.tools.development_support._source_root", lambda: tmp_path)
    rollback = rollback_tool_change(result["change_id"], allow_side_effects=True, verbose=False)
    assert rollback["status"] == "rolled_back"
    assert path.read_text() == MODULE
    assert all(
        not Path(artifact).exists() for artifact in result["artifacts"] if artifact != str(path)
    )


@pytest.mark.parametrize("optional", [False, True])
def test_generated_repair_keeps_name_updates_manifest_and_preserves_previous_module(
    tmp_path: Path,
    source: tuple[Path, AgentToolRegistry],
    monkeypatch: pytest.MonkeyPatch,
    optional: bool,
) -> None:
    fixture_path, _ = source
    directory = fixture_path.parent / "generated"
    directory.mkdir()
    package_source = Path(generated.__file__).read_text()
    (directory / "__init__.py").write_text(package_source)
    original = "from __future__ import annotations\n\n" + FUNCTION
    digest = hashlib.sha256(original.encode()).hexdigest()
    old_path = directory / f"magnitude_{digest[:16]}.py"
    old_path.write_text(original)
    specification = importlib.util.spec_from_file_location(
        f"dnakit.generated.{old_path.stem}", old_path
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    monkeypatch.setitem(sys.modules, specification.name, module)
    specification.loader.exec_module(module)
    specification_tool = _build_spec("dnakit.repair_fixture", "magnitude", module.magnitude)
    registry = AgentToolRegistry(
        (replace(specification_tool, name="dnakit.generated.magnitude", category="generated"),)
    )
    manifest = {
        "version": 1,
        "tools": {
            "magnitude": {
                "module": old_path.name,
                "sha256": digest,
                "effect": "local",
                "description": "Return fixture magnitude.",
                "audit_dir": "/old/machine/local/path",
                "contract": {
                    "parameters": registry.resolve(
                        "dnakit.generated.magnitude"
                    ).call_plan.input_schema,
                    "returns": {"type": "integer"},
                    "dependencies": [],
                    "effects": ["local"],
                    "description": "Return fixture magnitude.",
                    "success_criteria": ["Nonnegative magnitude"],
                },
            }
        },
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))
    action = _action(old_path, name="dnakit.generated.magnitude")
    action["target"]["kind"] = "generated"
    if optional:
        action["code"] = FUNCTION.replace(
            "verbose: bool = True", "verbose: bool = True, positive: bool = False"
        ).replace("return value", "return abs(value) if positive else value")
        action["tests"][1]["arguments"]["positive"] = True
    result = _repair(tmp_path, registry, action)
    assert result["status"] == "activated", result
    assert result["kind"] == "generated" and result["tool_name"] == "dnakit.generated.magnitude"
    assert old_path.read_text() == original
    new_path = Path(result["source_path"])
    assert new_path != old_path and action["code"] in new_path.read_text()
    updated = generated._entries(directory)["magnitude"]
    assert updated["module"] == new_path.name
    assert updated["sha256"] == hashlib.sha256(new_path.read_bytes()).hexdigest()
    assert "audit_dir" not in updated
    if optional:
        properties = updated["contract"]["parameters"]["properties"]
        assert properties["positive"] == {"type": "boolean", "default": False}
