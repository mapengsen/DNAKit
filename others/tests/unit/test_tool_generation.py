"""Executable generation, rejection, recovery and public-export regression tests."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from dnakit import generated
from dnakit.exceptions import ConfigurationError
from dnakit.tools.catalog import _build_spec
from dnakit.tools.development_support import DevelopmentJob, rollback_tool_change
from dnakit.tools.generation import generate_tool
from dnakit.tools.registry import AgentToolRegistry

_CODE = '''from __future__ import annotations

def scale_integer_total(values: list[int], factor: int = 2, *, verbose: bool = True) -> int:
    """Multiply the sum of integer observations by an integer factor."""
    return sum(values) * factor
'''


@pytest.fixture
def action() -> dict[str, Any]:
    return {
        "action": "generate_tool",
        "name": "scale_integer_total",
        "reason": "The provided empty registry lacks a scaled integer total operation.",
        "contract": {
            "description": "Multiply an integer observation total by a supplied factor.",
            "parameters": {
                "type": "object",
                "properties": {
                    "values": {"type": "array", "items": {"type": "integer"}},
                    "factor": {"type": "integer", "default": 2},
                },
                "required": ["values"],
                "additionalProperties": False,
            },
            "returns": {"type": "integer"},
            "dependencies": [],
            "effects": ["local"],
            "success_criteria": ["Empty observations return zero; integer factors scale totals."],
        },
        "code": _CODE,
        "tests": [
            {"arguments": {"values": []}, "expected": 0},
            {"arguments": {"values": [1, 2], "factor": 3}, "expected": 9},
        ],
    }


@pytest.fixture
def reviewer() -> Mock:
    return Mock(return_value={"approved": True, "issues": []})


def _generate(root: Path, action: dict[str, Any], reviewer: Mock) -> dict[str, Any]:
    return generate_tool(
        action, root=root, registry=AgentToolRegistry(()), ask=reviewer, verbose=False
    )


def test_activation_requires_real_tests_and_writes_audited_manifest(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated", result
    assert result["tool_name"] == "dnakit.generated.scale_integer_total"
    metadata = generated._entries(tmp_path / "src/dnakit/generated")[action["name"]]
    module = tmp_path / "src/dnakit/generated" / metadata["module"]
    assert "def scale_integer_total(" in module.read_text()
    assert '__all__ = ["scale_integer_total"]' in module.read_text()
    assert "audit_dir" not in metadata
    assert result["validation"]["quality"]["checks"]["mypy"]["passed"] is True
    assert len(list((tmp_path / "others/tests/generated").glob("test_*.py"))) == 1
    assert len(list((tmp_path / "docs/llm_tools/generated").glob("*.md"))) == 2
    assert (tmp_path / "registry/contributions" / f"{result['change_id']}.md").is_file()
    assert metadata["sha256"] == hashlib.sha256(module.read_bytes()).hexdigest()
    assert metadata["effect"] == "local"
    assert metadata["contract"] == action["contract"]
    assert result["validation"]["tests"][-1]["observed"] == 9
    assert "Test 2/2 passed" in result["validation"]["stdout"]
    assert Path(result["audit_dir"], "changes.diff").exists()
    assert json.loads(Path(result["audit_dir"], "cases-01.json").read_text()) == action["tests"]
    usage = Path(result["audit_dir"], "usage-01.md").read_text()
    assert "from dnakit.generated import scale_integer_total" in usage
    assert "call_dnakit_tool" in usage and "dnakit.generated.scale_integer_total" in usage
    assert reviewer.call_count == 1


def test_replay_does_not_repeat_review_or_execution(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = _generate(tmp_path, action, reviewer)
    assert original["status"] == "activated", original
    monkeypatch.setattr(
        "dnakit.tools.generation._run_tests", Mock(side_effect=AssertionError("must reuse receipt"))
    )
    reviewer.side_effect = AssertionError("must not repeat paid review")
    assert _generate(tmp_path, action, reviewer) == original


def test_failed_tests_leave_public_package_untouched_and_allow_corrected_revision(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    broken = {**action, "code": _CODE.replace("sum(values) * factor", "sum(values) + factor")}
    result = _generate(tmp_path, broken, reviewer)
    assert result["status"] == "rejected"
    assert result["validation"]["status"] == "failed"
    assert "output mismatch" in result["validation"]["error"]
    assert not (tmp_path / "src/dnakit/generated/manifest.json").exists()
    fixed = _generate(tmp_path, action, reviewer)
    assert fixed["status"] == "activated", fixed
    state = json.loads(Path(fixed["audit_dir"], "state.json").read_text())
    assert state["attempt"] == 2


def test_review_rejection_is_durable_without_executing_code(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    reviewer.return_value = {"approved": False, "issues": ["Insufficient domain justification"]}
    runner = Mock(side_effect=AssertionError("Rejected code must not execute"))
    monkeypatch.setattr("dnakit.tools.generation._run_tests", runner)
    first = _generate(tmp_path, action, reviewer)
    assert first["status"] == "rejected"
    assert _generate(tmp_path, action, reviewer) == first
    assert reviewer.call_count == 1
    runner.assert_not_called()


@pytest.mark.parametrize(
    "code",
    [
        "open('unexpected', 'w').write('bad')\n" + _CODE,
        _CODE.replace("verbose: bool = True", "debug: bool = True"),
        _CODE + "\ndef second_public() -> int:\n    return 1\n",
        _CODE.replace("from __future__ import annotations", "import sys"),
        _CODE.replace("list[int]", """list["__import__('os').system('echo unsafe')"]"""),
        _CODE.replace("from __future__ import annotations", "import pathlib"),
        _CODE.replace("sum(values) * factor", "eval('1')"),
        _CODE.replace("    return sum", "    writer = open\n    return sum"),
    ],
)
def test_static_rejections_never_request_paid_review(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, code: str
) -> None:
    result = _generate(tmp_path, {**action, "code": code}, reviewer)
    assert result["status"] == "rejected"
    reviewer.assert_not_called()
    assert not (tmp_path / "src/dnakit/generated").exists()


@pytest.mark.parametrize("field", ["tests", "contract"])
def test_validation_data_is_frozen_across_revisions(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, field: str
) -> None:
    reviewer.return_value = {"approved": False, "issues": ["Fix implementation"]}
    _generate(tmp_path, action, reviewer)
    changed = copy.deepcopy(action)
    if field == "tests":
        changed["tests"][1]["expected"] = 10
    else:
        changed["contract"]["description"] = "Changed goal after seeing failures"
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_FIXED_SCENARIOS"):
        _generate(tmp_path, changed, reviewer)


def test_existing_public_tool_cannot_be_replaced(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    activated = _generate(tmp_path, action, reviewer)
    assert activated["status"] == "activated", activated
    original = (tmp_path / "src/dnakit/generated/manifest.json").read_bytes()
    result = _generate(tmp_path, {**action, "reason": "Try replacing an active tool"}, reviewer)
    assert result["status"] == "rejected"
    assert "already exists" in str(result["reason"])
    assert (tmp_path / "src/dnakit/generated/manifest.json").read_bytes() == original
    assert reviewer.call_count == 1


def test_concurrent_manifest_edit_is_preserved(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    target = tmp_path / "src/dnakit/generated/manifest.json"

    def review(*args: object) -> dict[str, object]:
        target.parent.mkdir(parents=True)
        target.write_text('{"version": 1, "tools": {}, "other_writer": true}')
        return {"approved": True, "issues": []}

    reviewer.side_effect = review
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_CONFLICT"):
        _generate(tmp_path, action, reviewer)
    assert json.loads(target.read_text())["other_writer"] is True
    assert list(target.parent.glob("*.py")) == []


def test_activation_recovery_adds_tool_name_without_reexecuting(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    activated = _generate(tmp_path, action, reviewer)
    assert activated["status"] == "activated", activated
    state_path = Path(activated["audit_dir"], "state.json")
    state = json.loads(state_path.read_text())
    state["phase"] = "activating"
    state.pop("result")
    state_path.write_text(json.dumps(state))
    monkeypatch.setattr(
        "dnakit.tools.generation._run_tests",
        Mock(side_effect=AssertionError("must recover writes")),
    )
    reviewer.side_effect = AssertionError("must recover paid review")
    restored = _generate(tmp_path, action, reviewer)
    assert restored["status"] == "activated"
    assert restored["tool_name"] == activated["tool_name"]


def test_interrupted_validation_is_not_repeated(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    job = DevelopmentJob(tmp_path, "generation", action["name"], action, False)
    job.phase("validating")
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "rejected"
    assert result["code"] == "AGENT_DEVELOPMENT_VALIDATION_UNCERTAIN"
    reviewer.assert_not_called()


def test_completed_validation_is_reused_after_interruption_before_activation(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    with monkeypatch.context() as interrupted:
        interrupted.setattr(
            DevelopmentJob, "activate", Mock(side_effect=InterruptedError("before activation"))
        )
        with pytest.raises(InterruptedError, match="before activation"):
            _generate(tmp_path, action, reviewer)
    reviewer.side_effect = AssertionError("completed review must be reused")
    monkeypatch.setattr(
        "dnakit.tools.generation._run_tests",
        Mock(side_effect=AssertionError("completed validation must be reused")),
    )
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated"
    assert result["validation"]["tests"][-1]["observed"] == 9


def test_existing_native_function_alias_is_not_shadowed(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    def native_total(values: list[int], factor: int = 2) -> int:
        return sum(values) * factor

    registry = AgentToolRegistry((_build_spec("dnakit.ops", action["name"], native_total),))
    result = generate_tool(action, root=tmp_path, registry=registry, ask=reviewer, verbose=False)
    assert result["status"] == "rejected"
    assert "already exists" in str(result["reason"])
    reviewer.assert_not_called()


def test_three_revisions_are_the_maximum(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    reviewer.return_value = {"approved": False, "issues": ["Need stronger justification"]}
    for revision in range(3):
        assert (
            _generate(tmp_path, {**action, "code": _CODE + "\n" * revision}, reviewer)["status"]
            == "rejected"
        )
    with pytest.raises(ConfigurationError, match="AGENT_DEVELOPMENT_ROUND_LIMIT"):
        _generate(tmp_path, {**action, "code": _CODE + "\n" * 3}, reviewer)
    assert reviewer.call_count == 3


def test_third_revision_can_resume_completed_validation_without_a_fourth_attempt(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    reviewer.return_value = {"approved": False, "issues": ["Revise implementation"]}
    for revision in (1, 2):
        candidate = {**action, "code": _CODE + f"\n# revision {revision}\n"}
        assert _generate(tmp_path, candidate, reviewer)["status"] == "rejected"
    reviewer.return_value = {"approved": True, "issues": []}
    candidate = {**action, "code": _CODE + "\n# revision 3\n"}
    with monkeypatch.context() as interrupted:
        interrupted.setattr(
            DevelopmentJob, "activate", Mock(side_effect=InterruptedError("before activation"))
        )
        with pytest.raises(InterruptedError):
            _generate(tmp_path, candidate, reviewer)
    reviewer.side_effect = AssertionError("third revision already passed review")
    monkeypatch.setattr(
        "dnakit.tools.generation._run_tests",
        Mock(side_effect=AssertionError("third revision already passed validation")),
    )
    result = _generate(tmp_path, candidate, reviewer)
    assert result["status"] == "activated"
    assert json.loads(Path(result["audit_dir"], "state.json").read_text())["attempt"] == 3


def test_public_exports_refresh_and_rollback_removes_generated_function(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated", result
    try:
        with monkeypatch.context() as isolated:
            isolated.setattr(generated, "_DIRECTORY", tmp_path / "src/dnakit/generated")
            isolated.setattr("dnakit.tools.development_support._source_root", lambda: tmp_path)
            generated._refresh()
            assert generated.__all__ == [action["name"]]
            function = getattr(generated, action["name"])
            assert function([2, 3], factor=4, verbose=False) == 20
            module_name = function.__module__
            assert module_name in sys.modules
            rollback = rollback_tool_change(
                result["change_id"], allow_side_effects=True, verbose=False
            )
            assert rollback["status"] == "rolled_back"
            assert generated.__all__ == []
            assert not hasattr(generated, action["name"])
            assert module_name not in sys.modules
            assert not Path(result["artifacts"][0]).exists()
    finally:
        generated._refresh()


def test_modified_generated_source_is_not_imported(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated", result
    Path(result["artifacts"][0]).write_text(
        "raise AssertionError('must never execute tampered code')"
    )
    try:
        with monkeypatch.context() as isolated:
            isolated.setattr(generated, "_DIRECTORY", tmp_path / "src/dnakit/generated")
            with pytest.raises(ConfigurationError, match="GENERATED_MODULE_INVALID"):
                generated._refresh()
    finally:
        generated._refresh()


def test_test_timeout_prevents_activation(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("dnakit.tools.generation._MAX_TEST_SECONDS", 0.3)
    code = _CODE.replace("    return sum(values) * factor", "    while True:\n        pass")
    result = _generate(tmp_path, {**action, "code": code}, reviewer)
    assert result["status"] == "rejected"
    assert result["validation"]["status"] == "timeout"
    assert not (tmp_path / "src/dnakit/generated/manifest.json").exists()


def test_output_flood_prevents_activation(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("dnakit.tools.generation._MAX_LOG_BYTES", 128)
    code = _CODE.replace("    return sum", "    print('x' * 10000)\n    return sum")
    result = _generate(tmp_path, {**action, "code": code}, reviewer)
    assert result["status"] == "rejected"
    assert result["validation"]["status"] == "output_limit"


def test_missing_dependency_is_reported_without_installation(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    action["contract"]["dependencies"] = ["dnakit_dependency_that_does_not_exist"]
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "rejected"
    assert "not installed" in result["reason"]
    reviewer.assert_not_called()


@pytest.mark.parametrize("statement", ["import rich", "from rich.console import Console"])
def test_third_party_imports_must_stay_inside_functions(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, statement: str
) -> None:
    action["contract"]["dependencies"] = ["rich"]
    action["code"] = _CODE.replace("from __future__ import annotations", statement)
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "rejected"
    assert "imported inside functions" in result["reason"]
    reviewer.assert_not_called()


def test_missing_optional_dependency_does_not_break_generated_namespace_import(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    dependencies = tmp_path / "dependencies"
    dependencies.mkdir()
    module_path = dependencies / "optional_generated_scaler.py"
    module_path.write_text("def scale(value, factor):\n    return value * factor\n")
    monkeypatch.syspath_prepend(str(dependencies))
    monkeypatch.setenv("PYTHONPATH", str(dependencies))
    action["contract"]["dependencies"] = ["optional_generated_scaler"]
    action["contract"]["effects"] = ["external_process"]
    action["code"] = _CODE.replace(
        "    return sum(values) * factor",
        "    import optional_generated_scaler  # type: ignore[import-untyped]\n"
        "    return int(optional_generated_scaler.scale(sum(values), factor))",
    )
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated", result
    module_path.unlink()
    try:
        with monkeypatch.context() as isolated:
            isolated.setattr(generated, "_DIRECTORY", tmp_path / "src/dnakit/generated")
            generated._refresh()
            function = getattr(generated, action["name"])
            with pytest.raises(ModuleNotFoundError, match="optional_generated_scaler"):
                function([2], verbose=False)
    finally:
        generated._refresh()


def test_runtime_output_must_also_match_the_python_annotation(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Exercise the execution validator independently of the earlier static type check.
    monkeypatch.setattr(
        "dnakit.tools.generation.check_python_quality", lambda code, **kwargs: (code, {})
    )
    code = _CODE.replace(") -> int:", ") -> str:")
    result = _generate(tmp_path, {**action, "code": code}, reviewer)
    assert result["status"] == "rejected"
    assert "Output schema differs" in result["validation"]["error"]


@pytest.mark.parametrize("change", ["constraint", "default", "input_type", "generic_output"])
def test_contract_schema_must_match_actual_callable(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, change: str
) -> None:
    if change == "constraint":
        action["contract"]["parameters"]["properties"]["factor"]["minimum"] = 0
    elif change == "default":
        action["contract"]["parameters"]["properties"]["factor"]["default"] = 5
    elif change == "input_type":
        action["contract"]["parameters"]["properties"]["factor"]["type"] = "number"
    else:
        action["contract"]["returns"] = {}
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "rejected"
    assert "schema differs" in result["validation"]["error"]


def test_descriptions_and_omitted_generic_verbose_do_not_break_schema_comparison(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    properties = action["contract"]["parameters"]["properties"]
    properties["factor"]["description"] = "Multiply the observation sum by this integer."
    properties["factor"]["examples"] = [2, 3]
    action["contract"]["parameters"].pop("additionalProperties")
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "activated", result


def test_parameter_named_description_keeps_its_schema_constraints(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock
) -> None:
    parameters = action["contract"]["parameters"]
    parameters["properties"]["description"] = parameters["properties"].pop("values")
    parameters["properties"]["description"]["maxItems"] = 10
    parameters["required"] = ["description"]
    action["code"] = _CODE.replace("values", "description")
    for test in action["tests"]:
        test["arguments"]["description"] = test["arguments"].pop("values")
    result = _generate(tmp_path, action, reviewer)
    assert result["status"] == "rejected"
    assert "Input schema differs" in result["validation"]["error"]


def test_manifest_cannot_downgrade_a_declared_effect(tmp_path: Path) -> None:
    directory = tmp_path / "generated"
    directory.mkdir()
    manifest = {
        "version": 1,
        "tools": {
            "file_reader": {
                "module": "file_reader_" + "a" * 16 + ".py",
                "sha256": "a" * 64,
                "effect": "local",
                "contract": {"effects": ["filesystem"]},
            }
        },
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ConfigurationError, match="GENERATED_MANIFEST_INVALID"):
        generated._entries(directory)


@pytest.mark.parametrize("name", ["../escape", "_private", "class", "globals", "len", "9invalid"])
def test_invalid_names_fail_before_creating_job(
    tmp_path: Path, action: dict[str, Any], reviewer: Mock, name: str
) -> None:
    with pytest.raises(ConfigurationError, match="GENERATED_TOOL_INVALID"):
        _generate(tmp_path, {**action, "name": name}, reviewer)
    assert list(tmp_path.iterdir()) == []
