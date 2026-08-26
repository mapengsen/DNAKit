"""Tests for strict, bounded, resumable local workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError
from dnakit.workflows import WorkflowProgress, load_manifest, load_workflow, run_workflow


def _step(
    step_id: str,
    operation: str,
    input_ref: str,
    params: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        "id": step_id,
        "operation": operation,
        "input": input_ref,
        "params": params or {},
    }


def _workflow(
    input_name: str,
    steps: list[dict[str, object]],
    *,
    output_dir: str = "workflow-output",
    error_policy: str = "raise",
    overwrite: bool = False,
    limits: dict[str, object] | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "dnakit-workflow-v1",
        "run_id": "unit-run",
        "input": {"path": input_name, "format": "fasta", "alphabet": "iupac"},
        "output_dir": output_dir,
        "seed": 17,
        "error_policy": error_policy,
        "overwrite": overwrite,
        "steps": steps,
    }
    if limits is not None:
        payload["limits"] = limits
    return payload


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_yaml_schema_resolves_paths_and_rejects_unknown_fields(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.yaml"
    config.write_text(
        yaml.safe_dump(
            _workflow(
                source.name,
                [
                    _step(
                        "normalize",
                        "normalize",
                        "input",
                        {
                            "keep_ambiguous": True,
                            "keep_u": False,
                            "keep_other": False,
                        },
                    )
                ],
            ),
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    loaded = load_workflow(config)

    assert loaded.input_path == source
    assert loaded.output_dir == tmp_path / "workflow-output"
    assert loaded.spec.sha256 == loaded.spec.sha256
    assert loaded.spec.steps[0].operation == "normalize"

    payload = _workflow(source.name, [_step("normalize", "normalize", "input")])
    payload["typo"] = True
    _write_json(config.with_suffix(".json"), payload)
    with pytest.raises(ConfigurationError) as error:
        load_workflow(config.with_suffix(".json"))
    assert error.value.code == "UNKNOWN_WORKFLOW_FIELD"


@pytest.mark.parametrize(
    ("steps", "code"),
    [
        ([_step("bad", "shell", "input", {"command": "echo no"})], "UNKNOWN_WORKFLOW_OPERATION"),
        (
            [_step("bad", "write", "input", {"target": "../escape.fasta"})],
            "UNSAFE_WORKFLOW_OUTPUT_PATH",
        ),
        (
            [_step("bad", "write", "input", {"target": "..\\escape.fasta"})],
            "UNSAFE_WORKFLOW_OUTPUT_PATH",
        ),
        ([_step("later", "normalize", "missing")], "UNKNOWN_WORKFLOW_INPUT_REFERENCE"),
        (
            [
                _step("one", "write", "input", {"target": "same.fasta"}),
                _step("two", "write", "one", {"target": "same.fasta"}),
            ],
            "DUPLICATE_WORKFLOW_OUTPUT",
        ),
        (
            [_step("large", "fingerprint", "input", {"k": 10})],
            "WORKFLOW_FINGERPRINT_LIMIT",
        ),
    ],
)
def test_schema_rejects_unsafe_or_unbounded_steps(
    tmp_path: Path, steps: list[dict[str, object]], code: str
) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(config, _workflow(source.name, steps))

    with pytest.raises(ConfigurationError) as error:
        load_workflow(config)

    assert error.value.code == code


def test_schema_rejects_duplicate_json_and_yaml_keys(tmp_path: Path) -> None:
    duplicate_json = tmp_path / "duplicate.json"
    duplicate_json.write_text(
        '{"schema_version":"dnakit-workflow-v1","schema_version":"again"}',
        encoding="utf-8",
    )
    duplicate_yaml = tmp_path / "duplicate.yaml"
    duplicate_yaml.write_text(
        "schema_version: dnakit-workflow-v1\nschema_version: again\n",
        encoding="utf-8",
    )

    for path in (duplicate_json, duplicate_yaml):
        with pytest.raises(InputFormatError) as error:
            load_workflow(path)
        assert error.value.code == "DUPLICATE_WORKFLOW_FIELD"


def test_yaml_aliases_and_unsafe_tags_are_not_loaded(tmp_path: Path) -> None:
    alias = tmp_path / "alias.yaml"
    alias.write_text("value: &shared 1\ncopy: *shared\n", encoding="utf-8")
    unsafe_tag = tmp_path / "tag.yaml"
    unsafe_tag.write_text("!!python/object/apply:os.system ['echo forbidden']\n", encoding="utf-8")

    with pytest.raises(InputFormatError) as alias_error:
        load_workflow(alias)
    assert alias_error.value.code == "WORKFLOW_YAML_ALIAS_DISABLED"
    with pytest.raises(InputFormatError) as tag_error:
        load_workflow(unsafe_tag)
    assert tag_error.value.code == "INVALID_WORKFLOW_CONFIG"


def test_workflow_config_depth_and_byte_limits_are_structured(tmp_path: Path) -> None:
    deep = tmp_path / "deep.json"
    deep.write_text("[" * 1_200 + "]" * 1_200, encoding="utf-8")
    oversized = tmp_path / "oversized.json"
    oversized.write_bytes(b" " * 1_000_001)

    with pytest.raises(InputFormatError) as deep_error:
        load_workflow(deep)
    assert deep_error.value.code == "INVALID_WORKFLOW_CONFIG"
    with pytest.raises(InputFormatError) as size_error:
        load_workflow(oversized)
    assert size_error.value.code == "WORKFLOW_CONFIG_SIZE_LIMIT"


def test_end_to_end_workflow_records_all_steps_and_artifacts(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n>two\nACGT\n>three\nGGGG\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    steps = [
        _step("normalized", "normalize", "input"),
        _step("validated", "validate", "normalized", {"sequence_length": 4}),
        _step(
            "described",
            "descriptors",
            "validated",
            {"metrics": ["length", "gc", "entropy"]},
        ),
        _step("fingerprinted", "fingerprint", "described", {"k": 2}),
        _step("deduplicated", "deduplicate", "fingerprinted"),
        _step(
            "partitioned",
            "split",
            "deduplicated",
            {"method": "random", "ratios": {"train": 0.5, "test": 0.5}},
        ),
        _step(
            "written",
            "write",
            "partitioned",
            {"target": "splits/train.fasta", "subset": "train"},
        ),
        _step(
            "reported",
            "report",
            "written",
            {"target": "report.html", "subset": "test", "title": "Local workflow"},
        ),
    ]
    _write_json(config, _workflow(source.name, steps))
    events: list[WorkflowProgress] = []

    result = run_workflow(config, progress=events.append)

    assert result.status == "succeeded"
    assert result.resumed is False
    assert [step.status for step in result.steps] == ["succeeded"] * len(steps)
    assert (tmp_path / "workflow-output" / "splits" / "train.fasta").is_file()
    report = tmp_path / "workflow-output" / "report.html"
    assert report.is_file()
    assert "Local workflow" in report.read_text(encoding="utf-8")
    assert {event.status for event in events} >= {"started", "succeeded"}

    manifest = load_manifest(tmp_path / "workflow-output" / "run-manifest.json")
    assert manifest["status"] == "succeeded"
    assert manifest["seed"] == 17
    assert manifest["command"] == [
        "dnakit",
        "workflow",
        str(config.absolute()),
    ]
    outputs = manifest["outputs"]
    assert isinstance(outputs, list)
    assert len(outputs) == 2
    inputs = manifest["inputs"]
    assert isinstance(inputs, list)
    assert len(inputs) == 2
    assert all(item["sha256"] for item in inputs)
    resolved = manifest["resolved_config"]
    assert isinstance(resolved, dict)
    assert resolved["schema_version"] == "dnakit-workflow-execution-v1"
    execution = resolved["execution"]
    assert isinstance(execution, dict)
    audited_steps = execution["steps"]
    assert isinstance(audited_steps, list)
    assert [item["id"] for item in audited_steps] == [step["id"] for step in steps]
    assert all(item["duration_seconds"] >= 0 for item in audited_steps)
    assert audited_steps[-1]["artifacts"][0]["sha256"]


def test_dry_run_is_read_only_and_emits_planned_progress(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [
                _step("described", "descriptors", "input"),
                _step("written", "write", "described", {"target": "records.json"}),
            ],
        ),
    )
    events: list[WorkflowProgress] = []

    result = run_workflow(config, dry_run=True, progress=events.append)

    assert result.status == "dry-run"
    assert result.manifest_path is None
    assert [step.status for step in result.steps] == ["planned", "planned"]
    assert [event.status for event in events] == ["planned", "planned"]
    assert not (tmp_path / "workflow-output").exists()


def test_dry_run_and_resume_are_mutually_exclusive(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(config, _workflow(source.name, [_step("normalized", "normalize", "input")]))

    with pytest.raises(ConfigurationError) as error:
        run_workflow(config, dry_run=True, resume=True)

    assert error.value.code == "INVALID_WORKFLOW_MODE"


def test_dry_run_rejects_target_symlink_without_writing(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    output = tmp_path / "workflow-output"
    output.mkdir()
    (output / ".dnakit-workflow-output-v1").write_text(
        "dnakit-workflow-execution-v1\n", encoding="utf-8"
    )
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "linked").symlink_to(outside, target_is_directory=True)
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [_step("written", "write", "input", {"target": "linked/out.fasta"})],
        ),
    )

    with pytest.raises(ConfigurationError) as error:
        run_workflow(config, dry_run=True)

    assert error.value.code == "UNSAFE_WORKFLOW_OUTPUT_PATH"
    assert not (outside / "out.fasta").exists()


def test_progress_callback_failure_does_not_change_scientific_result(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(config, _workflow(source.name, [_step("normalized", "normalize", "input")]))

    def broken_callback(event: WorkflowProgress) -> None:
        del event
        raise RuntimeError("display failed")

    result = run_workflow(config, progress=broken_callback)

    assert result.status == "succeeded"
    assert result.steps[0].status == "succeeded"


def test_resume_skips_only_checksum_verified_outputs(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n>two\nGGGG\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [
                _step("described", "descriptors", "input"),
                _step("written", "write", "described", {"target": "records.fasta"}),
                _step("reported", "report", "written", {"target": "report.html"}),
            ],
        ),
    )
    first = run_workflow(config)
    assert first.status == "succeeded"

    resumed = run_workflow(config, resume=True)

    assert resumed.status == "succeeded"
    assert resumed.resumed is True
    assert [step.status for step in resumed.steps] == ["succeeded", "skipped", "skipped"]

    report = tmp_path / "workflow-output" / "report.html"
    report.write_text("tampered", encoding="utf-8")
    with pytest.raises(ConfigurationError) as error:
        run_workflow(config, resume=True)
    assert error.value.code == "WORKFLOW_RESUME_ARTIFACT_MISMATCH"


def test_resume_rejects_changed_input_before_skipping(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(source.name, [_step("written", "write", "input", {"target": "out.fasta"})]),
    )
    run_workflow(config)
    source.write_text(">one\nAAAA\n", encoding="utf-8")

    with pytest.raises(ConfigurationError) as error:
        run_workflow(config, resume=True)

    assert error.value.code == "WORKFLOW_RESUME_MISMATCH"


def test_collect_policy_runs_independent_branch_and_records_failure(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">ambiguous\nACNT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [
                _step("failed_descriptor", "descriptors", "input", {"metrics": ["gc"]}),
                _step("blocked_report", "report", "failed_descriptor", {"target": "blocked.html"}),
                _step("independent_write", "write", "input", {"target": "input.fasta"}),
            ],
            error_policy="collect",
        ),
    )

    result = run_workflow(config)

    assert result.status == "failed"
    assert [step.status for step in result.steps] == ["failed", "blocked", "succeeded"]
    assert (tmp_path / "workflow-output" / "input.fasta").is_file()
    manifest = load_manifest(tmp_path / "workflow-output" / "run-manifest.json")
    assert manifest["status"] == "failed"
    assert manifest["issues"]


def test_raise_policy_persists_failed_manifest_before_reraising(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">ambiguous\nACNT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [_step("failed_descriptor", "descriptors", "input", {"metrics": ["gc"]})],
        ),
    )

    with pytest.raises(DNAKitError):
        run_workflow(config)

    manifest = load_manifest(tmp_path / "workflow-output" / "run-manifest.json")
    assert manifest["status"] == "failed"
    resolved = manifest["resolved_config"]
    assert isinstance(resolved, dict)
    execution = resolved["execution"]
    assert isinstance(execution, dict)
    steps = execution["steps"]
    assert isinstance(steps, list)
    assert steps[0]["status"] == "failed"


def test_existing_unmarked_output_directory_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    output = tmp_path / "workflow-output"
    output.mkdir()
    (output / "user-data.txt").write_text("keep", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(config, _workflow(source.name, [_step("normalize", "normalize", "input")]))

    with pytest.raises(ConfigurationError) as error:
        run_workflow(config, dry_run=True)

    assert error.value.code == "UNSAFE_WORKFLOW_OUTPUT_DIR"
    assert (output / "user-data.txt").read_text(encoding="utf-8") == "keep"


def test_input_resource_failure_obeys_collect_policy(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [_step("normalize", "normalize", "input")],
            error_policy="collect",
            limits={"max_sequence_length": 3},
        ),
    )

    result = run_workflow(config)

    assert result.status == "failed"
    assert result.steps == ()
    manifest = load_manifest(tmp_path / "workflow-output" / "run-manifest.json")
    assert manifest["status"] == "failed"


def test_input_physical_byte_limit_applies_before_execution(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(
        config,
        _workflow(
            source.name,
            [_step("normalize", "normalize", "input")],
            limits={"max_input_bytes": 5},
        ),
    )

    with pytest.raises(InputFormatError) as error:
        run_workflow(config)
    assert error.value.code == "ARTIFACT_SIZE_LIMIT"
    assert not (tmp_path / "workflow-output").exists()


def test_result_objects_are_json_serializable(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    _write_json(config, _workflow(source.name, [_step("normalize", "normalize", "input")]))

    result = run_workflow(config)
    payload: dict[str, Any] = result.to_dict()

    assert json.loads(json.dumps(payload))["run_id"] == "unit-run"
