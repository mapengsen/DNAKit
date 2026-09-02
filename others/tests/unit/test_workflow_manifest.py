"""Tests for local workflow manifest persistence."""

import json
from pathlib import Path

import pytest

from dnakit.core import Issue, IssueSeverity
from dnakit.workflows import RunManifestBuilder, artifact_from_path, load_manifest, save_manifest


def test_manifest_builder_records_artifacts_seed_status_and_issues(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    input_artifact = artifact_from_path(source, media_type="text/x-fasta", schema_version="1")
    builder = RunManifestBuilder(
        "run-1",
        ("dnakit", "validate", str(source)),
        {"allow_empty": False},
        seed=7,
        seed_derivation="master-seed+input-index",
    )
    builder.add_input(input_artifact)
    builder.add_issue(Issue("NOTICE", IssueSeverity.INFO, "Checked input."))

    manifest = builder.build(status="succeeded", finished_at="2026-01-01T00:00:00+00:00")
    target = tmp_path / "run.json"
    artifact = save_manifest(manifest, target)
    restored = load_manifest(target)

    assert restored["run_id"] == "run-1"
    assert restored["seed"] == 7
    assert restored["status"] == "succeeded"
    assert artifact.sha256
    json.dumps(restored, sort_keys=True)


def test_manifest_writer_refuses_overwrite_and_cleans_temporary_file(tmp_path: Path) -> None:
    manifest = RunManifestBuilder("run", ("dnakit", "info"), {}).build()
    target = tmp_path / "run.json"
    save_manifest(manifest, target)

    with pytest.raises(FileExistsError):
        save_manifest(manifest, target)

    assert len(list(tmp_path.glob(".run.json.*.tmp"))) == 0


def test_artifact_rejects_directory_and_manifest_rejects_invalid_json(tmp_path: Path) -> None:
    with pytest.raises(Exception) as artifact_error:
        artifact_from_path(tmp_path, media_type="text/plain", schema_version="1")
    assert getattr(artifact_error.value, "code", None) == "INVALID_ARTIFACT_PATH"

    invalid = tmp_path / "invalid.json"
    invalid.write_text("[]", encoding="utf-8")
    with pytest.raises(Exception) as manifest_error:
        load_manifest(invalid)
    assert getattr(manifest_error.value, "code", None) == "INVALID_RUN_MANIFEST"


def test_manifest_loader_rejects_invalid_status_and_symlink_artifacts(tmp_path: Path) -> None:
    manifest = RunManifestBuilder("run", ("dnakit", "info"), {}).build()
    target = tmp_path / "run.json"
    save_manifest(manifest, target)
    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["status"] = "unknown"
    target.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(Exception) as status_error:
        load_manifest(target)
    assert getattr(status_error.value, "code", None) == "INVALID_RUN_MANIFEST"

    source = tmp_path / "source.txt"
    source.write_text("data", encoding="utf-8")
    link = tmp_path / "link.txt"
    link.symlink_to(source)
    with pytest.raises(Exception) as symlink_error:
        artifact_from_path(link, media_type="text/plain", schema_version="1")
    assert getattr(symlink_error.value, "code", None) == "INVALID_ARTIFACT_PATH"


def test_manifest_loader_reports_excessive_nesting(tmp_path: Path) -> None:
    target = tmp_path / "deep.json"
    target.write_text("[" * 1_200 + "0" + "]" * 1_200, encoding="utf-8")

    with pytest.raises(Exception) as error:
        load_manifest(target)
    assert getattr(error.value, "code", None) == "RUN_MANIFEST_STRUCTURE_LIMIT"


def test_manifest_loader_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    target = tmp_path / "duplicate.json"
    target.write_text('{"run_id":"a","run_id":"b"}', encoding="utf-8")

    with pytest.raises(Exception) as error:
        load_manifest(target)
    assert getattr(error.value, "code", None) == "INVALID_RUN_MANIFEST"


def test_manifest_writer_rejects_structurally_deep_payload_before_writing(tmp_path: Path) -> None:
    nested: object = "leaf"
    for _ in range(100):
        nested = [nested]
    target = tmp_path / "deep.json"

    with pytest.raises(Exception) as error:
        manifest = RunManifestBuilder("run", ("dnakit", "info"), {"nested": nested}).build()
        save_manifest(manifest, target)
    assert getattr(error.value, "code", None) == "JSON_STRUCTURE_LIMIT"
    assert not target.exists()
