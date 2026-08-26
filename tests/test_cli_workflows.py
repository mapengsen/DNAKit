"""End-to-end tests for the local DNAKit CLI workflows."""

from __future__ import annotations

import importlib
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from dnakit.cli.app import app

runner = CliRunner()


def test_cli_analysis_commands_emit_json() -> None:
    commands = (
        ["describe", "ACGTACGT"],
        ["fingerprint", "ACGTACGT", "--kind", "minhash", "--k", "2"],
        ["search", "TTACGT", "ACG"],
        ["orfs", "ATGAAATAA"],
        ["compare", "ACGT", "ACGA", "--method", "hamming"],
    )

    for command in commands:
        result = runner.invoke(app, command)
        assert result.exit_code == 0, result.output
        json.loads(result.stdout)


def test_cli_describe_defaults_to_240_fields_and_retains_compact_mode() -> None:
    complete = runner.invoke(app, ["describe", "ACGT"])
    compact = runner.invoke(app, ["describe", "ACGT", "--compact"])

    assert complete.exit_code == 0, complete.output
    complete_payload = json.loads(complete.stdout)
    assert complete_payload["schema_version"] == "descriptor_schema_v1"
    assert len(complete_payload["values"]) == 240
    assert compact.exit_code == 0, compact.output
    assert set(json.loads(compact.stdout)) == {
        "base_composition",
        "gc_at",
        "complexity",
        "repeat",
    }


def test_cli_normalize_and_validate_emit_json() -> None:
    normalized = runner.invoke(app, ["normalize", " acgu "])
    validated = runner.invoke(app, ["validate", "ACGT", "--sequence-length", "4"])
    invalid_length = runner.invoke(app, ["validate", "ACGT", "--sequence-length", "5"])

    assert normalized.exit_code == 0, normalized.output
    assert json.loads(normalized.stdout)["normalized_parts"] == ["ACG"]
    assert validated.exit_code == 0, validated.output
    assert json.loads(validated.stdout)["is_valid"] is True
    assert invalid_length.exit_code == 2
    invalid_payload = json.loads(invalid_length.stdout)
    assert invalid_payload["is_valid"] is False
    assert invalid_payload["issues"][0]["code"] == "STD_SEQUENCE_LENGTH_MISMATCH"


def test_cli_download_genome_emits_manifest_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class _FakeDownloadResult:
        def to_dict(self) -> dict[str, object]:
            return {"accession": "GCF_000001405.40", "fasta_path": "reference.fna"}

    def fake_download(query: str, output_dir: Path, **kwargs: object) -> _FakeDownloadResult:
        assert query == "hg38"
        assert output_dir == tmp_path
        assert kwargs.get("progress") is None
        return _FakeDownloadResult()

    cli_app = importlib.import_module("dnakit.cli.app")
    monkeypatch.setattr(cli_app, "download_genome", fake_download)
    result = runner.invoke(app, ["download-genome", "hg38", str(tmp_path), "--no-progress"])

    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["accession"] == "GCF_000001405.40"


def test_cli_deletes_other_characters_by_default_and_can_keep_them() -> None:
    deleted = runner.invoke(app, ["normalize", "ACX"])
    retained = runner.invoke(app, ["normalize", "ACX", "--keep-other"])

    assert deleted.exit_code == 0, deleted.output
    deleted_payload = json.loads(deleted.stdout)
    assert deleted_payload["normalized_parts"] == ["AC"]
    assert deleted_payload["changes"][-1]["operation"] == "delete_other"

    assert retained.exit_code == 2
    retained_payload = json.loads(retained.stdout)
    assert retained_payload["sequence"] is None
    assert retained_payload["invalid_symbols"][0]["symbol"] == "X"


def test_cli_convert_streams_fasta_to_json(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    target = tmp_path / "output.json"
    source.write_text(">one\nACGT\n>two\nNN\n", encoding="utf-8")

    result = runner.invoke(
        app,
        ["convert", str(source), str(target), "--no-progress"],
    )

    assert result.exit_code == 0, result.output
    assert [item["id"] for item in json.loads(target.read_text(encoding="utf-8"))] == [
        "one",
        "two",
    ]
    assert json.loads(result.stdout)["record_count"] == 2


def test_cli_report_writes_bounded_self_contained_html(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    target = tmp_path / "report.html"
    source.write_text(">one example\nACGT\n>two\nNN\n", encoding="utf-8")

    result = runner.invoke(
        app,
        ["report", str(source), str(target), "--title", "Local DNA report"],
    )

    assert result.exit_code == 0, result.output
    document = target.read_text(encoding="utf-8")
    assert document.startswith("<!doctype html>")
    assert "Local DNA report" in document
    assert "https://" not in document
    assert "record_count" in document
    payload = json.loads(result.stdout)
    assert payload["format"] == "html"
    assert payload["target_artifact"]["media_type"] == "text/html"

    refused = runner.invoke(app, ["report", str(source), str(target)])
    assert refused.exit_code == 2
    assert "Local DNA report" in target.read_text(encoding="utf-8")


def test_cli_report_enforces_record_limit_before_rendering(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    target = tmp_path / "report.html"
    source.write_text(">one\nA\n>two\nT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        ["report", str(source), str(target), "--max-records", "1"],
    )

    assert result.exit_code == 2
    assert not target.exists()


def test_cli_workflow_dry_run_uses_unified_entrypoint(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nACGT\n", encoding="utf-8")
    config = tmp_path / "workflow.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": "dnakit-workflow-v1",
                "run_id": "unified-cli",
                "input": {"path": "input.fasta", "format": "fasta"},
                "output_dir": "workflow-output",
                "steps": [
                    {
                        "id": "described",
                        "operation": "descriptors",
                        "input": "input",
                        "params": {"metrics": ["length"]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    result = runner.invoke(app, ["workflow", str(config), "--dry-run", "--no-progress"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)
    assert payload["status"] == "dry-run"
    assert payload["run_id"] == "unified-cli"
    assert not (tmp_path / "workflow-output").exists()


def test_cli_deduplicate_reverse_complement(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    target = tmp_path / "deduplicated.fasta"
    source.write_text(">forward\nAAGC\n>reverse\nGCTT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "deduplicate",
            str(source),
            str(target),
            "--equivalence",
            "reverse_complement",
        ],
    )

    assert result.exit_code == 0, result.output
    assert target.read_text(encoding="utf-8") == ">forward\nAAGC\n"
    assert json.loads(result.stdout)["removed_count"] == 1


def test_cli_split_writes_subsets_and_manifest(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    output_dir = tmp_path / "split"
    source.write_text(
        "".join(f">r{index}\nACGT\n" for index in range(4)),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(output_dir),
            "--ratios",
            "train=0.5,test=0.5",
            "--seed",
            "7",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "train.fasta").exists()
    assert (output_dir / "test.fasta").exists()
    manifest = json.loads((output_dir / "assignments.json").read_text(encoding="utf-8"))
    assert manifest["counts"] == {"test": 2, "train": 2}


def test_cli_refuses_to_overwrite_convert_output(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    target = tmp_path / "output.fasta"
    source.write_text(">one\nA\n", encoding="utf-8")
    target.write_text("keep", encoding="utf-8")

    result = runner.invoke(app, ["convert", str(source), str(target), "--no-progress"])

    assert result.exit_code == 2
    assert target.read_text(encoding="utf-8") == "keep"


def test_cli_split_rejects_unsafe_output_suffix(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    output_dir = tmp_path / "split"
    source.write_text(">one\nA\n>two\nT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(output_dir),
            "--ratios",
            "train=0.5,test=0.5",
            "--output-format",
            "../outside",
        ],
    )

    assert result.exit_code == 2
    assert not (tmp_path / "outside").exists()
    assert list(output_dir.glob(".*.tmp")) == []


@pytest.mark.parametrize(
    "unsafe_name",
    ["../escaped", "/tmp/escaped", "nested/name", r"nested\name", ".", "two words"],
)
def test_cli_split_rejects_unsafe_or_traversing_names(tmp_path: Path, unsafe_name: str) -> None:
    source = tmp_path / "input.fasta"
    output_dir = tmp_path / "split"
    source.write_text(">one\nA\n>two\nT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(output_dir),
            "--ratios",
            f"{unsafe_name}=0.5,train=0.5",
        ],
    )

    assert result.exit_code == 2
    assert not (tmp_path / "escaped.fasta").exists()
    assert not (output_dir / "assignments.json").exists()


def test_cli_split_rejects_duplicate_names(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    source.write_text(">one\nA\n>two\nT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(tmp_path / "split"),
            "--ratios",
            "train=0.5,train=0.5",
        ],
    )

    assert result.exit_code == 2
    assert not (tmp_path / "split" / "assignments.json").exists()


def test_cli_split_rejects_manifest_name_collision(tmp_path: Path) -> None:
    source = tmp_path / "input.fasta"
    output_dir = tmp_path / "split"
    source.write_text(">one\nA\n>two\nT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(output_dir),
            "--ratios",
            "assignments=0.5,test=0.5",
            "--output-format",
            "json",
            "--overwrite",
        ],
    )

    assert result.exit_code == 2
    assert list(output_dir.iterdir()) == []


def test_cli_split_writer_failure_leaves_no_partial_files(tmp_path: Path) -> None:
    source = tmp_path / "input.json"
    output_dir = tmp_path / "split"
    source.write_text(
        json.dumps(
            [
                {
                    "id": "with-quality",
                    "sequence": "A",
                    "letter_annotations": {"phred_quality": [40]},
                },
                {"id": "without-quality", "sequence": "T"},
            ]
        ),
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "split",
            str(source),
            str(output_dir),
            "--ratios",
            "train=0.5,test=0.5",
            "--output-format",
            "fastq",
            "--seed",
            "0",
        ],
    )

    assert result.exit_code == 2
    assert list(output_dir.iterdir()) == []
