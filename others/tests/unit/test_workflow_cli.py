"""Tests for the standalone local workflow CLI."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dnakit.cli.workflow import main


def _config(tmp_path: Path) -> Path:
    (tmp_path / "input.fasta").write_text(">one\nACGT\n", encoding="utf-8")
    path = tmp_path / "workflow.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "dnakit-workflow-v1",
                "run_id": "cli-run",
                "input": {"path": "input.fasta", "format": "fasta"},
                "output_dir": "output",
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
    return path


def test_workflow_cli_dry_run_outputs_json_without_writes(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)

    code = main(["run", str(config), "--dry-run", "--no-progress"])

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert code == 0
    assert payload["status"] == "dry-run"
    assert not (tmp_path / "output").exists()


def test_workflow_cli_reports_invalid_workflow(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["unknown"] = True
    config.write_text(json.dumps(payload), encoding="utf-8")

    code = main(["run", str(config), "--no-progress"])

    captured = capsys.readouterr()
    assert code == 2
    assert "UNKNOWN_WORKFLOW_FIELD" in captured.err
