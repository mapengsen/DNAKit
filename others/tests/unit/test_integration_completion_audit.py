"""Regression tests for completion record integrity, not scientific validation.

Version note: exercises count/row agreement, nested files, export targets and
complete non-skipped callable mappings using isolated synthetic records.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import types
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/integration/audit_completion.py"
)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("completion_audit_under_test", SCRIPT)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    base = tmp_path / "coverage"
    for folder in ("integration/packs", "integration/disposition", "coverage_gap/manifests"):
        (base / folder).mkdir(parents=True)
    module = tmp_path / "src/dnakit/area/functions.py"
    module.parent.mkdir(parents=True)
    module.write_text("def tool(*, verbose=True): return {}" + chr(10))
    (tmp_path / "nested.md").write_text("Nested documentation")
    (tmp_path / "nested_test.py").write_text("# Test fixture")
    row = {
        "group": "G1",
        "source": "tool",
        "operation": "operation",
        "status": "implemented",
        "dnakit_function": "dnakit.area.functions.tool",
        "reason": "Direct operation",
    }
    manifest = {
        "modules": [
            {
                "module": "dnakit.area.functions",
                "exports": ["tool"],
                "docs": ["nested.md"],
                "tests": ["nested_test.py"],
            }
        ],
        "operation_count": 1,
        "status_counts": {"implemented": 1, "skipped": 0},
        "operation_rows": [dict(row)],
    }
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    monkeypatch.setattr(audit, "BASE", base)
    monkeypatch.setattr(audit, "PACKS", (1,))

    def write(rows=None, data=None):
        (base / "integration/packs/P01.json").write_text(
            json.dumps({"items": [{"group": "G1", "source": "tool", "operation": "operation"}]})
        )
        with (base / "integration/disposition/p01.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerows(rows if rows is not None else [row])
        (base / "coverage_gap/manifests/integration_p01.json").write_text(
            json.dumps(data if data is not None else manifest)
        )

    write()
    return audit, row, manifest, write


def test_complete_nested_manifest_passes(fixture):
    audit, _row, _manifest, _write = fixture
    report = audit.audit()
    assert report["record_integrity_passed"]
    assert report["recorded_items"] == 1
    assert report["callables_static_checked"] == 1
    assert report["callables_import_checked"] == 0
    assert "do not prove scientific agreement" in report["scientific_validation"]


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("operation_count", 2, "operation_count differs"),
        ("operation_count", True, "operation_count differs"),
        ("status_counts", {"implemented": 0, "merged": 1}, "status_counts differs"),
        ("status_counts", {"implemented": True}, "status_counts differs"),
        ("operation_rows", [], "operation groups differ"),
        (
            "operation_rows",
            [{"group": "G1", "status": "merged", "dnakit_function": "dnakit.area.functions.tool"}],
            "status differs",
        ),
    ],
)
def test_manifest_count_and_operation_rows_match_csv(fixture, field, value, message):
    audit, _row, manifest, write = fixture
    manifest[field] = value
    write()
    assert any(message in error for error in audit.audit()["errors"])


@pytest.mark.parametrize("key", ["docs", "tests", "supporting_files"])
def test_nested_module_files_are_checked(fixture, key):
    audit, _row, manifest, write = fixture
    manifest["modules"][0][key] = ["missing.file"]
    write()
    assert any("missing " + key + " file" in error for error in audit.audit()["errors"])


def test_multiple_module_levels_are_traversed(fixture):
    audit, _row, manifest, write = fixture
    child = manifest["modules"][0]
    manifest["modules"] = [{"modules": [child]}]
    child["tests"] = ["missing.deep.test.py"]
    write()
    report = audit.audit()
    assert report["callables_static_checked"] == 1
    assert any(
        "modules[0]/modules[0]" in error and "missing tests" in error for error in report["errors"]
    )


def test_truncated_csv_is_reported_without_crashing(fixture):
    audit, _row, _manifest, _write = fixture
    path = audit.BASE / "integration/disposition/p01.csv"
    path.write_text(
        "group,source,operation,status,dnakit_function,reason" + chr(10) + "G1,tool" + chr(10)
    )
    assert any("missing or excess columns" in error for error in audit.audit()["errors"])


def test_export_name_must_exist_in_source(fixture):
    audit, _row, manifest, write = fixture
    manifest["modules"][0]["exports"].append("absent")
    write()
    assert any("target name missing" in error for error in audit.audit()["errors"])


def test_new_mapping_must_be_in_manifest_exports(fixture):
    audit, _row, manifest, write = fixture
    manifest["modules"][0]["exports"] = []
    write()
    assert any("absent from manifest exports" in error for error in audit.audit()["errors"])


@pytest.mark.parametrize("target", ["", "tool", "dnakit.area.functions.tool()"])
def test_non_skipped_mapping_requires_one_dotted_target(fixture, target):
    audit, row, manifest, write = fixture
    row["dnakit_function"] = target
    manifest.pop("operation_rows")
    write()
    assert any("without a valid dotted callable" in error for error in audit.audit()["errors"])


def test_already_covered_target_is_checked_without_new_export(fixture):
    audit, row, manifest, write = fixture
    row["status"] = "already_covered"
    manifest["modules"][0]["exports"] = []
    manifest["status_counts"] = {"already_covered": 1}
    manifest.pop("operation_rows")
    write()
    assert audit.audit()["record_integrity_passed"]


def test_import_check_rejects_noncallable_and_missing_keyword_verbose(fixture, monkeypatch):
    audit, _row, _manifest, _write = fixture
    for value in (42, lambda verbose=True: None):
        monkeypatch.setattr(
            audit.importlib,
            "import_module",
            lambda name, value=value: types.SimpleNamespace(tool=value),
        )
        report = audit.audit(check_imports=True)
        assert not report["record_integrity_passed"]
        assert report["callables_import_checked"] == 0

    def valid(*, verbose=True):
        return {}

    monkeypatch.setattr(
        audit.importlib, "import_module", lambda name: types.SimpleNamespace(tool=valid)
    )
    assert audit.audit(check_imports=True)["callables_import_checked"] == 1


def test_duplicate_csv_groups_fail(fixture):
    audit, row, _manifest, write = fixture
    write(rows=[row, row])
    assert any("duplicate group IDs" in error for error in audit.audit()["errors"])


def test_skipped_reason_and_manifest_row_are_checked(fixture):
    audit, row, manifest, write = fixture
    row.update(status="skipped", dnakit_function="", reason="")
    manifest.pop("operation_rows")
    manifest["status_counts"] = {"skipped": 1}
    manifest["skipped"] = [{"group": "G1", "reason": "different reason"}]
    write()
    errors = audit.audit()["errors"]
    assert any("skipped without reason" in error for error in errors)
    assert any("reason differs" in error for error in errors)
