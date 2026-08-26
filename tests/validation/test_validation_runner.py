from __future__ import annotations

import json
from pathlib import Path

import pytest
from validation.run_validation import build_report, write_report


def test_validation_report_has_manual_checks_and_explicit_nupack_boundary() -> None:
    report = build_report(include_optional=False, show_progress=False)

    assert report["schema_version"] == "dnakit.validation.v1"
    assert report["summary"]["fail"] == 0
    assert report["summary"]["pass"] >= 6
    assert report["prohibited_backend_audit"] == {
        "primer3_automatic_discovery_attempted": False,
        "primer3_installation_attempted": False,
        "primer3_import_attempted": False,
        "primer3_call_attempted": False,
        "nupack_installation_attempted": False,
        "nupack_probe_attempted": False,
        "nupack_import_attempted": False,
        "nupack_call_attempted": False,
        "note": (
            "Primer3 and NUPACK are deliberately outside this runner. These booleans "
            "describe runner behavior, not a scan of the environment."
        ),
    }
    assert {check["id"] for check in report["checks"]} >= {
        "MANUAL-001",
        "MANUAL-005",
        "ALGORITHM-001",
    }


def test_optional_validation_does_not_discover_or_execute_primer3() -> None:
    report = build_report(include_optional=True, show_progress=False)
    check_ids = {check["id"] for check in report["checks"]}

    assert not any(check_id.startswith("PRIMER3-") for check_id in check_ids)
    assert report["environment"]["primer3_cli"] == "not automatically discovered or executed"


def test_installed_optional_backends_pass_only_comparable_checks() -> None:
    report = build_report(include_optional=True, show_progress=False)
    checks = {check["id"]: check for check in report["checks"]}
    if "BIOPYTHON-000" in checks:
        pytest.skip("Biopython is not installed")

    assert report["summary"]["fail"] == 0
    assert checks["BIOPYTHON-001"]["status"] == "pass"
    assert checks["BIOPYTHON-002"]["status"] == "pass"
    assert checks["BIOPYTHON-003"]["status"] == "pass"
    assert checks["BIOPYTHON-004"]["status"] == "pass"
    assert checks["BIOCLUSTER-001"]["status"] == "pass"


def test_validation_write_refuses_accidental_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "validation.json"
    write_report({"ok": True}, output)
    assert json.loads(output.read_text(encoding="utf-8")) == {"ok": True}
    with pytest.raises(FileExistsError):
        write_report({"ok": True}, output)
