"""Checks the disposition of the first 30 kipoi precondition-skip rows.

No row in this batch performs a DNA or model-prediction operation, so nothing is implemented.
These tests check the recorded decisions: each row names a real object in the installed kipoi,
every skip reason starts with an allowed label, and the manifest matches the CSV.
"""

from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path

import pytest

pytest.importorskip("kipoi")

REPO = Path(__file__).resolve().parents[3]
COVERAGE = REPO / "paper/write_paper" / "figure" / "fig2" / "task_coverage" / "coverage_gap"
CSV_PATH = COVERAGE / "waves" / "batch_retry_F_kipoi_a_disposition.csv"
MANIFEST_PATH = COVERAGE / "manifests" / "kipoi_precondition_a.json"
ALLOWED_LABELS = ("(a)", "(b)", "(c)", "(d)")


def _rows() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _resolve(qualname: str) -> object:
    parts = qualname.split(".")
    for cut in range(len(parts), 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:cut]))
        except ImportError:
            continue
        try:
            for attr in parts[cut:]:
                obj = getattr(obj, attr)
        except AttributeError:
            return None
        return obj
    return None


def test_csv_has_header_and_thirty_skipped_kipoi_rows() -> None:
    rows = _rows()
    assert len(rows) == 30
    assert set(rows[0]) == {"package", "qualname", "work_type", "status", "reason", "module"}
    assert {row["package"] for row in rows} == {"kipoi"}
    assert {row["status"] for row in rows} == {"skipped"}
    assert {row["module"] for row in rows} == {""}


def test_every_reason_starts_with_an_allowed_label() -> None:
    for row in _rows():
        assert row["reason"].startswith(ALLOWED_LABELS), row["qualname"]


def test_every_row_resolves_in_installed_kipoi() -> None:
    missing = [row["qualname"] for row in _rows() if _resolve(row["qualname"]) is None]
    assert missing == []


def test_label_counts_match_the_recorded_decisions() -> None:
    labels = [row["reason"][:3] for row in _rows()]
    assert labels.count("(b)") == 2
    assert labels.count("(c)") == 14
    assert labels.count("(d)") == 14


def test_manifest_matches_csv() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["module"] is None
    assert manifest["exports"] == []
    assert manifest["new_public_names"] == 0
    expected = [{"operation": row["qualname"], "reason": row["reason"]} for row in _rows()]
    assert manifest["skipped"] == expected
