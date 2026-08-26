from __future__ import annotations

import json
from pathlib import Path

import pytest
from benchmarks.benchmark_core import run_benchmarks, write_report


def test_small_benchmark_records_parameters_environment_and_samples() -> None:
    report = run_benchmarks(
        sizes=(32,),
        repeats=2,
        warmups=0,
        seed=7,
        tasks=("construct", "gc_content"),
        show_progress=False,
    )

    assert report["schema_version"] == "dnakit.benchmark.v1"
    assert report["parameters"]["seed"] == 7
    assert report["environment"]["dnakit"]
    assert len(report["cases"]) == 2
    assert all(len(case["samples"]) == 2 for case in report["cases"])
    assert all(case["duration_ns"]["minimum"] >= 0 for case in report["cases"])
    assert all(case["peak_tracemalloc_bytes"]["maximum"] >= 0 for case in report["cases"])
    assert {case["implementation"] for case in report["cases"]} == {"dnakit"}
    assert report["source_code_metrics"]["dnakit"]["entries"]


def test_biopython_comparison_records_matched_tasks_versions_and_source_metrics() -> None:
    pytest.importorskip("Bio")

    report = run_benchmarks(
        sizes=(16, 32),
        repeats=1,
        warmups=0,
        seed=9,
        tasks=("construct", "gc_content", "reverse_complement", "minhash"),
        implementations=("dnakit", "biopython"),
        show_progress=False,
    )

    assert report["environment"]["dnakit"]
    assert report["parameters"]["implementations"] == ["dnakit", "biopython"]
    assert len(report["cases"]) == 14
    assert {case["implementation"] for case in report["cases"]} == {
        "dnakit",
        "biopython",
    }
    assert report["parameters"]["comparison_task_policy"]["biopython"] == [
        "construct",
        "gc_content",
        "reverse_complement",
    ]
    assert report["source_code_metrics"]["biopython"]["total_nonblank_noncomment_lines"] > 0


def test_benchmark_work_limits_are_enforced() -> None:
    with pytest.raises(ValueError, match="safety ceiling"):
        run_benchmarks(
            sizes=(200_001,),
            repeats=1,
            warmups=0,
            seed=0,
            tasks=("construct",),
            show_progress=False,
        )


def test_report_write_is_atomic_and_refuses_unrequested_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "report.json"
    report = {"finite": 1.0}

    write_report(report, output)
    assert json.loads(output.read_text(encoding="utf-8")) == report
    with pytest.raises(FileExistsError, match="already exists"):
        write_report(report, output)
    write_report({"finite": 2.0}, output, force=True)
    assert json.loads(output.read_text(encoding="utf-8")) == {"finite": 2.0}
