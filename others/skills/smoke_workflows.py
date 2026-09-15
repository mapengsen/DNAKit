"""Run small offline fixtures on installed skills and record the exact evidence.

This checks three scripts and their packaged resources, not complete scientific workflows.
Use Python 3.11+ with DNAKit installed. All generated files stay under temp/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from dnakit.skills import get_skill, install_skill

ROOT = Path(__file__).resolve().parents[2]


def fingerprint(entry: dict) -> str:
    return hashlib.sha256(json.dumps(entry["files"], sort_keys=True).encode()).hexdigest()


def check_fasta_curator(python: str, directory: Path, run: Path) -> list[dict]:
    fixture = run / "duplicates.fasta"
    fixture.write_text(
        ">same id\naacg\n>same id\nTTTT\n>same_sequence\nAACG\n>unique\nGGTT\n",
        encoding="utf-8",
    )
    output = run / "curated.fasta"
    args = [
        python,
        "-B",
        str(directory / "scripts" / "curate_fasta.py"),
        str(fixture),
        "--output",
        str(output),
        "--deduplicate",
        "both",
    ]
    completed = subprocess.run(args, cwd=run, check=True, capture_output=True, text=True)
    report = json.loads(completed.stdout)
    assert output.read_text() == ">same_id\nAACG\n>unique\nGGTT\n"
    assert (report["input_records"], report["retained_records"], report["removed_records"]) == (
        4,
        2,
        2,
    )
    assert report["duplicates_by_reason"] == {"id": 1, "sequence": 1}
    with (run / "curated.fasta.mapping.tsv").open() as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["duplicate_reason"] for row in mapping] == [
        "",
        "duplicate_id",
        "duplicate_sequence",
        "",
    ]
    before = {path: path.read_bytes() for path in run.glob("curated.fasta*")}
    repeated = subprocess.run(args, cwd=run, capture_output=True, text=True)
    assert repeated.returncode != 0 and "refusing to overwrite" in repeated.stderr
    assert all(path.read_bytes() == data for path, data in before.items())
    fixture.write_text(">empty_sequence\n", encoding="utf-8")
    invalid_output = run / "invalid.fasta"
    invalid_args = [str(invalid_output) if arg == str(output) else arg for arg in args]
    invalid = subprocess.run(invalid_args, cwd=run, capture_output=True, text=True)
    assert invalid.returncode != 0 and "has no sequence" in invalid.stderr
    assert not invalid_output.exists()
    return [
        {
            "case": "four-record-fasta-id-and-sequence-deduplication",
            "input_records": 4,
            "retained_records": 2,
            "removed_by_id": 1,
            "removed_by_sequence": 1,
            "mapping_and_report": "verified",
        },
        {"case": "existing-output-preserved", "result": "refused overwrite"},
        {"case": "empty-sequence-rejected", "result": "no output created"},
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python", default=sys.executable, help="Python 3.11+ for upstream scripts"
    )
    python = parser.parse_args().python
    version = json.loads(
        subprocess.check_output(
            [python, "-c", "import sys,json; print(json.dumps(list(sys.version_info[:3])))"],
            text=True,
        )
    )
    if tuple(version) < (3, 11):
        raise SystemExit(
            "Use Python 3.11+ to match the genomic-coordinates compatibility declaration"
        )
    work = ROOT / "temp" / "dnakit-skills-work"
    work.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="offline-examples-", dir=work))
    installed = run / "skills"
    records = []
    for identifier in (
        "genomic-coordinates",
        "sample-qc-triage",
        "bio-fasta-database-curator",
    ):
        result = install_skill(identifier, destination=installed)
        directory = Path(result["path"])
        cases = []
        if identifier == "genomic-coordinates":
            script = directory / "scripts" / "convert_coords.py"
            for source, target, start, end, expected in (
                ("bed", "gff", 999, 1000, (1000, 1000, 1)),
                ("gff", "bed", 1000, 1000, (999, 1000, 1)),
                ("bed", "gff", 0, 10, (1, 10, 10)),
            ):
                args = [
                    python,
                    "-B",
                    str(script),
                    "--from",
                    source,
                    "--to",
                    target,
                    "chr1",
                    str(start),
                    str(end),
                    "--format",
                    "json",
                ]
                completed = subprocess.run(
                    args, cwd=run, check=True, capture_output=True, text=True
                )
                rows = json.loads(completed.stdout)
                row = rows[0]
                assert (row["start"], row["end"], row["length"]) == expected, row
                cases.append(
                    {
                        "case": f"{source}-{target}-{start}-{end}",
                        "expected": list(expected),
                        "actual": row,
                    }
                )
        elif identifier == "bio-fasta-database-curator":
            cases = check_fasta_curator(python, directory, run)
        else:
            fixture = run / "qc.csv"
            fixture.write_text(
                "sample_id,batch,total_reads,mapped_pct,duplicate_pct,mitochondri"
                "al_pct,contamination_pct,complexity_score\n"
                "clean,A,100000,99,10,1,0.1,0.95\n"
                "contaminated,A,100000,99,10,1,8,0.95\n",
                encoding="utf-8",
            )
            output = run / "qc-results"
            args = [
                python,
                "-B",
                str(directory / "sample_qc_triage.py"),
                "--input",
                str(fixture),
                "--output",
                str(output),
            ]
            subprocess.run(args, cwd=run, check=True, capture_output=True, text=True)
            data = json.loads((output / "result.json").read_text())
            assert data["summary"]["sample_count"] == 2
            assert data["summary"]["flagged_count"] == 1
            assert data["samples"][0]["issues"] == []
            assert data["samples"][1]["issues"] == ["contamination"]
            assert (output / "report.md").is_file()
            assert (output / "tables" / "sample_flags.csv").is_file()
            cases.append(
                {
                    "case": "two-synthetic-samples",
                    "sample_count": 2,
                    "flagged_count": 1,
                    "expected_flag": "contamination",
                    "note": "Operational QC threshold fixture, not a clinical evaluation.",
                }
            )
        entry = get_skill(identifier)
        records.append(
            {
                "id": identifier,
                "upstream_revision": entry["upstream_revision"],
                "resource_manifest_sha256": fingerprint(entry),
                "execution": "smoke_passed",
                "cases": cases,
            }
        )
    report = {
        "schema_version": 1,
        "date": "2026-09-12",
        "python": ".".join(map(str, version)),
        "scope": (
            "Offline subprocess fixtures after installing complete selected "
            "folders outside the source layout. No LLM, API, large dataset, "
            "or full-pipeline validation."
        ),
        "records": records,
    }
    (work / "offline-validation.json").write_text(json.dumps(report, indent=2) + "\n")
    (ROOT / "registry" / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    path = ROOT / "registry" / "catalog.json"
    catalog = json.loads(path.read_text())
    for entry in catalog["skills"]:
        record = next((record for record in records if record["id"] == entry["id"]), None)
        if record and fingerprint(entry) == record["resource_manifest_sha256"]:
            entry["validation"]["execution"] = "smoke_passed"
            entry["validation"]["evidence"] = "registry/validation.json"
            entry["validation"]["scope"] = report["scope"]
    path.write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n")
    print(
        json.dumps(
            {"validated": [x["id"] for x in records], "output_directory": str(run)}, indent=2
        )
    )


if __name__ == "__main__":
    main()
