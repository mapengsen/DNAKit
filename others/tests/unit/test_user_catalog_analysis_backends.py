"""Optional genuine backend smoke tests on bounded synthetic biological fixtures.

Set DNAKIT_ANALYSIS_TOOLS to the directory containing separately installed real
executables. These tests exercise upstream programs, never substitute outputs,
and assess basic invariants rather than accuracy on real biological datasets.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import random
from pathlib import Path

import pytest

from dnakit.integrations.execution import run_analysis

CODONS = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


def _backend(name: str) -> str:
    directory = os.environ.get("DNAKIT_ANALYSIS_TOOLS")
    if not directory:
        pytest.skip("Set DNAKIT_ANALYSIS_TOOLS to separately installed genuine backends")
    path = Path(directory) / name
    if not path.is_file() or not os.access(path, os.X_OK):
        pytest.skip(f"Real optional executable is not installed: {name}")
    return str(path)


def _run(tmp_path: Path, analysis: str, program: str, inputs: dict, **kwargs) -> dict:
    result = run_analysis(
        analysis,
        inputs,
        str(tmp_path / analysis),
        executable_paths={program: _backend(program)},
        threads=1,
        timeout_seconds=120,
        max_artifact_bytes=100_000_000,
        verbose=False,
        **kwargs,
    )
    assert result["status"] == "completed"
    # Persistent evidence associates commands with exact synthetic inputs and
    # outputs; no hashes are interpreted as scientific validity.
    output = Path(result["output_directory"])
    evidence = {
        "analysis": analysis,
        "project_id": result["project_id"],
        "fixture": "deterministic synthetic smoke; scientific accuracy not assessed",
        "commands": [s["arguments"] for s in result["steps"]],
        "input_sha256": {
            name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
            for name, path in inputs.items()
            if Path(path).is_file()
        },
        "output_sha256": {
            str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in output.rglob("*")
            if path.is_file() and path.name != "dnakit-run.json"
        },
    }
    (output / "smoke-evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    return result


@pytest.fixture
def proteins(tmp_path: Path) -> tuple[str, Path, Path, Path]:
    rng = random.Random(335)
    sequence = "M" + "".join(rng.choice(tuple(CODONS)) for _ in range(239))
    query = tmp_path / "query protein.fa"
    reference = tmp_path / "reference protein.fa"
    genome = tmp_path / "genome sequence.fa"
    query.write_text(f">query\n{sequence}\n")
    reference.write_text(f">reference\n{sequence}\n")
    dna = "".join(CODONS[aa] for aa in sequence)
    genome.write_text(f">chromosome\n{'ACGT' * 25}{dna}TAA{'TGCA' * 25}\n")
    return sequence, query, reference, genome


@pytest.mark.parametrize("mode", ["blastp", "blastx"])
def test_real_diamond_exact_match(tmp_path: Path, proteins: tuple, mode: str) -> None:
    sequence, query, reference, genome = proteins
    result = _run(
        tmp_path,
        f"diamond-{mode}",
        "diamond",
        {
            "query": str(query if mode == "blastp" else genome),
            "reference_proteins": str(reference),
        },
    )
    rows = [
        line.split("\t")
        for line in (Path(result["output_directory"]) / "hits.tsv").read_text().splitlines()
    ]
    assert rows and all(len(row) == 12 for row in rows)
    assert any(float(row[2]) == 100 and int(row[3]) >= 0.95 * len(sequence) for row in rows)
    assert rows[0][1] == "reference"


def test_real_miniprot_coding_alignment_coordinates(tmp_path: Path, proteins: tuple) -> None:
    sequence, query, _, genome = proteins
    result = _run(
        tmp_path, "miniprot-align", "miniprot", {"genome": str(genome), "proteins": str(query)}
    )
    lines = (Path(result["output_directory"]) / "alignment.gff3").read_text().splitlines()
    features = [line.split("\t") for line in lines if not line.startswith("#")]
    cds = [row for row in features if row[2] == "CDS"]
    assert cds
    assert all(row[0] == "chromosome" and row[6] == "+" for row in cds)
    assert min(int(row[3]) for row in cds) == 101
    stops = [row for row in features if row[2] == "stop_codon"]
    assert len(stops) == 1
    # miniprot GFF includes the terminal stop codon in CDS coordinates.
    assert sum(int(row[4]) - int(row[3]) + 1 for row in cds) == 3 * (len(sequence) + 1)
    assert int(stops[0][3]) == 101 + 3 * len(sequence)


@pytest.mark.parametrize("analysis", ["vsearch-dereplicate", "vsearch-cluster"])
def test_real_vsearch_duplicates_and_cluster_membership(tmp_path: Path, analysis: str) -> None:
    rng = random.Random(913)
    first = "".join(rng.choice("ACGT") for _ in range(240))
    second = "".join(rng.choice("ACGT") for _ in range(240))
    sequences = tmp_path / "nucleotide reads.fa"
    sequences.write_text(f">a\n{first}\n>a_copy\n{first}\n>b\n{second}\n")
    result = _run(tmp_path, analysis, "vsearch", {"sequences": str(sequences)})
    output = Path(result["output_directory"])
    fasta = output / (
        "dereplicated.fasta" if analysis.endswith("dereplicate") else "centroids.fasta"
    )
    headers = [line for line in fasta.read_text().splitlines() if line.startswith(">")]
    assert len(headers) == 2
    if analysis.endswith("dereplicate"):
        assert sorted(int(line.split("size=")[1].rstrip(";")) for line in headers) == [1, 2]
    rows = [line.split("\t") for line in (output / "clusters.uc").read_text().splitlines()]
    assert sum(row[0] in ("S", "H") for row in rows) == 3


def test_real_mmseqs_protein_search(tmp_path: Path, proteins: tuple) -> None:
    sequence, query, reference, _ = proteins
    result = _run(
        tmp_path,
        "mmseqs-protein-search",
        "mmseqs",
        {"query": str(query), "reference": str(reference)},
    )
    rows = [
        line.split("\t")
        for line in (Path(result["output_directory"]) / "hits.tsv").read_text().splitlines()
    ]
    assert rows
    assert any(row[:2] == ["query", "reference"] and int(row[3]) == len(sequence) for row in rows)


def test_real_mmseqs_protein_cluster(tmp_path: Path, proteins: tuple) -> None:
    sequence, _, _, _ = proteins
    fasta = tmp_path / "proteins to cluster.fa"
    fasta.write_text(f">a\n{sequence}\n>b\n{sequence}\n")
    result = _run(tmp_path, "mmseqs-protein-cluster", "mmseqs", {"proteins": str(fasta)})
    rows = [
        line.split("\t")
        for line in (Path(result["output_directory"]) / "clusters_cluster.tsv")
        .read_text()
        .splitlines()
    ]
    assert {row[1] for row in rows} == {"a", "b"}
    assert len({row[0] for row in rows}) == 1


def test_real_drugz_unpaired_finite_gene_scores(tmp_path: Path) -> None:
    rng = random.Random(335)
    counts = tmp_path / "guide counts.tsv"
    with counts.open("w") as handle:
        handle.write("sgRNA\tgene\tcontrol1\tcontrol2\ttreated1\ttreated2\n")
        for index in range(1200):
            values = [rng.randint(80, 1500) for _ in range(4)]
            handle.write(f"guide{index}\tgene{index // 24}\t" + "\t".join(map(str, values)) + "\n")
    result = _run(
        tmp_path,
        "drugz-unpaired",
        "drugz.py",
        {"counts": str(counts)},
        parameters={
            "control": "control1,control2",
            "treatment": "treated1,treated2",
            "half_window_size": 50,
        },
    )
    with (Path(result["output_directory"]) / "drugz.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 50
    assert {row["GENE"] for row in rows} == {f"gene{index}" for index in range(50)}
    for row in rows:
        assert math.isfinite(float(row["normZ"]))
        for column in ("pval_synth", "pval_supp", "fdr_synth", "fdr_supp"):
            assert 0 <= float(row[column]) <= 1
