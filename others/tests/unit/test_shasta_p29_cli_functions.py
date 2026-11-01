from __future__ import annotations

import json
import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit.assembly import shasta_p29_cli_functions as sh
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

SHASTA = os.environ.get("DNAKIT_SHASTA") or shutil.which("shasta")
needs_shasta = pytest.mark.skipif(SHASTA is None, reason="Shasta is not installed")

FAKE = """#!{python}
import sys, pathlib, json
args = sys.argv[1:]
if args == ["--version"]:
    print("Shasta Release 0.14.0"); sys.exit(0)
run = pathlib.Path(args[args.index("--assemblyDirectory") + 1])
run.mkdir()
pathlib.Path(run / "args.json").write_text(json.dumps(args))
(run / "Assembly.fasta").write_text(">0 length 8\\nACGTACGT\\n")
(run / "AssemblySummary.json").write_text(json.dumps({{"Assembled segments": {{"Number of segments assembled": 1,
    "Total assembled segment length": 8, "Assembled segments N50": 8}}}}))
"""


@pytest.fixture()
def fake(tmp_path: Path) -> str:
    path = tmp_path / "fake_shasta"
    path.write_text(FAKE.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def _reads(path: Path, count: int = 450, error: float = 0.01, seed: int = 1) -> Path:
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(150_000))
    table = str.maketrans("ACGT", "TGCA")
    with path.open("w") as out:
        for index in range(count):
            length = rng.randint(12_000, 25_000)
            start = rng.randint(0, len(genome) - length)
            read = genome[start : start + length]
            if rng.random() < 0.5:
                read = read.translate(table)[::-1]
            mutated = []
            for base in read:
                x = rng.random()
                if x < error / 3:
                    continue
                if x < 2 * error / 3:
                    mutated.append(base + rng.choice("ACGT"))
                elif x < error:
                    mutated.append(rng.choice("ACGT"))
                else:
                    mutated.append(base)
            out.write(f">r{index}\n{''.join(mutated)}\n")
    return path


def test_arguments_are_mapped_to_shasta_flags(tmp_path: Path, fake: str) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">a\nACGT\n")
    run = tmp_path / "run"
    result = sh.shasta_assemble(
        [str(reads)], str(run), config="Nanopore-ncm23-May2024", assembly_mode=3, align_method=4,
        read_graph_creation_method=4, strand_separation_method=0, palindromic_flagging=False, simplify_max_length=[10, 100],
        detangle_method=2, write_reads_by_segment=True, minhash_all_pairs=True, flag_inconsistent_alignments=False,
        options={"Align.maxSkip": 40, "MarkerGraph.maxDistance": 25}, threads=2, shasta_path=fake, verbose=False,
    )
    args = json.loads((run / "args.json").read_text())
    assert result["status"] == "completed" and result["version"] == "0.14.0"
    assert result["summary"]["assembled_bases"] == 8 and result["assembly_fasta"].endswith("Assembly.fasta")
    for flag, value in [("--Assembly.mode", "3"), ("--Align.alignMethod", "4"), ("--ReadGraph.creationMethod", "4"),
                        ("--ReadGraph.strandSeparationMethod", "0"), ("--MarkerGraph.simplifyMaxLength", "10,100"),
                        ("--Assembly.detangleMethod", "2"), ("--Align.maxSkip", "40"), ("--MarkerGraph.maxDistance", "25"),
                        ("--threads", "2"), ("--config", "Nanopore-ncm23-May2024")]:
        assert args[args.index(flag) + 1] == value
    for flag in ("--Reads.palindromicReads.skipFlagging", "--Assembly.writeReadsByAssembledSegment", "--MinHash.allPairs"):
        assert flag in args
    assert "--ReadGraph.flagInconsistentAlignments" not in args and "--Reads.minReadLength" not in args


def test_unset_parameters_are_not_passed(tmp_path: Path, fake: str) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">a\nACGT\n")
    sh.shasta_assemble([str(reads)], str(tmp_path / "run"), shasta_path=fake, verbose=False)
    args = json.loads((tmp_path / "run" / "args.json").read_text())
    assert not [a for a in args if a.startswith("--") and "." in a]


def test_invalid_inputs_are_refused(tmp_path: Path, fake: str) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">a\nACGT\n")
    kwargs = {"shasta_path": fake, "verbose": False}
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "r1"), options={"Not.anOption": 1}, **kwargs)
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "r2"), options={"Align.maxSkip": True}, **kwargs)
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "r3"), align_method=2, **kwargs)
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "r4"), align_method=4, options={"Align.alignMethod": 3}, **kwargs)
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "r5"), config="--evil", **kwargs)
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(tmp_path / "missing.fa")], str(tmp_path / "r6"), **kwargs)
    (tmp_path / "exists").mkdir()
    with pytest.raises(ConfigurationError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "exists"), **kwargs)


def test_missing_program(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">a\nACGT\n")
    monkeypatch.delenv("DNAKIT_SHASTA", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    if (Path(sys.executable).parent / "shasta").is_file():
        pytest.skip("a shasta executable sits next to the interpreter")
    with pytest.raises(BackendUnavailableError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "run"), verbose=False)
    with pytest.raises(BackendUnavailableError):
        sh.shasta_assemble([str(reads)], str(tmp_path / "run2"), shasta_path=str(tmp_path / "nope"), verbose=False)


def _synthetic_run(path: Path) -> Path:
    path.mkdir()
    (path / "ReadGraphComponents.csv").write_text("Component,RepresentingRead,OrientedReadCount\n0,0-0,108\n1,0-1,9\n")
    (path / "ReadSummary.csv").write_text("Id,Name,Palindromic,Chimeric,\n0,r0,No,Yes,\n1,r1,Yes,No,\n")
    (path / "Assembly.fasta").write_text(">0 length 6\nGGCCAT\n>1 length 4\nACGT\n")
    (path / "Assembly-BothStrands-NoSequence.gfa").write_text("H\tVN:Z:1.0\nS\t0\t*\tLN:i:10\nS\t1\t*\tLN:i:30\nL\t0\t+\t1\t+\t0M\n")
    (path / "Assembly.gfa").write_text("H\tVN:Z:1.0\nS\t0\tACGTAC\n")
    return path


def test_read_run_table(tmp_path: Path) -> None:
    run = _synthetic_run(tmp_path / "run")
    table = sh.shasta_read_run_table(str(run), "read_graph_components", verbose=False)
    assert table["row_count"] == 2 and table["rows"][0] == {"Component": 0, "RepresentingRead": "0-0", "OrientedReadCount": 108}
    summary = sh.shasta_read_run_table(str(run), "read_summary", max_rows=1, verbose=False)
    assert summary["columns"] == ["Id", "Name", "Palindromic", "Chimeric"] and summary["truncated"] is True
    with pytest.raises(ConfigurationError):
        sh.shasta_read_run_table(str(run), "journeys", verbose=False)
    with pytest.raises(ConfigurationError):
        sh.shasta_read_run_table(str(run), "nothing", verbose=False)


def test_read_assembly_output(tmp_path: Path) -> None:
    run = _synthetic_run(tmp_path / "run")
    fasta = sh.shasta_read_assembly_output(str(run), verbose=False)
    assert fasta["record_count"] == 2 and fasta["total_length"] == 10 and fasta["n50"] == 6
    assert fasta["records"][0] == {"name": "0", "length": 6, "gc_fraction": 0.666667}
    gfa = sh.shasta_read_assembly_output(str(run), output="gfa_no_sequence", verbose=False)
    assert gfa["record_count"] == 2 and gfa["total_length"] == 40 and gfa["link_count"] == 1 and gfa["n50"] == 30
    seq = sh.shasta_read_assembly_output(str(run), output="gfa", include_sequences=True, verbose=False)
    assert seq["records"][0]["sequence"] == "ACGTAC"
    with pytest.raises(ConfigurationError):
        sh.shasta_read_assembly_output(str(run), file_name="../x", verbose=False)
    with pytest.raises(ConfigurationError):
        sh.shasta_read_assembly_output(str(run), output="gfa_both_strands", verbose=False)


@pytest.fixture(scope="module")
def reads_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _reads(tmp_path_factory.mktemp("shasta_reads") / "reads.fasta")


@needs_shasta
def test_real_assembly_matches_direct_call(tmp_path: Path, reads_file: Path) -> None:
    result = sh.shasta_assemble([str(reads_file)], str(tmp_path / "run"), min_read_length=5000, threads=4, shasta_path=SHASTA, verbose=False)
    assert result["status"] == "completed" and result["version"].startswith("0.14")
    direct = tmp_path / "direct"
    subprocess.run([SHASTA, "--input", str(reads_file), "--config", "Nanopore-May2022", "--assemblyDirectory", str(direct),
                    "--threads", "4", "--Reads.minReadLength", "5000"], check=True, capture_output=True)
    expected = json.loads((direct / "AssemblySummary.json").read_text())["Assembled segments"]
    # Shasta's multithreaded runs can differ by a base between two runs, so compare with a 1% tolerance
    assert abs(result["summary"]["assembled_bases"] - expected["Total assembled segment length"]) <= 0.01 * expected["Total assembled segment length"]
    assert expected["Total assembled segment length"] > 0
    fasta = sh.shasta_read_assembly_output(result["assembly_directory"], verbose=False)
    assert fasta["total_length"] == result["summary"]["assembled_bases"] and fasta["n50"] == result["summary"]["segment_n50"]
    gfa = sh.shasta_read_assembly_output(result["assembly_directory"], output="gfa_both_strands", verbose=False)
    assert gfa["record_count"] == 2 * fasta["record_count"]
    reads = sh.shasta_read_run_table(result["assembly_directory"], "read_summary", verbose=False)
    assert reads["row_count"] == result["summary"]["reads_used"] == 450
    components = sh.shasta_read_run_table(result["assembly_directory"], "read_graph_components", verbose=False)
    assert components["row_count"] >= 1


@needs_shasta
@pytest.mark.parametrize(
    "parameters",
    [
        {"align_method": 1, "read_graph_creation_method": 2},
        {"align_method": 4},
        {"read_graph_creation_method": 3},
        {"read_graph_creation_method": 4, "strand_separation_method": 0},
        {"detangle_method": 1, "write_reads_by_segment": True},
        {"strand_separation_method": 2, "assembly_mode": 2},
        {"palindromic_flagging": False, "flag_inconsistent_alignments": True, "cleanup_duplicate_markers": True},
    ],
)
def test_real_stage_selections(tmp_path: Path, reads_file: Path, parameters: dict) -> None:
    result = sh.shasta_assemble([str(reads_file)], str(tmp_path / "run"), min_read_length=5000, threads=4, shasta_path=SHASTA, verbose=False, **parameters)
    assert result["status"] == "completed", result["log_tail"]
    if parameters.get("write_reads_by_segment"):
        assert sh.shasta_read_run_table(result["assembly_directory"], "reads_by_segment", verbose=False)["row_count"] > 0


@needs_shasta
def test_real_mode3_and_bad_option(tmp_path: Path, reads_file: Path) -> None:
    mode3 = _reads(tmp_path / "reads_hq.fasta", count=800, error=0.002, seed=1)
    result = sh.shasta_assemble([str(mode3)], str(tmp_path / "m3"), config="Nanopore-ncm23-May2024", assembly_mode=3,
                                min_read_length=5000, threads=4, shasta_path=SHASTA, verbose=False)
    assert result["return_code"] == 0
    assert (tmp_path / "m3" / "AnchorCoverageHistogram.csv").is_file()
    assert sh.shasta_read_run_table(result["assembly_directory"], "anchor_coverage_histogram", verbose=False)["row_count"] > 0
    failed = sh.shasta_assemble([str(reads_file)], str(tmp_path / "bad"), read_graph_creation_method=4, min_read_length=5000,
                                shasta_path=SHASTA, verbose=False)
    assert failed["status"] == "failed" and "strandSeparationMethod 0" in failed["log_tail"]
