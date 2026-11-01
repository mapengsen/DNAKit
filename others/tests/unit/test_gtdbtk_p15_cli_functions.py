"""GTDB-Tk facade tests.

GTDB-Tk and its reference package (tens of gigabytes) are not available in the test environment, so a small stand-in
program that follows GTDB-Tk's command-line interface and output file names is used to test argument assembly, output
discovery and table parsing. ``trim_msa`` with a mask file is also checked against the real GTDB-Tk package when it is
installed.
"""

from __future__ import annotations

import json
import os
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.microbiome import gtdbtk_p15_cli_functions as gt

STUB = r'''#!{python}
import json, os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("gtdbtk: version 2.8.0")
    sys.exit(0)
log = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calls.jsonl")
with open(log, "a") as handle:
    handle.write(json.dumps({{"args": args, "data": os.environ.get("GTDBTK_DATA_PATH")}}) + "\n")
command = args[0]
def option(name):
    return args[args.index(name) + 1] if name in args else None
out = option("--out_dir")
prefix = option("--prefix") or "gtdbtk"
def put(relative, text):
    path = os.path.join(out, relative)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(text)
if command == "identify":
    put(f"identify/{{prefix}}.bac120.markers_summary.tsv", "name\tnumber_unique_genes\tnumber_multiple_genes\ng1\t100\t2\n")
    put(f"identify/{{prefix}}.translation_table_summary.tsv", "accession\ttranslation_table\ng1\t11\n")
elif command == "align":
    put(f"align/{{prefix}}.bac120.user_msa.fasta.gz", "x")
    put(f"align/{{prefix}}.bac120.filtered.tsv", "Genome ID\tFilter Reason\ng2\tToo few AA\n")
elif command == "classify":
    put(f"classify/{{prefix}}.bac120.summary.tsv", "user_genome\tclassification\tfastani_ani\ng1\td__Bacteria;p__X\t97.5\n")
    put(f"classify/{{prefix}}.bac120.classify.tree.1.tree", "(a,b);")
elif command == "de_novo_wf":
    put(f"{{prefix}}.bac120.decorated.tree", "(a,b);")
elif command == "ani_rep":
    put(f"{{prefix}}.ani_summary.tsv", "user_genome\tref_genome\tani\ng1\tGCF_1\t99.1\n")
elif command == "trim_msa":
    open(option("--output"), "w").write(">s1\nAC\n>s2\nGT\n")
elif command == "export_msa":
    open(option("--output"), "w").write(">r1\nACGT\n")
'''


@pytest.fixture()
def stub(tmp_path: Path) -> Path:
    path = tmp_path / "bin" / "gtdbtk"
    path.parent.mkdir()
    path.write_text(STUB.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def calls(stub: Path) -> list[dict[str, object]]:
    return [json.loads(line) for line in (stub.parent / "calls.jsonl").read_text().splitlines()]


def genomes(tmp_path: Path) -> Path:
    folder = tmp_path / "genomes"
    folder.mkdir(exist_ok=True)
    (folder / "g1.fna").write_text(">c\nACGT\n")
    return folder


def test_missing_program_has_install_hint(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        gt.gtdbtk_export_reference_msa(str(tmp_path / "o.fa"), database=str(tmp_path), gtdbtk_path=str(tmp_path / "missing"), verbose=False)


def test_database_is_required(tmp_path: Path, stub: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GTDBTK_DATA_PATH", raising=False)
    with pytest.raises(ConfigurationError):
        gt.gtdbtk_export_reference_msa(str(tmp_path / "o.fa"), gtdbtk_path=str(stub), verbose=False)


def test_genome_input_is_exclusive(tmp_path: Path, stub: Path) -> None:
    with pytest.raises(ConfigurationError):
        gt.gtdbtk_identify_markers(str(tmp_path / "out"), database=str(tmp_path), gtdbtk_path=str(stub), verbose=False)
    with pytest.raises(ConfigurationError):
        gt.gtdbtk_align_markers(
            str(tmp_path / "o2"), str(tmp_path), custom_msa_filters=True, skip_trimming=True, database=str(tmp_path),
            gtdbtk_path=str(stub), verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gt.gtdbtk_trim_msa(str(tmp_path / "m.fa"), str(tmp_path / "t.fa"), gtdbtk_path=str(stub), verbose=False)


def test_identify_and_align_commands(tmp_path: Path, stub: Path) -> None:
    database = tmp_path / "db"
    database.mkdir()
    identify = gt.gtdbtk_identify_markers(
        str(tmp_path / "id"), genome_dir=str(genomes(tmp_path)), database=str(database), prefix="run", cpus=2, genes=True,
        gtdbtk_path=str(stub), verbose=False,
    )
    assert identify["status"] == "completed" and identify["version"] == "2.8.0"
    assert identify["markers_summary"]["identify/run.bac120.markers_summary.tsv"][0]["number_unique_genes"] == 100
    first = calls(stub)[0]
    assert first["data"] == str(database)
    assert first["args"][0] == "identify" and "--genes" in first["args"] and first["args"][first["args"].index("--cpus") + 1] == "2"
    align = gt.gtdbtk_align_markers(
        str(tmp_path / "al"), str(tmp_path / "id"), database=str(database), prefix="run", skip_gtdb_refs=True, taxa_filter="p__X",
        gtdbtk_path=str(stub), verbose=False,
    )
    assert align["files"]["*.user_msa.fasta*"][0].endswith("run.bac120.user_msa.fasta.gz")
    assert align["filtered_genomes"]["align/run.bac120.filtered.tsv"][0]["Filter Reason"] == "Too few AA"
    args = calls(stub)[1]["args"]
    assert args[0] == "align" and "--skip_gtdb_refs" in args and args[args.index("--taxa_filter") + 1] == "p__X"


def test_classify_de_novo_and_ani(tmp_path: Path, stub: Path) -> None:
    database = tmp_path / "db"
    database.mkdir()
    classify = gt.gtdbtk_classify_placement(
        str(tmp_path / "cl"), str(tmp_path), genome_dir=str(genomes(tmp_path)), database=str(database), full_tree=True, min_af=0.5,
        gtdbtk_path=str(stub), verbose=False,
    )
    assert classify["summary"]["classify/gtdbtk.bac120.summary.tsv"][0]["classification"] == "d__Bacteria;p__X"
    args = calls(stub)[-1]["args"]
    assert args[0] == "classify" and "--full_tree" in args and args[args.index("--min_af") + 1] == "0.5"
    denovo = gt.gtdbtk_de_novo_workflow(
        str(tmp_path / "dn"), "p__Patescibacteria", batchfile=str(_batch(tmp_path)), database=str(database), prot_model="LG", gamma=True,
        gtdbtk_path=str(stub), verbose=False,
    )
    assert denovo["trees"]["*.decorated.tree"][0].endswith("gtdbtk.bac120.decorated.tree")
    args = calls(stub)[-1]["args"]
    assert args[0] == "de_novo_wf" and "--bacteria" in args and args[args.index("--prot_model") + 1] == "LG" and "--gamma" in args
    ani = gt.gtdbtk_ani_to_representatives(
        str(tmp_path / "an"), genome_dir=str(genomes(tmp_path)), database=str(database), gtdbtk_path=str(stub), verbose=False
    )
    assert ani["ani_summary"][0]["ani"] == 99.1


def _batch(tmp_path: Path) -> Path:
    path = tmp_path / "batch.tsv"
    path.write_text(f"{genomes(tmp_path)}/g1.fna\tg1\n")
    return path


def test_trim_and_export_msa(tmp_path: Path, stub: Path) -> None:
    mask = tmp_path / "mask.txt"
    mask.write_text("10")
    msa = tmp_path / "msa.fa"
    msa.write_text(">s1\nAC\n>s2\nGT\n")
    trimmed = gt.gtdbtk_trim_msa(str(msa), str(tmp_path / "trim.fa"), mask_file=str(mask), database=str(tmp_path), gtdbtk_path=str(stub), verbose=False)
    assert trimmed["sequence_count"] == 2
    exported = gt.gtdbtk_export_reference_msa(str(tmp_path / "ref.fa"), domain="arc", database=str(tmp_path), gtdbtk_path=str(stub), verbose=False)
    assert exported["sequence_count"] == 1 and calls(stub)[-1]["args"][:3] == ["export_msa", "--domain", "arc"]


REAL = shutil.which("gtdbtk") is not None and os.environ.get("GTDBTK_DATA_PATH")


@pytest.mark.skipif(not REAL, reason="GTDB-Tk with a reference package is not installed")
def test_real_trim_msa_with_mask_file(tmp_path: Path) -> None:
    mask = tmp_path / "mask.txt"
    mask.write_text("1010")
    msa = tmp_path / "msa.fa"
    msa.write_text(">s1\nACGT\n>s2\nTTGG\n")
    result = gt.gtdbtk_trim_msa(str(msa), str(tmp_path / "t.fa"), mask_file=str(mask), verbose=False)
    assert (tmp_path / "t.fa").read_text() == ">s1\nAG\n>s2\nTG\n"
    assert result["status"] == "completed"
