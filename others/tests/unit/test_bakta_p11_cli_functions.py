"""Tests for dnakit.annotation.bakta_p11_cli_functions: bakta_io / bakta_plot run for real when installed; the full
annotation (which needs a ~1.5 GB database) is tested with a stand-in executable that records its arguments."""

from __future__ import annotations

import json
import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.annotation import bakta_p11_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _json(path: Path) -> Path:
    rng = random.Random(2)
    seq = "".join(rng.choice("ACGT") for _ in range(3000))
    aa = "MKTAYIAKQRQISFVKSHFSRQ"
    cds = {"sequence": "c1", "strand": "+", "nt": seq[100:166], "aa": aa, "aa_digest": None, "start_type": "ATG",
           "rbs_motif": None, "frame": 2, "functions": [], "genes": []}
    data = {
        "genome": {"genus": "E", "species": "coli", "strain": "", "complete": False, "gram": "?", "translation_table": 11},
        "stats": {"no_sequences": 1, "size": 3000, "gc": 0.5, "n_ratio": 0.0, "n50": 3000, "n90": 3000, "coding_ratio": 0.1},
        "features": [
            {**cds, "type": "cds", "start": 101, "stop": 166, "gene": "abcD", "product": "test protein", "locus": "T_00005",
             "db_xrefs": ["SO:0001217", "COG:J"]},
            {**cds, "type": "cds", "start": 301, "stop": 366, "strand": "-", "hypothetical": True, "product": "hypothetical protein",
             "locus": "T_00010", "seq_stats": {"molecular_weight": 2500.0, "isoelectric_point": 9.1}, "db_xrefs": ["SO:0001217"]},
            {"type": "tRNA", "sequence": "c1", "start": 501, "stop": 576, "strand": "+", "gene": "tRNA-Ala", "product": "tRNA-Ala(tgc)",
             "locus": "T_00015", "nt": seq[500:576], "amino_acid": "Ala", "anti_codon": "tgc", "score": 80.5, "db_xrefs": ["SO:0001272"]},
            {"type": "gap", "sequence": "c1", "start": 901, "stop": 950, "strand": "?", "length": 50, "id": "gap1",
             "product": "assembly gap", "db_xrefs": ["SO:0000730"]},
        ],
        "sequences": [{"id": "c1", "description": "", "nt": seq, "length": 3000, "complete": False, "type": "contig",
                       "topology": "linear", "name": None, "orig_id": "c1", "orig_description": ""}],
        "version": {"bakta": "1.12.2", "db": {"type": "light", "version": "6.0", "doi": "10.5281/zenodo.0"}},
        "run": {"start": "2026-10-20 10:00:00", "end": "2026-10-20 10:01:00"},
    }
    path.write_text(json.dumps(data))
    return path


@pytest.mark.skipif(shutil.which("bakta_io") is None, reason="bakta_io (pip install bakta) is not installed")
def test_export_and_plot(tmp_path: Path) -> None:
    source = _json(tmp_path / "demo.json")
    exported = m.export_bakta_annotation(source, tmp_path / "e", verbose=False)
    for key in ("tsv", "gff3", "gbff", "embl", "fna", "ffn", "faa", "hypotheticals.tsv", "inference.tsv", "txt"):
        assert key in exported["outputs"]
    assert exported["summary"]["CDSs"] == 2 and exported["summary"]["tRNAs"] == 1
    plot = m.plot_bakta_genome(source, tmp_path / "p", plot_type="cog", size=4, dpi=150, verbose=False)
    assert Path(plot["outputs"]["png"]).is_file() and Path(plot["outputs"]["svg"]).is_file()


def test_missing_programs_and_validation(tmp_path: Path) -> None:
    source = _json(tmp_path / "demo.json")
    with pytest.raises(BackendUnavailableError):
        m.export_bakta_annotation(source, tmp_path / "e", bakta_io_path=str(tmp_path / "none"), verbose=False)
    with pytest.raises(ConfigurationError):
        m.plot_bakta_genome(source, tmp_path / "p", plot_type="bad", bakta_plot_path=sys.executable, verbose=False)


@pytest.mark.skipif(sys.platform == "win32", reason="stand-in executable is a POSIX script")
def test_annotate_features_arguments(tmp_path: Path) -> None:
    fake = tmp_path / "bakta"
    fake.write_text(
        f"#!{sys.executable}\nimport sys, json, pathlib\n"
        "a = sys.argv[1:]\n"
        "if a == ['--version']:\n    print('bakta 1.12.2'); sys.exit(0)\n"
        "out = pathlib.Path(a[a.index('--output') + 1]); out.mkdir(); p = a[a.index('--prefix') + 1]\n"
        "(out / (p + '.json')).write_text(json.dumps({'features': [{'type': 'tRNA'}, {'type': 'cds', 'hypothetical': True}]}))\n"
        "(out / (p + '.tsv')).write_text('x'); (out / (p + '.gff3')).write_text('x')\n"
        "(out / 'args.txt').write_text(' '.join(a))\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    genome = tmp_path / "g.fa"
    genome.write_text(">c\nACGT\n")
    db = tmp_path / "db"
    db.mkdir()
    result = m.annotate_bakta_features(
        genome, db, tmp_path / "o", features=["cds", "trna"], genus="Esch", translation_table=11, threads=2, bakta_path=str(fake), verbose=False
    )
    args = (Path(result["annotation_directory"]) / "args.txt").read_text().split()
    for flag in ("--skip-crispr", "--skip-sorf", "--skip-pseudo", "--skip-gap", "--skip-ori", "--skip-rrna", "--skip-plot"):
        assert flag in args
    assert "--skip-cds" not in args and "--skip-trna" not in args
    assert args[args.index("--genus") + 1] == "Esch" and result["feature_counts"] == {"tRNA": 1, "cds": 1, "hypothetical": 1}
    with pytest.raises(ConfigurationError):
        m.annotate_bakta_features(genome, db, tmp_path / "o2", features=["pseudogenes"], bakta_path=str(fake), verbose=False)
