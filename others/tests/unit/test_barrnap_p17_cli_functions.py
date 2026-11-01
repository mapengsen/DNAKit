"""Tests for dnakit.annotation.barrnap_p17_cli_functions (Barrnap ``--mrna`` feature prediction).

Command assembly, GFF3 parsing and failure handling are tested with a small stand-in ``barrnap`` script written into
tmp_path (it prints GFF3 on standard output and its arguments on standard error, as Barrnap separates them). The test
marked ``requires_barrnap`` runs the real program and skips when Barrnap or one of its tools is missing. No test uses the
network.
"""

from __future__ import annotations

import inspect
import json
import shutil
import stat
import sys
from pathlib import Path

import pytest

import dnakit.annotation.barrnap_p17_cli_functions as barrnap_module
from dnakit.annotation.barrnap_p17_cli_functions import predict_mrna_features_barrnap
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"

requires_barrnap = pytest.mark.skipif(
    any(shutil.which(name) is None for name in ("barrnap", "pyrodigal", "diamond", "transterm", "ziggypep")),
    reason="barrnap 1.x with pyrodigal, diamond, transterm and ziggypep is not on PATH",
)

FAKE_BARRNAP = """#!{python}
import json, sys

args = sys.argv[1:]
if args == ["--version"]:
    print("barrnap 1.10.5")
    sys.exit(0)
print("[12:00:00] fake barrnap", file=sys.stderr)
print("ARGS " + json.dumps(args), file=sys.stderr)
if "--fail" in sys.argv[0]:
    sys.exit(2)
print("##gff-version 3")
print("ctg1\\tbarrnap:1.10.5\\tmRNA\\t90\\t1000\\t.\\t+\\t.\\tproduct=messenger RNA with RBS")
print("ctg1\\tbarrnap:1.10.5\\tRBS\\t90\\t95\\t.\\t+\\t.\\tproduct=ribosome binding site GGAG/GAGG")
print("ctg1\\tdiamond:2.1.9\\tCDS\\t100\\t1000\\t1e-50\\t+\\t0\\tName=dxs;product=1-deoxy-D-xylulose-5-phosphate synthase;Dbxref=sP:P0DP")
print("ctg1\\tziggypep:0.3\\tsig_pep\\t100\\t160\\t.\\t+\\t0\\tproduct=signal peptide MKKLLV")
print("ctg2\\ttransterm:2.0\\tterminator\\t2000\\t2030\\t99\\t-\\t.\\tproduct=Rho-independent terminator")
"""


def _fake_barrnap(tmp_path: Path, name: str = "barrnap") -> str:
    path = tmp_path / "fake_bin" / name
    path.parent.mkdir(exist_ok=True)
    path.write_text(FAKE_BARRNAP.replace("{python}", sys.executable), encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def _fasta(tmp_path: Path, name: str = "contigs.fa") -> str:
    path = tmp_path / name
    path.write_text(">ctg1\nACGTACGTACGT\n>ctg2\nTTGACCGGTTAA\n", encoding="utf-8")
    return str(path)


def _arguments(result: dict) -> list[str]:
    line = next(line for line in result["log_tail"].splitlines() if line.startswith("ARGS "))
    return json.loads(line[5:])


def test_public_function_has_keyword_only_verbose() -> None:
    assert barrnap_module.__all__ == ["predict_mrna_features_barrnap"]
    parameter = inspect.signature(predict_mrna_features_barrnap).parameters["verbose"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is True


def test_command_assembly_and_parsed_features(tmp_path: Path) -> None:
    output = tmp_path / "mrna.gff3"
    result = predict_mrna_features_barrnap(
        [_fasta(tmp_path)],
        str(output),
        kingdom="arc",
        evalue=1e-5,
        threads=2,
        add_ids=True,
        barrnap_path=_fake_barrnap(tmp_path),
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["version"] == "1.10.5"
    arguments = _arguments(result)
    assert "--mrna" in arguments and "--no-rrna" in arguments and "--addids" in arguments
    assert arguments[arguments.index("--kingdom") + 1] == "arc"
    assert arguments[arguments.index("--gcode") + 1] == "auto"
    assert arguments[arguments.index("--evalue") + 1] == "1e-05"
    assert arguments[arguments.index("--threads") + 1] == "2"
    assert arguments[-1] == str(tmp_path / "contigs.fa")
    assert result["gff_file"] == str(output)
    assert output.read_text().startswith("##gff-version 3")
    assert result["feature_counts"] == {"mRNA": 1, "RBS": 1, "CDS": 1, "sig_pep": 1, "terminator": 1}
    assert result["feature_count"] == 5 and result["features_truncated"] is False
    cds = next(feature for feature in result["features"] if feature["type"] == "CDS")
    assert (cds["seqid"], cds["start"], cds["end"], cds["strand"], cds["phase"]) == ("ctg1", 100, 1000, "+", "0")
    assert cds["attributes"]["Name"] == "dxs"
    assert not list(tmp_path.glob("*.stderr"))
    json.dumps(result)


def test_include_rrna_genetic_code_and_database_dir(tmp_path: Path) -> None:
    database = tmp_path / "db"
    database.mkdir()
    result = predict_mrna_features_barrnap(
        [_fasta(tmp_path)],
        str(tmp_path / "all.gff3"),
        genetic_code=4,
        include_rrna=True,
        database_dir=str(database),
        max_features=2,
        barrnap_path=_fake_barrnap(tmp_path),
        verbose=False,
    )
    arguments = _arguments(result)
    assert "--rrna" in arguments and "--no-rrna" not in arguments
    assert arguments[arguments.index("--gcode") + 1] == "4"
    assert arguments[arguments.index("--dbdir") + 1] == str(database)
    assert len(result["features"]) == 2 and result["feature_count"] == 5 and result["features_truncated"] is True


def test_program_failure_is_reported_as_failed(tmp_path: Path) -> None:
    result = predict_mrna_features_barrnap(
        [_fasta(tmp_path)],
        str(tmp_path / "x.gff3"),
        barrnap_path=_fake_barrnap(tmp_path, "barrnap--fail"),
        verbose=False,
    )
    assert result["status"] == "failed" and result["return_code"] != 0
    assert result["features"] == [] and result["feature_count"] == 0
    assert "fake barrnap" in result["log_tail"]
    assert not list(tmp_path.glob("*.stderr"))


@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"kingdom": "euk"}, "kingdom"),
        ({"genetic_code": 0}, "genetic_code"),
        ({"genetic_code": "11"}, "genetic_code"),
        ({"evalue": 0.0}, "evalue"),
        ({"threads": 0}, "threads"),
        ({"max_features": -1}, "max_features"),
        ({"database_dir": "/nonexistent/db"}, "database_dir"),
        ({"sequence_files": []}, "sequence_files"),
        ({"sequence_files": "contigs.fa"}, "sequence_files"),
        ({"sequence_files": ["/nonexistent/contigs.fa"]}, "sequence_files"),
    ],
)
def test_invalid_arguments_are_refused(tmp_path: Path, override: dict, message: str) -> None:
    arguments = {"sequence_files": [_fasta(tmp_path)], "output_gff": str(tmp_path / "out.gff3")} | override
    with pytest.raises(ConfigurationError, match=message):
        predict_mrna_features_barrnap(**arguments, barrnap_path=_fake_barrnap(tmp_path), verbose=False)
    assert not (tmp_path / "out.gff3").exists()


def test_existing_output_is_never_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "out.gff3"
    output.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="must not exist"):
        predict_mrna_features_barrnap(
            [_fasta(tmp_path)], str(output), barrnap_path=_fake_barrnap(tmp_path), verbose=False
        )
    assert output.read_text() == "keep"


def test_missing_barrnap_raises_install_hint(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python"))
    with pytest.raises(BackendUnavailableError) as raised:
        predict_mrna_features_barrnap([_fasta(tmp_path)], str(tmp_path / "out.gff3"), verbose=False)
    assert "barrnap" in (raised.value.hint or "").lower()
    assert not (tmp_path / "out.gff3").exists()


# ---------------------------------------------------------------------------
# Real Barrnap (skipped without Barrnap 1.x and its tools)
# ---------------------------------------------------------------------------


@requires_barrnap
def test_real_barrnap_predicts_genes_on_the_shared_genome(tmp_path: Path) -> None:
    genome = DATA / "genome.fa"
    if not genome.is_file():
        pytest.skip("shared consistency/data/genome.fa is missing")
    result = predict_mrna_features_barrnap([str(genome)], str(tmp_path / "mrna.gff3"), verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert {"mRNA", "CDS"} <= set(result["feature_counts"])
    assert result["feature_counts"]["mRNA"] == result["feature_counts"]["CDS"]
