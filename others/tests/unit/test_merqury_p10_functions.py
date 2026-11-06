"""Tests for dnakit.assembly.merqury_p10_functions (Merqury scripts as typed facades).

Most tests build a stand-in Merqury folder and stand-in programs in tmp_path: the scripts record their arguments and create the
files the real scripts create. Tests marked ``requires_merqury`` run the real Merqury and Meryl on a tiny genome and skip when
they are not installed (``conda install -c bioconda merqury``).
"""

from __future__ import annotations

import inspect
import random
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from dnakit.assembly import merqury_p10_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _executable(path: Path, body: str) -> None:
    path.write_text("#!/bin/bash\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


@pytest.fixture
def fake(tmp_path: Path, monkeypatch) -> Path:
    home = tmp_path / "merqury"
    for directory in ("util", "eval", "build", "trio", "plot"):
        (home / directory).mkdir(parents=True)
    (home / "util" / "util.sh").write_text("", encoding="utf-8")
    log = 'echo "$0 $@" >> "$(dirname "$0")/../calls.log"\n'
    scripts = {
        "best_k.sh": 'echo "genome: $1"; echo 20.5\n',
        "eval/per_seq_qv.sh": 'printf "c1\\t4\\t1986\\t38.7\\t0.0001\\n" > "$3.qv"\n',
        "eval/false_duplications.sh": 'printf "hist\\tcutoff\\t1\\t2\\t3\\t4\\t>4\\tdup(>1)\\tall\\tdup%%\\n"; printf "%s\\t5\\t90\\t10\\t0\\t0\\t0\\t10\\t100\\t10\\n" "$1"\n',
        "build/intersect.sh": 'mkdir -p "$1_and_$2"; ln -s "$1_and_$2" "$3.meryl"\n',
        "build/filt.sh": "b=${1%.meryl}; echo 3 > $b.filt; : > $b.hist; : > $b.hist.ploidy; mkdir $b.gt3.meryl\n",
        "plot/to_hist_for_plotting.sh": 'printf "kmer\\tkmer_multiplicity\\tCount\\n"; printf "%s\\t1\\t5\\n" "$2"\n',
        "plot/plot_spectra_cn.R": 'for i in "$@"; do last=$i; done; : > "$(printf "%s" "$out").ln.png"\n',
    }
    for relative, body in scripts.items():
        _executable(home / relative, log + body)
    bin_directory = tmp_path / "bin"
    bin_directory.mkdir()
    for name in ("meryl", "meryl-lookup", "java", "awk", "samtools", "wigToBigWig", "Rscript"):
        if name == "Rscript":
            # runs the R script as bash, which is enough for the stand-in plot script
            _executable(
                bin_directory / name,
                'script=$1; shift; while [ $# -gt 0 ]; do [ "$1" = -o ] && out=$2; shift; done; : > $out.ln.png\n',
            )
        elif name == "awk":
            _executable(bin_directory / name, "exit 0\n")
        else:
            _executable(bin_directory / name, 'echo "$0 1.4.1"\n')
    monkeypatch.setenv("PATH", f"{bin_directory}:/usr/bin:/bin")
    return home


def _meryl(tmp_path: Path, name: str) -> str:
    path = tmp_path / f"{name}.meryl"
    path.mkdir()
    return str(path)


def test_best_k_parses_the_script_output(fake):
    result = m.merqury_best_k(3_000_000_000, merqury_home=str(fake), verbose=False)
    assert result["status"] == "completed" and result["k"] == 20.5 and result["k_rounded_up"] == 21
    assert (
        result["operation"] == "merqury_best_k"
        and result["provenance"]["tested_version"] == "1.4.1"
    )


def test_per_sequence_qv_links_inputs_and_parses_the_table(fake, tmp_path):
    fasta = tmp_path / "asm.fa"
    fasta.write_text(">c1\nACGT\n", encoding="utf-8")
    reads = _meryl(tmp_path, "reads")
    result = m.merqury_per_sequence_qv(
        str(fasta), reads, str(tmp_path / "out"), merqury_home=str(fake), verbose=False
    )
    assert result["status"] == "completed" and result["row_count"] == 1
    assert result["columns"][3] == "qv" and result["rows"][0][0] == "c1"
    assert (tmp_path / "out" / "asm.fa").is_symlink()  # the input folder is never written to
    assert not (tmp_path / "asm.0.meryl").exists()


def test_false_duplications_parses_counts(fake, tmp_path):
    hist = tmp_path / "x.hist"
    hist.write_text("1\t1\t5\n", encoding="utf-8")
    result = m.merqury_false_duplications(str(hist), merqury_home=str(fake), verbose=False)
    assert (
        result["cutoff"] == 5 and result["copies"]["1"] == 90 and result["duplicated_percent"] == 10
    )


def test_combine_kmer_dbs_checks_database_counts(fake, tmp_path):
    first, second = _meryl(tmp_path, "a"), _meryl(tmp_path, "b")
    result = m.merqury_combine_kmer_dbs(
        "intersect",
        [first, second],
        str(tmp_path / "o"),
        "both",
        merqury_home=str(fake),
        verbose=False,
    )
    assert result["status"] == "completed" and result["combine_operation"] == "intersect"
    with pytest.raises(ConfigurationError, match="exactly two"):
        m.merqury_combine_kmer_dbs(
            "difference", [first], str(tmp_path / "o2"), "x", merqury_home=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError, match="needs the k-mer size"):
        m.merqury_combine_kmer_dbs(
            "union_sum",
            [first, second],
            str(tmp_path / "o3"),
            "x",
            merqury_home=str(fake),
            verbose=False,
        )
    assert not (tmp_path / "o2").exists() and not (tmp_path / "o3").exists()


def test_filter_solid_kmers_reports_the_threshold(fake, tmp_path):
    result = m.merqury_filter_solid_kmers(
        _meryl(tmp_path, "reads"), str(tmp_path / "f"), merqury_home=str(fake), verbose=False
    )
    assert result["status"] == "completed" and result["threshold"] == 3
    assert result["solid_database"].endswith("reads.gt3.meryl")


def test_histogram_table_streams_to_the_file(fake, tmp_path):
    result = m.merqury_kmer_histogram_table(
        {"reads": _meryl(tmp_path, "reads")},
        str(tmp_path / "h.tsv"),
        merqury_home=str(fake),
        verbose=False,
    )
    assert result["columns"] == ["kmer", "kmer_multiplicity", "Count"] and result["rows"] == [
        ["reads", "1", "5"]
    ]


def test_inputs_with_the_same_name_are_refused(fake, tmp_path):
    other = tmp_path / "other"
    other.mkdir()
    (other / "reads.meryl").mkdir()
    with pytest.raises(ConfigurationError, match="same file name"):
        m.merqury_combine_kmer_dbs(
            "intersect",
            [_meryl(tmp_path, "reads"), str(other / "reads.meryl")],
            str(tmp_path / "o"),
            "x",
            merqury_home=str(fake),
            verbose=False,
        )


def test_bad_arguments_are_rejected(fake, tmp_path):
    fasta = tmp_path / "asm.txt"
    fasta.write_text(">c1\nA\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="named"):
        m.merqury_per_sequence_qv(
            str(fasta),
            _meryl(tmp_path, "r"),
            str(tmp_path / "o"),
            merqury_home=str(fake),
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="must be one of"):
        m.merqury_kmer_multiplicity_track(
            str(fasta),
            _meryl(tmp_path, "r2"),
            str(tmp_path / "o"),
            source="both",
            merqury_home=str(fake),
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="at least one"):
        m.merqury_plot_block_n(str(tmp_path / "o"), "x", merqury_home=str(fake), verbose=False)
    assert not (tmp_path / "o").exists()


def test_missing_merqury_and_missing_programs_raise_with_hints(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.delenv("MERQURY", raising=False)
    with pytest.raises(BackendUnavailableError) as error:
        m.merqury_best_k(1000, verbose=False)
    assert "Merqury" in str(error.value)


def test_missing_program_leaves_no_output_directory(fake, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    fasta = tmp_path / "asm.fa"
    fasta.write_text(">c1\nA\n", encoding="utf-8")
    with pytest.raises(BackendUnavailableError, match=r"meryl|bash"):
        m.merqury_per_sequence_qv(
            str(fasta),
            _meryl(tmp_path, "r"),
            str(tmp_path / "o"),
            merqury_home=str(fake),
            verbose=False,
        )
    assert not (tmp_path / "o").exists()


def test_public_functions_carry_the_verbosity_parameter():
    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True


# ---------------------------------------------------------------------------
# Real Merqury
# ---------------------------------------------------------------------------


def _merqury_installed() -> bool:
    if shutil.which("meryl") is None or shutil.which("merqury.sh") is None:
        return False
    try:
        m._home(None)
    except BackendUnavailableError:
        return False
    return True


requires_merqury = pytest.mark.skipif(
    not _merqury_installed(), reason="Merqury and Meryl are not installed"
)


@requires_merqury
def test_real_per_sequence_qv_and_set_operations(tmp_path):
    rng = random.Random(3)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    (tmp_path / "asm.fa").write_text(
        f">c1\n{genome[:2000]}\n>c2\n{genome[2000:]}\n", encoding="utf-8"
    )
    reads = []
    for index in range(600):
        start = rng.randint(0, 2900)
        reads.append(f">r{index}\n{genome[start : start + 100]}\n")
    (tmp_path / "reads.fa").write_text("".join(reads), encoding="utf-8")
    for name in ("asm", "reads"):
        subprocess.run(
            ["meryl", "k=15", "count", "output", f"{name}.meryl", f"{name}.fa"],
            cwd=tmp_path,
            check=True,
            capture_output=True,
        )
    qv = m.merqury_per_sequence_qv(
        str(tmp_path / "asm.fa"), str(tmp_path / "reads.meryl"), str(tmp_path / "qv"), verbose=False
    )
    assert qv["status"] == "completed" and [row[0] for row in qv["rows"]] == ["c1", "c2"]
    assert all(float(row[3]) > 20 or row[3] == "+inf" for row in qv["rows"])
    union = m.merqury_combine_kmer_dbs(
        "union_sum",
        [str(tmp_path / "asm.meryl"), str(tmp_path / "reads.meryl")],
        str(tmp_path / "union"),
        "all",
        k=15,
        verbose=False,
    )
    assert union["status"] == "completed" and Path(union["database"]).is_dir()
    best = m.merqury_best_k(3_000_000_000, verbose=False)
    assert best["k_rounded_up"] == 21
