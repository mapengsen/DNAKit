"""Tests of the seqtk facades: argument validation, argument vectors (stand-in program) and real seqtk runs."""

from __future__ import annotations

import json
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import seqtk_p12_cli_functions as seqtk

FASTA = ">a first\nACGTNNNNNNACGTAAGG\n>b\nacgtRYacgtSSKKggg\n>c\nCCCTAACCCTAACCCTAACCCTAA\n"
FASTQ = (
    "@r1/1 comment one\nACGTACGTAC\n+\nIIIIIIII##\n"
    "@r1/2 comment two\nTTGGCCAATT\n+\nIIIIIIIIII\n"
    "@r2/1\nGGGGAAAACC\n+\n5555555555\n"
    "@r3/1\nACACACACAC\n+\nIIIIIIIIII\n"
    "@r3/2\nTTTTAAAATT\n+\nIIIIIIIIII\n"
)
REAL = shutil.which("seqtk") is not None
needs_seqtk = pytest.mark.skipif(not REAL, reason="seqtk is not installed")


@pytest.fixture()
def files(tmp_path: Path) -> dict[str, Path]:
    fasta = tmp_path / "in.fa"
    fastq = tmp_path / "in.fq"
    fasta.write_text(FASTA, encoding="utf-8")
    fastq.write_text(FASTQ, encoding="utf-8")
    return {"fa": fasta, "fq": fastq, "dir": tmp_path}


def _direct(*arguments: str) -> str:
    return subprocess.run(["seqtk", *arguments], capture_output=True, text=True, check=True).stdout


def _stand_in(tmp_path: Path, stdout: str = "") -> tuple[Path, Path]:
    """A fake ``seqtk`` that records its argument vector and prints ``stdout`` (and a version line)."""
    record = tmp_path / "record.jsonl"
    script = tmp_path / "bin" / "seqtk"
    script.parent.mkdir()
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        f"open({str(record)!r}, 'a').write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if len(sys.argv) == 1:\n"
        "    sys.stderr.write('Usage: seqtk <command>\\nVersion: 9.9-test\\n'); sys.exit(1)\n"
        f"sys.stdout.write({stdout!r})\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script, record


def _first_call(record: Path) -> list[str]:
    return json.loads(record.read_text(encoding="utf-8").splitlines()[0])


# ---------------------------------------------------------------------------
# Validation and missing program
# ---------------------------------------------------------------------------


def test_public_names_have_verbose_and_are_exported() -> None:
    import inspect

    assert len(seqtk.__all__) == 16
    for name in seqtk.__all__:
        parameter = inspect.signature(getattr(seqtk, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True


def test_missing_program_raises_backend_unavailable(files: dict[str, Path]) -> None:
    with pytest.raises(BackendUnavailableError):
        seqtk.seqtk_size(
            str(files["fa"]),
            str(files["dir"] / "o.txt"),
            seqtk_path="/nonexistent/seqtk",
            verbose=False,
        )


@pytest.mark.parametrize(
    "keywords",
    [
        {"mate": 3},
        {"mate": True},
        {"mask_quality_below": 0},
        {"mask_character": "ab"},
        {"mask_character": " "},
        {"line_length": 0},
        {"drop_comments": "yes"},
    ],
)
def test_seq_rejects_bad_arguments(files: dict[str, Path], keywords: dict[str, object]) -> None:
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_seq(
            str(files["fq"]),
            str(files["dir"] / "o.fq"),
            seqtk_path="/nonexistent/seqtk",
            verbose=False,
            **keywords,  # type: ignore[arg-type]
        )


def test_inputs_and_outputs_are_checked(files: dict[str, Path]) -> None:
    with pytest.raises(ConfigurationError, match="existing file"):
        seqtk.seqtk_hpc(str(files["dir"] / "missing.fa"), str(files["dir"] / "o.fa"), verbose=False)
    existing = files["dir"] / "exists.txt"
    existing.write_text("x")
    with pytest.raises(ConfigurationError, match="must not exist"):
        seqtk.seqtk_hpc(str(files["fa"]), str(existing), verbose=False)
    with pytest.raises(ConfigurationError, match="folder"):
        seqtk.seqtk_hpc(str(files["fa"]), str(files["dir"] / "nodir" / "o.fa"), verbose=False)


@pytest.mark.parametrize("kmer", ["", "ACGN", "ACGTACGTACGTA", "-ACG", 5])
def test_kfreq_rejects_bad_kmers(files: dict[str, Path], kmer: object) -> None:
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_kfreq(kmer, str(files["fa"]), str(files["dir"] / "o.txt"), verbose=False)  # type: ignore[arg-type]


def test_other_parameters_are_validated(files: dict[str, Path]) -> None:
    out = str(files["dir"] / "o.txt")
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_gc(str(files["fa"]), out, min_fraction=0.0, verbose=False)
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_gc(str(files["fa"]), out, min_fraction=1.5, verbose=False)
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_telo(str(files["fa"]), out, motif="CCCTAN", verbose=False)
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_hety(str(files["fa"]), out, window_size=10, start_positions=20, verbose=False)
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_fqchk(str(files["fq"]), out, quality_threshold=-1, verbose=False)
    with pytest.raises(ConfigurationError):
        seqtk.seqtk_cutn(str(files["fa"]), out, gaps_only="no", verbose=False)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Argument vectors with a stand-in program
# ---------------------------------------------------------------------------


def test_seq_builds_the_expected_argument_vector(files: dict[str, Path]) -> None:
    script, record = _stand_in(files["dir"], stdout="@x\nA\n+\nI\n")
    out = files["dir"] / "o.fq"
    result = seqtk.seqtk_seq(
        str(files["fq"]),
        str(out),
        drop_comments=True,
        mask_quality_below=20,
        mask_character="N",
        quality_shift=64,
        uppercase=True,
        mate=2,
        line_length=60,
        seqtk_path=str(script),
        verbose=False,
    )
    assert _first_call(record) == [
        "seq",
        "-C",
        "-q",
        "20",
        "-n",
        "N",
        "-Q",
        "64",
        "-U",
        "-2",
        "-l",
        "60",
        str(files["fq"]),
    ]
    assert result["status"] == "completed" and result["record_count"] == 1
    assert result["outputs"] == [str(out)] and result["operation"] == "seqtk_seq"
    assert result["version"] == "9.9-test"
    assert json.dumps(result)


def test_table_parsers_with_a_stand_in(files: dict[str, Path]) -> None:
    script, record = _stand_in(files["dir"], stdout="s1\t0\t10\t4\ns2\t5\t9\t3\n")
    result = seqtk.seqtk_gc(
        str(files["fa"]),
        str(files["dir"] / "gc.txt"),
        at_rich=True,
        min_fraction=0.7,
        min_length=5,
        x_dropoff=8,
        seqtk_path=str(script),
        verbose=False,
    )
    assert _first_call(record) == ["gc", "-w", "-f", "0.7", "-l", "5", "-x", "8", str(files["fa"])]
    assert result["records"] == [
        {"seqid": "s1", "start": 0, "end": 10, "hit_bases": 4},
        {"seqid": "s2", "start": 5, "end": 9, "hit_bases": 3},
    ]
    assert result["record_count"] == 2


def test_max_records_limits_rows_but_not_the_count(files: dict[str, Path]) -> None:
    script, _ = _stand_in(files["dir"], stdout="a\t1\t2\nb\t3\t4\nc\t5\t6\n")
    result = seqtk.seqtk_gap(
        str(files["fa"]),
        str(files["dir"] / "g.txt"),
        max_records=1,
        seqtk_path=str(script),
        verbose=False,
    )
    assert len(result["records"]) == 1 and result["record_count"] == 3


def test_failed_program_raises(files: dict[str, Path]) -> None:
    script = files["dir"] / "fail"
    script.write_text(f"#!{sys.executable}\nimport sys\nsys.stderr.write('boom')\nsys.exit(3)\n")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    from dnakit.exceptions import BackendExecutionError

    with pytest.raises(BackendExecutionError, match="exited with code 3"):
        seqtk.seqtk_hpc(
            str(files["fa"]), str(files["dir"] / "o.fa"), seqtk_path=str(script), verbose=False
        )


# ---------------------------------------------------------------------------
# Real seqtk runs: results must equal the program's own output
# ---------------------------------------------------------------------------


@needs_seqtk
def test_seq_matches_seqtk(files: dict[str, Path]) -> None:
    out = files["dir"] / "seq.fq"
    result = seqtk.seqtk_seq(
        str(files["fq"]),
        str(out),
        drop_comments=True,
        mask_quality_below=30,
        uppercase=False,
        verbose=False,
    )
    assert out.read_text() == _direct("seq", "-C", "-q", "30", str(files["fq"]))
    assert "comment" not in out.read_text() and "ACGTACGTac" in out.read_text()
    assert result["record_count"] == 5 and result["version"]


@needs_seqtk
@pytest.mark.parametrize("mate,expected", [(1, ["r1/1", "r2/1", "r3/2"]), (2, ["r1/2", "r3/1"])])
def test_seq_selects_mates(files: dict[str, Path], mate: int, expected: list[str]) -> None:
    out = files["dir"] / f"mate{mate}.fq"
    seqtk.seqtk_seq(str(files["fq"]), str(out), mate=mate, verbose=False)
    names = [line.split()[0][1:] for line in out.read_text().splitlines()[::4]]
    assert names == expected


@needs_seqtk
def test_seq_uppercase(files: dict[str, Path]) -> None:
    out = files["dir"] / "up.fa"
    seqtk.seqtk_seq(str(files["fa"]), str(out), uppercase=True, verbose=False)
    sequences = [line for line in out.read_text().splitlines() if not line.startswith(">")]
    assert sequences and all(line == line.upper() for line in sequences)
    assert out.read_text() == _direct("seq", "-U", str(files["fa"]))


@needs_seqtk
def test_size_and_fqchk(files: dict[str, Path]) -> None:
    size = seqtk.seqtk_size(str(files["fq"]), str(files["dir"] / "size.txt"), verbose=False)
    assert (size["sequence_count"], size["base_count"]) == (5, 50)
    empty = files["dir"] / "empty.fa"
    empty.write_text("")
    zero = seqtk.seqtk_size(str(empty), str(files["dir"] / "zero.txt"), verbose=False)
    assert (zero["sequence_count"], zero["base_count"]) == (0, 0)
    chk = seqtk.seqtk_fqchk(str(files["fq"]), str(files["dir"] / "chk.txt"), verbose=False)
    assert chk["summary"]["min_len"] == 10 and chk["rows"][0]["POS"] == "ALL"
    assert chk["rows"][0]["#bases"] == 50 and chk["row_count"] == 11
    assert chk["columns"][-2:] == ["%low", "%high"]


@needs_seqtk
def test_dropse_hpc_and_listhet(files: dict[str, Path]) -> None:
    out = files["dir"] / "pairs.fq"
    result = seqtk.seqtk_dropse(str(files["fq"]), str(out), verbose=False)
    assert result["record_count"] == 4  # r2/1 has no mate
    assert "@r2/1" not in out.read_text()
    hpc = seqtk.seqtk_hpc(str(files["fa"]), str(files["dir"] / "hpc.fa"), verbose=False)
    assert hpc["record_count"] == 3
    assert "ACGTNACGTAG\n" in (files["dir"] / "hpc.fa").read_text()
    het = seqtk.seqtk_listhet(str(files["fa"]), str(files["dir"] / "het.txt"), verbose=False)
    assert het["records"][0] == {"seqid": "b", "position": 5, "base": "R"}
    assert [row["base"] for row in het["records"]] == ["R", "Y", "S", "S", "K", "K"]


@needs_seqtk
def test_gap_cutn_and_gc(files: dict[str, Path]) -> None:
    gap = seqtk.seqtk_gap(
        str(files["fa"]), str(files["dir"] / "gap.txt"), min_length=4, verbose=False
    )
    assert gap["records"] == [
        {"seqid": "a", "start": 4, "end": 10},
        {"seqid": "b", "start": 10, "end": 14},
    ]
    cut = seqtk.seqtk_cutn(
        str(files["fa"]), str(files["dir"] / "cut.fa"), min_n_tract=3, verbose=False
    )
    assert cut["record_count"] == 4  # a is cut into two pieces, b and c are kept
    gaps = seqtk.seqtk_cutn(
        str(files["fa"]),
        str(files["dir"] / "cutg.txt"),
        min_n_tract=3,
        gaps_only=True,
        verbose=False,
    )
    assert gaps["records"] == [{"seqid": "a", "start": 4, "end": 10}]
    gc = seqtk.seqtk_gc(str(files["fa"]), str(files["dir"] / "gc.txt"), min_length=4, verbose=False)
    assert gc["output_file"] and gc["record_count"] == len(
        (files["dir"] / "gc.txt").read_text().splitlines()
    )


@needs_seqtk
def test_mutfa_mergefa_famask_randbase(files: dict[str, Path]) -> None:
    snp = files["dir"] / "snp.txt"
    snp.write_text("a\t1\tA\tT\nc\t2\tC\tG\n")
    out = files["dir"] / "mut.fa"
    seqtk.seqtk_mutfa(str(files["fa"]), str(snp), str(out), verbose=False)
    assert out.read_text() == _direct("mutfa", str(files["fa"]), str(snp))
    assert out.read_text().splitlines()[1].startswith("TCGT") and out.read_text().splitlines()[
        5
    ].startswith("CGCT")
    other = files["dir"] / "other.fa"
    other.write_text(
        ">a\nACGTNNNNNNACGTAAGG\n>b\nacgtRYacgtSSKKggg\n>c\nCCCTTACCTTAACCTTAACCTTAA\n"
    )
    merged = files["dir"] / "merged.fa"
    seqtk.seqtk_mergefa(str(files["fa"]), str(other), str(merged), verbose=False)
    assert merged.read_text() == _direct("mergefa", "-q", "0", str(files["fa"]), str(other))
    mask = files["dir"] / "mask.fa"
    mask.write_text(">a\nXXXXXXXXXXXXXXXXXX\n>b\nNNNNNNNNNNNNNNNNN\n>c\nXXXXXXXXXXXXXXXXXXXXXXXX\n")
    masked = files["dir"] / "masked.fa"
    seqtk.seqtk_famask(str(files["fa"]), str(mask), str(masked), verbose=False)
    assert masked.read_text() == _direct("famask", str(files["fa"]), str(mask))
    rand = files["dir"] / "rand.fa"
    seqtk.seqtk_randbase(str(files["fa"]), str(rand), verbose=False)
    body = "".join(line for line in rand.read_text().splitlines() if not line.startswith(">"))
    assert not set(body) & set("RYSKrysk") and len(body) == 18 + 17 + 24


@needs_seqtk
def test_hety_kfreq_telo(files: dict[str, Path]) -> None:
    het = seqtk.seqtk_hety(
        str(files["fa"]),
        str(files["dir"] / "hety.txt"),
        window_size=10,
        start_positions=1,
        verbose=False,
    )
    assert het["record_count"] == len((files["dir"] / "hety.txt").read_text().splitlines())
    assert all(
        set(row) == {"seqid", "start", "end", "hets_per_window", "called_sites", "het_sites"}
        for row in het["records"]
    )
    kfreq = seqtk.seqtk_kfreq(
        "CCCTAA", str(files["fa"]), str(files["dir"] / "k.txt"), verbose=False
    )
    by_id = {row["seqid"]: row for row in kfreq["records"]}
    assert (
        by_id["c"]["exact_count"] == 4
        and by_id["c"]["strand"] == "+"
        and by_id["c"]["length"] == 24
    )
    telo = seqtk.seqtk_telo(
        str(files["fa"]), str(files["dir"] / "telo.txt"), min_score=10, verbose=False
    )
    assert telo["records"][0]["seqid"] == "c" and telo["records"][0]["start"] == 0
    assert telo["total_bases"] == 18 + 17 + 24 and telo["telomere_bases"] >= 18
