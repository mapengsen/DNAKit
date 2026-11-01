"""Tests for dnakit.sequence.pyfaidx_gap2_functions (faidx fetch, write, transform, mask and region
parsing).
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("pyfaidx")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequence import pyfaidx_gap2_functions as m

DATA = (
    Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENES = DATA / "genes.fa"
_ENV_BIN = Path(sys.executable).parent
_HAS_FAIDX = (_ENV_BIN / "faidx").is_file()
pytestmark = pytest.mark.skipif(
    not (_HAS_FAIDX and GENES.exists()),
    reason="the faidx program (pyfaidx CLI) or the shared genes.fa is missing",
)

SMALL_FASTA = (
    ">seq1 first record\nACGTNNacgtACGTAAAACCCCGGGGTTTT\nGGCCAATT\n>seq2\nGGGGCCCCAAAATTTT\n"
)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@pytest.fixture(autouse=True)
def _faidx_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif line.strip():
            records[name].append(line.strip())
    return {key: "".join(parts) for key, parts in records.items()}


def _genes_copy(tmp_path: Path) -> Path:
    """The shared FASTA is copied first because faidx writes a .fai index next to its input."""
    copy = tmp_path / "genes.fa"
    shutil.copyfile(GENES, copy)
    return copy


def _small(tmp_path: Path) -> Path:
    path = tmp_path / "small.fa"
    path.write_text(SMALL_FASTA, encoding="utf-8")
    return path


def test_wrap_sequence_lines_end_with_newline() -> None:
    assert m.faidx_wrap_sequence("ACGTACGT", 3, verbose=False) == ["ACG\n", "TAC\n", "GT\n"]
    assert m.faidx_wrap_sequence("", 5, verbose=False) == []
    with pytest.raises(ConfigurationError):
        m.faidx_wrap_sequence("ACGT", 0, verbose=False)


def test_ucsc_region_is_zero_based_half_open() -> None:
    assert m.faidx_parse_ucsc_region("chr1:11-20", verbose=False) == {
        "name": "chr1",
        "start": 10,
        "end": 20,
    }
    assert m.faidx_parse_ucsc_region("chr1", verbose=False) == {
        "name": "chr1",
        "start": None,
        "end": None,
    }


def test_bed_entry_comment_and_malformed_lines() -> None:
    assert m.faidx_parse_bed_entry("chr1\t5\t9\tname\n", verbose=False) == {
        "name": "chr1",
        "start": 5,
        "end": 9,
    }
    assert m.faidx_parse_bed_entry("# header", verbose=False) == {
        "name": None,
        "start": None,
        "end": None,
    }
    with pytest.raises(ConfigurationError):
        m.faidx_parse_bed_entry("chr1\t5", verbose=False)
    with pytest.raises(ConfigurationError):
        m.faidx_parse_bed_entry("", verbose=False)


def test_split_regions_from_strings_and_bed(tmp_path: Path) -> None:
    records = m.faidx_split_regions(["cds000:1-30", "cds001"], verbose=False)
    assert records == [
        {"name": "cds000", "start": 0, "end": 30},
        {"name": "cds001", "start": None, "end": None},
    ]
    bed = tmp_path / "regions.bed"
    bed.write_text("# comment\ncds000\t10\t20\ncds002\t0\t5\n", encoding="utf-8")
    assert m.faidx_split_regions(bed=str(bed), verbose=False) == [
        {"name": "cds000", "start": 10, "end": 20},
        {"name": "cds002", "start": 0, "end": 5},
    ]
    with pytest.raises(ConfigurationError):
        m.faidx_split_regions(["cds000"], bed=str(bed), verbose=False)


def test_fetch_regions_match_the_fasta_contents(tmp_path: Path) -> None:
    fasta = _genes_copy(tmp_path)
    expected = _read_fasta(fasta)
    out = tmp_path / "fetched.fa"
    result = m.faidx_fetch_sequences(
        str(fasta),
        str(out),
        regions=["cds000:1-30", "cds001:101-200", "cds002"],
        verbose=False,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["version"] is not None
    fetched = out.read_text(encoding="utf-8")
    assert fetched.splitlines()[0] == ">cds000:1-30"
    parsed = _read_fasta(out)
    assert parsed["cds000:1-30"] == expected["cds000"][:30]
    assert parsed["cds001:101-200"] == expected["cds001"][100:200]
    assert parsed["cds002"] == expected["cds002"]


def test_fetch_complement_and_missing_region(tmp_path: Path) -> None:
    fasta = _small(tmp_path)
    out = tmp_path / "complement.fa"
    result = m.faidx_fetch_sequences(
        str(fasta), str(out), regions=["seq2:1-8"], complement=True, verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _read_fasta(out) == {"seq2:1-8": "GGGGCCCC".translate(COMPLEMENT)}
    empty = tmp_path / "missing.fa"
    missing = m.faidx_fetch_sequences(str(fasta), str(empty), regions=["nope:1-5"], verbose=False)
    assert missing["status"] == "completed"
    assert empty.read_text(encoding="utf-8") == ""


def test_fetch_refuses_conflicting_or_existing_outputs(tmp_path: Path) -> None:
    fasta = _small(tmp_path)
    with pytest.raises(ConfigurationError):
        m.faidx_fetch_sequences(
            str(fasta), str(tmp_path / "a.fa"), no_names=True, long_names=True, verbose=False
        )
    existing = tmp_path / "exists.fa"
    existing.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.faidx_fetch_sequences(str(fasta), str(existing), verbose=False)
    with pytest.raises(ConfigurationError):
        m.faidx_fetch_sequences(str(tmp_path / "absent.fa"), str(tmp_path / "b.fa"), verbose=False)


def test_write_with_name_regex_and_size_range(tmp_path: Path) -> None:
    fasta = _genes_copy(tmp_path)
    expected = _read_fasta(fasta)
    out = tmp_path / "selected.fa"
    result = m.faidx_write_sequences(
        str(fasta),
        output=str(out),
        regex=r"cds00[12]$",
        verbose=False,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert set(_read_fasta(out)) == {"cds001", "cds002"}
    inverted = tmp_path / "inverted.fa"
    m.faidx_write_sequences(
        str(fasta), output=str(inverted), regex="cds00", invert_match=True, verbose=False
    )
    assert set(_read_fasta(inverted)) == {name for name in expected if not name.startswith("cds00")}
    sized = tmp_path / "sized.fa"
    m.faidx_write_sequences(str(fasta), output=str(sized), size_range=(1000, 1600), verbose=False)
    kept = _read_fasta(sized)
    assert kept and all(1000 <= len(expected[name]) <= 1600 for name in kept)
    assert kept == {name: seq for name, seq in expected.items() if 1000 <= len(seq) <= 1600}


def test_write_whole_fasta_keeps_every_record(tmp_path: Path) -> None:
    """Regression: writing with faidx -o lost the tail of larger outputs, so the whole file is
    compared.
    """
    fasta = _genes_copy(tmp_path)
    expected = _read_fasta(fasta)
    out = tmp_path / "all.fa"
    result = m.faidx_write_sequences(str(fasta), output=str(out), verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert _read_fasta(out) == expected
    assert len(expected) == 40


def test_write_split_directory_one_file_per_sequence(tmp_path: Path) -> None:
    fasta = _genes_copy(tmp_path)
    split = tmp_path / "split"
    result = m.faidx_write_sequences(
        str(fasta), split_directory=str(split), size_range=(1000, 1600), verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    files = sorted(Path(p) for p in result["outputs"])
    assert files and all(p.parent == split for p in files)
    expected = _read_fasta(fasta)
    for path in files:
        name = path.name[: -len(".fa")]
        assert _read_fasta(path) == {name: expected[name]}
    with pytest.raises(ConfigurationError):
        m.faidx_write_sequences(
            str(fasta), output=str(tmp_path / "x.fa"), split_directory=str(split), verbose=False
        )


def test_transform_chromsizes_nucleotide_and_bed(tmp_path: Path) -> None:
    fasta = _genes_copy(tmp_path)
    expected = _read_fasta(fasta)
    sizes = tmp_path / "sizes.tsv"
    result = m.faidx_transform_sequences(
        str(fasta), str(sizes), "chromsizes", regions=["cds000", "cds001"], verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    rows = [line.split("\t") for line in sizes.read_text(encoding="utf-8").splitlines()]
    assert rows == [
        ["cds000", str(len(expected["cds000"]))],
        ["cds001", str(len(expected["cds001"]))],
    ]

    counts = tmp_path / "counts.tsv"
    nuc = m.faidx_transform_sequences(
        str(fasta), str(counts), "nucleotide", regions=["cds000"], verbose=False
    )
    assert nuc["status"] == "completed", nuc["log_tail"]
    lines = counts.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[:3] == ["name", "start", "end"]
    seq = expected["cds000"].upper()
    fields = lines[1].split("\t")
    assert fields[0] == "cds000"
    assert [int(v) for v in fields[3:8]] == [seq.count(base) for base in "ATCGN"]

    coords = tmp_path / "coords.bed"
    bed = m.faidx_transform_sequences(
        str(fasta), str(coords), "bed", regions=["cds002:11-40"], verbose=False
    )
    assert bed["status"] == "completed", bed["log_tail"]
    assert coords.read_text(encoding="utf-8") == "cds002\t10\t40\n"

    transposed = tmp_path / "transposed.tsv"
    tr = m.faidx_transform_sequences(
        str(fasta), str(transposed), "transposed", regions=["cds001:1-12"], verbose=False
    )
    assert tr["status"] == "completed", tr["log_tail"]
    assert transposed.read_text(encoding="utf-8") == f"cds001\t1\t12\t{expected['cds001'][:12]}\n"

    with pytest.raises(ConfigurationError):
        m.faidx_transform_sequences(str(fasta), str(tmp_path / "bad.txt"), "protein", verbose=False)


def test_mask_copies_and_leaves_input_unchanged(tmp_path: Path) -> None:
    fasta = _genes_copy(tmp_path)
    before = fasta.read_bytes()
    expected = _read_fasta(fasta)
    out = tmp_path / "masked.fa"
    result = m.faidx_mask_sequences(
        str(fasta), str(out), regions=["cds000:1-10", "cds001:5-8"], verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    assert fasta.read_bytes() == before
    masked = _read_fasta(out)
    assert masked["cds000"] == "N" * 10 + expected["cds000"][10:]
    assert masked["cds001"] == expected["cds001"][:4] + "N" * 4 + expected["cds001"][8:]
    assert masked["cds002"] == expected["cds002"]
    with pytest.raises(ConfigurationError):
        m.faidx_mask_sequences(str(fasta), str(tmp_path / "none.fa"), verbose=False)
    with pytest.raises(ConfigurationError):
        m.faidx_mask_sequences(
            str(fasta),
            str(tmp_path / "two.fa"),
            regions=["cds000"],
            default_seq="NN",
            verbose=False,
        )


def test_default_seq_rule_and_dash_mask(tmp_path: Path) -> None:
    assert m.faidx_check_default_seq(None, verbose=False) is None
    assert m.faidx_check_default_seq("N", verbose=False) == "N"
    for bad in ("", "NN"):
        with pytest.raises(ConfigurationError):
            m.faidx_check_default_seq(bad, verbose=False)
    fasta = _genes_copy(tmp_path)
    out = tmp_path / "dash.fa"
    result = m.faidx_mask_sequences(
        str(fasta), str(out), regions=["cds000:1-5"], default_seq="-", verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _read_fasta(out)["cds000"].startswith("-----")


def test_missing_program_raises_with_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fasta = _small(tmp_path)
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    with pytest.raises(BackendUnavailableError) as excinfo:
        m.faidx_fetch_sequences(
            str(fasta),
            str(tmp_path / "out.fa"),
            regions=["seq1:1-5"],
            faidx_path=str(empty / "faidx"),
            verbose=False,
        )
    assert "pyfaidx" in str(excinfo.value.hint or "")
