from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import seqmagick_p15_cli_functions as sm

needs_seqmagick = pytest.mark.skipif(shutil.which("seqmagick") is None, reason="seqmagick is not installed")

ALIGNMENT = ">s1\nAC--GTACGT\n>s2\nAC--GTTCGT\n>s3\nACGTGTACGA\n>s4\n----------\n"
READS = (
    "@r1\nCAGATGACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIII\n"
    "@r2\nTTTTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIII\n"
    "@r3\nCAGATGACNTACGTACGT\n+\nIIIIIIIIIIIIIIIIII\n"
)


def write(path: Path, text: str) -> str:
    path.write_text(text, encoding="utf-8")
    return str(path)


def records(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            result[name] = ""
        else:
            result[name] += line
    return result


def test_missing_program_has_install_hint(tmp_path: Path) -> None:
    source = write(tmp_path / "a.fa", ALIGNMENT)
    with pytest.raises(BackendUnavailableError):
        sm.seqmagick_slice_residues(source, str(tmp_path / "o.fa"), cut="1:3", seqmagick_path=str(tmp_path / "missing"), verbose=False)


def test_argument_validation(tmp_path: Path) -> None:
    source = write(tmp_path / "a.fa", ALIGNMENT)
    with pytest.raises(ConfigurationError):
        sm.seqmagick_slice_residues(source, str(tmp_path / "o.fa"), verbose=False)  # no range at all
    with pytest.raises(ConfigurationError):
        sm.seqmagick_slice_residues(source, str(tmp_path / "o.fa"), cut="1-3", verbose=False)
    with pytest.raises(ConfigurationError):
        sm.seqmagick_clean_alignment(source, str(tmp_path / "o.fa"), verbose=False)
    with pytest.raises(ConfigurationError):
        sm.seqmagick_clean_alignment(source, str(tmp_path / "o.fa"), squeeze=True, squeeze_threshold=0.5, verbose=False)
    with pytest.raises(ConfigurationError):
        sm.seqmagick_select_sequences(source, str(tmp_path / "o.fa"), verbose=False)
    with pytest.raises(ConfigurationError):
        sm.seqmagick_quality_filter(source, str(tmp_path / "o.fq"), primer="ACGT!", verbose=False)


@needs_seqmagick
def test_cut_drop_and_mask_ranges(tmp_path: Path) -> None:
    source = write(tmp_path / "a.fa", ALIGNMENT)
    cut = sm.seqmagick_slice_residues(source, str(tmp_path / "cut.fa"), cut="3:6", verbose=False)
    assert records(tmp_path / "cut.fa")["s1"] == "--GT"
    assert cut["output_record_count"] == 4 and cut["version"]
    sm.seqmagick_slice_residues(source, str(tmp_path / "drop.fa"), drop="1:2,9:10", verbose=False)
    assert records(tmp_path / "drop.fa")["s3"] == "GTGTAC"
    sm.seqmagick_slice_residues(source, str(tmp_path / "mask.fa"), mask="1:3", verbose=False)
    assert records(tmp_path / "mask.fa")["s3"] == "---TGTACGA"


@needs_seqmagick
def test_squeeze_dash_gap_and_prune_empty(tmp_path: Path) -> None:
    source = write(tmp_path / "a.fa", ALIGNMENT.replace("AC--GTACGT", "AC..GTACGT"))
    cleaned = sm.seqmagick_clean_alignment(source, str(tmp_path / "c.fa"), dash_gap=True, prune_empty=True, verbose=False)
    assert "s4" not in records(tmp_path / "c.fa")
    assert records(tmp_path / "c.fa")["s1"] == "AC--GTACGT"
    assert cleaned["output_record_count"] == 3
    plain = write(tmp_path / "plain.fa", ALIGNMENT)
    sm.seqmagick_clean_alignment(plain, str(tmp_path / "q.fa"), squeeze_threshold=0.75, verbose=False)
    assert len(records(tmp_path / "q.fa")["s3"]) == 8


@needs_seqmagick
def test_select_by_id_list_and_sequence_pattern(tmp_path: Path) -> None:
    source = write(tmp_path / "a.fa", ALIGNMENT)
    ids = write(tmp_path / "ids.txt", "s1\ns3\n")
    keep = sm.seqmagick_select_sequences(source, str(tmp_path / "keep.fa"), include_ids_file=ids, verbose=False)
    assert sorted(records(tmp_path / "keep.fa")) == ["s1", "s3"] and keep["output_record_count"] == 2
    sm.seqmagick_select_sequences(source, str(tmp_path / "drop.fa"), exclude_ids_file=ids, verbose=False)
    assert sorted(records(tmp_path / "drop.fa")) == ["s2", "s4"]
    sm.seqmagick_select_sequences(source, str(tmp_path / "pat.fa"), seq_pattern_include="GTAC", verbose=False)
    assert sorted(records(tmp_path / "pat.fa")) == ["s1", "s3"]


@needs_seqmagick
def test_quality_filter_ambiguous_and_primer(tmp_path: Path) -> None:
    reads = write(tmp_path / "r.fq", READS)
    dropped = sm.seqmagick_quality_filter(
        reads, str(tmp_path / "d.fq"), min_length=5, max_length=100, ambiguous_action="drop", verbose=False
    )
    assert dropped["output_record_count"] == 2  # r3 contains an N
    truncated = sm.seqmagick_quality_filter(
        reads, str(tmp_path / "t.fq"), min_length=5, max_length=100, ambiguous_action="truncate", verbose=False
    )
    lines = (tmp_path / "t.fq").read_text().splitlines()
    assert lines[1::4][2] == "CAGATGAC"  # r3 cut before the N
    assert truncated["report"][-1]["name"].startswith("Ambiguous Base")
    primed = sm.seqmagick_quality_filter(reads, str(tmp_path / "p.fq"), min_length=5, max_length=100, primer="CAGATG", verbose=False)
    assert primed["output_record_count"] == 2
    assert (tmp_path / "p.fq").read_text().splitlines()[1] == "ACGTACGTACGT"  # primer removed
    assert not (tmp_path / "p.fq.primer_only.csv").exists()


@needs_seqmagick
def test_info_table(tmp_path: Path) -> None:
    alignment = write(tmp_path / "a.fa", ALIGNMENT)
    reads = write(tmp_path / "r.fq", READS)
    result = sm.seqmagick_info([alignment, reads], str(tmp_path / "info.tsv"), verbose=False)
    by_name = {Path(row["name"]).name: row for row in result["files"]}
    assert by_name["a.fa"]["alignment"] == "TRUE" and by_name["a.fa"]["num_seqs"] == 4
    assert by_name["r.fq"]["alignment"] == "TRUE" and by_name["r.fq"]["min_len"] == 18
    csv = sm.seqmagick_info([alignment], str(tmp_path / "info.csv"), table_format="csv", verbose=False)
    assert csv["files"][0]["max_len"] == 10
