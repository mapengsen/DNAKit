"""Contract tests for the SeqKit command-line facades (no network; SeqKit tests skip without seqkit)."""

from __future__ import annotations

import inspect
import json
import re
import shutil
from pathlib import Path

import pytest

import dnakit.fasta_fastq.seqkit_cli_functions as sk
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

requires_seqkit = pytest.mark.skipif(
    shutil.which("seqkit") is None, reason="the seqkit executable is not on PATH"
)

FASTA = ">chr1\nACGTACGTNNACGTTTGCA\n>chr2\nGGGGCCCCAAAATTTT\n>chr3_a\nACGTAC\n"
FASTQ = "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nTTTTGGGGCC\n+\nIIIIIIIIII\n"
BED = "chr1\t4\t10\tf1\t0\t+\n"


def _write(path: Path, text: str) -> str:
    path.write_text(text, encoding="utf-8")
    return str(path)


def _fasta_ids(path: Path) -> list[str]:
    return [line[1:] for line in path.read_text().splitlines() if line.startswith(">")]


# ---------------------------------------------------------------------------
# Contract without SeqKit
# ---------------------------------------------------------------------------


def test_public_functions_are_keyword_only_verbose() -> None:
    assert len(sk.__all__) == 21
    for name in sk.__all__:
        parameter = inspect.signature(getattr(sk, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_missing_seqkit_on_path_raises_with_install_hint(tmp_path: Path, monkeypatch) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(BackendUnavailableError) as raised:
        sk.seqkit_range_records(source, str(tmp_path / "out.fa"), start=1, end=2, verbose=False)
    assert raised.value.hint and "seqkit" in raised.value.hint


def test_explicit_missing_seqkit_path_raises(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    with pytest.raises(BackendUnavailableError):
        sk.seqkit_range_records(
            source,
            str(tmp_path / "out.fa"),
            start=1,
            end=2,
            seqkit_path=str(tmp_path / "no-such-seqkit"),
            verbose=False,
        )


def test_only_whitelisted_options_are_accepted() -> None:
    sk._check_options("duplicate_records", ["--times=2"])
    with pytest.raises(ConfigurationError):
        sk._check_options("duplicate_records", ["--times=2", "--threads=8"])
    with pytest.raises(ConfigurationError):
        sk._check_options("duplicate_records", ["-s"])


def test_existing_output_is_never_overwritten(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    existing = _write(tmp_path / "out.fa", "keep me\n")
    with pytest.raises(ConfigurationError):
        sk.seqkit_duplicate_records(source, existing, times=2, verbose=False)
    assert Path(existing).read_text() == "keep me\n"


@pytest.mark.parametrize(
    "call",
    [
        lambda d: sk.seqkit_duplicate_records(d["in"], d["out"], times=0, verbose=False),
        lambda d: sk.seqkit_split_by_region(d["in"], d["dir"], region="1-5", verbose=False),
        lambda d: sk.seqkit_grep_motif(d["in"], d["out"], motifs=["AC,GT"], verbose=False),
        lambda d: sk.seqkit_concat_by_id([d["in"]], d["out"], fill="N", verbose=False),
        lambda d: sk.seqkit_flanking_sequences(
            d["in"], d["bed"], d["out"], up_stream=2, down_stream=2, only_flank=True, verbose=False
        ),
        lambda d: sk.seqkit_flanking_sequences(
            d["in"], d["bed"], d["out"], up_stream=2, feature_types=["gene"], verbose=False
        ),
        lambda d: sk.seqkit_fish_search(
            d["in"], d["in"], d["out"], alignment_params="4,-4", verbose=False
        ),
    ],
)
def test_invalid_arguments_are_rejected_before_running(tmp_path: Path, call) -> None:
    data = {
        "in": _write(tmp_path / "in.fa", FASTA),
        "bed": _write(tmp_path / "f.bed", BED),
        "out": str(tmp_path / "out.fa"),
        "dir": str(tmp_path / "split"),
    }
    with pytest.raises(ConfigurationError):
        call(data)


# ---------------------------------------------------------------------------
# SeqKit-backed behaviour
# ---------------------------------------------------------------------------


@requires_seqkit
def test_range_and_duplicate_report_counts_and_version(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    first = sk.seqkit_range_records(
        source, str(tmp_path / "range.fa"), start=1, end=2, verbose=False
    )
    assert first["status"] == "completed" and first["backend"] == "seqkit"
    assert first["record_count"] == 2
    assert _fasta_ids(tmp_path / "range.fa") == ["chr1", "chr2"]
    assert re.match(r"\d+\.\d+", first["version"] or "")
    json.dumps(first)

    last = sk.seqkit_range_records(
        source, str(tmp_path / "last.fa"), start=-1, end=-1, verbose=False
    )
    assert _fasta_ids(tmp_path / "last.fa") == ["chr3_a"]
    assert last["record_count"] == 1

    doubled = sk.seqkit_duplicate_records(source, str(tmp_path / "dup.fa"), times=2, verbose=False)
    assert doubled["record_count"] == 6


@requires_seqkit
def test_count_and_motif_filter(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    counted = sk.seqkit_count_matching_records(source, patterns=["chr1", "chr3_a"], verbose=False)
    assert counted["count"] == 2
    inverted = sk.seqkit_count_matching_records(
        source, patterns=["chr1", "chr3_a"], invert_match=True, verbose=False
    )
    assert inverted["count"] == 1
    motif = sk.seqkit_grep_motif(
        source, str(tmp_path / "motif.fa"), motifs=["ACGTAC"], verbose=False
    )
    assert motif["record_count"] == 2
    assert _fasta_ids(tmp_path / "motif.fa") == ["chr1", "chr3_a"]


@requires_seqkit
def test_statistics_digest_and_quality_parse(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    stats = sk.seqkit_extended_stats([source], str(tmp_path / "stats.tsv"), verbose=False)
    row = stats["records"][0]
    assert row["format"] == "FASTA" and row["num_seqs"] == 3 and "N50" in row
    digest = sk.seqkit_sequence_digest([source], str(tmp_path / "sum.tsv"), verbose=False)
    assert digest["records"][0]["digest"].startswith("seqkit.")
    fastq = _write(tmp_path / "reads.fq", FASTQ)
    quality = sk.seqkit_mean_read_quality(fastq, str(tmp_path / "mq.tsv"), verbose=False)
    assert [row["id"] for row in quality["records"]] == ["r1", "r2"]
    assert quality["records"][0]["avg_qual"] == pytest.approx(40, abs=0.01)
    histogram = sk.seqkit_read_feature_histogram(fastq, "ReadLen", bins=2, verbose=False)
    assert sum(row["Count"] for row in histogram["records"]) == 2


@requires_seqkit
def test_flank_split_concat_rename_replace_and_merge(tmp_path: Path) -> None:
    source = _write(tmp_path / "in.fa", FASTA)
    bed = _write(tmp_path / "f.bed", BED)
    flank = sk.seqkit_flanking_sequences(
        source, bed, str(tmp_path / "flank.fa"), up_stream=2, only_flank=True, verbose=False
    )
    sequence = "".join(
        line
        for line in (tmp_path / "flank.fa").read_text().splitlines()
        if not line.startswith(">")
    )
    assert len(sequence) == 2 and flank["record_count"] == 1
    assert Path(flank["index_file"]).name == "in.fa.seqkit.fai"

    split = sk.seqkit_split_by_region(source, str(tmp_path / "split"), region="1:2", verbose=False)
    assert len(split["records"]) == 2
    assert sum(row["record_count"] for row in split["records"]) == 3

    left = _write(tmp_path / "a.fa", ">chr1\nAAAA\n")
    right = _write(tmp_path / "b.fa", ">chr1\nCCCC\n")
    sk.seqkit_concat_by_id([left, right], str(tmp_path / "cat.fa"), verbose=False)
    assert "AAAACCCC" in (tmp_path / "cat.fa").read_text()

    dups = _write(tmp_path / "dup_ids.fa", ">a\nAC\n>a\nGT\n")
    sk.seqkit_rename_duplicate_ids(dups, str(tmp_path / "renamed.fa"), verbose=False)
    assert _fasta_ids(tmp_path / "renamed.fa") == ["a", "a_2"]

    sk.seqkit_replace_sequence_regex(
        source, str(tmp_path / "rep.fa"), pattern="N+", replacement="-", verbose=False
    )
    assert "ACGTACGT-ACGTTTGCA" in (tmp_path / "rep.fa").read_text()

    names = _write(
        tmp_path / "windows.txt",
        "chr1_sliding:1-1000\nchr1_sliding:501-1500\nchr2_sliding:1-1000\n",
    )
    merged = sk.seqkit_merge_slides_to_bed(
        names, str(tmp_path / "merged.bed"), min_overlap=499, verbose=False
    )
    assert merged["records"] == [
        {"chrom": "chr1", "start": 0, "end": 1500},
        {"chrom": "chr2", "start": 0, "end": 1000},
    ]


@requires_seqkit
def test_fish_and_fasta_for_fastq(tmp_path: Path) -> None:
    # SeqKit 2.14.0 panics on N in fish targets, so the target here is N-free.
    target = _write(tmp_path / "target.fa", ">chr1\nACGTACGTTTGCAGGA\n>chr2\nGGGGCCCCAAAATTTT\n")
    query = _write(tmp_path / "q.fa", ">q\nACGTACGT\n")
    hits = sk.seqkit_fish_search(target, query, str(tmp_path / "fish.tsv"), verbose=False)
    assert hits["records"] and any(row["Ref"] == "chr1" for row in hits["records"])
    fastq = _write(tmp_path / "reads.fq", FASTQ)
    # fa2fq needs the same ID and the same full sequence: r1 matches, a different ID or sequence does not.
    same = _write(tmp_path / "same.fa", ">r1\nACGTACGTAC\n>r2\nACGTACGTAC\n")
    matched = sk.seqkit_fastq_records_for_fasta(fastq, same, str(tmp_path / "fq.fq"), verbose=False)
    assert matched["record_count"] == 1


@requires_seqkit
def test_histogram_can_save_image(tmp_path: Path) -> None:
    fastq = _write(tmp_path / "reads.fq", FASTQ)
    image = tmp_path / "hist.png"
    result = sk.seqkit_read_feature_histogram(
        fastq, "GC", bins=2, image_path=str(image), verbose=False
    )
    assert image.is_file() and result["outputs"]["image"] == str(image)


@requires_seqkit
def test_bam_facades(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as handle:
        for index, start in enumerate((10, 100, 300)):
            segment = pysam.AlignedSegment()
            segment.query_name = f"r{index}"
            segment.query_sequence = "ACGT" * (5 + index)
            segment.query_qualities = pysam.qualitystring_to_array(
                "I" * len(segment.query_sequence)
            )
            segment.reference_id = 0
            segment.reference_start = start
            segment.cigarstring = f"{len(segment.query_sequence)}M"
            segment.mapping_quality = 60
            handle.write(segment)
    pysam.index(str(bam_path))
    stats = sk.seqkit_bam_alignment_stats(str(bam_path), verbose=False)
    assert stats["records"][0]["TotalRecords"] == 3
    indexed = sk.seqkit_bam_alignment_stats(str(bam_path), index_based=True, verbose=False)
    assert indexed["records"][0]["TotalRec"] == 3
    histogram = sk.seqkit_bam_feature_histogram(str(bam_path), "ReadLen", bins=2, verbose=False)
    assert sum(row["Count"] for row in histogram["records"]) == 3
    bundles = sk.seqkit_bam_bundles(str(bam_path), str(tmp_path / "bundles"), verbose=False)
    assert bundles["record_count"] == 3
    top = sk.seqkit_bam_top_records(
        str(bam_path), str(tmp_path / "top.bam"), field="ReadLen", top_size=2, verbose=False
    )
    with pysam.AlignmentFile(top["outputs"]["bam"], "rb") as handle:
        assert sum(1 for _ in handle) == 2
