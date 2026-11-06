from __future__ import annotations

import array
import random
import re

import pytest

pysam = pytest.importorskip("pysam")

from dnakit.alignment import pysam_p36_alignment_functions as m
from dnakit.exceptions import BackendExecutionError, ConfigurationError

REF = "".join(random.Random(1).choice("ACGT") for _ in range(1000))


@pytest.fixture(scope="module")
def bam(tmp_path_factory):
    folder = tmp_path_factory.mktemp("p36bam")
    rng = random.Random(2)
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )

    def make(name, flag, pos, cigar, seq=None, mate=None, tags=()):
        read = pysam.AlignedSegment(header)
        read.query_name, read.flag, read.reference_id, read.reference_start = name, flag, 0, pos
        read.mapping_quality = 60
        read.cigarstring = cigar
        length = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar))
        read.query_sequence = seq or "".join(rng.choice("ACGT") for _ in range(length))
        read.query_qualities = pysam.qualitystring_to_array("I" * length)
        if mate is not None:
            read.next_reference_id, read.next_reference_start, read.template_length = 0, mate, 200
        for tag in tags:
            read.set_tag(*tag)
        return read

    mismatch = list(REF[600:640])
    mismatch[10] = "A" if REF[610] != "A" else "C"
    reads = [
        make("p1", 99, 100, "50M", REF[100:150], mate=250),
        make("p1", 147, 250, "50M", REF[250:300], mate=100),
        make("s1", 0, 400, "20M100N20M", REF[400:420] + REF[520:540]),
        make(
            "m1",
            0,
            600,
            "40M",
            "".join(mismatch),
            tags=[("MD", "10" + REF[610] + "29", "Z"), ("NM", 1, "i")],
        ),
        make(
            "mod",
            0,
            700,
            "40M",
            "ACGTCGACGT" * 4,
            tags=[("MM", "C+m,0,1;", "Z"), ("ML", array.array("B", [200, 100]))],
        ),
        make("ins", 0, 800, "10M2I5D10M3S"),
    ]
    reads.sort(key=lambda r: r.reference_start)
    path = folder / "t.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(str(path))
    return path


def test_mate_alignments_follow_pysam(bam):
    result = m.fetch_mate_alignments(bam, "chr1", verbose=False)
    assert result["count"] == 2
    assert {p["mate"]["flag"] for p in result["pairs"]} == {"99", "147"}
    with pysam.AlignmentFile(str(bam)) as handle:
        first = next(r for r in handle.fetch("chr1") if r.query_name == "p1" and not r.is_reverse)
        assert result["pairs"][0]["mate"] == handle.mate(first).to_dict()
    assert "s1" in result["without_mate"]


def test_count_base_coverage_equals_pysam_and_validates(bam):
    result = m.count_base_coverage(bam, "chr1", 100, 130, quality_threshold=15, verbose=False)
    with pysam.AlignmentFile(str(bam)) as handle:
        expected = handle.count_coverage("chr1", 100, 130, quality_threshold=15)
    assert [result[b] for b in "ACGT"] == [list(map(int, col)) for col in expected]
    assert sum(sum(result[b]) for b in "ACGT") == 30
    nofilter = m.count_base_coverage(bam, "chr1", 100, 130, read_callback="nofilter", verbose=False)
    assert nofilter["read_callback"] == "nofilter"
    with pytest.raises(ConfigurationError):
        m.count_base_coverage(bam, "chr1", 130, 100, verbose=False)
    with pytest.raises(ConfigurationError):
        m.count_base_coverage(bam, "chr1", 0, 10, read_callback="x", verbose=False)


def test_spliced_junctions_fast_and_slow_agree(bam):
    fast = m.count_spliced_junctions(bam, verbose=False)
    slow = m.count_spliced_junctions(bam, "chr1", method="slow", verbose=False)
    assert (
        fast["junctions"]
        == slow["junctions"]
        == [{"intron_start": 420, "intron_end": 520, "supporting_reads": 1}]
    )
    assert m.count_spliced_junctions(bam, "chr1", 0, 300, verbose=False)["junction_count"] == 0


def test_pileup_columns_report_aligned_reads(bam):
    result = m.pileup_columns(
        bam, "chr1", 120, 125, truncate=True, include_read_names=True, verbose=False
    )
    assert [c["reference_pos"] for c in result["columns"]] == [120, 121, 122, 123, 124]
    assert all(c["query_names"] == ["p1"] and c["num_aligned"] == 1 for c in result["columns"])
    capped = m.pileup_columns(bam, "chr1", 100, 150, truncate=True, max_columns=10, verbose=False)
    assert capped["truncated"] and capped["column_count"] == 10
    strict = m.pileup_columns(
        bam, "chr1", 100, 110, truncate=True, min_base_quality=60, verbose=False
    )
    assert all(c["num_aligned"] == 0 for c in strict["columns"])
    assert all(c["nsegments"] == 1 for c in strict["columns"])


def test_find_alignments_by_name_uses_name_index(bam):
    result = m.find_alignments_by_name(bam, ["p1", "nobody"], verbose=False)
    assert len(result["alignments"]["p1"]) == 2 and result["alignments"]["nobody"] == []
    assert result["matched_names"] == 1


def test_header_conversions_round_trip(bam, tmp_path):
    info = m.sam_header_to_dict(bam, verbose=False)
    assert info["header"]["SQ"] == [{"SN": "chr1", "LN": 1000}]
    out = tmp_path / "h.sam"
    written = m.write_sam_header(out, header_dict=info["header"], verbose=False)
    assert out.read_text() == info["header_text"] == written["header_text"]
    refs = m.write_sam_header(
        tmp_path / "r.sam", reference_names=["a", "b"], reference_lengths=[10, 20], verbose=False
    )
    assert refs["header"]["SQ"] == [{"SN": "a", "LN": 10}, {"SN": "b", "LN": 20}]
    assert refs["header_text"] == "@SQ\tSN:a\tLN:10\n@SQ\tSN:b\tLN:20\n"
    with pytest.raises(ConfigurationError):
        m.write_sam_header(out, header_dict=info["header"], verbose=False)  # exists
    with pytest.raises(ConfigurationError):
        m.write_sam_header(
            tmp_path / "x.sam", reference_names=["a"], reference_lengths=[1, 2], verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.write_sam_header(tmp_path / "y.sam", verbose=False)


def test_alignment_records_round_trip_through_dicts(bam, tmp_path):
    dumped = m.alignment_records_to_dicts(bam, verbose=False)
    assert dumped["count"] == 6 and dumped["records"][0]["name"] == "p1"
    out = tmp_path / "copy.bam"
    result = m.write_alignments_from_dicts(dumped["records"], dumped["header"], out, verbose=False)
    assert result["count"] == 6
    with pysam.AlignmentFile(str(bam)) as a, pysam.AlignmentFile(str(out)) as b:
        assert [r.to_dict() for r in a.fetch(until_eof=True)] == [
            r.to_dict() for r in b.fetch(until_eof=True)
        ]
    with pytest.raises(ConfigurationError):
        m.write_alignments_from_dicts([], dumped["header"], tmp_path / "bad.txt", verbose=False)


def test_per_read_geometry_matches_pysam(bam):
    with pysam.AlignmentFile(str(bam)) as handle:
        reads = list(handle.fetch("chr1"))
    pairs = m.read_aligned_pairs(bam, "chr1", matches_only=True, verbose=False)["records"]
    blocks = m.read_aligned_blocks(bam, "chr1", verbose=False)["records"]
    stats = m.read_cigar_stats(bam, "chr1", verbose=False)["records"]
    positions = m.read_reference_positions(bam, "chr1", verbose=False)["records"]
    lengths = m.infer_read_lengths(bam, "chr1", verbose=False)["records"]
    for read, p, b, s, pos, ln in zip(
        reads, pairs, blocks, stats, positions, lengths, strict=False
    ):
        assert p["pairs"] == [list(x) for x in read.get_aligned_pairs(matches_only=True)]
        assert b["blocks"] == [list(x) for x in read.get_blocks()]
        assert s["bases"]["M"] == read.get_cigar_stats()[0][0]
        assert pos["positions"] == read.get_reference_positions()
        assert (
            ln["query_length"] == read.infer_query_length()
            and ln["read_length"] == read.infer_read_length()
        )
    spliced = next(r for r in blocks if r["query_name"] == "s1")
    assert spliced["blocks"] == [[400, 420], [520, 540]]
    indel = next(r for r in stats if r["query_name"] == "ins")
    assert indel["bases"]["I"] == 2 and indel["bases"]["D"] == 5 and indel["bases"]["S"] == 3
    soft = next(r for r in lengths if r["query_name"] == "ins")
    assert soft["query_length"] == 25 and soft["stored_length"] == 25
    full = m.read_reference_positions(bam, "chr1", full_length=True, verbose=False)["records"]
    assert None in next(r for r in full if r["query_name"] == "ins")["positions"]


def test_reference_overlap_and_truncation(bam):
    result = m.read_reference_overlap(bam, "chr1", 100, 120, verbose=False)
    assert result["records"] == [{"query_name": "p1", "flag": 99, "overlap_bases": 20}]
    capped = m.read_aligned_blocks(bam, max_records=2, verbose=False)
    assert capped["truncated"] and capped["count"] == 2
    with pytest.raises(ConfigurationError):
        m.read_reference_overlap(bam, "chr1", None, None, verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_aligned_pairs(bam, None, 5, 10, verbose=False)


def test_reference_sequence_from_md_tag(bam):
    result = m.reconstruct_read_reference_sequence(bam, "chr1", verbose=False)
    assert [r["query_name"] for r in result["records"]] == ["m1"]
    assert result["records"][0]["reference_sequence"].upper() == REF[600:640]
    assert result["records"][0]["reference_sequence"][10] == REF[610].lower()
    assert "p1" in result["skipped_without_md"]


def test_base_modifications_both_orientations(bam):
    result = m.read_base_modifications(bam, "chr1", verbose=False)["records"]
    mod = next(r for r in result if r["query_name"] == "mod")["modifications"]
    assert mod == [{"base": "C", "strand": 0, "modification": "m", "calls": [[1, 200], [7, 100]]}]
    forward = m.read_base_modifications(bam, "chr1", forward=True, verbose=False)["records"]
    assert next(r for r in forward if r["query_name"] == "mod")["modifications"] == mod
    assert next(r for r in result if r["query_name"] == "p1")["modifications"] == []


def test_missing_file_and_bad_arguments(tmp_path, bam):
    with pytest.raises(ConfigurationError):
        m.read_aligned_blocks(tmp_path / "none.bam", verbose=False)
    with pytest.raises(BackendExecutionError):
        m.read_aligned_blocks(bam, "chrX", verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_aligned_blocks(bam, max_records=0, verbose=False)
    with pytest.raises(ConfigurationError):
        m.find_alignments_by_name(bam, "p1", verbose=False)
