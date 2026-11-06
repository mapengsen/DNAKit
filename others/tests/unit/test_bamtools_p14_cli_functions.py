"""Tests for the BamTools facades; program-dependent tests skip when BamTools or pysam is missing."""

from __future__ import annotations

import random
import shutil

import pytest

from dnakit.alignment import bamtools_p14_cli_functions as bt
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pysam = pytest.importorskip("pysam")
needs_bamtools = pytest.mark.skipif(
    shutil.which("bamtools") is None, reason="BamTools is not installed"
)


@pytest.fixture()
def bam(tmp_path):
    """60 reads on two references: 40 mapped (half reverse, MAPQ 60/10, tag XX), 20 unmapped, sorted and indexed."""
    rng = random.Random(3)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 5000}, {"SN": "chr2", "LN": 3000}],
    }
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for index in range(60):
            read = pysam.AlignedSegment(out.header)
            read.query_name = f"r{index}"
            read.query_sequence = "".join(rng.choice("ACGT") for _ in range(50))
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            if index < 40:
                read.reference_id = index % 2
                read.reference_start = 100 + 10 * index
                read.cigar = [(0, 50)]
                read.flag = 16 if index % 4 == 1 else 0
                read.mapping_quality = 60 if index % 3 else 10
                read.set_tag("XX", "a" if index % 2 else "b")
            else:
                read.flag = 4
                read.reference_id = -1
            out.write(read)
    path = tmp_path / "reads.bam"
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def test_missing_program_raises_install_hint(bam, tmp_path):
    with pytest.raises(BackendUnavailableError):
        bt.bamtools_count_alignments([str(bam)], bamtools_path=str(tmp_path / "absent"))


def test_argument_validation(bam, tmp_path):
    with pytest.raises(ConfigurationError):
        bt.bamtools_filter_alignments([str(bam)], str(tmp_path / "o.bam"))
    with pytest.raises(ConfigurationError):
        bt.bamtools_filter_alignments([str(bam)], str(tmp_path / "o.bam"), map_quality="50; rm")
    with pytest.raises(ConfigurationError):
        bt.bamtools_split_alignments(str(bam), str(tmp_path / "s"), by="tag")
    with pytest.raises(ConfigurationError):
        bt.bamtools_resolve_pairs(
            str(bam), str(tmp_path / "st"), mode="makeStats", output=str(tmp_path / "o.bam")
        )


@needs_bamtools
def test_count_filter_and_sample(bam, tmp_path):
    assert bt.bamtools_count_alignments([str(bam)], verbose=False)["count"] == 60
    assert bt.bamtools_count_alignments([str(bam)], region="chr1", verbose=False)["count"] == 20
    filtered = bt.bamtools_filter_alignments(
        [str(bam)],
        str(tmp_path / "f.bam"),
        map_quality=">=60",
        is_reverse_strand=False,
        is_mapped=True,
        verbose=False,
    )
    assert filtered["status"] == "completed"
    expected = sum(
        1
        for read in pysam.AlignmentFile(str(bam))
        if not read.is_unmapped and read.mapping_quality >= 60 and not read.is_reverse
    )
    assert (
        bt.bamtools_count_alignments([str(tmp_path / "f.bam")], verbose=False)["count"] == expected
    )
    sample = bt.bamtools_sample_alignments(
        [str(bam)], str(tmp_path / "r.bam"), count=15, seed=1, verbose=False
    )
    assert sample["status"] == "completed"
    assert bt.bamtools_count_alignments([str(tmp_path / "r.bam")], verbose=False)["count"] == 15


@needs_bamtools
def test_convert_json_and_yaml(bam, tmp_path):
    import json

    assert (
        bt.bamtools_convert_alignments(
            [str(bam)], str(tmp_path / "a.json"), output_format="json", region="chr2", verbose=False
        )["status"]
        == "completed"
    )
    records = [
        json.loads(line) for line in (tmp_path / "a.json").read_text().splitlines() if line.strip()
    ]
    assert len(records) == 20 and records[0]["reference"] == "chr2"
    assert (
        bt.bamtools_convert_alignments(
            [str(bam)], str(tmp_path / "a.yaml"), output_format="yaml", region="chr2", verbose=False
        )["status"]
        == "completed"
    )
    assert (tmp_path / "a.yaml").read_text().count("RefName: chr2") == 20


@needs_bamtools
def test_split_variants(bam, tmp_path):
    mapped = bt.bamtools_split_alignments(str(bam), str(tmp_path / "m"), by="mapped", verbose=False)
    assert sorted(path.rsplit(".", 2)[-2] for path in mapped["outputs"]) == ["MAPPED", "UNMAPPED"]
    paired = bt.bamtools_split_alignments(str(bam), str(tmp_path / "p"), by="paired", verbose=False)
    assert paired["outputs"][0].endswith("SINGLE_END.bam")
    reference = bt.bamtools_split_alignments(
        str(bam), str(tmp_path / "r"), by="reference", reference_prefix="REF-", verbose=False
    )
    assert [path.split(".")[-2] for path in reference["outputs"]] == [
        "REF-chr1",
        "REF-chr2",
        "REF-unmapped",
    ]
    tag = bt.bamtools_split_alignments(
        str(bam), str(tmp_path / "t"), by="tag", tag="XX", verbose=False
    )
    assert len(tag["outputs"]) == 2
    assert bt.bamtools_count_alignments([tag["outputs"][0]], verbose=False)["count"] == 20


@needs_bamtools
def test_index_revert_and_resolve(bam, tmp_path):
    index = bt.bamtools_index_bam_bti(str(bam), verbose=False)
    assert index["outputs"][0].endswith(".bti")
    bt.bamtools_revert_alignments([str(bam)], str(tmp_path / "rv.bam"), verbose=False)
    assert bt.bamtools_count_alignments([str(tmp_path / "rv.bam")], verbose=False)["count"] == 60
    resolved = bt.bamtools_resolve_pairs(
        str(bam), str(tmp_path / "st.txt"), output=str(tmp_path / "rs.bam"), verbose=False
    )
    assert resolved["status"] == "completed" and resolved["mode"] == "twoPass"
    stats = bt.bamtools_resolve_pairs(
        str(bam), str(tmp_path / "st2.txt"), mode="makeStats", verbose=False
    )
    assert stats["status"] == "completed"
    marked = bt.bamtools_resolve_pairs(
        str(bam),
        str(tmp_path / "st2.txt"),
        output=str(tmp_path / "rs2.bam"),
        mode="markPairs",
        verbose=False,
    )
    assert marked["status"] == "completed"
