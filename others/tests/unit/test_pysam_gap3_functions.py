"""Tests for dnakit.bioinformatics.pysam_gap3_functions (padded alignments and copy-number calling).

Inputs are written to pytest's tmp_path; nothing is downloaded.
"""

from __future__ import annotations

import inspect
import random
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError

pysam = pytest.importorskip("pysam")

from dnakit.bioinformatics import pysam_gap3_functions as facades

PADDED_REFERENCE = "ACG*TTA*CG"
READ_SEQUENCE = "ACGTTACG"


@pytest.fixture()
def padded_inputs(tmp_path: Path) -> dict[str, str]:
    reference = tmp_path / "padded.fa"
    reference.write_text(f">chr1\n{PADDED_REFERENCE}\n", encoding="utf-8")
    pysam.faidx(str(reference))
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(PADDED_REFERENCE)}],
    }
    alignment = tmp_path / "padded.bam"
    read = pysam.AlignedSegment()
    read.query_name = "r1"
    read.query_sequence = READ_SEQUENCE
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigar = [(0, len(READ_SEQUENCE))]
    read.query_qualities = pysam.qualitystring_to_array("I" * len(READ_SEQUENCE))
    with pysam.AlignmentFile(str(alignment), "wb", header=header) as handle:
        handle.write(read)
    return {"reference": str(reference), "alignment": str(alignment), "root": str(tmp_path)}


@pytest.fixture()
def cnv_vcf(tmp_path: Path) -> Path:
    rng = random.Random(3)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=10000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B allele frequency">',
        '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R ratio">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    for index in range(60):
        position = 100 * (index + 1)
        if 20 <= index < 40:
            genotype, baf, lrr = "0/0", 0.0, -0.6
        else:
            genotype, baf, lrr = (
                "0/1",
                rng.choice([0.5, 0.48, 0.52]),
                round(rng.uniform(-0.05, 0.05), 3),
            )
        lines.append(f"chr1\t{position}\t.\tA\tG\t50\tPASS\t.\tGT:BAF:LRR\t{genotype}:{baf}:{lrr}")
    path = tmp_path / "cnv.vcf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_depad_alignments_removes_pad_columns(padded_inputs):
    out = Path(padded_inputs["root"]) / "depad.bam"
    result = facades.depad_alignments(
        padded_inputs["alignment"], padded_inputs["reference"], str(out), verbose=False
    )
    assert result["status"] == "completed"
    assert result["output_format"] == "BAM"
    assert Path(result["outputs"][0]).is_file()
    [record] = result["records"]
    assert record["cigar"] == "3M1I3M1I"
    assert record["sequence"] == READ_SEQUENCE
    assert record["reference_start"] == 1


def test_depad_alignments_can_write_sam(padded_inputs):
    out = Path(padded_inputs["root"]) / "depad.sam"
    result = facades.depad_alignments(
        padded_inputs["alignment"],
        padded_inputs["reference"],
        str(out),
        output_sam=True,
        verbose=False,
    )
    assert result["output_format"] == "SAM"
    assert out.read_text(encoding="utf-8").startswith("@HD")


def test_pad_to_unpad_alignments_matches_depad(padded_inputs):
    root = Path(padded_inputs["root"])
    depad = facades.depad_alignments(
        padded_inputs["alignment"], padded_inputs["reference"], str(root / "d.bam"), verbose=False
    )
    legacy = facades.pad_to_unpad_alignments(
        padded_inputs["alignment"], padded_inputs["reference"], str(root / "p.bam"), verbose=False
    )
    assert legacy["records"] == depad["records"]
    assert legacy["command"][1] == "pad2unpad"


def test_depad_rejects_existing_output_and_missing_inputs(padded_inputs):
    root = Path(padded_inputs["root"])
    existing = root / "exists.bam"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        facades.depad_alignments(
            padded_inputs["alignment"], padded_inputs["reference"], str(existing), verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.depad_alignments(
            str(root / "none.bam"), padded_inputs["reference"], str(root / "o.bam"), verbose=False
        )


def test_depad_reports_samtools_failure(padded_inputs, tmp_path):
    bad_reference = tmp_path / "plain.fa"
    bad_reference.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    with pytest.raises(BackendExecutionError):
        facades.depad_alignments(
            padded_inputs["alignment"], str(bad_reference), str(tmp_path / "out.bam"), verbose=False
        )


def test_call_copy_number_states_reports_sites_and_segments(cnv_vcf, tmp_path):
    result = facades.call_copy_number_states(
        str(cnv_vcf), str(tmp_path / "cnv_out"), query_sample="S1", verbose=False
    )
    assert result["status"] == "completed"
    assert len(result["sites"]) == 60
    assert result["sites"][0]["chromosome"] == "chr1"
    assert all(len(site["probabilities"]) == 4 for site in result["sites"])
    assert result["segments"] and result["segments"][0]["chromosome"] == "chr1"
    assert {"cn", "dat", "summary"} <= set(result["outputs"])
    assert all(Path(path).is_file() for path in result["outputs"].values())


def test_call_copy_number_states_needs_new_output_directory(cnv_vcf, tmp_path):
    existing = tmp_path / "already"
    existing.mkdir()
    with pytest.raises(ConfigurationError):
        facades.call_copy_number_states(
            str(cnv_vcf), str(existing), query_sample="S1", verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.call_copy_number_states(
            str(cnv_vcf), str(tmp_path / "x"), query_sample="-S1", verbose=False
        )


def test_missing_pysam_gives_install_hint(monkeypatch, padded_inputs, tmp_path):
    import importlib

    real = importlib.import_module

    def blocked(name, package=None):
        if name == "pysam":
            raise ImportError("blocked for test")
        return real(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", blocked)
    from dnakit.exceptions import BackendUnavailableError

    with pytest.raises(BackendUnavailableError) as info:
        facades.depad_alignments(
            padded_inputs["alignment"],
            padded_inputs["reference"],
            str(tmp_path / "o.bam"),
            verbose=False,
        )
    assert "bio-reads" in (info.value.hint or "")
