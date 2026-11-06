"""Tests for the fastp facade; program-dependent tests skip when fastp is not installed."""

from __future__ import annotations

import random
import shutil

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import fastp_p14_cli_functions as fp

needs_fastp = pytest.mark.skipif(shutil.which("fastp") is None, reason="fastp is not installed")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
ADAPTER = "AGATCGGAAGAGC"


@pytest.fixture()
def pairs(tmp_path):
    rng = random.Random(9)
    first, second = [], []
    for index in range(300):
        fragment = "".join(rng.choice("ACGT") for _ in range(120))
        read1 = fragment[:100]
        read2 = fragment[::-1].translate(COMPLEMENT)[:100]
        if index < 100:  # short inserts: reads run into the adapter
            insert = fragment[:70]
            read1 = (insert + ADAPTER + "ACGTACGTAC" * 3)[:100]
            read2 = (insert[::-1].translate(COMPLEMENT) + ADAPTER + "ACGTACGTAC" * 3)[:100]
        first.append(f"@p{index}/1\n{read1}\n+\n{'I' * 100}\n")
        second.append(f"@p{index}/2\n{read2}\n+\n{'I' * 100}\n")
    (tmp_path / "r1.fq").write_text("".join(first))
    (tmp_path / "r2.fq").write_text("".join(second))
    (tmp_path / "adapters.fa").write_text(f">illumina\n{ADAPTER}\n")
    return tmp_path


def test_validation_and_missing_program(pairs):
    with pytest.raises(BackendUnavailableError):
        fp.fastp_process_reads(
            str(pairs / "r1.fq"),
            str(pairs / "o1.fq"),
            str(pairs / "r.json"),
            fastp_path=str(pairs / "nofastp"),
        )
    with pytest.raises(ConfigurationError):
        fp.fastp_process_reads(
            str(pairs / "r1.fq"),
            str(pairs / "o1.fq"),
            str(pairs / "r.json"),
            read2=str(pairs / "r2.fq"),
        )
    with pytest.raises(ConfigurationError):
        fp.fastp_process_reads(
            str(pairs / "r1.fq"), str(pairs / "o1.fq"), str(pairs / "r.json"), correction=True
        )
    with pytest.raises(ConfigurationError):
        fp.fastp_process_reads(
            str(pairs / "r1.fq"),
            str(pairs / "o1.fq"),
            str(pairs / "r.json"),
            umi=True,
            umi_location="read1",
        )
    with pytest.raises(ConfigurationError):
        fp.fastp_process_reads(
            str(pairs / "r1.fq"), str(pairs / "o1.fq"), str(pairs / "r.json"), umi_length=4
        )


@needs_fastp
def test_adapter_fasta_trimming_and_reports(pairs):
    result = fp.fastp_process_reads(
        str(pairs / "r1.fq"),
        str(pairs / "o1.fq"),
        str(pairs / "r.json"),
        read2=str(pairs / "r2.fq"),
        output_read2=str(pairs / "o2.fq"),
        html_report=str(pairs / "r.html"),
        adapter_fasta=str(pairs / "adapters.fa"),
        threads=2,
        verbose=False,
    )
    assert result["status"] == "completed" and result["version"]
    assert result["summary"]["adapter_cutting"]["adapter_trimmed_reads"] >= 150
    lengths = {len(line) for line in (pairs / "o1.fq").read_text().split("\n")[1::4] if line}
    assert min(lengths) <= 75 and (pairs / "r.html").is_file()


@needs_fastp
def test_low_complexity_filter_removes_homopolymers(pairs):
    low = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    with (pairs / "r1.fq").open("a") as handle:
        handle.write(f"@low\n{low}\n+\n{'I' * len(low)}\n")
    result = fp.fastp_process_reads(
        str(pairs / "r1.fq"),
        str(pairs / "o1.fq"),
        str(pairs / "r.json"),
        low_complexity_filter=True,
        complexity_threshold=30,
        trim_poly_x=False,
        verbose=False,
    )
    assert result["summary"]["filtering_result"]["low_complexity_reads"] == 1
    assert "@low" not in (pairs / "o1.fq").read_text()


@needs_fastp
def test_correction_umi_overrepresentation_and_index_filter(pairs):
    (pairs / "index.txt").write_text("ACGTAC\n")
    result = fp.fastp_process_reads(
        str(pairs / "r1.fq"),
        str(pairs / "o1.fq"),
        str(pairs / "r.json"),
        read2=str(pairs / "r2.fq"),
        output_read2=str(pairs / "o2.fq"),
        correction=True,
        umi=True,
        umi_location="read1",
        umi_length=4,
        umi_prefix="UMI",
        overrepresentation_analysis=True,
        verbose=False,
    )
    assert result["status"] == "completed"
    names = [line for line in (pairs / "o1.fq").read_text().split("\n")[::4] if line]
    assert names and all(":UMI_" in name for name in names)
    sequence = (pairs / "o1.fq").read_text().split("\n")[1]
    assert len(sequence) <= 96  # the 4 UMI bases are removed from the read


@needs_fastp
def test_index_barcode_filter_runs(pairs):
    (pairs / "index.txt").write_text("ACGTAC\n")
    result = fp.fastp_process_reads(
        str(pairs / "r1.fq"),
        str(pairs / "o1.fq"),
        str(pairs / "r.json"),
        filter_by_index1=str(pairs / "index.txt"),
        index_filter_threshold=1,
        verbose=False,
    )
    assert result["status"] == "completed" and "--filter_by_index1" in result["command"]
