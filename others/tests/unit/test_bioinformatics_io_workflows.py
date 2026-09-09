"""Roundtrips, coordinates, molecule geometry, and actual subprocess contracts."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.bioinformatics import alignment, cloning, reads, workflows
from dnakit.bioinformatics._common import ExecutionConfig, run_to_file
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
    InputFormatError,
)


def test_multiline_fastq_phred_roundtrip_and_malformed_rejection(tmp_path: Path) -> None:
    path = tmp_path / "reads.fq"
    path.write_text("@a/1 sample\nAC\nGT\n+\n@@\n++\n@b/1\nTT\n+\nII\n")
    records = tuple(reads.read_fastx(path))
    assert records[0].sequence == "ACGT"
    assert records[0].qualities == (31, 31, 10, 10)
    converted = tmp_path / "converted.fq.gz"
    assert reads.write_fastx(records, converted, phred_offset=64) == 2
    assert tuple(reads.read_fastx(converted, phred_offset=64)) == records
    malformed = tmp_path / "bad.fq"
    malformed.write_text("@a\nACGT\n+\nII\n")
    with pytest.raises(InputFormatError):
        tuple(reads.read_fastx(malformed))
    with pytest.raises(FileExistsError):
        reads.write_fastx(records, converted)


def test_pairing_out_of_order_orphans_duplicates_and_bounds() -> None:
    a = [reads.Read("a/1", "AC"), reads.Read("b/1", "AC"), reads.Read("c/1", "AC")]
    b = [reads.Read("b/2", "GT"), reads.Read("a/2", "GT"), reads.Read("d/2", "GT")]
    pairs = {pair.id: pair for pair in reads.pair_reads(a, b)}
    assert pairs["a"].first == a[0] and pairs["a"].second == b[1]
    assert pairs["c"].second is None and pairs["d"].first is None
    with pytest.raises(InputFormatError):
        tuple(reads.pair_reads([a[0], a[0]], b))
    with pytest.raises(ConfigurationError):
        tuple(reads.pair_reads(a, b, max_pending=1))


def test_recovery_does_not_invent_quality_or_bases(tmp_path: Path) -> None:
    source = tmp_path / "broken.fq"
    source.write_text("junk\n@a\nAC\n+\nII\n@bad\nCCC\n+\nI\n@b\nGT\n+\n@@\n")
    target = tmp_path / "fixed.fq"
    result = reads.repair_fastq(source, target)
    assert result.recovered_reads == 2
    assert result.discarded_lines == 5
    assert [record.id for record in reads.read_fastx(target)] == ["a", "b"]


def test_genomic_interval_coordinates_reference_and_reverse_strand(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nAACCGGTT\n")
    reads.index_reference(reference)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\tsource\texon\t2\t5\t.\t-\t.\tgene_id "gene"; transcript_id "tx";\n')
    (region,) = reads.read_regions(gtf, format="gtf")
    assert (region.start, region.end) == (1, 5)
    with reads.IndexedReference(reference) as handle:
        record = handle.fetch(region)
        assert record.sequence == "CGGT"
        assert record.id == "tx"
        with pytest.raises(ConfigurationError):
            handle.fetch(reads.Region("chr1", 0, 10))


def test_bam_coverage_handles_deletion_and_filters_duplicates(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    np = pytest.importorskip("numpy")
    target = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 8}]}
    with pysam.AlignmentFile(str(target), "wb", header=header) as handle:
        for flag in (0, 1024):
            read = pysam.AlignedSegment()
            read.query_name = "read" + str(flag)
            read.flag = flag
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.query_sequence = "ACGT"
            read.query_qualities = pysam.qualitystring_to_array("IIII")
            read.cigar = [(0, 2), (2, 1), (0, 2)]
            handle.write(read)
    pysam.index(str(target))
    result = reads.coverage_consensus(target, reads.Region("chr1", 0, 6))
    np.testing.assert_array_equal(result["depth"], [1, 1, 0, 1, 1, 0])
    assert result["consensus"] == "ACNGTN"
    assert reads.alignment_summary(target)["duplicates"] == 1


def test_duplex_staggered_cut_ligation_and_orientation() -> None:
    original = cloning.Duplex.from_sequence("GGAATTCC")
    left, right = original.cut(2, 6)
    assert left.ends["right"]["sequence"] == "AATT"
    assert right.ends["left"]["sequence"] == "AATT"
    assert left.ligate(right) == original
    assert left.reverse_complement().reverse_complement() == left
    with pytest.raises(ConfigurationError):
        cloning.Duplex("AAAA", "AAAA")
    with pytest.raises(ConfigurationError):
        cloning.Duplex("AAAA", "TTT", circular=True)
    plan = cloning.CloningPlan()
    plan.add("template", original)
    product = plan.apply("digest", lambda molecule: molecule.cut(2, 6), parents=["template"])
    assert product.parents == ("template",)
    assert product.product == (left, right)


def test_native_dag_failure_propagation_and_region_boolean_queries() -> None:
    def fail() -> None:
        raise RuntimeError("expected failure")

    results = workflows.run_workflow(
        [
            workflows.Step("failure", fail),
            workflows.Step("blocked", lambda value: value, {"value": "failure"}),
            workflows.Step("independent", lambda: 7),
        ],
        fail_fast=False,
    )
    assert {result.name: result.status for result in results} == {
        "failure": "failed",
        "blocked": "blocked",
        "independent": "completed",
    }
    a = [reads.Region("chr1", 0, 10)]
    b = [reads.Region("chr1", 2, 4), reads.Region("chr1", 6, 12)]
    difference = workflows.query_regions(a, b, operation="difference")
    assert [(item.start, item.end) for item in difference] == [(0, 2), (4, 6)]
    union = workflows.query_regions(a, b, operation="union")
    assert [(item.start, item.end) for item in union] == [(0, 12)]


def test_alignment_format_roundtrips(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    source = alignment.MultipleAlignment.from_sequences({"a": "AC-G", "b": "ACTG"})
    for format in ("fasta", "clustal", "stockholm", "phylip-relaxed"):
        path = tmp_path / format
        alignment.write_alignment(source, path, format=format)
        assert alignment.read_alignment(path, format=format) == source
    assert source.distances()[0, 1] == 0
    assert source.distances(gap_policy="state")[0, 1] == 0.25


def test_external_execution_separates_data_and_diagnostics_and_checks_exit(tmp_path: Path) -> None:
    output = tmp_path / "out.txt"
    log = run_to_file(
        sys.executable,
        ["-c", "import sys; print('DATA'); print('LOG', file=sys.stderr)"],
        output,
        ExecutionConfig(),
    )
    assert output.read_text() == "DATA\n"
    assert log == "LOG\n"
    with pytest.raises(BackendExecutionError):
        run_to_file(
            sys.executable,
            ["-c", "raise SystemExit(2)"],
            tmp_path / "failed.txt",
            ExecutionConfig(),
        )
    with pytest.raises(BackendTimeoutError):
        run_to_file(
            sys.executable,
            ["-c", "import time; time.sleep(10)"],
            tmp_path / "timeout.txt",
            ExecutionConfig(timeout_seconds=0.05),
        )
    with pytest.raises(BackendExecutionError):
        run_to_file(
            sys.executable,
            ["-c", "print('x'*10000)"],
            tmp_path / "large.txt",
            ExecutionConfig(max_artifact_bytes=10),
        )


def test_mafft_process_contract_and_id_restoration(tmp_path: Path) -> None:
    pytest.importorskip("Bio")
    program = tmp_path / "fake_mafft"
    program.write_text(
        f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\n"
        "print(Path(sys.argv[-1]).read_text(), end='')\n"
        "print('diagnostic', file=sys.stderr)\n"
    )
    program.chmod(0o755)
    result = alignment.align_sequences(
        {"my_a": "ACGT", "my_b": "ACTT"}, tmp_path / "alignment", program=str(program)
    )
    assert result.ids == ("my_a", "my_b")
    assert result.sequences == ("ACGT", "ACTT")
    assert (tmp_path / "alignment" / "alignment.log").read_text() == "diagnostic\n"


def test_native_sanger_mapping_and_chromatogram(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    pytest.importorskip("matplotlib")
    matplotlib = pytest.importorskip("matplotlib")

    matplotlib.use("Agg")
    from dnakit.bioinformatics.cloning import plot_chromatogram, sanger_consensus

    sequence = "ACGTTGCAAGTC"
    reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    traces = [
        {"id": "forward", "sequence": sequence, "qualities": (40,) * len(sequence)},
        {"id": "reverse", "sequence": reverse, "qualities": (40,) * len(sequence)},
    ]
    result = sanger_consensus(traces, sequence)
    assert result["consensus"] == sequence
    np.testing.assert_array_equal(result["depth"], 2)
    assert result["mappings"][1]["strand"] == "-"
    trace = {
        "sequence": "AC",
        "peak_positions": [1, 3],
        "channels": {base: [0, 1, 0, 2, 0] for base in "ACGT"},
    }
    axes = plot_chromatogram(trace)
    axes.figure.savefig(tmp_path / "trace.png")
    assert (tmp_path / "trace.png").stat().st_size > 0
