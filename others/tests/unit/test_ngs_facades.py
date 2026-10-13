"""JSON and Agent contracts for the formal NGS command facades."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.ngs import (
    compress_index_vcf_htslib,
    map_reads,
    mark_duplicates_picard,
    phase_variants_whatshap,
    sample_sequences_seqtk,
    sort_bam_bamtools,
)
from dnakit.tools import default_tool_registry


def _fake_bwa(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "if 'mem' in sys.argv:\n"
        "    print('@r\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACGT\\tIIII')\n"
        "    print('diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_picard(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "args = dict(item.split('=', 1) for item in sys.argv[2:] if '=' in item)\n"
        "Path(args['O']).write_bytes(b'BAM')\n"
        "Path(args['M']).write_text('## metrics\\n', encoding='ascii')\n"
        "print('picard diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_seqtk(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "print('@sample\\nACGT\\n+\\nIIII')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_bgzip(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "print(Path(sys.argv[sys.argv.index('-c') + 1]).read_text(), end='')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_tabix(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[-1] + '.tbi').write_bytes(b'INDEX')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_whatshap(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[sys.argv.index('--output') + 1]).write_text('##fileformat=VCFv4.2\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_bamtools(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[sys.argv.index('-out') + 1]).write_bytes(b'SORTED_BAM')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_map_reads_facade_serializes_existing_ngs_backend(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fq"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    reads.write_text("@r\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_bwa(tmp_path / "fake-bwa")

    result = map_reads(
        reference,
        reads,
        tmp_path / "mapped",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["operation"] == "map_reads"
    assert result["backend"] == "bwa"
    artifacts = result["result"]["artifacts"]
    assert "aligned.sam" in artifacts
    assert "aligned.sam.stderr" in artifacts
    assert "@r" in Path(artifacts["aligned.sam"]).read_text(encoding="ascii")


def test_map_reads_facade_supports_bwa_mem2_alias(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fq"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    reads.write_text("@r\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_bwa(tmp_path / "fake-bwa-mem2")

    result = map_reads(
        reference,
        reads,
        tmp_path / "mapped",
        method="bwa-mem2",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "bwa-mem2"
    assert "aligned.sam" in result["result"]["artifacts"]


def test_picard_mark_duplicates_facade_writes_bam_and_metrics(tmp_path: Path) -> None:
    source = tmp_path / "aligned.bam"
    source.write_bytes(b"BAM")
    executable = _fake_picard(tmp_path / "fake-picard")

    result = mark_duplicates_picard(
        source,
        tmp_path / "deduplicated",
        remove_duplicates=True,
        max_records_in_ram=1_000,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "marked_duplicates.bam" in result["result"]["artifacts"]
    assert "duplicate_metrics.txt" in result["result"]["artifacts"]


def test_seqtk_sampling_facade_captures_sampled_reads(tmp_path: Path) -> None:
    source = tmp_path / "reads.fq"
    source.write_text("@r\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_seqtk(tmp_path / "fake-seqtk")

    result = sample_sequences_seqtk(
        source,
        tmp_path / "sampled",
        sample_size=0.5,
        seed=7,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "sampled.fastx" in result["result"]["artifacts"]


def test_htslib_facade_compresses_and_indexes_vcf(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    bgzip = _fake_bgzip(tmp_path / "fake-bgzip")
    tabix = _fake_tabix(tmp_path / "fake-tabix")

    result = compress_index_vcf_htslib(
        source,
        tmp_path / "compressed",
        bgzip_program=str(bgzip),
        tabix_program=str(tabix),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "variants.vcf.gz" in result["result"]["artifacts"]
    assert "variants.vcf.gz.tbi" in result["result"]["artifacts"]


def test_whatshap_facade_writes_phased_vcf(tmp_path: Path) -> None:
    variants = tmp_path / "variants.vcf"
    alignment = tmp_path / "reads.bam"
    reference = tmp_path / "reference.fa"
    variants.write_text("##fileformat=VCFv4.2\n", encoding="ascii")
    alignment.write_bytes(b"BAM")
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    executable = _fake_whatshap(tmp_path / "fake-whatshap")

    result = phase_variants_whatshap(
        variants,
        alignment,
        reference,
        tmp_path / "phased",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "phased.vcf" in result["result"]["artifacts"]


def test_bamtools_facade_writes_sorted_bam(tmp_path: Path) -> None:
    source = tmp_path / "source.bam"
    source.write_bytes(b"BAM")
    executable = _fake_bamtools(tmp_path / "fake-bamtools")

    result = sort_bam_bamtools(
        source,
        tmp_path / "sorted",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "sorted.bam" in result["result"]["artifacts"]


def test_ngs_facades_are_registered_as_confirmed_external_tools() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.ngs.map_reads",
        "dnakit.ngs.prepare_bam",
        "dnakit.ngs.mark_duplicates_picard",
        "dnakit.ngs.call_variants",
        "dnakit.ngs.compress_index_vcf_htslib",
        "dnakit.ngs.phase_variants_whatshap",
        "dnakit.ngs.sample_sequences_seqtk",
        "dnakit.ngs.assemble_reads",
        "dnakit.ngs.run_dna_seq_workflow",
        "dnakit.ngs.run_rna_seq_workflow",
        "dnakit.ngs.call_chip_peaks",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True


def test_map_reads_rejects_existing_output_directory(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fq"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    reads.write_text("@r\nACGT\n+\nIIII\n", encoding="ascii")
    output = tmp_path / "mapped"
    output.mkdir()
    with pytest.raises(ConfigurationError, match="already exists"):
        map_reads(reference, reads, output, verbose=False)


def test_chip_peak_facade_validates_backend_inputs_before_execution(tmp_path: Path) -> None:
    treatment = tmp_path / "treatment.bam"
    control = tmp_path / "control.bam"
    treatment.write_bytes(b"BAM")
    control.write_bytes(b"BAM")
    from dnakit.ngs import call_chip_peaks

    with pytest.raises(ConfigurationError, match="qvalue"):
        call_chip_peaks(
            treatment,
            control,
            tmp_path / "peaks",
            genome_size=1000,
            qvalue=1.0,
            verbose=False,
        )
