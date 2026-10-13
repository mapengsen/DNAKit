"""JSON-safe facades for caller-installed NGS command-line backends.

The implementation reuses :mod:`dnakit.bioinformatics.ngs`, while this module
keeps the public Agent/MCP contract explicit: paths and scalar options in,
serialized stage metadata and artifact paths out.  BWA, BWA-MEM2, Bowtie2, SAMtools,
BCFtools, SPAdes, HISAT2, StringTie, featureCounts and MACS3 remain separate
executables and are never installed or invoked through a shell by DNAKit.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import ngs as _backend
from dnakit.bioinformatics._common import (
    ExecutionConfig,
    artifacts,
    input_path,
    run,
    run_to_file,
    workspace,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def _config(
    threads: int,
    timeout_seconds: float,
    max_output_bytes: int,
    max_artifact_bytes: int,
) -> ExecutionConfig:
    return ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )


def _result(operation: str, backend: str, value: object, *, verbose: bool) -> dict[str, object]:
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": operation,
            "backend": backend,
            "result": value,
            "provenance": {
                "implementation": "dnakit.bioinformatics.ngs",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def map_reads(
    reference: str | Path,
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    method: Literal["bwa", "bwa-mem2", "bowtie2"] = "bwa",
    program: str | None = None,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Map single- or paired-end reads with BWA, BWA-MEM2 or Bowtie2 into a new directory."""

    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    result = _backend.map_reads(
        reference,
        read1,
        output_directory,
        second=read2,
        method=method,
        program=program,
        config=config,
        verbose=verbose,
    )
    return _result("map_reads", method, result, verbose=verbose)


def prepare_bam(
    source: str | Path,
    output_directory: str | Path,
    *,
    mark_duplicates: bool = True,
    program: str = "samtools",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Sort, optionally mark duplicates, index and summarize a SAM/BAM input."""

    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    result = _backend.prepare_bam(
        source,
        output_directory,
        mark_duplicates=mark_duplicates,
        program=program,
        config=config,
        verbose=verbose,
    )
    return _result("prepare_bam", "samtools", result, verbose=verbose)


def mark_duplicates_picard(
    source: str | Path,
    output_directory: str | Path,
    *,
    remove_duplicates: bool = False,
    assume_sorted: bool = True,
    max_records_in_ram: int | None = None,
    program: str = "picard",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Mark PCR/optical duplicates with caller-installed Picard MarkDuplicates."""

    source_path = input_path(source)
    if not isinstance(remove_duplicates, bool) or not isinstance(assume_sorted, bool):
        raise ConfigurationError("remove_duplicates and assume_sorted must be booleans.")
    if max_records_in_ram is not None and (
        isinstance(max_records_in_ram, bool)
        or not isinstance(max_records_in_ram, int)
        or not 1_000 <= max_records_in_ram <= 100_000_000
    ):
        raise ConfigurationError("max_records_in_ram must be an integer in [1000, 100000000].")
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    directory = workspace(output_directory)
    output_bam = directory / "marked_duplicates.bam"
    metrics = directory / "duplicate_metrics.txt"
    arguments = [
        "MarkDuplicates",
        f"I={source_path}",
        f"O={output_bam}",
        f"M={metrics}",
        f"REMOVE_DUPLICATES={'true' if remove_duplicates else 'false'}",
        f"ASSUME_SORTED={'true' if assume_sorted else 'false'}",
    ]
    if max_records_in_ram is not None:
        arguments.append(f"MAX_RECORDS_IN_RAM={max_records_in_ram}")
    log = run(program, arguments, directory, config, outputs=(output_bam, metrics))
    return _result(
        "mark_duplicates_picard",
        "picard",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": output_bam,
            "metrics": metrics,
            "remove_duplicates": remove_duplicates,
            "assume_sorted": assume_sorted,
            "max_records_in_ram": max_records_in_ram,
            "log": log,
        },
        verbose=verbose,
    )


def call_variants(
    bam: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    ploidy: int = 2,
    program: str = "bcftools",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call and normalize variants from an indexed BAM into a new VCF directory."""

    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    result = _backend.call_variants(
        bam,
        reference,
        output_directory,
        ploidy=ploidy,
        program=program,
        config=config,
        verbose=verbose,
    )
    return _result("call_variants", "bcftools", result, verbose=verbose)


def assemble_reads(
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    nanopore: str | Path | None = None,
    program: str = "spades.py",
    memory_gb: int = 8,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assemble short reads or a hybrid short/long-read set with SPAdes."""

    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    result = _backend.assemble_reads(
        read1,
        output_directory,
        second=read2,
        nanopore=nanopore,
        program=program,
        memory_gb=memory_gb,
        config=config,
        verbose=verbose,
    )
    return _result("assemble_reads", "spades", result, verbose=verbose)


def run_dna_seq_workflow(
    reference: str | Path,
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    ploidy: int = 2,
    qc: bool = True,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Run the bounded QC -> mapping -> BAM -> variant workflow."""

    result = _backend.dna_seq_workflow(
        reference,
        read1,
        output_directory,
        second=read2,
        ploidy=ploidy,
        qc=qc,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result("run_dna_seq_workflow", "bwa-samtools-bcftools", result, verbose=verbose)


def run_rna_seq_workflow(
    reference: str | Path,
    annotation: str | Path,
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    strandedness: Literal[0, 1, 2] = 0,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Run the bounded HISAT2 -> BAM -> StringTie -> featureCounts workflow."""

    result = _backend.rna_seq_workflow(
        reference,
        annotation,
        read1,
        output_directory,
        second=read2,
        strandedness=strandedness,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result(
        "run_rna_seq_workflow", "hisat2-stringtie-featureCounts", result, verbose=verbose
    )


def call_chip_peaks(
    treatment: str | Path,
    control: str | Path,
    output_directory: str | Path,
    *,
    genome_size: int,
    paired_end: bool = False,
    qvalue: float = 0.05,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call ChIP-seq peaks with MACS3 and collect treatment coverage."""

    result = _backend.chip_seq_workflow(
        treatment,
        control,
        output_directory,
        genome_size=genome_size,
        paired_end=paired_end,
        qvalue=qvalue,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result("call_chip_peaks", "macs3-samtools", result, verbose=verbose)


def sample_sequences_seqtk(
    source: str | Path,
    output_directory: str | Path,
    *,
    sample_size: int | float,
    seed: int = 1,
    program: str = "seqtk",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Sample FASTA/FASTQ records with caller-installed seqtk."""

    source_path = input_path(source)
    if isinstance(sample_size, bool) or not isinstance(sample_size, (int, float)):
        raise ConfigurationError("sample_size must be a positive count or fraction.")
    if sample_size <= 0:
        raise ConfigurationError("sample_size must be positive.")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ConfigurationError("seed must be a non-negative integer.")
    if not isinstance(program, str) or not program.strip() or any(
        char in program for char in "\x00\r\n"
    ):
        raise ConfigurationError("program must be a non-empty executable name.")
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    directory = workspace(output_directory)
    target = directory / "sampled.fastx"
    arguments = ["sample", f"-s{seed}", str(source_path), str(sample_size)]
    log = run_to_file(program.strip(), arguments, target, config)
    return _result(
        "sample_sequences_seqtk",
        "seqtk",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "sample_size": sample_size,
            "seed": seed,
            "log": log,
        },
        verbose=verbose,
    )


def compress_index_vcf_htslib(
    source: str | Path,
    output_directory: str | Path,
    *,
    bgzip_program: str = "bgzip",
    tabix_program: str = "tabix",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Compress and index a VCF with caller-installed HTSlib utilities."""

    source_path = input_path(source)
    for value, name in ((bgzip_program, "bgzip_program"), (tabix_program, "tabix_program")):
        if not isinstance(value, str) or not value.strip() or any(char in value for char in "\x00\r\n"):
            raise ConfigurationError(f"{name} must be a non-empty executable name.")
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    directory = workspace(output_directory)
    compressed = directory / "variants.vcf.gz"
    log = run_to_file(bgzip_program.strip(), ("-c", str(source_path)), compressed, config)
    log += run(tabix_program.strip(), ("-p", "vcf", str(compressed)), directory, config)
    return _result(
        "compress_index_vcf_htslib",
        "htslib",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": compressed,
            "index": Path(str(compressed) + ".tbi"),
            "log": log,
        },
        verbose=verbose,
    )


def phase_variants_whatshap(
    variants: str | Path,
    alignment: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    program: str = "whatshap",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Phase a VCF against caller-installed read alignments with WhatsHap."""

    variants_path = input_path(variants)
    alignment_path = input_path(alignment)
    reference_path = input_path(reference)
    if not isinstance(program, str) or not program.strip() or any(
        char in program for char in "\x00\r\n"
    ):
        raise ConfigurationError("program must be a non-empty executable name.")
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    directory = workspace(output_directory)
    target = directory / "phased.vcf"
    arguments = [
        "phase",
        "--reference",
        str(reference_path),
        "--output",
        str(target),
        str(variants_path),
        str(alignment_path),
    ]
    log = run(program.strip(), arguments, directory, config, outputs=(target,))
    return _result(
        "phase_variants_whatshap",
        "whatshap",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "variants": variants_path,
            "alignment": alignment_path,
            "reference": reference_path,
            "output": target,
            "log": log,
        },
        verbose=verbose,
    )


def sort_bam_bamtools(
    source: str | Path,
    output_directory: str | Path,
    *,
    program: str = "bamtools",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Sort a BAM file with caller-installed BamTools."""

    source_path = input_path(source)
    if not isinstance(program, str) or not program.strip() or any(
        char in program for char in "\x00\r\n"
    ):
        raise ConfigurationError("program must be a non-empty executable name.")
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    directory = workspace(output_directory)
    target = directory / "sorted.bam"
    arguments = ["sort", "-in", str(source_path), "-out", str(target)]
    log = run(program.strip(), arguments, directory, config, outputs=(target,))
    return _result(
        "sort_bam_bamtools",
        "bamtools",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "log": log,
        },
        verbose=verbose,
    )


__all__ = [
    "assemble_reads",
    "call_chip_peaks",
    "call_variants",
    "compress_index_vcf_htslib",
    "mark_duplicates_picard",
    "map_reads",
    "prepare_bam",
    "phase_variants_whatshap",
    "run_dna_seq_workflow",
    "run_rna_seq_workflow",
    "sample_sequences_seqtk",
    "sort_bam_bamtools",
]

_register_public_verbose(globals(), tuple(__all__))
