"""Python orchestration of specialist mapping, assembly and NGS algorithms.

Programs are separately installed. Each operation writes into a new directory;
input data are not modified. Commands never use a shell or implicit downloads.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path

from dnakit.bioinformatics._common import (
    ExecutionConfig,
    ToolResult,
    artifacts,
    input_path,
    positive,
    run,
    run_to_file,
    workspace,
)
from dnakit.exceptions import ConfigurationError


def map_reads(
    reference: str | Path,
    first: str | Path,
    output_directory: str | Path,
    *,
    second: str | Path | None = None,
    method: str = "bwa",
    program: str | None = None,
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """Build private BWA/Bowtie2 index and map single/paired reads to SAM."""
    if method not in {"bwa", "bowtie2"}:
        raise ConfigurationError("method must be bwa or bowtie2.")
    fasta, one = input_path(reference), input_path(first)
    two = input_path(second) if second is not None else None
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    prefix, sam = directory / "reference", directory / "aligned.sam"
    if method == "bwa":
        log = run(program or "bwa", ("index", "-p", str(prefix), str(fasta)), directory, config)
        args = ["mem", "-t", str(config.threads), str(prefix), str(one)]
        if two is not None:
            args.append(str(two))
        log += run_to_file(program or "bwa", args, sam, config)
    else:
        log = run(
            "bowtie2-build",
            ("--threads", str(config.threads), str(fasta), str(prefix)),
            directory,
            config,
        )
        args = ["-x", str(prefix), "-p", str(config.threads), "-S", str(sam)]
        args.extend(("-U", str(one)) if two is None else ("-1", str(one), "-2", str(two)))
        log += run(program or "bowtie2", args, directory, config, outputs=(sam,))
    return ToolResult(method, directory, artifacts(directory), log)


def prepare_bam(
    source: str | Path,
    output_directory: str | Path,
    *,
    mark_duplicates: bool = True,
    program: str = "samtools",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """Name sort -> fixmate -m -> coordinate sort -> optional markdup -> index -> stats.

    The duplicate-marking path is intended for paired-end alignments with mate
    information. Set mark_duplicates=False for simple coordinate sort/index.
    """
    source = input_path(source)
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    target = directory / "aligned.bam"
    logs = []
    sort_input = source
    if mark_duplicates:
        names, fixed = directory / "names.bam", directory / "fixmate.bam"
        logs.append(
            run(
                program,
                ("sort", "-n", "-@", str(config.threads), "-o", str(names), str(source)),
                directory,
                config,
                outputs=(names,),
            )
        )
        logs.append(
            run(
                program,
                ("fixmate", "-m", str(names), str(fixed)),
                directory,
                config,
                outputs=(fixed,),
            )
        )
        sort_input = fixed
    sorted_path = directory / "sorted.bam" if mark_duplicates else target
    logs.append(
        run(
            program,
            ("sort", "-@", str(config.threads), "-o", str(sorted_path), str(sort_input)),
            directory,
            config,
            outputs=(sorted_path,),
        )
    )
    if mark_duplicates:
        logs.append(
            run(
                program,
                ("markdup", str(sorted_path), str(target)),
                directory,
                config,
                outputs=(target,),
            )
        )
    logs.append(run(program, ("index", str(target)), directory, config))
    logs.append(run_to_file(program, ("stats", str(target)), directory / "stats.txt", config))
    logs.append(run_to_file(program, ("flagstat", str(target)), directory / "flagstat.txt", config))
    return ToolResult("samtools", directory, artifacts(directory), "\n".join(logs))


def call_variants(
    bam: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    ploidy: int = 2,
    program: str = "bcftools",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """bcftools mpileup -> multiallelic call -> reference normalization -> VCF index.

    Requires a coordinate-sorted indexed BAM and an indexed reference FASTA.
    Ploidy is explicit; biological filtering remains an explicit downstream choice.
    """
    source, fasta = input_path(bam), input_path(reference)
    positive(ploidy, "ploidy")
    if not Path(str(fasta) + ".fai").is_file():
        raise ConfigurationError("Reference .fai is required for variant calling.")
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    likelihoods, calls, vcf = (
        directory / "likelihoods.bcf",
        directory / "calls.bcf",
        directory / "variants.vcf.gz",
    )
    logs = [
        run(
            program,
            ("mpileup", "-Ob", "-f", str(fasta), "-o", str(likelihoods), str(source)),
            directory,
            config,
            outputs=(likelihoods,),
        )
    ]
    logs.append(
        run(
            program,
            ("call", "-mv", "--ploidy", str(ploidy), "-Ob", "-o", str(calls), str(likelihoods)),
            directory,
            config,
            outputs=(calls,),
        )
    )
    logs.append(
        run(
            program,
            ("norm", "-f", str(fasta), "-Oz", "-o", str(vcf), str(calls)),
            directory,
            config,
            outputs=(vcf,),
        )
    )
    logs.append(run(program, ("index", "-t", str(vcf)), directory, config))
    return ToolResult("bcftools", directory, artifacts(directory), "\n".join(logs))


def assemble_reads(
    first: str | Path,
    output_directory: str | Path,
    *,
    second: str | Path | None = None,
    nanopore: str | Path | None = None,
    program: str = "spades.py",
    memory_gb: int = 8,
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """SPAdes short-read or hybrid de novo assembly, with an explicit memory budget."""
    one = input_path(first)
    two = input_path(second) if second is not None else None
    long_reads = input_path(nanopore) if nanopore is not None else None
    positive(memory_gb, "memory_gb")
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    args = ["-o", str(directory / "assembly"), "-t", str(config.threads), "-m", str(memory_gb)]
    args.extend(("-s", str(one)) if two is None else ("-1", str(one), "-2", str(two)))
    if long_reads is not None:
        args.extend(("--nanopore", str(long_reads)))
    log = run(program, args, directory, config)
    return ToolResult("spades", directory, artifacts(directory), log)


def quality_control(
    reads: Iterable[str | Path],
    output_directory: str | Path,
    *,
    program: str = "fastqc",
    config: ExecutionConfig | None = None,
) -> ToolResult:
    sources = tuple(str(input_path(path)) for path in reads)
    if not sources:
        raise ConfigurationError("Provide at least one read file.")
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    log = run(
        program,
        ("--threads", str(config.threads), "--outdir", str(directory), *sources),
        directory,
        config,
    )
    return ToolResult("fastqc", directory, artifacts(directory), log)


def dna_seq_workflow(
    reference: str | Path,
    first: str | Path,
    output_directory: str | Path,
    *,
    second: str | Path | None = None,
    ploidy: int = 2,
    qc: bool = True,
    config: ExecutionConfig | None = None,
) -> dict[str, ToolResult]:
    """QC -> BWA mapping -> sorted indexed BAM -> variant calling; retain every stage."""
    directory = workspace(output_directory)
    results = {}
    if qc:
        results["qc"] = quality_control(
            [first] + ([second] if second is not None else []), directory / "qc", config=config
        )
    results["mapping"] = map_reads(
        reference, first, directory / "mapping", second=second, config=config
    )
    results["bam"] = prepare_bam(
        results["mapping"].artifacts["aligned.sam"],
        directory / "bam",
        mark_duplicates=second is not None,
        config=config,
    )
    results["variants"] = call_variants(
        results["bam"].artifacts["aligned.bam"],
        reference,
        directory / "variants",
        ploidy=ploidy,
        config=config,
    )
    return results


def rna_seq_workflow(
    reference: str | Path,
    annotation: str | Path,
    first: str | Path,
    output_directory: str | Path,
    *,
    second: str | Path | None = None,
    strandedness: int = 0,
    config: ExecutionConfig | None = None,
) -> dict[str, ToolResult]:
    """HISAT2 splice-aware mapping -> BAM -> StringTie assembly -> featureCounts.

    strandedness: 0 unstranded, 1 forward, 2 reverse. Differential expression
    requires biological replicate samples and is a separate statistical step.
    """
    if strandedness not in {0, 1, 2}:
        raise ConfigurationError("strandedness must be 0, 1 or 2.")
    fasta, gtf, one = input_path(reference), input_path(annotation), input_path(first)
    two = input_path(second) if second is not None else None
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    prefix, sam = directory / "reference", directory / "rna.sam"
    log = run(
        "hisat2-build", ("-p", str(config.threads), str(fasta), str(prefix)), directory, config
    )
    args = ["-x", str(prefix), "-p", str(config.threads), "--dta", "-S", str(sam)]
    args.extend(("-U", str(one)) if two is None else ("-1", str(one), "-2", str(two)))
    if strandedness:
        args.extend(
            (
                "--rna-strandness",
                ("F" if strandedness == 1 else "R")
                if two is None
                else ("FR" if strandedness == 1 else "RF"),
            )
        )
    log += run("hisat2", args, directory, config, outputs=(sam,))
    bam = prepare_bam(sam, directory / "bam", mark_duplicates=False, config=config)
    bam_path = bam.artifacts["aligned.bam"]
    transcript = directory / "transcripts.gtf"
    stringtie_args = [
        str(bam_path),
        "-G",
        str(gtf),
        "-o",
        str(transcript),
        "-p",
        str(config.threads),
    ]
    if strandedness:
        stringtie_args.append("--fr" if strandedness == 1 else "--rf")
    log += run("stringtie", stringtie_args, directory, config, outputs=(transcript,))
    counts = directory / "counts.tsv"
    count_args = [
        "-a",
        str(gtf),
        "-o",
        str(counts),
        "-T",
        str(config.threads),
        "-s",
        str(strandedness),
    ]
    if two is not None:
        count_args.extend(("-p", "--countReadPairs"))
    log += run("featureCounts", (*count_args, str(bam_path)), directory, config, outputs=(counts,))
    return {
        "bam": bam,
        "rna": ToolResult("hisat2-stringtie-featureCounts", directory, artifacts(directory), log),
    }


def chip_seq_workflow(
    treatment: str | Path,
    control: str | Path,
    output_directory: str | Path,
    *,
    genome_size: int,
    paired_end: bool = False,
    qvalue: float = 0.05,
    config: ExecutionConfig | None = None,
) -> ToolResult:
    """MACS3 peak calling from prepared treatment/control BAMs, plus samtools coverage."""
    treatment, control = input_path(treatment), input_path(control)
    positive(genome_size, "genome_size")
    if not math.isfinite(qvalue) or not 0 < qvalue < 1:
        raise ConfigurationError("qvalue must be in (0,1).")
    config = config or ExecutionConfig()
    directory = workspace(output_directory)
    log = run(
        "macs3",
        (
            "callpeak",
            "-t",
            str(treatment),
            "-c",
            str(control),
            "-f",
            "BAMPE" if paired_end else "BAM",
            "-g",
            str(genome_size),
            "-q",
            str(qvalue),
            "-n",
            "peaks",
            "--outdir",
            str(directory),
        ),
        directory,
        config,
    )
    log += run_to_file("samtools", ("coverage", str(treatment)), directory / "coverage.tsv", config)
    return ToolResult("macs3", directory, artifacts(directory), log)
