"""JSON-safe facades for caller-installed alignment command-line backends.

The wrappers expose a small, fixed operation set for MAFFT/MUSCLE, NCBI BLAST+,
minimap2, BEDTools and BEDOPS. They reuse DNAKit's bounded subprocess runner, never
invoke a shell, and never install or download an upstream executable.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import alignment as _backend
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
                "implementation": "dnakit.alignment.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def align_sequences(
    sequences: Mapping[str, str],
    output_directory: str | Path,
    *,
    method: Literal["mafft", "muscle", "clustalo"] = "mafft",
    program: str | None = None,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Create a multiple alignment with a caller-installed aligner."""

    result = _backend.align_sequences(
        sequences,
        output_directory,
        method=method,
        program=program,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result("align_sequences", method, result, verbose=verbose)


def build_blast_database(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    program: str = "makeblastdb",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Build a nucleotide BLAST database in a fresh output directory."""

    result = _backend.build_blast_database(
        fasta,
        output_directory,
        program=program,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result("build_blast_database", "ncbi-blast+", result, verbose=verbose)


def local_blast(
    query: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    evalue: float = 1e-5,
    program: str = "blastn",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Search a local nucleotide BLAST database and write BLAST XML."""

    result = _backend.local_blast(
        query,
        database,
        output_directory,
        evalue=evalue,
        program=program,
        config=_config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes),
        verbose=verbose,
    )
    return _result("local_blast", "ncbi-blast+", result, verbose=verbose)


def map_long_reads(
    reference: str | Path,
    reads: str | Path,
    output_directory: str | Path,
    *,
    preset: Literal["map-ont", "map-pb", "sr", "asm5", "asm10", "asm20"] = "map-ont",
    output_format: Literal["sam", "paf"] = "sam",
    program: str = "minimap2",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Map reads to a reference with minimap2 and capture SAM or PAF output."""

    reference_path, reads_path = input_path(reference), input_path(reads)
    directory = workspace(output_directory)
    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    target = directory / ("alignments.sam" if output_format == "sam" else "alignments.paf")
    arguments = []
    if output_format == "sam":
        arguments.append("-a")
    arguments.extend(
        ("-x", preset, "-t", str(config.threads), str(reference_path), str(reads_path))
    )
    log = run_to_file(program, arguments, target, config)
    return _result(
        "map_long_reads",
        "minimap2",
        {
            "backend": "minimap2",
            "directory": directory,
            "artifacts": artifacts(directory),
            "log": log,
            "preset": preset,
            "output_format": output_format,
        },
        verbose=verbose,
    )


def intersect_bed(
    first: str | Path,
    second: str | Path,
    output_directory: str | Path,
    *,
    program: str = "bedtools",
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Write the BEDTools ``intersect -wa -wb`` result for two local BED files."""

    first_path, second_path = input_path(first), input_path(second)
    directory = workspace(output_directory)
    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    target = directory / "intersection.bed"
    log = run_to_file(
        program,
        ("intersect", "-a", str(first_path), "-b", str(second_path), "-wa", "-wb"),
        target,
        config,
    )
    return _result(
        "intersect_bed",
        "bedtools",
        {
            "backend": "bedtools",
            "directory": directory,
            "artifacts": artifacts(directory),
            "log": log,
            "operation": "intersect -wa -wb",
        },
        verbose=verbose,
    )


def sort_bed_bedops(
    source: str | Path,
    output_directory: str | Path,
    *,
    program: str = "sort-bed",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Sort a BED file with caller-installed BEDOPS ``sort-bed``."""

    source_path = input_path(source)
    directory = workspace(output_directory)
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    target = directory / "sorted.bed"
    log = run_to_file(program, (str(source_path),), target, config)
    return _result(
        "sort_bed_bedops",
        "bedops",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "log": log,
        },
        verbose=verbose,
    )


def search_diamond_blastp(
    query: str | Path,
    reference_proteins: str | Path,
    output_directory: str | Path,
    *,
    evalue: float = 1e-5,
    max_target_sequences: int = 25,
    threads: int = 1,
    program: str = "diamond",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Search protein queries against caller-provided proteins with DIAMOND.

    The facade builds a local ``.dmnd`` database in the fresh output directory
    and writes BLAST tabular output. DIAMOND and the biological interpretation
    of hits remain caller responsibilities.
    """

    return _diamond_search(
        query,
        reference_proteins,
        output_directory,
        translated=False,
        evalue=evalue,
        max_target_sequences=max_target_sequences,
        threads=threads,
        program=program,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
        verbose=verbose,
    )


def search_diamond_blastx(
    query: str | Path,
    reference_proteins: str | Path,
    output_directory: str | Path,
    *,
    evalue: float = 1e-5,
    max_target_sequences: int = 25,
    threads: int = 1,
    program: str = "diamond",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Search translated nucleotide queries against proteins with DIAMOND."""

    return _diamond_search(
        query,
        reference_proteins,
        output_directory,
        translated=True,
        evalue=evalue,
        max_target_sequences=max_target_sequences,
        threads=threads,
        program=program,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
        verbose=verbose,
    )


def _diamond_search(
    query: str | Path,
    reference_proteins: str | Path,
    output_directory: str | Path,
    *,
    translated: bool,
    evalue: float,
    max_target_sequences: int,
    threads: int,
    program: str,
    timeout_seconds: float,
    max_output_bytes: int,
    max_artifact_bytes: int,
    verbose: bool,
) -> dict[str, object]:
    query_path = input_path(query)
    reference_path = input_path(reference_proteins)
    if isinstance(evalue, bool) or not isinstance(evalue, (int, float)) or not 1e-100 <= evalue <= 1000:
        raise ConfigurationError("evalue must be a number in [1e-100, 1000].")
    if isinstance(max_target_sequences, bool) or not isinstance(max_target_sequences, int) or not 1 <= max_target_sequences <= 10_000:
        raise ConfigurationError("max_target_sequences must be an integer in [1, 10000].")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    directory = workspace(output_directory)
    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    database_prefix = directory / "reference"
    database = directory / "reference.dmnd"
    hits = directory / "hits.tsv"
    run(
        program,
        ("makedb", "--in", str(reference_path), "--db", str(database_prefix), "--threads", str(threads)),
        directory,
        config,
        outputs=(database,),
    )
    mode = "blastx" if translated else "blastp"
    run(
        program,
        (
            mode,
            "--query", str(query_path),
            "--db", str(database),
            "--out", str(hits),
            "--outfmt", "6",
            "--evalue", str(evalue),
            "--max-target-seqs", str(max_target_sequences),
            "--threads", str(threads),
        ),
        directory,
        config,
        outputs=(hits,),
    )
    return _result(
        f"search_diamond_{mode}",
        "diamond",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "query": query_path,
            "reference_proteins": reference_path,
            "database": database,
            "hits": hits,
            "evalue": evalue,
            "max_target_sequences": max_target_sequences,
            "threads": threads,
        },
        verbose=verbose,
    )


def search_mmseqs_protein(
    query: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    evalue: float = 1e-5,
    memory_megabytes: int = 4096,
    threads: int = 1,
    program: str = "mmseqs",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Search protein FASTA files with MMseqs2 ``easy-search``."""

    query_path, reference_path = input_path(query), input_path(reference)
    _validate_mmseqs_options(evalue, memory_megabytes, threads)
    directory = workspace(output_directory)
    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    hits = directory / "hits.tsv"
    arguments = (
        "easy-search", str(query_path), str(reference_path), str(hits), str(directory / "tmp"),
        "--search-type", "1", "--threads", str(threads), "-e", str(evalue),
        "--split-memory-limit", f"{memory_megabytes}M", "--remove-tmp-files", "1",
        "--createdb-mode", "0",
    )
    run(program, arguments, directory, config, outputs=(hits,))
    return _result(
        "search_mmseqs_protein",
        "mmseqs2",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "query": query_path,
            "reference": reference_path,
            "hits": hits,
            "evalue": evalue,
            "memory_megabytes": memory_megabytes,
            "threads": threads,
        },
        verbose=verbose,
    )


def cluster_mmseqs_protein(
    proteins: str | Path,
    output_directory: str | Path,
    *,
    identity: float = 0.5,
    coverage: float = 0.8,
    memory_megabytes: int = 4096,
    threads: int = 1,
    program: str = "mmseqs",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Cluster proteins with MMseqs2 ``easy-cluster``."""

    proteins_path = input_path(proteins)
    _validate_mmseqs_options(1e-5, memory_megabytes, threads)
    if isinstance(identity, bool) or not isinstance(identity, (int, float)) or not 0.1 <= identity <= 1:
        raise ConfigurationError("identity must be a number in [0.1, 1].")
    if isinstance(coverage, bool) or not isinstance(coverage, (int, float)) or not 0.01 <= coverage <= 1:
        raise ConfigurationError("coverage must be a number in [0.01, 1].")
    directory = workspace(output_directory)
    config = _config(threads, timeout_seconds, max_output_bytes, max_artifact_bytes)
    prefix = directory / "clusters"
    outputs = (
        Path(str(prefix) + "_cluster.tsv"),
        Path(str(prefix) + "_rep_seq.fasta"),
        Path(str(prefix) + "_all_seqs.fasta"),
    )
    arguments = (
        "easy-cluster", str(proteins_path), str(prefix), str(directory / "tmp"),
        "--min-seq-id", str(identity), "-c", str(coverage), "--cov-mode", "1",
        "--threads", str(threads), "--split-memory-limit", f"{memory_megabytes}M",
        "--remove-tmp-files", "1", "--createdb-mode", "0",
    )
    run(program, arguments, directory, config, outputs=outputs)
    return _result(
        "cluster_mmseqs_protein",
        "mmseqs2",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "proteins": proteins_path,
            "cluster_table": outputs[0],
            "representatives": outputs[1],
            "all_sequences": outputs[2],
            "identity": identity,
            "coverage": coverage,
            "memory_megabytes": memory_megabytes,
            "threads": threads,
        },
        verbose=verbose,
    )


def _validate_mmseqs_options(evalue: float, memory_megabytes: int, threads: int) -> None:
    if isinstance(evalue, bool) or not isinstance(evalue, (int, float)) or not 1e-100 <= evalue <= 1000:
        raise ConfigurationError("evalue must be a number in [1e-100, 1000].")
    if isinstance(memory_megabytes, bool) or not isinstance(memory_megabytes, int) or not 128 <= memory_megabytes <= 32_768:
        raise ConfigurationError("memory_megabytes must be an integer in [128, 32768].")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
def consensus_spoa(
    source: str | Path,
    output_directory: str | Path,
    *,
    algorithm: Literal[0, 1, 2] = 2,
    strand_ambiguous: bool = False,
    program: str = "spoa",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Generate a consensus FASTA with caller-installed SPOA."""

    source_path = input_path(source)
    if isinstance(algorithm, bool) or not isinstance(algorithm, int) or algorithm not in {0, 1, 2}:
        raise ConfigurationError("algorithm must be one of 0 (local), 1 (global), or 2 (semi-global).")
    if not isinstance(strand_ambiguous, bool):
        raise ConfigurationError("strand_ambiguous must be a boolean.")
    directory = workspace(output_directory)
    config = _config(1, timeout_seconds, max_output_bytes, max_artifact_bytes)
    target = directory / "consensus.fasta"
    arguments = ["-l", str(algorithm), "-r", "0"]
    if strand_ambiguous:
        arguments.append("-s")
    arguments.append(str(source_path))
    log = run_to_file(program, arguments, target, config)
    return _result(
        "consensus_spoa",
        "spoa",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "algorithm": algorithm,
            "strand_ambiguous": strand_ambiguous,
            "log": log,
        },
        verbose=verbose,
    )


__all__ = [
    "align_sequences",
    "build_blast_database",
    "cluster_mmseqs_protein",
    "consensus_spoa",
    "intersect_bed",
    "local_blast",
    "map_long_reads",
    "search_diamond_blastp",
    "search_diamond_blastx",
    "search_mmseqs_protein",
    "sort_bed_bedops",
]

_register_public_verbose(globals(), tuple(__all__))
