"""Bounded facades for caller-installed microbial classification tools."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, input_path, run, workspace
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def _program(value: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(char in value for char in "\x00\r\n"):
        raise ConfigurationError("program must be a non-empty executable name without newlines.")
    return value.strip()


def _threads(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    return value


def _directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {path}")
    return path


def classify_kraken2(
    reads: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    confidence: float = 0.0,
    threads: int = 1,
    program: str = "kraken2",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Classify FASTA/FASTQ reads against a caller-prepared Kraken2 database."""

    reads_path = input_path(reads)
    database_path = Path(database).expanduser().resolve()
    if not database_path.is_dir():
        raise ConfigurationError(f"database directory does not exist: {database_path}")
    if isinstance(confidence, bool) or not isinstance(confidence, (int, float)) or not 0 <= confidence <= 1:
        raise ConfigurationError("confidence must be a number in [0, 1].")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    if not isinstance(program, str) or not program.strip() or any(char in program for char in "\x00\r\n"):
        raise ConfigurationError("program must be a non-empty executable name without newlines.")
    directory = workspace(output_directory)
    report = directory / "report.tsv"
    assignments = directory / "assignments.tsv"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program.strip(),
        (
            "--db", str(database_path), "--confidence", str(confidence), "--threads", str(threads),
            "--report", str(report), "--output", str(assignments), str(reads_path),
        ),
        directory,
        config,
        outputs=(report, assignments),
    )
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "classify_kraken2",
            "backend": "kraken2",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "database": database_path,
                "report": report,
                "assignments": assignments,
                "confidence": confidence,
                "threads": threads,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def estimate_bracken_abundance(
    report: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    read_length: int = 100,
    level: str = "S",
    threshold: int = 10,
    program: str = "bracken",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Estimate taxon abundances from a Kraken2 report with Bracken."""

    report_path = input_path(report)
    database_path = _directory(database, "database")
    if isinstance(read_length, bool) or not isinstance(read_length, int) or not 1 <= read_length <= 10_000:
        raise ConfigurationError("read_length must be an integer in [1, 10000].")
    if level not in {"D", "P", "C", "O", "F", "G", "S"}:
        raise ConfigurationError("level must be one of D, P, C, O, F, G or S.")
    if isinstance(threshold, bool) or not isinstance(threshold, int) or not 0 <= threshold <= 1_000_000:
        raise ConfigurationError("threshold must be an integer in [0, 1000000].")
    program = _program(program)
    directory = workspace(output_directory)
    abundance = directory / "bracken.tsv"
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program,
        (
            "-d", str(database_path), "-i", str(report_path), "-o", str(abundance),
            "-r", str(read_length), "-l", level, "-t", str(threshold),
        ),
        directory,
        config,
        outputs=(abundance,),
    )
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "estimate_bracken_abundance",
            "backend": "bracken",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "report": report_path,
                "database": database_path,
                "abundance": abundance,
                "read_length": read_length,
                "level": level,
                "threshold": threshold,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def classify_centrifuge(
    index_prefix: str | Path,
    reads: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    program: str = "centrifuge",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Classify reads with a caller-prepared Centrifuge index."""

    reads_path = input_path(reads)
    prefix = Path(index_prefix).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"index parent directory does not exist: {prefix.parent}")
    if not prefix.name or any(char in prefix.name for char in "\x00\r\n"):
        raise ConfigurationError("index_prefix must be a valid index stem.")
    threads = _threads(threads)
    program = _program(program)
    directory = workspace(output_directory)
    assignments = directory / "assignments.tsv"
    report = directory / "report.tsv"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program,
        (
            "-x", str(prefix), "-U", str(reads_path), "-S", str(assignments),
            "--report-file", str(report), "-p", str(threads),
        ),
        directory,
        config,
        outputs=(assignments, report),
    )
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "classify_centrifuge",
            "backend": "centrifuge",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "index_prefix": prefix,
                "reads": reads_path,
                "assignments": assignments,
                "report": report,
                "threads": threads,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def classify_kaiju(
    reads: str | Path,
    protein_database: str | Path,
    taxonomy_nodes: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    program: str = "kaiju",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Classify reads with a caller-prepared Kaiju protein database and taxonomy."""

    reads_path = input_path(reads)
    database_path = input_path(protein_database)
    taxonomy_path = input_path(taxonomy_nodes)
    threads = _threads(threads)
    program = _program(program)
    directory = workspace(output_directory)
    assignments = directory / "assignments.tsv"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program,
        (
            "-t", str(taxonomy_path), "-f", str(database_path), "-i", str(reads_path),
            "-o", str(assignments), "-z", str(threads),
        ),
        directory,
        config,
        outputs=(assignments,),
    )
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "classify_kaiju",
            "backend": "kaiju",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "protein_database": database_path,
                "taxonomy_nodes": taxonomy_path,
                "assignments": assignments,
                "threads": threads,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def classify_gtdbtk(
    genomes: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    extension: str = "fna",
    threads: int = 1,
    program: str = "gtdbtk",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Classify bacterial/archaeal genome bins with caller-installed GTDB-Tk.

    GTDB-Tk reads its reference package from ``GTDBTK_DATA_PATH``.  DNAKit
    sets that variable for this bounded process from ``database`` and invokes
    the documented ``classify_wf`` command without downloading or modifying
    the reference package.
    """

    genomes_path = _directory(genomes, "genomes")
    database_path = _directory(database, "database")
    if not isinstance(extension, str) or not extension or not all(char.isalnum() or char == "." for char in extension):
        raise ConfigurationError("extension must contain only letters, digits or periods.")
    threads = _threads(threads)
    program = _program(program)
    directory = workspace(output_directory)
    result_directory = directory / "gtdbtk"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program,
        (
            "classify_wf", "--genome_dir", str(genomes_path), "--out_dir", str(result_directory),
            "--cpus", str(threads), "--extension", extension,
        ),
        directory,
        config,
        environment={"GTDBTK_DATA_PATH": str(database_path)},
    )
    summaries = tuple(result_directory.rglob("*.summary.tsv"))
    if not summaries:
        raise BackendExecutionError("GTDB-Tk did not produce a *.summary.tsv classification report.")
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "classify_gtdbtk",
            "backend": "gtdbtk",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "genomes": genomes_path,
                "database": database_path,
                "classification_directory": result_directory,
                "summaries": summaries,
                "extension": extension,
                "threads": threads,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def profile_humann(
    reads: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    program: str = "humann",
    nucleotide_database: str | Path | None = None,
    protein_database: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Profile metagenomic reads with caller-installed HUMAnN 3."""

    reads_path = input_path(reads)
    threads = _threads(threads)
    program = _program(program)
    nucleotide_path = _directory(nucleotide_database, "nucleotide_database") if nucleotide_database is not None else None
    protein_path = _directory(protein_database, "protein_database") if protein_database is not None else None
    directory = workspace(output_directory)
    result_directory = directory / "humann"
    result_directory.mkdir()
    stem = reads_path.name.split(".", 1)[0]
    outputs = tuple(result_directory / f"{stem}{suffix}" for suffix in ("_genefamilies.tsv", "_pathabundance.tsv", "_pathcoverage.tsv"))
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    arguments = ["--input", str(reads_path), "--output", str(result_directory), "--threads", str(threads)]
    if nucleotide_path is not None:
        arguments.extend(("--nucleotide-database", str(nucleotide_path)))
    if protein_path is not None:
        arguments.extend(("--protein-database", str(protein_path)))
    run(program, arguments, directory, config, outputs=outputs)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "profile_humann",
            "backend": "humann",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "output_directory": result_directory,
                "gene_families": outputs[0],
                "pathway_abundance": outputs[1],
                "pathway_coverage": outputs[2],
                "threads": threads,
            },
            "provenance": {
                "implementation": "dnakit.microbiome.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )
__all__ = [
    "classify_centrifuge",
    "classify_gtdbtk",
    "classify_kaiju",
    "classify_kraken2",
    "estimate_bracken_abundance",
    "profile_humann",
]

_register_public_verbose(globals(), tuple(__all__))
