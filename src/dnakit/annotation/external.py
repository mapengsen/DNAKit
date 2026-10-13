"""Bounded facades for caller-installed genome annotation backends."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, input_path, run, workspace
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def _program(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(char in value for char in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty executable name without newlines.")
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


def _simple_name(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(
        char in value for char in ("/", "\\", "\x00", "\n", "\r")
    ):
        raise ConfigurationError(f"{name} must be a simple output name.")
    return value.strip()


def _result(operation: str, backend: str, value: dict[str, object], *, verbose: bool) -> dict[str, object]:
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": operation,
            "backend": backend,
            "result": value,
            "provenance": {
                "implementation": "dnakit.annotation.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def annotate_prokka(
    assembly: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    prefix: str = "annotation",
    program: str = "prokka",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Annotate a prokaryotic assembly with caller-installed Prokka."""

    assembly_path = input_path(assembly)
    program = _program(program, "program")
    threads = _threads(threads)
    prefix = _simple_name(prefix, "prefix")
    directory = workspace(output_directory)
    annotation = directory / "annotation"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    gff = annotation / f"{prefix}.gff"
    faa = annotation / f"{prefix}.faa"
    run(
        program,
        (
            "--outdir", str(annotation), "--prefix", prefix,
            "--cpus", str(threads), str(assembly_path),
        ),
        directory,
        config,
        outputs=(gff, faa),
    )
    return _result(
        "annotate_prokka",
        "prokka",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "assembly": assembly_path,
            "annotation_directory": annotation,
            "gff": gff,
            "proteins": faa,
            "threads": threads,
            "prefix": prefix,
        },
        verbose=verbose,
    )


def annotate_bakta(
    genome: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    min_contig_length: int = 200,
    translation_table: int = 11,
    threads: int = 1,
    prefix: str = "annotation",
    program: str = "bakta",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Annotate a bacterial genome with a caller-prepared Bakta database."""

    genome_path = input_path(genome)
    database_path = _directory(database, "database")
    program = _program(program, "program")
    threads = _threads(threads)
    if isinstance(min_contig_length, bool) or not isinstance(min_contig_length, int) or not 1 <= min_contig_length <= 1_000_000:
        raise ConfigurationError("min_contig_length must be an integer in [1, 1000000].")
    if translation_table not in {4, 11, 25}:
        raise ConfigurationError("translation_table must be one of 4, 11 or 25.")
    prefix = _simple_name(prefix, "prefix")
    directory = workspace(output_directory)
    annotation = directory / "bakta"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    outputs = tuple(annotation / f"{prefix}{suffix}" for suffix in (".gff3", ".tsv", ".json", ".fna"))
    run(
        program,
        (
            "--db", str(database_path), "--output", str(annotation), "--prefix", prefix,
            "--threads", str(threads), "--min-contig-length", str(min_contig_length),
            "--translation-table", str(translation_table), "--skip-plot", str(genome_path),
        ),
        directory,
        config,
        outputs=outputs,
    )
    return _result(
        "annotate_bakta",
        "bakta",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "genome": genome_path,
            "database": database_path,
            "annotation_directory": annotation,
            "gff3": outputs[0],
            "table": outputs[1],
            "json": outputs[2],
            "translated_sequences": outputs[3],
            "min_contig_length": min_contig_length,
            "translation_table": translation_table,
            "threads": threads,
            "prefix": prefix,
        },
        verbose=verbose,
    )


def predict_checkm2(
    bins: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    extension: str = "fna",
    threads: int = 1,
    program: str = "checkm2",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Estimate bin completeness and contamination with caller-installed CheckM2."""

    bins_path = _directory(bins, "bins")
    database_path = input_path(database)
    program = _program(program, "program")
    threads = _threads(threads)
    if not isinstance(extension, str) or not extension or not all(char.isalnum() or char == "." for char in extension):
        raise ConfigurationError("extension must contain only letters, digits or periods.")
    directory = workspace(output_directory)
    result_directory = directory / "checkm2"
    report = result_directory / "quality_report.tsv"
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    run(
        program,
        (
            "predict", "--input", str(bins_path), "--output-directory", str(result_directory),
            "--database_path", str(database_path), "--extension", extension,
            "--threads", str(threads), "--lowmem", "--remove_intermediates",
        ),
        directory,
        config,
        outputs=(report,),
    )
    return _result(
        "predict_checkm2",
        "checkm2",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "bins": bins_path,
            "database": database_path,
            "report": report,
            "extension": extension,
            "threads": threads,
        },
        verbose=verbose,
    )


__all__ = ["annotate_bakta", "annotate_prokka", "predict_checkm2"]

_register_public_verbose(globals(), tuple(__all__))
