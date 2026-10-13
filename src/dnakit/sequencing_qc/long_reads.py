"""JSON-safe facades for caller-installed long-read QC programs.

The wrappers expose a small, stable option surface for NanoFilt, Filtlong and
NanoPlot.  They do not install backends, download reference data or invoke a
shell; each backend executable must be supplied by the caller or available on
``PATH``.
"""

from __future__ import annotations

import math
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
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


def filter_nanopore_reads(
    source: str | Path,
    output_directory: str | Path,
    *,
    min_quality: float | None = None,
    min_length: int | None = None,
    headcrop: int | None = None,
    tailcrop: int | None = None,
    program: str = "NanoFilt",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Filter and trim long-read FASTQ records with caller-installed NanoFilt.

    NanoFilt reads an uncompressed FASTQ path and writes filtered FASTQ to
    stdout.  The adapter captures that stream as ``filtered.fastq`` and keeps
    stderr as a diagnostic artifact.
    """

    source_path = input_path(source)
    _validate_program(program, "program")
    _validate_number(min_quality, "min_quality", minimum=0.0, maximum=100.0)
    _validate_integer(min_length, "min_length", minimum=0, maximum=10_000_000_000)
    _validate_integer(headcrop, "headcrop", minimum=0, maximum=10_000_000_000)
    _validate_integer(tailcrop, "tailcrop", minimum=0, maximum=10_000_000_000)
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / "filtered.fastq"
    arguments: list[str] = []
    if min_quality is not None:
        arguments.extend(("-q", str(min_quality)))
    if min_length is not None:
        arguments.extend(("-l", str(min_length)))
    if headcrop is not None:
        arguments.extend(("--headcrop", str(headcrop)))
    if tailcrop is not None:
        arguments.extend(("--tailcrop", str(tailcrop)))
    arguments.append(str(source_path))
    log = run_to_file(program, arguments, target, config)
    return _envelope(
        "filter_nanopore_reads",
        "nanofilt",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "log": log,
        },
        verbose=verbose,
    )


def filter_long_reads(
    source: str | Path,
    output_directory: str | Path,
    *,
    min_length: int | None = None,
    keep_percent: float | None = None,
    target_bases: int | None = None,
    program: str = "filtlong",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Filter long reads with caller-installed Filtlong.

    The wrapper captures Filtlong's stdout in ``filtered.fastq`` and exposes
    only numeric thresholds.  Unit-suffixed and reference-driven Filtlong
    modes remain available through the upstream CLI rather than this bounded
    facade.
    """

    source_path = input_path(source)
    _validate_program(program, "program")
    _validate_integer(min_length, "min_length", minimum=0, maximum=10_000_000_000)
    _validate_number(keep_percent, "keep_percent", minimum=0.0, maximum=100.0)
    _validate_integer(target_bases, "target_bases", minimum=0, maximum=10_000_000_000_000)
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / "filtered.fastq"
    arguments: list[str] = []
    if min_length is not None:
        arguments.extend(("--min_length", str(min_length)))
    if keep_percent is not None:
        arguments.extend(("--keep_percent", str(keep_percent)))
    if target_bases is not None:
        arguments.extend(("--target_bases", str(target_bases)))
    arguments.append(str(source_path))
    log = run_to_file(program, arguments, target, config)
    return _envelope(
        "filter_long_reads",
        "filtlong",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "log": log,
        },
        verbose=verbose,
    )


def summarize_long_reads(
    source: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    min_length: int | None = None,
    min_quality: float | None = None,
    program: str = "NanoPlot",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Create bounded long-read statistics with caller-installed NanoPlot.

    NanoPlot is run in FASTQ mode with static plots disabled.  The expected
    ``NanoStats.txt`` and ``NanoPlot-report.html`` artifacts are checked, while
    the returned envelope retains every file generated in the fresh directory.
    """

    source_path = input_path(source)
    _validate_program(program, "program")
    _validate_integer(threads, "threads", minimum=1, maximum=256)
    _validate_integer(min_length, "min_length", minimum=0, maximum=10_000_000_000)
    _validate_number(min_quality, "min_quality", minimum=0.0, maximum=100.0)
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    stats = directory / "NanoStats.txt"
    report = directory / "NanoPlot-report.html"
    arguments = [
        "--fastq",
        str(source_path),
        "--outdir",
        str(directory),
        "--threads",
        str(threads),
        "--no_static",
    ]
    if min_length is not None:
        arguments.extend(("--minlength", str(min_length)))
    if min_quality is not None:
        arguments.extend(("--minqual", str(min_quality)))
    log = run(program, arguments, directory, config, outputs=(stats, report))
    _check_artifact_budget(directory, config.max_artifact_bytes)
    return _envelope(
        "summarize_long_reads",
        "nanoplot",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "stats": stats,
            "report": report,
            "log": log,
        },
        verbose=verbose,
    )


def trim_porechop(
    source: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    no_split: bool = False,
    program: str = "porechop",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Trim Oxford Nanopore adapters with caller-installed Porechop.

    Porechop is an archived upstream project.  This boundary keeps only input,
    output, thread and no-split controls and records that the caller supplied
    the executable.
    """

    source_path = input_path(source)
    _validate_program(program, "program")
    _validate_integer(threads, "threads", minimum=1, maximum=256)
    if not isinstance(no_split, bool):
        raise ConfigurationError("no_split must be a boolean.")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / "trimmed.fastq"
    arguments = ["-i", str(source_path), "-o", str(target), "-t", str(threads)]
    if no_split:
        arguments.append("--no_split")
    log = run(program, arguments, directory, config, outputs=(target,))
    return _envelope(
        "trim_porechop",
        "porechop",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": target,
            "threads": threads,
            "no_split": no_split,
            "log": log,
        },
        verbose=verbose,
    )


def _envelope(
    operation: str,
    backend: str,
    result: dict[str, object],
    *,
    verbose: bool,
) -> dict[str, object]:
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": operation,
            "backend": backend,
            "result": result,
            "provenance": {
                "implementation": "dnakit.sequencing_qc.long_reads",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def _validate_program(value: str, name: str) -> None:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > 4096
        or any(char in value for char in "\x00\r\n")
    ):
        raise ConfigurationError(f"{name} must be a bounded executable name or path.")


def _validate_integer(value: int | None, name: str, *, minimum: int, maximum: int) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")


def _validate_number(
    value: float | None,
    name: str,
    *,
    minimum: float,
    maximum: float,
) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    if not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")


def _check_artifact_budget(directory: Path, maximum: int) -> None:
    total = sum(path.stat().st_size for path in directory.rglob("*") if path.is_file())
    if total > maximum:
        raise ConfigurationError(
            f"Generated NanoPlot artifacts exceed max_artifact_bytes ({total} > {maximum})."
        )


__all__ = [
    "filter_long_reads",
    "filter_nanopore_reads",
    "summarize_long_reads",
    "trim_porechop",
]

_register_public_verbose(globals(), tuple(__all__))
