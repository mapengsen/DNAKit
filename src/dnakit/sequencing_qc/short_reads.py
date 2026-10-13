"""Bounded facades for caller-installed short-read QC programs."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, input_path, run, workspace
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def trim_trimmomatic(
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    adapter_fasta: str | Path | None = None,
    leading: int | None = None,
    trailing: int | None = None,
    min_length: int | None = None,
    threads: int = 1,
    program: str = "trimmomatic",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Trim Illumina reads with caller-installed Trimmomatic.

    The wrapper supports a small fixed subset of Trimmomatic's steps.  Paired
    input retains both paired and unpaired output files; single-end input
    produces one trimmed FASTQ.
    """

    read1_path = input_path(read1)
    read2_path = input_path(read2) if read2 is not None else None
    if read2_path is not None and read2_path == read1_path:
        raise ConfigurationError("read1 and read2 must be different files.")
    _validate_integer(threads, "threads", minimum=1, maximum=256)
    _validate_integer(leading, "leading", minimum=0, maximum=93)
    _validate_integer(trailing, "trailing", minimum=0, maximum=93)
    _validate_integer(min_length, "min_length", minimum=0, maximum=10_000_000)
    adapter_path = input_path(adapter_fasta) if adapter_fasta is not None else None
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )


    directory = workspace(output_directory)
    arguments: list[str] = []
    outputs: list[Path]
    if read2_path is None:
        target = directory / "trimmed.fastq.gz"
        arguments.extend(("SE", "-threads", str(threads), "-phred33", str(read1_path), str(target)))
        outputs = [target]
    else:
        paired1 = directory / "trimmed_R1.paired.fastq.gz"
        unpaired1 = directory / "trimmed_R1.unpaired.fastq.gz"
        paired2 = directory / "trimmed_R2.paired.fastq.gz"
        unpaired2 = directory / "trimmed_R2.unpaired.fastq.gz"
        arguments.extend(
            (
                "PE",
                "-threads",
                str(threads),
                "-phred33",
                str(read1_path),
                str(read2_path),
                str(paired1),
                str(unpaired1),
                str(paired2),
                str(unpaired2),
            )
        )
        outputs = [paired1, unpaired1, paired2, unpaired2]
    steps: list[str] = []
    if adapter_path is not None:
        steps.append(f"ILLUMINACLIP:{adapter_path}:2:30:10")
    if leading is not None:
        steps.append(f"LEADING:{leading}")
    if trailing is not None:
        steps.append(f"TRAILING:{trailing}")
    if min_length is not None:
        steps.append(f"MINLEN:{min_length}")
    if not steps:
        raise ConfigurationError("At least one trimming step is required.")
    arguments.extend(steps)
    log = run(program, arguments, directory, config, outputs=outputs)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "trim_trimmomatic",
            "backend": "trimmomatic",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "read1": read1_path,
                "read2": read2_path,
                "outputs": outputs,
                "steps": steps,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.sequencing_qc.short_reads",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def run_afterqc(
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    qc_only: bool = False,
    trim_front: int | None = None,
    trim_tail: int | None = None,
    qualified_quality_phred: int | None = None,
    unqualified_base_limit: int | None = None,
    poly_size_limit: int | None = None,
    n_base_limit: int | None = None,
    seq_len_req: int | None = None,
    gzip_output: bool = False,
    compression: int = 2,
    program: str = "after.py",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Run a bounded AfterQC filtering and reporting job.

    AfterQC writes ``good``, ``bad`` and ``QC`` directories in the working
    directory.  The facade exposes the common read-pair and quality controls;
    callers still provide the upstream ``after.py`` executable.
    """

    read1_path = input_path(read1)
    read2_path = input_path(read2) if read2 is not None else None
    if read2_path is not None and read2_path == read1_path:
        raise ConfigurationError("read1 and read2 must be different files.")
    _validate_program(program, "program")
    _validate_bool(qc_only, "qc_only")
    _validate_bool(gzip_output, "gzip_output")
    _validate_integer(trim_front, "trim_front", minimum=-1, maximum=10_000_000)
    _validate_integer(trim_tail, "trim_tail", minimum=-1, maximum=10_000_000)
    _validate_integer(
        qualified_quality_phred,
        "qualified_quality_phred",
        minimum=0,
        maximum=93,
    )
    _validate_integer(
        unqualified_base_limit,
        "unqualified_base_limit",
        minimum=0,
        maximum=10_000_000,
    )
    _validate_integer(poly_size_limit, "poly_size_limit", minimum=0, maximum=10_000_000)
    _validate_integer(n_base_limit, "n_base_limit", minimum=0, maximum=10_000_000)
    _validate_integer(seq_len_req, "seq_len_req", minimum=0, maximum=10_000_000)
    _validate_integer(compression, "compression", minimum=0, maximum=9)
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    good = directory / "good"
    bad = directory / "bad"
    qc = directory / "QC"
    arguments = ["-1", str(read1_path)]
    if read2_path is not None:
        arguments.extend(("-2", str(read2_path)))
    arguments.extend(("-g", str(good), "-b", str(bad)))
    if qc_only:
        arguments.append("--qc_only")
    for flag, value in (
        ("-f", trim_front),
        ("-t", trim_tail),
        ("-q", qualified_quality_phred),
        ("-u", unqualified_base_limit),
        ("-p", poly_size_limit),
        ("-n", n_base_limit),
        ("-s", seq_len_req),
    ):
        if value is not None:
            arguments.extend((flag, str(value)))
    if gzip_output:
        arguments.append("--gzip")
    arguments.extend(("--compression", str(compression)))
    log = run(program, arguments, directory, config)
    if not qc.is_dir():
        raise BackendExecutionError(f"AfterQC did not produce the expected QC directory: {qc}")
    if not qc_only and (not good.is_dir() or not bad.is_dir()):
        raise BackendExecutionError(
            "AfterQC did not produce the expected good and bad output directories."
        )
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "run_afterqc",
            "backend": "afterqc",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "read1": read1_path,
                "read2": read2_path,
                "good": good,
                "bad": bad,
                "qc": qc,
                "qc_only": qc_only,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.sequencing_qc.short_reads",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def _validate_integer(value: int | None, name: str, *, minimum: int, maximum: int) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")


def _validate_bool(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")


def _validate_program(value: str, name: str) -> None:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > 4096
        or any(char in value for char in "\x00\r\n")
    ):
        raise ConfigurationError(f"{name} must be a bounded executable name or path.")


__all__ = ["run_afterqc", "trim_trimmomatic"]

_register_public_verbose(globals(), tuple(__all__))
