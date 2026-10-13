"""JSON-safe facades for caller-installed genome assembly programs."""

from __future__ import annotations

from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import (
    ExecutionConfig,
    artifacts,
    input_path,
    run,
    workspace,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def _input_directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {path}")
    return path


def _validate_program(program: str) -> str:
    if not isinstance(program, str) or not program.strip() or any(
        character in program for character in ("\x00", "\n", "\r")
    ):
        raise ConfigurationError("program must be a non-empty executable name without newlines.")
    return program.strip()


def assemble_shasta(
    reads: str | Path,
    output_directory: str | Path,
    *,
    config: str = "Nanopore-May2022",
    threads: int | None = None,
    program: str = "shasta",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assemble caller-provided long reads with Shasta.

    The facade exposes the documented input/config/output boundary and returns
    ``Assembly.fasta`` and ``Assembly.gfa``. Shasta, its configuration files,
    and read indexes remain caller-installed dependencies.
    """

    reads_path = input_path(reads)
    if not isinstance(config, str) or not config.strip() or len(config) > 200:
        raise ConfigurationError("config must be a non-empty string of at most 200 characters.")
    if any(character.isspace() for character in config) or any(
        character in config for character in (";", "|", "&", "\x00")
    ):
        raise ConfigurationError("config must be one Shasta configuration name without shell characters.")
    program = _validate_program(program)
    if threads is not None and (
        isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256
    ):
        raise ConfigurationError("threads must be an integer in [1, 256].")
    config_limits = ExecutionConfig(
        threads=threads or 1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    fasta = directory / "Assembly.fasta"
    graph = directory / "Assembly.gfa"
    arguments = [
        "--input",
        str(reads_path),
        "--config",
        config.strip(),
        "--assemblyDirectory",
        str(directory),
    ]
    if threads is not None:
        arguments.extend(("--threads", str(threads)))
    log = run(program, arguments, directory, config_limits, outputs=(fasta, graph))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "assemble_shasta",
            "backend": "shasta",
            "result": {
                "directory": directory,
                "assembly_fasta": fasta,
                "assembly_gfa": graph,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "config": config.strip(),
                "threads": threads,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def evaluate_merqury(
    read_database: str | Path,
    assembly: str | Path,
    output_directory: str | Path,
    *,
    program: str = "merqury.sh",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Evaluate a genome assembly with a caller-installed Merqury workflow.

    ``read_database`` is an existing Meryl database directory. The facade runs
    the documented ``merqury.sh read-db.meryl assembly.fasta prefix`` boundary
    and requires the QV and completeness reports.
    """

    read_database_path = _input_directory(read_database, "read_database")
    assembly_path = input_path(assembly)
    program = _validate_program(program)
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    prefix = directory / "merqury"
    qv = Path(str(prefix) + ".qv")
    completeness = Path(str(prefix) + ".completeness.stats")
    arguments = [str(read_database_path), str(assembly_path), str(prefix)]
    log = run(program, arguments, directory, config, outputs=(qv, completeness))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "evaluate_merqury",
            "backend": "merqury",
            "result": {
                "directory": directory,
                "qv": qv,
                "completeness": completeness,
                "artifacts": artifacts(directory),
                "read_database": read_database_path,
                "assembly": assembly_path,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def run_busco(
    input_data: str | Path,
    lineage: str,
    output_directory: str | Path,
    *,
    mode: str = "genome",
    threads: int = 1,
    output_name: str = "busco_run",
    program: str = "busco",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assess genome, transcriptome or protein completeness with BUSCO."""

    input_path_value = input_path(input_data)
    _validate_program(program)
    if not isinstance(lineage, str) or not lineage.strip() or len(lineage) > 256:
        raise ConfigurationError("lineage must be a non-empty string of at most 256 characters.")
    if any(character.isspace() for character in lineage) or any(
        character in lineage for character in (";", "|", "&", "\x00")
    ):
        raise ConfigurationError("lineage must be a BUSCO dataset name or path without shell characters.")
    if mode not in {"genome", "transcriptome", "proteins"}:
        raise ConfigurationError("mode must be genome, transcriptome or proteins.")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    if not isinstance(output_name, str) or not output_name.strip() or any(
        character in output_name for character in ("/", "\\", "\x00", "\n", "\r")
    ):
        raise ConfigurationError("output_name must be a simple directory name.")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    arguments = [
        "-i", str(input_path_value), "-l", lineage.strip(), "-m", mode,
        "-o", output_name.strip(), "--out_path", str(directory), "-c", str(threads),
    ]
    log = run(program, arguments, directory, config)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "run_busco",
            "backend": "busco",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "input": input_path_value,
                "lineage": lineage.strip(),
                "mode": mode,
                "output_name": output_name.strip(),
                "threads": threads,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def assemble_flye(
    reads: str | Path,
    output_directory: str | Path,
    *,
    genome_size: str,
    platform: str = "nano-hq",
    threads: int = 1,
    program: str = "flye",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assemble long reads with caller-installed Flye."""

    reads_path = input_path(reads)
    program = _validate_program(program)
    if not isinstance(genome_size, str) or not genome_size.strip() or len(genome_size) > 32:
        raise ConfigurationError("genome_size must be a non-empty Flye size expression.")
    if any(character.isspace() for character in genome_size) or any(
        character in genome_size for character in (";", "|", "&", "\x00")
    ):
        raise ConfigurationError("genome_size must not contain shell characters or whitespace.")
    if platform not in {"nano-raw", "nano-hq", "pacbio-raw", "pacbio-hifi", "pacbio-corr"}:
        raise ConfigurationError("platform must be a supported Flye read mode.")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    assembly = directory / "assembly.fasta"
    arguments = [f"--{platform}", str(reads_path), "--genome-size", genome_size.strip(), "--out-dir", str(directory), "--threads", str(threads)]
    log = run(program, arguments, directory, config, outputs=(assembly,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "assemble_flye",
            "backend": "flye",
            "result": {
                "directory": directory,
                "assembly": assembly,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "genome_size": genome_size.strip(),
                "platform": platform,
                "threads": threads,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def assemble_megahit(
    read1: str | Path,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    threads: int = 1,
    min_contig_length: int | None = None,
    program: str = "megahit",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assemble short reads with caller-installed MEGAHIT."""

    read1_path = input_path(read1)
    read2_path = input_path(read2) if read2 is not None else None
    program = _validate_program(program)
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    if min_contig_length is not None and (
        isinstance(min_contig_length, bool)
        or not isinstance(min_contig_length, int)
        or not 1 <= min_contig_length <= 10_000_000
    ):
        raise ConfigurationError("min_contig_length must be an integer in [1, 10000000].")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    megahit_directory = directory / "megahit"
    assembly = megahit_directory / "final.contigs.fa"
    arguments = ["-o", str(megahit_directory), "-t", str(threads)]
    if min_contig_length is not None:
        arguments.extend(("--min-contig-len", str(min_contig_length)))
    if read2_path is None:
        arguments.extend(("-r", str(read1_path)))
    else:
        arguments.extend(("-1", str(read1_path), "-2", str(read2_path)))
    log = run(program, arguments, directory, config, outputs=(assembly,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "assemble_megahit",
            "backend": "megahit",
            "result": {
                "directory": directory,
                "assembly": assembly,
                "artifacts": artifacts(directory),
                "read1": read1_path,
                "read2": read2_path,
                "threads": threads,
                "min_contig_length": min_contig_length,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def assemble_hifiasm(
    reads: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    prefix: str = "assembly",
    program: str = "hifiasm",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Build a HiFi assembly graph with caller-installed hifiasm."""

    reads_path = input_path(reads)
    program = _validate_program(program)
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    if not isinstance(prefix, str) or not prefix.strip() or any(
        character in prefix for character in ("/", "\\", "\x00", "\n", "\r")
    ):
        raise ConfigurationError("prefix must be a simple output stem.")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )

    directory = workspace(output_directory)
    output_prefix = directory / prefix.strip()
    graph = Path(str(output_prefix) + ".bp.p_ctg.gfa")
    arguments = ["-o", str(output_prefix), "-t", str(threads), str(reads_path)]
    log = run(program, arguments, directory, config, outputs=(graph,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "assemble_hifiasm",
            "backend": "hifiasm",
            "result": {
                "directory": directory,
                "primary_contig_graph": graph,
                "artifacts": artifacts(directory),
                "reads": reads_path,
                "threads": threads,
                "prefix": prefix.strip(),
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def assemble_unicycler(
    read1: str | Path | None,
    output_directory: str | Path,
    *,
    read2: str | Path | None = None,
    long_reads: str | Path | None = None,
    mode: str = "normal",
    threads: int = 1,
    program: str = "unicycler",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Assemble bacterial reads with caller-installed Unicycler.

    At least short-read ``read1`` or ``long_reads`` must be supplied.  If
    ``read2`` is supplied it is paired with ``read1``; Unicycler and all
    assemblers it invokes remain caller-installed dependencies.
    """

    if read1 is None and long_reads is None:
        raise ConfigurationError("read1 or long_reads must be supplied.")
    read1_path = input_path(read1) if read1 is not None else None
    read2_path = input_path(read2) if read2 is not None else None
    long_reads_path = input_path(long_reads) if long_reads is not None else None
    if read2_path is not None and read1_path is None:
        raise ConfigurationError("read2 requires read1.")
    if mode not in {"conservative", "normal", "bold"}:
        raise ConfigurationError("mode must be conservative, normal or bold.")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    program = _validate_program(program)
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    result_directory = directory / "unicycler"
    assembly = result_directory / "assembly.fasta"
    arguments: list[str] = []
    if read1_path is not None:
        arguments.extend(("-1", str(read1_path)))
    if read2_path is not None:
        arguments.extend(("-2", str(read2_path)))
    if long_reads_path is not None:
        arguments.extend(("-l", str(long_reads_path)))
    arguments.extend(("-o", str(result_directory), "-m", mode, "-t", str(threads)))
    log = run(program, arguments, directory, config, outputs=(assembly,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "assemble_unicycler",
            "backend": "unicycler",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "assembly": assembly,
                "read1": read1_path,
                "read2": read2_path,
                "long_reads": long_reads_path,
                "mode": mode,
                "threads": threads,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.assembly.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )
__all__ = [
    "assemble_flye",
    "assemble_hifiasm",
    "assemble_megahit",
    "assemble_shasta",
    "assemble_unicycler",
    "evaluate_merqury",
    "run_busco",
]

_register_public_verbose(globals(), tuple(__all__))
