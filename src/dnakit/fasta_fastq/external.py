"""Bounded facade for caller-installed SeqMagick conversion workflows."""

from __future__ import annotations

import re
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, input_path, run, workspace
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent

_FORMATS = {
    "afa": "afa",
    "aln": "aln",
    "clustal": "aln",
    "emboss": "needle",
    "fasta": "fasta",
    "fastq": "fastq",
    "genbank": "gb",
    "nexus": "nex",
    "phylip": "phy",
    "phylip-relaxed": "phyx",
    "qual": "qual",
    "sff-trim": "sff",
    "stockholm": "sto",
}
_ALPHABETS = {"dna", "dna-ambiguous", "rna", "rna-ambiguous", "protein"}


def convert_seqmagick(
    source: str | Path,
    output_directory: str | Path,
    *,
    output_format: str = "fasta",
    input_format: str | None = None,
    alphabet: str | None = None,
    line_wrap: int | None = None,
    ungap: bool = False,
    reverse: bool = False,
    reverse_complement: bool = False,
    upper: bool = False,
    lower: bool = False,
    min_length: int | None = None,
    max_length: int | None = None,
    head: int | None = None,
    tail: int | None = None,
    deduplicate_sequences: bool = False,
    program: str = "seqmagick",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Convert and apply bounded SeqMagick sequence transformations.

    The caller supplies the GPL-3.0 ``seqmagick`` executable.  Format names
    and transformations are allow-listed so a JSON/MCP caller cannot inject
    arbitrary command-line options or overwrite an existing output directory.
    """

    source_path = input_path(source)
    _validate_program(program, "program")
    output_format = _validate_format(output_format, "output_format")
    if input_format is not None:
        input_format = _validate_format(input_format, "input_format")
    if alphabet is not None and alphabet not in _ALPHABETS:
        raise ConfigurationError(f"alphabet must be one of {sorted(_ALPHABETS)}.")
    if output_format == "nexus" and alphabet is None:
        raise ConfigurationError("alphabet is required when output_format is nexus.")
    _validate_bool(ungap, "ungap")
    _validate_bool(reverse, "reverse")
    _validate_bool(reverse_complement, "reverse_complement")
    _validate_bool(upper, "upper")
    _validate_bool(lower, "lower")
    _validate_bool(deduplicate_sequences, "deduplicate_sequences")
    if reverse and reverse_complement:
        raise ConfigurationError("reverse and reverse_complement cannot both be true.")
    if upper and lower:
        raise ConfigurationError("upper and lower cannot both be true.")
    _validate_integer(line_wrap, "line_wrap", minimum=0, maximum=10_000_000)
    _validate_integer(min_length, "min_length", minimum=0, maximum=10_000_000_000)
    _validate_integer(max_length, "max_length", minimum=0, maximum=10_000_000_000)
    _validate_integer(head, "head", minimum=0, maximum=10_000_000)
    _validate_integer(tail, "tail", minimum=0, maximum=10_000_000)
    if min_length is not None and max_length is not None and min_length > max_length:
        raise ConfigurationError("min_length cannot exceed max_length.")
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )

    directory = workspace(output_directory)
    suffix = _FORMATS.get(output_format, "out")
    target = directory / f"converted.{suffix}"
    arguments = ["convert"]
    if line_wrap is not None:
        arguments.extend(("--line-wrap", str(line_wrap)))
    for flag, enabled in (
        ("--ungap", ungap),
        ("--reverse", reverse),
        ("--reverse-complement", reverse_complement),
        ("--upper", upper),
        ("--lower", lower),
        ("--deduplicate-sequences", deduplicate_sequences),
    ):
        if enabled:
            arguments.append(flag)
    if min_length is not None:
        arguments.extend(("--min-length", str(min_length)))
    if max_length is not None:
        arguments.extend(("--max-length", str(max_length)))
    if head is not None:
        arguments.extend(("--head", str(head)))
    if tail is not None:
        arguments.extend(("--tail", str(tail)))
    if input_format is not None:
        arguments.extend(("--input-format", input_format))
    if output_format is not None:
        arguments.extend(("--output-format", output_format))
    if alphabet is not None:
        arguments.extend(("--alphabet", alphabet))
    arguments.extend((str(source_path), str(target)))
    log = run(program, arguments, directory, config, outputs=(target,))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "convert_seqmagick",
            "backend": "seqmagick",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "source": source_path,
                "output": target,
                "output_format": output_format,
                "input_format": input_format,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.fasta_fastq.external",
                "execution": "caller-installed GPL-3.0 executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )



def cluster_vsearch(
    source: str | Path,
    output_directory: str | Path,
    *,
    identity: float = 0.97,
    threads: int = 1,
    program: str = "vsearch",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Cluster FASTA sequences with caller-installed VSEARCH."""

    source_path = input_path(source)
    _validate_program(program, "program")
    if isinstance(identity, bool) or not isinstance(identity, (int, float)) or not 0 < identity <= 1:
        raise ConfigurationError("identity must be a number in (0, 1].")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    centroids = directory / "centroids.fasta"
    uc = directory / "clusters.uc"
    arguments = [
        "--cluster_fast",
        str(source_path),
        "--id",
        str(identity),
        "--threads",
        str(threads),
        "--centroids",
        str(centroids),
        "--uc",
        str(uc),
    ]
    log = run(program, arguments, directory, config, outputs=(centroids, uc))
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "cluster_vsearch",
            "backend": "vsearch",
            "result": {
                "directory": directory,
                "artifacts": artifacts(directory),
                "source": source_path,
                "centroids": centroids,
                "clusters": uc,
                "identity": identity,
                "threads": threads,
                "log": log,
            },
            "provenance": {
                "implementation": "dnakit.fasta_fastq.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )

def _validate_format(value: str, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > 64
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value) is None
    ):
        raise ConfigurationError(f"{name} must be a bounded SeqMagick format name.")
    return value.casefold()


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


__all__ = ["cluster_vsearch", "convert_seqmagick"]

_register_public_verbose(globals(), tuple(__all__))
