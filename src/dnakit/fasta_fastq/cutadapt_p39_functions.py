"""Fixed-option cutadapt read filters, header editing and demultiplexing (P39)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.exceptions import ConfigurationError


def _run(
    source: str,
    output: str,
    options: list[str],
    mate: str | None,
    mate_output: str | None,
    executable_path: str | None,
) -> dict[str, Any]:
    reads = cli.existing_file(source, "source")
    target = cli.new_file(output, "output")
    arguments = [*options, "-o", str(target)]
    expected = [target]
    if (mate is None) != (mate_output is None):
        raise ConfigurationError("Supply both mate and mate_output, or neither.")
    if mate_output is not None:
        target2 = cli.new_file(mate_output, "mate_output")
        if target2 == target:
            raise ConfigurationError("Mate outputs must be distinct.")
        arguments += ["-p", str(target2)]
        expected.append(target2)
    arguments.append(str(reads))
    if mate is not None:
        arguments.append(str(cli.existing_file(mate, "mate")))
    return cli.run(
        "cutadapt",
        arguments,
        cwd=target.parent,
        expected=expected,
        explicit=executable_path,
        hint='python -m pip install "dnakit[bio-cutadapt]"',
    )


def cutadapt_filter_reads(
    source: str,
    output: str,
    *,
    max_n: float | None = None,
    max_expected_errors: float | None = None,
    max_average_error_rate: float | None = None,
    adapter: str | None = None,
    adapter_filter: str | None = None,
    discard_casava: bool = False,
    mate: str | None = None,
    mate_output: str | None = None,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter reads by N content, expected errors, average error rate, adapter presence or CASAVA chastity.

    adapter_filter is trimmed (discard adapter-containing reads) or untrimmed
    (discard reads without adapters). Other thresholds preserve upstream units.
    """
    arguments = []
    for flag, value in (
        ("--max-n", max_n),
        ("--max-ee", max_expected_errors),
        ("--max-aer", max_average_error_rate),
    ):
        if value is not None:
            if value < 0:
                raise ConfigurationError(f"{flag} must be nonnegative.")
            arguments += [flag, str(value)]
    if adapter is not None:
        arguments += ["-a", adapter]
    if adapter_filter is not None:
        if adapter_filter not in ("trimmed", "untrimmed") or adapter is None:
            raise ConfigurationError("adapter_filter requires an adapter and trimmed/untrimmed.")
        arguments.append("--discard-" + adapter_filter)
    if discard_casava:
        arguments.append("--discard-casava")
    return _run(source, output, arguments, mate, mate_output, executable_path)


def cutadapt_rewrite_read_names(
    source: str,
    output: str,
    *,
    rename: str | None = None,
    length_tag: str | None = None,
    strip_suffix: str | None = None,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rewrite read headers using a rename template, length tag or stripped suffix."""
    arguments = []
    for flag, value in (
        ("--rename", rename),
        ("--length-tag", length_tag),
        ("--strip-suffix", strip_suffix),
    ):
        if value is not None:
            arguments += [flag, value]
    if not arguments:
        raise ConfigurationError("Select at least one header rewriting operation.")
    return _run(source, output, arguments, None, None, executable_path)


def cutadapt_clamp_negative_qualities(
    source: str,
    output: str,
    *,
    quality_base: int = 33,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clamp negative Phred qualities to zero with cutadapt --zero-cap."""
    if quality_base not in (33, 64):
        raise ConfigurationError("quality_base must be 33 or 64.")
    return _run(
        source,
        output,
        ["--zero-cap", "--quality-base", str(quality_base)],
        None,
        None,
        executable_path,
    )


def cutadapt_interleaved_reads(
    source: str,
    output: str,
    *,
    adapter_first: str | None = None,
    adapter_second: str | None = None,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read and write interleaved paired-end FASTQ, optionally trimming each mate's 3' adapter."""
    options = ["--interleaved"]
    for flag, value in (("-a", adapter_first), ("-A", adapter_second)):
        if value is not None:
            options += [flag, value]
    return _run(source, output, options, None, None, executable_path)


def cutadapt_demultiplex_reads(
    source: str,
    output_dir: str,
    adapters: Mapping[str, str],
    *,
    mate: str | None = None,
    mate_adapters: Mapping[str, str] | None = None,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Single-index or combinatorial dual-index demultiplexing into a new directory.

    Adapter names are simple filenames. Supplying mate_adapters selects
    combinatorial {name1}-{name2} templates; otherwise both mates use {name}.
    """
    import re

    reads = cli.existing_file(source, "source")
    directory = Path(output_dir).expanduser().resolve()
    if directory.exists() or not directory.parent.is_dir():
        raise ConfigurationError("output_dir must be new and its parent must exist.")
    if not adapters or any(
        not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in [*adapters, *(mate_adapters or {})]
    ):
        raise ConfigurationError("Provide adapters with simple alphanumeric names.")
    if mate_adapters and not mate:
        raise ConfigurationError("mate_adapters requires paired reads.")
    options = []
    for flag, values in (("-g", adapters), ("-G", mate_adapters or {})):
        for name, sequence in values.items():
            options += [flag, name + "=" + sequence]
    template = "{name1}-{name2}" if mate_adapters else "{name}"
    options += ["-o", str(directory / (template + (".1.fastq" if mate else ".fastq")))]
    if mate:
        options += ["-p", str(directory / (template + ".2.fastq"))]
    options.append(str(reads))
    if mate:
        options.append(str(cli.existing_file(mate, "mate")))
    # Resolve first so a missing executable leaves no output directory.
    executable = cli.program(
        "cutadapt", explicit=executable_path, hint='python -m pip install "dnakit[bio-cutadapt]"'
    )
    directory.mkdir()
    result = cli.run(
        "cutadapt",
        options,
        cwd=directory,
        expected=[],
        explicit=executable,
        hint='python -m pip install "dnakit[bio-cutadapt]"',
    )
    result["outputs"] = [str(p) for p in sorted(directory.glob("*.fastq"))]
    return result


__all__ = [
    "cutadapt_clamp_negative_qualities",
    "cutadapt_demultiplex_reads",
    "cutadapt_filter_reads",
    "cutadapt_interleaved_reads",
    "cutadapt_rewrite_read_names",
]

_register_public_verbose(globals(), tuple(__all__))
