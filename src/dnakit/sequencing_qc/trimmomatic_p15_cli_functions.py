"""Trimmomatic command-line facades: read trimming and filtering steps, quality re-encoding and mate pairing.

Trimmomatic (0.39 is the tested version) is a Java program. The user installs it (``conda install -c bioconda trimmomatic`` gives a
``trimmomatic`` launcher) and either puts that launcher on PATH or passes ``trimmomatic_jar`` (the path of the
``trimmomatic-0.39.jar`` file; ``java`` must be on PATH). Each function runs the program with a fixed argument whitelist and no
shell, writes into a new output folder and returns a JSON-serialisable dictionary (status, backend, version, command,
outputs, parsed read-survival summary, log tail).

Adapter clipping, ``LEADING``, ``TRAILING`` and ``MINLEN`` for short fixed subsets are in
``dnakit.sequencing_qc.short_reads.trim_trimmomatic``; this module adds the remaining steps and the Pairomatic tool.
"""

from __future__ import annotations

import dataclasses
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="trimmomatic",
    executable="trimmomatic",
    hint=(
        "Install Trimmomatic 0.39 or later (for example `conda install -c bioconda trimmomatic`) and put the trimmomatic launcher on "
        "PATH, or pass trimmomatic_jar (the .jar file; java must be on PATH)."
    ),
    tested_version="0.39",
    citation="https://doi.org/10.1093/bioinformatics/btu170",
    version_arguments=("-version",),
    version_pattern=r"(\d+\.\d+(?:\.\d+)?)",
)
_JAVA_HINT = "Install a Java runtime (java must be on PATH) to run a Trimmomatic .jar file."
_SUMMARY = re.compile(r"^([A-Za-z][A-Za-z ]*?):\s*(\d+(?:\.\d+)?)\s*$", re.MULTILINE)
_BASES = re.compile(r"^[ACGTNacgtn]+$")


def _runner(jar: str | None) -> tuple[_cli.Tool, str | None, tuple[str, ...]]:
    """Launcher mode (``trimmomatic`` on PATH) or jar mode (``java -jar <jar>``)."""
    if jar is None:
        return _TOOL, None, ()
    path = _cli.existing_file(jar, "trimmomatic_jar")
    tool = dataclasses.replace(_TOOL, executable="java", hint=_JAVA_HINT, version_arguments=("-jar", str(path), "-version"))
    return tool, "java", ("-jar", str(path))


def _summary(path: Path) -> dict[str, Any]:
    """The ``key: number`` lines of Trimmomatic's ``-summary`` file (for example ``surviving_reads``, ``dropped_read_percent``)."""
    if not path.is_file():
        return {}
    return {
        key.strip().lower().replace(" ", "_"): float(value) if "." in value else int(value)
        for key, value in _SUMMARY.findall(path.read_text(encoding="utf-8", errors="replace"))
    }


def trimmomatic_process_reads(
    read1: str,
    output_directory: str,
    *,
    read2: str | None = None,
    input_phred: int | None = None,
    convert_quality: int | None = None,
    crop: int | None = None,
    head_crop: int | None = None,
    tail_crop: int | None = None,
    leading: int | None = None,
    trailing: int | None = None,
    sliding_window: Sequence[float] | None = None,
    max_info: Sequence[float] | None = None,
    average_quality: int | None = None,
    base_count: Sequence[Any] | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    threads: int = 1,
    trimmomatic_jar: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim and filter single-end or paired-end FASTQ reads with Trimmomatic steps (``SLIDINGWINDOW``, ``MAXINFO``, ``MAXLEN``, ...).

    Steps run in this fixed order, which is the order Trimmomatic's manual recommends: ``convert_quality`` (``TOPHRED33`` /
    ``TOPHRED64``: give 33 or 64), ``crop`` (``CROP``), ``head_crop`` (``HEADCROP``), ``tail_crop`` (``TAILCROP``), ``leading`` and
    ``trailing`` (quality cut-offs), ``sliding_window`` (``[window_size, required_quality]``: cut when the window average falls
    below the quality), ``max_info`` (``[target_length, strictness]``: adaptive trimming that balances read length and
    error rate, strictness 0 to 1), ``average_quality`` (``AVGQUAL``: drop reads with a lower mean quality), ``base_count``
    (``BASECOUNT``: ``[bases, minimum]`` or ``[bases, minimum, maximum]``, drops reads whose count of those bases is outside
    the range), ``min_length`` (``MINLEN``) and ``max_length`` (``MAXLEN``). At least one step is required. ``input_phred``
    (33 or 64) sets the input quality encoding (default: Trimmomatic detects it). With ``read2`` the reads are processed as
    pairs and four files are written (paired and unpaired for each mate), otherwise one. Output goes to the new folder
    ``output_directory``; ``survival`` holds the counts from Trimmomatic's ``-summary`` (input, surviving, dropped reads).
    """
    first = _cli.existing_file(read1, "read1")
    second = _cli.existing_file(read2, "read2") if read2 is not None else None
    if second is not None and second == first:
        raise ConfigurationError("read1 and read2 must be different files.")
    steps: list[str] = []
    if convert_quality is not None:
        steps.append(f"TOPHRED{_cli.choice(str(convert_quality), ('33', '64'), 'convert_quality')}")
    for name, value in (("CROP", crop), ("HEADCROP", head_crop), ("TAILCROP", tail_crop), ("LEADING", leading), ("TRAILING", trailing)):
        if value is not None:
            steps.append(f"{name}:{_cli.integer(value, name.lower(), minimum=0)}")
    if sliding_window is not None:
        if len(sliding_window) != 2:
            raise ConfigurationError("sliding_window must be [window_size, required_quality].")
        steps.append(
            f"SLIDINGWINDOW:{_cli.integer(int(sliding_window[0]), 'sliding_window window_size', minimum=1)}:"
            f"{format(_cli.number(sliding_window[1], 'sliding_window quality', minimum=0.0), '.12g')}"
        )
    if max_info is not None:
        if len(max_info) != 2:
            raise ConfigurationError("max_info must be [target_length, strictness].")
        steps.append(
            f"MAXINFO:{_cli.integer(int(max_info[0]), 'max_info target_length', minimum=1)}:"
            f"{format(_cli.number(max_info[1], 'max_info strictness', minimum=0.0, maximum=1.0), '.12g')}"
        )
    if average_quality is not None:
        steps.append(f"AVGQUAL:{_cli.integer(average_quality, 'average_quality', minimum=0)}")
    if base_count is not None:
        if len(base_count) not in (2, 3):
            raise ConfigurationError("base_count must be [bases, minimum] or [bases, minimum, maximum].")
        bases = base_count[0]
        if not isinstance(bases, str) or not _BASES.match(bases):
            raise ConfigurationError("base_count bases must be a string of nucleotide letters such as 'N' or 'GC'.")
        counts = [str(_cli.integer(int(item), "base_count", minimum=0)) for item in base_count[1:]]
        steps.append("BASECOUNT:" + ":".join([bases.upper(), *counts]))
    if min_length is not None:
        steps.append(f"MINLEN:{_cli.integer(min_length, 'min_length', minimum=0)}")
    if max_length is not None:
        steps.append(f"MAXLEN:{_cli.integer(max_length, 'max_length', minimum=1)}")
    if not steps:
        raise ConfigurationError("give at least one trimming or filtering step.")
    out = _cli.new_directory(output_directory, "output_directory")
    summary = out / "summary.txt"
    arguments = ["PE" if second is not None else "SE", "-threads", str(_cli.threads_value(threads))]
    if input_phred is not None:
        arguments.append(f"-phred{_cli.choice(str(input_phred), ('33', '64'), 'input_phred')}")
    arguments += ["-summary", str(summary)]
    if second is None:
        outputs = [out / "trimmed.fastq.gz"]
        arguments += [str(first), str(outputs[0])]
    else:
        outputs = [out / "R1.paired.fastq.gz", out / "R1.unpaired.fastq.gz", out / "R2.paired.fastq.gz", out / "R2.unpaired.fastq.gz"]
        arguments += [str(first), str(second), *map(str, outputs)]
    arguments += steps
    tool, name, prefix = _runner(trimmomatic_jar)
    record = _cli.invoke(
        tool, arguments, cwd=out, expected=[*outputs, summary], executable_name=name, prefix=prefix, timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record, "trim_and_filter_reads_with_trimmomatic_steps", output_directory=str(out), steps=steps, survival=_summary(summary),
        read_outputs={path.name: str(path) for path in outputs},
    )


def trimmomatic_pair_reads(
    read1: str,
    read2: str,
    output_directory: str,
    *,
    delimiter: str | None = None,
    trimmomatic_jar: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pair FASTQ mates by read name and split off the reads without a mate (Trimmomatic's ``Pairomatic`` tool).

    Reads of ``read1`` and ``read2`` with the same name are written as pairs to ``R1.paired.fastq`` and ``R2.paired.fastq``;
    reads without a mate go to ``R1.unpaired.fastq`` and ``R2.unpaired.fastq``. By default whole names are compared; give
    ``delimiter`` (one character, ``'/'`` for ``name/1`` and ``name/2``, ``' '`` for ``name 1:N:0`` and ``name 2:N:0``) to compare
    only the part before it. The reads that occur in both files must be in the same order in both, otherwise Pairomatic refuses. Output goes to the new folder ``output_directory``; the result gives the read count of every file.
    Pairomatic is part of the Trimmomatic jar, so this function needs ``trimmomatic_jar``.
    """
    first = _cli.existing_file(read1, "read1")
    second = _cli.existing_file(read2, "read2")
    if first == second:
        raise ConfigurationError("read1 and read2 must be different files.")
    jar = _cli.existing_file(trimmomatic_jar, "trimmomatic_jar") if trimmomatic_jar is not None else None
    if jar is None:
        raise ConfigurationError("trimmomatic_jar is required: Pairomatic is started from the Trimmomatic jar file.")
    out = _cli.new_directory(output_directory, "output_directory")
    outputs = [out / "R1.paired.fastq", out / "R1.unpaired.fastq", out / "R2.paired.fastq", out / "R2.unpaired.fastq"]
    arguments = ["-cp", str(jar), "org.usadellab.trimmomatic.Pairomatic"]
    if delimiter is not None:
        if not isinstance(delimiter, str) or len(delimiter) != 1 or delimiter in "\n\r\x00":
            raise ConfigurationError("delimiter must be a single character such as '/' or ' '.")
        arguments += ["-delim", delimiter]
    arguments += [str(first), str(second), *map(str, outputs)]
    tool = dataclasses.replace(_TOOL, executable="java", hint=_JAVA_HINT, version_arguments=("-jar", str(jar), "-version"))
    record = _cli.invoke(tool, arguments, cwd=out, expected=outputs, executable_name="java", timeout_seconds=timeout_seconds)
    return _cli.finish(
        record, "pair_fastq_mates_by_read_name", output_directory=str(out),
        read_counts={path.name: _cli.count_fastx_records(path) for path in outputs},
        read_outputs={path.name: str(path) for path in outputs},
    )


__all__ = ["trimmomatic_pair_reads", "trimmomatic_process_reads"]
_register_public_verbose(globals())
