"""BEDOPS command-line facades: fixed argument whitelists, no shell, bounded execution.

Each public function runs one BEDOPS program (2.4.42 is the tested version) with an explicit argument
list: ``bedops`` (set operations), ``bedextract``, ``closest-features``, ``bedmap``, ``sort-bed``,
``starch``/``unstarch``/``starchcat``/``starchstrip`` and ``convert2bed``. The user installs BEDOPS
(conda, apt or the release archive) and either puts its ``bin`` folder on PATH or passes
``bedops_bin_dir``. BEDOPS programs write their result to standard output; the facades stream it to an
explicit new output file (nothing is overwritten) and return a JSON-serialisable dictionary with the
status, backend, version, command, outputs, a small parsed preview and a bounded log tail.

Every BEDOPS input file must be sorted as ``sort-bed`` sorts it; the facades do not sort for you
(``error_check=True`` makes BEDOPS verify the order). ``convert2bed`` for ``sam`` and ``bam`` input
calls ``samtools``, which must be on PATH.
"""

from __future__ import annotations

import contextlib
import math
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
    DNAKitError,
)

_TESTED_VERSION = "2.4.42"
_CITATION = "https://doi.org/10.1093/bioinformatics/bts277"
_INSTALL_HINT = (
    "Install BEDOPS 2.4.42 or later (for example `conda install -c bioconda bedops`, or the release "
    "archive from https://github.com/bedops/bedops/releases), then put its bin folder on PATH or pass "
    "bedops_bin_dir."
)
_SAMTOOLS_HINT = "convert2bed needs samtools on PATH to read SAM and BAM (for example `conda install -c bioconda samtools`)."
_TIMEOUT_SECONDS = 3600.0
_VERSION_TIMEOUT_SECONDS = 60.0
_MAX_LOG_BYTES = 5_000_000
_MAX_CAPTURE_BYTES = 50_000_000
_MAX_ARTIFACT_BYTES = 2_000_000_000
_LOG_TAIL_CHARS = 4000
_PREVIEW_LINES = 10
_VERSION_LINE = re.compile(r"version:\s*(\d+(?:\.\d+)+)")
_OVERLAP_AMOUNT = re.compile(r"^\d+%?$")
_RANGE = re.compile(r"^[+-]?\d+(?::[+-]?\d+)?$")
_MEMORY = re.compile(r"^\d+[KMGkmg]?$")

SET_OPERATIONS = ("element-of", "not-element-of", "symmdiff", "everything")
SCORE_STATISTICS = (
    "cv",
    "kth",
    "mad",
    "max",
    "mean",
    "median",
    "min",
    "stdev",
    "sum",
    "tmean",
    "variance",
    "wmean",
    "count",
    "bases",
    "bases-uniq",
    "bases-uniq-f",
    "indicator",
)
ECHO_FIELDS = (
    "echo",
    "echo-map",
    "echo-map-id",
    "echo-map-id-uniq",
    "echo-map-range",
    "echo-map-score",
    "echo-map-size",
    "echo-overlap-size",
    "echo-ref-name",
    "echo-ref-row-id",
    "echo-ref-size",
    "max-element",
    "max-element-rand",
    "min-element",
    "min-element-rand",
)
OVERLAP_MODES = (
    "bp-ovr",
    "exact",
    "fraction-both",
    "fraction-either",
    "fraction-map",
    "fraction-ref",
    "range",
)
UNSTARCH_STATISTICS = (
    "elements",
    "elements-max-string-length",
    "bases",
    "bases-uniq",
    "has-duplicate",
    "has-nested",
    "list-chromosomes",
    "list-json",
)
CONVERT_INPUT_FORMATS = ("bam", "gff", "gtf", "gvf", "psl", "rmsk", "sam", "vcf", "wig")
CONVERT_OUTPUT_FORMATS = ("bed", "starch")
# Boolean convert2bed flags the facade accepts per input format (from ``convert2bed --help-<fmt>``).
_CONVERT_FLAGS: dict[str, frozenset[str]] = {
    "bam": frozenset({"all-reads", "keep-header", "split", "split-with-deletions", "reduced"}),
    "sam": frozenset({"all-reads", "keep-header", "split", "split-with-deletions", "reduced"}),
    "gff": frozenset({"keep-header"}),
    "gtf": frozenset(),
    "gvf": frozenset({"keep-header"}),
    "psl": frozenset({"keep-header", "split"}),
    "rmsk": frozenset({"keep-header"}),
    "vcf": frozenset({"do-not-split", "snvs", "insertions", "deletions", "keep-header"}),
    "wig": frozenset({"keep-header", "zero-indexed"}),
}
# Attribute keys convert2bed accepts for ``--attribute-key`` (from ``convert2bed --help-gff`` / ``--help-gtf``).
_ATTRIBUTE_KEYS: dict[str, tuple[str, ...]] = {
    "gff": (
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_id",
        "transcript_name",
        "transcript_type",
        "exon_id",
        "exon_number",
        "ccdsid",
        "protein_id",
        "ont",
        "hgnc_id",
        "mgi_id",
        "havana_gene",
        "havana_transcript",
    ),
    "gtf": (
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_id",
        "transcript_name",
        "exon_id",
        "havana_gene",
        "havana_transcript",
    ),
}
_VERSIONS: dict[str, str | None] = {}


# ---------------------------------------------------------------------------
# Validation and path handling
# ---------------------------------------------------------------------------


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _input(value: object, name: str = "input_path") -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    return _gap_cli.existing_file(value, name)


def _inputs(values: object, name: str, *, minimum: int) -> tuple[Path, ...]:
    if isinstance(values, (str, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    if len(values) < minimum:
        raise ConfigurationError(f"{name} needs at least {minimum} file(s); got {len(values)}.")
    return tuple(_input(value, name) for value in values)


def _output(value: object, name: str = "output_path") -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    return _gap_cli.new_file(value, name)


def _chromosome(value: object, name: str = "chromosome") -> str:
    return _gap_cli.text(value, name)


def _delimiter(value: object, name: str = "delimiter") -> str:
    """A delimiter may be whitespace (a tab), but it can never be read as a command-line flag."""
    if not isinstance(value, str) or not value or "\x00" in value or "\n" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL or newline.")
    if value.startswith("-"):
        raise ConfigurationError(f"{name} must not start with '-'.")
    return value


def _range(value: object) -> str:
    text = str(value)
    if not isinstance(value, (str, int)) or isinstance(value, bool) or not _RANGE.match(text):
        raise ConfigurationError(
            "range must look like '10', '-5:20' or '+3:-2' (bp added to start:end)."
        )
    return text


def _bin_dir(value: object) -> str | None:
    if value is None or value == "":
        return None
    path = Path(str(value)).expanduser()
    if not path.is_dir():
        raise BackendUnavailableError(f"bedops_bin_dir is not a folder: {path}", hint=_INSTALL_HINT)
    return str(path.resolve())


def _executable(tool: str, bin_dir: str | None) -> str:
    explicit = None if bin_dir is None else str(Path(bin_dir) / tool)
    return _gap_cli.program(tool, explicit=explicit, hint=_INSTALL_HINT)


def _version(executable: str) -> str | None:
    if executable not in _VERSIONS:
        try:
            result = execute_bounded_command(
                executable,
                ("--version",),
                backend_id=Path(executable).name,
                cwd=Path(executable).parent,
                timeout_seconds=_VERSION_TIMEOUT_SECONDS,
                max_output_bytes=100_000,
                verbose=False,
            )
            match = _VERSION_LINE.search(result.output or "")
            _VERSIONS[executable] = match.group(1) if match else None
        except DNAKitError:
            _VERSIONS[executable] = None
    return _VERSIONS[executable]


def _common_flags(
    *, chromosome: str | None, error_check: bool, header: bool, range_value: str | int | None = None
) -> list[str]:
    flags: list[str] = []
    if chromosome is not None:
        flags += ["--chrom", _chromosome(chromosome)]
    if range_value is not None:
        flags += ["--range", _range(range_value)]
    if _flag(error_check, "error_check"):
        flags.append("--ec")
    if _flag(header, "header"):
        flags.append("--header")
    return flags


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Completed:
    """Result of one bounded BEDOPS run."""

    tool: str
    return_code: int
    log: str
    executable: str
    version: str | None
    command: tuple[str, ...]
    stdout_path: Path | None


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


def _run(
    tool: str,
    arguments: Sequence[str],
    *,
    bin_dir: str | None,
    stdout: Path | None,
    stdin: Path | None = None,
    capture: bool = False,
    accept: Sequence[int] = (0,),
    version_tool: str | None = None,
) -> _Completed:
    """Run one fixed BEDOPS program without a shell.

    Standard output goes to ``stdout`` (a new file the facade owns), or to a temporary file that
    ``capture=True`` makes the caller read through ``_Completed.stdout_path`` before it disappears.
    """
    executable = _executable(tool, bin_dir)
    version = _version(_executable(version_tool, bin_dir) if version_tool else executable)
    command = (executable, *arguments)
    created = stdout is not None
    with tempfile.TemporaryDirectory(prefix="dnakit-bedops-") as scratch:
        scratch_dir = Path(scratch)
        target = stdout if stdout is not None else scratch_dir / "stdout"
        diagnostic = scratch_dir / "stderr"
        limit = _MAX_ARTIFACT_BYTES if created else _MAX_CAPTURE_BYTES
        input_handle = stdin.open("rb") if stdin is not None else None
        try:
            with target.open("xb") as out_handle, diagnostic.open("xb") as err_handle:
                process = subprocess.Popen(
                    command,
                    stdin=input_handle if input_handle is not None else subprocess.DEVNULL,
                    stdout=out_handle,
                    stderr=err_handle,
                    cwd=scratch_dir,
                    shell=False,
                    start_new_session=os.name == "posix",
                )
                deadline = time.monotonic() + _TIMEOUT_SECONDS
                try:
                    while process.poll() is None:
                        if time.monotonic() > deadline:
                            raise BackendTimeoutError(f"{tool} exceeded its time limit.")
                        if (
                            target.stat().st_size > limit
                            or diagnostic.stat().st_size > _MAX_LOG_BYTES
                        ):
                            raise BackendExecutionError(f"{tool} exceeded its output size limit.")
                        time.sleep(0.02)
                finally:
                    if process.poll() is None:
                        if os.name == "posix":
                            with contextlib.suppress(ProcessLookupError):
                                os.killpg(process.pid, signal.SIGKILL)
                        else:
                            process.kill()
                    process.wait()
        except BaseException:
            if created and stdout is not None:
                with contextlib.suppress(OSError):
                    stdout.unlink()
            raise
        finally:
            if input_handle is not None:
                input_handle.close()
        log = diagnostic.read_text(encoding="utf-8", errors="replace")
        if process.returncode not in accept:
            if created and stdout is not None:
                with contextlib.suppress(OSError):
                    stdout.unlink()
            raise BackendExecutionError(
                f"BEDOPS {tool} exited with code {process.returncode}.",
                code="BEDOPS_EXECUTION_FAILED",
                context={"tool": tool, "log_tail": _tail(log)},
            )
        kept: Path | None = stdout
        if capture:
            descriptor, copy_name = tempfile.mkstemp(prefix="dnakit-bedops-out-")
            os.close(descriptor)
            copy = Path(copy_name)
            shutil.copyfile(target, copy)
            kept = copy
    return _Completed(
        tool,
        int(process.returncode),
        log,
        executable,
        version,
        (tool, *arguments),
        kept,
    )


def _captured_text(run: _Completed) -> str:
    """Read (and delete) the standard output captured by ``_run(..., capture=True)``."""
    if run.stdout_path is None:
        return ""
    try:
        return run.stdout_path.read_text(encoding="utf-8", errors="replace")
    finally:
        with contextlib.suppress(OSError):
            run.stdout_path.unlink()


def _result(
    operation: str,
    run: _Completed,
    *,
    outputs: dict[str, Path] | None = None,
    **extra: object,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": "bedops",
        "operation": operation,
        "program": run.tool,
        "version": run.version,
        "command": list(run.command),
        "outputs": {name: str(path) for name, path in (outputs or {}).items()},
        "log_tail": _tail(run.log),
        "provenance": {
            "tool": "BEDOPS",
            "version": run.version,
            "tested_version": _TESTED_VERSION,
            "executable": run.executable,
            "citation": _CITATION,
        },
    }
    result.update(extra)
    return result


def _line_summary(path: Path) -> dict[str, Any]:
    """Line count and the first lines of a text output (BED rows, bedmap columns)."""
    count = 0
    preview: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            count += 1
            if len(preview) < _PREVIEW_LINES:
                preview.append(line.rstrip("\n"))
    return {"line_count": count, "preview": preview}


def _number(value: str) -> object:
    try:
        return int(value)
    except ValueError:
        pass
    try:
        number = float(value)
    except ValueError:
        return value
    return number if math.isfinite(number) else value


# ---------------------------------------------------------------------------
# Interval set operations, bedextract, closest-features, sort-bed
# ---------------------------------------------------------------------------


def bedops_set_operation(
    input_paths: Sequence[str],
    output_path: str,
    *,
    operation: str = "element-of",
    overlap: str | int | None = None,
    range: str | int | None = None,
    chromosome: str | None = None,
    error_check: bool = False,
    header: bool = False,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set operations of ``bedops`` that keep the original BED columns of the selected elements.

    ``operation`` selects the BEDOPS operation: ``element-of`` (``--element-of``: reference elements that
    overlap elements of the other files), ``not-element-of`` (``--not-element-of``: reference elements
    that do not), ``symmdiff`` (``--symmdiff``: symmetric difference of the interval sets) or
    ``everything`` (``--everything``: union of all elements with every column kept). The first path is the
    reference file; ``element-of``, ``not-element-of`` and ``symmdiff`` need at least two files,
    ``everything`` at least one. ``overlap`` is the required overlap for ``element-of`` and
    ``not-element-of``: a number of bp (``1``) or a percentage of the reference element (``"50%"``);
    BEDOPS uses 100% when it is omitted. ``range`` pads the non-reference intervals (``"10"`` or
    ``"-5:20"``). ``bedextract <query> <target>`` is the same as ``element-of`` with ``overlap=1``.
    All inputs must be sorted as ``sort-bed`` sorts them.
    """
    choice = _gap_cli.choice(operation, SET_OPERATIONS, "operation")
    minimum = 1 if choice == "everything" else 2
    inputs = _inputs(input_paths, "input_paths", minimum=minimum)
    target = _output(output_path)
    flags = _common_flags(
        chromosome=chromosome, error_check=error_check, header=header, range_value=range
    )
    operation_args = [f"--{choice}"]
    if overlap is not None:
        if choice not in ("element-of", "not-element-of"):
            raise ConfigurationError("overlap only applies to element-of and not-element-of.")
        text = str(overlap)
        if isinstance(overlap, bool) or not _OVERLAP_AMOUNT.match(text):
            raise ConfigurationError(
                "overlap must be a number of bp such as 1 or a percentage such as '50%'."
            )
        operation_args.append(text)
    arguments = [*flags, *operation_args, *(str(path) for path in inputs)]
    run = _run("bedops", arguments, bin_dir=_bin_dir(bedops_bin_dir), stdout=target)
    return _result(
        f"bedops_{choice.replace('-', '_')}",
        run,
        outputs={"output": target},
        **_line_summary(target),
    )


def bedextract_list_chromosomes(
    input_path: str,
    *,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the unique chromosome names of a sorted BED or Starch file (``bedextract --list-chr``)."""
    source = _input(input_path)
    run = _run(
        "bedextract",
        ["--list-chr", str(source)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=None,
        capture=True,
    )
    chromosomes = [line for line in _captured_text(run).splitlines() if line]
    return _result(
        "bedextract_list_chromosomes",
        run,
        chromosomes=chromosomes,
        chromosome_count=len(chromosomes),
    )


def bedextract_chromosome(
    input_path: str,
    chromosome: str,
    output_path: str,
    *,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write all rows of one chromosome from a sorted BED or Starch file (``bedextract <chromosome> <file>``)."""
    source = _input(input_path)
    name = _chromosome(chromosome)
    target = _output(output_path)
    run = _run(
        "bedextract",
        [name, str(source)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=target,
    )
    return _result(
        "bedextract_chromosome",
        run,
        outputs={"output": target},
        chromosome=name,
        **_line_summary(target),
    )


def bedops_closest_features(
    input_path: str,
    query_path: str,
    output_path: str,
    *,
    closest: bool = False,
    distances: bool = False,
    no_overlaps: bool = False,
    no_reference: bool = False,
    no_query: bool = False,
    delimiter: str = "|",
    chromosome: str | None = None,
    error_check: bool = False,
    header: bool = False,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """For each element of ``input_path`` report the nearest query elements on its left and right (``closest-features``).

    By default every output row is the input element, the nearest element to its left and the nearest
    element to its right, separated by ``delimiter`` (``|``); ``NA`` marks a missing neighbour. An
    overlapping query element wins over non-overlapping ones. ``closest=True`` keeps only the closer
    neighbour (ties go left), ``distances=True`` adds the signed distances, ``no_overlaps`` ignores
    overlapping query elements, ``no_reference`` / ``no_query`` drop the echoed input / query columns.
    Both files must be sorted as ``sort-bed`` sorts them.
    """
    source = _input(input_path, "input_path")
    query = _input(query_path, "query_path")
    target = _output(output_path)
    flags = _common_flags(chromosome=chromosome, error_check=error_check, header=header)
    for enabled, name, option in (
        (closest, "closest", "--closest"),
        (distances, "distances", "--dist"),
        (no_overlaps, "no_overlaps", "--no-overlaps"),
        (no_reference, "no_reference", "--no-ref"),
        (no_query, "no_query", "--no-query"),
    ):
        if _flag(enabled, name):
            flags.append(option)
    flags += ["--delim", _delimiter(delimiter)]
    run = _run(
        "closest-features",
        [*flags, str(source), str(query)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=target,
    )
    return _result(
        "bedops_closest_features", run, outputs={"output": target}, **_line_summary(target)
    )


def bedops_sort_bed_filter(
    input_paths: Sequence[str],
    output_path: str,
    *,
    mode: str = "unique",
    max_memory: str | None = None,
    temp_dir: str | None = None,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort BED files and keep only unique (``--unique``) or only repeated (``--duplicates``) elements (``sort-bed``).

    ``mode="unique"`` prints each distinct element once, like ``sort -u``; ``mode="duplicates"`` prints
    only the elements that occur more than once, like ``uniq -d``. ``max_memory`` (``"8G"``, ``"8000M"``)
    and ``temp_dir`` are the ``--max-mem`` and ``--tmpdir`` options for large inputs.
    """
    choice = _gap_cli.choice(mode, ("unique", "duplicates"), "mode")
    inputs = _inputs(input_paths, "input_paths", minimum=1)
    target = _output(output_path)
    arguments: list[str] = []
    if max_memory is not None:
        if not isinstance(max_memory, str) or not _MEMORY.match(max_memory):
            raise ConfigurationError("max_memory must look like '8G', '8000M' or '8000000000'.")
        arguments += ["--max-mem", max_memory]
    if temp_dir is not None:
        folder = Path(str(temp_dir)).expanduser()
        if not folder.is_dir():
            raise ConfigurationError(f"temp_dir must be an existing folder: {folder}")
        arguments += ["--tmpdir", str(folder.resolve())]
    arguments.append(f"--{choice}" if choice == "duplicates" else "--unique")
    arguments += [str(path) for path in inputs]
    run = _run("sort-bed", arguments, bin_dir=_bin_dir(bedops_bin_dir), stdout=target)
    return _result(f"sort_bed_{choice}", run, outputs={"output": target}, **_line_summary(target))


# ---------------------------------------------------------------------------
# bedmap
# ---------------------------------------------------------------------------


def _bedmap_flags(
    *,
    overlap_mode: str | None,
    overlap_value: float | None,
    delimiter: str | None,
    multi_delimiter: str | None,
    precision: int | None,
    scientific: bool,
    skip_unmapped: bool,
    unmapped_value: str | None,
    chromosome: str | None,
    error_check: bool,
    header: bool,
) -> list[str]:
    flags = _common_flags(chromosome=chromosome, error_check=error_check, header=header)
    if delimiter is not None:
        flags += ["--delim", _delimiter(delimiter)]
    if multi_delimiter is not None:
        flags += ["--multidelim", _delimiter(multi_delimiter, "multi_delimiter")]
    if precision is not None:
        flags += ["--prec", str(_gap_cli.integer(precision, "precision", minimum=0))]
    if _flag(scientific, "scientific"):
        flags.append("--sci")
    if _flag(skip_unmapped, "skip_unmapped"):
        flags.append("--skip-unmapped")
    if unmapped_value is not None:
        flags += ["--unmapped-val", _delimiter(unmapped_value, "unmapped_value")]
    if overlap_mode is None:
        if overlap_value is not None:
            raise ConfigurationError("overlap_value needs overlap_mode.")
        return flags
    mode = _gap_cli.choice(overlap_mode, OVERLAP_MODES, "overlap_mode")
    if mode == "exact":
        if overlap_value is not None:
            raise ConfigurationError("overlap_mode='exact' takes no overlap_value.")
        return [*flags, "--exact"]
    if overlap_value is None:
        raise ConfigurationError(f"overlap_mode={mode!r} needs overlap_value.")
    if mode == "bp-ovr":
        amount: float = _gap_cli.integer(overlap_value, "overlap_value", minimum=1)
    elif mode == "range":
        amount = _gap_cli.integer(overlap_value, "overlap_value", minimum=0)
    else:
        amount = _gap_cli.number(overlap_value, "overlap_value")
        if not 0 < amount <= 1:
            raise ConfigurationError("overlap_value must be in (0, 1] for fraction modes.")
    return [*flags, f"--{mode}", str(amount)]


def _bedmap_inputs(reference_path: object, map_path: object) -> list[str]:
    reference = _input(reference_path, "reference_path")
    if map_path is None:
        return [str(reference)]
    return [str(reference), str(_input(map_path, "map_path"))]


def bedmap_score_statistic(
    reference_path: str,
    output_path: str,
    *,
    map_path: str | None = None,
    statistic: str = "mean",
    kth_fraction: float | None = None,
    mad_multiplier: float | None = None,
    trim_low: float | None = None,
    trim_high: float | None = None,
    overlap_mode: str | None = None,
    overlap_value: float | None = None,
    delimiter: str | None = None,
    precision: int | None = None,
    scientific: bool = False,
    skip_unmapped: bool = False,
    unmapped_value: str | None = None,
    chromosome: str | None = None,
    error_check: bool = False,
    header: bool = False,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """One numeric summary of the map elements overlapping each reference element (``bedmap``).

    ``statistic`` is the bedmap operation, one value per reference row: score statistics ``mean``,
    ``median``, ``sum``, ``min``, ``max``, ``variance``, ``stdev``, ``cv`` (stdev / mean), ``mad``
    (median absolute deviation, scaled by ``mad_multiplier``, default 1), ``kth`` (generalized median:
    the ``kth_fraction`` quantile, 0 < value <= 1, required), ``tmean`` (trimmed mean that ignores the
    bottom ``trim_low`` and top ``trim_high`` fractions, both required), ``wmean`` (mean weighted by the
    coverage of the reference element); and overlap measures ``count`` (number of overlapping map
    elements), ``bases`` (total overlapping bases), ``bases-uniq`` / ``bases-uniq-f`` (distinct
    reference bases covered, as a count / fraction) and ``indicator`` (1 when any element overlaps).
    The score statistics read the 5th column of ``map_path``. ``map_path`` omitted means the reference
    file is also the map file. ``overlap_mode`` is one of ``bp-ovr``, ``range``, ``exact``,
    ``fraction-ref``, ``fraction-map``, ``fraction-both``, ``fraction-either`` with ``overlap_value``
    (BEDOPS default: ``bp-ovr`` 1). Both files must be sorted as ``sort-bed`` sorts them. The result
    holds ``values`` (parsed first rows) and the output file with every row.
    """
    choice = _gap_cli.choice(statistic, SCORE_STATISTICS, "statistic")
    target = _output(output_path)
    operation: list[str] = [f"--{choice}"]
    given = {
        "kth_fraction": kth_fraction,
        "mad_multiplier": mad_multiplier,
        "trim_low": trim_low,
        "trim_high": trim_high,
    }
    allowed = {
        "kth": ("kth_fraction",),
        "mad": ("mad_multiplier",),
        "tmean": ("trim_low", "trim_high"),
    }
    for name, value in given.items():
        if value is not None and name not in allowed.get(choice, ()):
            raise ConfigurationError(f"{name} does not apply to statistic {choice!r}.")
    if choice == "kth":
        if kth_fraction is None:
            raise ConfigurationError("statistic 'kth' needs kth_fraction in (0, 1].")
        fraction = _gap_cli.number(kth_fraction, "kth_fraction")
        if not 0 < fraction <= 1:
            raise ConfigurationError("kth_fraction must be in (0, 1].")
        operation.append(str(fraction))
    elif choice == "mad" and mad_multiplier is not None:
        multiplier = _gap_cli.number(mad_multiplier, "mad_multiplier")
        if multiplier <= 0:
            raise ConfigurationError("mad_multiplier must be positive.")
        operation.append(str(multiplier))
    elif choice == "tmean":
        if trim_low is None or trim_high is None:
            raise ConfigurationError("statistic 'tmean' needs trim_low and trim_high.")
        low = _gap_cli.number(trim_low, "trim_low", minimum=0, maximum=1)
        high = _gap_cli.number(trim_high, "trim_high", minimum=0, maximum=1)
        if low + high > 1:
            raise ConfigurationError("trim_low + trim_high must not exceed 1.")
        operation += [str(low), str(high)]
    flags = _bedmap_flags(
        overlap_mode=overlap_mode,
        overlap_value=overlap_value,
        delimiter=delimiter,
        multi_delimiter=None,
        precision=precision,
        scientific=scientific,
        skip_unmapped=skip_unmapped,
        unmapped_value=unmapped_value,
        chromosome=chromosome,
        error_check=error_check,
        header=header,
    )
    run = _run(
        "bedmap",
        [*flags, *operation, *_bedmap_inputs(reference_path, map_path)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=target,
    )
    summary = _line_summary(target)
    values = [_number(line) for line in summary["preview"]]
    return _result(
        f"bedmap_{choice.replace('-', '_')}",
        run,
        outputs={"output": target},
        statistic=choice,
        values=values,
        **summary,
    )


def bedmap_echo_elements(
    reference_path: str,
    output_path: str,
    *,
    map_path: str | None = None,
    field: str = "echo-map",
    overlap_mode: str | None = None,
    overlap_value: float | None = None,
    delimiter: str | None = None,
    multi_delimiter: str | None = None,
    precision: int | None = None,
    scientific: bool = False,
    skip_unmapped: bool = False,
    unmapped_value: str | None = None,
    chromosome: str | None = None,
    error_check: bool = False,
    header: bool = False,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print, for each reference element, the reference row or the overlapping map elements (``bedmap`` echo operations).

    ``field`` is the bedmap operation: ``echo`` (the reference row), ``echo-map`` (all overlapping map
    elements), ``echo-map-id`` / ``echo-map-id-uniq`` (their IDs from column 4, all / unique),
    ``echo-map-range`` (the coordinate range they span), ``echo-map-score`` (column 5 scores),
    ``echo-map-size`` (their full lengths), ``echo-overlap-size`` (the lengths of the overlaps),
    ``echo-ref-name`` (``chrom:start-end``), ``echo-ref-row-id`` (``id-<line number>``), ``echo-ref-size``
    (reference element length), and the highest / lowest scoring overlapping element ``max-element``,
    ``min-element`` (``max-element-rand`` and ``min-element-rand`` choose randomly among ties).
    ``multi_delimiter`` separates multiple values within one column (BEDOPS default ``;``).
    Other parameters are as in ``bedmap_score_statistic``. Both files must be sorted as ``sort-bed``.
    """
    choice = _gap_cli.choice(field, ECHO_FIELDS, "field")
    target = _output(output_path)
    flags = _bedmap_flags(
        overlap_mode=overlap_mode,
        overlap_value=overlap_value,
        delimiter=delimiter,
        multi_delimiter=multi_delimiter,
        precision=precision,
        scientific=scientific,
        skip_unmapped=skip_unmapped,
        unmapped_value=unmapped_value,
        chromosome=chromosome,
        error_check=error_check,
        header=header,
    )
    run = _run(
        "bedmap",
        [*flags, f"--{choice}", *_bedmap_inputs(reference_path, map_path)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=target,
    )
    return _result(
        f"bedmap_{choice.replace('-', '_')}",
        run,
        outputs={"output": target},
        field=choice,
        **_line_summary(target),
    )


# ---------------------------------------------------------------------------
# Starch archives
# ---------------------------------------------------------------------------


def _compression(value: object) -> str:
    return _gap_cli.choice(value, ("bzip2", "gzip"), "compression")


def _note(value: object) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value or "\n" in value:
        raise ConfigurationError("note must be a non-empty single-line string.")
    return value


def _starch_elements(archive: Path, bin_dir: str | None) -> int | None:
    run = _run("unstarch", ["--elements", str(archive)], bin_dir=bin_dir, stdout=None, capture=True)
    text = _captured_text(run).strip()
    return int(text) if text.isdigit() else None


def starch_compress(
    input_path: str,
    output_path: str,
    *,
    compression: str = "bzip2",
    note: str | None = None,
    omit_signature: bool = False,
    header: bool = False,
    unique_tag: str | None = None,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compress a sorted BED file into a BEDOPS Starch archive (``starch``).

    ``compression`` is ``bzip2`` (smaller, the default) or ``gzip`` (extracts faster); ``note`` is
    stored in the archive metadata; ``omit_signature`` skips the per-chromosome integrity signature;
    ``header`` accepts UCSC track, SAM or VCF headers; ``unique_tag`` is the optional identifier of the
    transformed data. The BED input must be sorted as ``sort-bed`` sorts it.
    """
    source = _input(input_path)
    target = _output(output_path)
    arguments = [f"--{_compression(compression)}"]
    if note is not None:
        arguments.append(f"--note={_note(note)}")
    if _flag(omit_signature, "omit_signature"):
        arguments.append("--omit-signature")
    if _flag(header, "header"):
        arguments.append("--header")
    if unique_tag is not None:
        arguments.append(_gap_cli.text(unique_tag, "unique_tag"))
    bin_dir = _bin_dir(bedops_bin_dir)
    run = _run("starch", [*arguments, str(source)], bin_dir=bin_dir, stdout=target)
    return _result(
        "starch_compress",
        run,
        outputs={"archive": target},
        element_count=_starch_elements(target, bin_dir),
        size_bytes=target.stat().st_size,
    )


def unstarch_to_bed(
    input_path: str,
    output_path: str,
    *,
    chromosome: str | None = None,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decompress a Starch archive to BED, optionally one chromosome only (``unstarch [<chromosome>]``)."""
    source = _input(input_path)
    target = _output(output_path)
    arguments = [] if chromosome is None else [_chromosome(chromosome)]
    run = _run(
        "unstarch",
        [*arguments, str(source)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=target,
    )
    return _result("unstarch_to_bed", run, outputs={"output": target}, **_line_summary(target))


def unstarch_summary(
    input_path: str,
    *,
    statistic: str = "elements",
    chromosome: str | None = None,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one metadata value from a Starch archive without decompressing it (``unstarch``).

    ``statistic`` selects the value: ``elements`` (element count), ``bases`` (total bases),
    ``bases-uniq`` (distinct bases), ``has-duplicate`` (True when a duplicate element exists),
    ``has-nested`` (True when a nested element exists), ``list-chromosomes`` (chromosome names),
    ``elements-max-string-length`` (longest element line) and ``list-json`` (the full archive metadata
    as a dictionary). With ``chromosome`` the value is restricted to that chromosome. The value is
    returned as ``value``.
    """
    choice = _gap_cli.choice(statistic, UNSTARCH_STATISTICS, "statistic")
    source = _input(input_path)
    flag = {"list-chromosomes": "--list-chr"}.get(choice, f"--{choice}")
    arguments = [] if chromosome is None else [_chromosome(chromosome)]
    run = _run(
        "unstarch",
        [*arguments, flag, str(source)],
        bin_dir=_bin_dir(bedops_bin_dir),
        stdout=None,
        capture=True,
    )
    text = _captured_text(run).strip()
    value: object
    if choice == "list-chromosomes":
        value = [line for line in text.splitlines() if line]
    elif choice in ("has-duplicate", "has-nested"):
        value = text == "1"
    elif choice == "list-json":
        import json

        try:
            value = json.loads(text)
        except ValueError as exc:
            raise BackendExecutionError(
                "unstarch --list-json did not return JSON.",
                code="BEDOPS_OUTPUT_UNREADABLE",
                context={"output_head": text[:200]},
            ) from exc
    else:
        value = _number(text)
    return _result("unstarch_summary", run, statistic=choice, chromosome=chromosome, value=value)


def starchcat_archives(
    input_paths: Sequence[str],
    output_path: str,
    *,
    compression: str = "bzip2",
    note: str | None = None,
    omit_signature: bool = False,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge Starch archives into one multiset union archive (``starchcat``).

    One input only re-compresses the archive with another backend or updates its note. ``compression``,
    ``note`` and ``omit_signature`` are as in ``starch_compress``. The inputs must be sorted,
    header-less archives.
    """
    inputs = _inputs(input_paths, "input_paths", minimum=1)
    target = _output(output_path)
    arguments = [f"--{_compression(compression)}"]
    if note is not None:
        arguments.append(f"--note={_note(note)}")
    if _flag(omit_signature, "omit_signature"):
        arguments.append("--omit-signature")
    bin_dir = _bin_dir(bedops_bin_dir)
    run = _run(
        "starchcat",
        [*arguments, *(str(path) for path in inputs)],
        bin_dir=bin_dir,
        stdout=target,
    )
    return _result(
        "starchcat_archives",
        run,
        outputs={"archive": target},
        element_count=_starch_elements(target, bin_dir),
        size_bytes=target.stat().st_size,
    )


def starchstrip_chromosomes(
    input_path: str,
    output_path: str,
    chromosomes: Sequence[str],
    *,
    mode: str = "include",
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep (``mode="include"``) or drop (``mode="exclude"``) chromosomes of a Starch archive (``starchstrip``).

    The selected chromosome streams are copied without recompression; chromosomes missing from the
    archive are ignored. BEDOPS 2.4.42 ``starchstrip --exclude`` writes an archive whose data ``unstarch``
    cannot read when the first chromosome is removed, so for ``mode="exclude"`` the facade lists the
    archive's chromosomes and passes the complementary ``--include`` list, which gives the intended
    archive (the result reports this in ``strategy``).
    """
    choice = _gap_cli.choice(mode, ("include", "exclude"), "mode")
    source = _input(input_path)
    target = _output(output_path)
    if isinstance(chromosomes, str) or not isinstance(chromosomes, Sequence) or not chromosomes:
        raise ConfigurationError("chromosomes must be a non-empty list of names.")
    names = [_chromosome(name, "chromosomes") for name in chromosomes]
    if any("," in name for name in names):
        raise ConfigurationError("chromosome names must not contain commas.")
    bin_dir = _bin_dir(bedops_bin_dir)
    selected = names
    strategy = "starchstrip --include"
    if choice == "exclude":
        listing = _run(
            "unstarch", ["--list-chr", str(source)], bin_dir=bin_dir, stdout=None, capture=True
        )
        present = [line for line in _captured_text(listing).splitlines() if line]
        selected = [name for name in present if name not in names]
        if not selected:
            raise ConfigurationError("Excluding these chromosomes would leave an empty archive.")
        strategy = "starchstrip --include of the complementary chromosome list (see docstring)"
    run = _run(
        "starchstrip",
        ["--include", ",".join(selected), str(source)],
        bin_dir=bin_dir,
        stdout=target,
    )
    chromosome_run = _run(
        "unstarch", ["--list-chr", str(target)], bin_dir=bin_dir, stdout=None, capture=True
    )
    kept = [line for line in _captured_text(chromosome_run).splitlines() if line]
    return _result(
        "starchstrip_chromosomes",
        run,
        outputs={"archive": target},
        mode=choice,
        strategy=strategy,
        chromosomes_kept=kept,
        element_count=_starch_elements(target, bin_dir),
    )


# ---------------------------------------------------------------------------
# convert2bed
# ---------------------------------------------------------------------------


def convert_to_bed(
    input_path: str,
    output_path: str,
    *,
    input_format: str,
    output_format: str = "bed",
    format_options: Sequence[str] | None = None,
    attribute_key: str | None = None,
    sort_output: bool = True,
    max_memory: str | None = None,
    temp_dir: str | None = None,
    starch_compression: str = "bzip2",
    starch_note: str | None = None,
    bedops_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a genomic file to sorted BED or to a Starch archive (``convert2bed``).

    ``input_format`` is one of ``bam``, ``sam``, ``gff`` (GFF3), ``gtf``, ``gvf``, ``psl``, ``rmsk``
    (RepeatMasker .out), ``vcf`` or ``wig``; ``output_format`` is ``bed`` or ``starch``. This one function
    covers the single-format programs ``bam2bed``, ``sam2bed``, ``gff2bed``, ``gtf2bed``, ``gvf2bed``,
    ``psl2bed``, ``rmsk2bed``, ``vcf2bed``, ``wig2bed`` and their ``*2starch`` variants (the cluster
    wrappers ``bam2starch_sge`` / ``_slurm`` / ``_gnuParallel`` only distribute the same conversion).
    ``format_options`` are boolean convert2bed switches valid for the input format, for example
    ``keep-header`` (all), ``split``, ``all-reads``, ``reduced`` (sam/bam), ``snvs``, ``insertions``,
    ``deletions``, ``do-not-split`` (vcf), ``zero-indexed`` (wig); ``attribute_key`` selects the
    GENCODE attribute (``gene_id`` by default, ``transcript_id``, ``gene_name``, ...) copied to the BED
    name column for ``gff`` / ``gtf``. ``sort_output=False`` skips the final
    sort-bed (BED output only). ``sam`` and ``bam`` input needs ``samtools`` on PATH.
    """
    fmt = _gap_cli.choice(input_format, CONVERT_INPUT_FORMATS, "input_format")
    out_fmt = _gap_cli.choice(output_format, CONVERT_OUTPUT_FORMATS, "output_format")
    source = _input(input_path)
    target = _output(output_path)
    arguments = [f"--input={fmt}", f"--output={out_fmt}"]
    if format_options is not None:
        if isinstance(format_options, str) or not isinstance(format_options, Sequence):
            raise ConfigurationError("format_options must be a list of option names.")
        for option_name in format_options:
            name = str(option_name).removeprefix("--")
            if name not in _CONVERT_FLAGS[fmt]:
                raise ConfigurationError(
                    f"format option {option_name!r} is not allowed for {fmt}; allowed: "
                    f"{sorted(_CONVERT_FLAGS[fmt])}."
                )
            arguments.append(f"--{name}")
    if attribute_key is not None:
        if fmt not in ("gff", "gtf"):
            raise ConfigurationError("attribute_key only applies to gff and gtf input.")
        key = _gap_cli.choice(attribute_key, _ATTRIBUTE_KEYS[fmt], "attribute_key")
        arguments.append(f"--attribute-key={key}")
    if not _flag(sort_output, "sort_output"):
        if out_fmt == "starch":
            raise ConfigurationError(
                "sort_output=False cannot be combined with output_format='starch'."
            )
        arguments.append("--do-not-sort")
    if max_memory is not None:
        if not isinstance(max_memory, str) or not _MEMORY.match(max_memory):
            raise ConfigurationError("max_memory must look like '8G', '8000M' or '8000000000'.")
        arguments.append(f"--max-mem={max_memory}")
    if temp_dir is not None:
        folder = Path(str(temp_dir)).expanduser()
        if not folder.is_dir():
            raise ConfigurationError(f"temp_dir must be an existing folder: {folder}")
        arguments.append(f"--sort-tmpdir={folder.resolve()}")
    if out_fmt == "starch":
        arguments.append(f"--starch-{_compression(starch_compression)}")
        if starch_note is not None:
            arguments.append(f"--starch-note={_note(starch_note)}")
    elif starch_note is not None or starch_compression != "bzip2":
        raise ConfigurationError(
            "starch_compression and starch_note only apply to output_format='starch'."
        )
    if fmt in ("sam", "bam") and shutil.which("samtools") is None:
        raise BackendUnavailableError("samtools was not found on PATH.", hint=_SAMTOOLS_HINT)
    bin_dir = _bin_dir(bedops_bin_dir)
    run = _run("convert2bed", arguments, bin_dir=bin_dir, stdout=target, stdin=source)
    if out_fmt == "starch":
        return _result(
            f"convert_{fmt}_to_starch",
            run,
            outputs={"archive": target},
            input_format=fmt,
            output_format=out_fmt,
            element_count=_starch_elements(target, bin_dir),
            size_bytes=target.stat().st_size,
        )
    return _result(
        f"convert_{fmt}_to_bed",
        run,
        outputs={"output": target},
        input_format=fmt,
        output_format=out_fmt,
        **_line_summary(target),
    )


__all__ = [
    "bedextract_chromosome",
    "bedextract_list_chromosomes",
    "bedmap_echo_elements",
    "bedmap_score_statistic",
    "bedops_closest_features",
    "bedops_set_operation",
    "bedops_sort_bed_filter",
    "convert_to_bed",
    "starch_compress",
    "starchcat_archives",
    "starchstrip_chromosomes",
    "unstarch_summary",
    "unstarch_to_bed",
]

_register_public_verbose(globals(), tuple(__all__))
