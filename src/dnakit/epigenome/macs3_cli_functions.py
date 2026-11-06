"""MACS3 command-line facades: fixed option sets, no shell, bounded execution.

Each public function runs one MACS3 subcommand (3.0.5 is the tested version) through
``dnakit.backends.external.execute_bounded_command``. The user installs MACS3 (``pip install macs3``
or conda) and either puts ``macs3`` on PATH or passes ``macs3_path``. Options are explicit keyword
arguments that map onto a fixed flag list, so no other MACS3 option can be passed. Outputs are
written to explicit paths or to a new output directory; standard output is never the only result.
Results are JSON-serialisable dictionaries with status, backend, version, command, outputs, parsed
counts and a bounded log tail.
"""

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_PROGRAM = "macs3"
_TESTED_VERSION = "3.0.5"
_CITATION = "https://doi.org/10.1186/gb-2008-9-9-r137"
_INSTALL_HINT = (
    "Install MACS3 3.0.5 or later (python -m pip install macs3, or conda install -c bioconda macs3), "
    "then put macs3 on PATH or pass macs3_path."
)
_TIMEOUT_SECONDS = 3600.0
_VERSION_TIMEOUT_SECONDS = 60.0
_MAX_OUTPUT_BYTES = 5_000_000
_MAX_ARTIFACT_BYTES = 2_000_000_000
_LOG_TAIL_CHARS = 4000
_VERSION_LINE = re.compile(r"macs3\s+(\S+)")
_FRAGMENT_LENGTH = re.compile(r"predicted fragment length is (\d+) bps")
_TAG_SIZE = re.compile(r"tag size is determined as (\d+) bps")
_FORMATS = (
    "AUTO",
    "BAM",
    "SAM",
    "BED",
    "ELAND",
    "ELANDMULTI",
    "ELANDEXPORT",
    "BOWTIE",
    "BAMPE",
    "BEDPE",
    "FRAG",
)
_PEAK_FORMATS = (
    "AUTO",
    "BAM",
    "SAM",
    "BED",
    "ELAND",
    "ELANDMULTI",
    "ELANDEXPORT",
    "BOWTIE",
    "BAMPE",
    "BEDPE",
)
_BDG_METHODS = ("ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR", "max")
_OPT_METHODS = ("multiply", "add", "p2q", "max", "min")
_COMBINE_METHODS = ("fisher", "max", "mean")


# ---------------------------------------------------------------------------
# Argument checks
# ---------------------------------------------------------------------------


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def _int(
    value: object, name: str, *, minimum: int | None = None, maximum: int | None = None
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def _number(
    value: object, name: str, *, minimum: float | None = None, maximum: float | None = None
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value != value:
        raise ConfigurationError(f"{name} must be a number.")
    number = float(value)
    if minimum is not None and number < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and number > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _input(value: object, name: str = "input_path") -> Path:
    path = Path(_text(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}")
    return path


def _inputs(values: object, name: str) -> tuple[Path, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of file paths.")
    return tuple(_input(value, name) for value in values)


def _output(value: object, name: str = "output_path") -> Path:
    """A new file: its parent must exist and the file itself must not, so nothing is overwritten."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _new_directory(value: object, name: str = "outdir") -> Path:
    """A directory that MACS3 creates; it must not exist yet."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output directory already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Parent of the output directory does not exist: {path.parent}")
    return path


def _name(value: object, name: str = "name") -> str:
    text = _text(value, name)
    if "/" in text or "\\" in text or text.startswith("-") or len(text) > 200:
        raise ConfigurationError(
            f"{name} must be a plain name without slashes and not starting with '-'."
        )
    return text


def _genome_size(value: object) -> str:
    """MACS3 accepts the shortcuts hs, mm, ce, dm, or an integer number of effective bases."""
    if isinstance(value, bool):
        raise ConfigurationError(
            "genome_size must be a shortcut (hs, mm, ce, dm) or a positive integer."
        )
    if isinstance(value, int):
        return str(_int(value, "genome_size", minimum=1))
    if isinstance(value, str) and value in {"hs", "mm", "ce", "dm"}:
        return value
    if isinstance(value, str) and value.isdigit() and int(value) > 0:
        return value
    raise ConfigurationError("genome_size must be hs, mm, ce, dm or a positive integer.")


def _mfold(value: object) -> tuple[int, int]:
    if isinstance(value, str) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError("mfold must be two integers (lower and upper fold enrichment).")
    low = _int(value[0], "mfold lower", minimum=1)
    high = _int(value[1], "mfold upper", minimum=low)
    return low, high


def _keep_duplicates(value: object) -> str:
    if isinstance(value, bool):
        raise ConfigurationError("keep_duplicates must be 'auto', 'all' or a positive integer.")
    if isinstance(value, int):
        return str(_int(value, "keep_duplicates", minimum=1))
    if isinstance(value, str) and value in {"auto", "all"}:
        return value
    raise ConfigurationError("keep_duplicates must be 'auto', 'all' or a positive integer.")


def _floats(values: object, name: str) -> tuple[float, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of numbers.")
    return tuple(_number(value, name) for value in values)


def _opt(flag: str, value: object) -> list[str]:
    """One option with its value, or nothing when the value is None."""
    return [] if value is None else [flag, str(value)]


def _switch(flag: str, enabled: bool) -> list[str]:
    return [flag] if enabled else []


# ---------------------------------------------------------------------------
# Executable, execution and results
# ---------------------------------------------------------------------------


def _resolve_executable(macs3_path: str | None) -> str:
    """Find macs3 on PATH, or use an explicit path. Never a hard-coded location."""
    if macs3_path is None or macs3_path == "":
        found = shutil.which(_PROGRAM)
        if found is None:
            raise BackendUnavailableError(
                "External executable 'macs3' was not found on PATH.",
                hint=_INSTALL_HINT,
                code="MACS3_UNAVAILABLE",
            )
        return str(Path(found).absolute())
    if os.sep not in macs3_path:
        found = shutil.which(macs3_path)
        if found is None:
            raise BackendUnavailableError(
                f"External executable {macs3_path!r} was not found on PATH.",
                hint=_INSTALL_HINT,
                code="MACS3_UNAVAILABLE",
            )
        return str(Path(found).absolute())
    candidate = Path(macs3_path).expanduser()
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise BackendUnavailableError(
            f"MACS3 executable is missing or not executable: {macs3_path}",
            hint=_INSTALL_HINT,
            code="MACS3_UNAVAILABLE",
        )
    return str(candidate.absolute())


@dataclass(frozen=True)
class _Completed:
    text: str
    version: str | None
    executable: str
    command: tuple[str, ...]


def _version(executable: str, cwd: Path, verbose: bool) -> str | None:
    result = execute_bounded_command(
        executable,
        ("--version",),
        backend_id=_PROGRAM,
        cwd=cwd,
        timeout_seconds=_VERSION_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_OUTPUT_BYTES,
        verbose=verbose,
    )
    match = _VERSION_LINE.search(result.output)
    return match.group(1) if match else None


def _run(
    subcommand: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    monitored: Sequence[Path],
    macs3_path: str | None,
    verbose: bool,
) -> _Completed:
    """Run one MACS3 subcommand with a fixed argument list. No shell."""
    executable = _resolve_executable(macs3_path)
    version = _version(executable, cwd, verbose)
    argv = (subcommand, *arguments)
    result = execute_bounded_command(
        executable,
        argv,
        backend_id=_PROGRAM,
        cwd=cwd,
        timeout_seconds=_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_OUTPUT_BYTES,
        monitored_output_paths=tuple(monitored),
        max_monitored_output_bytes=_MAX_ARTIFACT_BYTES,
        # MACS3 iterates over Python sets and dicts, so its output order (and a seeded sample) changes with
        # the process hash seed. A fixed seed keeps the order of every output the same from run to run.
        environment_overrides={"PYTHONHASHSEED": "0"},
        verbose=verbose,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"MACS3 {subcommand} exited with code {result.return_code}.",
            code="MACS3_EXECUTION_FAILED",
            context={"subcommand": subcommand, "log_tail": _tail(result.output)},
        )
    missing = [str(path) for path in monitored if not path.is_file()]
    if missing:
        raise BackendExecutionError(
            f"MACS3 {subcommand} did not write the expected files: {missing}",
            code="MACS3_OUTPUT_MISSING",
            context={"subcommand": subcommand},
        )
    return _Completed(result.output, version, executable, (_PROGRAM, *argv))


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


def _count_rows(path: Path) -> int:
    """Data rows of a BED, bedGraph, narrowPeak or VCF file (track, header and comment lines skipped)."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(
            1 for line in handle if line.strip() and not line.startswith(("#", "track", "browser"))
        )


def _directory_files(directory: Path) -> dict[str, Path]:
    return {path.name: path for path in sorted(directory.iterdir()) if path.is_file()}


def _result(
    subcommand: str,
    run: _Completed,
    *,
    outputs: Mapping[str, Path],
    counts: Mapping[str, object] | None = None,
    **extra: object,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": _PROGRAM,
        "operation": subcommand,
        "version": run.version,
        "command": list(run.command),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "counts": dict(counts or {}),
        "log_tail": _tail(run.text),
        "provenance": {
            "tool": "MACS3",
            "tested_version": _TESTED_VERSION,
            "citation": _CITATION,
            "executable": run.executable,
        },
    }
    result.update(extra)
    return result


def _fragment_info(text: str) -> dict[str, int | None]:
    """Tag size and predicted fragment length as MACS3 logs them (None when the model was not built)."""
    tag = _TAG_SIZE.findall(text)
    fragment = _FRAGMENT_LENGTH.findall(text)
    return {
        "tag_size": int(float(tag[-1])) if tag else None,
        "predicted_fragment_length": int(fragment[-1]) if fragment else None,
    }


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def macs3_call_peaks(
    treatment_paths: Sequence[str],
    outdir: str,
    name: str,
    *,
    control_paths: Sequence[str] | None = None,
    input_format: str = "AUTO",
    genome_size: str | int = "hs",
    tag_size: int | None = None,
    qvalue: float = 0.05,
    pvalue: float | None = None,
    keep_duplicates: str | int = 1,
    extsize: int | None = None,
    shift: int | None = None,
    no_model: bool = False,
    mfold: Sequence[int] = (5, 50),
    bandwidth: int = 300,
    call_summits: bool = False,
    broad: bool = False,
    broad_cutoff: float = 0.1,
    min_length: int | None = None,
    max_gap: int | None = None,
    save_bedgraph: bool = False,
    seed: int | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call peaks from alignment files (``macs3 callpeak``): narrowPeak or broadPeak regions.

    Writes ``<name>_peaks.narrowPeak`` (or ``broadPeak`` with ``broad``), ``<name>_peaks.xls`` and,
    when requested, the summit BED and bedGraph tracks, into the new directory ``outdir``. Give
    ``pvalue`` to use a p-value cutoff instead of ``qvalue``. ``no_model`` with ``extsize`` skips
    the shift model, which is needed when MACS3 cannot find enough paired peaks to build it.
    """
    treatment = _inputs(treatment_paths, "treatment_paths")
    controls = _inputs(control_paths, "control_paths") if control_paths else ()
    directory = _new_directory(outdir)
    label = _name(name)
    fmt = _choice(input_format, _FORMATS, "input_format")
    options: list[str] = ["-t", *map(str, treatment)]
    if controls:
        options += ["-c", *map(str, controls)]
    options += ["-f", fmt, "-g", _genome_size(genome_size), "-n", label, "--outdir", str(directory)]
    if pvalue is None:
        options += ["-q", str(_number(qvalue, "qvalue", minimum=1e-12, maximum=1.0))]
    else:
        options += ["-p", str(_number(pvalue, "pvalue", minimum=1e-12, maximum=1.0))]
    options += _opt("-s", None if tag_size is None else _int(tag_size, "tag_size", minimum=1))
    options += ["--keep-dup", _keep_duplicates(keep_duplicates)]
    options += _opt("--extsize", None if extsize is None else _int(extsize, "extsize", minimum=1))
    options += _opt(
        "--shift", None if shift is None else _int(shift, "shift", minimum=-100000, maximum=100000)
    )
    options += _switch("--nomodel", _flag(no_model, "no_model"))
    low, high = _mfold(mfold)
    options += ["-m", str(low), str(high), "--bw", str(_int(bandwidth, "bandwidth", minimum=1))]
    options += _switch("--call-summits", _flag(call_summits, "call_summits"))
    options += _switch("--broad", _flag(broad, "broad"))
    if broad:
        options += [
            "--broad-cutoff",
            str(_number(broad_cutoff, "broad_cutoff", minimum=1e-12, maximum=1.0)),
        ]
    options += _opt(
        "--min-length", None if min_length is None else _int(min_length, "min_length", minimum=1)
    )
    options += _opt("--max-gap", None if max_gap is None else _int(max_gap, "max_gap", minimum=0))
    options += _switch("-B", _flag(save_bedgraph, "save_bedgraph"))
    options += _opt("--seed", None if seed is None else _int(seed, "seed", minimum=0))

    peak_name = "broadPeak" if broad else "narrowPeak"
    monitored = [directory / f"{label}_peaks.{peak_name}", directory / f"{label}_peaks.xls"]
    if call_summits:
        monitored.append(directory / f"{label}_summits.bed")
    if save_bedgraph:
        monitored += [
            directory / f"{label}_treat_pileup.bdg",
            directory / f"{label}_control_lambda.bdg",
        ]
    run = _run(
        "callpeak",
        options,
        cwd=directory.parent,
        monitored=monitored,
        macs3_path=macs3_path,
        verbose=verbose,
    )
    files = _directory_files(directory)
    counts = {"peak_count": _count_rows(directory / f"{label}_peaks.{peak_name}")}
    if call_summits:
        counts["summit_count"] = _count_rows(directory / f"{label}_summits.bed")
    info = _fragment_info(run.text)
    if info["predicted_fragment_length"] is None and no_model and extsize is not None:
        info["predicted_fragment_length"] = int(extsize)
    return _result(
        "callpeak",
        run,
        outputs={name: directory / name for name in files},
        counts=counts,
        peak_format=peak_name,
        model=info,
        mode="nomodel" if no_model else "model",
    )


def macs3_predict_fragment_length(
    input_paths: Sequence[str],
    outdir: str,
    *,
    input_format: str = "AUTO",
    genome_size: str | int = "hs",
    mfold: Sequence[int] = (5, 50),
    bandwidth: int = 300,
    d_min: int | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the fragment length from alignments (``macs3 predictd``) and write its R model script.

    The predicted length is parsed from MACS3's log. MACS3 exits with status 0 when it cannot build the model
    (fewer than 100 paired peaks on the two strands), so this facade raises BackendExecutionError in that case
    instead of returning a result without a length.
    """
    inputs = _inputs(input_paths, "input_paths")
    directory = _new_directory(outdir)
    fmt = _choice(input_format, _PEAK_FORMATS, "input_format")
    low, high = _mfold(mfold)
    options = [
        "-i",
        *map(str, inputs),
        "-f",
        fmt,
        "-g",
        _genome_size(genome_size),
        "-m",
        str(low),
        str(high),
    ]
    options += ["--bw", str(_int(bandwidth, "bandwidth", minimum=1)), "--outdir", str(directory)]
    options += _opt("--d-min", None if d_min is None else _int(d_min, "d_min", minimum=1))
    run = _run(
        "predictd",
        options,
        cwd=directory.parent,
        monitored=(),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    info = _fragment_info(run.text)
    if info["predicted_fragment_length"] is None:
        raise BackendExecutionError(
            "MACS3 predictd did not build a fragment model; widen mfold or use a larger sample.",
            code="MACS3_MODEL_NOT_BUILT",
            context={"subcommand": "predictd", "log_tail": _tail(run.text)},
        )
    files = _directory_files(directory)
    return _result(
        "predictd",
        run,
        outputs={name: directory / name for name in files},
        counts={"file_count": len(files)},
        model=info,
    )


def macs3_pileup(
    input_paths: Sequence[str],
    output_path: str,
    *,
    extsize: int,
    input_format: str = "AUTO",
    both_directions: bool = False,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up single-end reads (or paired-end fragments) into a bedGraph, extended to ``extsize`` bp.

    ``macs3 pileup`` with ``--extsize``: each read is extended from its 5' end to ``extsize`` bases
    downstream (both directions with ``both_directions``).
    """
    inputs = _inputs(input_paths, "input_paths")
    target = _output(output_path)
    fmt = _choice(input_format, _FORMATS, "input_format")
    options = [
        "-i",
        *map(str, inputs),
        "-f",
        fmt,
        "--extsize",
        str(_int(extsize, "extsize", minimum=1)),
    ]
    options += _switch("-B", _flag(both_directions, "both_directions"))
    options += ["-o", str(target)]
    run = _run(
        "pileup",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "pileup", run, outputs={"bedgraph": target}, counts={"bedgraph_rows": _count_rows(target)}
    )


def macs3_filter_duplicates(
    input_paths: Sequence[str],
    output_path: str,
    *,
    genome_size: str | int = "hs",
    input_format: str = "AUTO",
    keep_duplicates: str | int = "auto",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove duplicate reads at the same position and strand, and write a BED file (``macs3 filterdup``).

    ``keep_duplicates`` is MACS3's own default ``auto`` (a binomial test at p = 1e-5 sets the number of
    tags kept per position), ``all`` keeps every read, and an integer keeps at most that many.
    """
    inputs = _inputs(input_paths, "input_paths")
    target = _output(output_path)
    fmt = _choice(input_format, _PEAK_FORMATS, "input_format")
    options = [
        "-i",
        *map(str, inputs),
        "-f",
        fmt,
        "-g",
        _genome_size(genome_size),
        "--keep-dup",
        _keep_duplicates(keep_duplicates),
        "-o",
        str(target),
    ]
    run = _run(
        "filterdup",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "filterdup", run, outputs={"bed": target}, counts={"kept_reads": _count_rows(target)}
    )


def macs3_random_sample(
    input_paths: Sequence[str],
    output_path: str,
    *,
    percentage: float | None = None,
    number: int | None = None,
    input_format: str = "AUTO",
    seed: int | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly keep a percentage or a number of reads or fragments and write them as BED/BEDPE."""
    if (percentage is None) == (number is None):
        raise ConfigurationError("Give exactly one of percentage or number.")
    inputs = _inputs(input_paths, "input_paths")
    target = _output(output_path)
    fmt = _choice(input_format, _PEAK_FORMATS, "input_format")
    options = ["-i", *map(str, inputs), "-f", fmt, "-o", str(target)]
    if percentage is not None:
        options += ["-p", str(_number(percentage, "percentage", minimum=0.0, maximum=100.0))]
    else:
        options += ["-n", str(_int(number, "number", minimum=1))]
    options += _opt("--seed", None if seed is None else _int(seed, "seed", minimum=0))
    run = _run(
        "randsample",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "randsample", run, outputs={"sample": target}, counts={"sampled_reads": _count_rows(target)}
    )


def macs3_refine_summits(
    candidate_peaks: str,
    input_paths: Sequence[str],
    output_path: str,
    *,
    input_format: str = "BAM",
    cutoff: float = 0.0,
    window_size: int = 100,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Move each candidate peak's summit to the position supported by the 5' read positions (``refinepeak``)."""
    candidate = _input(candidate_peaks, "candidate_peaks")
    inputs = _inputs(input_paths, "input_paths")
    target = _output(output_path)
    fmt = _choice(input_format, _PEAK_FORMATS, "input_format")
    options = [
        "-b",
        str(candidate),
        "-i",
        *map(str, inputs),
        "-f",
        fmt,
        "-c",
        str(_number(cutoff, "cutoff", minimum=0.0)),
        "-w",
        str(_int(window_size, "window_size", minimum=1)),
        "-o",
        str(target),
    ]
    run = _run(
        "refinepeak",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "refinepeak", run, outputs={"bed": target}, counts={"refined_peaks": _count_rows(target)}
    )


def macs3_call_variants(
    peak_bed: str,
    treatment_bam: str,
    output_path: str,
    *,
    control_bam: str | None = None,
    min_quality: int = 30,
    max_duplicates: int = 1,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call genotypes at SNP sites inside peaks from allele-specific read counts (``callvar``).

    ``min_quality`` (-Q) and ``max_duplicates`` (-D) are the base-quality and duplicate filters of MACS3.
    """
    peaks = _input(peak_bed, "peak_bed")
    treatment = _input(treatment_bam, "treatment_bam")
    target = _output(output_path)
    options = [
        "-b",
        str(peaks),
        "-t",
        str(treatment),
        "-o",
        str(target),
        "-Q",
        str(_int(min_quality, "min_quality", minimum=0)),
        "-D",
        str(_int(max_duplicates, "max_duplicates", minimum=1)),
    ]
    if control_bam is not None:
        options += ["-c", str(_input(control_bam, "control_bam"))]
    run = _run(
        "callvar",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "callvar", run, outputs={"variants": target}, counts={"variant_rows": _count_rows(target)}
    )


def macs3_combine_replicates(
    bedgraph_paths: Sequence[str],
    output_path: str,
    *,
    method: str = "fisher",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the MACS3 score bedGraphs of replicates with Fisher, max or mean (``cmbreps``)."""
    inputs = _inputs(bedgraph_paths, "bedgraph_paths")
    target = _output(output_path)
    chosen = _choice(method, _COMBINE_METHODS, "method")
    options = ["-i", *map(str, inputs), "-m", chosen, "-o", str(target)]
    run = _run(
        "cmbreps",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "cmbreps", run, outputs={"bedgraph": target}, counts={"bedgraph_rows": _count_rows(target)}
    )


def macs3_compare_bedgraphs(
    treatment_bedgraph: str,
    control_bedgraph: str,
    output_path: str,
    *,
    method: str = "ppois",
    scaling_factor: float = 1.0,
    pseudocount: float = 0.0,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score each bin of a treatment bedGraph against a control bedGraph (``bdgcmp``).

    ``method`` is one of ppois, qpois, subtract, logFE, FE, logLR, slogLR or max, as in MACS3.
    """
    treatment = _input(treatment_bedgraph, "treatment_bedgraph")
    control = _input(control_bedgraph, "control_bedgraph")
    target = _output(output_path)
    chosen = _choice(method, _BDG_METHODS, "method")
    options = [
        "-t",
        str(treatment),
        "-c",
        str(control),
        "-m",
        chosen,
        "-S",
        str(_number(scaling_factor, "scaling_factor", minimum=1e-12)),
        "-p",
        str(_number(pseudocount, "pseudocount", minimum=0.0)),
        "-o",
        str(target),
    ]
    run = _run(
        "bdgcmp",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "bdgcmp", run, outputs={"bedgraph": target}, counts={"bedgraph_rows": _count_rows(target)}
    )


def macs3_operate_bedgraph(
    input_bedgraph: str,
    output_path: str,
    *,
    method: str,
    extra_params: Sequence[float] = (),
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multiply, add, take -log10 p to q, or take the max or min of each bedGraph score (``bdgopt``)."""
    source = _input(input_bedgraph, "input_bedgraph")
    target = _output(output_path)
    chosen = _choice(method, _OPT_METHODS, "method")
    params = _floats(extra_params, "extra_params")
    if chosen != "p2q" and len(params) != 1:
        raise ConfigurationError(
            f"method {chosen!r} needs exactly one value in extra_params (MACS3 -p)."
        )
    if chosen == "p2q" and params:
        raise ConfigurationError("method 'p2q' takes no extra_params.")
    options = ["-i", str(source), "-m", chosen, "-o", str(target)]
    if params:
        options += ["-p", *map(str, params)]
    run = _run(
        "bdgopt",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "bdgopt", run, outputs={"bedgraph": target}, counts={"bedgraph_rows": _count_rows(target)}
    )


def macs3_call_bedgraph_peaks(
    input_bedgraph: str,
    output_path: str,
    *,
    cutoff: float,
    min_length: int | None = None,
    max_gap: int | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call peaks from a score bedGraph: regions above ``cutoff`` (``bdgpeakcall``)."""
    source = _input(input_bedgraph, "input_bedgraph")
    target = _output(output_path)
    options = ["-i", str(source), "-c", str(_number(cutoff, "cutoff")), "-o", str(target)]
    options += _opt("-l", None if min_length is None else _int(min_length, "min_length", minimum=1))
    options += _opt("-g", None if max_gap is None else _int(max_gap, "max_gap", minimum=0))
    run = _run(
        "bdgpeakcall",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "bdgpeakcall", run, outputs={"bed": target}, counts={"peak_count": _count_rows(target)}
    )


def macs3_call_bedgraph_broad_peaks(
    input_bedgraph: str,
    output_path: str,
    *,
    cutoff_peak: float,
    cutoff_link: float,
    min_length: int = 200,
    lvl1_max_gap: int = 50,
    lvl2_max_gap: int = 400,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call nested broad peaks from a score bedGraph with a strong and a weak cutoff (``bdgbroadcall``)."""
    source = _input(input_bedgraph, "input_bedgraph")
    target = _output(output_path)
    peak = _number(cutoff_peak, "cutoff_peak")
    link = _number(cutoff_link, "cutoff_link")
    if link > peak:
        raise ConfigurationError("cutoff_link must not be larger than cutoff_peak.")
    options = [
        "-i",
        str(source),
        "-c",
        str(peak),
        "-C",
        str(link),
        "-l",
        str(_int(min_length, "min_length", minimum=1)),
        "-g",
        str(_int(lvl1_max_gap, "lvl1_max_gap", minimum=0)),
        "-G",
        str(_int(lvl2_max_gap, "lvl2_max_gap", minimum=0)),
        "-o",
        str(target),
    ]
    run = _run(
        "bdgbroadcall",
        options,
        cwd=target.parent,
        monitored=(target,),
        macs3_path=macs3_path,
        verbose=verbose,
    )
    return _result(
        "bdgbroadcall",
        run,
        outputs={"bed": target},
        counts={"broad_region_count": _count_rows(target)},
    )


def macs3_diff_bedgraphs(
    treatment1_bedgraph: str,
    control1_bedgraph: str,
    treatment2_bedgraph: str,
    control2_bedgraph: str,
    output_paths: Sequence[str],
    *,
    cutoff: float = 2.0,
    min_length: int = 200,
    max_gap: int = 50,
    depth1: float | None = None,
    depth2: float | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential peaks between two conditions from four MACS3 bedGraphs (``bdgdiff``).

    ``output_paths`` are three new files, in the order MACS3 writes them: condition 1 specific,
    condition 2 specific and common regions. ``depth1``/``depth2`` are sequencing depths in millions of
    non-redundant reads, used to scale the two conditions when they differ.
    """
    t1 = _input(treatment1_bedgraph, "treatment1_bedgraph")
    c1 = _input(control1_bedgraph, "control1_bedgraph")
    t2 = _input(treatment2_bedgraph, "treatment2_bedgraph")
    c2 = _input(control2_bedgraph, "control2_bedgraph")
    if (
        isinstance(output_paths, str)
        or not isinstance(output_paths, Sequence)
        or len(output_paths) != 3
    ):
        raise ConfigurationError(
            "output_paths must list three new files: condition 1, condition 2, common."
        )
    targets = [_output(path, "output_paths") for path in output_paths]
    if len({str(path) for path in targets}) != 3:
        raise ConfigurationError("output_paths must be three different files.")
    options = [
        "--t1",
        str(t1),
        "--c1",
        str(c1),
        "--t2",
        str(t2),
        "--c2",
        str(c2),
        "-C",
        str(_number(cutoff, "cutoff")),
        "-l",
        str(_int(min_length, "min_length", minimum=1)),
        "-g",
        str(_int(max_gap, "max_gap", minimum=0)),
        "-o",
        *map(str, targets),
    ]
    options += _opt("--d1", None if depth1 is None else _number(depth1, "depth1", minimum=1e-12))
    options += _opt("--d2", None if depth2 is None else _number(depth2, "depth2", minimum=1e-12))
    run = _run(
        "bdgdiff",
        options,
        cwd=targets[0].parent,
        monitored=targets,
        macs3_path=macs3_path,
        verbose=verbose,
    )
    labels = ("condition1_specific", "condition2_specific", "common")
    counts = {
        f"{label}_regions": _count_rows(path) for label, path in zip(labels, targets, strict=False)
    }
    return _result("bdgdiff", run, outputs=dict(zip(labels, targets, strict=False)), counts=counts)


def macs3_hmmratac(
    input_paths: Sequence[str],
    outdir: str,
    name: str,
    *,
    input_format: str = "BAMPE",
    lower: int | None = None,
    upper: int | None = None,
    prescan_cutoff: float | None = None,
    min_length: int | None = None,
    blacklist: str | None = None,
    save_digested: bool = False,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call open chromatin regions from paired-end ATAC-seq fragments with HMMRATAC (``macs3 hmmratac``).

    Input must be paired-end (BAMPE, BEDPE or FRAG). Writes ``<name>_peaks.narrowPeak`` (open regions),
    ``<name>_peaks.gappedPeak`` and the summit BED into the new directory ``outdir``. ``lower`` and ``upper``
    are the fragment-length bounds (in bp) for the nucleosome-free and mononucleosome signal, as in MACS3.
    """
    inputs = _inputs(input_paths, "input_paths")
    directory = _new_directory(outdir)
    label = _name(name)
    fmt = _choice(input_format, ("BAMPE", "BEDPE", "FRAG"), "input_format")
    options = ["-i", *map(str, inputs), "-f", fmt, "-n", label, "--outdir", str(directory)]
    options += _opt("-l", None if lower is None else _int(lower, "lower", minimum=1))
    options += _opt("-u", None if upper is None else _int(upper, "upper", minimum=1))
    options += _opt(
        "-c",
        None if prescan_cutoff is None else _number(prescan_cutoff, "prescan_cutoff", minimum=0.0),
    )
    options += _opt(
        "--minlen", None if min_length is None else _int(min_length, "min_length", minimum=1)
    )
    if blacklist is not None:
        options += ["-e", str(_input(blacklist, "blacklist"))]
    options += _switch("--save-digested", _flag(save_digested, "save_digested"))
    monitored = [directory / f"{label}_peaks.narrowPeak", directory / f"{label}_summits.bed"]
    run = _run(
        "hmmratac",
        options,
        cwd=directory.parent,
        monitored=monitored,
        macs3_path=macs3_path,
        verbose=verbose,
    )
    files = _directory_files(directory)
    counts = {"open_region_count": _count_rows(directory / f"{label}_peaks.narrowPeak")}
    return _result(
        "hmmratac", run, outputs={name: directory / name for name in files}, counts=counts
    )


__all__ = [
    "macs3_call_bedgraph_broad_peaks",
    "macs3_call_bedgraph_peaks",
    "macs3_call_peaks",
    "macs3_call_variants",
    "macs3_combine_replicates",
    "macs3_compare_bedgraphs",
    "macs3_diff_bedgraphs",
    "macs3_filter_duplicates",
    "macs3_hmmratac",
    "macs3_operate_bedgraph",
    "macs3_pileup",
    "macs3_predict_fragment_length",
    "macs3_random_sample",
    "macs3_refine_summits",
]

_register_public_verbose(globals(), tuple(__all__))
