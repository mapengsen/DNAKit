"""Bounded adapters for reviewed FASTQ command-line tools.

The adapters intentionally expose a small, typed option surface instead of
passing arbitrary command-line arguments.  They use the caller's installed
executable and never install packages, download databases, or invoke a shell.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
from collections.abc import Mapping
from pathlib import Path

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


def _mapping(value: object, name: str) -> dict[str, object]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping.")
    if any(not isinstance(key, str) for key in value):
        raise ConfigurationError(f"{name} keys must be strings.")
    return dict(value)


def _reject_unknown(options: Mapping[str, object], allowed: set[str], backend: str) -> None:
    unknown = sorted(set(options) - allowed)
    if unknown:
        raise ConfigurationError(
            f"Unsupported {backend} options: {unknown}.",
            code="EXTERNAL_OPTION_NOT_ALLOWLISTED",
        )


def _adapter(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 4096
        or value.startswith(("-", "file:"))
        or any(c in value for c in "\x00\r\n")
    ):
        raise ConfigurationError(f"{name} requires bounded inline adapter specifications.")
    return value


def _path(value: object, name: str, *, input_file: bool) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError(f"{name} must be a non-empty NUL-free path.")
    original = Path(value).expanduser().absolute()
    if not input_file and (original.exists() or original.is_symlink()):
        raise ConfigurationError(f"Refusing to overwrite existing output: {original}.")
    path = original.resolve()
    if input_file and not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}.")
    if not input_file:
        if path.exists() and not path.is_file():
            raise ConfigurationError(f"Output path is not a regular file: {path}.")
        if not path.parent.is_dir():
            raise ConfigurationError(f"Output directory does not exist: {path.parent}.")
    return path


def _read_pair(value: object) -> tuple[Path, Path | None]:
    if isinstance(value, Mapping):
        read1 = value.get("read1", value.get("r1", value.get("input1")))
        read2 = value.get("read2", value.get("r2", value.get("input2")))
        if read1 is None:
            source = value.get("input", value.get("source"))
            if source is None or read2 is not None:
                raise ConfigurationError(
                    "Input mapping requires read1/r1 or a single input/source path."
                )
            return _path(source, "input", input_file=True), None
        return _path(read1, "read1", input_file=True), (
            _path(read2, "read2", input_file=True) if read2 is not None else None
        )
    return _path(value, "input", input_file=True), None


def _input_files(value: object) -> list[Path]:
    if isinstance(value, Mapping):
        value = value.get("files", value.get("inputs", value.get("input", value.get("source"))))
    if isinstance(value, (str, Path)):
        values = [value]
    elif isinstance(value, (list, tuple)):
        values = list(value)
    else:
        raise ConfigurationError("Input must be a FASTQ path or a list of paths.")
    if not values or len(values) > 8:
        raise ConfigurationError("FastQC accepts 1 to 8 files per call (16 monitored reports).")
    return [_path(item, f"input[{index}]", input_file=True) for index, item in enumerate(values)]


def _directory(options: Mapping[str, object], name: str = "output_directory") -> Path:
    value = options.get(name)
    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError(f"options.{name} must be a non-empty directory path.")
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path}.")
    return path


def _file_name(value: object, name: str, default: str) -> str:
    resolved = default if value is None else value
    if (
        not isinstance(resolved, str)
        or not resolved
        or Path(resolved).name != resolved
        or resolved.startswith("-")
        or not resolved.endswith(".html")
        or any(c in resolved for c in "\x00\r\n")
    ):
        raise ConfigurationError(f"{name} must be a single .html filename.")
    return resolved


def _executable(options: Mapping[str, object], name: str) -> Path:
    value = options.get("executable_path", options.get("executable"))
    if value is not None and (
        not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value)
    ):
        raise ConfigurationError("executable_path must be a non-empty path.")
    resolved = shutil.which(name) if value is None else str(value)
    if not resolved:
        raise BackendUnavailableError(
            f"{name} executable is unavailable.",
            context={"backend_id": name, "executable": name},
            hint=f"Install {name} separately or pass options.executable_path.",
        )
    path = Path(resolved).expanduser().absolute()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise BackendUnavailableError(
            f"{name} executable is unavailable.",
            context={"backend_id": name, "executable": str(path)},
        )
    return path


def _integer(options: Mapping[str, object], name: str, *, minimum: int, maximum: int) -> int | None:
    value = options.get(name)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(
    options: Mapping[str, object], name: str, *, minimum: float, maximum: float
) -> float | None:
    value = options.get(name)
    if value is None or isinstance(value, bool) or not isinstance(value, (int, float)):
        if value is None:
            return None
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    if not minimum <= value <= maximum or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return float(value)


def _boolean(options: Mapping[str, object], name: str) -> bool:
    value = options.get(name, False)
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _output_paths(input_value: object, options: Mapping[str, object]) -> tuple[Path, Path | None]:
    _, read2 = _read_pair(input_value)
    output1 = options.get("output1", options.get("output", options.get("output_r1")))
    if output1 is None:
        raise ConfigurationError("options.output or options.output1 is required.")
    output2_value = options.get("output2", options.get("output_r2"))
    if read2 is not None and output2_value is None:
        raise ConfigurationError("Paired input requires options.output2/output_r2.")
    if read2 is None and output2_value is not None:
        raise ConfigurationError("output2 requires paired input.")
    output2 = (
        _path(output2_value, "output2", input_file=False) if output2_value is not None else None
    )
    output1_path = _path(output1, "output", input_file=False)
    if output2 is not None and output2.parent != output1_path.parent:
        raise ConfigurationError("All output files must share one existing parent directory.")
    if any(path.exists() for path in (output1_path, output2) if path is not None):
        raise ConfigurationError("Refusing to overwrite an existing output file.")
    return output1_path, output2


def _report_path(options: Mapping[str, object], output: Path, key: str, suffix: str) -> Path | None:
    value = options.get(key)
    if value is None:
        return None
    path = _path(value, key, input_file=False)
    if path.parent != output.parent:
        raise ConfigurationError("Report files must share the output parent directory.")
    if path.exists():
        raise ConfigurationError("Refusing to overwrite an existing report file.")
    if path.suffix.casefold() != suffix:
        raise ConfigurationError(f"{key} must use the {suffix} suffix.")
    return path


def _run(
    executable: Path,
    arguments: list[str],
    outputs: list[Path],
    *,
    backend_id: str,
    options: Mapping[str, object],
    verbose: bool,
    allow_empty: tuple[Path, ...] = (),
) -> dict[str, object]:
    timeout = _number(options, "timeout_seconds", minimum=0.01, maximum=86_400.0) or 3600.0
    max_output = _integer(options, "max_output_bytes", minimum=1, maximum=100_000_000) or 1_000_000
    if len(set(outputs)) != len(outputs):
        raise ConfigurationError("Output paths must be distinct.")
    for path in outputs:
        if path.exists() or path.is_symlink():
            raise ConfigurationError(f"Refusing to overwrite existing output: {path}.")
    max_artifacts = (
        _integer(options, "max_artifact_bytes", minimum=1, maximum=10_000_000_000) or 100_000_000
    )
    result = execute_bounded_command(
        executable,
        arguments,
        backend_id=backend_id,
        cwd=outputs[0].parent,
        timeout_seconds=timeout,
        max_output_bytes=max_output,
        monitored_output_paths=outputs,
        max_monitored_output_bytes=max_artifacts,
        verbose=verbose,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"{backend_id} exited with a non-zero status.",
            code="BACKEND_NONZERO_EXIT",
            context={
                "backend_id": backend_id,
                "return_code": result.return_code,
                "output": result.output[-4000:],
            },
        )
    for path in outputs:
        if path.is_symlink() or not path.is_file():
            raise BackendExecutionError(f"{backend_id} did not produce a regular output: {path}.")
        if path not in allow_empty and path.stat().st_size == 0:
            raise BackendExecutionError(f"{backend_id} produced an empty report: {path}.")
        if path.suffix == ".json":
            try:
                with path.open(encoding="utf-8") as handle:
                    report = json.load(handle)
                if not isinstance(report, dict):
                    raise ValueError("report must be an object")
            except (OSError, UnicodeError, ValueError) as exc:
                raise BackendExecutionError(f"Invalid {backend_id} JSON report: {path}.") from exc
    return {
        "status": "completed",
        "backend": backend_id,
        "executable": str(executable),
        "arguments": arguments,
        "outputs": [str(path) for path in outputs],
        "backend_output": result.output[-4000:],
        "provenance": {"implementation": "bounded_external_adapter", "backend": backend_id},
    }


def run_cutadapt(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> dict[str, object]:
    """Trim FASTQ adapters with a caller-installed Cutadapt executable.

    ``input`` is a FASTQ path or ``{"read1": ..., "read2": ...}``. Required
    output paths and adapter specifications are supplied in ``options``.
    """

    opts = _mapping(options, "options")
    _reject_unknown(
        opts,
        {
            "executable_path",
            "executable",
            "output",
            "output1",
            "output_r1",
            "output2",
            "output_r2",
            "adapters",
            "adapters2",
            "adapter_5prime_r2",
            "adapter_anywhere_r2",
            "adapter_5prime",
            "adapter_anywhere",
            "cores",
            "minimum_overlap",
            "minimum_length",
            "error_rate",
            "quality_cutoff",
            "discard_untrimmed",
            "json_report",
            "timeout_seconds",
            "max_output_bytes",
            "max_artifact_bytes",
        },
        "Cutadapt",
    )
    read1, read2 = _read_pair(input)
    output1, output2 = _output_paths(input, opts)
    executable = _executable(opts, "cutadapt")
    arguments: list[str] = []
    for key, flag in (("adapters", "-a"), ("adapters2", "-A")):
        adapters = opts.get(key, [])
        if not isinstance(adapters, list) or len(adapters) > 32:
            raise ConfigurationError(f"{key} must be a list of up to 32 adapter strings.")
        if key == "adapters2" and adapters and read2 is None:
            raise ConfigurationError("adapters2 requires paired input.")
        for adapter in adapters:
            arguments.extend([flag, _adapter(adapter, key)])
    for key, flag in (
        ("adapter_5prime", "-g"),
        ("adapter_anywhere", "-b"),
        ("adapter_5prime_r2", "-G"),
        ("adapter_anywhere_r2", "-B"),
    ):
        value = opts.get(key)
        if value is not None:
            if key.endswith("_r2") and read2 is None:
                raise ConfigurationError(f"{key} requires paired input.")
            arguments.extend([flag, _adapter(value, key)])
    if not arguments:
        raise ConfigurationError("At least one adapter specification is required.")
    arguments.extend(["-o", str(output1)])
    if output2 is not None:
        arguments.extend(["-p", str(output2)])
    for key, flag, low, high in (
        ("cores", "-j", 1, 256),
        ("minimum_overlap", "-O", 1, 10_000),
        ("minimum_length", "-m", 0, 10_000_000),
    ):
        value = _integer(opts, key, minimum=low, maximum=high)
        if value is not None:
            arguments.extend([flag, str(value)])
    error_rate = _number(opts, "error_rate", minimum=0.0, maximum=0.999999)
    if error_rate is not None:
        arguments.extend(["-e", str(error_rate)])
    quality_cutoff = _integer(opts, "quality_cutoff", minimum=0, maximum=93)
    if quality_cutoff is not None:
        arguments.extend(["-q", str(quality_cutoff)])
    if _boolean(opts, "discard_untrimmed"):
        arguments.append("--discard-untrimmed")
    report = _report_path(opts, output1, "json_report", ".json")
    if report is not None:
        arguments.extend(["--json", str(report)])
    arguments.extend([str(read1), *([str(read2)] if read2 is not None else [])])
    outputs = [
        output1,
        *([output2] if output2 is not None else []),
        *([report] if report is not None else []),
    ]
    result = _run(
        executable,
        arguments,
        outputs,
        backend_id="cutadapt",
        options=opts,
        verbose=verbose,
        allow_empty=tuple(path for path in (output1, output2) if path is not None),
    )
    result.update(
        {
            "input": [str(read1), *([str(read2)] if read2 is not None else [])],
            "citation_doi": "10.14806/ej.17.1.200",
            "license": "MIT",
        }
    )
    return result


def run_fastp(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> dict[str, object]:
    """Preprocess FASTQ reads with a caller-installed fastp executable."""

    opts = _mapping(options, "options")
    _reject_unknown(
        opts,
        {
            "executable_path",
            "executable",
            "output",
            "output1",
            "output_r1",
            "output2",
            "output_r2",
            "thread",
            "qualified_quality_phred",
            "unqualified_percent_limit",
            "n_base_limit",
            "average_qual",
            "length_required",
            "disable_quality_filtering",
            "disable_length_filtering",
            "trim_poly_g",
            "trim_poly_x",
            "cut_front",
            "cut_tail",
            "cut_right",
            "adapter_sequence",
            "adapter_sequence_r2",
            "detect_adapter_for_pe",
            "dont_eval_duplication",
            "json_report",
            "html_report",
            "timeout_seconds",
            "max_output_bytes",
            "max_artifact_bytes",
        },
        "fastp",
    )
    read1, read2 = _read_pair(input)
    output1, output2 = _output_paths(input, opts)
    executable = _executable(opts, "fastp")
    arguments = ["-i", str(read1), "-o", str(output1)]
    if read2 is not None:
        arguments.extend(["-I", str(read2), "-O", str(output2)])
    for key, flag, low, high in (
        ("thread", "-w", 1, 256),
        ("qualified_quality_phred", "-q", 0, 93),
        ("unqualified_percent_limit", "-u", 0, 100),
        ("n_base_limit", "-n", 0, 10_000_000),
        ("average_qual", "-e", 0, 93),
        ("length_required", "-l", 0, 10_000_000),
    ):
        value = _integer(opts, key, minimum=low, maximum=high)
        if value is not None:
            arguments.extend([flag, str(value)])
    for key, flag in (
        ("disable_quality_filtering", "-Q"),
        ("disable_length_filtering", "-L"),
        ("trim_poly_g", "-g"),
        ("trim_poly_x", "-x"),
        ("detect_adapter_for_pe", "--detect_adapter_for_pe"),
        ("dont_eval_duplication", "--dont_eval_duplication"),
    ):
        if _boolean(opts, key):
            arguments.append(flag)
    for key, flag in (
        ("cut_front", "--cut_front"),
        ("cut_tail", "--cut_tail"),
        ("cut_right", "--cut_right"),
    ):
        if _boolean(opts, key):
            arguments.append(flag)
    for key in ("adapter_sequence", "adapter_sequence_r2"):
        value = opts.get(key)
        if value is not None:
            if key.endswith("_r2") and read2 is None:
                raise ConfigurationError(f"{key} requires paired input.")
            sequence = _adapter(value, key)
            if re.fullmatch("[ACGTacgt]+", sequence) is None:
                raise ConfigurationError(f"{key} must contain only A/C/G/T.")
            arguments.extend(["--" + key, sequence])
    opts.setdefault("json_report", output1.with_name(output1.name + ".fastp.json"))
    opts.setdefault("html_report", output1.with_name(output1.name + ".fastp.html"))
    report_json = _report_path(opts, output1, "json_report", ".json")
    report_html = _report_path(opts, output1, "html_report", ".html")
    if report_json is None or report_html is None:
        raise ConfigurationError("fastp report paths cannot be null.")
    arguments.extend(["-j", str(report_json), "-h", str(report_html)])
    outputs = [
        output1,
        *([output2] if output2 is not None else []),
        *([report_json] if report_json is not None else []),
        *([report_html] if report_html is not None else []),
    ]
    result = _run(
        executable,
        arguments,
        outputs,
        backend_id="fastp",
        options=opts,
        verbose=verbose,
        allow_empty=tuple(path for path in (output1, output2) if path is not None),
    )
    result.update(
        {
            "input": [str(read1), *([str(read2)] if read2 is not None else [])],
            "citation_doi": "10.1002/imt2.70078",
            "license": "MIT",
        }
    )
    return result


def run_fastqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> dict[str, object]:
    """Run FastQC on caller-provided read files and monitor HTML/ZIP reports."""

    opts = _mapping(options, "options")
    _reject_unknown(
        opts,
        {
            "executable_path",
            "executable",
            "output_directory",
            "threads",
            "noextract",
            "nogroup",
            "quiet",
            "format",
            "contaminant",
            "adapter",
            "kmer",
            "timeout_seconds",
            "max_output_bytes",
            "max_artifact_bytes",
        },
        "FastQC",
    )
    inputs = _input_files(input)
    output_directory = _directory(opts)
    executable = _executable(opts, "fastqc")
    threads = _integer(opts, "threads", minimum=1, maximum=256)
    arguments: list[str] = ["--outdir", str(output_directory)]
    if threads is not None:
        arguments.extend(["--threads", str(threads)])
    for key, flag in (("noextract", "--noextract"), ("nogroup", "--nogroup"), ("quiet", "--quiet")):
        if _boolean(opts, key):
            arguments.append(flag)
    format_value = opts.get("format")
    if format_value is not None:
        if not isinstance(format_value, str) or format_value not in {"fastq", "bam"}:
            raise ConfigurationError("format must be 'fastq' or 'bam'.")
        arguments.extend(["--format", str(format_value)])
    for key, flag in (
        ("contaminant", "--contaminants"),
        ("adapter", "--adapters"),
    ):
        value = opts.get(key)
        if value is not None:
            if not isinstance(value, (str, Path)) or not str(value):
                raise ConfigurationError(f"{key} must be a non-empty path string.")
            arguments.extend([flag, str(_path(value, key, input_file=True))])
    kmer = _integer(opts, "kmer", minimum=2, maximum=10)
    if kmer is not None:
        arguments.extend(["--kmers", str(kmer)])
    # FastQC 0.12.x extracts ZIP data by default in noninteractive mode.
    if "noextract" not in opts:
        arguments.append("--noextract")
    elif not opts["noextract"] and any(output_directory.iterdir()):
        raise ConfigurationError("FastQC extraction requires an empty output directory.")
    arguments.extend(str(path) for path in inputs)
    outputs: list[Path] = []
    for path in inputs:
        base = re.sub(r"(?i)\.gz$", "", path.name)
        base = re.sub(r"(?i)\.(fastq|fq|bam|sam|txt)$", "", base)
        for suffix in ("_fastqc.html", "_fastqc.zip"):
            output = output_directory / f"{base}{suffix}"
            if output.exists():
                raise ConfigurationError(
                    f"Refusing to overwrite an existing output file: {output}."
                )
            outputs.append(output)
    result = _run(
        executable, arguments, outputs, backend_id="fastqc", options=opts, verbose=verbose
    )
    result.update(
        {
            "input": [str(path) for path in inputs],
            "citation": "FastQC documentation; no DOI supplied by the survey",
            "license": "GPL-3.0",
        }
    )
    return result


def run_multiqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> dict[str, object]:
    """Aggregate caller-provided QC inputs with a bounded MultiQC invocation."""

    opts = _mapping(options, "options")
    _reject_unknown(
        opts,
        {
            "executable_path",
            "executable",
            "output_directory",
            "filename",
            "title",
            "force",
            "interactive",
            "zip_data",
            "flat",
            "dirs",
            "dirs_depth",
            "config",
            "timeout_seconds",
            "max_output_bytes",
            "max_artifact_bytes",
        },
        "MultiQC",
    )
    raw_inputs = (
        input.get("inputs", input.get("files", input.get("input")))
        if isinstance(input, Mapping)
        else input
    )
    if isinstance(raw_inputs, (str, Path)):
        raw_values = [raw_inputs]
    elif isinstance(raw_inputs, (list, tuple)):
        raw_values = list(raw_inputs)
    else:
        raise ConfigurationError("Input must be a file or directory path, or a list of paths.")
    if not raw_values or len(raw_values) > 256:
        raise ConfigurationError("Input must contain between 1 and 256 paths.")
    inputs = []
    for index, value in enumerate(raw_values):
        if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
            raise ConfigurationError(f"input[{index}] must be a non-empty path.")
        path = Path(value).expanduser().resolve()
        if not path.exists() or not (path.is_file() or path.is_dir()):
            raise ConfigurationError(f"Input path does not exist: {path}.")
        inputs.append(path)
    output_directory = _directory(opts)
    if any(output_directory.iterdir()):
        raise ConfigurationError("MultiQC requires an empty output directory for reports and data.")
    if _boolean(opts, "interactive") and _boolean(opts, "flat"):
        raise ConfigurationError("interactive and flat cannot both be true.")
    executable = _executable(opts, "multiqc")
    filename = _file_name(opts.get("filename"), "filename", "multiqc_report.html")
    report = output_directory / filename
    if report.exists():
        raise ConfigurationError(f"Refusing to overwrite an existing report: {report}.")
    arguments = ["--outdir", str(output_directory), "--filename", filename]
    for key, flag in (
        ("force", "--force"),
        ("interactive", "--interactive"),
        ("zip_data", "--zip-data-dir"),
        ("flat", "--flat"),
        ("dirs", "--dirs"),
    ):
        if _boolean(opts, key):
            arguments.append(flag)
    dirs_depth = _integer(opts, "dirs_depth", minimum=-100, maximum=100)
    if dirs_depth is not None:
        arguments.extend(["--dirs-depth", str(dirs_depth)])
    title = opts.get("title")
    if title is not None:
        if not isinstance(title, str) or not title or "\x00" in title or len(title) > 256:
            raise ConfigurationError("title must be a bounded non-empty string.")
        arguments.extend(["--title", title])
    config = opts.get("config")
    if config is not None:
        arguments.extend(["--config", str(_path(config, "config", input_file=True))])
    arguments.extend(str(path) for path in inputs)
    result = _run(
        executable, arguments, [report], backend_id="multiqc", options=opts, verbose=verbose
    )
    result.update(
        {
            "input": [str(path) for path in inputs],
            "citation_doi": "10.1093/bioinformatics/btw354",
            "license": "GPL-3.0-or-later",
        }
    )
    return result


__all__ = ["run_cutadapt", "run_fastp", "run_fastqc", "run_multiqc"]
