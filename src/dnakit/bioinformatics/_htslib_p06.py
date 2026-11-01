"""Private helpers shared by the samtools / bcftools command-line facades of integration pack P06.

A facade runs one fixed ``samtools`` or ``bcftools`` command with an explicit argument list (no shell).
Unlike ``_gap_cli.run``, which merges standard error into standard output, these programs print their
data (reports, VCF text, histograms) on standard output and their messages on standard error, so the two
streams are captured separately: standard output goes straight into a new file or a temporary file,
standard error into a bounded log tail. The program is resolved by ``_gap_cli.program`` (explicit path,
PATH, interpreter directory). A bcftools plugin is found through ``BCFTOOLS_PLUGINS``; when the variable
is not set and a ``libexec/bcftools`` folder sits next to the executable's ``bin`` folder, that folder is
used, so a relocated installation keeps its plugins.
"""

from __future__ import annotations

import math
import os
import re
import signal
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

from dnakit.bioinformatics import _gap_cli as cli
from dnakit.exceptions import BackendExecutionError, BackendTimeoutError, ConfigurationError

TESTED_VERSION = "1.24"
CITATION = "https://doi.org/10.1093/gigascience/giab008"
TIMEOUT_SECONDS = 3600.0
LOG_TAIL_CHARS = 4000
MAX_TEXT_BYTES = 20_000_000
MAX_PARSED_ROWS = 10_000
OUTPUT_TYPES = ("v", "z", "b", "u")

HINTS: Mapping[str, str] = {
    "samtools": (
        "Install samtools 1.24 or later (for example `conda install -c bioconda samtools`, or a release "
        "from https://www.htslib.org/download/), then put it on PATH or pass samtools_path."
    ),
    "bcftools": (
        "Install bcftools 1.24 or later with its plugins (for example `conda install -c bioconda "
        "bcftools`, or a release from https://www.htslib.org/download/), then put it on PATH or pass "
        "bcftools_path. Plugins are found through BCFTOOLS_PLUGINS or the libexec/bcftools folder "
        "next to the installation."
    ),
}

_NUMBER = re.compile(r"^-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$")


@dataclass(frozen=True)
class Completed:
    """What one finished command produced."""

    tool: str
    executable: str
    version: str | None
    command: tuple[str, ...]
    return_code: int
    stdout_path: Path | None
    stdout: str
    stderr_tail: str
    stderr: str = ""


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def existing_files(values: object, name: str, *, minimum: int = 1) -> list[Path]:
    """A list of existing files (a single path string is accepted as a list of one)."""
    if isinstance(values, (str, Path)):
        values = [values]
    if isinstance(values, (bytes, bytearray)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    if len(values) < minimum:
        raise ConfigurationError(f"{name} needs at least {minimum} file(s).")
    return [cli.existing_file(value, name) for value in values]


def text_list(values: object, name: str) -> list[str]:
    """Comma-joined list of safe text items (samples, tags, regions)."""
    if isinstance(values, str):
        values = [values]
    if isinstance(values, (bytes, bytearray)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    items = [cli.text(value, name) for value in values]
    if any("," in item for item in items):
        raise ConfigurationError(f"{name} items must not contain commas.")
    return items


def selection(
    *,
    regions: str | None = None,
    targets: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    samples: Sequence[str] | None = None,
) -> list[str]:
    """The record/sample selection options every VCF/BCF-reading bcftools command shares."""
    if include is not None and exclude is not None:
        raise ConfigurationError("include and exclude are mutually exclusive.")
    arguments: list[str] = []
    if regions is not None:
        arguments += ["-r", cli.text(regions, "regions")]
    if targets is not None:
        arguments += ["-t", cli.text(targets, "targets")]
    if include is not None:
        arguments += ["-i", cli.text(include, "include")]
    if exclude is not None:
        arguments += ["-e", cli.text(exclude, "exclude")]
    if samples is not None:
        arguments += ["-s", ",".join(text_list(samples, "samples"))]
    return arguments


def output_type_for(path: Path, output_type: str | None) -> str:
    """Explicit ``output_type`` or the one the file name implies (.vcf.gz -> z, .bcf -> b, else v)."""
    if output_type is not None:
        return cli.choice(output_type, OUTPUT_TYPES, "output_type")
    name = path.name.lower()
    if name.endswith((".vcf.gz", ".vcf.bgz")):
        return "z"
    if name.endswith(".bcf"):
        return "b"
    return "v"


def version_tuple(version: str | None) -> tuple[int, ...]:
    return tuple(int(part) for part in (version or "0").split(".") if part.isdigit())


def require_version(done_version: str | None, minimum: str, feature: str) -> None:
    if version_tuple(done_version) < version_tuple(minimum):
        raise BackendExecutionError(
            f"{feature} needs bcftools {minimum} or later (found {done_version}).",
            code="HTSLIB_VERSION_TOO_OLD",
        )


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def executable_for(tool: str, explicit: str | None) -> str:
    return cli.program(tool, explicit=explicit, hint=HINTS[tool])


def version_of(executable: str) -> str | None:
    return cli._version_of(executable)


def _environment(executable: str, plugin: bool) -> dict[str, str] | None:
    if not plugin or "BCFTOOLS_PLUGINS" in os.environ:
        return None
    folder = Path(executable).parent.parent / "libexec" / "bcftools"
    return {"BCFTOOLS_PLUGINS": str(folder)} if folder.is_dir() else None


def _kill(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        process.kill()
    process.wait()


def invoke(
    tool: str,
    explicit: str | None,
    arguments: Sequence[str],
    *,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    plugin: bool = False,
    created: Sequence[Path] = (),
    accept_codes: Sequence[int] = (0,),
    timeout_seconds: float = TIMEOUT_SECONDS,
) -> Completed:
    """Run ``tool`` once without a shell.

    Standard output is written to the new file ``stdout_path`` when given (never overwritten), else
    kept (up to ``MAX_TEXT_BYTES``) in memory. A return code outside ``accept_codes`` raises
    ``BackendExecutionError`` with the standard-error tail, and removes the files in ``created`` and
    ``stdout_path`` that the failed run left behind.
    """
    executable = executable_for(tool, explicit)
    version = version_of(executable)
    command = tuple(str(argument) for argument in arguments)
    if any("\x00" in argument for argument in command):
        raise ConfigurationError("Arguments must not contain NUL characters.")
    environment = os.environ.copy()
    extra = _environment(executable, plugin)
    if extra is not None:
        environment.update(extra)
    environment.update({"LANG": "C", "LC_ALL": "C"})
    scratch: tempfile.TemporaryDirectory[str] | None = None
    if cwd is None:
        scratch = tempfile.TemporaryDirectory(prefix="dnakit-htslib-")
        cwd = Path(scratch.name)
    try:
        with tempfile.TemporaryFile() as error_file:
            out_handle: IO[bytes] = (
                stdout_path.open("xb") if stdout_path is not None else tempfile.TemporaryFile()  # noqa: SIM115
            )
            try:
                process = subprocess.Popen(
                    [executable, *command],
                    cwd=str(cwd),
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=out_handle,
                    stderr=error_file,
                    shell=False,
                    start_new_session=os.name == "posix",
                )
            except OSError as exc:
                out_handle.close()
                if stdout_path is not None:
                    stdout_path.unlink(missing_ok=True)
                raise BackendExecutionError(
                    f"{tool} could not be started: {type(exc).__name__}.",
                    code="HTSLIB_START_FAILED",
                ) from exc
            deadline = time.monotonic() + float(timeout_seconds)
            try:
                while process.poll() is None:
                    if time.monotonic() >= deadline:
                        _kill(process)
                        out_handle.close()
                        for leftover in (*created, *([stdout_path] if stdout_path else [])):
                            leftover.unlink(missing_ok=True)
                        raise BackendTimeoutError(f"{tool} exceeded its time limit.")
                    time.sleep(0.005)
            finally:
                if process.poll() is None:
                    _kill(process)
            return_code = int(process.returncode)
            text = ""
            if stdout_path is None:
                out_handle.flush()
                out_handle.seek(0)
                text = out_handle.read(MAX_TEXT_BYTES).decode("utf-8", errors="replace")
            out_handle.close()
            error_file.seek(0)
            stderr_text = error_file.read(MAX_TEXT_BYTES).decode("utf-8", errors="replace")
            stderr_tail = stderr_text[-LOG_TAIL_CHARS * 4 :]
    finally:
        if scratch is not None:
            scratch.cleanup()
    if return_code not in accept_codes:
        for leftover in (*created, *([stdout_path] if stdout_path else [])):
            leftover.unlink(missing_ok=True)
        raise BackendExecutionError(
            f"{tool} {command[0] if command else ''} exited with code {return_code}.",
            code="HTSLIB_EXECUTION_FAILED",
            context={"command": [tool, *command][:12], "log_tail": stderr_tail[-LOG_TAIL_CHARS:]},
        )
    return Completed(
        tool, executable, version, (tool, *command), return_code, stdout_path, text, stderr_tail, stderr_text
    )


def output_check(done: Completed, outputs: Mapping[str, Path]) -> None:
    missing = [str(path) for path in outputs.values() if not path.exists()]
    if missing:
        raise BackendExecutionError(
            f"{done.tool} did not write the expected output: {missing}",
            code="HTSLIB_OUTPUT_MISSING",
            context={"log_tail": done.stderr_tail[-LOG_TAIL_CHARS:]},
        )


def result(
    operation: str,
    done: Completed,
    *,
    outputs: Mapping[str, Path] | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """The common JSON-serialisable result dictionary of every P06 facade."""
    files = dict(outputs or {})
    output_check(done, files)
    payload: dict[str, Any] = {
        "status": "completed",
        "backend": done.tool,
        "operation": operation,
        "version": done.version,
        "command": list(done.command),
        "outputs": {name: str(path) for name, path in files.items()},
        "log_tail": done.stderr_tail[-LOG_TAIL_CHARS:],
        "provenance": {
            "tool": done.tool,
            "version": done.version,
            "tested_version": TESTED_VERSION,
            "executable": done.executable,
            "citation": CITATION,
        },
    }
    payload.update(extra)
    return payload


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def number(value: str) -> object:
    """int, float or the unchanged text."""
    if _NUMBER.match(value):
        try:
            return int(value)
        except ValueError:
            parsed = float(value)
            return parsed if math.isfinite(parsed) else value
    return value


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def tab_table(text: str, *, header: str = "auto") -> dict[str, Any]:
    """Split tab-separated report text into comment lines, an optional header and numeric rows.

    ``header`` is ``auto`` (the last ``#`` line whose first field starts with ``[1]`` or ``#`` + a
    column name counts as the header when the report has no other one), ``first`` (the first data line)
    or ``none``. Rows are capped at ``MAX_PARSED_ROWS``; ``truncated`` tells whether that happened.
    """
    comments: list[str] = []
    rows: list[list[object]] = []
    names: list[str] | None = None
    truncated = False
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            comments.append(line)
            continue
        fields = line.rstrip("\t").split("\t")
        if header == "first" and names is None:
            names = fields
            continue
        if len(rows) >= MAX_PARSED_ROWS:
            truncated = True
            continue
        rows.append([number(field) for field in fields])
    if header == "auto" and comments:
        last = comments[-1].lstrip("#").strip()
        candidate = [re.sub(r"^\[\d+\]", "", part).strip() for part in last.split("\t")]
        if len(candidate) > 1 and rows and len(candidate) == len(rows[0]):
            names = candidate
    return {"comments": comments, "columns": names, "rows": rows, "truncated": truncated}


def vcf_record_count(path: Path) -> int | None:
    """Number of records in a plain-text VCF, or None for compressed/binary output."""
    if path.suffix.lower() in {".gz", ".bgz", ".bcf"}:
        return None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def directory_files(directory: Path) -> list[str]:
    return sorted(path.name for path in directory.iterdir() if path.is_file())


# ---------------------------------------------------------------------------
# bcftools plugins
# ---------------------------------------------------------------------------


def plugin_command(
    plugin: str,
    inputs: Sequence[Path],
    *,
    general: Sequence[str] = (),
    options: Sequence[str] = (),
    separator: bool = True,
) -> list[str]:
    """``+plugin general... inputs -- options...`` (classic plugins) or ``+plugin options... inputs``.

    Plugins written with the classic interface take their own options after ``--``; the newer ones
    (prune, split, scatter, trio-stats, mendelian2, trio-dnm3, ...) read all options before the file.
    """
    if separator:
        return [f"+{plugin}", *general, *(str(path) for path in inputs), "--", *options]
    return [f"+{plugin}", *general, *options, *(str(path) for path in inputs)]


def run_vcf_plugin(
    operation: str,
    plugin: str,
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    options: Sequence[str] = (),
    general: Sequence[str] = (),
    separator: bool = True,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    extra_inputs: Sequence[Path] = (),
    **extra: Any,
) -> dict[str, Any]:
    """Run a plugin that writes a VCF/BCF (``-o``); returns the common result with a record count."""
    source = cli.existing_file(vcf_path, "vcf_path")
    target = cli.new_file(output_path, "output_path")
    kind = output_type_for(target, output_type)
    arguments = plugin_command(
        plugin,
        [source, *extra_inputs],
        general=[*general, "-O", kind, "-o", str(target)],
        options=options,
        separator=separator,
    )
    done = invoke("bcftools", bcftools_path, arguments, plugin=True, created=[target])
    return result(
        operation,
        done,
        outputs={"output": target},
        output_type=kind,
        records=vcf_record_count(target),
        messages=[line for line in done.stderr_tail.splitlines() if line.strip()][-8:],
        **extra,
    )


def run_text_plugin(
    operation: str,
    plugin: str,
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    options: Sequence[str] = (),
    general: Sequence[str] = (),
    separator: bool = True,
    header: str = "auto",
    bcftools_path: str | None = None,
    extra_inputs: Sequence[Path] = (),
    environment_cwd: Path | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """Run a plugin whose result is text on standard output; the text is saved to ``output_path``."""
    source = cli.existing_file(vcf_path, "vcf_path")
    target = cli.new_file(output_path, "output_path")
    arguments = plugin_command(
        plugin, [source, *extra_inputs], general=general, options=options, separator=separator
    )
    done = invoke(
        "bcftools", bcftools_path, arguments, plugin=True, stdout_path=target, cwd=environment_cwd
    )
    report = tab_table(read_text(target), header=header)
    return result(operation, done, outputs={"report": target}, report=report, **extra)
