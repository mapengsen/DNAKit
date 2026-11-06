"""cutadapt command-line operations as typed DNAKit facades: adapter and quality trimming of FASTQ reads.

The user installs cutadapt (``pip install cutadapt``, or ``conda install -c bioconda cutadapt``). Each function runs
the cutadapt program with an argument list built from explicit typed parameters mapped to fixed flags; no shell is
used; outputs are new explicit paths; results are JSON-serialisable dictionaries with status, backend, version,
command, outputs, the run report and a bounded log tail.
"""

from __future__ import annotations

import functools
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = (
    "Install cutadapt (tested with 5.2: python -m pip install cutadapt, or conda install -c bioconda cutadapt), "
    "or pass cutadapt_path."
)
_MAX_TAIL = 4000
_TIMEOUT_SECONDS = 3600.0
_ADAPTER = re.compile(r"^[ACGTNRYSWKMBDHVU]+$", re.IGNORECASE)
_QUALITY = re.compile(r"^(\d+,)?\d+$")
_REPORTS = ("full", "minimal")


def _program(cutadapt_path: str | None) -> str:
    """Resolve the cutadapt executable: explicit path, then PATH, then the interpreter's own bin folder."""
    if cutadapt_path is not None:
        candidate = Path(cutadapt_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"cutadapt_path does not exist: {candidate}", hint=_HINT)
        return str(candidate)
    found = shutil.which("cutadapt")
    if found is None:
        local = Path(sys.executable).parent / "cutadapt"
        if local.is_file():
            found = str(local)
    if found is None:
        raise BackendUnavailableError("cutadapt was not found on PATH.", hint=_HINT)
    return found


@functools.lru_cache(maxsize=8)
def _version(program: str) -> str:
    completed = subprocess.run(
        [program, "--version"],
        capture_output=True,
        text=True,
        timeout=120,
        shell=False,
        check=False,
    )
    return (completed.stdout or completed.stderr).strip()


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _adapter(value: str | None, name: str) -> str | None:
    if value is None:
        return None
    if not _ADAPTER.match(value):
        raise ConfigurationError(f"{name} must be a non-empty DNA sequence of IUPAC symbols.")
    return value.upper()


def _length(value: int | None, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer or None.")
    return value


def _trim(
    input_path: str | Path,
    output_path: str | Path,
    *,
    adapter_3p: str | None,
    adapter_5p: str | None,
    anywhere_adapter: str | None,
    error_rate: float,
    quality_cutoff: str | None,
    minimum_length: int | None,
    maximum_length: int | None,
    discard_untrimmed: bool,
    cores: int,
    report: str,
    json_path: str | Path | None,
    cutadapt_path: str | None,
) -> dict[str, Any]:
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    report_json = _new(json_path, "json_path") if json_path is not None else None
    if (
        isinstance(error_rate, bool)
        or not isinstance(error_rate, (int, float))
        or not 0 <= error_rate <= 1
    ):
        raise ConfigurationError("error_rate must be a number between 0 and 1.")
    if isinstance(cores, bool) or not isinstance(cores, int) or cores < 1:
        raise ConfigurationError("cores must be a positive integer.")
    if report not in _REPORTS:
        raise ConfigurationError(f"report must be one of {_REPORTS}.")
    if quality_cutoff is not None and not _QUALITY.match(quality_cutoff):
        raise ConfigurationError(
            "quality_cutoff must look like '20' or '10,20' (5' cutoff, 3' cutoff)."
        )
    minimum = _length(minimum_length, "minimum_length")
    maximum = _length(maximum_length, "maximum_length")
    if minimum is not None and maximum is not None and maximum < minimum:
        raise ConfigurationError("maximum_length must not be smaller than minimum_length.")

    args: list[str] = ["--report", report, "-e", str(float(error_rate)), "-j", str(cores)]
    for flag, value in (
        ("-a", _adapter(adapter_3p, "adapter_3p")),
        ("-g", _adapter(adapter_5p, "adapter_5p")),
        ("-b", _adapter(anywhere_adapter, "anywhere_adapter")),
    ):
        if value is not None:
            args += [flag, value]
    if quality_cutoff is not None:
        args += ["-q", quality_cutoff]
    if minimum is not None:
        args += ["-m", str(minimum)]
    if maximum is not None:
        args += ["-M", str(maximum)]
    if discard_untrimmed:
        args.append("--discard-untrimmed")
    if report_json is not None:
        args += ["--json", str(report_json)]
    args += ["-o", str(target), str(source)]

    program = _program(cutadapt_path)
    completed = subprocess.run(
        [program, *args],
        capture_output=True,
        text=True,
        timeout=_TIMEOUT_SECONDS,
        shell=False,
        check=False,
        cwd=target.parent,
    )
    outputs = [p for p in (target, report_json) if p is not None and p.is_file()]
    status = (
        "completed"
        if completed.returncode == 0 and len(outputs) == (2 if report_json else 1)
        else "failed"
    )
    result: dict[str, Any] = {
        "status": status,
        "backend": "cutadapt",
        "version": _version(program),
        "command": ["cutadapt", *args],
        "return_code": completed.returncode,
        "outputs": [str(p) for p in outputs],
        "report": completed.stdout,
        "log_tail": (completed.stderr or "")[-_MAX_TAIL:],
    }
    if report_json is not None and report_json.is_file():
        result["json"] = json.loads(report_json.read_text(encoding="utf-8"))
    return result


def cutadapt_trim(
    input_path: str,
    output_path: str,
    *,
    adapter_3p: str | None = None,
    adapter_5p: str | None = None,
    anywhere_adapter: str | None = None,
    error_rate: float = 0.1,
    quality_cutoff: str | None = None,
    minimum_length: int | None = None,
    maximum_length: int | None = None,
    discard_untrimmed: bool = False,
    cores: int = 1,
    report: str = "full",
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim 3' and 5' adapters and low-quality ends of one FASTQ file (cutadapt -a, -g, -b, -q, -m, -M, -e)."""
    return _trim(
        input_path,
        output_path,
        adapter_3p=adapter_3p,
        adapter_5p=adapter_5p,
        anywhere_adapter=anywhere_adapter,
        error_rate=error_rate,
        quality_cutoff=quality_cutoff,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        discard_untrimmed=discard_untrimmed,
        cores=cores,
        report=report,
        json_path=None,
        cutadapt_path=cutadapt_path,
    )


def cutadapt_full_report(
    input_path: str,
    output_path: str,
    *,
    adapter_3p: str | None = None,
    adapter_5p: str | None = None,
    anywhere_adapter: str | None = None,
    error_rate: float = 0.1,
    quality_cutoff: str | None = None,
    minimum_length: int | None = None,
    maximum_length: int | None = None,
    discard_untrimmed: bool = False,
    cores: int = 1,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a trimming job and return cutadapt's full text report (summary, adapter and length sections)."""
    return _trim(
        input_path,
        output_path,
        adapter_3p=adapter_3p,
        adapter_5p=adapter_5p,
        anywhere_adapter=anywhere_adapter,
        error_rate=error_rate,
        quality_cutoff=quality_cutoff,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        discard_untrimmed=discard_untrimmed,
        cores=cores,
        report="full",
        json_path=None,
        cutadapt_path=cutadapt_path,
    )


def cutadapt_minimal_report(
    input_path: str,
    output_path: str,
    *,
    adapter_3p: str | None = None,
    adapter_5p: str | None = None,
    anywhere_adapter: str | None = None,
    error_rate: float = 0.1,
    quality_cutoff: str | None = None,
    minimum_length: int | None = None,
    maximum_length: int | None = None,
    discard_untrimmed: bool = False,
    cores: int = 1,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a trimming job and return cutadapt's one-line tab-separated minimal report."""
    return _trim(
        input_path,
        output_path,
        adapter_3p=adapter_3p,
        adapter_5p=adapter_5p,
        anywhere_adapter=anywhere_adapter,
        error_rate=error_rate,
        quality_cutoff=quality_cutoff,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        discard_untrimmed=discard_untrimmed,
        cores=cores,
        report="minimal",
        json_path=None,
        cutadapt_path=cutadapt_path,
    )


def cutadapt_json_report(
    input_path: str,
    output_path: str,
    json_path: str,
    *,
    adapter_3p: str | None = None,
    adapter_5p: str | None = None,
    anywhere_adapter: str | None = None,
    error_rate: float = 0.1,
    quality_cutoff: str | None = None,
    minimum_length: int | None = None,
    maximum_length: int | None = None,
    discard_untrimmed: bool = False,
    cores: int = 1,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a trimming job and write cutadapt's machine-readable JSON report to a new file (returned as ``json``)."""
    return _trim(
        input_path,
        output_path,
        adapter_3p=adapter_3p,
        adapter_5p=adapter_5p,
        anywhere_adapter=anywhere_adapter,
        error_rate=error_rate,
        quality_cutoff=quality_cutoff,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        discard_untrimmed=discard_untrimmed,
        cores=cores,
        report="full",
        json_path=json_path,
        cutadapt_path=cutadapt_path,
    )


__all__ = [
    "cutadapt_full_report",
    "cutadapt_json_report",
    "cutadapt_minimal_report",
    "cutadapt_trim",
]
_register_public_verbose(globals())
