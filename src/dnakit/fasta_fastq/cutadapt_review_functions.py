"""cutadapt read modifications that cutadapt_trim does not expose, as command-line facades.

``cutadapt_trim`` covers 3' / 5' / anywhere adapter trimming, quality trimming and length filters. The cutadapt
options below are the other read operations of the same program: a fixed-length cut from either read end (``-u``),
shortening to a given length (``-l``), the action taken on an adapter match (``--action``, ``-n``, ``--rc``) and
removal of N bases from both read ends (``--trim-n``). The user installs cutadapt (``pip install cutadapt`` or
``conda install -c bioconda cutadapt``). Each function runs the cutadapt program with an argument list built from
typed parameters mapped to fixed flags; no shell is used; outputs are new explicit paths; results are
JSON-serialisable dictionaries with status, backend, version, command, outputs, the one-line report and a bounded
log tail.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq.cutadapt_gap_functions import _adapter, _existing, _new, _program, _version

_MAX_TAIL = 4000
_TIMEOUT_SECONDS = 3600.0
_ACTIONS = ("trim", "retain", "mask", "lowercase", "crop", "none")


def _integer(value: object, name: str, *, minimum: int | None = None, nonzero: bool = False) -> int:
    """Return value when it is an int that satisfies the bounds; raise ConfigurationError otherwise."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if nonzero and value == 0:
        raise ConfigurationError(f"{name} must not be 0.")
    return value


def _execute(source: Path, target: Path, args: list[str], cutadapt_path: str | None) -> dict[str, Any]:
    """Run cutadapt once with a fixed argument list and return the JSON-serialisable run record."""
    program = _program(cutadapt_path)
    full = ["--report", "minimal", *args, "-o", str(target), str(source)]
    completed = subprocess.run(
        [program, *full],
        capture_output=True,
        text=True,
        timeout=_TIMEOUT_SECONDS,
        shell=False,
        check=False,
        cwd=target.parent,
    )
    produced = target.is_file()
    return {
        "status": "completed" if completed.returncode == 0 and produced else "failed",
        "backend": "cutadapt",
        "version": _version(program),
        "command": ["cutadapt", *full],
        "return_code": completed.returncode,
        "outputs": [str(target)] if produced else [],
        "report": completed.stdout.strip(),
        "log_tail": (completed.stderr or "")[-_MAX_TAIL:],
    }


def cutadapt_unconditional_cut(
    input_path: str,
    output_path: str,
    *,
    cut_5p: int = 0,
    cut_3p: int = 0,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove a fixed number of bases from the 5' and/or the 3' end of every read (cutadapt -u N, then -u -N)."""
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    five = _integer(cut_5p, "cut_5p", minimum=0)
    three = _integer(cut_3p, "cut_3p", minimum=0)
    if five == 0 and three == 0:
        raise ConfigurationError("Give cut_5p or cut_3p (or both) as a positive number of bases.")
    args: list[str] = []
    if five:
        args += ["-u", str(five)]
    if three:
        args += ["-u", str(-three)]
    return _execute(source, target, args, cutadapt_path)


def cutadapt_shorten(
    input_path: str,
    output_path: str,
    *,
    length: int,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shorten every read to length bases: a positive length keeps the first bases, a negative one the last (cutadapt -l)."""
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    value = _integer(length, "length", nonzero=True)
    return _execute(source, target, ["-l", str(value)], cutadapt_path)


def cutadapt_adapter_action(
    input_path: str,
    output_path: str,
    *,
    adapter_3p: str | None = None,
    adapter_5p: str | None = None,
    anywhere_adapter: str | None = None,
    action: str = "trim",
    times: int = 1,
    reverse_complement: bool = False,
    error_rate: float = 0.1,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find adapters in every read and apply one action to each match: trim, retain, mask with N, lowercase, crop or none.

    ``reverse_complement`` (cutadapt --rc) also checks the reverse complement of each read; a read written in its
    reverse-complemented form gets the cutadapt suffix " rc" appended to its name.
    """
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    adapters = [
        (flag, _adapter(value, name))
        for flag, value, name in (
            ("-a", adapter_3p, "adapter_3p"),
            ("-g", adapter_5p, "adapter_5p"),
            ("-b", anywhere_adapter, "anywhere_adapter"),
        )
        if value is not None
    ]
    if not adapters:
        raise ConfigurationError("Give at least one of adapter_3p, adapter_5p or anywhere_adapter.")
    if action not in _ACTIONS:
        raise ConfigurationError(f"action must be one of {_ACTIONS}.")
    count = _integer(times, "times", minimum=1)
    if isinstance(error_rate, bool) or not isinstance(error_rate, (int, float)) or not 0 <= error_rate <= 1:
        raise ConfigurationError("error_rate must be a number between 0 and 1.")
    args: list[str] = ["-e", str(float(error_rate)), "--action", action, "-n", str(count)]
    if reverse_complement:
        args.append("--rc")
    for flag, value in adapters:
        args += [flag, value]
    return _execute(source, target, args, cutadapt_path)


def cutadapt_trim_n_ends(
    input_path: str,
    output_path: str,
    *,
    cutadapt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove N bases from both ends of every read (cutadapt --trim-n)."""
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    return _execute(source, target, ["--trim-n"], cutadapt_path)


__all__ = [
    "cutadapt_adapter_action",
    "cutadapt_shorten",
    "cutadapt_trim_n_ends",
    "cutadapt_unconditional_cut",
]
_register_public_verbose(globals())
