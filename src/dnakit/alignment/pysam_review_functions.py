"""samtools tview text view for alignments, run with the system samtools (pysam review group).

pysam bundles samtools 1.21, and that build has no ``tview`` subcommand, so ``pysam.tview`` fails on
every call. ``pysam_tview_text`` runs the samtools found on PATH (or an explicit path) with fixed
arguments, ``samtools tview -d T -p <chrom:pos> -w <width> <alignment> [reference]``, and writes the
text view to a new file. No shell is used. When samtools is missing, the function raises
``BackendUnavailableError`` with an install hint.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, input_path, run, run_to_file
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_POSITION = re.compile(r"^(?P<chrom>[^:\s]+):(?P<pos>[1-9][0-9]{0,11})$")
_MIN_WIDTH = 20
_MAX_WIDTH = 1000
_INSTALL_HINT = (
    "Install samtools (tested with 1.24; for example conda install -c bioconda samtools, or apt install "
    "samtools) and put it on PATH, or pass samtools='/full/path/to/samtools'. "
    "pysam's bundled samtools has no tview subcommand."
)


def _program(value: object) -> str:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError("samtools must be a program name or a path to the executable.")
    resolved = shutil.which(str(value))
    if resolved is None:
        raise BackendUnavailableError(
            f"samtools executable {str(value)!r} was not found.", hint=_INSTALL_HINT
        )
    return str(Path(resolved).resolve())


def _position(value: object) -> str:
    if not isinstance(value, str):
        raise ConfigurationError("position must be a string such as 'chr1:100'.")
    match = _POSITION.fullmatch(value.strip())
    if match is None:
        raise ConfigurationError(
            f"position must be 'chrom:pos' with a 1-based position; got {value!r}."
        )
    return f"{match['chrom']}:{match['pos']}"


def _width(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError("width must be an integer.")
    if not _MIN_WIDTH <= value <= _MAX_WIDTH:
        raise ConfigurationError(f"width must be between {_MIN_WIDTH} and {_MAX_WIDTH} columns.")
    return value


def _new_output(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("output must be a file path.")
    target = Path(value).expanduser().resolve()
    stderr_log = target.with_name(target.name + ".stderr")
    if target.exists() or stderr_log.exists():
        raise ConfigurationError(f"Output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def pysam_tview_text(
    alignment: str | Path,
    position: str,
    *,
    output: str | Path,
    reference: str | Path | None = None,
    width: int = 80,
    samtools: str | Path = "samtools",
    verbose: bool = True,
) -> dict[str, Any]:
    """Render one screen of a BAM or CRAM alignment as text with the system ``samtools tview``.

    ``position`` is ``chrom:pos`` (1-based) and sets the view centre. ``width`` is the number of
    columns (``-w``). The text view is written to ``output``, a new file, and its lines are returned
    as data. ``reference`` is an optional FASTA for mismatch display. The samtools version used is
    recorded in the result. A missing samtools raises ``BackendUnavailableError``.
    """

    source = input_path(alignment)
    ref = input_path(reference) if reference is not None else None
    locus = _position(position)
    columns = _width(width)
    target = _new_output(output)
    program = _program(samtools)
    version_line = run(program, ["--version"], target.parent, ExecutionConfig(timeout_seconds=60.0))
    version = version_line.splitlines()[0].strip() if version_line.strip() else None
    arguments = ["tview", "-d", "T", "-p", locus, "-w", str(columns), str(source)]
    if ref is not None:
        arguments.append(str(ref))
    log = run_to_file(program, arguments, target, ExecutionConfig(timeout_seconds=300.0))
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    return {
        "status": "ok",
        "backend": "samtools",
        "program": program,
        "version": version,
        "alignment": str(source),
        "reference": str(ref) if ref is not None else None,
        "position": locus,
        "width": columns,
        "output": str(target),
        "line_count": len(lines),
        "lines": lines,
        "log": log[-2000:],
    }


__all__ = ["pysam_tview_text"]

_register_public_verbose(globals())
