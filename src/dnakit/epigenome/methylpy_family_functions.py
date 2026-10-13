"""methylpy wrapper for the paired-end read conversion that the existing wrappers do not reach.

``convert_paired_fastq_to_bisulfite`` wraps ``methylpy.call_mc_pe.convert_reads_pe``. The package is
imported only when the function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-methylpy]"`` hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.methylpy_functions import (
    _MAX_BUFFER_LINES,
    _METHYLPY_EXTRA,
    _integer,
    _line_count,
    _new_path,
    _path,
    _provenance,
    _upstream_call,
)
from dnakit.exceptions import ConfigurationError


def convert_paired_fastq_to_bisulfite(
    fastq_path: str | Path,
    output_path: str | Path,
    *,
    is_read2: bool = False,
    buffer_line_number: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the reads of one mate file of a paired-end library for bisulfite alignment.

    Calls ``methylpy.call_mc_pe.convert_reads_pe``. For read 1 (``is_read2=False``) every C becomes
    T, as in ``convert_fastq_to_bisulfite``. For read 2 (``is_read2=True``) every G becomes A.
    In both cases the positions
    of the original C (read 1) or G (read 2) bases are encoded in the read name after a ``!``;
    spaces in the name are also replaced by ``!``. The output file must not exist yet.
    """

    source = _path(fastq_path, "fastq_path")
    target = _new_path(output_path, "output_path")
    if not isinstance(is_read2, bool):
        raise ConfigurationError("is_read2 must be a bool.")
    buffer_lines = _integer(buffer_line_number, "buffer_line_number", 4, _MAX_BUFFER_LINES)
    call_mc = require("methylpy.call_mc_pe", _METHYLPY_EXTRA)
    _upstream_call(
        call_mc.convert_reads_pe,
        str(source),
        str(target),
        is_read2=is_read2,
        buffer_line_number=buffer_lines,
    )
    reads = _line_count(target) // 4
    return {
        "status": "completed",
        "operation": "convert_paired_fastq_to_bisulfite",
        "input": str(source),
        "output": str(target),
        "is_read2": is_read2,
        "converted_base": "G to A" if is_read2 else "C to T",
        "reads": reads,
        "provenance": _provenance(
            "G to A conversion of every read 2 base, C to T for read 1; the original positions "
            "are stored in the read name",
            note="reverse-strand reads are converted when mapped, not here",
        ),
    }


__all__ = ["convert_paired_fastq_to_bisulfite"]

_register_public_verbose(globals(), tuple(__all__))
