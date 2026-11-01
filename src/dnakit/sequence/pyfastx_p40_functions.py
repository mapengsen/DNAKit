"""pyfastx FASTA/FASTQ operations as DNAKit functions (pack P40).

The functions call the pyfastx Python package in process on explicit local files (plain or gzip). pyfastx builds
a small ``.fxi`` index file next to the input on first use, exactly as pyfastx does. Outputs go to new paths and
nothing is overwritten. Flanking sequences, subsequence search and FASTQ statistics are already available as
``seqkit_flanking_sequences``, ``locate_sequence_positions`` and ``vsearch_fastq_stats`` / ``summarize_pyfastx``.
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-pyfastx"


def _fasta(path: str | Path) -> tuple[Any, Any, Path]:
    module = require("pyfastx", _EXTRA)
    source = _h.existing(path, "fasta")
    try:
        return module, module.Fasta(str(source)), source
    except (RuntimeError, ValueError, OSError) as exc:
        raise ConfigurationError(f"pyfastx could not open {source} as FASTA: {exc}") from exc


def pyfastx_count_sequences(
    fasta: str | Path,
    *,
    min_length: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the FASTA sequences that are at least ``min_length`` bp long (pyfastx Fasta.count)."""

    module, fa, source = _fasta(fasta)
    cutoff = _h.integer(min_length, "min_length", minimum=0)
    return {"status": "completed", "backend": "pyfastx", "version": module.__version__, "min_length": cutoff,
            "sequence_count": int(fa.count(cutoff)), "total_sequence_count": len(fa)}


def pyfastx_list_sequence_names(
    fasta: str | Path,
    *,
    sort_by: str | None = None,
    reverse: bool = False,
    max_names: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the sequence names of a FASTA file (pyfastx Fasta.keys).

    ``sort_by`` is None (file order), ``"id"``, ``"name"`` or ``"length"`` (pyfastx Identifier.sort).
    """

    module, fa, _source = _fasta(fasta)
    limit = _h.integer(max_names, "max_names", minimum=1, maximum=_h.MAX_ROWS * 10)
    keys = fa.keys()
    if sort_by is not None:
        keys = keys.sort(_h.choice(sort_by, ("id", "name", "length"), "sort_by"), reverse=bool(reverse))
    names = [str(name) for name in keys]
    return {"status": "completed", "backend": "pyfastx", "version": module.__version__,
            "name_count": len(names), "names": names[:limit], "truncated": len(names) > limit}


def pyfastx_nl_statistic(
    fasta: str | Path,
    *,
    quantile: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the N50/L50-style statistic of a FASTA file (pyfastx Fasta.nl).

    ``quantile`` 50 gives N50 and L50; 90 gives N90 and L90. Returns the N length and the L count.
    """

    module, fa, _source = _fasta(fasta)
    q = _h.integer(quantile, "quantile", minimum=1, maximum=99)
    n_length, l_count = fa.nl(q)
    return {"status": "completed", "backend": "pyfastx", "version": module.__version__, "quantile": q,
            "n_length": int(n_length), "l_count": int(l_count), "sequence_count": len(fa)}


def pyfastx_split_fastq(
    fastq: str | Path,
    output_dir: str | Path,
    *,
    file_count: int | None = None,
    reads_per_file: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a FASTQ file into several files (``pyfastx split``).

    Give ``file_count`` for N files of even size or ``reads_per_file`` for files with that many reads. Files are
    named ``<name>.<n><ext>`` (gzip input gives gzip output), as pyfastx does, inside the new ``output_dir``.
    """

    module = require("pyfastx", _EXTRA)
    source = _h.existing(fastq, "fastq")
    if (file_count is None) == (reads_per_file is None):
        raise ConfigurationError("Give exactly one of file_count and reads_per_file.")
    try:
        fq = module.Fastq(str(source))
    except (RuntimeError, ValueError, OSError) as exc:
        raise ConfigurationError(f"pyfastx could not open {source} as FASTQ: {exc}") from exc
    total = len(fq)
    if total == 0:
        raise ConfigurationError("The FASTQ file contains no reads.")
    if file_count is not None:
        parts = _h.integer(file_count, "file_count", minimum=1)
        per = math.ceil(total / parts)
    else:
        per = _h.integer(reads_per_file, "reads_per_file", minimum=1)
        parts = math.ceil(total / per)
    folder = _h.new_dir(output_dir, "output_dir")
    base, ext = Path(source.name).stem, Path(source.name).suffix
    suffix = ext
    if fq.is_gzip:
        base, inner = Path(base).stem, Path(base).suffix
        suffix = inner + ext
    digits = len(str(parts))
    outputs: list[str] = []
    counts: list[int] = []
    handle = None
    written = 0
    number = 0
    try:
        for read in fq:
            if written == 0:
                number += 1
                target = folder / f"{base}.{str(number).zfill(digits)}{suffix}"
                handle = gzip.open(target, "wt") if fq.is_gzip else open(target, "w")
                outputs.append(str(target))
            handle.write(read.raw)
            written += 1
            if written == per:
                handle.close()
                counts.append(written)
                written = 0
                handle = None
        if handle is not None:
            handle.close()
            counts.append(written)
    finally:
        if handle is not None and not handle.closed:
            handle.close()
    return {"status": "completed", "backend": "pyfastx", "version": module.__version__, "read_count": total,
            "file_count": len(outputs), "reads_per_file": counts, "outputs": outputs}


def pyfastx_extract_records(
    source: str | Path,
    output: str | Path,
    names: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract FASTA sequences or FASTQ reads by name into a new file (``pyfastx extract``).

    ``names`` are the sequence or read names (the first word of the header); records are written in the order of
    ``names`` exactly as stored in the input (line wrapping and FASTQ layout preserved). The ``--reverse-complement``
    and ``--out-fasta`` flags of the pyfastx command are accepted by its parser but ignored by pyfastx 2.3, so they
    are not offered here. Names that are not in the file raise a ConfigurationError and nothing is written.
    """

    module = require("pyfastx", _EXTRA)
    path = _h.existing(source, "source")
    target = _h.new_file(output, "output")
    if isinstance(names, str) or not names:
        raise ConfigurationError("names must be a non-empty list of names.")
    wanted = [_h.text(n, "names") for n in names]
    with _h.open_text(path) as handle:
        marker = handle.read(1)
    if marker not in (">", "@"):
        raise ConfigurationError(f"{path} is neither FASTA nor FASTQ.")
    is_fastq = marker == "@"
    try:
        fx = module.Fastq(str(path)) if is_fastq else module.Fasta(str(path))
    except (RuntimeError, ValueError, OSError) as exc:
        raise ConfigurationError(f"pyfastx could not open {path}: {exc}") from exc
    missing = [n for n in wanted if n not in fx]
    if missing:
        raise ConfigurationError(f"These names are not in the file: {', '.join(missing[:10])}")
    with open(target, "w") as out:
        for name in wanted:
            out.write(fx[name].raw)
    return {"status": "completed", "backend": "pyfastx", "version": module.__version__,
            "format": "fastq" if is_fastq else "fasta", "record_count": len(wanted), "outputs": [str(target)]}


__all__ = [
    "pyfastx_count_sequences",
    "pyfastx_extract_records",
    "pyfastx_list_sequence_names",
    "pyfastx_nl_statistic",
    "pyfastx_split_fastq",
]

_register_public_verbose(globals(), tuple(__all__))
