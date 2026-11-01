"""pyfaidx sequence-region operations as DNAKit functions: fetch, write, transform, mask and region
parsing.

The ``faidx`` program that pyfaidx installs runs the FASTA operations (``python -m pip install
pyfaidx``, or ``conda install -c bioconda pyfaidx``). Each function maps its keyword arguments to
fixed ``faidx`` flags; no shell is used, outputs are new explicit paths, and results are
JSON-serialisable dictionaries with status, backend, version, command, outputs and a bounded log
tail. The wrapper and the region parsers call pyfaidx's own functions in-process. pyfaidx is
imported only when a function runs; a missing package raises BackendUnavailableError with the
install hint.

Sequence output is streamed from faidx's standard output into the new file, not written with ``faidx
-o``: in pyfaidx 0.8.2 the ``-o`` handle is never closed, and the tail of a larger output can be
lost (reproduction in the manifest notes). faidx writes a ``.fai`` index next to a FASTA that has
none (or whose index is out of date unless ``no_rebuild`` is set). ``faidx_mask_sequences`` masks a
copy and never changes its input.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendTimeoutError, ConfigurationError

_HINT = (
    "Install pyfaidx 0.8 or later (python -m pip install pyfaidx, or conda "
    "install -c bioconda pyfaidx)."
)
_TRANSFORMS = ("bed", "chromsizes", "nucleotide", "transposed")
_TIMEOUT_SECONDS = 3600.0
_LOG_TAIL_CHARS = 4000
_VERSION = re.compile(r"\d+(?:\.\d+)+")


def _pyfaidx() -> Any:
    return require("pyfaidx", "bio-pyfaidx")


def _selection(regions: Sequence[str], bed: str | Path | None) -> tuple[list[str], list[str]]:
    """Option flags and positional arguments that select regions; a bed file and region strings are
    exclusive.
    """
    if bed is not None and len(regions) > 0:
        raise ConfigurationError(
            "Give either regions or a bed file, not both: faidx ignores the regions when -b is set."
        )
    if bed is not None:
        return ["-b", str(_cli.existing_file(bed, "bed"))], []
    return [], [_cli.text(region, "region") for region in regions]


def _output_flags(
    *,
    complement: bool,
    reverse: bool,
    auto_strand: bool,
    no_names: bool,
    long_names: bool,
    no_coords: bool,
    lazy: bool,
    default_seq: str | None,
    delimiter: str | None,
    no_rebuild: bool,
) -> list[str]:
    if no_names and long_names:
        raise ConfigurationError(
            "no_names and long_names are mutually exclusive (faidx -n and -f)."
        )
    flags: list[str] = []
    for enabled, flag in (
        (complement, "-c"),
        (reverse, "-r"),
        (auto_strand, "-y"),
        (no_names, "-n"),
        (long_names, "-f"),
        (no_coords, "-t"),
        (lazy, "-l"),
        (no_rebuild, "--no-rebuild"),
    ):
        if enabled:
            flags.append(flag)
    checked = faidx_check_default_seq(default_seq, verbose=False)
    if checked is not None:
        flags += ["-s", checked]
    if delimiter is not None:
        flags += ["-d", _cli.text(delimiter, "delimiter")]
    return flags


def faidx_check_default_seq(value: str | None, *, verbose: bool = True) -> str | None:
    """Check the default base for missing positions and masking: None, or exactly one character.

    This is the rule of pyfaidx.cli.check_seq_length (``--default-seq``). A dash is allowed; faidx
    accepts it as a value.
    """
    if value is None:
        return None
    if not isinstance(value, str) or len(value) != 1 or value == "\x00":
        raise ConfigurationError("default_seq must be a single character.")
    return value


def _size_range_flags(size_range: Sequence[int] | None) -> list[str]:
    if size_range is None:
        return []
    if len(size_range) != 2:
        raise ConfigurationError("size_range must be a pair (low, high).")
    low = _cli.integer(size_range[0], "size_range low", minimum=0)
    high = _cli.integer(size_range[1], "size_range high", minimum=0)
    return ["-a", f"{low},{high}"]


@lru_cache(maxsize=16)
def _version_of(executable: str) -> str | None:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            shell=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = _VERSION.search(completed.stdout or "")
    return match.group(0) if match else None


def _run(
    arguments: list[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    faidx_path: str | None,
) -> dict[str, Any]:
    """Run faidx for operations whose outputs are files written by faidx itself (split files, masked
    copy).
    """
    return _cli.run("faidx", arguments, cwd=cwd, expected=expected, explicit=faidx_path, hint=_HINT)


def _run_to_file(arguments: list[str], *, out: Path, faidx_path: str | None) -> dict[str, Any]:
    """Run faidx with its standard output streamed into the new file ``out`` (no ``-o``; see module
    notes).

    stderr goes to a temporary file whose tail is returned; the command runs without a shell and is
    killed after the time limit. ``status`` is completed only for exit code 0.
    """
    executable = _cli.program("faidx", explicit=faidx_path, hint=_HINT)
    with out.open("xb") as stdout, tempfile.TemporaryFile() as stderr:
        process = subprocess.Popen(
            [executable, *arguments],
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            cwd=str(out.parent),
            shell=False,
        )
        try:
            return_code = process.wait(timeout=_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired as exc:
            process.kill()
            process.wait()
            raise BackendTimeoutError("faidx exceeded its time limit.") from exc
        stderr.seek(0)
        log_tail = stderr.read()[-_LOG_TAIL_CHARS:].decode("utf-8", errors="replace")
    complete = return_code == 0 and out.is_file()
    return {
        "status": "completed" if complete else "failed",
        "backend": "faidx",
        "version": _version_of(executable),
        "command": ["faidx", *arguments],
        "return_code": return_code,
        "outputs": [str(out)] if out.is_file() else [],
        "log_tail": log_tail,
    }


def faidx_wrap_sequence(sequence: str, width: int, *, verbose: bool = True) -> list[str]:
    """Wrap a sequence into lines of ``width`` characters, each line ending with a newline
    (pyfaidx.wrap_sequence).
    """
    if not isinstance(sequence, str):
        raise ConfigurationError("sequence must be a string.")
    checked = _cli.integer(width, "width", minimum=1)
    return list(_pyfaidx().wrap_sequence(checked, sequence))


def faidx_parse_ucsc_region(region: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a UCSC region ``name:start-end`` (1-based, inclusive) into name, 0-based start and
    exclusive end.

    A bare name gives start and end None (the whole sequence). This is pyfaidx.ucsc_split.
    """
    name, start, end = _pyfaidx().ucsc_split(_cli.text(region, "region"))
    return {"name": name, "start": start, "end": end}


def faidx_parse_bed_entry(line: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse one BED line into chromosome, 0-based start and end; comments (#) give None values.

    This is pyfaidx.bed_split. A malformed line raises ConfigurationError.
    """
    if not isinstance(line, str):
        raise ConfigurationError("line must be a string.")
    try:
        name, start, end = _pyfaidx().bed_split(line)
    except (ValueError, IndexError) as exc:
        raise ConfigurationError(f"Malformed BED entry: {line.rstrip()!r} ({exc})") from exc
    return {"name": name, "start": start, "end": end}


def faidx_split_regions(
    regions: Sequence[str] = (),
    *,
    bed: str | Path | None = None,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """Turn region strings (``name:start-end``) or a BED file into name/start/end records, in input
    order.

    Comment lines of a BED file are dropped. With neither regions nor a bed file the result is empty
    (faidx then means "all sequences", which is applied by the FASTA functions, not here).
    """
    if bed is not None and len(regions) > 0:
        raise ConfigurationError("Give either regions or a bed file, not both.")
    if bed is not None:
        path = _cli.existing_file(bed, "bed")
        records: list[dict[str, Any]] = []
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                record = faidx_parse_bed_entry(line, verbose=False)
                if record["name"] is not None:
                    records.append(record)
        return records
    return [faidx_parse_ucsc_region(region, verbose=False) for region in regions]


def faidx_fetch_sequences(
    fasta: str | Path,
    output: str | Path,
    *,
    regions: Sequence[str] = (),
    bed: str | Path | None = None,
    complement: bool = False,
    reverse: bool = False,
    auto_strand: bool = False,
    no_names: bool = False,
    long_names: bool = False,
    no_coords: bool = False,
    lazy: bool = False,
    default_seq: str | None = None,
    delimiter: str | None = None,
    no_rebuild: bool = False,
    faidx_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract FASTA sequences for regions (``name:start-end`` or a BED file) into a new FASTA.

    Without regions every sequence is written. A region missing from the FASTA is reported by faidx
    and produces no record. ``auto_strand`` reverse-complements regions whose start is greater than
    their end.
    """
    source = _cli.existing_file(fasta, "fasta")
    out = _cli.new_file(output, "output")
    selection, positional = _selection(regions, bed)
    flags = _output_flags(
        complement=complement,
        reverse=reverse,
        auto_strand=auto_strand,
        no_names=no_names,
        long_names=long_names,
        no_coords=no_coords,
        lazy=lazy,
        default_seq=default_seq,
        delimiter=delimiter,
        no_rebuild=no_rebuild,
    )
    arguments = [*flags, *selection, str(source), *positional]
    return _run_to_file(arguments, out=out, faidx_path=faidx_path)


def faidx_write_sequences(
    fasta: str | Path,
    *,
    output: str | Path | None = None,
    split_directory: str | Path | None = None,
    regions: Sequence[str] = (),
    bed: str | Path | None = None,
    regex: str = ".*",
    invert_match: bool = False,
    size_range: Sequence[int] | None = None,
    complement: bool = False,
    reverse: bool = False,
    auto_strand: bool = False,
    no_names: bool = False,
    long_names: bool = False,
    no_coords: bool = False,
    lazy: bool = False,
    default_seq: str | None = None,
    delimiter: str | None = None,
    no_rebuild: bool = False,
    faidx_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write selected FASTA sequences by name regex and length range, to one file or one file per
    sequence (faidx).

    ``regex`` selects sequence names (``invert_match`` selects the others); ``size_range`` = (low,
    high) keeps sequences whose length, or region length, lies in the inclusive range. Give exactly
    one of ``output`` (a new file) or ``split_directory`` (a new directory that receives one file
    per selected sequence, named after the sequence and its coordinates). The result lists the files
    written.
    """
    if (output is None) == (split_directory is None):
        raise ConfigurationError("Give exactly one of output or split_directory.")
    source = _cli.existing_file(fasta, "fasta")
    selection, positional = _selection(regions, bed)
    flags = _output_flags(
        complement=complement,
        reverse=reverse,
        auto_strand=auto_strand,
        no_names=no_names,
        long_names=long_names,
        no_coords=no_coords,
        lazy=lazy,
        default_seq=default_seq,
        delimiter=delimiter,
        no_rebuild=no_rebuild,
    )
    filters = ["-g", _cli.text(regex, "regex")]
    if invert_match:
        filters.append("-v")
    filters += _size_range_flags(size_range)
    if split_directory is not None:
        directory = _cli.new_directory(split_directory, "split_directory")
        arguments = ["-x", *flags, *filters, *selection, str(source), *positional]
        result = _run(arguments, cwd=directory, expected=[], faidx_path=faidx_path)
        if result["status"] == "completed":
            result["outputs"] = sorted(str(path) for path in directory.iterdir() if path.is_file())
        return result
    if output is None:  # pragma: no cover - excluded by the exactly-one check above
        raise ConfigurationError("Give exactly one of output or split_directory.")
    out = _cli.new_file(output, "output")
    arguments = [*flags, *filters, *selection, str(source), *positional]
    return _run_to_file(arguments, out=out, faidx_path=faidx_path)


def faidx_transform_sequences(
    fasta: str | Path,
    output: str | Path,
    transform: str,
    *,
    regions: Sequence[str] = (),
    bed: str | Path | None = None,
    complement: bool = False,
    reverse: bool = False,
    delimiter: str | None = None,
    no_rebuild: bool = False,
    faidx_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe regions as a table instead of FASTA (faidx -i): bed, chromsizes, nucleotide or
    transposed.

    ``bed`` writes name, 0-based start, end; ``chromsizes`` writes name and length; ``nucleotide``
    writes counts of A, T, C, G, N and other characters (with a header line); ``transposed`` writes
    name, start, end and the sequence.
    """
    kind = _cli.choice(transform, _TRANSFORMS, "transform")
    source = _cli.existing_file(fasta, "fasta")
    out = _cli.new_file(output, "output")
    selection, positional = _selection(regions, bed)
    flags = _output_flags(
        complement=complement,
        reverse=reverse,
        auto_strand=False,
        no_names=False,
        long_names=False,
        no_coords=False,
        lazy=False,
        default_seq=None,
        delimiter=delimiter,
        no_rebuild=no_rebuild,
    )
    arguments = ["-i", kind, *flags, *selection, str(source), *positional]
    return _run_to_file(arguments, out=out, faidx_path=faidx_path)


def faidx_mask_sequences(
    fasta: str | Path,
    output: str | Path,
    *,
    regions: Sequence[str] = (),
    bed: str | Path | None = None,
    default_seq: str = "N",
    delimiter: str | None = None,
    no_rebuild: bool = False,
    faidx_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the bases of the given regions with ``default_seq`` (faidx -m) in a new copy of the
    FASTA.

    The input is copied to ``output`` first and only the copy is changed. Regions are required:
    faidx masks nothing without them. Lowercase masking (faidx -M) is not offered because pyfaidx
    0.8.2 fails on it (see the manifest).
    """
    source = _cli.existing_file(fasta, "fasta")
    out = _cli.new_file(output, "output")
    if bed is None and len(regions) == 0:
        raise ConfigurationError("faidx_mask_sequences needs at least one region or a bed file.")
    selection, positional = _selection(regions, bed)
    mask_seq = faidx_check_default_seq(default_seq, verbose=False)
    if mask_seq is None:
        raise ConfigurationError("default_seq is required to mask sequences.")
    flags = _output_flags(
        complement=False,
        reverse=False,
        auto_strand=False,
        no_names=False,
        long_names=False,
        no_coords=False,
        lazy=False,
        default_seq=mask_seq,
        delimiter=delimiter,
        no_rebuild=no_rebuild,
    )
    shutil.copyfile(source, out)
    arguments = ["-m", *flags, *selection, str(out), *positional]
    return _run(arguments, cwd=out.parent, expected=[out], faidx_path=faidx_path)


__all__ = [
    "faidx_check_default_seq",
    "faidx_fetch_sequences",
    "faidx_mask_sequences",
    "faidx_parse_bed_entry",
    "faidx_parse_ucsc_region",
    "faidx_split_regions",
    "faidx_transform_sequences",
    "faidx_wrap_sequence",
    "faidx_write_sequences",
]
_register_public_verbose(globals())
