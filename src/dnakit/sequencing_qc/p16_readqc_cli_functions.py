"""Long-read filtering, trimming and demultiplexing facades (pack P16): Filtlong, NanoFilt, Porechop.

``filter_long_reads_with_reference`` (Filtlong reference scoring, trimming and splitting), ``filter_nanopore_reads_by_gc_summary``
(NanoFilt GC and sequencing-summary filters), ``filter_nanopore_reads_by_start_time`` (NanoFilt ``scripts/timefilt.py``,
re-implemented in Python), ``trim_nanopore_adapters_discard_middle`` and ``demultiplex_nanopore_barcodes`` (Porechop).
The programs are run without a shell through ``dnakit.bioinformatics._p16_cli``; results are JSON-serialisable dictionaries.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_FILTLONG = cli.Tool("filtlong", "filtlong", "Install Filtlong (conda install -c bioconda filtlong) and put it on PATH or pass filtlong_path.",
                     "0.3.1", "https://github.com/rrwick/Filtlong", ("--version",), r"v(\d+(?:\.\d+)+)")
_NANOFILT = cli.Tool("nanofilt", "NanoFilt", "Install NanoFilt (pip install nanofilt or conda install -c bioconda nanofilt) and put it on PATH or pass nanofilt_path.",
                     "2.8.0", "https://doi.org/10.1093/bioinformatics/bty149", ("--version",), None)
_PORECHOP = cli.Tool("porechop", "porechop", "Install Porechop (conda install -c bioconda porechop) and put it on PATH or pass porechop_path.",
                     "0.2.4", "https://github.com/rrwick/Porechop", ("--version",), None)


def filter_long_reads_with_reference(
    reads: str,
    output_path: str,
    *,
    illumina_reads_1: str | None = None,
    illumina_reads_2: str | None = None,
    assembly: str | None = None,
    trim: bool = False,
    split: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    min_mean_quality: float | None = None,
    target_bases: int | None = None,
    keep_percent: float | None = None,
    filtlong_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score long reads by k-mer matches to a reference, then filter, trim or split them (``filtlong``).

    Exactly one kind of reference may be given: Illumina reads (``illumina_reads_1`` and optionally ``illumina_reads_2``,
    ``--short_1``/``--short_2``) or an ``assembly`` (``--assembly``); read quality is then judged by k-mer identity to it
    instead of the Phred scores. ``trim`` removes non-matching bases from both read ends (``--trim``); ``split`` cuts
    reads at runs of at least that many non-matching bases (``--split``). Both need a reference. Thresholds are the
    Filtlong ones: ``min_length``, ``max_length``, ``min_mean_quality``, ``target_bases``, ``keep_percent``. The kept reads
    are written to ``output_path`` (FASTQ).
    """
    source = cli.existing_file(reads, "reads")
    target = cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    references = [illumina_reads_1 is not None, assembly is not None]
    if all(references) or (illumina_reads_2 is not None and illumina_reads_1 is None):
        raise ConfigurationError("Give illumina_reads_1 (and optionally illumina_reads_2) or assembly, not both.")
    if illumina_reads_1 is not None:
        arguments += ["--short_1", str(cli.existing_file(illumina_reads_1, "illumina_reads_1"))]
    if illumina_reads_2 is not None:
        arguments += ["--short_2", str(cli.existing_file(illumina_reads_2, "illumina_reads_2"))]
    if assembly is not None:
        arguments += ["--assembly", str(cli.existing_file(assembly, "assembly"))]
    has_reference = any(references)
    if cli.flag(trim, "trim"):
        if not has_reference:
            raise ConfigurationError("trim needs a reference (illumina_reads_1 or assembly).")
        arguments.append("--trim")
    if split is not None:
        if not has_reference:
            raise ConfigurationError("split needs a reference (illumina_reads_1 or assembly).")
        arguments += ["--split", str(cli.integer(split, "split", minimum=1))]
    for flag, value, name, lower in (("--min_length", min_length, "min_length", 0), ("--max_length", max_length, "max_length", 1),
                                     ("--target_bases", target_bases, "target_bases", 1)):
        if value is not None:
            arguments += [flag, str(cli.integer(value, name, minimum=lower))]
    if min_mean_quality is not None:
        arguments += ["--min_mean_q", str(cli.number(min_mean_quality, "min_mean_quality", minimum=0.0, maximum=100.0))]
    if keep_percent is not None:
        arguments += ["--keep_percent", str(cli.number(keep_percent, "keep_percent", minimum=0.0, maximum=100.0))]
    arguments.append(str(source))
    record = cli.invoke(_FILTLONG, arguments, cwd=target.parent, stdout_to=target, explicit=filtlong_path)
    kept = sum(1 for index, _ in enumerate(cli.open_text(target)) if index % 4 == 0)
    return cli.finish(record, "filter_long_reads_with_reference", reads_kept=kept, output=str(target))


def filter_nanopore_reads_by_gc_summary(
    reads: str,
    output_path: str,
    *,
    min_gc: float | None = None,
    max_gc: float | None = None,
    sequencing_summary: str | None = None,
    min_quality: float | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    read_type: str | None = None,
    nanofilt_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter long reads by GC content or by the quality in a sequencing-summary file (``NanoFilt``).

    ``min_gc``/``max_gc`` are fractions in 0-1 (``--minGC``/``--maxGC``; NanoFilt ignores them when a summary is used).
    ``sequencing_summary`` (``--summary``) takes read quality from the ONT ``sequencing_summary.txt`` instead of the FASTQ;
    ``read_type`` is 1D, 2D or 1D2. ``min_quality``, ``min_length`` and ``max_length`` are the NanoFilt thresholds. The
    uncompressed FASTQ is read from ``reads`` and the kept reads are written to ``output_path``.
    """
    source = cli.existing_file(reads, "reads")
    target = cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    for flag, value, name in (("--minGC", min_gc, "min_gc"), ("--maxGC", max_gc, "max_gc")):
        if value is not None:
            arguments += [flag, str(cli.number(value, name, minimum=0.0, maximum=1.0))]
    if sequencing_summary is not None:
        arguments += ["--summary", str(cli.existing_file(sequencing_summary, "sequencing_summary"))]
    if read_type is not None:
        arguments += ["--readtype", cli.choice(read_type, ("1D", "2D", "1D2"), "read_type")]
    if min_quality is not None:
        arguments += ["-q", str(cli.number(min_quality, "min_quality", minimum=0.0, maximum=100.0))]
    if min_length is not None:
        arguments += ["-l", str(cli.integer(min_length, "min_length", minimum=0))]
    if max_length is not None:
        arguments += ["--maxlength", str(cli.integer(max_length, "max_length", minimum=1))]
    if not arguments:
        raise ConfigurationError("Give at least one filter (min_gc, max_gc, sequencing_summary, min_quality, ...).")
    arguments.append(str(source))
    record = cli.invoke(_NANOFILT, arguments, cwd=target.parent, stdout_to=target, explicit=nanofilt_path)
    kept = sum(1 for index, _ in enumerate(cli.open_text(target)) if index % 4 == 0)
    return cli.finish(record, "filter_nanopore_reads_by_gc_summary", reads_kept=kept, output=str(target))


def filter_nanopore_reads_by_start_time(
    reads_gz: str,
    output_path: str,
    *,
    time_from: str,
    time_to: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep reads whose sequencing start time lies in a window (NanoFilt ``scripts/timefilt.py``).

    ``reads_gz`` is a gzip-compressed FASTQ whose headers carry ``start_time=<ISO time>`` (as written by MinKNOW/Guppy);
    ``time_from``/``time_to`` are inclusive bounds in any format ``dateutil`` parses; give a time zone (for example
    ``2017-10-13T00:00:00Z``) because read times carry one and naive and aware times cannot be compared. Reads are written unchanged
    (the original rewrites them through Biopython, which gives the same four lines) to ``output_path``; reads
    without ``start_time`` raise an error, as in the script. Needs python-dateutil (extra ``nanopore-time``).
    """
    parser = require("dateutil.parser", "nanopore-time").parse
    source = cli.existing_file(reads_gz, "reads_gz")
    target = cli.new_file(output_path, "output_path")
    lower, upper = parser(cli.text(time_from, "time_from")), parser(cli.text(time_to, "time_to"))
    total = kept = 0
    with gzip.open(source, "rt", encoding="utf-8") as handle, target.open("x", encoding="utf-8") as out:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            if not header.startswith("@") or not quality:
                raise ConfigurationError(f"{source.name} is not a four-line FASTQ file.")
            fields = [item for item in header[1:].split() if item.startswith("start_time")]
            if not fields:
                raise ConfigurationError(f"Read {header.strip()!r} has no start_time field.")
            total += 1
            if lower <= parser(fields[0].split("=")[1]) <= upper:
                kept += 1
                out.write(header.rstrip("\n") + "\n" + sequence.rstrip("\n") + "\n+\n" + quality.rstrip("\n") + "\n")
    return {
        "status": "completed", "backend": "nanofilt-scripts", "operation": "filter_nanopore_reads_by_start_time",
        "reads_total": total, "reads_kept": kept, "outputs": [str(target)],
        "provenance": {"tool": "NanoFilt scripts/timefilt.py", "reimplemented_in": "dnakit"},
    }


def trim_nanopore_adapters_discard_middle(
    reads: str,
    output_path: str,
    *,
    discard_middle: bool = True,
    threads: int = 1,
    porechop_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim adapters and discard reads with an adapter in the middle (``porechop --discard_middle``).

    Reads with an internal adapter are chimeras; Porechop's default splits them, ``discard_middle=True`` (this
    function's default) throws them away instead. ``False`` reproduces Porechop's default behaviour (split).
    """
    source = cli.existing_file(reads, "reads")
    target = cli.new_file(output_path, "output_path")
    arguments = ["-i", str(source), "-o", str(target), "-t", str(cli.threads_value(threads))]
    if cli.flag(discard_middle, "discard_middle"):
        arguments.append("--discard_middle")
    record = cli.invoke(_PORECHOP, arguments, cwd=target.parent, expected=[target], explicit=porechop_path)
    return cli.finish(record, "trim_nanopore_adapters_discard_middle", discard_middle=discard_middle)


def demultiplex_nanopore_barcodes(
    reads: str,
    output_directory: str,
    *,
    require_two_barcodes: bool = False,
    barcode_threshold: float | None = None,
    discard_unassigned: bool = False,
    threads: int = 1,
    porechop_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin reads by native/rapid barcode while trimming adapters (``porechop --barcode_dir``).

    One FASTQ per barcode (``BC01.fastq`` ...) and ``none.fastq`` for unassigned reads is written to the new
    ``output_directory``. ``require_two_barcodes`` (``--require_two_barcodes``) demands a barcode at both ends;
    ``barcode_threshold`` is ``--barcode_threshold`` (percent identity, default 75); ``discard_unassigned`` is
    ``--discard_unassigned``.
    """
    source = cli.existing_file(reads, "reads")
    out_dir = Path(output_directory).expanduser().resolve()
    if out_dir.exists():
        raise ConfigurationError(f"output_directory must not exist yet: {out_dir}")
    arguments = ["-i", str(source), "-b", str(out_dir), "-t", str(cli.threads_value(threads))]
    if cli.flag(require_two_barcodes, "require_two_barcodes"):
        arguments.append("--require_two_barcodes")
    if barcode_threshold is not None:
        arguments += ["--barcode_threshold", str(cli.number(barcode_threshold, "barcode_threshold", minimum=0.0, maximum=100.0))]
    if cli.flag(discard_unassigned, "discard_unassigned"):
        arguments.append("--discard_unassigned")
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    record = cli.invoke(_PORECHOP, arguments, cwd=out_dir.parent, expected=[out_dir], explicit=porechop_path)
    bins = {p.name: sum(1 for i, _ in enumerate(p.open()) if i % 4 == 0) for p in sorted(out_dir.glob("*.fastq*"))}
    return cli.finish(record, "demultiplex_nanopore_barcodes", reads_per_bin=bins, output_directory=str(out_dir))


__all__ = [
    "demultiplex_nanopore_barcodes",
    "filter_long_reads_with_reference",
    "filter_nanopore_reads_by_gc_summary",
    "filter_nanopore_reads_by_start_time",
    "trim_nanopore_adapters_discard_middle",
]

_register_public_verbose(globals(), tuple(__all__))
