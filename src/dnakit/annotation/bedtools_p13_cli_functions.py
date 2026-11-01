"""BEDTools command-line facades: window, random, shuffle, sample, multiinter, unionbedg, pairtopair, multicov,
reldist, coverage (histogram or per-position depth), bamtobed -bedpe and split.

The user installs BEDTools (``conda install -c bioconda bedtools``, or ``apt install bedtools``); 2.27.1 is the tested
version. Each function runs one ``bedtools`` subcommand with an explicit typed signature mapped to a fixed set of flags;
no shell is used and the program is resolved from an explicit path, then PATH. BEDTools writes its result to standard
output, which is streamed to a **new** output file (an existing file is never overwritten). Results are JSON-serialisable
dictionaries: ``status``, ``backend``, ``version``, ``command``, ``return_code``, ``outputs``, ``log_tail`` and the
number of output records.

The defaults are BEDTools' own (window 1000 bp, 100 bp random intervals, 1,000,000 random or sampled records, overlap
fraction 1E-9 and so on). Inputs that BEDTools requires to be sorted (``multiinter``, ``unionbedg``) are passed through
unchanged, so sort them first.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p13_cli as _p13
from dnakit.exceptions import ConfigurationError

_PROGRAM = "bedtools"
_HINT = (
    "Install BEDTools 2.27.1 or later (for example `conda install -c bioconda bedtools` or "
    "`apt install bedtools`), then put it on PATH or pass bedtools_path."
)
_STRAND_MODES = (None, "same", "opposite")
_PAIR_TYPES = ("neither", "either", "both", "notboth")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fraction(value: float, name: str) -> float:
    return _cli.number(value, name, minimum=0.0, maximum=1.0)


def _flag(on: bool, flag: str) -> list[str]:
    if not isinstance(on, bool):
        raise ConfigurationError(f"{flag} option must be True or False.")
    return [flag] if on else []


def _strand(value: str | None, same: str, opposite: str) -> list[str]:
    if value not in _STRAND_MODES:
        raise ConfigurationError("strand must be None, 'same' or 'opposite'.")
    return [same] if value == "same" else [opposite] if value == "opposite" else []


def _existing_list(paths: Sequence[str], name: str, minimum: int) -> list[Path]:
    if isinstance(paths, str) or not isinstance(paths, Sequence) or len(paths) < minimum:
        raise ConfigurationError(f"{name} must be a list of at least {minimum} existing files.")
    return [_cli.existing_file(path, name) for path in paths]


def _names(names: Sequence[str] | None, count: int) -> list[str]:
    if names is None:
        return []
    if isinstance(names, str) or len(names) != count:
        raise ConfigurationError(f"names must give exactly one name per input file ({count}).")
    return ["-names", *[_cli.text(name, "names") for name in names]]


def _stream(
    subcommand: str,
    arguments: Sequence[str],
    output: Path,
    bedtools_path: str | None,
    timeout_seconds: float,
    *,
    skip_prefixes: Sequence[str] = (),
) -> dict[str, Any]:
    """Run ``bedtools <subcommand> <arguments>`` streaming standard output to ``output``."""
    result = _p13.run_to_file(
        _PROGRAM,
        [subcommand, *arguments],
        output=output,
        explicit=bedtools_path,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = f"bedtools {subcommand}"
    result["record_count"] = _p13.count_lines(output, skip_prefixes=skip_prefixes) if output.is_file() else 0
    return result


def _check_program(bedtools_path: str | None) -> None:
    """Resolve the program before any output file is created, so a missing program leaves nothing behind."""
    _cli.program(_PROGRAM, explicit=bedtools_path, hint=_HINT)


# ---------------------------------------------------------------------------
# window
# ---------------------------------------------------------------------------


def bedtools_window(
    a_path: str,
    b_path: str,
    output_path: str,
    *,
    window: int | None = None,
    left: int | None = None,
    right: int | None = None,
    strand_aware_window: bool = False,
    strand: str | None = None,
    report: str = "overlaps",
    header: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report the B features that fall in a window around each A interval (``bedtools window``).

    ``window`` (``-w``) adds the same number of bases on both sides, ``left``/``right`` (``-l``/``-r``) add an
    asymmetric window; with none of them BEDTools uses 1000 bp. ``strand_aware_window`` (``-sw``) defines left and
    right by the strand of A. ``strand`` is ``None``, ``"same"`` (``-sm``) or ``"opposite"`` (``-Sm``). ``report`` is
    ``"overlaps"`` (every A/B pair, the default), ``"unique"`` (``-u``, each A once when any B overlaps),
    ``"count"`` (``-c``, A plus the number of B hits; this is the "count B features within a window" operation) or
    ``"none"`` (``-v``, only A intervals without a B hit).
    """
    a = _cli.existing_file(a_path, "a_path")
    b = _cli.existing_file(b_path, "b_path")
    if window is not None and (left is not None or right is not None):
        raise ConfigurationError("Give either window, or left and/or right, not both.")
    mode = _cli.choice(report, ("overlaps", "unique", "count", "none"), "report")
    args = ["-a", str(a), "-b", str(b)]
    args += _cli.option("-w", None if window is None else _cli.integer(window, "window", minimum=0))
    args += _cli.option("-l", None if left is None else _cli.integer(left, "left", minimum=0))
    args += _cli.option("-r", None if right is None else _cli.integer(right, "right", minimum=0))
    args += _flag(strand_aware_window, "-sw") + _strand(strand, "-sm", "-Sm")
    args += {"overlaps": [], "unique": ["-u"], "count": ["-c"], "none": ["-v"]}[mode]
    args += _flag(header, "-header")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("window", args, out, bedtools_path, timeout_seconds)


# ---------------------------------------------------------------------------
# random, shuffle, sample
# ---------------------------------------------------------------------------


def bedtools_random(
    genome_path: str,
    output_path: str,
    *,
    length: int = 100,
    count: int = 1_000_000,
    seed: int | None = None,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate random intervals of a fixed length across a genome (``bedtools random``).

    ``genome_path`` is a two-column ``<chrom><TAB><size>`` file. ``length`` (``-l``, default 100) and ``count``
    (``-n``, default 1,000,000) are BEDTools' defaults; ``seed`` (``-seed``) makes the result reproducible.
    """
    genome = _cli.existing_file(genome_path, "genome_path")
    args = ["-g", str(genome), "-l", str(_cli.integer(length, "length", minimum=1))]
    args += ["-n", str(_cli.integer(count, "count", minimum=1))]
    args += _cli.option("-seed", None if seed is None else _cli.integer(seed, "seed", minimum=0))
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("random", args, out, bedtools_path, timeout_seconds)


def bedtools_shuffle(
    intervals_path: str,
    genome_path: str,
    output_path: str,
    *,
    exclude_path: str | None = None,
    include_path: str | None = None,
    same_chromosome: bool = False,
    seed: int | None = None,
    max_exclude_overlap: float | None = None,
    choose_chromosome_first: bool = False,
    bedpe: bool = False,
    max_tries: int | None = None,
    no_overlapping: bool = False,
    allow_beyond_chrom_end: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly permute the genomic locations of intervals (``bedtools shuffle``).

    ``exclude_path`` (``-excl``) lists regions where intervals must not land and ``include_path`` (``-incl``) the only
    regions where they may land (not combinable with ``exclude_path`` or ``max_exclude_overlap``).
    ``same_chromosome`` (``-chrom``) keeps each interval on its chromosome; ``choose_chromosome_first`` is
    ``-chromFirst``; ``max_exclude_overlap`` is ``-f``; ``bedpe`` reads a BEDPE file; ``max_tries`` is ``-maxTries``.
    """
    intervals = _cli.existing_file(intervals_path, "intervals_path")
    genome = _cli.existing_file(genome_path, "genome_path")
    if exclude_path is not None and include_path is not None:
        raise ConfigurationError("exclude_path and include_path cannot be combined.")
    if include_path is not None and max_exclude_overlap is not None:
        raise ConfigurationError("max_exclude_overlap cannot be used with include_path.")
    args = ["-i", str(intervals), "-g", str(genome)]
    if exclude_path is not None:
        args += ["-excl", str(_cli.existing_file(exclude_path, "exclude_path"))]
    if include_path is not None:
        args += ["-incl", str(_cli.existing_file(include_path, "include_path"))]
    args += _flag(same_chromosome, "-chrom")
    args += _cli.option("-seed", None if seed is None else _cli.integer(seed, "seed", minimum=0))
    args += _cli.option(
        "-f", None if max_exclude_overlap is None else _fraction(max_exclude_overlap, "max_exclude_overlap")
    )
    args += _flag(choose_chromosome_first, "-chromFirst") + _flag(bedpe, "-bedpe")
    args += _cli.option("-maxTries", None if max_tries is None else _cli.integer(max_tries, "max_tries", minimum=1))
    args += _flag(no_overlapping, "-noOverlapping") + _flag(allow_beyond_chrom_end, "-allowBeyondChromEnd")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("shuffle", args, out, bedtools_path, timeout_seconds)


def bedtools_sample(
    input_path: str,
    output_path: str,
    *,
    count: int = 1_000_000,
    seed: int | None = None,
    strand: str | None = None,
    header: bool = False,
    bam_as_bed: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Take a random sample of records from a BED/GFF/VCF or BAM file by reservoir sampling (``bedtools sample``).

    ``count`` is ``-n`` (default 1,000,000), ``seed`` is ``-seed``, ``strand`` is ``None``, ``"forward"`` or
    ``"reverse"`` (``-s``), ``bam_as_bed`` writes BAM input as BED (``-bed``).
    """
    source = _cli.existing_file(input_path, "input_path")
    args = ["-i", str(source), "-n", str(_cli.integer(count, "count", minimum=1))]
    args += _cli.option("-seed", None if seed is None else _cli.integer(seed, "seed", minimum=0))
    if strand is not None:
        args += ["-s", _cli.choice(strand, ("forward", "reverse"), "strand")]
    args += _flag(header, "-header") + _flag(bam_as_bed, "-bed")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("sample", args, out, bedtools_path, timeout_seconds)


# ---------------------------------------------------------------------------
# multiinter, unionbedg
# ---------------------------------------------------------------------------


def _multi_file_options(
    paths: Sequence[str],
    names: Sequence[str] | None,
    header: bool,
    genome_path: str | None,
    report_empty: bool,
    filler: str | None,
) -> list[str]:
    files = _existing_list(paths, "input_paths", 2)
    args = ["-i", *[str(path) for path in files], *_names(names, len(files))]
    args += _flag(header, "-header")
    if genome_path is not None:
        args += ["-g", str(_cli.existing_file(genome_path, "genome_path"))]
    if report_empty and genome_path is None:
        raise ConfigurationError("report_empty needs genome_path.")
    args += _flag(report_empty, "-empty")
    if filler is not None:
        args += ["-filler", _cli.text(filler, "filler")]
    return args


def bedtools_multiinter(
    input_paths: Sequence[str],
    output_path: str,
    *,
    names: Sequence[str] | None = None,
    header: bool = False,
    cluster: bool = False,
    genome_path: str | None = None,
    report_empty: bool = False,
    filler: str | None = None,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Identify the common intervals among several sorted BED/GFF/VCF files (``bedtools multiinter``).

    ``names`` gives one name per file (printed in the header), ``cluster`` is Ryan Layer's clustering (``-cluster``),
    ``genome_path`` with ``report_empty`` reports regions covered by no file (``-g``/``-empty``) and ``filler`` is the
    text for absent values (default ``0``). Every input must be sorted by chromosome and start.
    """
    args = _multi_file_options(input_paths, names, header, genome_path, report_empty, filler)
    args += _flag(cluster, "-cluster")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("multiinter", args, out, bedtools_path, timeout_seconds)


def bedtools_unionbedg(
    input_paths: Sequence[str],
    output_path: str,
    *,
    names: Sequence[str] | None = None,
    header: bool = False,
    genome_path: str | None = None,
    report_empty: bool = False,
    filler: str | None = None,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine several sorted bedGraph files into one table with a value column per file (``bedtools unionbedg``).

    ``names`` gives one name per file (printed in the header), ``genome_path`` with ``report_empty`` reports regions
    covered by no file and ``filler`` is the text for absent values (default ``0``).
    """
    args = _multi_file_options(input_paths, names, header, genome_path, report_empty, filler)
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("unionbedg", args, out, bedtools_path, timeout_seconds)


# ---------------------------------------------------------------------------
# pairtopair, multicov
# ---------------------------------------------------------------------------


def bedtools_pairtopair(
    a_path: str,
    b_path: str,
    output_path: str,
    *,
    min_overlap: float | None = None,
    overlap_type: str = "both",
    slop: int = 0,
    strand_aware_slop: bool = False,
    ignore_strand: bool = False,
    require_different_names: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report overlaps between two paired-end BED (BEDPE) files (``bedtools pairtopair``).

    ``overlap_type`` is ``"both"`` (default), ``"either"``, ``"neither"`` or ``"notboth"`` (``-type``); ``min_overlap``
    is ``-f`` (fraction of A, BEDTools default 1E-9); ``slop`` is added to each footprint of A (``-slop``) on both sides,
    or only on the strand-dependent side with ``strand_aware_slop`` (``-ss``); ``ignore_strand`` is ``-is`` and
    ``require_different_names`` is ``-rdn``.
    """
    a = _cli.existing_file(a_path, "a_path")
    b = _cli.existing_file(b_path, "b_path")
    args = ["-a", str(a), "-b", str(b), "-type", _cli.choice(overlap_type, _PAIR_TYPES, "overlap_type")]
    args += _cli.option("-f", None if min_overlap is None else _fraction(min_overlap, "min_overlap"))
    args += ["-slop", str(_cli.integer(slop, "slop", minimum=0))]
    args += _flag(strand_aware_slop, "-ss") + _flag(ignore_strand, "-is") + _flag(require_different_names, "-rdn")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("pairtopair", args, out, bedtools_path, timeout_seconds)


def bedtools_multicov(
    bam_paths: Sequence[str],
    intervals_path: str,
    output_path: str,
    *,
    split: bool = False,
    strand: str | None = None,
    min_overlap: float | None = None,
    reciprocal: bool = False,
    min_mapping_quality: int | None = None,
    include_duplicates: bool = False,
    include_failed_qc: bool = False,
    proper_pairs_only: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the alignments of several BAM files that overlap each interval (``bedtools multicov``).

    Each output line is the interval followed by one count per BAM file, in the order given. The BAM files must be
    indexed. ``strand`` is ``None``, ``"same"`` (``-s``) or ``"opposite"`` (``-S``); ``split`` treats split alignments
    as separate intervals; ``min_overlap`` (``-f``) and ``reciprocal`` (``-r``) set the overlap fraction;
    ``min_mapping_quality`` is ``-q``; ``include_duplicates`` is ``-D``, ``include_failed_qc`` is ``-F`` and
    ``proper_pairs_only`` is ``-p``.
    """
    bams = _existing_list(bam_paths, "bam_paths", 1)
    intervals = _cli.existing_file(intervals_path, "intervals_path")
    args = ["-bams", *[str(path) for path in bams], "-bed", str(intervals)]
    args += _flag(split, "-split") + _strand(strand, "-s", "-S")
    args += _cli.option("-f", None if min_overlap is None else _fraction(min_overlap, "min_overlap"))
    args += _flag(reciprocal, "-r")
    args += _cli.option(
        "-q", None if min_mapping_quality is None else _cli.integer(min_mapping_quality, "min_mapping_quality", minimum=0)
    )
    args += _flag(include_duplicates, "-D") + _flag(include_failed_qc, "-F") + _flag(proper_pairs_only, "-p")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("multicov", args, out, bedtools_path, timeout_seconds)


# ---------------------------------------------------------------------------
# reldist, coverage, bamtobed -bedpe
# ---------------------------------------------------------------------------


def _reldist_rows(path: Path) -> list[dict[str, Any]]:
    """Parse the summary table (``reldist``, ``count``, ``total``, ``fraction``) of ``bedtools reldist``."""
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) == 4 and fields[0].replace(".", "", 1).isdigit():
            rows.append(
                {
                    "reldist": float(fields[0]),
                    "count": int(fields[1]),
                    "total": int(fields[2]),
                    "fraction": float(fields[3]),
                }
            )
    return rows


def bedtools_reldist(
    a_path: str,
    b_path: str,
    output_path: str,
    *,
    detail: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distribution of the relative distances between two interval sets (``bedtools reldist``).

    By default the summary table (relative distance 0.00-0.50, interval count, total, fraction) is written and also
    returned in ``summary``. With ``detail=True`` (``-detail``) the output has the relative distance of each A interval
    to the closest B intervals instead, and ``summary`` is empty.
    """
    a = _cli.existing_file(a_path, "a_path")
    b = _cli.existing_file(b_path, "b_path")
    args = ["-a", str(a), "-b", str(b)] + _flag(detail, "-detail")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    result = _stream("reldist", args, out, bedtools_path, timeout_seconds)
    result["summary"] = [] if detail or result["status"] != "completed" else _reldist_rows(out)
    return result


def bedtools_coverage_profile(
    a_path: str,
    b_path: str,
    output_path: str,
    *,
    report: str = "histogram",
    strand: str | None = None,
    min_overlap_a: float | None = None,
    min_overlap_b: float | None = None,
    reciprocal: bool = False,
    either_fraction: bool = False,
    split: bool = False,
    sorted_input: bool = False,
    genome_path: str | None = None,
    header: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage of A intervals by B features as a histogram or as per-position depth (``bedtools coverage``).

    ``report="histogram"`` (``-hist``) writes, after each A interval, the depth, the number of bases at that depth, the
    size of A and the fraction of A at that depth, plus a final summary histogram of all intervals. ``report="depth"``
    (``-d``) writes the depth at every one-based position inside each A interval. ``strand`` is ``None``, ``"same"``
    (``-s``) or ``"opposite"`` (``-S``); ``min_overlap_a``/``min_overlap_b`` are ``-f``/``-F``; ``reciprocal`` is
    ``-r``; ``either_fraction`` is ``-e``; ``split`` treats split BED12 or BAM entries as separate intervals;
    ``sorted_input`` is ``-sorted`` (with ``genome_path`` as ``-g``).
    """
    a = _cli.existing_file(a_path, "a_path")
    b = _cli.existing_file(b_path, "b_path")
    mode = _cli.choice(report, ("histogram", "depth"), "report")
    args = ["-a", str(a), "-b", str(b), "-hist" if mode == "histogram" else "-d"]
    args += _strand(strand, "-s", "-S")
    args += _cli.option("-f", None if min_overlap_a is None else _fraction(min_overlap_a, "min_overlap_a"))
    args += _cli.option("-F", None if min_overlap_b is None else _fraction(min_overlap_b, "min_overlap_b"))
    args += _flag(reciprocal, "-r") + _flag(either_fraction, "-e") + _flag(split, "-split")
    args += _flag(sorted_input, "-sorted")
    if genome_path is not None:
        if not sorted_input:
            raise ConfigurationError("genome_path only applies together with sorted_input=True.")
        args += ["-g", str(_cli.existing_file(genome_path, "genome_path"))]
    args += _flag(header, "-header")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    result = _stream("coverage", args, out, bedtools_path, timeout_seconds)
    result["report"] = mode
    return result


def bedtools_bam_to_bedpe(
    bam_path: str,
    output_path: str,
    *,
    mate1_first: bool = False,
    edit_distance_score: bool = False,
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert paired-end BAM alignments to BEDPE (``bedtools bamtobed -bedpe``).

    The BAM must be grouped or sorted by query name. ``mate1_first`` (``-mate1``) always reports mate one as the first
    block; ``edit_distance_score`` (``-ed``) uses the summed NM edit distance of the two mates as the score instead of
    the smaller mapping quality.
    """
    bam = _cli.existing_file(bam_path, "bam_path")
    args = ["-i", str(bam), "-bedpe"] + _flag(mate1_first, "-mate1") + _flag(edit_distance_score, "-ed")
    _check_program(bedtools_path)
    out = _cli.new_file(output_path, "output_path")
    return _stream("bamtobed", args, out, bedtools_path, timeout_seconds)


# ---------------------------------------------------------------------------
# split
# ---------------------------------------------------------------------------


def bedtools_split(
    input_path: str,
    output_directory: str,
    number_of_files: int,
    *,
    prefix: str = "split",
    algorithm: str = "size",
    bedtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a BED file into ``number_of_files`` files inside a fresh directory (``bedtools split``).

    ``algorithm="size"`` (BEDTools' default) groups records so that every file holds about the same number of bases;
    ``algorithm="simple"`` (``-a simple``) routes records round-robin so that every file holds about the same number of
    records. The files are named ``<prefix>.00001.bed`` and so on and listed in ``outputs`` with their record counts in
    ``files``.
    """
    source = _cli.existing_file(input_path, "input_path")
    number = _cli.integer(number_of_files, "number_of_files", minimum=1)
    mode = _cli.choice(algorithm, ("size", "simple"), "algorithm")
    stem = _cli.text(prefix, "prefix")
    if "/" in stem or "\\" in stem:
        raise ConfigurationError("prefix must be a simple file stem.")
    _check_program(bedtools_path)
    directory = _cli.new_directory(output_directory, "output_directory")
    args = ["split", "-i", str(source), "-n", str(number), "-p", stem, "-a", mode]
    result = _p13.run(
        _PROGRAM, args, cwd=directory, expected=[], explicit=bedtools_path, hint=_HINT, timeout_seconds=timeout_seconds
    )
    files = sorted(directory.glob(f"{stem}.*.bed"))
    result["operation"] = "bedtools split"
    result["outputs"] = [str(path) for path in files]
    result["files"] = [{"file": path.name, "record_count": _p13.count_lines(path)} for path in files]
    if result["return_code"] == 0 and not files:
        result["status"] = "failed"
    return result


__all__ = [
    "bedtools_bam_to_bedpe",
    "bedtools_coverage_profile",
    "bedtools_multicov",
    "bedtools_multiinter",
    "bedtools_pairtopair",
    "bedtools_random",
    "bedtools_reldist",
    "bedtools_sample",
    "bedtools_shuffle",
    "bedtools_split",
    "bedtools_unionbedg",
    "bedtools_window",
]

_register_public_verbose(globals(), tuple(__all__))
