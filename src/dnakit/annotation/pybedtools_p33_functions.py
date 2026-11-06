"""pybedtools operations that DNAKit did not have yet (integration pack P33, part 1).

Conversions (BED12 to BED6, BAM to BED, BED or BEDPE to BAM, an interval table as records), file
format detection, interval transforms (slop, shift, windows, concatenation, truncation, removal of
invalid records, introns), BEDTools statistics (map, annotate, genome coverage, total coverage,
spacing, Fisher test, nucleotide content, pair-to-BED), randomisation statistics (overlap count or
base-pair distribution, significance summary, pairwise matrix), BAM tagging and report helpers
(feature track plot, IGV batch script, UCSC link page).

pybedtools drives the BEDTools programs, so BEDTools must be installed (tested with BEDTools 2.27.1
and pybedtools 0.12.1). Coordinates are 0-based half-open as in BED. Functions take file paths,
return a JSON-serialisable dict with the first ``max_rows`` records and, when ``output_path`` is
given, write the full result to that new file. Inputs of ``map``, ``spacing`` and ``fisher`` are
sorted first because BEDTools needs sorted input there. The randomisation functions take a seed
(shuffle ``i`` uses ``seed + i``, which equals pybedtools ``debug=True`` for ``seed=0``). pybedtools is
imported only when a function runs; a missing package raises ``BackendUnavailableError``.
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _p33_common as _c
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-pybedtools"
_MAP_OPERATIONS = (
    "sum",
    "min",
    "max",
    "absmin",
    "absmax",
    "mean",
    "median",
    "mode",
    "antimode",
    "stdev",
    "sstdev",
    "collapse",
    "distinct",
    "distinct_sort_num",
    "distinct_sort_num_desc",
    "count",
    "count_distinct",
    "first",
    "last",
)
_COVERAGE_REPORTS = ("histogram", "bedgraph", "bedgraph_all", "per_base", "per_base_nonzero")
_PAIR_TYPES = (
    "either",
    "neither",
    "xor",
    "both",
    "notboth",
    "ispan",
    "ospan",
    "notispan",
    "notospan",
)
_MAX_ITERATIONS = 10_000
_MAX_FILES = 50
_MAX_TRACK_INTERVALS = 20_000


def _pybedtools() -> Any:
    module = require("pybedtools", _EXTRA)
    if shutil.which("bedtools") is None:
        raise BackendUnavailableError(
            "bedtools was not found on PATH.",
            hint="Install BEDTools (conda install -c bioconda bedtools, or apt install bedtools).",
        )
    return module


def _provenance(pybedtools: Any) -> dict[str, Any]:
    try:
        bedtools = ".".join(str(part) for part in pybedtools.settings.bedtools_version)
    except Exception:  # pragma: no cover - attribute differs between pybedtools builds
        bedtools = None
    return {
        "package": "pybedtools",
        "package_version": _c.package_version("pybedtools"),
        "bedtools_version": bedtools,
    }


def _guard(operation: str, call: Any) -> Any:
    try:
        return call()
    except (ConfigurationError, BackendUnavailableError):
        raise
    except Exception as exc:  # pybedtools raises BEDToolsError with the BEDTools message
        raise ConfigurationError(f"pybedtools {operation} failed: {exc}") from exc


def _discard(*bedtools: Any) -> None:
    """Remove pybedtools temporary files that were created for this call."""

    for item in bedtools:
        name = getattr(item, "fn", None)
        if isinstance(name, str) and "pybedtools." in Path(name).name:
            Path(name).unlink(missing_ok=True)


def _text_result(
    operation: str,
    pybedtools: Any,
    out_text: str,
    output_path: Path | None,
    limit: int,
    **extra: Any,
) -> dict[str, Any]:
    if out_text and not out_text.endswith("\n"):
        out_text += "\n"
    lines = [line for line in out_text.splitlines() if line.strip()]
    if output_path is not None:
        output_path.write_text(out_text, encoding="utf-8")
    return {
        "status": "completed",
        "operation": operation,
        "n_records": len(lines),
        "records": [line.split("\t") for line in lines[:limit]],
        "truncated": len(lines) > limit,
        "output_path": str(output_path) if output_path else None,
        **extra,
        "provenance": _provenance(pybedtools),
    }


def _bam_result(
    operation: str,
    pybedtools: Any,
    produced: str,
    output_path: Path,
    **extra: Any,
) -> dict[str, Any]:
    shutil.copyfile(produced, output_path)
    count = None
    try:
        import pysam  # optional: only used to count the alignments written

        with pysam.AlignmentFile(str(output_path), "rb", check_sq=False) as handle:
            count = handle.count(until_eof=True)
    except Exception:  # pragma: no cover - the count is informative only
        count = None
    return {
        "status": "completed",
        "operation": operation,
        "output_path": str(output_path),
        "n_alignments": count,
        **extra,
        "provenance": _provenance(pybedtools),
    }


def _genome(value: object) -> Path:
    path = input_path(value)  # type: ignore[arg-type]
    _c.check_genome_file(path)
    return path


def _sorted(pybedtools: Any, path: Path) -> Any:
    return pybedtools.BedTool(str(path)).sort()


def convert_bed12_to_bed6_pybedtools(
    bed12_path: str | Path,
    *,
    number_names: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split every BED12 record into one BED6 record per block (BedTool.bed6, bedtools bed12tobed6).

    ``number_names=True`` writes the 1-based block number into the score column (``bedtools bed12tobed6 -n``).
    """

    path = input_path(bed12_path)
    number = _c.flag(number_names, "number_names")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    kwargs = {"n": True} if number else {}
    result = _guard("bed6", lambda: pb.BedTool(str(path)).bed6(**kwargs))
    try:
        return _text_result("convert_bed12_to_bed6_pybedtools", pb, str(result), out, limit)
    finally:
        _discard(result)


def bam_to_bed_pybedtools(
    bam_path: str | Path,
    *,
    bed12: bool = False,
    split: bool = False,
    split_deletions: bool = False,
    edit_distance: bool = False,
    tag: str | None = None,
    color: str | None = None,
    cigar: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert BAM alignments to BED, BED12 or split BED (BedTool.bam_to_bed, bedtools bamtobed).

    ``bed12`` writes BED12, ``split`` splits on N CIGAR operations, ``split_deletions`` also splits
    on D, ``edit_distance`` puts the NM tag in the score, ``tag`` uses another tag as score,
    ``color`` is an ``R,G,B`` string for the itemRgb column, ``cigar`` appends the CIGAR string.
    """

    path = input_path(bam_path)
    kwargs: dict[str, Any] = {}
    for name, value, key in (
        ("bed12", bed12, "bed12"),
        ("split", split, "split"),
        ("split_deletions", split_deletions, "splitD"),
        ("edit_distance", edit_distance, "ed"),
        ("cigar", cigar, "cigar"),
    ):
        if _c.flag(value, name):
            kwargs[key] = True
    if tag is not None:
        kwargs["tag"] = _c.label(tag, "tag")
    if color is not None:
        kwargs["color"] = _c.text(color, "color", 32)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("bam_to_bed", lambda: pb.BedTool(str(path)).bam_to_bed(**kwargs))
    try:
        return _text_result("bam_to_bed_pybedtools", pb, str(result), out, limit)
    finally:
        _discard(result)


def convert_intervals_to_bam_pybedtools(
    intervals_path: str | Path,
    genome_path: str | Path,
    output_path: str | Path,
    *,
    input_format: str = "bed",
    mapq: int = 255,
    bed12: bool = False,
    uncompressed: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write BED (``input_format="bed"``, BedTool.to_bam, bedtools bedtobam) or BEDPE
    (``input_format="bedpe"``, BedTool.bedpe_to_bam, bedpeToBam) intervals as an unaligned BAM.

    ``genome_path`` is a chromosome-size file, ``mapq`` the mapping quality written, ``bed12``
    treats BED12 blocks as CIGAR blocks (BED only) and ``uncompressed`` writes an uncompressed BAM.
    """

    path = input_path(intervals_path)
    genome = _genome(genome_path)
    fmt = _c.choice(input_format, ("bed", "bedpe"), "input_format")
    quality = _c.integer(mapq, "mapq", 0, 255)
    block = _c.flag(bed12, "bed12")
    raw = _c.flag(uncompressed, "uncompressed")
    if block and fmt == "bedpe":
        raise ConfigurationError("bed12 applies to BED input only.")
    out = _c.new_output(output_path)
    assert out is not None
    pb = _pybedtools()
    kwargs: dict[str, Any] = {"g": str(genome), "mapq": quality}
    if block:
        kwargs["bed12"] = True
    if raw:
        kwargs["ubam"] = True
    bt = pb.BedTool(str(path))
    method = bt.to_bam if fmt == "bed" else bt.bedpe_to_bam
    result = _guard("to_bam", lambda: method(**kwargs))
    try:
        return _bam_result(
            "convert_intervals_to_bam_pybedtools", pb, result.fn, out, input_format=fmt
        )
    finally:
        _discard(result)


def intervals_to_records_pybedtools(
    path: str | Path,
    *,
    disable_auto_names: bool = False,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse BED, GFF, VCF or SAM records into named columns (BedTool.to_dataframe).

    Column names follow pybedtools (``chrom start end name score strand thickStart thickEnd
    itemRgb blockCount blockSizes blockStarts`` for BED, the GFF and VCF names for those formats)
    unless ``disable_auto_names`` is true, which numbers the columns instead. The first
    ``max_rows`` records are returned as dicts.
    """

    file = input_path(path)
    off = _c.flag(disable_auto_names, "disable_auto_names")
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    frame = _guard(
        "to_dataframe",
        lambda: pb.BedTool(str(file)).to_dataframe(
            disable_auto_names=off, **({"header": None} if off else {})
        ),
    )
    records = [
        {str(key): _c.plain(value) for key, value in row.items()}
        for row in frame.head(limit).to_dict("records")
    ]
    return {
        "status": "completed",
        "operation": "intervals_to_records_pybedtools",
        "columns": [str(name) for name in frame.columns],
        "dtypes": {str(name): str(kind) for name, kind in frame.dtypes.items()},
        "n_records": len(frame),
        "records": records,
        "truncated": len(frame) > limit,
        "provenance": _provenance(pb),
    }


def list_tabix_contigs_pybedtools(
    path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the contig names of a bgzip-compressed, tabix-indexed file (BedTool.tabix_contigs).

    The ``.tbi`` (or ``.csi``) index must sit next to the file. Needs pysam.
    """

    file = input_path(path)
    pysam = require("pysam", "bio-reads")
    try:
        handle = pysam.TabixFile(str(file))
    except (OSError, ValueError):
        try:
            handle = pysam.TabixFile(str(file), index=str(file) + ".csi")
        except (OSError, ValueError) as exc:
            raise ConfigurationError(
                f"{file.name} has no tabix index (.tbi or .csi) next to it."
            ) from exc
    with handle:
        contigs = list(handle.contigs)
    return {
        "status": "completed",
        "operation": "list_tabix_contigs_pybedtools",
        "contigs": contigs,
        "n_contigs": len(contigs),
        "provenance": {"package": "pysam", "package_version": _c.package_version("pysam")},
    }


def detect_alignment_file_format_pybedtools(
    path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect BAM, CRAM, BGZF and gzip files from their magic bytes (pybedtools.helpers isBAM,
    isCRAM, isBGZIP and isGZIP).

    ``is_bgzip`` is true for BAM files too, because a BAM file is BGZF-compressed;
    ``format`` names the most specific match.
    """

    file = input_path(path)
    pb = _pybedtools()
    helpers = pb.helpers
    is_bam = bool(_guard("isBAM", lambda: helpers.isBAM(str(file))))
    is_cram = bool(_guard("isCRAM", lambda: helpers.isCRAM(str(file))))
    is_bgzip = bool(_guard("isBGZIP", lambda: helpers.isBGZIP(str(file))))
    is_gzip = bool(_guard("isGZIP", lambda: helpers.isGZIP(str(file))))
    if is_bam:
        kind = "bam"
    elif is_cram:
        kind = "cram"
    elif is_bgzip:
        kind = "bgzf"
    elif is_gzip:
        kind = "gzip"
    else:
        kind = "other"
    return {
        "status": "completed",
        "operation": "detect_alignment_file_format_pybedtools",
        "is_bam": is_bam,
        "is_cram": is_cram,
        "is_bgzip": is_bgzip,
        "is_gzip": is_gzip,
        "format": kind,
        "provenance": _provenance(pb),
    }


def interval_has_hits_pybedtools(
    intervals_path: str | Path,
    chrom: str,
    start: int,
    end: int,
    *,
    strand: str = ".",
    same_strand: bool = False,
    min_overlap_fraction: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether one query interval overlaps any interval of a file (BedTool.any_hits).

    ``same_strand`` requires equal strands, ``min_overlap_fraction`` is the fraction of the query
    that must overlap (0 means any overlap). Coordinates are BED style (0-based, half-open).
    """

    file = input_path(intervals_path)
    name = _c.label(chrom, "chrom")
    begin = _c.integer(start, "start", 0, 2_000_000_000)
    stop = _c.integer(end, "end", 0, 2_000_000_000)
    if stop < begin:
        raise ConfigurationError("end must not be smaller than start.")
    symbol = _c.choice(strand, ("+", "-", "."), "strand")
    same = _c.flag(same_strand, "same_strand")
    fraction = _c.number(min_overlap_fraction, "min_overlap_fraction", 0.0, 1.0)
    pb = _pybedtools()
    query = pb.create_interval_from_list([name, str(begin), str(stop), ".", "0", symbol])
    hit = _guard(
        "any_hits",
        lambda: pb.BedTool(str(file)).any_hits(query, same_strand=same, overlap=fraction),
    )
    return {
        "status": "completed",
        "operation": "interval_has_hits_pybedtools",
        "has_hits": bool(hit),
        "query": {"chrom": name, "start": begin, "end": stop, "strand": symbol},
        "provenance": _provenance(pb),
    }


def slop_pybedtools_intervals(
    bed_path: str | Path,
    genome_path: str | Path,
    *,
    both: float | None = None,
    left: float | None = None,
    right: float | None = None,
    strand_aware: bool = False,
    as_fraction: bool = False,
    header: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extend intervals by flanking bases, clipped to the chromosome ends (BedTool.slop).

    ``both`` extends both sides; ``left`` and ``right`` (given together) extend the two sides
    separately (upstream/downstream when ``strand_aware``). ``as_fraction`` reads the values as fractions of the interval length.
    """

    path = input_path(bed_path)
    genome = _genome(genome_path)
    kwargs: dict[str, Any] = {"g": str(genome)}
    fraction = _c.flag(as_fraction, "as_fraction")
    top = 100.0 if fraction else 2_000_000_000.0
    for key, value in (("b", both), ("l", left), ("r", right)):
        if value is not None:
            kwargs[key] = _c.number(value, {"b": "both", "l": "left", "r": "right"}[key], 0.0, top)
    if len(kwargs) == 1:
        raise ConfigurationError("Give both, or left and right.")
    if "b" in kwargs and ("l" in kwargs or "r" in kwargs):
        raise ConfigurationError("both cannot be combined with left or right.")
    if "b" not in kwargs and not ("l" in kwargs and "r" in kwargs):
        raise ConfigurationError("left and right must be given together.")
    if not fraction:
        for key in ("b", "l", "r"):
            if key in kwargs:
                kwargs[key] = int(kwargs[key])
    if _c.flag(strand_aware, "strand_aware"):
        kwargs["s"] = True
    if fraction:
        kwargs["pct"] = True
    if _c.flag(header, "header"):
        kwargs["header"] = True
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("slop", lambda: pb.BedTool(str(path)).slop(**kwargs))
    try:
        return _text_result("slop_pybedtools_intervals", pb, str(result), out, limit)
    finally:
        _discard(result)


def shift_pybedtools_intervals(
    bed_path: str | Path,
    genome_path: str | Path,
    *,
    shift: float | None = None,
    plus_shift: float | None = None,
    minus_shift: float | None = None,
    as_fraction: bool = False,
    header: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Move intervals along the chromosome, clipped to chromosome ends (BedTool.shift).

    ``shift`` moves every interval (negative values move left); ``plus_shift`` and ``minus_shift``
    move plus-strand and minus-strand intervals separately and override ``shift`` for that strand.
    ``as_fraction`` reads the shifts as fractions of the interval length.
    """

    path = input_path(bed_path)
    genome = _genome(genome_path)
    fraction = _c.flag(as_fraction, "as_fraction")
    kwargs: dict[str, Any] = {"g": str(genome)}
    for key, value, name in (
        ("s", shift, "shift"),
        ("p", plus_shift, "plus_shift"),
        ("m", minus_shift, "minus_shift"),
    ):
        if value is not None:
            checked = _c.number(value, name, -2_000_000_000.0, 2_000_000_000.0)
            kwargs[key] = checked if fraction else int(checked)
    if len(kwargs) == 1:
        raise ConfigurationError("Give shift, plus_shift or minus_shift.")
    if fraction:
        kwargs["pct"] = True
    if _c.flag(header, "header"):
        kwargs["header"] = True
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("shift", lambda: pb.BedTool(str(path)).shift(**kwargs))
    try:
        return _text_result("shift_pybedtools_intervals", pb, str(result), out, limit)
    finally:
        _discard(result)


def map_pybedtools_values(
    a_path: str | Path,
    b_path: str | Path,
    *,
    column: int = 5,
    operation: str = "sum",
    null: str = ".",
    delimiter: str = ",",
    min_overlap_a: float | None = None,
    min_overlap_b: float | None = None,
    reciprocal: bool = False,
    same_strand: bool = False,
    opposite_strand: bool = False,
    split: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a column of B over the B intervals that overlap each A interval (BedTool.map).

    ``operation`` is one of sum, min, max, absmin, absmax, mean, median, mode, antimode, stdev,
    sstdev, collapse, distinct, distinct_sort_num, distinct_sort_num_desc, count, count_distinct,
    first, last. A and B are sorted first; ``null`` is written for A intervals without overlap.
    """

    a = input_path(a_path)
    b = input_path(b_path)
    kwargs: dict[str, Any] = {
        "c": _c.integer(column, "column", 1, 10_000),
        "o": _c.choice(operation, _MAP_OPERATIONS, "operation"),
        "null": _c.text(null, "null", 32),
        "delim": _c.text(delimiter, "delimiter", 8),
    }
    if min_overlap_a is not None:
        kwargs["f"] = _c.number(min_overlap_a, "min_overlap_a", 0.0, 1.0)
    if min_overlap_b is not None:
        kwargs["F"] = _c.number(min_overlap_b, "min_overlap_b", 0.0, 1.0)
    for name, value, key in (
        ("reciprocal", reciprocal, "r"),
        ("same_strand", same_strand, "s"),
        ("opposite_strand", opposite_strand, "S"),
        ("split", split, "split"),
    ):
        if _c.flag(value, name):
            kwargs[key] = True
    if "s" in kwargs and "S" in kwargs:
        raise ConfigurationError("same_strand and opposite_strand exclude each other.")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    first, second = _sorted(pb, a), _sorted(pb, b)
    result = _guard("map", lambda: first.map(second, **kwargs))
    try:
        return _text_result("map_pybedtools_values", pb, str(result), out, limit)
    finally:
        _discard(result, first, second)


def annotate_pybedtools_intervals(
    a_path: str | Path,
    feature_files: Sequence[str | Path],
    *,
    names: Sequence[str] | None = None,
    counts: bool = False,
    both: bool = False,
    same_strand: bool = False,
    opposite_strand: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report the coverage of each A interval by features of several files (BedTool.annotate).

    One column group per feature file: fraction covered (default), overlap counts (``counts``) or
    both (``both``). ``names`` labels the files in the header line.
    """

    a = input_path(a_path)
    if not isinstance(feature_files, Sequence) or isinstance(feature_files, (str, Path)):
        raise ConfigurationError("feature_files must be a list of paths.")
    if not 1 <= len(feature_files) <= _MAX_FILES:
        raise ConfigurationError(f"feature_files must hold 1 to {_MAX_FILES} paths.")
    files = [str(input_path(item)) for item in feature_files]
    kwargs: dict[str, Any] = {"files": files}
    if names is not None:
        if not isinstance(names, Sequence) or isinstance(names, str) or len(names) != len(files):
            raise ConfigurationError("names must be a list with one label per feature file.")
        kwargs["names"] = [_c.label(item, "names") for item in names]
    for name, value, key in (
        ("counts", counts, "counts"),
        ("both", both, "both"),
        ("same_strand", same_strand, "s"),
        ("opposite_strand", opposite_strand, "S"),
    ):
        if _c.flag(value, name):
            kwargs[key] = True
    if "s" in kwargs and "S" in kwargs:
        raise ConfigurationError("same_strand and opposite_strand exclude each other.")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("annotate", lambda: pb.BedTool(str(a)).annotate(**kwargs))
    try:
        return _text_result("annotate_pybedtools_intervals", pb, str(result), out, limit)
    finally:
        _discard(result)


def genome_coverage_pybedtools(
    path: str | Path,
    genome_path: str | Path | None = None,
    *,
    report: str = "histogram",
    scale: float = 1.0,
    strand: str | None = None,
    five_prime: bool = False,
    three_prime: bool = False,
    split: bool = False,
    max_depth: int | None = None,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage of a genome by BED intervals or BAM alignments (BedTool.genome_coverage).

    Interval input is sorted first; BAM input must be coordinate-sorted.
    ``report``: ``histogram`` (depth histogram per chromosome and genome-wide, the default),
    ``bedgraph`` (covered runs), ``bedgraph_all`` (also zero-coverage runs), ``per_base`` (depth at
    each position) or ``per_base_nonzero`` (positions with non-zero depth). ``genome_path`` is
    required for BED input; BAM input takes the sizes from its header when omitted.
    """

    file = input_path(path)
    mode = _c.choice(report, _COVERAGE_REPORTS, "report")
    kwargs: dict[str, Any] = {
        {
            "histogram": "",
            "bedgraph": "bg",
            "bedgraph_all": "bga",
            "per_base": "d",
            "per_base_nonzero": "dz",
        }[mode]: True
    }
    kwargs.pop("", None)
    is_bam = file.suffix.lower() == ".bam"
    if genome_path is not None:
        kwargs["g"] = str(_genome(genome_path))
    elif not is_bam:
        raise ConfigurationError("genome_path is required for interval (non-BAM) input.")
    factor = _c.number(scale, "scale", 0.0, 1e12)
    if factor != 1.0:
        if mode == "histogram":
            raise ConfigurationError("scale needs a bedgraph or per-base report.")
        kwargs["scale"] = factor
    if strand is not None:
        kwargs["strand"] = _c.choice(strand, ("+", "-"), "strand")
    if _c.flag(five_prime, "five_prime"):
        kwargs["5"] = True
    if _c.flag(three_prime, "three_prime"):
        kwargs["3"] = True
    if "5" in kwargs and "3" in kwargs:
        raise ConfigurationError("five_prime and three_prime exclude each other.")
    if _c.flag(split, "split"):
        kwargs["split"] = True
    if max_depth is not None:
        kwargs["max"] = _c.integer(max_depth, "max_depth", 1, 2_000_000_000)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    source = pb.BedTool(str(file)) if is_bam else _sorted(pb, file)
    result = _guard("genome_coverage", lambda: source.genome_coverage(**kwargs))
    try:
        return _text_result("genome_coverage_pybedtools", pb, str(result), out, limit, report=mode)
    finally:
        _discard(result, source if not is_bam else None)


def total_coverage_pybedtools(
    path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the bases covered by at least one interval after merging (BedTool.total_coverage)."""

    file = input_path(path)
    pb = _pybedtools()
    bt = pb.BedTool(str(file))
    total = _guard("total_coverage", lambda: bt.total_coverage())
    return {
        "status": "completed",
        "operation": "total_coverage_pybedtools",
        "total_bases": int(total),
        "n_intervals": int(bt.count()),
        "provenance": _provenance(pb),
    }


def pair_to_bed_pybedtools(
    bedpe_path: str | Path,
    bed_path: str | Path,
    *,
    pair_type: str = "either",
    min_overlap_fraction: float | None = None,
    same_strand: bool = False,
    opposite_strand: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the BEDPE pairs that overlap BED intervals (BedTool.pair_to_bed, bedtools pairtobed).

    ``pair_type``: either, neither, xor, both, notboth, ispan, ospan, notispan or notospan.
    """

    pair = input_path(bedpe_path)
    bed = input_path(bed_path)
    kwargs: dict[str, Any] = {"b": str(bed), "type": _c.choice(pair_type, _PAIR_TYPES, "pair_type")}
    if min_overlap_fraction is not None:
        kwargs["f"] = _c.number(min_overlap_fraction, "min_overlap_fraction", 0.0, 1.0)
    if _c.flag(same_strand, "same_strand"):
        kwargs["s"] = True
    if _c.flag(opposite_strand, "opposite_strand"):
        kwargs["S"] = True
    if "s" in kwargs and "S" in kwargs:
        raise ConfigurationError("same_strand and opposite_strand exclude each other.")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("pair_to_bed", lambda: pb.BedTool(str(pair)).pair_to_bed(**kwargs))
    try:
        return _text_result("pair_to_bed_pybedtools", pb, str(result), out, limit)
    finally:
        _discard(result)


def cat_pybedtools_intervals(
    paths: Sequence[str | Path],
    *,
    postmerge: bool = True,
    force_truncate: bool = False,
    merge_distance: int = 0,
    same_strand: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate interval files and, by default, merge the result (BedTool.cat).

    With ``postmerge`` the files are concatenated, sorted and merged (``merge_distance`` is the
    maximum gap, ``same_strand`` merges per strand); the result has the first three columns only
    unless ``force_truncate`` is false and the columns agree. Without ``postmerge`` the files are
    only concatenated.
    """

    if not isinstance(paths, Sequence) or isinstance(paths, (str, Path)):
        raise ConfigurationError("paths must be a list of interval files.")
    if not 2 <= len(paths) <= _MAX_FILES:
        raise ConfigurationError(f"paths must hold 2 to {_MAX_FILES} files.")
    files = [str(input_path(item)) for item in paths]
    kwargs: dict[str, Any] = {
        "postmerge": _c.flag(postmerge, "postmerge"),
        "force_truncate": _c.flag(force_truncate, "force_truncate"),
    }
    distance = _c.integer(merge_distance, "merge_distance", 0, 2_000_000_000)
    if distance:
        kwargs["d"] = distance
    if _c.flag(same_strand, "same_strand"):
        kwargs["s"] = True
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    first = pb.BedTool(files[0])
    result = _guard("cat", lambda: first.cat(*files[1:], **kwargs))
    try:
        return _text_result("cat_pybedtools_intervals", pb, str(result), out, limit)
    finally:
        _discard(result)


def make_windows_pybedtools(
    source_path: str | Path,
    *,
    source_type: str = "genome",
    window_size: int | None = None,
    step_size: int | None = None,
    number_of_windows: int | None = None,
    reverse_numbering: bool = False,
    id_naming: str | None = None,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut a genome or BED intervals into fixed-size or fixed-number windows (BedTool.window_maker).

    ``source_type`` is ``genome`` (chromosome-size file) or ``bed``. Give ``window_size`` (with an
    optional ``step_size`` for sliding windows) or ``number_of_windows``. ``id_naming`` is None,
    ``winnum`` or ``srcwinnum`` (names the windows by number or source plus number).
    """

    file = input_path(source_path)
    kind = _c.choice(source_type, ("genome", "bed"), "source_type")
    if (window_size is None) == (number_of_windows is None):
        raise ConfigurationError("Give exactly one of window_size and number_of_windows.")
    kwargs: dict[str, Any] = {}
    if kind == "genome":
        kwargs["g"] = str(_genome(file))
    else:
        kwargs["b"] = str(file)
    if window_size is not None:
        kwargs["w"] = _c.integer(window_size, "window_size", 1, 2_000_000_000)
    if step_size is not None:
        if window_size is None:
            raise ConfigurationError("step_size needs window_size.")
        kwargs["s"] = _c.integer(step_size, "step_size", 1, 2_000_000_000)
    if number_of_windows is not None:
        kwargs["n"] = _c.integer(number_of_windows, "number_of_windows", 1, 10_000_000)
    if _c.flag(reverse_numbering, "reverse_numbering"):
        kwargs["reverse"] = True
    if id_naming is not None:
        kwargs["i"] = _c.choice(id_naming, ("winnum", "srcwinnum"), "id_naming")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard(
        "window_maker", lambda: pb.BedTool("chr1\t0\t1\n", from_string=True).window_maker(**kwargs)
    )
    try:
        return _text_result("make_windows_pybedtools", pb, str(result), out, limit)
    finally:
        _discard(result)


def spacing_pybedtools_intervals(
    bed_path: str | Path,
    *,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report the gap to the previous interval for every interval (BedTool.spacing).

    The input is sorted first. The first interval of each chromosome gets ``.`` and overlapping
    intervals get -1, as BEDTools writes them.
    """

    path = input_path(bed_path)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    sorted_bt = _sorted(pb, path)
    result = _guard("spacing", lambda: sorted_bt.spacing())
    try:
        return _text_result("spacing_pybedtools_intervals", pb, str(result), out, limit)
    finally:
        _discard(result, sorted_bt)


def fisher_pybedtools_overlap(
    a_path: str | Path,
    b_path: str | Path,
    genome_path: str | Path,
    *,
    min_overlap_a: float | None = None,
    reciprocal: bool = False,
    same_strand: bool = False,
    opposite_strand: bool = False,
    split: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fisher's exact test of whether two interval sets overlap more than chance (BedTool.fisher).

    Returns the 2x2 contingency table of counts and the left, right and two-tail p-values and the
    odds ratio that bedtools fisher prints. A and B are sorted first.
    """

    a = input_path(a_path)
    b = input_path(b_path)
    genome = _genome(genome_path)
    kwargs: dict[str, Any] = {"g": str(genome)}
    if min_overlap_a is not None:
        kwargs["f"] = _c.number(min_overlap_a, "min_overlap_a", 0.0, 1.0)
    for name, value, key in (
        ("reciprocal", reciprocal, "r"),
        ("same_strand", same_strand, "s"),
        ("opposite_strand", opposite_strand, "S"),
        ("split", split, "split"),
    ):
        if _c.flag(value, name):
            kwargs[key] = True
    if "s" in kwargs and "S" in kwargs:
        raise ConfigurationError("same_strand and opposite_strand exclude each other.")
    pb = _pybedtools()
    first, second = _sorted(pb, a), _sorted(pb, b)
    try:
        result = _guard("fisher", lambda: first.fisher(second, **kwargs))
    finally:
        _discard(first, second)
    return {
        "status": "completed",
        "operation": "fisher_pybedtools_overlap",
        "table": {key: dict(value) for key, value in result.table.items()},
        "left_tail": result.left_tail,
        "right_tail": result.right_tail,
        "two_tail": result.two_tail,
        "ratio": result.ratio,
        "report": result.text,
        "provenance": _provenance(pb),
    }


def introns_pybedtools(
    annotation_path: str | Path,
    *,
    gene: str = "gene",
    exon: str = "exon",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Derive intron and UTR intervals from GFF gene and exon features, or from BED12 (BedTool.introns).

    ``gene`` and ``exon`` are the GFF feature types. The result is BED6 with names of the gene and
    the feature kind (``intron``, ``utr5``, ``utr3``). pybedtools pairs every gene with all
    same-strand exons inside it, so overlapping genes can yield records with start greater than
    end; those are dropped and counted in ``n_invalid_dropped``.
    """

    path = input_path(annotation_path)
    gene_type = _c.label(gene, "gene")
    exon_type = _c.label(exon, "exon")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard(
        "introns", lambda: pb.BedTool(str(path)).introns(gene=gene_type, exon=exon_type)
    )
    try:
        raw = Path(result.fn).read_text(encoding="utf-8").splitlines()
        kept = []
        for line in raw:
            fields = line.split("\t")
            if len(fields) >= 3 and int(fields[1]) <= int(fields[2]):
                kept.append(line)
        return _text_result(
            "introns_pybedtools",
            pb,
            "\n".join(kept),
            out,
            limit,
            n_invalid_dropped=len(raw) - len(kept),
        )
    finally:
        _discard(result)


def truncate_pybedtools_to_chrom(
    bed_path: str | Path,
    genome_path: str | Path,
    *,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip intervals to the chromosome bounds of a size file (BedTool.truncate_to_chrom).

    Intervals beyond a chromosome end are cut, intervals on unknown chromosomes are dropped.
    """

    path = input_path(bed_path)
    sizes = _c.check_genome_file(_genome(genome_path))
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    bounds = {name: (0, length) for name, length in sizes.items()}
    result = _guard("truncate_to_chrom", lambda: pb.BedTool(str(path)).truncate_to_chrom(bounds))
    try:
        return _text_result("truncate_pybedtools_to_chrom", pb, str(result), out, limit)
    finally:
        _discard(result)


def remove_invalid_pybedtools_intervals(
    bed_path: str | Path,
    *,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop records with invalid coordinates (BedTool.remove_invalid).

    Removes records whose start is greater than the end, negative coordinates and malformed lines.
    ``n_removed`` counts the dropped data lines.
    """

    path = input_path(bed_path)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    total = sum(
        1
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(("#", "track", "browser"))
    )
    result = _guard("remove_invalid", lambda: pb.BedTool(str(path)).remove_invalid().saveas())
    try:
        text = str(result)
        kept = sum(1 for line in text.splitlines() if line.strip())
        return _text_result(
            "remove_invalid_pybedtools_intervals", pb, text, out, limit, n_removed=total - kept
        )
    finally:
        _discard(result)


def tag_bam_pybedtools(
    bam_path: str | Path,
    interval_files: Sequence[str | Path],
    output_path: str | Path,
    *,
    labels: Sequence[str] | None = None,
    tag: str = "YB",
    same_strand: bool = False,
    opposite_strand: bool = False,
    min_overlap_fraction: float | None = None,
    report: str = "labels",
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a BAM tag naming the interval files (or interval names) that each alignment overlaps
    (BedTool.tag_bam, bedtools tag).

    ``labels`` has one label per interval file (default: file names); ``report`` selects
    ``labels``, ``names`` (BED name column), ``scores`` (BED score column) or ``intervals``
    (chrom:start-end). The tagged alignments are written to the new BAM ``output_path``.
    """

    bam = input_path(bam_path)
    if not isinstance(interval_files, Sequence) or isinstance(interval_files, (str, Path)):
        raise ConfigurationError("interval_files must be a list of paths.")
    if not 1 <= len(interval_files) <= _MAX_FILES:
        raise ConfigurationError(f"interval_files must hold 1 to {_MAX_FILES} paths.")
    files = [str(input_path(item)) for item in interval_files]
    if labels is None:
        names = [Path(item).stem for item in files]
    else:
        if not isinstance(labels, Sequence) or isinstance(labels, str) or len(labels) != len(files):
            raise ConfigurationError("labels must hold one label per interval file.")
        names = [_c.label(item, "labels") for item in labels]
    kwargs: dict[str, Any] = {
        "files": files,
        "labels": names,
        "tag": _c.text(tag, "tag", 2),
    }
    mode = _c.choice(report, ("labels", "names", "scores", "intervals"), "report")
    if mode != "labels":
        kwargs[{"names": "names", "scores": "scores", "intervals": "intervals"}[mode]] = True
    if _c.flag(same_strand, "same_strand"):
        kwargs["s"] = True
    if _c.flag(opposite_strand, "opposite_strand"):
        kwargs["S"] = True
    if "s" in kwargs and "S" in kwargs:
        raise ConfigurationError("same_strand and opposite_strand exclude each other.")
    if min_overlap_fraction is not None:
        kwargs["f"] = _c.number(min_overlap_fraction, "min_overlap_fraction", 0.0, 1.0)
    out = _c.new_output(output_path)
    assert out is not None
    pb = _pybedtools()
    result = _guard("tag_bam", lambda: pb.BedTool(str(bam)).tag_bam(**kwargs))
    produced = getattr(result, "filename", None) or getattr(result, "fn", None)
    if not isinstance(produced, str):
        raise ConfigurationError("bedtools tag produced no BAM file.")
    try:
        return _bam_result("tag_bam_pybedtools", pb, produced, out, tag=kwargs["tag"], labels=names)
    finally:
        Path(produced).unlink(missing_ok=True)


def _shuffled_counts(
    pb: Any,
    a: Path,
    b: Path,
    genome: Path,
    iterations: int,
    metric: str,
    seed: int,
) -> list[int]:
    first = pb.BedTool(str(a))
    second = pb.BedTool(str(b))
    values: list[int] = []
    for index in range(iterations):
        shuffled = first.shuffle(g=str(genome), seed=seed + index)
        try:
            if metric == "count":
                joined = shuffled.intersect(second, u=True)
                values.append(len(joined))
            else:
                joined = shuffled.intersect(second)
                values.append(sum(len(feature) for feature in joined))
        finally:
            _discard(shuffled)
            if "joined" in locals():
                _discard(joined)
    return values


def random_overlap_distribution_pybedtools(
    a_path: str | Path,
    b_path: str | Path,
    genome_path: str | Path,
    *,
    iterations: int = 100,
    metric: str = "count",
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distribution of overlaps between randomly shuffled A and fixed B intervals
    (BedTool.randomintersection for ``metric="count"``, randomintersection_bp for ``"bp"``).

    ``metric="count"`` counts the shuffled A intervals that overlap B; ``metric="bp"`` sums the
    overlapping base pairs. Iteration ``i`` shuffles A with ``bedtools shuffle -seed seed+i``.
    Also returns the overlap of the unshuffled A for comparison.
    """

    a = input_path(a_path)
    b = input_path(b_path)
    genome = _genome(genome_path)
    kind = _c.choice(metric, ("count", "bp"), "metric")
    n = _c.integer(iterations, "iterations", 1, _MAX_ITERATIONS)
    start = _c.integer(seed, "seed", 0, 2**31 - 1 - _MAX_ITERATIONS)
    pb = _pybedtools()
    values = _guard(
        "randomintersection", lambda: _shuffled_counts(pb, a, b, genome, n, kind, start)
    )
    observed = _guard("intersect", lambda: _observed(pb, a, b, kind))
    return {
        "status": "completed",
        "operation": "random_overlap_distribution_pybedtools",
        "metric": kind,
        "iterations": n,
        "seed": start,
        "actual": observed,
        "distribution": values,
        "provenance": _provenance(pb),
    }


def _observed(pb: Any, a: Path, b: Path, metric: str) -> int:
    first, second = pb.BedTool(str(a)), pb.BedTool(str(b))
    if metric == "count":
        joined = first.intersect(second, u=True)
        value = len(joined)
    else:
        joined = first.intersect(second)
        value = sum(len(feature) for feature in joined)
    _discard(joined)
    return int(value)


def _summary(actual: int, distribution: Sequence[int], n_a: int, n_b: int) -> dict[str, Any]:
    import numpy as np

    data = np.array(list(distribution))
    median = float(np.median(data))
    n = float(len(data))

    def percentile_of_score(values: Any, score: int) -> float:
        values = np.array(values)
        size = len(values)
        if not np.any(values == score):
            values = np.append(values, score)
            positions = np.array(list(range(len(values))))
        else:
            positions = np.array(list(range(len(values)))) + 1.0
        values = np.sort(values)
        index = tuple([values == score])
        return float(np.mean(positions[index]) / size * 100.0)

    lower, upper = np.percentile(data, [2.5, 97.5])
    return {
        "actual": int(actual),
        "self": n_a,
        "other": n_b,
        "frac_randomized_above_actual": float(sum(data > actual) / n),
        "frac_randomized_below_actual": float(sum(data < actual) / n),
        "median_randomized": median,
        "normalized": (actual / median) if median else None,
        "percentile": percentile_of_score(data, actual),
        "lower_2.5th": float(lower),
        "upper_97.5th": float(upper),
    }


def random_overlap_summary_pybedtools(
    a_path: str | Path,
    b_path: str | Path,
    genome_path: str | Path,
    *,
    iterations: int = 100,
    seed: int = 0,
    include_distribution: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Empirical significance of the overlap of A with B by shuffling A (BedTool.randomstats).

    Returns the actual overlap count, the median of the randomised counts, the fractions of
    randomised counts above and below the actual one, the actual count divided by the median
    (``normalized``), the percentile of the actual count and the 2.5th/97.5th percentiles.
    """

    a = input_path(a_path)
    b = input_path(b_path)
    genome = _genome(genome_path)
    n = _c.integer(iterations, "iterations", 1, _MAX_ITERATIONS)
    start = _c.integer(seed, "seed", 0, 2**31 - 1 - _MAX_ITERATIONS)
    keep = _c.flag(include_distribution, "include_distribution")
    pb = _pybedtools()
    values = _guard("randomstats", lambda: _shuffled_counts(pb, a, b, genome, n, "count", start))
    actual = _guard("intersect", lambda: _observed(pb, a, b, "count"))
    summary = _summary(actual, values, len(pb.BedTool(str(a))), len(pb.BedTool(str(b))))
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "random_overlap_summary_pybedtools",
        "iterations": n,
        "seed": start,
        **summary,
        "provenance": _provenance(pb),
    }
    if keep:
        result["distribution"] = values
    return result


def random_overlap_matrix_pybedtools(
    bed_paths: Sequence[str | Path],
    genome_path: str | Path,
    *,
    iterations: int = 100,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise randomisation overlap statistics between several interval files
    (pybedtools.contrib.intersection_matrix.IntersectionMatrix.create_matrix).

    For every ordered pair (A, B) of files, A is shuffled ``iterations`` times and the summary of
    ``random_overlap_summary_pybedtools`` is computed; ``matrix`` is the nested mapping
    ``name_a -> name_b -> summary``. File names (without extension) must be unique.
    """

    if not isinstance(bed_paths, Sequence) or isinstance(bed_paths, (str, Path)):
        raise ConfigurationError("bed_paths must be a list of interval files.")
    if not 2 <= len(bed_paths) <= 12:
        raise ConfigurationError("bed_paths must hold 2 to 12 files.")
    files = [input_path(item) for item in bed_paths]
    names = [item.stem for item in files]
    if len(set(names)) != len(names):
        raise ConfigurationError("The interval files need unique names without extension.")
    genome = _genome(genome_path)
    n = _c.integer(iterations, "iterations", 1, _MAX_ITERATIONS)
    start = _c.integer(seed, "seed", 0, 2**31 - 1 - _MAX_ITERATIONS)
    pb = _pybedtools()
    matrix: dict[str, dict[str, Any]] = {}
    for first, name_a in zip(files, names, strict=False):
        matrix[name_a] = {}
        for second, name_b in zip(files, names, strict=False):
            values = _guard(
                "randomstats",
                lambda a=first, b=second: _shuffled_counts(pb, a, b, genome, n, "count", start),
            )
            actual = _guard("intersect", lambda a=first, b=second: _observed(pb, a, b, "count"))
            matrix[name_a][name_b] = _summary(
                actual, values, len(pb.BedTool(str(first))), len(pb.BedTool(str(second)))
            )
    return {
        "status": "completed",
        "operation": "random_overlap_matrix_pybedtools",
        "names": names,
        "iterations": n,
        "seed": start,
        "matrix": matrix,
        "provenance": _provenance(pb),
    }


def nucleotide_content_pybedtools(
    bed_path: str | Path,
    fasta_path: str | Path,
    *,
    strand_aware: bool = False,
    pattern: str | None = None,
    case_insensitive: bool = False,
    include_sequence: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide composition of the intervals in a FASTA file (BedTool.nucleotide_content,
    bedtools nuc).

    Returns the AT and GC fractions, base counts (A, C, G, T, N, other) and sequence length per
    interval; ``pattern`` adds a count of an IUPAC-free text pattern, ``include_sequence`` adds the
    sequence, ``strand_aware`` reverse-complements minus-strand intervals first.
    """

    bed = input_path(bed_path)
    fasta = input_path(fasta_path)
    kwargs: dict[str, Any] = {"fi": str(fasta)}
    if _c.flag(strand_aware, "strand_aware"):
        kwargs["s"] = True
    if pattern is not None:
        value = _c.text(pattern, "pattern", 256)
        if not value.replace("[", "").replace("]", "").isalpha():
            raise ConfigurationError("pattern must contain letters only.")
        kwargs["pattern"] = value
    if _c.flag(case_insensitive, "case_insensitive"):
        kwargs["C"] = True
    if _c.flag(include_sequence, "include_sequence"):
        kwargs["seq"] = True
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pb = _pybedtools()
    result = _guard("nucleotide_content", lambda: pb.BedTool(str(bed)).nucleotide_content(**kwargs))
    try:
        text = Path(result.fn).read_text(encoding="utf-8")
        lines = [line for line in text.splitlines() if line.strip()]
        header = [item for item in lines[0].lstrip("#").split("\t")] if lines else []
        body = "\n".join(line for line in lines if not line.startswith("#"))
        outcome = _text_result("nucleotide_content_pybedtools", pb, body, out, limit)
        outcome["header"] = header
        return outcome
    finally:
        _discard(result)


def _read_track(path: Path, chrom: str, start: int, end: int) -> list[tuple[int, int, str]]:
    found: list[tuple[int, int, str]] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] != chrom:
                continue
            try:
                left, right = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(f"{path.name}: start and end must be integers.") from exc
            if right > start and left < end:
                found.append((left, right, fields[3] if len(fields) > 3 else ""))
            if len(found) > _MAX_TRACK_INTERVALS:
                raise ConfigurationError(f"{path.name} has too many intervals in the region.")
    return found


def plot_interval_tracks_pybedtools(
    tracks: Mapping[str, str | Path],
    chrom: str,
    start: int,
    end: int,
    output_path: str | Path,
    *,
    label_intervals: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw several BED files as stacked feature tracks over one region and save a PNG
    (pybedtools.contrib.plotting.TrackCollection.plot).

    ``tracks`` maps a track label to a BED file; the track order is the mapping order. Returns the
    image path and the intervals drawn per track (clipped by nothing: full BED coordinates).
    """

    if not isinstance(tracks, Mapping) or not 1 <= len(tracks) <= 20:
        raise ConfigurationError("tracks must map 1 to 20 labels to BED files.")
    name = _c.label(chrom, "chrom")
    left = _c.integer(start, "start", 0, 2_000_000_000)
    right = _c.integer(end, "end", 1, 2_000_000_000)
    if right <= left:
        raise ConfigurationError("end must be greater than start.")
    show_names = _c.flag(label_intervals, "label_intervals")
    out = _c.new_output(output_path)
    assert out is not None
    drawn: dict[str, list[tuple[int, int, str]]] = {}
    for label, file in tracks.items():
        drawn[_c.label(label, "track label")] = _read_track(input_path(file), name, left, right)
    plt = require("matplotlib.figure", "bio-plots")
    backend = require("matplotlib.backends.backend_agg", "bio-plots")
    patches = require("matplotlib.patches", "bio-plots")
    figure = plt.Figure(figsize=(8, max(1.5, 0.8 * len(drawn) + 0.8)))
    backend.FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    for row, (_label, items) in enumerate(drawn.items()):
        for index, (lo, hi, feature) in enumerate(items):
            axes.add_patch(
                patches.Rectangle(
                    (lo, row + 0.1 + 0.2 * (index % 3)), max(hi - lo, 1), 0.15, color="#4477aa"
                )
            )
            if show_names and feature:
                axes.text(lo, row + 0.1 + 0.2 * (index % 3) + 0.17, feature, fontsize=6)
    axes.set_xlim(left, right)
    axes.set_ylim(0, len(drawn))
    axes.set_yticks([row + 0.5 for row in range(len(drawn))])
    axes.set_yticklabels(list(drawn))
    axes.set_xlabel(f"{name}:{left}-{right}")
    figure.tight_layout()
    figure.savefig(str(out), format="png", dpi=150)
    return {
        "status": "completed",
        "operation": "plot_interval_tracks_pybedtools",
        "image_path": str(out),
        "region": {"chrom": name, "start": left, "end": right},
        "tracks": {
            label: [{"start": lo, "end": hi, "name": feature} for lo, hi, feature in items]
            for label, items in drawn.items()
        },
        "provenance": {
            "package": "matplotlib",
            "package_version": _c.package_version("matplotlib"),
        },
    }


def build_igv_script_pybedtools(
    bed_path: str | Path,
    *,
    snapshot_dir: str = "/tmp/snapshots",
    image_format: str = "png",
    slop: int = 0,
    sort: str | None = None,
    collapse: bool = False,
    name_in_file: bool = False,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create an IGV batch script that takes a snapshot of every interval (BedTool.igv,
    bedtools igv).

    ``slop`` pads each locus, ``sort`` is one of position, strand, base, quality, sample, readGroup
    or None, ``collapse`` collapses the track and ``name_in_file`` names snapshots by the BED name.
    The script text is returned; nothing is run in IGV.
    """

    path = input_path(bed_path)
    kwargs: dict[str, Any] = {
        "path": _c.text(snapshot_dir, "snapshot_dir", 1024),
        "img": _c.choice(image_format, ("png", "eps", "svg"), "image_format"),
    }
    pad = _c.integer(slop, "slop", 0, 2_000_000_000)
    if pad:
        kwargs["slop"] = pad
    if sort is not None:
        kwargs["sort"] = _c.choice(
            sort, ("position", "strand", "base", "quality", "sample", "readGroup"), "sort"
        )
    if _c.flag(collapse, "collapse"):
        kwargs["clps"] = True
    if _c.flag(name_in_file, "name_in_file"):
        kwargs["name"] = True
    out = _c.new_output(output_path)
    pb = _pybedtools()
    result = _guard("igv", lambda: pb.BedTool(str(path)).igv(**kwargs))
    try:
        script = (
            Path(result.igv_script).read_text(encoding="utf-8")
            if hasattr(result, "igv_script")
            else str(result)
        )
    finally:
        _discard(result)
    if out is not None:
        out.write_text(script, encoding="utf-8")
    return {
        "status": "completed",
        "operation": "build_igv_script_pybedtools",
        "script": script,
        "n_commands": sum(1 for line in script.splitlines() if line.startswith("snapshot ")),
        "output_path": str(out) if out else None,
        "provenance": _provenance(pb),
    }


def build_ucsc_links_page_pybedtools(
    bed_path: str | Path,
    *,
    base_url: str = "http://genome.ucsc.edu",
    organism: str = "human",
    database: str = "hg19",
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write an HTML table that links each interval to its UCSC Genome Browser location
    (BedTool.links, bedtools links).

    The page is built offline; the links point at ``base_url``/cgi-bin/hgTracks for ``organism``
    and ``database``.
    """

    path = input_path(bed_path)
    url = _c.text(base_url, "base_url", 512)
    if not url.startswith(("http://", "https://")):
        raise ConfigurationError("base_url must start with http:// or https://.")
    org = _c.label(organism, "organism")
    db = _c.label(database, "database")
    out = _c.new_output(output_path)
    pb = _pybedtools()
    program = shutil.which("bedtools")
    assert program is not None
    try:
        done = subprocess.run(
            [program, "links", "-i", str(path), "-base", url, "-org", org, "-db", db],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ConfigurationError("bedtools links timed out.") from exc
    if done.returncode != 0:
        raise ConfigurationError(f"bedtools links failed: {done.stderr.strip()[:500]}")
    if out is not None:
        out.write_text(done.stdout, encoding="utf-8")
    return {
        "status": "completed",
        "operation": "build_ucsc_links_page_pybedtools",
        "html": done.stdout,
        "n_links": done.stdout.count("<a href="),
        "output_path": str(out) if out else None,
        "provenance": _provenance(pb),
    }


__all__ = [
    "annotate_pybedtools_intervals",
    "bam_to_bed_pybedtools",
    "build_igv_script_pybedtools",
    "build_ucsc_links_page_pybedtools",
    "cat_pybedtools_intervals",
    "convert_bed12_to_bed6_pybedtools",
    "convert_intervals_to_bam_pybedtools",
    "detect_alignment_file_format_pybedtools",
    "fisher_pybedtools_overlap",
    "genome_coverage_pybedtools",
    "interval_has_hits_pybedtools",
    "intervals_to_records_pybedtools",
    "introns_pybedtools",
    "list_tabix_contigs_pybedtools",
    "make_windows_pybedtools",
    "map_pybedtools_values",
    "nucleotide_content_pybedtools",
    "pair_to_bed_pybedtools",
    "plot_interval_tracks_pybedtools",
    "random_overlap_distribution_pybedtools",
    "random_overlap_matrix_pybedtools",
    "random_overlap_summary_pybedtools",
    "remove_invalid_pybedtools_intervals",
    "shift_pybedtools_intervals",
    "slop_pybedtools_intervals",
    "spacing_pybedtools_intervals",
    "tag_bam_pybedtools",
    "total_coverage_pybedtools",
    "truncate_pybedtools_to_chrom",
]

_register_public_verbose(globals(), tuple(__all__))
