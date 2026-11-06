"""pyranges1 operations that DNAKit did not have yet (integration pack P33, part 2).

Readers for WIG, ENCODE broadPeak and gappedPeak, Parquet and UCSC chain files; gene-annotation
helpers (tag test, one transcript per gene, UTR labelling); interval arithmetic (outer ranges,
resize, 5' and 3' ends, lengths, upstream and downstream flanks, windows, slicing, maximal
disjoint sets, cumulative length, strand repair, duplicate positions, summary); overlap joins and
the per-pair interval metrics and column combinations that build on a join; coordinate mapping
between local and global frames; liftover with a chain file; writers (Parquet, narrowPeak,
broadPeak, gappedPeak, Hi-C pairs, bedGraph, bigWig) and bigWig statistics over intervals; run-length
coverage.

These functions call pyranges1 methods that exist from pyranges1 1.5.0 on (tested with 1.5.0); with
an older release they raise ``BackendUnavailableError`` with the upgrade hint. ``bigwig_stats`` and
``coverage_runs`` also need ``pyBigWig``, ``pyarrow`` and ``runlength``. Coordinates are 0-based
half-open as in BED. Functions take file paths and return JSON-serialisable dicts with the first
``max_rows`` rows; ``output_path`` (a new file) receives the full table as TSV.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _p33_common as _c
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-pyranges1"
_INPUT_FORMATS = (
    "bed",
    "bedgraph",
    "gtf",
    "gff3",
    "narrowPeak",
    "broadPeak",
    "gappedPeak",
    "parquet",
    "pairs",
)
_STRAND = ("auto", "yes", "no")
_STRAND_BEHAVIOR = ("auto", "same", "opposite", "ignore")
_WRITE_FORMATS = ("parquet", "narrowPeak", "broadPeak", "gappedPeak", "pairs", "bedgraph")


def _pyranges() -> Any:
    module = require("pyranges1", _EXTRA)
    if not hasattr(module, "read_parquet") or not hasattr(module.PyRanges, "liftover_ranges"):
        raise BackendUnavailableError(
            "The installed pyranges1 is older than 1.5.0.",
            hint='Upgrade with: python -m pip install -U "pyranges1>=1.5"',
        )
    return module


def _provenance() -> dict[str, Any]:
    return {"package": "pyranges1", "package_version": _c.package_version("pyranges1")}


def _guard(operation: str, call: Any) -> Any:
    try:
        return call()
    except (ConfigurationError, BackendUnavailableError):
        raise
    except ImportError as exc:
        raise BackendUnavailableError(
            f"pyranges1 {operation} needs an extra package: {exc}",
            hint='Install with: python -m pip install "dnakit[bio-pyranges1]"',
        ) from exc
    except Exception as exc:
        raise ConfigurationError(f"pyranges1 {operation} failed: {exc}") from exc


def _strand(value: object) -> Any:
    mode = _c.choice(value, _STRAND, "use_strand")
    return {"auto": "auto", "yes": True, "no": False}[mode]


def _group(value: object, name: str = "group_by") -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return [_c.label(value, name)]
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        if not 1 <= len(value) <= 20:
            raise ConfigurationError(f"{name} must hold 1 to 20 column names.")
        return [_c.label(item, name) for item in value]
    raise ConfigurationError(f"{name} must be a column name or a list of column names.")


def _load(pr: Any, path: Path, fmt: str, *, duplicate_attr: bool = False) -> Any:
    readers = {
        "bed": pr.read_bed,
        "bedgraph": pr.read_bedgraph,
        "gff3": pr.read_gff3,
        "narrowPeak": pr.read_narrowPeak,
        "broadPeak": pr.read_broadPeak,
        "gappedPeak": pr.read_gappedPeak,
        "parquet": pr.read_parquet,
        "pairs": pr.read_pairs,
    }
    if fmt == "gtf":
        return _guard("read_gtf", lambda: pr.read_gtf(path, duplicate_attr=duplicate_attr))
    return _guard(f"read {fmt}", lambda: readers[fmt](path))


def _records(frame: Any, limit: int) -> list[dict[str, Any]]:
    return [
        {str(key): _c.plain(value) for key, value in row.items()}
        for row in frame.head(limit).to_dict("records")
    ]


def _frame_result(
    operation: str,
    frame: Any,
    limit: int,
    output_path: Path | None,
    **extra: Any,
) -> dict[str, Any]:
    if output_path is not None:
        frame.to_csv(str(output_path), sep="\t", index=False)
    return {
        "status": "completed",
        "operation": operation,
        "columns": [str(name) for name in frame.columns],
        "n_rows": len(frame),
        "rows": _records(frame, limit),
        "truncated": len(frame) > limit,
        "output_path": str(output_path) if output_path else None,
        **extra,
        "provenance": _provenance(),
    }


def _one_input(path: object, fmt: object) -> tuple[Any, Path, str]:
    file = input_path(path)  # type: ignore[arg-type]
    kind = _c.choice(fmt, _INPUT_FORMATS, "input_format")
    pr = _pyranges()
    return pr, file, kind


def read_pyranges1_intervals(
    path: str | Path,
    file_format: str,
    *,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read WIG (``wig``, pyranges1.read_wig), ENCODE broadPeak (``broadPeak``), gappedPeak
    (``gappedPeak``), Parquet interval tables (``parquet``) or UCSC chain files (``chain``,
    pyranges1.read_chain) into interval records.

    Rows carry ``Chromosome``, ``Start`` and ``End`` plus the format's own columns (for example
    ``Value`` for WIG, the peak statistics for the peak formats, chain header fields for chain).
    """

    file = input_path(path)
    kind = _c.choice(
        file_format, ("wig", "broadPeak", "gappedPeak", "parquet", "chain"), "file_format"
    )
    limit = _c.rows_limit(max_rows)
    pr = _pyranges()
    readers = {
        "wig": pr.read_wig,
        "broadPeak": pr.read_broadPeak,
        "gappedPeak": pr.read_gappedPeak,
        "parquet": pr.read_parquet,
        "chain": pr.read_chain,
    }
    gr = _guard(f"read_{kind}", lambda: readers[kind](file))
    return _frame_result("read_pyranges1_intervals", gr, limit, None, file_format=kind)


def has_pyranges1_tag(
    annotation_path: str | Path,
    tag: str,
    *,
    tag_column: str = "tag",
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test, for every row of a GTF annotation, whether it carries a tag such as
    ``Ensembl_canonical`` or ``MANE_Select`` (pyranges1.genes.has_tag).

    The GTF is read with ``duplicate_attr=True`` so that all tags of a row are kept. Whole tags are
    matched. Returns the flag per row (as ``has_tag``) and the number of rows that carry the tag.
    """

    file = input_path(annotation_path)
    name = _c.label(tag, "tag")
    column = _c.label(tag_column, "tag_column")
    limit = _c.rows_limit(max_rows)
    pr = _pyranges()
    from pyranges1 import genes

    gr = _load(pr, file, "gtf", duplicate_attr=True)
    flags = _guard("has_tag", lambda: genes.has_tag(gr, name, tag_col=column))
    frame = gr.copy()
    frame["has_tag"] = flags.to_numpy()
    return _frame_result(
        "has_pyranges1_tag", frame, limit, None, tag=name, n_tagged=int(flags.sum())
    )


def select_pyranges1_transcripts(
    annotation_path: str | Path,
    *,
    how: str = "canonical",
    gene_column: str = "gene_id",
    transcript_column: str = "transcript_id",
    tag_column: str = "tag",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep one representative transcript per gene from a GTF annotation
    (pyranges1.genes.select_transcripts).

    ``how``: ``canonical`` (the Ensembl_canonical transcript), ``mane`` (MANE_Select, else
    Ensembl_canonical) or ``longest`` (longest span). Genes without the tag fall back to the next
    rule. Rows of the chosen transcripts and rows without a transcript (gene lines) are kept.
    """

    file = input_path(annotation_path)
    rule = _c.choice(how, ("canonical", "mane", "longest"), "how")
    gene = _c.label(gene_column, "gene_column")
    transcript = _c.label(transcript_column, "transcript_column")
    tags = _c.label(tag_column, "tag_column")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pr = _pyranges()
    from pyranges1 import genes

    gr = _load(pr, file, "gtf", duplicate_attr=True)
    chosen = _guard(
        "select_transcripts",
        lambda: genes.select_transcripts(
            gr, rule, gene_col=gene, transcript_col=transcript, tag_col=tags
        ),
    )
    return _frame_result("select_pyranges1_transcripts", chosen, limit, out, how=rule)


def label_pyranges1_utrs(
    annotation_path: str | Path,
    *,
    transcript_column: str = "transcript_id",
    feature_column: str = "Feature",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Relabel unsplit ``UTR`` rows of a GTF as ``five_prime_utr`` or ``three_prime_utr``
    from the transcript's coding extent (pyranges1.genes.label_utrs).

    A UTR row wholly before the first CDS base (in transcription order) becomes the 5' UTR, wholly
    after the last CDS base the 3' UTR; UTR rows of transcripts without CDS keep ``UTR``.
    """

    file = input_path(annotation_path)
    transcript = _c.label(transcript_column, "transcript_column")
    feature = _c.label(feature_column, "feature_column")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    pr = _pyranges()
    from pyranges1 import genes

    gr = _load(pr, file, "gtf")
    labelled = _guard(
        "label_utrs",
        lambda: genes.label_utrs(gr, transcript_col=transcript, feature_col=feature),
    )
    counts = labelled[feature].astype(str).value_counts().to_dict()
    return _frame_result(
        "label_pyranges1_utrs",
        labelled,
        limit,
        out,
        feature_counts={str(key): int(value) for key, value in counts.items()},
    )


def outer_pyranges1_ranges(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    group_by: str | Sequence[str] | None = None,
    use_strand: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Outer boundaries of interval groups: the smallest interval that spans each group
    (PyRanges.outer_ranges; ``boundaries`` in pyranges 0.1).

    Groups are chromosome (and strand when ``use_strand`` is ``auto`` and strands exist) plus the
    columns in ``group_by`` (for example a transcript id).
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    groups = _group(group_by)
    strand = _strand(use_strand)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard("outer_ranges", lambda: gr.outer_ranges(group_by=groups, use_strand=strand))
    return _frame_result("outer_pyranges1_ranges", result, limit, out)


def resize_pyranges1_ranges(
    intervals_path: str | Path,
    width: int,
    *,
    fix: str = "center",
    input_format: str = "bed",
    use_strand: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resize intervals to a fixed width, keeping the start, end or centre fixed
    (PyRanges.resize_ranges).

    ``fix`` is ``start``, ``end`` or ``center``; with ``use_strand="auto"`` and strands present,
    ``start`` and ``end`` mean the 5' and 3' end of each interval.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    size = _c.integer(width, "width", 1, 2_000_000_000)
    anchor = _c.choice(fix, ("start", "end", "center"), "fix")
    strand = _strand(use_strand)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard("resize_ranges", lambda: gr.resize_ranges(size, anchor, use_strand=strand))
    return _frame_result("resize_pyranges1_ranges", result, limit, out)


def end_positions_pyranges1(
    intervals_path: str | Path,
    *,
    end: str = "five",
    input_format: str = "bed",
    group_by: str | Sequence[str] | None = None,
    extension: int = 0,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """5' end (``end="five"``, PyRanges.five_end) or 3' end (``end="three"``, PyRanges.three_end)
    of intervals, optionally of whole groups such as transcripts.

    Without ``group_by`` every interval is reduced to its end; with it, one end per group.
    ``extension`` widens the reported end by that many bases.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    which = _c.choice(end, ("five", "three"), "end")
    groups = _group(group_by)
    ext = _c.integer(extension, "extension", 0, 2_000_000_000)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    method = gr.five_end if which == "five" else gr.three_end
    result = _guard(f"{which}_end", lambda: method(group_by=groups, ext=ext))
    return _frame_result("end_positions_pyranges1", result, limit, out, end=which)


def join_pyranges1_overlaps(
    a_path: str | Path,
    b_path: str | Path,
    *,
    a_format: str = "bed",
    b_format: str = "bed",
    multiple: str = "all",
    strand_behavior: str = "auto",
    join_type: str = "inner",
    match_by: str | Sequence[str] | None = None,
    contained_intervals_only: bool = False,
    slack: int = 0,
    suffix: str = "_b",
    overlap_column: str | None = None,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join intervals of A with the B intervals they overlap, one output row per pair
    (PyRanges.join_overlaps; ``join`` in pyranges 0.1).

    ``join_type``: inner, left, right or outer. ``multiple``: all, first or last partner.
    ``slack`` widens A before testing. B columns get ``suffix``; ``overlap_column`` names an
    optional column with the overlap length.
    """

    first_path = input_path(a_path)
    second_path = input_path(b_path)
    pr = _pyranges()
    fmt_a = _c.choice(a_format, _INPUT_FORMATS, "a_format")
    fmt_b = _c.choice(b_format, _INPUT_FORMATS, "b_format")
    options: dict[str, Any] = {
        "multiple": _c.choice(multiple, ("all", "first", "last"), "multiple"),
        "strand_behavior": _c.choice(strand_behavior, _STRAND_BEHAVIOR, "strand_behavior"),
        "join_type": _c.choice(join_type, ("inner", "left", "right", "outer"), "join_type"),
        "match_by": _group(match_by, "match_by"),
        "contained_intervals_only": _c.flag(contained_intervals_only, "contained_intervals_only"),
        "slack": _c.integer(slack, "slack", 0, 2_000_000_000),
        "suffix": _c.label(suffix, "suffix"),
    }
    if overlap_column is not None:
        options["report_overlap_column"] = _c.label(overlap_column, "overlap_column")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    first = _load(pr, first_path, fmt_a)
    second = _load(pr, second_path, fmt_b)
    result = _guard("join_overlaps", lambda: first.join_overlaps(second, **options))
    return _frame_result("join_pyranges1_overlaps", result, limit, out)


def interval_lengths_pyranges1(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Length (End minus Start) of every interval (PyRanges.lengths).

    Returns the first ``max_rows`` lengths in input order plus count, minimum, maximum, mean and sum.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    lengths = _guard("lengths", lambda: gr.lengths())
    return {
        "status": "completed",
        "operation": "interval_lengths_pyranges1",
        "n_intervals": len(lengths),
        "lengths": [int(value) for value in lengths.head(limit).tolist()],
        "truncated": len(lengths) > limit,
        "min": _c.plain(lengths.min()) if len(lengths) else None,
        "max": _c.plain(lengths.max()) if len(lengths) else None,
        "mean": _c.plain(lengths.mean()) if len(lengths) else None,
        "sum": _c.plain(lengths.sum()) if len(lengths) else None,
        "provenance": _provenance(),
    }


def map_pyranges1_coordinates(
    query_path: str | Path,
    reference_path: str | Path,
    reference_column: str,
    *,
    direction: str = "to_global",
    query_format: str = "bed",
    reference_format: str = "gtf",
    local_column: str = "Chromosome",
    match_by: str | Sequence[str] | None = None,
    keep_id: bool = False,
    keep_chrom: bool = False,
    keep_loc: bool = False,
    pep_to_cds: bool = False,
    split_at_junctions: bool = True,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map intervals between a local reference frame (for example transcript coordinates) and
    global (genome) coordinates.

    ``direction="to_global"`` (PyRanges.map_to_global) takes query intervals whose
    ``local_column`` (default the chromosome field, holding the reference id) names rows of the
    reference and returns genome coordinates; ``keep_id``, ``keep_loc`` and ``pep_to_cds`` apply.
    ``direction="to_local"`` (PyRanges.map_to_local) takes genome intervals and returns coordinates
    local to the reference rows matched by ``match_by``; ``keep_chrom``, ``keep_loc`` and
    ``split_at_junctions`` apply. ``reference_column`` is the reference id column (``global_on``
    or ``ref_on``), for example ``transcript_id`` of a GTF.
    """

    query = input_path(query_path)
    reference = input_path(reference_path)
    way = _c.choice(direction, ("to_global", "to_local"), "direction")
    column = _c.label(reference_column, "reference_column")
    pr = _pyranges()
    fmt_q = _c.choice(query_format, _INPUT_FORMATS, "query_format")
    fmt_r = _c.choice(reference_format, _INPUT_FORMATS, "reference_format")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, query, fmt_q)
    ref = _load(pr, reference, fmt_r)
    if way == "to_global":
        result = _guard(
            "map_to_global",
            lambda: gr.map_to_global(
                ref,
                column,
                local_on=_c.label(local_column, "local_column"),
                keep_id=_c.flag(keep_id, "keep_id"),
                keep_loc=_c.flag(keep_loc, "keep_loc"),
                pep_to_cds=_c.flag(pep_to_cds, "pep_to_cds"),
            ),
        )
    else:
        result = _guard(
            "map_to_local",
            lambda: gr.map_to_local(
                ref,
                column,
                match_by=_group(match_by, "match_by"),
                keep_chrom=_c.flag(keep_chrom, "keep_chrom"),
                keep_loc=_c.flag(keep_loc, "keep_loc"),
                split_at_junctions=_c.flag(split_at_junctions, "split_at_junctions"),
            ),
        )
    return _frame_result("map_pyranges1_coordinates", result, limit, out, direction=way)


def max_disjoint_pyranges1_overlaps(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    use_strand: str = "auto",
    slack: int = 0,
    match_by: str | Sequence[str] | None = None,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Largest set of intervals that do not overlap each other (PyRanges.max_disjoint_overlaps).

    ``slack`` also separates intervals closer than that many bases; ``match_by`` restricts the
    disjointness to rows with equal values in the given columns.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    strand = _strand(use_strand)
    gap = _c.integer(slack, "slack", 0, 2_000_000_000)
    keys = _group(match_by, "match_by")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard(
        "max_disjoint_overlaps",
        lambda: gr.max_disjoint_overlaps(use_strand=strand, slack=gap, match_by=keys),
    )
    return _frame_result("max_disjoint_pyranges1_overlaps", result, limit, out)


def slice_pyranges1_ranges(
    intervals_path: str | Path,
    *,
    start: int = 0,
    end: int | None = None,
    input_format: str = "bed",
    group_by: str | Sequence[str] | None = None,
    use_strand: str = "auto",
    count_introns: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut intervals at start/end offsets, genomic or spliced (PyRanges.slice_ranges;
    ``spliced_subsequence`` in pyranges 0.1).

    Without ``group_by`` each interval is sliced alone (offsets are within the interval). With
    ``group_by`` (for example ``transcript_id``) the rows of a group are one spliced sequence:
    offsets count transcript bases, introns are skipped unless ``count_introns``. Negative offsets
    count from the end; minus-strand groups are counted from their 5' end when ``use_strand`` is
    ``auto``.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    begin = _c.integer(start, "start", -2_000_000_000, 2_000_000_000)
    stop = None if end is None else _c.integer(end, "end", -2_000_000_000, 2_000_000_000)
    groups = _group(group_by)
    strand = _strand(use_strand)
    counting = _c.flag(count_introns, "count_introns")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard(
        "slice_ranges",
        lambda: gr.slice_ranges(
            begin, stop, group_by=groups, use_strand=strand, count_introns=counting
        ),
    )
    return _frame_result("slice_pyranges1_ranges", result, limit, out)


def make_pyranges1_strand_valid(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace invalid strand values (anything but ``+`` and ``-``, such as ``.``) by ``+``
    (PyRanges.make_strand_valid).

    Returns the repaired intervals and how many strand values changed.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    if "Strand" not in gr.columns:
        raise ConfigurationError("The intervals have no Strand column.")
    before = gr["Strand"].astype(str)
    result = _guard("make_strand_valid", lambda: gr.make_strand_valid())
    changed = int((before.to_numpy() != result["Strand"].astype(str).to_numpy()).sum())
    return _frame_result("make_pyranges1_strand_valid", result, limit, out, n_changed=changed)


def summarize_pyranges1_intervals(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary statistics of interval lengths and of the merged coverage (PyRanges.summary).

    Returns the table pyranges1 prints: the number of intervals, total and merged (non-overlapping)
    bases, and the distribution of interval lengths.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    gr = _load(pr, file, kind)
    table = _guard("summary", lambda: gr.summary(return_df=True))
    return {
        "status": "completed",
        "operation": "summarize_pyranges1_intervals",
        "columns": [str(name) for name in table.columns],
        "index": [str(name) for name in table.index],
        "rows": _records(table.reset_index(), 1000),
        "provenance": _provenance(),
    }


def write_pyranges1_intervals(
    intervals_path: str | Path,
    output_path: str | Path,
    output_format: str,
    *,
    input_format: str = "bed",
    value_column: str = "Value",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write intervals as a Parquet table (``parquet``, PyRanges.to_parquet), ENCODE ``narrowPeak``,
    ``broadPeak`` or ``gappedPeak`` (PyRanges.to_narrowPeak, to_broadPeak, to_gappedPeak), a Hi-C
    ``pairs`` file (PyRanges.to_pairs) or a ``bedgraph`` file (PyRanges.to_bedgraph).

    The input must carry the columns the output format needs (for example the peak statistics for
    the peak formats, which a narrowPeak, broadPeak or gappedPeak input provides). For ``bedgraph``
    ``value_column`` names the signal column (``Value`` by default, as WIG and bigWig readers
    produce it). The new output file is written in the requested format.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    fmt = _c.choice(output_format, _WRITE_FORMATS, "output_format")
    column = _c.label(value_column, "value_column")
    out = _c.new_output(output_path)
    assert out is not None
    gr = _load(pr, file, kind)
    methods = {
        "parquet": lambda: gr.to_parquet(str(out)),
        "narrowPeak": lambda: gr.to_narrowPeak(str(out)),
        "broadPeak": lambda: gr.to_broadPeak(str(out)),
        "gappedPeak": lambda: gr.to_gappedPeak(str(out)),
        "pairs": lambda: gr.to_pairs(str(out)),
        "bedgraph": lambda: gr.to_bedgraph(str(out), value_col=column),
    }
    _guard(f"write {fmt}", methods[fmt])
    return {
        "status": "completed",
        "operation": "write_pyranges1_intervals",
        "output_format": fmt,
        "n_rows": len(gr),
        "output_path": str(out),
        "provenance": _provenance(),
    }


def bigwig_stats_pyranges1(
    intervals_path: str | Path,
    bigwig_path: str | Path,
    *,
    input_format: str = "bed",
    stats: str | Sequence[str] = "mean",
    missing: float | None = None,
    pad: int = 0,
    bins: int | None = None,
    use_strand: str = "auto",
    prefix: str = "",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise bigWig signal over each interval (PyRanges.bigwig_stats).

    ``stats`` is any of mean, max, min, sum, std (one name or a list); one column per statistic is
    added. ``missing`` is the value for bases without data (None leaves them out; 0.0 gives the
    mean over the whole interval), ``pad`` widens each interval, and ``bins`` splits each interval
    into equal bins, one column per bin and statistic (a regions x bins matrix). Needs pyBigWig,
    runlength and pyarrow.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    track = input_path(bigwig_path)
    wanted = [stats] if isinstance(stats, str) else list(stats)
    allowed = ("mean", "max", "min", "sum", "std")
    chosen = [_c.choice(item, allowed, "stats") for item in wanted]
    options: dict[str, Any] = {
        "stats": chosen,
        "pad": _c.integer(pad, "pad", 0, 2_000_000_000),
        "use_strand": _strand(use_strand),
        "prefix": "" if not prefix else _c.label(prefix, "prefix"),
    }
    if missing is not None:
        options["missing"] = _c.number(missing, "missing", -1e300, 1e300)
    if bins is not None:
        options["bins"] = _c.integer(bins, "bins", 1, 10_000)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard("bigwig_stats", lambda: gr.bigwig_stats(track, **options))
    return _frame_result("bigwig_stats_pyranges1", result, limit, out)


def coverage_to_bigwig_pyranges1(
    intervals_path: str | Path,
    output_path: str | Path,
    *,
    input_format: str = "bed",
    chromsizes_path: str | Path | None = None,
    value_column: str | None = None,
    divide: bool = False,
    rpm: bool = True,
    precomputed: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute coverage of intervals and write it as a bigWig file (PyRanges.to_bigwig).

    The score at each position is the number of intervals covering it (or the sum of
    ``value_column``); ``rpm`` scales to reads per million, ``divide`` (with ``value_column``)
    divides by plain coverage and takes log2, and ``precomputed`` writes the input ranges with their
    ``value_column`` as they are (they must not overlap). ``chromsizes_path`` is a two-column size
    file; without it the largest end per chromosome is used. Needs pyBigWig.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    out = _c.new_output(output_path)
    assert out is not None
    sizes = None
    if chromsizes_path is not None:
        sizes = _c.check_genome_file(input_path(chromsizes_path))
    column = None if value_column is None else _c.label(value_column, "value_column")
    options = {
        "divide": _c.flag(divide, "divide"),
        "rpm": _c.flag(rpm, "rpm"),
        "precomputed": _c.flag(precomputed, "precomputed"),
    }
    gr = _load(pr, file, kind)
    _guard(
        "to_bigwig",
        lambda: gr.to_bigwig(str(out), chromosome_sizes=sizes, value_col=column, **options),
    )
    return {
        "status": "completed",
        "operation": "coverage_to_bigwig_pyranges1",
        "n_rows": len(gr),
        "output_path": str(out),
        "provenance": _provenance(),
    }


def coverage_runs_pyranges1(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    value_column: str | None = None,
    use_strand: str = "auto",
    rpm: bool = False,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage (or the sum of ``value_column``) as run-length encoded sequences per chromosome and
    strand (PyRanges.to_rle).

    Returns one row per run (``Chromosome``, ``Strand``, ``start``, ``end``, ``value``; ``rpm``
    scales to reads per million). Needs the runlength package.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    column = None if value_column is None else _c.label(value_column, "value_column")
    strand = _strand(use_strand)
    scale = _c.flag(rpm, "rpm")
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    batch = _guard("to_rle", lambda: gr.to_rle(value_col=column, strand=strand, rpm=scale))
    frame = _guard("to_rle table", lambda: batch.to_pandas())
    return _frame_result("coverage_runs_pyranges1", frame, limit, None)


def flank_pyranges1_ranges(
    intervals_path: str | Path,
    length: int,
    *,
    side: str = "upstream",
    gap: int = 0,
    input_format: str = "bed",
    group_by: str | Sequence[str] | None = None,
    use_strand: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Regions of ``length`` bases upstream (``side="upstream"``, PyRanges.upstream) or downstream
    (``side="downstream"``, PyRanges.downstream) of each interval or interval group.

    ``gap`` leaves that many bases between the interval and the flank; with strands present
    (``use_strand="auto"``) upstream means 5' of the interval.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    size = _c.integer(length, "length", 1, 2_000_000_000)
    which = _c.choice(side, ("upstream", "downstream"), "side")
    spacing = _c.integer(gap, "gap", 0, 2_000_000_000)
    groups = _group(group_by)
    strand = _strand(use_strand)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    method = gr.upstream if which == "upstream" else gr.downstream
    result = _guard(which, lambda: method(size, spacing, group_by=groups, use_strand=strand))
    return _frame_result("flank_pyranges1_ranges", result, limit, out, side=which)


def window_pyranges1_ranges(
    intervals_path: str | Path,
    window_size: int,
    *,
    input_format: str = "bed",
    use_strand: str = "auto",
    group_by: str | Sequence[str] | None = None,
    add_window_id: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split intervals into windows of ``window_size`` anchored at the 5' end (PyRanges.window_ranges;
    ``window`` in pyranges 0.1).

    The last window of an interval may be shorter. ``add_window_id`` adds a window number column.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    size = _c.integer(window_size, "window_size", 1, 2_000_000_000)
    strand = _strand(use_strand)
    groups = _group(group_by)
    ident = _c.flag(add_window_id, "add_window_id")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard(
        "window_ranges",
        lambda: gr.window_ranges(size, use_strand=strand, group_by=groups, add_window_id=ident),
    )
    return _frame_result("window_pyranges1_ranges", result, limit, out)


def group_cumsum_pyranges1(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    group_by: str | Sequence[str] | None = None,
    use_strand: str = "auto",
    start_column: str | None = None,
    end_column: str | None = None,
    keep_order: bool = True,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cumulative interval length within each group, in transcription order for minus-strand groups
    (PyRanges.group_cumsum).

    Adds two columns (default ``CumsumStart`` and ``CumsumEnd``, or ``start_column`` and
    ``end_column``) giving each interval's start and end offset within the concatenated group, for
    example exon offsets within a transcript when ``group_by="transcript_id"``.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    groups = _group(group_by)
    strand = _strand(use_strand)
    first = None if start_column is None else _c.label(start_column, "start_column")
    last = None if end_column is None else _c.label(end_column, "end_column")
    order = _c.flag(keep_order, "keep_order")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    result = _guard(
        "group_cumsum",
        lambda: gr.group_cumsum(
            group_by=groups,
            use_strand=strand,
            cumsum_start_column=first,
            cumsum_end_column=last,
            keep_order=order,
        ),
    )
    return _frame_result("group_pyranges1_cumsum", result, limit, out)


def interval_metrics_pyranges1(
    a_path: str | Path,
    b_path: str | Path,
    *,
    a_format: str = "bed",
    b_format: str = "bed",
    metrics: str | Sequence[str] = "fraction",
    denominator: str = "first",
    strand_behavior: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Relationship metrics for every overlapping pair of A and B intervals (join_overlaps followed
    by PyRanges.compute_interval_metrics).

    ``metrics`` is one name or a list from pyranges1 (for example fraction, jaccard, distance); each
    becomes a column. ``denominator`` (first, second, union or min/max as pyranges1 documents) sets
    the reference length for ``fraction``.
    """

    first_path = input_path(a_path)
    second_path = input_path(b_path)
    pr = _pyranges()
    fmt_a = _c.choice(a_format, _INPUT_FORMATS, "a_format")
    fmt_b = _c.choice(b_format, _INPUT_FORMATS, "b_format")
    wanted = [metrics] if isinstance(metrics, str) else list(metrics)
    names = [_c.label(item, "metrics") for item in wanted]
    denom = _c.label(denominator, "denominator")
    behavior = _c.choice(strand_behavior, _STRAND_BEHAVIOR, "strand_behavior")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    first = _load(pr, first_path, fmt_a)
    second = _load(pr, second_path, fmt_b)
    joined = _guard("join_overlaps", lambda: first.join_overlaps(second, strand_behavior=behavior))
    result = _guard(
        "compute_interval_metrics", lambda: joined.compute_interval_metrics(names, denom=denom)
    )
    return _frame_result("interval_metrics_pyranges1", result, limit, out, metrics=names)


def combine_pyranges1_interval_columns(
    a_path: str | Path,
    b_path: str | Path,
    *,
    function: str = "intersect",
    a_format: str = "bed",
    b_format: str = "bed",
    drop_old_columns: bool = True,
    strand_behavior: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the start/end pairs of overlapping A and B intervals into a new start and end
    (join_overlaps followed by PyRanges.combine_interval_columns).

    ``function``: ``intersect`` (the overlap), ``union`` (the span covering both) or ``swap``
    (exchange the A and B coordinates). ``drop_old_columns`` removes the B coordinate columns.
    """

    first_path = input_path(a_path)
    second_path = input_path(b_path)
    pr = _pyranges()
    fmt_a = _c.choice(a_format, _INPUT_FORMATS, "a_format")
    fmt_b = _c.choice(b_format, _INPUT_FORMATS, "b_format")
    how = _c.choice(function, ("intersect", "union", "swap"), "function")
    drop = _c.flag(drop_old_columns, "drop_old_columns")
    behavior = _c.choice(strand_behavior, _STRAND_BEHAVIOR, "strand_behavior")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    first = _load(pr, first_path, fmt_a)
    second = _load(pr, second_path, fmt_b)
    joined = _guard("join_overlaps", lambda: first.join_overlaps(second, strand_behavior=behavior))
    result = _guard(
        "combine_interval_columns",
        lambda: joined.combine_interval_columns(how, drop_old_columns=drop),
    )
    return _frame_result("combine_pyranges1_interval_columns", result, limit, out, function=how)


def liftover_pyranges1_ranges(
    intervals_path: str | Path,
    chain_path: str | Path,
    *,
    input_format: str = "bed",
    min_match: float = 0.95,
    return_unmapped: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert intervals between assemblies with a UCSC chain file (PyRanges.liftover_ranges).

    ``min_match`` is the minimum fraction of an interval that must map. With ``return_unmapped``
    the intervals that failed to map are returned as well (``unmapped_rows``).
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    chain = input_path(chain_path)
    fraction = _c.number(min_match, "min_match", 0.0, 1.0)
    unmapped_wanted = _c.flag(return_unmapped, "return_unmapped")
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    outcome = _guard(
        "liftover_ranges",
        lambda: gr.liftover_ranges(chain, min_match=fraction, return_unmapped=unmapped_wanted),
    )
    extra: dict[str, Any] = {}
    if unmapped_wanted:
        mapped, unmapped = outcome
        extra = {"n_unmapped": len(unmapped), "unmapped_rows": _records(unmapped, limit)}
    else:
        mapped = outcome
    return _frame_result("liftover_pyranges1_ranges", mapped, limit, out, **extra)


def drop_pyranges1_duplicate_positions(
    intervals_path: str | Path,
    *,
    input_format: str = "bed",
    keep: str = "first",
    use_strand: str = "auto",
    output_path: str | Path | None = None,
    max_rows: int = _c.DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop rows that repeat a genomic position, comparing Chromosome, Start, End and (when
    present and ``use_strand`` is ``auto`` or ``yes``) Strand (pyranges 0.1
    ``drop_duplicate_positions``, here PyRanges.drop_duplicates on the position columns).

    ``keep``: ``first`` or ``last`` row of each duplicate set, or ``none`` to drop all rows of
    duplicated positions.
    """

    pr, file, kind = _one_input(intervals_path, input_format)
    which = _c.choice(keep, ("first", "last", "none"), "keep")
    strand = _strand(use_strand)
    out = _c.new_output(output_path)
    limit = _c.rows_limit(max_rows)
    gr = _load(pr, file, kind)
    columns = ["Chromosome", "Start", "End"]
    if "Strand" in gr.columns and strand is not False:
        columns.append("Strand")
    result = _guard(
        "drop_duplicates",
        lambda: gr.drop_duplicates(subset=columns, keep=False if which == "none" else which),
    )
    return _frame_result(
        "drop_pyranges1_duplicate_positions",
        result,
        limit,
        out,
        n_dropped=int(len(gr) - len(result)),
        compared_columns=columns,
    )


__all__ = [
    "bigwig_stats_pyranges1",
    "combine_pyranges1_interval_columns",
    "coverage_runs_pyranges1",
    "coverage_to_bigwig_pyranges1",
    "drop_pyranges1_duplicate_positions",
    "end_positions_pyranges1",
    "flank_pyranges1_ranges",
    "group_cumsum_pyranges1",
    "has_pyranges1_tag",
    "interval_lengths_pyranges1",
    "interval_metrics_pyranges1",
    "join_pyranges1_overlaps",
    "label_pyranges1_utrs",
    "liftover_pyranges1_ranges",
    "make_pyranges1_strand_valid",
    "map_pyranges1_coordinates",
    "max_disjoint_pyranges1_overlaps",
    "outer_pyranges1_ranges",
    "read_pyranges1_intervals",
    "resize_pyranges1_ranges",
    "select_pyranges1_transcripts",
    "slice_pyranges1_ranges",
    "summarize_pyranges1_intervals",
    "window_pyranges1_ranges",
    "write_pyranges1_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
