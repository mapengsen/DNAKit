"""pysam alignment-record operations (pack P36): mates, pileup columns, CIGAR geometry, name lookup, headers.

Each public function opens a local SAM, BAM or CRAM file with pysam and calls the pysam method that carries
the operation, with the caller's arguments and pysam's own defaults. Coordinates follow pysam: contig names,
0-based half-open ``start``/``end``. Region-based functions need an index (BAI, CSI or CRAI) when ``contig``
is given; with ``contig=None`` the whole file is read sequentially. Result lists are cut at ``max_records``
(the result then says ``truncated: true``). Output files are written only to new paths. pysam is imported
when a function runs; a missing package raises ``BackendUnavailableError`` with the ``bio-reads`` hint.
Ordinary tests use only small files built on the fly.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-reads"
_MAX_RECORDS = 1_000_000
_ERRORS = (
    ValueError,
    KeyError,
    OSError,
    TypeError,
    AttributeError,
    IndexError,
    NotImplementedError,
)


def _pysam() -> Any:
    return require("pysam", _EXTRA)


def _provenance(operation: str) -> dict[str, str]:
    return {
        "package": "pysam",
        "version": importlib.metadata.version("pysam"),
        "operation": operation,
        "coordinates": "0-based half-open (pysam)",
    }


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return str(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _optional_region(
    contig: object, start: object, end: object
) -> tuple[str | None, int | None, int | None]:
    if contig is None:
        if start is not None or end is not None:
            raise ConfigurationError("start and end need a contig.")
        return None, None, None
    chrom = _text(contig, "contig")
    first = None if start is None else _int(start, "start", 0, 2**31 - 1)
    last = None if end is None else _int(end, "end", 0, 2**31 - 1)
    if first is not None and last is not None and last < first:
        raise ConfigurationError("end must not be smaller than start.")
    return chrom, first, last


def _new_output(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}")
    return path


def _open(path: str | Path, name: str = "alignment_file") -> tuple[Any, Any]:
    pysam = _pysam()
    resolved = input_path(path)
    try:
        return pysam, pysam.AlignmentFile(str(resolved), check_sq=False)
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam could not open {name} {resolved}: {exc}") from exc


def _reads(
    handle: Any, contig: str | None, start: int | None, end: int | None, **kwargs: Any
) -> Iterator[Any]:
    try:
        if contig is None:
            return iter(handle.fetch(until_eof=True, **kwargs))
        return iter(handle.fetch(contig, start, end, **kwargs))
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam fetch failed: {exc}") from exc


def _limited(items: Iterator[Any], limit: int) -> tuple[list[Any], bool]:
    out: list[Any] = []
    for item in items:
        if len(out) >= limit:
            return out, True
        out.append(item)
    return out, False


def _plain(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return value.tolist()
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    return value


def _per_read(
    operation: str,
    alignment_file: str | Path,
    contig: str | None,
    start: int | None,
    end: int | None,
    max_records: int,
    extract: Any,
    **extra: Any,
) -> dict[str, Any]:
    chrom, first, last = _optional_region(contig, start, end)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    _pysam, handle = _open(alignment_file)
    with handle:
        reads, truncated = _limited(_reads(handle, chrom, first, last), limit)
        rows = []
        for read in reads:
            try:
                value = extract(read)
            except _ERRORS as exc:
                raise BackendExecutionError(
                    f"pysam {operation} failed on read {read.query_name}: {exc}"
                ) from exc
            rows.append({"query_name": read.query_name, "flag": int(read.flag), **value})
    return {
        "status": "completed",
        "operation": operation,
        "source": str(input_path(alignment_file)),
        "region": {"contig": chrom, "start": first, "end": last},
        "count": len(rows),
        "truncated": truncated,
        "records": rows,
        "provenance": _provenance(operation) | extra,
    }


def fetch_mate_alignments(
    alignment_file: str | Path,
    contig: str,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the mate alignment of every paired, mapped read in a region (``AlignmentFile.mate``).

    The file must be indexed and coordinate sorted. Reads that are unpaired or whose mate is unmapped are
    listed in ``without_mate``; every other read gets the mate's ``pysam`` record dictionary.
    """

    chrom, first, last = _optional_region(contig, start, end)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    _pysam, handle = _open(alignment_file)
    pairs: list[dict[str, Any]] = []
    without: list[str] = []
    truncated = False
    with handle:
        for read in _reads(handle, chrom, first, last):
            if not read.is_paired or read.mate_is_unmapped:
                without.append(read.query_name)
                continue
            if len(pairs) >= limit:
                truncated = True
                break
            try:
                mate = handle.mate(read)
            except _ERRORS as exc:
                raise BackendExecutionError(
                    f"pysam mate failed for {read.query_name}: {exc}"
                ) from exc
            pairs.append({"read": read.to_dict(), "mate": mate.to_dict()})
    return {
        "status": "completed",
        "operation": "fetch_mate_alignments",
        "region": {"contig": chrom, "start": first, "end": last},
        "count": len(pairs),
        "truncated": truncated,
        "pairs": pairs,
        "without_mate": without,
        "provenance": _provenance("AlignmentFile.mate"),
    }


def count_base_coverage(
    alignment_file: str | Path,
    contig: str,
    start: int,
    end: int,
    *,
    quality_threshold: int = 15,
    read_callback: str = "all",
    verbose: bool = True,
) -> dict[str, Any]:
    """Count A, C, G and T bases per position of a region (``AlignmentFile.count_coverage``).

    ``quality_threshold`` skips bases below this Phred quality (pysam default 15); ``read_callback`` is
    ``all`` (skip unmapped, secondary, QC-fail and duplicate reads) or ``nofilter``. The four lists have one
    count per position from ``start`` to ``end``.
    """

    chrom, first, last = _optional_region(contig, start, end)
    if first is None or last is None or last <= first:
        raise ConfigurationError("start and end are required and end must be greater than start.")
    threshold = _int(quality_threshold, "quality_threshold", 0, 255)
    callback = _choice(read_callback, "read_callback", ("all", "nofilter"))
    _pysam, handle = _open(alignment_file)
    with handle:
        try:
            a, c, g, t = handle.count_coverage(
                chrom, first, last, quality_threshold=threshold, read_callback=callback
            )
        except _ERRORS as exc:
            raise BackendExecutionError(f"pysam count_coverage failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "count_base_coverage",
        "region": {"contig": chrom, "start": first, "end": last},
        "quality_threshold": threshold,
        "read_callback": callback,
        "A": [int(v) for v in a],
        "C": [int(v) for v in c],
        "G": [int(v) for v in g],
        "T": [int(v) for v in t],
        "provenance": _provenance("AlignmentFile.count_coverage"),
    }


def count_spliced_junctions(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    method: str = "fast",
    verbose: bool = True,
) -> dict[str, Any]:
    """Count reads that support each intron (CIGAR ``N`` gap) in a region (``AlignmentFile.find_introns``).

    ``method`` is ``fast`` (``find_introns``) or ``slow`` (``find_introns_slow``, which collects the
    reads' aligned blocks); both give the same junctions. Each junction is the 0-based half-open
    interval of the intron with its read count.
    """

    chrom, first, last = _optional_region(contig, start, end)
    mode = _choice(method, "method", ("fast", "slow"))
    _pysam, handle = _open(alignment_file)
    with handle:
        reads = _reads(handle, chrom, first, last)
        try:
            counter = (
                handle.find_introns(reads) if mode == "fast" else handle.find_introns_slow(reads)
            )
        except _ERRORS as exc:
            raise BackendExecutionError(f"pysam find_introns failed: {exc}") from exc
    junctions = [
        {"intron_start": int(s), "intron_end": int(e), "supporting_reads": int(n)}
        for (s, e), n in sorted(counter.items())
    ]
    return {
        "status": "completed",
        "operation": "count_spliced_junctions",
        "region": {"contig": chrom, "start": first, "end": last},
        "method": mode,
        "junction_count": len(junctions),
        "junctions": junctions,
        "provenance": _provenance("AlignmentFile.find_introns"),
    }


def pileup_columns(
    alignment_file: str | Path,
    contig: str,
    start: int | None = None,
    end: int | None = None,
    *,
    min_base_quality: int = 13,
    min_mapping_quality: int = 0,
    stepper: str = "all",
    truncate: bool = False,
    max_depth: int = 8000,
    ignore_overlaps: bool = True,
    ignore_orphans: bool = True,
    include_read_names: bool = False,
    max_columns: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Iterate the pileup columns of a region and report the reads aligned at each (``AlignmentFile.pileup``).

    For every column: ``reference_pos`` (0-based), ``nsegments`` (reads covering it, before the base-quality
    filter), ``num_aligned`` (``PileupColumn.get_num_aligned`` - reads whose base passes
    ``min_base_quality``), and optionally the read names. ``stepper`` is ``all``, ``nofilter``, ``samtools``
    or ``bam``; ``truncate=False`` is pysam's default and may list columns just outside the region.
    """

    chrom, first, last = _optional_region(contig, start, end)
    quality = _int(min_base_quality, "min_base_quality", 0, 255)
    mapq = _int(min_mapping_quality, "min_mapping_quality", 0, 255)
    steps = _choice(stepper, "stepper", ("all", "nofilter", "samtools", "bam"))
    depth = _int(max_depth, "max_depth", 1, 10_000_000)
    limit = _int(max_columns, "max_columns", 1, _MAX_RECORDS)
    names = _flag(include_read_names, "include_read_names")
    _pysam, handle = _open(alignment_file)
    columns: list[dict[str, Any]] = []
    truncated = False
    with handle:
        try:
            for column in handle.pileup(
                chrom,
                first,
                last,
                truncate=_flag(truncate, "truncate"),
                stepper=steps,
                max_depth=depth,
                min_base_quality=quality,
                min_mapping_quality=mapq,
                ignore_overlaps=_flag(ignore_overlaps, "ignore_overlaps"),
                ignore_orphans=_flag(ignore_orphans, "ignore_orphans"),
            ):
                if len(columns) >= limit:
                    truncated = True
                    break
                row: dict[str, Any] = {
                    "reference_pos": int(column.reference_pos),
                    "nsegments": int(column.nsegments),
                    "num_aligned": int(column.get_num_aligned()),
                }
                if names:
                    row["query_names"] = list(column.get_query_names())
                columns.append(row)
        except _ERRORS as exc:
            raise BackendExecutionError(f"pysam pileup failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "pileup_columns",
        "region": {"contig": chrom, "start": first, "end": last},
        "min_base_quality": quality,
        "column_count": len(columns),
        "truncated": truncated,
        "columns": columns,
        "provenance": _provenance("AlignmentFile.pileup"),
    }


def find_alignments_by_name(
    alignment_file: str | Path,
    query_names: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up all alignments of the given read names with pysam's in-memory name index (``IndexedReads``).

    The whole file is read once to build the index, so it needs no BAM index and works on unsorted files.
    Names without an alignment map to an empty list.
    """

    if isinstance(query_names, (str, bytes)) or not isinstance(query_names, Sequence):
        raise ConfigurationError("query_names must be a list of read names.")
    wanted = [_text(name, "query_names item") for name in query_names]
    if not wanted or len(wanted) > 100_000:
        raise ConfigurationError("query_names must hold between 1 and 100000 names.")
    pysam, handle = _open(alignment_file)
    found: dict[str, list[dict[str, Any]]] = {}
    with handle:
        try:
            index = pysam.IndexedReads(handle)
            index.build()
            for name in wanted:
                try:
                    found[name] = [read.to_dict() for read in index.find(name)]
                except KeyError:
                    found[name] = []
        except _ERRORS as exc:
            raise BackendExecutionError(f"pysam IndexedReads failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "find_alignments_by_name",
        "alignments": found,
        "matched_names": sum(1 for rows in found.values() if rows),
        "provenance": _provenance("IndexedReads"),
    }


def sam_header_to_dict(alignment_file: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Return the header of a SAM, BAM or CRAM file as nested dictionaries (``AlignmentHeader.to_dict``)."""

    _pysam, handle = _open(alignment_file)
    with handle:
        header = handle.header.to_dict()
        text = str(handle.header)
    return {
        "status": "completed",
        "operation": "sam_header_to_dict",
        "header": _plain(header),
        "header_text": text,
        "provenance": _provenance("AlignmentHeader.to_dict"),
    }


def write_sam_header(
    output_path: str | Path,
    *,
    header_dict: Mapping[str, Any] | None = None,
    reference_names: Sequence[str] | None = None,
    reference_lengths: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a SAM header and write it as text (``AlignmentHeader.from_dict`` / ``from_references``).

    Give either ``header_dict`` (the nested form of ``sam_header_to_dict``: ``HD``, ``SQ`` list, ``RG`` list,
    ``PG`` list, ``CO`` list) or ``reference_names`` with ``reference_lengths`` (one ``@SQ`` line each).
    The result holds the header text and the round-tripped dictionary.
    """

    target = _new_output(output_path, "output_path")
    pysam = _pysam()
    try:
        if header_dict is not None:
            if reference_names is not None or reference_lengths is not None:
                raise ConfigurationError(
                    "Give header_dict or reference_names/reference_lengths, not both."
                )
            header = pysam.AlignmentHeader.from_dict(dict(header_dict))
            mode = "from_dict"
        else:
            if reference_names is None or reference_lengths is None:
                raise ConfigurationError(
                    "Give header_dict, or both reference_names and reference_lengths."
                )
            if len(reference_names) != len(reference_lengths) or not reference_names:
                raise ConfigurationError(
                    "reference_names and reference_lengths need the same length (> 0)."
                )
            names = [_text(n, "reference_names item") for n in reference_names]
            lengths = [_int(n, "reference_lengths item", 1, 2**31 - 1) for n in reference_lengths]
            header = pysam.AlignmentHeader.from_references(names, lengths)
            mode = "from_references"
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam header construction failed: {exc}") from exc
    text = str(header)
    target.write_text(text, encoding="utf-8")
    return {
        "status": "completed",
        "operation": "write_sam_header",
        "mode": mode,
        "output": str(target),
        "header_text": text,
        "header": _plain(header.to_dict()),
        "provenance": _provenance(f"AlignmentHeader.{mode}"),
    }


def alignment_records_to_dicts(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert alignment records to SAM-style dictionaries (``AlignedSegment.to_dict``).

    Each record holds ``name, flag, ref_name, ref_pos`` (1-based), ``map_quality, cigar, next_ref_name,
    next_ref_pos, length, seq, qual`` and ``tags`` - all as SAM text, so the output can be written back with
    ``write_alignments_from_dicts``.
    """

    chrom, first, last = _optional_region(contig, start, end)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    _pysam, handle = _open(alignment_file)
    with handle:
        reads, truncated = _limited(_reads(handle, chrom, first, last), limit)
        rows = [read.to_dict() for read in reads]
        header = _plain(handle.header.to_dict())
    return {
        "status": "completed",
        "operation": "alignment_records_to_dicts",
        "region": {"contig": chrom, "start": first, "end": last},
        "count": len(rows),
        "truncated": truncated,
        "header": header,
        "records": rows,
        "provenance": _provenance("AlignedSegment.to_dict"),
    }


def write_alignments_from_dicts(
    records: Sequence[Mapping[str, Any]],
    header: Mapping[str, Any],
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build alignment records from SAM-style dictionaries and write them (``AlignedSegment.from_dict``).

    ``header`` is the nested header dictionary; ``records`` use the keys of ``alignment_records_to_dicts``.
    The output format follows the suffix (``.sam`` or ``.bam``).
    """

    target = _new_output(output_path, "output_path")
    mode = {".sam": "w", ".bam": "wb"}.get(target.suffix.lower())
    if mode is None:
        raise ConfigurationError("output_path must end with .sam or .bam.")
    if isinstance(records, (str, bytes)) or not isinstance(records, Sequence):
        raise ConfigurationError("records must be a list of dictionaries.")
    if len(records) > _MAX_RECORDS:
        raise ConfigurationError(f"records holds more than {_MAX_RECORDS} entries.")
    pysam = _pysam()
    try:
        header_obj = pysam.AlignmentHeader.from_dict(dict(header))
        built = [
            pysam.AlignedSegment.from_dict(
                {k: str(v) if k != "tags" else v for k, v in dict(r).items()}, header_obj
            )
            for r in records
        ]
        with pysam.AlignmentFile(str(target), mode, header=header_obj) as out:
            for read in built:
                out.write(read)
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam could not build the alignments: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_alignments_from_dicts",
        "output": str(target),
        "count": len(built),
        "provenance": _provenance("AlignedSegment.from_dict"),
    }


def read_aligned_pairs(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    matches_only: bool = False,
    with_seq: bool = False,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the (read position, reference position) pairs of each read (``AlignedSegment.get_aligned_pairs``).

    Insertions and soft clips have ``None`` as the reference position, deletions and skipped regions ``None``
    as the read position, unless ``matches_only=True``. ``with_seq=True`` adds the reference base (needs an MD
    tag) as third element.
    """

    only = _flag(matches_only, "matches_only")
    seq = _flag(with_seq, "with_seq")
    return _per_read(
        "read_aligned_pairs",
        alignment_file,
        contig,
        start,
        end,
        max_records,
        lambda read: {"pairs": _plain(read.get_aligned_pairs(matches_only=only, with_seq=seq))},
        matches_only=only,
        with_seq=seq,
    )


def read_aligned_blocks(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the gapless aligned blocks of each read in reference coordinates (``AlignedSegment.get_blocks``).

    A block is a 0-based half-open ``[start, end)`` interval; deletions and introns split blocks.
    """

    return _per_read(
        "read_aligned_blocks",
        alignment_file,
        contig,
        start,
        end,
        max_records,
        lambda read: {"blocks": [list(b) for b in read.get_blocks()]},
    )


def read_cigar_stats(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the CIGAR of each read: bases and blocks per operation, plus the NM edit distance
    (``AlignedSegment.get_cigar_stats``).

    ``bases`` and ``blocks`` map the operation names M, I, D, N, S, H, P, =, X and NM to counts.
    """

    labels = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "NM"]

    def extract(read: Any) -> dict[str, Any]:
        bases, blocks = read.get_cigar_stats()
        return {
            "bases": {label: int(bases[i]) for i, label in enumerate(labels)},
            "blocks": {label: int(blocks[i]) for i, label in enumerate(labels)},
        }

    return _per_read("read_cigar_stats", alignment_file, contig, start, end, max_records, extract)


def read_reference_overlap(
    alignment_file: str | Path,
    contig: str,
    start: int,
    end: int,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count how many aligned bases of each read overlap ``[start, end)`` on the reference
    (``AlignedSegment.get_overlap``)."""

    chrom, first, last = _optional_region(contig, start, end)
    if first is None or last is None:
        raise ConfigurationError("start and end are required.")
    return _per_read(
        "read_reference_overlap",
        alignment_file,
        chrom,
        first,
        last,
        max_records,
        lambda read: {"overlap_bases": int(read.get_overlap(first, last) or 0)},
    )


def read_reference_positions(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    full_length: bool = False,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the 0-based reference positions covered by aligned read bases (``AlignedSegment.get_reference_positions``).

    ``full_length=True`` returns one entry per read base with ``None`` for insertions and clips.
    """

    full = _flag(full_length, "full_length")
    return _per_read(
        "read_reference_positions",
        alignment_file,
        contig,
        start,
        end,
        max_records,
        lambda read: {"positions": read.get_reference_positions(full_length=full)},
        full_length=full,
    )


def reconstruct_read_reference_sequence(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rebuild the reference sequence under each read from its MD tag (``AlignedSegment.get_reference_sequence``).

    Reference bases at mismatches are lower case, matches upper case, as pysam returns them. Reads without an MD tag are
    listed in ``skipped_without_md`` instead of failing.
    """

    chrom, first, last = _optional_region(contig, start, end)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    _pysam, handle = _open(alignment_file)
    rows: list[dict[str, Any]] = []
    skipped: list[str] = []
    with handle:
        reads, truncated = _limited(_reads(handle, chrom, first, last), limit)
        for read in reads:
            if read.is_unmapped or not read.has_tag("MD"):
                skipped.append(read.query_name)
                continue
            try:
                rows.append(
                    {
                        "query_name": read.query_name,
                        "reference_sequence": read.get_reference_sequence(),
                    }
                )
            except _ERRORS as exc:
                raise BackendExecutionError(f"pysam get_reference_sequence failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "reconstruct_read_reference_sequence",
        "region": {"contig": chrom, "start": first, "end": last},
        "count": len(rows),
        "truncated": truncated,
        "records": rows,
        "skipped_without_md": skipped,
        "provenance": _provenance("AlignedSegment.get_reference_sequence"),
    }


def infer_read_lengths(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer each read's length from its CIGAR (``infer_query_length`` and ``infer_read_length``).

    ``query_length`` counts M, I, S, =, X (without hard clips); ``read_length`` also counts hard-clipped
    bases (the original read length). ``stored_length`` is the length of the stored sequence (None if
    the sequence is ``*``).
    """

    def extract(read: Any) -> dict[str, Any]:
        return {
            "query_length": read.infer_query_length(),
            "read_length": read.infer_read_length(),
            "stored_length": read.query_length if read.query_sequence else None,
        }

    return _per_read("infer_read_lengths", alignment_file, contig, start, end, max_records, extract)


def read_base_modifications(
    alignment_file: str | Path,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    *,
    forward: bool = False,
    max_records: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode the base modifications of each read from its MM/ML tags (``modified_bases``).

    For every modification (original base, strand, code such as ``m``) the result lists
    ``[read position, probability 0-255]`` pairs. ``forward=True`` uses ``modified_bases_forward`` (positions
    in the original sequencing direction). Reads without MM/ML tags have an empty mapping.
    """

    use_forward = _flag(forward, "forward")

    def extract(read: Any) -> dict[str, Any]:
        raw = read.modified_bases_forward if use_forward else read.modified_bases
        return {
            "modifications": [
                {
                    "base": key[0],
                    "strand": int(key[1]),
                    "modification": key[2] if isinstance(key[2], str) else int(key[2]),
                    "calls": [[int(p), int(q)] for p, q in calls],
                }
                for key, calls in raw.items()
            ]
        }

    return _per_read(
        "read_base_modifications",
        alignment_file,
        contig,
        start,
        end,
        max_records,
        extract,
        forward=use_forward,
    )


__all__ = [
    "alignment_records_to_dicts",
    "count_base_coverage",
    "count_spliced_junctions",
    "fetch_mate_alignments",
    "find_alignments_by_name",
    "infer_read_lengths",
    "pileup_columns",
    "read_aligned_blocks",
    "read_aligned_pairs",
    "read_base_modifications",
    "read_cigar_stats",
    "read_reference_overlap",
    "read_reference_positions",
    "reconstruct_read_reference_sequence",
    "sam_header_to_dict",
    "write_alignments_from_dicts",
    "write_sam_header",
]
_register_public_verbose(globals())
