"""genomicranges interval operations as DNAKit functions (pack P40).

Input is an explicit local BED file (``chrom start end [name score strand]``, 0-based half-open, as everywhere in
DNAKit). The upstream package stores ranges 1-based and inclusive, so the functions shift the start by one on the
way in and back on the way out; returned coordinates are BED again. Results are JSON-serialisable dictionaries
with bounded record lists. genomicranges is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-genomicranges"
_STRAND = {"+": 1, "-": -1, ".": 0, "*": 0}
_STRAND_BACK = {1: "+", -1: "-", 0: "*"}


def _read_bed(path: str | Path, name: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with _h.open_text(_h.existing(path, name)) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t") if "\t" in line else line.split()
            if len(fields) < 3:
                raise ConfigurationError(f"{name} line {number}: a BED line needs chrom, start and end.")
            try:
                start, end = int(fields[2 - 1]), int(fields[3 - 1])
            except ValueError as exc:
                raise ConfigurationError(f"{name} line {number}: start and end must be integers.") from exc
            if start < 0 or end <= start:
                raise ConfigurationError(f"{name} line {number}: need 0 <= start < end.")
            strand = fields[5] if len(fields) > 5 else "*"
            if strand not in _STRAND:
                raise ConfigurationError(f"{name} line {number}: strand must be +, - or '.'.")
            score = None
            if len(fields) > 4:
                try:
                    score = float(fields[4])
                except ValueError:
                    score = None
            rows.append({"chrom": fields[0], "start": start, "end": end, "strand": strand, "score": score,
                         "name": fields[3] if len(fields) > 3 else None})
            if len(rows) > 1_000_000:
                raise ConfigurationError(f"{name} has more than 1,000,000 intervals.")
    if not rows:
        raise ConfigurationError(f"{name} contains no intervals.")
    return rows


def _granges(rows: list[dict[str, Any]], *, seqlengths: Mapping[str, int] | None = None,
             score: bool = False) -> Any:
    gr = require("genomicranges", _EXTRA)
    pd = require("pandas", _EXTRA)
    data: dict[str, Any] = {
        "seqnames": [r["chrom"] for r in rows],
        "starts": [r["start"] + 1 for r in rows],
        "ends": [r["end"] for r in rows],
        "strand": [r["strand"] if r["strand"] != "." else "*" for r in rows],
    }
    if score:
        if any(r["score"] is None for r in rows):
            raise ConfigurationError("Every interval needs a numeric score (BED column 5).")
        data["score"] = [r["score"] for r in rows]
    ranges = gr.GenomicRanges.from_pandas(pd.DataFrame(data))
    if seqlengths is not None:
        names = sorted({r["chrom"] for r in rows} | set(seqlengths))
        lengths = [int(seqlengths[n]) if n in seqlengths else None for n in names]
        ranges = ranges.set_seqinfo(gr.SeqInfo(names, lengths))
    return ranges


def _records(ranges: Any, limit: int) -> list[dict[str, Any]]:
    out = []
    names = list(ranges.get_seqnames())
    starts, ends, strands = ranges.get_start(), ranges.get_end(), ranges.get_strand()
    for i in range(min(len(names), limit)):
        out.append({"chrom": str(names[i]), "start": int(starts[i]) - 1, "end": int(ends[i]),
                    "strand": _STRAND_BACK.get(int(strands[i]), "*")})
    return out


def _version() -> str:
    from importlib.metadata import version

    return version("genomicranges")


def _limit(value: int) -> int:
    return _h.integer(value, "max_rows", minimum=1, maximum=_h.MAX_ROWS * 10)


def _result(operation: str, **fields: Any) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, "backend": "genomicranges",
            "version": _version(), "coordinates": "BED, 0-based half-open", **fields}


def genomicranges_restrict_intervals(
    source: str | Path,
    *,
    start: int | None = None,
    end: int | None = None,
    keep_all_ranges: bool = False,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Restrict intervals to a coordinate window (genomicranges GenomicRanges.restrict).

    ``start`` and ``end`` are BED coordinates of the window (0-based start, exclusive end). Intervals
    outside the window are dropped unless ``keep_all_ranges`` is True (they are then kept as empty).
    """

    limit = _limit(max_rows)
    ranges = _granges(_read_bed(source, "source"))
    lo = None if start is None else _h.integer(start, "start", minimum=0) + 1
    hi = None if end is None else _h.integer(end, "end", minimum=1)
    result = ranges.restrict(start=lo, end=hi, keep_all_ranges=bool(keep_all_ranges))
    records = _records(result, limit)
    return _result("restrict", interval_count=len(result), intervals=records, truncated=len(result) > limit)


def genomicranges_narrow_intervals(
    source: str | Path,
    *,
    start: int | None = None,
    width: int | None = None,
    end: int | None = None,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Narrow every interval by relative ``start``, ``width`` or ``end`` (GenomicRanges.narrow).

    As in the R/Bioconductor original, the numbers are 1-based offsets inside each interval (``start=2``
    drops the first base); give two of the three.
    """

    limit = _limit(max_rows)
    ranges = _granges(_read_bed(source, "source"))
    kwargs = {k: _h.integer(v, k, minimum=1) for k, v in (("start", start), ("width", width), ("end", end))
              if v is not None}
    if len(kwargs) != 2:
        raise ConfigurationError("Give exactly two of start, width and end.")
    result = ranges.narrow(**kwargs)
    return _result("narrow", interval_count=len(result), intervals=_records(result, limit),
                   truncated=len(result) > limit)


def genomicranges_out_of_bound_intervals(
    source: str | Path,
    seqlengths: Mapping[str, int],
    *,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find intervals that extend past the chromosome length (GenomicRanges.get_out_of_bound_index).

    ``seqlengths`` maps chromosome name to length; intervals on chromosomes without a length are never
    reported.
    """

    limit = _limit(max_rows)
    if not isinstance(seqlengths, Mapping) or not seqlengths:
        raise ConfigurationError("seqlengths must be a non-empty mapping of chromosome to length.")
    lengths = {str(k): _h.integer(v, f"seqlengths[{k}]", minimum=1) for k, v in seqlengths.items()}
    rows = _read_bed(source, "source")
    ranges = _granges(rows, seqlengths=lengths)
    indices = [int(i) for i in ranges.get_out_of_bound_index()]
    return _result("get_out_of_bound_index", out_of_bound_count=len(indices), indices=indices[:limit],
                   intervals=[_records(ranges, len(rows))[i] for i in indices[:limit]])


def genomicranges_span_intervals(
    source: str | Path,
    *,
    ignore_strand: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Span (minimum start to maximum end) of the intervals per chromosome and strand (GenomicRanges.range)."""

    ranges = _granges(_read_bed(source, "source"))
    result = ranges.range(ignore_strand=bool(ignore_strand))
    return _result("range", ignore_strand=bool(ignore_strand), interval_count=len(result),
                   intervals=_records(result, _h.MAX_ROWS))


def genomicranges_is_disjoint(
    source: str | Path,
    *,
    ignore_strand: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether no two intervals overlap on the same chromosome (and strand) (GenomicRanges.is_disjoint)."""

    rows = _read_bed(source, "source")
    flag = bool(_granges(rows).is_disjoint(ignore_strand=bool(ignore_strand)))
    return _result("is_disjoint", ignore_strand=bool(ignore_strand), interval_count=len(rows), is_disjoint=flag)


def genomicranges_disjoint_bins(
    source: str | Path,
    *,
    ignore_strand: bool = False,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign each interval to a bin so that intervals within one bin do not overlap (GenomicRanges.disjoint_bins)."""

    limit = _limit(max_rows)
    rows = _read_bed(source, "source")
    bins = [int(b) for b in _granges(rows).disjoint_bins(ignore_strand=bool(ignore_strand))]
    return _result("disjoint_bins", ignore_strand=bool(ignore_strand), bin_count=len(set(bins)),
                   bins=bins[:limit],
                   intervals=[{"chrom": r["chrom"], "start": r["start"], "end": r["end"], "bin": b}
                              for r, b in list(zip(rows, bins))[:limit]],
                   truncated=len(rows) > limit)


def genomicranges_interval_coverage(
    source: str | Path,
    *,
    shift: int = 0,
    width: int | None = None,
    weight: int = 1,
    ignore_strand: bool = True,
    max_runs: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-position coverage depth of the intervals on each chromosome (GenomicRanges.coverage).

    The depth is returned as run-length encoded BED segments ``[start, end) depth`` (zero-depth runs
    included), starting from position 0 of each chromosome.
    """

    limit = _limit(max_runs)
    ranges = _granges(_read_bed(source, "source"))
    kwargs: dict[str, Any] = {"shift": _h.integer(shift, "shift"), "weight": _h.integer(weight, "weight"),
                              "ignore_strand": bool(ignore_strand)}
    if width is not None:
        kwargs["width"] = _h.integer(width, "width", minimum=1)
    coverage = ranges.coverage(**kwargs)
    runs: list[dict[str, Any]] = []
    for chrom, values in coverage.items():
        begin = 0
        for i in range(len(values)):
            if i + 1 == len(values) or values[i + 1] != values[i]:
                if len(runs) < limit:
                    runs.append({"chrom": str(chrom), "start": begin, "end": i + 1, "depth": float(values[i])})
                begin = i + 1
    return _result("coverage", chromosome_count=len(coverage), runs=runs,
                   lengths={str(k): int(len(v)) for k, v in coverage.items()},
                   max_depth=max((float(v.max()) for v in coverage.values()), default=0.0))


def genomicranges_match_intervals(
    source: str | Path,
    query: str | Path,
    *,
    ignore_strand: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find, for every query interval, the identical interval in ``source`` (GenomicRanges.match).

    Returns for each query line the index of the exact match (same chromosome, start, end and, unless
    ``ignore_strand``, strand) or None.
    """

    subject = _granges(_read_bed(source, "source"))
    queries = _read_bed(query, "query")
    matches = subject.match(_granges(queries), ignore_strand=bool(ignore_strand))
    indices = [None if (x is None or int(x) < 0) else int(x) for x in matches]
    return _result("match", query_count=len(queries), matched_count=sum(i is not None for i in indices),
                   indices=indices)


def genomicranges_order_intervals(
    source: str | Path,
    *,
    decreasing: bool = False,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Indices that sort the intervals by chromosome, strand and position (GenomicRanges.order)."""

    limit = _limit(max_rows)
    ranges = _granges(_read_bed(source, "source"))
    order = [int(i) for i in ranges.order(decreasing=bool(decreasing))]
    return _result("order", decreasing=bool(decreasing), order=order[:limit], truncated=len(order) > limit)


def genomicranges_rank_intervals(
    source: str | Path,
    *,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rank of each interval in coordinate order (GenomicRanges.rank)."""

    limit = _limit(max_rows)
    ranks = [int(i) for i in _granges(_read_bed(source, "source")).rank()]
    return _result("rank", ranks=ranks[:limit], truncated=len(ranks) > limit)


def genomicranges_tile_intervals(
    source: str | Path,
    *,
    n: int | None = None,
    width: int | None = None,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tile every interval into ``n`` equal sub-intervals or into pieces of ``width`` (GenomicRanges.tile)."""

    limit = _limit(max_rows)
    if (n is None) == (width is None):
        raise ConfigurationError("Give exactly one of n and width.")
    kwargs = {"n": _h.integer(n, "n", minimum=1)} if n is not None else {
        "width": _h.integer(width, "width", minimum=1)}
    tiles = _granges(_read_bed(source, "source")).tile(**kwargs)
    pieces = [_records(t, _h.MAX_ROWS) for t in tiles][:limit]
    return _result("tile", interval_count=len(tiles), tiles=pieces, truncated=len(tiles) > limit)


def genomicranges_binned_average(
    source: str | Path,
    bins: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average the BED score (column 5) of ``source`` over each interval of ``bins`` (GenomicRanges.binned_average).

    ``bins`` is a BED file of bins; the result lists the mean score per bin (None when a bin overlaps no
    interval), computed by genomicranges (partial overlaps are weighted by overlap width).
    """

    ranges = _granges(_read_bed(source, "source"), score=True)
    bin_rows = _read_bed(bins, "bins")
    result = ranges.binned_average("score", _granges(bin_rows), outname="binned_average")
    values = result.get_mcols().get_column("binned_average")
    averages = [None if v is None else _h.plain(v) for v in values]
    return _result("binned_average", bin_count=len(bin_rows),
                   bins=[{"chrom": b["chrom"], "start": b["start"], "end": b["end"], "mean_score": a}
                         for b, a in zip(bin_rows, averages)])


__all__ = [
    "genomicranges_binned_average",
    "genomicranges_disjoint_bins",
    "genomicranges_interval_coverage",
    "genomicranges_is_disjoint",
    "genomicranges_match_intervals",
    "genomicranges_narrow_intervals",
    "genomicranges_order_intervals",
    "genomicranges_out_of_bound_intervals",
    "genomicranges_rank_intervals",
    "genomicranges_restrict_intervals",
    "genomicranges_span_intervals",
    "genomicranges_tile_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
