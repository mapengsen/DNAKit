"""Native helpers behind the epigenome candidate functions (BED/bedGraph/methylation counts, intervals, Fisher test)."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_RECORDS = 2_000_000
MAX_TEXT_BYTES = 200_000_000


# ----------------------------------------------------------------------------- option helpers
def opt_str(
    options: Mapping[str, object], key: str, default: str, allowed: Sequence[str] | None = None
) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"Option {key!r} must be a non-empty string.", code="EPIGENOME_OPTION_INVALID"
        )
    value = value.strip().lower()
    if allowed is not None and value not in allowed:
        raise ConfigurationError(
            f"Option {key!r}={value!r} is not supported.",
            code="EPIGENOME_OPTION_INVALID",
            hint=f"Use one of: {', '.join(allowed)}.",
        )
    return value


def opt_int(
    options: Mapping[str, object], key: str, default: int, minimum: int, maximum: int
) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r} must be an integer in [{minimum}, {maximum}].",
            code="EPIGENOME_OPTION_INVALID",
        )
    return value


def opt_float(
    options: Mapping[str, object], key: str, default: float, minimum: float, maximum: float
) -> float:
    value = options.get(key, default)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not minimum <= value <= maximum
    ):
        raise ConfigurationError(
            f"Option {key!r} must be a number in [{minimum}, {maximum}].",
            code="EPIGENOME_OPTION_INVALID",
        )
    return float(value)


def entity(options: Mapping[str, object], default: str, allowed: Sequence[str]) -> str:
    return opt_str(options, "entity", default, allowed)


# ----------------------------------------------------------------------------- text / record parsing
def _is_file(text: str) -> bool:
    if "\n" in text or "\t" in text or len(text) >= 4096:
        return False
    try:
        return Path(text).is_file()
    except OSError:  # e.g. file name too long: it is text, not a path
        return False


def read_text(value: object) -> str:
    if isinstance(value, Path) or (isinstance(value, str) and _is_file(value)):
        path = Path(value)
        if path.stat().st_size > MAX_TEXT_BYTES:
            raise ConfigurationError("Input file is too large.", code="EPIGENOME_INPUT_TOO_LARGE")
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    raise InputFormatError("Expected text, bytes or a path.", code="EPIGENOME_INPUT_INVALID")


def _split_lines(text: str) -> list[tuple[int, list[str]]]:
    rows = []
    for number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "track", "browser")):
            continue
        rows.append((number, line.rstrip("\r\n").split("\t") if "\t" in line else stripped.split()))
    return rows


def _int(value: object, what: str, where: object) -> int:
    if isinstance(value, bool):
        raise InputFormatError(
            f"{what} must be an integer ({where}).", code="EPIGENOME_RECORD_INVALID"
        )
    try:
        number = int(str(value)) if not isinstance(value, float) else int(value)
    except ValueError as exc:
        raise InputFormatError(
            f"{what} {value!r} is not an integer ({where}).", code="EPIGENOME_RECORD_INVALID"
        ) from exc
    if isinstance(value, float) and value != number:
        raise InputFormatError(
            f"{what} {value!r} is not an integer ({where}).", code="EPIGENOME_RECORD_INVALID"
        )
    return number


def _float(value: object, what: str, where: object) -> float:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise InputFormatError(
            f"{what} {value!r} is not a number ({where}).", code="EPIGENOME_RECORD_INVALID"
        ) from exc
    if not math.isfinite(number):
        raise InputFormatError(f"{what} must be finite ({where}).", code="EPIGENOME_RECORD_INVALID")
    return number


def _rows(value: object) -> list[tuple[object, object]]:
    """Return (location, row) pairs where a row is a list of cells or a mapping."""

    if isinstance(value, (str, bytes, Path)):
        return [(f"line {n}", cells) for n, cells in _split_lines(read_text(value))]
    if isinstance(value, Mapping):
        value = value.get("records", value.get("regions", [value]))
    if isinstance(value, Sequence):
        return [(f"record {i + 1}", row) for i, row in enumerate(value)]
    raise InputFormatError("Unsupported interval input.", code="EPIGENOME_INPUT_INVALID")


def _get(row: Mapping[str, object], *names: str) -> object:
    for name in names:
        if name in row:
            return row[name]
    return None


def parse_intervals(value: object, kind: str = "bed") -> list[dict[str, object]]:
    """BED (``kind="bed"``: chrom, start, end, name, score, strand) or bedGraph (``kind="bedgraph"``: + value).

    Accepts text/path, rows of cells, or mappings with chrom/start/end (and name/value/strand/score).
    Coordinates are 0-based half-open; end must be greater than start.
    """

    out: list[dict[str, object]] = []
    for where, row in _rows(value):
        if isinstance(row, Mapping):
            chrom = _get(row, "chrom", "chr", "seqid", "contig")
            start, end = _get(row, "start", "begin"), _get(row, "end", "stop")
            extra = dict(row)
            name, strand, score = _get(row, "name", "id"), _get(row, "strand"), _get(row, "score")
            val = _get(row, "value", "signal")
            record_extra = {
                k: v
                for k, v in extra.items()
                if k not in ("chrom", "chr", "seqid", "contig", "start", "begin", "end", "stop")
            }
        elif isinstance(row, Sequence) and not isinstance(row, (str, bytes)):
            cells = list(row)
            if len(cells) < (4 if kind == "bedgraph" else 3):
                raise InputFormatError(
                    f"Too few columns ({where}).", code="EPIGENOME_RECORD_INVALID"
                )
            chrom, start, end = cells[0], cells[1], cells[2]
            name = cells[3] if kind == "bed" and len(cells) > 3 else None
            score = cells[4] if kind == "bed" and len(cells) > 4 else None
            strand = cells[5] if kind == "bed" and len(cells) > 5 else None
            val = cells[3] if kind == "bedgraph" else None
            record_extra = {}
        else:
            raise InputFormatError(
                f"Unsupported record ({where}).", code="EPIGENOME_RECORD_INVALID"
            )
        if chrom is None or start is None or end is None:
            raise InputFormatError(
                f"Record needs chrom, start and end ({where}).", code="EPIGENOME_RECORD_INVALID"
            )
        s, e = _int(start, "start", where), _int(end, "end", where)
        if s < 0 or e <= s:
            raise InputFormatError(
                f"Invalid interval {chrom}:{s}-{e} ({where}).",
                code="EPIGENOME_RECORD_INVALID",
                hint="Intervals are 0-based half-open with end > start >= 0.",
            )
        record: dict[str, object] = {"chrom": str(chrom), "start": s, "end": e}
        if name is not None:
            record["name"] = str(name)
        if score is not None and str(score) != ".":
            record["score"] = _float(score, "score", where)
        if strand is not None:
            if str(strand) not in ("+", "-", "."):
                raise InputFormatError(
                    f"Strand must be +, - or . ({where}).", code="EPIGENOME_RECORD_INVALID"
                )
            record["strand"] = str(strand)
        if kind == "bedgraph":
            if val is None:
                raise InputFormatError(
                    f"bedGraph record needs a value ({where}).", code="EPIGENOME_RECORD_INVALID"
                )
            record["value"] = _float(val, "value", where)
        elif val is not None:
            record["value"] = _float(val, "value", where)
        for k, v in record_extra.items():
            record.setdefault(k, v)
        out.append(record)
        if len(out) > MAX_RECORDS:
            raise ConfigurationError("Too many records.", code="EPIGENOME_INPUT_TOO_LARGE")
    if not out:
        raise InputFormatError("No intervals found.", code="EPIGENOME_INPUT_EMPTY")
    return out


_METH = ("methylated", "meth", "count_methylated", "n_methylated", "mc", "c")
_UNMETH = ("unmethylated", "unmeth", "count_unmethylated", "n_unmethylated", "t")


def parse_methylation(value: object) -> list[dict[str, object]]:
    """Per-site methylation counts: Bismark ``.cov`` rows (chrom, start, end, percent, methylated, unmethylated),
    or mappings with chrom, start (end optional) and methylated plus unmethylated (or coverage) counts."""

    out: list[dict[str, object]] = []
    for where, row in _rows(value):
        if isinstance(row, Mapping):
            chrom, start = _get(row, "chrom", "chr", "seqid"), _get(row, "start", "pos", "position")
            m, u = _get(row, *_METH), _get(row, *_UNMETH)
            if m is not None and u is None and _get(row, "coverage", "cov", "depth") is not None:
                u = _int(_get(row, "coverage", "cov", "depth"), "coverage", where) - _int(
                    m, "methylated", where
                )
        else:
            cells = list(row)  # type: ignore[arg-type]
            if len(cells) < 6:
                raise InputFormatError(
                    f"Bismark coverage rows need 6 columns ({where}).",
                    code="EPIGENOME_RECORD_INVALID",
                )
            chrom, start, m, u = cells[0], cells[1], cells[4], cells[5]
        if chrom is None or start is None or m is None or u is None:
            raise InputFormatError(
                f"Methylation record needs chrom, start, methylated and unmethylated ({where}).",
                code="EPIGENOME_RECORD_INVALID",
            )
        mi, ui = _int(m, "methylated", where), _int(u, "unmethylated", where)
        if mi < 0 or ui < 0:
            raise InputFormatError(
                f"Counts must be non-negative ({where}).", code="EPIGENOME_RECORD_INVALID"
            )
        out.append(
            {
                "chrom": str(chrom),
                "start": _int(start, "start", where),
                "methylated": mi,
                "unmethylated": ui,
            }
        )
        if len(out) > MAX_RECORDS:
            raise ConfigurationError("Too many records.", code="EPIGENOME_INPUT_TOO_LARGE")
    if not out:
        raise InputFormatError("No methylation records found.", code="EPIGENOME_INPUT_EMPTY")
    return out


# ----------------------------------------------------------------------------- intervals
def merge_intervals(
    intervals: Sequence[Mapping[str, object]], distance: int = 0
) -> list[dict[str, object]]:
    """bedtools-style merge: intervals on one chromosome whose gap is <= ``distance`` (book-ended included) are joined."""

    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for iv in intervals:
        by_chrom.setdefault(str(iv["chrom"]), []).append((int(iv["start"]), int(iv["end"])))  # type: ignore[call-overload]
    out: list[dict[str, object]] = []
    for chrom in sorted(by_chrom):
        current: list[int] | None = None
        count = 0
        for s, e in sorted(by_chrom[chrom]):
            if current is not None and s - current[1] <= distance:
                current[1] = max(current[1], e)
                count += 1
            else:
                if current is not None:
                    out.append(
                        {"chrom": chrom, "start": current[0], "end": current[1], "n_merged": count}
                    )
                current, count = [s, e], 1
        if current is not None:
            out.append({"chrom": chrom, "start": current[0], "end": current[1], "n_merged": count})
    return out


def depth_segments(
    per_sample: Sequence[Sequence[Mapping[str, object]]], min_depth: int
) -> list[dict[str, object]]:
    """Maximal segments where at least ``min_depth`` samples have coverage (each sample counted once per base)."""

    events: dict[str, list[tuple[int, int]]] = {}
    for intervals in per_sample:
        for merged in merge_intervals(intervals, 0):  # one vote per sample and base
            events.setdefault(str(merged["chrom"]), []).extend(
                [(int(merged["start"]), 1), (int(merged["end"]), -1)]
            )  # type: ignore[call-overload]
    out: list[dict[str, object]] = []
    for chrom in sorted(events):
        depth = 0
        start = None
        points = sorted(events[chrom])
        i = 0
        while i < len(points):
            pos = points[i][0]
            while i < len(points) and points[i][0] == pos:
                depth += points[i][1]
                i += 1
            if depth >= min_depth and start is None:
                start = pos
            elif depth < min_depth and start is not None:
                out.append({"chrom": chrom, "start": start, "end": pos})
                start = None
    # join segments touching end-to-start (can happen when depth stays >= min_depth across a boundary)
    return [{**m, "n_merged": None} for m in _join_touching(out)]


def _join_touching(segments: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for seg in segments:
        if out and out[-1]["chrom"] == seg["chrom"] and out[-1]["end"] == seg["start"]:
            out[-1] = {**out[-1], "end": seg["end"]}
        else:
            out.append(dict(seg))
    return out


def overlap_bp(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def group_by_chrom(
    intervals: Sequence[Mapping[str, object]],
) -> dict[str, list[Mapping[str, object]]]:
    out: dict[str, list[Mapping[str, object]]] = {}
    for iv in intervals:
        out.setdefault(str(iv["chrom"]), []).append(iv)
    for chrom in out:
        out[chrom].sort(key=lambda r: (int(r["start"]), int(r["end"])))  # type: ignore[call-overload]
    return out


def interval_gap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """0 when the half-open intervals overlap or touch, else the number of bases between them."""

    if b_start >= a_end:
        return b_start - a_end
    if a_start >= b_end:
        return a_start - b_end
    return 0


# ----------------------------------------------------------------------------- statistics
def _log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact(a: int, b: int, c: int, d: int) -> dict[str, float | None]:
    """Two-sided Fisher exact test on [[a, b], [c, d]]: sum of hypergeometric probabilities <= P(observed) (1 + 1e-7).

    Same definition as ``scipy.stats.fisher_exact(..., alternative="two-sided")``; the odds ratio is the sample
    odds ratio ad/bc (None when bc = 0).
    """

    n1, n2, k = a + b, c + d, a + c
    n = n1 + n2
    if n == 0:
        return {"p_value": 1.0, "odds_ratio": None}
    low, high = max(0, k - n2), min(n1, k)
    if high - low > 5_000_000:
        raise ConfigurationError(
            "Counts are too large for an exact test.", code="EPIGENOME_COUNTS_TOO_LARGE"
        )
    denom = _log_choose(n, k)

    def log_pmf(x: int) -> float:
        return _log_choose(n1, x) + _log_choose(n2, k - x) - denom

    observed = log_pmf(a)
    cutoff = observed + math.log1p(1e-7)
    total = 0.0
    for x in range(low, high + 1):
        lp = log_pmf(x)
        if lp <= cutoff:
            total += math.exp(lp)
    odds = (a * d) / (b * c) if b * c else None
    return {"p_value": min(total, 1.0), "odds_ratio": odds}


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, pvalues[index] * n / rank)
        adjusted[index] = running
    return adjusted


# ----------------------------------------------------------------------------- enhancer-gene links
def parse_genes(value: object) -> list[dict[str, object]]:
    """Genes as BED-like rows (chrom, start, end, name, score, strand) or mappings with chrom/start/end/name/strand/tss."""

    genes = parse_intervals(value, "bed")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for gene, raw in zip(genes, value, strict=False):
            if isinstance(raw, Mapping) and raw.get("tss") is not None:
                gene["tss"] = _int(raw["tss"], "tss", gene.get("name", ""))
    for gene in genes:
        gene.setdefault("name", f"{gene['chrom']}:{gene['start']}-{gene['end']}")
    return genes


def gene_tss(gene: Mapping[str, object]) -> int:
    """TSS of a gene record: explicit ``tss``; else ``start`` for + / unstranded and ``end - 1`` for - strand."""

    if gene.get("tss") is not None:
        return int(gene["tss"])  # type: ignore[call-overload]
    return int(gene["end"]) - 1 if gene.get("strand") == "-" else int(gene["start"])  # type: ignore[call-overload]


def _direction(mid: int, tss: int, strand: object) -> str:
    """Position of the region relative to the gene's orientation: upstream, downstream or at_tss."""

    if mid == tss:
        return "at_tss"
    return "upstream" if (mid < tss) == (strand != "-") else "downstream"


def link_regions_to_genes(
    regions: Sequence[Mapping[str, object]],
    genes: Sequence[Mapping[str, object]],
    method: str,
    max_distance: int,
) -> list[dict[str, object]]:
    """Assign regions to genes by distance from the region midpoint ``(start + end) // 2`` to the gene TSS.

    ``nearest``: each region gets its single nearest same-chromosome TSS (ties: lexicographically smallest gene
    name) if within ``max_distance``; ``window``: every gene whose TSS lies within ``max_distance``.
    """

    by_chrom: dict[str, list[tuple[int, str, Mapping[str, object]]]] = {}
    for g in genes:
        by_chrom.setdefault(str(g["chrom"]), []).append((gene_tss(g), str(g.get("name", "")), g))
    out: list[dict[str, object]] = []
    for r in regions:
        mid = (int(r["start"]) + int(r["end"])) // 2  # type: ignore[call-overload]
        candidates = sorted(
            ((abs(tss - mid), name, tss, g) for tss, name, g in by_chrom.get(str(r["chrom"]), [])),
            key=lambda t: (t[0], t[1]),
        )
        chosen = (
            candidates[:1]
            if method == "nearest"
            else [c for c in candidates if c[0] <= max_distance]
        )
        links = [c for c in chosen if c[0] <= max_distance]
        out.append(
            {
                "region": {
                    "chrom": r["chrom"],
                    "start": r["start"],
                    "end": r["end"],
                    **({"name": r["name"]} if "name" in r else {}),
                },
                "midpoint": mid,
                "genes": [
                    {
                        "gene": name,
                        "tss": tss,
                        "distance": dist,
                        "direction": _direction(mid, tss, g.get("strand")),
                        "strand": g.get("strand", "."),
                    }
                    for dist, name, tss, g in links
                ],
            }
        )
    return out


# ----------------------------------------------------------------------------- SVG
def svg_signal_track(
    records: Sequence[Mapping[str, object]],
    chrom: str,
    start: int,
    end: int,
    width: int = 800,
    height: int = 120,
) -> str:
    rows = [
        r
        for r in records
        if r["chrom"] == chrom and int(r["end"]) > start and int(r["start"]) < end
    ]  # type: ignore[call-overload]
    top = max([float(r["value"]) for r in rows] + [0.0]) or 1.0  # type: ignore[arg-type]
    scale = (width - 40) / (end - start)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height + 30}">',
        f'<text x="20" y="12" font-size="10">{chrom}:{start}-{end} (max {top:g})</text>',
    ]
    for r in rows:
        s, e = max(int(r["start"]), start), min(int(r["end"]), end)  # type: ignore[call-overload]
        h = float(r["value"]) / top * height  # type: ignore[arg-type]
        parts.append(
            f'<rect class="bin" data-start="{r["start"]}" data-end="{r["end"]}" data-value="{r["value"]}" '
            f'x="{20 + (s - start) * scale:.2f}" y="{20 + height - h:.2f}" width="{max((e - s) * scale, 0.5):.2f}" '
            f'height="{h:.2f}" fill="#4e79a7"/>'
        )
    parts.append(
        f'<line x1="20" y1="{20 + height}" x2="{width - 20}" y2="{20 + height}" stroke="#333"/></svg>'
    )
    return "".join(parts)
