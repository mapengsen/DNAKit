"""Local candidate adapters for the epigenome capability family."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.epigenome import _impl_local as _impl
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.epigenome.external import extract_methylation_methyldackel

def align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Align epigenomic reads to a reference (not implemented: needs an aligner such as Bismark, Bowtie2 or BWA)."""
    return candidate_not_implemented(
        "epigenome", "align", input, options=options,
        requires="a read aligner (Bismark/bwa-meth for bisulfite, Bowtie2/BWA for ChIP/ATAC) and an indexed reference",
        message="Read alignment is not implemented natively; no aligner is wired into this function.",
    )

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate genomic regions (e.g. peaks) with overlapping and nearest features from a supplied annotation.

    ``input``: regions as BED text/path, rows ``[chrom, start, end, ...]`` or mappings with chrom/start/end
    (0-based half-open). ``options``: ``features`` (required: the annotation as BED with a name column, rows or
    mappings with ``name``); ``max_distance`` (optional: nearest features farther than this are not reported).
    Per region: ``overlapping`` feature names and the ``nearest`` feature (``name``, ``distance`` in bases, 0 when
    overlapping or touching, and ``side``: ``overlap``, ``upstream`` or ``downstream`` of the region on the forward
    strand). Features on other chromosomes are ignored. Quadratic search; at most 50 million region x feature pairs.
    """
    opts = dict(options or {})
    _impl.entity(opts, "peaks", ("peaks",))
    regions = _impl.parse_intervals(input)
    if "features" not in opts:
        raise ConfigurationError("annotate needs options['features'].", code="EPIGENOME_OPTION_INVALID",
                                 hint="Pass a BED with names: {'features': 'chr1\\t100\\t200\\tgeneA'}.")
    features = _impl.parse_intervals(opts["features"])
    if len(regions) * len(features) > 50_000_000:
        raise ConfigurationError("Too many region x feature pairs.", code="EPIGENOME_INPUT_TOO_LARGE")
    max_distance = _impl.opt_int(opts, "max_distance", -1, -1, 10**12)
    by_chrom = _impl.group_by_chrom(features)
    rows = []
    for r in regions:
        rs, re_ = int(r["start"]), int(r["end"])  # type: ignore[call-overload]
        overlapping, best = [], None
        for f in by_chrom.get(str(r["chrom"]), []):
            fs, fe = int(f["start"]), int(f["end"])  # type: ignore[call-overload]
            gap = _impl.interval_gap(rs, re_, fs, fe)
            label = str(f.get("name", f"{f['chrom']}:{fs}-{fe}"))
            if fs < re_ and fe > rs:
                overlapping.append(label)
            if best is None or gap < best[0]:
                best = (gap, label, "overlap" if fs < re_ and fe > rs else "downstream" if fs >= re_ else "upstream")
        nearest = None
        if best is not None and (max_distance < 0 or best[0] <= max_distance):
            nearest = {"name": best[1], "distance": best[0], "side": best[2]}
        rows.append({"chrom": r["chrom"], "start": rs, "end": re_, "overlapping": overlapping, "nearest": nearest})
    return candidate_ok("epigenome", "annotate", input, {"regions": rows}, options=opts)

def build_chromatin_state_segments(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Collapse per-bin chromatin-state calls into contiguous state segments.

    ``input``: bins as BED text/path, rows or mappings ``chrom, start, end, name`` where ``name`` is the state label
    (``state`` is accepted as an alias key in mappings), or ``{chrom: [state, state, ...]}`` with
    ``options["bin_size"]`` (bins start at 0 and are consecutive). ``options``: ``max_gap`` (default 0): consecutive
    bins of the same state on one chromosome are joined when the gap between them is <= max_gap bases.
    Result: ``segments`` (``chrom``, ``start``, ``end``, ``state``, ``n_bins``) sorted by position; overlapping bins
    of the same state are also joined. This segments given calls; it does not learn states (no ChromHMM).
    """
    opts = dict(options or {})
    max_gap = _impl.opt_int(opts, "max_gap", 0, 0, 10**9)
    if isinstance(input, Mapping) and input and all(isinstance(v, Sequence) and not isinstance(v, (str, bytes)) for v in input.values()) \
            and "records" not in input:
        size = _impl.opt_int(opts, "bin_size", 0, 1, 10**9) if "bin_size" in opts else 0
        if not size:
            raise ConfigurationError("A {chrom: [states]} input needs options['bin_size'].", code="EPIGENOME_OPTION_INVALID")
        bins = [{"chrom": str(c), "start": i * size, "end": (i + 1) * size, "name": str(s)}
                for c, states in input.items() for i, s in enumerate(states)]
    else:
        raw = []
        if isinstance(input, Sequence) and not isinstance(input, (str, bytes)):
            for row in input:
                raw.append({**row, "name": row.get("name", row.get("state"))} if isinstance(row, Mapping) else row)
        else:
            raw = input  # type: ignore[assignment]
        bins = _impl.parse_intervals(raw)
        if any("name" not in b or b["name"] is None for b in bins):
            raise InputFormatError("Every bin needs a state label (BED name column).", code="EPIGENOME_RECORD_INVALID")
    segments: list[dict[str, object]] = []
    for b in sorted(bins, key=lambda r: (str(r["chrom"]), int(r["start"]), int(r["end"]))):  # type: ignore[call-overload]
        last = segments[-1] if segments else None
        if last and last["chrom"] == b["chrom"] and last["state"] == b["name"] and int(b["start"]) - int(last["end"]) <= max_gap:  # type: ignore[call-overload]
            last["end"] = max(int(last["end"]), int(b["end"]))  # type: ignore[call-overload]
            last["n_bins"] = int(last["n_bins"]) + 1  # type: ignore[call-overload]
        else:
            segments.append({"chrom": b["chrom"], "start": b["start"], "end": b["end"], "state": b["name"], "n_bins": 1})
    return candidate_ok("epigenome", "build_chromatin_state_segments", input,
                        {"bin_count": len(bins), "segment_count": len(segments), "segments": segments}, options=opts)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate methylation levels from counts, CpG observed/expected ratios from DNA, or interval overlap.

    ``options["entity"]``:
    ``"methylation_level"`` (default): ``input`` is per-site counts (Bismark ``.cov`` text/path/rows
    ``chrom, start, end, percent, methylated, unmethylated``, or mappings with ``chrom``, ``start``, ``methylated``,
    ``unmethylated`` (or ``coverage``)). Level = methylated / (methylated + unmethylated). ``min_coverage``
    (default 1) drops sites with fewer reads; ``one_based`` (default false; set true for ``.cov``) shifts
    ``start`` to 0-based; ``regions`` (optional BED) aggregates sites inside each region (0-based half-open):
    ``weighted_level`` = sum(methylated) / sum(coverage) and ``mean_site_level`` = unweighted mean of site levels;
    without ``regions`` per-site levels are returned.
    ``"cpg_observed_expected"``: ``input`` is DNA (string, FASTA text/path or ``{id: seq}``). Reports C, G and CpG
    counts, GC fraction and observed/expected CpG = CpG * L / (C * G) (null if C*G = 0); with ``window`` (and
    ``step``, default = window) also per window.
    ``"interval_overlap"``: ``input`` and ``options["b"]`` are interval sets; after merging each set, reports
    ``overlap_bp``, ``union_bp`` and Jaccard = overlap/union. ``a_intervals_overlapping_b`` counts the original
    (unmerged) input intervals that share at least one base with any interval of ``b`` (half-open, so book-ended
    intervals do not touch); each input interval is counted once however many ``b`` intervals it spans.
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "methylation_level", ("methylation_level", "cpg_observed_expected", "interval_overlap"))
    if which == "methylation_level":
        sites = _impl.parse_methylation(input)
        min_cov = _impl.opt_int(opts, "min_coverage", 1, 1, 10**9)
        shift = 1 if opts.get("one_based", False) is True else 0
        kept = []
        for s in sites:
            cov = int(s["methylated"]) + int(s["unmethylated"])  # type: ignore[call-overload]
            if cov >= min_cov:
                kept.append({"chrom": s["chrom"], "pos": int(s["start"]) - shift, "coverage": cov,  # type: ignore[call-overload]
                             "methylated": s["methylated"], "level": int(s["methylated"]) / cov})  # type: ignore[call-overload]
        if "regions" in opts:
            rows = []
            for reg in _impl.parse_intervals(opts["regions"]):
                inside = [k for k in kept if k["chrom"] == reg["chrom"] and int(reg["start"]) <= int(k["pos"]) < int(reg["end"])]  # type: ignore[call-overload]
                cov = sum(int(k["coverage"]) for k in inside)  # type: ignore[call-overload]
                meth = sum(int(k["methylated"]) for k in inside)  # type: ignore[call-overload]
                rows.append({"chrom": reg["chrom"], "start": reg["start"], "end": reg["end"], "n_sites": len(inside),
                             "coverage": cov, "weighted_level": meth / cov if cov else None,
                             "mean_site_level": sum(float(k["level"]) for k in inside) / len(inside) if inside else None})  # type: ignore[arg-type]
            result: dict[str, object] = {"entity": which, "regions": rows}
        else:
            result = {"entity": which, "sites": kept, "site_count": len(kept), "dropped_low_coverage": len(sites) - len(kept)}
    elif which == "cpg_observed_expected":
        def stats(seq: str) -> dict[str, object]:
            c, g, cg = seq.count("C"), seq.count("G"), seq.count("CG")
            return {"length": len(seq), "c": c, "g": g, "cpg": cg, "gc_fraction": (c + g) / len(seq),
                    "observed_expected": cg * len(seq) / (c * g) if c and g else None}

        window = _impl.opt_int(opts, "window", 0, 0, 10**7) if "window" in opts else 0
        step = _impl.opt_int(opts, "step", window or 1, 1, 10**7)
        out = []
        for name, seq in _dna_records(input):
            entry: dict[str, object] = {"id": name, **stats(seq)}
            if window:
                entry["windows"] = [{"start": i, "end": i + window, **stats(seq[i : i + window])}
                                    for i in range(0, len(seq) - window + 1, step)]
            out.append(entry)
        result = {"entity": which, "sequences": out}
    else:
        a = _impl.merge_intervals(_impl.parse_intervals(input))
        if "b" not in opts:
            raise ConfigurationError("interval_overlap needs options['b'].", code="EPIGENOME_OPTION_INVALID")
        b = _impl.merge_intervals(_impl.parse_intervals(opts["b"]))
        a_by, b_by = _impl.group_by_chrom(a), _impl.group_by_chrom(b)
        inter = sum(_impl.overlap_bp((int(x["start"]), int(x["end"])), (int(y["start"]), int(y["end"])))  # type: ignore[call-overload]
                    for chrom, xs in a_by.items() for x in xs for y in b_by.get(chrom, []))
        a_bp = sum(int(x["end"]) - int(x["start"]) for x in a)  # type: ignore[call-overload]
        b_bp = sum(int(x["end"]) - int(x["start"]) for x in b)  # type: ignore[call-overload]
        touching = sum(1 for x in _impl.parse_intervals(input)
                       if any(int(x["start"]) < int(y["end"]) and int(y["start"]) < int(x["end"])  # type: ignore[call-overload]
                              for y in b_by.get(str(x["chrom"]), [])))
        result = {"entity": which, "overlap_bp": inter, "a_bp": a_bp, "b_bp": b_bp, "union_bp": a_bp + b_bp - inter,
                  "jaccard": inter / (a_bp + b_bp - inter) if a_bp + b_bp - inter else 0.0, "a_intervals_overlapping_b": touching}
    return candidate_ok("epigenome", "calculate", input, result, options=opts)

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Call consensus peaks across samples: regions covered by at least ``min_samples`` samples' peak sets.

    ``input``: peak sets as ``{sample: intervals}`` (each BED text/path, rows or mappings) or a list of interval
    sets. ``options``: ``entity`` (only ``"consensus_peaks"``); ``min_samples`` (default 2). Each sample's own
    overlapping peaks are merged first so a sample votes at most once per base; the result is every maximal segment
    whose vote count is >= ``min_samples`` (0-based half-open). This combines existing peak calls; it does not call
    peaks from signal (no MACS).
    """
    opts = dict(options or {})
    _impl.entity(opts, "consensus_peaks", ("consensus_peaks",))
    if isinstance(input, Mapping) and not {"records", "regions"} & set(input):
        sets = {str(k): _impl.parse_intervals(v) for k, v in input.items()}
    elif isinstance(input, Sequence) and not isinstance(input, (str, bytes)) and input and all(
        isinstance(s, (str, bytes, Path)) or (isinstance(s, Sequence) and s and not isinstance(s[0], (str, int, float))) for s in input
    ):
        sets = {f"sample{i + 1}": _impl.parse_intervals(v) for i, v in enumerate(input)}
    else:
        raise InputFormatError("Provide peak sets as {sample: intervals} or a list of interval sets.",
                               code="EPIGENOME_INPUT_INVALID")
    min_samples = _impl.opt_int(opts, "min_samples", 2, 1, len(sets))
    peaks = _impl.depth_segments(list(sets.values()), min_samples)
    for p in peaks:
        p.pop("n_merged", None)
    result = {"samples": list(sets), "min_samples": min_samples, "peak_count": len(peaks), "peaks": peaks}
    return candidate_ok("epigenome", "call", input, result, options=opts)

def find(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find differentially methylated regions by Fisher's exact test, or CpG islands in DNA.

    ``options["entity"]``:
    ``"differential_regions"`` (default): ``input`` is rows/mappings with ``chrom``, ``start``, ``end`` and pooled
    counts ``methylated_a``, ``unmethylated_a``, ``methylated_b``, ``unmethylated_b`` (group A vs B). If
    ``options["regions"]`` (BED) is given, ``input`` rows are single sites (``chrom``, ``start`` and the four counts)
    that are summed inside each region (0-based half-open). Per region: group levels, ``delta = level_b - level_a``,
    two-sided Fisher exact ``p_value`` (the scipy definition: sum of table probabilities <= the observed one), sample
    ``odds_ratio``, Benjamini-Hochberg ``q_value`` over the tested regions. ``min_coverage`` (default 10) per group
    excludes regions; ``alpha`` (default 0.05) and ``min_delta`` (default 0.0) set the ``significant`` flag.
    ``"cpg_islands"``: ``input`` is one DNA sequence; delegates to ``dnakit.patterns.find_cpg_islands`` with
    ``window_size`` 200, ``step`` 1, ``min_gc`` 0.5, ``min_observed_expected`` 0.6, ``min_region_length`` 200.
    """
    opts = dict(options or {})
    which = _impl.entity(opts, "differential_regions", ("differential_regions", "cpg_islands"))
    if which == "cpg_islands":
        from dnakit.core import DNASequence
        from dnakit.patterns import find_cpg_islands

        name, seq = _dna_records(input)[0]
        found = find_cpg_islands(
            DNASequence(seq), window_size=_impl.opt_int(opts, "window_size", 200, 2, 10**6), step=_impl.opt_int(opts, "step", 1, 1, 10**6),
            min_gc=_impl.opt_float(opts, "min_gc", 0.5, 0, 1), min_observed_expected=_impl.opt_float(opts, "min_observed_expected", 0.6, 0, 100),
            min_region_length=_impl.opt_int(opts, "min_region_length", 200, 1, 10**7), verbose=False,
        )
        islands = [{"start": h.coordinate_location.start, "end": h.coordinate_location.end, "length": h.length,
                    "observed_expected": h.score, "gc_fraction": dict(h.attributes)["gc_fraction"]} for h in found.hits]
        return candidate_ok("epigenome", "find", input, {"entity": which, "sequence": name, "islands": islands}, options=opts)
    keys = ("methylated_a", "unmethylated_a", "methylated_b", "unmethylated_b")
    rows_in = []
    for where, row in _impl._rows(input):
        if not isinstance(row, Mapping) or not all(k in row for k in keys) or "chrom" not in row or "start" not in row:
            raise InputFormatError(f"Differential-region rows need chrom, start, end and {', '.join(keys)} ({where}).",
                                   code="EPIGENOME_RECORD_INVALID")
        counts = [_impl._int(row[k], k, where) for k in keys]
        if min(counts) < 0:
            raise InputFormatError(f"Counts must be non-negative ({where}).", code="EPIGENOME_RECORD_INVALID")
        rows_in.append({"chrom": str(row["chrom"]), "start": _impl._int(row["start"], "start", where),
                        "end": _impl._int(row.get("end", int(row["start"]) + 1), "end", where), "counts": counts})
    if "regions" in opts:
        regions = _impl.parse_intervals(opts["regions"])
        merged = []
        for reg in regions:
            inside = [r for r in rows_in if r["chrom"] == reg["chrom"] and int(reg["start"]) <= r["start"] < int(reg["end"])]  # type: ignore[call-overload]
            merged.append({"chrom": reg["chrom"], "start": reg["start"], "end": reg["end"],
                           "counts": [sum(r["counts"][i] for r in inside) for i in range(4)]})  # type: ignore[index]
        rows_in = merged  # type: ignore[assignment]
    min_cov = _impl.opt_int(opts, "min_coverage", 10, 1, 10**9)
    alpha = _impl.opt_float(opts, "alpha", 0.05, 0, 1)
    min_delta = _impl.opt_float(opts, "min_delta", 0.0, 0, 1)
    tested = []
    for r in rows_in:
        ma, ua, mb, ub = r["counts"]  # type: ignore[misc]
        if ma + ua < min_cov or mb + ub < min_cov:
            continue
        test = _impl.fisher_exact(ma, ua, mb, ub)
        level_a, level_b = ma / (ma + ua), mb / (mb + ub)
        tested.append({"chrom": r["chrom"], "start": r["start"], "end": r["end"], "level_a": level_a, "level_b": level_b,
                       "delta": level_b - level_a, "coverage_a": ma + ua, "coverage_b": mb + ub,
                       "p_value": test["p_value"], "odds_ratio": test["odds_ratio"]})
    for row, q in zip(tested, _impl.benjamini_hochberg([t["p_value"] for t in tested])):  # type: ignore[misc]
        row["q_value"] = q
        row["significant"] = q <= alpha and abs(float(row["delta"])) >= min_delta  # type: ignore[arg-type]
    tested.sort(key=lambda t: (t["p_value"], t["chrom"], t["start"]))  # type: ignore[arg-type]
    result = {"entity": which, "tested": len(tested), "excluded_low_coverage": len(rows_in) - len(tested), "regions": tested}
    return candidate_ok("epigenome", "find", input, result, options=opts)

def link_enhancers_to_genes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign enhancers (or any regions) to genes by genomic distance to the gene transcription start site (TSS).

    ``input``: enhancers as BED text/path, rows or mappings (0-based half-open). ``options``: ``genes`` (required;
    mappings/rows with ``chrom``, ``start``, ``end``, ``name`` and optional ``strand`` or explicit ``tss``; TSS =
    ``start`` for + or unstranded genes, ``end - 1`` for - genes); ``method``: ``"nearest"`` (default: the single
    nearest TSS on the same chromosome) or ``"window"`` (every gene with a TSS within ``max_distance``);
    ``max_distance`` (default 1,000,000 bp; for nearest, farther genes are not linked). Distance is measured from
    the enhancer midpoint ``(start + end) // 2``; ties go to the alphabetically first gene. This is a proximity
    heuristic, not evidence of regulation (no Hi-C or correlation support is used).
    """
    opts = dict(options or {})
    if "genes" not in opts:
        raise ConfigurationError("link_enhancers_to_genes needs options['genes'].", code="EPIGENOME_OPTION_INVALID")
    method = _impl.opt_str(opts, "method", "nearest", ("nearest", "window"))
    max_distance = _impl.opt_int(opts, "max_distance", 1_000_000, 0, 10**12)
    regions = _impl.parse_intervals(input)
    genes = _impl.parse_genes(opts["genes"])
    links = _impl.link_regions_to_genes(regions, genes, method, max_distance)
    return candidate_ok("epigenome", "link_enhancers_to_genes", input,
                        {"method": method, "max_distance": max_distance, "links": links}, options=opts,
                        warnings=["Distance-based assignment is a proximity heuristic, not evidence of regulation."])

def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge overlapping or nearby intervals (peaks, regions) per chromosome, like ``bedtools merge``.

    ``input``: BED text/path, rows ``[chrom, start, end, ...]`` or mappings (0-based half-open). ``options``:
    ``entity`` (only ``"peaks"``); ``distance`` (default 0): intervals whose gap is <= distance are merged, so
    with 0 overlapping and book-ended intervals merge. Strand and names are ignored. Result: ``intervals`` sorted by
    chromosome name and start, each with ``n_merged`` (how many input intervals it combines).
    """
    opts = dict(options or {})
    _impl.entity(opts, "peaks", ("peaks",))
    distance = _impl.opt_int(opts, "distance", 0, 0, 10**12)
    merged = _impl.merge_intervals(_impl.parse_intervals(input), distance)
    return candidate_ok("epigenome", "merge", input, {"distance": distance, "intervals": merged,
                                                      "count": len(merged)}, options=opts)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise signal values (bedGraph records or a plain list of numbers).

    ``input``: bedGraph text/path, rows ``[chrom, start, end, value]``, mappings with ``value``, or a list of numbers.
    ``options``: ``entity`` (only ``"signal"``); ``method``: ``"cpm"`` (default; value / sum(values) * 1e6),
    ``"rpkm"`` (cpm / (interval length / 1000); needs intervals), ``"log2cpm"`` (log2(cpm + ``pseudocount``,
    default 1)), ``"zscore"`` ((v - mean) / sd; population sd, ``ddof`` option default 0) or ``"minmax"``
    ((v - min) / (max - min)). Values must be non-negative for the cpm family. Result: the records with ``value``
    replaced by the normalised value (``raw_value`` keeps the input) and the parameters used.
    """
    opts = dict(options or {})
    _impl.entity(opts, "signal", ("signal",))
    method = _impl.opt_str(opts, "method", "cpm", ("cpm", "rpkm", "log2cpm", "zscore", "minmax"))
    if isinstance(input, Sequence) and not isinstance(input, (str, bytes)) and input and all(
        isinstance(v, (int, float)) and not isinstance(v, bool) for v in input
    ):
        records: list[dict[str, object]] = [{"value": float(v)} for v in input]
        if method == "rpkm":
            raise ConfigurationError("rpkm needs interval lengths; pass bedGraph records.", code="EPIGENOME_OPTION_INVALID")
    else:
        records = _impl.parse_intervals(input, "bedgraph")
    values = [float(r["value"]) for r in records]  # type: ignore[arg-type]
    if not all(math.isfinite(v) for v in values):
        raise InputFormatError("Signal values must be finite.", code="EPIGENOME_RECORD_INVALID")
    if method in ("cpm", "rpkm", "log2cpm"):
        if min(values) < 0:
            raise InputFormatError("Signal values must be non-negative for this method.", code="EPIGENOME_RECORD_INVALID")
        total = sum(values)
        if total <= 0:
            raise InputFormatError("The signal sums to zero.", code="EPIGENOME_RECORD_INVALID")
        cpm = [v / total * 1e6 for v in values]
        if method == "cpm":
            new = cpm
        elif method == "rpkm":
            new = [c / ((int(r["end"]) - int(r["start"])) / 1000.0) for c, r in zip(cpm, records)]  # type: ignore[call-overload]
        else:
            pc = _impl.opt_float(opts, "pseudocount", 1.0, 1e-12, 1e12)
            new = [math.log2(c + pc) for c in cpm]
    elif method == "zscore":
        ddof = _impl.opt_int(opts, "ddof", 0, 0, 1)
        mean = sum(values) / len(values)
        var = sum((v - mean) ** 2 for v in values) / (len(values) - ddof) if len(values) > ddof else 0.0
        if var == 0:
            raise InputFormatError("Zero variance: z-scores are undefined.", code="EPIGENOME_RECORD_INVALID")
        new = [(v - mean) / math.sqrt(var) for v in values]
    else:
        lo, hi = min(values), max(values)
        if hi == lo:
            raise InputFormatError("All values are equal: min-max scaling is undefined.", code="EPIGENOME_RECORD_INVALID")
        new = [(v - lo) / (hi - lo) for v in values]
    out = [{**r, "raw_value": r["value"], "value": n} for r, n in zip(records, new)]
    return candidate_ok("epigenome", "normalize", input, {"method": method, "records": out}, options=opts)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a signal track (bedGraph bins as bars) for one genomic window as an SVG string.

    ``input``: bedGraph text/path, rows ``[chrom, start, end, value]`` or mappings. ``options``: ``entity`` (only
    ``"tracks"``); ``chrom``, ``start``, ``end`` (default: the first chromosome and the span of its bins);
    ``output`` (optional ``.svg`` path). Bars are ``<rect class="bin">`` elements with ``data-start``, ``data-end``
    and ``data-value`` attributes, heights scaled to the maximum value in the window. Result: ``{"svg": ...}``.
    """
    opts = dict(options or {})
    _impl.entity(opts, "tracks", ("tracks",))
    records = _impl.parse_intervals(input, "bedgraph")
    chrom = str(opts.get("chrom", records[0]["chrom"]))
    on_chrom = [r for r in records if r["chrom"] == chrom]
    if not on_chrom:
        raise InputFormatError(f"No records on {chrom!r}.", code="EPIGENOME_INPUT_EMPTY")
    start = _impl.opt_int(opts, "start", min(int(r["start"]) for r in on_chrom), 0, 10**12)  # type: ignore[call-overload]
    end = _impl.opt_int(opts, "end", max(int(r["end"]) for r in on_chrom), 1, 10**12)  # type: ignore[call-overload]
    if end <= start:
        raise ConfigurationError("end must be greater than start.", code="EPIGENOME_OPTION_INVALID")
    svg = _impl.svg_signal_track(records, chrom, start, end)
    written = None
    if opts.get("output") is not None:
        Path(str(opts["output"])).write_text(svg, encoding="utf-8")
        written = str(opts["output"])
    return candidate_ok("epigenome", "plot", input, {"svg": svg, "output_path": written}, options=opts)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict epigenomic signal from sequence (not implemented: needs a trained model such as Enformer or Basenji)."""
    return candidate_not_implemented(
        "epigenome", "predict", input, options=options,
        requires="a trained sequence-to-epigenome model and its weights (Enformer, Basenji, DeepSEA)",
        message="Epigenome prediction needs trained model weights that are not bundled.",
    )

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quality-check bisulfite methylation counts: coverage distribution, global methylation and conversion rate.

    ``input``: per-site counts (Bismark ``.cov`` text/path/rows or mappings with chrom, start, methylated,
    unmethylated; see ``calculate``). ``options``: ``entity`` (only ``"methylation"``); ``min_coverage`` (default 10);
    ``control_chroms`` (optional list of unmethylated control sequences such as lambda phage or chrL).
    Result: ``sites``, total/mean/median ``coverage``, ``fraction_sites_min_coverage``, ``global_methylation`` =
    sum(methylated) / sum(coverage) excluding control chromosomes, and, with controls,
    ``conversion_rate = 1 - methylated/coverage`` on them (the share of unmethylated cytosines read as T).
    """
    opts = dict(options or {})
    _impl.entity(opts, "methylation", ("methylation",))
    sites = _impl.parse_methylation(input)
    min_cov = _impl.opt_int(opts, "min_coverage", 10, 1, 10**9)
    controls = opts.get("control_chroms", [])
    if isinstance(controls, str) or not isinstance(controls, Sequence):
        raise ConfigurationError("control_chroms must be a list of chromosome names.", code="EPIGENOME_OPTION_INVALID")
    controls = {str(c) for c in controls}
    sample = [s for s in sites if s["chrom"] not in controls]
    control = [s for s in sites if s["chrom"] in controls]
    covs = sorted(int(s["methylated"]) + int(s["unmethylated"]) for s in sample)  # type: ignore[call-overload]
    if not covs:
        raise InputFormatError("No sites remain after removing control chromosomes.", code="EPIGENOME_INPUT_EMPTY")
    mid = len(covs) // 2
    median = float(covs[mid]) if len(covs) % 2 else (covs[mid - 1] + covs[mid]) / 2.0
    total_cov = sum(covs)
    meth = sum(int(s["methylated"]) for s in sample)  # type: ignore[call-overload]
    result: dict[str, object] = {
        "sites": len(sample), "total_coverage": total_cov, "mean_coverage": total_cov / len(covs), "median_coverage": median,
        "min_coverage": min_cov, "fraction_sites_min_coverage": sum(1 for c in covs if c >= min_cov) / len(covs),
        "global_methylation": meth / total_cov if total_cov else None,
    }
    if controls:
        ccov = sum(int(s["methylated"]) + int(s["unmethylated"]) for s in control)  # type: ignore[call-overload]
        cmeth = sum(int(s["methylated"]) for s in control)  # type: ignore[call-overload]
        result["control_coverage"] = ccov
        result["conversion_rate"] = 1 - cmeth / ccov if ccov else None
    return candidate_ok("epigenome", "qc", input, result, options=opts)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read BED, bedGraph or Bismark coverage (per-site methylation counts) text, a file path or in-memory rows.

    ``input``: text, a path, rows of cells or mappings. ``options``: ``format``: ``"bed"`` (default; chrom, start,
    end, optional name/score/strand), ``"bedgraph"`` (chrom, start, end, value) or ``"bismark_cov"`` (chrom,
    start, end, percent, methylated, unmethylated). ``track``/``browser``/``#`` lines and blanks are skipped;
    malformed lines raise an error with the line number. Result: ``format``, ``count`` and ``records``
    (0-based half-open for bed/bedgraph; ``.cov`` positions are returned as given, 1-based).
    """
    opts = dict(options or {})
    fmt = _impl.opt_str(opts, "format", "bed", ("bed", "bedgraph", "bismark_cov"))
    records = _impl.parse_methylation(input) if fmt == "bismark_cov" else _impl.parse_intervals(input, fmt)
    return candidate_ok("epigenome", "read", input, {"format": fmt, "count": len(records), "records": records}, options=opts)

def _dna_records(value: object) -> list[tuple[str, str]]:
    """DNA string, FASTA text/path or ``{id: sequence}`` as cleaned (id, ACGT sequence) pairs."""

    if isinstance(value, Mapping):
        items = [(str(k), str(v)) for k, v in value.items()]
    else:
        text = _impl.read_text(value)
        if text.lstrip().startswith(">"):
            items, name, chunks = [], None, []
            for line in text.splitlines():
                if line.startswith(">"):
                    if name is not None:
                        items.append((name, "".join(chunks)))
                    name, chunks = line[1:].split()[0] if line[1:].strip() else "seq", []
                else:
                    chunks.append(line.strip())
            if name is not None:
                items.append((name, "".join(chunks)))
        else:
            items = [("seq1", "".join(text.split()))]
    out = []
    for name, seq in items:
        seq = seq.upper().replace("U", "T")
        if not seq or set(seq) - set("ACGT"):
            raise InputFormatError(f"Sequence {name!r} is empty or has non-ACGT characters.", code="EPIGENOME_SEQUENCE_INVALID")
        out.append((name, seq))
    if not out or len(out) > 10_000 or sum(len(s) for _, s in out) > 50_000_000:
        raise ConfigurationError("Provide 1-10000 sequences of at most 50 Mb in total.", code="EPIGENOME_INPUT_TOO_LARGE")
    return out


__all__ = ['align', 'annotate', 'build_chromatin_state_segments', 'calculate', 'call', 'extract_methylation_methyldackel', 'find', 'link_enhancers_to_genes', 'merge', 'normalize', 'plot', 'predict', 'qc', 'read']

__candidate_adapter_operations__ = set(['align', 'annotate', 'build_chromatin_state_segments', 'calculate', 'call', 'find', 'link_enhancers_to_genes', 'merge', 'normalize', 'plot', 'predict', 'qc', 'read'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    balance_contact_matrix,
    build_contact_matrix,
    build_extended_tag_pileup,
    build_fragment_pileup,
    calculate_coverage,
    calculate_expected_contacts,
    calculate_insulation_score,
    calculate_non_conversion_rate,
    calculate_site_binomial_pvalues,
    call_fragment_peaks,
    count_bam_reads_per_bin,
    deduplicate_pairs,
    estimate_fragment_lengths,
    fetch_contact_matrix,
    filter_allc_sites,
    read_pairs_header,
    summarize_bigwig_bins,
    summarize_pairs,
)

__all__ = [*__all__, "balance_contact_matrix", "build_contact_matrix", "build_extended_tag_pileup", "build_fragment_pileup", "calculate_coverage", "calculate_expected_contacts", "calculate_insulation_score", "calculate_non_conversion_rate", "calculate_site_binomial_pvalues", "call_fragment_peaks", "count_bam_reads_per_bin", "deduplicate_pairs", "estimate_fragment_lengths", "fetch_contact_matrix", "filter_allc_sites", "read_pairs_header", "summarize_bigwig_bins", "summarize_pairs"]

# Full-pass exports (see docs/api/features)
from .MACS3_functions import (
    build_paired_fragment_pileup,
    build_shifted_read_pileup,
    calculate_allele_model_likelihood,
    calculate_genotype_quality,
    filter_summits_by_peakiness,
    filter_summits_by_valleys,
    find_signal_maxima,
    guess_tag_file_format,
    hmmratac_fragment_weights,
    max_duplicate_tags,
    merge_pileup_arrays,
    refine_peak_summit,
)

from .cooler_functions import (
    coarsen_contact_matrix,
    fetch_chrom_sizes,
    load_fasta_lengths,
    merge_contact_matrices,
    parse_genomic_region,
    query_bed_intervals,
    read_chrom_sizes,
    read_cooler_chroms,
    read_cooler_pixels,
    rename_cooler_chromosomes,
    zoomify_contact_matrix,
)

from .cooltools_functions import (
    align_windows_to_bins,
    assign_features_to_view,
    calculate_compartments,
    calculate_directionality_index,
    calculate_expected_trans,
    calculate_logbin_expected,
    calculate_saddle,
    calculate_trans_compartments,
    calculate_virtual_4c,
    dense_cis_eigenvectors,
    digitize_track,
    expected_from_dense_matrix,
    iterative_correction_dense,
    make_chromosome_view,
    observed_over_expected_matrix,
    pileup_features,
    read_expected_table,
    read_viewframe_bed,
    rearrange_contact_matrix,
    sample_contacts,
)

from .deeptools_functions import (
    bedgraph_to_bigwig,
    bigwig_chromosome_sizes,
    chop_regions_for_profile,
    estimate_bam_scale_factors,
    estimate_count_size_factors,
    fingerprint_jensen_shannon,
    gc_fraction_from_twobit,
    map_twobit_to_bam_chromosomes,
    match_chromosome_name,
    parse_deeptools_region,
)

from .methylpy_functions import (
    benjamini_hochberg_methylation_cutoffs,
    convert_fastq_to_bisulfite,
    count_bisulfite_basecalls,
    decode_read_c_positions,
    encode_read_c_positions,
    expand_methylation_context,
    filter_bam_by_mch_level,
    index_allc_table,
    read_allc_index,
    read_fasta_chromosome,
)

from .pairtools_functions import (
    assign_genomic_regions,
    assign_haplotype_phase,
    compute_contact_scaling,
    find_read_mismatches,
    find_restriction_fragment,
    make_pairs_header,
    normalize_scaling_curve,
    parse_sam_alignments,
    read_chrom_order,
    read_pairs_table,
    subset_pairs_header,
)

__all__ = [*__all__, "build_paired_fragment_pileup", "build_shifted_read_pileup", "calculate_allele_model_likelihood", "calculate_genotype_quality", "filter_summits_by_peakiness", "filter_summits_by_valleys", "find_signal_maxima", "guess_tag_file_format", "hmmratac_fragment_weights", "max_duplicate_tags", "merge_pileup_arrays", "refine_peak_summit", "coarsen_contact_matrix", "fetch_chrom_sizes", "load_fasta_lengths", "merge_contact_matrices", "parse_genomic_region", "query_bed_intervals", "read_chrom_sizes", "read_cooler_chroms", "read_cooler_pixels", "rename_cooler_chromosomes", "zoomify_contact_matrix", "align_windows_to_bins", "assign_features_to_view", "calculate_compartments", "calculate_directionality_index", "calculate_expected_trans", "calculate_logbin_expected", "calculate_saddle", "calculate_trans_compartments", "calculate_virtual_4c", "dense_cis_eigenvectors", "digitize_track", "expected_from_dense_matrix", "iterative_correction_dense", "make_chromosome_view", "observed_over_expected_matrix", "pileup_features", "read_expected_table", "read_viewframe_bed", "rearrange_contact_matrix", "sample_contacts", "bedgraph_to_bigwig", "bigwig_chromosome_sizes", "chop_regions_for_profile", "estimate_bam_scale_factors", "estimate_count_size_factors", "fingerprint_jensen_shannon", "gc_fraction_from_twobit", "map_twobit_to_bam_chromosomes", "match_chromosome_name", "parse_deeptools_region", "benjamini_hochberg_methylation_cutoffs", "convert_fastq_to_bisulfite", "count_bisulfite_basecalls", "decode_read_c_positions", "encode_read_c_positions", "expand_methylation_context", "filter_bam_by_mch_level", "index_allc_table", "read_allc_index", "read_fasta_chromosome", "assign_genomic_regions", "assign_haplotype_phase", "compute_contact_scaling", "find_read_mismatches", "find_restriction_fragment", "make_pairs_header", "normalize_scaling_curve", "parse_sam_alignments", "read_chrom_order", "read_pairs_table", "subset_pairs_header"]

# Review-pass exports (see docs/api/features)
from .MACS3_review_functions import (
    digest_hmmratac_fragments,
    extract_hmmratac_training_data,
    hmmratac_proba_to_bedgraph,
    hmmratac_states_path,
    load_hmmratac_model,
    predict_hmmratac_states,
    save_hmmratac_model,
    train_hmmratac_model,
    write_hmmratac_accessible_regions,
    write_hmmratac_states_bed,
)

from .cooler_review_functions import (
    create_single_cell_cooler,
    infer_cooler_binsize,
    list_cooler_groups,
    list_scool_cells,
)

from .cooltools_review_functions import (
    adaptive_coarsegrain_matrix,
    align_track_to_cooler_bins,
    assign_bins_to_regions,
    assign_paired_features_to_view,
    call_chromatin_dots,
    dense_trans_eigenvectors,
    expand_features_to_bins,
    expected_full_table,
    make_track_view,
    mask_track_bad_bins,
    observed_over_expected_pixels,
    pairwise_block_sums,
    pairwise_diagonal_sums,
    replace_cis_with_decoys,
    saddle_compartment_strength,
    smooth_contact_decay,
)

from .deeptools_review_functions import (
    bigwig_bedgraph_tiles,
    count_reads_per_gc_bins,
    fingerprint_synthetic_jsd,
    tabulate_gc_fragment_counts,
)

from .methylpy_review_functions import (
    add_cytosine_info_to_dmrs,
    add_methylation_levels_to_dmrs,
    collapse_dmr_windows,
    dmr_residual_cutoff,
    estimate_dmr_fdr_cutoff,
    filter_allc_by_pvalue_cutoffs,
    merge_allc_tables,
    merge_dms_to_dmrs,
    run_allc_binomial_test,
    split_allc_by_chromosome,
    split_allc_into_chunks,
    split_mpileup_chunks,
)

from .pairtools_review_functions import (
    chromsizes_from_sam_header,
    filter_pairs_by_coverage,
    mark_sam_duplicate,
    sam_to_pairs,
)

__all__ = [*__all__, "digest_hmmratac_fragments", "extract_hmmratac_training_data", "hmmratac_proba_to_bedgraph", "hmmratac_states_path", "load_hmmratac_model", "predict_hmmratac_states", "save_hmmratac_model", "train_hmmratac_model", "write_hmmratac_accessible_regions", "write_hmmratac_states_bed", "create_single_cell_cooler", "infer_cooler_binsize", "list_cooler_groups", "list_scool_cells", "adaptive_coarsegrain_matrix", "align_track_to_cooler_bins", "assign_bins_to_regions", "assign_paired_features_to_view", "call_chromatin_dots", "dense_trans_eigenvectors", "expand_features_to_bins", "expected_full_table", "make_track_view", "mask_track_bad_bins", "observed_over_expected_pixels", "pairwise_block_sums", "pairwise_diagonal_sums", "replace_cis_with_decoys", "saddle_compartment_strength", "smooth_contact_decay", "bigwig_bedgraph_tiles", "count_reads_per_gc_bins", "fingerprint_synthetic_jsd", "tabulate_gc_fragment_counts", "add_cytosine_info_to_dmrs", "add_methylation_levels_to_dmrs", "collapse_dmr_windows", "dmr_residual_cutoff", "estimate_dmr_fdr_cutoff", "filter_allc_by_pvalue_cutoffs", "merge_allc_tables", "merge_dms_to_dmrs", "run_allc_binomial_test", "split_allc_by_chromosome", "split_allc_into_chunks", "split_mpileup_chunks", "chromsizes_from_sam_header", "filter_pairs_by_coverage", "mark_sam_duplicate", "sam_to_pairs"]

# Family-pass exports (see docs/api/features)
from .MACS3_family_functions import (
    build_counted_fragment_pileup,
    build_fragment_end_pileup,
    build_hmmratac_fragment_pileup,
    calculate_allele_specific_genotype_quality,
    calculate_allele_specific_model_likelihood,
)

from .cooltools_family_functions import (
    compute_region_contact_scaling,
    pairwise_region_sums,
)

from .deeptools_family_functions import (
    chop_regions_from_middle,
)

from .methylpy_family_functions import (
    convert_paired_fastq_to_bisulfite,
)

from .pairtools_family_functions import (
    deduplicate_pairs_online,
    normalize_scaling_curve_unity,
)

__all__ = [*__all__, "build_counted_fragment_pileup", "build_fragment_end_pileup", "build_hmmratac_fragment_pileup", "calculate_allele_specific_genotype_quality", "calculate_allele_specific_model_likelihood", "compute_region_contact_scaling", "pairwise_region_sums", "chop_regions_from_middle", "convert_paired_fastq_to_bisulfite", "deduplicate_pairs_online", "normalize_scaling_curve_unity"]

_register_public_verbose(globals())
