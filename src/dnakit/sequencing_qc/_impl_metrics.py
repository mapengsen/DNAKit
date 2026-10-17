"""Native read-level QC metrics behind the ``dnakit.sequencing_qc`` candidate functions.

Every ``op_*`` function takes ``(value, options, warnings)`` and returns a JSON-serialisable
result. Definitions follow textbook / FastQC / Picard conventions and are stated in the public
docstrings. Qualities are Phred scores decoded with ``phred_offset`` (default 33).
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

from ._impl_reads import (
    MAX_PHRED,
    Rd,
    as_text,
    gc_fraction,
    hamming,
    hist_percentile,
    kth,
    lx,
    mean,
    need_quality,
    nx,
    opt_bool,
    opt_choice,
    opt_int,
    opt_num,
    p_to_q,
    parse_fasta_text,
    parse_pair,
    parse_reads,
    percentile,
    q_to_p,
    revcomp,
    round6,
    stdev,
)

Warnings = list[str]

TENX_CHEMISTRY = {"v1": (14, 10), "v2": (16, 10), "v3": (16, 12)}

DEFAULT_ADAPTERS = {
    "Illumina Universal Adapter": "AGATCGGAAGAG",
    "Illumina Small RNA 3' Adapter": "TGGAATTCTCGG",
    "Illumina Small RNA 5' Adapter": "GTTCAGAGTTCTACAGTCCGACGATC",
    "Nextera Transposase Sequence": "CTGTCTCTTATA",
    "SOLID Small RNA Adapter": "CGCCTTGGCCGT",
}
TENX_EXTRA_ADAPTERS = {
    "10x Template Switch Oligo": "AAGCAGTGGTATCAACGCAGAGTACATGGG",
    "poly-A": "AAAAAAAAAAAA",
}


# --------------------------------------------------------------------------- generic plumbing


def by_mate(value: object, options: Mapping[str, object], fn: Callable[[list[Rd]], dict[str, Any]]) -> dict[str, Any]:
    """Apply ``fn`` to a single read set, or to read1 and read2 of a pair separately."""

    first, second = parse_pair(value, options)
    if second is None:
        return fn(first)
    return {"paired": True, "read1": fn(first), "read2": fn(second)}


def tenx_regions(value: object, options: Mapping[str, object]) -> tuple[dict[str, list[Rd]], dict[str, Any]]:
    """Split 10x Chromium reads into barcode, UMI and RNA regions.

    ``value`` is ``{'read1': ..., 'read2': ...}`` where read 1 = cell barcode + UMI and read 2 = cDNA.
    ``options['chemistry']`` is v1 (14+10), v2 (16+10) or v3 (16+12, default); ``barcode_length`` and
    ``umi_length`` override it.
    """

    first, second = parse_pair(value, options)
    if second is None:
        raise InputFormatError(
            "10x functions need a read pair: {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>}.",
            code="PAIR_REQUIRED",
            hint="Pass read1/read2 FASTQ text, paths or read lists.",
        )
    chemistry = opt_choice(options, "chemistry", "v3", tuple(TENX_CHEMISTRY))
    default_bc, default_umi = TENX_CHEMISTRY[chemistry]
    bc_len = opt_int(options, "barcode_length", default_bc, minimum=1, maximum=64)
    umi_len = opt_int(options, "umi_length", default_umi, minimum=1, maximum=64)
    too_short = sum(1 for read in first if len(read.seq) < bc_len + umi_len)
    if too_short:
        raise InputFormatError(
            f"{too_short} read-1 sequences are shorter than barcode_length+umi_length={bc_len + umi_len}.",
            code="READ1_TOO_SHORT",
            hint="Check options['chemistry'] / barcode_length / umi_length.",
        )

    def cut(read: Rd, start: int, stop: int) -> Rd:
        return Rd(read.name, read.seq[start:stop], None if read.qual is None else read.qual[start:stop], read.comment)

    regions = {
        "barcode": [cut(r, 0, bc_len) for r in first],
        "umi": [cut(r, bc_len, bc_len + umi_len) for r in first],
        "rna_read": list(second),
    }
    meta = {"chemistry": chemistry, "barcode_length": bc_len, "umi_length": umi_len, "read_pairs": len(first)}
    return regions, meta


def _per_region(
    regions: Mapping[str, list[Rd]], fn: Callable[[list[Rd]], dict[str, Any]]
) -> dict[str, Any]:
    return {name: fn(reads) for name, reads in regions.items()}


# --------------------------------------------------------------------------- yield / length


def _yield(reads: list[Rd]) -> dict[str, Any]:
    bases = sum(len(r.seq) for r in reads)
    qbases = q20 = q30 = 0
    for r in reads:
        if r.qual:
            qbases += len(r.qual)
            q20 += sum(1 for q in r.qual if q >= 20)
            q30 += sum(1 for q in r.qual if q >= 30)
    gc = sum(r.seq.count("G") + r.seq.count("C") for r in reads)
    acgt = sum(sum(r.seq.count(b) for b in "ACGT") for r in reads)
    return {
        "reads": len(reads),
        "bases": bases,
        "mean_read_length": round6(bases / len(reads)) if reads else None,
        "n_bases": sum(r.seq.count("N") for r in reads),
        "gc_fraction": round6(gc / acgt) if acgt else None,
        "q20_bases": q20 if qbases else None,
        "q30_bases": q30 if qbases else None,
        "q20_fraction": round6(q20 / qbases) if qbases else None,
        "q30_fraction": round6(q30 / qbases) if qbases else None,
    }


def op_yield(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    genome = opt_num(options, "genome_size", None, minimum=1, maximum=1e12)
    if second is None:
        result: dict[str, Any] = _yield(first)
        total_bases = result["bases"]
    else:
        r1, r2 = _yield(first), _yield(second)
        total_bases = r1["bases"] + r2["bases"]
        result = {"paired": True, "read_pairs": len(first), "read1": r1, "read2": r2, "total_bases": total_bases}
    if genome is not None:
        result["genome_size"] = genome
        result["expected_coverage"] = round6(total_bases / genome)
    return result


def op_yield_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    barcodes = Counter(r.seq for r in regions["barcode"])
    per_barcode = sorted(barcodes.values())
    result: dict[str, Any] = {**meta, "regions": _per_region(regions, _yield)}
    result["total_bases_read1_read2"] = sum(v["bases"] for v in result["regions"].values())
    result["distinct_barcodes"] = len(barcodes)
    result["median_reads_per_barcode"] = round6(percentile(per_barcode, 0.5))
    whitelist = whitelist_option(options)
    if whitelist:
        mm = opt_int(options, "max_mismatches", 1, minimum=0, maximum=3)
        matched = match_to_whitelist(barcodes, whitelist, mm)
        valid = sum(c for b, c in barcodes.items() if isinstance(matched[b], tuple))
        result["valid_barcode_fraction"] = round6(valid / len(regions["barcode"]))
        result["valid_barcode_definition"] = f"Barcode equals a whitelist entry or has a unique whitelist neighbour within {mm} mismatches."
    return result


def _length_stats(reads: list[Rd], options: Mapping[str, object]) -> dict[str, Any]:
    lengths = [len(r.seq) for r in reads]
    counts = Counter(lengths)
    ordered = sorted(lengths)
    span = ordered[-1] - ordered[0] + 1
    default_width = 1 if span <= 1000 else math.ceil(span / 200)
    width = opt_int(options, "bin_width", default_width, minimum=1, maximum=1_000_000)
    binned: Counter[int] = Counter()
    for length, count in counts.items():
        binned[(length // width) * width] += count
    return {
        "reads": len(lengths),
        "bases": sum(lengths),
        "min": ordered[0],
        "max": ordered[-1],
        "mean": round6(sum(lengths) / len(lengths)),
        "median": round6(percentile(ordered, 0.5)),
        "stdev": round6(stdev(lengths)),
        "n50": nx(lengths, 0.5),
        "n90": nx(lengths, 0.9),
        "l50": lx(lengths, 0.5),
        "bin_width": width,
        "histogram": [{"start": s, "end": s + width - 1, "count": binned[s]} for s in sorted(binned)],
    }


def op_read_length(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return by_mate(value, options, lambda reads: _length_stats(reads, options))


# --------------------------------------------------------------------------- quality


def _qual_hist(reads: Sequence[Rd]) -> list[list[int]]:
    maxlen = max((len(r.qual) for r in reads if r.qual), default=0)
    hist = [[0] * (MAX_PHRED + 1) for _ in range(maxlen)]
    for r in reads:
        if r.qual:
            for index, q in enumerate(r.qual):
                hist[index][q] += 1
    return hist


def _per_position_quality(reads: list[Rd]) -> list[dict[str, Any]]:
    rows = []
    for index, counts in enumerate(_qual_hist(reads)):
        n = sum(counts)
        rows.append(
            {
                "position": index + 1,
                "bases": n,
                "mean": round6(sum(q * c for q, c in enumerate(counts)) / n),
                "median": round6(hist_percentile(counts, 0.5)),
                "q1": round6(hist_percentile(counts, 0.25)),
                "q3": round6(hist_percentile(counts, 0.75)),
                "p10": round6(hist_percentile(counts, 0.10)),
                "p90": round6(hist_percentile(counts, 0.90)),
                "min": kth(counts, 0),
                "max": kth(counts, n - 1),
            }
        )
    return rows


def _quality(reads: list[Rd]) -> dict[str, Any]:
    need_quality(reads, "quality statistics")
    qreads = [r for r in reads if r.qual]
    read_means = [sum(r.qual) / len(r.qual) for r in qreads if r.qual]
    dist = Counter(int(m + 0.5) for m in read_means)
    ordered = sorted(read_means)
    return {
        "reads_with_quality": len(qreads),
        "per_position": _per_position_quality(reads),
        "per_read_mean_quality": {
            "mean": round6(mean(read_means)),
            "median": round6(percentile(ordered, 0.5)),
            "mode": min(dist, key=lambda q: (-dist[q], q)) if dist else None,
            "distribution": {str(q): dist[q] for q in sorted(dist)},
        },
    }


def op_quality(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return by_mate(value, options, _quality)


def _qscore_summary(reads: list[Rd], read_method: str) -> dict[str, Any]:
    need_quality(reads, "Q-score summary")
    qreads = [r for r in reads if r.qual]
    counts = [0] * (MAX_PHRED + 1)
    read_scores: list[float] = []
    for r in qreads:
        assert r.qual is not None
        for q in r.qual:
            counts[q] += 1
        if not r.qual:
            continue
        if read_method == "probability":
            score = p_to_q(sum(q_to_p(q) for q in r.qual) / len(r.qual))
            read_scores.append(93.0 if score is None else score)
        else:
            read_scores.append(sum(r.qual) / len(r.qual))
    n = sum(counts)
    mean_error = sum(q_to_p(q) * c for q, c in enumerate(counts)) / n
    ordered = sorted(read_scores)
    thresholds = (7, 10, 15, 20, 30)
    return {
        "bases": n,
        "mean_q": round6(sum(q * c for q, c in enumerate(counts)) / n),
        "q_of_mean_error_probability": round6(p_to_q(mean_error)),
        "median_q": round6(hist_percentile(counts, 0.5)),
        "min_q": kth(counts, 0),
        "max_q": kth(counts, n - 1),
        "fraction_bases_at_least": {
            f"Q{t}": round6(sum(counts[t:]) / n) for t in (10, 20, 30, 40)
        },
        "read_q_method": read_method,
        "reads": len(read_scores),
        "mean_read_q": round6(mean(read_scores)),
        "median_read_q": round6(percentile(ordered, 0.5)),
        "reads_with_mean_q_at_least": {
            f"Q{t}": sum(1 for s in read_scores if s >= t) for t in thresholds
        },
    }


def op_qscore_summary(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    method = opt_choice(options, "read_q_method", "probability", ("probability", "arithmetic"))
    return by_mate(value, options, lambda reads: _qscore_summary(reads, method))


_ILLUMINA_FILTER = re.compile(r"^\d+:([YN]):\d+")


def _basecalling(reads: list[Rd], min_mean_q: float) -> dict[str, Any]:
    need_quality(reads, "basecalling QC")
    qreads = [r for r in reads if r.qual]
    n_bases = sum(len(r.seq) for r in reads)
    n_calls = sum(r.seq.count("N") for r in reads)
    read_means = [sum(r.qual) / len(r.qual) for r in qreads if r.qual]
    all_q = [q for r in qreads for q in (r.qual or ())]
    flags = [m.group(1) for r in reads if (m := _ILLUMINA_FILTER.match(r.comment))]
    result: dict[str, Any] = {
        "reads": len(reads),
        "bases": n_bases,
        "no_call_fraction": round6(n_calls / n_bases) if n_bases else None,
        "reads_with_no_calls": sum(1 for r in reads if "N" in r.seq),
        "mean_q": round6(mean(all_q)),
        "q30_fraction": round6(sum(1 for q in all_q if q >= 30) / len(all_q)),
        "expected_errors_per_read": round6(mean([sum(q_to_p(q) for q in r.qual or ()) for r in qreads])),
        "min_mean_q": min_mean_q,
        "reads_below_min_mean_q": sum(1 for m in read_means if m < min_mean_q),
        "fraction_reads_below_min_mean_q": round6(sum(1 for m in read_means if m < min_mean_q) / len(read_means)),
        "illumina_filter_flag_reads": len(flags),
    }
    if flags:
        result["filtered_reads"] = flags.count("Y")
        result["pass_filter_fraction"] = round6(flags.count("N") / len(flags))
    else:
        result["pass_filter_fraction"] = None
    return result


def op_basecalling_qc(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    min_q = opt_num(options, "min_mean_quality", 20.0, minimum=0, maximum=MAX_PHRED) or 0.0
    return by_mate(value, options, lambda reads: _basecalling(reads, min_q))


def _error_profile(reads: list[Rd]) -> dict[str, Any]:
    need_quality(reads, "error profile")
    hist = _qual_hist(reads)
    per_position = []
    cumulative = 0.0
    for index, counts in enumerate(hist):
        n = sum(counts)
        rate = sum(q_to_p(q) * c for q, c in enumerate(counts)) / n
        cumulative += rate
        per_position.append({"position": index + 1, "bases": n, "expected_error_rate": rate, "cumulative_expected_errors": cumulative})
    total_bases = sum(sum(c) for c in hist)
    expected = sum(row["expected_error_rate"] * row["bases"] for row in per_position)
    n_by_pos = [sum(1 for r in reads if len(r.seq) > i and r.seq[i] == "N") for i in range(len(hist))]
    for row, ns in zip(per_position, n_by_pos, strict=False):
        row["n_fraction"] = ns / row["bases"] if row["bases"] else 0.0
    return {
        "source": "phred_quality",
        "overall_expected_error_rate": expected / total_bases,
        "overall_q_equivalent": round6(p_to_q(expected / total_bases)),
        "expected_errors_per_read": expected / len([r for r in reads if r.qual]),
        "per_position": per_position,
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate.",
    }


def op_error_profile(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return by_mate(value, options, _error_profile)


# --------------------------------------------------------------------------- GC


def _norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _gc_bias(reads: list[Rd], reference_gc: float | None) -> dict[str, Any]:
    percents = []
    for r in reads:
        fraction = gc_fraction(r.seq)
        if fraction is not None:
            percents.append(fraction * 100.0)
    if not percents:
        raise InputFormatError("no read contains an unambiguous A/C/G/T base.", code="NO_ACGT", hint="Check the input sequences.")
    n = len(percents)
    bins = Counter(round(p) for p in percents)
    mu = sum(percents) / n
    sd = stdev(percents) or 0.0
    deviation = None
    theoretical: list[float] = []
    if sd > 0:
        theoretical = [n * (_norm_cdf((i + 0.5 - mu) / sd) - _norm_cdf((i - 0.5 - mu) / sd)) for i in range(101)]
        deviation = 100.0 * sum(abs(bins.get(i, 0) - theoretical[i]) for i in range(101)) / n
    maxlen = max(len(r.seq) for r in reads)
    by_position = []
    for i in range(maxlen):
        g = sum(1 for r in reads if len(r.seq) > i and r.seq[i] in "GC")
        a = sum(1 for r in reads if len(r.seq) > i and r.seq[i] in "ACGT")
        by_position.append(round6(100.0 * g / a) if a else None)
    result = {
        "reads_with_gc": n,
        "mean_gc_percent": round6(mu),
        "stdev_gc_percent": round6(sd),
        "mode_gc_percent": min(bins, key=lambda k: (-bins[k], k)),
        "histogram": {str(i): bins[i] for i in sorted(bins)},
        "theoretical_normal_counts": [round(t, 4) for t in theoretical] if theoretical else None,
        "deviation_from_normal_percent": round6(deviation),
        "gc_percent_by_position": by_position,
        "definition": "Per-read GC% (G+C over A+C+G+T) histogrammed in 1% bins; the theoretical curve is a normal "
        "distribution with the observed mean and standard deviation (FastQC convention).",
    }
    if reference_gc is not None:
        result["reference_gc_percent"] = reference_gc * 100.0
        result["gc_shift_vs_reference_percent"] = round6(mu - reference_gc * 100.0)
    return result


def op_gc_bias(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    ref = opt_num(options, "reference_gc", None, minimum=0.0, maximum=1.0)
    return by_mate(value, options, lambda reads: _gc_bias(reads, ref))


# --------------------------------------------------------------------------- adapters


def adapters_option(options: Mapping[str, object], extra: Mapping[str, str] | None = None) -> dict[str, str]:
    raw = options.get("adapters")
    adapters: dict[str, str] = {}
    if raw is None:
        adapters = dict(DEFAULT_ADAPTERS)
        if extra:
            adapters.update(extra)
        return adapters
    items: list[tuple[str, object]]
    if isinstance(raw, Mapping):
        items = [(str(k), v) for k, v in raw.items()]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
        items = [(str(v), v) for v in raw]
    elif isinstance(raw, str):
        items = [(raw, raw)]
    else:
        raise ConfigurationError("options['adapters'] must be a mapping name->sequence or a list of sequences.", code="INVALID_OPTION")
    for name, seq in items:
        if not isinstance(seq, str) or not re.fullmatch(r"[ACGTacgt]+", seq):
            raise ConfigurationError(f"adapter {name!r} must be an A/C/G/T string.", code="INVALID_OPTION")
        adapters[name] = seq.upper()
    if not adapters:
        raise ConfigurationError("options['adapters'] is empty.", code="INVALID_OPTION")
    return adapters


def find_first(seq: str, kmer: str, mismatches: int) -> int:
    """Index of the first occurrence of ``kmer`` in ``seq`` allowing Hamming ``mismatches``; -1 if none."""

    if mismatches == 0:
        return seq.find(kmer)
    k = len(kmer)
    for start in range(len(seq) - k + 1):
        bad = 0
        for a, b in zip(seq[start : start + k], kmer, strict=False):
            if a != b:
                bad += 1
                if bad > mismatches:
                    break
        else:
            return start
    return -1


def _adapter_content(reads: list[Rd], adapters: Mapping[str, str], k: int, mismatches: int) -> dict[str, Any]:
    maxlen = max(len(r.seq) for r in reads)
    n = len(reads)
    per_adapter: dict[str, Any] = {}
    any_hit = [False] * n
    for name, seq in adapters.items():
        kmer = seq[: min(k, len(seq))]
        first_hits = Counter()
        for index, r in enumerate(reads):
            pos = find_first(r.seq, kmer, mismatches)
            if pos >= 0:
                first_hits[pos + 1] += 1
                any_hit[index] = True
        running = 0
        cumulative = []
        for position in range(1, maxlen + 1):
            running += first_hits.get(position, 0)
            cumulative.append(round6(100.0 * running / n))
        per_adapter[name] = {
            "sequence": seq,
            "kmer_used": kmer,
            "reads_with_adapter": sum(first_hits.values()),
            "fraction_reads_with_adapter": round6(sum(first_hits.values()) / n),
            "max_cumulative_percent": cumulative[-1] if cumulative else 0.0,
            "cumulative_percent_by_position": cumulative,
        }
    return {
        "reads": n,
        "kmer_size": k,
        "max_mismatches": mismatches,
        "reads_with_any_adapter": sum(any_hit),
        "fraction_reads_with_any_adapter": round6(sum(any_hit) / n),
        "adapters": per_adapter,
        "definition": "A read hits an adapter when the first min(kmer_size, len(adapter)) bases of the adapter occur in the "
        "read (Hamming <= max_mismatches). The cumulative percentage at position p counts reads whose first hit "
        "starts at or before p (1-based), as in FastQC. Adapters truncated by the read end are not detected.",
    }


def op_adapter_content(value: object, options: Mapping[str, object], warnings: Warnings, *, tenx: bool = False) -> dict[str, Any]:
    k = opt_int(options, "kmer_size", 12, minimum=4, maximum=64)
    mm = opt_int(options, "max_mismatches", 0, minimum=0, maximum=3)
    adapters = adapters_option(options, TENX_EXTRA_ADAPTERS if tenx else None)
    if tenx:
        regions, meta = tenx_regions(value, options)
        return {**meta, "region": "rna_read (read 2)", **_adapter_content(regions["rna_read"], adapters, k, mm)}
    return by_mate(value, options, lambda reads: _adapter_content(reads, adapters, k, mm))


# --------------------------------------------------------------------------- duplication


def _dup_key(read: Rd, prefix: int | None, canonical: bool) -> str:
    seq = read.seq if prefix is None else read.seq[:prefix]
    if canonical:
        other = revcomp(read.seq)
        other = other if prefix is None else other[:prefix]
        return min(seq, other)
    return seq


def _keys(value: object, options: Mapping[str, object]) -> tuple[list[str], bool]:
    prefix = opt_int(options, "prefix_length", 0, minimum=0, maximum=100_000) or None
    canonical = opt_bool(options, "reverse_complement_aware", False)
    first, second = parse_pair(value, options)
    if second is None:
        return [_dup_key(r, prefix, canonical) for r in first], False
    return [_dup_key(a, prefix, canonical) + "|" + _dup_key(b, prefix, canonical) for a, b in zip(first, second, strict=False)], True


def op_duplicate_rate(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    keys, paired = _keys(value, options)
    counts = Counter(keys)
    n = len(keys)
    duplicates = n - len(counts)
    return {
        "paired": paired,
        "reads": n,
        "distinct_sequences": len(counts),
        "duplicate_reads": duplicates,
        "duplicate_rate": round6(duplicates / n),
        "unique_fraction": round6(len(counts) / n),
        "most_common": [{"sequence": s[:80], "count": c} for s, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:5]],
        "definition": "duplicate_rate = (reads - distinct sequences) / reads; sequences are compared exactly "
        "(read pairs by the concatenation of both mates).",
    }


_DUP_BINS = [(str(i), i, i) for i in range(1, 10)] + [
    (">10", 10, 49), (">50", 50, 99), (">100", 100, 499), (">500", 500, 999),
    (">1k", 1000, 4999), (">5k", 5000, 9999), (">10k+", 10000, 10**12),
]


def dup_levels(counts: Counter[str]) -> dict[str, Any]:
    n = sum(counts.values())
    reads_at: Counter[str] = Counter()
    distinct_at: Counter[str] = Counter()
    for count in counts.values():
        for label, low, high in _DUP_BINS:
            if low <= count <= high:
                reads_at[label] += count
                distinct_at[label] += 1
                break
    return {
        "reads": n,
        "distinct_sequences": len(counts),
        "percent_remaining_after_deduplication": round6(100.0 * len(counts) / n),
        "levels": [
            {
                "level": label,
                "distinct_sequences": distinct_at[label],
                "percent_of_total_reads": round6(100.0 * reads_at[label] / n),
                "percent_of_deduplicated": round6(100.0 * distinct_at[label] / len(counts)),
            }
            for label, _, _ in _DUP_BINS
        ],
        "definition": "Level = number of copies of an identical sequence; bins 1-9, >10 (10-49), >50, >100, >500, >1k, "
        ">5k, >10k+ as in FastQC. No extrapolation beyond the supplied reads is performed.",
    }


def op_duplication(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    keys, paired = _keys(value, options)
    return {"paired": paired, **dup_levels(Counter(keys))}


def op_duplication_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    prefix = opt_int(options, "cdna_prefix_length", 25, minimum=1, maximum=1000)
    keys = [f"{b.seq}|{u.seq}|{c.seq[:prefix]}" for b, u, c in zip(regions["barcode"], regions["umi"], regions["rna_read"], strict=False)]
    counts = Counter(keys)
    n = len(keys)
    result = dup_levels(counts)
    result.update(meta)
    result["duplicate_rate"] = round6((n - len(counts)) / n)
    result["sequencing_saturation_proxy"] = result["duplicate_rate"]
    result["definition"] = (
        f"Reads are duplicates when cell barcode, UMI and the first {prefix} bases of read 2 are identical. This is a "
        "proxy: Cell Ranger's sequencing saturation uses (barcode, UMI, gene) from an alignment, which is not available here."
    )
    return result


# --------------------------------------------------------------------------- complexity


def solve_library_size(total: int, unique: int) -> float | None:
    """Solve ``unique / X = 1 - exp(-total / X)`` for X (Lander-Waterman, as in Picard); ``None`` if unbounded."""

    if total <= 0 or unique <= 0 or unique > total:
        raise ConfigurationError("need 0 < unique_reads <= total_reads.", code="INVALID_OPTION")
    if unique == total:
        return None
    if total == 1:
        return None

    def f(x: float) -> float:
        return x * (1.0 - math.exp(-total / x)) - unique

    low = float(unique)  # f(unique) < 0 because x*(1-exp(-N/x)) < x
    high = float(unique) * 2.0
    while f(high) < 0:
        high *= 2.0
        if high > 1e18:
            return None
    for _ in range(200):
        mid = 0.5 * (low + high)
        if f(mid) < 0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def op_library_complexity(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    total = unique = None
    if isinstance(value, Mapping) and "total_reads" in value and "unique_reads" in value:
        total, unique = value["total_reads"], value["unique_reads"]
    elif "total_reads" in options and "unique_reads" in options:
        total, unique = options["total_reads"], options["unique_reads"]
    if total is not None:
        if isinstance(total, bool) or isinstance(unique, bool) or not isinstance(total, int) or not isinstance(unique, int):
            raise ConfigurationError("total_reads and unique_reads must be integers.", code="INVALID_OPTION")
        source = "counts"
    else:
        keys, _ = _keys(value, options)
        total, unique = len(keys), len(set(keys))
        source = "reads (exact-sequence duplicates)"
    size = solve_library_size(total, unique)
    result: dict[str, Any] = {
        "source": source,
        "total_reads": total,
        "unique_reads": unique,
        "duplicate_rate": round6((total - unique) / total),
    }
    if size is None:
        result["estimated_library_size"] = None
        warnings.append("Every read is distinct (or only one read given): the library size is unbounded by this data.")
    else:
        result["estimated_library_size"] = round6(size)
        result["predicted_unique_reads"] = {
            f"{f}x": round6(size * (1.0 - math.exp(-f * total / size))) for f in (0.5, 1, 2, 5, 10)
        }
    result["definition"] = "Lander-Waterman / Picard EstimateLibraryComplexity: find X with unique/X = 1 - exp(-total/X)."
    return result


# --------------------------------------------------------------------------- barcodes


def match_to_whitelist(
    observed: Mapping[str, int], whitelist: Sequence[str], max_mismatches: int
) -> dict[str, tuple[str, int] | str | None]:
    """Map each observed barcode to ``(whitelist_entry, distance)``, ``'ambiguous'`` or ``None``."""

    wl = set(whitelist)
    big = len(wl) > 5000
    if big and max_mismatches > 1:
        raise ConfigurationError(
            "max_mismatches > 1 is only supported for whitelists up to 5000 entries.", code="INVALID_OPTION"
        )
    out: dict[str, tuple[str, int] | str | None] = {}
    wl_list = sorted(wl)
    unmatched = sum(1 for obs in observed if obs not in wl)
    if not big and max_mismatches > 0 and unmatched * len(wl) > 20_000_000:
        raise ConfigurationError(
            f"{unmatched} distinct non-whitelist barcodes x {len(wl)} whitelist entries is too large for brute-force matching.",
            code="INPUT_TOO_LARGE",
            hint="Use max_mismatches=0 or <= 1 with a larger whitelist, or fewer reads.",
        )
    for obs in observed:
        if obs in wl:
            out[obs] = (obs, 0)
            continue
        if max_mismatches == 0:
            out[obs] = None
            continue
        if big:
            hits = set()
            for i, base in enumerate(obs):
                for alt in "ACGT":
                    if alt != base:
                        cand = obs[:i] + alt + obs[i + 1 :]
                        if cand in wl:
                            hits.add(cand)
            out[obs] = "ambiguous" if len(hits) > 1 else ((next(iter(hits)), 1) if hits else None)
            continue
        best = max_mismatches + 1
        winners: list[str] = []
        for cand in wl_list:
            if len(cand) != len(obs):
                continue
            d = hamming(obs, cand)
            if d < best:
                best, winners = d, [cand]
            elif d == best:
                winners.append(cand)
        if best > max_mismatches:
            out[obs] = None
        elif len(winners) > 1:
            out[obs] = "ambiguous"
        else:
            out[obs] = (winners[0], best)
    return out


def whitelist_option(options: Mapping[str, object]) -> list[str]:
    raw = options.get("whitelist")
    if raw is None:
        return []
    if isinstance(raw, str):
        text = as_text(raw) or ""
        raw = [line.split()[0] for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not isinstance(raw, Sequence) or isinstance(raw, (bytes, bytearray)):
        raise ConfigurationError("options['whitelist'] must be a list of barcodes or newline text.", code="INVALID_OPTION")
    if len(raw) > 100_000:
        raise ConfigurationError("whitelist has more than 100000 barcodes.", code="INVALID_OPTION")
    cleaned = []
    for item in raw:
        if not isinstance(item, str) or not re.fullmatch(r"[ACGTacgt]+", item):
            raise ConfigurationError(f"invalid whitelist barcode {item!r}.", code="INVALID_OPTION")
        cleaned.append(item.upper())
    return cleaned


def header_index(read: Rd) -> str:
    """Sample index at the end of an Illumina header comment (``1:N:0:ATCACG+TTAGGC``)."""

    text = read.comment or read.name
    token = text.split(":")[-1].strip()
    if not re.fullmatch(r"[ACGTN+]+", token):
        raise InputFormatError(
            f"read {read.name!r} has no index sequence at the end of its header ({text!r}).",
            code="NO_HEADER_INDEX",
            hint="Use barcode_location 'read_start'/'read_end' or reads with Illumina-style headers.",
        )
    return token


def _barcode_region(reads: Sequence[Rd], options: Mapping[str, object], default_len: int | None) -> tuple[list[str], list[tuple[int, ...] | None], str]:
    location = opt_choice(options, "barcode_location", "read_start", ("read_start", "read_end", "header"))
    if location == "header":
        return [header_index(r) for r in reads], [None] * len(reads), location
    length = opt_int(options, "barcode_length", default_len or 8, minimum=1, maximum=64)
    short = [r.name for r in reads if len(r.seq) < length]
    if short:
        raise InputFormatError(f"{len(short)} reads are shorter than barcode_length={length}.", code="READ_TOO_SHORT")
    if location == "read_start":
        return [r.seq[:length] for r in reads], [None if r.qual is None else r.qual[:length] for r in reads], location
    return [r.seq[-length:] for r in reads], [None if r.qual is None else r.qual[-length:] for r in reads], location


def _barcode_quality(reads: list[Rd], options: Mapping[str, object]) -> dict[str, Any]:
    barcodes, quals, location = _barcode_region(reads, options, None)
    n = len(reads)
    counts = Counter(barcodes)
    result: dict[str, Any] = {
        "reads": n,
        "barcode_source": location,
        "distinct_barcodes": len(counts),
        "reads_with_N_in_barcode": sum(1 for b in barcodes if "N" in b),
        "top_barcodes": [{"barcode": b, "reads": c} for b, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:10]],
        "barcode_entropy_bits": round6(-sum((c / n) * math.log2(c / n) for c in counts.values())),
    }
    scored = [q for q in quals if q]
    if scored:
        flat = [x for q in scored for x in q]
        length = max(len(q) for q in scored)
        result["mean_quality"] = round6(sum(flat) / len(flat))
        result["q30_fraction"] = round6(sum(1 for x in flat if x >= 30) / len(flat))
        result["mean_quality_by_position"] = [
            round6(sum(q[i] for q in scored) / len(scored)) for i in range(length)
        ]
    else:
        result["mean_quality"] = None
        result["q30_fraction"] = None
    whitelist = whitelist_option(options)
    if whitelist:
        mm = opt_int(options, "max_mismatches", 1, minimum=0, maximum=3)
        matched = match_to_whitelist(counts, whitelist, mm)
        exact = corrected = ambiguous = invalid = 0
        for barcode, c in counts.items():
            m = matched[barcode]
            if m is None:
                invalid += c
            elif m == "ambiguous":
                ambiguous += c
            else:
                assert isinstance(m, tuple)
                if m[1] == 0:
                    exact += c
                else:
                    corrected += c
        result["whitelist"] = {
            "size": len(whitelist),
            "max_mismatches": mm,
            "exact_match_fraction": round6(exact / n),
            "corrected_fraction": round6(corrected / n),
            "ambiguous_fraction": round6(ambiguous / n),
            "invalid_fraction": round6(invalid / n),
        }
    return result


def op_barcode_quality(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    return _barcode_quality(reads, options)


def op_barcode_quality_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    bc_reads = regions["barcode"]
    result = {**meta, "cell_barcode": _barcode_quality(bc_reads, {**options, "barcode_location": "read_start", "barcode_length": meta["barcode_length"]})}
    umi = regions["umi"]
    umi_summary: dict[str, Any] = {"distinct_umis": len({u.seq for u in umi}), "reads_with_N_in_umi": sum(1 for u in umi if "N" in u.seq)}
    flat = [q for u in umi if u.qual for q in u.qual]
    if flat:
        umi_summary["mean_quality"] = round6(sum(flat) / len(flat))
        umi_summary["q30_fraction"] = round6(sum(1 for q in flat if q >= 30) / len(flat))
    homopolymer = sum(1 for u in umi if len(set(u.seq)) == 1)
    umi_summary["homopolymer_umi_reads"] = homopolymer
    result["umi"] = umi_summary
    return result


def _parse_samples(options: Mapping[str, object]) -> dict[str, str]:
    raw = options.get("samples", options.get("barcodes"))
    if raw is None:
        raise ConfigurationError(
            "options['samples'] is required: a mapping sample name -> barcode (dual indexes as 'ACGT+TTAG').",
            code="MISSING_OPTION",
            hint="Example: options={'samples': {'S1': 'ACGTAC', 'S2': 'TGCATG'}}.",
        )
    items: list[tuple[str, object]]
    if isinstance(raw, Mapping):
        items = [(str(k), v) for k, v in raw.items()]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
        items = []
        for entry in raw:
            if not isinstance(entry, Mapping) or "sample" not in entry or "barcode" not in entry:
                raise ConfigurationError("samples list entries must be {'sample':..., 'barcode':...}.", code="INVALID_OPTION")
            items.append((str(entry["sample"]), entry["barcode"]))
    else:
        raise ConfigurationError("options['samples'] must be a mapping or a list.", code="INVALID_OPTION")
    samples: dict[str, str] = {}
    for name, barcode in items:
        if not isinstance(barcode, str) or not re.fullmatch(r"[ACGTacgt]+(\+[ACGTacgt]+)?", barcode):
            raise ConfigurationError(f"invalid barcode {barcode!r} for sample {name!r}.", code="INVALID_OPTION")
        if not re.fullmatch(r"[A-Za-z0-9_.\-]{1,64}", name):
            raise ConfigurationError(f"sample name {name!r} must match [A-Za-z0-9_.-]{{1,64}}.", code="INVALID_OPTION")
        samples[name] = barcode.upper()
    if not samples:
        raise ConfigurationError("options['samples'] is empty.", code="INVALID_OPTION")
    return samples


def _barcode_distance(observed: str, expected: str) -> int | None:
    a, b = observed.replace("+", ""), expected.replace("+", "")
    if len(a) != len(b):
        return None
    return hamming(a, b)


def assign_reads(reads: Sequence[Rd], options: Mapping[str, object]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    samples = _parse_samples(options)
    mm = opt_int(options, "max_mismatches", 1, minimum=0, maximum=4)
    gap = opt_int(options, "min_distance_gap", 1, minimum=0, maximum=8)
    lengths = {len(code.replace("+", "")) for code in samples.values()}
    default_len = lengths.pop() if len(lengths) == 1 else None
    barcodes, _, location = _barcode_region(reads, options, default_len)
    cache: dict[str, tuple[str | None, int | None, str]] = {}
    rows = []
    for read, barcode in zip(reads, barcodes, strict=False):
        if barcode not in cache:
            ranked = []
            for name, code in samples.items():
                d = _barcode_distance(barcode, code)
                if d is not None:
                    ranked.append((d, name))
            ranked.sort()
            if not ranked:
                raise InputFormatError(
                    f"observed barcode {barcode!r} has a different length from every sample barcode.",
                    code="BARCODE_LENGTH_MISMATCH", hint="Set barcode_length or fix the sample sheet.",
                )
            best = ranked[0]
            if best[0] > mm:
                cache[barcode] = (None, best[0], "unassigned")
            elif len(ranked) > 1 and ranked[1][0] - best[0] < gap:
                cache[barcode] = (None, best[0], "ambiguous")
            else:
                cache[barcode] = (best[1], best[0], "assigned")
        sample, dist, status = cache[barcode]
        rows.append({"id": read.name, "barcode": barcode, "sample": sample, "distance": dist, "status": status})
    counts = Counter(r["sample"] for r in rows if r["sample"])
    unassigned = sum(1 for r in rows if r["status"] == "unassigned")
    ambiguous = sum(1 for r in rows if r["status"] == "ambiguous")
    summary = {
        "reads": len(rows),
        "barcode_location": location,
        "max_mismatches": mm,
        "min_distance_gap": gap,
        "reads_per_sample": {name: counts.get(name, 0) for name in samples},
        "assigned": sum(counts.values()),
        "unassigned": unassigned,
        "ambiguous": ambiguous,
        "assigned_fraction": round6(sum(counts.values()) / len(rows)),
    }
    codes = [c.replace("+", "") for c in samples.values()]
    if len(codes) > 1:
        closest = min(hamming(a, b) for i, a in enumerate(codes) for b in codes[i + 1 :] if len(a) == len(b)) if any(
            len(a) == len(b) for i, a in enumerate(codes) for b in codes[i + 1 :]
        ) else None
        summary["min_distance_between_sample_barcodes"] = closest
    return rows, summary


def op_assign_barcodes(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    rows, summary = assign_reads(reads, options)
    limit = opt_int(options, "max_listed", 100_000, minimum=1, maximum=1_000_000)
    if len(rows) > limit:
        warnings.append(f"per-read assignments truncated to the first {limit} reads (see max_listed).")
    d = summary.get("min_distance_between_sample_barcodes")
    if d is not None and d <= 2 * summary["max_mismatches"]:
        warnings.append(
            f"sample barcodes differ by only {d} bases; with max_mismatches={summary['max_mismatches']} some reads may be ambiguous."
        )
    return {**summary, "assignments": rows[:limit]}


# --------------------------------------------------------------------------- insert size


_SAM_PAIRED, _SAM_UNMAPPED, _SAM_FIRST, _SAM_SECONDARY, _SAM_SUPP = 0x1, 0x4, 0x40, 0x100, 0x800


def overlap_insert_size(r1: str, r2: str, min_overlap: int, max_mismatch_rate: float) -> int | None:
    """Fragment length from the overlap of read 1 with the reverse complement of read 2.

    Offsets ``o`` (start of rc(read2) relative to read1) are tried from the largest overlap down; the insert
    size is ``o + len(read2)``. Returns ``None`` when no overlap of at least ``min_overlap`` bases agrees.
    """

    rc2 = revcomp(r2)
    l1, l2 = len(r1), len(r2)
    candidates = []
    for o in range(-(l2 - min_overlap), l1 - min_overlap + 1):
        a_start, b_start = max(o, 0), max(-o, 0)
        overlap = min(l1 - a_start, l2 - b_start)
        if overlap >= min_overlap:
            candidates.append((overlap, o, a_start, b_start))
    for overlap, o, a_start, b_start in sorted(candidates, key=lambda c: (-c[0], c[1])):
        budget = max_mismatch_rate * overlap
        mism = 0
        for x, y in zip(r1[a_start : a_start + overlap], rc2[b_start : b_start + overlap], strict=False):
            if x != y:
                mism += 1
                if mism > budget:
                    break
        else:
            return o + l2
    return None


def insert_sizes(value: object, options: Mapping[str, object], warnings: Warnings) -> tuple[list[int], dict[str, Any]]:
    """Collect fragment lengths from paired reads, SAM text, a 10x fragments table or numbers."""

    meta: dict[str, Any] = {}
    if isinstance(options.get("insert_sizes"), Sequence) and not isinstance(options.get("insert_sizes"), (str, bytes)):
        value = options["insert_sizes"]
    if isinstance(value, Mapping) and "insert_sizes" in value:
        value = value["insert_sizes"]
    if isinstance(value, Mapping) and "alignments" in value:
        value = value["alignments"]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value and all(
        isinstance(v, (int, float)) and not isinstance(v, bool) for v in value
    ):
        sizes = [int(v) for v in value]
        if any(v <= 0 for v in sizes):
            raise InputFormatError("insert sizes must be positive.", code="INVALID_INSERT_SIZE")
        meta["source"] = "numbers"
        return sizes, meta
    if isinstance(value, Mapping) or (
        isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value and not isinstance(value[0], str)
    ):
        first, second = parse_pair(value, options)
        if second is None:
            raise InputFormatError("insert size from reads needs a read pair (read1/read2).", code="PAIR_REQUIRED")
        min_ov = opt_int(options, "min_overlap", 10, minimum=3, maximum=1000)
        rate = opt_num(options, "max_mismatch_rate", 0.1, minimum=0.0, maximum=0.5) or 0.0
        sizes = []
        no_overlap = 0
        pair_limit = opt_int(options, "max_pairs", 200_000, minimum=1, maximum=5_000_000)
        if len(first) > pair_limit:
            warnings.append(f"only the first {pair_limit} of {len(first)} pairs were used for overlap sizing (max_pairs).")
        for a, b in zip(first[:pair_limit], second[:pair_limit], strict=False):
            size = overlap_insert_size(a.seq, b.seq, min_ov, rate)
            if size is None:
                no_overlap += 1
            else:
                sizes.append(size)
        meta.update({"source": "read-pair overlap", "pairs": min(len(first), pair_limit), "pairs_without_overlap": no_overlap, "min_overlap": min_ov, "max_mismatch_rate": rate})
        if no_overlap:
            warnings.append(f"{no_overlap} pairs have no detectable overlap (insert longer than both reads combined) and are excluded.")
        return sizes, meta
    text = as_text(value) if not isinstance(value, Sequence) or isinstance(value, str) else "\n".join(str(v) for v in value)
    if text is None:
        raise InputFormatError("unsupported insert-size input.", code="INVALID_INPUT")
    sizes = []
    mode = None
    skipped = 0
    for line in text.splitlines():
        if not line.strip() or line.startswith(("@", "#", "track", "browser")):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 11 and fields[1].isdigit() and fields[3].lstrip("-").isdigit():  # SAM
            mode = "SAM TLEN"
            flag, tlen = int(fields[1]), int(fields[8])
            if not flag & _SAM_PAIRED or flag & (_SAM_UNMAPPED | _SAM_SECONDARY | _SAM_SUPP):
                continue
            if not flag & _SAM_FIRST:  # count each pair once, via the first mate
                continue
            if options.get("proper_pairs_only", True) and not flag & 0x2:
                skipped += 1
                continue
            if tlen != 0:
                sizes.append(abs(tlen))
        elif len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():  # 10x fragments / BED
            mode = "fragments (end - start)"
            start, end = int(fields[1]), int(fields[2])
            weight = int(fields[4]) if len(fields) >= 5 and fields[4].isdigit() else 1
            if end <= start:
                raise InputFormatError(f"fragment end <= start: {line[:80]!r}", code="INVALID_FRAGMENT")
            sizes.extend([end - start] * min(weight, 1_000_000))
        else:
            raise InputFormatError(f"unrecognised line (not SAM, fragments or BED): {line[:80]!r}", code="INVALID_INPUT")
        if len(sizes) > 10_000_000:
            raise InputFormatError("too many fragments.", code="INPUT_TOO_LARGE")
    if mode is None:
        raise InputFormatError("no alignments or fragments found.", code="EMPTY_INPUT")
    meta["source"] = mode
    if skipped:
        meta["non_proper_pairs_skipped"] = skipped
    return sizes, meta


def _size_summary(sizes: list[int], max_insert: int) -> dict[str, Any]:
    kept = [s for s in sizes if s <= max_insert]
    if not kept:
        raise InputFormatError(
            "no insert size could be determined: no read pair overlaps (or every size exceeds max_insert_size).",
            code="NO_INSERT_SIZES",
            hint="Overlap-based sizing needs fragments shorter than read 1 + read 2; use alignments (SAM TLEN) for longer inserts.",
        )
    ordered = sorted(kept)
    counts = Counter(kept)
    median = percentile(ordered, 0.5)
    assert median is not None
    mad = percentile(sorted(abs(s - median) for s in kept), 0.5)
    return {
        "count": len(kept),
        "excluded_above_max": len(sizes) - len(kept),
        "min": ordered[0],
        "max": ordered[-1],
        "mean": round6(sum(kept) / len(kept)),
        "median": median,
        "stdev": round6(stdev(kept)),
        "median_absolute_deviation": mad,
        "mode": min(counts, key=lambda s: (-counts[s], s)),
        "percentiles": {f"p{p}": round6(percentile(ordered, p / 100)) for p in (5, 25, 75, 95)},
    }


def op_insert_size(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    sizes, meta = insert_sizes(value, options, warnings)
    return {**meta, **_size_summary(sizes, opt_int(options, "max_insert_size", 100_000, minimum=1, maximum=10_000_000))}


def op_insert_size_distribution(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    sizes, meta = insert_sizes(value, options, warnings)
    max_insert = opt_int(options, "max_insert_size", 100_000, minimum=1, maximum=10_000_000)
    summary = _size_summary(sizes, max_insert)
    kept = [s for s in sizes if s <= max_insert]
    span = summary["max"] - summary["min"] + 1
    width = opt_int(options, "bin_width", 1 if span <= 2000 else math.ceil(span / 400), minimum=1, maximum=1_000_000)
    binned: Counter[int] = Counter((s // width) * width for s in kept)
    running = 0
    rows = []
    for start in sorted(binned):
        running += binned[start]
        rows.append({"start": start, "end": start + width - 1, "count": binned[start], "fraction": round6(binned[start] / len(kept)), "cumulative_fraction": round6(running / len(kept))})
    return {**meta, **summary, "bin_width": width, "histogram": rows}


# --------------------------------------------------------------------------- k-mer screens


def _canonical_kmers(seq: str, k: int) -> list[str]:
    out = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer or not set(kmer) <= set("ACGT"):
            continue
        other = revcomp(kmer)
        out.append(kmer if kmer <= other else other)
    return out


def _reference_sets(options: Mapping[str, object], k: int) -> dict[str, set[str]]:
    raw = options.get("references")
    if raw is None:
        raise ConfigurationError(
            "options['references'] is required: a mapping name -> reference sequence (or FASTA text).",
            code="MISSING_OPTION",
            hint="Example: options={'references': {'human': 'ACGT...', 'phiX': 'GAGTTTTATCGC...'}}.",
        )
    sets: dict[str, set[str]] = {}
    if isinstance(raw, Mapping):
        items = [(str(n), s) for n, s in raw.items()]
    else:
        records = parse_reads(raw, {}, name="references")
        items = [(r.name, r.seq) for r in records]
    for name, seq in items:
        if not isinstance(seq, str) or len(seq) > 50_000_000:
            raise ConfigurationError(f"reference {name!r} must be a sequence string (<= 50 Mb).", code="INVALID_OPTION")
        sets[name] = set(_canonical_kmers(seq.upper(), k))
        if not sets[name]:
            raise ConfigurationError(f"reference {name!r} yields no {k}-mers.", code="INVALID_OPTION", hint="Use a longer reference or a smaller kmer_size.")
    if not sets:
        raise ConfigurationError("options['references'] is empty.", code="INVALID_OPTION")
    return sets


def _contamination(reads: list[Rd], options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    k = opt_int(options, "kmer_size", 21, minimum=4, maximum=64)
    min_fraction = opt_num(options, "min_kmer_fraction", 0.5, minimum=0.01, maximum=1.0) or 0.5
    refs = _reference_sets(options, k)
    expected = options.get("expected_reference")
    if expected is not None and expected not in refs:
        raise ConfigurationError(f"expected_reference {expected!r} is not one of the references.", code="INVALID_OPTION")
    hits_by_ref: Counter[str] = Counter()
    only_other = no_hit = multi = skipped = 0
    for r in reads:
        kmers = _canonical_kmers(r.seq, k)
        if not kmers:
            skipped += 1
            continue
        hit = [name for name, kset in refs.items() if sum(1 for x in kmers if x in kset) / len(kmers) >= min_fraction]
        for name in hit:
            hits_by_ref[name] += 1
        if not hit:
            no_hit += 1
        elif len(hit) > 1:
            multi += 1
        if hit and expected is not None and expected not in hit:
            only_other += 1
    scored = len(reads) - skipped
    if skipped:
        warnings.append(f"{skipped} reads have no usable {k}-mer (too short or ambiguous) and were skipped.")
    if scored == 0:
        raise InputFormatError(f"no read is long enough for kmer_size={k}.", code="READ_TOO_SHORT")
    result: dict[str, Any] = {
        "reads_screened": scored,
        "kmer_size": k,
        "min_kmer_fraction": min_fraction,
        "hits": {name: {"reads": hits_by_ref[name], "fraction": round6(hits_by_ref[name] / scored)} for name in refs},
        "reads_without_hit": no_hit,
        "reads_with_multiple_hits": multi,
        "definition": "A read hits a reference when at least min_kmer_fraction of its canonical k-mers occur in that reference "
        "(exact k-mer containment, both strands). It is a k-mer screen, not an alignment.",
    }
    if expected is not None:
        result["expected_reference"] = expected
        result["contamination_fraction"] = round6(only_other / scored)
        result["contamination_definition"] = "Fraction of screened reads that hit at least one non-expected reference and not the expected one."
    return result


def op_contamination(value: object, options: Mapping[str, object], warnings: Warnings, *, tenx: bool = False) -> dict[str, Any]:
    if tenx:
        regions, meta = tenx_regions(value, options)
        return {**meta, "region": "rna_read (read 2)", **_contamination(regions["rna_read"], options, warnings)}
    return by_mate(value, options, lambda reads: _contamination(reads, options, warnings))


def op_cross_sample_contamination(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    if not isinstance(value, Mapping) or len(value) < 2:
        raise InputFormatError(
            "input must be a mapping sample name -> reads with at least two samples.",
            code="INVALID_INPUT",
            hint="Example: {'A': fastq_text_a, 'B': fastq_text_b}.",
        )
    k = opt_int(options, "kmer_size", 21, minimum=4, maximum=64)
    threshold = opt_num(options, "threshold", 0.1, minimum=0.0, maximum=1.0) or 0.0
    sets: dict[str, set[str]] = {}
    for name, reads_value in value.items():
        reads = parse_reads(reads_value, options, name=str(name))
        kmers = {x for r in reads for x in _canonical_kmers(r.seq, k)}
        if not kmers:
            raise InputFormatError(f"sample {name!r} has no usable {k}-mers.", code="READ_TOO_SHORT")
        sets[str(name)] = kmers
    names = list(sets)
    matrix: dict[str, dict[str, float]] = {a: {} for a in names}
    flagged = []
    for a in names:
        for b in names:
            if a == b:
                continue
            containment = len(sets[a] & sets[b]) / len(sets[a])
            matrix[a][b] = round6(containment) or 0.0
            if containment >= threshold:
                flagged.append({"sample": a, "found_in": b, "containment": round6(containment)})
    return {
        "samples": {n: len(s) for n, s in sets.items()},
        "kmer_size": k,
        "threshold": threshold,
        "containment": matrix,
        "flagged_pairs": sorted(flagged, key=lambda f: (-f["containment"], f["sample"], f["found_in"])),
        "definition": "containment[a][b] = |K(a) intersect K(b)| / |K(a)| over canonical k-mers. Heuristic: samples that are "
        "expected to be biologically unrelated should share almost no k-mers, so high containment suggests cross-"
        "contamination or index swapping. Related samples (same species, same amplicon) share k-mers legitimately.",
    }


# --------------------------------------------------------------------------- index hopping


def op_index_hopping(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    raw = options.get("expected_pairs", options.get("samples"))
    if raw is None:
        raise ConfigurationError(
            "options['samples'] is required: a mapping sample -> 'i7+i5' dual-index string.",
            code="MISSING_OPTION",
            hint="Example: options={'samples': {'S1': 'ACGTAC+TTGGCA', 'S2': 'GGTTAA+CCAATG'}}.",
        )
    samples = _parse_samples({"samples": raw})
    if not all("+" in code for code in samples.values()):
        raise ConfigurationError("index-hopping needs dual indexes written as 'i7+i5'.", code="INVALID_OPTION")
    i7s = {c.split("+")[0] for c in samples.values()}
    i5s = {c.split("+")[1] for c in samples.values()}
    expected = set(samples.values())
    counts: Counter[str] = Counter()
    for r in reads:
        token = header_index(r)
        if "+" not in token:
            raise InputFormatError(f"read {r.name!r} has a single index ({token}); dual indexes expected.", code="NO_HEADER_INDEX")
        counts[token] += 1
    n = len(reads)
    valid = sum(c for t, c in counts.items() if t in expected)
    hopped = {t: c for t, c in counts.items() if t not in expected and t.split("+")[0] in i7s and t.split("+")[1] in i5s}
    other = sum(c for t, c in counts.items() if t not in expected and t not in hopped)
    hopped_total = sum(hopped.values())
    return {
        "reads": n,
        "expected_combination_reads": valid,
        "index_hopped_reads": hopped_total,
        "unknown_index_reads": other,
        "index_hopping_rate": round6(hopped_total / (valid + hopped_total)) if valid + hopped_total else None,
        "unknown_index_fraction": round6(other / n),
        "hopped_combinations": [{"index": t, "reads": c} for t, c in sorted(hopped.items(), key=lambda kv: (-kv[1], kv[0]))[:20]],
        "reads_per_expected_sample": {name: counts.get(code, 0) for name, code in samples.items()},
        "definition": "A read is index-hopped when its i7 and i5 are each used in the sample sheet but the pair is not. "
        "rate = hopped / (hopped + expected-pair reads). Reads whose indexes are not in the sheet are reported separately; "
        "exact index matching, no mismatch tolerance.",
    }


# --------------------------------------------------------------------------- UMI


UMI_RE = re.compile(r"^[ACGTN]+$")


def _umi_from_read(read: Rd, options: Mapping[str, object]) -> str:
    sep = options.get("umi_separator", "_")
    if not isinstance(sep, str) or not sep:
        raise ConfigurationError("umi_separator must be a non-empty string.", code="INVALID_OPTION")
    token = read.name.rsplit(sep, 1)[-1] if sep in read.name else ""
    if not UMI_RE.match(token):
        raise InputFormatError(
            f"read {read.name!r} has no UMI after the last {sep!r} in its name.",
            code="NO_UMI",
            hint="Name reads like 'READ_ACGTACGT' (umi_tools convention), or run extract_umi_sequences first.",
        )
    return token


def group_umis(counts: Mapping[str, int], threshold: int, method: str) -> list[list[str]]:
    """Group UMIs by the UMI-tools 'cluster', 'adjacency' or 'directional' network (first member = representative)."""

    umis = sorted(counts, key=lambda u: (-counts[u], u))
    lengths = {len(u) for u in umis}
    if len(lengths) > 1:
        raise InputFormatError("all UMIs must have the same length.", code="UMI_LENGTH_MISMATCH")
    adjacency: dict[str, list[str]] = {u: [] for u in umis}
    for i, a in enumerate(umis):
        for b in umis[i + 1 :]:
            if sum(x != y for x, y in zip(a, b, strict=False)) <= threshold:
                if method == "directional":
                    if counts[a] >= 2 * counts[b] - 1:
                        adjacency[a].append(b)
                    if counts[b] >= 2 * counts[a] - 1:
                        adjacency[b].append(a)
                else:
                    adjacency[a].append(b)
                    adjacency[b].append(a)
    groups: list[list[str]] = []
    seen: set[str] = set()
    if method == "adjacency":
        components = []
        for u in umis:
            if u in seen:
                continue
            stack, comp = [u], []
            seen.add(u)
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adjacency[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            components.append(sorted(comp, key=lambda z: (-counts[z], z)))
        for comp in components:
            observed: set[str] = set()
            for u in comp:
                if u in observed:
                    continue
                group = [u] + [x for x in adjacency[u] if x not in observed and x != u]
                observed.update(group)
                groups.append(group)
        return groups
    for u in umis:  # cluster (undirected) and directional (BFS over directed edges)
        if u in seen:
            continue
        queue, comp = [u], []
        seen.add(u)
        while queue:
            x = queue.pop(0)
            comp.append(x)
            for y in sorted(adjacency[x], key=lambda z: (-counts[z], z)):
                if y not in seen:
                    seen.add(y)
                    queue.append(y)
        groups.append(comp)
    return groups


def _umi_counts(value: object, options: Mapping[str, object]) -> Counter[str]:
    if isinstance(value, Mapping) and value and all(isinstance(v, int) and not isinstance(v, bool) for v in value.values()) and not any(
        k in value for k in ("reads", "records")
    ):
        counts = Counter({str(k): int(v) for k, v in value.items()})
    elif isinstance(value, Mapping) and "umis" in value:
        return _umi_counts(value["umis"], options)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value and all(
        isinstance(v, str) and UMI_RE.match(v.upper()) for v in value
    ):
        counts = Counter(v.upper() for v in value)
    else:
        reads = parse_reads(value, options)
        counts = Counter(_umi_from_read(r, options) for r in reads)
    if not counts:
        raise InputFormatError("no UMIs found.", code="EMPTY_INPUT")
    for umi, c in counts.items():
        if not UMI_RE.match(umi) or c < 1:
            raise InputFormatError(f"invalid UMI or count: {umi!r}: {c}.", code="INVALID_UMI")
    if len(counts) > 5_000:
        raise InputFormatError("more than 5000 distinct UMIs (the network is quadratic); reduce the input.", code="INPUT_TOO_LARGE")
    return counts


def op_correct_umi_errors(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    counts = _umi_counts(value, options)
    method = opt_choice(options, "method", "directional", ("directional", "adjacency", "cluster"))
    threshold = opt_int(options, "threshold", 1, minimum=0, maximum=4)
    groups = group_umis(counts, threshold, method)
    mapping = {member: group[0] for group in groups for member in group}
    return {
        "method": method,
        "threshold": threshold,
        "distinct_umis_before": len(counts),
        "distinct_umis_after": len(groups),
        "reads": sum(counts.values()),
        "groups": [
            {"representative": g[0], "members": [*g[1:]], "reads": sum(counts[m] for m in g), "member_counts": {m: counts[m] for m in g}}
            for g in groups
        ],
        "corrections": {u: r for u, r in mapping.items() if u != r},
        "definition": "UMI-tools networks over Hamming distance <= threshold. directional: edge a->b when count(a) >= 2*count(b)-1, "
        "groups grown by breadth-first search from the most abundant UMI; adjacency: most abundant UMI plus its neighbours, repeated "
        "within each connected component; cluster: connected components. The representative is the most abundant UMI of a group.",
    }


def _consensus(members: Sequence[Rd]) -> tuple[str, str | None]:
    n = len(members)
    maxlen = max(len(m.seq) for m in members)
    seq, qual = [], []
    for i in range(maxlen):
        covered = [m for m in members if len(m.seq) > i]
        if len(covered) * 2 < n:
            break
        weight: Counter[str] = Counter()
        for m in covered:
            weight[m.seq[i]] += m.qual[i] if m.qual else 1
        best = max(weight, key=lambda b: (weight[b], -ord(b)))
        seq.append(best)
        if all(m.qual for m in covered):
            agree = sum(m.qual[i] for m in covered if m.seq[i] == best)  # type: ignore[index]
            disagree = sum(m.qual[i] for m in covered if m.seq[i] != best)  # type: ignore[index]
            qual.append(chr(33 + max(0, min(MAX_PHRED, agree - disagree))))
    return "".join(seq), ("".join(qual) if len(qual) == len(seq) and qual else None)


def op_collapse_umi_reads(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    group_by = opt_choice(options, "group_by", "umi_and_start", ("umi", "umi_and_start"))
    start_len = opt_int(options, "start_length", 10, minimum=1, maximum=1000)
    correct = opt_bool(options, "correct_errors", False)
    threshold = opt_int(options, "threshold", 1, minimum=0, maximum=4)
    umis = [_umi_from_read(r, options) for r in reads]
    mapping: dict[str, str] = {}
    if correct:
        counts = Counter(umis)
        mapping = {m: g[0] for g in group_umis(counts, threshold, "directional") for m in g}
    groups: dict[tuple[str, str], list[Rd]] = {}
    for r, u in zip(reads, umis, strict=False):
        umi = mapping.get(u, u)
        key = (umi, r.seq[:start_len] if group_by == "umi_and_start" else "")
        groups.setdefault(key, []).append(r)
    rows = []
    for (umi, start), members in sorted(groups.items()):
        seq, qual = _consensus(members)
        rows.append(
            {"umi": umi, "start": start or None, "reads": len(members), "read_names": [m.name for m in members][:20],
             "consensus_sequence": seq, "consensus_quality": qual}
        )
    return {
        "reads_in": len(reads),
        "molecules_out": len(rows),
        "collapse_ratio": round6(len(reads) / len(rows)),
        "group_by": group_by,
        "start_length": start_len if group_by == "umi_and_start" else None,
        "error_corrected": correct,
        "groups": rows,
        "definition": "Reads are grouped by (UMI, first start_length bases) because no alignment position is available; "
        "the consensus is a per-position vote weighted by Phred quality over positions covered by at least half of the group. "
        "Consensus quality = clip(sum Q of agreeing - sum Q of disagreeing, 0, 93).",
    }


def op_extract_umi(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    pattern = options.get("pattern")
    sep = options.get("umi_separator", "_")
    if not isinstance(sep, str) or not sep:
        raise ConfigurationError("umi_separator must be a non-empty string.", code="INVALID_OPTION")
    append = opt_bool(options, "append_to_name", True)
    location = opt_choice(options, "umi_location", "5prime", ("5prime", "3prime"))
    if pattern is not None:
        if not isinstance(pattern, str) or not re.fullmatch(r"[NXC]+", pattern):
            raise ConfigurationError(
                "pattern must consist of N (UMI base), C (cell-barcode base) and X (discarded base).",
                code="INVALID_OPTION", hint="Example: pattern='CCCCCCNNNNNNXX'.",
            )
        location = "5prime"
        length = len(pattern)
    else:
        length = opt_int(options, "umi_length", 8, minimum=1, maximum=100)
    out = []
    for r in reads:
        if len(r.seq) <= length:
            raise InputFormatError(f"read {r.name!r} ({len(r.seq)} bases) is not longer than the UMI region ({length}).", code="READ_TOO_SHORT")
        if location == "5prime":
            region, rest = r.seq[:length], r.seq[length:]
            rq = None if r.qual is None else r.qual[length:]
        else:
            region, rest = r.seq[-length:], r.seq[:-length]
            rq = None if r.qual is None else r.qual[:-length]
        if pattern:
            umi = "".join(b for b, p in zip(region, pattern, strict=False) if p == "N")
            cell = "".join(b for b, p in zip(region, pattern, strict=False) if p == "C")
        else:
            umi, cell = region, ""
        name = f"{r.name}{sep}{umi}" if append else r.name
        out.append(
            {"id": name, "umi": umi, "cell_barcode": cell or None, "sequence": rest,
             "quality": None if rq is None else "".join(chr(q + 33) for q in rq)}
        )
    counts = Counter(o["umi"] for o in out)
    return {
        "reads": len(out),
        "umi_location": location,
        "umi_region_length": length,
        "distinct_umis": len(counts),
        "records": out,
        "definition": "The first (5prime) or last (3prime) umi_length bases, or the N/C/X pattern positions at the 5' end, are "
        "removed from each read and appended to its name as '<name><umi_separator><UMI>' (umi_tools convention).",
    }


# --------------------------------------------------------------------------- error model with reference


def _seed_index(ref: str, k: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i in range(len(ref) - k + 1):
        index.setdefault(ref[i : i + k], []).append(i)
    return index


def op_read_error_model(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    reference = options.get("reference")
    if reference is None:
        need_quality(reads, "a quality-based error model (or pass options['reference'] for an observed model)")
        model = _error_profile(reads)
        model["note"] = "No reference supplied: the model is derived from reported Phred scores only. " + model["note"]
        return {"model": model}
    if not isinstance(reference, str):
        raise ConfigurationError("options['reference'] must be a sequence string.", code="INVALID_OPTION")
    ref_text = as_text(reference) or reference
    if ref_text.lstrip().startswith(">"):
        ref_text = "".join(r.seq for r in parse_fasta_text(ref_text, 1_000_000))
    ref = "".join(ref_text.split()).upper()
    if not 20 <= len(ref) <= 5_000_000:
        raise ConfigurationError("reference must be 20 to 5,000,000 bases.", code="INVALID_OPTION")
    if len(reads) > 200_000:
        raise InputFormatError("at most 200000 reads are aligned for the error model.", code="INPUT_TOO_LARGE")
    k = opt_int(options, "seed_length", 11, minimum=6, maximum=32)
    max_rate = opt_num(options, "max_mismatch_fraction", 0.2, minimum=0.0, maximum=0.5) or 0.0
    index = _seed_index(ref, k)
    mismatches_at: Counter[int] = Counter()
    bases_at: Counter[int] = Counter()
    substitutions: Counter[tuple[str, str]] = Counter()
    by_q_err: Counter[int] = Counter()
    by_q_tot: Counter[int] = Counter()
    aligned = 0
    for r in reads:
        best = None
        for strand, seq in (("+", r.seq), ("-", revcomp(r.seq))):
            if len(seq) < k:
                continue
            votes: Counter[int] = Counter()
            for s in range(0, len(seq) - k + 1, max(1, k // 2)):
                for pos in index.get(seq[s : s + k], ()):
                    votes[pos - s] += 1
            for diag, _ in votes.most_common(3):
                if diag < 0 or diag + len(seq) > len(ref):
                    continue
                mism = sum(1 for a, b in zip(seq, ref[diag : diag + len(seq)], strict=False) if a != b and a != "N")
                if best is None or mism < best[0]:
                    best = (mism, strand, diag, seq)
        if best is None or best[0] > max_rate * len(r.seq):
            continue
        _, strand, diag, seq = best
        aligned += 1
        for j, (a, b) in enumerate(zip(seq, ref[diag : diag + len(seq)], strict=False)):
            rj = j if strand == "+" else len(seq) - 1 - j  # position in the original read
            if a == "N":
                continue
            bases_at[rj] += 1
            wrong = a != b
            if wrong:
                mismatches_at[rj] += 1
                substitutions[(b, a)] += 1
            if r.qual:
                by_q_tot[r.qual[rj]] += 1
                by_q_err[r.qual[rj]] += wrong
    if aligned == 0:
        raise InputFormatError(
            "no read could be placed on the reference (ungapped seed-and-extend).",
            code="NO_ALIGNMENTS", hint="Check the reference, the seed_length and max_mismatch_fraction.",
        )
    total = sum(bases_at.values())
    errors = sum(mismatches_at.values())
    return {
        "model": {
            "source": "ungapped_reference_alignment",
            "reads": len(reads),
            "reads_aligned": aligned,
            "bases_compared": total,
            "mismatches": errors,
            "overall_mismatch_rate": round6(errors / total),
            "per_position_mismatch_rate": [
                {"position": p + 1, "bases": bases_at[p], "mismatch_rate": round6(mismatches_at[p] / bases_at[p])}
                for p in sorted(bases_at)
            ],
            "substitution_counts": {f"{a}>{b}": c for (a, b), c in sorted(substitutions.items())},
            "empirical_error_by_reported_q": {
                str(q): {"bases": by_q_tot[q], "empirical_error_rate": round6(by_q_err[q] / by_q_tot[q]),
                         "reported_error_rate": round6(q_to_p(q))}
                for q in sorted(by_q_tot)
            },
            "definition": "Each read is placed on the reference by k-mer seed voting (both strands, top 3 diagonals) and "
            "scored by Hamming distance without gaps; reads above max_mismatch_fraction are not counted. Insertions and "
            "deletions are not modelled, so indel-rich reads are excluded or inflate the substitution rate.",
        }
    }


# --------------------------------------------------------------------------- platform


_PLATFORM_RULES: list[tuple[str, re.Pattern[str]]] = [
    ("oxford_nanopore", re.compile(r"(^|\s)(runid|ch|start_time|read)=|^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")),
    ("pacbio", re.compile(r"^m\d+[A-Za-z0-9]*_\d{6}_\d{6}/\d+/|^[A-Za-z0-9_.\-]+/\d+/(ccs|\d+_\d+)(/|$)")),
    ("illumina", re.compile(r"^[A-Za-z0-9_\-]+:\d+:[A-Za-z0-9\-]+:\d+:\d+:\d+:\d+")),
    ("mgi_dnbseq", re.compile(r"^[A-Z]{1,3}\d{8,12}[A-Z]?L\d+C\d+R\d+")),
    ("ion_torrent", re.compile(r"^[A-Z0-9]{5}:\d+:\d+(\s|/|$)")),
]
_SRA = re.compile(r"^[SED]RR\d+\.\d+")


def _platform(reads: list[Rd]) -> dict[str, Any]:
    n = len(reads)
    votes: Counter[str] = Counter()
    for r in reads:
        header = f"{r.name} {r.comment}".strip()
        for platform, pattern in _PLATFORM_RULES:
            if pattern.search(header):
                votes[platform] += 1
                break
        else:
            votes["sra_renamed" if _SRA.match(r.name) else "unknown"] += 1
    best, best_n = max(votes.items(), key=lambda kv: (kv[1], kv[0]))
    lengths = sorted(len(r.seq) for r in reads)
    evidence: dict[str, Any] = {}
    first = reads[0]
    if best == "illumina":
        parts = first.name.split(":")
        evidence = {"instrument_id": parts[0], "run_number": parts[1], "flowcell": parts[2], "lane": parts[3]}
    uniform = lengths[0] == lengths[-1]
    return {
        "platform": best if best_n / n >= 0.5 else "unknown",
        "confidence_fraction": round6(best_n / n),
        "read_name_votes": dict(votes),
        "read_length": {"min": lengths[0], "median": round6(percentile(lengths, 0.5)), "max": lengths[-1], "uniform": uniform},
        "evidence": evidence,
        "definition": "Rule-based: the read-name (header) format of each read is matched against known layouts (Illumina "
        "instrument:run:flowcell:lane:tile:x:y, Oxford Nanopore runid=/ch=/UUID, PacBio movie/zmw/ccs, MGI/BGI "
        "L#C#R#, Ion Torrent). The platform is the majority layout when it covers >= 50% of reads. SRA-renamed headers "
        "(SRRxxxx.n) carry no platform information. Instrument model is not inferred.",
    }


def op_detect_platform(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return by_mate(value, options, _platform)


def op_detect_platform_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    if second is None:
        raise InputFormatError("detect_platform_tenx needs read1 and read2.", code="PAIR_REQUIRED")
    result = _platform(second)
    lengths = Counter(len(r.seq) for r in first)
    modal, _ = max(lengths.items(), key=lambda kv: (kv[1], kv[0]))
    layouts = {
        24: ("14 bp barcode + 10 bp UMI", "Chromium 3' v1"),
        26: ("16 bp barcode + 10 bp UMI", "Chromium 3' v2 / 5' v1-v2"),
        28: ("16 bp barcode + 12 bp UMI", "Chromium 3' v3 / v3.1 / 5' HT"),
    }
    layout = layouts.get(modal)
    result["tenx"] = {
        "read1_modal_length": modal,
        "layout": layout[0] if layout else None,
        "compatible_chemistries": layout[1] if layout else None,
        "note": "Chemistry is inferred only from the modal read-1 length (24/26/28). Cell Ranger's whitelist-based "
        "autodetection needs the 10x barcode whitelists and is not performed." + ("" if layout else " Read-1 length matches no known layout."),
    }
    return result


# --------------------------------------------------------------------------- 10x wrappers


def op_tenx_per_region(
    value: object, options: Mapping[str, object], fn: Callable[[list[Rd]], dict[str, Any]]
) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    return {**meta, "regions": _per_region(regions, fn)}


def op_quality_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return op_tenx_per_region(value, options, _quality)


def op_basecalling_qc_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    min_q = opt_num(options, "min_mean_quality", 20.0, minimum=0, maximum=MAX_PHRED) or 0.0
    return op_tenx_per_region(value, options, lambda reads: _basecalling(reads, min_q))


def op_error_profile_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return op_tenx_per_region(value, options, _error_profile)


def op_gc_bias_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    ref = opt_num(options, "reference_gc", None, minimum=0.0, maximum=1.0)
    return {**meta, "region": "rna_read (read 2)", **_gc_bias(regions["rna_read"], ref)}


def op_read_length_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    if second is None:
        raise InputFormatError("10x functions need a read pair: {'read1': ..., 'read2': ...}.", code="PAIR_REQUIRED")
    return {"read_pairs": len(first), "read1": _length_stats(first, options), "read2": _length_stats(second, options)}


# --------------------------------------------------------------------------- thin variants and demultiplexing


def op_adapter_content_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return op_adapter_content(value, options, warnings, tenx=True)


def op_contamination_tenx(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    return op_contamination(value, options, warnings, tenx=True)


def op_demultiplex(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    reads = parse_reads(value, options)
    samples = _parse_samples(options)
    if "unassigned" in samples:
        raise ConfigurationError("'unassigned' is reserved and cannot be a sample name.", code="INVALID_OPTION")
    rows, summary = assign_reads(reads, options)
    location = summary["barcode_location"]
    trim = opt_bool(options, "trim_barcode", True) and location != "header"
    limit = opt_int(options, "max_listed", 1000, minimum=1, maximum=1_000_000)
    groups: dict[str, list[dict[str, Any]]] = {name: [] for name in samples}
    groups["unassigned"] = []
    for read, row in zip(reads, rows, strict=False):
        key = row["sample"] or "unassigned"
        blen = len(row["barcode"])
        seq, qual = read.seq, read.qual
        if trim and key != "unassigned":
            seq, qual = (seq[blen:], None if qual is None else qual[blen:]) if location == "read_start" else (
                seq[:-blen], None if qual is None else qual[:-blen]
            )
        groups[key].append(
            {"id": read.name, "sequence": seq, "quality": None if qual is None else "".join(chr(q + 33) for q in qual)}
        )
    written: dict[str, str] = {}
    out_dir = options.get("output_directory")
    if out_dir is not None:
        directory = Path(str(out_dir))
        if not directory.is_dir():
            raise ConfigurationError(f"output_directory {directory} does not exist.", code="INVALID_OPTION")
        fastq = all(rec["quality"] is not None for recs in groups.values() for rec in recs)
        for name, recs in groups.items():
            target = directory / f"{name}.{'fastq' if fastq else 'fasta'}"
            if target.exists():
                raise ConfigurationError(f"{target} already exists.", code="OUTPUT_EXISTS", hint="Use an empty output_directory.")
            if fastq:
                text = "".join(f"@{r['id']}\n{r['sequence']}\n+\n{r['quality']}\n" for r in recs)
            else:
                text = "".join(f">{r['id']}\n{r['sequence']}\n" for r in recs)
            target.write_text(text, encoding="utf-8")
            written[name] = str(target)
    if any(len(v) > limit for v in groups.values()):
        warnings.append(f"returned reads per sample truncated to the first {limit} (see max_listed); counts are complete.")
    return {
        **{k: v for k, v in summary.items()},
        "barcode_trimmed": trim,
        "samples": {name: {"reads": len(recs), "records": recs[:limit]} for name, recs in groups.items()},
        "written_files": written or None,
    }
