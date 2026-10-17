"""Native analyses behind the ``variants`` candidate functions (statistics, genotype-level methods).

Pure standard library. Every function takes ``(input, options)`` and returns a JSON-friendly dict; a
``warnings`` key (list) is picked up by the caller. Invalid input raises a DNAKit exception.
"""

from __future__ import annotations

import itertools
import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from html import escape
from typing import Any

from ._impl_vcf import (
    VariantSet,
    Reference,
    call_ad,
    call_gt,
    chrom_key,
    diploid_class,
    fail,
    hwe_exact_pvalue,
    is_called,
    is_symbolic,
    load_variant_set,
    number_or_none,
    option_bool,
    option_choice,
    option_int,
    option_number,
    parse_gt,
    record_view,
    reference_from_options,
    split_multiallelic_set,
    to_vcf_text,
    trim_alleles,
    write_output,
)


def site_name(rec: Mapping[str, Any]) -> str:
    """Stable site label ``chrom:pos:ref:alt1,alt2``."""

    return f"{rec['chrom']}:{rec['pos']}:{rec['ref']}:{','.join(rec['alts'])}"


def pick_sample(vs: VariantSet, options: Mapping[str, object], key: str = "sample") -> int:
    """Index of the sample named by ``options[key]`` (default: the first sample)."""

    if not vs.samples:
        fail("the input has no sample genotype columns.", "VARIANTS_NO_SAMPLES", "Provide a VCF with sample columns or records with genotypes.")
    name = options.get(key)
    if name is None:
        return 0
    if name not in vs.samples:
        fail(f"sample {name!r} not found; available: {vs.samples[:10]}.", "VARIANTS_UNKNOWN_SAMPLE", f"Set options.{key} to one of the sample names.", config=True)
    return vs.samples.index(str(name))


# ---------------------------------------------------------------------------
# binomial tails
# ---------------------------------------------------------------------------
def _log_pmf(i: int, n: int, p: float) -> float:
    return (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
            + i * math.log(p) + (n - i) * math.log1p(-p))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""

    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    total = 0.0
    for i in range(k, n + 1):
        term = math.exp(_log_pmf(i, n, p))
        total += term
        if i > n * p and term < total * 1e-17:
            break
    return min(1.0, total)


def binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p)."""

    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0
    total = 0.0
    for i in range(k, -1, -1):
        term = math.exp(_log_pmf(i, n, p))
        total += term
        if i < n * p and term < total * 1e-17:
            break
    return min(1.0, total)


# ---------------------------------------------------------------------------
# genotype likelihood
# ---------------------------------------------------------------------------
def _likelihood(ref_reads: Sequence[float], alt_reads: Sequence[float], ploidy: int) -> dict[str, Any]:
    """Genotype likelihoods from per-read error probabilities (ref_reads/alt_reads hold the error rate of each read)."""

    logs = []
    for k in range(ploidy + 1):
        frac = k / ploidy
        total = 0.0
        for e in alt_reads:  # read supports ALT
            p = frac * (1 - e) + (1 - frac) * e
            total += math.log10(p) if p > 0 else -math.inf
        for e in ref_reads:  # read supports REF
            p = frac * e + (1 - frac) * (1 - e)
            total += math.log10(p) if p > 0 else -math.inf
        logs.append(total)
    best = max(logs)
    pl = [None if math.isinf(v) else int(round(-10 * (v - best))) for v in logs]
    finite = sorted(v for v in pl if v is not None)
    labels = ["/".join(["0"] * (ploidy - k) + ["1"] * k) if ploidy > 1 else str(k) for k in range(ploidy + 1)]
    best_index = logs.index(best)
    return {
        "genotypes": labels, "log10_likelihoods": [None if math.isinf(v) else v for v in logs],
        "phred_scaled": pl, "best_genotype": labels[best_index],
        "genotype_quality": min(99, finite[1]) if len(finite) > 1 else 99,
        "ref_reads": len(ref_reads), "alt_reads": len(alt_reads),
    }


def genotype_likelihood(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    error = option_number(options, "error_rate", 0.01, minimum=1e-9, maximum=0.5)
    ploidy = option_int(options, "ploidy", 2, minimum=1, maximum=6)
    assert error is not None and ploidy is not None

    def from_counts(item: Mapping[str, Any]) -> dict[str, Any]:
        if "reads" in item:
            reads = item["reads"]
            if not isinstance(reads, Sequence) or isinstance(reads, (str, bytes)) or len(reads) > 1_000_000:
                fail("reads must be a list.", "VARIANTS_BAD_INPUT", "Use [{'allele': 'alt', 'base_quality': 30}, ...].")
            ref_e: list[float] = []
            alt_e: list[float] = []
            for read in reads:
                if not isinstance(read, Mapping):
                    fail("each read must be a mapping.", "VARIANTS_BAD_INPUT", "Use {'allele': 'ref'|'alt', 'base_quality': q}.")
                allele = read.get("allele")
                if allele not in {"ref", "alt", 0, 1, "REF", "ALT"}:
                    fail(f"read allele {allele!r} must be 'ref' or 'alt'.", "VARIANTS_BAD_INPUT", "Use 'ref' or 'alt'.")
                quality = number_or_none(read.get("base_quality"))
                e = error if quality is None else min(0.5, 10 ** (-quality / 10))
                (alt_e if allele in {"alt", 1, "ALT"} else ref_e).append(max(e, 1e-12))
            return _likelihood(ref_e, alt_e, ploidy)
        r, a = item.get("ref_count"), item.get("alt_count")
        for name, v in (("ref_count", r), ("alt_count", a)):
            if isinstance(v, bool) or not isinstance(v, int) or v < 0 or v > 10_000_000:
                fail(f"{name} must be a non-negative integer.", "VARIANTS_BAD_INPUT", "Pass ref_count and alt_count (read counts).")
        return _likelihood([error] * r, [error] * a, ploidy)  # type: ignore[operator]

    def counts_form(item: object) -> bool:
        return isinstance(item, Mapping) and any(k in item for k in ("ref_count", "alt_count", "reads"))

    if counts_form(value):
        return {"model": _MODEL_NOTE, "error_rate": error, "ploidy": ploidy, "result": from_counts(value)}  # type: ignore[arg-type]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value and all(counts_form(x) for x in value):
        return {"model": _MODEL_NOTE, "error_rate": error, "ploidy": ploidy, "results": [from_counts(x) for x in value]}
    vs = split_multiallelic_set(load_variant_set(value))
    rows, warnings = [], []
    for rec in vs.records:
        for name, call in zip(vs.samples, rec["calls"]):
            ad = call_ad(call)
            if ad is None or len(ad) != 2:
                warnings.append(f"{rec['chrom']}:{rec['pos']} {name}: no usable AD, skipped")
                continue
            row = _likelihood([error] * ad[0], [error] * ad[1], ploidy)
            row.update({"site": site_name(rec), "sample": name})
            rows.append(row)
    if not rows:
        fail("no AD-bearing genotype calls found.", "VARIANTS_BAD_INPUT", "Provide AD (allelic depth) per sample, or ref_count/alt_count.")
    return {"model": _MODEL_NOTE, "error_rate": error, "ploidy": ploidy, "results": rows, "warnings": warnings[:50]}


_MODEL_NOTE = ("each read independently supports ALT with probability (k/m)(1-e)+(1-k/m)e for a genotype with k ALT copies "
               "of m; likelihoods are log10, PL = -10*(log10 L - max log10 L), GQ = second-smallest PL capped at 99")


# ---------------------------------------------------------------------------
# HWE
# ---------------------------------------------------------------------------
def hwe_pvalue(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    mid_p = option_bool(options, "mid_p", False)

    def stats(hom_ref: int, het: int, hom_alt: int) -> dict[str, Any]:
        n = hom_ref + het + hom_alt
        p = (2 * hom_ref + het) / (2 * n) if n else None
        expected_het = 2 * p * (1 - p) * n if p is not None else None
        return {
            "n_hom_ref": hom_ref, "n_het": het, "n_hom_alt": hom_alt, "n_samples": n,
            "observed_het_fraction": het / n if n else None, "expected_het_fraction": expected_het / n if n else None,
            "ref_allele_frequency": p, "p_value": hwe_exact_pvalue(hom_ref, het, hom_alt, mid_p=mid_p),
        }

    counts = value
    if isinstance(value, Mapping) and "counts" in value:
        counts = value["counts"]
    if isinstance(counts, Mapping) and {"hom_ref", "het", "hom_alt"} <= set(counts):
        counts = [counts["hom_ref"], counts["het"], counts["hom_alt"]]
    if isinstance(counts, Sequence) and not isinstance(counts, (str, bytes)) and len(counts) == 3 and all(
            isinstance(c, int) and not isinstance(c, bool) for c in counts):
        return {"method": _HWE_NOTE, "mid_p": mid_p, **stats(*counts)}
    vs = load_variant_set(value)
    if not vs.samples:
        fail("HWE needs genotypes.", "VARIANTS_NO_SAMPLES", "Provide sample genotypes or three genotype counts [AA, Aa, aa].")
    rows, warnings = [], []
    for rec in vs.records:
        if len(rec["alts"]) != 1:
            warnings.append(f"{site_name(rec)}: multi-allelic site skipped (split it first)")
            continue
        counts_map = {"hom_ref": 0, "het": 0, "hom_alt": 0}
        other = 0
        for call in rec["calls"]:
            cls = diploid_class(call_gt(call)[0])
            if cls in counts_map:
                counts_map[cls] += 1
            elif cls == "other":
                other += 1
        if other:
            warnings.append(f"{site_name(rec)}: {other} non-diploid calls ignored")
        row = stats(counts_map["hom_ref"], counts_map["het"], counts_map["hom_alt"])
        row["site"] = site_name(rec)
        rows.append(row)
    return {"method": _HWE_NOTE, "mid_p": mid_p, "sites": rows, "warnings": warnings[:50]}


_HWE_NOTE = "exact test of Hardy-Weinberg equilibrium (Wigginton, Cutler & Abecasis 2005), two-sided"


# ---------------------------------------------------------------------------
# LD
# ---------------------------------------------------------------------------
def ld_matrix(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    statistic = option_choice(options, "statistic", "r2", ("r2", "r", "d_prime"))
    cap = option_int(options, "max_variants", 200, minimum=2, maximum=1000)
    assert cap is not None
    vs = load_variant_set(value)
    if not vs.samples:
        fail("LD needs sample genotypes.", "VARIANTS_NO_SAMPLES", "Provide at least two samples with genotypes.")
    warnings: list[str] = []
    recs = [r for r in vs.records if len(r["alts"]) == 1]
    if len(recs) != len(vs.records):
        warnings.append(f"{len(vs.records) - len(recs)} multi-allelic sites skipped")
    if len(recs) < 2:
        fail("LD needs at least two biallelic variants.", "VARIANTS_BAD_INPUT", "Provide two or more biallelic sites.")
    if len(recs) > cap:
        fail(f"{len(recs)} variants exceed options.max_variants={cap}.", "VARIANTS_TOO_MANY_RECORDS", "Subset the region or raise max_variants (<=1000).", config=True)
    parsed = [[call_gt(c) for c in r["calls"]] for r in recs]
    all_phased = all(ph for row in parsed for al, ph in row if is_called(al) and len(al) == 2 and al[0] != al[1])
    phased_mode = option_choice(options, "phase", "auto", ("auto", "phased", "unphased"))
    use_haps = (phased_mode == "phased") or (phased_mode == "auto" and all_phased)
    if statistic == "d_prime" and not use_haps:
        fail("D' needs phased genotypes.", "VARIANTS_NEEDS_PHASE", "Provide phased GT (0|1) or choose statistic='r2'/'r'.", config=True)

    def dosage(row: list[tuple[list[int | None], bool]], i: int) -> int | None:
        al = row[i][0]
        if len(al) != 2 or not is_called(al):
            return None
        return sum(1 for a in al if a)

    size = len(recs)
    matrix: list[list[float | None]] = [[None] * size for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            val: float | None
            if use_haps:
                pairs = []
                for s in range(len(vs.samples)):
                    ai, aj = parsed[i][s][0], parsed[j][s][0]
                    if len(ai) == 2 and len(aj) == 2 and is_called(ai) and is_called(aj):
                        pairs.append((ai[0] != 0, aj[0] != 0))
                        pairs.append((ai[1] != 0, aj[1] != 0))
                n = len(pairs)
                if n < 2:
                    val = None
                else:
                    pa = sum(x for x, _ in pairs) / n
                    pb = sum(y for _, y in pairs) / n
                    pab = sum(x and y for x, y in pairs) / n
                    d = pab - pa * pb
                    denom = pa * (1 - pa) * pb * (1 - pb)
                    if statistic == "d_prime":
                        dmax = min(pa * (1 - pb), (1 - pa) * pb) if d >= 0 else min(pa * pb, (1 - pa) * (1 - pb))
                        val = abs(d) / dmax if dmax > 0 else None
                    elif denom <= 0:
                        val = None
                    else:
                        r = d / math.sqrt(denom)
                        val = r * r if statistic == "r2" else r
            else:
                xs, ys = [], []
                for s in range(len(vs.samples)):
                    a, b = dosage(parsed[i], s), dosage(parsed[j], s)
                    if a is not None and b is not None:
                        xs.append(a)
                        ys.append(b)
                n = len(xs)
                if n < 2:
                    val = None
                else:
                    mx, my = sum(xs) / n, sum(ys) / n
                    sxx = sum((x - mx) ** 2 for x in xs)
                    syy = sum((y - my) ** 2 for y in ys)
                    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                    if sxx <= 0 or syy <= 0:
                        val = None
                    else:
                        r = sxy / math.sqrt(sxx * syy)
                        val = r * r if statistic == "r2" else r
            matrix[i][j] = matrix[j][i] = val
    return {
        "statistic": statistic,
        "method": "haplotype-based (phased)" if use_haps else "squared/plain Pearson correlation of ALT-allele dosages (unphased, Rogers-Huff style r)",
        "variants": [site_name(r) for r in recs], "matrix": matrix, "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# runs of homozygosity
# ---------------------------------------------------------------------------
def runs_of_homozygosity(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    min_snps = option_int(options, "min_snps", 25, minimum=1)
    min_len = option_int(options, "min_length_bp", 1_000_000, minimum=0)
    max_hets = option_int(options, "max_heterozygotes", 1, minimum=0)
    max_missing = option_int(options, "max_missing", 2, minimum=0)
    max_gap = option_int(options, "max_gap_bp", 1_000_000, minimum=1)
    assert None not in (min_snps, min_len, max_hets, max_missing, max_gap)
    vs = load_variant_set(value)
    idx = pick_sample(vs, options)
    by_chrom: dict[str, list[tuple[int, str]]] = {}
    for rec in vs.records:
        cls = diploid_class(call_gt(rec["calls"][idx])[0])
        if cls == "hom_alt" or cls == "hom_ref":
            kind = "hom"
        elif cls == "het":
            kind = "het"
        else:
            kind = "missing"
        by_chrom.setdefault(rec["chrom"], []).append((rec["pos"], kind))
    runs = []
    for chrom in sorted(by_chrom, key=chrom_key):
        sites = sorted(by_chrom[chrom])
        i, n = 0, len(sites)
        while i < n:
            if sites[i][1] != "hom":
                i += 1
                continue
            hets = miss = 0
            last_hom = i
            j = i + 1
            while j < n:
                if sites[j][0] - sites[j - 1][0] > max_gap:  # type: ignore[operator]
                    break
                kind = sites[j][1]
                if kind == "het":
                    if hets + 1 > max_hets:  # type: ignore[operator]
                        break
                    hets += 1
                elif kind == "missing":
                    if miss + 1 > max_missing:  # type: ignore[operator]
                        break
                    miss += 1
                else:
                    last_hom = j
                j += 1
            seg = sites[i:last_hom + 1]
            length = seg[-1][0] - seg[0][0] + 1
            if len(seg) >= min_snps and length >= min_len:  # type: ignore[operator]
                runs.append({"chrom": chrom, "start": seg[0][0], "end": seg[-1][0], "length_bp": length, "n_snps": len(seg),
                             "n_heterozygous": sum(1 for _, k in seg if k == "het"), "n_missing": sum(1 for _, k in seg if k == "missing")})
            i = max(j, last_hom + 1)
    total = sum(r["length_bp"] for r in runs)
    return {"sample": vs.samples[idx], "method": "greedy left-to-right runs of homozygous genotypes tolerating a bounded number of "
            "heterozygous/missing calls (consecutive-SNP definition, as in PLINK --homozyg-het style runs)",
            "parameters": {"min_snps": min_snps, "min_length_bp": min_len, "max_heterozygotes": max_hets, "max_missing": max_missing, "max_gap_bp": max_gap},
            "runs": runs, "n_runs": len(runs), "total_roh_bp": total}


# ---------------------------------------------------------------------------
# identity by descent / kinship (KING-robust)
# ---------------------------------------------------------------------------
def identity_by_descent(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    n = len(vs.samples)
    if n < 2:
        fail("IBD/kinship needs at least two samples.", "VARIANTS_NO_SAMPLES", "Provide a multi-sample VCF.")
    if n > 150:
        fail(f"{n} samples exceed the 150-sample limit.", "VARIANTS_TOO_MANY_SAMPLES", "Compute kinship for subsets.")
    if len(vs.records) * n * (n - 1) // 2 > 60_000_000:
        fail("input too large for pairwise kinship.", "VARIANTS_TOO_MANY_RECORDS", "Thin the variants or reduce the sample count.")
    classes: list[list[str]] = []
    for rec in vs.records:
        if len(rec["alts"]) != 1:
            continue
        classes.append([diploid_class(call_gt(c)[0]) for c in rec["calls"]])
    if not classes:
        fail("no biallelic sites available.", "VARIANTS_BAD_INPUT", "Provide biallelic diploid genotypes.")
    pairs = []
    for i, j in itertools.combinations(range(n), 2):
        het_het = opposite = het_i = het_j = used = 0
        for row in classes:
            a, b = row[i], row[j]
            if a not in {"hom_ref", "het", "hom_alt"} or b not in {"hom_ref", "het", "hom_alt"}:
                continue
            used += 1
            het_i += a == "het"
            het_j += b == "het"
            if a == "het" and b == "het":
                het_het += 1
            elif {a, b} == {"hom_ref", "hom_alt"}:
                opposite += 1
        denom = het_i + het_j
        kinship = (het_het - 2 * opposite) / denom if denom else None
        pairs.append({
            "sample_1": vs.samples[i], "sample_2": vs.samples[j], "n_sites": used, "n_het_het": het_het, "n_opposite_homozygous": opposite,
            "ibs0_fraction": opposite / used if used else None,
            "kinship": kinship, "pi_hat_approx": 2 * kinship if kinship is not None else None, "relationship": _degree(kinship),
        })
    return {"method": "KING-robust kinship estimator (Manichaikul et al. 2010): (N_AaAa - 2 N_AAaa) / (N_Aa(i) + N_Aa(j)) over sites "
            "called in both samples; pi_hat_approx = 2 * kinship", "pairs": pairs}


def _degree(kinship: float | None) -> str | None:
    if kinship is None:
        return None
    if kinship > 0.3535:
        return "duplicate_or_monozygotic_twin"
    if kinship > 0.1768:
        return "first_degree"
    if kinship > 0.0884:
        return "second_degree"
    if kinship > 0.0442:
        return "third_degree"
    return "unrelated"


# ---------------------------------------------------------------------------
# haplotype frequency (EM)
# ---------------------------------------------------------------------------
def haplotype_frequency(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    max_iter = option_int(options, "max_iterations", 1000, minimum=1, maximum=100_000)
    tol = option_number(options, "tolerance", 1e-10, minimum=1e-15, maximum=1e-2)
    assert max_iter is not None and tol is not None
    warnings: list[str] = []
    loci: list[str] | None = None
    phased_haps: list[tuple[int, ...]] | None = None
    individuals: list[list[tuple[int, int] | None]] = []

    def gt_pair(text: object) -> tuple[int, int] | None:
        al, _ = parse_gt(text)
        if len(al) != 2 or not is_called(al):
            return None
        return (al[0], al[1])  # type: ignore[return-value]

    if isinstance(value, Mapping) and "haplotypes" in value:
        haps = value["haplotypes"]
        if not isinstance(haps, Sequence) or isinstance(haps, (str, bytes)):
            fail("haplotypes must be a list.", "VARIANTS_BAD_INPUT", "Use ['010', '011'] or [[0,1,0], ...].")
        phased_haps = []
        for h in haps:
            alleles = tuple(int(c) for c in h) if isinstance(h, str) else tuple(int(c) for c in h)
            phased_haps.append(alleles)
    elif isinstance(value, Mapping) and "genotypes" in value and not any(k in value for k in ("records", "variants", "vcf")):
        rows = value["genotypes"]
        if not isinstance(rows, Sequence) or any(not isinstance(r, Sequence) or isinstance(r, (str, bytes)) for r in rows):
            fail("genotypes must be a list (individuals) of lists (loci) of GT strings.", "VARIANTS_BAD_INPUT", "Use [['0/1','0/0'], ['1/1','0/1']].")
        individuals = [[gt_pair(g) for g in r] for r in rows]
        loci = [f"locus_{k + 1}" for k in range(len(rows[0]))] if rows else []
    elif (isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value
          and all(isinstance(r, Sequence) and not isinstance(r, (str, bytes)) for r in value)):
        individuals = [[gt_pair(g) for g in r] for r in value]
        loci = [f"locus_{k + 1}" for k in range(len(value[0]))]
    else:
        vs = load_variant_set(value)
        if not vs.samples:
            fail("haplotype frequencies need genotypes.", "VARIANTS_NO_SAMPLES", "Provide sample genotypes at two or more loci.")
        recs = [r for r in vs.records if len(r["alts"]) == 1]
        if len(recs) < len(vs.records):
            warnings.append("multi-allelic sites skipped")
        loci = [site_name(r) for r in recs]
        parsed = [[call_gt(r["calls"][s]) for r in recs] for s in range(len(vs.samples))]
        everything_phased = all(
            ph or (len(al) == 2 and is_called(al) and al[0] == al[1]) for row in parsed for al, ph in row if is_called(al))
        if option_bool(options, "use_phase", True) and everything_phased and all(is_called(al) for row in parsed for al, _ in row):
            phased_haps = []
            for row in parsed:
                for h in (0, 1):
                    phased_haps.append(tuple(al[h] for al, _ in row))  # type: ignore[arg-type]
        else:
            individuals = [[gt_pair_from(al) for al, _ in row] for row in parsed]
    if phased_haps is not None:
        if not phased_haps or len({len(h) for h in phased_haps}) != 1:
            fail("haplotypes must be a non-empty list of equal length.", "VARIANTS_BAD_INPUT", "Give one allele per locus for every haplotype.")
        counts = Counter(phased_haps)
        total = sum(counts.values())
        out = [{"haplotype": "-".join(map(str, h)), "alleles": list(h), "count": c, "frequency": c / total}
               for h, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
        return {"method": "direct counting of phased haplotypes", "n_chromosomes": total, "loci": loci, "haplotypes": out, "warnings": warnings}
    usable = [row for row in individuals if row and all(g is not None for g in row)]
    skipped = len(individuals) - len(usable)
    if skipped:
        warnings.append(f"{skipped} individuals with missing genotypes were excluded")
    if not usable:
        fail("no individual has complete genotypes.", "VARIANTS_BAD_INPUT", "Provide at least one individual genotyped at every locus.")
    n_loci = len(usable[0])
    if any(len(r) != n_loci for r in usable):
        fail("individuals have different numbers of loci.", "VARIANTS_BAD_INPUT", "Give every individual a genotype at every locus.")
    if n_loci < 2 or n_loci > 12:
        fail("EM haplotype estimation needs between 2 and 12 loci.", "VARIANTS_BAD_INPUT", "Provide 2-12 loci.")
    resolutions: list[list[tuple[tuple[int, ...], tuple[int, ...]]]] = []
    for row in usable:
        het_idx = [k for k, g in enumerate(row) if g[0] != g[1]]  # type: ignore[index]
        if len(het_idx) > 10:
            fail("an individual has more than 10 heterozygous loci.", "VARIANTS_BAD_INPUT", "Reduce the number of loci per individual.")
        options_list = []
        for mask in range(1 << max(0, len(het_idx) - 1)):
            h1, h2 = [], []
            for k, g in enumerate(row):
                a, b = g  # type: ignore[misc]
                if a == b:
                    h1.append(a)
                    h2.append(b)
                else:
                    pos = het_idx.index(k)
                    flip = pos > 0 and (mask >> (pos - 1)) & 1
                    h1.append(b if flip else a)
                    h2.append(a if flip else b)
            options_list.append((tuple(h1), tuple(h2)))
        resolutions.append(options_list)
    haplotypes = sorted({h for opts in resolutions for pair in opts for h in pair})
    freq = {h: 1.0 / len(haplotypes) for h in haplotypes}
    iterations = 0
    for iterations in range(1, max_iter + 1):
        counts_acc = dict.fromkeys(haplotypes, 0.0)
        for opts in resolutions:
            weights = [freq[a] * freq[b] for a, b in opts]
            total = sum(weights)
            for (a, b), w in zip(opts, weights):
                share = w / total
                counts_acc[a] += share
                counts_acc[b] += share
        new = {h: c / (2 * len(usable)) for h, c in counts_acc.items()}
        delta = max(abs(new[h] - freq[h]) for h in haplotypes)
        freq = new
        if delta < tol:
            break
    out = [{"haplotype": "-".join(map(str, h)), "alleles": list(h), "count": f * 2 * len(usable), "frequency": f}
           for h, f in sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))]
    return {"method": "EM algorithm for haplotype frequencies from unphased genotypes (Excoffier & Slatkin 1995; "
            "Hardy-Weinberg assumed, uniform start)", "n_individuals": len(usable), "n_chromosomes": 2 * len(usable),
            "iterations": iterations, "loci": loci, "haplotypes": out, "warnings": warnings}


def gt_pair_from(alleles: Sequence[int | None]) -> tuple[int, int] | None:
    if len(alleles) != 2 or not is_called(alleles):
        return None
    return (alleles[0], alleles[1])  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# copy-number profile from depth
# ---------------------------------------------------------------------------
def copy_number_profile(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    purity = option_number(options, "purity", 1.0, minimum=0.01, maximum=1.0)
    ploidy = option_number(options, "ploidy", 2.0, minimum=1.0, maximum=8.0)
    assert purity is not None and ploidy is not None
    exclude = options.get("median_exclude_chroms", ["X", "Y", "M", "MT", "chrX", "chrY", "chrM", "chrMT"])
    if isinstance(exclude, str):
        exclude = [exclude]
    bins_in = value["bins"] if isinstance(value, Mapping) and "bins" in value else value
    if not isinstance(bins_in, Sequence) or isinstance(bins_in, (str, bytes)) or not bins_in or len(bins_in) > 1_000_000:
        fail("input must be a non-empty list of depth bins.", "VARIANTS_BAD_INPUT", "Use [{'chrom','start','end','depth'[,'control_depth']}, ...].")
    bins = []
    for index, item in enumerate(bins_in):
        if isinstance(item, Mapping):
            chrom, start, end = item.get("chrom"), item.get("start"), item.get("end")
            depth, control = item.get("depth"), item.get("control_depth", item.get("normal_depth"))
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes)) and len(item) in (4, 5):
            chrom, start, end, depth = item[:4]
            control = item[4] if len(item) == 5 else None
        else:
            fail(f"bin {index} is not a mapping or [chrom,start,end,depth].", "VARIANTS_BAD_INPUT", "Provide bin mappings.")
        if chrom in (None, "") or not isinstance(start, int) or not isinstance(end, int) or end < start:
            fail(f"bin {index} needs chrom and integer start <= end.", "VARIANTS_BAD_INPUT", "Provide chrom, start, end.")
        d = number_or_none(depth)
        c = number_or_none(control)
        if d is None or d < 0 or (control is not None and (c is None or c < 0)):
            fail(f"bin {index} needs a non-negative numeric depth.", "VARIANTS_BAD_INPUT", "Provide numeric depth (and control_depth if used).")
        bins.append((str(chrom), start, end, d, c))
    paired = [b[4] is not None for b in bins]
    if any(paired) and not all(paired):
        fail("control_depth must be given for all bins or none.", "VARIANTS_BAD_INPUT", "Provide control_depth consistently.")
    use_control = all(paired)
    keep = [b for b in bins if b[0] not in set(exclude) and b[3] > 0 and (not use_control or b[4] > 0)]  # type: ignore[operator]
    if not keep:
        fail("no positive-depth bins to normalise against.", "VARIANTS_BAD_INPUT", "Provide bins with positive depth.")

    def median(values: list[float]) -> float:
        s = sorted(values)
        m = len(s) // 2
        return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2

    med = median([b[3] for b in keep])
    cmed = median([b[4] for b in keep]) if use_control else None  # type: ignore[misc]
    rows = []
    scale = purity * ploidy + 2 * (1 - purity)
    for chrom, start, end, depth, control in bins:
        ratio = (depth / med) / ((control / cmed) if use_control and control else 1.0) if (not use_control or control) else None
        if ratio is None:
            cn = None
            log2 = None
        else:
            cn = max(0.0, (ratio * scale - 2 * (1 - purity)) / purity)
            log2 = math.log2(ratio) if ratio > 0 else None
        rows.append({"chrom": chrom, "start": start, "end": end, "depth": depth, "ratio": ratio, "log2_ratio": log2,
                     "copy_number": cn, "rounded_copy_number": None if cn is None else int(round(cn))})
    baseline = int(round(ploidy))
    segments: list[dict[str, Any]] = []
    for row in rows:
        state = row["rounded_copy_number"]
        last = segments[-1] if segments else None
        if last and last["chrom"] == row["chrom"] and last["copy_number"] == state and state is not None:
            last["end"] = row["end"]
            last["_ratios"].append(row["ratio"])
            last["n_bins"] += 1
        else:
            segments.append({"chrom": row["chrom"], "start": row["start"], "end": row["end"], "copy_number": state,
                             "n_bins": 1, "_ratios": [row["ratio"]]})
    for seg in segments:
        ratios = [r for r in seg.pop("_ratios") if r is not None]
        mean_ratio = sum(ratios) / len(ratios) if ratios else None
        seg["mean_ratio"] = mean_ratio
        seg["mean_log2_ratio"] = math.log2(mean_ratio) if mean_ratio and mean_ratio > 0 else None
        cn = seg["copy_number"]
        seg["state"] = (None if cn is None else "homozygous_deletion" if cn == 0 else "loss" if cn < baseline else
                        "neutral" if cn == baseline else "gain" if cn < baseline + 3 else "amplification")
    return {
        "method": "per-bin depth ratio to the median (or to the control depth ratio), copy number = "
                  "(ratio*(purity*ploidy+2*(1-purity)) - 2*(1-purity))/purity; adjacent bins with the same rounded "
                  "copy number are merged into segments (no statistical segmentation)",
        "purity": purity, "ploidy": ploidy, "median_depth": med, "control_median_depth": cmed,
        "bins": rows, "segments": segments,
    }


# ---------------------------------------------------------------------------
# low-frequency and mosaic variant detection from allelic depths
# ---------------------------------------------------------------------------
def _allelic_rows(value: object, options: Mapping[str, object]) -> tuple[VariantSet, list[dict[str, Any]], list[str]]:
    vs = split_multiallelic_set(load_variant_set(value))
    if not vs.samples:
        fail("allelic-depth analysis needs sample columns with AD.", "VARIANTS_NO_SAMPLES", "Provide FORMAT/AD for at least one sample.")
    wanted = options.get("samples")
    if wanted is None and options.get("sample") is not None:
        wanted = [options["sample"]]
    if isinstance(wanted, str):
        wanted = [wanted]
    names = list(vs.samples) if wanted is None else list(wanted)
    for name in names:
        if name not in vs.samples:
            fail(f"sample {name!r} not found.", "VARIANTS_UNKNOWN_SAMPLE", "Use names from the VCF header.", config=True)
    rows, warnings = [], []
    for rec in vs.records:
        for name in names:
            call = rec["calls"][vs.samples.index(name)]
            ad = call_ad(call)
            if ad is None or len(ad) != 2:
                continue
            depth = sum(ad)
            rows.append({"rec": rec, "sample": name, "alt": ad[1], "depth": depth, "ref": ad[0], "gt": call.get("GT")})
    if not rows:
        warnings.append("no calls with AD were found")
    return vs, rows, warnings


def low_frequency_variants(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    min_vaf = option_number(options, "min_vaf", 0.01, minimum=0, maximum=1)
    max_vaf = option_number(options, "max_vaf", 0.2, minimum=0, maximum=1)
    min_depth = option_int(options, "min_depth", 100, minimum=1)
    min_alt = option_int(options, "min_alt_count", 3, minimum=1)
    error = option_number(options, "error_rate", 0.001, minimum=1e-9, maximum=0.5)
    alpha = option_number(options, "alpha", 0.01, minimum=0, maximum=1)
    assert None not in (min_vaf, max_vaf, min_depth, min_alt, error, alpha)
    _, rows, warnings = _allelic_rows(value, options)
    hits, tested = [], 0
    for row in rows:
        if row["depth"] < min_depth:  # type: ignore[operator]
            continue
        tested += 1
        vaf = row["alt"] / row["depth"]
        p = binom_sf(row["alt"], row["depth"], error)  # type: ignore[arg-type]
        if min_vaf <= vaf <= max_vaf and row["alt"] >= min_alt and p <= alpha:  # type: ignore[operator]
            rec = row["rec"]
            hits.append({"site": site_name(rec), "chrom": rec["chrom"], "pos": rec["pos"], "ref": rec["ref"], "alt": rec["alts"][0],
                         "sample": row["sample"], "alt_count": row["alt"], "depth": row["depth"], "vaf": vaf, "p_value_vs_error": p})
    return {"method": "calls where min_vaf <= alt/depth <= max_vaf, alt count >= min_alt_count, depth >= min_depth and the one-sided "
            "binomial P(X >= alt | depth, error_rate) <= alpha (heuristic, no strand/context error model)",
            "parameters": {"min_vaf": min_vaf, "max_vaf": max_vaf, "min_depth": min_depth, "min_alt_count": min_alt, "error_rate": error, "alpha": alpha},
            "n_tested": tested, "n_detected": len(hits), "variants": hits, "warnings": warnings}


def mosaic_variants(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    min_vaf = option_number(options, "min_vaf", 0.03, minimum=0, maximum=1)
    max_vaf = option_number(options, "max_vaf", 0.4, minimum=0, maximum=1)
    min_depth = option_int(options, "min_depth", 20, minimum=1)
    min_alt = option_int(options, "min_alt_count", 3, minimum=1)
    error = option_number(options, "error_rate", 0.002, minimum=1e-9, maximum=0.5)
    alpha = option_number(options, "alpha", 0.01, minimum=0, maximum=1)
    assert None not in (min_vaf, max_vaf, min_depth, min_alt, error, alpha)
    _, rows, warnings = _allelic_rows(value, options)
    hits, tested = [], 0
    for row in rows:
        if row["depth"] < min_depth:  # type: ignore[operator]
            continue
        tested += 1
        vaf = row["alt"] / row["depth"]
        p_err = binom_sf(row["alt"], row["depth"], error)  # type: ignore[arg-type]
        p_het = binom_cdf(row["alt"], row["depth"], 0.5)  # type: ignore[arg-type]
        if min_vaf <= vaf <= max_vaf and row["alt"] >= min_alt and p_err <= alpha and p_het <= alpha:  # type: ignore[operator]
            rec = row["rec"]
            hits.append({"site": site_name(rec), "chrom": rec["chrom"], "pos": rec["pos"], "ref": rec["ref"], "alt": rec["alts"][0],
                         "sample": row["sample"], "alt_count": row["alt"], "depth": row["depth"], "vaf": vaf,
                         "p_value_vs_error": p_err, "p_value_vs_heterozygous": p_het})
    return {"method": "mosaic candidates: ALT fraction in [min_vaf, max_vaf], enough ALT reads, and the ALT count is both significantly "
            "above sequencing error (one-sided binomial upper tail) and significantly below a germline heterozygous 0.5 "
            "(one-sided binomial lower tail), each at level alpha (heuristic)",
            "parameters": {"min_vaf": min_vaf, "max_vaf": max_vaf, "min_depth": min_depth, "min_alt_count": min_alt, "error_rate": error, "alpha": alpha},
            "n_tested": tested, "n_detected": len(hits), "variants": hits, "warnings": warnings}


# ---------------------------------------------------------------------------
# repeat expansions in sequences
# ---------------------------------------------------------------------------
def repeat_expansions(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    motifs = options.get("motif", options.get("motifs"))
    if isinstance(motifs, str):
        motifs = [motifs]
    if not isinstance(motifs, Sequence) or not motifs or any(not isinstance(m, str) or not re.fullmatch(r"[ACGTacgt]{1,12}", m) for m in motifs):
        fail("options.motif (a DNA repeat unit of 1-12 bases, or a list) is required.", "VARIANTS_BAD_OPTION", "Example: options={'motif': 'CAG', 'expansion_threshold': 36}.", config=True)
    threshold = option_int(options, "expansion_threshold", None, minimum=1)
    both = option_bool(options, "both_strands", False)
    if isinstance(value, Mapping):
        seqs = dict(value.get("sequences", value)) if isinstance(value.get("sequences", value), Mapping) else None
        if seqs is None or any(not isinstance(v, str) for v in seqs.values()):
            fail("sequence mapping must map names to sequences.", "VARIANTS_BAD_INPUT", "Use {'name': 'ACGT...'}.")
    elif isinstance(value, str):
        if value.lstrip().startswith(">"):
            seqs = {}
            name = None
            for line in value.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0] if line[1:].strip() else f"seq{len(seqs) + 1}"
                    seqs[name] = ""
                elif name is not None:
                    seqs[name] += line.strip()
        else:
            seqs = {"sequence_1": value}
    elif isinstance(value, Sequence) and all(isinstance(v, str) for v in value):
        seqs = {f"sequence_{i + 1}": v for i, v in enumerate(value)}
    else:
        fail("unsupported input.", "VARIANTS_BAD_INPUT", "Pass a DNA string, FASTA text, a list of strings or {name: sequence}.")
    comp = str.maketrans("ACGT", "TGCA")
    results = []
    for name, seq in seqs.items():
        seq = re.sub(r"\s+", "", seq).upper()
        if len(seq) > 50_000_000:
            fail("sequence too long.", "VARIANTS_BAD_INPUT", "Pass the locus window, not a chromosome.")
        for motif in motifs:
            motif = motif.upper()
            strands = [(motif, "+")] + ([(motif.translate(comp)[::-1], "-")] if both else [])
            for unit, strand in strands:
                best = None
                for m in re.finditer(rf"(?:{unit})+", seq):
                    units = len(m.group(0)) // len(unit)
                    if best is None or units > best[0]:
                        best = (units, m.start() + 1, m.end())
                units = best[0] if best else 0
                results.append({"sequence": name, "motif": unit, "strand": strand, "longest_pure_run_units": units,
                                "start": best[1] if best else None, "end": best[2] if best else None,
                                "expansion": None if threshold is None else units >= threshold})
    return {"method": "longest uninterrupted tandem run of the repeat unit in each sequence (exact matches, every phase; interruptions "
            "end a run); 'expansion' is true when the run has at least expansion_threshold units", "expansion_threshold": threshold, "results": results}


# ---------------------------------------------------------------------------
# mobile element insertion calls in an SV VCF
# ---------------------------------------------------------------------------
_ME_CLASSES = {"ALU": "Alu", "L1": "LINE-1", "LINE1": "LINE-1", "SVA": "SVA", "HERVK": "HERV-K", "HERV": "HERV", "LINE": "LINE-1"}


def mobile_element_insertions(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    wanted = options.get("families")
    if isinstance(wanted, str):
        wanted = [wanted]
    hits = []
    for rec in vs.records:
        for alt in rec["alts"]:
            family = None
            if alt.startswith("<INS:ME") or alt.startswith("<INS"):
                parts = alt.strip("<>").split(":")
                if len(parts) >= 3 and parts[1].upper() == "ME":
                    family = parts[2].upper()
                elif len(parts) == 2 and parts[1].upper() in _ME_CLASSES:
                    family = parts[1].upper()
            elif alt.startswith("<") and alt.strip("<>").upper() in _ME_CLASSES:
                family = alt.strip("<>").upper()
            info = rec["info"]
            svtype = str(info.get("SVTYPE", "")).upper()
            if family is None and svtype in _ME_CLASSES:
                family = svtype
            if family is None and svtype in {"INS", "MEI"} and isinstance(info.get("MEINFO"), str):
                family = info["MEINFO"].split(",")[0].upper()
            if family is None and isinstance(info.get("MEINFO"), str):
                family = info["MEINFO"].split(",")[0].upper()
            if family is None:
                continue
            label = _ME_CLASSES.get(family, family)
            if wanted is not None and label not in wanted and family not in wanted:
                continue
            hits.append({"chrom": rec["chrom"], "pos": rec["pos"], "id": rec["id"], "alt": alt, "family": label,
                         "svlen": info.get("SVLEN"), "filter": rec["filter"], "qual": rec["qual"]})
    return {"method": "selection of mobile-element-insertion records from an SV VCF: ALT <INS:ME:FAMILY>, INFO/SVTYPE in {ALU,L1,SVA,...} "
            "or INFO/MEINFO; no read-level evidence is examined", "n_records": len(vs.records), "n_insertions": len(hits), "insertions": hits}


# ---------------------------------------------------------------------------
# density plot
# ---------------------------------------------------------------------------
def density_plot(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    bin_size = option_int(options, "bin_size", 1_000_000, minimum=1)
    assert bin_size is not None
    vs = load_variant_set(value)
    chroms = options.get("chroms")
    if isinstance(chroms, str):
        chroms = [chroms]
    lengths = options.get("chrom_lengths") or {}
    if not isinstance(lengths, Mapping):
        fail("options.chrom_lengths must map contig to length.", "VARIANTS_BAD_OPTION", "Use {'chr1': 248956422}.", config=True)
    counts: dict[str, Counter[int]] = {}
    for rec in vs.records:
        if chroms is not None and rec["chrom"] not in chroms:
            continue
        counts.setdefault(rec["chrom"], Counter())[(rec["pos"] - 1) // bin_size] += 1
    if not counts:
        fail("no variants to plot.", "VARIANTS_BAD_INPUT", "Provide variants (and check options.chroms).")
    bins = []
    for chrom in sorted(counts, key=chrom_key):
        last = max(counts[chrom])
        if chrom in lengths:
            last = max(last, (int(lengths[chrom]) - 1) // bin_size)
        for index in range(last + 1):
            c = counts[chrom].get(index, 0)
            bins.append({"chrom": chrom, "start": index * bin_size + 1, "end": (index + 1) * bin_size, "count": c,
                         "density_per_mb": c * 1_000_000 / bin_size})
        if len(bins) > 20_000:
            fail("more than 20000 bins; increase options.bin_size.", "VARIANTS_BAD_OPTION", "Use a larger bin_size.", config=True)
    peak = max(b["count"] for b in bins) or 1
    width, height, pad = max(320, min(2400, 6 * len(bins) + 80)), 220, 30
    barw = (width - 2 * pad) / len(bins)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             f'<text x="{pad}" y="16" font-size="12">Variant density ({bin_size} bp bins, max {peak} per bin)</text>']
    previous = None
    for i, b in enumerate(bins):
        h = (height - 2 * pad) * b["count"] / peak
        x = pad + i * barw
        parts.append(f'<rect x="{x:.2f}" y="{height - pad - h:.2f}" width="{max(barw - 0.5, 0.5):.2f}" height="{h:.2f}" fill="#4477aa"/>')
        if b["chrom"] != previous:
            parts.append(f'<text x="{x:.2f}" y="{height - 10}" font-size="9">{escape(b["chrom"])}</text>')
            previous = b["chrom"]
    parts.append("</svg>")
    svg = "".join(parts)
    path = write_output({**options, "output": options.get("svg_output")}, svg) if options.get("svg_output") else None
    return {"bin_size": bin_size, "bins": bins, "svg": svg, "output_path": path, "max_count": peak}


# ---------------------------------------------------------------------------
# coding effect on a supplied CDS
# ---------------------------------------------------------------------------
_CODONS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASES = "TCAG"
_AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}
_CODON_TABLE = {a + b + c: _CODONS[16 * i + 4 * j + k]
                for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def translate_cds(seq: str) -> str:
    """Translate with the standard genetic code up to and including the first stop; 'X' for ambiguous codons."""

    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = _CODON_TABLE.get(seq[i:i + 3], "X")
        out.append(aa)
        if aa == "*":
            break
    return "".join(out)


def _p3(aas: str) -> str:
    return "".join(_AA3.get(a, "Xaa") for a in aas)


def _coding_effect(cds: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    ref = "" if ref in {"-", "."} else ref.upper()
    alt = "" if alt in {"-", "."} else alt.upper()
    offset, tref, talt = trim_alleles(ref, alt)
    p = pos + offset
    if not tref and not talt:
        fail("REF and ALT are identical.", "VARIANTS_BAD_RECORD", "Provide differing alleles.")
    if p < 1 or p - 1 + len(tref) > len(cds):
        fail(f"variant at CDS position {pos} lies outside the CDS (length {len(cds)}).", "VARIANTS_OUT_OF_RANGE", "Use 1-based positions inside the coding sequence.")
    if cds[p - 1:p - 1 + len(tref)] != tref:
        fail(f"REF {tref} does not match the CDS at position {p} (CDS has {cds[p - 1:p - 1 + len(tref)]}).", "VARIANTS_REF_MISMATCH",
             "Give REF/ALT on the coding strand, 1-based CDS positions.")
    mut = cds[:p - 1] + talt + cds[p - 1 + len(tref):]
    ref_prot, alt_prot = translate_cds(cds), translate_cds(mut)
    first_codon = (p - 1) // 3 + 1
    result: dict[str, Any] = {"cds_position": p, "ref": tref, "alt": talt}
    delta = len(talt) - len(tref)
    if delta == 0:
        last_codon = (p - 1 + len(tref) - 1) // 3 + 1
        changes = []
        for c in range(first_codon, last_codon + 1):
            ra = _CODON_TABLE.get(cds[3 * c - 3:3 * c], "X")
            aa_alt = _CODON_TABLE.get(mut[3 * c - 3:3 * c], "X")
            changes.append((c, ra, aa_alt, cds[3 * c - 3:3 * c], mut[3 * c - 3:3 * c]))
        kinds = []
        notes = []
        for c, ra, aa_alt, cr, ca in changes:
            if c == 1 and ra == "M" and aa_alt != "M":
                kinds.append("start_lost")
                notes.append("Met1?")
            elif ra == aa_alt:
                kinds.append("synonymous_variant")
                notes.append(f"{_p3(ra)}{c}=")
            elif aa_alt == "*":
                kinds.append("stop_gained")
                notes.append(f"{_p3(ra)}{c}Ter")
            elif ra == "*":
                kinds.append("stop_lost")
                notes.append(f"Ter{c}{_p3(aa_alt)}ext*?")
            else:
                kinds.append("missense_variant")
                notes.append(f"{_p3(ra)}{c}{_p3(aa_alt)}")
        order = ["start_lost", "stop_gained", "stop_lost", "missense_variant", "synonymous_variant"]
        result["consequence"] = min(kinds, key=order.index)
        result["codons"] = [{"codon_index": c, "ref_codon": cr, "alt_codon": ca, "ref_aa": ra, "alt_aa": aa_alt} for c, ra, aa_alt, cr, ca in changes]
        result["protein_position"] = first_codon
        result["hgvs_p"] = "p." + notes[0] if len(notes) == 1 else "p.[" + ";".join(notes) + "]"
        return result
    if delta % 3:
        i = 0
        while i < min(len(ref_prot), len(alt_prot)) and ref_prot[i] == alt_prot[i]:
            i += 1
        ref_aa = ref_prot[i] if i < len(ref_prot) else "X"
        if first_codon == 1 and not alt_prot.startswith("M"):
            result.update(consequence="start_lost", protein_position=1, hgvs_p="p.Met1?")
            return result
        if i < len(alt_prot) and alt_prot[i] == "*":
            result.update(consequence="stop_gained", protein_position=i + 1, hgvs_p=f"p.{_p3(ref_aa)}{i + 1}Ter")
            return result
        alt_aa = alt_prot[i] if i < len(alt_prot) else "X"
        ter = f"Ter{len(alt_prot) - i}" if alt_prot.endswith("*") else "Ter?"
        result.update(consequence="frameshift_variant", protein_position=i + 1, hgvs_p=f"p.{_p3(ref_aa)}{i + 1}{_p3(alt_aa)}fs{ter}")
        return result
    # in-frame
    if first_codon == 1 and not alt_prot.startswith("M"):
        result.update(consequence="start_lost", protein_position=1, hgvs_p="p.Met1?")
        return result
    a = 0
    while a < min(len(ref_prot), len(alt_prot)) and ref_prot[a] == alt_prot[a]:
        a += 1
    rb, ab = len(ref_prot), len(alt_prot)
    while rb > a and ab > a and ref_prot[rb - 1] == alt_prot[ab - 1]:
        rb -= 1
        ab -= 1
    ref_seg, alt_seg = ref_prot[a:rb], alt_prot[a:ab]
    if "*" in alt_seg and alt_seg.index("*") < len(alt_seg):
        star = a + alt_seg.index("*")
        ref_aa = ref_prot[star] if star < len(ref_prot) else "X"
        result.update(consequence="stop_gained", protein_position=star + 1, hgvs_p=f"p.{_p3(ref_aa)}{star + 1}Ter")
        return result
    if not tref and talt:
        kind = "inframe_insertion"
    elif tref and not talt:
        kind = "inframe_deletion"
    else:
        kind = "protein_altering_variant"
    if not ref_seg and alt_seg:
        hgvs = f"p.{_p3(ref_prot[a - 1])}{a}_{_p3(ref_prot[a])}{a + 1}ins{_p3(alt_seg)}" if 0 < a < len(ref_prot) else "p.?"
    elif ref_seg and not alt_seg:
        hgvs = f"p.{_p3(ref_seg[0])}{a + 1}del" if len(ref_seg) == 1 else f"p.{_p3(ref_seg[0])}{a + 1}_{_p3(ref_seg[-1])}{rb}del"
    elif ref_seg:
        left = f"{_p3(ref_seg[0])}{a + 1}" if len(ref_seg) == 1 else f"{_p3(ref_seg[0])}{a + 1}_{_p3(ref_seg[-1])}{rb}"
        hgvs = f"p.{left}delins{_p3(alt_seg)}"
    else:
        hgvs = "p.="
        kind = "synonymous_variant"
    result.update(consequence=kind, protein_position=a + 1, hgvs_p=hgvs)
    return result


def coding_effect(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    cds = options.get("cds")
    items: Any = value
    if isinstance(value, Mapping):
        cds = value.get("cds", cds)
        items = value.get("variants", value if "cds_position" in value else None)
    if not isinstance(cds, str) or not re.fullmatch(r"[ACGTacgt]+", cds) or len(cds) > 1_000_000:
        fail("a coding sequence (ACGT only) is required as input.cds or options.cds.", "VARIANTS_BAD_INPUT",
             "Pass {'cds': 'ATG...', 'variants': [{'cds_position': 5, 'ref': 'C', 'alt': 'T'}]}.")
    cds = cds.upper()
    if isinstance(items, Mapping):
        items = [items]
    if not isinstance(items, Sequence) or isinstance(items, (str, bytes)) or not items or len(items) > 10_000:
        fail("variants must be a non-empty list.", "VARIANTS_BAD_INPUT", "Use [{'cds_position', 'ref', 'alt'}] or 'c.5C>T' strings.")
    rows = []
    for item in items:
        if isinstance(item, str):
            m = re.fullmatch(r"c\.(\d+)([ACGTacgt]+)>([ACGTacgt]+)", item.strip())
            if not m:
                fail(f"cannot parse {item!r}.", "VARIANTS_BAD_INPUT", "Use 'c.5C>T' or a mapping.")
            item = {"cds_position": int(m.group(1)), "ref": m.group(2), "alt": m.group(3)}
        if not isinstance(item, Mapping) or not isinstance(item.get("cds_position"), int) or isinstance(item.get("cds_position"), bool):
            fail("each variant needs an integer cds_position, ref and alt.", "VARIANTS_BAD_INPUT", "Use {'cds_position': 5, 'ref': 'C', 'alt': 'T'}.")
        ref, alt = item.get("ref"), item.get("alt")
        if not isinstance(ref, str) or not isinstance(alt, str):
            fail("ref and alt must be strings.", "VARIANTS_BAD_INPUT", "Use '-' or '' for an empty allele.")
        rows.append(_coding_effect(cds, item["cds_position"], ref, alt))
    warnings = []
    if len(cds) % 3:
        warnings.append("CDS length is not a multiple of 3")
    if not cds.startswith("ATG"):
        warnings.append("CDS does not start with ATG")
    return {"method": "translation of the reference and mutated CDS with the standard genetic code; REF/ALT are on the coding strand, "
            "positions are 1-based CDS coordinates; protein HGVS is not 3'-normalised", "cds_length": len(cds),
            "protein_length": len(translate_cds(cds).rstrip("*")), "effects": rows, "warnings": warnings}


# ---------------------------------------------------------------------------
# construct haplotypes from phased genotypes
# ---------------------------------------------------------------------------
def construct_haplotypes(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    if not vs.samples:
        fail("haplotype construction needs phased sample genotypes.", "VARIANTS_NO_SAMPLES", "Provide GT with | separators.")
    wanted = options.get("samples")
    if wanted is None and options.get("sample") is not None:
        wanted = [options["sample"]]
    if isinstance(wanted, str):
        wanted = [wanted]
    names = list(vs.samples) if wanted is None else list(wanted)
    on_unphased = option_choice(options, "on_unphased", "error", ("error", "skip"))
    ref = reference_from_options(options)
    warnings: list[str] = []
    result: dict[str, Any] = {}
    for name in names:
        if name not in vs.samples:
            fail(f"sample {name!r} not found.", "VARIANTS_UNKNOWN_SAMPLE", "Use names from the VCF header.", config=True)
        idx = vs.samples.index(name)
        per_chrom: dict[str, list[tuple[dict[str, Any], list[str]]]] = {}
        for rec in vs.records:
            alleles, phased = call_gt(rec["calls"][idx])
            if len(alleles) != 2 or not is_called(alleles):
                warnings.append(f"{name} {site_name(rec)}: missing/non-diploid genotype skipped")
                continue
            if not phased and alleles[0] != alleles[1]:
                if on_unphased == "error":
                    fail(f"unphased heterozygous genotype for {name} at {rec['chrom']}:{rec['pos']}.", "VARIANTS_NEEDS_PHASE",
                         "Phase the genotypes first or set options.on_unphased='skip'.", config=True)
                warnings.append(f"{name} {site_name(rec)}: unphased heterozygote skipped")
                continue
            pool = [rec["ref"], *rec["alts"]]
            if any(a >= len(pool) for a in alleles):  # type: ignore[operator]
                fail(f"genotype allele out of range at {rec['chrom']}:{rec['pos']}.", "VARIANTS_BAD_GENOTYPE", "Check GT indices against ALT.")
            per_chrom.setdefault(rec["chrom"], []).append((rec, [pool[a] for a in alleles]))  # type: ignore[index]
        chrom_out = {}
        for chrom, items in per_chrom.items():
            items.sort(key=lambda t: t[0]["pos"])
            haps = {}
            for h in (0, 1):
                allele_rows = [{"pos": rec["pos"], "ref": rec["ref"], "allele": al[h]} for rec, al in items]
                entry: dict[str, Any] = {"alleles": allele_rows}
                if ref is not None:
                    entry["sequence"] = _apply_to_reference(ref, chrom, [(rec["pos"], rec["ref"], al[h]) for rec, al in items], options, warnings)
                haps[str(h + 1)] = entry
            chrom_out[chrom] = haps
        result[name] = chrom_out
    return {"method": "per-sample haplotype alleles taken from the phased GT (first allele = haplotype 1); optional reference sequence with the "
            "haplotype's alleles substituted over the window", "haplotypes": result, "warnings": warnings[:100]}


def _apply_to_reference(ref: Reference, chrom: str, edits: list[tuple[int, str, str]], options: Mapping[str, object], warnings: list[str]) -> str | None:
    lo = min(p for p, _, _ in edits)
    hi = max(p + len(r) - 1 for p, r, _ in edits)
    start = option_int(options, "window_start", lo, minimum=1)
    end = option_int(options, "window_end", hi, minimum=1)
    assert start is not None and end is not None
    if start > lo or end < hi:
        fail("window does not cover all variants.", "VARIANTS_BAD_OPTION", "Use a window covering every variant.", config=True)
    seq = ref.fetch(chrom, start, end - start + 1)
    if seq is None:
        warnings.append(f"{chrom}: window {start}-{end} is outside the reference; sequence omitted")
        return None
    pieces = []
    cursor = start
    for pos, r, a in sorted(edits):
        if is_symbolic(a):
            warnings.append(f"{chrom}:{pos}: symbolic allele {a} ignored")
            continue
        if pos < cursor:
            fail(f"overlapping variants at {chrom}:{pos}.", "VARIANTS_OVERLAP", "Remove or resolve overlapping records.")
        actual = seq[pos - start:pos - start + len(r)]
        if actual.upper() != r.upper():
            fail(f"REF {r} at {chrom}:{pos} does not match the reference ({actual}).", "VARIANTS_REF_MISMATCH", "Use the matching reference.")
        pieces.append(seq[cursor - start:pos - start])
        pieces.append(a.upper())
        cursor = pos + len(r)
    pieces.append(seq[cursor - start:])
    return "".join(pieces)


# ---------------------------------------------------------------------------
# naive imputation and trio phasing
# ---------------------------------------------------------------------------
def impute_missing(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    if not vs.samples:
        fail("imputation needs sample genotypes.", "VARIANTS_NO_SAMPLES", "Provide genotypes with some missing calls.")
    imputed = unresolved = 0
    for rec in vs.records:
        counts: Counter[str] = Counter()
        for call in rec["calls"]:
            alleles, _ = call_gt(call)
            if is_called(alleles):
                counts["/".join(str(a) for a in sorted(alleles))] += 1  # type: ignore[type-var]
        if not counts:
            unresolved += sum(1 for c in rec["calls"] if not is_called(call_gt(c)[0]))
            continue
        top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        for call in rec["calls"]:
            if not is_called(call_gt(call)[0]):
                call["GT"] = top
                call["IMPUTED"] = "1"
                imputed += 1
        if "IMPUTED" not in rec["format"] and any("IMPUTED" in c for c in rec["calls"]):
            rec["format"].append("IMPUTED")
    text = to_vcf_text(vs)
    path = write_output(options, text)
    return {"method": "naive per-site modal-genotype imputation: each missing call receives the most frequent called genotype at the same "
            "site (ties broken by genotype string); it does not use linkage disequilibrium or a reference panel and is NOT a "
            "substitute for Beagle/IMPUTE-style imputation", "imputed_genotypes": imputed, "sites_without_any_call": unresolved,
            "records": [record_view(r, vs.samples) for r in vs.records], "vcf": text, "output_path": path}


def phase_by_transmission(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    trio = options.get("trio") or {k: options.get(k) for k in ("child", "father", "mother")}
    if not isinstance(trio, Mapping) or not all(isinstance(trio.get(k), str) for k in ("child", "father", "mother")):
        fail("options.trio={'child','father','mother'} sample names are required.", "VARIANTS_BAD_OPTION",
             "Phasing is by transmission from parents; statistical phasing needs Beagle/SHAPEIT.", config=True)
    vs = load_variant_set(value)
    for key in ("child", "father", "mother"):
        if trio[key] not in vs.samples:
            fail(f"{key} sample {trio[key]!r} not found.", "VARIANTS_UNKNOWN_SAMPLE", "Use names from the VCF header.", config=True)
    ci, fi, mi = (vs.samples.index(trio[k]) for k in ("child", "father", "mother"))
    stats = Counter()
    for rec in vs.records:
        child, _ = call_gt(rec["calls"][ci])
        if len(child) != 2 or not is_called(child):
            stats["child_missing"] += 1
            continue
        if child[0] == child[1]:
            rec["calls"][ci]["GT"] = f"{child[0]}|{child[1]}"
            stats["homozygous"] += 1
            continue
        fa, _ = call_gt(rec["calls"][fi])
        ma, _ = call_gt(rec["calls"][mi])
        fset = set(fa) if is_called(fa) else None
        mset = set(ma) if is_called(ma) else None
        x, y = child
        ok_xy = (fset is None or x in fset) and (mset is None or y in mset)
        ok_yx = (fset is None or y in fset) and (mset is None or x in mset)
        if ok_xy and not ok_yx:
            rec["calls"][ci]["GT"] = f"{x}|{y}"
            stats["phased"] += 1
        elif ok_yx and not ok_xy:
            rec["calls"][ci]["GT"] = f"{y}|{x}"
            stats["phased"] += 1
        elif ok_xy and ok_yx:
            stats["ambiguous"] += 1
        else:
            stats["mendelian_error"] += 1
    text = to_vcf_text(vs)
    path = write_output(options, text)
    return {"method": "trio phasing by transmission: the child's heterozygous genotype is written paternal|maternal when exactly one "
            "assignment of its two alleles is compatible with the parents' alleles; ambiguous, inconsistent and uninformative sites stay unphased",
            "trio": dict(trio), "counts": dict(stats), "vcf": text, "output_path": path}
