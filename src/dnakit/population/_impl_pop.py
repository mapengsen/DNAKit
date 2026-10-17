"""Native population-genetics statistics on dosage genotype matrices (allele frequencies, HWE, Fst,
LD, GRM, simulation).

Genotypes: ``{'genotypes': samples x variants alternate-allele dosages, 'sample_ids',
'variant_ids', 'populations'}`` (see ``dnakit.gwas._impl_common.parse_genotype_matrix``).
Heterozygosity-based statistics need hard calls (0/1/2).
"""

from __future__ import annotations

import math
import random
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.gwas import _impl_stats as st
from dnakit.gwas._impl_common import (
    CandidateUnavailable,
    nice_ticks,
    opt_bool,
    opt_choice,
    opt_number,
    parse_genotype_matrix,
    serialize_records,
    svg_escape,
    svg_wrap,
    write_text_file,
)

MAX_VARIANTS_PAIRWISE = 300


def _genotypes(
    source: object, *, min_samples: int = 2, need_populations: bool = False
) -> dict[str, Any]:
    if source is None or isinstance(source, (str, bytes)):
        raise InputFormatError(
            "input must be a genotype matrix or a mapping with a 'genotypes' key.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'genotypes': [[0,1,2],...], 'populations': [...]}.",
        )
    g = parse_genotype_matrix(source, min_samples=min_samples)
    if need_populations and "populations" not in g:
        raise InputFormatError(
            "a 'populations' list (one label per sample) is required.",
            code="MISSING_POPULATIONS",
            hint="Add 'populations': ['A', 'A', 'B', ...] to the input mapping.",
        )
    return g


def _hard_calls(g: Mapping[str, Any]) -> list[list[int | None]]:
    out = []
    for row in g["matrix"]:
        new: list[int | None] = []
        for v in row:
            if v is None:
                new.append(None)
            elif abs(v - round(v)) > 1e-9:
                raise InputFormatError(
                    "this statistic needs hard genotype calls (dosages 0, 1 or 2), but fractional "
                    "dosages were given.",
                    code="FRACTIONAL_DOSAGE",
                    hint="Round imputed dosages to best-guess genotypes first.",
                )
            else:
                new.append(round(v))
        out.append(new)
    return out


def _groups(g: Mapping[str, Any]) -> OrderedDict[str, list[int]]:
    groups: OrderedDict[str, list[int]] = OrderedDict()
    pops = g.get("populations")
    if pops is None:
        groups["all"] = list(range(g["n_samples"]))
    else:
        for i, p in enumerate(pops):
            groups.setdefault(str(p), []).append(i)
    return groups


# ------------------------------------------------------------- allele frequency, heterozygosity, LD
def _allele_frequency(g: Mapping[str, Any], rows: Sequence[int]) -> list[dict[str, Any]]:
    out = []
    for j in range(g["n_variants"]):
        vals = [g["matrix"][i][j] for i in rows if g["matrix"][i][j] is not None]
        n = len(vals)
        p = sum(vals) / (2.0 * n) if n else None
        out.append(
            {
                "variant_id": g["variants"][j],
                "n_called": n,
                "alt_allele_frequency": p,
                "maf": None if p is None else min(p, 1 - p),
                "missing_rate": 1 - n / len(rows),
            }
        )
    return out


def _heterozygosity(g: Mapping[str, Any]) -> list[dict[str, Any]]:
    hard = _hard_calls(g)
    out = []
    for j in range(g["n_variants"]):
        vals = [row[j] for row in hard if row[j] is not None]
        n = len(vals)
        if n == 0:
            out.append(
                {"variant_id": g["variants"][j], "n_called": 0, "observed_heterozygosity": None}
            )
            continue
        p = sum(vals) / (2.0 * n)
        ho = sum(1 for v in vals if v == 1) / n
        he = 2 * p * (1 - p)
        he_unbiased = he * 2 * n / (2 * n - 1) if n > 0 and 2 * n > 1 else None
        out.append(
            {
                "variant_id": g["variants"][j],
                "n_called": n,
                "alt_allele_frequency": p,
                "observed_heterozygosity": ho,
                "expected_heterozygosity": he,
                "expected_heterozygosity_unbiased": he_unbiased,
                "f_statistic": (1 - ho / he) if he > 0 else None,
            }
        )
    return out


def _ld_pairs(g: Mapping[str, Any], limit: int) -> list[dict[str, Any]]:
    m = g["n_variants"]
    if m > limit:
        raise ConfigurationError(
            f"LD is computed for at most {limit} variants at a time (got {m}).",
            code="TOO_MANY_VARIANTS",
            hint="Subset the variants or raise options['max_variants'].",
        )
    cols = [[row[j] for row in g["matrix"]] for j in range(m)]
    out = []
    for a in range(m):
        for b in range(a + 1, m):
            pairs = [
                (x, y)
                for x, y in zip(cols[a], cols[b], strict=False)
                if x is not None and y is not None
            ]
            if len(pairs) < 3:
                continue
            xs, ys = [p[0] for p in pairs], [p[1] for p in pairs]
            n = len(pairs)
            mx, my = sum(xs) / n, sum(ys) / n
            cov = math.fsum((x - mx) * (y - my) for x, y in pairs) / n
            vx = math.fsum((x - mx) ** 2 for x in xs) / n
            vy = math.fsum((y - my) ** 2 for y in ys) / n
            if vx <= 0 or vy <= 0:
                continue
            pa, pb = mx / 2.0, my / 2.0
            delta = cov / 2.0
            dmax = (
                min(pa * (1 - pb), (1 - pa) * pb)
                if delta > 0
                else min(pa * pb, (1 - pa) * (1 - pb))
            )
            out.append(
                {
                    "variant_1": g["variants"][a],
                    "variant_2": g["variants"][b],
                    "n": n,
                    "r2": cov * cov / (vx * vy),
                    "d_composite": delta,
                    "d_prime": (abs(delta) / dmax) if dmax > 0 else None,
                }
            )
    return out


def calculate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    statistic = opt_choice(
        opts, "statistic", "allele_frequency", ("allele_frequency", "heterozygosity", "ld")
    )
    g = _genotypes(source)
    if statistic == "allele_frequency":
        result: dict[str, Any] = {
            "statistic": statistic,
            "overall": _allele_frequency(g, range(g["n_samples"])),
        }
        if "populations" in g:
            result["by_population"] = {
                name: _allele_frequency(g, rows) for name, rows in _groups(g).items()
            }
        return result, []
    if statistic == "heterozygosity":
        return {
            "statistic": statistic,
            "variants": _heterozygosity(g),
            "definition": (
                "Ho = fraction of called samples with dosage 1; He = 2pq from sample frequency; F "
                "= 1 - Ho/He"
            ),
        }, []
    limit = opt_number(
        opts, "max_variants", 100, minimum=2, maximum=MAX_VARIANTS_PAIRWISE, integer=True
    )
    return {
        "statistic": statistic,
        "pairs": _ld_pairs(g, limit),
        "definition": (
            "composite (genotype-dosage) LD: r2 = squared Pearson correlation of dosages; D = "
            "cov(dosage_a, dosage_b)/2; "
            "D' = |D| / Dmax with Dmax = min(pA qB, qA pB) if D > 0 else min(pA pB, qA qB)"
        ),
    }, []


# --------------------------------------------------------------------------- Fst
def _wc_components(counts: Sequence[tuple[int, float, float]]) -> tuple[float, float, float] | None:
    """Weir & Cockerham (1984) variance components (a, b, c) for one biallelic locus.

    ``counts`` holds (n_i diploid individuals, p_i alternate-allele frequency, h_i heterozygote
    frequency) per population.
    """

    pops = [c for c in counts if c[0] > 0]
    r = len(pops)
    if r < 2:
        return None
    total = sum(n for n, _, _ in pops)
    nbar = total / r
    if nbar <= 1:
        return None
    nc = (total - sum(n * n for n, _, _ in pops) / total) / (r - 1)
    pbar = sum(n * p for n, p, _ in pops) / total
    s2 = sum(n * (p - pbar) ** 2 for n, p, _ in pops) / ((r - 1) * nbar)
    hbar = sum(n * h for n, _, h in pops) / total
    a = (nbar / nc) * (
        s2 - (1.0 / (nbar - 1.0)) * (pbar * (1 - pbar) - (r - 1) / r * s2 - hbar / 4.0)
    )
    b = (nbar / (nbar - 1.0)) * (
        pbar * (1 - pbar) - (r - 1) / r * s2 - (2 * nbar - 1) / (4 * nbar) * hbar
    )
    c = hbar / 2.0
    return a, b, c


def _pop_stats(
    hard: Sequence[Sequence[int | None]], rows: Sequence[int], j: int
) -> tuple[int, float, float]:
    vals = [hard[i][j] for i in rows if hard[i][j] is not None]
    n = len(vals)
    if n == 0:
        return 0, 0.0, 0.0
    return n, sum(vals) / (2.0 * n), sum(1 for v in vals if v == 1) / n


def _fst_wc(
    hard: Sequence[Sequence[int | None]], groups: Sequence[Sequence[int]], m: int
) -> tuple[float | None, list[float | None]]:
    per_locus: list[float | None] = []
    sa = sabc = 0.0
    for j in range(m):
        comp = _wc_components([_pop_stats(hard, rows, j) for rows in groups])
        if comp is None:
            per_locus.append(None)
            continue
        a, b, c = comp
        if a + b + c == 0:
            per_locus.append(None)
            continue
        per_locus.append(a / (a + b + c))
        sa += a
        sabc += a + b + c
    return (sa / sabc if sabc else None), per_locus


def _fst_hudson(
    hard: Sequence[Sequence[int | None]], rows1: Sequence[int], rows2: Sequence[int], m: int
) -> tuple[float | None, list[float | None]]:
    per_locus: list[float | None] = []
    sn = sd = 0.0
    for j in range(m):
        n1, p1, _ = _pop_stats(hard, rows1, j)
        n2, p2, _ = _pop_stats(hard, rows2, j)
        h1, h2 = 2 * n1, 2 * n2
        if h1 < 2 or h2 < 2:
            per_locus.append(None)
            continue
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (h1 - 1) - p2 * (1 - p2) / (h2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        if den == 0:
            per_locus.append(None)
            continue
        per_locus.append(num / den)
        sn += num
        sd += den
    return (sn / sd if sd else None), per_locus


def compare(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    method = opt_choice(opts, "method", "weir_cockerham", ("weir_cockerham", "hudson"))
    g = _genotypes(source, need_populations=True)
    hard = _hard_calls(g)
    groups = _groups(g)
    names = list(groups)
    if len(names) < 2:
        raise InputFormatError(
            "at least two populations are needed to compare them.",
            code="TOO_FEW_POPULATIONS",
            hint="Provide at least two distinct population labels.",
        )
    m = g["n_variants"]
    result: dict[str, Any] = {
        "method": method,
        "populations": {n: len(r) for n, r in groups.items()},
    }
    pairwise: dict[str, float | None] = {}
    if method == "hudson":
        for i in range(len(names)):
            for k in range(i + 1, len(names)):
                value, per = _fst_hudson(hard, groups[names[i]], groups[names[k]], m)
                pairwise[f"{names[i]}|{names[k]}"] = value
                if len(names) == 2:
                    result["fst"], result["per_variant"] = (
                        value,
                        dict(zip(g["variants"], per, strict=False)),
                    )
        result["definition"] = (
            "Hudson Fst (Bhatia et al. 2013): ratio of sums over variants of "
            "[(p1-p2)^2 - p1q1/(n1-1) - p2q2/(n2-1)] / [p1 q2 + p2 q1], n = haploid sample sizes"
        )
    else:
        value, per = _fst_wc(hard, list(groups.values()), m)
        result["fst"], result["per_variant"] = value, dict(zip(g["variants"], per, strict=False))
        for i in range(len(names)):
            for k in range(i + 1, len(names)):
                pairwise[f"{names[i]}|{names[k]}"] = _fst_wc(
                    hard, [groups[names[i]], groups[names[k]]], m
                )[0]
        result["definition"] = (
            "Weir & Cockerham (1984) theta = sum(a) / sum(a + b + c) over variants (ratio of "
            "sums), "
            "diploid genotype-based components"
        )
    result["pairwise"] = pairwise
    return result, []


# ------------------------------------------------------------------ estimate (GRM, inbreeding, PCA)
def estimate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    quantity = opt_choice(opts, "quantity", "grm", ("grm", "kinship", "inbreeding", "pca"))
    g = _genotypes(source, min_samples=2)
    if quantity in ("grm", "kinship", "pca"):
        try:
            from dnakit.exceptions import BackendUnavailableError
            from dnakit.genomics.breeding import genomic_relationship_matrix
            from dnakit.genomics.population import genotype_pca
        except ImportError as exc:  # pragma: no cover
            raise CandidateUnavailable(
                "pip install numpy", f"the genomics helpers could not be imported: {exc}"
            ) from exc
        try:
            if quantity == "pca":
                components = opt_number(opts, "components", 2, minimum=1, integer=True)
                res = genotype_pca(
                    g["matrix"], components=components, sample_ids=g["samples"], verbose=False
                )
                return {"quantity": "pca", **res}, []
            res = genomic_relationship_matrix(g["matrix"], verbose=False)
        except BackendUnavailableError as exc:
            raise CandidateUnavailable(
                "pip install numpy", f"numpy is required for this estimate: {exc}"
            ) from exc
        if quantity == "kinship":
            res = dict(res)
            res["relationship_matrix"] = [
                [v / 2.0 for v in row] for row in res["relationship_matrix"]
            ]
            res["method"] = "vanraden_method_1 / 2 (coancestry)"
        return {"quantity": quantity, "sample_ids": g["samples"], **res}, []
    hard = _hard_calls(g)
    n = g["n_samples"]
    # per-marker sample frequencies and expected homozygosity, as in PLINK --het
    exp_hom: list[float | None] = []
    for j in range(g["n_variants"]):
        vals = [row[j] for row in hard if row[j] is not None]
        k = len(vals)
        if k < 2:
            exp_hom.append(None)
            continue
        p = sum(vals) / (2.0 * k)
        exp_hom.append(1.0 - 2 * p * (1 - p) * (2 * k) / (2 * k - 1))
    out = []
    for i in range(n):
        obs = e = nn = 0.0
        for j, row_val in enumerate(hard[i]):
            if row_val is None or exp_hom[j] is None:
                continue
            nn += 1
            e += exp_hom[j]
            obs += 1 if row_val in (0, 2) else 0
        f = (obs - e) / (nn - e) if nn - e > 0 else None
        out.append(
            {
                "sample_id": g["samples"][i],
                "n_markers": int(nn),
                "observed_homozygous": int(obs),
                "expected_homozygous": e,
                "f": f,
            }
        )
    return {
        "quantity": "inbreeding",
        "samples": out,
        "definition": (
            "F = (O_hom - E_hom) / (N - E_hom) with E_hom = sum(1 - 2pq*2n/(2n-1)) over called "
            "markers (PLINK --het method of moments)"
        ),
    }, []


# --------------------------------------------------------------------------- tests
def hwe_exact_p(n_aa: int, n_ab: int, n_bb: int) -> float:
    """Exact two-sided Hardy-Weinberg p-value (Wigginton, Cutler & Abecasis 2005): sum of all
    configurations as or less likely.
    """

    n = n_aa + n_ab + n_bb
    if n == 0:
        return 1.0
    n_a = 2 * n_aa + n_ab
    n_b = 2 * n_bb + n_ab
    rare = min(n_a, n_b)
    probs = [0.0] * (rare + 1)
    mid = rare * (n_a + n_b - rare) // (n_a + n_b) if (n_a + n_b) else 0
    if (mid % 2) != (rare % 2):
        mid += 1
    probs[mid] = 1.0
    total = 1.0
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets > 1:
        probs[curr_hets - 2] = (
            probs[curr_hets]
            * curr_hets
            * (curr_hets - 1.0)
            / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        )
        total += probs[curr_hets - 2]
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = (
            probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        total += probs[curr_hets + 2]
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
    observed = n_ab
    p_obs = probs[observed]
    p_sum = sum(p for p in probs if p <= p_obs * (1 + 1e-12))
    return min(1.0, p_sum / total)


def hwe_chisq_p(n_aa: int, n_ab: int, n_bb: int) -> tuple[float | None, float | None]:
    n = n_aa + n_ab + n_bb
    if n == 0:
        return None, None
    p = (2 * n_bb + n_ab) / (
        2.0 * n
    )  # alternate-allele frequency; counts are (hom ref, het, hom alt)
    exp = [n * (1 - p) ** 2, 2 * n * p * (1 - p), n * p * p]
    if min(exp) == 0:
        return None, None
    obs = [n_aa, n_ab, n_bb]
    chi = sum((o - e) ** 2 / e for o, e in zip(obs, exp, strict=False))
    return chi, st.chi2_sf(chi, 1)


def test(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    which = opt_choice(opts, "test", "hwe", ("hwe", "differentiation"))
    g = _genotypes(source, need_populations=(which == "differentiation"))
    hard = _hard_calls(g)
    records: list[dict[str, Any]] = []
    if which == "hwe":
        method = opt_choice(opts, "method", "exact", ("exact", "chisq"))
        for j in range(g["n_variants"]):
            vals = [row[j] for row in hard if row[j] is not None]
            counts = (vals.count(0), vals.count(1), vals.count(2))
            if method == "exact":
                p, chi = hwe_exact_p(*counts), None
            else:
                chi, p = hwe_chisq_p(*counts)
            records.append(
                {
                    "variant_id": g["variants"][j],
                    "n_hom_ref": counts[0],
                    "n_het": counts[1],
                    "n_hom_alt": counts[2],
                    "chi2": chi,
                    "p": p,
                }
            )
        definition = (
            "exact test of Wigginton et al. (2005): p = sum of probabilities of heterozygote "
            "counts no more likely than observed"
            if method == "exact"
            else "Pearson chi-square (1 df) of genotype counts against HWE expectations from the "
            "sample frequency"
        )
    else:
        groups = _groups(g)
        if len(groups) < 2:
            raise InputFormatError(
                "at least two populations are needed.",
                code="TOO_FEW_POPULATIONS",
                hint="Provide >= 2 population labels.",
            )
        for j in range(g["n_variants"]):
            table = []
            for rows in groups.values():
                vals = [hard[i][j] for i in rows if hard[i][j] is not None]
                alt = sum(vals)
                table.append((2 * len(vals) - alt, alt))
            col_ref, col_alt = sum(t[0] for t in table), sum(t[1] for t in table)
            total = col_ref + col_alt
            if col_ref == 0 or col_alt == 0 or any(sum(t) == 0 for t in table):
                records.append(
                    {"variant_id": g["variants"][j], "chi2": None, "df": len(table) - 1, "p": None}
                )
                continue
            chi = 0.0
            for ref, alt in table:
                size = ref + alt
                for obs, col in ((ref, col_ref), (alt, col_alt)):
                    expected = size * col / total
                    chi += (obs - expected) ** 2 / expected
            records.append(
                {
                    "variant_id": g["variants"][j],
                    "chi2": chi,
                    "df": len(table) - 1,
                    "p": st.chi2_sf(chi, len(table) - 1),
                }
            )
        definition = (
            "Pearson chi-square test of homogeneity of allele counts across populations "
            "(2 x k table, df = k - 1)"
        )
    tested = [r for r in records if r["p"] is not None]
    bonf = st.bonferroni([r["p"] for r in tested])
    bh = st.benjamini_hochberg([r["p"] for r in tested])
    for r, a, b in zip(tested, bonf, bh, strict=False):
        r["p_bonferroni"], r["q_bh"] = a, b
    return {
        "test": which,
        "records": records,
        "n_tested": len(tested),
        "definition": definition,
    }, []


# --------------------------------------------------------------------------- export
def export(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    fmt = opt_choice(opts, "format", "tsv", ("tsv", "csv", "ped_map", "frq"))
    g = _genotypes(source, min_samples=1)
    path = opts.get("path", opts.get("output"))
    overwrite = opt_bool(opts, "overwrite", False)
    result: dict[str, Any] = {
        "format": fmt,
        "n_samples": g["n_samples"],
        "n_variants": g["n_variants"],
    }
    if fmt in ("tsv", "csv"):
        delim = "\t" if fmt == "tsv" else ","

        def cell(v: float | None) -> float | int | None:
            return int(v) if v is not None and v == int(v) else v

        recs = [
            dict({"sample_id": s}, **{v: cell(row[j]) for j, v in enumerate(g["variants"])})
            for s, row in zip(g["samples"], g["matrix"], strict=False)
        ]
        text = serialize_records(recs, ["sample_id", *g["variants"]], delim)
        result["text"] = text
        if path is not None:
            result["path"] = write_text_file(path, text, overwrite=overwrite)
    elif fmt == "frq":
        freq = _allele_frequency(g, range(g["n_samples"]))
        recs = [
            {
                "SNP": r["variant_id"],
                "A1": (g.get("alt") or ["ALT"] * g["n_variants"])[k],
                "A2": (g.get("ref") or ["REF"] * g["n_variants"])[k],
                "MAF": r["maf"],
                "NCHROBS": 2 * r["n_called"],
            }
            for k, r in enumerate(freq)
        ]
        text = serialize_records(recs, ["SNP", "A1", "A2", "MAF", "NCHROBS"], "\t")
        result["text"] = text
        if path is not None:
            result["path"] = write_text_file(path, text, overwrite=overwrite)
    else:
        hard = _hard_calls(g)
        ref = g.get("ref") or ["A"] * g["n_variants"]
        alt = g.get("alt") or ["G"] * g["n_variants"]
        chrom = g.get("chrom") or ["0"] * g["n_variants"]
        pos = g.get("pos") or [0] * g["n_variants"]
        pops = g.get("populations") or ["1"] * g["n_samples"]
        ped_lines = []
        for i, sample in enumerate(g["samples"]):
            cells = []
            for j in range(g["n_variants"]):
                d = hard[i][j]
                cells.append(
                    "0 0"
                    if d is None
                    else {
                        0: f"{ref[j]} {ref[j]}",
                        1: f"{ref[j]} {alt[j]}",
                        2: f"{alt[j]} {alt[j]}",
                    }[d]
                )
            ped_lines.append(" ".join([str(pops[i]), sample, "0", "0", "0", "-9", *cells]))
        map_lines = [f"{chrom[j]}\t{g['variants'][j]}\t0\t{pos[j]}" for j in range(g["n_variants"])]
        result["ped"], result["map"] = "\n".join(ped_lines) + "\n", "\n".join(map_lines) + "\n"
        if path is not None:
            result["path"] = [
                write_text_file(f"{path}.ped", result["ped"], overwrite=overwrite),
                write_text_file(f"{path}.map", result["map"], overwrite=overwrite),
            ]
    return result, []


# --------------------------------------------------------------------------- plot
def _axes(
    width: int,
    height: int,
    left: int,
    top: int,
    right: int,
    bottom: int,
    xr: tuple[float, float],
    yr: tuple[float, float],
    xlabel: str,
    ylabel: str,
) -> list[str]:
    pw, ph = width - left - right, height - top - bottom
    body = [f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="black"/>']
    for t in nice_ticks(xr[0], xr[1], 5):
        x = left + pw * (t - xr[0]) / ((xr[1] - xr[0]) or 1)
        body.append(
            f'<text x="{x:.1f}" y="{height - bottom + 14}" font-size="10" '
            f'text-anchor="middle">{t:g}</text>'
        )
    for t in nice_ticks(yr[0], yr[1], 5):
        y = top + ph * (1 - (t - yr[0]) / ((yr[1] - yr[0]) or 1))
        body.append(
            f'<text x="{left - 6}" y="{y + 3:.1f}" font-size="10" text-anchor="end">{t:g}</text>'
        )
    body.append(
        f'<text x="{left + pw / 2}" y="{height - 8}" font-size="12" '
        f'text-anchor="middle">{svg_escape(xlabel)}</text>'
    )
    body.append(
        f'<text x="14" y="{top + ph / 2}" font-size="12" text-anchor="middle" '
        f'transform="rotate(-90 14 {top + ph / 2})">{svg_escape(ylabel)}</text>'
    )
    return body


def plot(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    kind = opt_choice(opts, "kind", "maf_histogram", ("maf_histogram", "pca", "fst"))
    width = int(opt_number(opts, "width", 640, minimum=300, maximum=3000, integer=True))
    height = int(opt_number(opts, "height", 420, minimum=200, maximum=3000, integer=True))
    left, top, right, bottom = 60, 30, 20, 50
    pw, ph = width - left - right, height - top - bottom
    g = _genotypes(source)
    body: list[str] = []
    if kind == "maf_histogram":
        bins = int(opt_number(opts, "bins", 10, minimum=2, maximum=100, integer=True))
        mafs = [
            r["maf"] for r in _allele_frequency(g, range(g["n_samples"])) if r["maf"] is not None
        ]
        counts = [0] * bins
        for v in mafs:
            counts[min(bins - 1, int(v / 0.5 * bins))] += 1
        top_count = max([*counts, 1])
        body += _axes(
            width,
            height,
            left,
            top,
            right,
            bottom,
            (0.0, 0.5),
            (0.0, float(top_count)),
            "Minor allele frequency",
            "Variants",
        )
        for b, c in enumerate(counts):
            h = ph * c / top_count
            body.append(
                f'<rect x="{left + pw * b / bins:.1f}" y="{top + ph - h:.1f}" '
                f'width="{pw / bins - 1:.1f}" height="{h:.1f}" '
                f'fill="#1f4e79"><title>{c}</title></rect>'
            )
        data: dict[str, Any] = {
            "bin_counts": counts,
            "bin_width": 0.5 / bins,
            "n_variants": len(mafs),
        }
    elif kind == "fst":
        fst_result, _ = compare(source, {"method": opts.get("method", "weir_cockerham")})
        per = [v for v in fst_result["per_variant"].values()] if "per_variant" in fst_result else []
        vals = [v for v in per if v is not None]
        if not vals:
            raise InputFormatError(
                "no per-variant Fst values could be computed.",
                code="NO_FST",
                hint="Provide >= 2 populations with >= 2 samples.",
            )
        lo, hi = min(0.0, min(vals)), max(vals) * 1.05 or 1.0
        body += _axes(
            width,
            height,
            left,
            top,
            right,
            bottom,
            (0.0, float(len(per))),
            (lo, hi),
            "Variant index",
            "Fst",
        )
        for idx, v in enumerate(per):
            if v is not None:
                body.append(
                    f'<circle cx="{left + pw * (idx + 0.5) / len(per):.1f}" '
                    f'cy="{top + ph * (1 - (v - lo) / (hi - lo)):.1f}" r="2.5" fill="#1f4e79"/>'
                )
        data = {"n_points": len(vals)}
    else:
        est, _ = estimate(source, {"quantity": "pca", "components": 2})
        scores = est["scores"]
        pops = g.get("populations") or ["all"] * g["n_samples"]
        names = list(dict.fromkeys(map(str, pops)))
        palette = ["#1f4e79", "#c0392b", "#27ae60", "#8e44ad", "#e69f00", "#16a085", "#7f8c8d"]
        xs, ys = [s[0] for s in scores], [s[1] for s in scores]
        xr = (min(xs), max(xs) if max(xs) > min(xs) else min(xs) + 1)
        yr = (min(ys), max(ys) if max(ys) > min(ys) else min(ys) + 1)
        body += _axes(width, height, left, top, right, bottom, xr, yr, "PC1", "PC2")
        for (x, y), pop in zip([s[:2] for s in scores], pops, strict=False):
            color = palette[names.index(str(pop)) % len(palette)]
            body.append(
                f'<circle cx="{left + pw * (x - xr[0]) / (xr[1] - xr[0]):.1f}" '
                f'cy="{top + ph * (1 - (y - yr[0]) / (yr[1] - yr[0])):.1f}" r="3.5" '
                f'fill="{color}"><title>{svg_escape(pop)}</title></circle>'
            )
        data = {"n_points": len(scores), "populations": names}
    svg = svg_wrap(width, height, body, f"population {kind} plot")
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, svg, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"kind": kind, "svg": svg, "path": written, **data}, []


# --------------------------------------------------------------------------- simulate
def simulate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    params: dict[str, Any] = dict(opts)
    if isinstance(source, Mapping):
        params.update(source)
    elif source is not None and not isinstance(source, (str, int, float)):
        raise InputFormatError(
            "input must be None or a mapping of simulation parameters.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'model': 'balding_nichols', 'fst': 0.05, 'seed': 1}.",
        )
    model = opt_choice(params, "model", "balding_nichols", ("balding_nichols", "wright_fisher"))
    seed = opt_number(params, "seed", 0, integer=True)
    rng = random.Random(seed)
    if model == "balding_nichols":
        k = opt_number(params, "n_populations", 2, minimum=1, maximum=20, integer=True)
        n = opt_number(
            params, "n_samples_per_population", 50, minimum=1, maximum=5000, integer=True
        )
        m = opt_number(params, "n_variants", 100, minimum=1, maximum=50_000, integer=True)
        fst = opt_number(params, "fst", 0.1, minimum=0, maximum=1, exclusive=True)
        maf_lo = opt_number(params, "maf_min", 0.1, minimum=0.001, maximum=0.5)
        if k * n * m > 5_000_000:
            raise ConfigurationError(
                "n_populations * n_samples_per_population * n_variants exceeds 5,000,000.",
                code="INPUT_TOO_LARGE",
                hint="Simulate fewer samples or variants.",
            )
        ancestral = [rng.uniform(maf_lo, 1 - maf_lo) for _ in range(m)]
        scale = (1 - fst) / fst
        matrix, pops = [], []
        pop_freqs = []
        for pidx in range(k):
            freqs = [rng.betavariate(p * scale, (1 - p) * scale) for p in ancestral]
            pop_freqs.append(freqs)
            for _ in range(n):
                matrix.append([(rng.random() < f) + (rng.random() < f) for f in freqs])
                pops.append(f"pop{pidx + 1}")
        ids = [f"var{j + 1}" for j in range(m)]
        return {
            "model": "balding_nichols",
            "genotypes": {
                "genotypes": matrix,
                "populations": pops,
                "variant_ids": ids,
                "sample_ids": [f"s{i + 1}" for i in range(len(matrix))],
            },
            "ancestral_frequencies": ancestral,
            "population_frequencies": pop_freqs,
            "fst": fst,
            "seed": seed,
            "definition": (
                "Balding-Nichols model: p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F) per population, "
                "genotypes ~ Binomial(2, p_k)"
            ),
        }, []
    size = opt_number(params, "population_size", 100, minimum=1, maximum=1_000_000, integer=True)
    gens = opt_number(params, "generations", 50, minimum=1, maximum=100_000, integer=True)
    p0 = opt_number(params, "initial_frequency", 0.5, minimum=0, maximum=1)
    reps = opt_number(params, "replicates", 100, minimum=1, maximum=10_000, integer=True)
    if size * gens * reps > 200_000_000:
        raise ConfigurationError(
            "population_size * generations * replicates is too large.",
            code="INPUT_TOO_LARGE",
            hint="Reduce population_size, generations or replicates.",
        )
    finals = []
    trajectories = []
    for _ in range(reps):
        p = p0
        path = [p]
        for _ in range(gens):
            if 0.0 < p < 1.0:
                p = sum(1 for _ in range(2 * size) if rng.random() < p) / (2.0 * size)
            path.append(p)
        finals.append(p)
        trajectories.append(path)
    return {
        "model": "wright_fisher",
        "population_size": size,
        "generations": gens,
        "initial_frequency": p0,
        "final_frequencies": finals,
        "mean_final_frequency": st.mean(finals),
        "variance_final_frequency": st.variance(finals) if reps > 1 else None,
        "expected_variance": p0 * (1 - p0) * (1 - (1 - 1.0 / (2 * size)) ** gens),
        "trajectories": trajectories[:20],
        "seed": seed,
        "definition": (
            "neutral Wright-Fisher drift of a biallelic locus in 2N gene copies: p' = "
            "Binomial(2N, p) / 2N; "
            "E[Var(p_t)] = p0 q0 (1 - (1 - 1/2N)^t)"
        ),
    }, []
