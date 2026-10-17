"""Forensic and kinship statistics from allele frequencies (product rule, NRC II 4.10 theta
correction, paternity and kinship likelihood ratios).  Genotypes are lists of two allele labels;
frequencies are ``{locus: {allele: freq}}``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.gwas._impl_common import (
    load_table,
    nice_ticks,
    opt_bool,
    opt_choice,
    opt_number,
    svg_escape,
    svg_wrap,
    to_number,
    write_text_file,
)

# Cotterman coefficients (k0, k1, k2): probabilities that two non-inbred relatives share
# 0, 1 or 2 alleles identical by descent.
RELATIONSHIPS: dict[str, tuple[float, float, float]] = {
    "unrelated": (1.0, 0.0, 0.0),
    "parent_child": (0.0, 1.0, 0.0),
    "full_sibling": (0.25, 0.5, 0.25),
    "half_sibling": (0.5, 0.5, 0.0),
    "grandparent_grandchild": (0.5, 0.5, 0.0),
    "uncle_nephew": (0.5, 0.5, 0.0),
    "first_cousin": (0.75, 0.25, 0.0),
    "identical_twin": (0.0, 0.0, 1.0),
}

ENFSI_SCALE = (
    (1e6, "extremely strong support"),
    (1e4, "very strong support"),
    (1e3, "strong support"),
    (1e2, "moderately strong support"),
    (10.0, "moderate support"),
    (1.0, "weak or limited support"),
)


# --------------------------------------------------------------------------- parsing
def allele_label(value: object) -> str:
    if isinstance(value, bool) or value is None:
        raise InputFormatError(
            f"invalid allele {value!r}.",
            code="INVALID_ALLELE",
            hint="Alleles are repeat numbers or labels such as 12 or 9.3.",
        )
    if isinstance(value, (int, float)):
        number = float(value)
        return str(int(number)) if number == int(number) else repr(number)
    text = str(value).strip()
    if not text:
        raise InputFormatError(
            "empty allele label.", code="INVALID_ALLELE", hint="Remove empty alleles."
        )
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number == int(number) else repr(number)


def parse_profile(value: object, name: str = "profile") -> dict[str, tuple[str, ...]]:
    """Profile as {locus: (allele1, allele2)}; a locus given with one allele is treated as
    homozygous (and reported by validate).
    """

    if isinstance(value, Mapping) and "records" in value:
        value = value["records"]
    if isinstance(value, Mapping):
        out: dict[str, tuple[str, ...]] = {}
        for locus, alleles in value.items():
            if isinstance(alleles, (str, int, float)) and not isinstance(alleles, bool):
                alleles = [alleles]
            if not isinstance(alleles, Sequence) or not 1 <= len(alleles) <= 2:
                raise InputFormatError(
                    f"{name}: locus {locus!r} needs one or two alleles (single-source profile).",
                    code="INVALID_PROFILE",
                    hint="Use {'D3S1358': [15, 16], ...}; mixtures are not supported.",
                )
            labels = tuple(allele_label(a) for a in alleles)
            out[str(locus)] = labels if len(labels) == 2 else (labels[0], labels[0])
        if not out:
            raise InputFormatError(
                f"{name} is empty.", code="INVALID_PROFILE", hint="Provide at least one locus."
            )
        return out
    records = load_table(value, what=name)
    out = {}
    for rec in records:
        norm = {str("".join(ch for ch in k.lower() if ch.isalnum())): v for k, v in rec.items()}
        locus = norm.get("locus", norm.get("marker"))
        a1, a2 = norm.get("allele1", norm.get("a1")), norm.get("allele2", norm.get("a2"))
        if locus is None or a1 is None:
            raise InputFormatError(
                f"{name}: each row needs locus and allele1 (allele2).",
                code="INVALID_PROFILE",
                hint="Columns: locus, allele1, allele2.",
            )
        out[str(locus)] = (allele_label(a1), allele_label(a2 if a2 is not None else a1))
    return out


def parse_frequencies(value: object, name: str = "frequencies") -> dict[str, dict[str, float]]:
    if (
        isinstance(value, Mapping)
        and value
        and all(isinstance(v, Mapping) for v in value.values())
        and "records" not in value
    ):
        out: dict[str, dict[str, float]] = {}
        for locus, table in value.items():
            freqs: dict[str, float] = {}
            for allele, f in table.items():
                try:
                    number = to_number(f)
                except ValueError as exc:
                    raise InputFormatError(
                        f"{name}: frequency {f!r} is not numeric.",
                        code="INVALID_FREQUENCY",
                        hint="Use numbers in (0, 1].",
                    ) from exc
                if number is None or not 0.0 < number <= 1.0:
                    raise InputFormatError(
                        f"{name}: frequency of {allele!r} at {locus!r} must be in (0, 1].",
                        code="INVALID_FREQUENCY",
                        hint="Use numbers in (0, 1].",
                    )
                freqs[allele_label(allele)] = number
            out[str(locus)] = freqs
        return out
    out = {}
    for rec in load_table(value, what=name):
        norm = {str("".join(ch for ch in k.lower() if ch.isalnum())): v for k, v in rec.items()}
        locus, allele, freq = (
            norm.get("locus", norm.get("marker")),
            norm.get("allele"),
            norm.get("frequency", norm.get("freq")),
        )
        if locus is None or allele is None or freq is None:
            raise InputFormatError(
                f"{name}: each row needs locus, allele and frequency.",
                code="INVALID_FREQUENCY",
                hint="Columns: locus, allele, frequency.",
            )
        number = to_number(freq)
        if number is None or not 0.0 < number <= 1.0:
            raise InputFormatError(
                f"{name}: frequency must be in (0, 1].",
                code="INVALID_FREQUENCY",
                hint="Use numbers in (0, 1].",
            )
        out.setdefault(str(locus), {})[allele_label(allele)] = number
    if not out:
        raise InputFormatError(
            f"{name} is empty.",
            code="INVALID_FREQUENCY",
            hint="Provide locus/allele/frequency rows.",
        )
    return out


class _Freq:
    def __init__(self, table: Mapping[str, Mapping[str, float]], minimum: float | None) -> None:
        self.table, self.minimum = table, minimum

    def get(self, locus: str, allele: str) -> float:
        if locus not in self.table:
            raise ConfigurationError(
                f"no allele frequencies for locus {locus!r}.",
                code="MISSING_FREQUENCIES",
                hint="Add the locus to the frequency table or remove it from the profile.",
            )
        f = self.table[locus].get(allele)
        if f is None:
            if self.minimum is None:
                raise ConfigurationError(
                    f"allele {allele!r} at locus {locus!r} is not in the frequency table.",
                    code="UNSEEN_ALLELE",
                    hint=(
                        "Add it, or set options['min_frequency'] (e.g. 5/(2N), NRC II 4.10 "
                        "minimum)."
                    ),
                )
            return self.minimum
        return f

    def distribution(self, locus: str) -> dict[str, float]:
        return dict(self.table[locus])


def _freq(source: Mapping[str, Any], opts: Mapping[str, Any]) -> _Freq:
    raw = source.get("frequencies", opts.get("frequencies"))
    if raw is None:
        raise InputFormatError(
            "allele 'frequencies' are required.",
            code="MISSING_FREQUENCIES",
            hint="Add {'locus': {'allele': frequency}}.",
        )
    return _Freq(
        parse_frequencies(raw),
        opt_number(opts, "min_frequency", None, minimum=0, maximum=1, exclusive=True),
    )


def _source(value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise InputFormatError(
            "input must be a mapping with profiles and 'frequencies'.",
            code="INVALID_INPUT_TYPE",
            hint="E.g. {'profile': {...}, 'frequencies': {...}}.",
        )
    return value


# --------------------------------------------------------------------------- genotype probabilities
def genotype_probability(
    f: _Freq, locus: str, genotype: Sequence[str], theta: float = 0.0
) -> float:
    """Population frequency of an unordered genotype: p^2 / 2pq, or the NRC II 4.10
    theta-corrected match probability.
    """

    a, b = genotype
    pa, pb = f.get(locus, a), f.get(locus, b)
    if theta == 0.0:
        return pa * pa if a == b else 2 * pa * pb
    denom = (1 + theta) * (1 + 2 * theta)
    if a == b:
        return (2 * theta + (1 - theta) * pa) * (3 * theta + (1 - theta) * pa) / denom
    return 2 * (theta + (1 - theta) * pa) * (theta + (1 - theta) * pb) / denom


def _same(g1: Sequence[str], g2: Sequence[str]) -> bool:
    return sorted(g1) == sorted(g2)


def _log10(x: float) -> float | None:
    return math.log10(x) if x > 0 else None


def verbal_scale(lr: float) -> str:
    if lr == 1.0 or lr <= 0.0 or math.isnan(lr):
        return (
            "no support for either proposition"
            if lr == 1.0
            else "not applicable (LR = 0: exclusion)"
        )
    value = lr if lr > 1 else 1.0 / lr
    side = (
        "the first (prosecution/relationship) proposition"
        if lr > 1
        else "the alternative proposition"
    )
    for bound, label in ENFSI_SCALE:
        if value > bound:
            return f"{label} for {side}"
    return f"weak or limited support for {side}"


# --------------------------------------------------------------------------- statistics
def random_match_probability(source: Mapping[str, Any], opts: Mapping[str, Any]) -> dict[str, Any]:
    f = _freq(source, opts)
    theta = opt_number(opts, "theta", 0.0, minimum=0, maximum=0.5)
    profile = parse_profile(source.get("profile"), "profile") if "profile" in source else None
    if profile is None:
        raise InputFormatError(
            "'profile' is required.",
            code="INVALID_PROFILE",
            hint="Add {'profile': {'locus': [a1, a2]}}.",
        )
    loci = []
    log_total = 0.0
    for locus, genotype in profile.items():
        p = genotype_probability(f, locus, genotype, theta)
        loci.append(
            {
                "locus": locus,
                "genotype": list(genotype),
                "genotype_probability": p,
                "log10": _log10(p),
            }
        )
        log_total += math.log10(p)
    return {
        "quantity": "random_match_probability",
        "theta": theta,
        "per_locus": loci,
        "combined": 10**log_total,
        "combined_log10": log_total,
        "one_in": 10 ** (-log_total),
        "formula": (
            "product rule: homozygote p^2, heterozygote 2pq; with theta>0 the NRC II 4.10 formulas "
            "[2θ+(1-θ)p][3θ+(1-θ)p]/[(1+θ)(1+2θ)] and 2[θ+(1-θ)pi][θ+(1-θ)pj]/[(1+θ)(1+2θ)]"
        ),
    }


def likelihood_ratio(source: Mapping[str, Any], opts: Mapping[str, Any]) -> dict[str, Any]:
    f = _freq(source, opts)
    theta = opt_number(opts, "theta", 0.0, minimum=0, maximum=0.5)
    evidence, suspect = (
        parse_profile(source.get("evidence"), "evidence"),
        parse_profile(source.get("suspect"), "suspect"),
    )
    shared = [loc for loc in evidence if loc in suspect]
    if not shared:
        raise InputFormatError(
            "evidence and suspect profiles share no locus.",
            code="NO_COMMON_LOCI",
            hint="Use the same locus names.",
        )
    loci = []
    log_total, excluded = 0.0, False
    for locus in shared:
        if _same(evidence[locus], suspect[locus]):
            p = genotype_probability(f, locus, evidence[locus], theta)
            lr = 1.0 / p
            loci.append({"locus": locus, "match": True, "lr": lr, "log10_lr": math.log10(lr)})
            log_total += math.log10(lr)
        else:
            excluded = True
            loci.append({"locus": locus, "match": False, "lr": 0.0, "log10_lr": None})
    combined = 0.0 if excluded else 10**log_total
    return {
        "quantity": "likelihood_ratio",
        "theta": theta,
        "per_locus": loci,
        "combined": combined,
        "combined_log10": None if excluded else log_total,
        "excluded": excluded,
        "verbal_equivalent": verbal_scale(combined),
        "formula": (
            "Hp: the suspect is the source of the evidence; Hd: an unrelated random person is the "
            "source. "
            "Per locus LR = 1/P(genotype) when the genotypes match, 0 otherwise; combined = "
            "product "
            "over loci"
        ),
    }


def _transmission(genotype: Sequence[str] | None, dist: Mapping[str, float]) -> dict[str, float]:
    if genotype is None:
        return dict(dist)
    out: dict[str, float] = {}
    for a in genotype:
        out[a] = out.get(a, 0.0) + 0.5
    return out


def _child_probability(
    child: Sequence[str], t_mother: Mapping[str, float], t_father: Mapping[str, float]
) -> float:
    total = 0.0
    for x, px in t_mother.items():
        for y, py in t_father.items():
            if sorted((x, y)) == sorted(child):
                total += px * py
    return total


def paternity(source: Mapping[str, Any], opts: Mapping[str, Any]) -> dict[str, Any]:
    if opt_number(opts, "theta", 0.0, minimum=0, maximum=0.5) != 0.0:
        raise ConfigurationError(
            "theta is only supported for random_match_probability and likelihood_ratio.",
            code="UNSUPPORTED_OPTION",
            hint="Paternity indices use the plain product rule.",
        )
    f = _freq(source, opts)
    child = parse_profile(source.get("child"), "child")
    father = parse_profile(source.get("alleged_father"), "alleged_father")
    mother = (
        parse_profile(source.get("mother"), "mother") if source.get("mother") is not None else None
    )
    prior = opt_number(opts, "prior", 0.5, minimum=0, maximum=1, exclusive=True)
    loci = []
    log_total, exclusions = 0.0, 0
    for locus in child:
        if locus not in father or (mother is not None and locus not in mother):
            continue
        dist = {
            a: f.get(locus, a)
            for a in set(child[locus])
            | set(father[locus])
            | (set(mother[locus]) if mother else set())
        }
        # all alleles need frequencies; a random parent transmits by the full locus distribution
        full = f.distribution(locus)
        for a, v in dist.items():
            full.setdefault(a, v)
        t_father = _transmission(father[locus], full)
        t_random = _transmission(None, full)
        t_mother = _transmission(mother[locus], full) if mother is not None else t_random
        num = _child_probability(child[locus], t_mother, t_father)
        den = _child_probability(child[locus], t_mother, t_random)
        if den == 0.0:
            raise InputFormatError(
                f"the child's genotype at {locus} is inconsistent with the mother's.",
                code="MOTHER_CHILD_INCONSISTENT",
                hint=(
                    "Check the profiles (a mutation would need a mutation model, which is not "
                    "implemented)."
                ),
            )
        pi = num / den
        if pi == 0.0:
            exclusions += 1
        else:
            log_total += math.log10(pi)
        loci.append(
            {"locus": locus, "numerator": num, "denominator": den, "pi": pi, "log10_pi": _log10(pi)}
        )
    if not loci:
        raise InputFormatError(
            "no locus is shared by the child, alleged father (and mother).",
            code="NO_COMMON_LOCI",
            hint="Use the same locus names.",
        )
    cpi = 0.0 if exclusions else 10**log_total
    probability = cpi * prior / (cpi * prior + (1 - prior)) if cpi > 0 else 0.0
    return {
        "quantity": "paternity_index",
        "design": "trio" if mother is not None else "duo",
        "per_locus": loci,
        "combined": cpi,
        "combined_log10": None if exclusions else log_total,
        "exclusions": exclusions,
        "prior": prior,
        "probability_of_paternity": probability,
        "verbal_equivalent": verbal_scale(cpi),
        "formula": (
            "PI = P(child | mother, alleged father) / P(child | mother, random man); each parent "
            "transmits one of its two alleles "
            "with probability 1/2 (a random man transmits allele i with frequency p_i); CPI = "
            "product "
            "of PI; "
            "W = CPI*prior/(CPI*prior + 1 - prior). No mutation model: one exclusion gives CPI = "
            "0; "
            "theta = 0."
        ),
    }


def _conditional(
    f: _Freq, locus: str, g1: Sequence[str], g2: Sequence[str], k: tuple[float, float, float]
) -> float:
    """P(g2 | g1, relationship with IBD probabilities k) at one locus."""

    dist = f.distribution(locus)
    for a in (*g1, *g2):
        dist.setdefault(a, f.get(locus, a))
    p0 = genotype_probability(f, locus, g2)
    p1 = 0.0
    for x in g1:
        for y, q in dist.items():
            if sorted((x, y)) == sorted(g2):
                p1 += 0.5 * q
    p2 = 1.0 if _same(g1, g2) else 0.0
    return k[0] * p0 + k[1] * p1 + k[2] * p2


def kinship(
    source: Mapping[str, Any], opts: Mapping[str, Any], relationship: str | None = None
) -> dict[str, Any]:
    if opt_number(opts, "theta", 0.0, minimum=0, maximum=0.5) != 0.0:
        raise ConfigurationError(
            "theta is only supported for random_match_probability and likelihood_ratio.",
            code="UNSUPPORTED_OPTION",
            hint="Kinship indices use the plain product rule.",
        )
    f = _freq(source, opts)
    a, b = (
        parse_profile(source.get("profile_1"), "profile_1"),
        parse_profile(source.get("profile_2"), "profile_2"),
    )
    ibd = opts.get("ibd")
    if ibd is not None:
        if (
            isinstance(ibd, str)
            or not isinstance(ibd, Sequence)
            or len(ibd) != 3
            or any(not isinstance(v, (int, float)) or v < 0 for v in ibd)
            or abs(sum(ibd) - 1) > 1e-9
        ):
            raise ConfigurationError(
                "ibd must be [k0, k1, k2], non-negative and summing to 1.",
                code="INVALID_OPTION",
                hint="E.g. [0.25, 0.5, 0.25] for full siblings.",
            )
        k, name = (float(ibd[0]), float(ibd[1]), float(ibd[2])), "custom"
    else:
        name = relationship or opt_choice(
            opts, "relationship", "full_sibling", tuple(RELATIONSHIPS)
        )
        k = RELATIONSHIPS[name]
    shared = [loc for loc in a if loc in b]
    if not shared:
        raise InputFormatError(
            "the two profiles share no locus.",
            code="NO_COMMON_LOCI",
            hint="Use the same locus names.",
        )
    loci, log_total = [], 0.0
    zero = False
    for locus in shared:
        num = _conditional(f, locus, a[locus], b[locus], k)
        den = genotype_probability(f, locus, b[locus])
        lr = num / den
        loci.append({"locus": locus, "lr": lr, "log10_lr": _log10(lr)})
        if lr == 0.0:
            zero = True
        else:
            log_total += math.log10(lr)
    combined = 0.0 if zero else 10**log_total
    return {
        "quantity": "kinship_index",
        "relationship": name,
        "ibd": list(k),
        "per_locus": loci,
        "combined": combined,
        "combined_log10": None if zero else log_total,
        "verbal_equivalent": verbal_scale(combined),
        "formula": (
            "LR = P(g2 | g1, relationship) / P(g2) with P(g2 | g1) = k0 P(g2) + k1 P(g2 | one IBD "
            "allele) + k2 [g2 = g1]; "
            "product over loci; theta = 0, no mutation"
        ),
    }


QUANTITIES = ("random_match_probability", "likelihood_ratio", "paternity_index", "kinship_index")


def calculate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    quantity = opt_choice(opts, "quantity", "random_match_probability", QUANTITIES)
    src = _source(source)
    return _by_quantity(quantity, src, opts), []


def _by_quantity(quantity: str, src: Mapping[str, Any], opts: Mapping[str, Any]) -> dict[str, Any]:
    if quantity == "random_match_probability":
        return random_match_probability(src, opts)
    if quantity == "likelihood_ratio":
        return likelihood_ratio(src, opts)
    if quantity == "paternity_index":
        return paternity(src, opts)
    return kinship(src, opts)


# --------------------------------------------------------------------------- compare / infer
def compare(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _source(source)
    a, b = (
        parse_profile(src.get("profile_1"), "profile_1"),
        parse_profile(src.get("profile_2"), "profile_2"),
    )
    shared = [loc for loc in a if loc in b]
    if not shared:
        raise InputFormatError(
            "the two profiles share no locus.",
            code="NO_COMMON_LOCI",
            hint="Use the same locus names.",
        )
    rows = []
    counts = {"identical": 0, "one_allele_shared": 0, "no_allele_shared": 0}
    for locus in shared:
        n_shared = _shared_alleles(a[locus], b[locus])
        status = (
            "identical"
            if _same(a[locus], b[locus])
            else ("one_allele_shared" if n_shared == 1 else "no_allele_shared")
        )
        if status == "identical":
            n_shared = 2
        counts[status] += 1
        rows.append(
            {
                "locus": locus,
                "profile_1": list(a[locus]),
                "profile_2": list(b[locus]),
                "status": status,
                "alleles_shared": n_shared,
            }
        )
    result: dict[str, Any] = {
        "per_locus": rows,
        "n_loci_compared": len(shared),
        **counts,
        "loci_only_in_profile_1": [x for x in a if x not in b],
        "loci_only_in_profile_2": [x for x in b if x not in a],
        "consistent_with_same_source": counts["identical"] == len(shared),
        "mean_alleles_shared": sum(r["alleles_shared"] for r in rows) / len(rows),
    }
    if "frequencies" in src or "frequencies" in opts:
        result["random_match_probability_profile_1"] = random_match_probability(
            {
                "profile": {k: a[k] for k in shared},
                "frequencies": src.get("frequencies", opts.get("frequencies")),
            },
            opts,
        )["combined"]
    return result, []


def _shared_alleles(g1: Sequence[str], g2: Sequence[str]) -> int:
    pool = list(g2)
    n = 0
    for x in g1:
        if x in pool:
            pool.remove(x)
            n += 1
    return n


def infer(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _source(source)
    names = opts.get(
        "relationships",
        ["unrelated", "parent_child", "full_sibling", "half_sibling", "first_cousin"],
    )
    if (
        isinstance(names, str)
        or not isinstance(names, Sequence)
        or not names
        or any(n not in RELATIONSHIPS for n in names)
    ):
        raise ConfigurationError(
            f"relationships must be a list drawn from {sorted(RELATIONSHIPS)}.",
            code="INVALID_OPTION",
            hint="E.g. ['unrelated', 'full_sibling'].",
        )
    priors = opts.get("priors")
    if priors is None:
        prior = {n: 1.0 / len(names) for n in names}
    else:
        if (
            not isinstance(priors, Mapping)
            or set(priors) != set(names)
            or any(not isinstance(v, (int, float)) or v < 0 for v in priors.values())
            or sum(priors.values()) <= 0
        ):
            raise ConfigurationError(
                "priors must map every listed relationship to a non-negative number.",
                code="INVALID_OPTION",
                hint="E.g. {'unrelated': 0.9, 'full_sibling': 0.1}.",
            )
        total = float(sum(priors.values()))
        prior = {n: float(priors[n]) / total for n in names}
    rows = []
    for name in names:
        res = kinship(src, opts, relationship=name)
        rows.append(
            {
                "relationship": name,
                "likelihood_ratio_vs_unrelated": res["combined"],
                "log10_lr": res["combined_log10"],
                "prior": prior[name],
            }
        )
    norm = sum(r["likelihood_ratio_vs_unrelated"] * r["prior"] for r in rows)
    if norm == 0:
        raise InputFormatError(
            "the profiles are inconsistent with every candidate relationship.",
            code="NO_CONSISTENT_RELATIONSHIP",
            hint="Add candidate relationships (including 'unrelated').",
        )
    for r in rows:
        r["posterior"] = r["likelihood_ratio_vs_unrelated"] * r["prior"] / norm
    rows.sort(key=lambda r: -r["posterior"])
    groups: dict[tuple[float, float, float], list[str]] = {}
    for name in names:
        groups.setdefault(RELATIONSHIPS[name], []).append(name)
    warnings = [
        f"relationships {', '.join(group)} have identical IBD probabilities (k0, k1, k2 = {', '.join(f'{x:g}' for x in k)}) and are mathematically "
        "indistinguishable from autosomal genotypes: their likelihood ratios are equal and their posteriors differ only through the priors; "
        "most_probable is one of them by list order, not by evidence"
        for k, group in groups.items()
        if len(group) > 1
    ]
    return {
        "relationships": rows,
        "most_probable": rows[0]["relationship"],
        "method": (
            "posterior(r) ∝ prior(r) * LR(r vs unrelated), LR from the Cotterman IBD model (theta "
            "= 0, no mutation); priors default to uniform"
        ),
    }, warnings


# --------------------------------------------------------------------------- read / validate
def read(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    kind = opt_choice(opts, "kind", "auto", ("auto", "profile", "frequencies"))
    if (
        isinstance(source, Mapping)
        and kind == "auto"
        and ("profile" in source or "frequencies" in source)
    ):
        out: dict[str, Any] = {}
        if "profile" in source:
            out["profile"] = {k: list(v) for k, v in parse_profile(source["profile"]).items()}
        if "frequencies" in source:
            out["frequencies"] = parse_frequencies(source["frequencies"])
        return out, []
    records = load_table(source, what="forensic table")
    columns = {"".join(ch for ch in k.lower() if ch.isalnum()) for k in records[0]}
    if kind == "auto":
        kind = "frequencies" if columns & {"frequency", "freq"} else "profile"
    if kind == "frequencies":
        freqs = parse_frequencies(records)
        return {"kind": "frequencies", "frequencies": freqs, "n_loci": len(freqs)}, []
    profile = parse_profile(records)
    return {
        "kind": "profile",
        "profile": {k: list(v) for k, v in profile.items()},
        "n_loci": len(profile),
    }, []


def validate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _source(source)
    tolerance = opt_number(opts, "tolerance", 0.02, minimum=0, maximum=0.5)
    issues: list[dict[str, str]] = []

    def add(severity: str, where: str, message: str) -> None:
        issues.append({"severity": severity, "where": where, "message": message})

    freqs = None
    if "frequencies" in src:
        freqs = parse_frequencies(src["frequencies"])
        for locus, table in freqs.items():
            total = sum(table.values())
            if abs(total - 1.0) > tolerance:
                add(
                    "error",
                    locus,
                    f"allele frequencies sum to {total:.4f}, not 1 (tolerance {tolerance})",
                )
    for key in [k for k in src if k not in ("frequencies",)]:
        value = src[key]
        raw = value
        if not isinstance(raw, (Mapping, str, Sequence)) or isinstance(raw, bytes):
            continue
        if isinstance(raw, Mapping):
            too_many = False
            for locus, alleles in raw.items():
                n = 1 if isinstance(alleles, (str, int, float)) else len(alleles)
                if n > 2:
                    too_many = True
                    add(
                        "error",
                        f"{key}:{locus}",
                        f"{n} alleles: a single-source profile has at most two",
                    )
                elif n == 1:
                    add(
                        "warning",
                        f"{key}:{locus}",
                        "one allele listed: treated as homozygous (possible dropout)",
                    )
            if too_many:  # the over-complete loci are already reported; check the remaining ones
                value = {
                    locus: a
                    for locus, a in raw.items()
                    if isinstance(a, (str, int, float)) or len(a) <= 2
                }
                if not value:
                    continue
        try:
            profile = parse_profile(value, key)
        except InputFormatError as exc:
            add("error", key, exc.message)
            continue
        if freqs is not None:
            for locus, genotype in profile.items():
                if locus not in freqs:
                    add("error", f"{key}:{locus}", "locus missing from the frequency table")
                    continue
                for allele in set(genotype):
                    if allele not in freqs[locus]:
                        add(
                            "warning",
                            f"{key}:{locus}",
                            f"allele {allele} not in the frequency table (needs "
                            "options['min_frequency'])",
                        )
    errors = [i for i in issues if i["severity"] == "error"]
    return {
        "valid": not errors,
        "n_errors": len(errors),
        "n_warnings": len(issues) - len(errors),
        "issues": issues,
    }, []


# --------------------------------------------------------------------------- report / plot
def _fmt(x: float | None) -> str:
    return "NA" if x is None else (f"{x:.4g}" if x != 0 else "0")


def report(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    quantity = opt_choice(opts, "quantity", "random_match_probability", QUANTITIES)
    res = _by_quantity(quantity, _source(source), opts)
    title = str(opts.get("title", "Forensic statistical report"))
    lines = [f"# {title}", "", f"Statistic: **{quantity.replace('_', ' ')}**", ""]
    if "design" in res:
        lines += [f"Design: {res['design']}", ""]
    lines += [f"Method: {res['formula']}", "", "| Locus | Value | log10 |", "|---|---|---|"]
    key = {
        "random_match_probability": "genotype_probability",
        "likelihood_ratio": "lr",
        "paternity_index": "pi",
        "kinship_index": "lr",
    }[quantity]
    for row in res["per_locus"]:
        lines.append(
            f"| {row['locus']} | {_fmt(row[key])} | "
            f"{_fmt(row.get('log10', row.get('log10_lr', row.get('log10_pi'))))} |"
        )
    lines += ["", f"Combined: **{_fmt(res['combined'])}** (log10 {_fmt(res['combined_log10'])})"]
    if quantity == "random_match_probability":
        lines.append(f"Random match probability about 1 in {_fmt(res['one_in'])}.")
    else:
        lines.append(f"Verbal equivalent (ENFSI 2015 scale): {res['verbal_equivalent']}.")
    if quantity == "paternity_index":
        lines.append(
            f"Probability of paternity (prior {res['prior']}): "
            f"{_fmt(res['probability_of_paternity'])}; exclusions: {res['exclusions']}."
        )
    lines += [
        "",
        "Limitations: single-source profiles, no mutation, drop-in/dropout or mixture model; "
        "allele frequencies as supplied. "
        "Interpretation requires a qualified examiner.",
    ]
    text = "\n".join(lines) + "\n"
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, text, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"format": "markdown", "text": text, "statistic": res, "path": written}, []


def plot(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    quantity = opt_choice(opts, "quantity", "random_match_probability", QUANTITIES)
    res = _by_quantity(quantity, _source(source), opts)
    logs = []
    for row in res["per_locus"]:
        value = (
            row.get("log10")
            if quantity == "random_match_probability"
            else row.get("log10_lr", row.get("log10_pi"))
        )
        logs.append((row["locus"], value))
    width, height = 720, 400
    left, top, right, bottom = 60, 30, 20, 90
    pw, ph = width - left - right, height - top - bottom
    finite = [v for _, v in logs if v is not None]
    lo, hi = min([*finite, 0.0]), max([*finite, 0.0])
    if hi == lo:
        hi = lo + 1.0
    body = [f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="black"/>']
    zero_y = top + ph * (1 - (0 - lo) / (hi - lo))
    body.append(
        f'<line x1="{left}" x2="{left + pw}" y1="{zero_y:.1f}" y2="{zero_y:.1f}" stroke="#7f7f7f"/>'
    )
    for t in nice_ticks(lo, hi, 5):
        y = top + ph * (1 - (t - lo) / (hi - lo))
        body.append(
            f'<text x="{left - 6}" y="{y + 3:.1f}" font-size="10" text-anchor="end">{t:g}</text>'
        )
    step = pw / len(logs)
    for idx, (locus, value) in enumerate(logs):
        x = left + idx * step
        if value is not None:
            y = top + ph * (1 - (value - lo) / (hi - lo))
            body.append(
                f'<rect x="{x + 2:.1f}" y="{min(y, zero_y):.1f}" width="{max(step - 4, 1):.1f}" '
                f'height="{abs(zero_y - y):.1f}" fill="#1f4e79"><title>{svg_escape(locus)}: '
                f"{value:.3f}</title></rect>"
            )
        body.append(
            f'<text x="{x + step / 2:.1f}" y="{height - bottom + 12}" font-size="9" '
            f'text-anchor="end" transform="rotate(-60 {x + step / 2:.1f} '
            f'{height - bottom + 12})">{svg_escape(locus)}</text>'
        )
    label = (
        "log10 genotype probability"
        if quantity == "random_match_probability"
        else "log10 likelihood ratio"
    )
    body.append(
        f'<text x="14" y="{top + ph / 2}" font-size="12" text-anchor="middle" '
        f'transform="rotate(-90 14 {top + ph / 2})">{label}</text>'
    )
    svg = svg_wrap(width, height, body, f"{quantity} per locus")
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, svg, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"quantity": quantity, "svg": svg, "n_loci": len(logs), "path": written}, []
