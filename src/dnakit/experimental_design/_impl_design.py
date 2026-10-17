"""Native experimental-design computations: power / sample size, Lander-Waterman coverage,
randomised layouts (CRD, RCBD, Latin square, factorial), group assignment, design checks, export
and Monte-Carlo power simulation.

All randomisation uses ``random.Random(seed)`` (default seed 0), so results are reproducible.
"""

from __future__ import annotations

import itertools
import json
import math
import random
from collections import Counter, OrderedDict
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.gwas import _impl_stats as st
from dnakit.gwas._impl_common import (
    load_table,
    opt_bool,
    opt_choice,
    opt_number,
    serialize_records,
    write_text_file,
)


def _params(source: object, opts: Mapping[str, Any]) -> dict[str, Any]:
    params: dict[str, Any] = dict(opts)
    if isinstance(source, Mapping):
        params.update(source)
    elif source is not None and not isinstance(source, (int, float, str)):
        raise InputFormatError(
            "input must be None or a mapping of parameters.",
            code="INVALID_INPUT_TYPE",
            hint="Pass a mapping, e.g. {'effect_size': 0.5}.",
        )
    return params


# --------------------------------------------------------------------------- power and sample size
TESTS = ("two_sample_t", "one_sample_t", "paired_t", "two_proportions", "chi_square")


def power_t(effect_size: float, n1: int, n2: int | None, alpha: float, alternative: str) -> float:
    """Exact power of the t-test (non-central t) for Cohen's d; ``n2`` None means a one-sample /
    paired test.
    """

    if n2 is None:
        df, delta = n1 - 1, abs(effect_size) * math.sqrt(n1)
    else:
        df, delta = n1 + n2 - 2, abs(effect_size) * math.sqrt(n1 * n2 / (n1 + n2))
    if alternative == "two-sided":
        return st.nct_power(st.t_ppf(1 - alpha / 2, df), df, delta, two_sided=True)
    return st.nct_power(st.t_ppf(1 - alpha, df), df, delta, two_sided=False)


def power_proportions(
    p1: float, p2: float, n1: float, ratio: float, alpha: float, alternative: str
) -> float:
    """Normal-approximation power for two independent proportions with pooled-variance null (as R
    power.prop.test).
    """

    k = ratio
    pbar = (p1 + k * p2) / (1 + k)
    z = st.normal_ppf(1 - alpha / 2) if alternative == "two-sided" else st.normal_ppf(1 - alpha)
    num = abs(p1 - p2) * math.sqrt(n1) - z * math.sqrt((1 + 1 / k) * pbar * (1 - pbar))
    return st.normal_cdf(num / math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / k))


def n_proportions(
    p1: float, p2: float, power: float, ratio: float, alpha: float, alternative: str
) -> float:
    k = ratio
    pbar = (p1 + k * p2) / (1 + k)
    za = st.normal_ppf(1 - alpha / 2) if alternative == "two-sided" else st.normal_ppf(1 - alpha)
    zb = st.normal_ppf(power)
    return (
        za * math.sqrt((1 + 1 / k) * pbar * (1 - pbar))
        + zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / k)
    ) ** 2 / (p1 - p2) ** 2


def power_chi_square(w: float, df: int, n: float, alpha: float) -> float:
    return st.ncx2_sf(st.chi2_ppf(1 - alpha, df), df, n * w * w)


def _smallest_n(power_fn: Any, target: float, start: int = 2, limit: int = 10_000_000) -> int:
    lo, hi = start, start
    while power_fn(hi) < target:
        lo, hi = hi, hi * 2
        if hi > limit:
            raise ConfigurationError(
                "the required sample size exceeds 10,000,000.",
                code="EFFECT_TOO_SMALL",
                hint="Use a larger effect size.",
            )
    if power_fn(lo) >= target:
        return lo
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if power_fn(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi


def calculate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    p = _params(source, opts)
    quantity = opt_choice(
        p, "quantity", "sample_size", ("sample_size", "power", "coverage", "replicates")
    )
    if quantity == "coverage":
        return lander_waterman(p), []
    if quantity == "replicates":
        return replicate_count(p), []
    test = opt_choice(p, "test", "two_sample_t", TESTS)
    alpha = opt_number(p, "alpha", 0.05, minimum=0, maximum=1, exclusive=True)
    alt = opt_choice(p, "alternative", "two-sided", ("two-sided", "larger", "smaller"))
    alt = "two-sided" if alt == "two-sided" else "one-sided"
    target = opt_number(p, "power", 0.8, minimum=0, maximum=1, exclusive=True)
    ratio = opt_number(p, "ratio", 1.0, minimum=0, exclusive=True)
    n_given = opt_number(p, "n", None, minimum=2, integer=True)
    result: dict[str, Any] = {
        "quantity": quantity,
        "test": test,
        "alpha": alpha,
        "alternative": alt,
    }
    if test in ("two_sample_t", "one_sample_t", "paired_t"):
        d = opt_number(p, "effect_size", None, minimum=0, exclusive=True)
        if d is None:
            raise InputFormatError(
                "'effect_size' (Cohen's d = mean difference / SD) is required.",
                code="MISSING_VALUES",
                hint="E.g. {'effect_size': 0.5}.",
            )
        two = test == "two_sample_t"
        n2_of = (lambda n: max(2, math.ceil(ratio * n))) if two else (lambda n: None)
        if quantity == "power":
            if n_given is None:
                raise InputFormatError(
                    "'n' (per group, or pairs/subjects) is required to compute power.",
                    code="MISSING_VALUES",
                    hint="Add 'n'.",
                )
            result.update(n=n_given, power=power_t(d, n_given, n2_of(n_given), alpha, alt))
        else:
            n = _smallest_n(lambda m: power_t(d, m, n2_of(m), alpha, alt), target)
            result.update(
                target_power=target,
                n_per_group=n if two else None,
                n=n,
                n_group2=n2_of(n),
                total=(n + n2_of(n)) if two else n,
                achieved_power=power_t(d, n, n2_of(n), alpha, alt),
            )
        result["method"] = (
            "exact non-central t power; delta = d*sqrt(n1*n2/(n1+n2)) (two-sample) or d*sqrt(n) "
            "(one-sample/paired), df = n1+n2-2 or n-1"
        )
    elif test == "two_proportions":
        p1, p2 = (
            opt_number(p, "p1", None, minimum=0, maximum=1, exclusive=True),
            opt_number(p, "p2", None, minimum=0, maximum=1, exclusive=True),
        )
        if p1 is None or p2 is None or p1 == p2:
            raise InputFormatError(
                "'p1' and 'p2' (different, in (0,1)) are required.",
                code="MISSING_VALUES",
                hint="E.g. {'p1': 0.5, 'p2': 0.75}.",
            )
        if quantity == "power":
            if n_given is None:
                raise InputFormatError(
                    "'n' (per group) is required to compute power.",
                    code="MISSING_VALUES",
                    hint="Add 'n'.",
                )
            result.update(n=n_given, power=power_proportions(p1, p2, n_given, ratio, alpha, alt))
        else:
            exact = n_proportions(p1, p2, target, ratio, alpha, alt)
            n = math.ceil(exact - 1e-9)
            result.update(
                target_power=target,
                n_per_group=n,
                n_group2=math.ceil(ratio * n),
                n_unrounded=exact,
                achieved_power=power_proportions(p1, p2, n, ratio, alpha, alt),
            )
        result["method"] = (
            "normal approximation with pooled null variance (R power.prop.test): n = "
            "(z_a*sqrt((1+1/k) pbar qbar) + "
            "z_b*sqrt(p1 q1 + p2 q2 / k))^2 / (p1-p2)^2"
        )
    else:
        w = opt_number(p, "effect_size", None, minimum=0, exclusive=True)
        df = opt_number(p, "df", None, minimum=1, integer=True)
        if w is None or df is None:
            raise InputFormatError(
                "'effect_size' (Cohen's w) and 'df' are required for the chi-square test.",
                code="MISSING_VALUES",
                hint="E.g. {'effect_size': 0.3, 'df': 1}.",
            )
        if quantity == "power":
            if n_given is None:
                raise InputFormatError(
                    "'n' (total observations) is required to compute power.",
                    code="MISSING_VALUES",
                    hint="Add 'n'.",
                )
            result.update(n=n_given, power=power_chi_square(w, df, n_given, alpha))
        else:
            n = _smallest_n(lambda m: power_chi_square(w, df, m, alpha), target, start=2)
            result.update(
                target_power=target, n_total=n, achieved_power=power_chi_square(w, df, n, alpha)
            )
        result["method"] = (
            "non-central chi-square with lambda = n*w^2 and the central chi-square critical value "
            "at 1-alpha"
        )
        result["alternative"] = "upper tail (chi-square tests are one-sided by construction)"
    return result, []


def lander_waterman(p: Mapping[str, Any]) -> dict[str, Any]:
    """Lander & Waterman (1988) random-clone coverage statistics (also the Poisson model for
    sequencing depth).
    """

    genome = opt_number(p, "genome_size", None, minimum=0, exclusive=True)
    length = opt_number(p, "read_length", None, minimum=0, exclusive=True)
    reads = opt_number(p, "n_reads", None, minimum=0, exclusive=True)
    cov = opt_number(p, "coverage", None, minimum=0, exclusive=True)
    target_fraction = opt_number(p, "target_fraction", None, minimum=0, maximum=1, exclusive=True)
    min_depth = opt_number(p, "min_depth", 1, minimum=1, integer=True)
    theta = opt_number(p, "overlap_fraction", 0.0, minimum=0, maximum=1)
    if cov is None:
        if genome is None or length is None or reads is None:
            raise InputFormatError(
                "give coverage, or genome_size with read_length and n_reads.",
                code="MISSING_VALUES",
                hint="E.g. {'genome_size': 5e6, 'read_length': 150, 'n_reads': 1e6}.",
            )
        cov = length * reads / genome
    result: dict[str, Any] = {
        "quantity": "coverage",
        "coverage": cov,
        "fraction_covered": 1 - math.exp(-cov),
        "fraction_covered_min_depth": 1
        - sum(math.exp(-cov) * cov**i / math.factorial(i) for i in range(min_depth)),
        "min_depth": min_depth,
    }
    if genome is not None:
        result["expected_uncovered_bases"] = genome * math.exp(-cov)
    if length is not None and reads is not None:
        sigma = 1 - theta
        result["expected_contigs"] = reads * math.exp(-cov * sigma)
        result["expected_contig_length"] = length * ((math.exp(cov * sigma) - 1) / cov + 1 - sigma)
    if target_fraction is not None:
        needed = -math.log(1 - target_fraction)
        result["coverage_needed_for_target_fraction"] = needed
        if genome is not None and length is not None:
            result["reads_needed_for_target_fraction"] = math.ceil(needed * genome / length)
    result["method"] = (
        "Lander-Waterman: C = L*N/G, covered fraction 1 - exp(-C), expected contigs "
        "N*exp(-C*(1-theta)), "
        "contig length L*((exp(C(1-theta)) - 1)/C + 1 - theta); depth >= k fraction from the "
        "Poisson(C) distribution"
    )
    return result


def replicate_count(p: Mapping[str, Any]) -> dict[str, Any]:
    sd = opt_number(p, "sd", None, minimum=0, exclusive=True)
    margin = opt_number(p, "margin", None, minimum=0, exclusive=True)
    conf = opt_number(p, "confidence", 0.95, minimum=0, maximum=1, exclusive=True)
    if sd is None or margin is None:
        raise InputFormatError(
            "'sd' and 'margin' (half-width of the confidence interval) are required.",
            code="MISSING_VALUES",
            hint="E.g. {'sd': 10, 'margin': 5}.",
        )
    z = st.normal_ppf(1 - (1 - conf) / 2)
    n_normal = math.ceil((z * sd / margin) ** 2 - 1e-12)
    n = max(2, n_normal)
    while st.t_ppf(1 - (1 - conf) / 2, n - 1) * sd / math.sqrt(n) > margin:
        n += 1
    while n > 2 and st.t_ppf(1 - (1 - conf) / 2, n - 2) * sd / math.sqrt(n - 1) <= margin:
        n -= 1
    return {
        "quantity": "replicates",
        "replicates": n,
        "replicates_normal_approximation": n_normal,
        "sd": sd,
        "margin": margin,
        "confidence": conf,
        "method": (
            "smallest n with t(1-(1-conf)/2, n-1) * sd / sqrt(n) <= margin (the z version is "
            "reported as replicates_normal_approximation)"
        ),
    }


# --------------------------------------------------------------------------- layouts
def _treatments(p: Mapping[str, Any]) -> list[str]:
    t = p.get("treatments")
    if t is None:
        k = opt_number(p, "n_treatments", None, minimum=2, maximum=1000, integer=True)
        if k is None:
            raise InputFormatError(
                "give 'treatments' (list) or 'n_treatments'.",
                code="MISSING_VALUES",
                hint="E.g. {'treatments': ['A', 'B', 'C']}.",
            )
        return [f"T{i + 1}" for i in range(k)]
    if isinstance(t, str) or not isinstance(t, Sequence) or len(t) < 2:
        raise InputFormatError(
            "'treatments' must be a list of at least two labels.",
            code="INVALID_TREATMENTS",
            hint="E.g. ['A', 'B', 'C'].",
        )
    labels = [str(x) for x in t]
    if len(set(labels)) != len(labels):
        raise InputFormatError(
            "treatment labels must be unique.", code="INVALID_TREATMENTS", hint="Remove duplicates."
        )
    return labels


def _factorial_levels(factors: object) -> OrderedDict[str, list[str]]:
    if not isinstance(factors, Mapping) or not factors:
        raise InputFormatError(
            "'factors' must map factor names to level lists.",
            code="INVALID_FACTORS",
            hint="E.g. {'A': [0, 1], 'B': ['x', 'y']}.",
        )
    out: OrderedDict[str, list[str]] = OrderedDict()
    for name, levels in factors.items():
        if (
            isinstance(levels, str)
            or not isinstance(levels, Sequence)
            or len(levels) < 2
            or len({str(v) for v in levels}) != len(levels)
        ):
            raise InputFormatError(
                f"factor {name!r} needs at least two unique levels.",
                code="INVALID_FACTORS",
                hint="E.g. 'A': [0, 1].",
            )
        out[str(name)] = [str(v) for v in levels]
    return out


def design(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    p = _params(source, opts)
    kind = opt_choice(p, "design_type", "crd", ("crd", "rcbd", "latin_square", "factorial"))
    seed = opt_number(p, "seed", 0, integer=True)
    rng = random.Random(seed)
    layout: list[dict[str, Any]] = []
    if kind == "latin_square":
        t = _treatments(p)
        k = len(t)
        rows, cols, labels = list(range(k)), list(range(k)), list(range(k))
        rng.shuffle(rows)
        rng.shuffle(cols)
        rng.shuffle(labels)
        plot = 1
        for r in range(k):
            for c in range(k):
                layout.append(
                    {
                        "plot": plot,
                        "row": r + 1,
                        "column": c + 1,
                        "treatment": t[labels[(rows[r] + cols[c]) % k]],
                    }
                )
                plot += 1
        definition = "cyclic Latin square with rows, columns and treatment labels randomly permuted"
    else:
        if kind == "factorial":
            levels = _factorial_levels(p.get("factors"))
            names = list(levels)
            t = [
                "|".join(f"{n}={v}" for n, v in zip(names, combo, strict=False))
                for combo in itertools.product(*levels.values())
            ]
            combos = list(itertools.product(*levels.values()))
        else:
            t = _treatments(p)
            combos = None
        reps = opt_number(p, "replicates", 2, minimum=1, maximum=1000, integer=True)
        if len(t) * reps > 100_000:
            raise ConfigurationError(
                "the design would have more than 100,000 plots.",
                code="INPUT_TOO_LARGE",
                hint="Reduce treatments or replicates.",
            )
        index = {lab: i for i, lab in enumerate(t)}
        if kind == "rcbd" or (kind == "factorial" and p.get("blocks") is not None):
            blocks = opt_number(p, "blocks", reps, minimum=1, maximum=1000, integer=True)
            plot = 1
            for b in range(blocks):
                order = list(t)
                rng.shuffle(order)
                for lab in order:
                    row: dict[str, Any] = {"plot": plot, "block": b + 1, "treatment": lab}
                    if combos is not None:
                        row.update(dict(zip(levels, combos[index[lab]], strict=False)))
                    layout.append(row)
                    plot += 1
            definition = (
                "each block holds every treatment once in random order (randomised complete block)"
            )
        else:
            units = [lab for lab in t for _ in range(reps)]
            rng.shuffle(units)
            for plot, lab in enumerate(units, start=1):
                row = {"plot": plot, "treatment": lab}
                if combos is not None:
                    row.update(dict(zip(levels, combos[index[lab]], strict=False)))
                layout.append(row)
            definition = (
                "treatments x replicates randomly allocated to plots (completely randomised)"
            )
    counts = Counter(r["treatment"] for r in layout)
    return {
        "design_type": kind,
        "seed": seed,
        "n_plots": len(layout),
        "layout": layout,
        "replications": dict(counts),
        "definition": definition,
    }, []


# --------------------------------------------------------------------------- generate
def generate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    p = _params(source, opts)
    kind = opt_choice(
        p, "kind", "full_factorial", ("full_factorial", "fractional_factorial", "latin_square")
    )
    if kind == "full_factorial":
        levels = _factorial_levels(p.get("factors"))
        total = math.prod(len(v) for v in levels.values())
        if total > 100_000:
            raise ConfigurationError(
                "the full factorial would have more than 100,000 runs.",
                code="INPUT_TOO_LARGE",
                hint="Use fewer factors/levels.",
            )
        names = list(levels)
        runs = []
        for i in range(total):
            idx, row = i, {"run": i + 1}
            for name in names:  # first factor varies fastest (standard / Yates order)
                k = len(levels[name])
                row[name] = levels[name][idx % k]
                idx //= k
            runs.append(row)
        return {
            "kind": kind,
            "n_runs": total,
            "runs": runs,
            "order": "standard (Yates): first factor varies fastest",
        }, []
    if kind == "latin_square":
        t = _treatments(p)
        k = len(t)
        return {
            "kind": kind,
            "square": [[t[(r + c) % k] for c in range(k)] for r in range(k)],
            "order": "cyclic, unrandomised",
        }, []
    n = opt_number(p, "n_factors", None, minimum=2, maximum=12, integer=True)
    gens = p.get("generators", [])
    if n is None or isinstance(gens, str) or not isinstance(gens, Sequence):
        raise InputFormatError(
            "give n_factors and generators for a two-level fractional factorial.",
            code="MISSING_VALUES",
            hint="E.g. {'n_factors': 4, 'generators': ['D=ABC']} (a 2^(4-1) design).",
        )
    letters = [chr(ord("A") + i) for i in range(n)]
    parsed = []
    for g in gens:
        try:
            left, right = str(g).replace(" ", "").upper().split("=")
        except ValueError as exc:
            raise InputFormatError(
                f"invalid generator {g!r}.", code="INVALID_GENERATOR", hint="Use e.g. 'D=ABC'."
            ) from exc
        if left not in letters or not right or any(ch not in letters or ch == left for ch in right):
            raise InputFormatError(
                f"invalid generator {g!r}.",
                code="INVALID_GENERATOR",
                hint=f"Use letters from {''.join(letters)}.",
            )
        parsed.append((left, right))
    generated = {left for left, _ in parsed}
    if len(generated) != len(parsed):
        raise InputFormatError(
            "each generated factor may appear once.",
            code="INVALID_GENERATOR",
            hint="Remove duplicates.",
        )
    base = [x for x in letters if x not in generated]
    if any(ch in generated for _, right in parsed for ch in right):
        raise InputFormatError(
            "generators must be expressed in base factors only.",
            code="INVALID_GENERATOR",
            hint="Use base factors on the right-hand side.",
        )
    runs = []
    for i in range(2 ** len(base)):
        row = {"run": i + 1}
        coded = {}
        for pos, name in enumerate(base):
            coded[name] = -1 if (i >> pos) & 1 == 0 else 1
        for left, right in parsed:
            coded[left] = math.prod(coded[ch] for ch in right)
        row.update({x: coded[x] for x in letters})
        runs.append(row)
    return {
        "kind": kind,
        "n_runs": len(runs),
        "runs": runs,
        "generators": [f"{a}={b}" for a, b in parsed],
        "order": (
            "standard order of the base factors (first varies fastest); generated columns are "
            "products of base columns"
        ),
    }, []


# --------------------------------------------------------------------------- randomize / assign
def randomize(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    seed = opt_number(opts, "seed", 0, integer=True)
    if (
        isinstance(source, Sequence)
        and not isinstance(source, (str, bytes))
        and source
        and not isinstance(source[0], Mapping)
    ):
        runs: list[dict[str, Any]] = [{"item": x} for x in source]
    else:
        runs = load_table(source, what="runs")
    if len(runs) > 1_000_000 or not runs:
        raise InputFormatError(
            "provide between 1 and 1,000,000 runs.",
            code="INVALID_INPUT_SIZE",
            hint="Pass a list of runs or a table.",
        )
    rng = random.Random(seed)
    within = opts.get("within")
    if within is not None:
        if any(within not in r for r in runs):
            raise ConfigurationError(
                f"column {within!r} is missing from some runs.",
                code="INVALID_OPTION",
                hint="Use a column present in every run.",
            )
        groups: OrderedDict[Any, list[dict[str, Any]]] = OrderedDict()
        for r in runs:
            groups.setdefault(r[within], []).append(r)
        ordered = []
        for members in groups.values():
            members = list(members)
            rng.shuffle(members)
            ordered.extend(members)
    else:
        ordered = list(runs)
        rng.shuffle(ordered)
    out = [dict(r, run_order=i + 1) for i, r in enumerate(ordered)]
    return {
        "seed": seed,
        "n_runs": len(out),
        "within": within,
        "runs": out,
        "method": (
            "Fisher-Yates shuffle (Python random.Random(seed)); with `within`, runs are shuffled "
            "inside each group and groups keep first-appearance order"
        ),
    }, []


def assign(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    seed = opt_number(opts, "seed", 0, integer=True)
    if isinstance(source, int) and not isinstance(source, bool):
        if not 1 <= source <= 1_000_000:
            raise InputFormatError(
                "the number of units must be between 1 and 1,000,000.",
                code="INVALID_INPUT_SIZE",
                hint="Pass a positive integer or a list.",
            )
        units: list[dict[str, Any]] = [{"unit": f"u{i + 1}"} for i in range(source)]
    elif (
        isinstance(source, Sequence)
        and not isinstance(source, (str, bytes))
        and source
        and not isinstance(source[0], Mapping)
    ):
        units = [{"unit": str(x)} for x in source]
    else:
        units = load_table(source, what="units")
        if not units:
            raise InputFormatError(
                "no units given.", code="EMPTY_TABLE", hint="Provide units to assign."
            )
    groups = opts.get("groups")
    if groups is None:
        k = opt_number(opts, "n_groups", 2, minimum=2, maximum=100, integer=True)
        groups = [f"G{i + 1}" for i in range(k)]
    if (
        isinstance(groups, str)
        or not isinstance(groups, Sequence)
        or len(groups) < 2
        or len(set(map(str, groups))) != len(groups)
    ):
        raise ConfigurationError(
            "groups must be a list of at least two unique names.",
            code="INVALID_OPTION",
            hint="E.g. ['control', 'treated'].",
        )
    groups = [str(g) for g in groups]
    ratio = opts.get("ratio", [1] * len(groups))
    if (
        isinstance(ratio, str)
        or not isinstance(ratio, Sequence)
        or len(ratio) != len(groups)
        or any(isinstance(r, bool) or not isinstance(r, int) or r < 1 for r in ratio)
    ):
        raise ConfigurationError(
            "ratio must be a list of positive integers, one per group.",
            code="INVALID_OPTION",
            hint="E.g. [1, 2] for twice as many in the second group.",
        )
    pattern = [g for g, r in zip(groups, ratio, strict=False) for _ in range(r)]
    strat = opts.get("stratify_by")
    if strat is not None and any(strat not in u for u in units):
        raise ConfigurationError(
            f"column {strat!r} is missing from some units.",
            code="INVALID_OPTION",
            hint="Use a column present in every unit.",
        )
    rng = random.Random(seed)
    strata: OrderedDict[Any, list[int]] = OrderedDict()
    for i, u in enumerate(units):
        strata.setdefault(u[strat] if strat is not None else None, []).append(i)
    assigned: dict[int, str] = {}
    for members in strata.values():
        rng.shuffle(members)
        for start in range(0, len(members), len(pattern)):
            block = list(pattern)
            rng.shuffle(block)
            for idx, label in zip(members[start : start + len(pattern)], block, strict=False):
                assigned[idx] = label
    rows = [dict(u, group=assigned[i]) for i, u in enumerate(units)]
    counts = Counter(r["group"] for r in rows)
    return {
        "seed": seed,
        "n_units": len(rows),
        "assignments": rows,
        "group_sizes": {g: counts.get(g, 0) for g in groups},
        "method": (
            "permuted-block randomisation: within each stratum, consecutive blocks of sum(ratio) "
            "units receive a random permutation of the group pattern, "
            "so group sizes differ from the target ratio by less than one block"
        ),
    }, []


# ------------------------------------------------------------------------ check / validate / export
def _design_table(source: object) -> list[dict[str, Any]]:
    # accept the whole envelope returned by design() ({'status', 'result', ...}) or its inner result mapping
    if isinstance(source, Mapping) and "status" in source and "result" in source:
        if source["status"] != "ok":
            raise InputFormatError(
                f"the design envelope has status {source['status']!r}, not 'ok'; there is no layout to use.",
                code="INVALID_DESIGN_ENVELOPE",
                hint="Pass the envelope of a successful design() call, its 'result' mapping, or the layout records.",
            )
        source = source["result"]
    if isinstance(source, Mapping) and isinstance(source.get("layout"), Sequence):
        return load_table(source["layout"], what="design")
    return load_table(source, what="design")


def check(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows = _design_table(source)
    tcol = str(opts.get("treatment_column", "treatment"))
    if not rows or any(tcol not in r for r in rows):
        raise InputFormatError(
            f"every row needs a '{tcol}' value.",
            code="MISSING_COLUMNS",
            hint="Set options['treatment_column'] if it has another name.",
        )
    reps = Counter(str(r[tcol]) for r in rows)
    findings: dict[str, Any] = {
        "n_plots": len(rows),
        "n_treatments": len(reps),
        "replications": dict(reps),
        "balanced": len(set(reps.values())) == 1,
    }
    issues = []
    if not findings["balanced"]:
        issues.append("treatments have unequal replication")
    if any(c < 2 for c in reps.values()):
        issues.append("a treatment is not replicated (no error df for that treatment)")
    if "block" in rows[0]:
        per_block: dict[Any, Counter] = {}
        for r in rows:
            per_block.setdefault(r["block"], Counter())[str(r[tcol])] += 1
        complete = all(
            set(c) == set(reps) and len(set(c.values())) == 1 for c in per_block.values()
        )
        findings["n_blocks"] = len(per_block)
        findings["complete_blocks"] = complete
        if not complete:
            issues.append("not every block contains every treatment equally often")
    if "row" in rows[0] and "column" in rows[0]:
        ok = True
        for key in ("row", "column"):
            groups: dict[Any, list[str]] = {}
            for r in rows:
                groups.setdefault(r[key], []).append(str(r[tcol]))
            ok = ok and all(sorted(v) == sorted(set(reps)) for v in groups.values())
        findings["latin_square"] = ok and len({r["row"] for r in rows}) == len(reps)
        if not findings["latin_square"]:
            issues.append(
                "treatments do not appear exactly once per row and column (not a Latin square)"
            )
    factor_cols = opts.get("factor_columns")
    if factor_cols is not None:
        if (
            isinstance(factor_cols, str)
            or not isinstance(factor_cols, Sequence)
            or any(c not in rows[0] for c in factor_cols)
        ):
            raise ConfigurationError(
                "factor_columns must name numeric (e.g. -1/+1) columns of the design.",
                code="INVALID_OPTION",
                hint="E.g. ['A', 'B'].",
            )
        try:
            cols = {c: [float(r[c]) for r in rows] for c in factor_cols}
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "factor_columns must be numeric.",
                code="INVALID_OPTION",
                hint="Use coded -1/+1 columns.",
            ) from exc
        dots = {}
        for a, b in itertools.combinations(factor_cols, 2):
            dots[f"{a}*{b}"] = math.fsum(x * y for x, y in zip(cols[a], cols[b], strict=False))
        findings["column_dot_products"] = dots
        findings["orthogonal"] = all(abs(v) < 1e-9 for v in dots.values())
        if not findings["orthogonal"]:
            issues.append("factor columns are not orthogonal")
    findings["issues"] = issues
    findings["ok"] = not issues
    return findings, []


def validate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows = _design_table(source)
    required = opts.get("required_columns", ["treatment"])
    if isinstance(required, str) or not isinstance(required, Sequence):
        raise ConfigurationError(
            "required_columns must be a list of column names.",
            code="INVALID_OPTION",
            hint="E.g. ['plot', 'treatment'].",
        )
    errors = []
    if not rows:
        errors.append("the design has no rows")
    else:
        for col in required:
            missing = sum(1 for r in rows if r.get(col) is None or str(r.get(col)).strip() == "")
            if col not in rows[0]:
                errors.append(f"required column '{col}' is missing")
            elif missing:
                errors.append(f"column '{col}' has {missing} empty value(s)")
        unit = opts.get("unit_column", "plot")
        if unit in rows[0]:
            ids = [r[unit] for r in rows]
            if len(set(ids)) != len(ids):
                errors.append(f"unit column '{unit}' has duplicated values")
        widths = {len(r) for r in rows}
        if len(widths) > 1:
            errors.append("rows have different numbers of columns")
    return {"valid": not errors, "n_rows": len(rows), "errors": errors}, []


def export(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows = _design_table(source)
    if not rows:
        raise InputFormatError(
            "the design has no rows.", code="EMPTY_TABLE", hint="Provide a layout."
        )
    fmt = opt_choice(opts, "format", "csv", ("csv", "tsv", "json"))
    columns = list(opts.get("columns", rows[0].keys()))
    if fmt == "json":
        text = json.dumps([{c: r.get(c) for c in columns} for r in rows], indent=2)
    else:
        text = serialize_records(rows, columns, "," if fmt == "csv" else "\t", missing="")
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, text, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {
        "format": fmt,
        "n_rows": len(rows),
        "columns": columns,
        "text": text,
        "path": written,
    }, []


# --------------------------------------------------------------------------- simulate
def simulate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    p = _params(source, opts)
    test = opt_choice(p, "test", "two_sample_t", ("two_sample_t", "one_sample_t"))
    n = opt_number(p, "n", None, minimum=2, maximum=100_000, integer=True)
    d = opt_number(p, "effect_size", None, minimum=0)
    if n is None or d is None:
        raise InputFormatError(
            "'n' (per group) and 'effect_size' (Cohen's d) are required.",
            code="MISSING_VALUES",
            hint="E.g. {'n': 20, 'effect_size': 0.8}.",
        )
    alpha = opt_number(p, "alpha", 0.05, minimum=0, maximum=1, exclusive=True)
    sims = opt_number(p, "n_simulations", 2000, minimum=10, maximum=200_000, integer=True)
    seed = opt_number(p, "seed", 0, integer=True)
    if sims * n > 20_000_000:
        raise ConfigurationError(
            "n_simulations * n is too large.",
            code="INPUT_TOO_LARGE",
            hint="Reduce n_simulations or n.",
        )
    rng = random.Random(seed)
    hits = 0
    for _ in range(sims):
        if test == "two_sample_t":
            a = [rng.gauss(0.0, 1.0) for _ in range(n)]
            b = [rng.gauss(d, 1.0) for _ in range(n)]
            ma, mb = st.mean(a), st.mean(b)
            sp2 = (st.variance(a) + st.variance(b)) / 2.0
            t = (mb - ma) / math.sqrt(sp2 * 2.0 / n)
            df = 2 * n - 2
        else:
            x = [rng.gauss(d, 1.0) for _ in range(n)]
            t = st.mean(x) / math.sqrt(st.variance(x) / n)
            df = n - 1
        if st.t_two_sided_p(t, df) < alpha:
            hits += 1
    power = hits / sims
    exact = (
        power_t(d, n, n if test == "two_sample_t" else None, alpha, "two-sided") if d > 0 else alpha
    )
    return {
        "test": test,
        "n": n,
        "effect_size": d,
        "alpha": alpha,
        "n_simulations": sims,
        "seed": seed,
        "empirical_power": power,
        "standard_error": math.sqrt(power * (1 - power) / sims),
        "analytic_power": exact,
        "method": (
            "Monte Carlo: normal data with unit SD, two-sided pooled t-test; compare with the "
            "exact non-central-t power"
        ),
    }, []
