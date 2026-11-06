"""Individual-level and LD-based GWAS computations: association scans, simulation, conditional/joint
analysis, polygenic scores and SNP heritability / genetic correlation estimators.
"""

from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

from . import _impl_stats as st
from ._impl_common import (
    CandidateUnavailable,
    opt_bool,
    opt_choice,
    opt_number,
    parse_genotype_matrix,
    require_option,
    to_number,
)
from ._impl_sumstats import (
    _need,
    chrom_order,
    load_sumstats,
)

_EXTERNAL_PROGRAMS = {
    "plink": "plink2",
    "plink2": "plink2",
    "plink1.9": "plink",
    "regenie": "regenie",
    "gemma": "gemma",
    "gcta": "gcta64",
    "gcta64": "gcta64",
}


def _numpy() -> Any:
    try:
        import numpy

        return numpy
    except ImportError:  # pragma: no cover - exercised only without numpy
        return None


# --------------------------------------------------------------------------- correlation helpers
def column_r(a: Sequence[float | None], b: Sequence[float | None]) -> float | None:
    """Pearson correlation over the samples observed in both columns (None if < 3 or no
    variance).
    """

    pairs = [(x, y) for x, y in zip(a, b, strict=False) if x is not None and y is not None]
    if len(pairs) < 3:
        return None
    r = st.pearson([p[0] for p in pairs], [p[1] for p in pairs])
    return None if math.isnan(r) else r


def correlation_matrix_from_options(
    opts: Mapping[str, Any], wanted: Sequence[str] | None = None, limit: int = 500
) -> tuple[list[str], list[list[float]]]:
    """LD correlation matrix from ``options['ld']`` ({'variant_ids', 'matrix'}) or
    ``options['genotypes']``.
    """

    ld = opts.get("ld")
    if ld is not None:
        if not isinstance(ld, Mapping):
            raise ConfigurationError(
                "options['ld'] must be {'variant_ids': [...], 'matrix': [[...]]}.",
                code="INVALID_OPTION",
                hint="Provide a square correlation (r) matrix with matching variant ids.",
            )
        ids = ld.get("variant_ids", ld.get("ids"))
        matrix = ld.get("matrix")
        if (
            not isinstance(ids, Sequence)
            or not isinstance(matrix, Sequence)
            or len(ids) != len(matrix)
            or any(not isinstance(r, Sequence) or len(r) != len(ids) for r in matrix)
        ):
            raise InputFormatError(
                "the LD matrix must be square with one variant id per row.",
                code="INVALID_LD",
                hint=(
                    "Use {'variant_ids': [...], 'matrix': [[1, r12], [r12, 1]]} holding "
                    "correlations r."
                ),
            )
        ids = [str(i) for i in ids]
        mat = [[float(v) for v in row] for row in matrix]
        for i in range(len(ids)):
            for j in range(len(ids)):
                if abs(mat[i][j] - mat[j][i]) > 1e-6 or abs(mat[i][j]) > 1 + 1e-9:
                    raise InputFormatError(
                        "the LD matrix must be symmetric with entries in [-1, 1].",
                        code="INVALID_LD",
                        hint="Provide the signed correlation r, not r^2.",
                    )
        return ids, mat
    genotypes = opts.get("genotypes")
    if genotypes is None:
        raise ConfigurationError(
            "a linkage-disequilibrium source is required.",
            code="MISSING_OPTION",
            hint=(
                "Pass options['ld'] = {'variant_ids', 'matrix'} or options['genotypes'] "
                "(reference panel)."
            ),
        )
    g = parse_genotype_matrix(genotypes, what="LD reference genotypes", min_samples=3)
    ids = g["variants"]
    keep = [i for i, v in enumerate(ids) if wanted is None or v in set(wanted)]
    if len(keep) > limit:
        raise ConfigurationError(
            f"{len(keep)} variants exceed the LD limit of {limit}.",
            code="TOO_MANY_VARIANTS",
            hint="Restrict to a locus or pass fewer variants.",
        )
    cols = {j: [row[j] for row in g["matrix"]] for j in keep}
    mat = [[1.0 if i == j else (column_r(cols[i], cols[j]) or 0.0) for j in keep] for i in keep]
    return [ids[i] for i in keep], mat


def correlation_lookup(
    opts: Mapping[str, Any], lead_id: str, ids: Sequence[str]
) -> dict[str, float]:
    """r^2 between ``lead_id`` and each of ``ids`` (used by the locus plot data)."""

    pool = list(dict.fromkeys([lead_id, *ids]))
    allowed, mat = correlation_matrix_from_options(opts, pool)
    if lead_id not in allowed:
        return {}
    k = allowed.index(lead_id)
    return {v: mat[k][j] ** 2 for j, v in enumerate(allowed)}


# --------------------------------------------------------------------------- association scan
def _binary(values: Sequence[float]) -> bool:
    return set(values) <= {0.0, 1.0} and len(set(values)) == 2


def _phenotype_vector(
    source: Mapping[str, Any], n: int, key: str = "phenotype"
) -> list[float | None]:
    raw = source.get(key)
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or len(raw) != n:
        raise InputFormatError(
            f"'{key}' must be a list with one value per sample ({n}).",
            code="INVALID_PHENOTYPE",
            hint=f"Add '{key}': [...] with None for missing values.",
        )
    out = []
    for v in raw:
        try:
            out.append(to_number(v))
        except ValueError as exc:
            raise InputFormatError(
                f"{key} value {v!r} is not numeric.",
                code="INVALID_PHENOTYPE",
                hint="Phenotypes must be numbers (0/1 for binary traits) or None.",
            ) from exc
    return out


def run_association(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    backend = str(opts.get("backend", "native")).lower()
    if backend != "native":
        import shutil

        program = _EXTERNAL_PROGRAMS.get(backend)
        if program is None:
            raise ConfigurationError(
                f"unknown GWAS backend {backend!r}.",
                code="INVALID_OPTION",
                hint="Use 'native' or one of: " + ", ".join(sorted(_EXTERNAL_PROGRAMS)) + ".",
            )
        if shutil.which(program) is None:
            raise CandidateUnavailable(
                program,
                f"external backend '{backend}' needs the '{program}' program on PATH, which was "
                "not found.",
            )
        raise CandidateUnavailable(
            program,
            f"'{program}' is installed but DNAKit has no command-line adapter for '{backend}' "
            "in this release; use backend='native' for the built-in association scan.",
        )
    if not isinstance(source, Mapping):
        raise InputFormatError(
            "input must be a mapping with 'genotypes' and 'phenotype'.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'genotypes': [[0,1,2],...], 'phenotype': [...], 'variant_ids': [...]}.",
        )
    g = parse_genotype_matrix(source, min_samples=5)
    n, m = g["n_samples"], g["n_variants"]
    y_all = _phenotype_vector(source, n)
    cov_raw = source.get("covariates", opts.get("covariates"))
    covariates: list[list[float]] | None = None
    if cov_raw is not None:
        if not isinstance(cov_raw, Sequence) or len(cov_raw) != n:
            raise InputFormatError(
                "covariates must be a list with one row per sample.",
                code="INVALID_COVARIATES",
                hint="Use a samples x covariates matrix.",
            )
        covariates = [
            [to_number(v) if to_number(v) is not None else math.nan for v in row] for row in cov_raw
        ]
    k = len(covariates[0]) if covariates else 0
    max_work = 100_000_000
    if n * m * (k + 2) ** 2 > max_work:
        raise ConfigurationError(
            "the scan is too large for the built-in implementation.",
            code="INPUT_TOO_LARGE",
            hint="Reduce samples/variants/covariates or use an external GWAS program.",
        )
    model = opt_choice(opts, "model", "auto", ("auto", "linear", "logistic"))
    maf_min = opt_number(opts, "maf_min", 0.0, minimum=0, maximum=0.5)
    valid_rows = [
        i
        for i in range(n)
        if y_all[i] is not None and not (covariates and any(math.isnan(v) for v in covariates[i]))
    ]
    observed = [y_all[i] for i in valid_rows]
    if len(valid_rows) < 5:
        raise InputFormatError(
            "fewer than 5 samples have a phenotype (and covariates).",
            code="TOO_FEW_SAMPLES",
            hint="Provide more phenotyped samples.",
        )
    if model == "auto":
        model = "logistic" if _binary(observed) else "linear"  # type: ignore[arg-type]
    if model == "logistic" and not _binary(observed):  # type: ignore[arg-type]
        raise InputFormatError(
            "logistic regression needs a 0/1 phenotype with both classes present.",
            code="INVALID_PHENOTYPE",
            hint="Recode cases as 1 and controls as 0, or use model='linear'.",
        )
    chrom, pos = g.get("chrom"), g.get("pos")
    alt, ref = g.get("alt"), g.get("ref")
    records = []
    for j in range(m):
        rows = [i for i in valid_rows if g["matrix"][i][j] is not None]
        dose = [g["matrix"][i][j] for i in rows]
        y = [y_all[i] for i in rows]
        rec: dict[str, Any] = {
            "variant_id": g["variants"][j],
            "chrom": None if chrom is None else str(chrom[j]),
            "pos": None if pos is None else int(pos[j]),
            "effect_allele": None if alt is None else str(alt[j]),
            "other_allele": None if ref is None else str(ref[j]),
            "beta": None,
            "se": None,
            "z": None,
            "p": None,
            "eaf": None,
            "n": len(rows),
            "status": "ok",
        }
        if len(rows) >= 1:
            rec["eaf"] = st.mean(dose) / 2.0
        if len(rows) < k + 4 or st.variance(dose) <= 1e-12:
            rec["status"] = "monomorphic_or_too_few_samples"
        elif min(rec["eaf"], 1 - rec["eaf"]) < maf_min:
            rec["status"] = "maf_below_threshold"
        else:
            try:
                if model == "linear":
                    if covariates:
                        design = [[1.0, *covariates[i], dose[t]] for t, i in enumerate(rows)]
                        fit = st.ols_multiple(design, y)
                        beta, se, df = fit["beta"][-1], fit["se"][-1], fit["df"]
                        t = beta / se if se else math.inf
                        p = st.t_two_sided_p(t, df)
                    else:
                        fit = st.simple_regression(dose, y)
                        beta, se, t, p = fit["beta"], fit["se"], fit["t"], fit["p"]
                    rec.update(beta=beta, se=se, z=t, p=p)
                else:
                    if len(set(y)) < 2:
                        raise ArithmeticError("one phenotype class only")
                    design = [
                        [1.0, *(covariates[i] if covariates else []), dose[t]]
                        for t, i in enumerate(rows)
                    ]
                    fit = st.logistic_regression(design, y)
                    beta, se = fit["beta"][-1], fit["se"][-1]
                    z = beta / se
                    rec.update(beta=beta, se=se, z=z, p=st.two_sided_p_from_z(z))
            except ArithmeticError:
                rec["status"] = "model_failed"
        records.append(rec)
    return {
        "model": model,
        "n_samples": len(valid_rows),
        "n_variants": m,
        "records": records,
        "test": "OLS t-test on the dosage coefficient (df = n - covariates - 2)"
        if model == "linear"
        else "logistic regression Wald z-test on the dosage coefficient (log OR per effect allele)",
    }, []


# --------------------------------------------------------------------------- simulation
def simulate_gwas_data(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    params: dict[str, Any] = dict(opts)
    if isinstance(source, Mapping):
        params.update(source)
    elif source is not None and not isinstance(source, (str, int, float)):
        raise InputFormatError(
            "input must be None or a mapping of simulation parameters.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'n_samples': 200, 'n_variants': 50, 'seed': 1}.",
        )
    n = opt_number(params, "n_samples", 500, minimum=10, maximum=100_000, integer=True)
    m = opt_number(params, "n_variants", 200, minimum=2, maximum=100_000, integer=True)
    n_causal = opt_number(params, "n_causal", min(5, m), minimum=0, maximum=m, integer=True)
    h2 = opt_number(params, "heritability", 0.5, minimum=0, maximum=1)
    seed = opt_number(params, "seed", 0, integer=True)
    maf_lo = opt_number(params, "maf_min", 0.05, minimum=0.001, maximum=0.5)
    maf_hi = opt_number(params, "maf_max", 0.5, minimum=maf_lo, maximum=0.5)
    n_chrom = opt_number(params, "n_chromosomes", 22, minimum=1, maximum=22, integer=True)
    if n * m > 5_000_000:
        raise ConfigurationError(
            "n_samples * n_variants exceeds 5,000,000.",
            code="INPUT_TOO_LARGE",
            hint="Simulate fewer samples or variants.",
        )
    rng = random.Random(seed)
    freqs = [rng.uniform(maf_lo, maf_hi) for _ in range(m)]
    matrix = [[(rng.random() < f) + (rng.random() < f) for f in freqs] for _ in range(n)]
    causal = sorted(rng.sample(range(m), n_causal))
    effects = {j: rng.gauss(0.0, math.sqrt(h2 / n_causal)) for j in causal} if n_causal else {}
    genetic = [
        sum(
            b * (matrix[i][j] - 2 * freqs[j]) / math.sqrt(2 * freqs[j] * (1 - freqs[j]))
            for j, b in effects.items()
        )
        for i in range(n)
    ]
    phenotype = [gv + rng.gauss(0.0, math.sqrt(1.0 - h2)) for gv in genetic]
    per_chrom = math.ceil(m / n_chrom)
    chrom = [str(j // per_chrom + 1) for j in range(m)]
    pos = [(j % per_chrom + 1) * 10_000 for j in range(m)]
    ids = [f"snp{j + 1}" for j in range(m)]
    data = {
        "genotypes": matrix,
        "sample_ids": [f"s{i + 1}" for i in range(n)],
        "variant_ids": ids,
        "chrom": chrom,
        "pos": pos,
        "alt": ["A"] * m,
        "ref": ["G"] * m,
    }
    scan, _ = run_association({**data, "phenotype": phenotype}, {"model": "linear"})
    return {
        "genotypes": data,
        "phenotype": phenotype,
        "true_effects": {ids[j]: b for j, b in effects.items()},
        "allele_frequencies": dict(zip(ids, freqs, strict=False)),
        "summary_statistics": scan["records"],
        "seed": seed,
        "model": (
            "dosage ~ Binomial(2, f_j); phenotype = sum_j b_j * standardised dosage_j + N(0, 1 - "
            "h2), b_j ~ N(0, h2 / n_causal)"
        ),
    }, []


# --------------------------------------------------------------------- joint / conditional analysis
def joint_ols(b: Sequence[float], R: list[list[float]], n: float) -> dict[str, Any]:
    """Exact OLS from marginal standardised effects ``b`` (= correlations with y) and the
    predictor correlation matrix.
    """

    k = len(b)
    df = n - k - 1
    if df < 1:
        raise ConfigurationError(
            "not enough samples for the number of predictors in the joint model.",
            code="TOO_FEW_SAMPLES",
            hint="Use a larger sample size n or fewer SNPs.",
        )
    try:
        inv = st.invert(R)
    except ArithmeticError as exc:
        raise ConfigurationError(
            "the LD matrix of the selected SNPs is singular (perfect LD).",
            code="SINGULAR_LD",
            hint="Drop one of each pair of perfectly correlated variants.",
        ) from exc
    beta = st.matvec(inv, b)
    r2 = sum(x * y for x, y in zip(b, beta, strict=False))
    sigma = (1.0 - r2) / df
    se = [math.sqrt(max(sigma * inv[i][i], 0.0)) for i in range(k)]
    t = [be / s if s > 0 else math.inf for be, s in zip(beta, se, strict=False)]
    return {
        "beta": beta,
        "se": se,
        "t": t,
        "p": [st.t_two_sided_p(x, df) for x in t],
        "r2": r2,
        "df": df,
        "partial_variance_factor": [1.0 / inv[i][i] for i in range(k)],
    }


def _marginal_b(z: float, n: float) -> float:
    return z / math.sqrt(n - 2 + z * z)


def _sample_size(rows: Sequence[Mapping[str, Any]], opts: Mapping[str, Any]) -> float:
    n = opt_number(opts, "n", None, minimum=4)
    if n is None:
        values = [r["n"] for r in rows if r["n"]]
        if not values:
            raise ConfigurationError(
                "the sample size is unknown.",
                code="MISSING_OPTION",
                hint="Add an n column to the table or pass options['n'].",
            )
        n = st.median(values)
    return float(n)


def _z_by_id(rows: Sequence[Mapping[str, Any]]) -> dict[str, Mapping[str, Any]]:
    return {r["variant_id"]: r for r in rows if r["variant_id"] and r["z"] is not None}


def conditional_gwas(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    given = require_option(
        opts, "condition_on", "list the variant ids whose effects are conditioned out"
    )
    if isinstance(given, str) or not isinstance(given, Sequence) or not given or len(given) > 20:
        raise ConfigurationError(
            "condition_on must be a list of 1-20 variant ids.",
            code="INVALID_OPTION",
            hint="Pass e.g. options['condition_on']=['rs1', 'rs2'].",
        )
    cond = [str(c) for c in given]
    table = _z_by_id(rows)
    ids, R = correlation_matrix_from_options(opts)
    index = {v: i for i, v in enumerate(ids)}
    unknown = [c for c in cond if c not in table or c not in index]
    if unknown:
        raise ConfigurationError(
            f"conditioning variants missing from the table or the LD matrix: {unknown}",
            code="UNKNOWN_VARIANT",
            hint="Conditioning SNPs need a z-score and an LD row.",
        )
    n = _sample_size(rows, opts)
    cset = [index[c] for c in cond]
    out = []
    for v in ids:
        if v in cond or v not in table:
            continue
        j = index[v]
        sel = [*cset, j]
        sub = [[R[a][b] for b in sel] for a in sel]
        bvec = [_marginal_b(table[ids[a]]["z"], n) for a in sel]
        fit = joint_ols(bvec, sub, n)
        out.append(
            {
                "variant_id": v,
                "marginal_z": table[v]["z"],
                "marginal_p": table[v]["p"],
                "conditional_beta_std": fit["beta"][-1],
                "conditional_se_std": fit["se"][-1],
                "conditional_t": fit["t"][-1],
                "conditional_p": fit["p"][-1],
                "r2_with_conditioning_set": 1.0 - fit["partial_variance_factor"][-1],
            }
        )
    return {
        "conditioned_on": cond,
        "n": n,
        "records": out,
        "method": (
            "exact joint OLS from marginal correlations and LD (COJO-style): standardised "
            "phenotype and genotypes, "
            "marginal r from z via r = z/sqrt(n-2+z^2), t-test with df = n - |C| - 2"
        ),
    }, notes


def joint_gwas_analysis(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    table = _z_by_id(rows)
    ids, R = correlation_matrix_from_options(opts)
    chosen = opts.get("snps")
    if chosen is None:
        chosen = [v for v in ids if v in table]
    elif isinstance(chosen, str) or not isinstance(chosen, Sequence):
        raise ConfigurationError(
            "snps must be a list of variant ids.",
            code="INVALID_OPTION",
            hint="Pass options['snps']=['rs1','rs2'].",
        )
    chosen = [str(c) for c in chosen]
    unknown = [c for c in chosen if c not in table or c not in ids]
    if unknown or not chosen:
        raise ConfigurationError(
            f"SNPs missing from the table or the LD matrix: {unknown or 'none selected'}",
            code="UNKNOWN_VARIANT",
            hint="Selected SNPs need a z-score and an LD row.",
        )
    if len(chosen) > 200:
        raise ConfigurationError(
            "at most 200 SNPs can be modelled jointly.",
            code="TOO_MANY_VARIANTS",
            hint="Select fewer SNPs.",
        )
    n = _sample_size(rows, opts)
    index = {v: i for i, v in enumerate(ids)}
    sel = [index[c] for c in chosen]
    sub = [[R[a][b] for b in sel] for a in sel]
    bvec = [_marginal_b(table[c]["z"], n) for c in chosen]
    fit = joint_ols(bvec, sub, n)
    recs = [
        {
            "variant_id": c,
            "marginal_z": table[c]["z"],
            "joint_beta_std": fit["beta"][i],
            "joint_se_std": fit["se"][i],
            "joint_t": fit["t"][i],
            "joint_p": fit["p"][i],
        }
        for i, c in enumerate(chosen)
    ]
    return {
        "n": n,
        "records": recs,
        "model_r2": fit["r2"],
        "df": fit["df"],
        "method": (
            "exact multiple regression of the standardised phenotype on the selected standardised "
            "genotypes, "
            "computed from marginal correlations and LD (beta = R^-1 r, se^2 = (1-R2)/(n-k-1) * "
            "diag(R^-1))"
        ),
    }, notes


# --------------------------------------------------------------------------- polygenic scores
def _score_pheno(source: object) -> tuple[list[float], list[float]]:
    if not isinstance(source, Mapping) or "score" not in source or "phenotype" not in source:
        raise InputFormatError(
            "input must be a mapping with 'score' and 'phenotype' lists.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'score': [...], 'phenotype': [...]}.",
        )
    s, p = source["score"], source["phenotype"]
    if (
        not isinstance(s, Sequence)
        or not isinstance(p, Sequence)
        or isinstance(s, str)
        or len(s) != len(p)
    ):
        raise InputFormatError(
            "'score' and 'phenotype' must be lists of equal length.",
            code="INVALID_INPUT_TYPE",
            hint="Provide one score and one phenotype per individual.",
        )
    pairs = []
    for a, b in zip(s, p, strict=False):
        try:
            x, y = to_number(a), to_number(b)
        except ValueError as exc:
            raise InputFormatError(
                "scores and phenotypes must be numeric.",
                code="INVALID_VALUE",
                hint="Use numbers or None.",
            ) from exc
        if x is not None and y is not None:
            pairs.append((x, y))
    if len(pairs) < 10:
        raise InputFormatError(
            "at least 10 individuals with both score and phenotype are needed.",
            code="TOO_FEW_SAMPLES",
            hint="Provide more individuals.",
        )
    return [x for x, _ in pairs], [y for _, y in pairs]


def auc_mann_whitney(scores: Sequence[float], labels: Sequence[float]) -> float:
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for t in range(i, j + 1):
            ranks[order[t]] = (i + j) / 2.0 + 1.0
        i = j + 1
    cases = [i for i, y in enumerate(labels) if y == 1]
    n1, n0 = len(cases), len(labels) - len(cases)
    return (sum(ranks[i] for i in cases) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def evaluate_prs(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    score, pheno = _score_pheno(source)
    trait = opt_choice(opts, "trait", "auto", ("auto", "quantitative", "binary"))
    if trait == "auto":
        trait = "binary" if _binary(pheno) else "quantitative"
    n = len(score)
    mu, sd = st.mean(score), math.sqrt(st.variance(score))
    if sd == 0:
        raise InputFormatError(
            "the polygenic score has no variance.",
            code="CONSTANT_SCORE",
            hint="Check how the score was computed.",
        )
    std = [(s - mu) / sd for s in score]
    if trait == "quantitative":
        fit = st.simple_regression(std, pheno)
        r = st.pearson(score, pheno)
        return {
            "trait": "quantitative",
            "n": n,
            "r": r,
            "r2": r * r,
            "slope_per_sd": fit["beta"],
            "slope_se": fit["se"],
            "p": fit["p"],
        }, []
    if not _binary(pheno):
        raise InputFormatError(
            "a binary trait needs 0/1 phenotypes with both classes present.",
            code="INVALID_PHENOTYPE",
            hint="Recode cases as 1 and controls as 0.",
        )
    n1 = int(sum(pheno))
    result: dict[str, Any] = {
        "trait": "binary",
        "n": n,
        "n_cases": n1,
        "n_controls": n - n1,
        "auc": auc_mann_whitney(score, pheno),
    }
    try:
        fit = st.logistic_regression([[1.0, x] for x in std], pheno)
        ll0 = n1 * math.log(n1 / n) + (n - n1) * math.log((n - n1) / n)
        beta, se = fit["beta"][1], fit["se"][1]
        z = beta / se
        result.update(
            odds_ratio_per_sd=math.exp(beta),
            or_ci95=[
                math.exp(beta - 1.959963984540054 * se),
                math.exp(beta + 1.959963984540054 * se),
            ],
            log_or=beta,
            se=se,
            p=st.two_sided_p_from_z(z),
            nagelkerke_r2=(1 - math.exp(2.0 / n * (ll0 - fit["loglik"])))
            / (1 - math.exp(2.0 * ll0 / n)),
        )
    except ArithmeticError:
        result["note"] = (
            "logistic regression did not converge (perfect separation); only the AUC is reported"
        )
    return result, []


def calibrate_prs_risk(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    score, pheno = _score_pheno(source)
    if not _binary(pheno):
        raise InputFormatError(
            "risk calibration needs a 0/1 phenotype with both classes present.",
            code="INVALID_PHENOTYPE",
            hint="Recode cases as 1 and controls as 0.",
        )
    n = len(score)
    mu, sd = st.mean(score), math.sqrt(st.variance(score))
    if sd == 0:
        raise InputFormatError(
            "the polygenic score has no variance.",
            code="CONSTANT_SCORE",
            hint="Check how the score was computed.",
        )
    bins = opt_number(opts, "n_bins", 10, minimum=2, maximum=100, integer=True)
    prevalence = opt_number(opts, "prevalence", None, minimum=0, maximum=1, exclusive=True)
    try:
        fit = st.logistic_regression([[1.0, (x - mu) / sd] for x in score], pheno)
    except ArithmeticError as exc:
        raise InputFormatError(
            "logistic calibration failed (perfect separation or singular data).",
            code="CALIBRATION_FAILED",
            hint="Check the score and phenotype.",
        ) from exc
    intercept, slope = fit["beta"]
    sample_prev = sum(pheno) / n
    if prevalence is not None:
        intercept += math.log(prevalence / (1 - prevalence)) - math.log(
            sample_prev / (1 - sample_prev)
        )

    def risk(x: float) -> float:
        eta = intercept + slope * (x - mu) / sd
        return 1.0 / (1.0 + math.exp(-eta)) if eta >= 0 else math.exp(eta) / (1.0 + math.exp(eta))

    predicted = [risk(x) for x in score]
    order = sorted(range(n), key=lambda i: predicted[i])
    table = []
    for b in range(bins):
        chunk = order[b * n // bins : (b + 1) * n // bins]
        if chunk:
            table.append(
                {
                    "bin": b + 1,
                    "n": len(chunk),
                    "mean_predicted": st.mean([predicted[i] for i in chunk]),
                    "observed": st.mean([pheno[i] for i in chunk]),
                }
            )
    new = opts.get("new_scores")
    new_risks = None
    if new is not None:
        if isinstance(new, str) or not isinstance(new, Sequence):
            raise ConfigurationError(
                "new_scores must be a list of numbers.",
                code="INVALID_OPTION",
                hint="Pass options['new_scores']=[...].",
            )
        new_risks = [risk(float(x)) for x in new]
    return {
        "intercept": intercept,
        "slope_per_sd": slope,
        "score_mean": mu,
        "score_sd": sd,
        "sample_prevalence": sample_prev,
        "prevalence_adjusted": prevalence is not None,
        "predicted_risk": predicted,
        "calibration_table": table,
        "brier_score": st.mean([(p - y) ** 2 for p, y in zip(predicted, pheno, strict=False)]),
        "mean_predicted": st.mean(predicted),
        "new_score_risk": new_risks,
        "method": (
            "logistic recalibration risk = expit(a + b * standardised score); with `prevalence`, "
            "a is shifted by "
            "logit(prevalence) - logit(sample prevalence)"
        ),
    }, []


def clump_and_threshold_prs(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    if not isinstance(source, Mapping) or "summary_statistics" not in source:
        raise InputFormatError(
            "input must be a mapping with 'summary_statistics' and 'genotypes'.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'summary_statistics': table, 'genotypes': [[...]], 'variant_ids': [...]}.",
        )
    rows, _, notes = load_sumstats(source["summary_statistics"], opts)
    _need(rows, "beta", "p", "chrom", "pos")
    rest = {k: v for k, v in source.items() if k != "summary_statistics"}
    target = parse_genotype_matrix(
        rest["genotypes"]
        if isinstance(rest.get("genotypes"), Mapping) and "genotypes" in rest["genotypes"]
        else rest,
        what="target genotypes",
    )
    ref_src = opts.get("ld_reference")
    reference = (
        parse_genotype_matrix(ref_src, what="LD reference genotypes", min_samples=3)
        if ref_src is not None
        else target
    )
    thresholds = opts.get("p_thresholds", [5e-8, 1e-5, 1e-3, 0.05, 0.5, 1.0])
    if (
        isinstance(thresholds, str)
        or not isinstance(thresholds, Sequence)
        or not thresholds
        or any(
            isinstance(t, bool) or not isinstance(t, (int, float)) or not 0 < t <= 1
            for t in thresholds
        )
    ):
        raise ConfigurationError(
            "p_thresholds must be a non-empty list of numbers in (0, 1].",
            code="INVALID_OPTION",
            hint="Use e.g. [5e-8, 1e-5, 0.05, 1].",
        )
    r2_max = opt_number(opts, "clump_r2", 0.1, minimum=0, maximum=1)
    window = opt_number(opts, "clump_window", 250_000, minimum=0, integer=True)
    t_index = {v: j for j, v in enumerate(target["variants"])}
    r_index = {v: j for j, v in enumerate(reference["variants"])}
    usable = [
        r
        for r in rows
        if r["variant_id"] in t_index
        and r["variant_id"] in r_index
        and r["beta"] is not None
        and r["p"] is not None
        and r["chrom"]
        and r["pos"] is not None
    ]
    if not usable:
        raise InputFormatError(
            "no summary-statistic variant matches the genotype variant ids.",
            code="NO_OVERLAP",
            hint="variant_id values must be shared between the two inputs.",
        )
    max_p = max(thresholds)
    cand = sorted(
        (r for r in usable if r["p"] <= max_p),
        key=lambda r: (r["p"], chrom_order(r["chrom"]), r["pos"]),
    )
    cols: dict[str, list[float | None]] = {}

    def column(v: str) -> list[float | None]:
        if v not in cols:
            cols[v] = [row[r_index[v]] for row in reference["matrix"]]
        return cols[v]

    taken: set[str] = set()
    clumps = []
    for r in cand:
        if r["variant_id"] in taken:
            continue
        taken.add(r["variant_id"])
        members = []
        for o in cand:
            if (
                o["variant_id"] in taken
                or o["chrom"] != r["chrom"]
                or abs(o["pos"] - r["pos"]) > window
            ):
                continue
            rr = column_r(column(r["variant_id"]), column(o["variant_id"]))
            if rr is not None and rr * rr >= r2_max:
                taken.add(o["variant_id"])
                members.append(o["variant_id"])
        clumps.append(
            {"index": r["variant_id"], "p": r["p"], "beta": r["beta"], "members": members}
        )
    by_id = {r["variant_id"]: r for r in usable}
    means = {}
    for c in clumps:
        vals = [
            row[t_index[c["index"]]]
            for row in target["matrix"]
            if row[t_index[c["index"]]] is not None
        ]
        means[c["index"]] = st.mean(vals) if vals else 0.0
    results = []
    for t in sorted(thresholds):
        chosen = [c for c in clumps if c["p"] <= t]
        scores = {}
        for i, sample in enumerate(target["samples"]):
            total = 0.0
            for c in chosen:
                d = target["matrix"][i][t_index[c["index"]]]
                total += by_id[c["index"]]["beta"] * (means[c["index"]] if d is None else d)
            scores[sample] = total
        results.append(
            {
                "p_threshold": t,
                "n_variants": len(chosen),
                "variants": [c["index"] for c in chosen],
                "scores": scores,
            }
        )
    return {
        "clumps": clumps,
        "thresholds": results,
        "clump_r2": r2_max,
        "clump_window": window,
        "method": (
            "PLINK-style greedy clumping by p-value (r^2 from the LD reference, same chromosome, "
            "within window bp), "
            "then score = sum(beta * effect-allele dosage) over index variants with p <= "
            "threshold; "
            "missing dosages imputed with the variant mean. Effect alleles are assumed to match "
            "the "
            "dosage coding."
        ),
    }, notes


# --------------------------------------------------------------------------- GRM, heritability
def vanraden_grm(matrix: Sequence[Sequence[float | None]]) -> dict[str, Any]:
    """VanRaden (2008) method-1 GRM with frequencies estimated from the sample; monomorphic
    markers dropped.
    """

    n, m = len(matrix), len(matrix[0])
    means = []
    for j in range(m):
        vals = [row[j] for row in matrix if row[j] is not None]
        means.append(st.mean(vals) if vals else None)
    keep = [j for j in range(m) if means[j] is not None and 0.0 < means[j] < 2.0]
    if not keep:
        raise InputFormatError(
            "no polymorphic markers remain.",
            code="NO_VARIABLE_MARKERS",
            hint="Provide variable markers.",
        )
    freqs = [means[j] / 2.0 for j in keep]
    denom = 2.0 * sum(p * (1 - p) for p in freqs)
    np = _numpy()
    if np is not None:
        z = np.array(
            [
                [(matrix[i][j] if matrix[i][j] is not None else means[j]) - means[j] for j in keep]
                for i in range(n)
            ],
            dtype=float,
        )
        g = (z @ z.T / denom).tolist()
    else:
        if n * n * len(keep) / 2 > 3e7:
            raise ConfigurationError(
                "the GRM is too large for the pure-Python path.",
                code="INPUT_TOO_LARGE",
                hint="Install numpy (pip install numpy) or reduce samples/markers.",
            )
        z = [
            [(matrix[i][j] if matrix[i][j] is not None else means[j]) - means[j] for j in keep]
            for i in range(n)
        ]
        g = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for k in range(i, n):
                value = math.fsum(a * b for a, b in zip(z[i], z[k], strict=False)) / denom
                g[i][k] = g[k][i] = value
    return {"matrix": g, "n_markers": len(keep), "denominator": denom}


def _he_slope(a_values: Sequence[float], c_values: Sequence[float]) -> float:
    ma, mc = st.mean(a_values), st.mean(c_values)
    sxx = math.fsum((a - ma) ** 2 for a in a_values)
    if sxx == 0:
        raise InputFormatError(
            "the relationship matrix has no off-diagonal variance.",
            code="DEGENERATE_GRM",
            hint="The samples look identical or too few markers were used.",
        )
    return math.fsum((a - ma) * (c - mc) for a, c in zip(a_values, c_values, strict=False)) / sxx


def _standardize_vec(values: Sequence[float]) -> list[float]:
    mu, sd = st.mean(values), math.sqrt(st.variance(values))
    if sd == 0:
        raise InputFormatError(
            "the phenotype has no variance.",
            code="CONSTANT_PHENOTYPE",
            hint="Check the phenotype vector.",
        )
    return [(v - mu) / sd for v in values]


def _individual_inputs(
    source: Mapping[str, Any], keys: Sequence[str]
) -> tuple[dict[str, Any], list[list[float]]]:
    g = parse_genotype_matrix(source, min_samples=10, min_variants=2)
    vecs = [_phenotype_vector(source, g["n_samples"], k) for k in keys]
    rows = [i for i in range(g["n_samples"]) if all(v[i] is not None for v in vecs)]
    if len(rows) < 10:
        raise InputFormatError(
            "fewer than 10 samples have all phenotypes.",
            code="TOO_FEW_SAMPLES",
            hint="Provide more samples.",
        )
    g["matrix"] = [g["matrix"][i] for i in rows]
    return g, [[v[i] for i in rows] for v in vecs]  # type: ignore[misc]


def _he_pairs(
    grm: Sequence[Sequence[float]], y1: Sequence[float], y2: Sequence[float] | None = None
) -> tuple[list[float], list[float]]:
    n = len(grm)
    a, c = [], []
    for i in range(n):
        for j in range(i + 1, n):
            a.append(grm[i][j])
            c.append(y1[i] * y1[j] if y2 is None else 0.5 * (y1[i] * y2[j] + y1[j] * y2[i]))
    return a, c


def _ldsc_table(rows: Sequence[Mapping[str, Any]], opts: Mapping[str, Any]) -> list[dict[str, Any]]:
    scores = opts.get("ld_scores")
    out = []
    for r in rows:
        ell = r.get("ld_score")
        if ell is None and isinstance(scores, Mapping):
            ell = scores.get(r["variant_id"])
        if ell is not None and r["z"] is not None:
            out.append({**r, "ld_score": float(ell)})
    if len(out) < 10:
        raise ConfigurationError(
            "LD-score regression needs at least 10 variants with a z-score and an LD score.",
            code="MISSING_LD_SCORES",
            hint="Add an 'ldscore' column or pass options['ld_scores'] = {variant_id: score}.",
        )
    return out


def _ldsc_fit(
    ell: Sequence[float], y: Sequence[float], fix_intercept: bool, target_intercept: float
) -> tuple[float, float]:
    if fix_intercept:
        slope = math.fsum(
            ld * (v - target_intercept) for ld, v in zip(ell, y, strict=False)
        ) / math.fsum(ld * ld for ld in ell)
        return slope, target_intercept
    ml, my = st.mean(ell), st.mean(y)
    sll = math.fsum((ld - ml) ** 2 for ld in ell)
    if sll == 0:
        raise InputFormatError(
            "LD scores have no variance.",
            code="DEGENERATE_LD_SCORES",
            hint="Provide per-variant LD scores.",
        )
    slope = math.fsum((ld - ml) * (v - my) for ld, v in zip(ell, y, strict=False)) / sll
    return slope, my - slope * ml


_NO_SE_WARNING = (
    "no standard errors or confidence intervals are computed (unweighted regression without a jackknife); "
    "treat the estimates as point estimates only"
)


def _h2_range_warnings(value: float, label: str) -> list[str]:
    if value < 0:
        return [
            f"{label} = {value:.4g} is negative: sampling noise (small sample or weak signal); the estimate is kept unclamped and is not a valid heritability"
        ]
    if value > 1:
        return [
            f"{label} = {value:.4g} exceeds 1: sampling noise or model misfit; the estimate is kept unclamped and is not a valid heritability"
        ]
    return []


def estimate_snp_heritability(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    method = opt_choice(opts, "method", "auto", ("auto", "haseman_elston", "ldsc"))
    if method == "auto":
        method = (
            "haseman_elston" if isinstance(source, Mapping) and "genotypes" in source else "ldsc"
        )
    if method == "haseman_elston":
        if not isinstance(source, Mapping):
            raise InputFormatError(
                "Haseman-Elston needs a mapping with 'genotypes' and 'phenotype'.",
                code="INVALID_INPUT_TYPE",
                hint="Use {'genotypes': [[...]], 'phenotype': [...]}.",
            )
        g, (y,) = _individual_inputs(source, ["phenotype"])
        grm = vanraden_grm(g["matrix"])
        a, c = _he_pairs(grm["matrix"], _standardize_vec(y))
        h2 = _he_slope(a, c)
        he_warnings = _h2_range_warnings(h2, "h2")
        return {
            "method": "haseman_elston_cross_product",
            "h2": h2,
            "n_samples": len(y),
            "n_markers": grm["n_markers"],
            "definition": (
                "slope of y_i*y_j on VanRaden GRM_ij over all pairs i<j (with intercept), y "
                "scaled to unit variance"
            ),
        }, he_warnings
    rows, _, notes = load_sumstats(source, opts)
    table = _ldsc_table(rows, opts)
    n = _sample_size(table, opts)
    m = opt_number(opts, "m", len(table), minimum=1)
    fix = opt_bool(opts, "fix_intercept", False)
    ell = [r["ld_score"] for r in table]
    chi = [r["z"] ** 2 for r in table]
    slope, intercept = _ldsc_fit(ell, chi, fix, 1.0)
    ldsc_h2 = slope * m / n
    ldsc_warnings = [_NO_SE_WARNING, *_h2_range_warnings(ldsc_h2, "h2")]
    return {
        "method": "ldsc_unweighted",
        "h2": ldsc_h2,
        "intercept": intercept,
        "slope": slope,
        "n": n,
        "m": m,
        "n_variants": len(table),
        "definition": (
            "OLS of chi2 = z^2 on LD score l_j: E[chi2] = N*h2/M*l_j + intercept; h2 = slope*M/N. "
            "Unweighted, no jackknife standard error."
        ),
    }, notes + ldsc_warnings


def _rg_warnings(rg: float | None) -> list[str]:
    if rg is not None and abs(rg) > 1:
        return [
            f"rg = {rg:.4g} is outside [-1, 1]: the ratio estimator is not bounded and small samples or noisy covariance can push it past 1; "
            "the value is reported unclamped and carries large uncertainty"
        ]
    return []


def estimate_genetic_correlation(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    if not isinstance(source, Mapping):
        raise InputFormatError(
            "input must be a mapping.",
            code="INVALID_INPUT_TYPE",
            hint=(
                "Use {'genotypes', 'phenotype_1', 'phenotype_2'} or "
                "{'trait1': table, 'trait2': table}."
            ),
        )
    if "genotypes" in source:
        g, (y1, y2) = _individual_inputs(source, ["phenotype_1", "phenotype_2"])
        grm = vanraden_grm(g["matrix"])
        z1, z2 = _standardize_vec(y1), _standardize_vec(y2)
        a, c11 = _he_pairs(grm["matrix"], z1)
        _, c22 = _he_pairs(grm["matrix"], z2)
        _, c12 = _he_pairs(grm["matrix"], z1, z2)
        h1, h2, cov = _he_slope(a, c11), _he_slope(a, c22), _he_slope(a, c12)
        rg = cov / math.sqrt(h1 * h2) if h1 > 0 and h2 > 0 else None
        he_warnings = (
            _rg_warnings(rg)
            + _h2_range_warnings(h1, "h2_trait1")
            + _h2_range_warnings(h2, "h2_trait2")
        )
        return {
            "method": "haseman_elston_cross_product",
            "rg": rg,
            "genetic_covariance_std": cov,
            "h2_trait1": h1,
            "h2_trait2": h2,
            "n_samples": len(y1),
            "definition": (
                "rg = slope(0.5*(y1_i*y2_j+y1_j*y2_i) ~ GRM_ij) / sqrt(h2_1*h2_2), phenotypes "
                "scaled to unit variance"
            ),
        }, (
            [] if rg is not None else ["a heritability estimate is not positive; rg is undefined"]
        ) + he_warnings
    if "trait1" not in source or "trait2" not in source:
        raise InputFormatError(
            "provide genotypes with phenotype_1/phenotype_2, or trait1 and trait2 summary "
            "statistics.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'trait1': table, 'trait2': table} with LD scores.",
        )
    r1, _, notes = load_sumstats(source["trait1"], opts)
    r2, _, _ = load_sumstats(source["trait2"], opts)
    t1 = _ldsc_table(r1, opts)
    second = {r["variant_id"]: r for r in r2 if r["variant_id"] and r["z"] is not None}
    paired = []
    for r in t1:
        o = second.get(r["variant_id"])
        if o is None:
            continue
        sign = 1.0
        if r["effect_allele"] and o["effect_allele"]:
            if (r["effect_allele"], r["other_allele"]) == (o["effect_allele"], o["other_allele"]):
                sign = 1.0
            elif (r["effect_allele"], r["other_allele"]) == (o["other_allele"], o["effect_allele"]):
                sign = -1.0
            else:
                continue
        paired.append((r["ld_score"], r["z"], sign * o["z"], r["n"], o["n"]))
    if len(paired) < 10:
        raise InputFormatError(
            "fewer than 10 variants are shared between the two traits.",
            code="NO_OVERLAP",
            hint="Check variant ids and alleles.",
        )
    n1 = float(opt_number(opts, "n1", None, minimum=1) or st.median([p[3] for p in paired if p[3]]))
    n2 = float(opt_number(opts, "n2", None, minimum=1) or st.median([p[4] for p in paired if p[4]]))
    m = opt_number(opts, "m", len(paired), minimum=1)
    fix = opt_bool(opts, "fix_intercept", False)
    ell = [p[0] for p in paired]
    s1, _ = _ldsc_fit(ell, [p[1] ** 2 for p in paired], fix, 1.0)
    s2, _ = _ldsc_fit(ell, [p[2] ** 2 for p in paired], fix, 1.0)
    sx, ix = _ldsc_fit(ell, [p[1] * p[2] for p in paired], fix, 0.0)
    h1, h2 = s1 * m / n1, s2 * m / n2
    rho = sx * m / math.sqrt(n1 * n2)
    warnings = [_NO_SE_WARNING]
    rg = None
    if h1 > 0 and h2 > 0:
        rg = rho / math.sqrt(h1 * h2)
        warnings.extend(_rg_warnings(rg))
    else:
        warnings.append("a heritability estimate is not positive; rg is undefined")
    warnings.extend(_h2_range_warnings(h1, "h2_trait1") + _h2_range_warnings(h2, "h2_trait2"))
    return {
        "method": "ldsc_unweighted",
        "rg": rg,
        "genetic_covariance": rho,
        "h2_trait1": h1,
        "h2_trait2": h2,
        "cross_intercept": ix,
        "n_variants": len(paired),
        "definition": (
            "rg = rho_g / sqrt(h1*h2); rho_g = slope(z1*z2 ~ LD score)*M/sqrt(N1*N2); unweighted "
            "OLS"
        ),
    }, notes + warnings
