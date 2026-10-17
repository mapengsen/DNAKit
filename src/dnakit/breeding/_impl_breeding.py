"""Native quantitative-genetics and breeding computations (breeder's equation, pedigree
relationships, GBLUP/ridge, selection indices, optimum-contribution selection and a small forward
breeding simulation).
"""

from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.gwas import _impl_stats as st
from dnakit.gwas._impl_common import (
    CandidateUnavailable,
    load_table,
    nice_ticks,
    opt_bool,
    opt_choice,
    opt_number,
    parse_genotype_matrix,
    svg_escape,
    svg_wrap,
    to_number,
    write_text_file,
)

_UNKNOWN_TEXT = {"", "0", "NA", "na", "NaN", "nan", ".", "None"}


def _unknown(value: object) -> bool:
    return value is None or str(value).strip() in _UNKNOWN_TEXT


def _mapping(source: object, what: str) -> Mapping[str, Any]:
    if not isinstance(source, Mapping):
        raise InputFormatError(
            f"input must be a mapping ({what}).", code="INVALID_INPUT_TYPE", hint=f"Pass {what}."
        )
    return source


def _num(source: Mapping[str, Any], opts: Mapping[str, Any], key: str, **kw: Any) -> Any:
    merged = dict(opts)
    merged.update({k: v for k, v in source.items() if k == key})
    return opt_number(merged, key, kw.pop("default", None), **kw)


# --------------------------------------------------------------------------- selection intensity
def selection_intensity(proportion: float) -> float:
    """Mean standardised selection differential of truncation selection of the top ``proportion``
    (infinite population).
    """

    x = st.normal_ppf(1.0 - proportion)
    return math.exp(-x * x / 2.0) / math.sqrt(2.0 * math.pi) / proportion


# --------------------------------------------------------------------------- pedigree relationship
def pedigree_relationship(
    records: Sequence[Mapping[str, Any]], limit: int = 2000
) -> dict[str, Any]:
    """Additive relationship matrix A by Henderson's tabular method; F_i = A_ii - 1."""

    if len(records) > limit:
        raise ConfigurationError(
            f"the pedigree has {len(records)} individuals; the limit is {limit}.",
            code="INPUT_TOO_LARGE",
            hint="Use a smaller pedigree or a dedicated pedigree package.",
        )
    parents: dict[str, tuple[str | None, str | None]] = {}
    for rec in records:
        ident = rec.get("id", rec.get("animal", rec.get("individual")))
        if _unknown(ident):
            raise InputFormatError(
                "every pedigree record needs an id.",
                code="INVALID_PEDIGREE",
                hint="Columns: id, sire, dam.",
            )
        ident = str(ident)
        if ident in parents:
            raise InputFormatError(
                f"duplicate individual {ident!r} in the pedigree.",
                code="INVALID_PEDIGREE",
                hint="Ids must be unique.",
            )
        sire, dam = rec.get("sire", rec.get("father")), rec.get("dam", rec.get("mother"))
        parents[ident] = (
            None if _unknown(sire) else str(sire),
            None if _unknown(dam) else str(dam),
        )
    order: list[str] = []
    state: dict[str, int] = {}

    def visit(ident: str) -> None:
        if state.get(ident) == 2:
            return
        if state.get(ident) == 1:
            raise InputFormatError(
                f"the pedigree contains a cycle through {ident!r}.",
                code="INVALID_PEDIGREE",
                hint="An individual cannot be its own ancestor.",
            )
        state[ident] = 1
        for parent in parents[ident]:
            if parent is not None:
                if parent not in parents:
                    raise InputFormatError(
                        f"parent {parent!r} of {ident!r} has no pedigree record.",
                        code="INVALID_PEDIGREE",
                        hint="List founders with unknown (empty or 0) parents.",
                    )
                visit(parent)
        state[ident] = 2
        order.append(ident)

    for ident in parents:
        visit(ident)
    index = {ident: i for i, ident in enumerate(order)}
    n = len(order)
    a = [[0.0] * n for _ in range(n)]
    for i, ident in enumerate(order):
        s, d = parents[ident]
        si, di = (index[s] if s else None), (index[d] if d else None)
        for j in range(i):
            v = 0.0
            if si is not None:
                v += a[si][j]
            if di is not None:
                v += a[di][j]
            a[i][j] = a[j][i] = 0.5 * v
        a[i][i] = 1.0 + (0.5 * a[si][di] if si is not None and di is not None else 0.0)
    return {
        "ids": order,
        "matrix": a,
        "inbreeding": {ident: a[i][i] - 1.0 for i, ident in enumerate(order)},
    }


# --------------------------------------------------------------------------- calculate
def calculate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    quantity = opt_choice(
        opts,
        "quantity",
        "breeders_equation",
        ("breeders_equation", "heritability", "pedigree_relationship", "genomic_relationship"),
    )
    if quantity == "pedigree_relationship":
        records = load_table(
            source.get("pedigree")
            if isinstance(source, Mapping) and "pedigree" in source
            else source,
            what="pedigree",
        )
        return {
            "quantity": quantity,
            **pedigree_relationship(records),
            "definition": (
                "Henderson tabular method: A_ij = (A_sj + A_dj)/2, A_ii = 1 + A_sd/2; inbreeding "
                "F = A_ii - 1"
            ),
        }, []
    if quantity == "genomic_relationship":
        try:
            from dnakit.exceptions import BackendUnavailableError
            from dnakit.genomics.breeding import genomic_relationship_matrix

            g = parse_genotype_matrix(source, min_samples=2)
            res = genomic_relationship_matrix(g["matrix"], verbose=False)
        except BackendUnavailableError as exc:
            raise CandidateUnavailable(
                "pip install numpy", f"numpy is required for the genomic relationship matrix: {exc}"
            ) from exc
        return {"quantity": quantity, "sample_ids": g["samples"], **res}, []
    src = _mapping(source, "a mapping of the quantities")
    if quantity == "heritability":
        if "parent" in src and "offspring" in src:
            x = [to_number(v) for v in src["parent"]]
            y = [to_number(v) for v in src["offspring"]]
            pairs = [(a, b) for a, b in zip(x, y, strict=False) if a is not None and b is not None]
            if len(pairs) < 3:
                raise InputFormatError(
                    "at least 3 parent-offspring pairs are needed.",
                    code="TOO_FEW_SAMPLES",
                    hint="Provide more pairs.",
                )
            kind = opt_choice(opts, "relationship", "midparent", ("midparent", "single_parent"))
            fit = st.simple_regression([p[0] for p in pairs], [p[1] for p in pairs])
            h2 = fit["beta"] if kind == "midparent" else 2.0 * fit["beta"]
            return {
                "quantity": quantity,
                "method": f"parent-offspring regression ({kind})",
                "h2": h2,
                "slope": fit["beta"],
                "slope_se": fit["se"],
                "n_pairs": len(pairs),
                "definition": (
                    "h2 = b (offspring on mid-parent) or 2b (offspring on a single parent)"
                ),
            }, []
        if "response" in src and "selection_differential" in src:
            s = _num(src, opts, "selection_differential", minimum=0, exclusive=True)
            r = _num(src, opts, "response")
            return {
                "quantity": quantity,
                "method": "realised heritability",
                "h2": r / s,
                "definition": "h2 = R / S",
            }, []
        if "genetic_variance" in src and "phenotypic_variance" in src:
            vg = _num(src, opts, "genetic_variance", minimum=0)
            vp = _num(src, opts, "phenotypic_variance", minimum=0, exclusive=True)
            if vg > vp:
                raise InputFormatError(
                    "the genetic variance cannot exceed the phenotypic variance.",
                    code="INVALID_VALUE",
                    hint="Check the variances.",
                )
            return {
                "quantity": quantity,
                "method": "variance ratio",
                "h2": vg / vp,
                "definition": "h2 = Vg / Vp (narrow-sense when Vg is additive)",
            }, []
        raise InputFormatError(
            "provide parent/offspring lists, response with selection_differential, or "
            "genetic_variance with phenotypic_variance.",
            code="MISSING_VALUES",
            hint="E.g. {'response': 2.0, 'selection_differential': 5.0}.",
        )
    h2 = _num(src, opts, "heritability", minimum=0, maximum=1)
    if h2 is None:
        raise InputFormatError(
            "'heritability' is required for the breeder's equation.",
            code="MISSING_VALUES",
            hint="E.g. {'heritability': 0.4, 'selection_differential': 5}.",
        )
    sd = _num(src, opts, "selection_differential")
    intensity = _num(src, opts, "selection_intensity", minimum=0)
    proportion = _num(src, opts, "proportion_selected", minimum=0, maximum=1, exclusive=True)
    sigma = _num(src, opts, "phenotypic_sd", minimum=0, exclusive=True)
    interval = _num(src, opts, "generation_interval", minimum=0, exclusive=True)
    if proportion is not None and intensity is None:
        intensity = selection_intensity(proportion)
    if sd is None:
        if intensity is None or sigma is None:
            raise InputFormatError(
                "give selection_differential, or (selection_intensity or proportion_selected) "
                "with phenotypic_sd.",
                code="MISSING_VALUES",
                hint="S = i * sigma_p.",
            )
        sd = intensity * sigma
    elif intensity is None and sigma is not None:
        intensity = sd / sigma
    response = h2 * sd
    result: dict[str, Any] = {
        "quantity": "breeders_equation",
        "heritability": h2,
        "selection_differential": sd,
        "selection_intensity": intensity,
        "response": response,
        "definition": (
            "R = h2 * S, S = i * sigma_p; truncation intensity i = phi(z)/p with z = Phi^-1(1-p) "
            "(infinite population, normal trait)"
        ),
    }
    if interval is not None:
        result["response_per_year"] = response / interval
    return result, []


# --------------------------------------------------------------------------- fit / predict
def _training(source: Mapping[str, Any]) -> tuple[dict[str, Any], list[float]]:
    g = (
        parse_genotype_matrix(source.get("train_genotypes", source), min_samples=3)
        if "train_genotypes" in source
        else parse_genotype_matrix(source, min_samples=3)
    )
    key = "phenotypes" if "phenotypes" in source else "phenotype"
    raw = source.get(key)
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or len(raw) != g["n_samples"]:
        raise InputFormatError(
            "'phenotypes' must be a list with one value per training sample.",
            code="INVALID_PHENOTYPE",
            hint="Add 'phenotypes': [...] (None allowed; those samples are dropped).",
        )
    values = []
    for v in raw:
        try:
            values.append(to_number(v))
        except ValueError as exc:
            raise InputFormatError(
                f"phenotype {v!r} is not numeric.",
                code="INVALID_PHENOTYPE",
                hint="Use numbers or None.",
            ) from exc
    rows = [i for i, v in enumerate(values) if v is not None]
    if len(rows) < 3:
        raise InputFormatError(
            "fewer than 3 training samples have a phenotype.",
            code="TOO_FEW_SAMPLES",
            hint="Provide more phenotyped samples.",
        )
    g["matrix"] = [g["matrix"][i] for i in rows]
    g["samples"] = [g["samples"][i] for i in rows]
    g["n_samples"] = len(rows)
    return g, [values[i] for i in rows]  # type: ignore[misc]


def _ridge_lambda(g: Mapping[str, Any], opts: Mapping[str, Any]) -> tuple[float, str]:
    lam = opt_number(opts, "ridge_lambda", None, minimum=0, exclusive=True)
    if lam is not None:
        return lam, "ridge_lambda"
    h2 = opt_number(opts, "h2", None, minimum=0, maximum=1, exclusive=True)
    if h2 is None:
        raise ConfigurationError(
            "give options['ridge_lambda'] or options['h2'] to set the shrinkage.",
            code="MISSING_OPTION",
            hint="lambda = 2*sum(p*q)*(1-h2)/h2 when only the heritability is known.",
        )
    total = 0.0
    for j in range(g["n_variants"]):
        vals = [row[j] for row in g["matrix"] if row[j] is not None]
        if vals:
            mean = sum(vals) / len(vals)
            if 0.0 < mean < 2.0:
                p = mean / 2.0
                total += 2 * p * (1 - p)
    return total * (1 - h2) / h2, "h2"


def _genomic_prediction(train: Any, y: Any, new: Any, lam: float) -> dict[str, Any]:
    try:
        from dnakit.exceptions import BackendUnavailableError
        from dnakit.genomics.breeding import genomic_prediction

        return genomic_prediction(train, y, new, ridge_lambda=lam, verbose=False)
    except BackendUnavailableError as exc:
        raise CandidateUnavailable(
            "pip install numpy", f"numpy is required for ridge/GBLUP genomic prediction: {exc}"
        ) from exc


def fit(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _mapping(source, "{'genotypes': ..., 'phenotypes': [...]}")
    g, y = _training(src)
    lam, origin = _ridge_lambda(g, opts)
    res = _genomic_prediction(g["matrix"], y, g["matrix"], lam)
    fitted = res["training_predictions"]
    r = st.pearson(fitted, y) if len(set(fitted)) > 1 else None
    return {
        "model": "rrblup_ridge",
        "ridge_lambda": lam,
        "ridge_lambda_source": origin,
        "intercept": res["intercept"],
        "raw_intercept": res["raw_intercept"],
        "marker_effects": res["marker_effects"],
        "allele_frequencies": res["allele_frequencies"],
        "n_markers_used": res["n_markers_used"],
        "n_training_samples": res["n_training_samples"],
        "training_predictions": fitted,
        "training_breeding_values": res["breeding_values"],
        "training_correlation": r,
        "sample_ids": g["samples"],
        "definition": (
            "ridge regression of the phenotype on centred dosages (SNP-BLUP), objective "
            "sum((y-mu-Zb)^2) + lambda*sum(b^2); "
            "equivalent to GBLUP with G = ZZ'/(2*sum(pq)) when lambda = 2*sum(pq)*(1-h2)/h2. "
            "Training "
            "correlation is optimistic."
        ),
    }, ["training_correlation is computed on the training data and overstates prediction accuracy"]


def predict(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _mapping(
        source,
        "{'model': <fit result>, 'genotypes': ...} or "
        "{'train_genotypes', 'phenotypes', 'genotypes'}",
    )
    if "model" in src:
        model = src["model"]
        if (
            not isinstance(model, Mapping)
            or "marker_effects" not in model
            or "raw_intercept" not in model
            or "allele_frequencies" not in model
        ):
            raise InputFormatError(
                "'model' must be a result of breeding.fit.",
                code="INVALID_MODEL",
                hint="Pass the 'result' of breeding.fit.",
            )
        new = parse_genotype_matrix(
            src.get("genotypes", src), what="genotypes to predict", min_samples=1
        )
        effects, freqs = model["marker_effects"], model["allele_frequencies"]
        if len(effects) != new["n_variants"]:
            raise InputFormatError(
                f"the model has {len(effects)} markers but the genotypes have {new['n_variants']}.",
                code="MARKER_MISMATCH",
                hint="Use the same markers in the same order as in the training data.",
            )
        preds = []
        for row in new["matrix"]:
            total = float(model["raw_intercept"])
            for x, b, f in zip(row, effects, freqs, strict=False):
                if b:
                    total += b * (2.0 * f if x is None else x)
            preds.append(total)
        return {
            "sample_ids": new["samples"],
            "predictions": preds,
            "breeding_values": [p - model["intercept"] for p in preds],
            "definition": (
                "prediction = raw_intercept + sum_j b_j * dosage_j (missing dosage -> 2*training "
                "frequency)"
            ),
        }, []
    if "train_genotypes" not in src:
        raise InputFormatError(
            "provide a fitted 'model' or 'train_genotypes' with 'phenotypes'.",
            code="MISSING_VALUES",
            hint="See the function documentation.",
        )
    g, y = _training(src)
    new = parse_genotype_matrix(src.get("genotypes"), what="genotypes to predict", min_samples=1)
    lam, origin = _ridge_lambda(g, opts)
    res = _genomic_prediction(g["matrix"], y, new["matrix"], lam)
    return {
        "sample_ids": new["samples"],
        "predictions": res["predictions"],
        "breeding_values": res["breeding_values"],
        "ridge_lambda": lam,
        "ridge_lambda_source": origin,
        "n_markers_used": res["n_markers_used"],
    }, []


# --------------------------------------------------------------------------- select
def select(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    method = opt_choice(opts, "method", "truncation", ("truncation", "index"))
    records = load_table(source, what="selection candidates")
    if not records:
        raise InputFormatError(
            "no candidates given.",
            code="EMPTY_TABLE",
            hint="Provide records with an id and trait values.",
        )
    ids = [str(r.get("id", i + 1)) for i, r in enumerate(records)]
    if len(set(ids)) != len(ids):
        raise InputFormatError(
            "candidate ids must be unique.",
            code="DUPLICATE_IDS",
            hint="Give every candidate a unique id.",
        )
    n_select = opt_number(opts, "n_select", None, minimum=1, integer=True)
    proportion = opt_number(opts, "proportion", None, minimum=0, maximum=1, exclusive=True)
    if (n_select is None) == (proportion is None):
        raise ConfigurationError(
            "give exactly one of options['n_select'] or options['proportion'].",
            code="MISSING_OPTION",
            hint="E.g. options={'n_select': 10} or {'proportion': 0.1}.",
        )
    k = n_select if n_select is not None else max(1, round(proportion * len(ids)))
    if k > len(ids):
        raise ConfigurationError(
            f"cannot select {k} of {len(ids)} candidates.",
            code="INVALID_OPTION",
            hint="Lower n_select.",
        )
    direction = 1.0 if opt_bool(opts, "higher_is_better", True) else -1.0
    extra: dict[str, Any] = {}
    if method == "truncation":
        trait = str(opts.get("trait", "ebv"))
        values = []
        for r in records:
            try:
                v = to_number(r.get(trait))
            except ValueError as exc:
                raise InputFormatError(
                    f"column '{trait}' must be numeric.", code="INVALID_VALUE", hint="Use numbers."
                ) from exc
            if v is None:
                raise InputFormatError(
                    f"candidate {r.get('id')!r} has no value for '{trait}'.",
                    code="MISSING_VALUES",
                    hint=f"Fill '{trait}' or choose options['trait'].",
                )
            values.append(v)
        scores = [direction * v for v in values]
        extra["trait"] = trait
    else:
        traits = opts.get("traits")
        weights = opts.get("economic_weights")
        gm, pm = opts.get("genetic_covariance"), opts.get("phenotypic_covariance")
        if (
            not isinstance(traits, Sequence)
            or isinstance(traits, str)
            or not traits
            or weights is None
            or gm is None
            or pm is None
        ):
            raise ConfigurationError(
                "index selection needs options traits, economic_weights, genetic_covariance and "
                "phenotypic_covariance.",
                code="MISSING_OPTION",
                hint="Smith-Hazel: b = P^-1 G a.",
            )
        t = len(traits)
        try:
            a = [float(w) for w in weights]
            gmat = [[float(v) for v in row] for row in gm]
            pmat = [[float(v) for v in row] for row in pm]
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "weights and covariance matrices must be numeric.",
                code="INVALID_OPTION",
                hint="Provide numbers.",
            ) from exc
        if len(a) != t or len(gmat) != t or len(pmat) != t or any(len(r) != t for r in gmat + pmat):
            raise ConfigurationError(
                "weights and covariance matrices must match the number of traits.",
                code="INVALID_OPTION",
                hint=f"Use {t} traits.",
            )
        try:
            ga = st.matvec(gmat, a)
            b = st.solve_linear(pmat, ga)
        except ArithmeticError as exc:
            raise ConfigurationError(
                "the phenotypic covariance matrix is singular.",
                code="SINGULAR_MATRIX",
                hint="Check the traits are not collinear.",
            ) from exc
        cols = [[to_number(r.get(c)) for r in records] for c in traits]
        if any(v is None for col in cols for v in col):
            raise InputFormatError(
                "every candidate needs a value for every index trait.",
                code="MISSING_VALUES",
                hint="Fill the trait columns.",
            )
        values = [sum(b[c] * cols[c][i] for c in range(t)) for i in range(len(ids))]
        scores = [direction * v for v in values]
        var_i = sum(bi * pi for bi, pi in zip(b, st.matvec(pmat, b), strict=False))
        var_h = sum(ai * gi for ai, gi in zip(a, ga, strict=False))
        extra = {
            "index_weights": b,
            "traits": list(traits),
            "index_variance": var_i,
            "accuracy": math.sqrt(var_i / var_h) if var_h > 0 else None,
        }
    order = sorted(range(len(ids)), key=lambda i: (-scores[i], i))
    chosen = order[:k]
    return {
        "method": method,
        "n_select": k,
        "selected": [ids[i] for i in chosen],
        "ranking": [
            {"id": ids[i], "score": values[i], "rank": rank + 1} for rank, i in enumerate(order)
        ],
        "selection_differential": st.mean([values[i] for i in chosen]) - st.mean(values),
        **extra,
        "definition": (
            "truncation: keep the k best by the trait value; index: Smith-Hazel index I = b'x "
            "with b = P^-1 G a, "
            "accuracy r_IH = sqrt(b'Pb / a'Ga)"
        ),
    }, []


# --------------------------------------------------------------------------- optimize
def optimize(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _mapping(
        source,
        "{'candidates': table(id, ebv), 'relationship': {'ids', 'matrix'} or 'pedigree': table}",
    )
    cands = load_table(src.get("candidates"), what="candidates")
    trait = str(opts.get("trait", "ebv"))
    if "relationship" in src:
        rel = src["relationship"]
        if not isinstance(rel, Mapping) or "ids" not in rel or "matrix" not in rel:
            raise InputFormatError(
                "'relationship' must be {'ids': [...], 'matrix': [[...]]}.",
                code="INVALID_RELATIONSHIP",
                hint="Provide a square matrix.",
            )
        rids, rmat = (
            [str(i) for i in rel["ids"]],
            [[float(v) for v in row] for row in rel["matrix"]],
        )
        if len(rids) != len(rmat) or any(len(r) != len(rids) for r in rmat):
            raise InputFormatError(
                "the relationship matrix must be square and match its ids.",
                code="INVALID_RELATIONSHIP",
                hint="Provide a square matrix.",
            )
    elif "pedigree" in src:
        ped = pedigree_relationship(load_table(src["pedigree"], what="pedigree"))
        rids, rmat = ped["ids"], ped["matrix"]
    else:
        raise InputFormatError(
            "provide 'relationship' ({'ids','matrix'}) or a 'pedigree'.",
            code="MISSING_VALUES",
            hint="Add the relationship source.",
        )
    pos = {i: k for k, i in enumerate(rids)}
    ids, u = [], []
    for r in cands:
        ident = str(r.get("id"))
        if ident not in pos:
            raise InputFormatError(
                f"candidate {ident!r} is missing from the relationship matrix.",
                code="UNKNOWN_CANDIDATE",
                hint="Include all candidates.",
            )
        try:
            v = to_number(r.get(trait))
        except ValueError as exc:
            raise InputFormatError(
                f"column '{trait}' must be numeric.", code="INVALID_VALUE", hint="Use numbers."
            ) from exc
        if v is None:
            raise InputFormatError(
                f"candidate {ident!r} has no '{trait}'.",
                code="MISSING_VALUES",
                hint=f"Fill '{trait}'.",
            )
        ids.append(ident)
        u.append(v)
    n = len(ids)
    if n < 2 or n > 300:
        raise ConfigurationError(
            "optimum contribution selection supports 2-300 candidates.",
            code="INVALID_INPUT_SIZE",
            hint="Pre-select candidates first.",
        )
    penalty = opt_number(opts, "penalty", 1.0, minimum=0, exclusive=True)
    a = [[rmat[pos[ids[i]]][pos[ids[j]]] for j in range(n)] for i in range(n)]

    def solve_active(active: list[int]) -> list[float]:
        sub = [[a[i][j] for j in active] for i in active]
        inv = st.invert(sub)
        ones = [1.0] * len(active)
        ainv_u = st.matvec(inv, [u[i] for i in active])
        ainv_1 = st.matvec(inv, ones)
        nu = (sum(ainv_u) - penalty) / sum(ainv_1)
        return [(x - nu * y) / penalty for x, y in zip(ainv_u, ainv_1, strict=False)] + [nu]

    active = list(range(n))
    solution = [0.0] * n
    nu = 0.0
    try:
        for _ in range(20 * n):
            sol = solve_active(active)
            nu = sol.pop()
            negatives = [(c, i) for c, i in zip(sol, active, strict=False) if c < -1e-12]
            if negatives:
                drop = min(negatives)[1]
                active.remove(drop)
                if not active:
                    raise ArithmeticError("no active candidates")
                continue
            solution = [0.0] * n
            for c, i in zip(sol, active, strict=False):
                solution[i] = max(c, 0.0)
            grads = [
                (u[i] - penalty * sum(a[i][j] * solution[j] for j in range(n)) - nu, i)
                for i in range(n)
                if i not in active
            ]
            violators = [(gv, i) for gv, i in grads if gv > 1e-9]
            if not violators:
                break
            active.append(max(violators)[1])
            active.sort()
        else:
            raise ArithmeticError("active-set iteration did not converge")
    except ArithmeticError as exc:
        raise ConfigurationError(
            f"optimum contribution selection failed: {exc}.",
            code="OPTIMIZATION_FAILED",
            hint=(
                "The relationship matrix may be singular; add a pedigree or blend the GRM with "
                "the identity."
            ),
        ) from exc
    total = sum(solution)
    solution = [c / total for c in solution]
    quad = sum(solution[i] * a[i][j] * solution[j] for i in range(n) for j in range(n))
    return {
        "contributions": {ids[i]: solution[i] for i in range(n) if solution[i] > 1e-12},
        "mean_ebv": sum(c * v for c, v in zip(solution, u, strict=False)),
        "mean_relationship": quad,
        "group_coancestry": quad / 2.0,
        "penalty": penalty,
        "objective": sum(c * v for c, v in zip(solution, u, strict=False)) - penalty / 2.0 * quad,
        "definition": (
            "maximise c'u - (penalty/2) c'Ac subject to sum(c) = 1, c >= 0 (convex QP solved by "
            "an active-set method; "
            "KKT conditions verified at the solution)"
        ),
    }, []


# --------------------------------------------------------------------------- simulate
def simulate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    params: dict[str, Any] = dict(opts)
    if isinstance(source, Mapping):
        params.update(source)
    elif source is not None and not isinstance(source, (str, int, float)):
        raise InputFormatError(
            "input must be None or a mapping of simulation parameters.",
            code="INVALID_INPUT_TYPE",
            hint="Use {'generations': 5, 'seed': 1}.",
        )
    n = opt_number(params, "n_individuals", 1000, minimum=10, maximum=50_000, integer=True)
    loci = opt_number(params, "n_loci", 100, minimum=1, maximum=2000, integer=True)
    h2 = opt_number(params, "heritability", 0.5, minimum=0, maximum=1, exclusive=True)
    prop = opt_number(params, "selection_proportion", 0.2, minimum=0, maximum=1, exclusive=True)
    gens = opt_number(params, "generations", 5, minimum=1, maximum=100, integer=True)
    seed = opt_number(params, "seed", 0, integer=True)
    if n * loci * gens > 20_000_000:
        raise ConfigurationError(
            "n_individuals * n_loci * generations exceeds 20,000,000.",
            code="INPUT_TOO_LARGE",
            hint="Reduce population size, loci or generations.",
        )
    rng = random.Random(seed)
    effects = [rng.gauss(0.0, 1.0) for _ in range(loci)]
    full = (1 << loci) - 1
    pop = [(rng.getrandbits(loci), rng.getrandbits(loci)) for _ in range(n)]

    def genetic(ind: tuple[int, int]) -> float:
        h1, h2_ = ind
        return sum(e * (((h1 >> j) & 1) + ((h2_ >> j) & 1)) for j, e in enumerate(effects))

    g_values = [genetic(i) for i in pop]
    base_var = st.variance(g_values)
    ve = base_var * (1 - h2) / h2
    sd_e = math.sqrt(ve)
    rows = []
    for gen in range(gens + 1):
        phen = [g + rng.gauss(0.0, sd_e) for g in g_values]
        vg, vp = st.variance(g_values), st.variance(phen)
        order = sorted(range(n), key=lambda i: -phen[i])
        k = max(2, round(prop * n))
        chosen = order[:k]
        sel_diff = st.mean([phen[i] for i in chosen]) - st.mean(phen)
        row = {
            "generation": gen,
            "mean_genetic_value": st.mean(g_values),
            "genetic_variance": vg,
            "phenotypic_variance": vp,
            "realised_heritability_in_generation": vg / vp,
            "selection_differential": sel_diff,
            "predicted_response": (vg / vp) * sel_diff,
        }
        rows.append(row)
        if gen == gens:
            break
        new_pop = []
        for _ in range(n):
            p1, p2 = rng.choice(chosen), rng.choice(chosen)
            while p2 == p1 and k > 1:
                p2 = rng.choice(chosen)
            kids = []
            for parent in (p1, p2):
                mask = rng.getrandbits(loci)
                a, b = pop[parent]
                kids.append((a & mask) | (b & (full ^ mask)))
            new_pop.append((kids[0], kids[1]))
        pop = new_pop
        g_values = [genetic(i) for i in pop]
    for idx in range(1, len(rows)):
        rows[idx]["observed_response"] = (
            rows[idx]["mean_genetic_value"] - rows[idx - 1]["mean_genetic_value"]
        )
    return {
        "generations": rows,
        "parameters": {
            "n_individuals": n,
            "n_loci": loci,
            "heritability": h2,
            "selection_proportion": prop,
            "environmental_variance": ve,
            "seed": seed,
        },
        "definition": (
            "unlinked additive biallelic loci (base frequency 0.5, effects ~ N(0,1)); phenotype = "
            "genetic value + N(0, Ve) with Ve fixed "
            "from the base-population variance so base h2 = heritability; truncation selection of "
            "the "
            "best phenotypes, random mating "
            "of selected parents (no selfing), free recombination. predicted_response = (Vg/Vp)*S "
            "for "
            "the parent generation."
        ),
    }, []


# --------------------------------------------------------------------------- plot
def plot(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    src = _mapping(source, "{'series': {name: [values]}} or a breeding.simulate result")
    if "series" in src and isinstance(src["series"], Mapping):
        series = {str(k): [float(v) for v in vals] for k, vals in src["series"].items()}
    elif "generations" in src and isinstance(src["generations"], Sequence):
        key = str(opts.get("metric", "mean_genetic_value"))
        try:
            series = {key: [float(r[key]) for r in src["generations"]]}
        except (KeyError, TypeError, ValueError) as exc:
            raise InputFormatError(
                f"the simulation rows have no numeric '{key}'.",
                code="INVALID_INPUT",
                hint="Set options['metric'].",
            ) from exc
    else:
        raise InputFormatError(
            "input needs a 'series' mapping or simulation 'generations'.",
            code="INVALID_INPUT",
            hint="See the function documentation.",
        )
    if not series or any(len(v) < 1 for v in series.values()):
        raise InputFormatError(
            "series must be non-empty.", code="INVALID_INPUT", hint="Provide values per generation."
        )
    width, height = 640, 400
    left, top, right, bottom = 60, 30, 20, 50
    pw, ph = width - left - right, height - top - bottom
    flat = [v for vals in series.values() for v in vals]
    lo, hi = min(flat), max(flat)
    if hi == lo:
        hi = lo + 1.0
    longest = max(len(v) for v in series.values())
    body = [f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="black"/>']
    for t in nice_ticks(lo, hi, 5):
        y = top + ph * (1 - (t - lo) / (hi - lo))
        body.append(
            f'<text x="{left - 6}" y="{y + 3:.1f}" font-size="10" text-anchor="end">{t:g}</text>'
        )
    for g in range(longest):
        x = left + pw * (g / max(1, longest - 1))
        body.append(
            f'<text x="{x:.1f}" y="{height - bottom + 14}" font-size="10" '
            f'text-anchor="middle">{g}</text>'
        )
    palette = ["#1f4e79", "#c0392b", "#27ae60", "#8e44ad", "#e69f00"]
    for idx, (name, vals) in enumerate(series.items()):
        pts = " ".join(
            f"{left + pw * (g / max(1, longest - 1)):.1f},"
            f"{top + ph * (1 - (v - lo) / (hi - lo)):.1f}"
            for g, v in enumerate(vals)
        )
        color = palette[idx % len(palette)]
        body.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2" '
            f'points="{pts}"><title>{svg_escape(name)}</title></polyline>'
        )
        body.append(
            f'<text x="{left + 8}" y="{top + 14 + 14 * idx}" font-size="11" '
            f'fill="{color}">{svg_escape(name)}</text>'
        )
    body.append(
        f'<text x="{left + pw / 2}" y="{height - 8}" font-size="12" '
        'text-anchor="middle">Generation</text>'
    )
    svg = svg_wrap(width, height, body, "breeding progress")
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, svg, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"svg": svg, "series": {k: len(v) for k, v in series.items()}, "path": written}, []
