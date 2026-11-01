"""scikit-bio compositional data analysis and differential abundance (pack P20, part 4).

Compositions are lists of rows (one row per sample, one column per part). Count tables follow the
same layout with optional ``sample_ids`` and ``feature_ids``; metadata is a mapping of column name
to one value per sample. Every function imports scikit-bio only when it runs.
"""

from __future__ import annotations

import io
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _skbio_p20_common as common
from dnakit.exceptions import ConfigurationError

_TRANSFORMS = ("clr", "rclr", "ilr", "alr")
_POSTHOC = ("global", "pairwise", "dunnett", "trend")
_P_ADJUST = (
    "holm",
    "bonferroni",
    "sidak",
    "holm-sidak",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
    "none",
)
_MAX_ITER = 10_000


def _composition(value: object, name: str, *, positive: bool = False) -> Any:
    array = common.matrix(value, name, nonnegative=True)
    if positive and (array <= 0).any():
        raise ConfigurationError(f"{name} must be strictly positive (use zero replacement first).")
    return array


def _vector(value: object, name: str) -> Any:
    import numpy as np

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of numbers.")
    try:
        array = np.array(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a list of numbers.") from exc
    if array.ndim != 1 or not np.isfinite(array).all() or (array < 0).any():
        raise ConfigurationError(f"{name} must be a one-dimensional list of non-negative numbers.")
    return array


def _run(label: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return function(*args, **kwargs)
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio {label} failed: {exc}") from exc


def _result(operation: str, source: str, **values: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"status": "completed", "operation": operation}
    out.update({k: common.jsonable(v) for k, v in values.items()})
    out["provenance"] = common.provenance(method=source)
    return out


def _p_adjust(value: str) -> str:
    return common.choice(value, _P_ADJUST, "p_adjust")


def closure_compositions(
    matrix: Sequence[Sequence[float]], *, verbose: bool = True
) -> dict[str, Any]:
    """Scale every row so that it sums to one (``skbio.stats.composition.closure``)."""

    array = _composition(matrix, "matrix")
    if (array.sum(axis=1) == 0).any():
        raise ConfigurationError("every row of matrix needs a positive sum.")
    closure = common.skbio().stats.composition.closure
    return _result(
        "closure_compositions",
        "skbio.stats.composition.closure",
        compositions=_run("closure", closure, array),
    )


def replace_zeros_multiplicative(
    matrix: Sequence[Sequence[float]], *, delta: float | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Replace zeros by a small value and shrink the other parts (multiplicative replacement).

    ``delta`` is the replacement value (default: scikit-bio's ``1/N**2`` rule);
    rows are closed first.
    """

    array = _composition(matrix, "matrix")
    value = None if delta is None else common.number(delta, "delta", minimum=1e-300, maximum=1.0)
    function = common.skbio().stats.composition.multi_replace
    return _result(
        "replace_zeros_multiplicative",
        "skbio.stats.composition.multi_replace",
        compositions=_run("multi_replace", function, array, delta=value),
    )


def perturb_compositions(
    x: Sequence[float], y: Sequence[float], *, inverse: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Aitchison perturbation of composition ``x`` by ``y``.

    ``inverse=True`` perturbs by the inverse of ``y``.

    Wraps ``skbio.stats.composition.perturb`` and ``perturb_inv``; the result is closed to sum 1.
    """

    left, right = _vector(x, "x"), _vector(y, "y")
    if len(left) != len(right):
        raise ConfigurationError("x and y must have the same length.")
    module = common.skbio().stats.composition
    name = "perturb_inv" if common.flag(inverse, "inverse") else "perturb"
    return _result(
        "perturb_compositions",
        f"skbio.stats.composition.{name}",
        inverse=inverse,
        composition=_run(name, getattr(module, name), left, right),
    )


def power_composition(
    x: Sequence[float], exponent: float, *, verbose: bool = True
) -> dict[str, Any]:
    """Raise a composition to a power (Aitchison scalar multiplication) and close it."""

    vector = _vector(x, "x")
    power = common.number(exponent, "exponent")
    module = common.skbio().stats.composition
    return _result(
        "power_composition",
        "skbio.stats.composition.power",
        composition=_run("power", module.power, vector, power),
    )


def aitchison_inner_product(
    x: Sequence[float], y: Sequence[float], *, verbose: bool = True
) -> dict[str, Any]:
    """Aitchison inner product of two compositions (``skbio.stats.composition.inner``)."""

    left, right = _vector(x, "x"), _vector(y, "y")
    if len(left) != len(right) or (left <= 0).any() or (right <= 0).any():
        raise ConfigurationError("x and y must be positive and of equal length.")
    module = common.skbio().stats.composition
    return _result(
        "aitchison_inner_product",
        "skbio.stats.composition.inner",
        inner_product=_run("inner", module.inner, left, right),
    )


def logratio_transform(
    matrix: Sequence[Sequence[float]],
    method: str,
    *,
    inverse: bool = False,
    basis: Sequence[Sequence[float]] | None = None,
    reference_index: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Log-ratio transform (or its inverse) of compositions.

    ``method`` is ``clr`` (centred), ``rclr`` (robust centred; zeros become missing values),
    ``ilr`` (isometric; optional orthonormal ``basis``, default Gram-Schmidt basis) or ``alr``
    (additive; ``reference_index`` is the denominator part). With ``inverse=True`` the matrix is in
    transformed space and the result is the composition (closed to sum 1). Data must be strictly
    positive except for ``rclr``.
    """

    chosen = common.choice(method, _TRANSFORMS, "method")
    flag = common.flag(inverse, "inverse")
    array = common.matrix(matrix, "matrix", nonnegative=not flag and chosen != "ilr")
    module = common.skbio().stats.composition
    name = chosen + ("_inv" if flag else "")
    kwargs: dict[str, Any] = {}
    if chosen == "ilr" and basis is not None:
        kwargs["basis"] = common.matrix(basis, "basis")
    if chosen == "alr":
        kwargs["ref_idx"] = common.integer(
            reference_index, "reference_index", minimum=0, maximum=array.shape[1] - 1
        )
    if not flag and chosen in ("clr", "ilr", "alr") and (array <= 0).any():
        raise ConfigurationError(
            "matrix must be strictly positive; use replace_zeros_multiplicative first."
        )
    values = _run(name, getattr(module, name), array, **kwargs)
    return _result(
        "logratio_transform",
        f"skbio.stats.composition.{name}",
        method=chosen,
        inverse=flag,
        values=values,
    )


def centralize_compositions(
    matrix: Sequence[Sequence[float]], *, verbose: bool = True
) -> dict[str, Any]:
    """Centre compositions on their geometric mean (``skbio.stats.composition.centralize``)."""

    array = _composition(matrix, "matrix", positive=True)
    module = common.skbio().stats.composition
    return _result(
        "centralize_compositions",
        "skbio.stats.composition.centralize",
        compositions=_run("centralize", module.centralize, array),
    )


def variance_log_ratio(
    counts: Sequence[Sequence[float]],
    *,
    feature_ids: Sequence[str] | None = None,
    feature_x: str | None = None,
    feature_y: str | None = None,
    ddof: int = 1,
    robust: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variance of the log-ratio between features, a proportionality measure.

    With ``feature_x`` and ``feature_y`` the single value ``vlr`` is returned; otherwise the
    matrix of all feature pairs (``pairwise_vlr``). ``robust`` uses the median absolute deviation.
    """

    table = common.table_frame(counts, None, feature_ids)
    degrees = common.integer(ddof, "ddof", minimum=0, maximum=1)
    rob = common.flag(robust, "robust")
    module = common.skbio().stats.composition
    if (feature_x is None) != (feature_y is None):
        raise ConfigurationError("give both feature_x and feature_y, or neither.")
    if feature_x is not None:
        for name in (feature_x, feature_y):
            if name not in table.columns:
                raise ConfigurationError(f"feature {name!r} is not in feature_ids.")
        value = _run(
            "vlr",
            module.vlr,
            table[feature_x].to_numpy(),
            table[feature_y].to_numpy(),
            ddof=degrees,
            robust=rob,
        )
        return _result("variance_log_ratio", "skbio.stats.composition.vlr", vlr=value)
    matrix = _run(
        "pairwise_vlr",
        module.pairwise_vlr,
        table.to_numpy(),
        ids=list(table.columns),
        ddof=degrees,
        robust=rob,
    )
    return _result(
        "variance_log_ratio",
        "skbio.stats.composition.pairwise_vlr",
        ids=list(matrix.ids),
        matrix=matrix.data,
    )


def ilr_basis(
    *,
    sbp: Sequence[Sequence[int]] | None = None,
    newick: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an orthonormal ILR basis from a sequential binary partition or a tree.

    ``sbp`` has one row per balance with +1/-1/0 entries over the parts
    (``sbp_basis``); ``newick`` is a tree with named internal nodes whose tips are the parts
    (``tree_basis``; returns the node names of the balances).
    """

    if (sbp is None) == (newick is None):
        raise ConfigurationError("give exactly one of sbp and newick.")
    module = common.skbio().stats.composition
    if sbp is not None:
        import numpy as np

        array = common.matrix(sbp, "sbp")
        if not np.allclose(array, np.round(array)):
            raise ConfigurationError("sbp may only contain -1, 0 and 1.")
        basis = _run("sbp_basis", module.sbp_basis, np.round(array).astype(int))
        return _result("ilr_basis", "skbio.stats.composition.sbp_basis", basis=basis)
    tree = _run(
        "tree parse",
        common.skbio().TreeNode.read,
        io.StringIO(common.text(newick, "newick", maximum=5_000_000)),
    )
    sparse, names = _run("tree_basis", module.tree_basis, tree)
    return _result(
        "ilr_basis",
        "skbio.stats.composition.tree_basis",
        basis=sparse.toarray(),
        balance_nodes=[str(n) for n in names],
    )


def _metadata_frame(columns: object, sample_ids: Sequence[str], name: str = "metadata") -> Any:
    return common.dataframe(columns, sample_ids, name)


def ancom_test(
    counts: Sequence[Sequence[float]],
    grouping: Sequence[str],
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    alpha: float = 0.05,
    tau: float = 0.02,
    theta: float = 0.1,
    p_adjust: str = "holm",
    significance_test: str = "f_oneway",
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential abundance by ANCOM (``skbio.stats.composition.ancom``).

    ``counts`` must be strictly positive (add a pseudo-count first). ``significance_test`` is a
    SciPy test name such as ``f_oneway``, ``kruskal`` or ``ttest_ind``. Returns W and the
    significance per feature.
    """

    table = common.table_frame(counts, sample_ids, feature_ids)
    if (table.to_numpy() <= 0).any():
        raise ConfigurationError("ANCOM needs strictly positive counts; add a pseudo-count.")
    groups = common.grouping_list(grouping, len(table))
    import pandas as pd

    series = pd.Series(groups, index=table.index)
    module = common.skbio().stats.composition
    test = common.choice(
        significance_test,
        ("f_oneway", "kruskal", "mannwhitneyu", "ttest_ind", "ttest_rel", "wilcoxon"),
        "significance_test",
    )
    frame, _percentiles = _run(
        "ANCOM",
        module.ancom,
        table,
        series,
        alpha=common.number(alpha, "alpha", minimum=0, maximum=1),
        tau=common.number(tau, "tau", minimum=0, maximum=1),
        theta=common.number(theta, "theta", minimum=0, maximum=1),
        p_adjust=_p_adjust(p_adjust),
        sig_test=test,
    )
    return _result(
        "ancom_test",
        "skbio.stats.composition.ancom",
        features=list(frame.index),
        W=frame["W"].tolist(),
        significant=frame["Signif"].tolist(),
    )


def _bc_inputs(
    counts: object,
    sample_ids: object,
    feature_ids: object,
    metadata: object,
    formula: object,
) -> tuple[Any, Any, str]:
    table = common.table_frame(counts, sample_ids, feature_ids)
    meta = _metadata_frame(metadata, list(table.index))
    return table, meta, common.text(formula, "formula", maximum=1000)


def differential_abundance_ancombc(
    counts: Sequence[Sequence[float]],
    metadata: Mapping[str, Sequence[Any]],
    formula: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    max_iter: int = 100,
    tol: float = 1e-5,
    alpha: float = 0.05,
    p_adjust: str = "holm",
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential abundance by ANCOM-BC (bias-corrected, ``skbio.stats.composition.ancombc``).

    ``formula`` is a Patsy model over the ``metadata`` columns, for example ``"group"`` or
    ``"C(group) + age"``. Returns log fold change, standard error, W, p- and q-values per
    feature and covariate.
    """

    table, meta, model = _bc_inputs(counts, sample_ids, feature_ids, metadata, formula)
    module = common.skbio().stats.composition
    fitted = _run(
        "ANCOM-BC",
        module.ancombc,
        table,
        meta,
        model,
        max_iter=common.integer(max_iter, "max_iter", minimum=1, maximum=_MAX_ITER),
        tol=common.number(tol, "tol", minimum=1e-12, maximum=1.0),
        alpha=common.number(alpha, "alpha", minimum=0, maximum=1),
        p_adjust=_p_adjust(p_adjust),
    )
    return {
        "status": "completed",
        "operation": "differential_abundance_ancombc",
        "formula": model,
        **common.frame_result(fitted.result),
        "provenance": common.provenance(method="skbio.stats.composition.ancombc"),
    }


def differential_abundance_ancombc2(
    counts: Sequence[Sequence[float]],
    metadata: Mapping[str, Sequence[Any]],
    formula: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    pseudocount: float = 0.0,
    var_quantile: float = 0.05,
    max_iter: int = 100,
    tol: float = 1e-5,
    alpha: float = 0.05,
    p_adjust: str = "holm",
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential abundance by ANCOM-BC2 (``skbio.stats.composition.ancombc2``).

    Same inputs as :func:`differential_abundance_ancombc`; ``pseudocount`` is added to the counts
    and ``var_quantile`` regularises small variances.
    """

    table, meta, model = _bc_inputs(counts, sample_ids, feature_ids, metadata, formula)
    module = common.skbio().stats.composition
    fitted = _run(
        "ANCOM-BC2",
        module.ancombc2,
        table,
        meta,
        model,
        pseudocount=common.number(pseudocount, "pseudocount", minimum=0),
        var_quantile=common.number(var_quantile, "var_quantile", minimum=0, maximum=1),
        max_iter=common.integer(max_iter, "max_iter", minimum=1, maximum=_MAX_ITER),
        tol=common.number(tol, "tol", minimum=1e-12, maximum=1.0),
        alpha=common.number(alpha, "alpha", minimum=0, maximum=1),
        p_adjust=_p_adjust(p_adjust),
    )
    return {
        "status": "completed",
        "operation": "differential_abundance_ancombc2",
        "formula": model,
        **common.frame_result(fitted.result),
        "provenance": common.provenance(method="skbio.stats.composition.ancombc2"),
    }


def ancombc_posthoc_test(
    counts: Sequence[Sequence[float]],
    metadata: Mapping[str, Sequence[Any]],
    formula: str,
    grouping: str,
    test: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    pseudocount: float = 0.0,
    alpha: float = 0.05,
    p_adjust: str = "holm",
    bootstraps: int = 100,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit ANCOM-BC2 with a grouping term and run a post-hoc group test.

    ``test`` is ``global`` (any difference between groups), ``pairwise`` (directional pairwise
    differences), ``dunnett`` (each group against the reference, ``bootstraps`` and ``seed``) or
    ``trend`` (monotone trend over ordered groups). ``grouping`` must be a term of ``formula``.
    """

    chosen = common.choice(test, _POSTHOC, "test")
    table, meta, model = _bc_inputs(counts, sample_ids, feature_ids, metadata, formula)
    group = common.text(grouping, "grouping", maximum=256)
    module = common.skbio().stats.composition
    fitted = _run(
        "ANCOM-BC2",
        module.ancombc2,
        table,
        meta,
        model,
        grouping=group,
        pseudocount=common.number(pseudocount, "pseudocount", minimum=0),
        alpha=common.number(alpha, "alpha", minimum=0, maximum=1),
        p_adjust=_p_adjust(p_adjust),
    )
    kwargs: dict[str, Any] = {}
    if chosen in ("dunnett", "trend"):
        kwargs = {
            "bootstraps": common.integer(bootstraps, "bootstraps", minimum=1, maximum=100_000),
            "seed": seed,
        }
    frame = _run(
        f"{chosen} test",
        getattr(fitted, f"{chosen}_test" if chosen != "global" else "global_test"),
        **kwargs,
    )
    return {
        "status": "completed",
        "operation": "ancombc_posthoc_test",
        "test": chosen,
        **common.frame_result(frame),
        "provenance": common.provenance(
            method=f"skbio.stats.composition.ANCOMBCResult.{chosen}_test"
        ),
    }


def detect_structural_zeros(
    counts: Sequence[Sequence[float]],
    grouping: Sequence[str],
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    negative_lower_bound: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag features that are structurally absent from a whole group (``struc_zero``).

    Returns, for every feature, one True/False per group.
    """

    table = common.table_frame(counts, sample_ids, feature_ids)
    groups = common.grouping_list(grouping, len(table))
    meta = _metadata_frame({"group": groups}, list(table.index))
    module = common.skbio().stats.composition
    frame = _run(
        "struc_zero",
        module.struc_zero,
        table,
        meta,
        "group",
        neg_lb=common.flag(negative_lower_bound, "negative_lower_bound"),
    )
    return _result(
        "detect_structural_zeros",
        "skbio.stats.composition.struc_zero",
        groups=[str(c) for c in frame.columns],
        features=[str(i) for i in frame.index],
        structural_zero=frame.to_numpy().tolist(),
    )


def dirmult_ttest_test(
    counts: Sequence[Sequence[float]],
    grouping: Sequence[str],
    *,
    treatment: str | None = None,
    reference: str | None = None,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    pseudocount: float = 0.5,
    draws: int = 128,
    p_adjust: str = "holm",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential abundance with a Dirichlet-multinomial model and a t-test between two groups.

    ``treatment`` and ``reference`` name the two groups (``skbio.stats.composition.dirmult_ttest``);
    ``draws`` is the number of posterior draws.
    """

    table = common.table_frame(counts, sample_ids, feature_ids)
    import pandas as pd

    series = pd.Series(common.grouping_list(grouping, len(table)), index=table.index)
    module = common.skbio().stats.composition
    frame = _run(
        "dirmult_ttest",
        module.dirmult_ttest,
        table,
        series,
        treatment=treatment,
        reference=reference,
        pseudocount=common.number(pseudocount, "pseudocount", minimum=0),
        draws=common.integer(draws, "draws", minimum=1, maximum=10_000),
        p_adjust=_p_adjust(p_adjust),
        seed=seed,
    )
    return {
        "status": "completed",
        "operation": "dirmult_ttest_test",
        **common.frame_result(frame),
        "provenance": common.provenance(method="skbio.stats.composition.dirmult_ttest"),
    }


def dirmult_mixed_model(
    counts: Sequence[Sequence[float]],
    metadata: Mapping[str, Sequence[Any]],
    formula: str,
    grouping: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    pseudocount: float = 0.5,
    draws: int = 128,
    p_adjust: str = "holm",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential abundance with a Dirichlet-multinomial linear mixed model.

    ``formula`` holds the fixed effects, ``grouping`` is the metadata column that defines the
    random-effect groups (for example the subject), as in ``skbio.stats.composition.dirmult_lme``.
    """

    table, meta, model = _bc_inputs(counts, sample_ids, feature_ids, metadata, formula)
    group = common.text(grouping, "grouping", maximum=256)
    module = common.skbio().stats.composition
    frame = _run(
        "dirmult_lme",
        module.dirmult_lme,
        table,
        meta,
        model,
        group,
        pseudocount=common.number(pseudocount, "pseudocount", minimum=0),
        draws=common.integer(draws, "draws", minimum=1, maximum=10_000),
        p_adjust=_p_adjust(p_adjust),
        seed=seed,
    )
    return {
        "status": "completed",
        "operation": "dirmult_mixed_model",
        **common.frame_result(frame),
        "provenance": common.provenance(method="skbio.stats.composition.dirmult_lme"),
    }


__all__ = [
    "aitchison_inner_product",
    "ancom_test",
    "ancombc_posthoc_test",
    "centralize_compositions",
    "closure_compositions",
    "detect_structural_zeros",
    "differential_abundance_ancombc",
    "differential_abundance_ancombc2",
    "dirmult_mixed_model",
    "dirmult_ttest_test",
    "ilr_basis",
    "logratio_transform",
    "perturb_compositions",
    "power_composition",
    "replace_zeros_multiplicative",
    "variance_log_ratio",
]

_register_public_verbose(globals(), tuple(__all__))
