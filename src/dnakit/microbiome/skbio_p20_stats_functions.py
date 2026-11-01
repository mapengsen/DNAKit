"""scikit-bio distance-matrix tests, ordination methods and alpha diversity (pack P20, part 3).

Count tables are lists of rows (one row per sample, one column per feature) with optional
``sample_ids`` and ``feature_ids``. Distance matrices are square lists with ``ids``. Every
function imports scikit-bio only when it runs.
"""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _skbio_p20_common as common
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_MAX_PERMUTATIONS = 99_999
_ALPHA_METRICS = (
    "ace",
    "chao1",
    "chao1_ci",
    "berger_parker_d",
    "brillouin_d",
    "dominance",
    "doubles",
    "inv_simpson",
    "esty_ci",
    "faith_pd",
    "fisher_alpha",
    "gini_index",
    "goods_coverage",
    "heip_e",
    "hill",
    "kempton_taylor_q",
    "lladser_pe",
    "lladser_ci",
)
_ORDINATION_SCALINGS = (1, 2)


def _run(label: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except ConfigurationError:
        raise
    except Exception as exc:  # scikit-bio raises ValueError and its own error classes
        raise ConfigurationError(f"scikit-bio {label} failed: {exc}") from exc


def _series_result(
    label: str, series: Any, extra: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    result = {
        "status": "completed",
        "operation": label,
        "method": str(series["method name"]),
        "test_statistic_name": str(series["test statistic name"]),
        "test_statistic": common.jsonable(series["test statistic"]),
        "p_value": common.jsonable(series["p-value"]),
        "sample_size": int(series["sample size"]),
        "number_of_groups": int(series["number of groups"]),
        "permutations": int(series["number of permutations"]),
    }
    if extra:
        result.update(extra)
    return result


def anosim_test(
    distance_matrix: Sequence[Sequence[float]],
    grouping: Sequence[str],
    *,
    ids: Sequence[str] | None = None,
    permutations: int = 999,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Analysis of similarities (ANOSIM) between groups (``skbio.stats.distance.anosim``)."""

    matrix = common.distance_matrix(distance_matrix, ids)
    groups = common.grouping_list(grouping, len(matrix.ids))
    count = common.integer(permutations, "permutations", minimum=0, maximum=_MAX_PERMUTATIONS)
    distance = common.skbio().stats.distance
    series = _run("ANOSIM", distance.anosim, matrix, groups, permutations=count, seed=seed)
    result = _series_result("anosim_test", series)
    result["provenance"] = common.provenance(method="skbio.stats.distance.anosim")
    return result


def permdisp_test(
    distance_matrix: Sequence[Sequence[float]],
    grouping: Sequence[str],
    *,
    ids: Sequence[str] | None = None,
    test: str = "median",
    permutations: int = 999,
    dimensions: int = 10,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test homogeneity of group dispersions (PERMDISP, ``skbio.stats.distance.permdisp``).

    ``test`` is ``median`` or ``centroid``. ``dimensions`` is the number of principal
    coordinates used and is reduced to the number of samples minus one when larger.
    """

    matrix = common.distance_matrix(distance_matrix, ids)
    groups = common.grouping_list(grouping, len(matrix.ids))
    count = common.integer(permutations, "permutations", minimum=0, maximum=_MAX_PERMUTATIONS)
    chosen = common.choice(test, ("median", "centroid"), "test")
    dims = min(
        common.integer(dimensions, "dimensions", minimum=1, maximum=10_000), len(matrix.ids) - 1
    )
    distance = common.skbio().stats.distance
    series = _run(
        "PERMDISP",
        distance.permdisp,
        matrix,
        groups,
        test=chosen,
        permutations=count,
        dimensions=dims,
        seed=seed,
    )
    result = _series_result("permdisp_test", series, {"test": chosen, "dimensions": dims})
    result["provenance"] = common.provenance(method="skbio.stats.distance.permdisp")
    return result


def mantel_test(
    distance_matrices: Sequence[Sequence[Sequence[float]]],
    *,
    ids: Sequence[str] | None = None,
    labels: Sequence[str] | None = None,
    method: str = "pearson",
    permutations: int = 999,
    alternative: str = "two-sided",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mantel test of the correlation between distance matrices.

    Two matrices give one Mantel test (``skbio.stats.distance.mantel``: coefficient, p-value, n);
    three or more give all pairwise tests (``pwmantel``). ``method`` is ``pearson``,
    ``spearman`` or ``kendalltau``. All matrices share the same ``ids``.
    """

    if not isinstance(distance_matrices, Sequence) or isinstance(distance_matrices, (str, bytes)):  # type: ignore[unreachable]
        raise ConfigurationError("distance_matrices must be a list of square matrices.")
    if not 2 <= len(distance_matrices) <= 50:
        raise ConfigurationError("distance_matrices must hold 2 to 50 matrices.")
    matrices = [
        common.distance_matrix(m, ids, f"distance_matrices[{i}]")
        for i, m in enumerate(distance_matrices)
    ]
    chosen = common.choice(method, ("pearson", "spearman", "kendalltau"), "method")
    count = common.integer(permutations, "permutations", minimum=0, maximum=_MAX_PERMUTATIONS)
    side = common.choice(alternative, ("two-sided", "greater", "less"), "alternative")
    distance = common.skbio().stats.distance
    result: dict[str, Any] = {"status": "completed", "operation": "mantel_test", "method": chosen}
    if len(matrices) == 2:
        coefficient, p_value, n = _run(
            "Mantel",
            distance.mantel,
            matrices[0],
            matrices[1],
            method=chosen,
            permutations=count,
            alternative=side,
            seed=seed,
        )
        result.update(
            {
                "coefficient": common.jsonable(coefficient),
                "p_value": common.jsonable(p_value),
                "n": int(n),
            }
        )
        used = "skbio.stats.distance.mantel"
    else:
        names = common.labels(labels, len(matrices), "labels", "M")
        frame = _run(
            "pairwise Mantel",
            distance.pwmantel,
            matrices,
            labels=names,
            method=chosen,
            permutations=count,
            alternative=side,
            seed=seed,
        )
        result.update(common.frame_result(frame))
        used = "skbio.stats.distance.pwmantel"
    result["provenance"] = common.provenance(method=used)
    return result


def bioenv_search(
    distance_matrix: Sequence[Sequence[float]],
    variables: Mapping[str, Sequence[float]],
    *,
    ids: Sequence[str] | None = None,
    columns: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """BIO-ENV search for the variable subset best correlated with a distance matrix.

    ``variables`` maps a numeric variable name to one value per sample
    (``skbio.stats.distance.bioenv``); subsets are scored by Euclidean distance.
    """

    matrix = common.distance_matrix(distance_matrix, ids)
    frame = common.dataframe(variables, [str(i) for i in matrix.ids], "variables")
    chosen = None if columns is None else [str(c) for c in columns]
    table = _run("BIO-ENV", common.skbio().stats.distance.bioenv, matrix, frame, columns=chosen)
    return {
        "status": "completed",
        "operation": "bioenv_search",
        **common.frame_result(table),
        "provenance": common.provenance(method="skbio.stats.distance.bioenv"),
    }


# ---------------------------------------------------------------------------- ordination
def pcoa_biplot_projection(
    distance_matrix: Sequence[Sequence[float]],
    descriptors: Mapping[str, Sequence[float]],
    *,
    ids: Sequence[str] | None = None,
    dimensions: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run PCoA on a distance matrix and project descriptor variables into it (``pcoa_biplot``).

    ``descriptors`` maps a variable name to one numeric value per sample. The result holds the
    sample coordinates and the projected variable coordinates in ``features``.
    """

    matrix = common.distance_matrix(distance_matrix, ids)
    frame = common.dataframe(descriptors, [str(i) for i in matrix.ids], "descriptors")
    dims = min(
        common.integer(dimensions, "dimensions", minimum=1, maximum=1000), len(matrix.ids) - 1
    )
    ordination = common.skbio().stats.ordination
    pcoa = _run("PCoA", ordination.pcoa, matrix, number_of_dimensions=dims)
    projected = _run("PCoA biplot", ordination.pcoa_biplot, pcoa, frame)
    return {
        "status": "completed",
        "operation": "pcoa_biplot_projection",
        "ordination": common.ordination_to_dict(projected),
        "provenance": common.provenance(method="skbio.stats.ordination.pcoa_biplot"),
    }


def correspondence_analysis(
    counts: Sequence[Sequence[float]],
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    scaling: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Correspondence analysis (CA) of a non-negative count table (``skbio.stats.ordination.ca``).

    ``scaling`` 1 preserves distances between samples, 2 between features.
    """

    frame = common.table_frame(counts, sample_ids, feature_ids)
    scale = common.integer(scaling, "scaling", minimum=1, maximum=2)
    result = _run("CA", common.skbio().stats.ordination.ca, frame, scaling=scale)
    return {
        "status": "completed",
        "operation": "correspondence_analysis",
        "scaling": scale,
        "ordination": common.ordination_to_dict(result),
        "provenance": common.provenance(method="skbio.stats.ordination.ca"),
    }


def constrained_ordination(
    response: Sequence[Sequence[float]],
    constraints: Mapping[str, Sequence[float]],
    *,
    method: str = "rda",
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    scaling: int = 1,
    scale_response: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Constrained ordination of a sample-by-feature table by environmental variables.

    ``method`` is ``rda`` (redundancy analysis, ``skbio.stats.ordination.rda``; optional
    ``scale_response`` standardises the columns) or ``cca`` (canonical correspondence analysis,
    ``cca``; the table must be non-negative). ``constraints`` maps a variable to one
    value per sample.
    """

    chosen = common.choice(method, ("rda", "cca"), "method")
    frame = common.table_frame(response, sample_ids, feature_ids, nonnegative=chosen == "cca")
    env = common.dataframe(constraints, list(frame.index), "constraints")
    scale = common.integer(scaling, "scaling", minimum=1, maximum=2)
    ordination = common.skbio().stats.ordination
    if chosen == "rda":
        result = _run(
            "RDA",
            ordination.rda,
            frame,
            env,
            scale_Y=common.flag(scale_response, "scale_response"),
            scaling=scale,
        )
    else:
        result = _run("CCA", ordination.cca, frame, env, scaling=scale)
    return {
        "status": "completed",
        "operation": "constrained_ordination",
        "method": chosen,
        "scaling": scale,
        "ordination": common.ordination_to_dict(result),
        "provenance": common.provenance(method=f"skbio.stats.ordination.{chosen}"),
    }


def mmvec_analysis(
    x_counts: Sequence[Sequence[float]],
    y_counts: Sequence[Sequence[float]],
    *,
    sample_ids: Sequence[str] | None = None,
    x_ids: Sequence[str] | None = None,
    y_ids: Sequence[str] | None = None,
    dimensions: int = 3,
    optimizer: str = "lbfgs",
    max_iter: int = 1000,
    learning_rate: float = 0.001,
    seed: int | None = None,
    predict_x: Sequence[Sequence[float]] | None = None,
    score_x: Sequence[Sequence[float]] | None = None,
    score_y: Sequence[Sequence[float]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit an MMvec microbe-metabolite co-occurrence model (``skbio.stats.ordination.mmvec``).

    ``x_counts`` (for example microbes) and ``y_counts`` (metabolites) share their samples. The
    result gives the X-by-Y ``ranks`` (log conditional probabilities), the embeddings and
    whether the optimiser converged. ``predict_x`` (new X table) additionally returns predicted
    Y proportions (``MMvecResult.predict``); ``score_x`` with ``score_y`` returns the held-out
    Q-squared (``MMvecResult.score``). Needs PyTorch.
    """

    x_table = common.table_frame(x_counts, sample_ids, x_ids)
    y_table = common.table_frame(y_counts, list(x_table.index), y_ids)
    if x_table.shape[0] != y_table.shape[0]:
        raise ConfigurationError("x_counts and y_counts must have the same number of samples.")
    dims = common.integer(dimensions, "dimensions", minimum=1, maximum=100)
    kind = common.choice(optimizer, ("lbfgs", "adam"), "optimizer")
    iterations = common.integer(max_iter, "max_iter", minimum=1, maximum=100_000)
    rate = common.number(learning_rate, "learning_rate", minimum=1e-9, maximum=10.0)
    require("torch", common.EXTRA)
    package = common.skbio()
    fitted = _run(
        "MMvec",
        package.stats.ordination.mmvec,
        x_table,
        y_table,
        dimensions=dims,
        optimizer=kind,
        max_iter=iterations,
        learning_rate=rate,
        seed=seed,
        verbose=False,
    )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "mmvec_analysis",
        "x_ids": list(x_table.columns),
        "y_ids": list(y_table.columns),
        "ranks": common.jsonable(fitted.ranks),
        "convergence": common.jsonable(fitted.convergence),
    }
    if predict_x is not None:
        new = common.table_frame(predict_x, None, list(x_table.columns))
        result["predicted_y_proportions"] = common.jsonable(
            _run("MMvec predict", fitted.predict, new)
        )
    if (score_x is None) != (score_y is None):
        raise ConfigurationError("score_x and score_y must be given together.")
    if score_x is not None:
        held_x = common.table_frame(score_x, None, list(x_table.columns))
        held_y = common.table_frame(score_y, list(held_x.index), list(y_table.columns))
        result["q_squared"] = common.jsonable(_run("MMvec score", fitted.score, held_x, held_y))
    result["provenance"] = common.provenance(method="skbio.stats.ordination.mmvec")
    return result


def plot_ordination(
    distance_matrix: Sequence[Sequence[float]],
    output_path: str,
    *,
    ids: Sequence[str] | None = None,
    grouping: Sequence[str] | None = None,
    dimensions: int = 3,
    axes: Sequence[int] = (0, 1),
    title: str = "",
    verbose: bool = True,
) -> dict[str, Any]:
    """Run PCoA on a distance matrix and draw the sample scatter plot (``OrdinationResults.plot``).

    ``grouping`` colours the samples by group; ``axes`` selects the two axes (0-based). The
    image (PNG or SVG, by extension) is written to ``output_path`` and the coordinates are returned.
    """

    matrix = common.distance_matrix(distance_matrix, ids)
    dims = min(
        common.integer(dimensions, "dimensions", minimum=2, maximum=1000), len(matrix.ids) - 1
    )
    chosen = [common.integer(a, "axes item", minimum=0, maximum=dims - 1) for a in axes]
    if len(chosen) != 2 or chosen[0] == chosen[1]:
        raise ConfigurationError("axes must name two different axes.")
    target = common.output_path(output_path)
    if target.suffix.lower() not in (".png", ".svg", ".pdf"):
        raise ConfigurationError("output_path must end with .png, .svg or .pdf.")
    pcoa = _run("PCoA", common.skbio().stats.ordination.pcoa, matrix, number_of_dimensions=dims)
    frame = None
    column = None
    if grouping is not None:
        frame = common.dataframe(
            {"group": common.grouping_list(grouping, len(matrix.ids))},
            [str(i) for i in matrix.ids],
            "grouping",
        )
        column = "group"
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = _run(
        "ordination plot", pcoa.plot, df=frame, column=column, axes=chosen, title=str(title)
    )
    figure.savefig(str(target), dpi=150, bbox_inches="tight")
    plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_ordination",
        "image": str(target),
        "axes": chosen,
        "ordination": common.ordination_to_dict(pcoa),
        "provenance": common.provenance(method="skbio.OrdinationResults.plot"),
    }


# ---------------------------------------------------------------------------- alpha diversity
def alpha_diversity_metric(
    counts: Sequence[Sequence[float]],
    metric: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    tree_newick: str | None = None,
    order: float = 2,
    bias_corrected: bool = True,
    finite: bool = False,
    rare_threshold: int = 10,
    zscore: float = 1.96,
    lower_quantile: float = 0.25,
    upper_quantile: float = 0.75,
    r: int = 10,
    alpha: float = 0.95,
    f: int = 10,
    ci_type: str = "ULCL",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute one scikit-bio alpha-diversity metric for every sample of a count table.

    ``metric`` is one of ace, chao1, chao1_ci, berger_parker_d, brillouin_d, dominance,
    doubles, inv_simpson (also ``enspie``), esty_ci, faith_pd, fisher_alpha, gini_index,
    goods_coverage, heip_e, hill, kempton_taylor_q, lladser_pe, lladser_ci. Counts are
    integers. Parameters apply to the metrics that define them and keep scikit-bio defaults:
    ``order`` (hill), ``bias_corrected`` (chao1, chao1_ci), ``finite`` (dominance,
    inv_simpson), ``rare_threshold`` (ace), ``zscore`` (chao1_ci), the two quantiles
    (kempton_taylor_q), ``r``, ``alpha``, ``f``, ``ci_type`` and ``seed`` (lladser_*;
    scikit-bio 0.7.4 does not make these two metrics reproducible from ``seed``).
    ``faith_pd`` needs ``tree_newick`` whose tips are the ``feature_ids``. Interval metrics return
    ``[lower, upper]`` for each sample.
    """

    import numpy as np

    name = common.choice("inv_simpson" if metric == "enspie" else metric, _ALPHA_METRICS, "metric")
    array = common.matrix(counts, "counts", nonnegative=True)
    if not np.allclose(array, np.round(array)):
        raise ConfigurationError("counts must be integers.")
    table = array.astype(int)
    samples = common.labels(sample_ids, table.shape[0], "sample_ids", "S")
    features = common.labels(feature_ids, table.shape[1], "feature_ids", "F")
    module = common.skbio().diversity.alpha
    function = getattr(module, name)
    kwargs: dict[str, Any] = {}
    if name == "hill":
        kwargs["order"] = common.number(order, "order", minimum=0)
    elif name in ("chao1", "chao1_ci"):
        kwargs["bias_corrected"] = common.flag(bias_corrected, "bias_corrected")
        if name == "chao1_ci":
            kwargs["zscore"] = common.number(zscore, "zscore", minimum=0)
    elif name in ("dominance", "inv_simpson"):
        kwargs["finite"] = common.flag(finite, "finite")
    elif name == "ace":
        kwargs["rare_threshold"] = common.integer(
            rare_threshold, "rare_threshold", minimum=1, maximum=10**9
        )
    elif name == "kempton_taylor_q":
        kwargs["lower_quantile"] = common.number(
            lower_quantile, "lower_quantile", minimum=0, maximum=1
        )
        kwargs["upper_quantile"] = common.number(
            upper_quantile, "upper_quantile", minimum=0, maximum=1
        )
    elif name in ("lladser_pe", "lladser_ci"):
        kwargs["r"] = common.integer(r, "r", minimum=1, maximum=10_000)
        kwargs["seed"] = seed
        if name == "lladser_ci":
            kwargs["alpha"] = common.number(alpha, "alpha", minimum=0, maximum=1)
            kwargs["f"] = common.integer(f, "f", minimum=1, maximum=10_000)
            kwargs["ci_type"] = common.choice(ci_type, ("ULCL", "ULCU", "U", "L"), "ci_type")
    elif name == "faith_pd":
        if tree_newick is None:
            raise ConfigurationError("faith_pd needs tree_newick.")
        tree = common.skbio().TreeNode.read(
            io.StringIO(common.text(tree_newick, "tree_newick", maximum=5_000_000))
        )
        kwargs["taxa"] = features
        kwargs["tree"] = tree
    values = []
    for row in table:
        values.append(_run(name, function, row.tolist(), **kwargs))
    return {
        "status": "completed",
        "operation": "alpha_diversity_metric",
        "metric": name,
        "sample_ids": samples,
        "values": common.jsonable(values),
        "provenance": common.provenance(method=f"skbio.diversity.alpha.{name}"),
    }


__all__ = [
    "alpha_diversity_metric",
    "anosim_test",
    "bioenv_search",
    "constrained_ordination",
    "correspondence_analysis",
    "mantel_test",
    "mmvec_analysis",
    "pcoa_biplot_projection",
    "permdisp_test",
    "plot_ordination",
]

_register_public_verbose(globals(), tuple(__all__))
