"""scikit-bio diversity, resampling, table augmentation and metadata functions (pack P21).

Every function takes plain JSON-friendly inputs (lists, mappings, Newick strings) and returns a
JSON-serialisable dict. scikit-bio is imported only when a function runs. Random operations take
an integer ``seed`` so that results are reproducible.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-skibio"
_PACKAGE = "scikit-bio"
_MAX_SAMPLES = 2_000
_MAX_FEATURES = 20_000
_MAX_TEXT = 1_000_000
_MAX_ITEMS = 1_000_000

ALPHA_METRICS = (
    "margalef",
    "mcintosh_d",
    "mcintosh_e",
    "menhinick",
    "michaelis_menten_fit",
    "sobs",
    "osd",
    "phydiv",
    "pielou_e",
    "renyi",
    "robbins",
    "shannon",
    "simpson",
    "simpson_e",
    "singles",
    "strong",
    "tsallis",
)
_GRADIENT_ALGORITHMS = ("average", "trajectory", "first_difference", "window_difference")


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _int(value: object, name: str, *, minimum: int, maximum: int = 10**9) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _seed(value: object) -> int | None:
    if value is None:
        return None
    return _int(value, "seed", minimum=0, maximum=2**32 - 1)


def _number(value: object, name: str, *, minimum: float | None = None, exclusive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if minimum is not None and (value <= minimum if exclusive else value < minimum):
        raise ConfigurationError(f"{name} must be {'greater than' if exclusive else 'at least'} {minimum}.")
    return float(value)


def _text(value: object, name: str, *, maximum: int = _MAX_TEXT) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {maximum} characters.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _vector(values: object, name: str, *, integer: bool = False, minimum: int = 1) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of numbers.")
    if not minimum <= len(values) <= _MAX_FEATURES:
        raise ConfigurationError(f"{name} must hold {minimum} to {_MAX_FEATURES} values.")
    out: list[float] = []
    for item in values:
        value = _number(item, f"{name} item", minimum=0.0)
        if integer and value != int(value):
            raise ConfigurationError(f"{name} must hold whole-number counts.")
        out.append(int(value) if integer else value)
    return out  # type: ignore[return-value]


def _matrix(rows: object, name: str, *, integer: bool = False, min_rows: int = 1) -> list[list[float]]:
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise ConfigurationError(f"{name} must be a list of rows (samples x features).")
    if not min_rows <= len(rows) <= _MAX_SAMPLES:
        raise ConfigurationError(f"{name} must hold {min_rows} to {_MAX_SAMPLES} rows.")
    matrix = [_vector(row, f"{name} row", integer=integer) for row in rows]
    if len({len(row) for row in matrix}) != 1:
        raise ConfigurationError(f"{name} rows must all have the same length.")
    return matrix


def _ids(values: object, count: int, name: str = "ids") -> list[str]:
    if values is None:
        return [f"S{i + 1}" for i in range(count)]
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or len(values) != count:
        raise ConfigurationError(f"{name} must be a list of {count} strings.")
    out = [_text(v, f"{name} item", maximum=256) for v in values]
    if len(set(out)) != count:
        raise ConfigurationError(f"{name} must be unique.")
    return out


def _tree(newick: object) -> Any:
    text = _text(newick, "tree_newick")
    tree_module = require("skbio.tree", _EXTRA)
    import io

    try:
        return tree_module.TreeNode.read(io.StringIO(text if text.rstrip().endswith(";") else text + ";"))
    except Exception as exc:  # noqa: BLE001 - parser errors vary
        raise ConfigurationError(f"tree_newick could not be parsed: {exc}") from exc


def _plain(value: Any) -> Any:
    """Convert numpy scalars and arrays to JSON-friendly Python values."""

    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    return value


def _labels(labels: object, rows: int) -> Any:
    if labels is None:
        return None
    import numpy as np

    if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence) or len(labels) != rows:
        raise ConfigurationError(f"labels must be a list with one class index per sample ({rows}).")
    out = []
    for item in labels:
        out.append(_int(item, "labels item", minimum=0, maximum=10_000))
    return np.asarray(out)


def _augmented(operation: str, matrix: Any, labels: Any, **fields: object) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "n_synthetic_rows": len(matrix),
        "matrix": _plain(matrix),
        "labels": None if labels is None else _plain(labels),
        **fields,
        "provenance": _provenance(),
    }


# ---------------------------------------------------------------------------- alpha diversity
def calculate_alpha_diversity(
    counts: Sequence[float],
    *,
    metric: str = "shannon",
    base: float | None = None,
    exp: bool = False,
    order: float = 2.0,
    finite: bool = False,
    num_repeats: int = 1,
    taxa: Sequence[str] | None = None,
    tree_newick: str | None = None,
    rooted: bool | None = None,
    weight: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute one scikit-bio alpha-diversity metric for one sample count vector.

    ``metric`` selects the function of ``skbio.diversity.alpha``: margalef, mcintosh_d,
    mcintosh_e, menhinick, michaelis_menten_fit (uses ``num_repeats``; scikit-bio subsamples randomly and offers no seed, so the value varies between calls), sobs (observed
    features), osd (returns observed, singletons, doubletons), phydiv (needs ``taxa`` and
    ``tree_newick``; uses ``rooted`` and ``weight``), pielou_e (``base``), renyi (``order``,
    ``base``), robbins, shannon (``base``, ``exp``), simpson (``finite``), simpson_e, singles,
    strong, tsallis (``order``). Options that a metric does not use are ignored. Counts must be
    whole numbers.
    """

    chosen = _choice(metric, ALPHA_METRICS, "metric")
    vector = _vector(counts, "counts", integer=True)
    log_base = None if base is None else _number(base, "base", minimum=1.0, exclusive=True)
    exp_flag = _flag(exp, "exp")
    order_value = _number(order, "order", minimum=0.0)
    finite_flag = _flag(finite, "finite")
    repeats = _int(num_repeats, "num_repeats", minimum=1, maximum=1000)
    alpha = require("skbio.diversity.alpha", _EXTRA)
    func = getattr(alpha, chosen)
    kwargs: dict[str, Any] = {}
    if chosen in {"pielou_e", "renyi"}:
        kwargs["base"] = log_base
    if chosen == "shannon":
        kwargs.update(base=log_base, exp=exp_flag)
    if chosen in {"renyi", "tsallis"}:
        kwargs["order"] = order_value
    if chosen == "simpson":
        kwargs["finite"] = finite_flag
    if chosen == "michaelis_menten_fit":
        kwargs["num_repeats"] = repeats
    extra: dict[str, Any] = {}
    if chosen == "phydiv":
        if taxa is None or tree_newick is None:
            raise ConfigurationError("phydiv needs taxa and tree_newick.")
        extra["taxa"] = _ids(taxa, len(vector), "taxa")
        extra["tree"] = _tree(tree_newick)
        kwargs["rooted"] = rooted
        kwargs["weight"] = _flag(weight, "weight")
    try:
        if chosen == "phydiv":
            value = func(vector, extra["taxa"], extra["tree"], **kwargs)
        else:
            value = func(vector, **kwargs)
    except (ValueError, TypeError, RuntimeError) as exc:
        raise ConfigurationError(f"scikit-bio {chosen} failed: {exc}") from exc
    if chosen == "osd":
        value = dict(zip(("observed", "singles", "doubles"), (int(v) for v in value)))
    else:
        value = float(value)
    return {
        "status": "completed",
        "operation": "calculate_alpha_diversity",
        "metric": chosen,
        "value": value,
        "n_features": len(vector),
        "total_count": sum(vector),
        "provenance": _provenance(method=f"skbio.diversity.alpha.{chosen}"),
    }


def calculate_unifrac_distance(
    counts: Sequence[Sequence[float]],
    taxa: Sequence[str],
    tree_newick: str,
    *,
    ids: Sequence[str] | None = None,
    weighted: bool = False,
    normalized: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the pairwise UniFrac distance matrix of samples on a phylogenetic tree.

    ``weighted=False`` gives unweighted UniFrac (``beta.unweighted_unifrac``); ``weighted=True``
    gives weighted UniFrac (``beta.weighted_unifrac``, with ``normalized`` selecting the normalised
    form). ``counts`` is samples x taxa; ``taxa`` are the tip names of the tree, one per column.
    """

    matrix = _matrix(counts, "counts", min_rows=2)
    labels = _ids(ids, len(matrix))
    names = _ids(taxa, len(matrix[0]), "taxa")
    weighted_flag = _flag(weighted, "weighted")
    norm_flag = _flag(normalized, "normalized")
    tree = _tree(tree_newick)
    diversity = require("skbio.diversity", _EXTRA)
    metric = "weighted_unifrac" if weighted_flag else "unweighted_unifrac"
    kwargs: dict[str, Any] = {"taxa": names, "tree": tree}
    if weighted_flag:
        kwargs["normalized"] = norm_flag
    try:
        distance = diversity.beta_diversity(metric, matrix, ids=labels, **kwargs)
    except (ValueError, TypeError, RuntimeError) as exc:
        raise ConfigurationError(f"scikit-bio {metric} failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_unifrac_distance",
        "metric": metric,
        "normalized": norm_flag if weighted_flag else None,
        "ids": labels,
        "distance_matrix": distance.data.tolist(),
        "provenance": _provenance(method=f"skbio.diversity.beta_diversity({metric})"),
    }


# ---------------------------------------------------------------------------- resampling
def rarefy_count_vector(
    counts: Sequence[int],
    depth: int,
    *,
    replace: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly subsample a count vector to ``depth`` observations (``skbio.stats.subsample_counts``).

    ``replace=True`` draws with replacement. Without replacement ``depth`` cannot exceed the total.
    """

    vector = _vector(counts, "counts", integer=True)
    n = _int(depth, "depth", minimum=0, maximum=10**9)
    replace_flag = _flag(replace, "replace")
    seed_value = _seed(seed)
    stats = require("skbio.stats", _EXTRA)
    import numpy as np

    try:
        result = stats.subsample_counts(np.asarray(vector, dtype=int), n, replace=replace_flag, seed=seed_value)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio subsample_counts failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "rarefy_count_vector",
        "depth": n,
        "replace": replace_flag,
        "seed": seed_value,
        "counts": result.tolist(),
        "provenance": _provenance(method="skbio.stats.subsample_counts"),
    }


def subsample_items_per_bin(
    items: Sequence[str],
    maximum: int,
    *,
    minimum: int = 1,
    bins: Sequence[str] | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly keep up to ``maximum`` items per bin, dropping bins with fewer than ``minimum``.

    Uses ``skbio.stats.isubsample``. ``bins`` gives the bin label of each item (same length);
    without it all items share one bin. Returns ``selected`` as [bin, item] pairs.
    """

    if isinstance(items, (str, bytes)) or not isinstance(items, Sequence) or len(items) > _MAX_ITEMS:
        raise ConfigurationError("items must be a list of strings.")
    values = [_text(item, "items item", maximum=10_000) for item in items]
    top = _int(maximum, "maximum", minimum=1)
    low = _int(minimum, "minimum", minimum=1)
    if low > top:
        raise ConfigurationError("minimum cannot exceed maximum.")
    seed_value = _seed(seed)
    bin_labels = None
    if bins is not None:
        bin_labels = _ids_like(bins, len(values), "bins")
    stats = require("skbio.stats", _EXTRA)
    if bin_labels is None:
        stream: Any = values
        bin_f = None
    else:
        stream = list(zip(bin_labels, values))

        def bin_f(pair: tuple[str, str]) -> str:
            return pair[0]

    selected = []
    for bin_name, item in stats.isubsample(stream, top, minimum=low, bin_f=bin_f, seed=seed_value):
        selected.append([bin_name, item[1] if bin_labels is not None else item])
    return {
        "status": "completed",
        "operation": "subsample_items_per_bin",
        "maximum": top,
        "minimum": low,
        "seed": seed_value,
        "selected": [[str(b), str(i)] for b, i in selected],
        "provenance": _provenance(method="skbio.stats.isubsample"),
    }


def _ids_like(values: object, count: int, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or len(values) != count:
        raise ConfigurationError(f"{name} must be a list of {count} strings.")
    return [_text(v, f"{name} item", maximum=256) for v in values]


# ---------------------------------------------------------------------------- cospeciation
def run_hommola_cospeciation_test(
    host_distances: Sequence[Sequence[float]],
    parasite_distances: Sequence[Sequence[float]],
    interaction: Sequence[Sequence[int]],
    *,
    host_ids: Sequence[str] | None = None,
    parasite_ids: Sequence[str] | None = None,
    permutations: int = 999,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the Hommola et al. (2009) host-parasite cospeciation permutation test.

    ``interaction`` is a parasites x hosts 0/1 matrix (rows parasites, columns hosts, as scikit-bio expects). Returns the correlation, the permutation
    p-value and the number of interaction edges (``skbio.stats.evolve.hommola_cospeciation``).
    """

    host = _matrix(host_distances, "host_distances", min_rows=3)
    para = _matrix(parasite_distances, "parasite_distances", min_rows=3)
    link = _matrix(interaction, "interaction")
    if len(host) != len(host[0]) or len(para) != len(para[0]):
        raise ConfigurationError("distance matrices must be square.")
    if len(link) != len(para) or len(link[0]) != len(host):
        raise ConfigurationError("interaction must be parasites x hosts.")
    host_labels = _ids(host_ids, len(host), "host_ids")
    para_labels = _ids(parasite_ids, len(para), "parasite_ids")
    perms = _int(permutations, "permutations", minimum=0, maximum=100_000)
    seed_value = _seed(seed)
    skbio = require("skbio", _EXTRA)
    evolve = require("skbio.stats.evolve", _EXTRA)
    import numpy as np

    try:
        host_dm = skbio.DistanceMatrix(np.asarray(host), host_labels)
        para_dm = skbio.DistanceMatrix(np.asarray(para), para_labels)
        r, p, _null = evolve.hommola_cospeciation(
            host_dm, para_dm, np.asarray(link, dtype=int), permutations=perms, seed=seed_value
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio hommola_cospeciation failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "run_hommola_cospeciation_test",
        "correlation": float(r),
        "p_value": float(p),
        "interaction_edges": int(sum(sum(row) for row in link)),
        "permutations": perms,
        "seed": seed_value,
        "provenance": _provenance(method="skbio.stats.evolve.hommola_cospeciation"),
    }


# ---------------------------------------------------------------------------- gradient ANOVA
def calculate_gradient_anova(
    coordinates: Mapping[str, Sequence[float]],
    proportion_explained: Sequence[float],
    metadata: Mapping[str, Mapping[str, Any]],
    *,
    algorithm: str = "average",
    trajectory_categories: Sequence[str] | None = None,
    sort_category: str | None = None,
    axes: int = 3,
    weighted: bool = False,
    window_size: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute gradient-ANOVA trajectories of ordination coordinates by metadata category.

    ``algorithm`` selects the scikit-bio class: ``average`` (RMS average of the group; AverageGradientANOVA),
    ``trajectory`` (RMS trajectory length; TrajectoryGradientANOVA), ``first_difference``
    (FirstDifferenceGradientANOVA) or ``window_difference`` (WindowDifferenceGradientANOVA, uses
    ``window_size``). ``coordinates`` maps sample id to coordinates, ``metadata`` maps sample id to
    {category: value}; ``axes=0`` uses all axes.
    """

    chosen = _choice(algorithm, _GRADIENT_ALGORITHMS, "algorithm")
    if not isinstance(coordinates, Mapping) or not 2 <= len(coordinates) <= _MAX_SAMPLES:
        raise ConfigurationError("coordinates must map 2 to 2000 sample ids to coordinate lists.")
    ids = [_text(k, "sample id", maximum=256) for k in coordinates]
    rows = [_vector_any(v, "coordinates row") for v in coordinates.values()]
    if len({len(r) for r in rows}) != 1:
        raise ConfigurationError("coordinate rows must have equal length.")
    prop = _vector_any(proportion_explained, "proportion_explained")
    if len(prop) != len(rows[0]):
        raise ConfigurationError("proportion_explained needs one value per coordinate axis.")
    if not isinstance(metadata, Mapping) or not metadata:
        raise ConfigurationError("metadata must map sample ids to category mappings.")
    axes_count = _int(axes, "axes", minimum=0, maximum=len(prop))
    weighted_flag = _flag(weighted, "weighted")
    window = _int(window_size, "window_size", minimum=1, maximum=1000)
    categories = None if trajectory_categories is None else [str(c) for c in trajectory_categories]
    gradient = require("skbio.stats.gradient", _EXTRA)
    pd = require("pandas", _EXTRA)
    import numpy as np

    coords = pd.DataFrame(np.asarray(rows), index=ids)
    meta = pd.DataFrame.from_dict({k: dict(v) for k, v in metadata.items()}, orient="index")
    classes = {
        "average": gradient.AverageGradientANOVA,
        "trajectory": gradient.TrajectoryGradientANOVA,
        "first_difference": gradient.FirstDifferenceGradientANOVA,
        "window_difference": gradient.WindowDifferenceGradientANOVA,
    }
    kwargs: dict[str, Any] = {
        "trajectory_categories": categories,
        "sort_category": sort_category,
        "axes": axes_count,
        "weighted": weighted_flag,
    }
    if chosen == "window_difference":
        kwargs["window_size"] = window
    try:
        result = classes[chosen](coords, np.asarray(prop), meta, **kwargs).get_trajectories()
    except (ValueError, TypeError, KeyError) as exc:
        raise ConfigurationError(f"scikit-bio gradient ANOVA failed: {exc}") from exc
    output = []
    for item in result.categories:
        groups = {}
        for group in item.groups or []:
            groups[str(group.name)] = {
                "trajectory": _plain(group.trajectory),
                "mean": float(group.mean),
                "info": _plain(group.info),
                "message": group.message,
            }
        output.append(
            {
                "category": item.category,
                "probability": None if item.probability is None or item.probability != item.probability else float(item.probability),
                "message": item.message,
                "groups": groups,
            }
        )
    return {
        "status": "completed",
        "operation": "calculate_gradient_anova",
        "algorithm": result.algorithm,
        "weighted": bool(result.weighted),
        "categories": output,
        "provenance": _provenance(method=f"skbio.stats.gradient ({chosen})"),
    }


def _vector_any(values: object, name: str) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not 1 <= len(values) <= 1000:
        raise ConfigurationError(f"{name} must be a list of 1 to 1000 numbers.")
    return [_number(v, f"{name} item") for v in values]


# ---------------------------------------------------------------------------- table augmentation
def _augment_inputs(table: object, n: int, labels: object, seed: object) -> tuple[Any, int, Any, int | None]:
    import numpy as np

    matrix = _matrix(table, "table", min_rows=2)
    count = _int(n, "n", minimum=1, maximum=100_000)
    return np.asarray(matrix, dtype=float), count, _labels(labels, len(matrix)), _seed(seed)


def mixup_table(
    table: Sequence[Sequence[float]],
    n: int,
    *,
    labels: Sequence[int] | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate ``n`` synthetic samples by vanilla mixup of random sample pairs (``skbio.table.mixup``)."""

    matrix, count, label_array, seed_value = _augment_inputs(table, n, labels, seed)
    shape = _number(alpha, "alpha", minimum=0.0, exclusive=True)
    skbio_table = require("skbio.table", _EXTRA)
    try:
        aug, aug_labels = skbio_table.mixup(
            matrix, count, labels=label_array, intra_class=_flag(intra_class, "intra_class"),
            alpha=shape, append=_flag(append, "append"), seed=seed_value,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio mixup failed: {exc}") from exc
    return _augmented("mixup_table", aug, aug_labels, alpha=shape, seed=seed_value)


def aitchison_mixup_table(
    table: Sequence[Sequence[float]],
    n: int,
    *,
    labels: Sequence[int] | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    normalize: bool = True,
    append: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate ``n`` synthetic compositions by mixup in Aitchison geometry (``skbio.table.aitchison_mixup``)."""

    matrix, count, label_array, seed_value = _augment_inputs(table, n, labels, seed)
    shape = _number(alpha, "alpha", minimum=0.0, exclusive=True)
    skbio_table = require("skbio.table", _EXTRA)
    try:
        aug, aug_labels = skbio_table.aitchison_mixup(
            matrix, count, labels=label_array, intra_class=_flag(intra_class, "intra_class"), alpha=shape,
            normalize=_flag(normalize, "normalize"), append=_flag(append, "append"), seed=seed_value,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio aitchison_mixup failed: {exc}") from exc
    return _augmented("aitchison_mixup_table", aug, aug_labels, alpha=shape, seed=seed_value)


def compositional_cutmix_table(
    table: Sequence[Sequence[float]],
    n: int,
    *,
    labels: Sequence[int] | None = None,
    normalize: bool = True,
    append: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate ``n`` synthetic compositions by compositional cutmix within classes (``skbio.table.compos_cutmix``)."""

    matrix, count, label_array, seed_value = _augment_inputs(table, n, labels, seed)
    skbio_table = require("skbio.table", _EXTRA)
    try:
        aug, aug_labels = skbio_table.compos_cutmix(
            matrix, count, labels=label_array, normalize=_flag(normalize, "normalize"),
            append=_flag(append, "append"), seed=seed_value,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio compos_cutmix failed: {exc}") from exc
    return _augmented("compositional_cutmix_table", aug, aug_labels, seed=seed_value)


def phylomix_table(
    table: Sequence[Sequence[float]],
    n: int,
    tree_newick: str,
    taxa: Sequence[str],
    *,
    labels: Sequence[int] | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate ``n`` synthetic samples by PhyloMix, mixing the taxa below a random tree node (``skbio.table.phylomix``)."""

    matrix, count, label_array, seed_value = _augment_inputs(table, n, labels, seed)
    names = _ids(taxa, matrix.shape[1], "taxa")
    tree = _tree(tree_newick)
    shape = _number(alpha, "alpha", minimum=0.0, exclusive=True)
    skbio_table = require("skbio.table", _EXTRA)
    try:
        aug, aug_labels = skbio_table.phylomix(
            matrix, count, tree, taxa=names, labels=label_array,
            intra_class=_flag(intra_class, "intra_class"), alpha=shape,
            append=_flag(append, "append"), seed=seed_value,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio phylomix failed: {exc}") from exc
    return _augmented("phylomix_table", aug, aug_labels, alpha=shape, seed=seed_value)


# ---------------------------------------------------------------------------- metadata
def _sample_metadata(table: object, name: str) -> Any:
    if not isinstance(table, Mapping) or not table or len(table) > _MAX_SAMPLES:
        raise ConfigurationError(f"{name} must map 1 to {_MAX_SAMPLES} sample ids to {{column: value}} mappings.")
    pd = require("pandas", _EXTRA)
    skbio = require("skbio", _EXTRA)
    rows = {}
    for sample, record in table.items():
        if not isinstance(record, Mapping) or not record:
            raise ConfigurationError(f"{name} entries must be non-empty column mappings.")
        rows[_text(sample, "sample id", maximum=256)] = dict(record)
    frame = pd.DataFrame.from_dict(rows, orient="index")
    frame.index.name = "id"
    try:
        return skbio.SampleMetadata(frame)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"{name} is not valid sample metadata: {exc}") from exc


def merge_sample_metadata(
    tables: Sequence[Mapping[str, Mapping[str, Any]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sample-metadata tables on their shared sample ids (``skbio.SampleMetadata.merge``).

    Each table maps sample id to {column: value}; ids not present in every table are dropped
    and column names must differ between tables.
    """

    if isinstance(tables, (str, bytes)) or not isinstance(tables, Sequence) or not 2 <= len(tables) <= 20:
        raise ConfigurationError("tables must be a list of 2 to 20 metadata mappings.")
    objects = [_sample_metadata(t, "tables item") for t in tables]
    try:
        merged = objects[0].merge(*objects[1:])
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"scikit-bio SampleMetadata.merge failed: {exc}") from exc
    frame = merged.to_dataframe()
    records = {
        str(idx): {c: _plain(v) for c, v in row.items()} for idx, row in frame.to_dict(orient="index").items()
    }
    return {
        "status": "completed",
        "operation": "merge_sample_metadata",
        "ids": [str(i) for i in frame.index],
        "columns": [str(c) for c in frame.columns],
        "metadata": records,
        "provenance": _provenance(method="skbio.SampleMetadata.merge"),
    }


def query_sample_ids(
    metadata: Mapping[str, Mapping[str, Any]],
    *,
    where: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the sample ids that satisfy an SQLite ``WHERE`` clause (``skbio.SampleMetadata.get_ids``).

    Columns are referenced by name, quoted with double quotes when needed. Without ``where`` all
    ids are returned.
    """

    obj = _sample_metadata(metadata, "metadata")
    clause = None if where is None else _text(where, "where", maximum=2_000)
    try:
        found = obj.get_ids(clause)
    except Exception as exc:  # noqa: BLE001 - sqlite errors are wrapped by scikit-bio
        raise ConfigurationError(f"scikit-bio SampleMetadata.get_ids failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "query_sample_ids",
        "where": clause,
        "ids": sorted(str(i) for i in found),
        "provenance": _provenance(method="skbio.SampleMetadata.get_ids"),
    }


def _interval_metadata(spec: object, name: str) -> Any:
    if not isinstance(spec, Mapping):
        raise ConfigurationError(f"{name} must be a mapping with upper_bound and intervals.")
    metadata_module = require("skbio.metadata", _EXTRA)
    upper = _int(spec.get("upper_bound"), f"{name}.upper_bound", minimum=0)
    intervals = spec.get("intervals", [])
    if isinstance(intervals, (str, bytes)) or not isinstance(intervals, Sequence) or len(intervals) > 100_000:
        raise ConfigurationError(f"{name}.intervals must be a list.")
    obj = metadata_module.IntervalMetadata(upper)
    for item in intervals:
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"{name}.intervals items must be mappings with bounds and metadata.")
        bounds = [tuple(_int(v, "bound", minimum=0) for v in pair) for pair in item.get("bounds", [])]
        if not bounds or any(len(b) != 2 for b in bounds):
            raise ConfigurationError("each interval needs bounds as a list of [start, end] pairs.")
        attrs = item.get("metadata") or {}
        if not isinstance(attrs, Mapping):
            raise ConfigurationError("interval metadata must be a mapping.")
        try:
            obj.add(bounds, metadata=dict(attrs))
        except (ValueError, TypeError) as exc:
            raise ConfigurationError(f"invalid interval {item!r}: {exc}") from exc
    return obj


def _interval_rows(obj: Any) -> list[dict[str, Any]]:
    return [
        {"bounds": [list(map(int, b)) for b in iv.bounds], "metadata": _plain(dict(iv.metadata))}
        for iv in obj.query()
    ]


def query_interval_metadata(
    interval_metadata: Mapping[str, Any],
    *,
    bounds: Sequence[Sequence[int]] | None = None,
    attributes: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find interval features that overlap given spans and carry given attributes (``IntervalMetadata.query``).

    ``interval_metadata`` is {"upper_bound": int, "intervals": [{"bounds": [[start, end]], "metadata": {...}}]}
    with half-open 0-based coordinates. A feature matches if it overlaps any of ``bounds`` and has
    every key/value of ``attributes``.
    """

    obj = _interval_metadata(interval_metadata, "interval_metadata")
    spans = None
    if bounds is not None:
        spans = [tuple(_int(v, "bounds item", minimum=0) for v in pair) for pair in bounds]
        if any(len(s) != 2 for s in spans):
            raise ConfigurationError("bounds must be a list of [start, end] pairs.")
    if attributes is not None and not isinstance(attributes, Mapping):
        raise ConfigurationError("attributes must be a mapping.")
    hits = list(obj.query(spans, None if attributes is None else dict(attributes)))
    return {
        "status": "completed",
        "operation": "query_interval_metadata",
        "n_matches": len(hits),
        "matches": [
            {"bounds": [list(map(int, b)) for b in iv.bounds], "metadata": _plain(dict(iv.metadata))}
            for iv in hits
        ],
        "provenance": _provenance(method="skbio.metadata.IntervalMetadata.query"),
    }


def concat_interval_metadata(
    interval_metadata: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate interval annotations of several sequences into one coordinate space (``IntervalMetadata.concat``).

    The coordinates of each later object are shifted by the summed upper bounds of the earlier
    ones. Input objects use the format of :func:`query_interval_metadata`.
    """

    if isinstance(interval_metadata, (str, bytes)) or not isinstance(interval_metadata, Sequence) \
            or not 1 <= len(interval_metadata) <= 1000:
        raise ConfigurationError("interval_metadata must be a list of 1 to 1000 objects.")
    objects = [_interval_metadata(spec, "interval_metadata item") for spec in interval_metadata]
    metadata_module = require("skbio.metadata", _EXTRA)
    merged = metadata_module.IntervalMetadata.concat(objects)
    return {
        "status": "completed",
        "operation": "concat_interval_metadata",
        "upper_bound": int(merged.upper_bound),
        "intervals": _interval_rows(merged),
        "provenance": _provenance(method="skbio.metadata.IntervalMetadata.concat"),
    }


# ---------------------------------------------------------------------------- embeddings
def ordinate_embedding_vectors(
    vectors: Mapping[str, Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Ordinate embedding vectors (for example sequence-model embeddings) by singular value decomposition.

    ``vectors`` maps an object id to its embedding. Runs ``skbio.embedding.embed_vec_to_ordination``
    and returns the sample coordinates, singular values and proportion explained.
    """

    if not isinstance(vectors, Mapping) or not 3 <= len(vectors) <= _MAX_SAMPLES:
        raise ConfigurationError("vectors must map 3 to 2000 ids to embedding lists.")
    ids = [_text(k, "vector id", maximum=256) for k in vectors]
    rows = [_vector_any(v, "vector") for v in vectors.values()]
    if len({len(r) for r in rows}) != 1:
        raise ConfigurationError("all vectors must have equal length.")
    embedding = require("skbio.embedding", _EXTRA)
    import numpy as np

    try:
        items = [embedding.EmbeddingVector(np.asarray(row, dtype=float).reshape(1, -1), [ident.encode('ascii', 'replace')]) for row, ident in zip(rows, ids)]
        result = embedding.embed_vec_to_ordination(items)
    except (ValueError, TypeError, AttributeError) as exc:
        raise ConfigurationError(f"scikit-bio embedding ordination failed: {exc}") from exc
    samples = result.samples
    return {
        "status": "completed",
        "operation": "ordinate_embedding_vectors",
        "ids": ids,
        "coordinates": samples.to_numpy().tolist(),
        "axes": [str(c) for c in samples.columns],
        "eigvals": np.asarray(result.eigvals).tolist(),
        "proportion_explained": np.asarray(result.proportion_explained).tolist(),
        "provenance": _provenance(method="skbio.embedding.embed_vec_to_ordination"),
    }


__all__ = [
    "calculate_alpha_diversity",
    "calculate_unifrac_distance",
    "rarefy_count_vector",
    "subsample_items_per_bin",
    "run_hommola_cospeciation_test",
    "calculate_gradient_anova",
    "mixup_table",
    "aitchison_mixup_table",
    "compositional_cutmix_table",
    "phylomix_table",
    "merge_sample_metadata",
    "query_sample_ids",
    "query_interval_metadata",
    "concat_interval_metadata",
    "ordinate_embedding_vectors",
]

_register_public_verbose(globals(), tuple(__all__))
