"""Embedding-based single-cell and spatial operations that need numpy: PCA, k-means, graphs, label transfer."""

from __future__ import annotations

import heapq
import math
from collections import Counter

from dnakit.exceptions import ConfigurationError, InputFormatError

from ._impl_de import log1p_columns
from ._impl_matrix import (
    CountMatrix,
    labels_for,
    numpy_module,
    opt_bool,
    opt_float,
    opt_int,
    parse_coordinates,
    parse_matrix,
    require_option,
)

Result = tuple[object, list[str]]
MAX_CELLS = 20_000  # k-means only needs cells x k distances
MAX_PAIRWISE_CELLS = 5_000  # operations that build a cells x cells distance matrix
MAX_CELLS_AUTO = 5_000


def _check_size(cm: CountMatrix, limit: int) -> None:
    if cm.n_cells > limit:
        raise InputFormatError(
            f"{cm.n_cells} cells exceed the {limit}-cell limit of this pairwise-distance operation.",
            code="INPUT_TOO_LARGE",
            hint="Subsample cells first.",
        )


def sq_distances(a, b):
    """Matrix of squared Euclidean distances between the rows of ``a`` and ``b`` (memory ~ len(a) * len(b))."""

    np = numpy_module()
    d = (a * a).sum(axis=1)[:, None] + (b * b).sum(axis=1)[None, :] - 2.0 * (a @ b.T)
    return np.maximum(d, 0.0)


def pca_embedding(cm: CountMatrix, options: dict[str, object]):
    """Cells x PCs from ln(1 + CP10k) values (optionally the ``n_top_genes`` highest-variance genes).

    Returns ``(points, explained_variance_ratio)``. PCs come from the SVD of the column-centred matrix;
    each PC's sign is fixed so that its largest-magnitude loading is positive. ``n_pcs=0`` skips PCA and
    returns the centred expression vectors.
    """

    np = numpy_module()
    log = log1p_columns(
        cm, opt_float(options, "target_sum", 1e4, 0.0, math.inf, exclusive_min=True)
    )
    x = np.array(log.values, dtype=float).T  # cells x genes
    n_top = opt_int(options, "n_top_genes", None, 1)
    if n_top is not None and n_top < x.shape[1]:
        keep = np.argsort(
            -x.var(axis=0, ddof=1) if x.shape[0] > 1 else -x.var(axis=0), kind="stable"
        )[:n_top]
        x = x[:, np.sort(keep)]
    x = x - x.mean(axis=0)
    if opt_bool(options, "scale", False):
        sd = x.std(axis=0, ddof=1) if x.shape[0] > 1 else np.ones(x.shape[1])
        sd[sd == 0] = 1.0
        x = x / sd
    n_pcs = opt_int(options, "n_pcs", 10, 0, 1000)
    if n_pcs == 0:
        return x, None
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    n_pcs = min(n_pcs, len(s))
    for k in range(n_pcs):
        if vt[k][np.argmax(np.abs(vt[k]))] < 0:
            vt[k] = -vt[k]
            u[:, k] = -u[:, k]
    points = u[:, :n_pcs] * s[:n_pcs]
    variance = s**2
    total = variance.sum()
    ratio = (variance[:n_pcs] / total) if total > 0 else np.zeros(n_pcs)
    return points, ratio


def kmeans(points, k: int, max_iter: int = 300):
    """Lloyd's k-means with deterministic farthest-point initialisation; labels are numbered by first appearance."""

    np = numpy_module()
    n = points.shape[0]
    centre0 = int(np.argmax(((points - points.mean(axis=0)) ** 2).sum(axis=1)))
    chosen = [centre0]
    nearest = ((points - points[centre0]) ** 2).sum(axis=1)
    while len(chosen) < k:
        nxt = int(np.argmax(nearest))
        if nearest[nxt] == 0:
            break  # fewer than k distinct points
        chosen.append(nxt)
        nearest = np.minimum(nearest, ((points - points[nxt]) ** 2).sum(axis=1))
    centres = points[chosen].copy()
    labels = np.full(n, -1, dtype=int)
    for _iteration in range(max_iter):
        distances = sq_distances(points, centres)
        new = distances.argmin(axis=1)
        if (new == labels).all():
            break
        labels = new
        for c in range(len(centres)):
            members = points[labels == c]
            if len(members):
                centres[c] = members.mean(axis=0)
    # renumber by first appearance
    mapping: dict[int, int] = {}
    for value in labels:
        mapping.setdefault(int(value), len(mapping))
    labels = np.array([mapping[int(v)] for v in labels])
    inertia = 0.0
    for c in range(len(mapping)):
        members = points[labels == c]
        inertia += float(((members - members.mean(axis=0)) ** 2).sum())
    return labels, inertia


def silhouette(points, labels) -> float:
    """Mean silhouette coefficient (Euclidean); a singleton cluster contributes 0."""

    np = numpy_module()
    dist = np.sqrt(sq_distances(points, points))
    scores = []
    ids = sorted(set(labels.tolist()))
    for i in range(len(labels)):
        own = labels == labels[i]
        if own.sum() <= 1:
            scores.append(0.0)
            continue
        a = dist[i][own].sum() / (own.sum() - 1)
        b = min(dist[i][labels == c].mean() for c in ids if c != labels[i])
        scores.append(0.0 if max(a, b) == 0 else (b - a) / max(a, b))
    return float(np.mean(scores))


def knn_indices(points, k: int):
    """Indices of the ``k`` nearest other points per row (ties broken by lower index)."""

    np = numpy_module()
    dist = sq_distances(points, points)
    np.fill_diagonal(dist, np.inf)
    return np.argsort(dist, axis=1, kind="stable")[:, :k]


def _choose_k(points, options: dict[str, object]) -> tuple[int, dict[str, float] | None]:
    n = points.shape[0]
    raw = options.get("n_clusters", "auto")
    if raw == "auto":
        if n > MAX_CELLS_AUTO:
            raise InputFormatError(
                "n_clusters='auto' supports up to 5000 cells.",
                code="INPUT_TOO_LARGE",
                hint="Set n_clusters explicitly.",
            )
        k_max = opt_int(options, "max_clusters", 10, 2, 100)
        candidates = list(range(2, min(k_max, n - 1) + 1))
        if not candidates:
            raise InputFormatError(
                "Too few cells to cluster.", code="TOO_FEW_CELLS", hint="Need at least 3 cells."
            )
        scores = {}
        for k in candidates:
            labels, _ = kmeans(points, k)
            scores[k] = silhouette(points, labels) if len(set(labels.tolist())) > 1 else -1.0
        best = max(candidates, key=lambda k: (scores[k], -k))
        return best, {str(k): v for k, v in scores.items()}
    k = opt_int({"n_clusters": raw}, "n_clusters", None, 1, 10_000)
    assert k is not None
    if k > n:
        raise ConfigurationError(
            f"n_clusters={k} exceeds the {n} cells.",
            code="INVALID_OPTION",
            hint="Use fewer clusters.",
        )
    return k, None


def _cluster_result(cm: CountMatrix, points, ratio, options, k, scores, extra_def: str) -> Result:
    labels, inertia = kmeans(points, k)
    names = [str(int(v)) for v in labels]
    sizes = Counter(names)
    result: dict[str, object] = {
        "definition": extra_def,
        "n_clusters": len(sizes),
        "labels": dict(zip(cm.cells, names, strict=False)),
        "cluster_sizes": dict(sorted(sizes.items(), key=lambda kv: int(kv[0]))),
        "inertia": inertia,
        "explained_variance_ratio": None if ratio is None else [float(v) for v in ratio],
        "silhouette_by_k": scores,
    }
    if opt_bool(options, "return_embedding", False):
        result["embedding"] = {
            cell: [float(v) for v in row] for cell, row in zip(cm.cells, points, strict=False)
        }
    return result, []


def op_cluster(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    _check_size(cm, MAX_CELLS)
    points, ratio = pca_embedding(cm, options)
    k, scores = _choose_k(points, options)
    return _cluster_result(
        cm,
        points,
        ratio,
        options,
        k,
        scores,
        "k-means (deterministic farthest-point start, Lloyd iterations) on the top PCs of ln(1+CP10k) values; k from options or the best mean silhouette",
    )


def spatial_smooth(points, coords, k: int, alpha: float):
    np = numpy_module()
    neighbours = knn_indices(np.array(coords, dtype=float), k)
    return (1.0 - alpha) * points + alpha * points[neighbours].mean(axis=1)


def op_cluster_spatial(value: object, options: dict[str, object]) -> Result:
    cm = parse_matrix(value, options)
    _check_size(cm, MAX_PAIRWISE_CELLS)
    coords = parse_coordinates(options.get("coordinates"), cm.cells, options)
    k_spatial = opt_int(options, "n_neighbors", 6, 1, 100)
    alpha = opt_float(options, "alpha", 0.5, 0.0, 1.0)
    if k_spatial >= cm.n_cells:  # type: ignore[operator]
        raise ConfigurationError(
            "options['n_neighbors'] must be smaller than the number of spots.",
            code="INVALID_OPTION",
            hint="Lower n_neighbors.",
        )
    points, ratio = pca_embedding(cm, options)
    smoothed = spatial_smooth(points, coords, k_spatial, alpha)  # type: ignore[arg-type]
    k, scores = _choose_k(smoothed, options)
    return _cluster_result(
        cm,
        smoothed,
        ratio,
        options,
        k,
        scores,
        f"spatial-domain heuristic: each spot's PC vector is replaced by (1-alpha)*own + alpha*mean of its {k_spatial} nearest spatial neighbours (alpha={alpha}), then k-means as in cluster()",
    )


def dijkstra(adjacency: dict[int, dict[int, float]], source: int, n: int) -> list[float]:
    distance = [math.inf] * n
    distance[source] = 0.0
    queue = [(0.0, source)]
    while queue:
        d, node = heapq.heappop(queue)
        if d > distance[node]:
            continue
        for other, weight in adjacency[node].items():
            candidate = d + weight
            if candidate < distance[other]:
                distance[other] = candidate
                heapq.heappush(queue, (candidate, other))
    return distance


def _trajectory(cm: CountMatrix, options: dict[str, object], spatial: bool) -> Result:
    _check_size(cm, MAX_PAIRWISE_CELLS)
    np = numpy_module()
    root_name = str(require_option(options, "root", "name of the root (start) cell"))
    index = cm.cell_index()
    if root_name not in index:
        raise InputFormatError(
            f"Root cell {root_name!r} is not in the matrix.",
            code="UNKNOWN_CELL",
            hint="Use a column name of the matrix.",
        )
    n = cm.n_cells
    k = opt_int(options, "n_neighbors", 15 if not spatial else 6, 1, 1000)
    if n < 2:
        raise InputFormatError(
            "A trajectory needs at least two cells.",
            code="TOO_FEW_CELLS",
            hint="Provide more cells.",
        )
    k = min(k, n - 1)  # type: ignore[type-var]
    points, _ratio = pca_embedding(cm, {**options, "n_pcs": options.get("n_pcs", 0)})
    expression = np.array(points, dtype=float)
    if spatial:
        coords = np.array(
            parse_coordinates(options.get("coordinates"), cm.cells, options), dtype=float
        )
        neighbours = knn_indices(coords, k)  # type: ignore[arg-type]
    else:
        neighbours = knn_indices(expression, k)  # type: ignore[arg-type]
    adjacency: dict[int, dict[int, float]] = {i: {} for i in range(n)}
    for i in range(n):
        for j in neighbours[i]:
            weight = float(np.linalg.norm(expression[i] - expression[int(j)]))
            for a, b in ((i, int(j)), (int(j), i)):
                if b not in adjacency[a] or weight < adjacency[a][b]:
                    adjacency[a][b] = weight
    distance = dijkstra(adjacency, index[root_name], n)
    finite = [d for d in distance if math.isfinite(d)]
    top = max(finite)
    warnings = []
    if len(finite) < n:
        warnings.append(
            f"{n - len(finite)} cells are not connected to the root in the kNN graph (pseudotime null); raise n_neighbors."
        )
    pseudotime = {
        cell: (d / top if top > 0 else 0.0) if math.isfinite(d) else None
        for cell, d in zip(cm.cells, distance, strict=False)
    }
    definition = (
        "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph built on spatial coordinates, with edge weights = Euclidean distance between the cells' expression vectors; divided by the maximum"
        if spatial
        else "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph of the cells' expression space (edge weight = Euclidean distance); divided by the maximum"
    )
    return {
        "definition": definition,
        "root": root_name,
        "n_neighbors": k,
        "pseudotime": pseudotime,
    }, warnings


def op_trajectory(value: object, options: dict[str, object]) -> Result:
    return _trajectory(parse_matrix(value, options), options, spatial=False)


def op_trajectory_spatial(value: object, options: dict[str, object]) -> Result:
    return _trajectory(parse_matrix(value, options), options, spatial=True)


def op_map_cells(value: object, options: dict[str, object]) -> Result:
    np = numpy_module()
    query = parse_matrix(value, options)
    if options.get("reference") is None:
        raise ConfigurationError(
            "options['reference'] (annotated reference matrix) is required.",
            code="MISSING_OPTION",
            hint="Pass genes x cells counts plus options['reference_labels'].",
        )
    reference = parse_matrix(
        options["reference"],
        {k: v for k, v in options.items() if k in {"delimiter", "max_matrix_elements"}},
    )
    ref_labels = labels_for(
        reference.cells,
        require_option(options, "reference_labels", "label of each reference cell"),
        "options['reference_labels']",
    )
    _check_size(query, MAX_PAIRWISE_CELLS)
    _check_size(reference, MAX_PAIRWISE_CELLS)
    ref_index = reference.gene_index()
    shared = [g for g in query.genes if g in ref_index]
    if len(shared) < 2:
        raise InputFormatError(
            f"Only {len(shared)} genes are shared between query and reference.",
            code="TOO_FEW_SHARED_GENES",
            hint="Use the same gene identifiers.",
        )
    target = opt_float(options, "target_sum", 1e4, 0.0, math.inf, exclusive_min=True)
    q_log, r_log = log1p_columns(query, target), log1p_columns(reference, target)
    q_idx = query.gene_index()
    q = np.array([q_log.values[q_idx[g]] for g in shared], dtype=float).T
    r = np.array([r_log.values[ref_index[g]] for g in shared], dtype=float).T
    k = min(opt_int(options, "n_neighbors", 5, 1, 1000), reference.n_cells)  # type: ignore[type-var]
    dist = np.sqrt(sq_distances(q, r))
    order = np.argsort(dist, axis=1, kind="stable")[:, :k]
    assignments: dict[str, dict[str, object]] = {}
    for i, cell in enumerate(query.cells):
        votes = Counter(ref_labels[int(j)] for j in order[i])
        top = max(votes.values())
        tied = {label for label, n in votes.items() if n == top}
        winner = next(ref_labels[int(j)] for j in order[i] if ref_labels[int(j)] in tied)
        assignments[cell] = {
            "label": winner,
            "confidence": top / k,
            "nearest_reference_cell": reference.cells[int(order[i][0])],
        }
    return {
        "definition": f"k-nearest-neighbour label transfer (k={k}): Euclidean distance on ln(1+CP10k) values of the {len(shared)} shared genes; majority vote, ties go to the label of the nearest tied neighbour; confidence = winning votes / k",
        "n_shared_genes": len(shared),
        "assignments": assignments,
    }, []
