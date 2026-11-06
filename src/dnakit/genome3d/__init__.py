"""Local candidate adapters for the genome3d capability family."""

from __future__ import annotations

import itertools
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.genome3d import _impl_local as _impl


def aggregate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Coarsen a Hi-C contact matrix by summing ``factor`` x ``factor`` blocks of bins.

    ``input``: a square symmetric matrix (list of rows, ``{"matrix": ...}``, dense/sparse text or a path; see
    ``read``). ``options``: ``entity`` (only ``"contact_matrix"``); ``factor`` (required integer >= 2): bin ``k`` of
    the result sums original bins ``k*factor .. k*factor+factor-1`` (the last bin may cover fewer when ``n`` is not a
    multiple of ``factor``); ``bin_size`` (optional bp; the new bin size is returned); ``output`` (optional
    dense-TSV path). Contacts are conserved: the total of the result equals the total of the input.
    """
    opts = dict(options or {})
    _impl.opt_str(opts, "entity", "contact_matrix", ("contact_matrix",))
    matrix = _impl.parse_matrix(input, opts)
    if "factor" not in opts:
        raise ConfigurationError(
            "aggregate needs options['factor'] (integer >= 2).", code="GENOME3D_OPTION_INVALID"
        )
    factor = _impl.opt_int(opts, "factor", 2, 2, len(matrix))
    coarse = _impl.coarsen(matrix, factor)
    result: dict[str, object] = {
        "factor": factor,
        "n_bins": len(coarse),
        "original_n_bins": len(matrix),
    }
    if opts.get("bin_size") is not None:
        result["bin_size"] = _impl.opt_int(opts, "bin_size", 1, 1, 10**9) * factor
    return candidate_ok(
        "genome3d", "aggregate", input, _impl.with_matrix(result, coarse, opts), options=opts
    )


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate the bins of a contact matrix with TAD boundaries (insulation minima) or A/B compartment labels.

    ``input``: square symmetric matrix (see ``read``). ``options["entity"]``:
    ``"tad"`` (default): insulation score with ``window`` (default 5 bins; see ``calculate``), then bin ``i`` is a
    boundary iff its score is below ``threshold`` (default 0.0 = below the genome-wide average) and is the minimum
    within ``boundary_radius`` bins (default = window) on both sides, ties going to the leftmost bin. ``tads`` are
    bin intervals ``[start_bin, end_bin)`` (end exclusive). With boundaries ``b_1 < ... < b_m`` in a matrix of ``n``
    bins, the interior domains are ``[b_k, b_k+1)`` (``edge: false``); with ``include_edge_domains`` (default true)
    the two edge domains are reported too, flagged ``edge: true``: ``[0, b_1)`` (omitted when ``b_1 == 0``) and
    ``[b_m, n)`` (the last boundary opens the final domain). A single boundary therefore gives up to two edge
    domains, and with no boundary at all the whole matrix ``[0, n)`` is one edge domain. Edge domains are truncated
    by the matrix end, so their size is not comparable to interior domains; set ``include_edge_domains`` false to
    get only the interior domains (empty with fewer than two boundaries). Boundary coordinates are bin indices
    (``bin_size`` optional gives bp). This is an insulation-minimum rule, not a TAD-caller with statistical support.
    ``"compartment"``: A/B labels from the sign of the leading eigenvector of the observed/expected correlation matrix
    (see ``calculate``); optional ``reference`` per-bin track (e.g. GC content) orients A as positive.
    Needs numpy (optional dependency).
    """
    opts = dict(options or {})
    which = _impl.opt_str(opts, "entity", "tad", ("tad", "compartment"))
    matrix = _impl.parse_matrix(input, opts)
    np = _impl.numpy_or_none()
    if np is None:
        return candidate_not_implemented(
            "genome3d",
            "annotate",
            input,
            options=opts,
            requires="pip install numpy",
            message="Insulation and compartment calling use numpy, which is not installed.",
        )
    if which == "compartment":
        ref = opts.get("reference")
        comp = _impl.compartments(
            np, matrix, ref if isinstance(ref, Sequence) and not isinstance(ref, str) else None
        )
        return candidate_ok("genome3d", "annotate", input, {"entity": which, **comp}, options=opts)
    window = _impl.opt_int(opts, "window", 5, 1, len(matrix))
    radius = _impl.opt_int(opts, "boundary_radius", window, 1, len(matrix))
    threshold = _impl.opt_float(opts, "threshold", 0.0, -1e6, 1e6)
    ins = _impl.insulation(np, matrix, window)
    cuts = _impl.boundaries(ins["score"], radius, threshold)  # type: ignore[arg-type]
    size = _impl.opt_int(opts, "bin_size", 0, 0, 10**9) if "bin_size" in opts else 0
    tads = [{"start_bin": a, "end_bin": b, "edge": False} for a, b in itertools.pairwise(cuts)]
    if _impl.opt_bool(opts, "include_edge_domains", True):
        n_bins = len(matrix)
        if not cuts:
            tads = [{"start_bin": 0, "end_bin": n_bins, "edge": True}]
        else:
            if cuts[0] > 0:
                tads.insert(0, {"start_bin": 0, "end_bin": cuts[0], "edge": True})
            tads.append({"start_bin": cuts[-1], "end_bin": n_bins, "edge": True})
    result: dict[str, object] = {
        "entity": which,
        "window": window,
        "boundary_radius": radius,
        "threshold": threshold,
        "boundary_bins": cuts,
        "tads": tads,
    }
    if size:
        result["boundary_bp"] = [c * size for c in cuts]
        result["bin_size"] = size
    return candidate_ok("genome3d", "annotate", input, result, options=opts)


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate distance-decay curves, insulation scores or A/B compartment eigenvectors from a contact matrix.

    ``input``: square symmetric non-negative matrix (list of rows, ``{"matrix": ...}``, dense/sparse text or path).
    ``options["entity"]``:
    ``"distance_decay"`` (default): P(s) = mean contact over all valid bin pairs on diagonal ``s = |i - j|``
    (empty bins, with zero row sum, are excluded); ``bin_size`` (optional) adds ``distance_bp``; ``s_min`` (default
    2) and ``s_max`` (default n/2) give the log10-log10 least-squares ``slope`` of the curve.
    ``"insulation"``: ``window`` w (default 5): S(i) = mean of M[a][b] for a in [i-w, i-1], b in [i, i+w-1] (the
    square on the diagonal whose corner is the start of bin i; defined for w <= i <= n-w, else null) and
    ``score = log2(S(i) / mean(S))``; low scores mark insulating positions.
    ``"compartment"``: leading eigenvector of the Pearson correlation matrix of the observed/expected matrix
    (expected = mean per diagonal); sign gives A (>= 0) / B (< 0); orientation is arbitrary unless a per-bin
    ``reference`` (e.g. GC) is supplied, otherwise the largest-magnitude entry is made positive.
    Needs numpy (optional dependency).
    """
    opts = dict(options or {})
    which = _impl.opt_str(
        opts, "entity", "distance_decay", ("distance_decay", "insulation", "compartment")
    )
    matrix = _impl.parse_matrix(input, opts)
    np = _impl.numpy_or_none()
    if np is None:
        return candidate_not_implemented(
            "genome3d",
            "calculate",
            input,
            options=opts,
            requires="pip install numpy",
            message="Decay, insulation and compartment calculations use numpy, which is not installed.",
        )
    n = len(matrix)
    if which == "distance_decay":
        size = _impl.opt_int(opts, "bin_size", 0, 0, 10**9) if "bin_size" in opts else 0
        curve = _impl.distance_decay(np, matrix, size or None)
        result: dict[str, object] = {
            "entity": which,
            "curve": curve,
            "scaling": _impl.scaling_slope(
                curve,
                _impl.opt_int(opts, "s_min", 2, 1, n),
                _impl.opt_int(opts, "s_max", max(n // 2, 2), 1, n),
            ),
        }
    elif which == "insulation":
        result = {
            "entity": which,
            **_impl.insulation(np, matrix, _impl.opt_int(opts, "window", 5, 1, n)),
        }
    else:
        ref = opts.get("reference")
        result = {
            "entity": which,
            **_impl.compartments(
                np, matrix, ref if isinstance(ref, Sequence) and not isinstance(ref, str) else None
            ),
        }
    return candidate_ok("genome3d", "calculate", input, result, options=opts)


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two contact matrices by correlation of their upper-triangle contacts.

    ``input``: ``{"a": matrix, "b": matrix}`` or a list of two square matrices of the same size (see ``read``).
    ``options``: ``entity`` (only ``"contact_matrix"``); ``min_distance`` (default 1): only bin pairs with
    ``j - i >= min_distance`` are used; ``max_distance`` (optional); ``distance_normalize`` (default false): divide each
    entry by the mean of its diagonal (observed/expected) before correlating. Reports Pearson and Spearman
    (average ranks) correlation, the number of pairs and the mean absolute difference. Pairs where either matrix has an
    empty bin are skipped. This is a plain correlation, not stratum-adjusted (HiCRep).
    Needs numpy (optional dependency).
    """
    opts = dict(options or {})
    _impl.opt_str(opts, "entity", "contact_matrix", ("contact_matrix",))
    a, b = _impl.parse_pair(input, opts)
    if len(a) != len(b):
        raise InputFormatError(
            "The two matrices must have the same number of bins.", code="GENOME3D_INPUT_INVALID"
        )
    np = _impl.numpy_or_none()
    if np is None:
        return candidate_not_implemented(
            "genome3d",
            "compare",
            input,
            options=opts,
            requires="pip install numpy",
            message="Matrix comparison uses numpy, which is not installed.",
        )
    n = len(a)
    lo = _impl.opt_int(opts, "min_distance", 1, 0, n)
    hi = _impl.opt_int(opts, "max_distance", n, lo, n)
    ma, mb = np.array(a, dtype=float), np.array(b, dtype=float)
    va, vb = ma.sum(axis=1) > 0, mb.sum(axis=1) > 0
    valid = va & vb
    if _impl.opt_bool(opts, "distance_normalize", False):
        for m in (ma, mb):
            for s in range(n):
                diag = np.diagonal(m, offset=s)
                ok = valid[: n - s] & valid[s:]
                mean = diag[ok].mean() if ok.any() else 0.0
                if mean > 0:
                    idx = np.arange(n - s)
                    m[idx, idx + s] = diag / mean
    xs, ys = [], []
    for s in range(lo, hi + 1):
        if s >= n:
            break
        ok = valid[: n - s] & valid[s:]
        xs.append(np.diagonal(ma, offset=s)[ok])
        ys.append(np.diagonal(mb, offset=s)[ok])
    x, y = (np.concatenate(xs), np.concatenate(ys)) if xs else (np.array([]), np.array([]))
    if x.size < 2:
        raise InputFormatError(
            "Fewer than two comparable bin pairs.", code="GENOME3D_INPUT_INVALID"
        )
    rx, ry = np.array(_impl.average_ranks(x.tolist())), np.array(_impl.average_ranks(y.tolist()))
    result = {
        "pairs": int(x.size),
        "pearson": _impl.pearson(np, x, y),
        "spearman": _impl.pearson(np, rx, ry),
        "mean_absolute_difference": float(np.abs(x - y).mean()),
    }
    return candidate_ok("genome3d", "compare", input, result, options=opts)


def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise a contact matrix as dense TSV or sparse triplet text and return the text.

    ``input``: square symmetric matrix (see ``read``; ``format`` there selects how text input is parsed).
    ``options``: ``out_format`` ``"dense"`` (default: one tab-separated row per bin) or ``"sparse"`` (``i<TAB>j<TAB>value`` for the non-zero upper triangle including the diagonal,
    0-based bins); ``output`` (optional path, written when given; existing files need ``overwrite=True``).
    Result: ``text`` and ``output_path``. ``write`` is the same serialisation but requires an output path.
    """
    opts = dict(options or {})
    fmt = _impl.opt_str(opts, "out_format", "dense", ("dense", "sparse"))
    matrix = _impl.parse_matrix(input, opts)
    text = _impl.to_dense_text(matrix) if fmt == "dense" else _impl.to_sparse_text(matrix)
    written = None
    if opts.get("output") is not None:
        path = Path(str(opts["output"]))
        if path.exists() and not _impl.opt_bool(opts, "overwrite", False):
            raise ConfigurationError(
                f"Output file {path} already exists.",
                code="GENOME3D_OUTPUT_EXISTS",
                hint="Set options['overwrite']=True.",
            )
        path.write_text(text, encoding="utf-8")
        written = str(path)
    return candidate_ok(
        "genome3d",
        "export",
        input,
        {"format": fmt, "text": text, "output_path": written, **_impl.matrix_info(matrix)},
        options=opts,
    )


def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise (balance) a Hi-C contact matrix by iterative correction (ICE) or vanilla coverage.

    ``input``: square symmetric non-negative matrix (see ``read``). ``options``: ``entity`` (only ``"contact_matrix"``);
    ``method``: ``"ice"`` (default; Imakaev et al. 2012) or ``"coverage"`` (M_ij / (r_i r_j) with r_i the row sum divided
    by the mean row sum). ICE parameters: ``max_iter`` (default 200), ``tol`` (default 1e-5, on max |s_i - 1| of
    the relative row sums), ``ignore_diagonals`` (default 0: number of near-diagonals zeroed first, 1 = main
    diagonal only), ``min_nonzero`` (default 1: bins with fewer non-zero entries are masked), ``target_row_sum`` (default
    1.0: valid rows of the result sum to this value, as in cooler). Result: balanced ``matrix`` (symmetric),
    ``bias`` (the matrix equals raw / (bias_i bias_j); null for masked bins), ``masked_bins``, ``iterations``,
    ``converged`` and ``max_row_sum_deviation``. Matrices above 1000 bins require ``output``.
    Needs numpy (optional dependency).
    """
    opts = dict(options or {})
    _impl.opt_str(opts, "entity", "contact_matrix", ("contact_matrix",))
    method = _impl.opt_str(opts, "method", "ice", ("ice", "coverage"))
    matrix = _impl.parse_matrix(input, opts)
    np = _impl.numpy_or_none()
    if np is None:
        return candidate_not_implemented(
            "genome3d",
            "normalize",
            input,
            options=opts,
            requires="pip install numpy",
            message="Matrix balancing uses numpy, which is not installed.",
        )
    if method == "ice":
        out = _impl.ice_balance(
            np,
            matrix,
            _impl.opt_int(opts, "max_iter", 200, 1, 100_000),
            _impl.opt_float(opts, "tol", 1e-5, 1e-15, 1.0),
            _impl.opt_int(opts, "ignore_diagonals", 0, 0, len(matrix)),
            _impl.opt_int(opts, "min_nonzero", 1, 1, len(matrix)),
            _impl.opt_float(opts, "target_row_sum", 1.0, 1e-12, 1e12),
        )
    else:
        out = _impl.coverage_normalize(np, matrix)
    balanced = out.pop("matrix")
    result: dict[str, object] = {"entity": "contact_matrix", "method": method, **out}
    return candidate_ok(
        "genome3d", "normalize", input, _impl.with_matrix(result, balanced, opts), options=opts
    )  # type: ignore[arg-type]


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a contact matrix as an SVG heatmap (white to red on log10(1 + contacts)).

    ``input``: square symmetric matrix (see ``read``). ``options``: ``entity`` (only ``"contact_matrix"``);
    ``max_cells`` (default 150): larger matrices are coarsened by summing blocks until at most this many cells per side
    remain; ``output`` (optional ``.svg`` path). Each cell is a ``<rect class="cell">`` with ``data-i``, ``data-j`` and
    ``data-value``; the root carries ``data-bins`` and ``data-factor``. Result: ``{"svg": ...}``.
    """
    opts = dict(options or {})
    _impl.opt_str(opts, "entity", "contact_matrix", ("contact_matrix",))
    matrix = _impl.parse_matrix(input, opts)
    svg = _impl.svg_heatmap(matrix, _impl.opt_int(opts, "max_cells", 150, 2, 400))
    written = None
    if opts.get("output") is not None:
        Path(str(opts["output"])).write_text(svg, encoding="utf-8")
        written = str(opts["output"])
    return candidate_ok(
        "genome3d", "plot", input, {"svg": svg, "output_path": written}, options=opts
    )


def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a Hi-C contact matrix from dense or sparse text, a text file path, or an in-memory matrix.

    ``input``: list of rows, ``{"matrix": ...}``, matrix text or a path. ``options``: ``format``: ``"dense"``
    (default; whitespace/comma/tab-separated rows, a non-numeric header row and/or label column are skipped) or
    ``"sparse"`` (lines ``bin_i bin_j value``, 0-based, every unordered pair listed once, mirrored into a symmetric
    matrix; ``n_bins`` optional). Values must be finite and non-negative and the matrix square and symmetric (up to
    1e-9). Binary ``.hic``/``.cool``/``.mcool`` files are rejected with a hint to use the cooler wrappers. Result:
    ``n_bins``, ``total_contacts``, ``nonzero_entries`` and ``matrix`` (matrices above 1000 bins need ``output``).
    """
    opts = dict(options or {})
    matrix = _impl.parse_matrix(input, opts)
    return candidate_ok(
        "genome3d",
        "read",
        input,
        _impl.with_matrix(_impl.matrix_info(matrix), matrix, opts),
        options=opts,
    )


def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Write a contact matrix to a dense TSV or sparse triplet file.

    ``input``: square symmetric matrix (see ``read``). ``options``: ``output`` (required file path); ``out_format``
    ``"dense"`` (default) or ``"sparse"`` (see ``export``); ``overwrite`` (default false: an existing file is an
    error). Result: ``output_path``, ``bytes_written``, ``format`` and the matrix summary.
    """
    opts = dict(options or {})
    if opts.get("output") is None:
        raise ConfigurationError(
            "write needs options['output'] (a file path).", code="GENOME3D_OPTION_INVALID"
        )
    fmt = _impl.opt_str(opts, "out_format", "dense", ("dense", "sparse"))
    matrix = _impl.parse_matrix(input, opts)
    path = Path(str(opts["output"]))
    if path.exists() and not _impl.opt_bool(opts, "overwrite", False):
        raise ConfigurationError(
            f"Output file {path} already exists.",
            code="GENOME3D_OUTPUT_EXISTS",
            hint="Set options['overwrite']=True.",
        )
    text = _impl.to_dense_text(matrix) if fmt == "dense" else _impl.to_sparse_text(matrix)
    path.write_text(text, encoding="utf-8")
    return candidate_ok(
        "genome3d",
        "write",
        input,
        {
            "format": fmt,
            "output_path": str(path),
            "bytes_written": len(text.encode("utf-8")),
            **_impl.matrix_info(matrix),
        },
        options=opts,
    )


__all__ = [
    "aggregate",
    "annotate",
    "calculate",
    "compare",
    "export",
    "normalize",
    "plot",
    "read",
    "write",
]

__candidate_adapter_operations__ = set(
    [
        "aggregate",
        "annotate",
        "calculate",
        "compare",
        "export",
        "normalize",
        "plot",
        "read",
        "write",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
