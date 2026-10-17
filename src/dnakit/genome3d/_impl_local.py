"""Native helpers behind the genome3d candidate functions (contact matrices: parse, ICE, decay, insulation, compartments)."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_BINS = 3_000
MAX_RETURN_BINS = 1_000
BINARY_SUFFIXES = (".hic", ".cool", ".mcool", ".scool")


# ----------------------------------------------------------------------------- option helpers
def opt_str(options: Mapping[str, object], key: str, default: str | None, allowed: Sequence[str] | None = None) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"Option {key!r} must be a non-empty string.", code="GENOME3D_OPTION_INVALID")
    value = value.strip().lower()
    if allowed is not None and value not in allowed:
        raise ConfigurationError(
            f"Option {key!r}={value!r} is not supported.", code="GENOME3D_OPTION_INVALID", hint=f"Use one of: {', '.join(allowed)}."
        )
    return value


def opt_int(options: Mapping[str, object], key: str, default: int, minimum: int, maximum: int) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"Option {key!r} must be an integer in [{minimum}, {maximum}].", code="GENOME3D_OPTION_INVALID")
    return value


def opt_float(options: Mapping[str, object], key: str, default: float, minimum: float, maximum: float) -> float:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(f"Option {key!r} must be a number in [{minimum}, {maximum}].", code="GENOME3D_OPTION_INVALID")
    return float(value)


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(f"Option {key!r} must be a boolean.", code="GENOME3D_OPTION_INVALID")
    return value


def numpy_or_none() -> Any:
    try:
        import numpy
    except ImportError:
        return None
    return numpy


# ----------------------------------------------------------------------------- parsing
def _is_file(text: str) -> bool:
    if "\n" in text or "\t" in text or len(text) >= 4096:
        return False
    try:
        return Path(text).is_file()
    except OSError:
        return False


def _read_text(value: object) -> str:
    if isinstance(value, (str, Path)) and str(value).lower().endswith(BINARY_SUFFIXES):
        raise ConfigurationError(
            "Binary .hic/.cool/.mcool files are not read by genome3d.",
            code="GENOME3D_BINARY_FORMAT",
            hint="Use the cooler/hic wrappers in dnakit.epigenome, or export a dense/sparse text matrix first.",
        )
    if isinstance(value, Path) or (isinstance(value, str) and _is_file(value)):
        path = Path(value)
        if path.stat().st_size > 200_000_000:
            raise ConfigurationError("Matrix file is too large.", code="GENOME3D_INPUT_TOO_LARGE")
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    raise InputFormatError("Expected matrix text, a path or an in-memory matrix.", code="GENOME3D_INPUT_INVALID")


def _num(token: object, where: str) -> float:
    if isinstance(token, bool):
        raise InputFormatError(f"Matrix values must be numbers ({where}).", code="GENOME3D_MATRIX_INVALID")
    try:
        value = float(token)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise InputFormatError(f"Non-numeric matrix value {token!r} ({where}).", code="GENOME3D_MATRIX_INVALID") from exc
    if math.isnan(value) or math.isinf(value) or value < 0:
        raise InputFormatError(f"Contacts must be finite and non-negative ({where}).", code="GENOME3D_MATRIX_INVALID")
    return value


def _check_square(matrix: list[list[float]], symmetric_tol: float = 1e-9) -> list[list[float]]:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise InputFormatError("A contact matrix must be square and non-empty.", code="GENOME3D_MATRIX_INVALID")
    if n > MAX_BINS:
        raise ConfigurationError(f"At most {MAX_BINS} bins are supported.", code="GENOME3D_INPUT_TOO_LARGE")
    for i in range(n):
        for j in range(i + 1, n):
            a, b = matrix[i][j], matrix[j][i]
            if abs(a - b) > symmetric_tol * max(1.0, abs(a), abs(b)):
                raise InputFormatError(f"The matrix is not symmetric at ({i}, {j}).", code="GENOME3D_MATRIX_ASYMMETRIC",
                                       hint="Provide a full symmetric matrix or sparse upper-triangle triplets.")
    return matrix


def parse_matrix(value: object, options: Mapping[str, object]) -> list[list[float]]:
    """Dense square matrix as lists; accepts a matrix, ``{"matrix": ...}`` or dense/sparse text or a path.

    Text: ``format="dense"`` (default; whitespace/comma/tab-separated rows; a non-numeric header row and/or label
    column are skipped) or ``format="sparse"`` (lines ``bin_i bin_j value``, 0-based, each unordered pair once;
    ``n_bins`` optional, default max index + 1; the matrix is mirrored).
    """

    if isinstance(value, Mapping):
        if "matrix" not in value:
            raise InputFormatError("Provide the matrix under the key 'matrix'.", code="GENOME3D_INPUT_INVALID")
        value = value["matrix"]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        matrix = [[_num(x, f"row {i}") for x in row] for i, row in enumerate(value)]  # type: ignore[union-attr]
        return _check_square(matrix)
    text = _read_text(value)
    fmt = opt_str(options, "format", "dense", ("dense", "sparse"))
    lines = [ln for ln in (raw.strip() for raw in text.splitlines()) if ln and not ln.startswith("#")]
    if not lines:
        raise InputFormatError("The matrix text is empty.", code="GENOME3D_MATRIX_INVALID")
    rows = [re.split(r"[\t, ]+", ln) for ln in lines]
    if fmt == "sparse":
        triples = []
        for k, cells in enumerate(rows, 1):
            if len(cells) != 3:
                raise InputFormatError(f"Sparse line {k} needs 'i j value'.", code="GENOME3D_MATRIX_INVALID")
            i, j = int(_num(cells[0], f"line {k}")), int(_num(cells[1], f"line {k}"))
            triples.append((i, j, _num(cells[2], f"line {k}")))
        n = opt_int(options, "n_bins", max(max(i, j) for i, j, _ in triples) + 1, 1, MAX_BINS)
        if any(i >= n or j >= n for i, j, _ in triples):
            raise InputFormatError("A bin index exceeds n_bins.", code="GENOME3D_MATRIX_INVALID")
        matrix = [[0.0] * n for _ in range(n)]
        seen: set[tuple[int, int]] = set()
        for i, j, v in triples:
            key = (min(i, j), max(i, j))
            if key in seen:
                raise InputFormatError(f"Pair {key} appears twice.", code="GENOME3D_MATRIX_INVALID")
            seen.add(key)
            matrix[i][j] = matrix[j][i] = v
        return matrix
    def numeric(cell: str) -> bool:
        try:
            float(cell)
            return True
        except ValueError:
            return False

    if not all(numeric(c) for c in rows[0]):
        rows = rows[1:]
    if rows and len(rows[0]) == len(rows) + 1 and not numeric(rows[0][0]):
        rows = [r[1:] for r in rows]
    matrix = [[_num(c, f"row {i}") for c in r] for i, r in enumerate(rows)]
    return _check_square(matrix)


def parse_pair(value: object, options: Mapping[str, object]) -> tuple[list[list[float]], list[list[float]]]:
    if isinstance(value, Mapping) and "a" in value and "b" in value:
        return parse_matrix(value["a"], options), parse_matrix(value["b"], options)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2:
        return parse_matrix(value[0], options), parse_matrix(value[1], options)
    raise InputFormatError("compare expects {'a': matrix, 'b': matrix} or a list of two matrices.", code="GENOME3D_INPUT_INVALID")


# ----------------------------------------------------------------------------- serialisation
def fmt_value(x: float) -> str:
    return str(int(x)) if x == int(x) and abs(x) < 1e15 else repr(float(x))


def to_dense_text(matrix: Sequence[Sequence[float]]) -> str:
    return "\n".join("\t".join(fmt_value(x) for x in row) for row in matrix) + "\n"


def to_sparse_text(matrix: Sequence[Sequence[float]]) -> str:
    """Upper triangle including the diagonal, non-zero entries only: ``i<TAB>j<TAB>value``."""

    n = len(matrix)
    return "".join(f"{i}\t{j}\t{fmt_value(matrix[i][j])}\n" for i in range(n) for j in range(i, n) if matrix[i][j] != 0)


def matrix_info(matrix: Sequence[Sequence[float]]) -> dict[str, object]:
    n = len(matrix)
    return {"n_bins": n, "total_contacts": sum(sum(r) for r in matrix) / 2 + sum(matrix[i][i] for i in range(n)) / 2,
            "nonzero_entries": sum(1 for r in matrix for x in r if x != 0)}


def with_matrix(result: dict[str, object], matrix: Sequence[Sequence[float]], options: Mapping[str, object],
                key: str = "matrix") -> dict[str, object]:
    """Attach a matrix to a result unless it is too large; optionally write it (dense TSV) to ``options['output']``."""

    if options.get("output") is not None:
        path = Path(str(options["output"]))
        if path.exists() and not opt_bool(options, "overwrite", False):
            raise ConfigurationError(f"Output file {path} already exists.", code="GENOME3D_OUTPUT_EXISTS",
                                     hint="Set options['overwrite']=True.")
        path.write_text(to_dense_text(matrix), encoding="utf-8")
        result["output_path"] = str(path)
    if len(matrix) <= MAX_RETURN_BINS:
        result[key] = [list(r) for r in matrix]
    else:
        if options.get("output") is None:
            raise ConfigurationError(f"Matrices above {MAX_RETURN_BINS} bins are only returned via options['output'].",
                                     code="GENOME3D_OUTPUT_REQUIRED")
        result[key + "_omitted"] = True
    return result


# ----------------------------------------------------------------------------- aggregation / normalisation
def coarsen(matrix: Sequence[Sequence[float]], factor: int) -> list[list[float]]:
    """Sum factor x factor blocks; the last block may be smaller when n is not a multiple of factor."""

    n = len(matrix)
    m = -(-n // factor)
    out = [[0.0] * m for _ in range(m)]
    for i in range(n):
        for j in range(n):
            out[i // factor][j // factor] += matrix[i][j]
    return out


def ice_balance(np: Any, matrix: Sequence[Sequence[float]], max_iter: int, tol: float, ignore_diagonals: int,
                min_nonzero: int, target_row_sum: float) -> dict[str, object]:
    """Iterative correction (Imakaev et al. 2012) of a symmetric non-negative matrix.

    Bins with fewer than ``min_nonzero`` non-zero entries (default 1) are masked. Each iteration computes the row sums
    s_i of the current matrix over valid bins, divides them by their mean, and divides M_ij by s_i s_j; the product of
    all s_i is the total bias b_i, so the result is M_ij / (b_i b_j). Iteration stops when max |s_i - 1| < tol.
    Finally the matrix is scaled so every valid row sums to ``target_row_sum``. Bias of masked bins is NaN (None).
    """

    m = np.array(matrix, dtype=float)
    n = len(m)
    if ignore_diagonals > 0:
        for d in range(-(ignore_diagonals - 1), ignore_diagonals):
            idx = np.arange(max(0, -d), min(n, n - d))
            m[idx, idx + d] = 0.0
    valid = (m > 0).sum(axis=1) >= min_nonzero
    if not valid.any():
        raise InputFormatError("All bins are empty.", code="GENOME3D_MATRIX_INVALID")
    m[~valid, :] = 0.0
    m[:, ~valid] = 0.0
    bias = np.ones(n)
    deviation = math.inf
    iterations = 0
    for iterations in range(1, max_iter + 1):
        sums = m.sum(axis=1)
        s = np.where(valid, sums / sums[valid].mean(), 1.0)
        deviation = float(np.abs(s[valid] - 1.0).max())
        if deviation < tol:
            break
        s = np.where(valid & (s > 0), s, 1.0)
        m = m / np.outer(s, s)
        bias = bias * s
    sums = m.sum(axis=1)
    scale = target_row_sum / sums[valid].mean()
    m = m * scale
    bias = bias / math.sqrt(scale)
    final_sums = m.sum(axis=1)
    return {
        "matrix": m.tolist(), "bias": [float(b) if v else None for b, v in zip(bias, valid)], "masked_bins": [int(i) for i in np.where(~valid)[0]],
        "iterations": iterations, "converged": deviation < tol, "max_row_sum_deviation": float(np.abs(final_sums[valid] - target_row_sum).max()),
    }


def coverage_normalize(np: Any, matrix: Sequence[Sequence[float]]) -> dict[str, object]:
    """Vanilla coverage: M_ij / (r_i r_j) * (mean r)^2 with r_i the row sum; empty bins stay 0."""

    m = np.array(matrix, dtype=float)
    r = m.sum(axis=1)
    valid = r > 0
    mean = r[valid].mean() if valid.any() else 1.0
    w = np.where(valid, r / mean, 1.0)
    out = m / np.outer(w, w)
    return {"matrix": out.tolist(), "bias": [float(x) if v else None for x, v in zip(w, valid)],
            "masked_bins": [int(i) for i in np.where(~valid)[0]]}


# ----------------------------------------------------------------------------- derived quantities
def masked_nan(np: Any, matrix: Sequence[Sequence[float]]) -> tuple[Any, Any]:
    m = np.array(matrix, dtype=float)
    valid = m.sum(axis=1) > 0
    m[~valid, :] = np.nan
    m[:, ~valid] = np.nan
    return m, valid


def distance_decay(np: Any, matrix: Sequence[Sequence[float]], bin_size: int | None) -> list[dict[str, object]]:
    """Mean contact per diagonal s = |i - j| over valid (non-empty) bin pairs; s = 0 is the main diagonal."""

    m, _ = masked_nan(np, matrix)
    n = len(m)
    out = []
    for s in range(n):
        values = np.diagonal(m, offset=s)
        values = values[~np.isnan(values)]
        row: dict[str, object] = {"distance_bins": s, "n_pairs": int(values.size),
                                  "mean_contact": float(values.mean()) if values.size else None}
        if bin_size:
            row["distance_bp"] = s * bin_size
        out.append(row)
    return out


def scaling_slope(decay: Sequence[Mapping[str, object]], s_min: int, s_max: int) -> dict[str, object]:
    """Least-squares slope of log10(mean_contact) vs log10(distance in bins) for s_min <= s <= s_max (s >= 1)."""

    pts = [(math.log10(int(r["distance_bins"])), math.log10(float(r["mean_contact"])))  # type: ignore[call-overload]
           for r in decay if s_min <= int(r["distance_bins"]) <= s_max and r["mean_contact"]]  # type: ignore[call-overload]
    if len(pts) < 2:
        return {"slope": None, "points": len(pts)}
    mx, my = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    return {"slope": sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx if sxx else None, "points": len(pts)}


def insulation(np: Any, matrix: Sequence[Sequence[float]], window: int) -> dict[str, object]:
    """Insulation score of each bin i (window w): S(i) = mean of M[a][b] for a in [i-w, i-1], b in [i, i+w-1].

    The square sits on the diagonal with its corner at the start of bin i (the boundary position), so i ranges
    over [w, n - w]; other bins are undefined (None). Masked bins are ignored in the mean. The score is
    log2(S(i) / mean over defined bins of S): negative values mean fewer contacts across i than average.
    """

    m, _ = masked_nan(np, matrix)
    n = len(m)
    raw: list[float | None] = [None] * n
    for i in range(window, n - window + 1):
        block = m[i - window : i, i : i + window]
        finite = block[~np.isnan(block)]
        if finite.size:
            raw[i] = float(finite.mean())
    defined = [v for v in raw if v is not None]
    mean_all = sum(defined) / len(defined) if defined else None
    score = [None if v is None or not mean_all or v <= 0 else math.log2(v / mean_all) for v in raw]
    return {"window": window, "mean_contact": raw, "score": score}


def boundaries(score: Sequence[float | None], radius: int, threshold: float) -> list[int]:
    """Bin i is a boundary iff score[i] < threshold and score[i] <= every defined score within +-radius, with ties
    resolved to the leftmost bin (a strictly smaller score is required on the left)."""

    out = []
    n = len(score)
    for i, v in enumerate(score):
        if v is None or v >= threshold:
            continue
        left = [score[j] for j in range(max(0, i - radius), i) if score[j] is not None]
        right = [score[j] for j in range(i + 1, min(n, i + radius + 1)) if score[j] is not None]
        if all(v < x for x in left) and all(v <= x for x in right):  # type: ignore[operator]
            out.append(i)
    return out


def compartments(np: Any, matrix: Sequence[Sequence[float]], reference: Sequence[float] | None) -> dict[str, object]:
    """A/B compartments from the sign of the leading eigenvector of the O/E Pearson correlation matrix.

    1) mask empty bins; 2) expected(s) = mean contact on diagonal s over valid pairs, O/E = M_ij / expected(|i-j|)
    (entries with expected 0 get 1); 3) Pearson correlation between the O/E rows of valid bins; 4) leading
    eigenvector (numpy.linalg.eigh); 5) orientation: with ``reference`` (per-bin track such as GC content) the sign is
    flipped so the eigenvector correlates positively with it (A = positive); without it the sign is arbitrary
    and is fixed so that the entry of largest magnitude is positive. A = eigenvector >= 0, B = < 0.
    """

    m, valid = masked_nan(np, matrix)
    n = len(m)
    expected = np.array([np.nanmean(np.diagonal(m, offset=s)) if np.any(~np.isnan(np.diagonal(m, offset=s))) else 1.0 for s in range(n)])
    idx = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    exp_matrix = expected[idx]
    exp_matrix = np.where(exp_matrix > 0, exp_matrix, np.nan)
    oe = m / exp_matrix
    oe = np.where(np.isnan(oe) & ~np.isnan(m), 1.0, oe)
    sub = oe[np.ix_(valid, valid)]
    if sub.shape[0] < 3:
        raise InputFormatError("Need at least 3 non-empty bins for compartments.", code="GENOME3D_MATRIX_INVALID")
    if np.isnan(sub).any():
        sub = np.where(np.isnan(sub), 1.0, sub)
    corr = np.corrcoef(sub)
    corr = np.nan_to_num(corr, nan=0.0)
    values, vectors = np.linalg.eigh(corr)
    leading = vectors[:, -1]
    if reference is not None:
        ref = np.array(reference, dtype=float)[valid]
        if ref.size != leading.size or np.std(ref) == 0 or np.std(leading) == 0:
            raise InputFormatError("reference must have one value per bin and non-zero variance.", code="GENOME3D_INPUT_INVALID")
        if np.corrcoef(leading, ref)[0, 1] < 0:
            leading = -leading
    elif leading[np.argmax(np.abs(leading))] < 0:
        leading = -leading
    full: list[float | None] = [None] * n
    for slot, value in zip(np.where(valid)[0], leading):
        full[int(slot)] = float(value)
    return {"eigenvector": full, "labels": [None if v is None else ("A" if v >= 0 else "B") for v in full],
            "eigenvalue": float(values[-1]), "explained_variance": float(values[-1] / values.sum()) if values.sum() else None,
            "oriented_by_reference": reference is not None}


# ----------------------------------------------------------------------------- comparison
def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def pearson(np: Any, x: Any, y: Any) -> float | None:
    if x.size < 2 or np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


# ----------------------------------------------------------------------------- SVG
def svg_heatmap(matrix: Sequence[Sequence[float]], max_cells: int = 150, cell: int = 5) -> str:
    n = len(matrix)
    factor = -(-n // max_cells)
    grid = coarsen(matrix, factor) if factor > 1 else [list(r) for r in matrix]
    m = len(grid)
    top = max(math.log10(1 + x) for r in grid for x in r) or 1.0
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{m * cell + 20}" height="{m * cell + 20}" data-bins="{n}" data-factor="{factor}">']
    for i in range(m):
        for j in range(m):
            level = math.log10(1 + grid[i][j]) / top
            g = int(255 * (1 - level))
            parts.append(f'<rect class="cell" data-i="{i}" data-j="{j}" data-value="{fmt_value(grid[i][j])}" x="{10 + j * cell}" '
                         f'y="{10 + i * cell}" width="{cell}" height="{cell}" fill="rgb(255,{g},{g})"/>')
    parts.append("</svg>")
    return "".join(parts)
