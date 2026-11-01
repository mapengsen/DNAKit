"""PLINK text genotypes and a single-marker linear-mixed-model (LMM) GWAS scan.

``read_plink_ped`` turns a PLINK 1 ``.ped``/``.map`` text pair into diploid dosages.
``mixed_model_gwas`` fits the null model ``y = mu + u + e`` with ``u ~ N(0, Vu K)`` by REML
once and then tests each marker with variance components held fixed (P3D/EMMAX). The REML
fit and the per-marker F test follow R rrBLUP ``mixed.solve``/``GWAS(P3D=TRUE)``; the
``coding='bglr'`` reader mode reproduces R BGLR ``read_ped``. These are independent
implementations of the published methods, checked against the R packages, not copied code.
NumPy and SciPy are imported lazily (``dnakit[bio-analysis]``).
"""

from __future__ import annotations

import math
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit._progress import Progress
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

__all__ = ["mixed_model_gwas", "read_plink_ped"]

_RRBLUP = "https://doi.org/10.3835/plantgenome2011.08.0024"
_P3D = "https://doi.org/10.1038/ng.546"
_EMMAX = "https://doi.org/10.1038/ng.548"
_PED_CODINGS = ("minor", "bglr")
_MISSING_MODES = ("exact", "rrblup", "mean")
_MAX_CELLS = 200_000_000


# --------------------------------------------------------------------------- PED reader
def _ped_paths(path: object, map_path: object) -> tuple[Path, Path]:
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise ConfigurationError("path must name a PLINK .ped file.")
    ped = Path(path).expanduser()
    if not ped.is_file():
        raise ConfigurationError(f"PLINK .ped file not found: {ped}")
    if map_path is None:
        mapped = ped.with_suffix(".map")
    elif isinstance(map_path, (str, Path)) and str(map_path).strip():
        mapped = Path(map_path).expanduser()
    else:
        raise ConfigurationError("map_path must be a path or None.")
    if not mapped.is_file():
        raise ConfigurationError(f"PLINK .map file not found: {mapped}")
    return ped, mapped


def _read_map(path: Path) -> tuple[list[str], list[str], list[int]]:
    chromosomes: list[str] = []
    variant_ids: list[str] = []
    positions: list[int] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) not in (3, 4):
                raise ConfigurationError(f"{path.name} line {number}: expected 3 or 4 columns.")
            try:
                position = int(fields[-1])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number}: base-pair position must be an integer."
                ) from exc
            chromosomes.append(fields[0])
            variant_ids.append(fields[1])
            positions.append(position)
    if not variant_ids:
        raise ConfigurationError(f"{path.name} lists no variants.")
    return chromosomes, variant_ids, positions


def _read_ped(np: Any, path: Path, n_variants: int, max_cells: int) -> tuple[list[list[str]], Any]:
    meta: list[list[str]] = []
    rows: list[list[str]] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 6 + 2 * n_variants:
                raise ConfigurationError(
                    f"{path.name} line {number}: expected 6 + 2 x {n_variants} columns "
                    f"(the .map lists {n_variants} variants), found {len(fields)}."
                )
            if (len(rows) + 1) * n_variants > max_cells:
                raise ConfigurationError("The genotype matrix exceeds max_cells.")
            meta.append(fields[:6])
            rows.append(fields[6:])
    if len(rows) < 2:
        raise ConfigurationError(f"{path.name} needs at least two samples.")
    alleles = np.asarray(rows, dtype=str).reshape(len(rows), n_variants, 2)
    return meta, alleles


def _bglr_dosages(np: Any, alleles: Any) -> tuple[Any, list[str | None], list[str | None]]:
    """BGLR read_ped: genotypes are coded against the first character of the FIRST sample's
    first allele; '0 0' is missing and any two different characters are heterozygous."""
    first = np.vectorize(lambda token: token[:1], otypes=["U1"])(alleles)
    one, two = first[:, :, 0], first[:, :, 1]
    reference = one[0, :]
    missing = (one == "0") & (two == "0")
    heterozygous = (one != two) & ~missing
    hom_reference = (one == two) & (one == reference)
    dosage = np.where(missing, np.nan, np.where(heterozygous, 1.0, np.where(hom_reference, 0.0, 2.0)))
    counted: list[str | None] = []
    for j in range(first.shape[1]):
        seen = [a for a in np.unique(first[:, j, :]) if a not in ("0", reference[j])]
        counted.append(seen[0] if len(seen) == 1 else None)
    return dosage, [None if r == "0" else str(r) for r in reference], counted


def _minor_dosages(
    np: Any, alleles: Any, missing_allele: str
) -> tuple[Any, list[str | None], list[str | None], list[int], int]:
    n, m, _ = alleles.shape
    missing = (alleles[:, :, 0] == missing_allele) | (alleles[:, :, 1] == missing_allele)
    half_missing = int(np.sum((alleles[:, :, 0] == missing_allele) ^ (alleles[:, :, 1] == missing_allele)))
    dosage = np.full((n, m), np.nan)
    counted: list[str | None] = []
    other: list[str | None] = []
    multiallelic: list[int] = []
    for j in range(m):
        observed = alleles[~missing[:, j], j, :]
        names, counts = np.unique(observed, return_counts=True)
        if len(names) > 2:
            multiallelic.append(j)
            counted.append(None)
            other.append(None)
            continue
        if len(names) == 0:
            counted.append(None)
            other.append(None)
            continue
        if len(names) == 1:
            counted.append(None)
            other.append(str(names[0]))
            dosage[~missing[:, j], j] = 0.0
            continue
        # minor = the less frequent allele; on a tie, the allele that sorts last
        minor_index = 0 if counts[0] < counts[1] else 1
        minor = str(names[minor_index])
        counted.append(minor)
        other.append(str(names[1 - minor_index]))
        copies = (alleles[:, j, 0] == minor).astype(float) + (alleles[:, j, 1] == minor)
        dosage[:, j] = np.where(missing[:, j], np.nan, copies)
    return dosage, counted, other, multiallelic, half_missing


def read_plink_ped(
    path: str | Path,
    *,
    map_path: str | Path | None = None,
    coding: str = "minor",
    missing_allele: str = "0",
    include_genotypes: bool = True,
    max_cells: int = _MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a PLINK 1 text ``.ped`` + ``.map`` pair into a samples x variants dosage matrix.

    ``path`` is the .ped file (whitespace-delimited: family ID, individual ID, father,
    mother, sex, phenotype, then two allele columns per variant); ``map_path`` defaults to
    the same name with ``.map`` (chromosome, variant ID, [cM,] base-pair position). Rows
    follow the .ped order; columns follow the .map order. Dosages are 0/1/2 copies of the
    counted allele, NaN (JSON null) when missing.

    ``coding='minor'`` (default) counts the less frequent observed allele per variant (on
    a tie, the allele that sorts last); a genotype with either allele equal to
    ``missing_allele`` is missing; variants with more than two alleles are all-missing and
    listed in ``multiallelic_indices``; single-allele variants are all 0.

    ``coding='bglr'`` reproduces R ``BGLR::read_ped`` followed by the usual recoding
    (code 2 -> NA, code 3 -> 2) used in rrBLUP/SNPRelate tutorials: each variant is coded
    against the first character of the FIRST sample's first allele (0 = homozygous for it,
    1 = two different characters, 2 = any other homozygote, missing = '0 0'). When the first
    sample is missing at a variant, both homozygotes become 2, so that variant looks
    monomorphic and is removed by a later MAF filter; ``bglr_first_sample_missing_indices``
    lists these variants. Use this mode only to reproduce results of R pipelines built on
    BGLR::read_ped; ``counted_alleles`` then names the allele counted where it is unique.

    ``include_genotypes=False`` returns only the metadata and counts (useful for large
    files over MCP). ``max_cells`` caps samples x variants.

    Returns JSON-compatible ``sample_ids`` (individual IDs), ``family_ids``,
    ``ped_phenotypes`` (column 6 as text), ``variant_ids``, ``chromosomes``, ``positions``,
    ``counted_alleles``, ``other_alleles``, per-variant ``missing_counts``,
    ``n_missing_genotypes`` and ``genotypes`` (rows = samples).

    Reference: PLINK 1 PED/MAP formats (https://www.cog-genomics.org/plink/1.9/formats).
    """
    if coding not in _PED_CODINGS:
        raise ConfigurationError(f"coding must be one of {', '.join(_PED_CODINGS)}.")
    if not isinstance(missing_allele, str) or len(missing_allele) == 0:
        raise ConfigurationError("missing_allele must be a nonempty string.")
    if coding == "bglr" and missing_allele != "0":
        raise ConfigurationError("coding='bglr' always treats '0' as the missing allele.")
    for name, flag in (("include_genotypes", include_genotypes), ("verbose", verbose)):
        if not isinstance(flag, bool):
            raise ConfigurationError(f"{name} must be boolean.")
    if isinstance(max_cells, bool) or not isinstance(max_cells, int) or max_cells < 1:
        raise ConfigurationError("max_cells must be a positive integer.")
    np = require("numpy", "bio-analysis")
    ped, mapped = _ped_paths(path, map_path)
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Read PLINK .ped", total=3)
        chromosomes, variant_ids, positions = _read_map(mapped)
        meta, alleles = _read_ped(np, ped, len(variant_ids), max_cells)
        progress.advance(task)
        multiallelic: list[int] = []
        half_missing = 0
        first_missing: list[int] = []
        if coding == "bglr":
            dosage, other, counted = _bglr_dosages(np, alleles)
            first_missing = [j for j, allele in enumerate(other) if allele is None]
        else:
            dosage, counted, other, multiallelic, half_missing = _minor_dosages(
                np, alleles, missing_allele
            )
        progress.advance(task)
        missing_counts = np.sum(np.isnan(dosage), axis=0)
        result: dict[str, Any] = {
            "method": "plink_ped_text",
            "coding": coding,
            "ped_path": str(ped),
            "map_path": str(mapped),
            "n_samples": int(dosage.shape[0]),
            "n_variants": int(dosage.shape[1]),
            "sample_ids": [row[1] for row in meta],
            "family_ids": [row[0] for row in meta],
            "ped_phenotypes": [row[5] for row in meta],
            "variant_ids": variant_ids,
            "chromosomes": chromosomes,
            "positions": positions,
            "counted_alleles": counted,
            "other_alleles": other,
            "missing_counts": missing_counts.astype(int).tolist(),
            "n_missing_genotypes": int(missing_counts.sum()),
            "multiallelic_indices": multiallelic,
            "half_missing_genotype_count": half_missing,
        }
        if coding == "bglr":
            result["bglr_first_sample_missing_indices"] = first_missing
        if include_genotypes:
            result["genotypes"] = [
                [None if math.isnan(value) else value for value in row] for row in dosage.tolist()
            ]
        progress.advance(task)
    return result


# --------------------------------------------------------------------------- LMM scan
def _matrix(np: Any, values: Any, name: str, *, allow_nan: bool) -> Any:
    try:
        if not isinstance(values, np.ndarray) or values.dtype.kind not in "iuf":
            objects = np.asarray(values, dtype=object)
            if any(
                value is not None
                and (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real))
                for value in objects.flat
            ):
                raise ValueError("only real numbers or null are allowed")
        data = np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"{name} must be numeric (null allowed for missing).") from exc
    if np.isinf(data).any() or (not allow_nan and np.isnan(data).any()):
        raise ConfigurationError(f"{name} must be finite" + (" or null." if allow_nan else "."))
    return data


def _brent_fmin(function: Any, lower: float, upper: float, tol: float) -> float:
    """Brent's one-dimensional minimizer, step for step as R ``optimize`` (Brent_fmin)."""
    c = (3.0 - math.sqrt(5.0)) * 0.5
    eps = math.sqrt(2.220446049250313e-16)
    a, b = lower, upper
    v = a + c * (b - a)
    w = x = v
    d = e = 0.0
    fx = function(x)
    fv = fw = fx
    tol3 = tol / 3.0
    while True:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2.0
        if abs(x - xm) <= t2 - (b - a) * 0.5:
            break
        p = q = r = 0.0
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            else:
                q = -q
            r = e
            e = d
        if abs(p) >= abs(q * 0.5 * r) or p <= q * (a - x) or p >= q * (b - x):
            e = (b - x) if x < xm else (a - x)
            d = c * e
        else:
            d = p / q
            u = x + d
            if u - a < t2 or b - u < t2:
                d = tol1 if x < xm else -tol1
        if abs(d) >= tol1:
            u = x + d
        elif d > 0.0:
            u = x + tol1
        else:
            u = x - tol1
        fu = function(u)
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x


def _reml_null(np: Any, y: Any, kinship: Any, bounds: tuple[float, float]) -> dict[str, Any]:
    """Intercept-only REML as rrBLUP mixed.solve(method='REML') with Z = I (eigen route)."""
    n = len(y)
    offset = math.sqrt(n)
    hb = kinship + offset * np.eye(n)
    values, vectors = np.linalg.eigh(hb)
    phi = values - offset
    if phi.min() < -1e-6:
        raise ConfigurationError("kinship is not positive semi-definite (eigenvalue < -1e-6).")
    s = np.eye(n) - np.full((n, n), 1.0 / n)
    shbs = s @ hb @ s
    s_values, s_vectors = np.linalg.eigh((shbs + shbs.T) / 2)
    order = np.argsort(s_values)[::-1][: n - 1]
    theta = s_values[order] - offset
    omega_sq = (s_vectors[:, order].T @ y) ** 2

    def objective(lam: float) -> float:
        return float((n - 1) * math.log(np.sum(omega_sq / (theta + lam))) + np.sum(np.log(theta + lam)))

    lam = _brent_fmin(objective, bounds[0], bounds[1], 2.220446049250313e-16**0.25)
    vu = float(np.sum(omega_sq / (theta + lam)) / (n - 1))
    hinv = (vectors / (phi + lam)) @ vectors.T
    return {"lambda": lam, "vu": vu, "ve": lam * vu, "objective": objective(lam), "hinv": hinv}


def _f_tests(np: Any, special: Any, hi_1: Any, hi_x: Any, hi_y: Any, x: Any, y: Any, n: int) -> tuple[Any, ...]:
    """Vectorised rrBLUP score.calc: GLS on [1, x], residual variance re-estimated, F(1, n-2)."""
    s11 = hi_1.sum(axis=0)
    s1x = hi_x.sum(axis=0)
    sxx = np.einsum("ij,ij->j", x, hi_x)
    s1y = hi_y.sum(axis=0)
    sxy = np.einsum("ij,ij->j", x, hi_y)
    syy = np.einsum("ij,ij->j", y, hi_y)
    det = s11 * sxx - s1x**2
    with np.errstate(divide="ignore", invalid="ignore"):
        b0 = (sxx * s1y - s1x * sxy) / det
        b1 = (s11 * sxy - s1x * s1y) / det
        v2 = n - 2
        rss = syy - 2 * (b0 * s1y + b1 * sxy) + (b0**2 * s11 + 2 * b0 * b1 * s1x + b1**2 * sxx)
        s2 = rss / v2
        winv11 = s11 / det
        f = b1**2 / (s2 * winv11)
        p = special.betainc(v2 / 2.0, 0.5, v2 / (v2 + f))
    return b1, np.sqrt(s2 * winv11), f, p, det


def mixed_model_gwas(
    genotypes: Any,
    phenotypes: Any,
    kinship: Any,
    *,
    missing_genotypes: str = "exact",
    variant_ids: list[str] | None = None,
    chromosomes: list[Any] | None = None,
    positions: list[int] | None = None,
    min_maf: float = 0.0,
    alpha: float = 0.05,
    reml_bounds: tuple[float, float] = (1e-9, 1e9),
    verbose: bool = True,
) -> dict[str, Any]:
    """Single-marker LMM GWAS with a kinship random effect (P3D / EMMAX, rrBLUP-compatible).

    ``genotypes``: samples x markers dosages in [0, 2] (NaN/null = missing), rows aligned
    with ``phenotypes`` (NaN/null = missing) and with the rows/columns of the symmetric,
    positive semi-definite ``kinship`` matrix (e.g. a GRM from genomic_relationship_matrix;
    set negative eigenvalues to zero first). Samples with a missing phenotype are dropped.

    Null model ``y = mu + u + e``, ``u ~ N(0, Vu K)``, ``e ~ N(0, Ve I)`` fitted once by
    REML, exactly as R rrBLUP ``mixed.solve`` (spectral form, ``lambda = Ve/Vu`` searched
    on ``reml_bounds`` with R's Brent optimizer and tolerance). With ``H = K + lambda I``
    held fixed, each marker x is tested by GLS of y on [1, x] with the residual variance
    re-estimated, an F(1, n-2) test as in rrBLUP ``GWAS(P3D=TRUE)``; P values are
    identical to that function's 10^-score. Effect sizes are per counted allele.

    ``missing_genotypes`` handles samples missing a marker:
      'exact'  - drop them for that marker; use the exact inverse of H on the remaining
                 samples (default);
      'rrblup' - drop them but take the sub-matrix of the full H^-1, exactly what rrBLUP
                 ``GWAS`` does (an approximation; use it to reproduce rrBLUP results);
      'mean'   - fill them with the marker mean (as GEMMA/EMMAX-style pipelines).
    Markers with minor-allele frequency below ``min_maf`` (among phenotyped samples with
    data) or without variation are reported with status and no P value. The Bonferroni
    threshold is ``alpha / n_tested``; ``lambda_gc`` is the median 1-df chi-square of the
    tested P values divided by 0.4549.

    Returns JSON-compatible variance components, ``records`` (one per marker, input
    order: variant_id, chromosome, position, n, maf, beta, se, f_statistic, p,
    minus_log10_p, status), ``n_tested``, ``bonferroni_threshold``, ``significant_indices``
    (sorted by P) and ``lambda_gc``.

    References: Endelman (2011) doi:10.3835/plantgenome2011.08.0024; Zhang et al. (2010)
    doi:10.1038/ng.546; Kang et al. (2010) doi:10.1038/ng.548.
    """
    if missing_genotypes not in _MISSING_MODES:
        raise ConfigurationError(f"missing_genotypes must be one of {', '.join(_MISSING_MODES)}.")
    for name, value, top in (("min_maf", min_maf, 0.5), ("alpha", alpha, 1.0)):
        if isinstance(value, bool) or not isinstance(value, Real) or not 0 <= float(value) <= top:
            raise ConfigurationError(f"{name} must be a number in [0, {top}].")
    if (
        not isinstance(reml_bounds, (tuple, list))
        or len(reml_bounds) != 2
        or not all(isinstance(b, Real) and not isinstance(b, bool) for b in reml_bounds)
        or not 0 < float(reml_bounds[0]) < float(reml_bounds[1])
    ):
        raise ConfigurationError("reml_bounds must be two positive numbers, lower < upper.")
    if not isinstance(verbose, bool):
        raise ConfigurationError("verbose must be boolean.")
    np = require("numpy", "bio-analysis")
    special = require("scipy.special", "bio-analysis")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("LMM GWAS scan", total=3)
        g = _matrix(np, genotypes, "genotypes", allow_nan=True)
        y_all = _matrix(np, phenotypes, "phenotypes", allow_nan=True)
        k_all = _matrix(np, kinship, "kinship", allow_nan=False)
        if g.ndim != 2 or y_all.ndim != 1 or g.shape[0] != y_all.shape[0]:
            raise ConfigurationError("genotypes (samples x markers) and phenotypes must align.")
        if k_all.shape != (g.shape[0], g.shape[0]):
            raise ConfigurationError("kinship must be samples x samples, aligned with genotypes.")
        if not np.allclose(k_all, k_all.T, atol=1e-8):
            raise ConfigurationError("kinship must be symmetric.")
        if ((g < 0) | (g > 2)).any():
            raise ConfigurationError("genotype dosages must lie in [0, 2].")
        m = g.shape[1]
        for name, labels in (("variant_ids", variant_ids), ("chromosomes", chromosomes), ("positions", positions)):
            if labels is not None and (not isinstance(labels, list) or len(labels) != m):
                raise ConfigurationError(f"{name} must be a list with one entry per marker.")
        keep = ~np.isnan(y_all)
        n = int(keep.sum())
        if n < 5:
            raise ConfigurationError("At least five phenotyped samples are required.")
        y = y_all[keep]
        k = k_all[np.ix_(keep, keep)]
        x = g[keep] - 1.0  # rrBLUP coding {-1, 0, 1}; the F test is shift invariant
        null = _reml_null(np, y, k, (float(reml_bounds[0]), float(reml_bounds[1])))
        hinv = null["hinv"]
        progress.advance(task)

        observed = ~np.isnan(x)
        counts = observed.sum(axis=0)
        means = np.divide(np.nansum(x + 1.0, axis=0), counts, out=np.zeros(m), where=counts > 0)
        maf = np.minimum(means / 2, 1 - means / 2)
        variable = np.array([
            counts[j] > 0 and np.nanmax(x[:, j]) > np.nanmin(x[:, j]) for j in range(m)
        ])
        testable = variable & (maf >= float(min_maf)) & (counts >= 4)
        n_used = np.where(missing_genotypes == "mean", n, counts)
        beta = np.full(m, np.nan)
        se = np.full(m, np.nan)
        fstat = np.full(m, np.nan)
        pval = np.full(m, np.nan)
        ones = np.ones(n)
        hi_1_full = hinv @ ones
        hi_y_full = hinv @ y
        complete = testable & (counts == n)
        if missing_genotypes == "mean":
            complete = testable.copy()
            x = np.where(observed, x, (means - 1.0)[None, :])
        cols = np.flatnonzero(complete)
        for start in range(0, len(cols), 4096):
            block = cols[start : start + 4096]
            xb = x[:, block]
            hi_x = hinv @ xb
            b1, s, f, p, det = _f_tests(
                np, special, np.repeat(hi_1_full[:, None], len(block), 1), hi_x,
                np.repeat(hi_y_full[:, None], len(block), 1), xb, np.repeat(y[:, None], len(block), 1), n,
            )
            ok = (det > 0) & np.isfinite(p)
            beta[block], se[block], fstat[block], pval[block] = (np.where(ok, v, np.nan) for v in (b1, s, f, p))
        for j in np.flatnonzero(testable & ~complete):
            obs = observed[:, j]
            xo, yo = x[obs, j], y[obs]
            sub = hinv[np.ix_(obs, obs)]
            rhs = np.column_stack([np.ones(len(yo)), xo, yo])
            hi_rhs = sub @ rhs
            if missing_genotypes == "exact":
                mis = ~obs
                cross = hinv[np.ix_(obs, mis)]
                hi_rhs = hi_rhs - cross @ np.linalg.solve(hinv[np.ix_(mis, mis)], cross.T @ rhs)
            b1, s, f, p, det = _f_tests(
                np, special, hi_rhs[:, [0]], hi_rhs[:, [1]], hi_rhs[:, [2]],
                xo[:, None], yo[:, None], len(yo),
            )
            if det[0] > 0 and math.isfinite(p[0]):
                beta[j], se[j], fstat[j], pval[j] = b1[0], s[0], f[0], p[0]
        progress.advance(task)

        tested = np.flatnonzero(np.isfinite(pval))
        n_tested = len(tested)
        threshold = float(alpha) / n_tested if n_tested else None
        significant = [] if threshold is None else sorted(
            (int(j) for j in tested if pval[j] < threshold), key=lambda j: pval[j]
        )
        stats = require("scipy.stats", "bio-analysis")
        lambda_gc = (
            float(np.median(stats.chi2.isf(pval[tested], 1)) / stats.chi2.ppf(0.5, 1)) if n_tested else None
        )

        def status(j: int) -> str:
            if math.isfinite(pval[j]):
                return "ok"
            if not variable[j]:
                return "monomorphic_or_all_missing"
            if maf[j] < float(min_maf):
                return "maf_below_threshold"
            return "too_few_samples_or_singular"

        def number(value: float) -> float | None:
            return float(value) if math.isfinite(value) else None

        records = [
            {
                "variant_id": None if variant_ids is None else str(variant_ids[j]),
                "chromosome": None if chromosomes is None else str(chromosomes[j]),
                "position": None if positions is None else int(positions[j]),
                "n": int(n_used[j]),
                "maf": number(maf[j]),
                "beta": number(beta[j]),
                "se": number(se[j]),
                "f_statistic": number(fstat[j]),
                "p": number(pval[j]),
                "minus_log10_p": number(-math.log10(pval[j])) if pval[j] > 0 else None,
                "status": status(j),
            }
            for j in range(m)
        ]
        result = {
            "method": "lmm_p3d_f_test",
            "parameters": {
                "missing_genotypes": missing_genotypes,
                "min_maf": float(min_maf),
                "alpha": float(alpha),
                "reml_bounds": [float(reml_bounds[0]), float(reml_bounds[1])],
                "variance_components": "REML once on the null model (P3D)",
                "test": "GLS F(1, n-2) on the marker effect, residual variance re-estimated",
            },
            "n_samples": int(g.shape[0]),
            "n_phenotyped": n,
            "n_markers": int(m),
            "n_tested": n_tested,
            "lambda": float(null["lambda"]),
            "vu": null["vu"],
            "ve": null["ve"],
            "heritability": null["vu"] / (null["vu"] + null["ve"]),
            "reml_objective": null["objective"],
            "bonferroni_threshold": threshold,
            "significant_indices": significant,
            "lambda_gc": lambda_gc,
            "records": records,
            "citations": [_RRBLUP, _P3D, _EMMAX],
        }
        progress.advance(task)
    return result
