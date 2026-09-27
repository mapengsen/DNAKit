"""Independent, dense diploid dosage PCA with explicit preprocessing rules.

NumPy is imported lazily through ``dnakit[bio-analysis]``. This is a PCA
baseline for an already prepared dosage matrix, not a PLINK replacement.

This version adds mean imputation, invariant-marker filtering, SVD-based PCA,
JSON-safe provenance metadata and configurable progress reporting.
"""

from __future__ import annotations

from numbers import Real
from typing import Any

from dnakit._progress import Progress
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

__all__ = ["genotype_pca"]


def _dosage_matrix(np: Any, genotypes: Any) -> Any:
    """Reject coercions of labels/flags while accepting JSON null and numeric NaN."""
    try:
        if not isinstance(genotypes, np.ndarray) or genotypes.dtype.kind not in "iuf":
            objects = np.asarray(genotypes, dtype=object)
            if any(
                value is not None
                and (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real))
                for value in objects.flat
            ):
                raise ValueError("only real numbers or null are allowed")
        data = np.asarray(genotypes, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(
            "genotypes must be a rectangular numeric matrix, without text or booleans."
        ) from exc
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] < 1:
        raise ConfigurationError("genotypes needs at least two samples and one marker column.")
    if np.isinf(data).any() or ((data < 0) | (data > 2)).any():
        raise ConfigurationError("genotypes must contain dosages in [0, 2], or None/NaN.")
    return data


def _sample_labels(sample_ids: Any, n_samples: int) -> list[str]:
    if sample_ids is None:
        return [f"sample_{index + 1}" for index in range(n_samples)]
    if (
        not isinstance(sample_ids, list)
        or len(sample_ids) != n_samples
        or any(not isinstance(label, str) or not label.strip() for label in sample_ids)
        or len(set(sample_ids)) != n_samples
    ):
        raise ConfigurationError("sample_ids must be a list of unique, nonempty sample names.")
    return list(sample_ids)


def genotype_pca(
    genotypes: list[list[float | None]],
    *,
    components: int = 2,
    sample_ids: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute standard PCA on frequency-standardized diploid allele dosages.

    Rows are samples and columns are aligned biallelic markers. Dosages may be
    fractional and must be in [0, 2]. None/NaN means missing; strings, booleans,
    complex values and infinity are rejected. Numeric NumPy matrices are also
    accepted. At least two samples and one variable marker are required.

    For each marker, estimate mean dosage mu using observed samples and set
    p = mu/2. Drop entirely missing or observation-invariant columns (including
    all-heterozygous columns); fill remaining missing entries with mu. Calculate
    X = (dosages - mu) / sqrt(2*p*(1-p)), then a thin SVD X = U*S*V.T.
    Scores are U*S, marker loadings are V, and covariance eigenvalues are
    S**2/(n_samples-1). Variance ratios use the full spectrum as denominator,
    even when only a subset of components is requested.

    The largest absolute loading in each component is assigned a positive
    sign (ties select the first retained marker). Equal eigenvalues may still
    rotate their non-unique subspace across numerical libraries. Components
    must be in 1..min(n_samples-1, n_markers_used); collinearity can produce
    zero-variance components. ``numerical_rank`` reports the nonzero rank.

    Results are JSON-compatible. Scores have shape samples x components;
    loadings have shape retained-markers x components, in ``marker_indices``
    order. Indices are zero-based. Means and frequencies cover ALL input
    columns, with null for wholly missing markers. Imputation counts concern
    retained markers only. Neither dosages nor sample labels are modified.

    This dense, in-memory implementation does not read PLINK files, prune LD,
    filter related samples, or provide PLINK-compatible eigenvalue scaling.
    Its covariance eigenvalues differ from eigenvalues of X*X.T/n_markers.
    Large cohorts should use an appropriate specialized analysis backend.
    """
    if isinstance(components, bool) or not isinstance(components, int) or components < 1:
        raise ConfigurationError("components must be a positive integer.")
    if not isinstance(verbose, bool):
        raise ConfigurationError("verbose must be boolean.")
    np = require("numpy", "bio-analysis")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Genotype PCA", total=4)
        data = _dosage_matrix(np, genotypes)
        labels = _sample_labels(sample_ids, data.shape[0])
        missing = np.isnan(data)
        observed_counts = np.sum(~missing, axis=0)
        absent = observed_counts == 0
        means = np.divide(
            np.nansum(data, axis=0),
            observed_counts,
            out=np.zeros(data.shape[1]),
            where=~absent,
        )
        minima = np.min(np.where(missing, np.inf, data), axis=0)
        maxima = np.max(np.where(missing, -np.inf, data), axis=0)
        invariant = ~absent & (minima == maxima)
        keep = ~(absent | invariant)
        n_markers_used = int(np.sum(keep))
        if not n_markers_used:
            raise ConfigurationError("No variable markers remain after filtering.")
        maximum_components = min(data.shape[0] - 1, n_markers_used)
        if components > maximum_components:
            raise ConfigurationError(
                f"components must be <= {maximum_components} after marker filtering."
            )
        progress.advance(task)

        frequencies = means[keep] / 2
        scales = np.sqrt(2 * frequencies * (1 - frequencies))
        if np.any(scales <= 0) or not np.isfinite(scales).all():
            raise ConfigurationError("Variable marker frequencies exceed floating-point precision.")
        selected_missing = missing[:, keep]
        imputed = np.where(selected_missing, means[keep], data[:, keep])
        standardized = (imputed - means[keep]) / scales
        if not np.isfinite(standardized).all():
            raise ConfigurationError("Genotype standardization produced non-finite values.")
        progress.advance(task)

        try:
            left, singular, right = np.linalg.svd(standardized, full_matrices=False)
        except np.linalg.LinAlgError as exc:
            raise ConfigurationError("Genotype PCA SVD did not converge.") from exc
        full_variances = singular**2 / (data.shape[0] - 1)
        total_variance = float(np.sum(full_variances))
        if not np.isfinite(total_variance) or total_variance <= 0:
            raise ConfigurationError("Genotype PCA requires finite, positive total variance.")
        loadings = right[:components].T.copy()
        scores = left[:, :components] * singular[:components]
        pivot_rows = np.argmax(np.abs(loadings), axis=0)
        signs = np.where(loadings[pivot_rows, np.arange(components)] < 0, -1, 1)
        loadings *= signs
        scores *= signs
        if not np.isfinite(loadings).all() or not np.isfinite(scores).all():
            raise ConfigurationError("Genotype PCA produced non-finite scores or loadings.")
        tolerance = float(singular[0] * max(standardized.shape) * np.finfo(float).eps)
        progress.advance(task)

        result = {
            "method": "diploid_dosage_pca",
            "parameters": {
                "ploidy": 2,
                "components": components,
                "missing": "marker_mean",
                "standardization": "(dosage - 2*p) / sqrt(2*p*(1-p))",
                "eigenvalue_denominator": "n_samples - 1",
                "sign_convention": "largest_absolute_loading_positive",
            },
            "sample_ids": labels,
            "n_samples": int(data.shape[0]),
            "n_markers": int(data.shape[1]),
            "n_markers_used": n_markers_used,
            "numerical_rank": int(np.sum(singular > tolerance)),
            "scores": scores.tolist(),
            "loadings": loadings.tolist(),
            "eigenvalues": full_variances[:components].tolist(),
            "explained_variance_ratio": (full_variances[:components] / total_variance).tolist(),
            "total_variance": total_variance,
            "marker_indices": np.flatnonzero(keep).tolist(),
            "dropped_all_missing_indices": np.flatnonzero(absent).tolist(),
            "dropped_invariant_indices": np.flatnonzero(invariant).tolist(),
            "marker_means": [
                None if is_absent else float(mean)
                for mean, is_absent in zip(means, absent, strict=True)
            ],
            "allele_frequencies": [
                None if is_absent else float(mean / 2)
                for mean, is_absent in zip(means, absent, strict=True)
            ],
            "imputed_value_count": int(np.sum(selected_missing)),
            "imputed_counts_per_marker": np.sum(selected_missing, axis=0).tolist(),
        }
        progress.advance(task)
    return result
