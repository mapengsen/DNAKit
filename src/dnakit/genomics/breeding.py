"""Diploid genomic relationships and ridge predictions for breeding research.

These are independent implementations of the published mathematical methods,
not copied rrBLUP code and not a reproduction of the BOLE agent framework.
NumPy is loaded only when a function is called (``dnakit[bio-analysis]``).
"""

from __future__ import annotations

import math
from numbers import Real
from typing import Any, cast

from dnakit._progress import Progress
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

__all__ = ["genomic_prediction", "genomic_relationship_matrix"]

_VANRADEN = "https://doi.org/10.3168/jds.2007-0980"
_ENDELMAN = "https://doi.org/10.3835/plantgenome2011.08.0024"
_RELATIONSHIP = "https://doi.org/10.1534/g3.112.004259"


def _settings(ploidy: int, missing: str, verbose: bool) -> None:
    if isinstance(ploidy, bool) or not isinstance(ploidy, int) or ploidy != 2:
        raise ConfigurationError("Only diploid dosages are supported; ploidy must be 2.")
    if missing not in ("mean", "error"):
        raise ConfigurationError("missing must be 'mean' or 'error'.")
    if not isinstance(verbose, bool):
        raise ConfigurationError("verbose must be boolean.")


def _numeric_array(np: Any, values: Any, name: str) -> Any:
    """Preserve JSON null as NaN without coercing text or flags into measurements."""
    try:
        # Keep the common NumPy numeric-matrix path vectorized. For JSON lists,
        # inspect original objects before coercion can turn True into 1 or text
        # such as "NaN" into an apparently legitimate missing measurement.
        if not isinstance(values, np.ndarray) or values.dtype.kind not in "iuf":
            objects = np.asarray(values, dtype=object)
            if any(
                value is not None
                and (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real))
                for value in objects.flat
            ):
                raise ValueError("only real numeric values or null are allowed")
        return np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(
            f"{name} must contain real numeric values or null, not text or booleans."
        ) from exc


def _dosages(np: Any, values: Any, name: str, *, min_rows: int, missing: str) -> Any:
    data = _numeric_array(np, values, name)
    if data.ndim != 2 or data.shape[0] < min_rows or data.shape[1] < 1:
        raise ConfigurationError(
            f"{name} must have at least {min_rows} rows and one marker column."
        )
    if np.isinf(data).any() or ((data < 0) | (data > 2)).any():
        raise ConfigurationError(f"{name} dosages must be finite values in [0, 2], or NaN.")
    if missing == "error" and np.isnan(data).any():
        raise ConfigurationError(f"{name} contains missing dosages with missing='error'.")
    return data


def _reference(np: Any, data: Any) -> tuple[Any, Any, Any, Any]:
    """Compute allele frequencies solely from the reference/training samples."""
    count = np.sum(~np.isnan(data), axis=0)
    means = np.divide(np.nansum(data, axis=0), count, out=np.zeros(data.shape[1]), where=count > 0)
    absent = count == 0
    monomorphic = ~absent & ((means == 0) | (means == 2))
    keep = ~(absent | monomorphic)
    if not np.any(keep):
        raise ConfigurationError("No observed polymorphic markers remain after filtering.")
    return means, keep, absent, monomorphic


def _center(np: Any, data: Any, means: Any, keep: Any) -> Any:
    selected = data[:, keep]
    return np.where(np.isnan(selected), means[keep], selected) - means[keep]


def _metadata(np: Any, means: Any, keep: Any, absent: Any, monomorphic: Any) -> dict[str, Any]:
    return {
        "allele_frequencies": [
            None if is_absent else float(mean / 2)
            for mean, is_absent in zip(means, absent, strict=True)
        ],
        "marker_indices": np.flatnonzero(keep).tolist(),
        "dropped_all_missing_indices": np.flatnonzero(absent).tolist(),
        "dropped_monomorphic_indices": np.flatnonzero(monomorphic).tolist(),
        "n_markers_used": int(np.sum(keep)),
    }


def genomic_relationship_matrix(
    genotypes: Any,
    *,
    ploidy: int = 2,
    missing: str = "mean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate the VanRaden method-1 additive genomic relationship matrix.

    Rows are samples (at least two), columns are biallelic markers, and entries
    are alternate-allele dosages in [0, 2]; fractional dosages are accepted.
    NaN/None (JSON null) is mean-imputed per marker with ``missing='mean'`` or rejected with
    ``missing='error'``. Entirely missing and allele-monomorphic (frequency 0 or
    1) columns are excluded and reported by zero-based indices. An all-
    heterozygous marker is allele-polymorphic and retained, as in rrBLUP A.mat.
    Text and booleans are not dosage measurements and are rejected.

    With frequencies estimated from these samples, Z = M - 2p and
    G = Z Z' / (2 sum(p (1-p))). This is an additive relationship matrix;
    the coancestry/kinship matrix under the usual convention is G/2. It is
    relative to this sample's frequencies, not an inferred base population.
    No shrinkage, pedigree blending, or polyploid generalization is applied.

    Returns a JSON-compatible dictionary including ``relationship_matrix``,
    full-width frequencies (None for entirely missing markers), and filtering
    metadata. Uses O(samples**2) memory for the resulting dense matrix.

    References: P. M. VanRaden (2008), doi:10.3168/jds.2007-0980;
    J. B. Endelman and J.-L. Jannink (2012), doi:10.1534/g3.112.004259.
    """
    _settings(ploidy, missing, verbose)
    np = require("numpy", "bio-analysis")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Genomic relationship matrix", total=3)
        data = _dosages(np, genotypes, "genotypes", min_rows=2, missing=missing)
        means, keep, absent, monomorphic = _reference(np, data)
        progress.advance(task)
        centered = _center(np, data, means, keep)
        frequencies = means[keep] / 2
        denominator = float(2 * np.sum(frequencies * (1 - frequencies)))
        relationship = (centered @ centered.T) / denominator
        if not np.isfinite(relationship).all():
            raise ConfigurationError("Relationship calculation produced non-finite values.")
        progress.advance(task)
        result = {
            "method": "vanraden_method_1",
            "ploidy": 2,
            "relationship_matrix": relationship.tolist(),
            "n_samples": int(data.shape[0]),
            "denominator": denominator,
            "imputed_value_count": int(np.isnan(data[:, keep]).sum()),
            **_metadata(np, means, keep, absent, monomorphic),
            "citations": [_VANRADEN, _RELATIONSHIP],
        }
        progress.advance(task)
    return result


def genomic_prediction(
    train_genotypes: Any,
    phenotypes: Any,
    predict_genotypes: Any,
    *,
    ridge_lambda: float,
    ploidy: int = 2,
    missing: str = "mean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit an intercept and fixed-lambda RR-BLUP marker effects, then predict.

    Both genotype matrices use sample rows and the SAME marker/allele column
    order, with diploid alternate-allele dosages in [0, 2]. Arrays carry no
    marker IDs: callers must harmonize both marker identity and effect alleles.
    Training requires at least two samples and one finite phenotype per row.
    Missing phenotypes are rejected; filter them and their genotype rows first.

    The objective is sum((y - mu - Z b)**2) + ridge_lambda * sum(b**2),
    with unpenalized mu and Z centered using only training-marker means.
    ``ridge_lambda`` must be finite and positive. It is the residual/marker
    effect variance ratio for this dosage coding, not an estimated heritability.
    Unlike R rrBLUP ``mixed.solve``, this function does not estimate lambda by
    REML. Choose it with training-only/nested validation outside this function.

    Missing dosages are replaced by training means (or rejected), and entirely
    missing or allele-monomorphic training columns are excluded. Prediction
    samples never affect filtering, frequencies, centering, or model fitting.
    A thin SVD solves the ridge objective without forming a marker-square matrix.

    Returns JSON-compatible ``predictions``, ``training_predictions``, full-width
    ``marker_effects`` (zero for excluded markers), and filtering metadata.
    ``intercept`` is the phenotype mean for CENTERED dosages; ``raw_intercept``
    applies to uncentered, imputed dosages. ``breeding_values`` are predictions
    minus the training phenotype mean. This is an additive research baseline,
    without environmental covariates, genotype-by-environment effects, or a
    prediction-accuracy claim.

    Reference: J. B. Endelman (2011), doi:10.3835/plantgenome2011.08.0024.
    """
    _settings(ploidy, missing, verbose)
    if isinstance(ridge_lambda, bool) or not isinstance(cast(object, ridge_lambda), Real):
        raise ConfigurationError("ridge_lambda must be finite and positive.")
    try:
        penalty = float(ridge_lambda)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError("ridge_lambda must be finite and positive.") from exc
    if not math.isfinite(penalty) or penalty <= 0:
        raise ConfigurationError("ridge_lambda must be finite and positive.")
    np = require("numpy", "bio-analysis")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Genomic ridge prediction", total=4)
        training = _dosages(np, train_genotypes, "train_genotypes", min_rows=2, missing=missing)
        prediction = _dosages(
            np, predict_genotypes, "predict_genotypes", min_rows=1, missing=missing
        )
        if prediction.shape[1] != training.shape[1]:
            raise ConfigurationError("Training and prediction matrices need the same columns.")
        response = _numeric_array(np, phenotypes, "phenotypes")
        if response.shape != (training.shape[0],) or not np.isfinite(response).all():
            raise ConfigurationError("One finite phenotype is required per training row.")
        means, keep, absent, monomorphic = _reference(np, training)
        progress.advance(task)
        centered = _center(np, training, means, keep)
        prediction_centered = _center(np, prediction, means, keep)
        intercept = float(response.mean())
        centered_response = response - intercept
        if not math.isfinite(intercept) or not np.isfinite(centered_response).all():
            raise ConfigurationError("Phenotype magnitudes overflow; rescale them first.")
        progress.advance(task)
        try:
            left, singular, right = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError as exc:
            raise ConfigurationError("Genomic prediction SVD did not converge.") from exc
        effects = right.T @ ((singular / (singular**2 + penalty)) * (left.T @ centered_response))
        fitted = intercept + centered @ effects
        breeding_values = prediction_centered @ effects
        predicted = intercept + breeding_values
        raw_intercept = float(intercept - means[keep] @ effects)
        if (
            not np.isfinite(effects).all()
            or not np.isfinite(fitted).all()
            or not np.isfinite(predicted).all()
            or not np.isfinite(breeding_values).all()
            or not math.isfinite(raw_intercept)
        ):
            raise ConfigurationError("Prediction produced non-finite values; rescale phenotypes.")
        progress.advance(task)
        full_effects = np.zeros(training.shape[1])
        full_effects[keep] = effects
        result = {
            "method": "rrblup_fixed_lambda",
            "ridge_lambda": penalty,
            "ploidy": 2,
            "intercept": intercept,
            "raw_intercept": raw_intercept,
            "marker_effects": full_effects.tolist(),
            "training_predictions": fitted.tolist(),
            "predictions": predicted.tolist(),
            "breeding_values": breeding_values.tolist(),
            "n_training_samples": int(training.shape[0]),
            "n_prediction_samples": int(prediction.shape[0]),
            "imputed_training_value_count": int(np.isnan(training[:, keep]).sum()),
            "imputed_prediction_value_count": int(np.isnan(prediction[:, keep]).sum()),
            **_metadata(np, means, keep, absent, monomorphic),
            "citations": [_ENDELMAN],
        }
        progress.advance(task)
    return result
