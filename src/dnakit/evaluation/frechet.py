"""Fréchet distance between two DNA representation distributions."""

from __future__ import annotations

import math
from typing import Any

from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.representations import RepresentationBackend, extract_representations

from ._shared import EvaluationInput, materialize_input, record_for, report, require_nonempty
from .config import FrechetDistanceConfig
from .results import EvaluationReport


def _require_numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:
        raise BackendUnavailableError(
            "Fréchet DNA distance requires NumPy and a representation backend.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
        ) from exc
    return np


def _l2_normalize(matrix: Any, np: Any) -> Any:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    zero_rows = np.flatnonzero(norms[:, 0] == 0)
    if len(zero_rows):
        raise ConfigurationError(
            "Cannot L2-normalize zero-valued model representations.",
            code="ZERO_REPRESENTATION_VECTOR",
            context={"row_indices": tuple(int(value) for value in zero_rows[:20])},
        )
    return matrix / norms


def _frechet_components(left: Any, right: Any, np: Any) -> dict[str, float]:
    """Compute Gaussian Fréchet terms through an exact sample-space identity.

    If ``A`` and ``B`` are centered sample matrices scaled by ``sqrt(n - 1)``,
    the covariance cross term equals the nuclear norm of ``A @ B.T``.  This
    avoids materializing dense feature-by-feature covariance matrices.
    """

    left_mean = np.mean(left, axis=0, dtype=np.float64)
    right_mean = np.mean(right, axis=0, dtype=np.float64)
    mean_delta = left_mean - right_mean
    mean_component = float(np.dot(mean_delta, mean_delta))

    left_scaled = (left - left_mean) / math.sqrt(int(left.shape[0]) - 1)
    right_scaled = (right - right_mean) / math.sqrt(int(right.shape[0]) - 1)
    left_covariance_trace = float(np.sum(left_scaled * left_scaled, dtype=np.float64))
    right_covariance_trace = float(np.sum(right_scaled * right_scaled, dtype=np.float64))
    try:
        cross_gram = left_scaled @ right_scaled.T
        singular_values = np.linalg.svd(cross_gram, compute_uv=False)
    except (ValueError, np.linalg.LinAlgError) as exc:
        raise BackendExecutionError(
            "Could not compute the Fréchet covariance cross term.",
            code="FRECHET_DECOMPOSITION_FAILED",
            context={"error_type": type(exc).__name__},
        ) from exc
    cross_nuclear_norm = float(np.sum(singular_values, dtype=np.float64))
    raw_covariance_component = (
        left_covariance_trace + right_covariance_trace - 2.0 * cross_nuclear_norm
    )
    tolerance = 1e-10 * max(
        1.0,
        left_covariance_trace + right_covariance_trace,
        cross_nuclear_norm,
    )
    if raw_covariance_component < -tolerance:
        raise BackendExecutionError(
            "The Fréchet covariance term is numerically inconsistent.",
            code="INVALID_FRECHET_COVARIANCE_TERM",
            context={
                "covariance_component": raw_covariance_component,
                "tolerance": tolerance,
            },
        )
    covariance_component = max(0.0, raw_covariance_component)
    distance = mean_component + covariance_component
    if not math.isfinite(distance):
        raise BackendExecutionError(
            "The Fréchet distance is not finite.",
            code="INVALID_FRECHET_DISTANCE",
        )
    return {
        "frechet_distance": distance,
        "mean_component": mean_component,
        "covariance_component": covariance_component,
        "left_covariance_trace": left_covariance_trace,
        "right_covariance_trace": right_covariance_trace,
        "cross_covariance_nuclear_norm": cross_nuclear_norm,
    }


def evaluate_frechet_distance(
    left: EvaluationInput,
    right: EvaluationInput,
    *,
    config: FrechetDistanceConfig | None = None,
    backend: RepresentationBackend | None = None,
) -> EvaluationReport:
    """Compare two DNA collections with a Gaussian Fréchet representation distance.

    The default representation model is LucaOne, matching ``DATA-027``.  A
    standard LucaOne backend still requires explicit ``allow_remote_code=True``
    after the caller reviews its checkpoint code.  Passing ``backend`` enables
    controlled custom encoders and tests without a checkpoint download.
    """

    resolved = FrechetDistanceConfig() if config is None else config
    if not isinstance(resolved, FrechetDistanceConfig):
        raise ConfigurationError(
            "config must be FrechetDistanceConfig or None.",
            code="INVALID_FRECHET_DISTANCE_CONFIG",
        )
    left_items = materialize_input(left, limits=resolved.limits)
    right_items = materialize_input(right, limits=resolved.limits)
    require_nonempty(left_items, "Fréchet DNA distance")
    require_nonempty(right_items, "Fréchet DNA distance")
    if len(left_items) < 2 or len(right_items) < 2:
        raise ConfigurationError(
            "Fréchet DNA distance requires at least two records in each collection.",
            code="FRECHET_MINIMUM_SAMPLE_SIZE",
            context={"left_count": len(left_items), "right_count": len(right_items)},
        )
    total_symbols = sum(
        item.sequence.symbol_length for item in left_items + right_items
    )
    if total_symbols > resolved.limits.max_total_symbols:
        raise ConfigurationError(
            "Combined Fréchet inputs exceed max_total_symbols.",
            code="EVALUATION_SYMBOL_LIMIT",
            context={
                "total_symbols": total_symbols,
                "max_total_symbols": resolved.limits.max_total_symbols,
            },
        )
    cross_gram_elements = len(left_items) * len(right_items)
    if cross_gram_elements > resolved.max_cross_gram_elements:
        raise ConfigurationError(
            "Fréchet DNA distance exceeds max_cross_gram_elements.",
            code="FRECHET_CROSS_GRAM_LIMIT",
            context={
                "cross_gram_elements": cross_gram_elements,
                "max_cross_gram_elements": resolved.max_cross_gram_elements,
            },
        )

    records = tuple(record_for(item) for item in left_items + right_items)
    representation_result = extract_representations(
        records,
        config=resolved.representation,
        backend=backend,
    )
    np = _require_numpy()
    matrix = np.asarray(representation_result.representations, dtype=np.float64)
    if resolved.normalize:
        matrix = _l2_normalize(matrix, np)
    left_matrix = matrix[: len(left_items)]
    right_matrix = matrix[len(left_items) :]
    components = _frechet_components(left_matrix, right_matrix, np)

    return report(
        name="frechet_dna_distance",
        method="gaussian-frechet-distance-over-dna-model-representations",
        version="eval-frechet-dna-v1",
        parameters={
            "formula": (
                "||mu_left-mu_right||_2^2 + trace(Sigma_left + Sigma_right "
                "- 2*(Sigma_left^(1/2)*Sigma_right*Sigma_left^(1/2))^(1/2))"
            ),
            "encoder": representation_result.model_name,
            "checkpoint_path": representation_result.checkpoint_path,
            "pooling": representation_result.pooling,
            "representation_config": {
                "model": representation_result.model_name,
                "checkpoint_path": representation_result.checkpoint_path,
                "model_source_path": resolved.representation.model_source_path,
                "pooling": representation_result.pooling,
                "ambiguity_policy": representation_result.ambiguity_policy,
                "device": resolved.representation.device,
                "dtype": resolved.representation.dtype,
                "batch_size": resolved.representation.batch_size,
                "max_length": resolved.representation.max_length,
                "allow_remote_code": resolved.representation.allow_remote_code,
            },
            "normalize": resolved.normalize,
            "sample_covariance_ddof": 1,
            "numerical_method": "exact sample-space cross-Gram nuclear norm",
            "lower_is_better": True,
            "max_cross_gram_elements": resolved.max_cross_gram_elements,
            "limits": resolved.limits,
            "inference": (
                "FCD-analogous descriptive distance; not ChemNet FCD and not an "
                "experimentally validated DNA quality score"
            ),
        },
        metrics={
            **components,
            "left_count": len(left_items),
            "right_count": len(right_items),
            "embedding_dimension": representation_result.embedding_dimension,
            "left_covariance_rank_upper_bound": min(
                len(left_items) - 1,
                representation_result.embedding_dimension,
            ),
            "right_covariance_rank_upper_bound": min(
                len(right_items) - 1,
                representation_result.embedding_dimension,
            ),
        },
    )


__all__ = ["evaluate_frechet_distance"]
