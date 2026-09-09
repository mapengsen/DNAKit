"""Reproducible numeric, sequence-composition and frozen-embedding ML inputs."""

from __future__ import annotations

from typing import Any

from dnakit.core import DNASequence
from dnakit.exceptions import InputFormatError
from dnakit.representations import (
    RepresentationBackend,
    RepresentationConfig,
    extract_representations,
)

from .data import PropertyData
from .ml_models import require_ml_module
from .schema import CustomPropertyTrainingConfig


def machine_learning_features(
    data: PropertyData,
    config: CustomPropertyTrainingConfig,
    backend: RepresentationBackend | None,
) -> tuple[Any, dict[str, object]]:
    """Extract features without fitting preprocessing on held-out rows."""
    np = require_ml_module("numpy")
    metadata: dict[str, object] = {
        "source": config.feature_source,
        "kmer_k": config.kmer_k,
        "schema_version": "dnakit.ml-features.v1",
    }
    names: tuple[str, ...]
    if config.feature_source == "numeric":
        names = data.inspection.mapping.feature_fields
        values = np.asarray(data.features, dtype=np.float64)
    elif config.feature_source == "sequence":
        from dnakit.fingerprints import kmer_fingerprint

        rows: list[list[float]] = []
        names = ()
        for sequence in data.sequences:
            fingerprint = kmer_fingerprint(
                DNASequence(sequence, alphabet="iupac"),
                k=config.kmer_k,
                mode="frequency",
                ambiguity_policy="ignore",
            )
            canonical = sum(sequence.count(base) for base in "ACGT")
            rows.append(
                [
                    float(len(sequence)),
                    (sequence.count("G") + sequence.count("C")) / canonical if canonical else 0.0,
                    1.0 - canonical / len(sequence),
                    *fingerprint.values,
                ]
            )
            names = (
                "length",
                "gc_fraction",
                "ambiguity_fraction",
                *(f"kmer_{word}" for word in fingerprint.feature_names),
            )
        values = np.asarray(rows, dtype=np.float64)
    else:
        from .training import _records

        representation = extract_representations(
            _records(data),
            config=RepresentationConfig(
                model=config.embedding_model,
                checkpoint_dir=config.checkpoint_dir,
                checkpoint_path=config.checkpoint_path,
                model_source_path=config.model_source_path,
                pooling=config.pooling,
                device=config.device,
                dtype=config.dtype,
                batch_size=config.batch_size,
                max_length=config.max_length,
                max_records=config.max_records,
                show_progress=config.show_progress,
                allow_remote_code=config.allow_remote_code,
                hf_token=config.hf_token,
            ),
            backend=backend,
        )
        values = np.asarray(representation.representations, dtype=np.float64)
        names = tuple(f"embedding_{index}" for index in range(representation.embedding_dimension))
        metadata.update(
            embedding_model=config.embedding_model,
            checkpoint_path=representation.checkpoint_path,
            model_source_path=config.model_source_path,
            pooling=config.pooling,
            max_length=config.max_length,
            dtype=config.dtype,
        )
    if (
        values.ndim != 2
        or values.shape != (data.record_count, len(names))
        or not names
        or not bool(np.isfinite(values).all())
    ):
        raise InputFormatError(
            "Machine-learning inputs must be a non-empty finite feature matrix.",
            code="INVALID_CUSTOM_FEATURE",
        )
    metadata["names"] = list(names)
    metadata["dimension"] = len(names)
    return values, metadata
