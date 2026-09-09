"""Supervised ML training and trusted artifact inference through the property API."""

from __future__ import annotations

import platform
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import BackendExecutionError, InputFormatError
from dnakit.representations import RepresentationBackend

from .data import DataReader, PropertyData, inspect_property_data, load_property_data
from .llm import SchemaLLMClient
from .ml_features import machine_learning_features
from .ml_models import create_estimator, get_machine_learning_model, require_ml_module
from .schema import (
    CustomPropertyPredictionResult,
    CustomPropertyTrainingConfig,
    CustomPropertyTrainingResult,
    TargetSchema,
)
from .splitting import split_property_data
from .targets import compute_metrics, decoded_prediction_rows, encode_targets, fit_target_schema


def _fit_estimators(
    values: Any, targets: Any, config: CustomPropertyTrainingConfig
) -> tuple[list[Any], list[str]]:
    np = require_ml_module("numpy")
    preprocessing = require_ml_module("sklearn.preprocessing")
    pipeline_type = require_ml_module("sklearn.pipeline").Pipeline
    spec = get_machine_learning_model(config.model)
    estimators: list[Any] = []
    notes: list[str] = []
    if spec.nonnegative_features and bool((values < 0).any()):
        raise InputFormatError(
            "This Naive Bayes model requires nonnegative features; choose another model or input.",
            code="NEGATIVE_MACHINE_LEARNING_FEATURE",
        )
    matrix = targets.reshape(-1, 1) if targets.ndim == 1 else targets
    for column in range(matrix.shape[1]):
        target = matrix[:, column]
        estimator = create_estimator(config)
        if config.task_type == "classification":
            target = target.astype(np.int64)
            classes, counts = np.unique(target, return_counts=True)
            if len(classes) == 1:
                estimator = require_ml_module("sklearn.dummy").DummyClassifier(
                    strategy="constant", constant=int(classes[0])
                )
                notes.append(
                    f"Target {column} is constant in training; a constant classifier was fitted."
                )
        scaler = (
            (
                preprocessing.MaxAbsScaler()
                if spec.nonnegative_features
                else preprocessing.StandardScaler()
            )
            if config.standardize_features
            else "passthrough"
        )
        pipeline = pipeline_type([("scale", scaler), ("model", estimator)])
        if config.task_type == "classification" and not hasattr(estimator, "predict_proba"):
            folds = min(3, int(counts.min()))
            if folds < 2:
                raise InputFormatError(
                    "Probability calibration requires at least two training samples per class.",
                    code="INSUFFICIENT_CALIBRATION_SAMPLES",
                )
            cv = require_ml_module("sklearn.model_selection").StratifiedKFold(
                n_splits=folds, shuffle=True, random_state=config.seed
            )
            # The scaler is inside CV, so calibration folds do not share fitted transforms.
            pipeline = require_ml_module("sklearn.calibration").CalibratedClassifierCV(
                estimator=pipeline, method="sigmoid", cv=cv, n_jobs=1
            )
            notes.append(
                f"Target {column} uses {folds}-fold sigmoid probability calibration "
                "on training rows."
            )
        try:
            pipeline.fit(values, target)
        except Exception as exc:
            raise BackendExecutionError(
                "The machine-learning estimator could not fit the training data.",
                code="MACHINE_LEARNING_FIT_FAILED",
                context={"model": config.model, "target_index": column, "reason": str(exc)},
            ) from exc
        estimators.append(pipeline)
    return estimators, notes


def _outputs(estimators: list[Any], values: Any, schema: TargetSchema) -> Any:
    """Adapt genuine probabilities to the existing logits-based reporting contract."""
    np = require_ml_module("numpy")
    if schema.task_type == "regression":
        result = np.column_stack([model.predict(values) for model in estimators])
    else:
        blocks: list[Any] = []
        for column, estimator in enumerate(estimators):
            classes = schema.class_names[column]
            probabilities = np.asarray(estimator.predict_proba(values), dtype=np.float64)
            aligned = np.zeros((len(values), len(classes)), dtype=np.float64)
            aligned[:, np.asarray(estimator.classes_, dtype=np.int64)] = probabilities
            if (
                not bool(np.isfinite(aligned).all())
                or bool((aligned < 0).any())
                or not bool(np.allclose(aligned.sum(axis=1), 1.0))
            ):
                raise BackendExecutionError(
                    "The classifier returned invalid probabilities.",
                    code="INVALID_CUSTOM_PREDICTION_OUTPUT",
                )
            aligned = np.clip(aligned, 1e-15, 1 - 1e-15)
            if schema.mode == "multilabel":
                positive = aligned[:, 1:2]
                blocks.append(np.log(positive) - np.log1p(-positive))
            else:
                blocks.append(np.log(aligned))
        result = np.concatenate(blocks, axis=1)
    if result.shape != (len(values), schema.output_dimension) or not bool(
        np.isfinite(result).all()
    ):
        raise BackendExecutionError(
            "The estimator returned invalid predictions.", code="INVALID_CUSTOM_PREDICTION_OUTPUT"
        )
    return result


def _versions(model: str) -> dict[str, str]:
    versions = {"python": platform.python_version()}
    modules = ["numpy", "scipy", "sklearn", "joblib"]
    dependency = get_machine_learning_model(model).dependency
    if dependency != "scikit-learn":
        modules.append(dependency)
    for name in modules:
        versions[name] = str(require_ml_module(name).__version__)
    return versions


def train_machine_learning_property(
    source: str | Path,
    *,
    config: CustomPropertyTrainingConfig,
    llm_client: SchemaLLMClient | None,
    data_reader: DataReader | None,
    representation_backend: RepresentationBackend | None,
) -> CustomPropertyTrainingResult:
    from .training import (
        _ARTIFACT_SCHEMA,
        _staged_output,
        _validated_output_dir,
        _write_json,
        _write_predictions,
        _write_splits,
    )

    # Fail early for an unavailable optional package or unsupported constructor parameter.
    create_estimator(config)
    output_dir = _validated_output_dir(config.output_dir, source=source, config=config)
    inspection = inspect_property_data(
        source,
        mapping=config.mapping,
        llm=config.llm,
        llm_client=llm_client,
        reader=data_reader,
        max_records=config.max_records,
        max_file_bytes=config.max_file_bytes,
    )
    data = load_property_data(
        source,
        inspection=inspection,
        reader=data_reader,
        max_records=config.max_records,
        max_file_bytes=config.max_file_bytes,
        missing_label_policy=config.missing_label_policy,
        require_sequences=config.feature_source != "numeric",
    )
    split = split_property_data(
        data,
        train_ratio=config.train_ratio,
        validation_ratio=config.validation_ratio,
        test_ratio=config.test_ratio,
        seed=config.seed,
        task_type=config.task_type,
        stratify=config.stratify,
    )
    schema = fit_target_schema(
        data.labels,
        data.label_names,
        task_type=config.task_type,
        classification_mode=config.classification_mode,
        train_indices=split.train,
    )
    targets = encode_targets(data.labels, schema=schema)
    values, feature_metadata = machine_learning_features(data, config, representation_backend)
    warnings = [*inspection.warnings, *split.warnings]
    if data.dropped_records:
        warnings.append(f"Dropped {data.dropped_records} records with missing labels.")
    if config.feature_source != "numeric" and len(set(data.sequences)) < data.record_count:
        warnings.append(
            "Duplicate sequences exist; verify that the split avoids train/test leakage."
        )
    manifest: dict[str, object] = {
        "schema_version": _ARTIFACT_SCHEMA,
        "model_name": config.model,
        "task_type": config.task_type,
        "finetune_method": "machine_learning",
        "target_schema": schema.to_dict(),
        "threshold": config.threshold,
        "input_dimension": values.shape[1],
        "output_dimension": schema.output_dimension,
        "features": feature_metadata,
        "data_mapping": inspection.mapping.to_dict(),
        "versions": _versions(config.model),
        "training": {
            "seed": config.seed,
            "workers": config.workers,
            "model_params": to_json_compatible(config.model_params),
            "standardize_features": config.standardize_features,
            "multioutput_strategy": "independent_estimator_per_target",
            "ratios": {
                "train": config.train_ratio,
                "validation": config.validation_ratio,
                "test": config.test_ratio,
            },
        },
    }
    artifacts = {
        "manifest": "model.json",
        "estimator": "estimator.joblib",
        "metrics": "metrics.json",
        "data_inspection": "data_inspection.json",
        "splits": "splits.csv",
        "test_predictions": "test_predictions.jsonl",
    }
    with _staged_output(output_dir, overwrite=config.overwrite) as staging:
        # Limit OpenMP/BLAS as well as estimators exposing n_jobs.
        with require_ml_module("threadpoolctl").threadpool_limits(limits=max(1, config.workers)):
            estimators, notes = _fit_estimators(
                values[list(split.train)], targets[list(split.train)], config
            )
            warnings.extend(notes)
            metrics: dict[str, object] = {}
            for name, indices in (("validation", split.validation), ("test", split.test)):
                outputs = _outputs(estimators, values[list(indices)], schema)
                metrics[name] = compute_metrics(
                    targets[list(indices)], outputs, schema=schema, threshold=config.threshold
                )
                if name == "test":
                    _write_predictions(
                        staging / "test_predictions.jsonl",
                        tuple(data.record_ids[index] for index in indices),
                        outputs,
                        schema=schema,
                        threshold=config.threshold,
                    )
        require_ml_module("joblib").dump(
            {"estimators": estimators}, staging / "estimator.joblib", compress=3
        )
        manifest["warnings"] = warnings
        _write_json(staging / "model.json", manifest)
        _write_json(staging / "metrics.json", metrics)
        _write_json(staging / "data_inspection.json", inspection.to_dict())
        _write_splits(staging / "splits.csv", data, split)
    return CustomPropertyTrainingResult(
        config.model,
        config.task_type,
        "machine_learning",
        str(output_dir),
        schema,
        FrozenDict(split.counts),
        FrozenDict(metrics),
        (),
        FrozenDict({name: str(output_dir / filename) for name, filename in artifacts.items()}),
        inspection,
        tuple(warnings),
    )


def predict_machine_learning_property(
    root: Path,
    manifest: Mapping[str, object],
    data: PropertyData,
    schema: TargetSchema,
    *,
    representation_backend: RepresentationBackend | None,
    batch_size: int,
    device: str,
    threshold: float | None,
    allow_remote_code: bool,
    show_progress: bool,
) -> CustomPropertyPredictionResult:
    features = manifest.get("features")
    if (
        not isinstance(features, Mapping)
        or features.get("schema_version") != "dnakit.ml-features.v1"
    ):
        raise InputFormatError(
            "Invalid saved feature schema.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    selected_threshold = manifest.get("threshold", 0.5) if threshold is None else threshold
    try:
        config = CustomPropertyTrainingConfig(
            task_type=schema.task_type,
            model=str(manifest["model_name"]),
            finetune_method="machine_learning",
            feature_source=cast(Any, features["source"]),
            embedding_model=str(features.get("embedding_model", "dnabert2")),
            mapping=data.inspection.mapping,
            kmer_k=cast(int, features.get("kmer_k", 3)),
            checkpoint_path=cast(Any, features.get("checkpoint_path")),
            model_source_path=cast(Any, features.get("model_source_path")),
            pooling=cast(Any, features.get("pooling", "mean")),
            max_length=cast(Any, features.get("max_length")),
            dtype=cast(Any, features.get("dtype", "auto")),
            batch_size=batch_size,
            device=device,
            threshold=cast(float, selected_threshold),
            allow_remote_code=allow_remote_code,
            show_progress=show_progress,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise InputFormatError(
            "Invalid saved ML configuration.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        ) from exc
    values, actual = machine_learning_features(data, config, representation_backend)
    if actual["names"] != features.get("names") or values.shape[1] != manifest.get(
        "input_dimension"
    ):
        raise InputFormatError(
            "Feature names/order or embedding dimension differ from training.",
            code="CUSTOM_FEATURE_SCHEMA_MISMATCH",
        )
    current_versions = _versions(config.model)
    if manifest.get("versions") != current_versions:
        raise InputFormatError(
            "Load this ML artifact with the Python and library versions recorded in model.json.",
            code="MACHINE_LEARNING_VERSION_MISMATCH",
            context={"saved": manifest.get("versions"), "current": current_versions},
        )
    try:
        bundle = require_ml_module("joblib").load(root / "estimator.joblib")
        estimators = bundle["estimators"]
        expected = 1 if schema.mode in {"binary", "multiclass"} else len(schema.label_names)
        if not isinstance(estimators, list) or len(estimators) != expected:
            raise ValueError("Invalid estimator count")
    except Exception as exc:
        raise InputFormatError(
            "Could not load the ML estimator artifact.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        ) from exc
    with require_ml_module("threadpoolctl").threadpool_limits(limits=1):
        outputs = _outputs(estimators, values, schema)
    predictions, probabilities = decoded_prediction_rows(
        data.record_ids, outputs, schema=schema, threshold=config.threshold
    )
    return CustomPropertyPredictionResult(
        str(root),
        data.record_ids,
        schema.label_names,
        tuple(FrozenDict(row) for row in predictions),
        tuple(FrozenDict(row) for row in probabilities),
    )
