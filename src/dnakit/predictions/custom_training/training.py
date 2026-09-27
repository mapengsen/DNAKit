"""Training, artifact persistence, and inference for user-defined properties."""

from __future__ import annotations

import csv
import json
import math
import os
import random
import shutil
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import closing, contextmanager, suppress
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

from dnakit._checkpoint_status import checkpoint_status
from dnakit._dependency_hints import dependency_hint
from dnakit._progress import track
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core import DNARecord, DNASequence
from dnakit.core._json import FrozenDict
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    InputFormatError,
)
from dnakit.io import TableSchema, export_table
from dnakit.representations import (
    RepresentationBackend,
    RepresentationConfig,
    ensure_model_checkpoint,
    extract_representations,
    get_embedding_model,
)

from ._loading import _PropertyBatchLoader
from .data import DataReader, PropertyData, inspect_property_data, load_property_data
from .encoder import TrainableSequenceEncoder
from .encoder_factory import create_trainable_encoder
from .llm import SchemaLLMClient
from .schema import (
    CustomPropertyPredictionResult,
    CustomPropertyTrainingConfig,
    CustomPropertyTrainingResult,
    DataMapping,
    LLMSchemaConfig,
    TargetSchema,
)
from .splitting import SplitAssignment, split_property_data
from .targets import (
    compute_metrics,
    decoded_prediction_rows,
    encode_targets,
    fit_target_schema,
)

_ARTIFACT_SCHEMA = "dnakit.custom-property-model.v1"


def _numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - supplied by the training extra.
        raise BackendUnavailableError(
            "Custom property training requires NumPy.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("training"),
        ) from exc
    return np


def _torch() -> Any:
    try:
        import torch  # type: ignore[import-not-found]
    except ImportError as exc:
        raise BackendUnavailableError(
            "Custom property training requires PyTorch.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("training"),
        ) from exc
    return torch


def _resolve_device(torch: Any, requested: str) -> Any:
    from dnakit.representations.backends import _resolve_torch_device

    return _resolve_torch_device(torch, requested)


def _seed_everything(torch: Any, seed: int) -> None:
    np = _numpy()
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if bool(torch.cuda.is_available()):
        torch.cuda.manual_seed_all(seed)
    with suppress(AttributeError, RuntimeError):
        torch.use_deterministic_algorithms(True, warn_only=True)


def _validated_output_dir(
    output_dir: str | os.PathLike[str],
    *,
    source: str | Path,
    config: CustomPropertyTrainingConfig,
) -> Path:
    requested = Path(output_dir).expanduser()
    if requested.is_symlink():
        raise ConfigurationError(
            "output_dir must not be a symbolic link.",
            code="UNSAFE_CUSTOM_OUTPUT_DIR",
            hint="Choose a new dedicated model-output subdirectory.",
        )
    resolved = requested.resolve()
    source_path = Path(source).expanduser().resolve()
    working_directory = Path.cwd().resolve()
    home_directory = Path.home().resolve()
    protected = {
        Path(resolved.anchor).resolve(),
        working_directory,
        *working_directory.parents,
        home_directory,
        *home_directory.parents,
    }
    protected_inputs = [source_path]
    for configured_path in (config.checkpoint_path, config.model_source_path):
        if configured_path is not None:
            protected_inputs.append(Path(configured_path).expanduser().resolve())
    if resolved in protected or any(
        resolved == input_path or resolved in input_path.parents for input_path in protected_inputs
    ):
        raise ConfigurationError(
            "output_dir must not be a protected directory or contain an input source.",
            code="UNSAFE_CUSTOM_OUTPUT_DIR",
            context={"output_dir": str(resolved)},
            hint="Choose a new dedicated model-output subdirectory.",
        )
    return resolved


def _records(data: PropertyData) -> tuple[DNARecord, ...]:
    return tuple(
        DNARecord(DNASequence(sequence, alphabet="iupac"), record_id)
        for record_id, sequence in zip(data.record_ids, data.sequences, strict=True)
    )


def _representation_config(
    config: CustomPropertyTrainingConfig,
    *,
    checkpoint_path: str | None = None,
    max_records: int | None = None,
) -> RepresentationConfig:
    return RepresentationConfig(
        model=config.model,
        checkpoint_dir=None if checkpoint_path is not None else config.checkpoint_dir,
        checkpoint_path=checkpoint_path or config.checkpoint_path,
        model_source_path=config.model_source_path,
        pooling=config.pooling,
        device=config.device,
        dtype=config.dtype,
        batch_size=config.batch_size,
        max_length=config.max_length,
        max_records=max_records or config.max_records,
        show_progress=config.show_progress,
        allow_remote_code=config.allow_remote_code,
        hf_token=config.hf_token,
    )


def _resolve_checkpoint(config: CustomPropertyTrainingConfig) -> str:
    if config.checkpoint_path is not None:
        path = Path(config.checkpoint_path).expanduser().resolve()
        if not path.is_dir():
            raise ConfigurationError(
                "checkpoint_path must be an existing checkpoint directory.",
                code="MODEL_CHECKPOINT_NOT_FOUND",
                context={"path": str(path)},
            )
        checkpoint_status("Using local checkpoint", path)
        return str(path)
    return ensure_model_checkpoint(
        config.model,
        checkpoint_dir=config.checkpoint_dir,
        hf_token=config.hf_token,
        show_progress=config.show_progress,
    ).path


def _head_loss(torch: Any, logits: Any, targets: Any, schema: TargetSchema) -> Any:
    functional = torch.nn.functional
    if schema.task_type == "regression":
        return functional.mse_loss(logits, targets.float())
    if schema.mode in {"binary", "multiclass"}:
        return functional.cross_entropy(logits, targets.long())
    if schema.mode == "multilabel":
        return functional.binary_cross_entropy_with_logits(logits, targets.float())
    losses: list[Any] = []
    start = 0
    for column, classes in enumerate(schema.class_names):
        end = start + len(classes)
        losses.append(functional.cross_entropy(logits[:, start:end], targets[:, column].long()))
        start = end
    return torch.stack(losses).mean()


def _batches(
    indices: Sequence[int],
    *,
    batch_size: int,
    shuffle: bool,
    seed: int,
) -> Iterator[tuple[int, ...]]:
    np = _numpy()
    values = np.asarray(indices, dtype=np.int64)
    if shuffle:
        values = np.random.default_rng(seed).permutation(values)
    for start in range(0, len(values), batch_size):
        yield tuple(int(value) for value in values[start : start + batch_size])


def _features(
    torch: Any,
    batch: Sequence[int],
    *,
    method: str,
    sequences: Sequence[str],
    representations: Any | None,
    encoder: TrainableSequenceEncoder | None,
    device: Any,
) -> Any:
    if method == "linear_probe":
        if representations is None:
            raise AssertionError("Linear probing requires frozen representations.")
        return torch.as_tensor(representations[list(batch)], dtype=torch.float32, device=device)
    if encoder is None:
        raise AssertionError("Full fine-tuning requires an encoder.")
    encoded = encoder.encode([sequences[index] for index in batch])
    if not isinstance(encoded, torch.Tensor) or encoded.ndim != 2:
        raise BackendExecutionError(
            "The trainable encoder must return a two-dimensional PyTorch tensor.",
            code="INVALID_CUSTOM_ENCODER_OUTPUT",
        )
    return encoded.float()


def _loaded_features(
    torch: Any,
    batch: Mapping[str, Any],
    *,
    encoder: TrainableSequenceEncoder | None,
    device: Any,
) -> Any:
    if encoder is None:
        return batch["representations"].to(device=device, dtype=torch.float32)
    encoded = encoder.encode(batch["sequence"])
    if (
        not isinstance(encoded, torch.Tensor)
        or encoded.ndim != 2
        or encoded.shape[0] != len(batch["sequence"])
    ):
        raise BackendExecutionError(
            "The trainable encoder must return one tensor row per sequence.",
            code="INVALID_CUSTOM_ENCODER_OUTPUT",
        )
    return encoded.float()


def _save_head(path: Path, head: Any) -> None:
    np = _numpy()
    weight = head.weight.detach().float().cpu().numpy()
    bias = head.bias.detach().float().cpu().numpy()
    np.savez_compressed(path, weight=weight, bias=bias)


def _load_head(
    torch: Any, path: Path, *, input_dimension: int, output_dimension: int, device: Any
) -> Any:
    np = _numpy()
    try:
        with np.load(path, allow_pickle=False) as payload:
            weight = np.asarray(payload["weight"], dtype=np.float32)
            bias = np.asarray(payload["bias"], dtype=np.float32)
    except (OSError, KeyError, ValueError) as exc:
        raise InputFormatError(
            "The saved custom-property head is invalid.",
            code="INVALID_CUSTOM_MODEL_ARTIFACT",
        ) from exc
    if weight.shape != (output_dimension, input_dimension) or bias.shape != (output_dimension,):
        raise InputFormatError(
            "The saved head dimensions do not match the model manifest.",
            code="INVALID_CUSTOM_MODEL_ARTIFACT",
            context={"weight_shape": tuple(int(value) for value in weight.shape)},
        )
    head = torch.nn.Linear(input_dimension, output_dimension).to(device)
    with torch.no_grad():
        head.weight.copy_(torch.from_numpy(weight).to(device))
        head.bias.copy_(torch.from_numpy(bias).to(device))
    return head


def _evaluate(
    torch: Any,
    indices: Sequence[int],
    *,
    head: Any,
    schema: TargetSchema,
    batch_loader: _PropertyBatchLoader,
    encoder: TrainableSequenceEncoder | None,
    device: Any,
) -> tuple[float, Any, Any]:
    head.eval()
    if encoder is not None:
        encoder.eval()
    logits: list[Any] = []
    truths: list[Any] = []
    weighted_loss = 0.0
    count = 0
    with torch.no_grad():
        for batch in batch_loader.batches(indices, shuffle=False, seed=0):
            features = _loaded_features(
                torch,
                batch,
                encoder=encoder,
                device=device,
            )
            batch_logits = head(features)
            targets = batch["targets"].to(device)
            loss = _head_loss(torch, batch_logits, targets, schema)
            weighted_loss += float(loss.detach().cpu()) * len(targets)
            count += len(targets)
            logits.append(batch_logits.detach().float().cpu().numpy())
            truths.append(targets.detach().cpu().numpy())
    np = _numpy()
    return weighted_loss / count, np.concatenate(truths), np.concatenate(logits)


def _metric_value(metrics: Mapping[str, object], name: str) -> float | None:
    section: Mapping[str, object] = metrics
    aggregate = metrics.get("aggregate")
    if isinstance(aggregate, Mapping):
        section = aggregate
    key = "matthews_correlation_coefficient" if name == "mcc" else name
    if key == "accuracy" and "subset_accuracy" in section:
        key = "subset_accuracy"
    value = section.get(key)
    return float(value) if isinstance(value, (int, float)) else None


def _train_head_and_encoder(
    split: SplitAssignment,
    schema: TargetSchema,
    config: CustomPropertyTrainingConfig,
    staging: Path,
    *,
    batch_loader: _PropertyBatchLoader,
    encoder: TrainableSequenceEncoder | None,
    input_dimension: int,
) -> tuple[Any, TrainableSequenceEncoder | None, tuple[dict[str, object], ...]]:
    torch = _torch()
    if config.monitor == "mcc" and schema.mode == "multilabel":
        raise ConfigurationError(
            "monitor='mcc' requires binary, multiclass, or multioutput classification.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
            hint="Use monitor='f1_macro' or 'validation_loss' for multilabel targets.",
        )
    device = encoder.device if encoder is not None else _resolve_device(torch, config.device)
    _seed_everything(torch, config.seed)
    head = torch.nn.Linear(input_dimension, schema.output_dimension).to(device)
    parameters = list(head.parameters())
    if encoder is not None:
        parameters.extend(
            parameter for parameter in encoder.parameters() if parameter.requires_grad
        )
        if any(parameter.dtype == torch.float16 for parameter in parameters):
            raise ConfigurationError(
                "Full fine-tuning requires float32 or bfloat16 optimizer parameters.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
                hint=(
                    "Set dtype='float32' or 'bfloat16'; "
                    "float16 AdamW can produce non-finite weights."
                ),
            )
    optimizer = torch.optim.AdamW(
        parameters,
        lr=config.resolved_learning_rate,
        betas=(float(config.beta1), float(config.beta2)),
        weight_decay=float(config.weight_decay),
    )
    monitor_mode = config.resolved_monitor_mode
    scheduler = (
        torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode=monitor_mode,
            patience=config.scheduler_patience,
            factor=float(config.scheduler_factor),
            # An absolute threshold also handles zero and negative MCC correctly.
            threshold=float(config.monitor_min_delta),
            threshold_mode="abs",
        )
        if config.lr_scheduler == "reduce_on_plateau"
        else None
    )
    history: list[dict[str, object]] = []
    best_value = float("inf") if monitor_mode == "min" else -float("inf")
    stale_epochs = 0
    head_path = staging / "head.npz"
    backbone_dir = staging / "backbone"
    for epoch in range(1, config.epochs + 1):
        learning_rate = float(optimizer.param_groups[0]["lr"])
        head.train()
        if encoder is not None:
            encoder.train()
        optimizer.zero_grad(set_to_none=True)
        total_loss = 0.0
        count = 0
        epoch_batches = batch_loader.batches(
            split.train,
            shuffle=True,
            seed=config.seed + epoch,
        )
        batch_count = (len(split.train) + config.batch_size - 1) // config.batch_size
        progress_batches = track(
            epoch_batches,
            description=f"Training {config.model}: epoch {epoch}/{config.epochs}",
            total=batch_count,
            disable=not config.show_progress,
        )
        for batch_index, batch in enumerate(progress_batches, start=1):
            features = _loaded_features(
                torch,
                batch,
                encoder=encoder,
                device=device,
            )
            logits = head(features)
            targets = batch["targets"].to(device)
            loss = _head_loss(torch, logits, targets, schema)
            if not bool(torch.isfinite(loss)):
                raise BackendExecutionError(
                    "Training produced a non-finite loss.",
                    code="NONFINITE_CUSTOM_TRAINING_LOSS",
                    context={"epoch": epoch, "batch": batch_index},
                )
            if encoder is not None and not features.requires_grad:
                raise BackendExecutionError(
                    "Full fine-tuning requires differentiable backbone representations.",
                    code="CUSTOM_ENCODER_GRADIENT_MISSING",
                )
            accumulation_start = (
                ((batch_index - 1) // config.gradient_accumulation_steps)
                * config.gradient_accumulation_steps
                * config.batch_size
            )
            accumulation_samples = min(
                config.gradient_accumulation_steps * config.batch_size,
                len(split.train) - accumulation_start,
            )
            (loss * (len(targets) / accumulation_samples)).backward()
            total_loss += float(loss.detach().cpu()) * len(targets)
            count += len(targets)
            should_step = (
                batch_index % config.gradient_accumulation_steps == 0 or batch_index == batch_count
            )
            if should_step:
                try:
                    torch.nn.utils.clip_grad_norm_(
                        parameters, float(config.max_grad_norm), error_if_nonfinite=True
                    )
                except RuntimeError as exc:
                    raise BackendExecutionError(
                        "Training produced non-finite gradients.",
                        code="NONFINITE_CUSTOM_TRAINING_GRADIENT",
                        context={"epoch": epoch, "batch": batch_index},
                    ) from exc
                optimizer.step()
                optimizer.zero_grad(set_to_none=True)
        train_metrics = {}
        if config.log_train_metrics:
            _, train_targets, train_logits = _evaluate(
                torch,
                split.train,
                head=head,
                schema=schema,
                batch_loader=batch_loader,
                encoder=encoder,
                device=device,
            )
            train_metrics = compute_metrics(
                train_targets,
                train_logits,
                schema=schema,
                threshold=config.threshold,
            )
        validation_loss, validation_metrics = None, {}
        if split.validation:
            validation_loss, validation_targets, validation_logits = _evaluate(
                torch,
                split.validation,
                head=head,
                schema=schema,
                batch_loader=batch_loader,
                encoder=encoder,
                device=device,
            )
            validation_metrics = compute_metrics(
                validation_targets,
                validation_logits,
                schema=schema,
                threshold=config.threshold,
            )
        monitor_value = (
            validation_loss
            if config.monitor == "validation_loss"
            else _metric_value(validation_metrics, config.monitor)
        )
        if not split.validation:
            monitor_value = total_loss / count
        if monitor_value is None or not math.isfinite(monitor_value):
            raise BackendExecutionError(
                "The monitored validation metric is unavailable or non-finite.",
                code="INVALID_CUSTOM_MONITOR_VALUE",
                context={"epoch": epoch, "monitor": config.monitor, "value": monitor_value},
            )
        if scheduler is not None:
            scheduler.step(monitor_value)
        primary_metric_name = (
            ("rmse" if schema.task_type == "regression" else "f1_macro")
            if config.monitor == "validation_loss"
            else config.monitor
        )
        history.append(
            {
                "epoch": epoch,
                "train_loss": total_loss / count,
                "validation_loss": validation_loss,
                "primary_metric": _metric_value(validation_metrics, primary_metric_name),
                "primary_metric_name": primary_metric_name,
                "monitor": config.monitor if split.validation else "train_loss",
                "train_mcc": _metric_value(train_metrics, "mcc"),
                "validation_mcc": _metric_value(validation_metrics, "mcc"),
                "monitor_value": monitor_value,
                "monitor_mode": monitor_mode if split.validation else "min",
                "learning_rate": learning_rate,
                "next_learning_rate": float(optimizer.param_groups[0]["lr"]),
            }
        )
        improved = (
            monitor_value < best_value - config.monitor_min_delta
            if monitor_mode == "min"
            else monitor_value > best_value + config.monitor_min_delta
        )
        if improved or not split.validation:
            best_value = monitor_value
            stale_epochs = 0
            _save_head(head_path, head)
            if encoder is not None:
                if backbone_dir.exists():
                    shutil.rmtree(backbone_dir)
                encoder.save_pretrained(backbone_dir)
        else:
            stale_epochs += 1
            if (
                config.early_stopping_patience > 0
                and stale_epochs >= config.early_stopping_patience
            ):
                break
    best_head = _load_head(
        torch,
        head_path,
        input_dimension=input_dimension,
        output_dimension=schema.output_dimension,
        device=device,
    )
    if encoder is not None:
        restore = getattr(encoder, "load_pretrained", None)
        if callable(restore):
            restore(backbone_dir)
    return best_head, encoder, tuple(history)


@contextmanager
def _staged_output(output_dir: Path, *, overwrite: bool) -> Iterator[Path]:
    parent = output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    if output_dir.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}-", dir=parent))
    backup: Path | None = None
    try:
        yield staging
        if output_dir.exists():
            backup = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}-backup-", dir=parent))
            backup.rmdir()
            os.replace(output_dir, backup)
        os.replace(staging, output_dir)
        if backup is not None:
            shutil.rmtree(backup)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        if backup is not None and backup.exists() and not output_dir.exists():
            os.replace(backup, output_dir)
        raise


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_splits(path: Path, data: PropertyData, split: SplitAssignment) -> None:
    assigned = {
        **{index: "train" for index in split.train},
        **{index: "validation" for index in split.validation},
        **{index: "test" for index in split.test},
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("source_row", "record_id", "split"))
        for index, record_id in enumerate(data.record_ids):
            writer.writerow((index, record_id, assigned[index]))


def _write_predictions(
    path: Path,
    record_ids: Sequence[str],
    logits: Any,
    *,
    schema: TargetSchema,
    threshold: float,
) -> None:
    predicted, probabilities = decoded_prediction_rows(
        record_ids,
        logits,
        schema=schema,
        threshold=threshold,
    )
    by_id = {str(row["record_id"]): row for row in probabilities}
    with path.open("w", encoding="utf-8") as handle:
        for row in predicted:
            record_id = str(row["record_id"])
            payload = dict(row)
            if record_id in by_id:
                payload["probabilities"] = {
                    key: value for key, value in by_id[record_id].items() if key != "record_id"
                }
            handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def _manifest(
    config: CustomPropertyTrainingConfig,
    data: PropertyData,
    schema: TargetSchema,
    *,
    checkpoint_path: str | None,
    input_dimension: int,
) -> dict[str, object]:
    return {
        "schema_version": _ARTIFACT_SCHEMA,
        "model_name": config.model,
        "task_type": config.task_type,
        "finetune_method": config.finetune_method,
        "target_schema": schema.to_dict(),
        "input_dimension": input_dimension,
        "output_dimension": schema.output_dimension,
        "threshold": config.threshold,
        "pooling": config.pooling,
        "dtype": config.dtype,
        "max_length": config.max_length,
        "checkpoint_path": checkpoint_path,
        "model_source_path": (
            os.fspath(config.model_source_path) if config.model_source_path is not None else None
        ),
        "requires_remote_code": get_embedding_model(config.model).trust_remote_code,
        "data_mapping": data.inspection.mapping.to_dict(),
        "training": {
            "seed": config.seed,
            "ratios": {
                "train": config.train_ratio,
                "validation": config.validation_ratio,
                "test": config.test_ratio,
            },
            "epochs_requested": config.epochs,
            "batch_size": config.batch_size,
            "batch_size_scope": "per_device",
            "workers": config.workers,
            "optimizer": "AdamW",
            "learning_rate": config.resolved_learning_rate,
            "beta1": config.beta1,
            "beta2": config.beta2,
            "weight_decay": config.weight_decay,
            "gradient_accumulation_steps": config.gradient_accumulation_steps,
            "early_stopping_patience": config.early_stopping_patience,
            "monitor": config.monitor,
            "monitor_mode": config.resolved_monitor_mode,
            "monitor_min_delta": config.monitor_min_delta,
            "lr_scheduler": config.lr_scheduler,
            "scheduler_patience": config.scheduler_patience,
            "scheduler_factor": config.scheduler_factor,
            "scheduler_threshold_mode": "abs",
            "log_train_metrics": config.log_train_metrics,
        },
    }


def train_custom_property_model(
    source: str | Path,
    *,
    config: CustomPropertyTrainingConfig,
    llm_client: SchemaLLMClient | None = None,
    data_reader: DataReader | None = None,
    representation_backend: RepresentationBackend | None = None,
    training_encoder: TrainableSequenceEncoder | None = None,
    verbose: bool = True,
) -> CustomPropertyTrainingResult:
    """Train and evaluate a user-defined sequence property model from one file.

    The default split is seeded 80/10/10.  A mapped split field takes precedence.
    Full fine-tuning updates every exposed backbone parameter; ``linear_probe``
    freezes the selected DNA model and trains only a linear head. Both methods
    default to batch_size=32 on the selected device, workers=6, and AdamW with
    lr=1e-4, betas=(0.9, 0.95), weight_decay=0.1. Use a __main__ guard in scripts
    with workers > 0; notebooks may use workers=0. This trainer is single-device.

    Set ``monitor='mcc', lr_scheduler='reduce_on_plateau'`` to drive scheduling,
    early stopping, and best-checkpoint selection with validation MCC (maximized
    in auto mode). ``history.json`` records the monitored value and each epoch's
    learning rate before and after the scheduler step; ``model.json`` records
    the configured optimizer, scheduler, and stopping settings.

    ``machine_learning`` selects a supervised estimator through ``config.model``
    and uses ``model_params`` instead of neural optimizer/epoch settings. Inputs
    can be sequence features, frozen embeddings, or explicit numeric columns.
    For full K-fold validation on an existing training file, set config.cv_folds
    (e.g. 10). All folds are evaluated, then the saved ML model is refitted on all
    source rows; split ratios are ignored in this mode. Pass independent test
    data only to predict_custom_properties afterwards. Fold assignments, OOF
    predictions and validation metrics are saved under config.output_dir.
    """

    if not isinstance(config, CustomPropertyTrainingConfig):
        raise ConfigurationError(
            "config must be CustomPropertyTrainingConfig.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
        )
    if config.finetune_method == "machine_learning":
        if training_encoder is not None:
            raise ConfigurationError(
                "training_encoder is only supported for full deep-learning fine-tuning.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        from .ml_training import train_machine_learning_property

        return train_machine_learning_property(
            source,
            config=config,
            llm_client=llm_client,
            data_reader=data_reader,
            representation_backend=representation_backend,
        )
    if training_encoder is not None and config.finetune_method != "full":
        raise ConfigurationError(
            "training_encoder can only be supplied for full fine-tuning.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
        )
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
    if not split.validation and (
        config.early_stopping_patience > 0 or config.lr_scheduler != "none"
    ):
        raise ConfigurationError(
            "Training without validation requires early_stopping_patience=0 "
            "and lr_scheduler='none'.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
        )
    schema = fit_target_schema(
        data.labels,
        data.label_names,
        task_type=config.task_type,
        classification_mode=config.classification_mode,
        train_indices=split.train,
    )
    encoded_targets = encode_targets(data.labels, schema=schema)
    representations = None
    _seed_everything(_torch(), config.seed)
    encoder = training_encoder
    checkpoint_path: str | None = None
    input_dimension: int
    if config.finetune_method == "linear_probe":
        representation_result = extract_representations(
            _records(data),
            config=_representation_config(config),
            backend=representation_backend,
        )
        representations = representation_result.representations
        checkpoint_path = representation_result.checkpoint_path or (
            str(Path(config.checkpoint_path).expanduser().resolve())
            if config.checkpoint_path is not None
            else None
        )
        input_dimension = representation_result.embedding_dimension
    else:
        checkpoint_path = (
            _resolve_checkpoint(config)
            if encoder is None
            else (os.fspath(config.checkpoint_path) if config.checkpoint_path is not None else None)
        )
        if encoder is None:
            encoder = create_trainable_encoder(
                get_embedding_model(config.model),
                cast(str, checkpoint_path),
                device=config.device,
                dtype=config.dtype,
                pooling=config.pooling,
                max_length=config.max_length,
                allow_remote_code=config.allow_remote_code,
                model_source_path=(
                    os.fspath(config.model_source_path)
                    if config.model_source_path is not None
                    else None
                ),
                seed=config.seed,
            )
        torch = _torch()
        encoder.eval()
        with torch.no_grad():
            sample = encoder.encode([data.sequences[split.train[0]]])
        if not isinstance(sample, torch.Tensor) or sample.ndim != 2 or sample.shape[0] != 1:
            raise BackendExecutionError(
                "The trainable encoder returned an invalid sample representation.",
                code="INVALID_CUSTOM_ENCODER_OUTPUT",
            )
        input_dimension = int(sample.shape[1])

    warnings = [*inspection.warnings, *split.warnings]
    if data.dropped_records:
        warnings.append(f"Dropped {data.dropped_records} records with missing labels.")
    sequence_duplicates = len(data.sequences) - len(set(data.sequences))
    if sequence_duplicates:
        warnings.append(
            f"Detected {sequence_duplicates} duplicate sequence rows; verify that random "
            "splitting does not create train/test leakage."
        )
    if config.finetune_method == "full" and get_embedding_model(config.model).loader in {
        "transformers",
        "janusdna",
        "evo2",
    }:
        default_limit = getattr(
            encoder,
            "max_length",
            config.max_length or min(get_embedding_model(config.model).chunk_length, 4_096),
        )
        if any(len(sequence) > default_limit for sequence in data.sequences):
            warnings.append(
                f"Input exceeds the fine-tuning limit ({default_limit}); sequences may be "
                "truncated by the encoder. Check max_length and the model's input units."
            )

    with (
        _staged_output(output_dir, overwrite=config.overwrite) as staging,
        closing(
            _PropertyBatchLoader(
                _torch(),
                data.sequences,
                encoded_targets,
                representations=representations,
                batch_size=config.batch_size,
                workers=config.workers,
                seed=config.seed,
            )
        ) as batch_loader,
    ):
        head, best_encoder, history = _train_head_and_encoder(
            split,
            schema,
            config,
            staging,
            batch_loader=batch_loader,
            encoder=encoder,
            input_dimension=input_dimension,
        )
        torch = _torch()
        device = best_encoder.device if best_encoder is not None else next(head.parameters()).device
        metrics = {}
        evaluation_outputs = {}
        for split_name, indices in (("validation", split.validation), ("test", split.test)):
            if not indices:
                continue
            loss, targets, logits = _evaluate(
                torch,
                indices,
                head=head,
                schema=schema,
                batch_loader=batch_loader,
                encoder=best_encoder,
                device=device,
            )
            metrics[split_name] = compute_metrics(
                targets,
                logits,
                schema=schema,
                threshold=config.threshold,
            )
            metrics[split_name]["loss"] = loss
            filename = f"{split_name}_predictions.jsonl"
            _write_predictions(
                staging / filename,
                tuple(data.record_ids[i] for i in indices),
                logits,
                schema=schema,
                threshold=config.threshold,
            )
            evaluation_outputs[f"{split_name}_predictions"] = str(output_dir / filename)
        manifest = _manifest(
            config,
            data,
            schema,
            checkpoint_path=checkpoint_path,
            input_dimension=input_dimension,
        )
        _write_json(staging / "model.json", manifest)
        _write_json(staging / "metrics.json", metrics)
        _write_json(staging / "history.json", list(history))
        _write_json(staging / "data_inspection.json", inspection.to_dict())
        _write_splits(staging / "splits.csv", data, split)

    artifact_names = {
        "manifest": str(output_dir / "model.json"),
        "head": str(output_dir / "head.npz"),
        "metrics": str(output_dir / "metrics.json"),
        "history": str(output_dir / "history.json"),
        "data_inspection": str(output_dir / "data_inspection.json"),
        "splits": str(output_dir / "splits.csv"),
        **evaluation_outputs,
    }
    if config.finetune_method == "full":
        artifact_names["backbone"] = str(output_dir / "backbone")
    return CustomPropertyTrainingResult(
        config.model,
        config.task_type,
        config.finetune_method,
        str(output_dir),
        schema,
        FrozenDict(split.counts),
        FrozenDict(metrics),
        tuple(FrozenDict(item) for item in history),
        FrozenDict(artifact_names),
        inspection,
        tuple(warnings),
    )


def _read_manifest(model_dir: Path) -> dict[str, object]:
    path = model_dir / "model.json"
    try:
        if path.stat().st_size > 10_000_000:
            raise InputFormatError(
                "The custom model manifest exceeds the size limit.",
                code="INVALID_CUSTOM_MODEL_ARTIFACT",
            )
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise InputFormatError(
            "Could not read the custom-property model manifest.",
            code="INVALID_CUSTOM_MODEL_ARTIFACT",
        ) from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != _ARTIFACT_SCHEMA:
        raise InputFormatError(
            "The custom-property model manifest version is unsupported.",
            code="INVALID_CUSTOM_MODEL_ARTIFACT",
        )
    return cast(dict[str, object], payload)


def _target_schema(value: object) -> TargetSchema:
    if not isinstance(value, Mapping):
        raise InputFormatError(
            "The saved target schema is invalid.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    try:
        return TargetSchema(
            task_type=cast(Any, value["task_type"]),
            mode=str(value["mode"]),
            label_names=tuple(str(item) for item in cast(Sequence[object], value["label_names"])),
            class_names=tuple(
                tuple(str(item) for item in cast(Sequence[object], classes))
                for classes in cast(Sequence[object], value.get("class_names", ()))
            ),
            regression_means=tuple(
                float(cast(Any, item))
                for item in cast(Sequence[object], value.get("regression_means", ()))
            ),
            regression_scales=tuple(
                float(cast(Any, item))
                for item in cast(Sequence[object], value.get("regression_scales", ()))
            ),
        )
    except (ConfigurationError, KeyError, TypeError, ValueError) as exc:
        raise InputFormatError(
            "The saved target schema is invalid.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        ) from exc


def _mapping_from_payload(value: object, *, include_labels: bool) -> DataMapping:
    if not isinstance(value, Mapping):
        raise InputFormatError(
            "The saved data mapping is invalid.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    payload = dict(value)
    payload["label_fields"] = tuple(payload.get("label_fields", ())) if include_labels else ()
    if not include_labels:
        payload["split_field"] = None
    payload["missing_values"] = tuple(payload.get("missing_values", ()))
    payload["sequence_transforms"] = tuple(payload.get("sequence_transforms", ()))
    try:
        return DataMapping(**payload)
    except (ConfigurationError, TypeError, ValueError) as exc:
        raise InputFormatError(
            "The saved data mapping is invalid.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        ) from exc


def _prediction_output(
    result: CustomPropertyPredictionResult,
    output_path: Path | None,
    *,
    overwrite: bool,
) -> CustomPropertyPredictionResult:
    if output_path is None:
        return result
    # Persist locally, returning a small path/count/schema result to the Agent.
    exported = export_table(
        result.predictions,
        output_path,
        format=cast(Any, output_path.suffix.lower().lstrip(".")),
        schema=TableSchema(columns=tuple(result.predictions[0])),
        overwrite=overwrite,
    )
    return replace(result, output_file=FrozenDict(exported.to_dict()))


def predict_custom_properties(
    source: str | Path,
    *,
    model_dir: str | Path,
    model_source_path: str | Path | None = None,
    mapping: DataMapping | None = None,
    llm: LLMSchemaConfig | None = None,
    llm_client: SchemaLLMClient | None = None,
    data_reader: DataReader | None = None,
    representation_backend: RepresentationBackend | None = None,
    batch_size: int = 8,
    device: str = "auto",
    threshold: float | None = None,
    allow_remote_code: bool = False,
    show_progress: bool = True,
    trust_model: bool = False,
    output_path: str | Path | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> CustomPropertyPredictionResult:
    """Run inference; ML joblib artifacts require explicit ``trust_model=True``.

    Only enable this for artifacts whose source you trust: joblib deserialization
    can execute Python code. Neural artifacts do not use this option.

    Set output_path to a .csv, .tsv, .json or .parquet file to save predicted labels
    locally. The result's to_dict() then returns only an output summary (path,
    row count and columns), avoiding all predictions in the Agent/LLM. CSV can go
    directly to visualization.plot_csv_distributions. Without output_path the
    original in-memory predictions and probabilities are returned.
    """

    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise ConfigurationError(
            "batch_size must be a positive integer.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
        )
    root = Path(model_dir).expanduser().resolve()
    destination = Path(output_path).expanduser().resolve() if output_path is not None else None
    if destination is not None:
        if destination.suffix.lower() not in {".csv", ".tsv", ".json", ".parquet"}:
            raise ConfigurationError(
                "Unsupported prediction output format.", code="INVALID_TABLE_FORMAT"
            )
        if destination in {Path(source).expanduser().resolve(), root, root / "model.json"}:
            raise ConfigurationError(
                "Prediction output must be separate from the input and saved model directory.",
                code="UNSAFE_CUSTOM_OUTPUT_PATH",
            )
        if destination.exists() and not overwrite:
            raise FileExistsError(f"Prediction output already exists: {destination}")
    manifest = _read_manifest(root)
    if model_source_path is not None:
        manifest["model_source_path"] = str(Path(model_source_path).expanduser().resolve())
    is_ml = manifest.get("finetune_method") == "machine_learning"
    if is_ml and trust_model is not True:
        raise ConfigurationError(
            "ML artifacts use joblib; set trust_model=True only for a trusted model source.",
            code="UNTRUSTED_MACHINE_LEARNING_ARTIFACT",
        )
    artifact = root / ("estimator.joblib" if is_ml else "head.npz")
    if not artifact.is_file():
        raise InputFormatError(
            "The trained property model artifact is missing; "
            "restore the training output directory.",
            code="INVALID_CUSTOM_MODEL_ARTIFACT",
            context={"path": str(artifact)},
        )
    checkpoint_status("Using local checkpoint", artifact)
    saved_features = manifest.get("features", {})
    if is_ml and not isinstance(saved_features, Mapping):
        raise InputFormatError(
            "Invalid saved feature schema.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    schema = _target_schema(manifest.get("target_schema"))
    stored_mapping = _mapping_from_payload(manifest.get("data_mapping"), include_labels=False)
    selected_mapping = mapping or (None if llm is not None else stored_mapping)
    inspection = inspect_property_data(
        source,
        mapping=selected_mapping,
        llm=llm,
        llm_client=llm_client,
        reader=data_reader,
    )
    if inspection.mapping.label_fields or inspection.mapping.split_field is not None:
        inspection = type(inspection)(
            inspection.source_path,
            inspection.detected_format,
            inspection.source_size_bytes,
            inspection.columns,
            inspection.estimated_records,
            _mapping_from_payload(inspection.mapping.to_dict(), include_labels=False),
            inspection.mapping_origin,
            inspection.recommendations,
            inspection.warnings,
            inspection.llm_model,
            inspection.shared_sample_records,
            inspection.shared_sample_characters,
        )
    data = load_property_data(
        source,
        inspection=inspection,
        reader=data_reader,
        require_labels=False,
        require_sequences=not is_ml
        or cast(Mapping[str, object], saved_features).get("source") != "numeric",
    )
    if is_ml:
        from .ml_training import predict_machine_learning_property

        prediction = predict_machine_learning_property(
            root,
            manifest,
            data,
            schema,
            representation_backend=representation_backend,
            batch_size=batch_size,
            device=device,
            threshold=threshold,
            allow_remote_code=allow_remote_code,
            show_progress=show_progress,
        )
        return _prediction_output(prediction, destination, overwrite=overwrite)
    torch = _torch()
    method = str(manifest.get("finetune_method"))
    model_name = str(manifest.get("model_name"))
    input_dimension = int(cast(Any, manifest.get("input_dimension")))
    selected_threshold = float(
        cast(Any, manifest.get("threshold", 0.5) if threshold is None else threshold)
    )
    if not 0 < selected_threshold < 1:
        raise ConfigurationError(
            "threshold must be in (0, 1).", code="INVALID_CUSTOM_PROPERTY_CONFIG"
        )
    representations = None
    encoder: TrainableSequenceEncoder | None = None
    if method == "linear_probe":
        checkpoint_value = manifest.get("checkpoint_path")
        checkpoint_path = str(checkpoint_value) if isinstance(checkpoint_value, str) else None
        representation_config = RepresentationConfig(
            model=model_name,
            checkpoint_path=(
                checkpoint_path if checkpoint_path and Path(checkpoint_path).is_dir() else None
            ),
            model_source_path=(
                str(manifest["model_source_path"])
                if isinstance(manifest.get("model_source_path"), str)
                else None
            ),
            pooling=str(manifest.get("pooling", "mean")),  # type: ignore[arg-type]
            dtype=str(manifest.get("dtype", "auto")),  # type: ignore[arg-type]
            max_length=cast(int | None, manifest.get("max_length")),
            batch_size=batch_size,
            device=device,
            max_records=max(1, data.record_count),
            allow_remote_code=allow_remote_code,
            show_progress=show_progress,
        )
        representation_result = extract_representations(
            _records(data),
            config=representation_config,
            backend=representation_backend,
        )
        representations = representation_result.representations
        runtime_device = _resolve_device(torch, device)
    elif method == "full":
        checkpoint_status("Using local checkpoint", root / "backbone")
        encoder = create_trainable_encoder(
            get_embedding_model(model_name),
            root / "backbone",
            device=device,
            dtype=str(manifest.get("dtype", "auto")),
            pooling=str(manifest.get("pooling", "mean")),
            max_length=cast(int | None, manifest.get("max_length")),
            allow_remote_code=allow_remote_code,
            model_source_path=cast(str | None, manifest.get("model_source_path")),
        )
        runtime_device = encoder.device
    else:
        raise InputFormatError(
            "The saved fine-tuning method is invalid.", code="INVALID_CUSTOM_MODEL_ARTIFACT"
        )
    head = _load_head(
        torch,
        root / "head.npz",
        input_dimension=input_dimension,
        output_dimension=schema.output_dimension,
        device=runtime_device,
    )
    head.eval()
    if encoder is not None:
        encoder.eval()
    output_batches: list[Any] = []
    with torch.no_grad():
        for batch in track(
            _batches(tuple(range(data.record_count)), batch_size=batch_size, shuffle=False, seed=0),
            description=f"Predicting {model_name} properties",
            total=(data.record_count + batch_size - 1) // batch_size,
            disable=not show_progress,
        ):
            features = _features(
                torch,
                batch,
                method=method,
                sequences=data.sequences,
                representations=representations,
                encoder=encoder,
                device=runtime_device,
            )
            output_batches.append(head(features).detach().float().cpu().numpy())
    logits = _numpy().concatenate(output_batches)
    prediction_rows, probability_rows = decoded_prediction_rows(
        data.record_ids,
        logits,
        schema=schema,
        threshold=selected_threshold,
    )
    prediction = CustomPropertyPredictionResult(
        str(root),
        data.record_ids,
        schema.label_names,
        tuple(FrozenDict(row) for row in prediction_rows),
        tuple(FrozenDict(row) for row in probability_rows),
    )
    return _prediction_output(prediction, destination, overwrite=overwrite)


__all__ = ["predict_custom_properties", "train_custom_property_model"]

_register_public_verbose(globals())
