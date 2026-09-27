# Custom Model Training and Prediction

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) Machine-learning classification and regression {#machine-learning-classification-and-regression}

- **Purpose:** Train a classifier or regressor from a labelled file and predict new samples.
- **Method:** Extract sequence features and fit a random forest.
- **API:** `dnakit.predictions.train_custom_property_model(source, config=...)`; `dnakit.predictions.predict_custom_properties(source, model_dir=...)`
- **Input:** CSV columns `sequence` and `target`; this example uses regression. Install `.[ml]` first. Only load your own or trusted models with `trust_model=True`.
- **Output:** Returns `CustomPropertyTrainingResult` with training artifacts and evaluation metrics.
- **Example code:**

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

result = train_custom_property_model(
    "my_labels.csv",
    config=CustomPropertyTrainingConfig(
        task_type="regression",
        model="random_forest",
        finetune_method="machine_learning",
        feature_source="sequence",
        mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
        output_dir="runs/random-forest",
        workers=0,
        seed=42,
    ),
)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir, trust_model=True,
)
print(predictions.predictions)
```

- **Example output:** Test metrics and per-record predictions; values depend on your data and training. Model files are saved under `result.output_dir`.

## 2) Custom deep-learning classification and regression {#custom-deep-learning-classification-and-regression}

- **Purpose:** Train a classifier or regressor from a labelled file and predict new samples.
- **Method:** All 11 registered models provide full fine-tuning (`full`) and frozen-backbone linear probing (`linear_probe`); this example uses DNABERT-2. See [models and fine-tuning methods](../../faq/tutorial_integrations.en.md#models-and-fine-tuning-methods) for dependencies, representation choices, and validation scope.
- **API:** `dnakit.predictions.train_custom_property_model(source, config=...)`; `dnakit.predictions.predict_custom_properties(source, model_dir=...)`
- **Input:** CSV columns `sequence` and `target`; this example uses regression. Install `.[training]` first. The first run may download weights; enable remote model code only for trusted sources.
- **Output:** Returns `CustomPropertyTrainingResult` with training artifacts and evaluation metrics.
- **Example code:**

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

result = train_custom_property_model(
    "my_labels.csv",
    config=CustomPropertyTrainingConfig(
        task_type="regression",
        model="dnabert2",
        finetune_method="full",
        allow_remote_code=True,
        mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
        output_dir="runs/my-property",
        workers=0,
        seed=42,
    ),
)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir, allow_remote_code=True,
)
print(predictions.predictions)
```

- **Example output:** Test metrics and per-record predictions; values depend on your data and training. Model files are saved under `result.output_dir`.

<span id="checkpoint-cache"></span>Checkpoint locations and automatic downloads: [FAQ](../../faq/tutorial_integrations.md#property-prediction).

## 3) Full fine-tuning settings for NT tasks

This configuration uses AdamW, ReduceLROnPlateau on validation MCC, and early stopping after five epochs without improvement:

```python
from dnakit.predictions import CustomPropertyTrainingConfig

config = CustomPropertyTrainingConfig(
    task_type="classification",
    model="enformer",
    finetune_method="full",
    learning_rate=1e-4,
    batch_size=64,
    beta1=0.9,
    beta2=0.95,
    weight_decay=0.1,
    gradient_accumulation_steps=1,
    lr_scheduler="reduce_on_plateau",
    scheduler_patience=1,
    scheduler_factor=0.5,
    monitor="mcc",
    monitor_mode="auto",
    early_stopping_patience=5,
    device="cuda:0",
    output_dir="runs/enformer-nt",
)
# Pass config to train_custom_property_model(source, config=config).
```

`batch_size=64` means 64 samples per batch on the selected device, with one optimizer update per batch when `gradient_accumulation_steps=1`. This trainer runs on one device; it does not launch distributed training across GPUs. The generic defaults remain batch size 32, no scheduler, and selection by validation loss.

Supported monitors are `mcc`, `f1_macro`, and `accuracy` for classification; `rmse` and `mae` for regression; and `validation_loss` for either task. Auto mode maximizes classification scores and minimizes losses/errors. Explicit `min` and `max` overrides are also available. MCC monitoring supports binary, multiclass, and multioutput classification; use `f1_macro` for multilabel targets. Only validation results drive these decisions.

The scheduler steps after each epoch's validation. Following [PyTorch's patience semantics](https://docs.pytorch.org/docs/stable/generated/torch.optim.lr_scheduler.ReduceLROnPlateau.html), patience 1 tolerates one bad epoch and multiplies the next epoch's learning rate by 0.5 after the second consecutive bad epoch. Early stopping has its own counter, which a learning-rate reduction does not reset. Scheduling, early stopping, and best-checkpoint selection share an absolute `monitor_min_delta` (default `1e-12`, scheduler `threshold_mode="abs"`), including for zero or negative MCC. Training restores the weights with the best monitored score.

`history.json` records `monitor_value`, `primary_metric_name`, the current epoch's `learning_rate`, and the scheduled `next_learning_rate`. The `training` section in `model.json` records the settings. `epochs` sets the training limit. The same options are exposed by `dnakit train-property`, including `--monitor`, `--lr-scheduler`, `--scheduler-patience`, `--scheduler-factor`, and `--early-stopping-patience`.

## 4) Fold training and refitting on all training data

Neural `full` and `linear_probe` training accept optional validation and test partitions. With `DataMapping(split_field="split")`, a source containing only `train` and `valid` rows trains a fold without creating a test partition. A supplied validation partition produces `validation_predictions.jsonl` with record IDs, predictions and class probabilities. `test_predictions.jsonl` is only produced when test rows exist.

Set `log_train_metrics=True` to evaluate training MCC after every epoch; this adds one forward pass over the training set. `history.json` contains `train_mcc` and `validation_mcc`, with `null` for unavailable or unevaluated values. Training loss is the sample-weighted mean of the training batches; validation loss and both MCC values are evaluated after the epoch.

Once validation has selected a configuration and epoch count, refit on all official training data:

```python
from dnakit.predictions import CustomPropertyTrainingConfig, train_custom_property_model

config = CustomPropertyTrainingConfig(
    task_type="classification", model="enformer", finetune_method="full",
    checkpoint_path="models/enformer", output_dir="runs/enformer-final",
    train_ratio=1.0, validation_ratio=0.0, test_ratio=0.0,
    epochs=8,  # Fixed beforehand using internal validation only.
    learning_rate=1e-4, batch_size=64, beta1=0.9, beta2=0.95, weight_decay=0.1,
    gradient_accumulation_steps=1,
    early_stopping_patience=0, lr_scheduler="none", log_train_metrics=True,
)
result = train_custom_property_model("official_train.csv", config=config)
```

Without validation, disable early stopping and validation-driven scheduling. Training runs for the fixed epoch count and saves the final epoch; logs use training loss to describe progress. Keep official test data separate and infer with `predict_custom_properties` only after model selection. The manifest records the actual split and configuration; returned metrics omit absent partitions. Machine-learning training without CV still requires three nonempty partitions. Neural cross-validation requires an explicit user-side fold loop; `cv_folds` is reserved for machine-learning mode.
