# Property prediction: deep learning and machine learning

The **14 added trained-model interfaces**, including AlphaGenome, ESM, Evo 2 and NVIDIA NIM, have individual parameters, returns and registered credentials in [trained model predictions](../native/models.md). CADD, AlphaMissense and SpliceAI annotation/scoring interfaces are in [variant analysis](../native/variants.md). The sections below cover DNAKit's training and prediction workflows.

If a prediction dependency is missing or cannot be imported, DNAKit prints installation guidance to stderr before raising the existing exception. The exception also retains the installation commands. Dependency error diagnostics remain visible with `verbose=False`; DNAKit does not install packages automatically.

From the **DNAKit project root**, install the extras for your workflow:

| Workflow | Installation command |
| --- | --- |
| Enformer, including its 27 classification tasks | `python -m pip install -e ".[neural,neural-enformer]"` |
| General PyTorch / Transformers backbones | `python -m pip install -e ".[neural]"` |
| Custom deep-learning training | `python -m pip install -e ".[training]"` |
| Basic machine learning | `python -m pip install -e ".[ml]"` |
| XGBoost / LightGBM / CatBoost | Use `.[xgboost]` / `.[lightgbm]` / `.[catboost]`, respectively |

For a published package, replace `-e ".[extras]"` with `"dnakit[extras]"`; no source directory is needed. Evo 2, Caduceus and source-based backends receive their own guidance. The Evo 2 extra requires Python 3.11–3.12; AlphaGenome, LucaOneTasks and JanusDNA also require their official environments.

Models can share an installed PyTorch or TensorFlow version when it satisfies all dependency constraints; the two frameworks are not interchangeable. These commands do not force upgrades or reinstallations of compatible packages. Missing or incompatible dependencies are resolved by pip. The current `neural` / `training` extras specify a PyTorch range of `>=2.5,<3`, rather than different exact versions for every model; third-party backends may impose tighter constraints. See [pip installation behavior](https://pip.pypa.io/en/stable/cli/pip_install/).

This page covers three paths:

1. **User-defined properties:** train classification or regression models from your own files, with full fine-tuning or linear probing.
2. **Official pretrained properties:** run inference with released task heads without retraining locally.
3. **Machine-learning properties:** train your own classification or regression model with scikit-learn or optional XGBoost, LightGBM, or CatBoost.

The existing [neural representations](08_fingerprints.md) API remains embedding-only. Custom training is introduced first; numbered sections 1–54 describe the 54 official direct-inference functions.

Each function below describes its purpose, prediction method, API, input, and output. Placeholders in example results indicate output structure only. Actual values, label counts, and track counts depend on the input, species, tissue filters, and checkpoint.

Model papers, official repositories, and weight sources are collected in [FAQ: references for deep-learning property prediction](../../faq.md#deep-learning-property-prediction-references).

## Input and output error handling

Direct prediction rejects empty sequences and invalid symbols with `record_id` in the sequence
error. Enformer benchmark inputs may be 1, 199, 200 or 201 bp; the default upper bound is
196608 bp. `max_length` may tighten but cannot relax that architecture limit. Inputs are never
silently truncated, and IUPAC symbols follow `ambiguity_policy`.

`PredictionOutput` rejects empty, non-finite, complex or unconvertible values and duplicate
output names. Direct classification/segmentation outputs must be probabilities in `[0,1]`.
Enformer classification additionally requires the class dimension to match and probabilities
to sum to one within `1e-5`. Regression, scores and tracks retain their own numeric semantics.
Custom backends must follow the selected task's output contract.

Direct prediction and embedding APIs report memory exhaustion, including wrapped backend OOM
errors, as `MODEL_OUT_OF_MEMORY`, with batch/device guidance instead of substitute predictions.
An unavailable explicit CUDA device or device index is rejected before model loading.

## Checkpoint reuse and automatic downloads {#checkpoint-cache}

Pretrained weights are checked locally before prediction: missing weights are downloaded only for the selected model/task, with start/completion messages and absolute local paths. Existing weights are reused. Set `verbose=False` to silence feedback. Enformer task weights default to `./ckpt/enformer-benchmarks/`; other backbones use `./ckpt/<model>/`. ML embedding backbones follow the same policy. User-trained estimators (`estimator.joblib`) and neural heads (`head.npz`) must be retained in `model_dir`; they have no registered public download source. LucaOneTasks uses its official downloader in the supplied source checkout.

For the 27 Enformer classification tasks, no checkpoint path is required. Use `checkpoint_dir` to select a cache directory: existing files are also checked in its `enformer-benchmarks/` and task-family (`nt-revised/` or `genomic-benchmarks/`) subdirectories; a missing file is downloaded directly into the selected directory. Alternatively, `checkpoint_path` selects an exact file, downloading the selected Enformer task there if absent. The two parameters are mutually exclusive.

Install with `python -m pip install "dnakit[neural,neural-enformer]"`. Google Drive downloads use `gdown`, included in these extras. Downloads are written to temporary files and checked before installation; failed or interrupted files are not reused as checkpoints. Prediction validates the task and model parameters. Read class meanings from `output.output_names`.

## Machine-learning classification and regression

Use the existing `CustomPropertyTrainingConfig`, `train_custom_property_model()` and `predict_custom_properties()` APIs with `finetune_method="machine_learning"`. Here `model` selects the estimator family; `embedding_model` separately selects an optional frozen DNA backbone.

Install ML and tabular readers without PyTorch:

```bash
python -m pip install -e ".[ml]"
```

Install the three optional boosting libraries:

```bash
python -m pip install -e ".[xgboost,lightgbm,catboost]"
```

Embedding inputs additionally require `.[neural]` and the selected backbone's extras. Missing optional libraries affect only their own algorithms and produce an installation hint.

### Available algorithms

`available_machine_learning_models()` lists all 37 families, supported tasks and dependencies without importing optional libraries.

| Tasks | `model` names |
| --- | --- |
| Classification and regression | `ridge`, `sgd`, `svm`, `linear_svm`, `nu_svm`, `knn`, `decision_tree`, `extra_tree`, `random_forest`, `extra_trees`, `gradient_boosting`, `hist_gradient_boosting`, `adaboost`, `bagging`, `gaussian_process`, `mlp`, `dummy`, `xgboost`, `lightgbm`, `catboost` |
| Classification only | `logistic_regression`, `perceptron`, `gaussian_nb`, `multinomial_nb`, `bernoulli_nb`, `complement_nb`, `lda`, `qda` |
| Regression only | `linear_regression`, `lasso`, `elastic_net`, `bayesian_ridge`, `ard`, `huber`, `ransac`, `theil_sen`, `kernel_ridge` |

`mlp` is scikit-learn's shallow neural network; `dummy` is a baseline. This interface covers supervised prediction, not clustering, anomaly detection or automatic model search. Individual estimators retain their sample-size, parameter and computational constraints; Gaussian processes, for example, are costly on large datasets.

### Feature sources

| `feature_source` | Behavior |
| --- | --- |
| `sequence` (default) | Sequence length, GC fraction, ambiguous-base fraction and fixed-order k-mer frequencies. `kmer_k=3`, allowed range 1–6; dimension is `3 + 4**kmer_k`. Ambiguous windows are skipped; GC uses only A/C/G/T in its denominator. Undefined frequencies become zero. No deep-learning checkpoint is required. |
| `embedding` | Frozen vectors from `embedding_model` (default `dnabert2`). Checkpoint, pooling, device and batching settings follow the representation API. |
| `numeric` | Explicit, ordered `DataMapping.feature_fields`. Set `sequence_field=None` for files without sequences. Missing, nonnumeric or infinite features fail validation; there is no automatic imputation. |

Numeric columns must be explicitly selected and must not overlap labels, ID, split or sequence fields. Prediction requires the same feature names and order. Avoid target-derived features or other leakage.

Train a random forest from sequence data:

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

config = CustomPropertyTrainingConfig(
    task_type="regression",  # or classification
    finetune_method="machine_learning",
    model="random_forest",
    model_params={"n_estimators": 300, "max_depth": 12},
    feature_source="sequence",
    kmer_k=3,
    mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
    output_dir="runs/random-forest",
    workers=6,
    seed=42,
)
result = train_custom_property_model("my_labels.csv", config=config)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir,
    trust_model=True,  # only for artifacts from a trusted source
)
print(predictions.predictions)
```

For numeric input use `feature_source="numeric"` and, for example, `DataMapping(sequence_field=None, feature_fields=("gc", "length"), label_fields=("target",))`. For embeddings use `feature_source="embedding", embedding_model="dnabert2"` and the checkpoint's required loading options, including `allow_remote_code` where applicable.

### Training and persistence

Defaults remain seed 42 and an 80%/10%/10% split; a mapped split column takes precedence. Target semantics and metrics match the deep-learning path. Multilabel, multioutput classification and multi-target regression fit an independent estimator per target. Constant binary targets use constant classifiers. Class predictions use the highest probability; multilabel predictions use `threshold`.

`model_params` supplies serializable estimator constructor parameters, not arbitrary estimator objects. `workers` controls default `n_jobs` and numerical-library threads (0 means 1); CatBoost uses `thread_count`. Explicit estimator parameters can override their parallelism and seeds. `standardize_features=True` fits scaling only on training rows. Nonnegative-only `multinomial_nb` and `complement_nb` use MaxAbsScaler; negative inputs are rejected instead of silently shifted.

Classifiers without `predict_proba` use sigmoid calibration with up to three stratified folds inside the training split. Each class needs at least two training samples. Scaling is also fitted inside each calibration fold. Validation/test rows are excluded, and calibration is recorded in result warnings. See the [scikit-learn calibration API](https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html).

The ML path does not use neural `epochs`, AdamW settings or gradient accumulation. Set iterations/tree counts through `model_params`. Validation is evaluation-only; there is no automatic tuning or best-model selection.

Regression targets use training-set mean/std normalization, as in the existing API; reports and predictions return original units. Classification `class_weight` / `class_weights` dictionaries use zero-based encoded class indices, ordered as saved in `model.json` under `target_schema.class_names`. Options such as `"balanced"` remain available where the estimator supports them.

Artifacts are `model.json`, `estimator.joblib`, `metrics.json`, `data_inspection.json`, `splits.csv` and `test_predictions.jsonl`. Fitted transforms are stored with estimators; the manifest records feature order, target schema, settings and Python/library versions. The returned `history` is empty because there is no neural training history.

Loading joblib can execute Python code, so ML inference requires explicit `trust_model=True` for artifacts you trained or otherwise trust. Reloading checks the recorded Python and dependency versions and requires the original environment if they differ. See [scikit-learn persistence guidance](https://scikit-learn.org/stable/model_persistence.html).

Edit `others/examples/train_machine_learning_properties.py` to use your own data, or use `--demo` for synthetic training and reload verification:

```bash
python others/examples/train_machine_learning_properties.py --demo
```

## User-defined deep-learning classification and regression

Install the training and common tabular-format dependencies first: `python -m pip install "dnakit[training]"`.

`train_custom_property_model()` accepts one local data file and maps an ID, a DNA/RNA sequence, one or more targets, and an optional existing split. The default train/validation/test ratios are 0.8/0.1/0.1, with seed 42. A mapped existing split column takes precedence. The caller must explicitly set `task_type` to `classification` or `regression`; DNAKit does not infer the scientific task merely from how values look.

Built-in readers cover CSV, TSV, other single-character delimiters, JSON, JSONL, YAML, XML, Excel, Parquet, SQLite, FASTA, FASTQ, and GenBank. Dot paths address nested records, and `key=value` FASTA/FASTQ description tokens become `header.key`. A caller can implement the `DataReader` protocol for any other format. An LLM may propose only a validated declarative `DataMapping`; DNAKit never executes model-generated Python, SQL, regular expressions, or shell commands.

### Label semantics

- One categorical label column: binary or multiclass classification.
- Multiple binary label columns: automatically detected as one-hot multiclass or multilabel classification; `classification_mode` can select the interpretation explicitly.
- Multiple nonbinary categorical columns: multioutput classification, with an independent class table for each column.
- One or more numeric label columns: single-target or multi-target regression.

Classification output includes accuracy, balanced accuracy, precision, recall, F1, Jaccard, MCC, Cohen's kappa, log loss, AUROC, average precision, confusion matrices, and per-class/per-label results. Multilabel output also includes subset accuracy, Hamming loss, LRAP, and coverage error. Regression output includes MAE, MSE, RMSE, R², explained variance, median absolute error, max error, bias, MAPE, SMAPE, Pearson, and Spearman metrics, both aggregate and per target. Undefined metrics are `null`; classification metrics such as F1 are not reported for regression.

### Models and fine-tuning methods

The default is `finetune_method="full"`. DNABERT-2, NT-v2, HyenaDNA, Caduceus, GROVER, LucaOne, and GENERator currently expose the portable Transformers full-training path. Every registered representation model supports a frozen `linear_probe`, including AlphaGenome, Enformer, Evo 2, and JanusDNA. Query `available_custom_training_models()` for the version-specific support matrix. Models that load checkpoint-provided code still require `allow_remote_code=True`.

### Training parameters (Python API)

Configure training directly with `CustomPropertyTrainingConfig`; no command-line arguments are needed. Full fine-tuning and linear probing share these defaults:

| Python parameter | Default | Meaning |
| --- | --- | --- |
| `workers` | `6` | CPU data-loading processes; `0` loads in the main process |
| `batch_size` | `32` | Samples per training device per batch, before gradient accumulation |
| `learning_rate` | `1e-4` | AdamW learning rate (lr) |
| `beta1` | `0.9` | AdamW first-moment decay |
| `beta2` | `0.95` | AdamW second-moment decay |
| `weight_decay` | `0.1` | AdamW weight decay |

`workers` fetch records and assemble CPU batches for training, validation, and testing; encoder/GPU operations remain in the main process. This setting does not parallelize initial file parsing or LLM requests. Changing the worker count does not change the seeded split or epoch sample order. The trainer currently runs on one selected `device` (for example `"cuda:0"`), not multi-GPU DDP; available GPU count does not multiply the batch size. `gradient_accumulation_steps` defaults to `1` and is also configurable through Python.

Scripts using `workers>0` need the `if __name__ == "__main__":` guard below. Use `workers=0` in notebooks or environments without subprocess support. Worker processes are released on completion, early stopping, or errors.

```python
import os

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    LLMSchemaConfig,
    train_custom_property_model,
)

def main():
    result = train_custom_property_model(
        "my_labels.xlsx",
        config=CustomPropertyTrainingConfig(
            task_type="regression",  # or "classification"
            model="dnabert2",
            finetune_method="full",  # or "linear_probe"
            output_dir="runs/my-property",
            workers=6,
            batch_size=32,
            learning_rate=1e-4,
            beta1=0.9,
            beta2=0.95,
            weight_decay=0.1,
            seed=42,
            train_ratio=0.8,
            validation_ratio=0.1,
            test_ratio=0.1,
            allow_remote_code=True,
            llm=LLMSchemaConfig(
                provider="deepseek",
                model="deepseek-v4-flash",
                api_key=os.environ["DEEPSEEK_API_KEY"],
            ),  # Omit llm when schema assistance is unnecessary.
        ),
    )
    print(result.metrics["test"])


if __name__ == "__main__":
    main()
```

The output directory contains `model.json`, `head.npz`, `metrics.json`, `history.json`, `data_inspection.json`, `splits.csv`, and `test_predictions.jsonl`; full fine-tuning also writes `backbone/`. The `training` section of `model.json` records the worker count, batch size, learning rate, both betas, and weight decay. `predict_custom_properties()` reloads those artifacts for later inference. The editable script `others/examples/train_custom_properties.py` can be run from an IDE after configuring the source file; LLM assistance is off by default in that script.

### LLM file analysis and privacy boundary

`LLMSchemaConfig(provider=..., model=...)` supports the following connection presets through OpenAI-compatible `POST /chat/completions`, without vendor SDK dependencies. Query `available_schema_llm_providers()` offline to list them. Except for the legacy DeepSeek default, explicitly select a model ID available in your account; DNAKit does not automatically substitute a model or provider.

| `provider` | Default `base_url` | API-key environment variable |
| --- | --- | --- |
| `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| `custom` | Required explicitly | `LLM_API_KEY` |

For example, use Google Gemini to analyze only field names and type structure:

```python
from dnakit.predictions import LLMSchemaConfig, inspect_property_data

llm = LLMSchemaConfig(
    provider="gemini",
    model="gemini-3.8-flash",  # Replace with an available text/JSON model if needed.
    include_sample_values=False,
)
inspection = inspect_property_data("my_labels.xlsx", llm=llm)
print(inspection.mapping.to_dict())
print(inspection.recommendations)
# Pass the same llm into CustomPropertyTrainingConfig when training.
```

Override `api_key_env` for your own environment-variable name or pass `api_key` in Python without committing it to source control. Explicit `base_url` overrides support regional or private gateways. HTTP redirects are refused so that keys and sampled profiles are not forwarded to another endpoint. Gemini uses its [official compatibility endpoint](https://ai.google.dev/gemini-api/docs/openai), not native `generateContent` or Vertex AI. Qwen defaults to the still-supported Beijing endpoint; a workspace-specific endpoint such as `https://<workspaceID>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` can also be used. Endpoints and API keys must match the region, as described in [Alibaba Cloud's documentation](https://help.aliyun.com/zh/model-studio/compatibility-of-openai-with-dashscope). For Kimi's international endpoint, set `base_url="https://api.moonshot.ai/v1"`, per its [quickstart](https://platform.kimi.ai/docs/overview). OpenRouter requires its complete model ID; check the chosen model's supported parameters and downstream provider/data policies in the [quickstart](https://openrouter.ai/docs/quickstart) and [structured-output documentation](https://openrouter.ai/docs/guides/features/structured-outputs).

Optional Python settings:

- `temperature=None` (default) omits temperature for models that reject custom sampling parameters.
- `token_limit_parameter="auto"` sends `max_completion_tokens` for OpenAI and `max_tokens` for other compatible services; either field can also be selected explicitly. This follows [OpenAI Docs](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create), which notes that the older `max_tokens` field is incompatible with some reasoning models.
- `json_mode=True` requests JSON Object output. Set it to `False` for a model without this parameter; the prompt and mapping validation still require a JSON object.
- `extra_body` passes model-specific fields, such as `{"enable_thinking": False}` for supported Qwen models or `{"thinking": {"type": "disabled"}}` for supported Kimi models. Check the [Qwen](https://help.aliyun.com/zh/model-studio/qwen-structured-output) or [Kimi parameter reference](https://platform.kimi.ai/docs/api/models-overview) before applying these options; DNAKit does not force thinking off for all models. Requests are non-streaming, so models requiring streamed thinking must disable thinking or be replaced with a compatible model.

`extra_body` cannot replace authentication, messages, model selection, token limits, or other core request fields. Truncated responses fail explicitly: adjust `max_output_tokens` or thinking settings yourself; DNAKit never increases the budget or retries another provider automatically. Presets describe interface support, not universal model compatibility. Provider contracts are covered by mocked tests, not paid live tests against every service.

`LLMSchemaConfig` uses OpenAI-compatible `POST /chat/completions` with JSON Output, including DeepSeek and compatible providers. By default, it shares at most 8 sampled records, 12,000 characters, and a structural profile—not the whole source. Even a file with fewer than 8 records retains at least one complete unshared record. Set `include_sample_values=False` to share only field/type structure. `DataInspection` records how many sample rows and characters were sent. API keys are excluded from request bodies, results, exception displays, and training artifacts; the Python example reads the key from an environment variable.

## 1) RNA-seq coverage and expression

- **Purpose:** Predict RNA-seq coverage and gene-expression tracks from DNA sequences.
- **Prediction method:** Run direct inference with the official AlphaGenome all-folds weights and the `RNA_SEQ` output. Short sequences are centered and padded with `N` to an officially supported context length; use `ontology_terms` to restrict tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="rna_seq")[required], backend[optional])`
- **Input:** `BiologicalSequence`; alternatively use `VariantContext` with `predict_variant_effects()` to obtain REF, ALT, and ALT−REF tracks. Human and mouse are supported.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="rna_seq",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of RNA-seq tracks>)
('position', 'track')
```

## 2) CAGE transcription initiation and expression

- **Purpose:** Predict CAGE signals describing capped-transcript initiation sites and expression levels.
- **Prediction method:** The official AlphaGenome `CAGE` head directly produces a position × track matrix; `ontology_terms` filters target tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="cage")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a single-SNV `VariantContext`. Sequences may be up to 1,048,576 bp; human and mouse are supported.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="cage",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of CAGE tracks>)
('position', 'track')
```

## 3) PRO-cap nascent transcription initiation

- **Purpose:** Predict PRO-cap nascent RNA 5′-end signals to locate active transcription initiation.
- **Prediction method:** Run the official AlphaGenome `PROCAP` output head directly. Output resolution is recorded in `metadata["resolution_bp"]`.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="procap")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Set `organism` and `ontology_terms` as needed.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="procap",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of PRO-cap tracks>)
<bp per output position>
```

## 4) ATAC-seq chromatin accessibility

- **Purpose:** Predict ATAC-seq open-chromatin tracks to identify potentially accessible regulatory regions.
- **Prediction method:** Run the official AlphaGenome `ATAC` output head directly, returning accessibility signals by position for different tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="atac")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a `VariantContext` containing exactly one base change. Use `ontology_terms` to restrict tracks.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="atac",
        ontology_terms=("UBERON:0002048",),
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of filtered ATAC-seq tracks>)
(<first two official track names>,)
```

## 5) DNase-seq chromatin accessibility

- **Purpose:** Predict DNase-seq hypersensitivity signals describing chromatin accessibility.
- **Prediction method:** Run the official AlphaGenome `DNASE` output head directly. Short inputs are centered and padded; padding lengths are recorded in result metadata.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="dnase")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Select `organism="human"` or `"mouse"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="dnase",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["input_length"])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of DNase-seq tracks>)
16384
```

## 6) Histone modifications

- **Purpose:** Predict histone ChIP-seq tracks describing histone modification signals across a sequence region.
- **Prediction method:** Run the official AlphaGenome `CHIP_HISTONE` output head directly, returning a position × histone-track matrix.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="chip_histone")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Select outputs by species and ontology term.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_histone",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of histone ChIP-seq tracks>)
('position', 'track')
```

## 7) Transcription factor binding

- **Purpose:** Predict transcription factor ChIP-seq binding tracks to locate potential binding signals.
- **Prediction method:** Run the official AlphaGenome `CHIP_TF` output head directly, returning a position × transcription-factor-track matrix.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="chip_tf")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Use `ontology_terms` to restrict tissues or cell types.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_tf",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **Example result (schematic structure):**

```text
(<number of positions or bins>, <number of TF ChIP-seq tracks>)
(<first two official track names>,)
```

## 8) Splice donor/acceptor sites

- **Purpose:** Predict splice donor and acceptor tracks at each sequence position.
- **Prediction method:** Run the official AlphaGenome `SPLICE_SITES` output head directly, retaining output resolution and track names.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_sites")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`; variant results additionally include an ALT−REF difference axis.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_sites",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example result (schematic structure):**

```text
(<number of positions>, <number of donor/acceptor tracks>)
('position', 'track')
```

## 9) Splice-site usage

- **Purpose:** Predict usage signals or usage-proportion tracks for candidate splice sites.
- **Prediction method:** Run the official AlphaGenome `SPLICE_SITE_USAGE` output head directly. `output_names` identifies the tracks.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_site_usage")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Specify the species, ontology terms, and official context length.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_site_usage",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **Example result (schematic structure):**

```text
(<number of positions>, <number of splice-site usage tracks>)
<bp per output position>
```

## 10) Splice junction counts

- **Purpose:** Predict connection-count tracks for candidate splice junctions.
- **Prediction method:** Run the official AlphaGenome `SPLICE_JUNCTIONS` output head directly. `junctions` in result metadata contains candidate junction identifiers.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_junctions")[required], backend[optional])`
- **Input:** Only `BiologicalSequence` is accepted. Candidate junction sets may change between alleles, so this task currently rejects `VariantContext`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_junctions",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
print(len(output.metadata.get("junctions", ())))
```

- **Example result (schematic structure):**

```text
(<number of candidate splice junctions>, <number of tracks>)
('junction', 'track')
<number of candidate splice junctions>
```

## 11) 3D DNA contact maps

- **Purpose:** Predict 3D DNA contact probability maps between chromatin regions.
- **Prediction method:** Run the official AlphaGenome `CONTACT_MAPS` output head directly, producing two position axes and one track axis.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="contact_maps")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`; variant results add a leading REF, ALT, ALT−REF axis.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="contact_maps",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example result (schematic structure):**

```text
(<number of bins>, <number of bins>, <number of contact-map tracks>)
('position_1', 'position_2', 'track')
```

## 12) Human regulatory tracks

- **Purpose:** Predict 5,313 human CAGE, DNase/ATAC, transcription factor ChIP, and histone ChIP tracks.
- **Prediction method:** Run public pretrained Enformer weights directly. By default, center and pad the input with `N` to 196,608 bp and return 896 bins of 128 bp.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="enformer", task="human_tracks")[required], backend[optional])`
- **Input:** `BiologicalSequence`; alternatively use `VariantContext` with `predict_variant_effects()` for REF, ALT, and ALT−REF. `max_length` must be a multiple of 128 in the range 114,688–393,216.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="human_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["bin_size_bp"])
```

- **Example result (default context):**

```text
(896, 5313)
128
```

## 13) Mouse regulatory tracks

- **Purpose:** Predict 1,643 mouse CAGE, DNase/ATAC, transcription factor ChIP, and histone ChIP tracks.
- **Prediction method:** Run Enformer’s mouse output head directly with a default context of 196,608 bp and tracks binned at 128 bp.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="enformer", task="mouse_tracks")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a single-SNV `VariantContext`. Short inputs are automatically centered with `N` padding at both ends.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("mouse-region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="mouse_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["organism"])
```

- **Example result (default context):**

```text
(896, 1643)
mouse
```

## 14) 14-class single-base genomic segmentation

- **Purpose:** Predict probabilities for 14 gene and regulatory-element classes at each base of human DNA.
- **Prediction method:** Use `InstaDeepAI/segment_nt`, an NT-v2 backbone with a trained U-Net segmentation head. Apply softmax to binary logits and return the probability that each feature is present.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation", allow_remote_code=True)[required], backend[optional])`
- **Input:** Human `BiologicalSequence` up to 30,000 bp. Sections 14.1–14.14 list the 14 classes.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 1000)],
    config=PropertyPredictionConfig(
        model="segmentnt",
        task="genomic_segmentation",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(len(output.output_names))
```

- **Example result:**

```text
(4000, 14)
14
```

The following 14 classes share this API. Each corresponds to a column in `output.values`; use `output.output_names` for the actual column names and order.

### 14.1) Protein-coding gene

Predict the probability that each base belongs to a protein-coding gene region.

### 14.2) Long noncoding RNA (lncRNA)

Predict the probability that each base belongs to a long noncoding rna (lncrna) region.

### 14.3) Exon

Predict the probability that each base belongs to a exon region.

### 14.4) Intron

Predict the probability that each base belongs to a intron region.

### 14.5) Splice donor

Predict the probability that each base belongs to a splice donor region.

### 14.6) Splice acceptor

Predict the probability that each base belongs to a splice acceptor region.

### 14.7) 5′ untranslated region (5′ UTR)

Predict the probability that each base belongs to a 5′ untranslated region (5′ utr) region.

### 14.8) 3′ untranslated region (3′ UTR)

Predict the probability that each base belongs to a 3′ untranslated region (3′ utr) region.

### 14.9) CTCF-bound region

Predict the probability that each base belongs to a ctcf-bound region region.

### 14.10) Polyadenylation signal (polyA signal)

Predict the probability that each base belongs to a polyadenylation signal (polya signal) region.

### 14.11) Tissue-specific enhancer

Predict the probability that each base belongs to a tissue-specific enhancer region.

### 14.12) Tissue-invariant enhancer

Predict the probability that each base belongs to a tissue-invariant enhancer region.

### 14.13) Tissue-specific promoter

Predict the probability that each base belongs to a tissue-specific promoter region.

### 14.14) Tissue-invariant promoter

Predict the probability that each base belongs to a tissue-invariant promoter region.

## 15) Long-context zero-shot variant effects

- **Purpose:** Calculate REF and ALT sequence likelihoods and an ALT−REF effect score for a single SNV.
- **Prediction method:** Use the official Evo 2 7B base model to calculate mean sequence log-likelihood over the complete REF and ALT contexts, then compute `alternate − reference`. This is a zero-shot ranking score, not a clinically calibrated pathogenicity probability.
- **API:** `dnakit.predictions.predict_variant_effects(inputs[required], config=PropertyPredictionConfig(model="evo2", task="variant_effect")[required], backend[optional])`
- **Input:** `VariantContext`; REF and ALT must have equal lengths and differ at exactly one A/C/G/T position.
- **Example code:**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="A" * 4096 + "C" + "G" * 4095,
    alternate_sequence="A" * 4096 + "T" + "G" * 4095,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="evo2",
        task="variant_effect",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example result:**

```text
('reference_likelihood', 'alternate_likelihood', 'alternate_minus_reference')
(3,)
```

## 16) Exon probability

- **Purpose:** Predict whether a candidate site belongs to an exon using forward and reverse contexts around it.
- **Prediction method:** Read 4,096-dimensional last-token embeddings for both contexts from Evo 2 `blocks.26`, concatenate them into 8,192 dimensions, and apply the officially trained MLP head. The classifier is not retrained.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="evo2", task="exon_probability", allow_remote_code=True)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, ordered as forward context then reverse context, both with `sequence_type="gene"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

contexts = BiologicalSequencePair(
    "site-1",
    BiologicalSequence("forward-1", "ACGT" * 256),
    BiologicalSequence("reverse-1", "ACGT" * 256),
)
result = predict_pair_properties(
    [contexts],
    config=PropertyPredictionConfig(
        model="evo2",
        task="exon_probability",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example result (label depends on probability and threshold):**

```text
('exon_probability',)
exon or non_exon
```

## 17) Allele-conditional variant effects

- **Purpose:** Compare next-token support for REF and ALT bases from the sequence upstream of a variant.
- **Prediction method:** Use the official GENERator v2 eukaryote 1.2B weights with up to 8,192 bp of upstream context. Remove leading `N`, left-truncate to a multiple of 6, and calculate `log(p_ref/(p_alt+1e-10))`. This is a zero-shot ranking score, not a clinically calibrated probability.
- **API:** `dnakit.predictions.predict_variant_effects(inputs[required], config=PropertyPredictionConfig(model="generator", task="variant_effect", allow_remote_code=True)[required], backend[optional])`
- **Input:** `VariantContext`; REF and ALT must have equal lengths and exactly one SNV. Only the reference context before the variant is used.
- **Example code:**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="ACGT" * 1024 + "C" + "G" * 255,
    alternate_sequence="ACGT" * 1024 + "T" + "G" * 255,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="generator",
        task="variant_effect",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example result:**

```text
('reference_probability', 'alternate_probability', 'log_ratio_score')
(3,)
```

## 18) Nucleic acid–protein central dogma relationship

- **Purpose:** Predict whether a gene sequence and a protein sequence have a CentralDogma relationship.
- **Prediction method:** Run the official LucaOneTasks `CentralDogma` binary-classification checkpoint through `src/predict_v1.py`, returning the positive-class probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="central_dogma", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, with `sequence_type="gene"` first and `sequence_type="protein"` second. `model_source_path` must point to the official LucaOneTasks source directory.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "central-dogma-1",
    BiologicalSequence("gene-1", "ATGGCC" * 100, "gene"),
    BiologicalSequence("protein-1", "MAMAPRTEINSTRING" * 10, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="central_dogma",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example result (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 19) SupKTax classification

- **Purpose:** Classify gene sequences using the official SupKTax checkpoint.
- **Prediction method:** Run the official LucaOneTasks `SupKTax` multiclass checkpoint directly, returning the `top_k` highest-probability labels from `label.txt`.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="supktax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. Interpret labels using the checkpoint’s `label.txt`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="supktax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **Example result (schematic structure):**

```text
(5,)
(<five highest-probability official labels>,)
```

## 20) GenusTax classification

- **Purpose:** Classify gene sequences using the official GenusTax checkpoint.
- **Prediction method:** Run the official LucaOneTasks `GenusTax` multiclass checkpoint directly, returning the Top-k labels by probability.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="genustax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. The LucaOneTasks README swaps the natural-language taxonomic levels for SupKTax/GenusTax, so use the actual `label.txt` to interpret results.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="genustax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["predicted_label"])
```

- **Example result (schematic structure):**

```text
(5,)
<highest-probability official label>
```

## 21) Species classification

- **Purpose:** Predict the species class of a gene sequence.
- **Prediction method:** Run the official LucaOneTasks `SpeciesTax` multiclass checkpoint directly, returning Top-k species labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="speciestax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`; use `top_k` to set the number of labels returned.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="speciestax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **Example result (schematic structure):**

```text
(5,)
(<five highest-probability species labels>,)
```

## 22) Prokaryotic protein subcellular localization

- **Purpose:** Predict subcellular localization classes for prokaryotic proteins.
- **Prediction method:** Run the official LucaOneTasks `ProtLoc` multiclass checkpoint directly, returning Top-k localization labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_location", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="protein"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [
        BiologicalSequence(
            "protein-1",
            "MAMAPRTEINSTRING" * 10,
            "protein",
        )
    ],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_location",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **Example result (schematic structure):**

```text
(5,)
(<five highest-probability subcellular localization labels>,)
```

## 23) Protein stability

- **Purpose:** Predict a stability regression value from a protein sequence under the official ProtStab task definition.
- **Prediction method:** Run the official LucaOneTasks `ProtStab` regression checkpoint through `predict_v1.py`. DNAKit does not fit a regressor.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_stability", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="protein"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [
        BiologicalSequence(
            "protein-1",
            "MAMAPRTEINSTRING" * 10,
            "protein",
        )
    ],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_stability",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example result (value depends on input and checkpoint):**

```text
('prediction',)
(1,)
```

## 24) ncRNA families

- **Purpose:** Predict the family of a noncoding RNA sequence.
- **Prediction method:** Run the official LucaOneTasks `ncRNAFam` multiclass checkpoint directly, returning Top-k family labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="ncrna_family", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. DNA or RNA alphabets are accepted; the official task pipeline handles `U`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_family",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **Example result (schematic structure):**

```text
(5,)
(<five highest-probability ncRNA family labels>,)
```

## 25) Influenza A antigenic relationships

- **Purpose:** Predict the antigenic relationship class between two influenza A gene sequences.
- **Prediction method:** Run the official LucaOneTasks `InfA` binary-classification checkpoint directly, returning the positive-class probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="influenza_antigenicity", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`; both sequences must have `sequence_type="gene"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "influenza-pair-1",
    BiologicalSequence("strain-a", "ATGGCC" * 100, "gene"),
    BiologicalSequence("strain-b", "ATGGCT" * 100, "gene"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="influenza_antigenicity",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example result (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 26) Protein–protein interactions

- **Purpose:** Predict whether two protein sequences have a protein–protein interaction.
- **Prediction method:** Run the official LucaOneTasks `PPI` binary-classification checkpoint directly, returning the interaction probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_interaction", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`; both sequences must have `sequence_type="protein"`.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "protein-pair-1",
    BiologicalSequence("protein-a", "MAMAPRTEINSTRING" * 10, "protein"),
    BiologicalSequence("protein-b", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example result (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 27) ncRNA–protein interactions

- **Purpose:** Predict whether an ncRNA and a protein sequence interact.
- **Prediction method:** Run the official LucaOneTasks `ncRPI` binary-classification checkpoint directly, returning the interaction probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="ncrna_protein_interaction", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, with `sequence_type="gene"` first and `sequence_type="protein"` second.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "ncrna-protein-1",
    BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene"),
    BiologicalSequence("protein-1", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example result (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 28) H2A.Z histone-variant-region classification

- **Purpose:** Classify whether a sequence corresponds to an H2A.Z (H2AFZ) ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h2afz"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H2AFZ.ckpt`; the default path is `./ckpt/enformer-benchmarks/H2AFZ.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 29) H3K27ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k27ac"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K27ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27ac.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27ac",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 30) H3K27me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k27me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K27me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27me3.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 31) H3K36me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K36me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k36me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K36me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K36me3.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k36me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 32) H3K4me1 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me1"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K4me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me1.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me1",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 33) H3K4me2 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me2 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me2"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K4me2.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me2.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me2",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 34) H3K4me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K4me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me3.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 35) H3K9ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k9ac"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K9ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9ac.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9ac",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 36) H3K9me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k9me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H3K9me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9me3.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **CUDA smoke result with the local sample checkpoint:**

```text
('0', '1')
[0.04569203 0.95430803]
1
```

## 37) H4K20me1 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H4K20me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h4k20me1"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `H4K20me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H4K20me1.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h4k20me1",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 38) Enhancer classification

- **Purpose:** Classify enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="enhancers"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `enhancers.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 39) Enhancer-type classification

- **Purpose:** Classify the enhancer category represented by the input sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="enhancers_types"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `enhancers_types.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers_types.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers_types",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 40) Promoter classification

- **Purpose:** Classify promoter versus non-promoter sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_all"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `promoter_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_all.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 41) Non-TATA promoter classification

- **Purpose:** Classify non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_no_tata"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `promoter_no_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_no_tata",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 42) TATA promoter classification

- **Purpose:** Classify TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_tata"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `promoter_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_tata.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_tata",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 43) Splice-acceptor-site classification

- **Purpose:** Classify splice acceptor versus non-acceptor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_acceptors"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `splice_sites_acceptors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_acceptors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 44) Splice-site-type classification

- **Purpose:** Classify no-splice, splice-acceptor, or splice-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_all"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `splice_sites_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 45) Splice-donor-site classification

- **Purpose:** Classify splice donor versus non-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_donors"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `splice_sites_donors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_donors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 46) Coding-region versus intergenic-region classification

- **Purpose:** Classify coding versus intergenic genomic sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="demo_coding_vs_intergenomic_seqs"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `demo_coding_vs_intergenomic_seqs.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_coding_vs_intergenomic_seqs",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 47) Human versus worm sequence classification

- **Purpose:** Classify a sequence as human or Caenorhabditis elegans.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="demo_human_or_worm"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `demo_human_or_worm.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_human_or_worm",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 48) Drosophila enhancer classification

- **Purpose:** Classify Drosophila enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="drosophila_enhancers_stark"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `drosophila_enhancers_stark.ckpt`; the default path is `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="drosophila_enhancers_stark",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 49) Mouse enhancer classification

- **Purpose:** Classify mouse Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="dummy_mouse_enhancers_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `dummy_mouse_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="dummy_mouse_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 50) Human Cohn enhancer classification

- **Purpose:** Classify human Cohn enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_enhancers_cohn"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `human_enhancers_cohn.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="human_enhancers_cohn",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 51) Human Ensembl enhancer classification

- **Purpose:** Classify human Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_enhancers_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `human_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="human_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 52) Human regulatory-element-type classification

- **Purpose:** Classify human enhancer, promoter, or open-chromatin-region sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_ensembl_regulatory"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `human_ensembl_regulatory.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="human_ensembl_regulatory",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 53) Human non-TATA promoter classification

- **Purpose:** Classify human non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_nontata_promoters"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `human_nontata_promoters.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", ("ACGT" * 63)[:251])],
    task="human_nontata_promoters",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **CUDA smoke result with the local sample checkpoint:**

```text
('negative', 'positive')
[3.8891116e-08 1.0000000e+00]
positive
```

## 54) Human open-chromatin-region classification

- **Purpose:** Classify human Ensembl open-chromatin region versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_ocr_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Automatically download (if absent; see [checkpoint details](#checkpoint-cache)) `human_ocr_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`; a custom directory or exact path can also be used.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 80)],
    task="human_ocr_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```
