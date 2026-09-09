# Property prediction: deep learning and machine learning

This page covers user-defined deep-learning training, supervised machine-learning training, and the 54 tasks for which model authors released a **trained prediction head or native pretrained output**. The existing [neural representations](08_fingerprints.md) API remains embedding-only.

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
prediction = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir,
    trust_model=True,  # only for artifacts from a trusted source
)
print(prediction.predictions)
```

For numeric input use `feature_source="numeric"` and, for example, `DataMapping(sequence_field=None, feature_fields=("gc", "length"), label_fields=("target",))`. For embeddings use `feature_source="embedding", embedding_model="dnabert2"` and the checkpoint's required loading options, including `allow_remote_code` where applicable.

### Training and persistence

Defaults remain seed 42 and an 80/10/10 split; a mapped split column takes precedence. Target semantics and metrics match the deep-learning path. Multilabel, multioutput classification and multi-target regression fit an independent estimator per target. Constant binary targets use constant classifiers. Class predictions use the highest probability; multilabel predictions use `threshold`.

`model_params` supplies serializable estimator constructor parameters, not arbitrary estimator objects. `workers` controls default `n_jobs` and numerical-library threads (0 means 1); CatBoost uses `thread_count`. Explicit estimator parameters can override their parallelism and seeds. `standardize_features=True` fits scaling only on training rows. Nonnegative-only `multinomial_nb` and `complement_nb` use MaxAbsScaler; negative inputs are rejected instead of silently shifted.

Classifiers without `predict_proba` use sigmoid calibration with up to three stratified folds inside the training split. Each class needs at least two training samples. Scaling is also fitted inside each calibration fold. Validation/test rows are excluded, and calibration is recorded in result warnings. See the [scikit-learn calibration API](https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html).

The ML path does not use neural `epochs`, AdamW settings or gradient accumulation. Set iterations/tree counts through `model_params`. Validation is evaluation-only; there is no automatic tuning or best-model selection.

Regression targets use training-set mean/std normalization, as in the existing API; reports and predictions return original units. Classification `class_weight` / `class_weights` dictionaries use zero-based encoded class indices, ordered as saved in `model.json` under `target_schema.class_names`. Options such as `"balanced"` remain available where the estimator supports them.

Artifacts are `model.json`, `estimator.joblib`, `metrics.json`, `data_inspection.json`, `splits.csv` and `test_predictions.jsonl`. Fitted transforms are stored with estimators; the manifest records feature order, target schema, settings and Python/library versions. The returned `history` is empty because there is no neural training history.

Loading joblib can execute Python code, so ML inference requires explicit `trust_model=True` for artifacts you trained or otherwise trust. Reloading checks the recorded Python and dependency versions and requires the original environment if they differ. See [scikit-learn persistence guidance](https://scikit-learn.org/stable/model_persistence.html).

Edit `others/examples/train_machine_learning_properties.py` to use your own data, or run synthetic training and reload verification:

```bash
python others/examples/train_machine_learning_properties.py --demo
```

## User-defined deep-learning classification and regression

Install the training and common tabular-format dependencies first: `python -m pip install "dnakit[training]"`.

`train_custom_property_model()` accepts one local data file and maps an ID, a DNA/RNA sequence, one or more targets, and an optional existing split. The default is a seeded 80/10/10 train/validation/test split with seed 42. The caller must explicitly select `task_type="classification"` or `task_type="regression"`; DNAKit does not infer the scientific task merely from how values look.

Built-in readers cover CSV, TSV, other single-character delimiters, JSON, JSONL, YAML, XML, Excel, Parquet, SQLite, FASTA, FASTQ, and GenBank. Dot paths address nested records, and `key=value` FASTA/FASTQ description tokens become `header.key`. A caller can implement the `DataReader` protocol for any other format. An LLM may propose only a validated declarative `DataMapping`; DNAKit never executes model-generated Python, SQL, regular expressions, or shell commands.

One categorical target becomes binary or multiclass classification. Multiple binary columns become one-hot multiclass or multilabel targets, while multiple general categorical columns become independent multioutput targets. One or more numeric columns become single- or multi-target regression. `classification_mode` can override automatic classification-shape detection.

Classification output includes accuracy, balanced accuracy, precision, recall, F1, Jaccard, MCC, Cohen's kappa, log loss, AUROC, average precision, confusion matrices, and per-class/per-label results. Multilabel output also includes subset accuracy, Hamming loss, LRAP, and coverage error. Regression output includes MAE, MSE, RMSE, R², explained variance, median absolute error, max error, bias, MAPE, SMAPE, Pearson, and Spearman metrics, both aggregate and per target. Undefined metrics are `null`; classification metrics such as F1 are not reported for regression.

The default is full fine-tuning. DNABERT-2, NT-v2, HyenaDNA, Caduceus, GROVER, LucaOne, and GENERator currently expose the portable Transformers full-training path. Every registered representation model supports a frozen `linear_probe`, including AlphaGenome, Enformer, Evo 2, and JanusDNA. Query `available_custom_training_models()` for the version-specific support matrix. Models that load checkpoint-provided code still require `allow_remote_code=True`.

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

Workers fetch records and assemble CPU batches for training, validation, and testing; encoder/GPU operations remain in the main process. This setting does not parallelize initial file parsing or LLM requests. Changing the worker count does not change the seeded split or epoch sample order. The trainer currently runs on one selected `device` (for example `"cuda:0"`), not multi-GPU DDP; available GPU count does not multiply the batch size. `gradient_accumulation_steps` defaults to `1` and is also configurable through Python.

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

Override `api_key_env` for your own environment-variable name or pass `api_key` in Python without committing it to source control. Explicit `base_url` overrides support regional or private gateways. HTTP redirects are refused so that keys and sampled profiles are not forwarded to another endpoint. Gemini uses its [official compatibility endpoint](https://ai.google.dev/gemini-api/docs/openai), not native `generateContent` or Vertex AI. Qwen defaults to the still-supported Beijing endpoint; workspace-specific endpoints and API keys must match the region, as described in [Alibaba Cloud's documentation](https://help.aliyun.com/zh/model-studio/compatibility-of-openai-with-dashscope). Kimi's international endpoint is `https://api.moonshot.ai/v1`, per its [quickstart](https://platform.kimi.ai/docs/overview). OpenRouter requires its complete model ID; check the chosen model's supported parameters and downstream provider/data policies in the [quickstart](https://openrouter.ai/docs/quickstart) and [structured-output documentation](https://openrouter.ai/docs/guides/features/structured-outputs).

Optional Python settings:

- `temperature=None` (default) omits temperature for models that reject custom sampling parameters.
- `token_limit_parameter="auto"` sends `max_completion_tokens` for OpenAI and `max_tokens` for other compatible services; either field can also be selected explicitly. This follows [OpenAI Docs](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create), which notes that the older field is incompatible with some reasoning models.
- `json_mode=True` requests JSON Object output. Set it to `False` for a model without this parameter; the prompt and mapping validation still require a JSON object.
- `extra_body` passes model-specific fields, such as `{"enable_thinking": False}` for supported Qwen models or `{"thinking": {"type": "disabled"}}` for supported Kimi models. Check the [Qwen](https://help.aliyun.com/zh/model-studio/qwen-structured-output) or [Kimi parameter reference](https://platform.kimi.ai/docs/api/models-overview) before applying these options; they are not supported by every model. Requests are non-streaming, so models requiring streamed thinking must disable thinking or be replaced with a compatible model.

Extra fields cannot replace authentication, messages, model selection, token limits, or other core request fields. Truncated responses fail explicitly: adjust `max_output_tokens` or thinking settings yourself; DNAKit never increases the budget or retries another provider automatically. Presets describe interface support, not universal model compatibility. Provider contracts are covered by mocked tests, not paid live tests against every service.

`LLMSchemaConfig` uses OpenAI-compatible `POST /chat/completions` with JSON Output, including DeepSeek and compatible providers. By default, it shares at most 8 sampled records, 12,000 characters, and a structural profile—not the whole source. Even a file with fewer than 8 records retains at least one complete unshared record; an unparsed text source shares at most half of its probe. Set `include_sample_values=False` to share only field/type structure. `DataInspection` records how many sample rows and characters were sent. API keys are excluded from results, exceptions, and manifests; the Python example reads the key from an environment variable.

When the fields are known, omit the LLM and use `DataMapping(sequence_field="sequence", label_fields=("label",), id_field="index")` through the Python API.

## Integrated direct tasks

| `model` | Direct `task` values | Output | Input |
| --- | --- | --- | --- |
| `alphagenome` | `atac`, `cage`, `dnase`, `rna_seq`, `chip_histone`, `chip_tf`, `splice_sites`, `splice_site_usage`, `splice_junctions`, `contact_maps`, `procap` | tissue/cell-type tracks; all except `splice_junctions` also support REF, ALT and ALT−REF arrays | DNA sequence or one-SNV context |
| `enformer` | `human_tracks`, `mouse_tracks`, plus 27 NT Revised/Genomic Benchmarks classifiers listed below | regulatory tracks or checkpoint-ordered class probabilities | DNA sequence or one-SNV context for tracks; DNA sequence for classifiers |
| `segmentnt` | `genomic_segmentation` | per-base probabilities for 14 genomic annotations | human DNA, up to 30 kb |
| `evo2` | `variant_effect`, `exon_probability` | zero-shot likelihood delta or trained exon probability | one-SNV context or forward/reverse contexts |
| `generator` | `variant_effect` | `p_ref`, `p_alt`, and `log(p_ref/(p_alt+1e-10))` | one-SNV context |
| `lucaone` | ten LucaOneTasks heads | classification probabilities/Top-k labels or regression | gene/protein sequences and pairs |

The papers, official repositories, and checkpoint sources used by this chapter are collected in [FAQ: references for deep-learning property prediction](../../faq.md#deep-learning-property-prediction-references).

AlphaGenome predicts RNA-seq, CAGE, PRO-cap, ATAC/DNase accessibility, histone and TF ChIP, splice sites/usage/junctions, and contact maps. `ontology_terms` limits tissue or cell-type tracks. DNAKit center-pads to an official supported length (16,384, 131,072, 524,288, or 1,048,576 bp) and records that padding.

SegmentNT returns these 14 features: `protein_coding_gene`, `lncRNA`, `exon`, `intron`, `splice_donor`, `splice_acceptor`, `5UTR`, `3UTR`, `CTCF-bound`, `polyA_signal`, tissue-specific/invariant enhancers, and tissue-specific/invariant promoters. This uses the trained `InstaDeepAI/segment_nt` checkpoint, not the NT-v2 base checkpoint.

LucaOneTasks exposes `central_dogma`, `supktax`, `genustax`, `speciestax`, `protein_location`, `protein_stability`, `ncrna_family`, `influenza_antigenicity`, `protein_interaction`, and `ncrna_protein_interaction`. DNAKit uses the exact timestamps and steps in the official `prediction.sh`. Because the upstream README swaps the natural-language descriptions of SupKTax and GenusTax, the checkpoint's actual `label.txt` remains authoritative.

## Unified API

Registry queries do not download weights:

```python
from dnakit.predictions import available_prediction_models, available_prediction_tasks

print(available_prediction_models())
print(available_prediction_tasks("evo2"))
```

Use `predict_properties()` or the typed convenience wrappers `predict_sequence_properties()`, `predict_pair_properties()`, and `predict_variant_effects()`. A result contains one read-only float32 array per input. `output_names` names the final axis, while metadata records axes, padding, resolution, thresholds, and score definitions. Use `to_dict(include_values=False)` to avoid expanding large tracks into JSON.

### Single-nucleotide genome segmentation

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
print(result.records[0].output.values.shape)  # (4000, 14)
```

### Variant scoring

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
    config=PropertyPredictionConfig(model="evo2", task="variant_effect"),
)
print(result.records[0].output.values)
```

Evo 2's likelihood delta and GENERator's log-ratio are continuous zero-shot ranking scores, not calibrated clinical pathogenicity probabilities.

### Released downstream task heads

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ACGT" * 100)],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="speciestax",
        model_source_path="/opt/LucaOneTasks",
        device="cuda:0",
    ),
)
print(result.records[0].output.metadata["predicted_label"])
```

This adapter executes only the official `src/predict_v1.py` inside the explicitly supplied checkout. Its first run may download LucaOne and task-head checkpoints.

## 28) H2A.Z histone-variant-region classification

- **Purpose:** Classify whether a sequence corresponds to an H2A.Z (H2AFZ) ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h2afz", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H2AFZ.ckpt`; the default path is `./ckpt/enformer-benchmarks/H2AFZ.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 29) H3K27ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k27ac", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K27ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27ac.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 30) H3K27me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k27me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K27me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27me3.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 31) H3K36me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K36me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k36me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K36me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K36me3.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 32) H3K4me1 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me1", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me1.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 33) H3K4me2 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me2 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me2", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me2.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me2.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 34) H3K4me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me3.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 35) H3K9ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k9ac", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K9ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9ac.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 36) H3K9me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k9me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K9me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9me3.ckpt`.
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
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h4k20me1", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H4K20me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H4K20me1.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 38) Enhancer classification

- **Purpose:** Classify enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="enhancers", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `enhancers.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="enhancers",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 39) Enhancer-type classification

- **Purpose:** Classify the enhancer category represented by the input sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="enhancers_types", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `enhancers_types.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers_types.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="enhancers_types",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 40) Promoter classification

- **Purpose:** Classify promoter versus non-promoter sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_all", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_all.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 41) Non-TATA promoter classification

- **Purpose:** Classify non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_no_tata", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_no_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 42) TATA promoter classification

- **Purpose:** Classify TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_tata", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_tata.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 43) Splice-acceptor-site classification

- **Purpose:** Classify splice acceptor versus non-acceptor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_acceptors", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_acceptors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_acceptors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 44) Splice-site-type classification

- **Purpose:** Classify no-splice, splice-acceptor, or splice-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_all", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 45) Splice-donor-site classification

- **Purpose:** Classify splice donor versus non-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_donors", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_donors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_donors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 46) Coding-region versus intergenic-region classification

- **Purpose:** Classify coding versus intergenic genomic sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="demo_coding_vs_intergenomic_seqs", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `demo_coding_vs_intergenomic_seqs.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 47) Human versus worm sequence classification

- **Purpose:** Classify a sequence as human or Caenorhabditis elegans.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="demo_human_or_worm", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `demo_human_or_worm.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 48) Drosophila enhancer classification

- **Purpose:** Classify Drosophila enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="drosophila_enhancers_stark", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `drosophila_enhancers_stark.ckpt`; the default path is `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="drosophila_enhancers_stark",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 49) Mouse enhancer classification

- **Purpose:** Classify mouse Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="dummy_mouse_enhancers_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `dummy_mouse_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 25)],
    task="dummy_mouse_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 50) Human Cohn enhancer classification

- **Purpose:** Classify human Cohn enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_enhancers_cohn", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_enhancers_cohn.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`.
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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 51) Human Ensembl enhancer classification

- **Purpose:** Classify human Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_enhancers_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="human_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 52) Human regulatory-element-type classification

- **Purpose:** Classify human enhancer, promoter, or open-chromatin-region sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_ensembl_regulatory", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_ensembl_regulatory.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="human_ensembl_regulatory",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 53) Human non-TATA promoter classification

- **Purpose:** Classify human non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_nontata_promoters", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_nontata_promoters.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`.
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
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_ocr_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_ocr_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="human_ocr_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```


## Dependencies and boundaries

| Model | Requirement |
| --- | --- |
| SegmentNT, GENERator | `dnakit[neural]` and explicit `allow_remote_code=True` after reviewing official code |
| Enformer | `dnakit[neural,neural-enformer]` |
| Evo 2 | Python 3.11–3.12, `dnakit[neural,neural-evo2]`, compatible GPU |
| AlphaGenome | Python ≥3.11, official research package/JAX, accepted model terms, generally H100-class GPU |
| LucaOneTasks | its separate official environment and source checkout |

`checkpoint_dir` changes the cache root; `checkpoint_path` selects prepared local weights. Prediction loops show progress unless `show_progress=False`.

The NT-v2 base model, DNABERT-2, HyenaDNA, Caduceus-Ph, JanusDNA, and CrossDNA still require a separately trained task head. This addition runs only the 27 supplied full-model task checkpoints; it does not turn those embedding models into zero-shot classifiers. ConvNova does not provide a reusable official pretrained task checkpoint. SegmentNT, the Evo 2 exon classifier, and LucaOneTasks are included because their trained downstream heads are separately released.
