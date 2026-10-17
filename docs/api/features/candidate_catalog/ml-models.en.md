# Machine learning models

This category contains **48** documented functions. 25 compute a result from ordinary inputs (`export_torchscript_model` additionally needs PyTorch). 5 compute only once something outside DNAKit is supplied and return `status: "not_implemented"` until then: `attribute`, `mutagenesis` and `explain` need a Python callable (a scoring function or a model) in `options`, `finetune` and `generate_sequence_embeddings` need local pretrained weights (`checkpoint_path`) or a representation backend. 18 always return `not_implemented`: the 12 architecture exports (`export_chromatin_model`, `export_cnn_model`, `export_expression_model`, `export_graph_model`, `export_long_context_model`, `export_pangenome_model`, `export_phenotype_model`, `export_protein_model`, `export_rnn_model`, `export_splice_model`, `export_transformer_model`, `export_variant_effect_model`), the two ONNX exports, `pretrain`, `quantize_model`, `run_cross_species_transfer` and `serve`. All functions are exported through the DNAKit namespace and the Agent/MCP registry, but `attribute`, `mutagenesis`, `explain` and `export_torchscript_model` need a Python callable or a torch module that a JSON call plan cannot carry, so they can only be used from Python, not through JSON/MCP. Examples below shorten long arrays and text with `...`.

## 1) attribute

- **Purpose:** Compute per-position in-silico mutagenesis attributions for one DNA sequence with a user-supplied scoring function.

- **Method:** For position i, `attribution[i] = score(sequence) - mean(score of the three single-base substitutions at i)`; the score is whatever the supplied callable returns. DNAKit ships no model: the scorer is your own Python function (or an object you wrap in one). Without a scorer the call returns `not_implemented`. **This function cannot be used through JSON/MCP**, because `options['scorer']` (or `options['model']`) must be a Python callable, which a JSON call plan cannot carry; call it from Python.

- **API:** `dnakit.ml_models.attribute(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one unambiguous DNA sequence (A/C/G/T; an ambiguous base raises `ML_AMBIGUOUS_SEQUENCE`). `options`: `scorer` or `model` (required, Python callable sequence -> number), `positions` (0-based positions to score; default all), `max_length` (2000).

- **Output:** The standard envelope. `result` has `sequence_length`, `method` (`in_silico_mutagenesis`), `definition`, `reference_score`, `positions` and `attributions` (one value per position). Without a scorer: `status: "not_implemented"` with `required_dependency` and `message`.

- **Example code:**

```python
import json
from dnakit.ml_models import attribute

def a_count(sequence):
    return float(sequence.count("A"))

result = attribute("ACG", options={"scorer": a_count})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "attribute",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 3,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "77ef74518ccd6437d9f93ad907bf0bd990f86c2611175f3472aab7852bb1d0e4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_length": 3,
    "method": "in_silico_mutagenesis",
    "definition": "attribution[i] = score(sequence) - mean(score(sequence with base i replaced by each other base))",
    "reference_score": 1.0,
    "positions": [0, 1, 2],
    "attributions": [1.0, -0.3333333333333333, -0.3333333333333333]
  }
}
```

## 2) build

- **Purpose:** Describe the untrained scikit-learn estimator that DNAKit would create for a model name and parameters.

- **Method:** Builds a `CustomPropertyTrainingConfig` from the options and instantiates the estimator through DNAKit's estimator factory (`dnakit.predictions`); nothing is fitted and no data is read. Requires scikit-learn (the estimator families `random_forest`, `logistic_regression`, ... are listed by `ml_models.search`). Unknown option names raise `UNKNOWN_ML_OPTION`.

- **API:** `dnakit.ml_models.build(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `None` or a mapping of the options below (merged under `options`). `options`: `task_type` (required, `classification` or `regression`), `model` (default `random_forest`), `model_params` (dict of estimator parameters), plus other `CustomPropertyTrainingConfig` fields.

- **Output:** The standard envelope. `result` has `model`, `task_type`, `estimator_class` (fully qualified class name), `params` (the complete parameter dictionary of the estimator) and `trained: false`.

- **Example code:**

```python
import json
from dnakit.ml_models import build

result = build(None, options={"task_type": "classification", "model": "random_forest",
                              "model_params": {"n_estimators": 10}})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "build",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "random_forest",
    "task_type": "classification",
    "estimator_class": "sklearn.ensemble._forest.RandomForestClassifier",
    "params": {
      "bootstrap": true,
      "ccp_alpha": 0.0,
      "class_weight": null,
      "criterion": "gini",
      "max_depth": null,
      "max_features": "sqrt",
      "max_leaf_nodes": null,
      "max_samples": null,
      "min_impurity_decrease": 0.0,
      "min_samples_leaf": 1,
      "min_samples_split": 2,
      "min_weight_fraction_leaf": 0.0,
      "monotonic_cst": null,
      "n_estimators": 10,
      "n_jobs": 1,
      "oob_score": false,
      "random_state": 42,
      "verbose": 0,
      "warm_start": false
    },
    "trained": false
  }
}
```

## 3) calculate_prediction_uncertainty

- **Purpose:** Compute predictive uncertainty: entropy and margin for class probabilities, or an entropy decomposition / variance for an ensemble.

- **Method:** Shannon entropy `H = -sum p log p` (natural log by default), least confidence `1 - max p` and margin `p(1st) - p(2nd)` for one probability vector per sample. For an ensemble of probability matrices: predictive entropy = entropy of the member-averaged probabilities, expected entropy = mean member entropy, mutual information = predictive - expected entropy (epistemic part). For an ensemble of regression predictions: mean, variance, standard deviation, min and max across members (population variance, or sample variance with `ddof=1`). No calibration is applied.

- **API:** `dnakit.ml_models.calculate_prediction_uncertainty(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping with one of `probabilities` (samples x classes, or members x samples x classes) or `predictions` (members x samples regression values). `options`: `log_base` (default e), `ddof` (default 0; use 1 for sample variance across members).

- **Output:** The standard envelope. `result` has `kind` (`classification`, `classification_ensemble` or `regression_ensemble`), `n` and `samples`. A classification sample has `predicted_class`, `confidence`, `entropy`, `least_confidence`, `margin`; an ensemble classification sample has `mean_probabilities`, `predicted_class`, `predictive_entropy`, `expected_entropy`, `mutual_information`; a regression ensemble sample has `mean`, `variance`, `std`, `min`, `max`.

- **Example code:**

```python
import json
from dnakit.ml_models import calculate_prediction_uncertainty

result = calculate_prediction_uncertainty({"probabilities": [[0.5, 0.5], [0.9, 0.1]]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "calculate_prediction_uncertainty",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "53d478214e363e76beca7009946ab2ce05820050a42638ac7076d894bf9b271b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "classification",
    "n": 2,
    "samples": [
      {
        "predicted_class": 0,
        "confidence": 0.5,
        "entropy": 0.6931471805599453,
        "least_confidence": 0.5,
        "margin": 0.0
      },
      {
        "predicted_class": 0,
        "confidence": 0.9,
        "entropy": 0.3250829733914482,
        "least_confidence": 0.09999999999999998,
        "margin": 0.8
      }
    ]
  }
}
```

## 4) calibrate

- **Purpose:** Fit a probability calibrator (Platt scaling or isotonic regression) on binary scores and return the calibrated probabilities and the expected calibration error.

- **Method:** Platt: `p = 1/(1+exp(A*s+B))` fitted by the Lin-Lin-Weng smoothed-target Newton method (agrees with scikit-learn's sigmoid calibration in the tests). Isotonic: pool-adjacent-violators on the scores sorted ascending, linear interpolation between blocks, clipped outside the fitted range. ECE is the bin-weighted mean of |mean label - mean score| over `n_bins` equal-width bins of [0, 1]; it is measured on the fitting data (resubstitution), so `ece_after_on_fit_data` is optimistic and `ece_before` is `null` when scores fall outside [0, 1].

- **API:** `dnakit.ml_models.calibrate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping of equally long columns `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`, a list of such records, or a CSV/TSV/JSON/JSONL path; the column names can be changed with `true_column`, `pred_column` and `score_column`; `y_score` is required. `options`: `method` (`platt` default, or `isotonic`), `positive_label` (default the larger label), `apply_to` (scores to transform; default the fitting scores), `n_bins` (10).

- **Output:** The standard envelope. `result` has `method`, `model` (the fitted calibrator: `type`, `method`, `classes`, `positive_label`, `parameters`), `n_fit`, `calibrated` (probabilities for `apply_to`), `ece_before` (ECE of the raw scores; `null` when any fitting score lies outside [0, 1], because ECE needs probabilities), `ece_after_on_fit_data` (always reported) and a `note`.

- **Example code:**

```python
import json
from dnakit.ml_models import calibrate

data = {"y_true": [0, 0, 1, 0, 1, 1], "y_score": [0.1, 0.2, 0.35, 0.4, 0.8, 0.9]}
result = calibrate(data, options={"method": "isotonic", "apply_to": [0.15, 0.6]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "calibrate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 74,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e550f9b569b2e26747538240c96c650fe28d96d5680df20f5c22dc7ab0e44a3d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "isotonic",
    "model": {
      "type": "dnakit.calibration",
      "version": 1,
      "method": "isotonic",
      "classes": ["0", "1"],
      "positive_label": "1",
      "parameters": {
        "blocks": [
          {"x_first": 0.1, "x_last": 0.2, "value": 0.0},
          {"x_first": 0.35, "x_last": 0.4, "value": 0.5},
          {"x_first": 0.8, "x_last": 0.9, "value": 1.0}
        ]
      }
    },
    "n_fit": 6,
    "calibrated": [0.0, 0.75],
    "ece_before": 0.27499999999999997,
    "ece_after_on_fit_data": 0.0,
    "note": "ECE is measured on the data used to fit the calibrator (resubstitution); use held-out data for an unbiased estimate."
  }
}
```

## 5) check

- **Purpose:** Check prediction columns for equal lengths, missing values, non-finite scores and label consistency before evaluation.

- **Method:** Reads the columns and reports problems instead of raising: length mismatches (`ML_LENGTH_MISMATCH`) and missing values (`ML_MISSING_VALUES`) make `valid` false and `n` null; otherwise non-finite scores (`NON_FINITE_SCORE`), non-numeric regression columns (`NON_NUMERIC`) and missing columns are errors, and scores outside [0, 1] or predicted labels never seen in `y_true` are warnings. The task is inferred (`auto`) as regression when `y_true` contains non-integer numbers.

- **API:** `dnakit.ml_models.check(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping of equally long columns `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`, a list of such records, or a CSV/TSV/JSON/JSONL path; the column names can be changed with `true_column`, `pred_column` and `score_column`. `options`: `task` (`auto` default, `classification` or `regression`).

- **Output:** The standard envelope. `result` has `valid` (no errors), `n`, `task`, `errors` (items with `code`, `message`) and `warnings` (list of messages).

- **Example code:**

```python
import json
from dnakit.ml_models import check

result = check({"y_true": [0, 1], "y_pred": [0, 2], "y_score": [0.2, 1.7]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "check",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "35d924e53360654aff2b31ce7dc6a74ac892654614a2275375cc9f572ad69679"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": true,
    "n": 2,
    "task": "classification",
    "errors": [],
    "warnings": ["y_score has values outside [0, 1]; they are treated as raw scores, not probabilities.", "y_pred contains labels never seen in y_true: ['2']."]
  }
}
```

## 6) encode

- **Purpose:** Encode DNA sequences as one-hot matrices, integer codes or k-mer vectors using the DNAKit fingerprint encoders.

- **Method:** `one_hot`: one row per base over the columns A, C, G, T. `integer`: A=0, C=1, G=2, T=3 (ambiguity handled by `ambiguity_policy`). `kmer`: counts, frequencies or presence of all 4^k k-mers in lexicographic A,C,G,T order, optionally canonical (k-mer merged with its reverse complement). Limits: `max_sequences` sequences per call and `max_output_length` bases per sequence.

- **API:** `dnakit.ml_models.encode(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a sequence string, FASTA text, a list of sequences, records, or a table/FASTA path. `options`: `method` (`one_hot` default, `integer`, `kmer`), `k` (3, kmer), `mode` (`count`, `frequency`, `binary`), `canonical` (`False`), `ambiguity_policy` (one_hot: `error` default, `zero`, `fractional`; integer: `iupac`, `sentinel`; kmer: `ignore`, `error`), `sequence_column` / `id_column` (`sequence` / `id`), `max_sequences` (1000), `max_output_length` (100000).

- **Output:** The standard envelope. `result` has `method`, `n_sequences`, `encodings` (per sequence: `id`, `length`, `values`) and, for one_hot and kmer, `feature_names`.

- **Example code:**

```python
import json
from dnakit.ml_models import encode

result = encode(["ACGT", "TTA"], options={"method": "integer"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "encode",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 15,
    "sequence_bases": 7,
    "lengths": [4, 3],
    "sha256": "c3ac31c902756f55c7f34b1b8c26595f71d5fe231a9f97f8bce75d5da9d16a84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "integer",
    "n_sequences": 2,
    "encodings": [
      {"id": "seq1", "length": 4, "values": [0, 1, 2, 3]},
      {"id": "seq2", "length": 3, "values": [3, 3, 0]}
    ]
  }
}
```

## 7) evaluate

- **Purpose:** Compute classification or regression metrics from predictions.

- **Method:** Classification: accuracy, balanced accuracy, macro/weighted precision, recall and F1, Matthews correlation (binary formula, Gorodkin's R_K for multiclass), Cohen's kappa, confusion matrix and per-class table; for binary tasks also precision, recall/sensitivity, specificity, F1, AUROC (rank statistic, average ranks for ties), AUPRC (average precision) and log loss when `y_score` is given. Regression: MAE, MSE, RMSE, R2, explained variance, Pearson and Spearman (undefined correlations are `null`). The task is `auto`: regression when `y_true` has non-integer numbers. The tests compare AUROC, average precision, MCC, F1, kappa, log loss, R2 and others against scikit-learn on fixed fixtures.

- **API:** `dnakit.ml_models.evaluate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping of equally long columns `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`, a list of such records, or a CSV/TSV/JSON/JSONL path; the column names can be changed with `true_column`, `pred_column` and `score_column`. `options`: `task` (`auto` default, `classification`, `regression`), `positive_label` (binary; default the larger label), `threshold` (0.5; used only when `y_pred` is absent and `y_score` is given), `classes` (label order), `metrics` (subset of metric names to return).

- **Output:** The standard envelope. `result` has `task`, `n` and `metrics`; classification adds `per_class` (precision, recall, f1, support) and `confusion_matrix` (`classes`, `matrix`).

- **Example code:**

```python
import json
from dnakit.ml_models import evaluate

data = {"y_true": [0, 0, 1, 1, 1, 0], "y_pred": [0, 1, 1, 1, 0, 0],
        "y_score": [0.1, 0.6, 0.8, 0.7, 0.4, 0.2]}
result = evaluate(data)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "evaluate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 103,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a0420a37d6f513ef723edff4da0739377894cd20eec327bfed038500efc105d7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "task": "classification",
    "n": 6,
    "metrics": {
      "n": 6,
      "accuracy": 0.6666666666666666,
      "balanced_accuracy": 0.6666666666666666,
      "precision_macro": 0.6666666666666666,
      "recall_macro": 0.6666666666666666,
      "f1_macro": 0.6666666666666666,
      "f1_weighted": 0.6666666666666666,
      "mcc": 0.3333333333333334,
      "cohen_kappa": 0.33333333333333326,
      "positive_label": "1",
      "precision": 0.6666666666666666,
      "recall": 0.6666666666666666,
      "sensitivity": 0.6666666666666666,
      "f1": 0.6666666666666666,
      "specificity": 0.6666666666666666,
      "tp": 2,
      "fp": 1,
      "tn": 2,
      "fn": 1,
      "auroc": 0.8888888888888888,
      "auprc": 0.9166666666666666,
      "log_loss": 0.45681733766221466
    },
    "per_class": {
      "0": {
        "precision": 0.6666666666666666,
        "recall": 0.6666666666666666,
        "f1": 0.6666666666666666,
        "support": 3
      },
      "1": {
        "precision": 0.6666666666666666,
        "recall": 0.6666666666666666,
        "f1": 0.6666666666666666,
        "support": 3
      }
    },
    "confusion_matrix": {
      "classes": ["0", "1"],
      "matrix": [
        [2, 1],
        [1, 2]
      ]
    }
  }
}
```

## 8) explain

- **Purpose:** Explain a model with permutation feature importance on tabular data.

- **Method:** Baseline score on the data (R2 for `regression`, accuracy for `classification`), then for each column the mean drop in that score over `n_repeats` random shuffles of the column (`importance_mean` = baseline - permuted score, with its standard deviation); the shuffling uses `random.Random(seed)`. The model is not supplied by DNAKit: `options['model']` must be a Python callable (or an object with `predict`) mapping a list of feature rows to predictions, otherwise the call returns `not_implemented`. **This function cannot be used through JSON/MCP**, because a callable cannot be carried in a JSON call plan; call it from Python.

- **API:** `dnakit.ml_models.explain(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'X': [[...]], 'y': [...], 'feature_names': [...]}` (`feature_names` optional; default `x0`, `x1`, ...). `options`: `model` (required, Python callable or object with `predict`; `options['scorer']` is accepted as an alias of `model`, and `model` wins when both are given), `task` (`regression` default, or `classification`), `n_repeats` (5), `seed` (42).

- **Output:** The standard envelope. `result` has `method` (`permutation_importance`), `metric` (`r2` or `accuracy`), `baseline_score`, `n_repeats` and `importances` (per feature: `feature`, `importance_mean`, `importance_std`). Without a model: `status: "not_implemented"` with `required_dependency` and `message`.

- **Example code:**

```python
import json
from dnakit.ml_models import explain

def model(rows):
    return [row[0] for row in rows]      # prediction depends only on feature "f1"

data = {"X": [[1, 5], [2, 3], [3, 8], [4, 1], [5, 9]], "y": [1, 2, 3, 4, 5],
        "feature_names": ["f1", "f2"]}
result = explain(data, options={"model": model, "task": "regression", "n_repeats": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "explain",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9f65a07f1c2350287e76f40e140f8bb4caf3f699469411bb5b635b93a9541e86"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "permutation_importance",
    "metric": "r2",
    "baseline_score": 1.0,
    "n_repeats": 3,
    "importances": [
      {"feature": "f1", "importance_mean": 2.2666666666666666, "importance_std": 0.33993463423951903},
      {"feature": "f2", "importance_mean": 0.0, "importance_std": 0.0}
    ]
  }
}
```

## 9) export_benchmark

- **Purpose:** Aggregate model-by-dataset benchmark results into per-dataset ranks and a model summary, and write them.

- **Method:** Within each dataset the models are ranked by value (standard competition ranking: ties share a rank, the next rank is skipped). The summary lists per model the number of datasets, mean value, mean rank and number of wins (rank 1), sorted by mean rank. The direction comes from `higher_is_better`, or is inferred from the metric name (lower is better for names containing loss, error, rmse, mae, mse, ece, nll, perplexity). A duplicated model/dataset pair is rejected with `ML_DUPLICATE_RESULT`.

- **API:** `dnakit.ml_models.export_benchmark(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: long-format records `{'model', 'dataset', 'value'}` or a table path. `options`: `model_column`, `dataset_column`, `value_column` (defaults `model`, `dataset`, `value`), `metric` (name used to decide the direction), `higher_is_better`, `format` (`json` default: per_dataset + summary; `tsv` / `csv` write the per-dataset table), `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `higher_is_better`, `format`, `per_dataset` (`dataset`, `model`, `value`, `rank`), `summary` (`model`, `n_datasets`, `mean_value`, `mean_rank`, `wins`) and the written-file fields `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_benchmark

rows = [
    {"model": "A", "dataset": "d1", "value": 0.9}, {"model": "B", "dataset": "d1", "value": 0.8},
    {"model": "A", "dataset": "d2", "value": 0.7}, {"model": "B", "dataset": "d2", "value": 0.8},
]
result = export_benchmark(rows)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_benchmark",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 188,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "57cac65e332b3518094499f796e4b10d0bdc704ae39dab85ee7ff68389cae3ac"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "higher_is_better": true,
    "format": "json",
    "per_dataset": [
      {"dataset": "d1", "model": "A", "value": 0.9, "rank": 1},
      {"dataset": "d1", "model": "B", "value": 0.8, "rank": 2},
      {"dataset": "d2", "model": "A", "value": 0.7, "rank": 2},
      {"dataset": "d2", "model": "B", "value": 0.8, "rank": 1}
    ],
    "summary": [
      {"model": "A", "n_datasets": 2, "mean_value": 0.8, "mean_rank": 1.5, "wins": 1},
      {"model": "B", "n_datasets": 2, "mean_value": 0.8, "mean_rank": 1.5, "wins": 1}
    ],
    "bytes": 655,
    "sha256": "065c29228ad19e6d40ec03896140ed39c702a8c6b6b97418fb8365bec4592984",
    "content": "{\n  \"per_dataset\": [\n    {\n      \"dataset\": \"d1\",\n      \"model\": \"A\",\n      \"value\": 0.9,\n      \"rank\": 1\n    },\n    {\n      \"dataset\": \"d1\",\n      \"model\": \"B\"...",
    "output_path": null
  }
}
```

## 10) export_calibration_model

- **Purpose:** Fit a Platt or isotonic calibrator on binary scores and write it as a JSON model file.

- **Method:** The same fit as `calibrate` (Platt: Lin-Lin-Weng smoothed-target Newton fit of `p = 1/(1+exp(A*s+B))`; isotonic: pool-adjacent-violators blocks with linear interpolation). The JSON file has `type: "dnakit.calibration"`, `version`, `method`, `classes`, `positive_label` and `parameters`.

- **API:** `dnakit.ml_models.export_calibration_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping of equally long columns `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`, a list of such records, or a CSV/TSV/JSON/JSONL path; the column names can be changed with `true_column`, `pred_column` and `score_column`; `y_score` is required. `options`: `method` (`platt` or `isotonic`), `positive_label`, `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `model` (the calibrator) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_calibration_model

data = {"y_true": [0, 0, 1, 1], "y_score": [0.1, 0.4, 0.6, 0.9]}
result = export_calibration_model(data, options={"method": "isotonic", "output": "calibration.json"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_calibration_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 57,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "93a0e7086bac066520556985b76dd7bc916a6e8aba4d95e30579692f1cbcf76b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "type": "dnakit.calibration",
      "version": 1,
      "method": "isotonic",
      "classes": ["0", "1"],
      "positive_label": "1",
      "parameters": {
        "blocks": [
          {"x_first": 0.1, "x_last": 0.4, "value": 0.0},
          {"x_first": 0.6, "x_last": 0.9, "value": 1.0}
        ]
      }
    },
    "bytes": 351,
    "sha256": "3379ac095b66564832d46e6927ce2c2bd4debb15a3b5cf058557aa74dc2e95bc",
    "output_path": "calibration.json"
  }
}
```

## 11) export_chromatin_model

- **Purpose:** Export a trained chromatin model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no chromatin model class to export and nothing is trained or serialised; a chromatin model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_chromatin_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_chromatin_model

result = export_chromatin_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_chromatin_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained chromatin model object and its architecture-specific serializer",
  "message": "DNAKit has no chromatin model class to export; nothing was trained or serialised."
}
```

## 12) export_cnn_model

- **Purpose:** Export a trained cnn model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no cnn model class to export and nothing is trained or serialised; a cnn model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_cnn_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_cnn_model

result = export_cnn_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_cnn_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained cnn model object and its architecture-specific serializer",
  "message": "DNAKit has no cnn model class to export; nothing was trained or serialised."
}
```

## 13) export_dataset

- **Purpose:** Write a dataset of records to CSV, TSV, JSON, JSONL or FASTA.

- **Method:** Reads the records and serialises them in the chosen format; FASTA uses the `id` column as header and the `sequence` column as sequence. The file is written atomically (temporary file then rename) and an existing file is never replaced unless `overwrite` is true.

- **API:** `dnakit.ml_models.export_dataset(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records, a column mapping, or a CSV/TSV/JSON/JSONL/FASTA path. `options`: `format` (from the `output` suffix, default `csv`), `columns` (subset and order), `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true). FASTA needs a `sequence` column.

- **Output:** The standard envelope. `result` has `format`, `n_records` and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_dataset

rows = [{"id": "s1", "sequence": "ACGT", "label": 1}, {"id": "s2", "sequence": "GG", "label": 0}]
result = export_dataset(rows, options={"format": "fasta"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_dataset",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 90,
    "sequence_bases": 6,
    "lengths": [4, 2],
    "sha256": "ed916cf32c3cd499c0f54338d0812cffb240170647ba4b99a3b73a7f49bb4bcb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "n_records": 2,
    "bytes": 16,
    "sha256": "d3dcb573654302517a57f0d7c34c2e0b1706be6d333f0001d397795dd305c803",
    "content": ">s1\nACGT\n>s2\nGG\n",
    "output_path": null
  }
}
```

## 14) export_expression_model

- **Purpose:** Export a trained expression model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no expression model class to export and nothing is trained or serialised; a expression model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_expression_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_expression_model

result = export_expression_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_expression_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained expression model object and its architecture-specific serializer",
  "message": "DNAKit has no expression model class to export; nothing was trained or serialised."
}
```

## 15) export_graph_model

- **Purpose:** Export a trained graph model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no graph model class to export and nothing is trained or serialised; a graph model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_graph_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_graph_model

result = export_graph_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_graph_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained graph model object and its architecture-specific serializer",
  "message": "DNAKit has no graph model class to export; nothing was trained or serialised."
}
```

## 16) export_interpretability

- **Purpose:** Export per-position attributions (supplied, or computed by in-silico mutagenesis with a scorer) as TSV or JSON.

- **Method:** With supplied `attributions` the values are written as they are (`method: "user_supplied"`). With one sequence and `options['scorer']` the attributions are computed as in `attribute` (score minus the mean score of the three substitutions, `method: "in_silico_mutagenesis"`). Only the supplied-attribution form works through JSON/MCP; the `scorer` form needs a Python callable and so only works from Python.

- **API:** `dnakit.ml_models.export_interpretability(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'sequence': 'ACGT', 'attributions': [...]}`, or one DNA sequence plus `options['scorer']` (callable sequence -> number). `options`: `format` (`tsv` or `json`), `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `method`, `format`, `n_positions` and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`. The TSV columns are `position`, `base`, `attribution`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_interpretability

result = export_interpretability({"sequence": "ACG", "attributions": [1.0, -0.25, 0.5]},
                                 options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_interpretability",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 54,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "86fc2c5469e64995c2214e7ff45d8407a644e95e93b328c62fdb9630ec18396c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "user_supplied",
    "format": "tsv",
    "n_positions": 3,
    "bytes": 50,
    "sha256": "a76e9bba121ed9b16af1330a282ef9e61756a14d8de97660a794b24e19e03b19",
    "content": "position\tbase\tattribution\n0\tA\t1\n1\tC\t-0.25\n2\tG\t0.5\n",
    "output_path": null
  }
}
```

## 17) export_kmer_embedding

- **Purpose:** Compute k-mer frequency (or count) vectors for sequences and write them as a TSV matrix or JSON.

- **Method:** Each sequence becomes a vector over all 4^k k-mers in A,C,G,T order (frequency = count / number of k-mers, or counts, or presence/absence), optionally canonical. This is the exact k-mer profile, not a learned embedding.

- **API:** `dnakit.ml_models.export_kmer_embedding(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: as in `encode` (sequence string, FASTA text, list, records, table/FASTA path). `options`: `k` (3), `mode` (`frequency` default, `count`, `binary`), `canonical` (`False`), `ambiguity_policy` (`ignore`), `format` (`tsv` or `json`), `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `k`, `mode`, `format`, `dimension`, `n_sequences` and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`. The TSV has a header `id` plus the k-mers, one row per sequence.

- **Example code:**

```python
import json
from dnakit.ml_models import export_kmer_embedding

result = export_kmer_embedding(["AAAA", "ACGT"], options={"k": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_kmer_embedding",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 16,
    "sequence_bases": 8,
    "lengths": [4, 4],
    "sha256": "debb33074aea8f8fcaa086cd144fb49d7cec8ac50769824ba1b3f6924f9c85a6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 1,
    "mode": "frequency",
    "format": "tsv",
    "dimension": 4,
    "n_sequences": 2,
    "bytes": 49,
    "sha256": "1fa9ecba797eaa9bbc97fa4cbb11d67aa7388d29655635d44848b95a5bb31151",
    "content": "id\tA\tC\tG\tT\nseq1\t1\t0\t0\t0\nseq2\t0.25\t0.25\t0.25\t0.25\n",
    "output_path": null
  }
}
```

## 18) export_leaderboard

- **Purpose:** Rank models by one metric and write the leaderboard as TSV, CSV, JSON or JSONL.

- **Method:** Models are sorted by the metric with standard competition ranking (equal values share a rank, the next rank is skipped). Direction: `higher_is_better`, or automatic (lower is better when the metric name contains loss, error, rmse, mae, mse, ece, nll or perplexity).

- **API:** `dnakit.ml_models.export_leaderboard(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records such as `[{'model': 'A', 'score': 0.9}, ...]` or a table path. `options`: `metric` (`score`), `model_column` (`model`), `higher_is_better` (auto), `format` (`tsv` default, `csv`, `json`, `jsonl`, or from the `output` suffix), `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `metric`, `higher_is_better`, `format`, `n`, `ranking` (`rank`, `model`, metric value) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_leaderboard

rows = [{"model": "a", "score": 0.9}, {"model": "b", "score": 0.8},
        {"model": "c", "score": 0.8}, {"model": "d", "score": 0.5}]
result = export_leaderboard(rows, options={"format": "json"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_leaderboard",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 120,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "47f9a912020c2e579d40c9984cb88103254860549918754debc660556057f9db"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "metric": "score",
    "higher_is_better": true,
    "format": "json",
    "n": 4,
    "ranking": [
      {"rank": 1, "model": "a", "score": 0.9},
      {"rank": 2, "model": "b", "score": 0.8},
      {"rank": 2, "model": "c", "score": 0.8},
      {"rank": 4, "model": "d", "score": 0.5}
    ],
    "bytes": 239,
    "sha256": "406c950e5890c02aedb274bc99661459a36734ea4d8a25a6995a96332151e510",
    "content": "[\n  {\n    \"rank\": 1,\n    \"model\": \"a\",\n    \"score\": 0.9\n  },\n  {\n    \"rank\": 2,\n    \"model\": \"b\",\n    \"score\": 0.8\n  },\n  {\n    \"rank\": 2,\n    \"model\": \"c\",\n   ...",
    "output_path": null
  }
}
```

## 19) export_long_context_model

- **Purpose:** Export a trained long context model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no long context model class to export and nothing is trained or serialised; a long context model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_long_context_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_long_context_model

result = export_long_context_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_long_context_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained long context model object and its architecture-specific serializer",
  "message": "DNAKit has no long context model class to export; nothing was trained or serialised."
}
```

## 20) export_model_onnx

- **Purpose:** Export a trained model to the ONNX format; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. ONNX export would need the `onnx` and `onnxruntime` packages (not installed here) to write and validate the graph, and a trained model object. For PyTorch modules use `export_torchscript_model`. The function is a duplicate entry point of `export_onnx_model`.

- **API:** `dnakit.ml_models.export_model_onnx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_model_onnx

result = export_model_onnx(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_model_onnx",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install onnx onnxruntime and a trained model object",
  "message": "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules."
}
```

## 21) export_onehot_embedding

- **Purpose:** One-hot encode sequences and write them as JSON or a long-format TSV.

- **Method:** Each base becomes a row over the columns A, C, G, T (a non-ACGT base follows `ambiguity_policy`: `error`, an all-zero row for `zero`, or for `fractional` the weight split equally over the bases of the IUPAC code, e.g. N = 0.25 each). The TSV has one row per position (`id`, `position`, `A`, `C`, `G`, `T`); the JSON maps each sequence id to its matrix. This is a fixed encoding, not a learned embedding.

- **API:** `dnakit.ml_models.export_onehot_embedding(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: as in `encode`. `options`: `format` (`json` or `tsv`; default from the `output` suffix), `ambiguity_policy` (`error` default, `zero`, `fractional`), `max_sequences`, `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `format`, `n_sequences`, `total_positions` and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_onehot_embedding

result = export_onehot_embedding([{"id": "s1", "sequence": "AC"}], options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_onehot_embedding",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 32,
    "sequence_bases": 2,
    "lengths": [2],
    "sha256": "6fa4b6fb7d80cb98f276d1517d2e3a9501748dd9111ad544e8a71cc26c91c3f4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "n_sequences": 1,
    "total_positions": 2,
    "bytes": 46,
    "sha256": "5aa43c64f3d9770a1af9bf6a57191c6ad88809a3b0de37e262ebb61c31df0919",
    "content": "id\tposition\tA\tC\tG\tT\ns1\t0\t1\t0\t0\t0\ns1\t1\t0\t1\t0\t0\n",
    "output_path": null
  }
}
```

## 22) export_onnx_model

- **Purpose:** Export a trained model to the ONNX format; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. ONNX export would need the `onnx` and `onnxruntime` packages (not installed here) to write and validate the graph, and a trained model object. For PyTorch modules use `export_torchscript_model`. The function is a duplicate entry point of `export_model_onnx`.

- **API:** `dnakit.ml_models.export_onnx_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_onnx_model

result = export_onnx_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_onnx_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install onnx onnxruntime and a trained model object",
  "message": "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules."
}
```

## 23) export_pangenome_model

- **Purpose:** Export a trained pangenome model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no pangenome model class to export and nothing is trained or serialised; a pangenome model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_pangenome_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_pangenome_model

result = export_pangenome_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_pangenome_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained pangenome model object and its architecture-specific serializer",
  "message": "DNAKit has no pangenome model class to export; nothing was trained or serialised."
}
```

## 24) export_phenotype_model

- **Purpose:** Export a trained phenotype model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no phenotype model class to export and nothing is trained or serialised; a phenotype model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_phenotype_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_phenotype_model

result = export_phenotype_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_phenotype_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained phenotype model object and its architecture-specific serializer",
  "message": "DNAKit has no phenotype model class to export; nothing was trained or serialised."
}
```

## 25) export_protein_model

- **Purpose:** Export a trained protein model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no protein model class to export and nothing is trained or serialised; a protein model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_protein_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_protein_model

result = export_protein_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_protein_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained protein model object and its architecture-specific serializer",
  "message": "DNAKit has no protein model class to export; nothing was trained or serialised."
}
```

## 26) export_rnn_model

- **Purpose:** Export a trained rnn model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no rnn model class to export and nothing is trained or serialised; a rnn model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_rnn_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_rnn_model

result = export_rnn_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_rnn_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained rnn model object and its architecture-specific serializer",
  "message": "DNAKit has no rnn model class to export; nothing was trained or serialised."
}
```

## 27) export_splice_model

- **Purpose:** Export a trained splice model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no splice model class to export and nothing is trained or serialised; a splice model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_splice_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_splice_model

result = export_splice_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_splice_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained splice model object and its architecture-specific serializer",
  "message": "DNAKit has no splice model class to export; nothing was trained or serialised."
}
```

## 28) export_test_split

- **Purpose:** Write the test part of a deterministic train/validation/test split of a dataset.

- **Method:** Reads the records and returns one part of a train/validation/test split. If every record has a split label (column `split_column`, default `split`; values train, validation|valid|val|dev, test; case-insensitive) that label is used (an unknown label raises `ML_UNKNOWN_SPLIT`). Otherwise a seeded random split is made with `random.Random(seed)`: the part sizes use largest-remainder allocation of `train_ratio`/`validation_ratio`/`test_ratio` (they must sum to 1, else `INVALID_SPLIT_RATIOS`); with `group_column` all records of a group stay in the same part. The same inputs and seed always give the same split, and the three exports together cover each record exactly once.

- **API:** `dnakit.ml_models.export_test_split(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records or a table path. `options`: `split_column` (`split`), `train_ratio` (0.8), `validation_ratio` (0.1), `test_ratio` (0.1), `seed` (42), `group_column` (none), `format` (`csv` default), `columns`, `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `split`, `format`, `n_records`, `counts` (records in each part), `row_indices` (0-based rows of the input in this part) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_test_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_test_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_test_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "test",
    "format": "csv",
    "n_records": 1,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 20,
    "sha256": "88bc5ac8ecfcc2658ce1f8c7c983864052ad5cfb27c4ecafbac3d77a5ae4f3cc",
    "content": "id,sequence\nr3,ACGT\n",
    "output_path": null,
    "row_indices": [3]
  }
}
```

## 29) export_tokenizer

- **Purpose:** Build a k-mer DNA tokenizer vocabulary and write it as JSON.

- **Method:** The vocabulary lists the special tokens first (ids 0, 1, ...) and then all 4^k k-mers in A,C,G,T order; `k` up to 8 is accepted (larger values raise a configuration error). The JSON has `type: "dnakit.kmer_tokenizer"`, the model (`k`, `stride`, `alphabet`, `vocab`), `special_tokens` and `unk_token`. No tokenizer is trained and the input is not used.

- **API:** `dnakit.ml_models.export_tokenizer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: ignored (may be `None`). `options`: `k` (3, at most 8), `stride` (1), `special_tokens` (`['[PAD]', '[UNK]', '[CLS]', '[SEP]', '[MASK]']`), `output` (JSON path; otherwise returned inline), `overwrite` (`False`).

- **Output:** The standard envelope. `result` has `k`, `stride`, `vocab_size`, `tokenizer` (the vocabulary object) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_tokenizer

result = export_tokenizer(None, options={"k": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_tokenizer",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 1,
    "stride": 1,
    "vocab_size": 9,
    "tokenizer": {
      "type": "dnakit.kmer_tokenizer",
      "version": 1,
      "model": {
        "type": "kmer",
        "k": 1,
        "stride": 1,
        "alphabet": "ACGT",
        "vocab": {"[PAD]": 0, "[UNK]": 1, "[CLS]": 2, "[SEP]": 3, "[MASK]": 4, "A": 5, "C": 6, "G": 7, "T": 8}
      },
      "special_tokens": ["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"],
      "unk_token": "[UNK]"
    },
    "bytes": 428,
    "sha256": "bbeb600527ff07e8930ef642b5d334cc49fa2ac734568f0c4705e81309079e97",
    "content": "{\n  \"type\": \"dnakit.kmer_tokenizer\",\n  \"version\": 1,\n  \"model\": {\n    \"type\": \"kmer\",\n    \"k\": 1,\n    \"stride\": 1,\n    \"alphabet\": \"ACGT\",\n    \"vocab\": {\n      ...",
    "output_path": null
  }
}
```

## 30) export_torchscript_model

- **Purpose:** Trace a PyTorch module with an example input and save it as a TorchScript archive, verifying the reloaded output.

- **Method:** `torch.jit.trace` in eval mode under `torch.no_grad()`, then `save`, `torch.jit.load` and a comparison of the original and reloaded outputs on the example input (the maximum absolute difference is reported). The module's previous training/eval mode is restored. Needs PyTorch (otherwise `not_implemented` with `required_dependency: "pip install torch"`). The input must be a `torch.nn.Module` object, so this function cannot be used through JSON/MCP; call it from Python. ONNX export is a separate, unimplemented function (`export_onnx_model`).

- **API:** `dnakit.ml_models.export_torchscript_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a `torch.nn.Module` (or `options['model']`). `options`: `example_input` (required, rectangular nested list of numbers, converted to a float32 tensor), `output` (required, `.pt` path), `overwrite` (`False`).

- **Output:** The standard envelope. `result` has `output_path`, `bytes`, `input_shape`, `roundtrip_max_abs_difference` and `torch_version`.

- **Example code:**

```python
import json
import torch
from dnakit.ml_models import export_torchscript_model

model = torch.nn.Linear(4, 2)
result = export_torchscript_model(model, options={"example_input": [[0.0, 1.0, 0.0, 0.0]],
                                                  "output": "linear.pt"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_torchscript_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 35,
    "lengths": [35],
    "sha256": "6eb7634af5ec8a7ccccd464481ae0dee1f0f1e9237132e1977d422a8a6cc68b9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "output_path": "linear.pt",
    "bytes": 3620,
    "input_shape": [1, 4],
    "roundtrip_max_abs_difference": 0.0,
    "torch_version": "2.6.0+cu118"
  }
}
```

## 31) export_train_split

- **Purpose:** Write the train part of a deterministic train/validation/test split of a dataset.

- **Method:** Reads the records and returns one part of a train/validation/test split. If every record has a split label (column `split_column`, default `split`; values train, validation|valid|val|dev, test; case-insensitive) that label is used (an unknown label raises `ML_UNKNOWN_SPLIT`). Otherwise a seeded random split is made with `random.Random(seed)`: the part sizes use largest-remainder allocation of `train_ratio`/`validation_ratio`/`test_ratio` (they must sum to 1, else `INVALID_SPLIT_RATIOS`); with `group_column` all records of a group stay in the same part. The same inputs and seed always give the same split, and the three exports together cover each record exactly once.

- **API:** `dnakit.ml_models.export_train_split(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records or a table path. `options`: `split_column` (`split`), `train_ratio` (0.8), `validation_ratio` (0.1), `test_ratio` (0.1), `seed` (42), `group_column` (none), `format` (`csv` default), `columns`, `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `split`, `format`, `n_records`, `counts` (records in each part), `row_indices` (0-based rows of the input in this part) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_train_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_train_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_train_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "train",
    "format": "csv",
    "n_records": 8,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 76,
    "sha256": "59391b1a49e78777e0b9a06b303fa4654b760ddbeab878620fd61bbfc93b6766",
    "content": "id,sequence\nr0,ACGT\nr1,ACGT\nr2,ACGT\nr4,ACGT\nr5,ACGT\nr6,ACGT\nr7,ACGT\nr9,ACGT\n",
    "output_path": null,
    "row_indices": [0, 1, 2, 4, 5, 6, "..."]
  }
}
```

## 32) export_transformer_model

- **Purpose:** Export a trained transformer model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no transformer model class to export and nothing is trained or serialised; a transformer model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_transformer_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_transformer_model

result = export_transformer_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_transformer_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained transformer model object and its architecture-specific serializer",
  "message": "DNAKit has no transformer model class to export; nothing was trained or serialised."
}
```

## 33) export_conformal_uncertainty_model

- **Purpose:** Fit a split-conformal uncertainty model on a held-out calibration set and write it as JSON.

- **Method:** Split conformal prediction: nonconformity scores are the absolute residual |y - y_hat| (regression) or 1 - p(true class) (classification). The stored quantile is the k-th smallest score with k = ceil((n+1)(1-alpha)); if k > n (calibration set too small) `rank` and `quantile` are `null`. The coverage statement stored in the file is the marginal guarantee 1 - alpha under exchangeability; DNAKit does not test exchangeability.

- **API:** `dnakit.ml_models.export_conformal_uncertainty_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'y_true': [...], 'y_pred': [...]}` (regression) or `{'y_true': [...], 'probabilities': [[...]], 'classes': [...]}` (classification). `options`: `alpha` (0.1), `output` (JSON path; otherwise returned inline), `overwrite` (`False`).

- **Output:** The standard envelope. `result` has `model` (`type: "dnakit.split_conformal"`, `task`, `alpha`, `n_calibration`, `rank`, `quantile`, `nonconformity`, `coverage_guarantee`) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_conformal_uncertainty_model

data = {"y_true": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], "y_pred": [0.0] * 9}
result = export_conformal_uncertainty_model(data, options={"alpha": 0.2})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_conformal_uncertainty_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 114,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "160c66dde6c7b5f182bfb495e54b0da7b0b0b19f194c5db8f96e4b6428bdb83e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "type": "dnakit.split_conformal",
      "version": 1,
      "task": "regression",
      "alpha": 0.2,
      "n_calibration": 9,
      "rank": 8,
      "quantile": 8.0,
      "nonconformity": "absolute residual |y - y_hat|",
      "coverage_guarantee": ">= 0.8 marginal coverage under exchangeability"
    },
    "bytes": 275,
    "sha256": "68e78ed9c4245f54b1334992778c392f8bb8643183305896ac3028e9e66eba1b",
    "content": "{\n  \"type\": \"dnakit.split_conformal\",\n  \"version\": 1,\n  \"task\": \"regression\",\n  \"alpha\": 0.2,\n  \"n_calibration\": 9,\n  \"rank\": 8,\n  \"quantile\": 8.0,\n  \"nonconfor...",
    "output_path": null
  }
}
```

## 34) export_validation_split

- **Purpose:** Write the validation part of a deterministic train/validation/test split of a dataset.

- **Method:** Reads the records and returns one part of a train/validation/test split. If every record has a split label (column `split_column`, default `split`; values train, validation|valid|val|dev, test; case-insensitive) that label is used (an unknown label raises `ML_UNKNOWN_SPLIT`). Otherwise a seeded random split is made with `random.Random(seed)`: the part sizes use largest-remainder allocation of `train_ratio`/`validation_ratio`/`test_ratio` (they must sum to 1, else `INVALID_SPLIT_RATIOS`); with `group_column` all records of a group stay in the same part. The same inputs and seed always give the same split, and the three exports together cover each record exactly once.

- **API:** `dnakit.ml_models.export_validation_split(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records or a table path. `options`: `split_column` (`split`), `train_ratio` (0.8), `validation_ratio` (0.1), `test_ratio` (0.1), `seed` (42), `group_column` (none), `format` (`csv` default), `columns`, `output` (file path; without it nothing is written and the text is returned inline), `overwrite` (`False`; an existing file is not replaced unless true).

- **Output:** The standard envelope. `result` has `split`, `format`, `n_records`, `counts` (records in each part), `row_indices` (0-based rows of the input in this part) and `bytes`, `sha256`, `output_path` and (when no `output` path is given) the text itself under `content`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_validation_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_validation_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_validation_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "validation",
    "format": "csv",
    "n_records": 1,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 20,
    "sha256": "f59306c2734ef4b9cdb82d3f09d0205dad106eae24bd55b63f9d94c0d6815b10",
    "content": "id,sequence\nr8,ACGT\n",
    "output_path": null,
    "row_indices": [8]
  }
}
```

## 35) export_variant_effect_model

- **Purpose:** Export a trained variant effect model to a portable file; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no variant effect model class to export and nothing is trained or serialised; a variant effect model would need its architecture-specific serializer.

- **API:** `dnakit.ml_models.export_variant_effect_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import export_variant_effect_model

result = export_variant_effect_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "export_variant_effect_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained variant effect model object and its architecture-specific serializer",
  "message": "DNAKit has no variant effect model class to export; nothing was trained or serialised."
}
```

## 36) finetune

- **Purpose:** Fine-tune a DNA model (linear probe or full) on a labelled file via `dnakit.predictions.train_custom_property_model`, given local pretrained weights.

- **Method:** A thin, checked delegation to the custom-property trainer. DNAKit does not download weights during this call: without `checkpoint_path` (local weights), `representation_backend` or `training_encoder` the call returns `not_implemented`. With them it trains a linear probe on the frozen representation (`finetune_method: "linear_probe"`) or fine-tunes the whole model (`"full"`), and writes the model directory. Other option names must be `CustomPropertyTrainingConfig` fields (otherwise `UNKNOWN_ML_OPTION`). Needs the torch-based neural runtime and the chosen model's own dependencies. Not run here with real weights.

- **API:** `dnakit.ml_models.finetune(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of a data file (CSV/TSV/JSON/FASTA). `options`: `checkpoint_path`, `representation_backend` or `training_encoder` (one is required), `task_type` (required), `output_dir` (required), `finetune_method` (`linear_probe` default, or `full`), `model`, and other `CustomPropertyTrainingConfig` fields.

- **Output:** Without weights, checkpoint or encoder: `status: "not_implemented"` with `required_dependency` (pretrained DNA foundation-model weights or a representation backend / training encoder) and `message`. Otherwise the standard envelope with `result.delegate` and `result.training` (model name, target schema, split counts, metrics, written artifacts).

- **Example code:**

```python
import json
from pathlib import Path
from dnakit.ml_models import finetune

Path("data.csv").write_text("sequence,label\nAAAACC,0\nGGGTTT,1\n")
# no checkpoint_path, representation_backend or training_encoder is given
result = finetune("data.csv", options={"task_type": "classification", "output_dir": "model_out"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "finetune",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 8,
    "sequence_bases": 8,
    "lengths": [8],
    "sha256": "803d8eb0fe6902e3880bee63c537f65f3917f39e743a36703b093591e7054c79"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pretrained DNA foundation-model weights (options['checkpoint_path']) or a representation_backend/training_encoder",
  "message": "Fine-tuning needs a pretrained model; no local checkpoint or encoder was supplied (DNAKit does not download during this call)."
}
```

## 37) generate_sequence_embeddings

- **Purpose:** Embed DNA sequences with a DNAKit representation model via `dnakit.representations.extract_representations`; needs a local checkpoint or a supplied backend.

- **Method:** Delegates to the representation extractor, which pools the DNA language model's hidden states (`pooling`). Without `checkpoint_path` (local weights) or `representation_backend` (object with `extract(sequences, show_progress)`) the call returns `not_implemented`; no weights are downloaded and there is no network access. This is not the k-mer profile (see `export_kmer_embedding` for that). Not run here with real weights.

- **API:** `dnakit.ml_models.generate_sequence_embeddings(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a sequence string, list, records or FASTA/table path. `options`: `checkpoint_path` or `representation_backend` (one is required), `model` (name from `ml_models.search`), `model_source_path`, `pooling`, `ambiguity_policy`, `device`, `dtype`, `batch_size`, `max_length`, `max_sequences` (10000), `sequence_column`, `id_column`.

- **Output:** Without a checkpoint or backend: `status: "not_implemented"` with `required_dependency` and `message`. Otherwise the standard envelope with `result.delegate`, `model`, `pooling`, `dimension`, `record_ids` and `embeddings` (one vector per sequence).

- **Example code:**

```python
import json
from dnakit.ml_models import generate_sequence_embeddings

result = generate_sequence_embeddings(["ACGTACGT", "GGGGCCCC"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "generate_sequence_embeddings",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 24,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "777448ae73191b73a052271ab7cb66335d35637751e950cb0b24eee1ba9c8c73"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pretrained DNA foundation-model weights (options['checkpoint_path']) or options['representation_backend']",
  "message": "Sequence embeddings come from a DNA language model; no local checkpoint or backend was supplied."
}
```

## 38) load

- **Purpose:** Read the manifest (`model.json`) of a model directory saved by `ml_models.train`, without deserialising any weights.

- **Method:** Parses `model.json` (at most 10 MB) and checks whether the artifact file (`estimator.joblib` for machine-learning models, `head.npz` otherwise) exists; the joblib estimator is never loaded, so nothing untrusted is executed. To use the model call `ml_models.predict` with `trust_model=True`.

- **API:** `dnakit.ml_models.load(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the model directory path. No `options`.

- **Output:** The standard envelope. `result` has `model_dir`, `manifest` (task, schema, feature description, training settings, versions), `artifact`, `artifact_present`, `files` (directory listing), `weights_loaded: false` and a `note`. A directory without `model.json` raises `ML_MODEL_NOT_FOUND`.

- **Example code:**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import load, train

train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
result = load("model_out")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "load",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 9,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "74c9778b599cd5ee1b6c3de1b6ea549fe2d68dfe49523bb0e51ad3d99f1100b2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model_dir": "model_out",
    "manifest": {
      "data_mapping": "...",
      "features": "...",
      "finetune_method": "machine_learning",
      "input_dimension": 19,
      "model_name": "logistic_regression",
      "output_dimension": 2,
      "schema_version": "dnakit.custom-property-model.v1",
      "target_schema": "...",
      "task_type": "classification",
      "threshold": 0.5,
      "training": "...",
      "versions": {"joblib": "1.5.3", "numpy": "2.2.6", "python": "3.12.13", "scipy": "1.17.1", "sklearn": "1.8.0"},
      "warnings": []
    },
    "artifact": "estimator.joblib",
    "artifact_present": true,
    "files": ["data_inspection.json", "estimator.joblib", "metrics.json", "model.json", "oof_predictions.jsonl", "splits.csv"],
    "weights_loaded": false,
    "note": "Only the JSON manifest was read; use ml_models.predict with trust_model=True to run the model."
  }
}
```

## 39) mask_dna_tokens

- **Purpose:** Mask a fraction of non-overlapping k-base tokens of DNA sequences for masked-language-model pre-processing.

- **Method:** Each sequence is cut into non-overlapping tokens of `token_size` bases (the last token may be shorter). The number of masked tokens is `round(mask_rate * n_tokens)`, at least 1 when `mask_rate > 0`; the masked positions are sampled without replacement with a random generator seeded with `seed + i` for sequence i, so the result is reproducible. No model is involved.

- **API:** `dnakit.ml_models.mask_dna_tokens(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: as in `encode` (sequence string, FASTA text, list, records, table/FASTA path). `options`: `token_size` (1), `mask_rate` (0.15), `mask_token` (`[MASK]`), `seed` (42), `sequence_column`, `id_column`.

- **Output:** The standard envelope. `result` has `token_size`, `mask_rate`, `mask_token`, `seed` and `sequences`; each sequence has `id`, `tokens`, `masked_tokens`, `masked_positions` and `labels` (original token at each masked position).

- **Example code:**

```python
import json
from dnakit.ml_models import mask_dna_tokens

result = mask_dna_tokens("ACGTACGTAC", options={"token_size": 3, "mask_rate": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "mask_dna_tokens",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "6612c0e673bc5b6636f9147556b461347642dcde65d5bc6b39e209d2b6b95af1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "token_size": 3,
    "mask_rate": 0.5,
    "mask_token": "[MASK]",
    "seed": 42,
    "sequences": [
      {
        "id": "seq1",
        "tokens": ["ACG", "TAC", "GTA", "C"],
        "masked_tokens": ["[MASK]", "TAC", "GTA", "[MASK]"],
        "masked_positions": [0, 3],
        "labels": {"0": "ACG", "3": "C"}
      }
    ]
  }
}
```

## 40) mutagenesis

- **Purpose:** Score every single-base substitution of a DNA sequence with a user-supplied scoring function (in-silico saturation mutagenesis).

- **Method:** For each selected position and each of the three other bases the substituted sequence is scored; `delta = variant score - reference score`, and per position the mean delta and the maximum |delta| are reported. DNAKit supplies no model: the scorer is your own Python function. Without a scorer the call returns `not_implemented`. **This function cannot be used through JSON/MCP**, because `options['scorer']` (or `options['model']`) must be a Python callable, which a JSON call plan cannot carry; call it from Python.

- **API:** `dnakit.ml_models.mutagenesis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one unambiguous DNA sequence (A/C/G/T; an ambiguous base raises `ML_AMBIGUOUS_SEQUENCE`). `options`: `scorer` or `model` (required, Python callable sequence -> number), `positions` (0-based; default all), `max_length` (2000).

- **Output:** The standard envelope. `result` has `sequence_length`, `n_variants`, `reference_score`, `effects` (per variant: `position`, `ref`, `alt`, `score`, `delta`) and `per_position` (`position`, `ref`, `mean_delta`, `max_abs_delta`). Without a scorer: `status: "not_implemented"` with `required_dependency` and `message`.

- **Example code:**

```python
import json
from dnakit.ml_models import mutagenesis

def gc_count(sequence):
    return float(sum(base in "GC" for base in sequence))

result = mutagenesis("AAA", options={"scorer": gc_count, "positions": [1]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "mutagenesis",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 3,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "cb1ad2119d8fafb69566510ee712661f9f14b83385006ef92aec47f523a38358"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_length": 3,
    "n_variants": 3,
    "reference_score": 0.0,
    "effects": [
      {"position": 1, "ref": "A", "alt": "C", "score": 1.0, "delta": 1.0},
      {"position": 1, "ref": "A", "alt": "G", "score": 1.0, "delta": 1.0},
      {"position": 1, "ref": "A", "alt": "T", "score": 0.0, "delta": 0.0}
    ],
    "per_position": [
      {"position": 1, "ref": "A", "mean_delta": 0.6666666666666666, "max_abs_delta": 1.0}
    ]
  }
}
```

## 41) predict

- **Purpose:** Predict with a model saved by `ml_models.train`, via `dnakit.predictions.predict_custom_properties`.

- **Method:** Loads the model directory and predicts for the records of the data file with the same feature definition used in training (class labels and class probabilities, or regression values). Classical models are stored as joblib files, which can execute code when loaded, so they are loaded only with `trust_model=True`: load only models you trust. Unknown option names raise `UNKNOWN_ML_OPTION`.

- **API:** `dnakit.ml_models.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of a data file with the same feature columns. `options`: `model_dir` (required), `trust_model` (must be `True` for joblib models), `mapping` (dict of `DataMapping` fields), `threshold`, `batch_size`, `device`, `output_path` (save predictions to csv/tsv/json/parquet), `overwrite`.

- **Output:** The standard envelope. `result` has `delegate` and `prediction` (`model_dir`, `record_ids`, `label_names`, `predictions` with the predicted label per record, and for classifiers `probabilities` per class).

- **Example code:**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import predict, train

train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
Path("new.csv").write_text("sequence\nAAAACC\nGGGTTT\n")
result = predict("new.csv", options={"model_dir": "model_out", "trust_model": True})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "predict",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 7,
    "sequence_bases": 6,
    "lengths": [6],
    "sha256": "4d4732a2cb66fff9c0c6d137ffa6452fd04bb9ff620556d7615ee7f746331aa2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "delegate": "dnakit.predictions.predict_custom_properties",
    "prediction": {
      "model_dir": "/work/model_out",
      "record_ids": ["record_1", "record_2"],
      "label_names": ["label"],
      "predictions": [
        {"record_id": "record_1", "label": "0"},
        {"record_id": "record_2", "label": "1"}
      ],
      "probabilities": [
        {
          "record_id": "record_1",
          "label": {"0": 0.9232336713763536, "1": 0.07676632862364648}
        },
        {
          "record_id": "record_2",
          "label": {"0": 0.07620323087640403, "1": 0.923796769123596}
        }
      ]
    }
  }
}
```

## 42) pretrain

- **Purpose:** Pre-train a DNA foundation model from scratch; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. Self-supervised pre-training needs a large-scale DNA corpus, a GPU training stack and a foundation-model architecture, none of which DNAKit provides.

- **API:** `dnakit.ml_models.pretrain(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import pretrain

result = pretrain(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "pretrain",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "large-scale DNA corpus, GPU training stack and a foundation-model architecture",
  "message": "Self-supervised pre-training of a DNA foundation model is not implemented in DNAKit."
}
```

## 43) quantize_model

- **Purpose:** Quantise a trained neural network to lower precision; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no neural model object to quantise; this would need a trained neural model and a quantisation backend (`torch.ao` or ONNX Runtime).

- **API:** `dnakit.ml_models.quantize_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import quantize_model

result = quantize_model(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "quantize_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained neural model object and a quantisation backend (torch.ao / ONNX Runtime)",
  "message": "Model quantisation is not implemented; DNAKit has no neural model object to quantise."
}
```

## 44) run_cross_species_transfer

- **Purpose:** Run a cross-species transfer-learning experiment (train on one species, evaluate on another); not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. It needs trained models and species-labelled data with a defined transfer protocol, which DNAKit does not provide as a pipeline. Do it manually with `train` on one species file, then `predict` and `evaluate` on the other.

- **API:** `dnakit.ml_models.run_cross_species_transfer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import run_cross_species_transfer

result = run_cross_species_transfer(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "run_cross_species_transfer",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained models and species-labelled data with a defined transfer protocol",
  "message": "Cross-species transfer needs a pretrained model and a species transfer protocol; DNAKit has no such pipeline (use train and predict on the two species files)."
}
```

## 45) search

- **Purpose:** Search the catalogue of models known to DNAKit: classical ML estimators, fine-tunable DNA models, embedding models and pretrained predictors.

- **Method:** Builds the catalogue from the registries of `dnakit.predictions` and `dnakit.representations`, and keeps the entries in which every word of the query occurs (case-insensitive) in name, kind or description, optionally filtered by `kind` and `task`. No model is loaded or run and nothing is downloaded.

- **API:** `dnakit.ml_models.search(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a query string (empty = all). `options`: `query` (alternative to `input`), `kind` (`machine_learning_estimator`, `dna_model_finetuning`, `dna_embedding_model` or `pretrained_prediction_model`), `task` (`classification`, `regression`, ...), `limit` (50).

- **Output:** The standard envelope. `result` has `query`, `kinds` (all catalogue kinds), `n_matches` and `matches` (per entry: `name`, `kind`, `tasks`, `description`).

- **Example code:**

```python
import json
from dnakit.ml_models import search

result = search("random forest")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "search",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 13,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dc25bdcdacc105874b43917c189e1b76b0dea5a4afe5e9c499e529ea33d48a7d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "query": "random forest",
    "kinds": ["dna_embedding_model", "dna_model_finetuning", "machine_learning_estimator", "pretrained_prediction_model"],
    "n_matches": 1,
    "matches": [
      {
        "name": "random_forest",
        "kind": "machine_learning_estimator",
        "tasks": ["classification", "regression"],
        "description": "scikit-learn estimator"
      }
    ]
  }
}
```

## 46) serve

- **Purpose:** Serve a trained model over a network endpoint; not implemented.

- **Method:** Not implemented in DNAKit: the call computes and writes nothing and always returns `status: "not_implemented"`. DNAKit has no model-serving runtime: the call does not start a server, open a port or load a model. It would need an HTTP/gRPC model server and a deployed model.

- **API:** `dnakit.ml_models.serve(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is accepted and only summarised in `input_summary`; any `options` mapping is accepted and ignored.

- **Output:** `status: "not_implemented"` in the standard envelope (`module`, `operation`, `provenance`, `input_summary`, `warnings`), plus `required_dependency` and `message` naming what is missing (shown below). There is no `result`.

- **Example code:**

```python
import json
from dnakit.ml_models import serve

result = serve(None)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "serve",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a model server runtime (HTTP/gRPC) and a deployed model",
  "message": "DNAKit has no model-serving runtime; ml_models.serve does not start a server."
}
```

## 47) train

- **Purpose:** Train a classical machine-learning property model (scikit-learn estimators) on a labelled sequence/feature file via `dnakit.predictions.train_custom_property_model`.

- **Method:** Delegates to the custom-property trainer with `finetune_method="machine_learning"`: k-mer (or numeric / embedding) features, an estimator chosen by `model` (see `ml_models.search`), optional cross-validation (`cv_folds`) and train/validation/test splits, and a refit on all training rows. It writes the model directory (`model.json`, `estimator.joblib`, `metrics.json`, `splits.csv`, ...). Only classical estimators are trained here; neural fine-tuning is `ml_models.finetune`. Needs scikit-learn. Option names must be `CustomPropertyTrainingConfig` fields (otherwise `UNKNOWN_ML_OPTION`).

- **API:** `dnakit.ml_models.train(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of a CSV/TSV/JSON/FASTA data file. `options`: `task_type` (required, `classification` or `regression`), `output_dir` (required, model directory to write), `model` (default `random_forest`), `model_params`, `feature_source` (`sequence` k-mer features, `numeric`, `embedding`), `kmer_k`, `mapping` (dict of `DataMapping` fields), `cv_folds`, `seed`, split ratios, `overwrite`, and other `CustomPropertyTrainingConfig` fields.

- **Output:** The standard envelope. `result` has `delegate` and `training` (`model_name`, `task_type`, `finetune_method`, `output_dir`, `target_schema`, `split_counts`, `metrics` with cross-validation folds and out-of-fold metrics, `artifacts`, `data_inspection`, `warnings`). `finetune_method` other than `machine_learning` returns `not_implemented`.

- **Example code:**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import train

result = train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "train",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 8,
    "sequence_bases": 7,
    "lengths": [7],
    "sha256": "803d8eb0fe6902e3880bee63c537f65f3917f39e743a36703b093591e7054c79"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "delegate": "dnakit.predictions.train_custom_property_model",
    "training": {
      "model_name": "logistic_regression",
      "task_type": "classification",
      "finetune_method": "machine_learning",
      "output_dir": "/work/model_out",
      "target_schema": "...",
      "split_counts": {"train": 8, "validation": 0, "test": 0},
      "metrics": {
        "cross_validation": {
          "n_splits": 2,
          "splitter": "StratifiedKFold",
          "seed": 42,
          "folds": "...",
          "out_of_fold": {"accuracy": 1.0, "...": "..."},
          "refit_record_count": 8
        }
      },
      "history": [],
      "artifacts": {
        "manifest": "/work/model_out/model.json",
        "estimator": "/work/model_out/estimator.joblib",
        "metrics": "/work/model_out/metrics.json",
        "data_inspection": "/work/model_out/data_inspection.json",
        "splits": "/work/model_out/splits.csv",
        "oof_predictions": "/work/model_out/oof_predictions.jsonl"
      },
      "data_inspection": "...",
      "warnings": []
    }
  }
}
```

## 48) validate

- **Purpose:** Validate a labelled DNA-sequence dataset for machine-learning use: alphabet, empty fields, duplicate ids, conflicting labels and split leakage.

- **Method:** Row checks: characters outside `alphabet` (`INVALID_CHARACTERS`; the message lists the alphabet in the order given, `ACGTN` by default), empty sequence or label, length outside `min_length`..`max_length`, duplicate ids (`DUPLICATE_ID`), and identical sequences (case-insensitive) with different labels (`CONFLICTING_LABELS`). When `split_column` is given, an identical sequence in different splits is reported as `SPLIT_LEAKAGE`. Repeated sequences with consistent labels are only a warning. At most 100 errors are listed.

- **API:** `dnakit.ml_models.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records (list of dicts) or a CSV/TSV/JSON/JSONL/FASTA path. `options`: `sequence_column` (`sequence`), `label_column` (`label`), `id_column` (`id`), `split_column` (none), `alphabet` (`ACGTN`), `min_length` (1), `max_length` (10000000).

- **Output:** The standard envelope. `result` has `valid` (no errors), `errors` and `warnings` (items with `row`, `code`, `message`) and `summary` (`n_records`, `label_counts`, `length_min`, `length_max`, `duplicate_sequence_rows`).

- **Example code:**

```python
import json
from dnakit.ml_models import validate

rows = [
    {"id": "a", "sequence": "ACGT", "label": 1, "split": "train"},
    {"id": "b", "sequence": "ACGT", "label": 1, "split": "test"},
    {"id": "c", "sequence": "AXGT", "label": 0, "split": "train"},
]
result = validate(rows, options={"split_column": "split"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "ml_models",
  "operation": "validate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 188,
    "sequence_bases": 12,
    "lengths": [4, 4, 4],
    "sha256": "7a9decb3d90a6d6d42525cd21ab09db13d04d403102ca0dd1daeaf89dfe74bba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "errors": [
      {
        "row": 2,
        "code": "INVALID_CHARACTERS",
        "message": "characters ['X'] are outside the alphabet ACGTN"
      },
      {
        "row": 0,
        "code": "SPLIT_LEAKAGE",
        "message": "identical sequence appears in several splits (rows [0, 1])"
      }
    ],
    "warnings": [
      {"row": null, "code": "DUPLICATE_SEQUENCES", "message": "1 rows repeat an earlier sequence"}
    ],
    "summary": {
      "n_records": 3,
      "label_counts": {"0": 1, "1": 2},
      "length_min": 4,
      "length_max": 4,
      "duplicate_sequence_rows": 1
    }
  }
}
```
