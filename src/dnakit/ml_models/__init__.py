"""Real local implementations behind the ml models capability family.

Metrics, calibration, uncertainty, dataset/split/leaderboard exports, sequence encodings, masking and
in-silico mutagenesis are computed natively; training and prediction delegate to
``dnakit.predictions``. Operations that need pretrained foundation-model weights or model objects DNAKit
does not provide report ``status: "not_implemented"``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError
from dnakit.ml_models import _impl_eval as _eval
from dnakit.ml_models import _impl_models as _models
from dnakit.ml_models import _impl_seq as _seq
from dnakit.ml_models._impl_common import Unavailable

_MODULE = "ml_models"


def _run(
    operation: str,
    input: object,
    options: Mapping[str, object] | None,
    implementation: Callable[..., object],
    **extra: object,
) -> object:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping or None.", code="INVALID_ML_OPTIONS")
    opts = dict(options or {})
    try:
        result = implementation(input, opts, **extra)
    except Unavailable as exc:
        return candidate_not_implemented(
            _MODULE, operation, input, options=opts, requires=exc.requires, message=exc.message
        )
    return candidate_ok(_MODULE, operation, input, result, options=opts)


def _unavailable(
    operation: str, input: object, options: Mapping[str, object] | None, requires: str, message: str
) -> object:
    return candidate_not_implemented(
        _MODULE, operation, input, options=dict(options or {}), requires=requires, message=message
    )


def attribute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute per-position in-silico mutagenesis attributions for one DNA sequence with a user-supplied scoring function.

    Definition: attribution[i] = score(sequence) - mean(score of the three single-base substitutions at i).
    Input: one DNA sequence. Required option ``scorer`` (or ``model``): a Python callable sequence -> number
    (otherwise not_implemented). Options: ``positions``, ``max_length`` (2000).
    """
    return _run("attribute", input, options, _seq.attribute)


def build(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Describe the untrained scikit-learn estimator that DNAKit would create for a model name and parameters.

    Input: None or a mapping of the options below. Options: ``task_type`` (required), ``model`` (default
    'random_forest'), ``model_params``. Returns the estimator class and its full parameter dictionary;
    nothing is fitted.
    """
    return _run("build", input, options, _models.build)


def calculate_prediction_uncertainty(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Compute predictive uncertainty: entropy/margin for class probabilities, or entropy decomposition / variance for an ensemble.

    Input: a mapping with one of ``probabilities`` (samples x classes: entropy, least confidence, margin;
    members x samples x classes: predictive entropy, expected entropy, mutual information = predictive -
    expected entropy) or ``predictions`` (members x samples regression values: mean, variance, std).
    Options: ``log_base`` (e), ``ddof`` (0; 1 for sample variance across ensemble members).
    """
    return _run("calculate_prediction_uncertainty", input, options, _eval.uncertainty)


def calibrate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Fit a probability calibrator (Platt scaling or isotonic regression) on binary scores and return the calibrated probabilities and the expected calibration error.

    Input: a mapping of equally long columns ``{'y_true': [...], 'y_pred': [...], 'y_score': [...]}``,
    a list of such records, or a CSV/TSV/JSON/JSONL path. Column names can be changed with the options
    ``true_column``, ``pred_column`` and ``score_column``. ``y_score`` is required.
    Options: ``method`` ('platt' default | 'isotonic'), ``positive_label`` (default the larger label),
    ``apply_to`` (scores to transform; default the fitting scores), ``n_bins`` (10, for ECE). Platt uses
    the Lin-Lin-Weng smoothed-target Newton fit ``p = 1/(1+exp(A*s+B))``; isotonic uses pool-adjacent-
    violators with linear interpolation and clipping. ECE is measured on the fitting data (resubstitution).
    ``ece_before`` is ``null`` when any fitting score lies outside [0, 1] (ECE needs probabilities); ``ece_after_on_fit_data``
    is always reported.
    """
    return _run("calibrate", input, options, _eval.calibrate)


def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Check prediction columns for equal lengths, missing values, non-finite scores and label consistency before evaluation.

    Input: a mapping of equally long columns ``{'y_true': [...], 'y_pred': [...], 'y_score': [...]}``,
    a list of such records, or a CSV/TSV/JSON/JSONL path. Column names can be changed with the options
    ``true_column``, ``pred_column`` and ``score_column``.
    Options: ``task`` ('auto' | 'classification' | 'regression'). Returns ``valid``, ``errors`` and
    ``warnings`` (for example scores outside [0, 1] or predicted labels never seen in y_true).
    """
    return _run("check", input, options, _eval.check_predictions)


def encode(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Encode DNA sequences as one-hot matrices, integer codes or k-mer vectors using the DNAKit fingerprint encoders.

    Input: a sequence string, FASTA text, a list of sequences, records, or a table/FASTA path. Options:
    ``method`` ('one_hot' | 'integer' | 'kmer'), ``k`` (3, kmer), ``mode`` ('count' | 'frequency' |
    'binary'), ``canonical`` (False), ``ambiguity_policy`` (one_hot 'error'|'zero'|'fractional'; integer
    'iupac'|'sentinel'; kmer 'ignore'|'error'), ``sequence_column`` / ``id_column``, ``max_sequences`` (1000),
    ``max_output_length`` (100000 per sequence).
    """
    return _run("encode", input, options, _seq.encode)


def evaluate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute classification or regression metrics (accuracy, balanced accuracy, precision, recall, F1, MCC, kappa, AUROC, AUPRC, log loss; MAE, MSE, RMSE, R2, Pearson, Spearman) from predictions.

    Input: a mapping of equally long columns ``{'y_true': [...], 'y_pred': [...], 'y_score': [...]}``,
    a list of such records, or a CSV/TSV/JSON/JSONL path. Column names can be changed with the options
    ``true_column``, ``pred_column`` and ``score_column``.
    Options: ``task`` ('auto' | 'classification' | 'regression'; auto = regression when y_true has
    non-integer numbers), ``positive_label`` (binary, default the larger label), ``threshold`` (0.5, used
    only when y_pred is absent and y_score is given), ``classes`` (label order), ``metrics`` (subset of
    metric names to return). AUROC/AUPRC (average precision) need a binary y_score. Multiclass MCC is
    Gorodkin's R_K; undefined correlations are None.
    """
    return _run("evaluate", input, options, _eval.evaluate)


def explain(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Explain a model with permutation feature importance on tabular data (drop in R2 or accuracy when one column is shuffled).

    Input: ``{'X': [[...]], 'y': [...], 'feature_names': [...]}``. Required option ``model``: a callable (or object
    with ``predict``) mapping a list of feature rows to predictions; otherwise not_implemented. ``options['scorer']`` is accepted as an alias of ``model`` (``model`` wins when both are given). Options:
    ``task`` ('regression' | 'classification'), ``n_repeats`` (5), ``seed`` (42).
    """
    return _run("explain", input, options, _seq.permutation_importance)


def export_benchmark(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Aggregate model-by-dataset benchmark results into per-dataset ranks and a model summary (mean value, mean rank, wins) and write them.

    Input: long-format records ``{'model', 'dataset', 'value'}`` or a table path. Options: ``model_column``,
    ``dataset_column``, ``value_column`` (defaults 'model', 'dataset', 'value'), ``metric`` (name used to
    decide direction), ``higher_is_better``, ``format`` ('json' default: per_dataset + summary; 'tsv'/'csv'
    write the per-dataset table), ``output``, ``overwrite``. Duplicate model/dataset pairs are rejected.
    """
    return _run("export_benchmark", input, options, _eval.benchmark)


def export_calibration_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Fit a Platt or isotonic calibrator on binary scores and write it as a JSON model file.

    Input: a mapping of equally long columns ``{'y_true': [...], 'y_pred': [...], 'y_score': [...]}``,
    a list of such records, or a CSV/TSV/JSON/JSONL path. Column names can be changed with the options
    ``true_column``, ``pred_column`` and ``score_column``. ``y_score`` is required.
    Options: ``method`` ('platt' | 'isotonic'), ``positive_label``, ``output`` (JSON path; without it the
    JSON text is returned under ``content``), ``overwrite`` (False).
    """
    return _run("export_calibration_model", input, options, _eval.export_calibration)


def export_chromatin_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained chromatin model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no chromatin model class to export; nothing was trained or serialised. Needs: a trained chromatin model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_chromatin_model",
        input,
        options,
        "a trained chromatin model object and its architecture-specific serializer",
        "DNAKit has no chromatin model class to export; nothing was trained or serialised.",
    )


def export_cnn_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained cnn model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no cnn model class to export; nothing was trained or serialised. Needs: a trained cnn model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_cnn_model",
        input,
        options,
        "a trained cnn model object and its architecture-specific serializer",
        "DNAKit has no cnn model class to export; nothing was trained or serialised.",
    )


def export_dataset(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Write a dataset of records to CSV, TSV, JSON, JSONL or FASTA.

    Input: records, a column mapping, or a CSV/TSV/JSON/JSONL/FASTA path. Options: ``format`` (from the
    ``output`` suffix, default 'csv'), ``columns`` (subset and order), ``output`` (path; otherwise the
    text is returned under ``content``), ``overwrite`` (False). FASTA needs a ``sequence`` column.
    """
    return _run("export_dataset", input, options, _eval.export_dataset)


def export_expression_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Export a trained expression model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no expression model class to export; nothing was trained or serialised. Needs: a trained expression model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_expression_model",
        input,
        options,
        "a trained expression model object and its architecture-specific serializer",
        "DNAKit has no expression model class to export; nothing was trained or serialised.",
    )


def export_graph_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained graph model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no graph model class to export; nothing was trained or serialised. Needs: a trained graph model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_graph_model",
        input,
        options,
        "a trained graph model object and its architecture-specific serializer",
        "DNAKit has no graph model class to export; nothing was trained or serialised.",
    )


def export_interpretability(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Export per-position attributions (supplied, or computed by in-silico mutagenesis with a scorer) as TSV or JSON.

    Input: ``{'sequence': 'ACGT', 'attributions': [...]}`` or one DNA sequence plus ``options['scorer']`` (a
    callable sequence -> number). Options: ``format`` ('tsv' | 'json'), ``output``, ``overwrite``.
    """
    return _run("export_interpretability", input, options, _seq.export_interpretability)


def export_kmer_embedding(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute k-mer frequency (or count) vectors for sequences and write them as a TSV matrix or JSON.

    Input as in encode. Options: ``k`` (3), ``mode`` ('frequency' default | 'count' | 'binary'),
    ``canonical`` (False), ``ambiguity_policy`` ('ignore'), ``format`` ('tsv' | 'json'), ``output``,
    ``overwrite``. This is the exact k-mer profile, not a learned embedding.
    """
    return _run("export_kmer_embedding", input, options, _seq.export_kmer)


def export_leaderboard(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Rank models by one metric (standard competition ranking) and write the leaderboard as TSV, CSV, JSON or JSONL.

    Input: records such as ``[{'model': 'A', 'score': 0.9}, ...]`` or a table path. Options: ``metric``
    ('score'), ``model_column`` ('model'), ``higher_is_better`` (auto: False for names containing loss, error,
    rmse, mae, mse, ece, nll, perplexity), ``format`` ('tsv' default, or from the ``output`` suffix),
    ``output`` (path; otherwise returned inline), ``overwrite`` (False).
    """
    return _run("export_leaderboard", input, options, _eval.leaderboard)


def export_long_context_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Export a trained long context model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no long context model class to export; nothing was trained or serialised. Needs: a trained long context model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_long_context_model",
        input,
        options,
        "a trained long context model object and its architecture-specific serializer",
        "DNAKit has no long context model class to export; nothing was trained or serialised.",
    )


def export_model_onnx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained model to the ONNX format; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules. Needs: pip install onnx onnxruntime and a trained model object.
    """
    return _unavailable(
        "export_model_onnx",
        input,
        options,
        "pip install onnx onnxruntime and a trained model object",
        "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules.",
    )


def export_onehot_embedding(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """One-hot encode sequences (columns A, C, G, T) and write them as JSON or long-format TSV (id, position, A, C, G, T).

    Input as in encode. Options: ``format`` ('json' | 'tsv', default from the ``output`` suffix),
    ``ambiguity_policy`` ('error' | 'zero' | 'fractional'), ``output``, ``overwrite``, ``max_sequences``.
    """
    return _run("export_onehot_embedding", input, options, _seq.export_onehot)


def export_onnx_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained model to the ONNX format; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules. Needs: pip install onnx onnxruntime and a trained model object.
    """
    return _unavailable(
        "export_onnx_model",
        input,
        options,
        "pip install onnx onnxruntime and a trained model object",
        "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules.",
    )


def export_pangenome_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained pangenome model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no pangenome model class to export; nothing was trained or serialised. Needs: a trained pangenome model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_pangenome_model",
        input,
        options,
        "a trained pangenome model object and its architecture-specific serializer",
        "DNAKit has no pangenome model class to export; nothing was trained or serialised.",
    )


def export_phenotype_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained phenotype model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no phenotype model class to export; nothing was trained or serialised. Needs: a trained phenotype model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_phenotype_model",
        input,
        options,
        "a trained phenotype model object and its architecture-specific serializer",
        "DNAKit has no phenotype model class to export; nothing was trained or serialised.",
    )


def export_protein_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained protein model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no protein model class to export; nothing was trained or serialised. Needs: a trained protein model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_protein_model",
        input,
        options,
        "a trained protein model object and its architecture-specific serializer",
        "DNAKit has no protein model class to export; nothing was trained or serialised.",
    )


def export_rnn_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained rnn model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no rnn model class to export; nothing was trained or serialised. Needs: a trained rnn model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_rnn_model",
        input,
        options,
        "a trained rnn model object and its architecture-specific serializer",
        "DNAKit has no rnn model class to export; nothing was trained or serialised.",
    )


def export_splice_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a trained splice model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no splice model class to export; nothing was trained or serialised. Needs: a trained splice model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_splice_model",
        input,
        options,
        "a trained splice model object and its architecture-specific serializer",
        "DNAKit has no splice model class to export; nothing was trained or serialised.",
    )


def export_test_split(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Write the test part of a deterministic train/validation/test split of a dataset.

    Input: records or a table path. If every record has a ``split`` value (column ``split_column``,
    default 'split'; train/validation|valid|val|dev/test) it is used. Otherwise a seeded random split is made:
    ``train_ratio`` 0.8, ``validation_ratio`` 0.1, ``test_ratio`` 0.1 (sum 1), ``seed`` 42, optional
    ``group_column`` keeping each group in one split; counts use largest-remainder allocation. Other options:
    ``format``, ``columns``, ``output``, ``overwrite`` as in export_dataset.
    """
    return _run("export_test_split", input, options, _eval.export_split, which="test")


def export_tokenizer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a k-mer DNA tokenizer vocabulary (special tokens, then all 4^k k-mers in A,C,G,T order) and write it as JSON.

    Input is ignored (may be None). Options: ``k`` (3, max 8), ``stride`` (1), ``special_tokens``
    (['[PAD]', '[UNK]', '[CLS]', '[SEP]', '[MASK]']), ``output`` (JSON path; otherwise returned), ``overwrite``.
    """
    return _run("export_tokenizer", input, options, lambda _value, opts: _seq.tokenizer(opts))


def export_torchscript_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Trace a PyTorch module with an example input and save it as a TorchScript archive, verifying the reloaded output.

    Input: a ``torch.nn.Module`` (or ``options['model']``). Required options: ``example_input`` (nested list
    of numbers) and ``output`` (.pt path). ``overwrite`` (False). Needs torch (otherwise not_implemented).
    Returns the file size and the maximum absolute difference between original and reloaded outputs.
    """
    return _run("export_torchscript_model", input, options, _models.export_torchscript)


def export_train_split(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Write the train part of a deterministic train/validation/test split of a dataset.

    Input: records or a table path. If every record has a ``split`` value (column ``split_column``,
    default 'split'; train/validation|valid|val|dev/test) it is used. Otherwise a seeded random split is made:
    ``train_ratio`` 0.8, ``validation_ratio`` 0.1, ``test_ratio`` 0.1 (sum 1), ``seed`` 42, optional
    ``group_column`` keeping each group in one split; counts use largest-remainder allocation. Other options:
    ``format``, ``columns``, ``output``, ``overwrite`` as in export_dataset.
    """
    return _run("export_train_split", input, options, _eval.export_split, which="train")


def export_transformer_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Export a trained transformer model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no transformer model class to export; nothing was trained or serialised. Needs: a trained transformer model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_transformer_model",
        input,
        options,
        "a trained transformer model object and its architecture-specific serializer",
        "DNAKit has no transformer model class to export; nothing was trained or serialised.",
    )


def export_conformal_uncertainty_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Fit a split-conformal uncertainty model on a held-out calibration set and write it as JSON.

    Input: a mapping ``{'y_true': [...], 'y_pred': [...]}`` (regression, score |y - y_hat|) or
    ``{'y_true': [...], 'probabilities': [[...]], 'classes': [...]}`` (classification, score 1 - p(true class)).
    Options: ``alpha`` (0.1), ``output`` (JSON path; otherwise returned inline), ``overwrite`` (False). The stored
    quantile is the ceil((n+1)(1-alpha))-th smallest score (null when the calibration set is too small).
    """
    return _run("export_conformal_uncertainty_model", input, options, _eval.export_uncertainty)


def export_validation_split(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Write the validation part of a deterministic train/validation/test split of a dataset.

    Input: records or a table path. If every record has a ``split`` value (column ``split_column``,
    default 'split'; train/validation|valid|val|dev/test) it is used. Otherwise a seeded random split is made:
    ``train_ratio`` 0.8, ``validation_ratio`` 0.1, ``test_ratio`` 0.1 (sum 1), ``seed`` 42, optional
    ``group_column`` keeping each group in one split; counts use largest-remainder allocation. Other options:
    ``format``, ``columns``, ``output``, ``overwrite`` as in export_dataset.
    """
    return _run("export_validation_split", input, options, _eval.export_split, which="validation")


def export_variant_effect_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Export a trained variant effect model to a portable file; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no variant effect model class to export; nothing was trained or serialised. Needs: a trained variant effect model object and its architecture-specific serializer.
    """
    return _unavailable(
        "export_variant_effect_model",
        input,
        options,
        "a trained variant effect model object and its architecture-specific serializer",
        "DNAKit has no variant effect model class to export; nothing was trained or serialised.",
    )


def finetune(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Fine-tune a DNA model (linear probe or full) on a labelled file via dnakit.predictions.train_custom_property_model, given local pretrained weights.

    Input: path of a data file. Required (else not_implemented): ``checkpoint_path`` (local weights),
    ``representation_backend`` or ``training_encoder``. Required options: ``task_type`` and ``output_dir``;
    ``finetune_method`` ('linear_probe' default | 'full'); other keys are ``CustomPropertyTrainingConfig`` fields.
    No weights are downloaded.
    """
    return _run("finetune", input, options, _models.finetune)


def generate_sequence_embeddings(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Embed DNA sequences with a DNAKit representation model via dnakit.representations.extract_representations (needs a local checkpoint or a supplied backend).

    Input: sequence string, list, records or FASTA/table path. Required (else not_implemented): a local
    ``checkpoint_path`` or ``representation_backend`` (object with ``extract(sequences, show_progress)``).
    Options: ``model`` (name from ml_models.search), ``pooling``, ``ambiguity_policy``, ``device``, ``dtype``,
    ``batch_size``, ``max_length``, ``max_sequences``. No network access.
    """
    return _run("generate_sequence_embeddings", input, options, _models.embeddings)


def load(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read the manifest (model.json) of a model directory saved by ml_models.train, without deserialising any weights.

    Input: the model directory path. Returns the manifest (task, schema, metrics), whether the artifact file
    exists and the directory listing.
    """
    return _run("load", input, options, _models.load)


def mask_dna_tokens(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Mask a fraction of non-overlapping k-base tokens of DNA sequences for masked-language-model pre-processing.

    Input as in encode. Options: ``token_size`` (1), ``mask_rate`` (0.15; the number of masked tokens is
    round(rate * n_tokens), at least 1 when rate > 0), ``mask_token`` ('[MASK]'), ``seed`` (42; sequence i uses
    seed + i). Returns tokens, masked tokens, masked positions and the original token at each position.
    """
    return _run("mask_dna_tokens", input, options, _seq.mask_tokens)


def mutagenesis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Score every single-base substitution of a DNA sequence with a user-supplied scoring function (in-silico saturation mutagenesis).

    Input: one DNA sequence (A/C/G/T). Required option ``scorer`` (or ``model``): a Python callable mapping a
    sequence string to a number; without it the call reports not_implemented. Options: ``positions`` (0-based,
    default all), ``max_length`` (2000). Returns the reference score, each variant's score and delta
    (variant - reference) and per-position mean/max |delta|.
    """
    return _run("mutagenesis", input, options, _seq.mutagenesis)


def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict with a model saved by ml_models.train via dnakit.predictions.predict_custom_properties.

    Input: path of a data file with the same feature columns. Required option ``model_dir``. Options:
    ``trust_model`` (must be True for joblib ML artifacts: only load models you trust), ``mapping``,
    ``threshold``, ``batch_size``, ``output_path`` (save predictions to csv/tsv/json/parquet), ``overwrite``.
    """
    return _run("predict", input, options, _models.predict)


def pretrain(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Pre-train a DNA foundation model from scratch; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: Self-supervised pre-training of a DNA foundation model is not implemented in DNAKit. Needs: large-scale DNA corpus, GPU training stack and a foundation-model architecture.
    """
    return _unavailable(
        "pretrain",
        input,
        options,
        "large-scale DNA corpus, GPU training stack and a foundation-model architecture",
        "Self-supervised pre-training of a DNA foundation model is not implemented in DNAKit.",
    )


def quantize_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quantise a trained neural network to lower precision; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: Model quantisation is not implemented; DNAKit has no neural model object to quantise. Needs: a trained neural model object and a quantisation backend (torch.ao / ONNX Runtime).
    """
    return _unavailable(
        "quantize_model",
        input,
        options,
        "a trained neural model object and a quantisation backend (torch.ao / ONNX Runtime)",
        "Model quantisation is not implemented; DNAKit has no neural model object to quantise.",
    )


def run_cross_species_transfer(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Run a cross-species transfer-learning experiment (train on one species, evaluate on another); DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: Cross-species transfer needs a pretrained model and a species transfer protocol; DNAKit has no such pipeline (use train and predict on the two species files). Needs: trained models and species-labelled data with a defined transfer protocol.
    """
    return _unavailable(
        "run_cross_species_transfer",
        input,
        options,
        "trained models and species-labelled data with a defined transfer protocol",
        "Cross-species transfer needs a pretrained model and a species transfer protocol; DNAKit has no such pipeline (use train and predict on the two species files).",
    )


def search(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Search the catalogue of models known to DNAKit: classical ML estimators, fine-tunable DNA models, embedding models and pretrained predictors.

    Input: a query string (all words must occur in name, kind or description; empty = all). Options: ``query``
    (alternative to input), ``kind`` ('machine_learning_estimator' | 'dna_model_finetuning' |
    'dna_embedding_model' | 'pretrained_prediction_model'), ``task`` ('classification' | 'regression' | ...),
    ``limit`` (50). Returns matching catalogue entries; no model is run.
    """
    return _run("search", input, options, _models.search)


def serve(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serve a trained model over a network endpoint; DNAKit does not implement this, so the call always returns status 'not_implemented' and computes nothing.

    Reason: DNAKit has no model-serving runtime; ml_models.serve does not start a server. Needs: a model server runtime (HTTP/gRPC) and a deployed model.
    """
    return _unavailable(
        "serve",
        input,
        options,
        "a model server runtime (HTTP/gRPC) and a deployed model",
        "DNAKit has no model-serving runtime; ml_models.serve does not start a server.",
    )


def train(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Train a classical machine-learning property model (scikit-learn estimators) on a labelled sequence/feature file via dnakit.predictions.train_custom_property_model.

    Input: path of a CSV/TSV/JSON/FASTA data file. Required options: ``task_type`` ('classification' |
    'regression') and ``output_dir`` (model directory to write). Other options are the fields of
    ``CustomPropertyTrainingConfig``: ``model`` (default 'random_forest'; see ml_models.search),
    ``model_params``, ``feature_source`` ('sequence' k-mer features | 'numeric' | 'embedding'), ``kmer_k``,
    ``mapping`` (dict of DataMapping fields), ``cv_folds``, ``seed``, split ratios, ``overwrite`` ...
    Only finetune_method='machine_learning' is accepted here; neural fine-tuning is ml_models.finetune.
    """
    return _run("train", input, options, _models.train)


def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a labelled DNA-sequence dataset for machine-learning use (alphabet, empty fields, duplicate ids, conflicting labels, split leakage).

    Input: records (list of dicts) or a CSV/TSV/JSON/JSONL/FASTA path. Options: ``sequence_column``
    ('sequence'), ``label_column`` ('label'), ``id_column`` ('id'), ``split_column`` (None; when set,
    identical sequences in different splits are reported as leakage), ``alphabet`` ('ACGTN'),
    ``min_length`` (1), ``max_length`` (10000000); the alphabet is reported in the order given, e.g. 'ACGTN'. Returns ``valid`` (no errors), up to 100 ``errors`` and
    ``warnings``, and a ``summary`` with label counts and length range.
    """
    return _run("validate", input, options, _eval.validate_dataset)


__all__ = [
    "attribute",
    "build",
    "calculate_prediction_uncertainty",
    "calibrate",
    "check",
    "encode",
    "evaluate",
    "explain",
    "export_benchmark",
    "export_calibration_model",
    "export_chromatin_model",
    "export_cnn_model",
    "export_conformal_uncertainty_model",
    "export_dataset",
    "export_expression_model",
    "export_graph_model",
    "export_interpretability",
    "export_kmer_embedding",
    "export_leaderboard",
    "export_long_context_model",
    "export_model_onnx",
    "export_onehot_embedding",
    "export_onnx_model",
    "export_pangenome_model",
    "export_phenotype_model",
    "export_protein_model",
    "export_rnn_model",
    "export_splice_model",
    "export_test_split",
    "export_tokenizer",
    "export_torchscript_model",
    "export_train_split",
    "export_transformer_model",
    "export_validation_split",
    "export_variant_effect_model",
    "finetune",
    "generate_sequence_embeddings",
    "load",
    "mask_dna_tokens",
    "mutagenesis",
    "predict",
    "pretrain",
    "quantize_model",
    "run_cross_species_transfer",
    "search",
    "serve",
    "train",
    "validate",
]

__candidate_adapter_operations__ = set(
    [
        "attribute",
        "build",
        "calculate_prediction_uncertainty",
        "calibrate",
        "check",
        "encode",
        "evaluate",
        "explain",
        "export_benchmark",
        "export_calibration_model",
        "export_chromatin_model",
        "export_cnn_model",
        "export_dataset",
        "export_expression_model",
        "export_graph_model",
        "export_interpretability",
        "export_kmer_embedding",
        "export_leaderboard",
        "export_long_context_model",
        "export_model_onnx",
        "export_onehot_embedding",
        "export_onnx_model",
        "export_pangenome_model",
        "export_phenotype_model",
        "export_protein_model",
        "export_rnn_model",
        "export_splice_model",
        "export_test_split",
        "export_tokenizer",
        "export_torchscript_model",
        "export_train_split",
        "export_transformer_model",
        "export_conformal_uncertainty_model",
        "export_validation_split",
        "export_variant_effect_model",
        "finetune",
        "generate_sequence_embeddings",
        "load",
        "mask_dna_tokens",
        "mutagenesis",
        "predict",
        "pretrain",
        "quantize_model",
        "run_cross_species_transfer",
        "search",
        "serve",
        "train",
        "validate",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
