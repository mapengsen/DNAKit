# gReLU sequence datasets, models and metrics

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: Training and checkpoint/model inference require caller data and weights. Ordinary tests compare local label transforms and metrics; they do not establish pretrained-model accuracy.

## `dnakit.predictions.grelu_p36_functions.augment_sequences`

**Purpose:** Generate the training augmentations of DNA sequences (and labels): shift (jitter), reverse complement

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def augment_sequences(sequences: Sequence[str], *, labels: Sequence[Sequence[Sequence[float]]] | None=None, seq_len: int | None=None, label_len: int | None=None, rc: bool=False, max_seq_shift: int=0, max_pair_shift: int=0, n_mutated_seqs: int=0, n_mutated_bases: int | None=None, protect: Sequence[int]=(), seed: int | None=None, mode: str='serial', n_draws: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Generate the training augmentations of DNA sequences (and labels): shift (jitter), reverse complement
and random mutation (gReLU ``Augmenter``, ``shift``, ``rc_seq``, ``rc_label``).

``rc=True`` adds the reverse-complement copy (and reverses the label along its length); ``max_seq_shift``
shifts the sequence window by up to that many bases either way; ``max_pair_shift`` shifts sequence and label
together; ``n_mutated_seqs`` with ``n_mutated_bases`` makes randomly mutated copies (``protect`` lists
positions to keep). Input sequences must be ``seq_len + 2*(max_seq_shift + max_pair_shift)`` bases long
(``seq_len`` defaults to the unshifted length) and ``labels`` (optional) are tasks x ``label_len +
2*max_pair_shift``. ``mode="serial"`` returns every combination for every sequence, ``mode="random"``
returns ``n_draws`` random draws per sequence. Each output record gives the source index, the augmented
sequence and, when labels were given, the augmented label.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.build_labeled_sequence_dataset`

**Purpose:** Build a labelled sequence dataset and report it, with selected items (gReLU ``LabeledSeqDataset``,

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def build_labeled_sequence_dataset(dataset: Mapping[str, Any], *, indices: Sequence[int]=(0,), verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build a labelled sequence dataset and report it, with selected items (gReLU ``LabeledSeqDataset``,
``DFSeqDataset``, ``BigWigSeqDataset``, ``AnnDataSeqDataset``).

``dataset`` is described in the module description: a ``table`` (DataFrame-style, one label column per
task), ``sequences`` with ``labels`` arrays, ``intervals`` with ``bigwigs`` coverage, or an ``anndata``
file. The result gives the number of sequences, tasks, augmented length and the decoded sequence and
label of every requested ``indices`` entry (indices run over augmented items).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.scan_motif_insertions`

**Purpose:** Insert each motif at every (or the given) position of background sequences (gReLU ``MotifScanDataset``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def scan_motif_insertions(sequences: Sequence[str], motifs: Sequence[str], *, positions: Sequence[int] | None=None, max_records: int=10000, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Insert each motif at every (or the given) position of background sequences (gReLU ``MotifScanDataset``).

Motifs overwrite the bases at the position, so the sequence length does not change. The result lists, for
every sequence, position and motif, the new sequence; pair it with a model prediction
(``predict_sequence_model``) to get a positional motif-effect scan.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.filter_cell_types`

**Purpose:** Drop cell types (tasks) composed of few cells from an AnnData file (gReLU ``filter_cells``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def filter_cell_types(anndata_path: str | Path, output_path: str | Path, *, cutoff: int=1000, count_key: str='n_cells', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Drop cell types (tasks) composed of few cells from an AnnData file (gReLU ``filter_cells``).

``anndata_path`` is an ``.h5ad`` file with intervals in ``.var`` and per-cell-type cell counts in
``.obs[count_key]``; cell types with fewer than ``cutoff`` cells are removed and the result is written to
the new ``output_path``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.build_sequence_model`

**Purpose:** Create an untrained gReLU sequence model and save it as a Lightning checkpoint (``LightningModel``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def build_sequence_model(architecture: str, n_tasks: int, output_path: str | Path, *, model_params: Mapping[str, Any] | None=None, train_params: Mapping[str, Any] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Create an untrained gReLU sequence model and save it as a Lightning checkpoint (``LightningModel``).

``architecture`` selects the network (gReLU ``grelu.model.models``): ``conv`` (ConvModel),
``dilated_conv`` (DilatedConvModel), ``conv_gru`` (ConvGRUModel), ``conv_transformer``
(ConvTransformerModel), ``conv_mlp`` (ConvMLPModel, needs ``seq_len`` in ``model_params``), ``borzoi``
(BorzoiModel), ``enformer`` (EnformerModel), and ``borzoi_pretrained`` / ``enformer_pretrained``
(BorzoiPretrainedModel / EnformerPretrainedModel, which download the published weights and are not run
by the offline tests). ``model_params`` carries the architecture's own arguments (channels, kernel sizes,
``crop_len``, ``final_pool_func`` ...); ``train_params`` the training settings (``task``
binary/multiclass/regression, ``loss``, ``lr``, ``batch_size``, ``max_epochs`` ...). The checkpoint is
the input of the training, evaluation, prediction and embedding functions.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.replace_model_head`

**Purpose:** Replace the output head of a model by a new randomly initialised head with ``n_tasks`` outputs

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def replace_model_head(checkpoint: str | Path, n_tasks: int, output_path: str | Path, *, final_pool_func: str='avg', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Replace the output head of a model by a new randomly initialised head with ``n_tasks`` outputs
(``LightningModel.change_head``) and save the result as a new checkpoint.

``final_pool_func`` (``avg``, ``max`` or ``attn``) is the pooling applied to the new head; the trunk
weights are kept, as the first step of transfer learning.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.train_sequence_model`

**Purpose:** Train a model on a dataset, or fine-tune a pretrained model (``LightningModel.train_on_dataset`` /

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def train_sequence_model(checkpoint: str | Path, train_dataset: Mapping[str, Any], val_dataset: Mapping[str, Any], output_dir: str | Path, *, tune: bool=False, freeze_embedding: bool=False, final_pool_func: str='avg', max_epochs: int | None=None, batch_size: int | None=None, learning_rate: float | None=None, device: str='cpu', seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Train a model on a dataset, or fine-tune a pretrained model (``LightningModel.train_on_dataset`` /
``tune_on_dataset``).

``checkpoint`` is a model from ``build_sequence_model`` (or an earlier training). With ``tune=False`` the
model is trained as it is. With ``tune=True`` a new head with the dataset's number of tasks is built
first (``final_pool_func``), optionally keeping the trunk fixed (``freeze_embedding``), and then trained.
``max_epochs``, ``batch_size`` and ``learning_rate`` override the checkpoint's training settings. Lightning
writes its checkpoints and logs under ``output_dir`` (new or empty); the best and last checkpoint paths
and the final metrics are returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.evaluate_sequence_model`

**Purpose:** Evaluate a model on a labelled test dataset and return one row of metrics per task

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def evaluate_sequence_model(checkpoint: str | Path, dataset: Mapping[str, Any], *, device: str='cpu', batch_size: int=256, output_checkpoint: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Evaluate a model on a labelled test dataset and return one row of metrics per task
(``LightningModel.test_on_dataset``).

The metrics follow the model's ``task`` setting (for example ``test_mse`` and ``test_pearson`` for
regression, AUROC/average precision/accuracy for binary tasks). ``output_checkpoint`` (optional, new
``.ckpt`` path) saves a copy of the model that stores the test data parameters and performance.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.predict_sequence_model`

**Purpose:** Predict with one model, or with several models as an ensemble, on DNA sequences

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def predict_sequence_model(checkpoints: str | Path | Sequence[str | Path], sequences: Sequence[str], *, model_names: Sequence[str] | None=None, rc: bool=False, max_seq_shift: int=0, device: str='cpu', batch_size: int=256, output_path: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Predict with one model, or with several models as an ensemble, on DNA sequences
(``LightningModel.predict_on_dataset`` with ``SeqDataset``; ``LightningModelEnsemble`` for several).

``checkpoints`` is one checkpoint path or a list. For a list the predictions of all models are
concatenated along the task axis in the given order and the tasks are prefixed with ``model_names``
(default ``model0``, ``model1`` ...). The output has shape sequences x tasks x length. ``rc`` and
``max_seq_shift`` average over the reverse-complement and shifted copies (gReLU's ``augment_aggfunc``,
default mean). ``output_path`` (``.npy``) saves large results instead of returning them.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.embed_sequences`

**Purpose:** Return the model's embeddings (the output of its convolutional/attention trunk, before the head) for

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def embed_sequences(checkpoint: str | Path, sequences: Sequence[str], *, device: str='cpu', batch_size: int=256, output_path: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Return the model's embeddings (the output of its convolutional/attention trunk, before the head) for
DNA sequences (``LightningModel.embed_on_dataset``).

The embedding array has shape sequences x channels x length (length after pooling and cropping).
``output_path`` (``.npy``) saves large results instead of returning them.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.compute_prediction_metric`

**Purpose:** Compute a gReLU evaluation metric on predictions and labels, both shaped sequences x tasks x length

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compute_prediction_metric(metric: str, predictions: Sequence[Sequence[Sequence[float]]], targets: Sequence[Sequence[Sequence[float]]], *, average: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Compute a gReLU evaluation metric on predictions and labels, both shaped sequences x tasks x length
(gReLU ``MSE`` and ``BestF1`` metric classes with their ``update`` and ``compute``).

``metric`` is ``mse`` (mean squared error per task) or ``best_f1`` (best F1 over all thresholds per task;
targets must be 0/1 and predictions probabilities). With ``average=True`` the mean over tasks is
returned in ``value``; ``per_task`` always lists the task values.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.transform_labels`

**Purpose:** Clip and/or transform label values (gReLU ``LabelTransform``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def transform_labels(labels: Sequence[Sequence[Sequence[float]]], *, min_clip: float | None=None, max_clip: float | None=None, transform_func: str | None=None, output_path: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Clip and/or transform label values (gReLU ``LabelTransform``).

``labels`` is sequences x tasks x length. Values are first clipped to ``[min_clip, max_clip]``, then
``transform_func`` (``log`` or ``log1p``) is applied.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.score_prediction_specificity`

**Purpose:** Cell-type specificity of model predictions: foreground against background tasks (gReLU ``Specificity``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def score_prediction_specificity(predictions: Sequence[Sequence[Sequence[float]]], on_tasks: Sequence[int], *, off_tasks: Sequence[int] | None=None, on_aggfunc: str='mean', off_aggfunc: str='mean', off_weight: float=1.0, off_thresh: float | None=None, positions: Sequence[int] | None=None, length_aggfunc: str='sum', compare_func: str='divide', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Cell-type specificity of model predictions: foreground against background tasks (gReLU ``Specificity``).

``predictions`` is sequences x tasks x length. The ``on_tasks`` (indices) and ``off_tasks`` (default: all
others) are each aggregated over tasks (``on_aggfunc`` / ``off_aggfunc``: sum, mean, min, max) after
``length_aggfunc`` over ``positions``; the background is multiplied by ``off_weight`` (only values above
``off_thresh`` when given) and combined with ``compare_func`` (``divide``, ``subtract`` or ``log2FC``).
One score per sequence is returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.score_sequence_patterns`

**Purpose:** Weighted count of (overlapping) occurrences of subsequences in each DNA sequence (gReLU ``PatternScore``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def score_sequence_patterns(sequences: Sequence[str], patterns: Sequence[str], *, weights: Sequence[float] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Weighted count of (overlapping) occurrences of subsequences in each DNA sequence (gReLU ``PatternScore``).

Patterns are regular expressions (for example ``ACGT`` or ``TA[CT]A``) matched on the forward strand with
overlaps allowed; ``weights`` (default 1 each) multiply the per-pattern counts, which are summed per
sequence.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.score_sequence_motifs`

**Purpose:** Weighted count of PWM motif occurrences (FIMO-style hits) in each DNA sequence (gReLU ``MotifScore``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def score_sequence_motifs(sequences: Sequence[str], motifs: Mapping[str, Sequence[Sequence[float]]] | str | Path, *, names: Sequence[str] | None=None, weights: Sequence[float] | None=None, pthresh: float=0.001, rc: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Weighted count of PWM motif occurrences (FIMO-style hits) in each DNA sequence (gReLU ``MotifScore``).

``motifs`` is a MEME file path or a mapping from motif name to a 4 x L matrix (rows A, C, G, T);
``names`` selects motifs from a MEME file; hits need a p-value at or below ``pthresh`` and both strands
are scanned when ``rc`` is true. The score of a sequence is the sum over motifs of
``weight * number of hits`` (weights default to 1). gReLU's ``MotifScore.forward`` returns hit counts per
(sequence, motif) row in value-count order; here the hits are aggregated per input sequence, which is
the documented intent of the score.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.validate_genomic_file`

**Purpose:** Check that a path is a valid bigWig or FASTA file by gReLU's rules (``check_bigwig`` / ``check_fasta``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def validate_genomic_file(path: str | Path, kind: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Check that a path is a valid bigWig or FASTA file by gReLU's rules (``check_bigwig`` / ``check_fasta``).

``kind`` is ``bigwig`` (extension ``.bw``, ``.bigwig`` or ``.bigWig`` and an existing file) or ``fasta``
(``.fa``, ``.fasta``, ``.fa.gz``, ``.fasta.gz``, ``.fa.bgz``, ``.fasta.bgz`` and an existing file). The
check is by name and existence, as in gReLU; it does not parse the content.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.grelu_p36_functions.get_grelu_resource_path`

**Purpose:** Locate a data file bundled with gReLU (``get_meme_file_path`` / ``get_blacklist_file``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def get_grelu_resource_path(kind: str, name: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Locate a data file bundled with gReLU (``get_meme_file_path`` / ``get_blacklist_file``).

``kind="motifs"`` with ``name`` ``hocomoco_v12``, ``hocomoco_v13`` or ``consensus`` (JASPAR 2024
consensus motifs) returns the MEME file; ``kind="blacklist"`` with ``name`` ``hg19``, ``hg38`` or
``mm10`` returns the ENCODE blacklist BED file.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.predictions.grelu_p36_functions import transform_labels
r = transform_labels([[[-1., 0., 2.]]], min_clip=0, max_clip=1, verbose=False)
print(r['labels'])
```

**Example output:**

```text
[[[0.0, 0.0, 1.0]]]
```
