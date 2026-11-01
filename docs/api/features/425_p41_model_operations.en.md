# P41 model operations

Public checks use Kipoi 0.8.6, genomic-benchmarks 1.0.0 and PyTorch 2.14.1. Kipoi and CNN checks compare real local small models with native upstream calls. Evo2 uses contract tests only. Real TensorFlow export and Evo2 pretrained inference remain NOT_RUN. Independent acceptance has not been performed.

## Kipoi activation and gradients

**Purpose:** Extract intermediate activations and input gradients.

**Method:** Load a local model with get_model and call predict_activation_on_batch or input_grad. input_dtype defaults to float32 for JSON arrays; None preserves supplied dtype. Native model restrictions on pre_nonlinearity still apply.

**API:** Module `dnakit.predictions.kipoi_p41_functions`:

- `kipoi_layer_activations(model, inputs, layer, *, source='dir', pre_nonlinearity=False, input_dtype='float32', verbose=True)`
- `kipoi_input_gradients(model, inputs, *, source='dir', filter_idx=None, avg_func=None, layer=None, final_layer=True, pre_nonlinearity=False, input_dtype='float32', verbose=True)`

**Input:** A trusted local model directory, arrays or nested dictionaries and a layer selector. Changing source may trigger upstream model lookup.

**Output:** JSON activations or gradients and Kipoi version. Native PyTorchModel comparisons use one-hot, uniform and all-zero input compositions.

## Bounded Kipoi dataloader

**Purpose:** Read bounded raw, train or predict batches.

**Method:** Call native iterators, disabling cycling for train. One additional batch detects truncation.

**API:** `dnakit.predictions.kipoi_p41_functions.kipoi_dataloader_batches(dataloader, *, source='dir', dataloader_args=None, mode='predict', batch_size=32, max_batches=1, verbose=True)`

**Input:** A trusted local loader and its arguments; batch_size 1–10000, max_batches 1–100.

**Output:** batches, batch_count and truncated. Direct checks cover all three modes and 1, 4 and 5 records.

## Prediction BED and per-base tracks

**Purpose:** Write BED6 with predictions, bedGraph or bigWig.

**Method:** Call native BedBatchWriter or region_write. Upstream sums channels of two-dimensional bedGraph data.

**API:** Module `dnakit.predictions.kipoi_p41_functions`:

- `kipoi_write_bed_predictions(ranges, predictions, output_path, *, header=True, verbose=True)`
- `kipoi_write_region_track(regions, values, output_path, *, output_format='bedgraph', chrom_sizes=None, verbose=True)`

**Input:** Parallel BED arrays chr/start/end/id/strand or regions with matching finite per-base values. bigWig requires chromosome sizes and sorted intervals. Output paths must be new.

**Output:** Actual output paths. Checks compare text contents and decoded bigWig values with upstream.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.predictions.kipoi_p41_functions import kipoi_write_region_track
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / 'track.bedgraph'
    result = kipoi_write_region_track(
        [{'chr': 'chr1', 'start': 0, 'end': 2}],
        [[0.1, 0.2]], str(path), verbose=False)
    print(result['status'], result['region_count'])
    print(len(path.read_text().splitlines()))
```

**Example output:**

```text
completed 1
2
```

## genomic-benchmarks CNN

**Purpose:** Train, evaluate and predict with a local DNA token CNN, or export an initialised TensorFlow architecture.

**Method:** Call native CNN.fit/test/forward on CPU. Training saves state_dict and loading uses weights_only=True. The upstream binary sigmoid and binary_cross_entropy_with_logits combination is preserved. TensorFlow export calls get_basic_cnn_model_v0 without training.

**API:** Module `dnakit.predictions.genomic_benchmarks_p41_functions`:

- `train_genomic_benchmark_cnn(token_ids, labels, config, output_checkpoint, *, epochs=1, batch_size=32, seed=0, verbose=True)`
- `evaluate_genomic_benchmark_cnn(token_ids, labels, config, checkpoint, *, batch_size=32, verbose=True)`
- `predict_genomic_benchmark_cnn(token_ids, config, checkpoint, *, verbose=True)`
- `build_genomic_benchmark_tf_cnn(number_of_classes, vocab_size, output_path, *, input_len=100, verbose=True)`

**Input:** Explicit vocabulary indices and config containing number_of_classes, vocab_size, embedding_dim and input_len.

**Output:** Checkpoint paths, native accuracy/F1/loss at printed precision, or scores. Fixed-seed checks cover mixed DNA tokens, homopolymers and multiclass inputs. These checks do not establish biological predictive quality. TensorFlow and tensorflow-addons are absent; real Keras export is unverified.

**Example code:**

```python
import tempfile
from pathlib import Path
import torch
from dnakit.predictions.genomic_benchmarks_p41_functions import train_genomic_benchmark_cnn, predict_genomic_benchmark_cnn
torch.set_num_threads(1)
config = {'number_of_classes': 2, 'vocab_size': 5, 'embedding_dim': 2, 'input_len': 100}
tokens = [[0] * 100, [1] * 100, [2] * 100]
with tempfile.TemporaryDirectory() as folder:
    checkpoint = str(Path(folder) / 'cnn.pth')
    trained = train_genomic_benchmark_cnn(tokens, [0, 1, 0], config, checkpoint, batch_size=3, verbose=False)
    predicted = predict_genomic_benchmark_cnn(tokens, config, checkpoint, verbose=False)
    print(trained['status'], predicted['status'])
    print(len(predicted['scores']), len(predicted['scores'][0]))
```

**Example output:**

```text
completed completed
3 1
```

## Local Evo2 inference

**Purpose:** Sequence likelihood, forward logits/embeddings and positional entropy or perplexity.

**Method:** Require explicit local weights without download. Preserve case. Scoring optionally averages reverse complements. Positional averaging occurs in entropy space before exponentiating for perplexity. Forward output has max_output_values bounds.

**API:** Module `dnakit.predictions.evo2_p41_functions`:

- `score_evo2_sequences(sequences, checkpoint, *, model_name='evo2_7b', batch_size=1, prepend_bos=False, reduce_method='mean', average_reverse_complement=False, verbose=True)`
- `forward_evo2_sequences(sequences, checkpoint, *, model_name='evo2_7b', layer_names=None, prepend_bos=False, device='cuda:0', max_output_values=1000000, verbose=True)`
- `evo2_positional_uncertainty(sequence, checkpoint, *, model_name='evo2_7b', metric='entropy', reverse_complement=False, prepend_bos=False, device='cuda:0', verbose=True)`

**Input:** DNA strings of 2–100000 A/C/G/T/N bases, local checkpoint and compatible Evo2/Vortex/CUDA.

**Output:** JSON arrays with version. Contracts cover flags, reverse-position alignment, output limits and missing checkpoints. Real pretrained numerical output remains unverified.

**Example code:**

```python
from dnakit.predictions.evo2_p41_functions import score_evo2_sequences
from dnakit.exceptions import ConfigurationError
try:
    score_evo2_sequences(['AC-X'], 'unused.ckpt', verbose=False)
except ConfigurationError:
    print('invalid DNA rejected before checkpoint loading')
```

**Example output:**

```text
invalid DNA rejected before checkpoint loading
```
