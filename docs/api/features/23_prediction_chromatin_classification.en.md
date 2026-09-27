# Epigenetic Marks and Open Chromatin

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) H2A.Z histone-variant-region classification {#28-h2az-histone-variant-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H2A.Z (H2AFZ) ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h2afz"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 2) H3K27ac histone-mark-region classification {#29-h3k27ac-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K27ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k27ac"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 3) H3K27me3 histone-mark-region classification {#30-h3k27me3-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K27me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k27me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 4) H3K36me3 histone-mark-region classification {#31-h3k36me3-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K36me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k36me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 5) H3K4me1 histone-mark-region classification {#32-h3k4me1-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K4me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me1"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 6) H3K4me2 histone-mark-region classification {#33-h3k4me2-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K4me2 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me2"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 7) H3K4me3 histone-mark-region classification {#34-h3k4me3-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K4me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k4me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 8) H3K9ac histone-mark-region classification {#35-h3k9ac-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K9ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k9ac"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 9) H3K9me3 histone-mark-region classification {#36-h3k9me3-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H3K9me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h3k9me3"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (previous local checkpoint CUDA smoke run):** Class names, probabilities and the predicted label; values depend on the sequence and checkpoint.

```text
('0', '1')
[0.04569203 0.95430803]
1
```

## 10) H4K20me1 histone-mark-region classification {#37-h4k20me1-histone-mark-region-classification}

- **Purpose:** Classify whether a sequence corresponds to an H4K20me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="h4k20me1"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 11) Human open-chromatin-region classification {#54-human-open-chromatin-region-classification}

- **Purpose:** Classify human Ensembl open-chromatin region versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_ocr_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```
