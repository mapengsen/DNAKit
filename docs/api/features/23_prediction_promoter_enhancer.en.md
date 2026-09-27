# Promoters and Enhancers

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) Enhancer classification {#38-enhancer-classification}

- **Purpose:** Classify enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="enhancers"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 2) Enhancer-type classification {#39-enhancer-type-classification}

- **Purpose:** Classify the enhancer category represented by the input sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="enhancers_types"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 3) Promoter classification {#40-promoter-classification}

- **Purpose:** Classify promoter versus non-promoter sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_all"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 4) Non-TATA promoter classification {#41-non-tata-promoter-classification}

- **Purpose:** Classify non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_no_tata"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 5) TATA promoter classification {#42-tata-promoter-classification}

- **Purpose:** Classify TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="promoter_tata"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 6) Drosophila enhancer classification {#48-drosophila-enhancer-classification}

- **Purpose:** Classify Drosophila enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="drosophila_enhancers_stark"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 7) Mouse enhancer classification {#49-mouse-enhancer-classification}

- **Purpose:** Classify mouse Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="dummy_mouse_enhancers_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 8) Human Cohn enhancer classification {#50-human-cohn-enhancer-classification}

- **Purpose:** Classify human Cohn enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_enhancers_cohn"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 9) Human Ensembl enhancer classification {#51-human-ensembl-enhancer-classification}

- **Purpose:** Classify human Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_enhancers_ensembl"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 10) Human regulatory-element-type classification {#52-human-regulatory-element-type-classification}

- **Purpose:** Classify human enhancer, promoter, or open-chromatin-region sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_ensembl_regulatory"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 11) Human non-TATA promoter classification {#53-human-non-tata-promoter-classification}

- **Purpose:** Classify human non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="human_nontata_promoters"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (previous local checkpoint CUDA smoke run):** Class names, probabilities and the predicted label; values depend on the sequence and checkpoint.

```text
('negative', 'positive')
[3.8891116e-08 1.0000000e+00]
positive
```
