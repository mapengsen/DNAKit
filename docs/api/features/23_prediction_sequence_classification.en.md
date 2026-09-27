# Taxonomy and Molecular Function Classification

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) SupKTax classification {#19-supktax-classification}

- **Purpose:** Classify gene sequences using the official SupKTax checkpoint.
- **Prediction method:** Run the official LucaOneTasks `SupKTax` multiclass checkpoint directly, returning the `top_k` highest-probability labels from `label.txt`.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="supktax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. Interpret labels using the checkpoint’s `label.txt`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(5,)
(<five highest-probability official labels>,)
```

## 2) GenusTax classification {#20-genustax-classification}

- **Purpose:** Classify gene sequences using the official GenusTax checkpoint.
- **Prediction method:** Run the official LucaOneTasks `GenusTax` multiclass checkpoint directly, returning the Top-k labels by probability.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="genustax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. The LucaOneTasks README swaps the natural-language taxonomic levels for SupKTax/GenusTax, so use the actual `label.txt` to interpret results.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(5,)
<highest-probability official label>
```

## 3) Species classification {#21-species-classification}

- **Purpose:** Predict the species class of a gene sequence.
- **Prediction method:** Run the official LucaOneTasks `SpeciesTax` multiclass checkpoint directly, returning Top-k species labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="speciestax", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`; use `top_k` to set the number of labels returned.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(5,)
(<five highest-probability species labels>,)
```

## 4) Prokaryotic protein subcellular localization {#22-prokaryotic-protein-subcellular-localization}

- **Purpose:** Predict subcellular localization classes for prokaryotic proteins.
- **Prediction method:** Run the official LucaOneTasks `ProtLoc` multiclass checkpoint directly, returning Top-k localization labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_location", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="protein"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(5,)
(<five highest-probability subcellular localization labels>,)
```

## 5) ncRNA families {#24-ncrna-families}

- **Purpose:** Predict the family of a noncoding RNA sequence.
- **Prediction method:** Run the official LucaOneTasks `ncRNAFam` multiclass checkpoint directly, returning Top-k family labels and probabilities.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="ncrna_family", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="gene"`. DNA or RNA alphabets are accepted; the official task pipeline handles `U`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(5,)
(<five highest-probability ncRNA family labels>,)
```

## 6) Human versus worm sequence classification {#47-human-versus-worm-sequence-classification}

- **Purpose:** Classify a sequence as human or Caenorhabditis elegans.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="demo_human_or_worm"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```
