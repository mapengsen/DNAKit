# Molecular Relationships and Interactions

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) Nucleic acid–protein central dogma relationship {#18-nucleic-acidprotein-central-dogma-relationship}

- **Purpose:** Predict whether a gene sequence and a protein sequence have a CentralDogma relationship.
- **Prediction method:** Run the official LucaOneTasks `CentralDogma` binary-classification checkpoint through `src/predict_v1.py`, returning the positive-class probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="central_dogma", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, with `sequence_type="gene"` first and `sequence_type="protein"` second. `model_source_path` must point to the official LucaOneTasks source directory.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 2) Influenza A antigenic relationships {#25-influenza-a-antigenic-relationships}

- **Purpose:** Predict the antigenic relationship class between two influenza A gene sequences.
- **Prediction method:** Run the official LucaOneTasks `InfA` binary-classification checkpoint directly, returning the positive-class probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="influenza_antigenicity", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`; both sequences must have `sequence_type="gene"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 3) Protein–protein interactions {#26-proteinprotein-interactions}

- **Purpose:** Predict whether two protein sequences have a protein–protein interaction.
- **Prediction method:** Run the official LucaOneTasks `PPI` binary-classification checkpoint directly, returning the interaction probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_interaction", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`; both sequences must have `sequence_type="protein"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```

## 4) ncRNA–protein interactions {#27-ncrnaprotein-interactions}

- **Purpose:** Predict whether an ncRNA and a protein sequence interact.
- **Prediction method:** Run the official LucaOneTasks `ncRPI` binary-classification checkpoint directly, returning the interaction probability and official label.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="ncrna_protein_interaction", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, with `sequence_type="gene"` first and `sequence_type="protein"` second.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (label determined by the official checkpoint):**

```text
('probability',)
<official predicted label>
```
