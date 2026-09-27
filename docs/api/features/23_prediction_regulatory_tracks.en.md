# Expression, Regulatory Tracks and Contact Maps

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) RNA-seq coverage and expression {#1-rna-seq-coverage-and-expression}

- **Purpose:** Predict RNA-seq coverage and gene-expression tracks from DNA sequences.
- **Prediction method:** Run direct inference with the official AlphaGenome all-folds weights and the `RNA_SEQ` output. Short sequences are centered and padded with `N` to an officially supported context length; use `ontology_terms` to restrict tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="rna_seq")[required], backend[optional])`
- **Input:** `BiologicalSequence`; alternatively use `VariantContext` with `predict_variant_effects()` to obtain REF, ALT, and ALT−REF tracks. Human and mouse are supported.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of RNA-seq tracks>)
('position', 'track')
```

## 2) CAGE transcription initiation and expression {#2-cage-transcription-initiation-and-expression}

- **Purpose:** Predict CAGE signals describing capped-transcript initiation sites and expression levels.
- **Prediction method:** The official AlphaGenome `CAGE` head directly produces a position × track matrix; `ontology_terms` filters target tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="cage")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a single-SNV `VariantContext`. Sequences may be up to 1,048,576 bp; human and mouse are supported.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of CAGE tracks>)
('position', 'track')
```

## 3) PRO-cap nascent transcription initiation {#3-pro-cap-nascent-transcription-initiation}

- **Purpose:** Predict PRO-cap nascent RNA 5′-end signals to locate active transcription initiation.
- **Prediction method:** Run the official AlphaGenome `PROCAP` output head directly. Output resolution is recorded in `metadata["resolution_bp"]`.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="procap")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Set `organism` and `ontology_terms` as needed.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of PRO-cap tracks>)
<bp per output position>
```

## 4) ATAC-seq chromatin accessibility {#4-atac-seq-chromatin-accessibility}

- **Purpose:** Predict ATAC-seq open-chromatin tracks to identify potentially accessible regulatory regions.
- **Prediction method:** Run the official AlphaGenome `ATAC` output head directly, returning accessibility signals by position for different tissues or cell types.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="atac")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a `VariantContext` containing exactly one base change. Use `ontology_terms` to restrict tracks.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of filtered ATAC-seq tracks>)
(<first two official track names>,)
```

## 5) DNase-seq chromatin accessibility {#5-dnase-seq-chromatin-accessibility}

- **Purpose:** Predict DNase-seq hypersensitivity signals describing chromatin accessibility.
- **Prediction method:** Run the official AlphaGenome `DNASE` output head directly. Short inputs are centered and padded; padding lengths are recorded in result metadata.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="dnase")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Select `organism="human"` or `"mouse"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of DNase-seq tracks>)
16384
```

## 6) Histone modifications {#6-histone-modifications}

- **Purpose:** Predict histone ChIP-seq tracks describing histone modification signals across a sequence region.
- **Prediction method:** Run the official AlphaGenome `CHIP_HISTONE` output head directly, returning a position × histone-track matrix.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="chip_histone")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Select outputs by species and ontology term.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of histone ChIP-seq tracks>)
('position', 'track')
```

## 7) Transcription factor binding {#7-transcription-factor-binding}

- **Purpose:** Predict transcription factor ChIP-seq binding tracks to locate potential binding signals.
- **Prediction method:** Run the official AlphaGenome `CHIP_TF` output head directly, returning a position × transcription-factor-track matrix.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="chip_tf")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Use `ontology_terms` to restrict tissues or cell types.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of positions or bins>, <number of TF ChIP-seq tracks>)
(<first two official track names>,)
```

## 8) 3D DNA contact maps {#11-3d-dna-contact-maps}

- **Purpose:** Predict 3D DNA contact probability maps between chromatin regions.
- **Prediction method:** Run the official AlphaGenome `CONTACT_MAPS` output head directly, producing two position axes and one track axis.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="contact_maps")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`; variant results add a leading REF, ALT, ALT−REF axis.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (schematic structure):**

```text
(<number of bins>, <number of bins>, <number of contact-map tracks>)
('position_1', 'position_2', 'track')
```

## 9) Human regulatory tracks {#12-human-regulatory-tracks}

- **Purpose:** Predict 5,313 human CAGE, DNase/ATAC, transcription factor ChIP, and histone ChIP tracks.
- **Prediction method:** Run public pretrained Enformer weights directly. By default, center and pad the input with `N` to 196,608 bp and return 896 bins of 128 bp.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="enformer", task="human_tracks")[required], backend[optional])`
- **Input:** `BiologicalSequence`; alternatively use `VariantContext` with `predict_variant_effects()` for REF, ALT, and ALT−REF. `max_length` must be a multiple of 128 in the range 114,688–393,216.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (default context):**

```text
(896, 5313)
128
```

## 10) Mouse regulatory tracks {#13-mouse-regulatory-tracks}

- **Purpose:** Predict 1,643 mouse CAGE, DNase/ATAC, transcription factor ChIP, and histone ChIP tracks.
- **Prediction method:** Run Enformer’s mouse output head directly with a default context of 196,608 bp and tracks binned at 128 bp.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="enformer", task="mouse_tracks")[required], backend[optional])`
- **Input:** `BiologicalSequence` or a single-SNV `VariantContext`. Short inputs are automatically centered with `N` padding at both ends.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
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

- **Example output (default context):**

```text
(896, 1643)
mouse
```

<!-- BEGIN GENERATED API: model-tracks -->

<span id="model-tracks"></span>

## 11) Predict multimodal genomic tracks with AlphaGenome {#dnakit.databases.alpha_genome_predict_interval}

- **Purpose:** Predict multimodal genomic tracks with AlphaGenome.
- **Method:** Use `AlphaGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alpha_genome_predict_interval(*, chromosome, start, end, output_types=None, ontology_terms=None, organism=None, sequence_length=None, verbose=True)`
- **Input:** Provide `chromosome`, `start`, `end`. Configure: `ALPHA_GENOME_API_KEY`.
- **Output:** Returns a mapping with `interval`, `tracks` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alpha_genome_predict_interval(
    chromosome='chr19',
    start=40991281,
    end=41018398,
    output_types=['RNA_SEQ', 'ATAC'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "interval": "<string>",
  "tracks": [
    "<object>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_regulatory_tracks.md#dnakit.databases.alpha_genome_predict_interval-parameters)

<!-- END GENERATED API: model-tracks -->
