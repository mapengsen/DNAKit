# Splicing and Genomic Regions

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) Splice donor/acceptor sites {#8-splice-donoracceptor-sites}

- **Purpose:** Predict splice donor and acceptor tracks at each sequence position.
- **Prediction method:** Run the official AlphaGenome `SPLICE_SITES` output head directly, retaining output resolution and track names.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_sites")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`; variant results additionally include an ALT−REF difference axis.
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
        task="splice_sites",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **Example output (schematic structure):**

```text
(<number of positions>, <number of donor/acceptor tracks>)
('position', 'track')
```

## 2) Splice-site usage {#9-splice-site-usage}

- **Purpose:** Predict usage signals or usage-proportion tracks for candidate splice sites.
- **Prediction method:** Run the official AlphaGenome `SPLICE_SITE_USAGE` output head directly. `output_names` identifies the tracks.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_site_usage")[required], backend[optional])`
- **Input:** `BiologicalSequence` or `VariantContext`. Specify the species, ontology terms, and official context length.
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
        task="splice_site_usage",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **Example output (schematic structure):**

```text
(<number of positions>, <number of splice-site usage tracks>)
<bp per output position>
```

## 3) Splice junction counts {#10-splice-junction-counts}

- **Purpose:** Predict connection-count tracks for candidate splice junctions.
- **Prediction method:** Run the official AlphaGenome `SPLICE_JUNCTIONS` output head directly. `junctions` in result metadata contains candidate junction identifiers.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="alphagenome", task="splice_junctions")[required], backend[optional])`
- **Input:** Only `BiologicalSequence` is accepted. Candidate junction sets may change between alleles, so this task currently rejects `VariantContext`.
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
        task="splice_junctions",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
print(len(output.metadata.get("junctions", ())))
```

- **Example output (schematic structure):**

```text
(<number of candidate splice junctions>, <number of tracks>)
('junction', 'track')
<number of candidate splice junctions>
```

## 4) 14-class single-base genomic segmentation {#14-14-class-single-base-genomic-segmentation}

- **Purpose:** Predict probabilities for 14 gene and regulatory-element classes at each base of human DNA.
- **Prediction method:** Use `InstaDeepAI/segment_nt`, an NT-v2 backbone with a trained U-Net segmentation head. Apply softmax to binary logits and return the probability that each feature is present.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation", allow_remote_code=True)[required], backend[optional])`
- **Input:** Human `BiologicalSequence` up to 30,000 bp. The table below lists the 14 output classes.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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
output = result.records[0].output
print(output.values.shape)
print(len(output.output_names))
```

- **Example output:**

```text
(4000, 14)
14
```

The following 14 classes share this API. Each corresponds to a column in `output.values`; use `output.output_names` for the actual column names and order.

| Class | Meaning |
| --- | --- |
| <span id="141-protein-coding-gene"></span>Protein-coding gene | Predict the probability that each base belongs to a protein-coding gene region. |
| <span id="142-long-noncoding-rna-lncrna"></span>Long noncoding RNA (lncRNA) | Predict the probability that each base belongs to a long noncoding rna (lncrna) region. |
| <span id="143-exon"></span>Exon | Predict the probability that each base belongs to a exon region. |
| <span id="144-intron"></span>Intron | Predict the probability that each base belongs to a intron region. |
| <span id="145-splice-donor"></span>Splice donor | Predict the probability that each base belongs to a splice donor region. |
| <span id="146-splice-acceptor"></span>Splice acceptor | Predict the probability that each base belongs to a splice acceptor region. |
| <span id="147-5-untranslated-region-5-utr"></span>5′ untranslated region (5′ UTR) | Predict the probability that each base belongs to a 5′ untranslated region (5′ utr) region. |
| <span id="148-3-untranslated-region-3-utr"></span>3′ untranslated region (3′ UTR) | Predict the probability that each base belongs to a 3′ untranslated region (3′ utr) region. |
| <span id="149-ctcf-bound-region"></span>CTCF-bound region | Predict the probability that each base belongs to a ctcf-bound region region. |
| <span id="1410-polyadenylation-signal-polya-signal"></span>Polyadenylation signal (polyA signal) | Predict the probability that each base belongs to a polyadenylation signal (polya signal) region. |
| <span id="1411-tissue-specific-enhancer"></span>Tissue-specific enhancer | Predict the probability that each base belongs to a tissue-specific enhancer region. |
| <span id="1412-tissue-invariant-enhancer"></span>Tissue-invariant enhancer | Predict the probability that each base belongs to a tissue-invariant enhancer region. |
| <span id="1413-tissue-specific-promoter"></span>Tissue-specific promoter | Predict the probability that each base belongs to a tissue-specific promoter region. |
| <span id="1414-tissue-invariant-promoter"></span>Tissue-invariant promoter | Predict the probability that each base belongs to a tissue-invariant promoter region. |

## 5) Exon probability {#16-exon-probability}

- **Purpose:** Predict whether a candidate site belongs to an exon using forward and reverse contexts around it.
- **Prediction method:** Read 4,096-dimensional last-token embeddings for both contexts from Evo 2 `blocks.26`, concatenate them into 8,192 dimensions, and apply the officially trained MLP head. The classifier is not retrained.
- **API:** `dnakit.predictions.predict_pair_properties(inputs[required], config=PropertyPredictionConfig(model="evo2", task="exon_probability", allow_remote_code=True)[required], backend[optional])`
- **Input:** `BiologicalSequencePair`, ordered as forward context then reverse context, both with `sequence_type="gene"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

contexts = BiologicalSequencePair(
    "site-1",
    BiologicalSequence("forward-1", "ACGT" * 256),
    BiologicalSequence("reverse-1", "ACGT" * 256),
)
result = predict_pair_properties(
    [contexts],
    config=PropertyPredictionConfig(
        model="evo2",
        task="exon_probability",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **Example output (label depends on probability and threshold):**

```text
('exon_probability',)
exon or non_exon
```

## 6) Splice-acceptor-site classification {#43-splice-acceptor-site-classification}

- **Purpose:** Classify splice acceptor versus non-acceptor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_acceptors"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_acceptors",
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

## 7) Splice-site-type classification {#44-splice-site-type-classification}

- **Purpose:** Classify no-splice, splice-acceptor, or splice-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_all"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_all",
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

## 8) Splice-donor-site classification {#45-splice-donor-site-classification}

- **Purpose:** Classify splice donor versus non-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="splice_sites_donors"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_donors",
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

## 9) Coding-region versus intergenic-region classification {#46-coding-region-versus-intergenic-region-classification}

- **Purpose:** Classify coding versus intergenic genomic sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs[required], task="demo_coding_vs_intergenomic_seqs"[required], checkpoint_dir or checkpoint_path[optional])`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

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

- **Example output (values vary with sequence and checkpoint):**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```
