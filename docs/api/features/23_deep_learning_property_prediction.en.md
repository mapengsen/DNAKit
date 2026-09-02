# Deep-learning property prediction

This page lists only tasks for which the authors released a **trained prediction head or a native pretrained output**. DNAKit runs those weights without fitting or fine-tuning. The existing [neural representations](08_fingerprints.md) API remains embedding-only.

## Integrated direct tasks

| `model` | Direct `task` values | Output | Input |
| --- | --- | --- | --- |
| `alphagenome` | `atac`, `cage`, `dnase`, `rna_seq`, `chip_histone`, `chip_tf`, `splice_sites`, `splice_site_usage`, `splice_junctions`, `contact_maps`, `procap` | tissue/cell-type tracks; all except `splice_junctions` also support REF, ALT and ALT−REF arrays | DNA sequence or one-SNV context |
| `enformer` | `human_tracks`, `mouse_tracks`, plus 27 NT Revised/Genomic Benchmarks classifiers listed below | regulatory tracks or checkpoint-ordered class probabilities | DNA sequence or one-SNV context for tracks; DNA sequence for classifiers |
| `segmentnt` | `genomic_segmentation` | per-base probabilities for 14 genomic annotations | human DNA, up to 30 kb |
| `evo2` | `variant_effect`, `exon_probability` | zero-shot likelihood delta or trained exon probability | one-SNV context or forward/reverse contexts |
| `generator` | `variant_effect` | `p_ref`, `p_alt`, and `log(p_ref/(p_alt+1e-10))` | one-SNV context |
| `lucaone` | ten LucaOneTasks heads | classification probabilities/Top-k labels or regression | gene/protein sequences and pairs |

The papers, official repositories, and checkpoint sources used by this chapter are collected in [FAQ: references for deep-learning property prediction](../../faq.md#deep-learning-property-prediction-references).

AlphaGenome predicts RNA-seq, CAGE, PRO-cap, ATAC/DNase accessibility, histone and TF ChIP, splice sites/usage/junctions, and contact maps. `ontology_terms` limits tissue or cell-type tracks. DNAKit center-pads to an official supported length (16,384, 131,072, 524,288, or 1,048,576 bp) and records that padding.

SegmentNT returns these 14 features: `protein_coding_gene`, `lncRNA`, `exon`, `intron`, `splice_donor`, `splice_acceptor`, `5UTR`, `3UTR`, `CTCF-bound`, `polyA_signal`, tissue-specific/invariant enhancers, and tissue-specific/invariant promoters. This uses the trained `InstaDeepAI/segment_nt` checkpoint, not the NT-v2 base checkpoint.

LucaOneTasks exposes `central_dogma`, `supktax`, `genustax`, `speciestax`, `protein_location`, `protein_stability`, `ncrna_family`, `influenza_antigenicity`, `protein_interaction`, and `ncrna_protein_interaction`. DNAKit uses the exact timestamps and steps in the official `prediction.sh`. Because the upstream README swaps the natural-language descriptions of SupKTax and GenusTax, the checkpoint's actual `label.txt` remains authoritative.

## Unified API

Registry queries do not download weights:

```python
from dnakit.predictions import available_prediction_models, available_prediction_tasks

print(available_prediction_models())
print(available_prediction_tasks("evo2"))
```

Use `predict_properties()` or the typed convenience wrappers `predict_sequence_properties()`, `predict_pair_properties()`, and `predict_variant_effects()`. A result contains one read-only float32 array per input. `output_names` names the final axis, while metadata records axes, padding, resolution, thresholds, and score definitions. Use `to_dict(include_values=False)` to avoid expanding large tracks into JSON.

### Single-nucleotide genome segmentation

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
print(result.records[0].output.values.shape)  # (4000, 14)
```

### Variant scoring

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="A" * 4096 + "C" + "G" * 4095,
    alternate_sequence="A" * 4096 + "T" + "G" * 4095,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(model="evo2", task="variant_effect"),
)
print(result.records[0].output.values)
```

Evo 2's likelihood delta and GENERator's log-ratio are continuous zero-shot ranking scores, not calibrated clinical pathogenicity probabilities.

### Released downstream task heads

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ACGT" * 100)],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="speciestax",
        model_source_path="/opt/LucaOneTasks",
        device="cuda:0",
    ),
)
print(result.records[0].output.metadata["predicted_label"])
```

This adapter executes only the official `src/predict_v1.py` inside the explicitly supplied checkout. Its first run may download LucaOne and task-head checkpoints.

## 28) H2A.Z histone-variant-region classification

- **Purpose:** Classify whether a sequence corresponds to an H2A.Z (H2AFZ) ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h2afz", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H2AFZ.ckpt`; the default path is `./ckpt/enformer-benchmarks/H2AFZ.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 29) H3K27ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k27ac", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K27ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27ac.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 30) H3K27me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K27me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k27me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K27me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K27me3.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 31) H3K36me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K36me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k36me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K36me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K36me3.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 32) H3K4me1 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me1", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me1.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 33) H3K4me2 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me2 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me2", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me2.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me2.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 34) H3K4me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K4me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k4me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K4me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K4me3.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 35) H3K9ac histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9ac ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k9ac", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K9ac.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9ac.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 36) H3K9me3 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H3K9me3 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h3k9me3", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H3K9me3.ckpt`; the default path is `./ckpt/enformer-benchmarks/H3K9me3.ckpt`.
- **Example:**

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

- **CUDA smoke result with the local sample checkpoint:**

```text
('0', '1')
[0.04569203 0.95430803]
1
```

## 37) H4K20me1 histone-mark-region classification

- **Purpose:** Classify whether a sequence corresponds to an H4K20me1 ChIP-seq peak.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="h4k20me1", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `H4K20me1.ckpt`; the default path is `./ckpt/enformer-benchmarks/H4K20me1.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 38) Enhancer classification

- **Purpose:** Classify enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="enhancers", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `enhancers.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="enhancers",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 39) Enhancer-type classification

- **Purpose:** Classify the enhancer category represented by the input sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="enhancers_types", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `enhancers_types.ckpt`; the default path is `./ckpt/enformer-benchmarks/enhancers_types.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="enhancers_types",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 40) Promoter classification

- **Purpose:** Classify promoter versus non-promoter sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_all", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_all.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 41) Non-TATA promoter classification

- **Purpose:** Classify non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_no_tata", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_no_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 42) TATA promoter classification

- **Purpose:** Classify TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="promoter_tata", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `promoter_tata.ckpt`; the default path is `./ckpt/enformer-benchmarks/promoter_tata.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 43) Splice-acceptor-site classification

- **Purpose:** Classify splice acceptor versus non-acceptor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_acceptors", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_acceptors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_acceptors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 44) Splice-site-type classification

- **Purpose:** Classify no-splice, splice-acceptor, or splice-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_all", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_all.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 45) Splice-donor-site classification

- **Purpose:** Classify splice donor versus non-donor sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="splice_sites_donors", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `splice_sites_donors.ckpt`; the default path is `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="splice_sites_donors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 46) Coding-region versus intergenic-region classification

- **Purpose:** Classify coding versus intergenic genomic sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="demo_coding_vs_intergenomic_seqs", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `demo_coding_vs_intergenomic_seqs.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 47) Human versus worm sequence classification

- **Purpose:** Classify a sequence as human or Caenorhabditis elegans.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="demo_human_or_worm", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `demo_human_or_worm.ckpt`; the default path is `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 48) Drosophila enhancer classification

- **Purpose:** Classify Drosophila enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="drosophila_enhancers_stark", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `drosophila_enhancers_stark.ckpt`; the default path is `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="drosophila_enhancers_stark",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 49) Mouse enhancer classification

- **Purpose:** Classify mouse Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="dummy_mouse_enhancers_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `dummy_mouse_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 25)],
    task="dummy_mouse_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 50) Human Cohn enhancer classification

- **Purpose:** Classify human Cohn enhancer versus non-enhancer sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_enhancers_cohn", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_enhancers_cohn.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`.
- **Example:**

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

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 51) Human Ensembl enhancer classification

- **Purpose:** Classify human Ensembl enhancer versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_enhancers_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_enhancers_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="human_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```

## 52) Human regulatory-element-type classification

- **Purpose:** Classify human enhancer, promoter, or open-chromatin-region sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 3 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_ensembl_regulatory", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_ensembl_regulatory.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="human_ensembl_regulatory",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>, <class-2 probability>)
<label with the highest checkpoint probability>
```

## 53) Human non-TATA promoter classification

- **Purpose:** Classify human non-TATA promoter versus negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_nontata_promoters", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_nontata_promoters.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`.
- **Example:**

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

- **CUDA smoke result with the local sample checkpoint:**

```text
('negative', 'positive')
[3.8891116e-08 1.0000000e+00]
positive
```

## 54) Human open-chromatin-region classification

- **Purpose:** Classify human Ensembl open-chromatin region versus matched negative sequence.
- **Method:** Load the fully fine-tuned backbone and classification head, valid-bin mean-pool the final embeddings, and return 2 softmax probabilities. No training or fine-tuning occurs at runtime.
- **API:** `dnakit.predictions.predict_enformer_benchmark(inputs, task="human_ocr_ensembl", checkpoint_dir=...)`
- **Input:** `BiologicalSequence` or `DNARecord`; A/C/G/T and IUPAC ambiguity symbols are accepted, and biological bases are never silently truncated.
- **Checkpoint:** Follow [the task checkpoint download instructions](25_deep_learning_checkpoint_download.md) to download `human_ocr_ensembl.ckpt`; the default path is `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`.
- **Example:**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="human_ocr_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **Example result structure:**

```text
(<class-0 probability>, <class-1 probability>)
<label with the highest checkpoint probability>
```


## Dependencies and boundaries

| Model | Requirement |
| --- | --- |
| SegmentNT, GENERator | `dnakit[neural]` and explicit `allow_remote_code=True` after reviewing official code |
| Enformer | `dnakit[neural,neural-enformer]` |
| Evo 2 | Python 3.11–3.12, `dnakit[neural,neural-evo2]`, compatible GPU |
| AlphaGenome | Python ≥3.11, official research package/JAX, accepted model terms, generally H100-class GPU |
| LucaOneTasks | its separate official environment and source checkout |

`checkpoint_dir` changes the cache root; `checkpoint_path` selects prepared local weights. Prediction loops show progress unless `show_progress=False`.

The NT-v2 base model, DNABERT-2, HyenaDNA, Caduceus-Ph, JanusDNA, and CrossDNA still require a separately trained task head. This addition runs only the 27 supplied full-model task checkpoints; it does not turn those embedding models into zero-shot classifiers. ConvNova does not provide a reusable official pretrained task checkpoint. SegmentNT, the Evo 2 exon classifier, and LucaOneTasks are included because their trained downstream heads are separately released.
