# Deep-learning property prediction

This page lists only tasks for which the authors released a **trained prediction head or a native pretrained output**. DNAKit runs those weights without fitting or fine-tuning. The existing [neural representations](08_fingerprints.md) API remains embedding-only.

## Integrated direct tasks

| `model` | Direct `task` values | Output | Input |
| --- | --- | --- | --- |
| `alphagenome` | `atac`, `cage`, `dnase`, `rna_seq`, `chip_histone`, `chip_tf`, `splice_sites`, `splice_site_usage`, `splice_junctions`, `contact_maps`, `procap` | tissue/cell-type tracks; all except `splice_junctions` also support REF, ALT and ALT−REF arrays | DNA sequence or one-SNV context |
| `enformer` | `human_tracks`, `mouse_tracks` | 5,313 human or 1,643 mouse CAGE, DNase/ATAC and TF/histone ChIP tracks in 128-bp bins | DNA sequence or one-SNV context |
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

## Dependencies and boundaries

| Model | Requirement |
| --- | --- |
| SegmentNT, GENERator | `dnakit[neural]` and explicit `allow_remote_code=True` after reviewing official code |
| Enformer | `dnakit[neural,neural-enformer]` |
| Evo 2 | Python 3.11–3.12, `dnakit[neural,neural-evo2]`, compatible GPU |
| AlphaGenome | Python ≥3.11, official research package/JAX, accepted model terms, generally H100-class GPU |
| LucaOneTasks | its separate official environment and source checkout |

`checkpoint_dir` changes the cache root; `checkpoint_path` selects prepared local weights. Prediction loops show progress unless `show_progress=False`.

The NT-v2 base model, DNABERT-2, HyenaDNA, Caduceus-Ph, JanusDNA, CrossDNA, and the benchmark classifiers built on them remain embedding/fine-tuning workflows. ConvNova does not provide a reusable official pretrained task checkpoint. SegmentNT, the Evo 2 exon classifier, and LucaOneTasks are included because their trained downstream heads are separately released.
