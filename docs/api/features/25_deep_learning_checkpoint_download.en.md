# Deep-learning task checkpoint download

This page centralizes download, filename, and local-placement instructions for the 27 sequence-classification task checkpoints. Download only the files for tasks you intend to run. See [Deep-learning property prediction](23_deep_learning_property_prediction.md) for task behavior and parameters.

## 1) Download entry

[Open the Google Drive checkpoint folder](https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing)

The weights are not included in the DNAKit wheel or sdist, and DNAKit does not automatically download the whole Google Drive folder. Preserve the original filenames; Linux filenames are case-sensitive.

## 2) Checkpoint file list

### NT Revised (18)

```text
H2AFZ.ckpt
H3K27ac.ckpt
H3K27me3.ckpt
H3K36me3.ckpt
H3K4me1.ckpt
H3K4me2.ckpt
H3K4me3.ckpt
H3K9ac.ckpt
H3K9me3.ckpt
H4K20me1.ckpt
enhancers.ckpt
enhancers_types.ckpt
promoter_all.ckpt
promoter_no_tata.ckpt
promoter_tata.ckpt
splice_sites_acceptors.ckpt
splice_sites_all.ckpt
splice_sites_donors.ckpt
```

### Genomic Benchmarks (9)

```text
demo_coding_vs_intergenomic_seqs.ckpt
demo_human_or_worm.ckpt
drosophila_enhancers_stark.ckpt
dummy_mouse_enhancers_ensembl.ckpt
human_enhancers_cohn.ckpt
human_enhancers_ensembl.ckpt
human_ensembl_regulatory.ckpt
human_nontata_promoters.ckpt
human_ocr_ensembl.ckpt
```

## 3) Local placement

### Default directory

The default lookup location is `./ckpt/enformer-benchmarks/` under the current working directory. For example:

```text
ckpt/
└── enformer-benchmarks/
    └── H2AFZ.ckpt
```

No checkpoint path is required in the prediction call:

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    device="cuda",
)
```

### Custom directory

A file may be placed directly in the custom directory or in its `enformer-benchmarks/`, `nt-revised/`, or `genomic-benchmarks/` subdirectory:

```python
result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_dir="/data/checkpoints",
    device="cuda",
)
```

### Exact file path

An exact checkpoint file can also be selected:

```python
result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_path="/data/models/H2AFZ.ckpt",
    device="cuda",
)
```

`checkpoint_dir` and `checkpoint_path` are mutually exclusive.

## 4) Install inference dependencies

```bash
python -m pip install "dnakit[neural,neural-enformer]"
```

## 5) Loading validation and output

DNAKit reads weights with `torch.load(weights_only=True, mmap=True)` and validates checkpoint format, task name, dataset family, class count, label mapping, classification-head shape, and every parameter key. Incomplete files, task mismatches, or incompatible parameter structures stop prediction.

Probabilities follow the label order stored inside the checkpoint. Read `output.output_names` and `output.metadata["predicted_label"]`; do not assume class-index meanings.
