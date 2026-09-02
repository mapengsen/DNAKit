# 深度学习任务 Checkpoint 下载

本页集中说明 27 个序列分类任务的 checkpoint 下载、文件名和本地放置方法。只需下载准备预测的任务文件，不需要下载全部权重。预测功能与参数说明见[深度学习性质预测](23_deep_learning_property_prediction.md)。

## 1) 下载入口

[打开 Google Drive checkpoint 文件夹](https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing)

这些权重不会包含在 DNAKit 的 wheel 或 sdist 中，DNAKit 也不会自动下载整个 Google Drive 文件夹。下载后请保留原始文件名；Linux 文件名区分大小写。

## 2) Checkpoint 文件清单

### NT Revised（18 个）

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

### Genomic Benchmarks（9 个）

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

## 3) 本地放置方式

### 默认目录

默认从当前工作目录下的 `./ckpt/enformer-benchmarks/` 查找。例如：

```text
ckpt/
└── enformer-benchmarks/
    └── H2AFZ.ckpt
```

预测时不需要传入 checkpoint 路径：

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    device="cuda",
)
```

### 自定义目录

文件可以直接放在自定义目录，也可以放在该目录下的 `enformer-benchmarks/`、`nt-revised/` 或 `genomic-benchmarks/` 子目录：

```python
result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_dir="/data/checkpoints",
    device="cuda",
)
```

### 精确文件路径

也可以直接指定文件：

```python
result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_path="/data/models/H2AFZ.ckpt",
    device="cuda",
)
```

`checkpoint_dir` 与 `checkpoint_path` 不能同时使用。

## 4) 安装推理依赖

```bash
python -m pip install "dnakit[neural,neural-enformer]"
```

## 5) 加载校验与输出

DNAKit 使用 `torch.load(weights_only=True, mmap=True)` 读取权重，并核对 checkpoint 格式、任务名、数据集家族、类别数、标签映射、分类头形状及全部参数键。文件不完整、文件名对应任务错误或参数结构不兼容时会停止预测。

输出概率按 checkpoint 内部标签顺序排列。请读取 `output.output_names` 和 `output.metadata["predicted_label"]`，不要自行假设类别索引。
