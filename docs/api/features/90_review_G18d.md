# G18d 复核：kipoi、genomic_benchmarks 与 gReLU

本页记录 G18d 组的复核结果。台账中状态为 `other` 或 `cap` 的行共 121 个：kipoi 86 个、genomic_benchmarks 3 个、grelu 32 个（31 个 other 与 1 个 cap）。复核后每一行都有新状态和一行原因。本组为 gReLU 新增 20 个包装函数，位于 `src/dnakit/predictions/grelu_review_functions.py`；kipoi 与 genomic_benchmarks 没有可离线运行的 DNA 级函数，对应模块 `kipoi_review_functions.py` 与 `genomic_benchmarks_review_functions.py` 不导出函数。

- kipoi 的 86 行全部是插件注册、YAML 与文件辅助、工作目录上下文、命令行参数与写出器等非 DNA 代码。`kipoi.list_sources` 与 `kipoi.config.list_sources` 会在首次调用时克隆 202 MB 的模型仓库，标记为 `needs-network`，测试中以模拟对象替代，从不实际运行。已有的 `list_kipoi_models` 与 `describe_kipoi_model` 不重复包装。
- genomic_benchmarks 的 3 行：标签顺序辅助函数（不含序列内容，且文件系统顺序不确定）标记为 `not-dna`；两个会删除生成数据目录的缓存函数标记为 `internal`。
- gReLU 的 32 行中，20 个可离线运行的 DNA 级函数已包装；`filter_cells` 为 `not-dna`；`update_ckpt`、`check_bigwig`、`check_fasta`、`filter_intervals`、`filter_obs` 与两组资源路径函数为 `internal`；两个 `duplicate-of` 指向已包装的同一上游对象。

## 复核结果

| 包 | 复核前 other/cap | 新状态 | 新增 wrapped |
| --- | --- | --- | --- |
| kipoi | 86 | internal 58，not-dna 24，needs-network 2，cli 2 | 0 |
| genomic_benchmarks | 3 | internal 2，not-dna 1 | 0 |
| grelu | 32（31 other + 1 cap） | wrapped 20（其中 cap 1），duplicate-of 2，internal 9，not-dna 1 | 20 |

## 上游问题与包装处理

- `grelu.io.bed.read_narrowpeak`：上游把第一行当作表头，无表头文件的第一个峰会被丢弃。包装函数默认 `has_header=False`。
- `grelu.variant.predict_variant_effects`：`check_reference` 旗标与同名函数冲突，上游检查不可用。包装函数用本地 FASTA 自行核对参考碱基，`pos` 按 1 起算。
- `grelu.interpret.simulate.marginalize_patterns`：原始输出模式与 `n_shuffles=0` 在当前 torch 上都报 concat 错误。包装函数要求 `compare_func` 为 subtract 或 divide，且 `n_shuffles` 至少为 1。
- `grelu.design.ledidi`：上游传入 0 维目标，ledidi 2.2.0 报 IndexError。包装函数直接调用 `ledidi.Ledidi`，目标长度为 1，损失相同。
- `grelu.data.preprocess.filter_coverage`：percentile 模式缺少分位数参数，无法运行；负样本随机抽取依赖全局种子。两者均不对外提供。
- `grelu.design.evolve`：模型必须在 `data_params['tasks']['name']` 中有任务名，否则包装函数拒绝执行。
- `grelu.interpret.modisco.run_modisco`：随机序列通常找不到模式，此时不写出任何文件，结果中 `patterns_found` 为 False。

## 未包装的其余行

- kipoi：`list_sources` 两行为 `needs-network`；其余为 `not-dna`、`internal` 与 `cli`（见上）。
- genomic_benchmarks：`labels_in_order` 为 `not-dna`；两个 `remove_dataset_from_disk` 为 `internal`。
- grelu：`filter_cells` 为 `not-dna`；`filter_intervals`、`filter_obs`、`update_ckpt`、`check_bigwig`、`check_fasta`，以及 `get_meme_file_path`、`get_blacklist_file` 的两组路径函数为 `internal`；`data.dataset.check_intervals` 与 `data.preprocess.get_chromosomes` 为 `duplicate-of`。

## 包装函数

测试文件为 `others/tests/unit/test_review_G18d.py`。每个函数都由真实上游调用或手工计算检验。下面各节的示例输出都来自实际运行。

## 1) `read_narrowpeak_intervals`

- **作用：** 读取 narrowPeak（BED6+4）峰文件，保留文件中的顺序，每个峰包含 chrom、start、end、summit、signal、p 值与 q 值。
- **处理方法：** 调用 `grelu.io.bed.read_narrowpeak`。上游默认把第一行当作表头，会悄悄丢掉无表头文件的第一个峰，因此包装函数默认 `has_header=False`，并通过 `header=None` 传给 pandas。少于十列的文件会被拒绝。
- **API：** `read_narrowpeak_intervals(path[必须], skiprows[可选], has_header[可选], verbose[可选])`
- **输入：** 本地 narrowPeak 文件路径。
- **输出：** 字典，含 `count`、`peaks`（每个峰的字段）与 `provenance`。

**示例文件：** `peaks.narrowPeak`

```text
chr1	100	200	p1	500	.	5.5	3.2	1.1	50
chr2	300	420	p2	80	+	2.0	-1	-1	60
```

**示例代码：**

```python
from dnakit.predictions import read_narrowpeak_intervals

result = read_narrowpeak_intervals("peaks.narrowPeak")
print(result["count"])
for peak in result["peaks"]:
    print(peak["chrom"], peak["start"], peak["end"], peak["summit"], peak["signal"])
```

**示例输出：**

```text
2
chr1 100 200 50 5.5
chr2 300 420 60 2.0
```

## 2) `expand_chromosome_names`

- **作用：** 把染色体简写（`autosomes`、`autosomesX`、`autosomesXY`）展开为 GRCh38 风格的染色体名，或原样检查显式列表。
- **处理方法：** 调用 `grelu.data.utils.get_chromosomes`。上游对未知简写会原样返回字符串，包装函数改为拒绝未知简写。
- **API：** `expand_chromosome_names(chroms[必须], verbose[可选])`
- **输入：** 简写字符串，或染色体名字符串列表。
- **输出：** 字典，含 `count`、`chromosomes` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import expand_chromosome_names

print(expand_chromosome_names("autosomesXY")["chromosomes"][-3:])
print(expand_chromosome_names(["chrM", "chr1"])["count"])
```

**示例输出：**

```text
['chr22', 'chrX', 'chrY']
2
```

## 3) `check_genomic_interval_table`

- **作用：** 检查一张区间表是否符合 gReLU 的基因组区间格式：前三列为 `chrom`、`start`、`end`，`start` 与 `end` 为整数。
- **处理方法：** 调用 `grelu.sequence.format.check_intervals`。行按给定顺序组成表格，列顺序错误或坐标为浮点数时 `valid` 为 False。
- **API：** `check_genomic_interval_table(intervals[必须], verbose[可选])`
- **输入：** 映射列表，每项含 `chrom`、`start`、`end`。
- **输出：** 字典，含 `valid`、`columns`、`dtypes`、`row_count` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import check_genomic_interval_table

good = [{"chrom": "chr1", "start": 1, "end": 5}]
bad = [{"chrom": "chr1", "start": 1.5, "end": 5}]
print(check_genomic_interval_table(good)["valid"])
print(check_genomic_interval_table(bad)["dtypes"]["start"])
print(check_genomic_interval_table(bad)["valid"])
```

**示例输出：**

```text
True
float64
False
```

## 4) `filter_interval_coverage`

- **作用：** 按多个任务上的覆盖度保留基因组区间：`cutoff` 模式保留不低于阈值的区间，`top` 模式保留覆盖度最高的前 n 个区间。
- **处理方法：** 调用 `grelu.data.preprocess.filter_coverage`。每个区间的值是各任务的聚合值（mean、sum、max 或 min）。上游的 percentile 模式有缺陷，负样本随机抽取需要随机种子，因此两者都不对外提供。
- **API：** `filter_interval_coverage(coverage[必须], intervals[必须], tasks[可选], aggfunc[可选], method[可选], cutoff[可选], verbose[可选])`
- **输入：** coverage 为每个任务一行、每个区间一列的矩阵；intervals 为区间映射列表。
- **输出：** 字典，含 `count`、`intervals`（保留区间的原始下标）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import filter_interval_coverage

intervals = [
    {"chrom": "chr1", "start": 0, "end": 5},
    {"chrom": "chr1", "start": 5, "end": 9},
    {"chrom": "chr2", "start": 1, "end": 4},
]
coverage = [[0.5, 3.0, 1.0]]
cutoff = filter_interval_coverage(coverage, intervals, cutoff=1.0)
print(cutoff["count"], [item["index"] for item in cutoff["intervals"]])
top = filter_interval_coverage(coverage, intervals, method="top", cutoff=1)
print([item["chrom"] + ":" + str(item["start"]) for item in top["intervals"]])
```

**示例输出：**

```text
2 [1, 2]
['chr1:5']
```

## 5) `add_negative_intervals`

- **作用：** 在正区间表之后追加阴性对照区间，阴性列的覆盖度为零，并用 `class` 字段区分正负样本。
- **处理方法：** 调用 `grelu.data.preprocess.add_negatives`，并使用 `negative_labels=0`。上游只有该取值可以正常运行。
- **API：** `add_negative_intervals(coverage[必须], intervals[必须], negatives[必须], tasks[可选], verbose[可选])`
- **输入：** coverage 为任务 × 正区间矩阵；intervals 为正区间；negatives 为阴性区间。
- **输出：** 字典，含 `count`、`positive_count`、`negative_count`、`intervals`、`coverage` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import add_negative_intervals

positives = [{"chrom": "chr1", "start": 0, "end": 5}, {"chrom": "chr1", "start": 5, "end": 9}]
negatives = [{"chrom": "chr3", "start": 0, "end": 7}]
result = add_negative_intervals([[1.0, 2.0]], positives, negatives)
print(result["count"], result["negative_count"])
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["class"])
print(result["coverage"])
```

**示例输出：**

```text
3 1
chr1 0 5 positive
chr1 5 9 positive
chr3 0 7 negative
[[1.0, 2.0, 0.0]]
```

## 6) `read_tomtom_matches`

- **作用：** 读取 TOMTOM 输出目录中所有 `tomtom.tsv` 文件，保留 q 值不超过阈值的基序匹配。
- **处理方法：** 调用 `grelu.io.read_tomtom`，并丢弃含缺失值的行（与上游一致）。结果按查询、靶标、方向和偏移排序；目录中没有 `tomtom.tsv` 时拒绝。
- **API：** `read_tomtom_matches(tomtom_dir[必须], qthresh[可选], verbose[可选])`
- **输入：** 包含 `tomtom.tsv` 文件的本地目录。
- **输出：** 字典，含 `count`、`matches`（每行的列）与 `provenance`。

**示例文件：** `tomtom.tsv`

```text
Target_ID	Optimal_offset	p-value	E-value	q-value	Overlap	Query_consensus	Target_consensus	Orientation
MotifA	0	0.001	0.01	0.02	8	ACGTACGT	ACGTACGA	+
MotifB	1	0.2	0.5	0.4	7	ACGTACGT	TTTTACGT	-
```

**示例代码：**

```python
from dnakit.predictions import read_tomtom_matches

result = read_tomtom_matches(".", qthresh=0.05)
print(result["count"])
print(result["matches"][0]["Target_ID"], result["matches"][0]["q-value"])
```

**示例输出：**

```text
1
MotifA 0.02
```

## 7) `read_modisco_motifs`

- **作用：** 从 TF-MoDISco 的 H5 报告读取基序，返回按两端修剪后的位置概率矩阵（PPM）。
- **处理方法：** 调用 `grelu.io.motifs.read_modisco_report`。`group` 选择 pos 或 neg，`names` 按名称选择并覆盖 `group`；修剪阈值 `trim_threshold` 默认 0.3。
- **API：** `read_modisco_motifs(h5_path[必须], group[可选], names[可选], trim_threshold[可选], verbose[可选])`
- **输入：** 本地 TF-MoDISco H5 文件路径。
- **输出：** 字典，含 `count`、`motifs`（每个基序的 `name`、`length` 与 `ppm`）与 `provenance`。

**示例代码：**

```python
import h5py
import numpy as np

with h5py.File("modisco_report.h5", "w") as handle:
    group = handle.create_group("pos_patterns").create_group("pattern_0")
    group.create_dataset("sequence", data=np.tile([4.0, 0.0, 0.0, 0.0], (6, 1)))
    group.create_dataset("contrib_scores", data=np.ones((6, 4)))
from dnakit.predictions import read_modisco_motifs

result = read_modisco_motifs("modisco_report.h5")
print(result["count"], result["motifs"][0]["name"], result["motifs"][0]["length"])
```

**示例输出：**

```text
1 pos_pattern_0 6
```

## 8) `score_motif_sites`

- **作用：** 为每个基序命中计算位点归因分数：命中区间内各碱基归因值的平均值。
- **处理方法：** 调用 `grelu.interpret.motifs.score_sites`。归因值与序列按列表顺序或同名键一一对应；命中区间必须落在序列之内。
- **API：** `score_motif_sites(hits[必须], sequences[必须], attributions[必须], verbose[可选])`
- **输入：** hits 为 `scan_motifs_with_pwms` 的命中列表；sequences 为序列；attributions 为每条序列一个 4 × L 矩阵。
- **输出：** 字典，含 `count`、`hits`（每个命中的 `site_attr_score`）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import score_motif_sites

sequences = {"s0": "ACGTAC"}
attributions = {"s0": [[1, 1, 1, 1, 0, 0], [0] * 6, [0] * 6, [0] * 6]}
hits = [{"motif": "m1", "sequence": "s0", "start": 0, "end": 4, "strand": "+"}]
result = score_motif_sites(hits, sequences, attributions)
print(result["hits"][0]["site_attr_score"])
```

**示例输出：**

```text
0.25
```

## 9) `score_motif_pwm_products`

- **作用：** 为每个基序命中计算 PWM 与归因的逐元素乘积之和，衡量基序形状与归因曲线的吻合程度。
- **处理方法：** 调用 `grelu.interpret.motifs.score_motifs`。负链命中使用反向互补的 PWM。PWM 长度必须等于命中区间长度。
- **API：** `score_motif_pwm_products(hits[必须], sequences[必须], attributions[必须], motifs[必须], verbose[可选])`
- **输入：** 与 `score_motif_sites` 相同，另加 motifs：基序名到 4 × L 概率矩阵的映射。
- **输出：** 字典，含 `count`、`hits`（每个命中的 `motif_attr_score`）与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import score_motif_pwm_products

pwm = [[1, 1], [0, 0], [0, 0], [0, 0]]
attributions = [[[2, 3, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]
hits = [{"motif": "m1", "sequence": "seq0", "start": 0, "end": 2, "strand": "+"}]
result = score_motif_pwm_products(hits, ["ACGT"], attributions, {"m1": pwm})
print(result["hits"][0]["motif_attr_score"])
```

**示例输出：**

```text
2.5
```

## 10) `run_modisco_motifs`

- **作用：** 用 TF-MoDISco 从模型的归因中发现序列模式。若找到模式，报告 `modisco_report.h5` 与标志图写入新的输出目录。
- **处理方法：** 调用 `grelu.interpret.modisco.run_modisco`。先计算指定任务的归因，再运行 TF-MoDISco；输出目录必须不存在或为空。随机序列通常找不到模式，此时不写出任何文件。
- **API：** `run_modisco_motifs(model[必须], sequences[必须], out_dir[必须], task_index[可选], window[可选], method[可选], n_shuffles[可选], seed[可选], verbose[可选])`
- **输入：** LightningModel 对象与等长 DNA 序列列表；out_dir 为输出目录。
- **输出：** 字典，含 `patterns_found`、`motifs`、`files` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
from dnakit.predictions import run_modisco_motifs

rng = random.Random(1)
sequences = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(40)]
result = run_modisco_motifs(model, sequences, out_dir="modisco_out", window=60, n_shuffles=2, seed=0)
print(result["patterns_found"], result["files"])
```

**示例输出：**

```text
False []
```

## 11) `get_sequence_attributions`

- **作用：** 计算一个任务上每条等长 DNA 序列的逐碱基归因，形状为 (B, 4, L)。
- **处理方法：** 调用 `grelu.interpret.score.get_attributions`，支持 deepshap、saliency、inputxgradient 与 integratedgradients。DeepSHAP 的二核苷酸打乱使用全局随机状态，包装函数会在调用前后恢复 NumPy、torch 与 Python 的随机状态。
- **API：** `get_sequence_attributions(model[必须], sequences[必须], task_index[可选], method[可选], n_shuffles[可选], seed[可选], correct_grad[可选], verbose[可选])`
- **输入：** LightningModel 对象与等长 DNA 序列。
- **输出：** 字典，含 `count`、`length`、`attributions` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import get_sequence_attributions

result = get_sequence_attributions(model, [SEQ], method="saliency")
print(result["count"], result["length"])
print([round(value, 4) for value in result["attributions"][0][0][:5]])
```

**示例输出：**

```text
1 60
[0.0029, 0.0025, -0.0017, -0.0013, -0.0002]
```

## 12) `in_silico_mutagenesis_scores`

- **作用：** 对一条序列做计算机模拟诱变：逐个位置、逐个替换碱基，计算某一任务预测值的变化。
- **处理方法：** 调用 `grelu.interpret.score.ISM_predict`。`compare_func` 取 subtract、divide 或 log2FC（突变体相对参考）；subtract 时参考碱基的分数为零。
- **API：** `in_silico_mutagenesis_scores(model[必须], sequence[必须], task_index[可选], compare_func[可选], start[可选], end[可选], batch_size[可选], verbose[可选])`
- **输入：** LightningModel 对象与一条 DNA 序列；start 与 end 限定位置范围（end 不含）。
- **输出：** 字典，含 `positions`、`bases` 与 `scores`（每个位置一行，列为 A、C、G、T）。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import in_silico_mutagenesis_scores

result = in_silico_mutagenesis_scores(model, SEQ, compare_func="subtract", start=0, end=3)
print(result["positions"], result["bases"])
for row in result["scores"]:
    print([round(value, 4) for value in row])
```

**示例输出：**

```text
3 ['A', 'C', 'G', 'T']
[0.0, -0.0029, -0.0032, -0.0023]
[0.0032, 0.0, 0.0023, 0.0009]
[0.0009, -0.0019, 0.0, -0.0026]
```

## 13) `get_attention_scores`

- **作用：** 返回 Transformer 模型对一条序列的注意力矩阵。
- **处理方法：** 调用 `grelu.interpret.score.get_attention_scores`。需要模型的 `n_transformers` 至少为 1；`block_idx` 为空时返回所有层。
- **API：** `get_attention_scores(model[必须], sequence[必须], block_idx[可选], verbose[可选])`
- **输入：** 含 Transformer 层的 LightningModel 对象与一条 DNA 序列。
- **输出：** 字典，含 `shape`、`attention` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
attention_model = LightningModel(
    model_params={
        "model_type": "ConvTransformerModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "n_transformers": 1,
        "n_heads": 2,
        "key_len": 4,
    },
    train_params={"task": "regression", "loss": "mse"},
)
attention_model.eval()
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import get_attention_scores

result = get_attention_scores(attention_model, SEQ, block_idx=0)
print(result["shape"])
print(round(sum(result["attention"][0][0]), 4))
```

**示例输出：**

```text
[2, 60, 60]
1.0
```

## 14) `marginalize_pattern_insertions`

- **作用：** 把 DNA 模式插入背景序列（可先做二核苷酸打乱），比较插入前后的预测值。
- **处理方法：** 调用 `grelu.interpret.simulate.marginalize_patterns`。`compare_func` 取 subtract 或 divide；`n_shuffles` 至少为 1。上游的原始输出模式与 `n_shuffles=0` 在本环境中会因 concat 错误失败，因此不对外提供。
- **API：** `marginalize_pattern_insertions(model[必须], sequences[必须], patterns[必须], task_index[可选], compare_func[可选], n_shuffles[可选], seed[可选], rc[可选], batch_size[可选], verbose[可选])`
- **输入：** LightningModel 对象、背景序列列表与模式列表；模式必须短于背景序列。
- **输出：** 字典，含 `shape`、`effects` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import marginalize_pattern_insertions

result = marginalize_pattern_insertions(model, [SEQ], ["GATTACA"], compare_func="subtract", n_shuffles=2, seed=0)
print(result["shape"])
values = [value for shuffle in result["effects"][0] for pattern in shuffle for task in pattern for value in task]
print([round(value, 4) for value in values])
```

**示例输出：**

```text
[1, 2, 1, 1, 1]
[-0.0051, -0.0008]
```

## 15) `marginalize_pattern_spacing`

- **作用：** 研究两个模式之间间距的影响：固定模式放在序列中央，移动模式按步长放在不同位置。
- **处理方法：** 调用 `grelu.interpret.simulate.marginalize_pattern_spacing`，返回间距列表与各间距的预测变化。
- **API：** `marginalize_pattern_spacing(model[必须], sequence[必须], fixed_pattern[必须], moving_pattern[必须], task_index[可选], stride[可选], n_shuffles[可选], seed[可选], rc[可选], compare_func[可选], batch_size[可选], verbose[可选])`
- **输入：** LightningModel 对象、一条背景序列与两个模式。
- **输出：** 字典，含 `spacings`、`shape`、`effects` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import marginalize_pattern_spacing

result = marginalize_pattern_spacing(model, SEQ, "GATTACA", "CCC", stride=20, n_shuffles=1, compare_func="subtract")
print(result["spacings"], result["shape"])
```

**示例输出：**

```text
[-26, -6, 14] [1, 1, 3, 1, 1]
```

## 16) `shuffle_tile_effects`

- **作用：** 把序列按固定长度切成片段，逐段做二核苷酸打乱，测量预测值的变化，用于寻找调控元件。
- **处理方法：** 调用 `grelu.interpret.simulate.shuffle_tiles`。片段起点每隔 `stride` 个碱基取一个（默认等于 `tile_len`）；`protect_center` 可保护中央区域不被打乱。
- **API：** `shuffle_tile_effects(model[必须], sequences[必须], tile_len[必须], task_index[可选], stride[可选], protect_center[可选], n_shuffles[可选], seed[可选], compare_func[可选], batch_size[可选], verbose[可选])`
- **输入：** LightningModel 对象与等长 DNA 序列列表。
- **输出：** 字典，含 `tile_count`、`tiles`（每段的 start 与 end）、`effects` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import shuffle_tile_effects

result = shuffle_tile_effects(model, [SEQ], 20, stride=20, n_shuffles=1, compare_func="subtract")
print(result["tile_count"], [(tile["start"], tile["end"]) for tile in result["tiles"]])
```

**示例输出：**

```text
3 [(0, 20), (20, 40), (40, 60)]
```

## 17) `grelu_predict_variant_effects`

- **作用：** 在本地参考基因组上预测每个 SNV 或插入缺失对某一任务的影响。
- **处理方法：** 调用 `grelu.variant.predict_variant_effects`。基因组为本地 FASTA，需有 `.fai` 与 `.sizes` 文件。上游的参考等位基因检查因旗标与函数同名而不可用，因此包装函数自行用 FASTA 核对参考碱基（`pos` 为 1 起算）。
- **API：** `grelu_predict_variant_effects(model[必须], variants[必须], genome_fasta[必须], task_index[可选], seq_len[可选], compare_func[可选], rc[可选], max_seq_shift[可选], check_reference[可选], batch_size[可选], verbose[可选])`
- **输入：** LightningModel 对象、变异映射列表（chrom、pos、ref、alt）与本地 FASTA 路径。
- **输出：** 字典，含 `count`、`variants`、`effects` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
rng = random.Random(9)
sequence = "".join(rng.choice("ACGT") for _ in range(300))
with open("genome.fa", "w") as handle:
    handle.write(">chr1\n" + sequence + "\n")
with open("genome.fa.fai", "w") as handle:
    handle.write("chr1\t300\t6\t300\t301\n")
with open("genome.fa.sizes", "w") as handle:
    handle.write("chr1\t300\n")
from dnakit.predictions import grelu_predict_variant_effects

ref = sequence[99]
alt = "A" if ref != "A" else "C"
result = grelu_predict_variant_effects(model, [{"chrom": "chr1", "pos": 100, "ref": ref, "alt": alt}], "genome.fa", seq_len=60)
print(result["count"], round(result["effects"][0][0], 4))
```

**示例输出：**

```text
1 0.9891
```

## 18) `grelu_marginalize_variants`

- **作用：** 对每个变异做边缘化检验：在打乱的背景序列上比较参考窗口与突变窗口的预测分布。
- **处理方法：** 调用 `grelu.variant.marginalize_variants`，基因组与变异规则同 `grelu_predict_variant_effects`。输出为上游的结果表。
- **API：** `grelu_marginalize_variants(model[必须], variants[必须], genome_fasta[必须], task_index[可选], seq_len[可选], n_shuffles[可选], seed[可选], compare_func[可选], rc[可选], max_seq_shift[可选], verbose[可选])`
- **输入：** LightningModel 对象、变异映射列表与本地 FASTA 路径。
- **输出：** 字典，含 `count`、`result`（上游结果表）与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
rng = random.Random(9)
sequence = "".join(rng.choice("ACGT") for _ in range(300))
with open("genome.fa", "w") as handle:
    handle.write(">chr1\n" + sequence + "\n")
with open("genome.fa.fai", "w") as handle:
    handle.write("chr1\t300\t6\t300\t301\n")
with open("genome.fa.sizes", "w") as handle:
    handle.write("chr1\t300\n")
from dnakit.predictions import grelu_marginalize_variants

ref = sequence[99]
alt = "A" if ref != "A" else "C"
result = grelu_marginalize_variants(model, [{"chrom": "chr1", "pos": 100, "ref": ref, "alt": alt}], "genome.fa", seq_len=60, n_shuffles=2, seed=0)
print(result["count"], sorted(result["result"]))
```

**示例输出：**

```text
1 ['effect_size', 'mean', 'pvalue', 'sd']
```

## 19) `design_sequences_evolve`

- **作用：** 对每条序列做贪心定向进化：逐步尝试单碱基替换（ism）或插入模式（pattern），保留提高某一任务预测值的候选。
- **处理方法：** 调用 `grelu.design.evolve`。模型必须在 `data_params['tasks']['name']` 中包含任务名，gReLU 在训练或加载模型时会设置它；缺少时包装函数拒绝执行。
- **API：** `design_sequences_evolve(sequences[必须], model[必须], task_index[可选], method[可选], patterns[可选], max_iter[可选], positions[可选], batch_size[可选], return_seqs[可选], verbose[可选])`
- **输入：** DNA 序列列表与 LightningModel 对象；positions 限定可改变的位置。
- **输出：** 字典，含 `count`、`records`（每轮每个候选的记录）与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import design_sequences_evolve

sequence = SEQ[:30]
result = design_sequences_evolve([sequence], model, max_iter=2, positions=[0, 1, 2])
changed = sorted({index for record in result["records"] for index, (a, b) in enumerate(zip(sequence, record["seq"])) if a != b})
print(result["count"], changed)
```

**示例输出：**

```text
19 [0, 1, 2]
```

## 20) `design_sequence_ledidi`

- **作用：** 用 Ledidi 设计一条序列的编辑，使某一任务的预测值升高；只允许改变指定位置。
- **处理方法：** 直接调用 `ledidi.Ledidi`，损失与 `grelu.design.ledidi` 相同（预测均值的负值）。上游把 0 维标量作为目标传入，ledidi 2.2.0 会因 `y_bar.shape[0]` 报 IndexError，因此包装函数传入长度为 1 的目标。`edit_penalty` 为 L1 编辑惩罚（上游默认 0.1）。全局 torch 随机状态会在调用后恢复。
- **API：** `design_sequence_ledidi(sequence[必须], model[必须], task_index[可选], max_iter[可选], batch_size[可选], positions[可选], lr[可选], edit_penalty[可选], seed[可选], verbose[可选])`
- **输入：** 一条 DNA 序列与 LightningModel 对象；positions 为允许改变的位置。
- **输出：** 字典，含 `designs`、`edit_counts`、`count` 与 `provenance`。

**示例代码：**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import design_sequence_ledidi

result = design_sequence_ledidi(SEQ, model, max_iter=30, positions=[0, 1, 2, 3, 4], edit_penalty=0.0)
print(result["count"], sorted(set(result["edit_counts"])))
changed = sorted({index for design in result["designs"] for index, (a, b) in enumerate(zip(SEQ, design)) if a != b})
print(changed)
```

**示例输出：**

```text
16 [3]
[1, 3, 4]
```
