# 表观基因组学

本类收录 **13** 项有文档的候选功能，其中 11 项在本地原生计算出真实结果（BED/bedGraph/Bismark 解析、区间合并与共识峰、注释、甲基化水平与质控、基于 Fisher 检验的差异甲基化、CpG 统计、归一化、绘图以及增强子-基因连接），2 项（`align`、`predict`）因需要比对程序或已训练的模型权重而始终返回 `status: "not_implemented"`。有些结果明确是启发式：`link_enhancers_to_genes` 是按到 TSS 距离的邻近规则，不是调控证据；`call` 和 `build_chromatin_state_segments` 只组合已给出的峰或状态判定（不调用峰，也不学习状态）。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) align

- **作用：** 把表观基因组读段比对到参考基因组。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算，也没有接入任何比对程序。亚硫酸氢盐数据需要 Bismark 或 bwa-meth，ChIP/ATAC 数据需要 Bowtie2 或 BWA，并且都需要带索引的参考基因组。不会编造结果。

- **API：** `dnakit.epigenome.align(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a read aligner (Bismark/bwa-meth for bisulfite, Bowtie2/BWA for ChIP/ATAC) and an indexed reference`，`message` 为 `Read alignment is not implemented natively; no aligner is wired into this function.`

- **示例代码：**

```python
import json
from dnakit.epigenome import align

result = align(["ACGTACGTACGTACGT"], options={"reference": "genome.fa"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "align",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "87cd4b17356c3dee4b13d54fba3beefea5678a7b651616a23b27c782da631af2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a read aligner (Bismark/bwa-meth for bisulfite, Bowtie2/BWA for ChIP/ATAC) and an indexed reference",
  "message": "Read alignment is not implemented natively; no aligner is wired into this function."
}
```

## 2) annotate

- **作用：** 用提供的注释，为基因组区域（如峰）标注重叠的特征和最近的特征。

- **处理方法：** 对每个区域扫描同一染色体上的特征（二次复杂度搜索，区域 x 特征对数至多 5000 万）。`feature.start < region.end` 且 `feature.end > region.start` 时重叠；`nearest` 是间隙最小的特征（重叠或相接时为 0，否则为两区间之间的碱基数；并列时取排序后第一个），`side` 为 `overlap`、`upstream`（特征在区域之前）或 `downstream`（在区域之后），均按正链坐标。其他染色体上的特征被忽略，也不使用基因链方向。注释就是提供的内容，不会去查询任何数据库。

- **API：** `dnakit.epigenome.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：待注释区域，区间：BED 文本/路径、行 `[chrom, start, end, ...]` 或带 `chrom`/`start`/`end` 的映射（0 起始半开，`end` > `start` >= 0；跳过 `track`、`browser`、`#` 行和空行；至多 200 万条记录）。`options`：`features`（必填：带名称列的 BED 注释、行或带 `name` 的映射）；`max_distance`（可选：更远的最近特征不报告）；`entity`（只支持 `"peaks"`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result.regions` 每个区域一项：`chrom`、`start`、`end`、`overlapping`（特征名称）和 `nearest`（`name`、`distance`、`side`；没有特征或超出 `max_distance` 时为 `null`）。

- **示例代码：**

```python
import json
from dnakit.epigenome import annotate

peaks = "chr1\t100\t200\tpeak1\nchr1\t1000\t1100\tpeak2\n"
features = "chr1\t150\t400\tgeneA\nchr1\t5000\t6000\tgeneB\n"
result = annotate(peaks, options={"features": features, "max_distance": 10000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1e0b18cc52249734849452cd1247dd17a4569b367d0b45b2f578d88f9e3f1d25"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "regions": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "overlapping": [
          "geneA"
        ],
        "nearest": {
          "name": "geneA",
          "distance": 0,
          "side": "overlap"
        }
      },
      {
        "chrom": "chr1",
        "start": 1000,
        "end": 1100,
        "overlapping": [],
        "nearest": {
          "name": "geneA",
          "distance": 600,
          "side": "upstream"
        }
      }
    ]
  }
}
```

## 3) build_chromatin_state_segments

- **作用：** 把逐 bin 的染色质状态判定合并为连续的状态片段。

- **处理方法：** bin 按染色体、起点、终点排序；同一染色体上状态相同的相邻 bin，与前一片段的间隙 <= `max_gap` 个碱基时合并（状态相同且重叠的 bin 也合并）。它只是对已有的状态判定做分段，不学习也不调用状态（不是 ChromHMM）。

- **API：** `dnakit.epigenome.build_chromatin_state_segments(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：bin 为 BED 文本/路径、行或映射 `chrom, start, end, name`，其中 `name` 是状态标签（映射中也接受别名键 `state`），或 `{chrom: [状态, 状态, ...]}` 并配合 `options["bin_size"]`（bin 从 0 开始连续排列）。`options`：`max_gap`（默认 0）、`bin_size`（`{chrom: [状态]}` 形式必填）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `bin_count`、`segment_count` 和按位置排序的 `segments`（`chrom`、`start`、`end`、`state`、`n_bins`）。

- **示例代码：**

```python
import json
from dnakit.epigenome import build_chromatin_state_segments

states = {"chr1": ["Quiescent", "Enhancer", "Enhancer", "Promoter", "Promoter", "Promoter"]}
result = build_chromatin_state_segments(states, options={"bin_size": 200})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "build_chromatin_state_segments",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 83,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b24877c056115a3c3ed7a4314d1dc9a56d57e546ddbf588405cab52d2608c84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_count": 6,
    "segment_count": 3,
    "segments": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 200,
        "state": "Quiescent",
        "n_bins": 1
      },
      {
        "chrom": "chr1",
        "start": 200,
        "end": 600,
        "state": "Enhancer",
        "n_bins": 2
      },
      {
        "chrom": "chr1",
        "start": 600,
        "end": 1200,
        "state": "Promoter",
        "n_bins": 3
      }
    ]
  }
}
```

## 4) calculate

- **作用：** 根据计数计算甲基化水平，根据 DNA 计算 CpG 观测/期望比，或计算两组区间的重叠（`options["entity"]`）。

- **处理方法：** `methylation_level`（默认）：每个位点的水平 = 甲基化 / (甲基化 + 未甲基化)；覆盖度低于 `min_coverage` 的位点被丢弃；`one_based` 把 Bismark 的 1 起始位置移为 0 起始。给定 `regions` 时，按 0 起始半开区间聚合区域内的位点：`weighted_level` = sum(甲基化) / sum(覆盖度)，`mean_site_level` = 位点水平的未加权平均。`cpg_observed_expected`：统计 C、G、CpG 数量、GC 比例和 CpG o/e = CpG * L / (C * G)（C*G = 0 时为 `null`），设置 `window` 时还按滑动窗口计算。`interval_overlap`：先合并每个集合，再按碱基数计算 `overlap_bp`、`union_bp` 和 Jaccard = 重叠 / 并集。`a_intervals_overlapping_b` 统计与 `b` 的任一区间至少共享一个碱基的原始（未合并）输入区间个数（半开区间，首尾相接不算相交）；每个输入区间只计一次，无论它跨越多少个 `b` 区间。

- **API：** `dnakit.epigenome.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `methylation_level`：`input` 为逐位点计数：Bismark `.cov` 文本/路径或行 `chrom, start, end, percent, methylated, unmethylated`，或带 `chrom`、`start`、`methylated` 和 `unmethylated`（或 `coverage`）的映射；`options`：`min_coverage`（默认 1）、`one_based`（默认 `false`；`.cov` 请用 `true`）、`regions`（可选 BED）。`cpg_observed_expected`：`input` 为 DNA（字符串、FASTA 文本/路径或 `{id: 序列}`，ACGT；至多 10000 条序列、共 50 Mb）；`options`：`window`（可选）、`step`（默认等于 window）。`interval_overlap`：`input` 和 `options["b"]`（必填）都是区间集合（BED 文本/路径、行或映射）。`options["entity"]` 取 `"methylation_level"`（默认）、`"cpg_observed_expected"` 或 `"interval_overlap"`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result.entity` 回显所选类型。无 `regions` 的 `methylation_level`：`sites`（`chrom`、`pos`、`coverage`、`methylated`、`level`）、`site_count`、`dropped_low_coverage`；有 `regions` 时：`regions`（`chrom`、`start`、`end`、`n_sites`、`coverage`、`weighted_level`、`mean_site_level`）。`cpg_observed_expected`：`sequences`（`id`、`length`、`c`、`g`、`cpg`、`gc_fraction`、`observed_expected`，要求时还有 `windows`）。`interval_overlap`：`overlap_bp`、`a_bp`、`b_bp`、`union_bp`、`jaccard`、`a_intervals_overlapping_b`。

- **示例代码：**

```python
import json
from dnakit.epigenome import calculate

cov = "chr1\t101\t101\t66.7\t2\t1\nchr1\t120\t120\t100\t5\t0\nchr1\t5000\t5000\t0\t0\t4\n"   # Bismark .cov rows
result = calculate(cov, options={"entity": "methylation_level", "one_based": True,
                                 "regions": "chr1\t0\t1000\tregionA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "22911cbbf5fca89bad46f6d2e400272575fc4bd02d8613cfdf73a18d901f0f50"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "methylation_level",
    "regions": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 1000,
        "n_sites": 2,
        "coverage": 8,
        "weighted_level": 0.875,
        "mean_site_level": 0.8333333333333333
      }
    ]
  }
}
```

## 5) call

- **作用：** 跨样本调用共识峰：被至少 `min_samples` 个样本的峰集覆盖的区域。

- **处理方法：** 先合并每个样本自身重叠的峰，使一个样本在每个碱基上至多投一票；结果是票数 >= `min_samples` 的每个最大片段（0 起始半开），通过对区间起止事件的扫描线计算。它只是组合已有的峰调用结果，不从信号调用峰（不是 MACS），也不使用任何统计。

- **API：** `dnakit.epigenome.call(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：峰集，可以是映射 `{样本: 区间}`（每个为 BED 文本/路径、行或映射）或区间集合列表。`options`：`entity`（只支持 `"consensus_peaks"`）；`min_samples`（默认 2，不超过样本数）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `samples`、`min_samples`、`peak_count` 和 `peaks`（`chrom`、`start`、`end`）。

- **示例代码：**

```python
import json
from dnakit.epigenome import call

peak_sets = {
    "rep1": "chr1\t100\t300\nchr1\t1000\t1200\n",
    "rep2": "chr1\t200\t400\n",
    "rep3": "chr1\t250\t350\nchr1\t1100\t1300\n",
}
result = call(peak_sets, options={"entity": "consensus_peaks", "min_samples": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "call",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 120,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6e5f5041ddc207500e99845fb2c45b47675903c48b38032e44d579799a510e75"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "rep1",
      "rep2",
      "rep3"
    ],
    "min_samples": 2,
    "peak_count": 2,
    "peaks": [
      {
        "chrom": "chr1",
        "start": 200,
        "end": 350
      },
      {
        "chrom": "chr1",
        "start": 1100,
        "end": 1200
      }
    ]
  }
}
```

## 6) find

- **作用：** 用 Fisher 精确检验寻找差异甲基化区域，或在 DNA 中寻找 CpG 岛（`options["entity"]`）。

- **处理方法：** `differential_regions`（默认）：对每个区域，A、B 两组的汇总计数给出各组水平、`delta` = level_b - level_a，并对 2x2 表 [[methylated_a, unmethylated_a], [methylated_b, unmethylated_b]] 做双侧 Fisher 精确检验（scipy 的定义：概率不大于观测表的所有表的概率之和），同时给出样本比值比 ad/bc（bc = 0 时为 `null`）以及对被检验区域的 Benjamini-Hochberg `q_value`。任一组覆盖度低于 `min_coverage` 的区域被排除；`significant` 表示 q <= `alpha` 且 |delta| >= `min_delta`。计数在重复和胞嘧啶之间合并，所以没有对重复间变异建模（不是 beta-二项式/DSS 类模型）。`cpg_islands`：调用 `dnakit.patterns.find_cpg_islands`，即滑动窗口规则（默认窗口 200、步长 1、GC >= 0.5、CpG 观测/期望 >= 0.6、区域长度 >= 200）。

- **API：** `dnakit.epigenome.find(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `differential_regions`：`input` 为带 `chrom`、`start`、`end`（可选）以及汇总计数 `methylated_a`、`unmethylated_a`、`methylated_b`、`unmethylated_b` 的行或映射；给定 `options["regions"]`（BED）时，输入行是单个位点，在每个区域内求和。`options`：`min_coverage`（默认 10，按组）、`alpha`（默认 0.05）、`min_delta`（默认 0.0）。`cpg_islands`：`input` 为一条 DNA 序列（字符串、FASTA 文本/路径或 `{id: 序列}`）；`options`：`window_size`（200）、`step`（1）、`min_gc`（0.5）、`min_observed_expected`（0.6）、`min_region_length`（200）。`options["entity"]` 取 `"differential_regions"`（默认）或 `"cpg_islands"`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`differential_regions`：`result` 含 `entity`、`tested`、`excluded_low_coverage` 和按 p 值排序的 `regions`（`chrom`、`start`、`end`、`level_a`、`level_b`、`delta`、`coverage_a`、`coverage_b`、`p_value`、`odds_ratio`、`q_value`、`significant`）。`cpg_islands`：`entity`、`sequence` 和 `islands`（`start`、`end`、`length`、`observed_expected`、`gc_fraction`）。

- **示例代码：**

```python
import json
from dnakit.epigenome import find

regions = [
    {"chrom": "chr1", "start": 100, "end": 200,
     "methylated_a": 40, "unmethylated_a": 60, "methylated_b": 80, "unmethylated_b": 20},
    {"chrom": "chr1", "start": 500, "end": 600,
     "methylated_a": 50, "unmethylated_a": 50, "methylated_b": 52, "unmethylated_b": 48},
]
result = find(regions, options={"entity": "differential_regions"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "find",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 258,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "759876ad6dfe28bf2acab3945ddd4e23a745b2f8098b355e76502f0a3913a0a7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "differential_regions",
    "tested": 2,
    "excluded_low_coverage": 0,
    "regions": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "level_a": 0.4,
        "level_b": 0.8,
        "delta": 0.4,
        "coverage_a": 100,
        "coverage_b": 100,
        "p_value": 1.0636025781775309e-08,
        "odds_ratio": 0.16666666666666666,
        "q_value": 2.1272051563550617e-08,
        "significant": true
      },
      {
        "chrom": "chr1",
        "start": 500,
        "end": 600,
        "level_a": 0.5,
        "level_b": 0.52,
        "delta": 0.020000000000000018,
        "coverage_a": 100,
        "coverage_b": 100,
        "p_value": 0.8875621593770635,
        "odds_ratio": 0.9230769230769231,
        "q_value": 0.8875621593770635,
        "significant": false
      }
    ]
  }
}
```

## 7) link_enhancers_to_genes

- **作用：** 按基因组距离，把增强子（或任意区域）分配给基因的转录起始位点（TSS）。

- **处理方法：** 基于距离的分配：距离为 |TSS - 中点|，中点 = 增强子的 (start + end) // 2。TSS 取显式给出的 `tss`，否则对 + 链或无链方向的基因取 `start`，对 - 链基因取 `end - 1`。`method="nearest"` 把每个增强子连到同一染色体上最近的单个 TSS（距离并列时取名称字母序在前的基因），前提是在 `max_distance` 之内；`method="window"` 连接 TSS 在 `max_distance` 之内的所有基因。`direction` 表示增强子相对基因方向位于上游、下游还是 TSS 处。这是邻近性启发式，不是调控证据：不使用 Hi-C、表达相关性或 TAD 边界，调用会附带相应警告。

- **API：** `dnakit.epigenome.link_enhancers_to_genes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：增强子，区间：BED 文本/路径、行 `[chrom, start, end, ...]` 或带 `chrom`/`start`/`end` 的映射（0 起始半开，`end` > `start` >= 0；跳过 `track`、`browser`、`#` 行和空行；至多 200 万条记录）。`options`：`genes`（必填：带 `chrom`、`start`、`end`、`name` 以及可选 `strand` 或显式 `tss` 的行或映射）；`method`（`"nearest"` 默认 | `"window"`）；`max_distance`（默认 1,000,000 bp；`nearest` 时更远的基因不连接）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `method`、`max_distance` 和 `links`（每个增强子一项：`region`、`midpoint` 和 `genes`，其中每个基因有 `gene`、`tss`、`distance`、`direction`、`strand`；范围内没有基因时列表为空）。

- **示例代码：**

```python
import json
from dnakit.epigenome import link_enhancers_to_genes

enhancers = "chr1\t10000\t10400\tenh1\nchr1\t90000\t90200\tenh2\n"
genes = [
    {"name": "GENE_A", "chrom": "chr1", "start": 12000, "end": 20000, "strand": "+"},
    {"name": "GENE_B", "chrom": "chr1", "start": 60000, "end": 91000, "strand": "-"},
]
result = link_enhancers_to_genes(enhancers, options={"genes": genes, "method": "nearest", "max_distance": 50000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "link_enhancers_to_genes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b35ce0cf833e4ae26ca3cce75f823f7c9a74d1755641cf0d0f0eb456361e590c"
  },
  "warnings": [
    "Distance-based assignment is a proximity heuristic, not evidence of regulation."
  ],
  "status": "ok",
  "result": {
    "method": "nearest",
    "max_distance": 50000,
    "links": [
      {
        "region": {
          "chrom": "chr1",
          "start": 10000,
          "end": 10400,
          "name": "enh1"
        },
        "midpoint": 10200,
        "genes": [
          {
            "gene": "GENE_A",
            "tss": 12000,
            "distance": 1800,
            "direction": "upstream",
            "strand": "+"
          }
        ]
      },
      {
        "region": {
          "chrom": "chr1",
          "start": 90000,
          "end": 90200,
          "name": "enh2"
        },
        "midpoint": 90100,
        "genes": [
          {
            "gene": "GENE_B",
            "tss": 90999,
            "distance": 899,
            "direction": "downstream",
            "strand": "-"
          }
        ]
      }
    ]
  }
}
```

## 8) merge

- **作用：** 按染色体合并重叠或相近的区间（峰、区域），类似 `bedtools merge`。

- **处理方法：** 区间按染色体排序，与当前合并区间的间隙 <= `distance` 时合并（因此 `distance` 为 0 时，重叠和首尾相接的区间都会合并）。忽略链方向和名称。输出按染色体名和起点排序，`n_merged` 为合并的输入区间数。

- **API：** `dnakit.epigenome.merge(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：区间：BED 文本/路径、行 `[chrom, start, end, ...]` 或带 `chrom`/`start`/`end` 的映射（0 起始半开，`end` > `start` >= 0；跳过 `track`、`browser`、`#` 行和空行；至多 200 万条记录）。`options`：`entity`（只支持 `"peaks"`）；`distance`（默认 0）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `distance`、`intervals`（`chrom`、`start`、`end`、`n_merged`）和 `count`。

- **示例代码：**

```python
import json
from dnakit.epigenome import merge

peaks = "chr1\t100\t200\nchr1\t150\t300\nchr1\t310\t400\nchr2\t5\t50\n"
result = merge(peaks, options={"distance": 20})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "merge",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 49,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c1640184d9b52456608b48132dcfc54c62503331851a4b0255849c89789202c3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "distance": 20,
    "intervals": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 400,
        "n_merged": 3
      },
      {
        "chrom": "chr2",
        "start": 5,
        "end": 50,
        "n_merged": 1
      }
    ],
    "count": 2
  }
}
```

## 9) normalize

- **作用：** 对信号值（bedGraph 记录或普通数字列表）做归一化。

- **处理方法：** `cpm`（默认）：值 / 总和 * 1e6；`rpkm`：cpm / (区间长度 / 1000)，需要区间；`log2cpm`：log2(cpm + `pseudocount`)；`zscore`：(v - 均值) / 标准差，默认总体标准差，`ddof` 为 1 时用样本标准差；`minmax`：(v - 最小值) / (最大值 - 最小值)。cpm 系列要求值非负且总和为正，z-score 要求方差不为 0，min-max 要求值不全相同，否则抛出输入错误。归一化只作用于提供的值，不做文库大小或 GC 校正。

- **API：** `dnakit.epigenome.normalize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：bedGraph 文本/路径、行 `[chrom, start, end, value]` 或带 `value`（或 `signal`）字段的映射，或数字列表（`rpkm` 不可用）。`options`：`entity`（只支持 `"signal"`）；`method`（`"cpm"` 默认 | `"rpkm"` | `"log2cpm"` | `"zscore"` | `"minmax"`）；`pseudocount`（默认 1，用于 `log2cpm`）；`ddof`（0 或 1，默认 0，用于 `zscore`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `method` 和 `records`（输入记录，其中 `value` 被替换为归一化后的值，`raw_value` 保留原值）。

- **示例代码：**

```python
import json
from dnakit.epigenome import normalize

signal = "chr1\t0\t100\t5\nchr1\t100\t200\t15\nchr1\t200\t400\t30\n"
result = normalize(signal, options={"method": "rpkm"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 45,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fa1a617de98e7c794ff977eaf8ac2571d21211bfe540f311c20841acfddc1a81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "rpkm",
    "records": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 100,
        "value": 1000000.0,
        "raw_value": 5.0
      },
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "value": 3000000.0,
        "raw_value": 15.0
      },
      {
        "chrom": "chr1",
        "start": 200,
        "end": 400,
        "value": 3000000.0,
        "raw_value": 30.0
      }
    ]
  }
}
```

## 10) plot

- **作用：** 把一个基因组窗口内的信号轨道（bedGraph bin 作为柱子）画成 SVG 字符串。

- **处理方法：** 手写的 SVG（不依赖绘图库）：与窗口重叠的每个 bedGraph bin 画成带 `data-start`、`data-end`、`data-value` 属性的 `<rect class="bin">`，在窗口内线性定位，高度按窗口内最大值缩放。窗口标签显示最大值。只绘制第一条染色体（或 `chrom` 指定的染色体）。

- **API：** `dnakit.epigenome.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：bedGraph 文本/路径、行 `[chrom, start, end, value]` 或带 `value`（或 `signal`）字段的映射。`options`：`entity`（只支持 `"tracks"`）；`chrom`、`start`、`end`（默认取第一条染色体及其 bin 的跨度）；`output`（可选 `.svg` 路径）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `svg`（SVG 文本）和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.epigenome import plot

signal = "chr1\t0\t100\t5\nchr1\t100\t200\t15\nchr1\t200\t400\t30\n"
result = plot(signal, options={"chrom": "chr1", "start": 0, "end": 400})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "epigenome",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 45,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fa1a617de98e7c794ff977eaf8ac2571d21211bfe540f311c20841acfddc1a81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"150\"><text x=\"20\" y=\"12\" font-size=\"10\">chr1:0-400 (max 30)</text><rec...",
    "output_path": null
  }
}
```

## 11) predict

- **作用：** 由序列预测表观基因组信号。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。由序列预测表观基因组需要 Enformer、Basenji、DeepSEA 之类的已训练模型及其权重，DNAKit 并未内置。不会编造结果。

- **API：** `dnakit.epigenome.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a trained sequence-to-epigenome model and its weights (Enformer, Basenji, DeepSEA)`，`message` 为 `Epigenome prediction needs trained model weights that are not bundled.`

- **示例代码：**

```python
import json
from dnakit.epigenome import predict

result = predict("ACGTACGTACGTACGT")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 16,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "cf573e65038d08ff910a3345642ffd1e8329844633c2dcb15964b324ebdba4d0"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained sequence-to-epigenome model and its weights (Enformer, Basenji, DeepSEA)",
  "message": "Epigenome prediction needs trained model weights that are not bundled."
}
```

## 12) qc

- **作用：** 检查亚硫酸氢盐甲基化计数的质量：覆盖度分布、整体甲基化水平和转化率。

- **处理方法：** 先把对照染色体（`control_chroms`，如未甲基化的 lambda 噬菌体或 chrL）与样本位点分开。对样本位点：位点数、总覆盖度、平均和中位覆盖度、覆盖度 >= `min_coverage` 的位点比例，以及 `global_methylation` = sum(甲基化) / sum(覆盖度)。给定对照时，`conversion_rate` = 1 - 对照位点上的甲基化 / 覆盖度，即未甲基化胞嘧啶被读作 T 的比例。只是汇总统计，除报告外不设阈值。

- **API：** `dnakit.epigenome.qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：逐位点计数：Bismark `.cov` 文本/路径或行 `chrom, start, end, percent, methylated, unmethylated`，或带 `chrom`、`start`、`methylated` 和 `unmethylated`（或 `coverage`）的映射。`options`：`entity`（只支持 `"methylation"`）；`min_coverage`（默认 10）；`control_chroms`（可选，对照染色体名称列表）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `sites`、`total_coverage`、`mean_coverage`、`median_coverage`、`min_coverage`、`fraction_sites_min_coverage`、`global_methylation`，给定对照时还有 `control_coverage` 和 `conversion_rate`。去掉对照后没有样本位点时抛出输入错误。

- **示例代码：**

```python
import json
from dnakit.epigenome import qc

cov = ("chr1\t101\t101\t66.7\t20\t10\nchr1\t120\t120\t100\t30\t0\nchr1\t5000\t5000\t0\t0\t4\n"
       "chrL\t10\t10\t1.0\t1\t99\n")
result = qc(cov, options={"min_coverage": 10, "control_chroms": ["chrL"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 87,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d643d3abac8dc0545d26854a692a340ad0c42126f53306d4fc3057691324ee32"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sites": 3,
    "total_coverage": 64,
    "mean_coverage": 21.333333333333332,
    "median_coverage": 30.0,
    "min_coverage": 10,
    "fraction_sites_min_coverage": 0.6666666666666666,
    "global_methylation": 0.78125,
    "control_coverage": 100,
    "conversion_rate": 0.99
  }
}
```

## 13) read

- **作用：** 读取 BED、bedGraph 或 Bismark 覆盖度（逐位点甲基化计数）的文本、文件路径或内存中的行。

- **处理方法：** 跳过 `track`、`browser`、`#` 行和空行；按制表符（或空白）切分。校验记录（整数坐标、`end` > `start` >= 0、有限数值、链方向为 `+`/`-`/`.`），格式错误的行会抛出带行号的错误。bed/bedgraph 坐标为 0 起始半开；`bismark_cov` 的位置按原样返回（1 起始），且只保留 chrom、start、methylated 和 unmethylated。

- **API：** `dnakit.epigenome.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本、路径、单元格行或映射。`options`：`format`（`"bed"` 默认：chrom、start、end，可选 name/score/strand | `"bedgraph"`：chrom、start、end、value | `"bismark_cov"`：chrom、start、end、percent、methylated、unmethylated）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"epigenome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `format`、`count` 和 `records`。

- **示例代码：**

```python
import json
from dnakit.epigenome import read

bedgraph = "track type=bedGraph\nchr1\t0\t100\t1.5\nchr1\t100\t200\t3\n"
result = read(bedgraph, options={"format": "bedgraph"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "epigenome",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50e7b9157f4acdeffe27e5b4d0ac53b9255932da665686b2a6158b7f2dad752d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "bedgraph",
    "count": 2,
    "records": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 100,
        "value": 1.5
      },
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "value": 3.0
      }
    ]
  }
}
```
