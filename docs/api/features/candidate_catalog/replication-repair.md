# DNA 复制、重组与修复

本类收录 **10** 项功能。10 项均以纯 Python 在本地实际计算并返回 `status: "ok"`（不依赖外部程序或可选库），没有返回 `not_implemented` 的函数。所用科学模型较为简单，其中一部分明确为启发式：基于 GC 偏斜的起点/终点估计、DnaA 盒密度、基于插入缺失特征的修复途径规则以及确定性的复制叉模型，都不能替代复制或修复的实验数据。这些函数通过正式的 DNAKit 命名空间导出，JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。下面的示例输出均为真实运行结果，较长的数组和字符串以 `...` 省略。

## 1) annotate

- **作用：** 按环状染色体的复制起点和终点，将基因组特征标注为前导链或滞后链方向。

- **处理方法：** 染色体按环状处理，坐标为 0 起始、半开区间。右复制子沿坐标增大方向从起点走向终点，左复制子从终点回到起点。特征（取其中点）位于右复制子的 + 链或左复制子的 - 链时，与复制叉同向（前导链），否则为迎面（滞后链）。起点和终点取自输入或选项中的 `origin`/`terminus`；未给出起点时，由 `sequence` 的累积 GC 偏斜估计（启发式，见 `calculate`），终点默认为起点的对跖点。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为映射，含 `features`（`{id, start, end, strand}` 列表，`strand` 取 `+` 或 `-`），以及 `length`（染色体长度）或 `sequence`（DNA 字符串）之一，另可含 `origin` 和 `terminus`。选项：`origin`、`terminus`（覆盖输入中的值；默认由 `sequence` 的 GC 偏斜估计起点，终点为起点的对跖点）。`options` 可省略。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `length`、`origin`、`terminus`、`origin_source`、两条复制子的长度、`features`（每项含 `replichore`、`orientation`、`leading_strand`、`distance_from_origin_along_replichore`）、`leading_strand_features`、`lagging_strand_features`、`leading_strand_fraction` 以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import annotate

chromosome = {
    "length": 1000,
    "origin": 0,
    "terminus": 500,
    "features": [
        {"id": "g1", "start": 100, "end": 200, "strand": "+"},
        {"id": "g2", "start": 300, "end": 400, "strand": "-"},
        {"id": "g3", "start": 700, "end": 800, "strand": "-"},
    ],
}
result = annotate(chromosome)
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 225,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "453b7bb62ffd20dd42f2b582bfd5ca17c525b4b2bd305999ddc7c9b05f6b1f77"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 1000,
    "origin": 0,
    "terminus": 500,
    "origin_source": "supplied",
    "right_replichore_length": 500,
    "left_replichore_length": 500,
    "features": [
      {
        "id": "g1",
        "start": 100,
        "end": 200,
        "strand": "+",
        "replichore": "right",
        "orientation": "co-directional (leading strand)",
        "leading_strand": true,
        "distance_from_origin_along_replichore": 150
      },
      {
        "id": "g2",
        "start": 300,
        "end": 400,
        "strand": "-",
        "replichore": "right",
        "orientation": "head-on (lagging strand)",
        "leading_strand": false,
        "distance_from_origin_along_replichore": 350
      },
      {
        "id": "g3",
        "start": 700,
        "end": 800,
        "strand": "-",
        "replichore": "left",
        "orientation": "co-directional (leading strand)",
        "leading_strand": true,
        "distance_from_origin_along_replichore": 250
      }
    ],
    "leading_strand_features": 2,
    "lagging_strand_features": 1,
    "leading_strand_fraction": 0.666667,
    "definition": "Circular chromosome; the right replichore runs origin -> terminus in increasing coordinates, the left replichore terminus -> origin. A feature (position = midpo..."
  }
}
```

## 2) build_recombination_map

- **作用：** 由有序标记构建遗传（重组）图谱：标记间的厘摩距离及每个区间的 cM/Mb 重组率。

- **处理方法：** 每个非首位标记到前一标记的距离取自给定的 `cm`（累积 cM）、重组率 `r`，或 `r = recombinants / meioses`。重组率按 Haldane 图谱函数 d = -1/2 ln(1 - 2r)（默认）、Kosambi 图谱函数 d = 1/4 ln((1 + 2r) / (1 - 2r)) 转换，或保持 d = r（`none`）；Haldane 与 Kosambi 要求 r 小于 0.5。距离乘以 100 得 cM。区间重组率为该区间 delta cM 除以区间长度（Mb）；结果还给出累积 cM、图谱总长和平均重组率（总 cM 除以跨度 Mb）。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.build_recombination_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为按位置排序、至少含两个标记的列表（或 `{"markers": [...]}`），每个标记形如 `{id, position}`（bp，非负整数）。除首个标记外，每个标记须给出 `cm`（累积 cM）、`recombination_fraction`（到前一标记的 r）或 `recombinants` 与 `meioses` 之一。首个标记可带起始 `cm`。选项：`map_function`（`haldane` 默认、`kosambi`、`none`），仅作用于基于 r 的标记。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `map_function`、`markers`、`total_cm`、`span_bp`、`mean_rate_cm_per_mb`、`map`（每个标记一行，含 `cumulative_cm`；第二个标记起另有 `interval_start`、`interval_end`、`recombination_fraction`、`delta_cm`、`rate_cm_per_mb`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import build_recombination_map

markers = [
    {"id": "m1", "position": 0},
    {"id": "m2", "position": 1_000_000, "recombination_fraction": 0.1},
    {"id": "m3", "position": 3_000_000, "recombinants": 5, "meioses": 50},
]
result = build_recombination_map(markers, options={"map_function": "kosambi"})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "build_recombination_map",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 164,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d84d990f59e801eaea2a20335fc868493341e19a259092cf7d02ba067d4f388d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "map_function": "kosambi",
    "markers": 3,
    "total_cm": 20.273256,
    "span_bp": 3000000,
    "mean_rate_cm_per_mb": 6.757752,
    "map": [
      {
        "id": "m1",
        "position": 0,
        "cumulative_cm": 0.0
      },
      {
        "id": "m2",
        "position": 1000000,
        "interval_start": 0,
        "interval_end": 1000000,
        "recombination_fraction": 0.1,
        "delta_cm": 10.136628,
        "cumulative_cm": 10.136628,
        "rate_cm_per_mb": 10.136628
      },
      {
        "id": "m3",
        "position": 3000000,
        "interval_start": 1000000,
        "interval_end": 3000000,
        "recombination_fraction": 0.1,
        "delta_cm": 10.136628,
        "cumulative_cm": 20.273255,
        "rate_cm_per_mb": 5.068314
      }
    ],
    "definition": "Interval rate = delta cM / (interval length in Mb). Distances from recombination fractions r use the Haldane (d = -1/2 ln(1-2r)) or Kosambi (d = 1/4 ln((1+2r)/(..."
  }
}
```

## 3) calculate

- **作用：** 计算 DNA 序列的窗口 GC（或 AT）偏斜与累积偏斜，并启发式估计复制起点和终点。

- **处理方法：** 窗口偏斜为 (G - C) / (G + C)（`skew_type` 为 `GC`）或 (A - T) / (A + T)（`AT`）；累积偏斜为窗口偏斜的累加和。启发式（Lobry 1996、Grigoriev 1998）：起点取累积偏斜最小的窗口中心，终点取最大的窗口中心。仅适用于存在链组成偏好的类细菌染色体，分辨率取决于窗口大小。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 DNA 字符串、单条记录的 FASTA 文本、文件路径或 `{"sequence": ...}`（最长 20 Mb）。选项：`window`（默认 min(1000, 序列长度)）、`step`（默认等于 `window`）、`skew_type`（`GC` 默认或 `AT`）、`include_windows`（默认 `True`）、`max_windows_returned`（默认 20000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `skew_type`、`length`、`window`、`step`、`windows`、`overall_skew`、`origin_estimate`（及 `origin_window`）、`terminus_estimate`（及 `terminus_window`）、`cumulative_skew_range`、`heuristic`、`definition` 以及 `window_table`（每个窗口的 `start`、`end`、`skew`、`cumulative_skew` 和碱基计数）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import calculate

sequence = "GGGC" * 125 + "CCCG" * 125
result = calculate(sequence, options={"window": 100})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 1000,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "3490f072d0775230a9eb8cb9e2bc13b5ea9c4ae2878fa5265873210e8ad3c480"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "skew_type": "GC",
    "length": 1000,
    "window": 100,
    "step": 100,
    "windows": 10,
    "overall_skew": 0.0,
    "origin_estimate": 950,
    "origin_window": [
      900,
      1000
    ],
    "terminus_estimate": 450,
    "terminus_window": [
      400,
      500
    ],
    "cumulative_skew_range": 2.5,
    "heuristic": "Origin = centre of the window with the minimum cumulative skew, terminus = maximum (Lobry 1996; Grigoriev 1998). Valid only for circular bacterial-like chromoso...",
    "definition": "skew = (G-C)/(G+C) per window; cumulative skew = running sum of the window skews.",
    "window_table": [
      {
        "start": 0,
        "end": 100,
        "skew": 0.5,
        "cumulative_skew": 0.5,
        "count_G": 75,
        "count_C": 25
      },
      {
        "start": 100,
        "end": 200,
        "skew": 0.5,
        "cumulative_skew": 1.0,
        "count_G": 75,
        "count_C": 25
      },
      "..."
    ]
  }
}
```

## 4) classify_repair_pathway

- **作用：** 依据参考序列与修复后序列，用插入缺失特征规则判断可能的 DNA 双链断裂修复途径（NHEJ、MMEJ、SSA 或 HDR）。

- **处理方法：** 这是启发式规则，并非对实际使用途径的测量。修复后序列与提供的 `donor` 相同时判为 HDR。纯缺失且两侧断点共享至少 `ssa_min_homology` 个相同碱基时判为 SSA，至少 `min_microhomology` 个时判为 MMEJ，否则为 NHEJ。长度不超过 `max_nhej_insertion` bp 的纯插入判为 NHEJ。碱基替换、复杂插入缺失和更长的插入为 `undetermined`，序列完全相同则为 `no_edit`。断点同源长度指缺失可滑动的范围：前缀匹配长度 + 后缀匹配长度 - 修复后序列长度。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.classify_repair_pathway(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为映射 `{"reference": 未编辑序列, "edited": 修复后序列}`，可含 `donor`（预期的 HDR 产物）。选项：`min_microhomology`（默认 2）、`ssa_min_homology`（默认 25）、`max_nhej_insertion`（默认 20）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `pathway`（`NHEJ`、`MMEJ`、`SSA`、`HDR`、`no_edit` 或 `undetermined`）、`reason`、`rule`（规则说明），以及序列变化字段 `reference_length`、`edited_length`、`deleted_bases`、`inserted_bases`、`first_difference`，缺失时还有 `microhomology_length`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import classify_repair_pathway

pair = {"reference": "GGATGCATGCCGG", "edited": "GGATGCCGG"}
result = classify_repair_pathway(pair)
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "classify_repair_pathway",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4ed8687da434749ff2110833c98a49bd59d81d81f13b2294402bfa1059c84aa9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reference_length": 13,
    "edited_length": 9,
    "deleted_bases": 4,
    "inserted_bases": 0,
    "first_difference": 6,
    "microhomology_length": 4,
    "pathway": "MMEJ",
    "reason": "deletion of 4 bp with 4 bp microhomology at the junction (>= 2): microhomology-mediated end joining signature.",
    "rule": "Heuristic indel-signature rule, not a measurement of the pathway used: edited == donor -> HDR; pure deletion whose two junctions share >= ssa_min_homology (25) ..."
  }
}
```

## 5) compare

- **作用：** 比较两条 DNA 序列的累积 GC（或 AT）偏斜曲线及起点估计。

- **处理方法：** 各序列的累积偏斜（同 `calculate`）按相对位置重采样为相同数量的点（`points`），再用 Pearson 相关系数比较两条曲线：+1 表示形状相同，-1 表示偏斜反向（例如互补序列）。同时比较相对起点位置（累积偏斜最小值，启发式）。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 `{"a": 序列, "b": 序列}`，或恰好含两条 DNA 字符串/FASTA 文本的列表。选项：`window`（默认序列长度/20，最大 1000）、`skew_type`（`GC` 默认或 `AT`）、`points`（默认 100）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `skew_type`、`sequences`（每条序列的 `length`、`windows`，起点和终点估计及其相对位置）、`cumulative_skew_correlation`、`origin_relative_position_difference` 以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import compare

seq_a = "GGGC" * 125 + "CCCG" * 125
seq_b = seq_a.translate(str.maketrans("ACGT", "TGCA"))
result = compare({"a": seq_a, "b": seq_b}, options={"window": 50, "points": 20})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 2018,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "426e3a41b1350557a532a1a49d3da3c8d80283e11497373df65cc2cdb90834ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "skew_type": "GC",
    "sequences": {
      "a": {
        "length": 1000,
        "windows": 20,
        "origin_estimate": 975,
        "origin_relative_position": 0.975,
        "terminus_estimate": 475,
        "terminus_relative_position": 0.475
      },
      "b": {
        "length": 1000,
        "windows": 20,
        "origin_estimate": 475,
        "origin_relative_position": 0.475,
        "terminus_estimate": 975,
        "terminus_relative_position": 0.975
      }
    },
    "cumulative_skew_correlation": -1.0,
    "origin_relative_position_difference": 0.5,
    "definition": "Cumulative GC skew of each sequence is resampled to 20 points over relative position and compared with Pearson correlation (+1 = identical shape, -1 = inverted ..."
  }
}
```

## 6) detect

- **作用：** 检测 DNA 序列中与错配修复相关的基序：GATC（Dam）位点以及完美的单至四核苷酸微卫星串联重复。

- **处理方法：** 精确字符串扫描。GATC 是大肠杆菌 MutHLS 用来区分亲本链与新生链的甲基化位点；它是回文序列，因此每个位点只记一个位置。微卫星串联重复为 1 至 4 bp 单元、拷贝数不低于最小值的完美重复，是错配修复缺陷细胞中发生扩增或收缩的位点。扫描不对不稳定性建模；GATC 位点密度与随机期望密度一并给出。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.detect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 DNA 字符串、单条记录的 FASTA 文本、文件路径或 `{"sequence": ...}`。选项：`min_copies`（单元长度 1-6 到最小拷贝数的映射，默认 `{1: 8, 2: 5, 3: 4, 4: 3}`）、`max_results`（默认 1000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `length`、`gatc_sites`（`count`、`density_per_kb`、`expected_random_per_kb`、`positions`、`note`）、`microsatellite_tracts`（`count`、`min_copies`、`tracts`——每项含 `start`、`end`、`unit`、`copies`、`length`——以及 `note`）和 `heuristic` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import detect

sequence = "AAGATCTT" + "A" * 10 + "G" + "CA" * 6 + "GATC" + "TTT"
result = detect(sequence)
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 38,
    "sequence_bases": 38,
    "lengths": [
      38
    ],
    "sha256": "b3adc1d06006d25e6bd4d7958addc957746ee5779f68906991301c079daef47d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 38,
    "gatc_sites": {
      "count": 2,
      "density_per_kb": 52.631579,
      "expected_random_per_kb": 3.9062,
      "positions": [
        2,
        31
      ],
      "note": "GATC is the Dam-methylation site used by E. coli MutHLS to tell the parental from the new strand; it is palindromic, so each site is one position on both strand..."
    },
    "microsatellite_tracts": {
      "count": 2,
      "min_copies": {
        "1": 8,
        "2": 5,
        "3": 4,
        "4": 3
      },
      "tracts": [
        {
          "start": 8,
          "end": 18,
          "unit": "A",
          "copies": 10,
          "length": 10
        },
        {
          "start": 19,
          "end": 31,
          "unit": "CA",
          "copies": 6,
          "length": 12
        }
      ],
      "note": "Perfect repeats of 1-4 bp units are the loci that mismatch-repair-deficient cells expand or contract (microsatellite instability); presence does not imply insta..."
    },
    "heuristic": "Exact string scan for the GATC motif and perfect mono- to tetranucleotide tandem repeats at or above the minimum copy numbers; no statistical model of instabili..."
  }
}
```

## 7) map

- **作用：** 以 S 期与 G1 期读深之比的中位数归一化 log2 值绘制基因组复制时间图。

- **处理方法：** 每个窗口的得分为 log2((S / S 中位数) / (G1 / G1 中位数))（沿用 Koren 等 2012 的做法）；得分越高表示复制越早。`pseudocount` 在求比值前加到每个经中位数归一化的读深上，窗口读深可能为零时必须设置。可对得分按 `smoothing_windows` 个窗口做滑动中位数平滑。高于 `center` 的窗口标为 `early`，低于的标为 `late`（其余为 `intermediate`），高于 `center` 的内部局部极大值窗口作为早复制峰报告。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为映射 `{"s_phase": [各窗口读深], "g1": [各窗口读深]}`，两者基于相同的定长窗口。选项：`pseudocount`（默认 0；窗口可能为空时必须大于 0）、`smoothing_windows`（奇数，默认 1 即不平滑）、`window_size`（bp，默认 1，用于坐标）、`center`（默认 0）、`max_windows_returned`（默认 20000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `windows`、`pseudocount`、`smoothing_windows`、`early_fraction`、`late_fraction`、`early_peaks`（`window`、`start`、`log2_s_over_g1`）、`window_table`（每个窗口的 `start`、`end`、`log2_s_over_g1`、`timing`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import map

depths = {"s_phase": [20, 10, 5, 10, 40, 10], "g1": [10, 10, 10, 10, 10, 10]}
result = map(depths, options={"window_size": 1000})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "map",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 68,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "22da7963a7d6b1d21d72ddc8d686c2e9f526341676ebc8e99db9f63470a79320"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "windows": 6,
    "pseudocount": 0.0,
    "smoothing_windows": 1,
    "early_fraction": 0.333333,
    "late_fraction": 0.166667,
    "early_peaks": [
      {
        "window": 4,
        "start": 4000,
        "log2_s_over_g1": 2.0
      }
    ],
    "window_table": [
      {
        "start": 0,
        "end": 1000,
        "log2_s_over_g1": 1.0,
        "timing": "early"
      },
      {
        "start": 1000,
        "end": 2000,
        "log2_s_over_g1": 0.0,
        "timing": "intermediate"
      },
      {
        "start": 2000,
        "end": 3000,
        "log2_s_over_g1": -1.0,
        "timing": "late"
      },
      "..."
    ],
    "definition": "Replication timing score = log2( (S-phase depth / median S depth) / (G1 depth / median G1 depth) ) per window (Koren et al. 2012 style): higher = replicated ear..."
  }
}
```

## 8) plot

- **作用：** 绘制 DNA 序列的累积 GC（或 AT）偏斜图，标出启发式估计的起点和终点，输出 SVG 字符串。

- **处理方法：** 曲线为按 `calculate` 方式计算的累积偏斜；起点和终点标记为累积偏斜的最小值和最大值（仅适用于类细菌染色体的启发式）。SVG 以文本返回，仅在给出 `output_path` 时写入文件。纯 Python 实现，不依赖绘图库。

- **API：** `dnakit.replication_repair.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 及偏斜选项同 `calculate`（`window`、`step`、`skew_type`）。额外选项：`output_path`（可选的 `.svg` 文件；已存在的文件除非 `overwrite=True` 否则不会被覆盖）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `svg`（SVG 文本）、`output_path`（未写文件时为 `null`）、`skew_type`、`windows`、`origin_estimate`、`terminus_estimate` 以及 `note`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import plot

sequence = "GGGC" * 125 + "CCCG" * 125
result = plot(sequence, options={"window": 100})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 1000,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "3490f072d0775230a9eb8cb9e2bc13b5ea9c4ae2878fa5265873210e8ad3c480"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"760\" height=\"380\" viewBox=\"0 0 760 380\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Cumulative GC skew</title...",
    "output_path": null,
    "skew_type": "GC",
    "windows": 10,
    "origin_estimate": 950,
    "terminus_estimate": 450,
    "note": "Origin/terminus markers are cumulative-skew extrema (heuristic)."
  }
}
```

## 9) predict

- **作用：** 依据 DnaA 盒（TTATNCACA）的密度预测候选细菌复制起点（oriC）区域。

- **处理方法：** 启发式：定位 DnaA 盒（共有序列 TTATNCACA，两条链，汉明距离不超过 `max_mismatches`），返回含盒数最多的互不重叠窗口，并给出 GC 偏斜起点估计（见 `calculate`）及各窗口到它的环状距离，供交叉验证。仅凭盒密度的预测能力较弱，且只适用于细菌。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.replication_repair.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 DNA 字符串、单条记录的 FASTA 文本、文件路径或 `{"sequence": ...}`。选项：`max_mismatches`（默认 0）、`window`（默认 500 bp）、`top_n`（默认 3）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `length`、`dnaa_box_consensus`、`max_mismatches`、`boxes_found`、`candidate_regions`（每个窗口的 `start`、`end`、`dnaa_boxes`、以 `[位置, 链]` 表示的 `box_positions` 以及 `circular_distance_to_gc_skew_origin`）、`gc_skew_origin_estimate` 和 `heuristic` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import predict

filler = "CGGC" * 150
sequence = filler[:600] + "TTATCCACAGG" * 5 + filler[600:]
result = predict(sequence, options={"window": 100, "top_n": 2})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "predict",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 655,
    "sequence_bases": 655,
    "lengths": [
      655
    ],
    "sha256": "d0db177af5706fbac1f08bbbfff162d69116eb8f2169590c43d9b7a5eb7a7e84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 655,
    "dnaa_box_consensus": "TTATNCACA",
    "max_mismatches": 0,
    "boxes_found": 5,
    "candidate_regions": [
      {
        "start": 600,
        "end": 655,
        "dnaa_boxes": 5,
        "box_positions": [
          [
            600,
            "+"
          ],
          [
            611,
            "+"
          ],
          [
            622,
            "+"
          ],
          [
            633,
            "+"
          ],
          "..."
        ],
        "circular_distance_to_gc_skew_origin": 57
      }
    ],
    "heuristic": "oriC is predicted as the non-overlapping windows (default 500 bp) with the most DnaA boxes (TTATNCACA on either strand, Hamming <= max_mismatches). DnaA-box den...",
    "gc_skew_origin_estimate": 570
  }
}
```

## 10) simulate

- **作用：** 模拟从一个或多个复制起点出发的复制叉推进，返回复制时间曲线。

- **处理方法：** 确定性的双向复制叉模型：每个起点在其 `firing_time` 激活，复制叉以 `fork_speed` 向两侧推进。位置 x 的复制时间为各起点 firing_time + |x - 起点| / fork_speed 的最小值（包含被邻近复制叉被动复制的情形；`circular` 时使用环状距离）。终止位点为相向复制叉相遇处。时间在窗口中心处计算，因此 `s_phase_duration` 是窗口中心时间的最大值，窗口较宽时可能小于终止时间。模型不包含随机激活、复制叉停滞或起点许可；随机起点（`n_origins`、`seed`）为均匀分布的位置，均在时间 0 激活。

- **API：** `dnakit.replication_repair.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为基因组长度（bp，整数或 `{"length": n}`）或 DNA 序列（取其长度）。选项：`origins`（位置列表或 `{position, firing_time}`；默认用 `seed` 随机生成 `n_origins` 个位置）、`n_origins`（1）、`seed`（0）、`fork_speed`（每时间单位的 bp 数，默认 1000）、`circular`（默认 `False`）、`window`（默认长度/100）、`time_points`（5）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `length`、`circular`、`fork_speed_bp_per_time`、`origins`、`s_phase_duration`、`mean_replication_time`、`termination_sites`（`position`、`time`）、`window`、`timing_profile`（每个窗口的 `start`、`end`、`time`）、`replication_curve`（`time`、`fraction_replicated`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.replication_repair import simulate

result = simulate(1000, options={"origins": [250, 750], "fork_speed": 10, "window": 100, "time_points": 3})
print(result)
```

- **示例输出：**

```json
{
  "module": "replication_repair",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 4,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "40510175845988f13f6162ed8526f0b09f73384467fa855e1e79b44a56562a58"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 1000,
    "circular": false,
    "fork_speed_bp_per_time": 10.0,
    "origins": [
      {
        "position": 250.0,
        "firing_time": 0.0
      },
      {
        "position": 750.0,
        "firing_time": 0.0
      }
    ],
    "s_phase_duration": 20.0,
    "mean_replication_time": 12.0,
    "termination_sites": [
      {
        "position": 500.0,
        "time": 25.0
      }
    ],
    "window": 100,
    "timing_profile": [
      {
        "start": 0,
        "end": 100,
        "time": 20.0
      },
      {
        "start": 100,
        "end": 200,
        "time": 10.0
      },
      {
        "start": 200,
        "end": 300,
        "time": 0.0
      },
      "..."
    ],
    "replication_curve": [
      {
        "time": 0.0,
        "fraction_replicated": 0.2
      },
      {
        "time": 6.666667,
        "fraction_replicated": 0.2
      },
      {
        "time": 13.333333,
        "fraction_replicated": 0.6
      },
      "..."
    ],
    "definition": "Deterministic bidirectional-fork model: all origins fire at their firing_time and forks move at fork_speed in both directions; the replication time of position ..."
  }
}
```
