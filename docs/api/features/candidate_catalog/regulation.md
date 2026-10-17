# 基因调控

本类收录 **11** 项有文档的候选功能，其中 10 项在本地原生计算出真实结果（PWM 基序扫描与打分、基序比较与校验、富集检验、变异效应代理、计算机模拟突变、基于位置的区域分类以及基于距离的增强子-基因连接），1 项（`predict`）因需要已训练的序列模型权重而始终返回 `status: "not_implemented"`。所有打分都明确是启发式：基序得分是 PWM 对数几率（log2(p / 背景) 之和），没有命中 p 值，也没有训练得到的模型；`link`/`annotate(entity="target_gene")` 只按到 TSS 的距离给区域分配基因；`classify` 只使用位置规则。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) annotate

- **作用：** 用 PWM 基序命中注释序列，或为区域标注最近的靶基因（`options["entity"]`）。

- **处理方法：** `motif`（默认；别名 `tf_binding_site`）：用每个 PWM 给每条序列的每个窗口打分。PWM 打分采用标准的对数几率模型：计数或概率矩阵（各列之和在 1e-3 内等于 1 时视为概率，否则视为计数；可用 `matrix_type` 强制指定）先转换为概率 p（计数：(计数 + `pseudocount`) / (列总数 + 4 `pseudocount`)，`pseudocount` 默认 0；概率下限取 1e-6，所以在 0.25 的背景下零值项约扣 17.9 bit），窗口得分为各位置 log2(p(碱基) / 背景(碱基)) 之和，`background` 默认每个碱基 0.25。相对得分为矩阵的 (得分 - 最小值) / (最大值 - 最小值)。两条链都会扫描（反链对反向互补后的窗口打分）；含 N 的窗口被跳过。这只是一个简单的亲和力代理：没有命中 p 值，没有全基因组背景模型，也没有训练得到的模型。窗口得分 >= `threshold`（绝对 bit 值）时为命中；未给 `threshold` 时，相对得分 >= `min_relative_score`（默认 0.8）为命中。命中按位置排序，每条序列至多保留 `max_hits` 个（`truncated` 表示是否触顶）。`target_gene`：等同于 `method="nearest"` 的 `link`。距离为 |TSS - 中点|，中点 = 区域的 (start + end) // 2；TSS 取显式给出的 `tss`，否则对 + 链或无链方向的基因取 `start`，对 - 链基因取 `end - 1`。`method="nearest"` 把每个区域连到同一染色体上最近的单个 TSS（距离并列时取名称字母序在前的基因），前提是在 `max_distance` 之内；`method="window"` 连接 TSS 在 `max_distance` 之内的所有基因。`direction` 是区域相对基因方向的位置（`upstream`、`downstream` 或 `at_tss`）。这是邻近性启发式，不是调控证据：不使用接触、表达或 TAD 信息，调用会附带相应警告。

- **API：** `dnakit.regulation.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `motif`：`input` 为DNA（字符串、FASTA 文本/路径、列表或 `{id: 序列}`；ACGT 与 N，U 视为 T）；`options`：`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）；`threshold`（绝对 bit 值）或 `min_relative_score`（默认 0.8）；`strand`（`"both"` 默认 | `"forward"` | `"reverse"`）；`max_hits`（默认 10000）。`target_gene`：`input` 为区域：BED 文本/路径、行或映射（0 起始半开）；`options`：`genes`（必填：带 `chrom`、`start`、`end`、`name` 以及可选 `strand` 或显式 `tss` 的行或映射）、`max_distance`（默认 1,000,000 bp）。`options["entity"]` 取 `"motif"`（默认）、`"tf_binding_site"` 或 `"target_gene"`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`motif` 时 `result` 含 `entity` 和 `sequences`（每条序列有 `id`、`length`、`hit_count`、`hits`、`truncated`；每个命中有 `motif`、正链 0 起始半开区间 `start`/`end`、`strand`、`score`（bit）、`relative_score` 和 `matched`，即按 motif 方向给出的匹配序列）。`target_gene` 时结果同 `link`：`result` 含 `method`、`max_distance` 和 `links`（每个区域一项：`region`、`midpoint` 和 `genes`，其中每个基因有 `gene`、`tss`、`distance`、`direction`、`strand`；范围内没有基因时列表为空）。

- **示例代码：**

```python
import json
from dnakit.regulation import annotate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
result = annotate("TTGATAAGGCTTATCAC",
                  options={"motifs": {"GATA": gata}, "pseudocount": 0.5, "min_relative_score": 0.9})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 17,
    "sequence_bases": 17,
    "lengths": [
      17
    ],
    "sha256": "ef32207549d8a6455d7b256ce50fb87992ed624845b1ae6279f906be3d504032"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "motif",
    "sequences": [
      {
        "id": "seq1",
        "length": 17,
        "hit_count": 2,
        "hits": [
          {
            "motif": "GATA",
            "start": 2,
            "end": 7,
            "strand": "+",
            "score": 8.897784945974175,
            "relative_score": 1.0,
            "matched": "GATAA"
          },
          {
            "motif": "GATA",
            "start": 10,
            "end": 15,
            "strand": "-",
            "score": 8.897784945974175,
            "relative_score": 1.0,
            "matched": "GATAA"
          }
        ],
        "truncated": false
      }
    ]
  }
}
```

## 2) calculate

- **作用：** 计算基于 PWM 的变异调控效应差值，或计算 motif 的信息含量（`options["entity"]`）。

- **处理方法：** `regulatory_effect_delta`（默认）：对每个 motif，分别计算参考序列和替代序列在所有窗口和两条链上的最佳对数几率得分，`delta` = alt_score - ref_score（bit，负值表示预测结合变弱）；`max_abs_delta` 是各 motif 上 |delta| 的最大值。PWM 打分采用标准的对数几率模型：计数或概率矩阵（各列之和在 1e-3 内等于 1 时视为概率，否则视为计数；可用 `matrix_type` 强制指定）先转换为概率 p（计数：(计数 + `pseudocount`) / (列总数 + 4 `pseudocount`)，`pseudocount` 默认 0；概率下限取 1e-6，所以在 0.25 的背景下零值项约扣 17.9 bit），窗口得分为各位置 log2(p(碱基) / 背景(碱基)) 之和，`background` 默认每个碱基 0.25。相对得分为矩阵的 (得分 - 最小值) / (最大值 - 最小值)。两条链都会扫描（反链对反向互补后的窗口打分）；含 N 的窗口被跳过。这只是一个简单的亲和力代理：没有命中 p 值，没有全基因组背景模型，也没有训练得到的模型。它是位置权重矩阵代理，不是训练得到的模型。`pwm_information_content`：逐位置信息含量 sum(p log2(p / 背景)) 及其总和（bit；背景均匀时，完全保守的一列为 2 bit）。

- **API：** `dnakit.regulation.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `regulatory_effect_delta`：`input` 为变异记录或列表，形式为 `{"id", "ref_sequence", "alt_sequence"}` 或 `{"id", "sequence", "position", "ref", "alt"}`（`position` 0 起始；`ref` 必须与序列一致，否则报错），至多 10000 个；`options`：`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）。`pwm_information_content`：`input` 为一个或多个 PWM（上述 `motifs` 形式）；`options`：`background`、`pseudocount`、`matrix_type`。`options["entity"]` 取 `"regulatory_effect_delta"`（默认）或 `"pwm_information_content"`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`regulatory_effect_delta`：`result` 含 `entity` 和 `variants`（每个变异有 `id`、`motifs`（含 `motif`、`ref_score`、`alt_score`、`delta`）和 `max_abs_delta`；没有可打分窗口的 motif 被省略），`warnings` 说明 PWM 对数几率只是简单代理。`pwm_information_content`：`result` 含 `entity` 和 `motifs`（`motif`、`length`、`per_position`、`total_bits`）。

- **示例代码：**

```python
import json
from dnakit.regulation import calculate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
variant = {"id": "v1", "sequence": "TTGATAAGG", "position": 4, "ref": "T", "alt": "C"}
result = calculate(variant, options={"entity": "regulatory_effect_delta", "motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 9,
    "lengths": [
      9
    ],
    "sha256": "d79c2457c84e2ff160d36cf926e7da44d70d976b7e8eee511dbfc530fac2bd4f"
  },
  "warnings": [
    "PWM log-odds is a simple binding-affinity proxy."
  ],
  "status": "ok",
  "result": {
    "entity": "regulatory_effect_delta",
    "variants": [
      {
        "id": "v1",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 3.6123827271119264,
            "delta": -5.285402218862249
          }
        ],
        "max_abs_delta": 5.285402218862249
      }
    ]
  }
}
```

## 3) classify

- **作用：** 按区域相对基因 TSS 的位置，把区域分为启动子、近端或远端（仅用位置规则）。

- **处理方法：** TSS 取显式给出的 `tss`，否则对 + 链或无链方向的基因取 `start`，对 - 链基因取 `end - 1`。启动子窗口对 + 链基因为 [TSS - `promoter_upstream`, TSS + `promoter_downstream`)，对 - 链基因为 [TSS - `promoter_downstream` + 1, TSS + `promoter_upstream` + 1)。与任一启动子窗口重叠的区域为 `promoter`；否则其中点 (start + end) // 2 与同一染色体上某个 TSS 的距离在 `proximal_distance` 以内为 `proximal`，其余为 `distal`（染色体上没有基因时也是 `distal`）。最近基因按中点到 TSS 的距离给出。阈值只是常规默认值，不使用染色质、保守性或表达证据（调用会附带警告）。

- **API：** `dnakit.regulation.classify(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：区域：BED 文本/路径、行或映射（0 起始半开）。`options`：`genes`（必填：带 `chrom`、`start`、`end`、`name` 以及可选 `strand` 或显式 `tss` 的行或映射）；`promoter_upstream`（默认 1000）、`promoter_downstream`（默认 100）、`proximal_distance`（默认 10000）；`entity`（只支持 `"promoter"`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.regions` 每个区域一项：`chrom`、`start`、`end`、`class`（`promoter`、`proximal` 或 `distal`）、`nearest_gene` 和 `distance_to_tss`（染色体上没有基因时两者为 `null`）。

- **示例代码：**

```python
import json
from dnakit.regulation import classify

genes = [{"name": "GENE_A", "chrom": "chr1", "start": 5000, "end": 9000, "strand": "+"}]
regions = "chr1\t4600\t4800\tpromoter_like\nchr1\t9500\t9700\tnear\nchr1\t80000\t80400\tfar\n"
result = classify(regions, options={"genes": genes})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "classify",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ae3771f589ac85420ee7b88632e87b7bebaffb31e60bde9d3553caf7ca010ee8"
  },
  "warnings": [
    "Classification uses genomic position only."
  ],
  "status": "ok",
  "result": {
    "regions": [
      {
        "chrom": "chr1",
        "start": 4600,
        "end": 4800,
        "class": "promoter",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 300
      },
      {
        "chrom": "chr1",
        "start": 9500,
        "end": 9700,
        "class": "proximal",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 4600
      },
      {
        "chrom": "chr1",
        "start": 80000,
        "end": 80400,
        "class": "distal",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 75200
      }
    ]
  }
}
```

## 4) compare

- **作用：** 通过最佳无空位比对（逐列 Pearson 相关系数的平均值）两两比较 PWM 基序。

- **处理方法：** 先把 motif 转成概率矩阵。对每一对 a、b，尝试 b 的所有偏移以及两种方向（正向和反向互补），要求至少 `min_overlap` 列对齐；一个比对的相似度是对齐各列 A/C/G/T 概率向量 Pearson 相关系数的平均值（方差为 0 的列按 0 计）。最佳比对给出 `similarity`（1.0 表示相同）、`offset`（b 的第一列相对 a 第一列的位置）、`orientation` 和 `overlap`。这只是简单的逐列相关得分，不是 Tomtom 的 p 值/E 值，也不给出显著性。

- **API：** `dnakit.regulation.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：两个或更多 motif（矩阵列表或 `{名称: 矩阵}`；矩阵为计数或概率，形式同 `annotate`）。`options`：`entity`（只支持 `"motif"`）；`min_overlap`（默认 3，不超过最短 motif 的长度）；`background`、`pseudocount`、`matrix_type` 同 `annotate`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.pairs` 每个无序 motif 对一项：`a`、`b`、`similarity`、`offset`、`orientation`（`forward` 或 `reverse_complement`）和 `overlap`。

- **示例代码：**

```python
import json
from dnakit.regulation import compare

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
ttatc = [[0, 0, 0, 10], [0, 0, 0, 10], [10, 0, 0, 0], [0, 0, 0, 10], [0, 10, 0, 0]]   # rows A, C, G, T per position
result = compare({"GATA": gata, "TTATC": ttatc}, options={"min_overlap": 4})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 189,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8983b08a91514321197f33dac0b3688b1e0ef8afdc615b50ecc871d94d24680e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": [
      {
        "a": "GATA",
        "b": "TTATC",
        "similarity": 0.9989409928538651,
        "offset": 0,
        "orientation": "reverse_complement",
        "overlap": 5
      }
    ]
  }
}
```

## 5) enrich

- **作用：** 检验 motif 在前景序列与背景序列中的富集（超几何或二项上尾检验）。

- **处理方法：** 检验单位是至少含一个 motif 命中的序列条数。`hypergeometric`（默认）：令 N = fg_total + bg_total、K = fg_hits + bg_hits、n = fg_total，p = P(X >= fg_hits)，X ~ 超几何(N, K, n)。`binomial`：p = P(X >= fg_hits)，X ~ 二项(fg_total, bg_hits / bg_total)。`fold_enrichment` = (fg_hits / fg_total) / (bg_hits / bg_total)（背景率为 0 时为 `null`），`q_value` 是对被检验 motif 的 Benjamini-Hochberg 校正。给定序列时，命中由 `annotate` 的 PWM 扫描（`threshold` 或 `min_relative_score`）判定，因此结果取决于该阈值；背景就是所提供的序列（不做 GC 或长度匹配）。

- **API：** `dnakit.regulation.enrich(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射，可以是计数 `{"foreground_hits", "foreground_total", "background_hits", "background_total"}`（非负整数；命中数不能超过总数，总数必须为正），也可以是序列 `{"foreground": [...], "background": [...]}`。`options`：`entity`（`"motif"` 默认 | `"tf_binding_site"`）；`test`（`"hypergeometric"` 默认 | `"binomial"`）；`name`（计数形式的标签，默认 `motif`）；给定序列时还需要`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）、`threshold`/`min_relative_score`、`strand`；每个 motif 各自检验。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `test` 和 `motifs`（每个 motif 有 `motif`、`foreground_hits`、`foreground_total`、`background_hits`、`background_total`、`fold_enrichment`、`p_value`、`q_value`）。

- **示例代码：**

```python
import json
from dnakit.regulation import enrich

counts = {"foreground_hits": 30, "foreground_total": 100, "background_hits": 10, "background_total": 100}
result = enrich(counts, options={"test": "hypergeometric", "name": "GATA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "enrich",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 96,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "74bb3c6e8fb50effdc77642b2707471638887481358940928a82775ea5670cf4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "hypergeometric",
    "motifs": [
      {
        "motif": "GATA",
        "foreground_hits": 30,
        "foreground_total": 100,
        "background_hits": 10,
        "background_total": 100,
        "fold_enrichment": 2.9999999999999996,
        "p_value": 0.000325205353801686,
        "q_value": 0.000325205353801686
      }
    ]
  }
}
```

## 6) link

- **作用：** 按区域中点到基因 TSS 的距离，把调控区域（增强子等）分配给基因。

- **处理方法：** 基于距离的分配。距离为 |TSS - 中点|，中点 = 区域的 (start + end) // 2；TSS 取显式给出的 `tss`，否则对 + 链或无链方向的基因取 `start`，对 - 链基因取 `end - 1`。`method="nearest"` 把每个区域连到同一染色体上最近的单个 TSS（距离并列时取名称字母序在前的基因），前提是在 `max_distance` 之内；`method="window"` 连接 TSS 在 `max_distance` 之内的所有基因。`direction` 是区域相对基因方向的位置（`upstream`、`downstream` 或 `at_tss`）。这是邻近性启发式，不是调控证据：不使用接触、表达或 TAD 信息，调用会附带相应警告。`entity` 的取值（`enhancer`、`promoter`、`silencer`、`insulator`、`target_gene`）只是标签，对所有取值规则相同。

- **API：** `dnakit.regulation.link(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：区域：BED 文本/路径、行或映射（0 起始半开）。`options`：`entity`（`"enhancer"` 默认 | `"promoter"` | `"silencer"` | `"insulator"` | `"target_gene"`）；`genes`（必填：带 `chrom`、`start`、`end`、`name` 以及可选 `strand` 或显式 `tss` 的行或映射）；`method`（`"nearest"` 默认 | `"window"`）；`max_distance`（默认 1,000,000 bp）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `method`、`max_distance` 和 `links`（每个区域一项：`region`、`midpoint` 和 `genes`，其中每个基因有 `gene`、`tss`、`distance`、`direction`、`strand`；范围内没有基因时列表为空）。

- **示例代码：**

```python
import json
from dnakit.regulation import link

enhancers = "chr1\t10000\t10400\tenh1\nchr1\t90000\t90200\tenh2\n"
genes = [
    {"name": "GENE_A", "chrom": "chr1", "start": 12000, "end": 20000, "strand": "+"},
    {"name": "GENE_B", "chrom": "chr1", "start": 60000, "end": 91000, "strand": "-"},
]
result = link(enhancers, options={"entity": "enhancer", "genes": genes, "method": "nearest", "max_distance": 50000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "link",
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

## 7) perform_in_silico_mutagenesis

- **作用：** 用 PWM 对序列做计算机模拟的饱和突变（每一种单碱基替换）。

- **处理方法：** 对每个位置 i 和每个替代碱基 b，效应 = best(突变后) - best(参考)，其中 best 是所有窗口和两条链上的最大 PWM 对数几率得分（打分方式同 `score`：各位置 log2(p(碱基) / 背景(碱基)) 之和，背景默认 0.25）。模型就是 PWM 本身，不是训练得到的神经网络，替换既没有降低得分最高的窗口、也没有产生更好的窗口时，效应为 0。序列短于 motif 时报错。

- **API：** `dnakit.regulation.perform_in_silico_mutagenesis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：一条 DNA 序列（字符串、FASTA 文本/路径或 `{"sequence": ...}`），至多 2000 个碱基，仅限 ACGT。`options`：`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `sequence`（id）和 `motifs`；每个 motif 有 `motif`、`reference_score` 和 `positions`，后者每个位置一项：`position`、`ref` 和 `effects`（各替代碱基到 bit 变化量的映射）。`warnings` 说明模型是 PWM。

- **示例代码：**

```python
import json
from dnakit.regulation import perform_in_silico_mutagenesis

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
result = perform_in_silico_mutagenesis("AGATAAC", options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "regulation",
  "operation": "perform_in_silico_mutagenesis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 7,
    "sequence_bases": 7,
    "lengths": [
      7
    ],
    "sha256": "06ce638d4a9efc08740e01618439ea14f498bb9df12324428d7885607a97afcd"
  },
  "warnings": [
    "The scoring model is a PWM, not a trained sequence model."
  ],
  "status": "ok",
  "result": {
    "sequence": "seq1",
    "motifs": [
      {
        "motif": "GATA",
        "reference_score": 8.897784945974175,
        "positions": [
          {
            "position": 0,
            "ref": "A",
            "effects": {
              "C": 0.0,
              "G": 0.0,
              "T": 0.0
            }
          },
          {
            "position": 1,
            "ref": "G",
            "effects": {
              "A": -3.5443205162238103,
              "C": -3.5443205162238103,
              "T": -3.5443205162238103
            }
          },
          {
            "position": 2,
            "ref": "A",
            "effects": {
              "C": -5.285402218862249,
              "G": -5.285402218862249,
              "T": -3.7004397181410926
            }
          },
          {
            "position": 3,
            "ref": "T",
            "effects": {
              "A": -3.7004397181410926,
              "C": -5.285402218862249,
              "G": -5.285402218862249
            }
          },
          "..."
        ]
      }
    ]
  }
}
```

## 8) predict

- **作用：** 由序列预测调控活性。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。由序列预测功能需要 Enformer、Borzoi、DeepSEA、Sei 之类的已训练模型及其权重，DNAKit 并未内置。需要基于 PWM 的代理时请用 `score`、`calculate` 或 `perform_in_silico_mutagenesis`。不会编造结果。

- **API：** `dnakit.regulation.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a trained regulatory sequence model and weights (Enformer, Borzoi, DeepSEA, Sei)`，`message` 为 `Regulatory-activity prediction needs trained model weights that are not bundled.`

- **示例代码：**

```python
import json
from dnakit.regulation import predict

result = predict("ACGTACGTACGTACGT")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
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
  "required_dependency": "a trained regulatory sequence model and weights (Enformer, Borzoi, DeepSEA, Sei)",
  "message": "Regulatory-activity prediction needs trained model weights that are not bundled."
}
```

## 9) prioritize

- **作用：** 按各 motif 上 PWM 得分变化的最大绝对值对变异排序（调控影响的代理指标）。

- **处理方法：** 每个变异的打分方式同 `calculate`（`regulatory_effect_delta`）：参考序列和替代序列在所有窗口与两条链上的最佳对数几率得分，`delta` = alt_score - ref_score（bit）。变异按 `max_abs_delta` 降序排列，并列按 id；没有可打分 motif 窗口的变异排在最后。排序只依据这个 PWM 代理指标，不使用保守性、染色质状态或训练模型（调用会附带警告）。

- **API：** `dnakit.regulation.prioritize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：变异记录或列表，形式为 `{"id", "ref_sequence", "alt_sequence"}` 或 `{"id", "sequence", "position", "ref", "alt"}`（`position` 0 起始；`ref` 必须与序列一致，否则报错），至多 10000 个。`options`：`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）；`entity`（只支持 `"regulatory_variant"`）；`top_n`（可选）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.ranked` 含 `calculate` 的变异条目（`id`、`motifs`、`max_abs_delta`）以及从 1 开始的 `rank`。

- **示例代码：**

```python
import json
from dnakit.regulation import prioritize

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
variants = [
    {"id": "in_motif", "sequence": "TTGATAAGG", "position": 4, "ref": "T", "alt": "C"},
    {"id": "outside", "sequence": "TTGATAAGG", "position": 8, "ref": "G", "alt": "A"},
]
result = prioritize(variants, options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "prioritize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 167,
    "sequence_bases": 18,
    "lengths": [
      9,
      9
    ],
    "sha256": "e31572e30c31d4dff3508cf3b45e845f40532b2bb237769b2d2ac7f3b605f6f8"
  },
  "warnings": [
    "Ranking uses a PWM score change proxy only."
  ],
  "status": "ok",
  "result": {
    "ranked": [
      {
        "id": "in_motif",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 3.6123827271119264,
            "delta": -5.285402218862249
          }
        ],
        "max_abs_delta": 5.285402218862249,
        "rank": 1
      },
      {
        "id": "outside",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 8.897784945974175,
            "delta": 0.0
          }
        ],
        "max_abs_delta": 0.0,
        "rank": 2
      }
    ]
  }
}
```

## 10) score

- **作用：** 用 PWM 基序给序列打分：每条序列、每个 motif 的最佳对数几率窗口。

- **处理方法：** 对允许的链上的每个窗口打分，报告得分最高的窗口及其位置。PWM 打分采用标准的对数几率模型：计数或概率矩阵（各列之和在 1e-3 内等于 1 时视为概率，否则视为计数；可用 `matrix_type` 强制指定）先转换为概率 p（计数：(计数 + `pseudocount`) / (列总数 + 4 `pseudocount`)，`pseudocount` 默认 0；概率下限取 1e-6，所以在 0.25 的背景下零值项约扣 17.9 bit），窗口得分为各位置 log2(p(碱基) / 背景(碱基)) 之和，`background` 默认每个碱基 0.25。相对得分为矩阵的 (得分 - 最小值) / (最大值 - 最小值)。两条链都会扫描（反链对反向互补后的窗口打分）；含 N 的窗口被跳过。这只是一个简单的亲和力代理：没有命中 p 值，没有全基因组背景模型，也没有训练得到的模型。最佳得分是窗口得分的最大值，会随序列变长而增大，它不是概率，也不代表显著性。

- **API：** `dnakit.regulation.score(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：DNA（字符串、FASTA 文本/路径、列表或 `{id: 序列}`；ACGT 与 N，U 视为 T）。`options`：`options["motifs"]`（必填；别名 `motif`）：一个或多个 PWM（列表或 `{名称: pwm}`），可以是计数或概率，形式为 `{"A": [...], "C": [...], "G": [...], "T": [...]}`、行 `[[a, c, g, t], ...]` 或 JASPAR 风格文本（至多 1000 个位置、500 个 motif）；`background`（A/C/G/T 映射，默认各 0.25）、`pseudocount`（默认 0）、`matrix_type`（`"auto"` 默认 | `"counts"` | `"probabilities"`）；`entity`（`"motif"` 默认 | `"tf_binding_site"`）；`strand`（`"both"` 默认 | `"forward"` | `"reverse"`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.sequences` 每条序列一项：`id` 和 `motifs`，其中每个 motif 有 `motif` 和 `best`（`start`、`end` 为正链 0 起始半开坐标，以及 `strand`、`score`（bit）、`relative_score`）；序列短于 motif 或只有含 N 的窗口时 `best` 为 `null`。

- **示例代码：**

```python
import json
from dnakit.regulation import score

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
sequences = {"s1": "CCCGATAAGGCC", "s2": "CCCCCCCCCCCC"}
result = score(sequences, options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "db0e3c9cf6a6c083e8cbd0b701de8166870199076e7dbe5f73b61e223ee0f361"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequences": [
      {
        "id": "s1",
        "motifs": [
          {
            "motif": "GATA",
            "best": {
              "start": 3,
              "end": 8,
              "strand": "+",
              "score": 8.897784945974175,
              "relative_score": 1.0
            }
          }
        ]
      },
      {
        "id": "s2",
        "motifs": [
          {
            "motif": "GATA",
            "best": {
              "start": 0,
              "end": 5,
              "strand": "-",
              "score": -12.16787507624152,
              "relative_score": 0.1440196391333088
            }
          }
        ]
      }
    ]
  }
}
```

## 11) validate

- **作用：** 校验 PWM 基序：形状、有限的非负数值、非零列，以及各列之和是否为 1。

- **处理方法：** 逐个解析矩阵（A/C/G/T 列表、Lx4 的行列表或 JASPAR 风格文本）。无法解析、含负数或非有限数值、或某列之和为零时矩阵无效（错误信息放入 `errors`）。各列之和在 `tolerance` 内等于 1 时 `matrix_type` 为 `probabilities`，否则为 `counts`。对有效矩阵，`total_information_bits` 为 PWM 信息含量，即各位置 sum_b p log2(p / 0.25) 之和，概率下限取 1e-6。这只是格式与合理性检查，不判断 motif 是否有生物学意义。

- **API：** `dnakit.regulation.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：一个或多个 PWM（列表或 `{名称: 矩阵}`），计数或概率，形式同 `annotate`。`options`：`entity`（只支持 `"motif"`）；`tolerance`（默认 1e-3，用于列和检验）；`background`、`pseudocount` 同 `annotate`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"regulation"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `motifs`（每个 motif 有 `motif`、`valid`、`errors`、`length`，有效矩阵还有 `matrix_type` 和 `total_information_bits`）和 `all_valid`。格式错误的矩阵以 `valid: false` 报告，而不是抛出错误。

- **示例代码：**

```python
import json
from dnakit.regulation import validate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
bad = [[1, 0, 0, 0], [0, -1, 0, 0]]
result = validate({"GATA": gata, "bad": bad})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "regulation",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 141,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a8ead67ae3934d9b11ccd1049efcc1a88aac7adb7a3ccf0433059478c62c867e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "motifs": [
      {
        "motif": "GATA",
        "valid": true,
        "errors": [],
        "length": 5,
        "matrix_type": "counts",
        "total_information_bits": 7.7240815230348225
      },
      {
        "motif": "bad",
        "valid": false,
        "errors": [
          "bad: negative or non-finite entries at positions [1]."
        ],
        "length": 2
      }
    ],
    "all_valid": false
  }
}
```
