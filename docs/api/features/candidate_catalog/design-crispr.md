# CRISPR、引物与 DNA 设计

本类收录 **19** 项有文档的候选功能，其中 18 项在本地原生计算出真实结果（使用标准库，部分调用 `dnakit.patterns` 和 `dnakit.molbio`），1 项（`predict_editing_outcomes`）因需要训练好的修复结果模型权重而始终返回 `status: "not_implemented"`。若干步骤明确是启发式而非经过验证的预测器：向导 `proxy_score` 是基于经验规则的代理分数，没有用编辑数据拟合；脱靶搜索是不带 CFD/MIT 评分的无空位 Hamming 搜索；引物检查只统计连续的 Watson-Crick 配对，没有自由能模型。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) assemble_dna_fragments

- **作用：** 按给定顺序，依据精确的末端重叠把 DNA 片段连接起来（计算机模拟的 Gibson/重叠拼接）。

- **处理方法：** 调用 `dnakit.molbio.simulate_assembly(method="gibson")`：相邻片段在长度落入 `[min_overlap, max_overlap]` 的最长精确后缀/前缀重叠处连接；`circularize` 为真时，最后一个片段还会与第一个片段首尾相接。只使用精确的序列重叠，不模拟反应动力学、二级结构和错配退火。相邻片段没有符合条件的重叠时，函数抛出 DNAKit 错误而不是返回产物。

- **API：** `dnakit.design_crispr.assemble_dna_fragments(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：按顺序排列的片段，可以是 DNA 字符串列表、映射 `{id: 序列}`（插入顺序即拼接顺序）或 FASTA 文本/路径；至少 2 个、至多 1000 个片段，字符限于 ACGT（U 视为 T）。`options`：`min_overlap`（默认 15）、`max_overlap`（默认 200）、`circularize`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `product`（拼接产物序列）、`length`、`circular` 和 `junctions`（每个连接处一项，含 `left`、`right`、`overlap`、`overlap_length`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import assemble_dna_fragments

overlap = "ACGTACGTACGTACGTAAAA"
fragments = {"left": "AAAAAC" + overlap, "right": overlap + "CCCCCT"}
result = assemble_dna_fragments(fragments, options={"min_overlap": 15})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "assemble_dna_fragments",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 77,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c9c510fcdf0c197ea95010181d5c1150e5924cb4cb3ea553a5a3e7054f95e81f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "product": "AAAAACACGTACGTACGTACGTAAAACCCCCT",
    "length": 32,
    "circular": false,
    "junctions": [
      {
        "left": "left",
        "right": "right",
        "overlap": "ACGTACGTACGTACGTAAAA",
        "overlap_length": 20
      }
    ]
  }
}
```

## 2) calculate_primer_tm

- **作用：** 计算一条或多条引物的解链温度（Tm）：短寡核苷酸用 Wallace 规则，较长的用 SantaLucia 最近邻模型，同时给出每条引物的长度和 GC 含量比例。

- **处理方法：** 每条引物先转大写并去掉空白；GC 比例 = (G + C) / 长度。`method` 为 `"wallace"`：Tm = 2 x (A + T) + 4 x (G + C)（摄氏度），是适用于约 14 nt 及以下寡核苷酸的经验规则，忽略盐浓度、寡核苷酸浓度和序列上下文，对较长的寡核苷酸会高估 Tm（长于 14 nt 的引物使用 Wallace 时会在 `warnings` 中提示）。`"nearest_neighbor"`：调用 `dnakit.thermodynamics.melting_temperature`（SantaLucia 1998 统一 DNA/DNA 参数，引物与其互补链完全匹配的双链，2-60 nt），使用给定的 `sodium_molar`、`potassium_molar` 和 `strand_concentration_molar`，结果保留 0.01 度；不含 Mg2+、dNTP、DMSO 校正，因此在含 Mg2+ 的 PCR 缓冲液中实际值通常更高。`"auto"`（默认）：长度 <= 13 用 Wallace，更长用最近邻。简并碱基（N、IUPAC 码）会被拒绝，因为 Tm 对它们没有定义；`U` 默认被拒绝，`rna` 为 true 时按 `T` 处理。

- **API：** `dnakit.design_crispr.calculate_primer_tm(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：字符串形式的引物（5'->3'）、`{"id": ..., "sequence": ...}`、映射 `{名称: 引物}`，或字符串 / `{"id", "sequence"}` 映射的列表；1 到 1000 条，每条至多 200 个碱基（含多行的字符串视为一条引物，空白会被去掉）。`options`：`method`（`"auto"` 默认 | `"wallace"` | `"nearest_neighbor"`）、`sodium_molar`（默认 0.05）、`potassium_molar`（默认 0.0）、`strand_concentration_molar`（默认 250e-9），单位均为 mol/L，仅用于最近邻 Tm；以及 `rna`（默认 `false`）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`（`"calculate_primer_tm"`）、`provenance`（`implementation`、`backend`）、`input_summary`（通用文本摘要：`record_count`、`text_length`、`sequence_bases`、`lengths`、`sha256`，按输入的文本形式计算）、`warnings` 和 `result`。`result.primers` 按输入顺序每条引物一项：`id`（未命名的为 `primer1`、`primer2`……）、`sequence`（清理后）、`length`、`gc_fraction`、`tm_celsius`（浮点数）和 `method`（`"wallace"` 或 `"nearest_neighbor"`，实际使用的方法）。只要有引物用了最近邻模型，就有 `result.conditions`（`sodium_molar`、`potassium_molar`、`strand_concentration_molar`）。序列为空或含简并碱基时抛出 `InputFormatError`（`CRISPR_SEQUENCE_INVALID`）；对 2-60 nt 之外的引物要求最近邻时抛出 `ConfigurationError`（`CRISPR_PRIMER_LENGTH`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import calculate_primer_tm

result = calculate_primer_tm({"short": "ATGCGTAC", "long": "GACGTTAGCATCAGGCTAAC"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "calculate_primer_tm",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d1cef6d866bd036122f644676770ae56c33d025666a8977ff9212d2c0653bfd5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "id": "short",
        "sequence": "ATGCGTAC",
        "length": 8,
        "gc_fraction": 0.5,
        "tm_celsius": 24.0,
        "method": "wallace"
      },
      {
        "id": "long",
        "sequence": "GACGTTAGCATCAGGCTAAC",
        "length": 20,
        "gc_fraction": 0.5,
        "tm_celsius": 52.86,
        "method": "nearest_neighbor"
      }
    ],
    "conditions": {
      "sodium_molar": 0.05,
      "potassium_molar": 0.0,
      "strand_concentration_molar": 2.5e-07
    }
  }
}
```

## 3) check_primer_dimer

- **作用：** 用最长的完全互补、无空位片段筛查引物的自身二聚体和引物间二聚体。

- **处理方法：** 每条引物与自身及其他每条引物比较。对引物 A、B，`longest_complementary_run` 是 A 与 B 的反向互补序列之间最长公共子串的长度（动态规划，O(len(A)*len(B))）；`three_prime_run_a`/`three_prime_run_b` 是一条引物 3' 端最长的、与另一条引物某部分反向互补的后缀。任一 3' 端互补长度超过 `max_3prime_run` 时标记 `three_prime_dimer`，最长互补片段超过 `max_run` 时标记 `long_complementary_run`。只统计连续的 Watson-Crick 配对，不是自由能（最近邻）计算，阈值为经验值。

- **API：** `dnakit.design_crispr.check_primer_dimer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：单条引物字符串（只检查自身二聚体）、引物列表（至多 100 条）或映射如 `{"forward": ..., "reverse": ...}`；引物仅限 ACGT，至多 200 nt。`options`：`max_3prime_run`（默认 3）、`max_run`（默认 6）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.pairs` 对每个无序引物对（含自身配对）给出一项：`a`、`b`、`self`、`longest_complementary_run`、`run_start_in_a`（0 起始，无则为 `null`）、`three_prime_run_a`、`three_prime_run_b`、`three_prime_dimer`、`long_complementary_run`。

- **示例代码：**

```python
import json
from dnakit.design_crispr import check_primer_dimer

result = check_primer_dimer({"forward": "AAAAGGGG", "reverse": "CCCCTTTT"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_dimer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 46,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5d740a084f018177281992c5b4d9cb8ee12cee423c0c05d77fe37c9bf5920439"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": [
      {
        "a": "forward",
        "b": "forward",
        "self": true,
        "longest_complementary_run": 0,
        "run_start_in_a": null,
        "three_prime_run_a": 0,
        "three_prime_run_b": 0,
        "three_prime_dimer": false,
        "long_complementary_run": false
      },
      {
        "a": "forward",
        "b": "reverse",
        "self": false,
        "longest_complementary_run": 8,
        "run_start_in_a": 0,
        "three_prime_run_a": 8,
        "three_prime_run_b": 8,
        "three_prime_dimer": true,
        "long_complementary_run": true
      },
      {
        "a": "reverse",
        "b": "reverse",
        "self": true,
        "longest_complementary_run": 0,
        "run_start_in_a": null,
        "three_prime_run_a": 0,
        "three_prime_run_b": 0,
        "three_prime_dimer": false,
        "long_complementary_run": false
      }
    ]
  }
}
```

## 4) check_primer_secondary_structure

- **作用：** 筛查引物的发夹结构：由不少于 `min_loop` 个碱基的环封闭的、最大化的完全 Watson-Crick 配对茎。

- **处理方法：** 枚举所有位置对 (i, j)：茎是嵌套的 A-T/C-G 配对序列 (seq[i+k], seq[j-k])；只报告最大化的茎（外侧配对不能再向外延伸，向内延伸直到环长小于 `min_loop`）。短于 `min_stem` 的茎被丢弃，发夹按茎长排序。不允许错配、凸起和 G-U 配对，也没有能量模型，所以这是按茎计数的筛查，不是热力学折叠（不等同于 NUPACK 或 mfold）。

- **API：** `dnakit.design_crispr.check_primer_secondary_structure(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：单条引物字符串、引物列表（至多 1000 条）或映射 `{名称: 引物}`；仅限 ACGT，至多 200 nt。`options`：`min_stem`（默认 4，范围 2-50）最小茎长（碱基对）；`min_loop`（默认 3，范围 0-50）最小未配对环长。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.primers` 每条引物一项：`primer`、`sequence`、`hairpins`（每个含 `stem_length`、`left_start`、`right_start`、`loop_length`、`stem_5prime`、`stem_3prime`、`involves_3prime_end`）、`longest_stem` 和 `three_prime_hairpin`。

- **示例代码：**

```python
import json
from dnakit.design_crispr import check_primer_secondary_structure

result = check_primer_secondary_structure({"p1": "GCGCAAAAGCGC", "p2": "ACACACACACAC"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_secondary_structure",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "82ba75d2df38faa9ec3b57f89e3989391a1df3dfe43ef1202550a994c32f5885"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "primer": "p1",
        "sequence": "GCGCAAAAGCGC",
        "hairpins": [
          {
            "stem_length": 4,
            "left_start": 0,
            "right_start": 8,
            "loop_length": 4,
            "stem_5prime": "GCGC",
            "stem_3prime": "GCGC",
            "involves_3prime_end": true
          }
        ],
        "longest_stem": 4,
        "three_prime_hairpin": true
      },
      {
        "primer": "p2",
        "sequence": "ACACACACACAC",
        "hairpins": [],
        "longest_stem": 0,
        "three_prime_hairpin": false
      }
    ]
  }
}
```

## 5) check_primer_specificity

- **作用：** 统计引物在每条模板序列两条链上的结合位点数（无空位，允许错配）。

- **处理方法：** 对每条模板调用 `dnakit.molbio.match_primer`：在两条链上做无空位搜索，允许最多 `max_mismatches` 个错配，且 3' 端最后 `strict_three_prime_bases` 个碱基不允许错配。汇总所有模板的所有命中；恰好找到 1 个位点时 `specific` 为真。不评估结合能和延伸能力，命中只表示序列匹配，不代表预测会扩增。

- **API：** `dnakit.design_crispr.check_primer_specificity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：单条引物字符串、引物列表或映射 `{名称: 引物}`（ACGT，至多 200 nt）。`options`：`templates`（必填：映射 `{id: 序列}`、序列列表或 FASTA 文本/路径；总长至多 20 Mb）、`max_mismatches`（默认 0，范围 0-10）、`strict_three_prime_bases`（默认 0，范围 0-50）。缺少 `templates` 时抛出配置错误。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result.primers` 每条引物一项：`primer`、`sequence`、`total_hits`、`specific` 和 `hits`（每个含 `template`、`strand`、正链 0 起始半开区间 `start`/`end`、自 5' 端起的 `mismatch_positions`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import check_primer_specificity

primer = "GAATTCGGGACG"
template = "TTT" + primer + "AAAA" + "CGTCCCGAATTC" + "CC"   # second site is on the reverse strand
result = check_primer_specificity(primer, options={"templates": {"t1": template}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_specificity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 12,
    "sequence_bases": 12,
    "lengths": [
      12
    ],
    "sha256": "7994e92a211084996bdf46d96ef74833b03dd460002c7633153f9953445df383"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "primer": "primer1",
        "sequence": "GAATTCGGGACG",
        "total_hits": 2,
        "specific": false,
        "hits": [
          {
            "template": "t1",
            "strand": "forward",
            "start": 3,
            "end": 15,
            "mismatch_positions": []
          },
          {
            "template": "t1",
            "strand": "reverse",
            "start": 19,
            "end": 31,
            "mismatch_positions": []
          }
        ]
      }
    ]
  }
}
```

## 6) design

- **作用：** 为目标序列设计 CRISPR 向导 RNA，或围绕切割位点设计 HDR 供体模板（`options["entity"]` 取 `"grna"` 或 `"donor_template"`）。

- **处理方法：** `grna`：用 `dnakit.patterns.scan_pam_candidates` 在两条链上枚举 PAM 相邻的向导（SpCas9 `NGG` 3' PAM、20 nt；SaCas9 `NNGRRT`、21 nt；AsCas12a `TTTV` 5' PAM、20 nt；或自定义 PAM），再按 GC 范围、同聚物长度、`TTTT`（RNA Pol III 终止子）以及需排除的基序/酶切位点（双链）过滤，按 `proxy_score` 再按位置排序。`proxy_score` 是基于经验规则的启发式分数：从 100 起，GC 不在 40-70 % 扣 15（不在 20-80 % 则改扣 40），含 `TTTT` 扣 40，A/C/G 同聚物长度 >= 4 扣 15；没有用编辑数据拟合，不是预测的效率。提供 `references` 时，对每个向导做 Hamming 脱靶计数（无空位、双链、要求 PAM）：`exact_sites`（完全匹配）、`near_sites`（1..`max_mismatches` 个错配）和 `off_target_sites` = near_sites + max(exact_sites - 1, 0)。SpCas9/SaCas9 的 `cut_site` 为距 PAM 3 bp 的边界，其他核酸酶为 `null`。`donor_template`：供体 = 左同源臂 + `insert` + 右同源臂，切割位点之后 `delete_length` 个碱基被删除；不考虑 HDR 效率，也不加 PAM 阻断突变。

- **API：** `dnakit.design_crispr.design(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：目标 DNA（字符串、FASTA 文本/路径或 `{"sequence": ...}`；单条序列，ACGT）。`options`：`entity`（`"grna"` 默认 | `"donor_template"`）。`grna` 模式：`nuclease`（`"SpCas9"` 默认、`"SaCas9"`、`"AsCas12a"`）或自定义 `pam`/`pam_side`（`"3prime"` | `"5prime"`）/`guide_length`；`strand`（`"both"` 默认 | `"forward"` | `"reverse"`）；`min_gc` 0.4、`max_gc` 0.7、`max_homopolymer` 4（更长的被拒绝）、`avoid_poly_t` `true`、`exclude_motifs`（IUPAC 列表）、`avoid_enzymes`（内置限制酶名称）、`max_guides` 100，可选 `references`（映射/FASTA）和 `max_mismatches`（默认 3）。`donor_template` 模式：`cut_site`（必填，0 起始边界）、`insert`（默认空）、`delete_length`（默认 0）、`left_arm`/`right_arm`（默认 60）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`grna` 模式下 `result` 含 `target`、`nuclease`、`candidates_scanned` 和 `guides`（每个含 `id`、`sequence`、`pam`、`strand`、正链 0 起始半开区间 `start`/`end`、`pam_start`、`cut_site`、`gc`、`longest_homopolymer`、`proxy_score`，提供 `references` 时还有 `exact_sites`、`near_sites`、`off_target_sites`）；`warnings` 声明 `proxy_score` 为启发式分数。`donor_template` 模式下 `result` 含 `target`、`donor`、`donor_length`、`left_arm`/`right_arm`（`start`、`end`、`sequence`）、`insert`、`deleted` 和 `edited_sequence`。同源臂超出目标序列时抛出配置错误。

- **示例代码：**

```python
import json
from dnakit.design_crispr import design

target = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
result = design(target, options={"nuclease": "SpCas9", "max_guides": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 46,
    "sequence_bases": 46,
    "lengths": [
      46
    ],
    "sha256": "899ca71f07ea07f2bf5addd556b05562c9a7750d20636de9e9f49275ea9b61f9"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "target": "seq1",
    "nuclease": "SpCas9",
    "candidates_scanned": 4,
    "guides": [
      {
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "pam": "AGG",
        "strand": "forward",
        "start": 0,
        "end": 20,
        "pam_start": 20,
        "cut_site": 17,
        "gc": 0.5,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g1"
      },
      {
        "sequence": "AGCTTCAAGGCCAGCTAGCT",
        "pam": "AGG",
        "strand": "forward",
        "start": 13,
        "end": 33,
        "pam_start": 33,
        "cut_site": 30,
        "gc": 0.55,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g2"
      },
      {
        "sequence": "GGCCAGCTAGCTAGGCTAAT",
        "pam": "CGG",
        "strand": "forward",
        "start": 21,
        "end": 41,
        "pam_start": 41,
        "cut_site": 38,
        "gc": 0.55,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g3"
      }
    ]
  }
}
```

## 7) encode_data_into_dna

- **作用：** 把字节编码成 DNA 字符串，使用固定的每碱基 2 比特编码（每字节 4 个碱基）。原始字节和文件按其真实字节编码。这只是一个最简的演示映射，不是可用于合成的 DNA 存储方案。

- **处理方法：** 输入先转为字节（见“输入”），每个字节写成 4 个碱基，高位的两比特在前，A = 00、C = 01、G = 10、T = 11（`"Hi!"` = 字节 0x48 0x69 0x21，得到 `CAGA CGGC AGAC`；字节 0x00 0x48 0xFF 0xB2 得到 `AAAA CAGA TTTT GTAG`）。该编码对 0-255 的每个字节值都是无损的：`dnakit.design_crispr._impl_local.dna_to_bytes` 是它的精确逆过程（要求只含 ACGT 且长度是 4 的倍数），共享运行时中的 `decode_dna_storage` 假定的也是同一映射（该运行时分支没有被任何公开模块暴露，且只解码为 UTF-8 文本；二进制数据请用 `dna_to_bytes`）。编码没有索引、头部、校验、纠错、引物，也不限制 GC 含量和均聚物长度，因此输出不能直接用于合成。`dnakit.genomics.encode_dna_storage` 与 `decode_dna_storage` 提供带完整性摘要的分段索引编码（需要可选依赖 `dnakit[dna-storage]`）。

- **API：** `dnakit.design_crispr.encode_data_into_dna(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`bytes`/`bytearray`/`memoryview`（按原样编码，任意字节值）、已存在文件的 `Path`（编码文件的字节）、`str`（其 UTF-8 字节，所以 `"é"` 得到 8 个碱基；设置 `options["source"] = "path"` 时字符串被当作文件路径），或映射/列表（键排序的 JSON 文本，`ensure_ascii=False`，UTF-8）。至多 1,000,000 字节（超出抛出 `ConfigurationError`，`CRISPR_INPUT_TOO_LARGE`）；文件不存在、类型不支持（如 `int`）或映射内容不能转 JSON 时抛出 `InputFormatError`。`options`：`source`（`"auto"` 默认 | `"text"` | `"path"`）；`"auto"` 和 `"text"` 下字符串始终是文本，不会被当作路径。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`（`"encode_data_into_dna"`）、`provenance`（`implementation`、`backend`）、`input_summary`（通用文本摘要：`record_count`、`text_length`、`sequence_bases`、`lengths`、`sha256`，按输入的文本形式计算）、`warnings` 和 `result`。`result` 是由 `A`、`C`、`G`、`T` 组成的字符串，长度为 4 x 字节数。注意 `input_summary` 是按文本描述输入的（原始字节显示为 UTF-8 解码结果，不是合法 UTF-8 时显示为 base64 的 JSON）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import encode_data_into_dna

result = encode_data_into_dna(bytes([0x00, 0x48, 0xFF, 0xB2]))   # raw bytes, encoded as they are
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "encode_data_into_dna",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 4,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d909adf4e65a93c9f454770461199a5eff63c53a49a20aeff3418ae597a85212"
  },
  "warnings": [],
  "status": "ok",
  "result": "AAAACAGATTTTGTAG"
}
```

## 8) evaluate_primer_pair

- **作用：** 评估一对正向/反向引物：两条引物各自的 Tm 及其差值、GC 含量比例、3' 端 GC 钳、自身互补和相互互补，并按文档化的阈值给出标记列表。

- **处理方法：** 每条引物按 `calculate_primer_tm` 计算 Tm（两条使用同一个 `method`；`auto` 下若一对引物分居 13 nt 两侧，会用两种不同方法并给出警告，此时 `delta_tm_celsius` 没有意义）。`delta_tm_celsius` = |Tm(正向) - Tm(反向)|。`gc_clamp` = 引物最后 5 个碱基（3' 端）中 G/C 的个数。`longest_self_complementary` = 引物中最长的连续片段，其反向互补序列也出现在同一条引物中（两个拷贝反向平行退火的自身二聚体情形；回文引物如 `GAATTC` 得到全长；发夹结构请用 `check_primer_secondary_structure`）。`longest_cross_complementary` = 正向引物中最长的连续片段，它是反向引物某一片段的反向互补。互补性是精确的连续 Watson-Crick 配对，不计算自由能，也不考虑错配和凸起（3' 端互补请用 `check_primer_dimer`）。标记（括号内为默认阈值，均可覆盖）：delta Tm > `max_delta_tm` [5.0 C] 时为 `tm_mismatch`；逐条引物的 `length_out_of_range`（超出 `min_length`-`max_length` [18-30]）、`gc_out_of_range`（超出 `min_gc`-`max_gc` [0.40-0.60]）、`gc_clamp_out_of_range`（超出 `min_gc_clamp`-`max_gc_clamp` [1-3]）、`self_complementarity`（自身互补长度 > `max_self_complementary` [4]）；以及相互互补长度 > `max_cross_complementary` [4] 时的 `cross_complementarity`。这些是常用的设计经验阈值，并不保证 PCR 成功。

- **API：** `dnakit.design_crispr.evaluate_primer_pair(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{"forward": ..., "reverse": ...}`（恰好这两个键）或两条引物的列表 `[forward, reverse]`，均为 5'->3' 书写的字符串（反向引物按订购序列给出，不要再取反向互补）；只能含 ACGT（`rna` 为 true 时可含 `U`）。其他输入抛出 `InputFormatError`。`options`：与 `calculate_primer_tm` 相同的 `method`、`rna`、`sodium_molar`、`potassium_molar`、`strand_concentration_molar`，以及阈值 `max_delta_tm`、`min_length`、`max_length`、`min_gc`、`max_gc`、`min_gc_clamp`、`max_gc_clamp`、`max_self_complementary`、`max_cross_complementary`。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`（`"evaluate_primer_pair"`）、`provenance`（`implementation`、`backend`）、`input_summary`（通用文本摘要：`record_count`、`text_length`、`sequence_bases`、`lengths`、`sha256`，按输入的文本形式计算）、`warnings` 和 `result`。`result` 含 `forward` 和 `reverse`（各有 `sequence`、`length`、`gc_fraction`、`tm_celsius`、`method`、`gc_clamp`、`longest_self_complementary`）、`delta_tm_celsius`、`longest_cross_complementary`、`thresholds`（实际使用的阈值）、`flags`（`{"flag", "primer", "value", "threshold"}` 的列表；`primer` 为 `"forward"`、`"reverse"`，成对级别的标记为 `null`）、`pass`（没有任何标记时为 true），计算过最近邻 Tm 时还有 `conditions`。

- **示例代码：**

```python
import json
from dnakit.design_crispr import evaluate_primer_pair

result = evaluate_primer_pair({"forward": "GACGTTAGCATCAGGCTAAC", "reverse": "CCATGTTCAGCTAGGATCAT"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "evaluate_primer_pair",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3fdd37d2399c743d84981ed823a66973cf44c6c25c083abdbdc8bdd2613dda0b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "forward": {
      "sequence": "GACGTTAGCATCAGGCTAAC",
      "length": 20,
      "gc_fraction": 0.5,
      "tm_celsius": 52.86,
      "method": "nearest_neighbor",
      "gc_clamp": 2,
      "longest_self_complementary": 6
    },
    "reverse": {
      "sequence": "CCATGTTCAGCTAGGATCAT",
      "length": 20,
      "gc_fraction": 0.45,
      "tm_celsius": 50.99,
      "method": "nearest_neighbor",
      "gc_clamp": 1,
      "longest_self_complementary": 4
    },
    "delta_tm_celsius": 1.87,
    "longest_cross_complementary": 4,
    "thresholds": {
      "max_delta_tm": 5.0,
      "min_length": 18,
      "max_length": 30,
      "min_gc": 0.4,
      "max_gc": 0.6,
      "min_gc_clamp": 1,
      "max_gc_clamp": 3,
      "max_self_complementary": 4,
      "max_cross_complementary": 4
    },
    "flags": [
      {
        "flag": "self_complementarity",
        "primer": "forward",
        "value": 6,
        "threshold": 4
      }
    ],
    "pass": false,
    "conditions": {
      "sodium_molar": 0.05,
      "potassium_molar": 0.0,
      "strand_concentration_molar": 2.5e-07
    }
  }
}
```

## 9) export

- **作用：** 把向导 RNA 序列化为 TSV、CSV、FASTA 或 JSON 文本，可选写入文件。

- **处理方法：** 先解析并校验向导（ACGT），再按所选格式输出。TSV/CSV 的列为 `id`、`sequence`，其后是各记录额外字段（例如 `design` 返回的 `start`、`strand`、`proxy_score`）的排序并集；缺少某字段的记录该单元格为空。FASTA 输出 `>id` 和序列；JSON 按排序键输出记录。只有 `overwrite` 为真时才覆盖已有输出文件。

- **API：** `dnakit.design_crispr.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表（可直接使用 `design` 的返回，也接受 `{"guides": [...]}`）。`options`：`format`（`"tsv"` 默认 | `"csv"` | `"fasta"` | `"json"`）；`output`（可选文件路径）；`overwrite`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `format`、`count`、`text`（序列化后的向导）和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import export

guides = [{"id": "g1", "sequence": "GACGTACGTCGTAAGCTTCA", "start": 0, "strand": "forward"},
          {"id": "g2", "sequence": "ACGTACGTTTTTACGTACGT"}]
result = export(guides, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 133,
    "sequence_bases": 40,
    "lengths": [
      20,
      20
    ],
    "sha256": "b99bcb771c0c55eaedd9ee50be13be0a112cdaf706aaa2e103450482ec3eca20"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "count": 2,
    "text": "id\tsequence\tstart\tstrand\ng1\tGACGTACGTCGTAAGCTTCA\t0\tforward\ng2\tACGTACGTTTTTACGTACGT\t\t\n",
    "output_path": null
  }
}
```

## 10) filter

- **作用：** 按 GC 含量、同聚物、poly-T、长度和需排除的基序过滤向导 RNA，并给出每个被拒向导的原因。

- **处理方法：** 每个向导按固定的经验设计规则检查：GC 比例低于 `min_gc`（`gc_below_min`）或高于 `max_gc`（`gc_above_max`）；最长同聚物长于 `max_homopolymer`（`homopolymer_run`）；含 `TTTT`，即 RNA Pol III 终止子（`poly_t`，在 `avoid_poly_t` 时）；任一链含需排除的基序或酶切位点（`excluded_motif:<名称>`）；长度不等于 `length`（`length`）。这些是常用的设计经验规则，不是编辑效率模型。

- **API：** `dnakit.design_crispr.filter(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表。`options`：`min_gc`（默认 0.4）、`max_gc`（0.7）、`max_homopolymer`（4：更长则不通过）、`avoid_poly_t`（`true`）、`exclude_motifs`（IUPAC 列表，双链）、`avoid_enzymes`（内置限制酶名称）、`length`（可选，要求的向导长度）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `kept` 和 `rejected`（每条记录有 `id`、`sequence`、`gc`、`longest_homopolymer` 以及违反规则列表 `failed`），以及 `kept_count` 和 `rejected_count`。

- **示例代码：**

```python
import json
from dnakit.design_crispr import filter as filter_guides

guides = {"ok": "GACGTACGTCGTAAGCTTCA", "at_rich": "ATATATATATATATATATAT",
          "poly_t": "GACTTTTCGACTGACTGACT", "ecori": "GACGAATTCGTAAGCTTCAG"}
result = filter_guides(guides, options={"avoid_enzymes": ["EcoRI"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "filter",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 132,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "02704497d75ca487b918a3590bdae63dbfa789f44f7116c23d2bba76448b6562"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kept": [
      {
        "id": "ok",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "gc": 0.5,
        "longest_homopolymer": 2,
        "failed": []
      }
    ],
    "rejected": [
      {
        "id": "at_rich",
        "sequence": "ATATATATATATATATATAT",
        "gc": 0.0,
        "longest_homopolymer": 1,
        "failed": [
          "gc_below_min"
        ]
      },
      {
        "id": "poly_t",
        "sequence": "GACTTTTCGACTGACTGACT",
        "gc": 0.45,
        "longest_homopolymer": 4,
        "failed": [
          "poly_t"
        ]
      },
      {
        "id": "ecori",
        "sequence": "GACGAATTCGTAAGCTTCAG",
        "gc": 0.45,
        "longest_homopolymer": 2,
        "failed": [
          "excluded_motif:EcoRI"
        ]
      }
    ],
    "kept_count": 1,
    "rejected_count": 3
  }
}
```

## 11) plot

- **作用：** 把向导 RNA 在目标序列上的位置画成 SVG 图（正链向导在线上方，反链向导在线下方）。

- **处理方法：** 手写的简单 SVG（不依赖绘图库）：一条与目标长度对应的基线，每个向导一个 `<rect class="guide">`，按 `start` 到 `end` 线性定位，正链向导为蓝色位于线上方，其余为红色位于线下方。每个矩形带有 `data-start`、`data-end` 和 `data-strand` 属性。图中只显示位置，不显示分数。

- **API：** `dnakit.design_crispr.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射 `{"length": N, "guides": [...]}` 或 `{"sequence": ..., "guides": [...]}`；每个向导需要 0 <= `start` < `end`（0 起始半开）和 `strand`（`forward`/`reverse`），例如 `design` 的 `result` 中的向导。缺少长度时取最大的 `end`。`options`：`output`（可选 `.svg` 路径）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `svg`（SVG 文本）和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import plot

guides = [{"start": 0, "end": 20, "strand": "forward"}, {"start": 30, "end": 50, "strand": "reverse"}]
result = plot({"length": 60, "guides": guides})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "design_crispr",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 119,
    "sequence_bases": 60,
    "lengths": [
      60
    ],
    "sha256": "8f8a17b9efbb550b409b4f24f119aadff3263fc6863ad7764bffe66e5c0b8b80"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"90\"><line x1=\"20\" y1=\"45\" x2=\"780\" y2=\"45\" stroke=\"#333\" stroke-width=\"2\"/><text x=\"20\" y=\"85\" font-...",
    "output_path": null
  }
}
```

## 12) predict_editing_outcomes

- **作用：** 预测 CRISPR 修复结果。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算，也不编造结果。插入缺失和修复结果预测需要 inDelphi、FORECasT、Lindel 之类的已训练模型及其权重，DNAKit 并未内置。只需要 HDR 几何结构而不涉及结果模型时，请用 `design`（`entity="donor_template"`）和 `simulate`。

- **API：** `dnakit.design_crispr.predict_editing_outcomes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。此时 `status` 为 `"not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），没有 `result`。`required_dependency` 为 `a trained repair-outcome model (inDelphi, FORECasT, Lindel) and its weights`，`message` 为 `Indel/repair outcome prediction needs trained model weights that are not bundled.`

- **示例代码：**

```python
import json
from dnakit.design_crispr import predict_editing_outcomes

result = predict_editing_outcomes("GACGTACGTCGTAAGCTTCAAGG")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "predict_editing_outcomes",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 23,
    "sequence_bases": 23,
    "lengths": [
      23
    ],
    "sha256": "e17f729834a93bd40c59a9c547a2621d3e85a9c76868193434768ff96cb3b1d2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained repair-outcome model (inDelphi, FORECasT, Lindel) and its weights",
  "message": "Indel/repair outcome prediction needs trained model weights that are not bundled."
}
```

## 13) rank

- **作用：** 按某个数值字段对向导 RNA 排序；默认的 `proxy_score` 是基于经验规则的启发式分数（见 `score`）。

- **处理方法：** 按所选字段排序，默认降序，`ascending` 为真则升序；并列时保持输入顺序（稳定排序），`top_n` 截取前若干项。`by` 为 `proxy_score` 且记录中没有该字段时，用与 `score` 相同的规则代理计算（从 100 起；GC 不在 40-70 % 扣 15，不在 20-80 % 改扣 40；含 `TTTT` 扣 40；A/C/G 同聚物 >= 4 扣 15）。`proxy_score` 是基于经验规则的启发式分数，没有用编辑实验数据拟合或验证，不是预测的编辑效率（需要效率预测请用 Rule Set 2、DeepSpCas9 等）。其他字段（例如 `design` 给出的 `off_target_sites`）在每条记录中都必须是数值。

- **API：** `dnakit.design_crispr.rank(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表。`options`：`by`（默认 `"proxy_score"`）、`ascending`（默认 `false`）、`top_n`（可选）。记录缺少数值型的 `by` 字段时抛出输入错误。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `by`、`ascending` 和 `ranked`（记录列表，使用时加入 `proxy_score`，每条带 1 起始的 `rank`）。`by` 为 `proxy_score` 时，`warnings` 声明其为启发式分数。

- **示例代码：**

```python
import json
from dnakit.design_crispr import rank

guides = ["GGGGGGGGGGGGGGGGGGGG", "GACGTACGTCGTAAGCTTCA", "ACGTACGTTTTTACGTACGT"]
result = rank(guides)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "rank",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 72,
    "sequence_bases": 60,
    "lengths": [
      20,
      20,
      20
    ],
    "sha256": "9bb60ae6e4e1ebb70a61019350afe0fb603f219157932ffb6cd98efba658b434"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "by": "proxy_score",
    "ascending": false,
    "ranked": [
      {
        "id": "guide2",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "proxy_score": 100.0,
        "rank": 1
      },
      {
        "id": "guide3",
        "sequence": "ACGTACGTTTTTACGTACGT",
        "proxy_score": 60.0,
        "rank": 2
      },
      {
        "id": "guide1",
        "sequence": "GGGGGGGGGGGGGGGGGGGG",
        "proxy_score": 45.0,
        "rank": 3
      }
    ]
  }
}
```

## 14) remove_restriction_sites

- **作用：** 通过最少的碱基替换去除 DNA 序列中的限制酶识别位点。

- **处理方法：** 调用 `dnakit.molbio.optimize_sequence_rules`（贪心、给定随机种子下结果确定的优化器），把识别序列作为双链上的禁用基序；编辑前后分别统计位点数（计数支持 IUPAC 模式）。`preserve_translation` 为真时输入必须是读框内的编码序列，且只做同义替换。在约束下无法去除的位点会保留，此时 `removed_all` 为假。不检查替换是否产生了未列出的新位点，也不评估对表达的影响。

- **API：** `dnakit.design_crispr.remove_restriction_sites(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：DNA 字符串、FASTA 文本/路径或 `{"sequence": ...}`（单条序列）。`options`：`enzymes`（内置目录中的酶名，如 `"EcoRI"`）和/或 `sites`（IUPAC 识别序列），至少给出一个；`preserve_translation`（默认 `false`）；`seed`（默认 0）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `sequence_id`、`sequence`（编辑后）、`sites_before`、`sites_after`、`removed_all` 和 `changes`（每项含 0 起始的 `position`、`before`、`after`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import remove_restriction_sites

result = remove_restriction_sites("TTTTGAATTCAAAACCCGGATCCTTTT", options={"enzymes": ["EcoRI", "BamHI"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "remove_restriction_sites",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 27,
    "sequence_bases": 27,
    "lengths": [
      27
    ],
    "sha256": "f370022c244baaa8f4faedb685da22ddd8021fff5f05b2aaa52b195cdc90afe7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_id": "seq1",
    "sequence": "TTTTGCATTCAAAACCCAGATCCTTTT",
    "sites_before": 2,
    "sites_after": 0,
    "removed_all": true,
    "changes": [
      {
        "position": 5,
        "before": "A",
        "after": "C"
      },
      {
        "position": 17,
        "before": "G",
        "after": "A"
      }
    ]
  }
}
```

## 15) run_crispr_screen_design

- **作用：** 设计汇集式 CRISPR 向导文库：为每条目标序列选出过滤后最好的向导。

- **处理方法：** 对每条目标序列，向导的枚举和过滤与 `design` 相同，按启发式 `proxy_score` 排序后贪心挑选，直到选满 `guides_per_target` 个；切割位点与已选向导相距小于 `min_cut_spacing` bp 的被跳过（无 Cas9 切割位点时用向导起点），提供 `references` 和 `max_off_target_sites` 时，非自身 Hamming 位点数超限的向导也被跳过。`proxy_score` 是基于经验规则的启发式分数，没有用编辑实验数据拟合或验证，不是预测的编辑效率（需要效率预测请用 Rule Set 2、DeepSpCas9 等）。不生成非靶向对照，也不使用基因或外显子注释。

- **API：** `dnakit.design_crispr.run_crispr_screen_design(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：目标序列，可以是映射 `{target_id: DNA}`、序列列表或 FASTA 文本/路径（例如外显子），至多 10000 条。`options`：`guides_per_target`（默认 4）；`design` 的枚举与过滤选项（`nuclease`/`pam`、`min_gc` 0.4、`max_gc` 0.7、`max_homopolymer` 4、`avoid_poly_t`、`exclude_motifs`、`avoid_enzymes`）；`min_cut_spacing`（默认 0）；可选 `references`、`max_mismatches`（默认 3）和 `max_off_target_sites`（默认不限）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `guides_per_target`、`target_count`、`guide_count`、`library`（每行含 `target`、`guide_id`、`design` 的向导字段，提供 `references` 时还有脱靶计数）和 `targets_short`（得到的向导数少于要求的目标，含 `requested` 和 `designed`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import run_crispr_screen_design

exon = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
result = run_crispr_screen_design({"exon1": exon}, options={"guides_per_target": 2, "min_cut_spacing": 5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "run_crispr_screen_design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7e584ed8cb48bf6be7db6302c06095c9fde603d3140b4f062b450b482e52a1fb"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "guides_per_target": 2,
    "target_count": 1,
    "guide_count": 2,
    "library": [
      {
        "target": "exon1",
        "guide_id": "exon1_g1",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "pam": "AGG",
        "strand": "forward",
        "start": 0,
        "end": 20,
        "pam_start": 20,
        "cut_site": 17,
        "gc": 0.5,
        "proxy_score": 100.0
      },
      {
        "target": "exon1",
        "guide_id": "exon1_g2",
        "sequence": "AGCTTCAAGGCCAGCTAGCT",
        "pam": "AGG",
        "strand": "forward",
        "start": 13,
        "end": 33,
        "pam_start": 33,
        "cut_site": 30,
        "gc": 0.55,
        "proxy_score": 100.0
      }
    ],
    "targets_short": []
  }
}
```

## 16) scan_grna_offtargets

- **作用：** 在提供的序列中，用精确的无空位 Hamming 搜索在两条链上查找向导 RNA 的脱靶位点。

- **处理方法：** 把每条参考序列两条链上与向导等长的每个窗口与向导比较；错配位置数 <= `max_mismatches` 且（设置了 `pam` 时）相邻 PAM（3' 或 5' 侧）符合 IUPAC 模式时报告为位点。向导自身的靶位点也会被报告（0 个错配）。不考虑凸起和插入缺失，也没有 CFD/MIT 式的位置加权评分，错配数是唯一的排序依据；搜索为暴力枚举，向导数 x 参考碱基数 x 2 超过 5000 万时拒绝执行。参考序列允许含简并碱基，参考中的简并碱基按错配计。

- **API：** `dnakit.design_crispr.scan_grna_offtargets(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表（至多 1000 条）。`options`：`references`（必填：映射 `{id: 序列}`、序列列表或 FASTA 文本/路径）、`max_mismatches`（默认 3，范围 0-6）、`pam`（IUPAC，默认 `"NGG"`；`null` 表示忽略 PAM）、`pam_side`（`"3prime"` 默认 | `"5prime"`）。缺少 `references` 时抛出配置错误。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `max_mismatches`、`pam` 和 `guides`（每个向导有 `id`、`sequence`、`hit_count`、`counts_by_mismatch` 和按错配数排序的 `hits`，每个命中含 `reference`、`strand`（`+`/`-`）、正链 0 起始半开区间 `start`/`end`、`mismatches`、自向导 5' 端起的 1 起始 `mismatch_positions` 和匹配的 `site`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import scan_grna_offtargets

guide = "GACGTACGTCGTAAGCTTCA"
one_mismatch = guide[:9] + "A" + guide[10:]
reference = guide + "TGG" + "TT" + one_mismatch + "CGG" + "TT" + one_mismatch + "AAA"   # last site has no NGG PAM
result = scan_grna_offtargets(guide, options={"references": {"chr1": reference}, "max_mismatches": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "scan_grna_offtargets",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 20,
    "lengths": [
      20
    ],
    "sha256": "b015b64a47743134d8b718d27eacebc4a00519b9304d5f0f27607012ff7efa0c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "max_mismatches": 3,
    "pam": "NGG",
    "guides": [
      {
        "id": "seq1",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "hit_count": 2,
        "counts_by_mismatch": {
          "0": 1,
          "1": 1
        },
        "hits": [
          {
            "reference": "chr1",
            "strand": "+",
            "start": 0,
            "end": 20,
            "mismatches": 0,
            "mismatch_positions": [],
            "site": "GACGTACGTCGTAAGCTTCA"
          },
          {
            "reference": "chr1",
            "strand": "+",
            "start": 25,
            "end": 45,
            "mismatches": 1,
            "mismatch_positions": [
              10
            ],
            "site": "GACGTACGTAGTAAGCTTCA"
          }
        ]
      }
    ]
  }
}
```

## 17) score

- **作用：** 用基于经验规则的启发式靶向代理分数（0 到 100）给向导 RNA 打分；它不是训练得到的效率模型。

- **处理方法：** 代理分数从 100 起：GC 比例不在 40-70 % 时扣 15（不在 20-80 % 则改扣 40），含 `TTTT`（RNA Pol III 终止子）扣 40，A/C/G 同聚物长度 >= 4 扣 15，最低为 0。它只体现常见的设计经验规则。`proxy_score` 是基于经验规则的启发式分数，没有用编辑实验数据拟合或验证，不是预测的编辑效率（需要效率预测请用 Rule Set 2、DeepSpCas9 等）。向导须为 10-40 nt 的 ACGT（U 视为 T）。

- **API：** `dnakit.design_crispr.score(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表。`options`：未使用。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `method`（`"heuristic_rule_proxy"`）和 `scores`（每个向导有 `id`、`sequence`、`proxy_score`、`gc`、`longest_homopolymer`、`poly_t`）；`warnings` 声明该分数为启发式。向导长度不在 10-40 nt 或含其他字符时抛出输入错误。

- **示例代码：**

```python
import json
from dnakit.design_crispr import score

result = score({"balanced": "GACGTACGTCGTAAGCTTCA", "poly_t": "ACGTACGTTTTTACGTACGT",
                "all_g": "GGGGGGGGGGGGGGGGGGGG"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 103,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6c8cb8de2dc6b2a2809628762c7b49bc5a2f66b8c7654fb54181ebb6d99e9f63"
  },
  "warnings": [
    "Heuristic rule-based proxy, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "method": "heuristic_rule_proxy",
    "scores": [
      {
        "id": "balanced",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "proxy_score": 100.0,
        "gc": 0.5,
        "longest_homopolymer": 2,
        "poly_t": false
      },
      {
        "id": "poly_t",
        "sequence": "ACGTACGTTTTTACGTACGT",
        "proxy_score": 60.0,
        "gc": 0.4,
        "longest_homopolymer": 5,
        "poly_t": true
      },
      {
        "id": "all_g",
        "sequence": "GGGGGGGGGGGGGGGGGGGG",
        "proxy_score": 45.0,
        "gc": 1.0,
        "longest_homopolymer": 20,
        "poly_t": false
      }
    ]
  }
}
```

## 18) simulate

- **作用：** 在计算机中模拟同源定向修复：把供体同源臂之间的片段换入目标序列。

- **处理方法：** 按给定的臂长把供体拆成左臂、载荷和右臂。两条臂都必须在目标序列正链上恰好出现一次，且左臂在右臂之前；两臂之间的目标序列被载荷替换。不模拟修复通路竞争、切割效率、NHEJ 插入缺失以及臂内突变，所以这是字符串替换而不是效率预测。臂不存在、重复或顺序错误时抛出输入错误。

- **API：** `dnakit.design_crispr.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：目标 DNA（字符串、FASTA 文本/路径或 `{"sequence": ...}`；单条序列）。`options`：`donor`（必填，左臂 + 载荷 + 右臂）、`left_arm` 和 `right_arm`（必填，臂长，须能放入供体）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `target`、`edited_sequence`、`region`（被替换区域的 0 起始半开 `start`/`end`）、`deleted`、`inserted` 和 `length_change`。

- **示例代码：**

```python
import json
from dnakit.design_crispr import simulate

target = "AAAACCCCTTGGGGTTTT"
donor = "CCCC" + "AA" + "GGGG"          # left arm + payload + right arm
result = simulate(target, options={"donor": donor, "left_arm": 4, "right_arm": 4})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 18,
    "lengths": [
      18
    ],
    "sha256": "82435c45a7de0432d66f77bce4656ba8124463ba82f4c4688556c6b902bf5651"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "target": "seq1",
    "edited_sequence": "AAAACCCCAAGGGGTTTT",
    "region": {
      "start": 8,
      "end": 10
    },
    "deleted": "TT",
    "inserted": "AA",
    "length_change": 0
  }
}
```

## 19) validate

- **作用：** 校验向导 RNA：字符、长度、（可选）在目标序列中的 PAM 相邻性，以及设计规则警告。

- **处理方法：** 以下错误使向导无效：含 ACGT 之外的字符（`invalid_characters`）、长度不等于期望长度（`length_<n>_expected_<m>`），以及给定 `target` 时向导在任一链上都没有与匹配 PAM 相邻的出现位置（`no_pam_adjacent_site_in_target`）。警告不使向导无效：GC 不在 40-70 %（`gc_outside_40_70`）、同聚物长于 4（`homopolymer_run_gt_4`）、含 `TTTT`（`poly_t`）。靶位点通过两条链上的精确字符串匹配查找，不容许错配，也不评估效率或脱靶。

- **API：** `dnakit.design_crispr.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：向导 RNA 可以是字符串列表、FASTA 文本/路径、映射 `{id: 向导序列}`，或带 `guide`/`sequence` 字段的记录列表。`options`：`length`（期望的向导长度，默认取核酸酶默认值，SpCas9 为 20）；`nuclease`（`"SpCas9"` 默认、`"SaCas9"`、`"AsCas12a"`）或与 `design` 相同的自定义 `pam`/`pam_side`/`guide_length`；`target`（可选 DNA；给定后向导必须在其任一链上与匹配的 PAM 相邻出现）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"design_crispr"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `nuclease`、`pam`、`expected_length`、`all_valid` 和 `guides`（每个向导有 `id`、`sequence`、`valid`、`errors`、`warnings` 以及 `target_sites`，其中有 `strand`、正链 0 起始半开区间 `start`/`end`、`pam`）。

- **示例代码：**

```python
import json
from dnakit.design_crispr import validate

target = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
guides = {"good": "GACGTACGTCGTAAGCTTCA", "short": "GACGTACGTCGTAAGCTTC", "no_pam": "CATCATCATCATCATCATCA"}
result = validate(guides, options={"target": target})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "design_crispr",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 98,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d888e1a0acbc03ba69ad873c3fd9569ce541e899e96b5a9b738889ecb86b1d1c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "nuclease": "SpCas9",
    "pam": "NGG",
    "expected_length": 20,
    "guides": [
      {
        "id": "good",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "valid": true,
        "errors": [],
        "warnings": [],
        "target_sites": [
          {
            "strand": "forward",
            "start": 0,
            "end": 20,
            "pam": "AGG"
          }
        ]
      },
      {
        "id": "short",
        "sequence": "GACGTACGTCGTAAGCTTC",
        "valid": false,
        "errors": [
          "length_19_expected_20"
        ],
        "warnings": [],
        "target_sites": []
      },
      {
        "id": "no_pam",
        "sequence": "CATCATCATCATCATCATCA",
        "valid": false,
        "errors": [
          "no_pam_adjacent_site_in_target"
        ],
        "warnings": [
          "gc_outside_40_70"
        ],
        "target_sites": []
      }
    ],
    "all_valid": false
  }
}
```
