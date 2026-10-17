# 实验设计

本类收录 **9** 项有文档的功能，全部由 DNAKit 使用标准库原生计算（正常调用均不会返回 `not_implemented`）：样本量、效能、覆盖度和重复数计算，随机化布局（CRD、RCBD、拉丁方、析因），标准顺序析因矩阵，运行顺序随机化，处理分配，设计检查与校验，CSV/TSV/JSON 导出，以及蒙特卡洛效能模拟。随机化可通过 `seed` 复现。这些函数通过 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) assign

- **作用：** 用（分层）置换区组随机化把实验单位随机分配到处理组。

- **处理方法：** 在每个层内（未设 `stratify_by` 时全部单位为一层），单位被打乱并切成长度为 sum(ratio) 的连续区组；每个区组得到分组模式（每个组名重复 `ratio` 次）的一个随机排列。因此每层内各组大小与目标比例的偏差小于一个区组。随机化使用 Python `random.Random(seed)`，可复现。

- **API：** `dnakit.experimental_design.assign(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：单位 id 列表、单位个数（整数，id 为 `u1`..`un`），或记录（字典，可含分层列）。`options`：`groups`（至少两个互不相同的名称；默认 `G1..Gk`，`n_groups` 默认 2）、`ratio`（正整数，每组一个；默认相等）、`stratify_by`（列名）、`seed`（0）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `seed`、`n_units`、`assignments`（每个单位记录加上 `group`）、`group_sizes` 与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import assign

units = [
    {"unit": "m1", "sex": "F"}, {"unit": "m2", "sex": "F"}, {"unit": "m3", "sex": "F"},
    {"unit": "m4", "sex": "F"}, {"unit": "m5", "sex": "M"}, {"unit": "m6", "sex": "M"},
    {"unit": "m7", "sex": "M"}, {"unit": "m8", "sex": "M"},
]
result = assign(units, options={"groups": ["control", "treated"], "stratify_by": "sex", "seed": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "experimental_design",
  "operation": "assign",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 8,
    "text_length": 224,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "004f3d304d047d6db323ca275724203cedd8d41f637060e2a047d422cc41db24"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 3,
    "n_units": 8,
    "assignments": [
      {
        "unit": "m1",
        "sex": "F",
        "group": "treated"
      },
      {
        "unit": "m2",
        "sex": "F",
        "group": "treated"
      },
      "..."
    ],
    "group_sizes": {
      "control": 4,
      "treated": 4
    },
    "method": "permuted-block randomisation: within each stratum, consecutive blocks of sum(ratio) units receive a random permutation o..."
  }
}
```

## 2) calculate

- **作用：** 为常见设计计算样本量、检验效能、测序覆盖度（Lander-Waterman）或重复数。

- **处理方法：** `quantity` 选择计算内容。`two_sample_t`、`one_sample_t`、`paired_t` 的 `sample_size` / `power` 使用精确的非中心 t 分布（delta = d*sqrt(n1*n2/(n1+n2)) 或 d*sqrt(n)；效应量为 Cohen's d）；`two_proportions` 使用与 R `power.prop.test` 相同的合并零假设方差正态近似；`chi_square` 使用非中心卡方（lambda = n*w^2，w 为 Cohen's w）与 1-alpha 处的中心卡方临界值。样本量为达到目标效能的最小 n（`larger` 与 `smaller` 均按单侧处理）。`coverage` 遵循 Lander & Waterman：C = L*N/G，覆盖比例 1 - exp(-C)，期望 contig 数 N*exp(-C*(1-theta))，期望 contig 长度 L*((exp(C(1-theta)) - 1)/C + 1 - theta)，深度 >= k 的比例取自 Poisson(C)。`replicates` 为使 t 置信区间半宽 t * sd / sqrt(n) <= `margin` 的最小 n（同时报告基于 z 的 n）。

- **API：** `dnakit.experimental_design.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 参数来自 `input` 映射与 `options`（名称冲突时以 `input` 为准）。`quantity` 为 `sample_size`（默认）、`power`、`coverage` 或 `replicates`。检验：`test` 为 `two_sample_t`（默认）、`one_sample_t`、`paired_t`（`effect_size` 为 Cohen's d）、`two_proportions`（`p1`、`p2`）或 `chi_square`（`effect_size` 为 Cohen's w，另需 `df`）；`alpha`（0.05）、`power`（0.8）、`alternative` 为 `two-sided`（默认）/`larger`/`smaller`、`ratio`（n2/n1，默认 1）、`n`（计算效能时必需）。覆盖度：`genome_size`、`read_length`、`n_reads` 或 `coverage`，可选 `target_fraction`、`min_depth`（1）、`overlap_fraction`（0）。重复数：`sd`、`margin`、`confidence`（0.95）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 随计算量而变。样本量：`quantity`、`test`、`alpha`、`alternative`、`target_power`、`n_per_group`、`n`、`n_group2`、`total`（t 检验）或 `n_unrounded`（比例）或 `n_total`（卡方）、`achieved_power` 与 `method`。效能：`n` 与 `power`。覆盖度：`coverage`、`fraction_covered`、`fraction_covered_min_depth`、`expected_uncovered_bases`、`expected_contigs`、`expected_contig_length`，给出 `target_fraction` 时另有 `coverage_needed_for_target_fraction` 与 `reads_needed_for_target_fraction`。重复数：`replicates`、`replicates_normal_approximation`、`sd`、`margin`、`confidence`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import calculate

result = calculate(
    {"effect_size": 0.5, "alpha": 0.05, "power": 0.8},
    options={"quantity": "sample_size", "test": "two_sample_t"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "experimental_design",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 49,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3228c308633c7dde297f035b87166f4969c8f503a7f49fa9c7d59a78fbd49aca"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "sample_size",
    "test": "two_sample_t",
    "alpha": 0.05,
    "alternative": "two-sided",
    "target_power": 0.8,
    "n_per_group": 64,
    "n": 64,
    "n_group2": 64,
    "total": 128,
    "achieved_power": 0.8014595579222591,
    "method": "exact non-central t power; delta = d*sqrt(n1*n2/(n1+n2)) (two-sample) or d*sqrt(n) (one-sample/paired), df = n1+n2-2 or n-1"
  }
}
```

## 3) check

- **作用：** 检查设计表的统计性质：重复是否均衡、区组是否完整、是否满足拉丁方性质，以及编码因子列是否正交。

- **处理方法：** 统计各处理的重复数（全部相等即均衡；重复少于两次的处理会被标记）。存在 `block` 列时检查每个区组是否均含各处理且次数相同；存在 `row` 与 `column` 列时检查每个处理在每行、每列各出现一次；对 `factor_columns` 计算数值列两两点积，全部为零（容差 1e-9）时判定正交。每项未通过的检查写入 `issues`；无问题时 `ok` 为真。

- **API：** `dnakit.experimental_design.check(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：设计记录、`design()` 返回的完整封装（含 `status` 和 `result` 的映射；使用其中的 `result`，`status` 不是 `ok` 时抛出错误码为 `INVALID_DESIGN_ENVELOPE` 的 `InputFormatError`），或 `design()` 内层的 `result` 映射（含 `layout` 的那个），需有 `treatment` 列（`options['treatment_column']`），可有 `block`，或 `row` 与 `column`。`options`：`treatment_column`（`treatment`）、`factor_columns`（数值列名列表，如 -1/+1 编码）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `n_plots`、`n_treatments`、`replications`、`balanced`；有 `block` 列时另有 `n_blocks` 与 `complete_blocks`，有 `row` 与 `column` 时另有 `latin_square`，有 `factor_columns` 时另有 `column_dot_products` 与 `orthogonal`；最后是 `issues`（消息列表）与 `ok`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import check, design

layout = design({"design_type": "rcbd", "treatments": ["A", "B", "C"], "replicates": 3, "seed": 1})
result = check(layout)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "experimental_design",
  "operation": "check",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 82,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50d77680785df2e143b723bc829268992fcc49786029d1d325a320e054846b10"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_plots": 9,
    "n_treatments": 3,
    "replications": {
      "B": 3,
      "C": 3,
      "A": 3
    },
    "balanced": true,
    "n_blocks": 3,
    "complete_blocks": true,
    "issues": [],
    "ok": true
  }
}
```

## 4) design

- **作用：** 生成随机化的试验布局：完全随机、随机完全区组、拉丁方或析因设计。

- **处理方法：** `crd`：每个处理重复 `replicates` 次，并将单位打乱后分配到小区。`rcbd`：每个区组含每个处理一次，顺序在区组内重新随机（区组数 = `blocks`，默认等于 `replicates`）。`latin_square`：对循环拉丁方的行、列和处理标签做随机置换。`factorial`：处理为各因子水平的全部组合（标签形如 `A=1|B=x`），按 CRD 布局，给出 `blocks` 时按 RCBD 布局，每个因子一列。打乱使用 Python `random.Random(seed)`；超过 100,000 个小区的布局会被拒绝。

- **API：** `dnakit.experimental_design.design(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射（与 `options` 合并）。`design_type` 为 `crd`（默认）、`rcbd`、`latin_square` 或 `factorial`；`treatments`（列表）或 `n_treatments`；`replicates`（2；CRD 为重复数，RCBD 在未给 `blocks` 时为区组数）；析因设计用 `factors`（{名称: 水平}）；`seed`（0）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `design_type`、`seed`、`n_plots`、`layout`（含 `plot`、`block` 或 `row`/`column`、`treatment` 的记录，析因设计另含各因子列）、`replications`（{处理: 次数}）与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import design

result = design({"design_type": "rcbd", "treatments": ["A", "B", "C"], "replicates": 2, "seed": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "experimental_design",
  "operation": "design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 82,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "37fea3b5d78613b96f698d7b9dbf24c5dd634e12891148cf0671e29b399d67a1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "design_type": "rcbd",
    "seed": 1,
    "n_plots": 6,
    "layout": [
      {
        "plot": 1,
        "block": 1,
        "treatment": "B"
      },
      {
        "plot": 2,
        "block": 1,
        "treatment": "C"
      },
      "..."
    ],
    "replications": {
      "B": 2,
      "C": 2,
      "A": 2
    },
    "definition": "each block holds every treatment once in random order (randomised complete block)"
  }
}
```

## 5) export

- **作用：** 把设计布局导出为 CSV、TSV 或 JSON 文本，并可写入文件。

- **处理方法：** 按所选列（默认取第一条记录的全部列）序列化设计记录；CSV/TSV 中的缺失值写为空，JSON 保留 `null`。仅在给出 `path` 时写文件，已存在的文件不会被覆盖，除非 `overwrite=True`。文本总会返回。

- **API：** `dnakit.experimental_design.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：设计记录、`design()` 返回的完整封装（含 `status` 和 `result` 的映射；使用其中的 `result`，`status` 不是 `ok` 时抛出错误码为 `INVALID_DESIGN_ENVELOPE` 的 `InputFormatError`），或 `design()` 内层的 `result` 映射（含 `layout` 的那个）。`options`：`format` 为 `csv`（默认）、`tsv` 或 `json`；`columns`（子集与顺序）；`path`（或 `output`）；`overwrite`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `format`、`n_rows`、`columns`、`text`（序列化后的表格）与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import design, export

layout = design({"design_type": "crd", "treatments": ["A", "B"], "replicates": 2, "seed": 1})
result = export(layout, options={"format": "csv", "columns": ["plot", "treatment"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "experimental_design",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d8dbc7378da676ba95efe489f195176d930d39dd9ba2373fd6e2742835265bd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "csv",
    "n_rows": 4,
    "columns": [
      "plot",
      "treatment"
    ],
    "text": "plot,treatment\n1,B\n2,A\n3,B\n4,A\n",
    "path": null
  }
}
```

## 6) generate

- **作用：** 生成标准顺序的设计矩阵：全因子、二水平部分因子或循环拉丁方。

- **处理方法：** `full_factorial`：按标准（Yates）顺序列出所有水平组合，第一个因子变化最快（至多 100,000 次运行）。`fractional_factorial`：`n_factors`（2 至 12，命名为 A、B、...）个因子的 2^(k-p) 设计，每个生成元（如 `D=ABC`）把生成列定义为基础列的乘积，编码为 -1/+1；基础因子按标准顺序排列。`latin_square`：循环方阵 t[(r+c) mod k]，不做随机化（需要随机化请用 `design`）。所生成的矩阵均未随机化。

- **API：** `dnakit.experimental_design.generate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射（与 `options` 合并）。`kind` 为 `full_factorial`（默认；`factors` = {名称: 水平}）、`fractional_factorial`（`n_factors` 与 `generators`，如 `['D=ABC']`，只能用基础因子表示）或 `latin_square`（`treatments`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `kind`；析因类另含 `n_runs`、`runs`（含 `run` 与每个因子一列的记录）、`order`，部分因子另含 `generators`；`latin_square` 则含 `square`（行列表）与 `order`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import generate

result = generate(
    {"kind": "fractional_factorial", "n_factors": 4, "generators": ["D=ABC"]}
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "experimental_design",
  "operation": "generate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 73,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "feb488336a79b02b2e5ed7e00157348fe2840f9536dcde128d34618844093248"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "fractional_factorial",
    "n_runs": 8,
    "runs": [
      {
        "run": 1,
        "A": -1,
        "B": -1,
        "C": -1,
        "D": -1
      },
      {
        "run": 2,
        "A": 1,
        "B": -1,
        "C": -1,
        "D": 1
      },
      "..."
    ],
    "generators": [
      "D=ABC"
    ],
    "order": "standard order of the base factors (first varies fastest); generated columns are products of base columns"
  }
}
```

## 7) randomize

- **作用：** 随机化一组运行的顺序（可在组内进行），并加上从 1 开始的 `run_order`。

- **处理方法：** 使用 Python `random.Random(seed)` 的 Fisher-Yates 洗牌。给出 `within` 时，在该列的每个组内分别洗牌，各组保持首次出现的顺序。接受 1 至 1,000,000 次运行。

- **API：** `dnakit.experimental_design.randomize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：条目列表（每项变为 `{'item': x}`）或记录。`options`：`seed`（0）、`within`（列名）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `seed`、`n_runs`、`within`、`runs`（按新顺序排列的记录，各带 `run_order`）与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import randomize

runs = [{"sample": "s1", "batch": "b1"}, {"sample": "s2", "batch": "b1"},
        {"sample": "s3", "batch": "b2"}, {"sample": "s4", "batch": "b2"}]
result = randomize(runs, options={"seed": 7, "within": "batch"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "experimental_design",
  "operation": "randomize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 132,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "455cf2af29029a45016dc21756ea96ad764597784664f0e65a9305c082fe8872"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 7,
    "n_runs": 4,
    "within": "batch",
    "runs": [
      {
        "sample": "s1",
        "batch": "b1",
        "run_order": 1
      },
      {
        "sample": "s2",
        "batch": "b1",
        "run_order": 2
      },
      "..."
    ],
    "method": "Fisher-Yates shuffle (Python random.Random(seed)); with `within`, runs are shuffled inside each group and groups keep first-appearance order"
  }
}
```

## 8) simulate

- **作用：** 通过蒙特卡洛模拟 t 检验实验来估计检验效能，并与精确效能对比。

- **处理方法：** 生成 `n_simulations` 份单位标准差的正态数据（第 1 组均值 0，第 2 组均值 = `effect_size`；或单样本均值 `effect_size`），做双侧合并方差 t 检验（`one_sample_t` 则为单样本 t 检验），以 p < alpha 的比例作为经验效能并给出其二项标准误，同时给出精确的非中心 t 效能。给定 `seed`（Python `random`）时结果确定。

- **API：** `dnakit.experimental_design.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射（与 `options` 合并），含 `n`（每组样本数）、`effect_size`（Cohen's d）、`test`（`two_sample_t`（默认）或 `one_sample_t`）、`alpha`（0.05）、`n_simulations`（2000）、`seed`（0）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `test`、`n`、`effect_size`、`alpha`、`n_simulations`、`seed`、`empirical_power`、`standard_error`、`analytic_power` 与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import simulate

result = simulate({"n": 20, "effect_size": 0.8, "n_simulations": 500, "seed": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "experimental_design",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 62,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dd66b43353d0d0b3ff488fa6d81c75710a1925dd884f65c2e3d613f7f6e6beca"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "two_sample_t",
    "n": 20,
    "effect_size": 0.8,
    "alpha": 0.05,
    "n_simulations": 500,
    "seed": 1,
    "empirical_power": 0.704,
    "standard_error": 0.020414896521902825,
    "analytic_power": 0.6934041966153844,
    "method": "Monte Carlo: normal data with unit SD, two-sided pooled t-test; compare with the exact non-central-t power"
  }
}
```

## 9) validate

- **作用：** 校验设计表的结构：必需列存在、无空单元格、小区编号唯一、各行一致。

- **处理方法：** 检查表非空、每个必需列存在且无空值、单位列（若存在）无重复、各行列数相同。内容问题收集在 `errors` 中并使 `valid` 为 false，不会抛出异常。

- **API：** `dnakit.experimental_design.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：设计记录、`design()` 返回的完整封装（含 `status` 和 `result` 的映射；使用其中的 `result`，`status` 不是 `ok` 时抛出错误码为 `INVALID_DESIGN_ENVELOPE` 的 `InputFormatError`），或 `design()` 内层的 `result` 映射（含 `layout` 的那个）。`options`：`required_columns`（默认 `['treatment']`）、`unit_column`（默认 `plot`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `valid`、`n_rows` 与 `errors`（消息列表）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.experimental_design import validate

rows = [{"plot": 1, "treatment": "A"}, {"plot": 2, "treatment": ""}, {"plot": 2, "treatment": "B"}]
result = validate(rows)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "experimental_design",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 92,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ad0afe16f1236e53b11be4318dc7efee802a144ac7090c559777a3a886c87510"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_rows": 3,
    "errors": [
      "column 'treatment' has 1 empty value(s)",
      "unit column 'plot' has duplicated values"
    ]
  }
}
```
