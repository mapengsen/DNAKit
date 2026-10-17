# 法医与生殖基因组学

本类收录 **9** 项有文档的功能。其中 **7** 项使用标准库，基于等位基因频率和单一来源 STR/SNP 基因型谱原生计算：`calculate`（随机匹配概率、似然比、父权指数、亲缘指数）、`compare`、`infer`（亲缘关系后验）、`plot`、`read`、`report` 与 `validate`。**2** 项在任何调用下都返回 `status: "not_implemented"`，因为需要经过验证的外部模型或软件：`deconvolute`（混合物概率基因分型）与 `predict`（表型和祖先预测）。本类没有突变、drop-in/dropout 或混合物模型。这些函数通过 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) calculate

- **作用：** 由等位基因频率计算法医统计量：随机匹配概率、似然比、父权指数或亲缘指数。

- **处理方法：** `random_match_probability`：对位点连乘的乘积法则，纯合子 p^2，杂合子 2pq；`theta` > 0 时采用 NRC II 建议 4.10 公式 [2θ+(1-θ)p][3θ+(1-θ)p]/[(1+θ)(1+2θ)] 与 2[θ+(1-θ)pi][θ+(1-θ)pj]/[(1+θ)(1+2θ)]。`likelihood_ratio`：Hp 为嫌疑人是来源，Hd 为无关随机个体是来源；匹配位点的 LR = 1/P(基因型)，不匹配时为 0，对位点连乘（支持 `theta`）。`paternity_index`：PI = P(子代 | 母亲, 待检父亲) / P(子代 | 母亲, 随机男性)，每个亲本以 1/2 概率传递任一等位基因，累计父权指数 CPI 为各 PI 之积，父权概率 W = CPI*prior/(CPI*prior + 1 - prior)；二联体（无母亲）或三联体；无突变模型，出现一次排除则 CPI = 0；theta = 0。`kinship_index`：LR = P(g2 | g1, 亲缘关系)/P(g2)，P(g2 | g1) = k0 P(g2) + k1 P(g2 | 一个 IBD 等位基因) + k2 [g2 = g1]，其中 (k0, k1, k2) 为该亲缘关系的 Cotterman IBD 概率（theta = 0，无突变）。文字等价说明采用 ENFSI 2015 量表。仅依赖标准库。

- **API：** `dnakit.forensic_reproductive.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `options['quantity']`：`random_match_probability`（默认；输入 `{'profile', 'frequencies'}`；`theta` 取 0 至 0.5）、`likelihood_ratio`（`{'evidence', 'suspect', 'frequencies'}`；`theta`）、`paternity_index`（`{'child', 'alleged_father', 'mother'（可选）, 'frequencies'}`；`prior` 默认 0.5）或 `kinship_index`（`{'profile_1', 'profile_2', 'frequencies'}`；`relationship` 为 `parent_child`、`full_sibling`（默认）、`half_sibling`、`grandparent_grandchild`、`uncle_nephew`、`first_cousin`、`identical_twin`，或自定义 `ibd=[k0, k1, k2]`，和为 1）。`frequencies` 也可放在 `options` 中。频率表中缺少的等位基因会报错，除非设置了 `min_frequency`。基因型谱为 `{位点: [等位基因1, 等位基因2]}`（只列一个等位基因时按纯合处理）；等位基因频率为 `{位点: {等位基因: 频率}}`，或含 `locus`、`allele`、`frequency` 列的记录。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 始终含 `quantity`、`per_locus`（每个位点一条记录）、`combined`、`combined_log10` 与 `formula`。随机匹配概率：`theta`，各位点的 `genotype_probability` 与 `log10`，以及 `one_in`。似然比：`theta`，各位点的 `match`、`lr`、`log10_lr`，以及 `excluded`、`verbal_equivalent`。父权指数：`design`（`duo` 或 `trio`），各位点的 `numerator`、`denominator`、`pi`、`log10_pi`，以及 `exclusions`、`prior`、`probability_of_paternity`、`verbal_equivalent`。亲缘指数：`relationship`、`ibd`，各位点的 `lr`、`log10_lr`，以及 `verbal_equivalent`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import calculate

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = calculate(
    {"profile": profile, "frequencies": frequencies},
    options={"quantity": "random_match_probability"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 264,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "604f5e370a70a7f4e8a5dad7a9412e6c919b171c7c709b6d9e46d955bf60fed1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "random_match_probability",
    "theta": 0.0,
    "per_locus": [
      {
        "locus": "D3S1358",
        "genotype": [
          "15",
          "16"
        ],
        "genotype_probability": 0.15,
        "log10": -0.8239087409443188
      },
      {
        "locus": "vWA",
        "genotype": [
          "17",
          "17"
        ],
        "genotype_probability": 0.12249999999999998,
        "log10": -0.9118639112994488
      },
      {
        "locus": "FGA",
        "genotype": [
          "22",
          "24"
        ],
        "genotype_probability": 0.12,
        "log10": -0.9208187539523752
      }
    ],
    "combined": 0.002205,
    "combined_log10": -2.6565914061961426,
    "one_in": 453.51473922902494,
    "formula": "product rule: homozygote p^2, heterozygote 2pq; with theta>0 the NRC II 4.10 formulas [2θ+(1-θ)p][3θ+(1-θ)p]/[(1+θ)(1+2θ)] and 2[θ+(1-θ)pi][θ+(1-θ)pj]/[(1+θ)(1+2θ)]"
  }
}
```

## 2) compare

- **作用：** 逐位点比较两个单一来源的基因型谱（完全相同、共享一个等位基因、无共享等位基因）。

- **处理方法：** 对两个谱共有的每个位点，把基因型作为无序对比较：完全相同（共享 2 个等位基因）、共享一个等位基因（按重数计）或无共享。所有被比较位点都相同时 `consistent_with_same_source` 为真；除非给出 `frequencies`（此时另加谱 1 在被比较位点上的随机匹配概率），否则不计算统计权重。

- **API：** `dnakit.forensic_reproductive.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'profile_1': {...}, 'profile_2': {...}}`，可选 `'frequencies'`（或 `options['frequencies']`）。基因型谱为 `{位点: [等位基因1, 等位基因2]}`（只列一个等位基因时按纯合处理）；等位基因频率为 `{位点: {等位基因: 频率}}`，或含 `locus`、`allele`、`frequency` 列的记录。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `per_locus`（`locus`、两个基因型、`status`、`alleles_shared`）、`n_loci_compared`、计数 `identical`、`one_allele_shared`、`no_allele_shared`、`loci_only_in_profile_1`、`loci_only_in_profile_2`、`consistent_with_same_source`、`mean_alleles_shared`，给出频率时另有 `random_match_probability_profile_1`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import compare

profile_1 = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
profile_2 = {"D3S1358": ["15", "16"], "vWA": ["17", "18"], "FGA": ["21", "23"]}
result = compare({"profile_1": profile_1, "profile_2": profile_2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 164,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9a2cd44364446204b1f0719418ea04084042426d9b10b0dcd968f5421393ffe4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "per_locus": [
      {
        "locus": "D3S1358",
        "profile_1": [
          "15",
          "16"
        ],
        "profile_2": [
          "15",
          "16"
        ],
        "status": "identical",
        "alleles_shared": 2
      },
      {
        "locus": "vWA",
        "profile_1": [
          "17",
          "17"
        ],
        "profile_2": [
          "17",
          "18"
        ],
        "status": "one_allele_shared",
        "alleles_shared": 1
      },
      {
        "locus": "FGA",
        "profile_1": [
          "22",
          "24"
        ],
        "profile_2": [
          "21",
          "23"
        ],
        "status": "no_allele_shared",
        "alleles_shared": 0
      }
    ],
    "n_loci_compared": 3,
    "identical": 1,
    "one_allele_shared": 1,
    "no_allele_shared": 1,
    "loci_only_in_profile_1": [],
    "loci_only_in_profile_2": [],
    "consistent_with_same_source": false,
    "mean_alleles_shared": 1.0
  }
}
```

## 3) deconvolute

- **作用：** 把 DNA 混合物解卷积为各贡献者的基因型。DNAKit 未实现该功能。

- **处理方法：** 不做任何计算。混合物解卷积需要经过验证的、带峰高以及 drop-in/dropout 与 stutter 参数的概率基因分型模型（如 EuroForMix、STRmix）；本函数只报告这一需求。

- **API：** `dnakit.forensic_reproductive.deconvolute(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意 `input` 与可选 `options`；它们仅在 `input_summary` 中被概括。

- **输出：** 不计算结果。正常调用返回 `status: "not_implemented"` 的封装（`provenance.implementation` 为 `not_implemented`），并附 `required_dependency` 与 `message` 说明缺少什么：需要经过验证的概率基因分型软件（如 EuroForMix、STRmix）以及峰高数据。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import deconvolute

mixture = {"D3S1358": ["15", "16", "17"], "vWA": ["16", "17", "18"]}
result = deconvolute({"mixture": mixture})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "deconvolute",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 71,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "249f54c134c4b9460798aba692cf46ed0ed2bfc4fff3969693bb051e5987f4ef"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated probabilistic genotyping software (e.g. EuroForMix, STRmix) and peak-height data",
  "message": "Mixture deconvolution needs a validated probabilistic genotyping model with peak heights, drop-in/dropout and stutter parameters; DNAKit does not implement it."
}
```

## 4) infer

- **作用：** 根据 IBD 模型的似然比推断两个基因型谱之间最可能的亲缘关系。

- **处理方法：** 对每个候选亲缘关系，由 Cotterman IBD 概率 (k0, k1, k2) 计算相对于 `unrelated` 的亲缘似然比（theta = 0，无突变），后验概率 posterior(r) 正比于 prior(r) * LR(r 对 unrelated)，并在候选间归一化；先验默认均匀。具有相同 (k0, k1, k2) 的关系（half_sibling、grandparent_grandchild、uncle_nephew）无法被该模型区分：它们的似然比始终相等（后验概率只因先验而不同），其中的 `most_probable` 反映的是列表顺序而不是证据。只要所评估的关系中有两个或更多共享相同的 (k0, k1, k2)，`warnings` 中就会有一条消息列出它们。若谱与所有候选关系都不一致则报错。

- **API：** `dnakit.forensic_reproductive.infer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'profile_1', 'profile_2', 'frequencies'}`。`options`：`relationships`（默认 `unrelated`、`parent_child`、`full_sibling`、`half_sibling`、`first_cousin`；可用 `calculate` 接受的任意名称）、`priors`（{亲缘关系: 先验}，每个所列关系一项）、`min_frequency`。基因型谱为 `{位点: [等位基因1, 等位基因2]}`（只列一个等位基因时按纯合处理）；等位基因频率为 `{位点: {等位基因: 频率}}`，或含 `locus`、`allele`、`frequency` 列的记录。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `relationships`（按后验概率排序；每项含 `relationship`、`likelihood_ratio_vs_unrelated`、`log10_lr`、`prior`、`posterior`）、`most_probable` 与 `method`；在所评估的关系中，每组数学上无法区分的关系在 `warnings` 中有一条消息。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import infer

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile_1 = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
profile_2 = {"D3S1358": ["15", "17"], "vWA": ["17", "18"], "FGA": ["24", "23"]}
result = infer(
    {"profile_1": profile_1, "profile_2": profile_2, "frequencies": frequencies},
    options={"relationships": ["unrelated", "parent_child", "full_sibling"]},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "infer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 348,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3adbb8c661a53b053f3e78105cc3481951e6d5846c902b6a1ba2d867d3e90554"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "relationships": [
      {
        "relationship": "parent_child",
        "likelihood_ratio_vs_unrelated": 1.7857142857142858,
        "log10_lr": 0.2518119729937996,
        "prior": 0.3333333333333333,
        "posterior": 0.5223636957231472
      },
      {
        "relationship": "unrelated",
        "likelihood_ratio_vs_unrelated": 1.0,
        "log10_lr": 0.0,
        "prior": 0.3333333333333333,
        "posterior": 0.29252366960496246
      },
      {
        "relationship": "full_sibling",
        "likelihood_ratio_vs_unrelated": 0.6328125000000001,
        "log10_lr": -0.19872495076921853,
        "prior": 0.3333333333333333,
        "posterior": 0.18511263467189032
      }
    ],
    "most_probable": "parent_child",
    "method": "posterior(r) ∝ prior(r) * LR(r vs unrelated), LR from the Cotterman IBD model (theta = 0, no mutation); priors default to uniform"
  }
}
```

## 5) plot

- **作用：** 把法医统计量的各位点 log10 基因型概率或似然比绘制为 SVG 柱状图。

- **处理方法：** 按 `calculate` 的方式计算统计量，在固定 720 x 400 的 SVG 上为每个位点画一根柱（`random_match_probability` 为 log10 基因型概率，其余为 log10 LR / PI），带零线、刻度和旋转的位点名。LR 为 0（排除）的位点不画柱。未给出 `path` 时不写文件。

- **API：** `dnakit.forensic_reproductive.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 与 `options` 同 `calculate`（`quantity` 选择统计量）；`path` 将 SVG 写入文件（已有文件不会被覆盖，除非 `overwrite=True`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `quantity`、`svg`（SVG 文本）、`n_loci` 与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import plot

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = plot(
    {"profile": profile, "frequencies": frequencies},
    options={"quantity": "random_match_probability"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "forensic_reproductive",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 264,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "604f5e370a70a7f4e8a5dad7a9412e6c919b171c7c709b6d9e46d955bf60fed1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "random_match_probability",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"400\" viewBox=\"0 0 720 400\"><title>random_match_probability p...",
    "n_loci": 3,
    "path": null
  }
}
```

## 6) predict

- **作用：** 由 DNA 谱预测表型或祖先成分。DNAKit 未实现该功能。

- **处理方法：** 不做任何计算。法医 DNA 表型与祖先预测需要经过验证的训练模型和参考群体（如 HIrisPlex-S 或参考面板），本库未内置；本函数只报告这一需求。

- **API：** `dnakit.forensic_reproductive.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意 `input` 与可选 `options`；它们仅在 `input_summary` 中被概括。

- **输出：** 不计算结果。正常调用返回 `status: "not_implemented"` 的封装（`provenance.implementation` 为 `not_implemented`），并附 `required_dependency` 与 `message` 说明缺少什么：需要经过验证的法医表型/祖先预测模型和参考数据（如 HIrisPlex-S、参考面板）。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import predict

result = predict({"profile": {"D3S1358": ["15", "16"], "vWA": ["17", "17"]}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "36756aa017f6c6033830b17c63bfe422f26f28d09aaa87be8526a8ab01acb794"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated forensic phenotype/ancestry prediction models and reference data (e.g. HIrisPlex-S, a reference panel)",
  "message": "Forensic DNA phenotyping and ancestry prediction require validated trained models and reference populations that are not bundled."
}
```

## 7) read

- **作用：** 把基因型谱或等位基因频率表（路径、分隔文本、记录或映射）读入 DNAKit 数据结构。

- **处理方法：** 分隔文本、文件、记录或列映射先作为表格加载。`kind='auto'` 时，含 `frequency`/`freq` 列的表为频率表，否则为基因型谱。基因型谱的列为 `locus`、`allele1`、`allele2`；频率表的列为 `locus`、`allele`、`frequency`。含 `profile` 和/或 `frequencies` 键的映射直接规范化（等位基因标签转为字符串）。仅依赖标准库。

- **API：** `dnakit.forensic_reproductive.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：路径、分隔文本、记录，或含 `profile` / `frequencies` 的映射。`options`：`kind` 为 `auto`（默认）、`profile` 或 `frequencies`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `kind`、`profile`（{位点: [等位基因1, 等位基因2]}）或 `frequencies`（{位点: {等位基因: 频率}}）以及 `n_loci`；映射输入则只返回规范化后的 `profile` 和/或 `frequencies`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import read

text = "locus,allele1,allele2\nD3S1358,15,16\nvWA,17,17\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "forensic_reproductive",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 46,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8491b2344053456425e5c6a64b5e7257c361b1691765ccfd9662f0422ca24a98"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "profile",
    "profile": {
      "D3S1358": [
        "15",
        "16"
      ],
      "vWA": [
        "17",
        "17"
      ]
    },
    "n_loci": 2
  }
}
```

## 8) report

- **作用：** 为法医统计量生成 Markdown 统计报告（各位点表、合并值、ENFSI 文字量表、局限性）。

- **处理方法：** 按 `calculate` 的方式计算统计量并格式化为 Markdown：标题、统计量与设计、所用公式、含 log10 的各位点数值表、合并值、以 1/N 表示的随机匹配概率或 ENFSI 2015 文字等价说明（父权指数另含父权概率与排除数），以及固定的局限性说明（单一来源谱、无突变、drop-in/dropout 或混合物模型）。它只是格式化辅助，不是经验证的案件报告。

- **API：** `dnakit.forensic_reproductive.report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 与 `options` 同 `calculate`（`quantity` 选择统计量）；`title`（默认 `Forensic statistical report`）；`path` 将 Markdown 写入文件（已有文件不会被覆盖，除非 `overwrite=True`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `format`（`markdown`）、`text`（报告正文）、`statistic`（`calculate` 的完整结果）与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import report

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
evidence = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = report(
    {"evidence": evidence, "suspect": evidence, "frequencies": frequencies},
    options={"quantity": "likelihood_ratio", "title": "Case 1"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "forensic_reproductive",
  "operation": "report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 345,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "95080364deae6a65d99264fa39751ecd91ce28574936d2073d6f7aca4effb96c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "text": "# Case 1\n\nStatistic: **likelihood ratio**\n\nMethod: Hp: the suspect is the source of the evidence; Hd: an unrelated rando...",
    "statistic": {
      "quantity": "likelihood_ratio",
      "theta": 0.0,
      "per_locus": [
        {
          "locus": "D3S1358",
          "match": true,
          "lr": 6.666666666666667,
          "log10_lr": 0.8239087409443188
        },
        {
          "locus": "vWA",
          "match": true,
          "lr": 8.163265306122451,
          "log10_lr": 0.9118639112994489
        },
        {
          "locus": "FGA",
          "match": true,
          "lr": 8.333333333333334,
          "log10_lr": 0.9208187539523752
        }
      ],
      "combined": 453.5147392290254,
      "combined_log10": 2.656591406196143,
      "excluded": false,
      "verbal_equivalent": "moderately strong support for the first (prosecution/relationship) proposition",
      "formula": "Hp: the suspect is the source of the evidence; Hd: an unrelated random person is the source. Per locus LR = 1/P(genotype) when the genotypes match, 0 otherwise; combined = product over loci"
    },
    "path": null
  }
}
```

## 9) validate

- **作用：** 校验法医基因型谱与等位基因频率表：频率之和、每个位点的等位基因数，以及频率表中缺失的位点或等位基因。

- **处理方法：** 错误：某位点的频率之和与 1 相差超过 `tolerance`；某位点有多于两个等位基因（混合物，而非单一来源谱）；某位点不在频率表中。警告：只列出一个等位基因（按纯合处理，可能有 dropout）；某等位基因不在频率表中（此时计算需要 `min_frequency`）。内容问题使 `valid` 为 false，而不是抛出异常。

- **API：** `dnakit.forensic_reproductive.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射，可含 `frequencies`、`profile`、`evidence`、`suspect`、`child` 等（每个谱类条目都会被检查）。`options`：`tolerance`（0.02）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `valid`、`n_errors`、`n_warnings` 与 `issues`（含 `severity`、`where`、`message` 的记录）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.forensic_reproductive import validate

frequencies = {"D3S1358": {"15": 0.5, "16": 0.3}, "vWA": {"17": 0.5, "18": 0.5}}
profile = {"D3S1358": ["15", "16", "17"], "vWA": ["17", "19"], "FGA": ["22", "24"]}
result = validate({"frequencies": frequencies, "profile": profile})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "forensic_reproductive",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 169,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1d22dcd2d8b8e87b5a0373b7cc09edab13b7cbd26fa968b9172c8e9b52a50855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_errors": 3,
    "n_warnings": 1,
    "issues": [
      {
        "severity": "error",
        "where": "D3S1358",
        "message": "allele frequencies sum to 0.8000, not 1 (tolerance 0.02)"
      },
      {
        "severity": "error",
        "where": "profile:D3S1358",
        "message": "3 alleles: a single-source profile has at most two"
      },
      "..."
    ]
  }
}
```
