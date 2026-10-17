# 育种分析

本类收录 **7** 项有文档的功能，全部由 DNAKit 原生计算（正常调用均不会返回 `not_implemented`）：育种者方程与遗传力计算、系谱与基因组亲缘关系、岭回归/GBLUP 拟合与预测、截断选择与 Smith-Hazel 选择、最优贡献优化、前向选择模拟，以及 SVG 进展图。`fit`、`predict` 和 `calculate(quantity='genomic_relationship')` 需要 numpy，缺失时返回 `status: "not_implemented"`。这些函数通过 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) calculate

- **作用：** 计算数量遗传学量：选择反应（育种者方程）、遗传力、系谱（加性）亲缘关系矩阵或 VanRaden 基因组亲缘关系矩阵。

- **处理方法：** 由 `options['quantity']` 选择计算量。`breeders_equation`（默认）：R = h2 * S，S = i * sigma_p；截断选择的选择强度 i = phi(z)/p，z = Phi^-1(1-p)（无限群体、正态性状）；给出 generation_interval 时另报告 `response_per_year` = R / 世代间隔。`heritability`：由亲子列表对子代做回归（对中亲值 h2 = b，对单亲 h2 = 2b），或实现遗传力 R / S，或 Vg / Vp。`pedigree_relationship`：Henderson 表格法，A_ij = (A_sj + A_dj)/2，A_ii = 1 + A_sd/2，近交系数 F = A_ii - 1。`genomic_relationship`：VanRaden 基因组关系矩阵（需要 numpy，缺失时返回 `not_implemented`）。

- **API：** `dnakit.breeding.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 与 `options`。`breeders_equation`：映射，含 `heritability`（0 到 1），以及 `selection_differential`，或 `selection_intensity` / `proportion_selected` 加 `phenotypic_sd`；可选 `generation_interval`（这些值也可放在 `options` 中）。`heritability`：含 `parent` 与 `offspring` 列表的映射（至少 3 对；`options['relationship']` 为 `midparent`（默认）或 `single_parent`），或 `response` 加 `selection_differential`，或 `genetic_variance` 加 `phenotypic_variance`。`pedigree_relationship`：含 `id`、`sire`、`dam` 的记录（或含 `pedigree` 的映射）。`genomic_relationship`：基因型矩阵（样本 x 标记的剂量）。`options['quantity']` 默认为 `breeders_equation`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`quantity` 为 `breeders_equation` 时，`result` 含 `quantity`、`heritability`、`selection_differential`、`selection_intensity`（无法推得时为 null）、`response`、`response_per_year`（仅在给出世代间隔时）与 `definition`。其他量分别返回 `h2`、`slope`、`slope_se`、`n_pairs`、`method`（遗传力），`ids`、`matrix` 与 `inbreeding`（{id: F}）（系谱），或 `sample_ids`、`relationship_matrix`、`method`、`n_markers_used`、`allele_frequencies` 及标记过滤信息（基因组）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import calculate

result = calculate(
    {"heritability": 0.5, "selection_differential": 2.0},
    options={"quantity": "breeders_equation"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "breeding",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 52,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d78fc459837b52d2ed27a19a773b13dd1a63f39d2106ec21454b9fa6df84268"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "breeders_equation",
    "heritability": 0.5,
    "selection_differential": 2.0,
    "selection_intensity": null,
    "response": 1.0,
    "definition": "R = h2 * S, S = i * sigma_p; truncation intensity i = phi(z)/p with z = Phi^-1(1-p) (infinite population, normal trait)"
  }
}
```

## 2) fit

- **作用：** 用标记剂量对表型拟合基因组岭回归模型（SNP-BLUP，等价于 GBLUP）。

- **处理方法：** 对中心化的剂量做岭回归，最小化 sum((y - mu - Zb)^2) + lambda * sum(b^2)；当 lambda = 2*sum(pq)*(1-h2)/h2 时等价于 G = ZZ'/(2*sum(pq)) 的 GBLUP，p 由训练样本估计。育种值为拟合值减去截距。训练相关系数在训练数据上计算，因而偏乐观（会给出警告）。需要 numpy，缺失时返回 `not_implemented`。

- **API：** `dnakit.breeding.fit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'genotypes': 样本 x 标记的剂量, 'phenotypes': [...]}`（表型为 null 的样本被丢弃）。`options`：`ridge_lambda`（对中心化剂量的惩罚，> 0），或 `h2`（0 到 1，此时 lambda = 2*sum(p*q)*(1-h2)/h2）；二者必须给出其一。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `model`（`rrblup_ridge`）、`ridge_lambda` 与 `ridge_lambda_source`、`intercept`、`raw_intercept`、`marker_effects`、`allele_frequencies`、`n_markers_used`、`n_training_samples`、`training_predictions`、`training_breeding_values`、`training_correlation`、`sample_ids` 与 `definition`。警告说明训练相关系数高估了准确度。该结果可作为 `model` 传给 `predict`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import fit

genotypes = [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 1, 0], [1, 0, 2, 1], [2, 1, 0, 2]]
phenotypes = [1.0, 2.1, 2.9, 0.8, 3.2, 4.0]
result = fit({"genotypes": genotypes, "phenotypes": phenotypes}, options={"ridge_lambda": 1.0})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "breeding",
  "operation": "fit",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 145,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eef533a99c44baec24b7857c098f6dfa38f9dd4eb4bcdf24fec13da360ff651f"
  },
  "warnings": [
    "training_correlation is computed on the training data and overstates prediction accuracy"
  ],
  "status": "ok",
  "result": {
    "model": "rrblup_ridge",
    "ridge_lambda": 1.0,
    "ridge_lambda_source": "ridge_lambda",
    "intercept": 2.3333333333333335,
    "raw_intercept": 1.2766513056835636,
    "marker_effects": [
      0.6798771121351764,
      -0.26359447004608283,
      0.12319508448540714,
      0.47327188940092196
    ],
    "allele_frequencies": [
      0.5,
      0.4166666666666667,
      0.5,
      0.5
    ],
    "n_markers_used": 4,
    "n_training_samples": 6,
    "training_predictions": [
      1.259447004608295,
      2.6394777265745013,
      3.2328725038402455,
      0.8726574500768052,
      2.6761904761904765,
      3.3193548387096774
    ],
    "training_breeding_values": [
      -1.0738863287250384,
      0.3061443932411677,
      0.8995391705069121,
      -1.4606758832565283,
      0.34285714285714286,
      0.9860215053763441
    ],
    "training_correlation": 0.9290475470728303,
    "sample_ids": [
      "sample1",
      "sample2",
      "sample3",
      "sample4",
      "sample5",
      "sample6"
    ],
    "definition": "ridge regression of the phenotype on centred dosages (SNP-BLUP), objective sum((y-mu-Zb)^2) + lambda*sum(b^2); equivalen..."
  }
}
```

## 3) optimize

- **作用：** 优化亲本贡献（最优贡献选择）：最大化平均育种值并扣除亲缘惩罚。

- **处理方法：** 用有效集法求解凸二次规划 max c'u - (penalty/2) c'Ac，约束 sum(c) = 1、c >= 0，其中 u 为候选个体的育种值，A 为亲缘关系矩阵；在解处检验 KKT 条件。群体共祖系数按 c'Ac/2 报告。亲缘关系矩阵可直接给出，或由系谱用 Henderson 表格法计算。

- **API：** `dnakit.breeding.optimize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'candidates': 含 id 与 ebv 的记录, 'relationship': {'ids': [...], 'matrix': [[...]]}}`，或 `{'candidates': ..., 'pedigree': 含 id/sire/dam 的记录}`。`options`：`penalty`（默认 1.0）、`trait`（候选记录中存放育种值的列，默认 `ebv`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `contributions`（{id: 贡献}，仅非零项）、`mean_ebv`、`mean_relationship`（c'Ac）、`group_coancestry`（c'Ac/2）、`penalty`、`objective` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import optimize

candidates = [{"id": "A", "ebv": 3.0}, {"id": "B", "ebv": 2.5}, {"id": "C", "ebv": 1.0}]
relationship = {"ids": ["A", "B", "C"], "matrix": [[1.0, 0.8, 0.0], [0.8, 1.0, 0.0], [0.0, 0.0, 1.0]]}
result = optimize({"candidates": candidates, "relationship": relationship}, options={"penalty": 10.0})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "breeding",
  "operation": "optimize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 196,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ecac1283d9753f119385ad279d1a07411bc14b7fc9f605e69b801c825d0af4bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "contributions": {
      "A": 0.43421052631578966,
      "B": 0.18421052631578932,
      "C": 0.3815789473684211
    },
    "mean_ebv": 2.1447368421052633,
    "mean_relationship": 0.49605263157894747,
    "group_coancestry": 0.24802631578947373,
    "penalty": 10.0,
    "objective": -0.3355263157894739,
    "definition": "maximise c'u - (penalty/2) c'Ac subject to sum(c) = 1, c >= 0 (convex QP solved by an active-set method; KKT conditions verified at the solution)"
  }
}
```

## 4) plot

- **作用：** 将育种进展（每代的数值）绘制为 SVG 折线图。

- **处理方法：** 生成固定 640 x 400 的 SVG：带边框的绘图区、取整的坐标刻度、每个序列一条折线（颜色在五色调色板中循环）、序列标签，横轴按世代序号标注。纵轴范围取所有序列的最小/最大值；常数序列使用单位范围。仅依赖标准库；未给出 `path` 时不写文件。

- **API：** `dnakit.breeding.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'series': {名称: [每代的数值]}}`，或 `breeding.simulate` 的结果（含 `generations` 的映射），此时由 `options['metric']`（默认 `mean_genetic_value`）选择绘制的列。`options`：`path`（将 SVG 写入该路径；已有文件不会被覆盖，除非 `overwrite=True`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `svg`（SVG 文本）、`series`（{名称: 点数}）与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import plot

result = plot({"series": {"mean_genetic_value": [0.0, 0.8, 1.5, 2.1]}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "breeding",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a3ac612626547dbf535b0080025c610c92cd77554db08194a08ac88dfaef3fd7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"400\" viewBox=\"0 0 640 400\"><title>breeding progress</title><...",
    "series": {
      "mean_genetic_value": 4
    },
    "path": null
  }
}
```

## 5) predict

- **作用：** 用岭回归/GBLUP 模型预测新个体的基因组估计育种值。

- **处理方法：** 给定已拟合的 `model`（`fit` 的 `result`）时，预测值 = raw_intercept + sum_j b_j * 剂量_j，育种值为预测值减去模型截距；缺失剂量以训练等位基因频率的 2 倍代替。未给 `model` 时，先在 `train_genotypes` 与 `phenotypes` 上拟合与 `fit` 相同的岭回归模型（需要 numpy，缺失时返回 `not_implemented`），缺失剂量以训练均值代替。

- **API：** `dnakit.breeding.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'model': <breeding.fit 的结果>, 'genotypes': 新个体剂量}`（标记及顺序需与训练一致），或 `{'train_genotypes': ..., 'phenotypes': [...], 'genotypes': 新个体剂量}`。`options`（仅重新拟合路径）：`ridge_lambda` 或 `h2`，含义同 `fit`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `sample_ids`、`predictions` 与 `breeding_values`；重新拟合路径另含 `ridge_lambda`、`ridge_lambda_source` 与 `n_markers_used`，使用模型的路径另含 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import predict

result = predict({
    "train_genotypes": [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 1, 0], [1, 0, 2, 1], [2, 1, 0, 2]],
    "phenotypes": [1.0, 2.1, 2.9, 0.8, 3.2, 4.0],
    "genotypes": [[1, 1, 1, 1], [2, 0, 2, 0]],
}, options={"ridge_lambda": 1.0})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "breeding",
  "operation": "predict",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 194,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b87bce67e2f2c7d47c89d18b80153af40dae06f8663c05579f7b6b60bb332f59"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sample_ids": [
      "sample1",
      "sample2"
    ],
    "predictions": [
      2.2894009216589866,
      2.8827956989247308
    ],
    "breeding_values": [
      -0.04393241167434713,
      0.5494623655913973
    ],
    "ridge_lambda": 1.0,
    "ridge_lambda_source": "ridge_lambda",
    "n_markers_used": 4
  }
}
```

## 6) select

- **作用：** 按性状/EBV 截断选择，或按 Smith-Hazel 选择指数挑选育种候选个体。

- **处理方法：** `truncation` 保留性状值最优的 k 个候选。`index` 计算 Smith-Hazel 指数 I = b'x，b = P^-1 G a（P 为表型协方差矩阵，G 为遗传协方差矩阵，a 为经济权重），准确度 r_IH = sqrt(b'Pb / a'Ga)；P 奇异时报错。选择差为入选个体均值减去全部候选均值。仅依赖标准库。

- **API：** `dnakit.breeding.select(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：含 `id`（唯一）与性状列的记录。`options`：`method` 为 `truncation`（默认）或 `index`；`trait`（截断选择所用的列，默认 `ebv`）；`index` 需要 `traits`、`economic_weights`、`genetic_covariance`、`phenotypic_covariance`；`n_select` 与 `proportion` 必须恰给一个；`higher_is_better`（默认 True）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `method`、`n_select`、`selected`（入选 id）、`ranking`（每个候选的 id、score、rank）、`selection_differential`，以及 `trait`（截断）或 `index_weights`、`traits`、`index_variance`、`accuracy`（指数），和 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import select

records = [{"id": "a", "ebv": 1.2}, {"id": "b", "ebv": 3.4}, {"id": "c", "ebv": 2.2}, {"id": "d", "ebv": 0.4}]
result = select(records, options={"method": "truncation", "trait": "ebv", "n_select": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "breeding",
  "operation": "select",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "92fceebe3688fb1f2ee7418794195107fa73df082d1b93eace96de08af4e055f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "truncation",
    "n_select": 2,
    "selected": [
      "b",
      "c"
    ],
    "ranking": [
      {
        "id": "b",
        "score": 3.4,
        "rank": 1
      },
      {
        "id": "c",
        "score": 2.2,
        "rank": 2
      },
      "..."
    ],
    "selection_differential": 0.9999999999999998,
    "trait": "ebv",
    "definition": "truncation: keep the k best by the trait value; index: Smith-Hazel index I = b'x with b = P^-1 G a, accuracy r_IH = sqrt(b'Pb / a'Ga)"
  }
}
```

## 7) simulate

- **作用：** 模拟多基因性状（不连锁加性位点）的轮回截断选择，并报告每代的选择反应。

- **处理方法：** 双等位加性位点，基础频率 0.5，效应 ~ N(0,1)；表型 = 遗传值 + N(0, Ve)，Ve 由基础群体遗传方差确定，使基础遗传力等于 `heritability`。每代保留表型最高的 `selection_proportion`，随机选择亲本（不自交），配子自由重组。每代报告平均遗传值、方差、选择差、育种者方程预测（Vg/Vp）* S（针对亲本世代）以及观察到的反应（平均遗传值的变化）。给定 `seed`（Python `random`）时结果确定；这是一个小型模型，不是 AlphaSimR 这类带连锁基因组的前向模拟器。限制：个体数 10 至 50,000，位点至多 2,000，世代至多 100。

- **API：** `dnakit.breeding.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：None 或下列选项的映射（优先于 `options`）。`options`：`n_individuals`（1000）、`n_loci`（100）、`heritability`（0.5，基础群体）、`selection_proportion`（0.2）、`generations`（5）、`seed`（0）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `generations`（世代数 + 1 行：`generation`、`mean_genetic_value`、`genetic_variance`、`phenotypic_variance`、`realised_heritability_in_generation`、`selection_differential`、`predicted_response`，从第 1 代起另有 `observed_response`）、`parameters` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.breeding import simulate

result = simulate(None, options={"n_individuals": 200, "n_loci": 50, "generations": 3, "seed": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "breeding",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "generations": [
      {
        "generation": 0,
        "mean_genetic_value": 5.112519418461102,
        "genetic_variance": 18.871869146667038,
        "phenotypic_variance": 37.95812172495276,
        "realised_heritability_in_generation": 0.49717605321501257,
        "selection_differential": 8.570761111862272,
        "predicted_response": 4.261177182644397
      },
      {
        "generation": 1,
        "mean_genetic_value": 8.822494073920149,
        "genetic_variance": 18.800333738727726,
        "phenotypic_variance": 43.15957790675255,
        "realised_heritability_in_generation": 0.43560050052728416,
        "selection_differential": 9.147702786858458,
        "predicted_response": 3.9847439126303765,
        "observed_response": 3.709974655459047
      },
      "..."
    ],
    "parameters": {
      "n_individuals": 200,
      "n_loci": 50,
      "heritability": 0.5,
      "selection_proportion": 0.2,
      "environmental_variance": 18.871869146667038,
      "seed": 1
    },
    "definition": "unlinked additive biallelic loci (base frequency 0.5, effects ~ N(0,1)); phenotype = genetic value + N(0, Ve) with Ve fi..."
  }
}
```
