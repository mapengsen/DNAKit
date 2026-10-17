# 群体遗传学

本类收录 **8** 项有文档的功能（`dnakit.population` 同时导出的可选 dadi/tskit 后端封装不在本页介绍）。其中 **7** 项原生计算结果：`calculate`（等位基因频率、杂合度、成对 LD）、`compare`（Fst）、`estimate`（GRM、亲缘系数、近交系数、PCA）、`export`、`plot`、`simulate`（Balding-Nichols 与 Wright-Fisher）和 `test`（Hardy-Weinberg、群体分化）；只有 `estimate`（grm、kinship、pca）和 `plot(kind='pca')` 需要 numpy，缺失时返回 `status: "not_implemented"`。**1** 项 `infer`（祖先与结构推断）在任何调用下都返回 `status: "not_implemented"`，因为它需要参考面板或模型拟合程序。这些函数通过 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) calculate

- **作用：** 由基因型矩阵计算逐变异的等位基因频率、杂合度或成对连锁不平衡。

- **处理方法：** `allele_frequency`：替代等位基因频率 = 剂量之和/(2 x 已分型样本数)，MAF = min(p, 1-p)，并给出已分型数和缺失率，报告全体，给出 `populations` 时另报各群体；允许分数（填补）剂量。`heterozygosity`（仅限硬判读）：Ho = 已分型样本中剂量为 1 的比例，He = 由样本频率得到的 2pq，无偏 He = 2pq x 2n/(2n-1)，F = 1 - Ho/He。`ld`：对所有变异对计算复合（基因型剂量）LD，r^2 为剂量 Pearson 相关系数的平方，D = cov(剂量_a, 剂量_b)/2，D' = |D|/Dmax，其中 D > 0 时 Dmax = min(pA qB, qA pB)，否则 Dmax = min(pA pB, qA qB)；完整样本少于 3 个或变异单态的变异对被跳过。仅依赖标准库。

- **API：** `dnakit.population.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 `options`：`statistic` 为 `allele_frequency`（默认）、`heterozygosity` 或 `ld`；`max_variants`（100，最大 300；用于 `ld`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `statistic`；`allele_frequency` 另含 `overall`（每个变异：`variant_id`、`n_called`、`alt_allele_frequency`、`maf`、`missing_rate`）与 `by_population`；`heterozygosity` 另含 `variants`（`n_called`、`alt_allele_frequency`、`observed_heterozygosity`、`expected_heterozygosity`、`expected_heterozygosity_unbiased`、`f_statistic`）与 `definition`；`ld` 另含 `pairs`（`variant_1`、`variant_2`、`n`、`r2`、`d_composite`、`d_prime`）与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import calculate

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = calculate(genotypes, options={"statistic": "allele_frequency"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "population",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "allele_frequency",
    "overall": [
      {
        "variant_id": "v1",
        "n_called": 6,
        "alt_allele_frequency": 0.5,
        "maf": 0.5,
        "missing_rate": 0.0
      },
      {
        "variant_id": "v2",
        "n_called": 6,
        "alt_allele_frequency": 0.5833333333333334,
        "maf": 0.41666666666666663,
        "missing_rate": 0.0
      },
      "..."
    ],
    "by_population": {
      "A": [
        {
          "variant_id": "v1",
          "n_called": 3,
          "alt_allele_frequency": 0.16666666666666666,
          "maf": 0.16666666666666666,
          "missing_rate": 0.0
        },
        {
          "variant_id": "v2",
          "n_called": 3,
          "alt_allele_frequency": 0.3333333333333333,
          "maf": 0.3333333333333333,
          "missing_rate": 0.0
        },
        "..."
      ],
      "B": [
        {
          "variant_id": "v1",
          "n_called": 3,
          "alt_allele_frequency": 0.8333333333333334,
          "maf": 0.16666666666666663,
          "missing_rate": 0.0
        },
        {
          "variant_id": "v2",
          "n_called": 3,
          "alt_allele_frequency": 0.8333333333333334,
          "maf": 0.16666666666666663,
          "missing_rate": 0.0
        },
        "..."
      ]
    }
  }
}
```

## 2) compare

- **作用：** 用 Fst（Weir & Cockerham 1984 或 Hudson）比较群体，给出全局、逐变异和两两结果。

- **处理方法：** `weir_cockerham`（默认）：对每个双等位变异计算 Weir & Cockerham (1984) 的二倍体方差分量 a、b、c，theta = sum(a)/sum(a+b+c)（变异间的比值之和），对所有群体联合以及每对群体分别计算。`hudson`：Bhatia 等 (2013) 的两两估计量，即各变异 [(p1-p2)^2 - p1q1/(n1-1) - p2q2/(n2-1)] / [p1 q2 + p2 q1] 的求和比，n 为单倍体样本量；仅恰有两个群体时给出逐变异表。估计量未定义的变异记为 null 并被排除出比值。需要硬判读（0/1/2）和至少两个群体。

- **API：** `dnakit.population.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 必须给出 `populations` 列表（每个样本一个标签）。`options`：`method` 为 `weir_cockerham`（默认）或 `hudson`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `method`、`populations`（{名称: 样本数}）、`fst` 与 `per_variant`（{variant_id: Fst}；Hudson 仅在两个群体时给出）、`pairwise`（{`popA|popB`: Fst}）与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import compare

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = compare(genotypes, options={"method": "weir_cockerham"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "population",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "weir_cockerham",
    "populations": {
      "A": 3,
      "B": 3
    },
    "fst": 0.459016393442623,
    "per_variant": {
      "v1": 0.5384615384615384,
      "v2": 0.31818181818181834,
      "v3": 0.5384615384615384,
      "v4": 0.5384615384615384,
      "v5": 0.3181818181818181
    },
    "definition": "Weir & Cockerham (1984) theta = sum(a) / sum(a + b + c) over variants (ratio of sums), diploid genotype-based components",
    "pairwise": {
      "A|B": 0.459016393442623
    }
  }
}
```

## 3) estimate

- **作用：** 估计基因组亲缘关系（VanRaden GRM 或亲缘系数）、逐样本近交系数 F，或基因型主成分。

- **处理方法：** `grm`（默认）：VanRaden 方法 1 的基因组关系矩阵。`kinship`：同一矩阵除以 2（共祖系数）。`pca`：对频率标准化剂量 X = (剂量 - 2p)/sqrt(2p(1-p)) 做薄 SVD，得分为 U*S，协方差特征值为 S^2/(n-1)；全缺失或不变的标记被剔除，缺失剂量以标记均值填补。这三者使用 `dnakit.genomics` 中基于 numpy 的辅助函数，缺少 numpy 时返回 `not_implemented`。`inbreeding`：PLINK `--het` 的矩估计，逐样本 F = (O_hom - E_hom)/(N - E_hom)，E_hom = 对已分型标记求和 1 - 2pq x 2n/(2n-1)（硬判读，标准库）。

- **API：** `dnakit.population.estimate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 `options`：`quantity` 为 `grm`（默认）、`kinship`、`inbreeding` 或 `pca`；`components`（2；用于 `pca`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `quantity`；`grm` 与 `kinship` 另含 `sample_ids`、`relationship_matrix`、`method`、`n_markers_used`、`allele_frequencies` 及标记过滤信息；`pca` 另含 `scores`、`loadings`、`eigenvalues`、`explained_variance_ratio`、`sample_ids`、`parameters` 及标记信息；`inbreeding` 另含 `samples`（`sample_id`、`n_markers`、`observed_homozygous`、`expected_homozygous`、`f`）与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import estimate

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = estimate(genotypes, options={"quantity": "grm"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "population",
  "operation": "estimate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "grm",
    "sample_ids": [
      "s1",
      "s2",
      "s3",
      "s4",
      "s5",
      "s6"
    ],
    "method": "vanraden_method_1",
    "ploidy": 2,
    "relationship_matrix": [
      [
        1.2359550561797752,
        0.8988764044943819,
        0.49438202247191015,
        -1.2584269662921348,
        -0.9213483146067415,
        -0.449438202247191
      ],
      [
        0.8988764044943819,
        1.3707865168539324,
        0.15730337078651685,
        -1.1910112359550562,
        -0.8539325842696629,
        -0.38202247191011235
      ],
      "..."
    ],
    "n_samples": 6,
    "denominator": 2.4722222222222223,
    "imputed_value_count": 0,
    "allele_frequencies": [
      0.5,
      0.5833333333333334,
      0.5,
      0.5,
      0.4166666666666667
    ],
    "marker_indices": [
      0,
      1,
      2,
      3,
      4
    ],
    "dropped_all_missing_indices": [],
    "dropped_monomorphic_indices": [],
    "n_markers_used": 5,
    "citations": [
      "https://doi.org/10.3168/jds.2007-0980",
      "https://doi.org/10.1534/g3.112.004259"
    ]
  }
}
```

## 4) export

- **作用：** 把基因型矩阵导出为 TSV/CSV、PLINK 文本 PED+MAP，或等位基因频率表。

- **处理方法：** `tsv`/`csv`：每个样本一行，含 `sample_id` 及每个变异一列。`ped_map`（硬判读）：PLINK 文本文件，基因型单元格为 `ref ref`、`ref alt`、`alt alt`，缺失为 `0 0`，群体标签作为家系 ID；MAP 含染色体、变异 ID、0 和位置；输入中的 `ref`/`alt`/`chrom`/`pos` 列表给出等位基因名（否则用 `A`/`G` 和染色体 0）。`frq`：列为 `SNP`、`A1`、`A2`、`MAF`、`NCHROBS`（2 x 已分型样本数），输入没有等位基因名时用 `ALT`/`REF`。文本总会返回；仅在给出 `path` 时写文件，已有文件不会被覆盖，除非 `overwrite=True`。

- **API：** `dnakit.population.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 `options`：`format` 为 `tsv`（默认）、`csv`、`ped_map` 或 `frq`；`path`（`ped_map` 时为文件名主干，写出 `<path>.ped` 与 `<path>.map`）；`overwrite`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `format`、`n_samples`、`n_variants`，`text`（tsv、csv、frq）或 `ped` 与 `map`（ped_map），以及 `path`（仅在写了文件时）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import export

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = export(genotypes, options={"format": "frq"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "population",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "frq",
    "n_samples": 6,
    "n_variants": 5,
    "text": "SNP\tA1\tA2\tMAF\tNCHROBS\nv1\tALT\tREF\t0.5\t12\nv2\tALT\tREF\t0.41666666666666663\t12\nv3\tALT\tREF\t0.5\t12\nv4\tALT\tREF\t0.5\t12\nv5\tALT\tREF\t0.4166666666666667\t12\n"
  }
}
```

## 5) infer

- **作用：** 推断遗传祖先成分或群体结构。DNAKit 未实现该功能。

- **处理方法：** 不做任何计算。祖先或结构推断需要带标签的参考基因型（如 1000 Genomes）或 ADMIXTURE/STRUCTURE 类的模型拟合程序，本库均未内置。消息中提示可用 `estimate(quantity='pca')` 对自己的样本做 PCA。

- **API：** `dnakit.population.infer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意 `input` 与可选 `options`；它们仅在 `input_summary` 中被概括。

- **输出：** 不计算结果。正常调用返回 `status: "not_implemented"` 的封装（`provenance.implementation` 为 `not_implemented`），并附 `required_dependency` 与 `message` 说明缺少什么：需要祖先参考面板（如 1000 Genomes）或 ADMIXTURE/STRUCTURE 类程序。

- **示例代码：**

```python
import json
from dnakit.population import infer

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = infer(genotypes)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "population",
  "operation": "infer",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "ancestry reference panel (e.g. 1000 Genomes) or an ADMIXTURE/STRUCTURE-type program",
  "message": "Ancestry/structure inference needs labelled reference genotypes or a model-fitting backend, neither of which is bundled; use estimate(quantity='pca') for a PCA of your own samples."
}
```

## 6) plot

- **作用：** 把群体遗传学图（MAF 直方图、PCA 散点图或逐变异 Fst）绘制为 SVG 字符串。

- **处理方法：** `maf_histogram`（默认）：将 MAF 在 [0, 0.5] 上分入 `bins` 个等宽区间。`pca`：`estimate(quantity='pca')` 的前两个主成分，点按 `populations` 着色（需要 numpy，缺失时返回 `not_implemented`）。`fst`：`compare` 给出的逐变异 Fst（`method` 同其含义；需要 `populations`），按变异序号画散点。坐标轴带取整刻度，用纯标准库生成 SVG；未给出 `path` 时不写文件。

- **API：** `dnakit.population.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 `options`：`kind` 为 `maf_histogram`（默认）、`pca` 或 `fst`；`bins`（10）；`method`（用于 `fst`，默认 `weir_cockerham`）；`width`（640）、`height`（420）；`path`（已有文件不会被覆盖，除非 `overwrite=True`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `kind`、`svg`、`path`（未写文件时为 null），以及各类型特有的字段：`bin_counts`、`bin_width`、`n_variants`（直方图）、`n_points`（fst），或 `n_points` 与 `populations`（pca）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import plot

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = plot(genotypes, options={"kind": "maf_histogram", "bins": 5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "population",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "maf_histogram",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"420\" viewBox=\"0 0 640 420\"><title>population maf_histogram p...",
    "path": null,
    "bin_counts": [
      0,
      0,
      0,
      0,
      5
    ],
    "bin_width": 0.1,
    "n_variants": 5
  }
}
```

## 7) simulate

- **作用：** 模拟群体基因型（Balding-Nichols 结构化群体）或中性 Wright-Fisher 等位基因频率漂变。

- **处理方法：** `balding_nichols`（默认）：每个变异的祖先频率 p 在 [maf_min, 1-maf_min] 上均匀抽取；每个群体抽取 p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F)，每个样本的基因型为 Binomial(2, p_k)。`wright_fisher`：每次重复中，起始频率为 `initial_frequency` 的双等位位点漂变 `generations` 代，p' = Binomial(2N, p)/2N（固定或丢失后停止）；同时给出期望方差 p0 q0 (1 - (1 - 1/2N)^t) 与观察方差。给定 `seed`（Python `random`）时结果确定。它不是 msprime 或 SLiM 那样的溯祖或前向模拟器。

- **API：** `dnakit.population.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：None 或选项的映射（优先于 `options`）。`model` 为 `balding_nichols`（默认）：`n_populations`（2）、`n_samples_per_population`（50）、`n_variants`（100）、`fst`（0.1）、`maf_min`（0.1）。`model` 为 `wright_fisher`：`population_size`（100，二倍体）、`generations`（50）、`initial_frequency`（0.5）、`replicates`（100）。两者均有 `seed`（0）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`balding_nichols` 的 `result` 含 `model`、`genotypes`（可直接用于其他函数的基因型映射：`genotypes`、`populations`、`variant_ids`、`sample_ids`）、`ancestral_frequencies`、`population_frequencies`、`fst`、`seed` 与 `definition`。`wright_fisher` 的 `result` 含 `model`、`population_size`、`generations`、`initial_frequency`、`final_frequencies`、`mean_final_frequency`、`variance_final_frequency`、`expected_variance`、`trajectories`（前 20 次重复）、`seed` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import simulate

result = simulate(
    None,
    options={"model": "balding_nichols", "n_populations": 2, "n_samples_per_population": 3,
             "n_variants": 4, "fst": 0.1, "seed": 1},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "population",
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
    "model": "balding_nichols",
    "genotypes": {
      "genotypes": [
        [
          0,
          2,
          2,
          0
        ],
        [
          0,
          1,
          2,
          0
        ],
        "..."
      ],
      "populations": [
        "pop1",
        "pop1",
        "pop1",
        "pop2",
        "pop2",
        "pop2"
      ],
      "variant_ids": [
        "var1",
        "var2",
        "var3",
        "var4"
      ],
      "sample_ids": [
        "s1",
        "s2",
        "s3",
        "s4",
        "s5",
        "s6"
      ]
    },
    "ancestral_frequencies": [
      0.207491395289921,
      0.7779469895497861,
      0.7110196951812913,
      0.30405522059153733
    ],
    "population_frequencies": [
      [
        0.17901345124622997,
        0.8619574841214532,
        0.910923497312446,
        0.20385620033598056
      ],
      [
        0.39824481339193973,
        0.926723126702485,
        0.8023538517295546,
        0.3022080001595772
      ]
    ],
    "fst": 0.1,
    "seed": 1,
    "definition": "Balding-Nichols model: p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F) per population, genotypes ~ Binomial(2, p_k)"
  }
}
```

## 8) test

- **作用：** 逐变异检验 Hardy-Weinberg 平衡（精确检验或卡方检验），或检验群体间的等位基因频率分化。

- **处理方法：** `hwe` 且 `method='exact'`（默认）：Wigginton、Cutler 与 Abecasis (2005) 的精确检验，p 为杂合子数概率不大于观察值的所有构型的概率之和。`method='chisq'`：基因型计数相对于样本频率期望的 Pearson 卡方（1 自由度）（期望计数为 0 时结果为 null）。`differentiation`：对各群体的 2 x k 等位基因计数表做 Pearson 卡方同质性检验（自由度 k - 1；变异单态时为 null）。每个被检验的变异另给出 Bonferroni 校正 p（p x m）和 Benjamini-Hochberg q 值。需要硬判读。

- **API：** `dnakit.population.test(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因型输入为 `{'genotypes': 样本 x 变异的替代等位基因剂量（0/1/2；null 表示缺失）, 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}`，也可只给矩阵（此时样本名为 `sample1`..，变异名为 `variant1`..）。 `differentiation` 必须给出 `populations`。`options`：`test` 为 `hwe`（默认）或 `differentiation`；`method` 为 `exact`（默认）或 `chisq`（用于 `hwe`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `test`、`records`（每个变异：`variant_id`，HWE 的基因型计数 `n_hom_ref`、`n_het`、`n_hom_alt` 或分化检验的 `df`，以及 `chi2`、`p`、`p_bonferroni`、`q_bh`）、`n_tested` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.population import test

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = test(genotypes, options={"test": "hwe", "method": "exact"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "population",
  "operation": "test",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "hwe",
    "records": [
      {
        "variant_id": "v1",
        "n_hom_ref": 2,
        "n_het": 2,
        "n_hom_alt": 2,
        "chi2": null,
        "p": 0.4805194805194805,
        "p_bonferroni": 1.0,
        "q_bh": 0.8008658008658008
      },
      {
        "variant_id": "v2",
        "n_hom_ref": 1,
        "n_het": 3,
        "n_hom_alt": 2,
        "chi2": null,
        "p": 1.0,
        "p_bonferroni": 1.0,
        "q_bh": 1.0
      },
      "..."
    ],
    "n_tested": 5,
    "definition": "exact test of Wigginton et al. (2005): p = sum of probabilities of heterozygote counts no more likely than observed"
  }
}
```
