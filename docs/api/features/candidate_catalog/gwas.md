# GWAS 与多基因风险评分

本类收录 **29** 项有文档的功能。这 29 项都由 DNAKit 原生计算，正常调用均不会返回 `not_implemented`：汇总统计处理（`read`、`standardize`、`validate`、`filter`、`harmonize`、`write`、`export`）、质控与绘图数据（`calculate`、`apply_genomic_control`、`compute_qq_data`、`compute_manhattan_data`、`compute_locus_plot_data`、`plot`）、位点定义与精细定位（`select_lead_snps`、`merge_nearby_significant_hits`、`define_association_loci`、`annotate_lead_snps`、`build_credible_set`、`conditional_gwas`、`joint_gwas_analysis`）、荟萃分析（`fixed_effect_meta_analysis`、`random_effect_meta_analysis`）、遗传力与遗传相关（`estimate_snp_heritability`、`estimate_genetic_correlation`）、多基因评分（`clump_and_threshold_prs`、`evaluate_prs`、`calibrate_prs_risk`），以及原生关联扫描和模拟（`run`、`simulate_gwas_data`）。numpy 是可选的（仅用于加速基因组关系矩阵）。唯一的例外是把 `backend` 设为外部程序（PLINK、REGENIE、GEMMA、GCTA）的 `run`：它返回 `status: "not_implemented"`，因为 DNAKit 不驱动外部程序。这些函数通过 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) annotate_lead_snps

- **作用：** 用用户提供的基因表，为全基因组显著的 lead SNP 注释最近基因和重叠基因。

- **处理方法：** lead SNP 可由 `leads` 给出，或按 `select_lead_snps` 的方式选出（`p_threshold`、`window`）。对每个 lead，与同一染色体上每个基因的距离在基因内为 0，否则为到最近基因边界的距离；报告重叠基因、最近基因（仅当在 `max_distance` 以内）以及 `flank` bp 范围内的基因。不查询任何数据库，也未内置基因模型。

- **API：** `dnakit.gwas.annotate_lead_snps(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`genes`（必需；含 `name`、`chrom`、`start`、`end` 与可选 `strand` 的记录）、`leads`（显式的变异 id；默认以 `p_threshold` = 5e-8 与 `window` = 500000 选出）、`flank`（100000）、`max_distance`（1000000）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `annotations`（每个 lead：标准化后的记录字段、`neglog10p`、`n_nearby_significant`、`overlapping_genes`、`nearest_gene`、`distance_to_nearest_gene`、`relation`、`genes_within_flank`）、`n_leads` 与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import annotate_lead_snps

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
genes = [
    {"name": "GENE1", "chrom": "1", "start": 90000, "end": 120000, "strand": "+"},
    {"name": "GENE2", "chrom": "1", "start": 880000, "end": 940000, "strand": "-"},
    {"name": "GENE3", "chrom": "2", "start": 500000, "end": 520000, "strand": "+"},
]
result = annotate_lead_snps(rows, options={"genes": genes, "flank": 100000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "annotate_lead_snps",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "annotations": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0,
        "overlapping_genes": [
          "GENE1"
        ],
        "nearest_gene": "GENE1",
        "distance_to_nearest_gene": 0,
        "relation": "within_gene",
        "genes_within_flank": [
          "GENE1"
        ]
      },
      {
        "variant_id": "rs3",
        "chrom": "1",
        "pos": 900000,
        "effect_allele": "G",
        "other_allele": "A",
        "beta": -0.15,
        "se": 0.025,
        "p": 2e-09,
        "eaf": 0.4,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0,
        "overlapping_genes": [
          "GENE2"
        ],
        "nearest_gene": "GENE2",
        "distance_to_nearest_gene": 0,
        "relation": "within_gene",
        "genes_within_flank": [
          "GENE2"
        ]
      }
    ],
    "n_leads": 2,
    "method": "nearest gene by interval distance from the supplied gene table"
  }
}
```

## 2) apply_genomic_control

- **作用：** 应用基因组控制：把卡方统计量除以膨胀因子 lambda_GC 并重新计算 p 值。

- **处理方法：** lambda_GC = median(chi2)/0.4549364，其中 chi2 为每个 p 值对应的 1 自由度卡方统计量，0.4549364 是 1 自由度卡方分布的中位数。使用的因子为 lambda_GC，或给定的 `options['lambda_gc']`；`only_inflated=True`（默认）时取 max(1, lambda)，不会向下校正。每条记录增加 `p_gc`（由 chi2/lambda 得到）、`se_gc` = se*sqrt(lambda) 与 `z_gc` = z/sqrt(lambda)。

- **API：** `dnakit.gwas.apply_genomic_control(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `p`（或 `z`，或 `beta` 加 `se`）。`options`：`lambda_gc`（覆盖估计值）、`only_inflated`（True）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `lambda_gc`（估计值）、`lambda_used`、`records`（加入 `p_gc`、`se_gc`、`z_gc` 的标准化记录）与 `n_records`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import apply_genomic_control

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = apply_genomic_control(rows, options={"lambda_gc": 1.2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "apply_genomic_control",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "lambda_gc": 26.65064698571636,
    "lambda_used": 1.2,
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0,
        "p_gc": 4.3695935095779474e-08,
        "se_gc": 0.02190890230020664,
        "z_gc": 5.477225575051661
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5,
        "p_gc": 3.994693946558859e-05,
        "se_gc": 0.02190890230020664,
        "z_gc": 4.107919181288747
      },
      "..."
    ],
    "n_records": 6
  }
}
```

## 3) build_credible_set

- **作用：** 用 Wakefield 近似 Bayes 因子（单一因果变异）由逐变异的 beta 与 se 构建精细定位可信集。

- **处理方法：** log ABF = 0.5*(ln(V/(V+W)) + W/(V+W)*z^2)，V = se^2，z = beta/se，W = prior_sd^2。在均匀先验下后验包含概率 PIP = ABF/sum(ABF)；可信集是 PIP 之和不小于 `coverage` 的、PIP 最高的最小变异集合。方法假定区域内恰有一个因果变异，且不考虑 LD。

- **API：** `dnakit.gwas.build_credible_set(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 表中必须是同一个位点的变异，且含 `beta` 与 `se`。`options`：`prior_sd`（0.15）、`coverage`（0.95）、`region`（`chrom:start-end`，用于截取子集）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `credible_set`（变异 id）、`set_size`、`cumulative_pip`、`coverage_target`、`prior_sd`、`variants`（每个变异的 `log_abf` 与 `pip`，按 PIP 排序）与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import build_credible_set

locus = [
    {"variant_id": "v1", "chrom": "1", "pos": 1000, "beta": 0.30, "se": 0.05},
    {"variant_id": "v2", "chrom": "1", "pos": 1100, "beta": 0.28, "se": 0.05},
    {"variant_id": "v3", "chrom": "1", "pos": 1200, "beta": 0.10, "se": 0.05},
    {"variant_id": "v4", "chrom": "1", "pos": 1300, "beta": 0.02, "se": 0.05},
    {"variant_id": "v5", "chrom": "1", "pos": 1400, "beta": -0.05, "se": 0.05},
]
result = build_credible_set(locus, options={"prior_sd": 0.15, "coverage": 0.95})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "build_credible_set",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 374,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9edc672e99be6afc5d0c531ed01d24ab9a6f80d910fe1821d4838015578a5c21"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "credible_set": [
      "v1",
      "v2"
    ],
    "set_size": 2,
    "cumulative_pip": 0.9999992874121282,
    "coverage_target": 0.95,
    "prior_sd": 0.15,
    "variants": [
      {
        "variant_id": "v1",
        "chrom": "1",
        "pos": 1000,
        "log_abf": 15.048707453502972,
        "pip": 0.8897307255825493
      },
      {
        "variant_id": "v2",
        "chrom": "1",
        "pos": 1100,
        "log_abf": 12.960707453502978,
        "pip": 0.11026856182957885
      },
      "..."
    ],
    "method": "Wakefield approximate Bayes factors, single causal variant per region, uniform prior; log ABF = 0.5*(ln(V/(V+W)) + W/(V+W)*z^2) with V = se^2, W = prior_sd^2"
  }
}
```

## 4) calculate

- **作用：** 计算 GWAS 汇总统计的汇总量；默认给出含基因组膨胀因子 lambda_GC 的表摘要。

- **处理方法：** `summary`：变异数与 p 值个数、最小 p 及其变异、lambda_GC = median(chi2)/0.4549364、chi2（由 p 推得）的均值与中位数、p <= `significance` 的变异数及染色体数。`lambda_gc` 与 `n_significant` 仅给出对应值；`z_scores`、`neglog10p`（-log10 p）和 `maf`（min(eaf, 1-eaf)）给出逐变异的值；`adjusted_p` 给出 Bonferroni 校正 p（p x m，上限为 1）和 Benjamini-Hochberg q 值。

- **API：** `dnakit.gwas.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`statistic` 为 `summary`（默认）、`lambda_gc`、`z_scores`、`neglog10p`、`maf`、`n_significant` 或 `adjusted_p`；`significance`（5e-8）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `statistic`；`summary` 另含 `n_variants`、`n_with_p`、`min_p`、`lead_variant`、`lambda_gc`、`mean_chi2`、`median_chi2`、`significance_threshold`、`n_significant`、`n_chromosomes`；`lambda_gc` 另含 `value` 与 `median_chi2`；`n_significant` 另含 `threshold` 与 `value`；逐变异统计另含 `values`（`adjusted_p` 还有 `n_tests`）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import calculate

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = calculate(rows, options={"statistic": "summary"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "summary",
    "n_variants": 6,
    "n_with_p": 6,
    "min_p": 2e-09,
    "lead_variant": "rs1",
    "lambda_gc": 26.65064698571636,
    "mean_chi2": 16.24053827366549,
    "median_chi2": 12.124350013504214,
    "significance_threshold": 5e-08,
    "n_significant": 2,
    "n_chromosomes": 2
  }
}
```

## 5) calibrate_prs_risk

- **作用：** 用 logistic 重校准模型把多基因评分校准为绝对患病风险。

- **处理方法：** 对 0/1 表型关于标准化评分 z 做 logistic 回归，拟合 risk = expit(a + b*z)（含缺失值的个体被丢弃；至少需要 10 个）。给出 `prevalence` 时，截距平移 logit(prevalence) - logit(样本患病率)，以抵消病例对照抽样的影响。报告 Brier 评分，以及按评分排序分成 `n_bins` 个区间的平均预测风险与观察风险的校准表。

- **API：** `dnakit.gwas.calibrate_prs_risk(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'score': [...], 'phenotype': [0/1, ...]}`。`options`：`prevalence`（总体患病率，取值于 (0,1)）、`n_bins`（10）、`new_scores`（要转换为风险的评分列表）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `intercept`、`slope_per_sd`、`score_mean`、`score_sd`、`sample_prevalence`、`prevalence_adjusted`、`predicted_risk`、`calibration_table`（`bin`、`n`、`mean_predicted`、`observed`）、`brier_score`、`mean_predicted`、`new_score_risk`（给出 `new_scores` 时）与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import calibrate_prs_risk

score = [-1.2, -0.8, -0.5, -0.3, 0.0, 0.1, 0.3, 0.6, 0.9, 1.4, -1.0, -0.2, 0.2, 0.7, 1.1, 1.6]
phenotype = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1]
result = calibrate_prs_risk(
    {"score": score, "phenotype": phenotype},
    options={"n_bins": 4, "new_scores": [-1.0, 0.0, 1.0]},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "calibrate_prs_risk",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "bf1d3cbdd5f0dd2548bfac81709c02e63edd6840e4f64c0e340776a98b2507e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "intercept": -0.46097257162881095,
    "slope_per_sd": 2.360136952448344,
    "score_mean": 0.18125,
    "score_sd": 0.835239486614468,
    "sample_prevalence": 0.4375,
    "prevalence_adjusted": false,
    "predicted_risk": [
      0.012567774392076724,
      0.037917350382972455,
      0.08424795828956246,
      "..."
    ],
    "calibration_table": [
      {
        "bin": 1,
        "n": 4,
        "mean_predicted": 0.039159842928464444,
        "observed": 0.0
      },
      {
        "bin": 2,
        "n": 4,
        "mean_predicted": 0.23107124190436826,
        "observed": 0.5
      },
      "..."
    ],
    "brier_score": 0.13582393605542953,
    "mean_predicted": 0.4375,
    "new_score_risk": [
      0.021906288649246153,
      0.27425689620543464,
      0.8644285709560183
    ],
    "method": "logistic recalibration risk = expit(a + b * standardised score); with `prevalence`, a is shifted by logit(prevalence) - logit(sample prevalence)"
  }
}
```

## 6) clump_and_threshold_prs

- **作用：** 由汇总统计和剂量基因型计算聚类加阈值（C+T）多基因评分。

- **处理方法：** PLINK 风格的贪心聚类：按 p 值升序取变异，每个成为索引变异，并吸收同一染色体上 `clump_window` bp 内、与其 r^2（由 LD 参考基因型计算）不小于 `clump_r2` 且尚未归入聚类的变异。对每个 p 值阈值，样本评分为 p <= 阈值的索引变异上 sum(beta*剂量)，缺失剂量以该变异均值代替。LD 参考默认就是目标基因型本身。

- **API：** `dnakit.gwas.clump_and_threshold_prs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'summary_statistics': 含 variant_id、chrom、pos、beta、p 的表, 'genotypes': 样本 x 变异的效应等位基因剂量 [0, 2], 'sample_ids': [...], 'variant_ids': [...]}`。`options`：`p_thresholds`（`[5e-8, 1e-5, 1e-3, 0.05, 0.5, 1]`）、`clump_r2`（0.1）、`clump_window`（250000）、`ld_reference`（基因型映射；默认即目标基因型）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `clumps`（`index`、`p`、`beta`、`members`）、`thresholds`（每个阈值：`p_threshold`、`n_variants`、`variants` 与 `scores` {sample_id: 评分}）、`clump_r2`、`clump_window` 与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import clump_and_threshold_prs

summary = [
    {"variant_id": "v1", "chrom": "1", "pos": 1000, "beta": 0.30, "p": 1e-9},
    {"variant_id": "v2", "chrom": "1", "pos": 2000, "beta": 0.25, "p": 1e-7},
    {"variant_id": "v3", "chrom": "1", "pos": 3000, "beta": -0.20, "p": 1e-4},
]
genotypes = [[0, 0, 2], [1, 1, 1], [2, 2, 0], [1, 1, 0], [0, 0, 1], [2, 2, 2]]
result = clump_and_threshold_prs(
    {"summary_statistics": summary, "genotypes": genotypes,
     "sample_ids": ["a", "b", "c", "d", "e", "f"], "variant_ids": ["v1", "v2", "v3"]},
    options={"p_thresholds": [1e-8, 1e-3], "clump_r2": 0.5},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "clump_and_threshold_prs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 411,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3c2cc5467f66341cafe87e4b0959583a1cf91dd85b0704c13b015e61777908d7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "clumps": [
      {
        "index": "v1",
        "p": 1e-09,
        "beta": 0.3,
        "members": [
          "v2"
        ]
      },
      {
        "index": "v3",
        "p": 0.0001,
        "beta": -0.2,
        "members": []
      }
    ],
    "thresholds": [
      {
        "p_threshold": 1e-08,
        "n_variants": 1,
        "variants": [
          "v1"
        ],
        "scores": {
          "a": 0.0,
          "b": 0.3,
          "c": 0.6,
          "d": 0.3,
          "e": 0.0,
          "f": 0.6
        }
      },
      {
        "p_threshold": 0.001,
        "n_variants": 2,
        "variants": [
          "v1",
          "v3"
        ],
        "scores": {
          "a": -0.4,
          "b": 0.09999999999999998,
          "c": 0.6,
          "d": 0.3,
          "e": -0.2,
          "f": 0.19999999999999996
        }
      }
    ],
    "clump_r2": 0.5,
    "clump_window": 250000,
    "method": "PLINK-style greedy clumping by p-value (r^2 from the LD reference, same chromosome, within window bp), then score = sum(..."
  }
}
```

## 7) compute_locus_plot_data

- **作用：** 计算以 lead 变异为中心的区域（LocusZoom 风格）关联图的点。

- **处理方法：** 区域为 lead 所在染色体上 lead 位置 +/- `flank` bp（起点不小于 1）；区域内每个变异成为一个点，纵坐标为 -log10 p。若提供 LD 矩阵或参考基因型，则计算与 lead 的 r^2（剂量 Pearson 相关系数的平方，或给定 r 的平方），并把每个点归入区间（0-0.2、0.2-0.4、0.4-0.6、0.6-0.8、0.8-1、`lead`，无 LD 信息时为 `unknown`）。不会虚构 LD。

- **API：** `dnakit.gwas.compute_locus_plot_data(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`lead`（变异 id；默认 p 最小者）、`flank`（250000）、`ld`（`{'variant_ids': [...], 'matrix': [[...]]}`，相关系数 r）或 `genotypes`（参考剂量），用于按 r^2 着色。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `lead`（`variant_id`、`chrom`、`pos`、`p`）、`region`（`chrom`、`start`、`end`）与 `points`（`variant_id`、`pos`、`neglog10p`、`r2`、`ld_band`）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import compute_locus_plot_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_locus_plot_data(rows, options={"lead": "rs1", "flank": 100000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "compute_locus_plot_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "lead": {
      "variant_id": "rs1",
      "chrom": "1",
      "pos": 100000,
      "p": 2e-09
    },
    "region": {
      "chrom": "1",
      "start": 1,
      "end": 200000
    },
    "points": [
      {
        "variant_id": "rs1",
        "pos": 100000,
        "neglog10p": 8.698970004336019,
        "r2": 1.0,
        "ld_band": "lead"
      },
      {
        "variant_id": "rs2",
        "pos": 150000,
        "neglog10p": 5.167491087293763,
        "r2": null,
        "ld_band": "unknown"
      }
    ]
  }
}
```

## 8) compute_manhattan_data

- **作用：** 由汇总统计计算曼哈顿图坐标（-log10 p 与累计基因组 x 位置）。

- **处理方法：** x = 染色体偏移 + pos，染色体偏移为前面各染色体观测到的最大位置之和（染色体顺序为 1..22、X、Y、MT）。y = -log10 p。弱于 `min_neglog10p` 的点被丢弃；超过 `max_points` 时，最显著的一半点原样保留，其余抽稀。同时报告全基因组（5e-8）与提示性（1e-5）阈值和各染色体标签中心。

- **API：** `dnakit.gwas.compute_manhattan_data(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `chrom`、`pos` 与 `p`。`options`：`min_neglog10p`（0）、`max_points`（200000）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `points`（`variant_id`、`chrom`、`pos`、`x`、`neglog10p`）、`chromosomes`（`chrom`、`offset`、`length`、`center`）、`total_length`、`truncated`、`thresholds`（`genome_wide`、`suggestive`）与 `x_definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import compute_manhattan_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_manhattan_data(rows)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "compute_manhattan_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "points": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "x": 100000,
        "neglog10p": 8.698970004336019
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "x": 150000,
        "neglog10p": 5.167491087293763
      },
      "..."
    ],
    "chromosomes": [
      {
        "chrom": "1",
        "offset": 0,
        "length": 900000,
        "center": 450000.0
      },
      {
        "chrom": "2",
        "offset": 900000,
        "length": 90000,
        "center": 945000.0
      }
    ],
    "total_length": 990000,
    "truncated": false,
    "thresholds": {
      "genome_wide": {
        "p": 5e-08,
        "neglog10p": 7.301029995663981
      },
      "suggestive": {
        "p": 1e-05,
        "neglog10p": 5.0
      }
    },
    "x_definition": "cumulative base-pair position: chromosome offset (sum of preceding chromosomes' maximum observed positions) + pos"
  }
}
```

## 9) compute_qq_data

- **作用：** 计算 QQ 图数据（观察的与期望的 -log10 p）和 lambda_GC。

- **处理方法：** p 值升序排列、最小者秩为 1 时，期望值 = -log10((rank - 0.5)/n)，观察值 = -log10 p。超过 `max_points` 时，最显著的一半点原样保留，其余抽稀。lambda_GC = median(chi2)/0.4549364。

- **API：** `dnakit.gwas.compute_qq_data(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `p`。`options`：`max_points`（50000）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `n`、`points`（`rank`、`expected`、`observed`）、`lambda_gc`、`downsampled` 与 `expected_definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import compute_qq_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_qq_data(rows)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "compute_qq_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 6,
    "points": [
      {
        "rank": 6,
        "expected": 0.0377885608893998,
        "observed": 0.2076083105017461
      },
      {
        "rank": 5,
        "expected": 0.12493873660829993,
        "observed": 0.49894073778224846
      },
      "..."
    ],
    "lambda_gc": 26.65064698571636,
    "downsampled": false,
    "expected_definition": "-log10((rank - 0.5) / n) with rank 1 = smallest p"
  }
}
```

## 10) conditional_gwas

- **作用：** 由汇总统计和 LD 做近似条件关联分析（COJO 风格，精确的 OLS 代数）。

- **处理方法：** 每个变异的边际相关为 r = z/sqrt(n - 2 + z^2)。对条件集 C 之外的每个变异，其标准化效应、标准误、t 统计量与 p 值来自联合 OLS 模型 {C + 该变异}，由边际相关和 LD 相关矩阵精确算出，t 检验自由度为 n - |C| - 2。同时报告 `r2_with_conditioning_set`。结果的质量取决于所给的 LD 矩阵（假定为同样本 LD 参考）。

- **API：** `dnakit.gwas.conditional_gwas(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `z`（或 `beta` 与 `se`），最好有 `n`。`options`：`condition_on`（必需；1 至 20 个变异 id）、`ld`（`{'variant_ids': [...], 'matrix': [[...]]}`，相关系数矩阵 r）或 `genotypes`（参考剂量）、`n`（样本量；默认取 `n` 列的中位数）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `conditioned_on`、`n`、`records`（每个其他变异：`variant_id`、`marginal_z`、`marginal_p`、`conditional_beta_std`、`conditional_se_std`、`conditional_t`、`conditional_p`、`r2_with_conditioning_set`）与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import conditional_gwas

summary = [
    {"variant_id": "v1", "beta": 0.30, "se": 0.05, "n": 1000},
    {"variant_id": "v2", "beta": 0.27, "se": 0.05, "n": 1000},
    {"variant_id": "v3", "beta": 0.05, "se": 0.05, "n": 1000},
]
ld = {"variant_ids": ["v1", "v2", "v3"],
      "matrix": [[1.0, 0.9, 0.1], [0.9, 1.0, 0.1], [0.1, 0.1, 1.0]]}
result = conditional_gwas(summary, options={"condition_on": ["v1"], "ld": ld})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "conditional_gwas",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 176,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "332f55a26240696b32413b21615aea7711c83f5056381b5d1ae715c26ff00db2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "conditioned_on": [
      "v1"
    ],
    "n": 1000.0,
    "records": [
      {
        "variant_id": "v2",
        "marginal_z": 5.4,
        "marginal_p": 6.664089697085715e-08,
        "conditional_beta_std": 0.0029379654824893686,
        "conditional_se_std": 0.07138060961532508,
        "conditional_t": 0.04115915370185633,
        "conditional_p": 0.9671772573822426,
        "r2_with_conditioning_set": 0.81
      },
      {
        "variant_id": "v3",
        "marginal_z": 1.0,
        "marginal_p": 0.31731050786291415,
        "conditional_beta_std": 0.013110593385134293,
        "conditional_se_std": 0.031268103211206585,
        "conditional_t": 0.41929608894329784,
        "conditional_p": 0.6750900785228465,
        "r2_with_conditioning_set": 0.010000000000000009
      }
    ],
    "method": "exact joint OLS from marginal correlations and LD (COJO-style): standardised phenotype and genotypes, marginal r from z via r = z/sqrt(n-2+z^2), t-test with df = n - |C| - 2"
  }
}
```

## 11) define_association_loci

- **作用：** 把独立关联位点定义为 lead SNP +/- flank，并合并重叠窗口。

- **处理方法：** 按 `select_lead_snps` 的方式选出 lead SNP（`p_threshold`、`window`）；每个 lead 定义区间 [pos - flank, pos + flank]（起点不小于 1），同一染色体上重叠的区间合并为一个位点。每个位点报告其 lead、最强变异、所含变异数和显著变异数。

- **API：** `dnakit.gwas.define_association_loci(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `chrom`、`pos` 与 `p`。`options`：`p_threshold`（5e-8）、`window`（500000，用于选 lead）、`flank`（500000）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `loci`（`chrom`、`start`、`end`、`leads`、`lead_variant`、`lead_p`、`n_variants`、`n_significant`）、`n_loci` 与 `locus_definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import define_association_loci

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = define_association_loci(rows, options={"p_threshold": 5e-8, "flank": 100000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "define_association_loci",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "loci": [
      {
        "chrom": "1",
        "start": 1,
        "end": 200000,
        "leads": [
          "rs1"
        ],
        "lead_variant": "rs1",
        "lead_p": 2e-09,
        "n_variants": 2,
        "n_significant": 1
      },
      {
        "chrom": "1",
        "start": 800000,
        "end": 1000000,
        "leads": [
          "rs3"
        ],
        "lead_variant": "rs3",
        "lead_p": 2e-09,
        "n_variants": 1,
        "n_significant": 1
      }
    ],
    "n_loci": 2,
    "locus_definition": "lead variant +/- flank bp (clamped at 1), overlapping loci on one chromosome merged"
  }
}
```

## 12) estimate_genetic_correlation

- **作用：** 估计两个性状之间基于 SNP 的遗传相关。

- **处理方法：** 使用个体数据时，在 VanRaden GRM 上做跨性状 Haseman-Elston 回归：表型缩放为单位方差后，h2_1 和 h2_2 为 y_i*y_j 对 GRM_ij 在所有 i<j 对上的斜率（含截距），协方差为 0.5*(y1_i*y2_j + y1_j*y2_i) 的斜率，rg = 协方差/sqrt(h2_1*h2_2)。使用两张汇总统计表时，采用未加权的跨性状 LD 得分回归：各性状的 h2 来自 chi2 对 LD 得分的 OLS，rho_g 为 z1*z2 对 LD 得分的斜率乘以 M/sqrt(N1*N2)（性状 2 的等位基因与性状 1 对齐），rg = rho_g/sqrt(h1*h2)。任一遗传力估计不为正时 rg 为 `null`（并给出警告），估计值不受限于 [-1, 1]（它是两个带噪声估计量的比值，小样本下可能得到 1.008 这样的值）：数值原样保留，不做截断，且只要 |rg| > 1 就会添加一条警告。单个性状的 h2 估计小于 0 或大于 1 时也会添加警告，汇总统计路径还总会添加一条“未计算标准误”的警告。不给出标准误。GRM 在安装了 numpy 时使用 numpy，否则走纯 Python 路径（乘加次数超过约 3e7 时拒绝）。

- **API：** `dnakit.gwas.estimate_genetic_correlation(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'genotypes': 样本 x 变异的剂量, 'phenotype_1': [...], 'phenotype_2': [...]}`，或带 `z` 与 LD 得分的 `{'trait1': 表, 'trait2': 表}`（`ldscore` 列或 `options['ld_scores']`；至少 10 个变异）。`options`：`ld_scores`（{variant_id: 得分}）、`n1`、`n2`、`m`、`fix_intercept`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。个体数据的 `result` 含 `method`（`haseman_elston_cross_product`）、`rg`、`genetic_covariance_std`、`h2_trait1`、`h2_trait2`、`n_samples` 与 `definition`；汇总统计的 `result` 含 `method`（`ldsc_unweighted`）、`rg`、`genetic_covariance`、`h2_trait1`、`h2_trait2`、`cross_intercept`、`n_variants` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import estimate_genetic_correlation

genotypes = [
    [0, 1, 2, 0, 1, 2], [1, 0, 1, 2, 2, 0], [2, 2, 0, 1, 0, 1], [0, 1, 1, 0, 1, 2],
    [1, 2, 0, 2, 0, 1], [2, 0, 2, 1, 1, 0], [0, 2, 1, 1, 2, 2], [1, 1, 2, 0, 0, 1],
    [2, 1, 0, 2, 1, 0], [0, 0, 2, 2, 1, 1], [1, 2, 1, 0, 2, 0], [2, 0, 0, 1, 2, 2],
]
phenotype_1 = [1.2, 0.3, -0.8, 0.9, -0.1, 0.4, 1.0, 0.2, -1.1, 0.6, -0.4, -0.2]
phenotype_2 = [0.9, 0.5, -0.6, 0.7, -0.3, 0.1, 1.1, 0.0, -0.9, 0.8, -0.2, -0.5]
result = estimate_genetic_correlation(
    {"genotypes": genotypes, "phenotype_1": phenotype_1, "phenotype_2": phenotype_2}
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "estimate_genetic_correlation",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 419,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d01b512cd38fea620a4c1050f752fb747559dbe156281ace6fa12cd20548bcd"
  },
  "warnings": ["rg = 1.008 is outside [-1, 1]: the ratio estimator is not bounded and small samples or noisy covariance can push it past 1; the value is reported unclamped and carries large uncertainty"],
  "status": "ok",
  "result": {
    "method": "haseman_elston_cross_product",
    "rg": 1.0079706026145199,
    "genetic_covariance_std": 0.8772411283478185,
    "h2_trait1": 0.9015109250540984,
    "h2_trait2": 0.8401778798532548,
    "n_samples": 12,
    "definition": "rg = slope(0.5*(y1_i*y2_j+y1_j*y2_i) ~ GRM_ij) / sqrt(h2_1*h2_2), phenotypes scaled to unit variance"
  }
}
```

## 13) estimate_snp_heritability

- **作用：** 用 Haseman-Elston 回归（个体数据）或未加权 LD 得分回归（汇总数据）估计 SNP 遗传力。

- **处理方法：** `haseman_elston`：表型缩放为单位方差后，h2 为 y_i*y_j 对 VanRaden GRM_ij 在所有 i<j 对上的斜率（含截距）；单态标记被剔除。`ldsc`：chi2 = z^2 对 LD 得分 l_j 做 OLS，E[chi2] = N*h2/M*l_j + 截距，故 h2 = 斜率*M/N（截距自由，`fix_intercept` 时固定为 1）；需要至少 10 个同时有 z 与 LD 得分的变异。不给出标准误或 jackknife，小样本下估计可能为负或大于 1。数值原样保留，不做截断；只要 h2 < 0 或 h2 > 1 就会添加一条警告，汇总统计（`ldsc`）路径还总会添加一条“未计算标准误”的警告。GRM 在安装了 numpy 时使用 numpy，否则走纯 Python 路径（乘加次数超过约 3e7 时拒绝）。

- **API：** `dnakit.gwas.estimate_snp_heritability(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'genotypes': 样本 x 变异的剂量, 'phenotype': [...]}`，或带 `z` 与 LD 得分的汇总统计表（`ld_score`/`ldscore` 列或 `options['ld_scores']`）。`options`：`method` 为 `auto`（默认；有 `genotypes` 时用 Haseman-Elston，否则用 `ldsc`）、`haseman_elston` 或 `ldsc`；`ld_scores`（{variant_id: 得分}）、`n`（默认取 `n` 列中位数）、`m`（默认取所用变异数）、`fix_intercept`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。Haseman-Elston：`result` 含 `method`（`haseman_elston_cross_product`）、`h2`、`n_samples`、`n_markers` 与 `definition`。LD 得分回归：`method`（`ldsc_unweighted`）、`h2`、`intercept`、`slope`、`n`、`m`、`n_variants` 与 `definition`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import estimate_snp_heritability

genotypes = [
    [0, 1, 2, 0, 1, 2], [1, 0, 1, 2, 2, 0], [2, 2, 0, 1, 0, 1], [0, 1, 1, 0, 1, 2],
    [1, 2, 0, 2, 0, 1], [2, 0, 2, 1, 1, 0], [0, 2, 1, 1, 2, 2], [1, 1, 2, 0, 0, 1],
    [2, 1, 0, 2, 1, 0], [0, 0, 2, 2, 1, 1], [1, 2, 1, 0, 2, 0], [2, 0, 0, 1, 2, 2],
]
phenotype = [1.2, 0.3, -0.8, 0.9, -0.1, 0.4, 1.0, 0.2, -1.1, 0.6, -0.4, -0.2]
result = estimate_snp_heritability({"genotypes": genotypes, "phenotype": phenotype})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "estimate_snp_heritability",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 335,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f045666eb6204e967b3cc93c9d1e87ae7425f0fbf3754fb48f58fbbff9b878bf"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "haseman_elston_cross_product",
    "h2": 0.9015109250540984,
    "n_samples": 12,
    "n_markers": 6,
    "definition": "slope of y_i*y_j on VanRaden GRM_ij over all pairs i<j (with intercept), y scaled to unit variance"
  }
}
```

## 14) evaluate_prs

- **作用：** 对照表型评估多基因评分：定量性状给出 R^2 与斜率，二分类性状给出 AUC 与比值比。

- **处理方法：** 定量性状：Pearson r、r^2，以及表型对标准化评分（每 SD）的 OLS 斜率及 t 检验 p 值。二分类性状：由 Mann-Whitney 秩统计量得到 AUC（并列取平均秩），并由对标准化评分的 logistic 回归得到每 SD 的比值比、Wald 95% 置信区间与 p 值，以及 Nagelkerke R^2；logistic 拟合不收敛（完全分离）时只返回 AUC 和说明。含缺失值的个体被丢弃（至少需要 10 个）。

- **API：** `dnakit.gwas.evaluate_prs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'score': [...], 'phenotype': [...]}`。`options`：`trait` 为 `auto`（默认；表型为含两类的 0/1 时按二分类）、`quantitative` 或 `binary`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。定量性状：`result` 含 `trait`、`n`、`r`、`r2`、`slope_per_sd`、`slope_se`、`p`。二分类：`trait`、`n`、`n_cases`、`n_controls`、`auc`、`odds_ratio_per_sd`、`or_ci95`、`log_or`、`se`、`p`、`nagelkerke_r2`（logistic 拟合失败时为 `note`）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import evaluate_prs

score = [-1.2, -0.8, -0.5, -0.3, 0.0, 0.1, 0.3, 0.6, 0.9, 1.4, -1.0, -0.2, 0.2, 0.7, 1.1, 1.6]
phenotype = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1]
result = evaluate_prs({"score": score, "phenotype": phenotype}, options={"trait": "binary"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "evaluate_prs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "bf1d3cbdd5f0dd2548bfac81709c02e63edd6840e4f64c0e340776a98b2507e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "trait": "binary",
    "n": 16,
    "n_cases": 7,
    "n_controls": 9,
    "auc": 0.8888888888888888,
    "odds_ratio_per_sd": 10.592402008491819,
    "or_ci95": [
      1.1190821428442463,
      100.25982545332894
    ],
    "log_or": 2.360136952448344,
    "se": 1.1467701121929976,
    "p": 0.03958309732067795,
    "nagelkerke_r2": 0.5816456123759354
  }
}
```

## 15) export

- **作用：** 把 GWAS 汇总统计序列化为另一种文件格式（tsv、csv、json、jsonl、bed、ldsc、metal）。

- **处理方法：** 先标准化表格，再按所选布局写出：`tsv`/`csv`（全部非空的规范列，或 `columns`）、`json`、`jsonl`、`bed`（带 chr 前缀的染色体、0 起点的 start、end、变异 id，以及作为得分的 -log10 p）、`ldsc`（munge 后的 `SNP A1 A2 N Z` 布局）或 `metal`（`MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value`）。文本总会返回；仅在给出 `path` 时写文件（`.gz` 会被 gzip 压缩），已有文件不会被覆盖，除非 `overwrite=True`。

- **API：** `dnakit.gwas.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`format` 为 `tsv`（默认）、`csv`、`json`、`jsonl`、`bed`、`ldsc` 或 `metal`；`columns`（tsv/csv 的列子集）；`path`（或 `output`）；`overwrite`（False）。`bed` 需要 `chrom` 与 `pos`；`ldsc` 需要 `variant_id`、两个等位基因和 `z`；`metal` 需要 `variant_id`、两个等位基因、`beta` 与 `se`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `format`、`n_records`、`text`（序列化后的表格）与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import export

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = export(rows[:3], options={"format": "metal"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 374,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "879d97b1e5ae6c404c3f251cf85e86d9b0381c2aa31d5343b9285bc0814d8299"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "metal",
    "n_records": 3,
    "text": "MarkerName\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP-value\nrs1\tA\tG\t0.3\t0.12\t0.02\t2e-09\nrs2\tC\tT\t0.25\t0.09\t0.02\t6.8e-06\nrs3\tG\tA\t0.4\t-0.15\t0.025\t2e-09\n",
    "path": null
  }
}
```

## 16) filter

- **作用：** 按 p 值、MAF、样本量、INFO、染色体、区域或变异类型过滤 GWAS 汇总统计的行。

- **处理方法：** 每一行按固定顺序（`require_p`、`max_p`、`min_maf`（MAF = min(eaf, 1-eaf)）、`min_n`、`min_info`、`chromosomes`、`region`、`snp_only`、`remove_palindromic`）检验各个已启用的过滤条件；在已启用条件上有缺失值的行被移除，被移除的行只记在第一个适用的原因之下。

- **API：** `dnakit.gwas.filter(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`max_p`、`min_maf`、`min_n`、`min_info`、`chromosomes`（列表）、`region`（`chrom:start-end`）、`snp_only`（False）、`remove_palindromic`（False）、`require_p`（False）。未给出的过滤条件不起作用。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `records`（保留的标准化行）、`n_input`、`n_kept` 与 `removed`（{原因: 计数}，原因如 `p_above_threshold`、`maf_below_threshold`、`n_below_threshold`、`info_below_threshold`、`chromosome_excluded`、`outside_region`、`not_a_snp`、`palindromic`、`missing_p`）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import filter

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = filter(rows, options={"max_p": 0.05, "min_maf": 0.2, "chromosomes": ["1", "2"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "filter",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5
      },
      "..."
    ],
    "n_input": 6,
    "n_kept": 4,
    "removed": {
      "p_above_threshold": 2
    }
  }
}
```

## 17) fixed_effect_meta_analysis

- **作用：** 用逆方差加权的固定效应荟萃分析合并各研究中逐变异的效应。

- **处理方法：** 对每个变异（按 `variant_id` 或 `chrom:pos` 匹配），权重 w = 1/se^2，合并 beta = sum(w*beta)/sum(w)，se = 1/sqrt(sum(w))，z = beta/se 并给出双侧正态 p 值，Cochran Q = sum(w*(beta - 合并值)^2) 及其 k-1 自由度卡方 p 值，以及 I^2 = max(0, (Q - df)/Q)。各研究的等位基因与第一个报告该变异的研究对齐（等位基因互换时取负效应；等位基因不匹配的研究被跳过，不处理链翻转）。出现在少于 `min_studies` 个研究中的变异被略去。

- **API：** `dnakit.gwas.fixed_effect_meta_analysis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：字典 {研究名: 表} 或表的列表（每张表含 `variant_id` 或 `chrom` + `pos`、`beta`、`se` 和等位基因），或带 `study` 列的单张表。`options`：`min_studies`（2）、`study_column`（`study`）。至少需要两个研究。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `records`（每个变异：`variant_id`、`chrom`、`pos`、`effect_allele`、`other_allele`、`n_studies`、`beta`、`se`、`z`、`p`、`q`、`q_p`、`i2`）、`n_variants` 与 `model`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import fixed_effect_meta_analysis

study_a = [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.10, "se": 0.03},
           {"variant_id": "rs2", "effect_allele": "C", "other_allele": "T", "beta": 0.02, "se": 0.03}]
study_b = [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.14, "se": 0.04},
           {"variant_id": "rs2", "effect_allele": "T", "other_allele": "C", "beta": -0.05, "se": 0.04}]
result = fixed_effect_meta_analysis({"A": study_a, "B": study_b})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "fixed_effect_meta_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 382,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3113603d7296e5791a4642b96f6fd8c8a8dd6407e46ab7d703c2c4f3b959d454"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": null,
        "pos": null,
        "effect_allele": "A",
        "other_allele": "G",
        "n_studies": 2,
        "beta": 0.11440000000000002,
        "se": 0.024,
        "z": 4.7666666666666675,
        "p": 1.872988003352597e-06,
        "q": 0.6400000000000003,
        "q_p": 0.4237107971667936,
        "i2": 0.0
      },
      {
        "variant_id": "rs2",
        "chrom": null,
        "pos": null,
        "effect_allele": "C",
        "other_allele": "T",
        "n_studies": 2,
        "beta": 0.0308,
        "se": 0.024,
        "z": 1.2833333333333334,
        "p": 0.19937531925904434,
        "q": 0.36,
        "q_p": 0.5485062355001472,
        "i2": 0.0
      }
    ],
    "n_variants": 2,
    "model": "fixed effect, inverse-variance weights w=1/se^2; Cochran's Q and I^2 reported"
  }
}
```

## 18) harmonize

- **作用：** 把汇总统计表的效应等位基因与参考表对齐（等位基因互换、链翻转与回文位点处理）。

- **处理方法：** 变异按 `chrom:pos`（或 id；`match_on`='auto' 在每行都有位置时使用位置）与参考表匹配。若效应等位基因是参考的另一个等位基因，则互换等位基因，beta 与 z 取负，eaf 变为 1-eaf；位于反链的等位基因取互补（`strand_flip`、`strand_flip_swap`）。非回文的不匹配位点被丢弃。回文（A/T、C/G）SNP 有歧义：`palindromes='drop_ambiguous'` 时，若目标或参考的频率与 0.5 相差在 `palindrome_eaf_tolerance` 以内则丢弃，否则通过比较频率来判定；`assume_forward_strand` 时直接保留，不推断链。每个保留的行带有 `harmonize_action`。

- **API：** `dnakit.gwas.harmonize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'target': 表, 'reference': 表}`，或一张表加 `options['reference']`；两者都需要效应等位基因和另一等位基因。`options`：`match_on` 为 `auto`（默认）、`id` 或 `position`；`palindromes` 为 `drop_ambiguous`（默认）或 `assume_forward_strand`；`palindrome_eaf_tolerance`（0.08）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `records`（带 `harmonize_action` 的已对齐行：`none`、`swap`、`strand_flip`、`strand_flip_swap`）、`n_target`、`n_harmonized`、`actions`（计数，包含丢弃原因 `not_in_reference`、`missing_alleles`、`allele_mismatch`、`ambiguous_palindrome`）与 `matched_on`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import harmonize

target = [
    {"variant_id": "rs1", "chrom": "1", "pos": 100, "effect_allele": "A", "other_allele": "G", "beta": 0.10, "eaf": 0.30},
    {"variant_id": "rs2", "chrom": "1", "pos": 200, "effect_allele": "T", "other_allele": "C", "beta": 0.20, "eaf": 0.20},
    {"variant_id": "rs3", "chrom": "1", "pos": 300, "effect_allele": "A", "other_allele": "T", "beta": 0.05, "eaf": 0.49},
]
reference = [
    {"variant_id": "rs1", "chrom": "1", "pos": 100, "effect_allele": "A", "other_allele": "G"},
    {"variant_id": "rs2", "chrom": "1", "pos": 200, "effect_allele": "C", "other_allele": "T"},
    {"variant_id": "rs3", "chrom": "1", "pos": 300, "effect_allele": "A", "other_allele": "T"},
]
result = harmonize({"target": target, "reference": reference})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "harmonize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 656,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "cf206a3f2cac9b758aa78197bc7656674394c05b3e6cc4872f65ad7220010816"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.1,
        "se": null,
        "p": null,
        "eaf": 0.3,
        "n": null,
        "z": null,
        "harmonize_action": "none"
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 200,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": -0.2,
        "se": null,
        "p": null,
        "eaf": 0.8,
        "n": null,
        "z": null,
        "harmonize_action": "swap"
      }
    ],
    "n_target": 3,
    "n_harmonized": 2,
    "actions": {
      "none": 1,
      "swap": 1,
      "ambiguous_palindrome": 1
    },
    "matched_on": "position"
  }
}
```

## 19) joint_gwas_analysis

- **作用：** 由汇总统计和 LD 运行多 SNP 联合关联模型（COJO 联合风格）。

- **处理方法：** 由边际相关 r = z/sqrt(n - 2 + z^2) 和 LD 相关矩阵 R 精确计算标准化表型对所选标准化基因型的多元回归：beta = R^-1 r，模型 R^2 = r'beta，se^2 = (1 - R^2)/(n - k - 1) * diag(R^-1)，t 检验自由度为 n - k - 1。LD 矩阵奇异（完全 LD）时报错。结果的质量取决于所给的 LD 矩阵（假定为同样本 LD 参考）。

- **API：** `dnakit.gwas.joint_gwas_analysis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `z`（或 `beta` 与 `se`），最好有 `n`。`options`：`ld`（`{'variant_ids': [...], 'matrix': [[...]]}`，相关系数矩阵 r）或 `genotypes`（参考剂量）、`snps`（联合建模的变异 id，至多 200 个；默认取 LD 矩阵中且在表内的全部变异）、`n`（默认取 `n` 列中位数）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `n`、`records`（每个变异：`variant_id`、`marginal_z`、`joint_beta_std`、`joint_se_std`、`joint_t`、`joint_p`）、`model_r2`、`df` 与 `method`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import joint_gwas_analysis

summary = [
    {"variant_id": "v1", "beta": 0.30, "se": 0.05, "n": 1000},
    {"variant_id": "v2", "beta": 0.27, "se": 0.05, "n": 1000},
    {"variant_id": "v3", "beta": 0.05, "se": 0.05, "n": 1000},
]
ld = {"variant_ids": ["v1", "v2", "v3"],
      "matrix": [[1.0, 0.9, 0.1], [0.9, 1.0, 0.1], [0.1, 0.1, 1.0]]}
result = joint_gwas_analysis(summary, options={"ld": ld})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "joint_gwas_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 176,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "332f55a26240696b32413b21615aea7711c83f5056381b5d1ae715c26ff00db2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 1000.0,
    "records": [
      {
        "variant_id": "v1",
        "marginal_z": 5.999999999999999,
        "joint_beta_std": 0.1832581196160083,
        "joint_se_std": 0.07142915304024819,
        "joint_t": 2.565592784122021,
        "joint_p": 0.01044535198931147
      },
      {
        "variant_id": "v2",
        "marginal_z": 5.4,
        "joint_beta_std": 0.002249129961401298,
        "joint_se_std": 0.07142915304024819,
        "joint_t": 0.03148756307013721,
        "joint_p": 0.9748870188760195
      },
      {
        "variant_id": "v3",
        "marginal_z": 1.0,
        "joint_beta_std": 0.013087874900675695,
        "joint_se_std": 0.03129209958977394,
        "joint_t": 0.4182485378818342,
        "joint_p": 0.6758555722527371
      }
    ],
    "model_r2": 0.03498737698481886,
    "df": 996.0,
    "method": "exact multiple regression of the standardised phenotype on the selected standardised genotypes, computed from marginal correlations and LD (beta = R^-1 r, se^2 = (1-R2)/(n-k-1) * diag(R^-1))"
  }
}
```

## 20) merge_nearby_significant_hits

- **作用：** 按碱基对距离的单连锁规则把邻近的全基因组显著变异合并为命中区域。

- **处理方法：** p <= `p_threshold` 的变异按染色体和位置排序；同一染色体上相距不超过 `distance` bp 的连续显著变异归入同一区域。每个区域报告起止位置、命中数和 lead（p 最小）变异。

- **API：** `dnakit.gwas.merge_nearby_significant_hits(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `chrom`、`pos` 与 `p`。`options`：`p_threshold`（5e-8）、`distance`（500000）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `groups`（`chrom`、`start`、`end`、`n_hits`、`lead_variant`、`lead_p`）、`n_groups`、`n_hits`、`distance` 与 `rule`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import merge_nearby_significant_hits

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = merge_nearby_significant_hits(rows, options={"p_threshold": 1e-5, "distance": 100000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "merge_nearby_significant_hits",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "groups": [
      {
        "chrom": "1",
        "start": 100000,
        "end": 150000,
        "n_hits": 2,
        "lead_variant": "rs1",
        "lead_p": 2e-09
      },
      {
        "chrom": "1",
        "start": 900000,
        "end": 900000,
        "n_hits": 1,
        "lead_variant": "rs3",
        "lead_p": 2e-09
      }
    ],
    "n_groups": 2,
    "n_hits": 3,
    "distance": 100000,
    "rule": "single-linkage: consecutive significant variants closer than or equal to `distance` bp are merged"
  }
}
```

## 21) plot

- **作用：** 把 GWAS p 值的曼哈顿图或 QQ 图绘制为 SVG 字符串。

- **处理方法：** 点来自 `compute_manhattan_data` 或 `compute_qq_data`。曼哈顿图：染色体交替着色，用虚线标出全基因组（5e-8）与提示性（1e-5）阈值，并在标签中心标注染色体。QQ 图：观察与期望 -log10 p 的散点及 y = x 直线。使用纯标准库生成 SVG；未给出 `path` 时不写文件。

- **API：** `dnakit.gwas.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`kind` 为 `manhattan`（默认；需要 `chrom`、`pos`、`p`）或 `qq`（需要 `p`）；`width`（900）、`height`（450）、`max_points`（20000）、`path`（已有文件不会被覆盖，除非 `overwrite=True`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `kind`、`svg`（SVG 文本）、`n_points_drawn` 与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import plot

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = plot(rows, options={"kind": "manhattan"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "manhattan",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"450\" viewBox=\"0 0 900 450\"><title>GWAS manhattan plot</title...",
    "n_points_drawn": 6,
    "path": null
  }
}
```

## 22) random_effect_meta_analysis

- **作用：** 用 DerSimonian-Laird 随机效应荟萃分析合并各研究中逐变异的效应。

- **处理方法：** 输入处理与等位基因对齐同 `fixed_effect_meta_analysis`。令 w = 1/se^2，S1 = sum(w)，S2 = sum(w^2)，Cochran Q 的自由度为 k-1，则 tau^2 = max(0, (Q - df)/(S1 - S2/S1))；合并效应使用权重 1/(se^2 + tau^2)，se = 1/sqrt(这些权重之和)，p 值采用正态理论。Q、其 p 值和 I^2 与固定效应模型相同。研究数很少时 tau^2 估计不稳。

- **API：** `dnakit.gwas.random_effect_meta_analysis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 同 `fixed_effect_meta_analysis`：字典 {研究名: 表}、表的列表，或带 `study` 列的单张表；`options`：`min_studies`（2）、`study_column`（`study`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `records`（每个变异：`variant_id`、`chrom`、`pos`、`effect_allele`、`other_allele`、`n_studies`、`beta`、`se`、`z`、`p`、`tau2`、`q`、`q_p`、`i2`）、`n_variants` 与 `model`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import random_effect_meta_analysis

studies = {
    "A": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.10, "se": 0.03}],
    "B": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.30, "se": 0.04}],
    "C": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.05, "se": 0.05}],
}
result = random_effect_meta_analysis(studies)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "random_effect_meta_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 295,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "85abbef3f8debc7c12e8a2bc7be9bb7a242a534e0d4eee0f114d34a321caef77"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": null,
        "pos": null,
        "effect_allele": "A",
        "other_allele": "G",
        "n_studies": 3,
        "beta": 0.15101036270181903,
        "se": 0.07331952361816618,
        "z": 2.05962007456911,
        "p": 0.03943487558282798,
        "tau2": 0.014487000000000002,
        "q": 20.838751625487642,
        "q_p": 2.9848504713101197e-05,
        "i2": 0.9040249609984399
      }
    ],
    "n_variants": 1,
    "model": "DerSimonian-Laird random effects: tau^2 = max(0, (Q-df)/(S1 - S2/S1)), weights 1/(se^2+tau^2), normal-theory p"
  }
}
```

## 23) read

- **作用：** 把 GWAS 汇总统计表（路径、.gz 文件、分隔文本、记录或列）读入记录。

- **处理方法：** 自动检测分隔符（制表符、逗号、分号或空白），`NA`/空单元格变为 null，数字文本转为数值。默认保留原列名；`standardize=True` 时列被重命名为规范列名（见 `standardize`）。超过 `max_rows` 行的表会被拒绝。

- **API：** `dnakit.gwas.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文件路径（允许 `.gz`）、分隔文本、记录列表或列字典。`options`：`delimiter`、`max_rows`（500000）、`standardize`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `columns`、`n_records` 与 `records`；`standardize=True` 时另含 `column_map`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import read

text = "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\nrs1\t1\t100000\tA\tG\t0.12\t0.02\t2e-9\nrs2\t1\t150000\tC\tT\t0.09\t0.02\tNA\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 89,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ce0498fb00c9e38118484def9a109a673112328311015fc7e7c629675152e6dd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "columns": [
      "SNP",
      "CHR",
      "BP",
      "..."
    ],
    "n_records": 2,
    "records": [
      {
        "SNP": "rs1",
        "CHR": 1,
        "BP": 100000,
        "A1": "A",
        "A2": "G",
        "BETA": 0.12,
        "SE": 0.02,
        "P": 2e-09
      },
      {
        "SNP": "rs2",
        "CHR": 1,
        "BP": 150000,
        "A1": "C",
        "A2": "T",
        "BETA": 0.09,
        "SE": 0.02,
        "P": null
      }
    ]
  }
}
```

## 24) run

- **作用：** 对剂量基因型运行原生的单变异关联扫描（线性或 logistic 回归）。

- **处理方法：** 对每个变异，含截距和协变量：定量表型用 OLS（对剂量系数做 t 检验，自由度 = n - 协变量数 - 2），0/1 表型用 logistic 回归（Wald z 检验）；`model='auto'` 在表型为含两类的 0/1 时选 logistic。基因型或表型缺失的样本按变异丢弃。单态变异或样本过少的 `status` 为 `monomorphic_or_too_few_samples`，低于 `maf_min` 的为 `maf_below_threshold`，拟合失败的为 `model_failed`。这是普通回归扫描：除所给协变量外，不校正亲缘关系或群体结构。不会驱动外部程序（PLINK、REGENIE、GEMMA、GCTA）：把 `backend` 设为其中之一时，调用返回 `not_implemented`，指出缺少的程序；若程序已安装，则说明本版本没有对应的命令行适配器。

- **API：** `dnakit.gwas.run(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'genotypes': 样本 x 变异的剂量 [0, 2]（null 表示缺失）, 'phenotype': [...], 'covariates': 样本 x k（可选）, 'variant_ids', 'chrom', 'pos', 'alt', 'ref', 'sample_ids'}`（至少 5 个有表型的样本）。`options`：`model` 为 `auto`（默认）、`linear` 或 `logistic`；`maf_min`（0）；`backend` 为 `native`（默认），或 `plink`、`plink2`、`plink1.9`、`regenie`、`gemma`、`gcta`（这些返回 `not_implemented`）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `model`（`linear` 或 `logistic`）、`n_samples`、`n_variants`、`records`（每个变异：`variant_id`、`chrom`、`pos`、`effect_allele`、`other_allele`、`beta`、`se`、`z`、`p`、`eaf`、`n`、`status`）与 `test`。指定外部 `backend` 时，封装的 `status` 为 `"not_implemented"`，并带 `required_dependency`（程序名，如 `plink2`）与 `message`；未知的后端名会抛出错误。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import run

genotypes = [[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 1], [1, 0, 2], [2, 1, 0], [0, 0, 1], [1, 2, 2]]
phenotype = [0.1, 0.9, 1.8, 0.2, 1.1, 2.0, 0.0, 1.2]
result = run({"genotypes": genotypes, "phenotype": phenotype,
              "variant_ids": ["v1", "v2", "v3"], "chrom": ["1", "1", "1"], "pos": [100, 200, 300]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "run",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 243,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7e3665d19858b0c6bbf219dda78749a79278644f55aaaebfdafbd34aedfbebcf"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "linear",
    "n_samples": 8,
    "n_variants": 3,
    "records": [
      {
        "variant_id": "v1",
        "chrom": "1",
        "pos": 100,
        "effect_allele": null,
        "other_allele": null,
        "beta": 0.905128205128205,
        "se": 0.05695155982549665,
        "z": 15.892948461843325,
        "p": 3.938267710360089e-06,
        "eaf": 0.4375,
        "n": 8,
        "status": "ok"
      },
      {
        "variant_id": "v2",
        "chrom": "1",
        "pos": 200,
        "effect_allele": null,
        "other_allele": null,
        "beta": -0.12051282051282054,
        "se": 0.3706288819493145,
        "z": -0.32515766142936836,
        "p": 0.7561008935006274,
        "eaf": 0.4375,
        "n": 8,
        "status": "ok"
      },
      {
        "variant_id": "v3",
        "chrom": "1",
        "pos": 300,
        "effect_allele": null,
        "other_allele": null,
        "beta": -0.2897435897435897,
        "se": 0.35467510976489625,
        "z": -0.8169267641466362,
        "p": 0.4451878633327042,
        "eaf": 0.5625,
        "n": 8,
        "status": "ok"
      }
    ],
    "test": "OLS t-test on the dosage coefficient (df = n - covariates - 2)"
  }
}
```

## 25) select_lead_snps

- **作用：** 通过对全基因组显著变异做基于距离的贪心剪枝来选择独立 lead SNP。

- **处理方法：** p <= `p_threshold` 的变异按 p 升序取用；每个成为 lead，并移除同一染色体上 `window` bp 内的其他候选（其数量记为 `n_nearby_significant`）。不使用 LD。lead 按基因组顺序返回；`max_leads` 可提前终止选择。

- **API：** `dnakit.gwas.select_lead_snps(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 需要 `chrom`、`pos` 与 `p`。`options`：`p_threshold`（5e-8）、`window`（500000）、`max_leads`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `leads`（`variant_id`、`chrom`、`pos`、等位基因、`beta`、`se`、`p`、`eaf`、`neglog10p`、`n_nearby_significant`）、`n_leads`、`p_threshold` 与 `window`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import select_lead_snps

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = select_lead_snps(rows, options={"p_threshold": 5e-8, "window": 500000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "select_lead_snps",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "leads": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0
      },
      {
        "variant_id": "rs3",
        "chrom": "1",
        "pos": 900000,
        "effect_allele": "G",
        "other_allele": "A",
        "beta": -0.15,
        "se": 0.025,
        "p": 2e-09,
        "eaf": 0.4,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0
      }
    ],
    "n_leads": 2,
    "p_threshold": 5e-08,
    "window": 500000
  }
}
```

## 26) simulate_gwas_data

- **作用：** 模拟基因型、定量表型及由此得到的关联汇总统计（给定种子时结果确定）。

- **处理方法：** 等位基因频率 f_j 在 [maf_min, maf_max] 上均匀分布，剂量服从 Binomial(2, f_j)（变异间不连锁）。随机选出的 `n_causal` 个变异对标准化剂量的效应为 b_j ~ N(0, heritability/n_causal)；表型 = sum_j b_j * 标准化剂量_j + N(0, 1 - heritability)。变异以 10 kb 间隔分布在 `n_chromosomes` 条染色体上，alt 等位基因为 `A`，ref 为 `G`。汇总统计为 `run` 对模拟数据做的线性回归扫描。使用 Python `random`，给定 `seed` 时结果确定；没有 LD 结构。

- **API：** `dnakit.gwas.simulate_gwas_data(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：None 或下列选项的映射（优先于 `options`）。`options`：`n_samples`（500）、`n_variants`（200）、`n_causal`（5）、`heritability`（0.5）、`seed`（0）、`maf_min`（0.05）、`maf_max`（0.5）、`n_chromosomes`（22）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `genotypes`（可用于 `run` 及其他函数的基因型映射：`genotypes`、`sample_ids`、`variant_ids`、`chrom`、`pos`、`alt`、`ref`）、`phenotype`、`true_effects`（{variant_id: 效应}）、`allele_frequencies`、`summary_statistics`（`run` 的记录）、`seed` 与 `model`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import simulate_gwas_data

result = simulate_gwas_data(None, options={"n_samples": 50, "n_variants": 10, "n_causal": 2, "seed": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "simulate_gwas_data",
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
    "genotypes": {
      "genotypes": [
        [
          0,
          1,
          0,
          "..."
        ],
        [
          0,
          1,
          2,
          "..."
        ],
        "..."
      ],
      "sample_ids": [
        "s1",
        "s2",
        "s3",
        "..."
      ],
      "variant_ids": [
        "snp1",
        "snp2",
        "snp3",
        "..."
      ],
      "chrom": [
        "1",
        "2",
        "3",
        "..."
      ],
      "pos": [
        10000,
        10000,
        10000,
        "..."
      ],
      "alt": [
        "A",
        "A",
        "A",
        "..."
      ],
      "ref": [
        "G",
        "G",
        "G",
        "..."
      ]
    },
    "phenotype": [
      -0.04701902898502952,
      0.9976611520850311,
      -0.618585574379673,
      "..."
    ],
    "true_effects": {
      "snp5": -0.132523395903016,
      "snp8": -0.3711891596172627
    },
    "allele_frequencies": {
      "snp1": 0.11046390985058055,
      "snp2": 0.4313451816217547,
      "snp3": 0.3936985785394763,
      "snp4": 0.16478106158273975,
      "snp5": 0.2729457891913734,
      "snp6": 0.2522709791549322,
      "snp7": 0.3432168377252433,
      "snp8": 0.40492550801098093,
      "snp9": 0.0922368140484057,
      "snp10": 0.06275636443490284
    },
    "summary_statistics": [
      {
        "variant_id": "snp1",
        "chrom": "1",
        "pos": 10000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": -0.06099515839356878,
        "se": 0.18813268750275483,
        "z": -0.3242135069838708,
        "p": 0.7471850214511193,
        "eaf": 0.14,
        "n": 50,
        "status": "ok"
      },
      {
        "variant_id": "snp2",
        "chrom": "2",
        "pos": 10000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": -0.0013527033652641059,
        "se": 0.14615643899646016,
        "z": -0.00925517462352013,
        "p": 0.9926539046237707,
        "eaf": 0.41,
        "n": 50,
        "status": "ok"
      },
      "..."
    ],
    "seed": 1,
    "model": "dosage ~ Binomial(2, f_j); phenotype = sum_j b_j * standardised dosage_j + N(0, 1 - h2), b_j ~ N(0, h2 / n_causal)"
  }
}
```

## 27) standardize

- **作用：** 把 GWAS 汇总统计的列标准化为规范列名，并推出缺失的统计量。

- **处理方法：** 可识别的别名被重命名为 `variant_id`、`chrom`（无 `chr` 前缀；23 -> X，24 -> Y，26 -> MT）、`pos`、`effect_allele` 与 `other_allele`（大写）、`beta`、`se`、`p`、`eaf`、`n`、`z`；存在 `maf`、`info`、`ld_score` 列时予以保留。缺失值会被推出：由比值比得 beta = ln(OR)，z = beta/se，p 由 z 得到（双侧正态），se 由 beta 与 p 得到，缺少时 `variant_id` = `chrom:pos:other:effect`。无法识别的列会被丢弃并给出警告。

- **API：** `dnakit.gwas.standardize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`max_rows`（500000）、`delimiter`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `records`（具有规范列的行）、`column_map`（{原列名: 规范列名}）、`columns` 与 `n_records`；`warnings` 列出被忽略的列。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import standardize

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = standardize(rows[:2])
print(json.dumps(result, indent=2))
```

- **示例输出：**

较长的数组和字符串以 `...` 缩写。

```json
{
  "module": "gwas",
  "operation": "standardize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 249,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a530a3275efbe926d438abe1d12407274b4c4cbfd0513d44cad410c43b350a28"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5
      }
    ],
    "column_map": {
      "SNP": "variant_id",
      "CHR": "chrom",
      "BP": "pos",
      "A1": "effect_allele",
      "A2": "other_allele",
      "BETA": "beta",
      "SE": "se",
      "P": "p",
      "EAF": "eaf",
      "N": "n"
    },
    "columns": [
      "variant_id",
      "chrom",
      "pos",
      "..."
    ],
    "n_records": 2
  }
}
```

## 28) validate

- **作用：** 校验 GWAS 汇总统计表：必需列、取值范围、等位基因和重复变异。

- **处理方法：** 标准化之后报告错误（无变异标识或 `chrom` + `pos`；没有 `p`、`z` 或 `beta` 加 `se`；p 不在 [0, 1]；se <= 0；eaf 不在 [0, 1]；n <= 0；beta/se/z 非有限值；效应等位基因与另一等位基因相同）和警告（非 A/C/G/T 字符串的等位基因、重复的变异 id、缺少等位基因或 eaf 列）。内容问题使 `valid` 为 false，而不是抛出异常；无法读取的表则会抛出异常。

- **API：** `dnakit.gwas.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`max_rows`（500000）、`delimiter`。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `valid`、`n_records`、`errors`（消息列表）、`warnings`（消息列表）与 `column_map`。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import validate

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100, "A1": "A", "A2": "G", "BETA": 0.1, "SE": 0.02, "P": 1.5},
    {"SNP": "rs1", "CHR": 1, "BP": 200, "A1": "C", "A2": "C", "BETA": 0.2, "SE": -0.02, "P": 0.01},
]
result = validate(rows)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 190,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f4830ef0886039ab4da0ddc10686680f90e3be2a0f403baedda6e82657802523"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_records": 2,
    "errors": [
      "1 p-values outside [0, 1]",
      "1 standard errors that are not > 0",
      "1 variants have identical effect and other alleles"
    ],
    "warnings": [
      "1 duplicated variant identifiers",
      "no effect-allele frequency column"
    ],
    "column_map": {
      "SNP": "variant_id",
      "CHR": "chrom",
      "BP": "pos",
      "A1": "effect_allele",
      "A2": "other_allele",
      "BETA": "beta",
      "SE": "se",
      "P": "p"
    }
  }
}
```

## 29) write

- **作用：** 把 GWAS 汇总统计写成使用规范列名的分隔表（tsv 或 csv）。

- **处理方法：** 先标准化表格，再按非空的规范列（或 `columns`）序列化。文本总会返回；仅在给出 `path` 时写文件（`.gz` 会被 gzip 压缩），已有文件不会被覆盖，除非 `overwrite=True`。格式取自 `format`，否则取决于路径后缀（`.csv` 为 csv，其余为 tsv）。

- **API：** `dnakit.gwas.write(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 汇总统计：字典记录列表、等长列的字典，或分隔文本/路径（允许 `.gz`）。列名按不区分大小写和标点的方式匹配 `variant_id`（SNP、rsid）、`chrom`（CHR）、`pos`（BP、POS）、`effect_allele`（A1、EA、ALT）、`other_allele`（A2、NEA、REF）、`beta`、`se`、`p`（P、PVAL）、`eaf`（EAF、FREQ）、`n`、`z` 与 `OR`（转换为 beta = ln OR）；缺失的 `z`、`p`、`se`、`variant_id` 会尽可能推出。 `options`：`path`（或 `output`）、`format`（`tsv` 或 `csv`）、`columns`（子集，规范列名）、`overwrite`（False）。

- **输出：** 返回 DNAKit 封装：`status`（`ok`）、`module`、`operation`、`provenance`（`implementation: native`）、`input_summary`、`warnings` 与 `result`。`result` 含 `format`、`columns`、`n_records`、`text` 与 `path`（未写文件时为 null）。输入格式错误或选项无效时抛出 DNAKit 异常。

- **示例代码：**

```python
import json
from dnakit.gwas import write

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = write(rows[:2], options={"format": "tsv", "columns": ["variant_id", "chrom", "pos", "beta", "p"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "gwas",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 249,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a530a3275efbe926d438abe1d12407274b4c4cbfd0513d44cad410c43b350a28"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "columns": [
      "variant_id",
      "chrom",
      "pos",
      "beta",
      "p"
    ],
    "n_records": 2,
    "text": "variant_id\tchrom\tpos\tbeta\tp\nrs1\t1\t100000\t0.12\t2e-09\nrs2\t1\t150000\t0.09\t6.8e-06\n",
    "path": null
  }
}
```
