# 可选 Hail 后端：MatrixTable 导入、质控与群体分析

这一组函数是 `dnakit.variants` 中基于 [hail](https://hail.is/) 的可选适配器。它们把 VCF 导入为 Hail MatrixTable，计算位点与样本质控指标、Hardy-Weinberg 精确检验、连锁不平衡剪枝、基因型 PCA，以及无协变量的线性或逻辑回归。函数只处理显式给出的 VCF、MatrixTable 目录或表型文件，不访问网络，不调用变异检测程序，也不读取 BGEN、PLINK 或 gVCF 输入。

依赖是可选的，只在调用函数时导入。第一次调用会启动本地 Spark 会话，需要 Java；本页示例在 Java 17 下运行，Hail 会打印 Java 11 兼容性提示，该提示输出到 stderr。

## 示例数据

下面的示例都使用这两个文件。VCF 有 4 个样本、3 个双等位位点和 1 个多等位位点；多等位位点在需要双等位位点的分析中被跳过并计数。数据只用于演示接口，统计值不代表真实群体。

**示例文件：** `sample.vcf`

```text
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	10	rs1	A	G	50	PASS	.	GT	0/1	1/1	0/0	./.
chr1	20	rs2	C	T	.	PASS	.	GT	0|1	0/0	0/1	1/1
chr1	30	rs3	G	C	99	PASS	.	GT	0/0	0/1	1/1	0/1
chr2	5	rs4	A	C,T	60	PASS	.	GT	0/1	1/2	./.	0/0
```

**示例文件：** `phenotype.tsv`

```text
sample	y	case
S1	1.0	0
S2	2.5	1
S3	3.0	0
S4	2.0	1
```

## 1) `dnakit.variants.import_vcf_matrix_table()`

- **作用：** 把 VCF 导入为 Hail MatrixTable，并写入一个新目录。
- **处理方法：** 调用 `hl.import_vcf`，参考基因组由 `reference_genome` 指定（`GRCh37` 或 `GRCh38`）。输出路径必须不存在，写入时不覆盖已有内容。
- **API：** `dnakit.variants.import_vcf_matrix_table(source, output_path, *, reference_genome="GRCh38")`。
- **输入：** VCF 文件路径（普通或 bgzip 压缩）；`output_path` 为尚不存在的路径。
- **输出：** `n_variants`、`n_samples`、写出的 `path` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import import_vcf_matrix_table

target = Path(tempfile.mkdtemp()) / "sample.mt"
result = import_vcf_matrix_table("sample.vcf", target)
print(result["n_variants"], result["n_samples"])
```

- **示例输出：**

```text
4 4
```

## 2) `dnakit.variants.summarize_vcf_matrix_table()`

- **作用：** 汇总 VCF 或 MatrixTable 目录中的位点数、双等位与 SNV 数量、染色体分布和样本名。
- **处理方法：** 用 `hl.agg.count_where` 与 `hl.agg.counter` 聚合；染色体按名称排序，样本名按字典序排序。染色体数超过 `max_contigs` 或样本数超过 `max_samples` 时报错。
- **API：** `dnakit.variants.summarize_vcf_matrix_table(source, *, reference_genome="GRCh38", max_contigs=100, max_samples=100000)`。
- **输入：** VCF 文件路径，或由 `import_vcf_matrix_table()` 写出的 MatrixTable 目录。
- **输出：** `n_variants`、`n_samples`、`biallelic_variants`、`multiallelic_variants`、`snv_variants`、`contigs` 和 `samples`。
- **示例代码：**

```python
from dnakit.variants import summarize_vcf_matrix_table

result = summarize_vcf_matrix_table("sample.vcf")
print(result["biallelic_variants"], result["multiallelic_variants"], result["contigs"])
```

- **示例输出：**

```text
3 1 [{'contig': 'chr1', 'variants': 3}, {'contig': 'chr2', 'variants': 1}]
```

## 3) `dnakit.variants.calculate_variant_qc()`

- **作用：** 计算每个双等位位点的 call rate、等位基因频率和杂合、纯合计数。
- **处理方法：** 先去掉多等位位点，再调用 `hl.variant_qc`。等位基因频率为 `AC / AN`，缺失基因型不计入 `AN`。Hardy-Weinberg 检验单独由 `calculate_hardy_weinberg_pvalues()` 完成。
- **API：** `dnakit.variants.calculate_variant_qc(source, *, reference_genome="GRCh38", max_variants=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录；`max_variants` 为双等位位点数的上限。
- **输出：** `variants` 列表（每项含 `chrom`、`pos`、`ref`、`alt`、`n_called`、`call_rate`、`allele_count`、`allele_number`、`allele_frequency`、`n_het`、`n_hom_alt`）和 `skipped_multiallelic`。
- **示例代码：**

```python
from dnakit.variants import calculate_variant_qc

result = calculate_variant_qc("sample.vcf")
first = result["variants"][0]
print(first["pos"], first["call_rate"], first["allele_frequency"], first["n_het"])
```

- **示例输出：**

```text
10 0.75 0.5 1
```

## 4) `dnakit.variants.calculate_sample_qc()`

- **作用：** 计算每个样本的 call rate、基因型类别计数和转换与颠换比。
- **处理方法：** 对全部位点调用 `hl.sample_qc`。转换与颠换按等位基因拷贝计数，因此纯合替代基因型计两次。没有 SNV 时 `r_ti_tv` 为 `None`。
- **API：** `dnakit.variants.calculate_sample_qc(source, *, reference_genome="GRCh38", max_samples=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录。
- **输出：** `samples` 列表（按样本名排序，每项含 `sample`、`call_rate`、`n_called`、`n_hom_ref`、`n_het`、`n_hom_var`、`r_het_hom_var`、`n_transition`、`n_transversion`、`r_ti_tv`）。
- **示例代码：**

```python
from dnakit.variants import calculate_sample_qc

result = calculate_sample_qc("sample.vcf")
print([(row["sample"], row["call_rate"]) for row in result["samples"]])
```

- **示例输出：**

```text
[('S1', 1.0), ('S2', 1.0), ('S3', 0.75), ('S4', 0.75)]
```

## 5) `dnakit.variants.calculate_hardy_weinberg_pvalues()`

- **作用：** 对每个双等位位点检验 Hardy-Weinberg 平衡。
- **处理方法：** 统计纯合参考、杂合和纯合替代的基因型数，再调用 `hl.hardy_weinberg_test`。默认是双侧精确检验，使用中点 p 值校正（Levene–Haldane 分布）；`one_sided=True` 时检验杂合过剩。
- **API：** `dnakit.variants.calculate_hardy_weinberg_pvalues(source, *, reference_genome="GRCh38", one_sided=False, max_variants=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录；`one_sided` 为布尔值。
- **输出：** `variants` 列表（每项含基因型计数、`het_freq_hwe`、`p_value`）、`alternative`（`two_sided` 或 `excess_heterozygosity`）和 `skipped_multiallelic`。
- **示例代码：**

```python
from dnakit.variants import calculate_hardy_weinberg_pvalues

result = calculate_hardy_weinberg_pvalues("sample.vcf")
print([(row["pos"], round(row["p_value"], 4)) for row in result["variants"]])
```

- **示例输出：**

```text
[(10, 0.7), (20, 0.6571), (30, 0.6571)]
```

## 6) `dnakit.variants.prune_linkage_disequilibrium()`

- **作用：** 在连锁不平衡窗口内去掉相关的位点，返回保留的位点集合。
- **处理方法：** 调用 `hl.ld_prune`，只使用双等位位点的二倍体基因型，忽略相位。窗口内两两 r² 严格小于 `r2` 的位点被保留，结果是极大独立集合，不一定是最大集合。
- **API：** `dnakit.variants.prune_linkage_disequilibrium(source, *, reference_genome="GRCh38", r2=0.2, bp_window_size=1000000, max_variants=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录；`r2` 在 (0, 1] 内；`bp_window_size` 为正整数（bp）。
- **输出：** `n_input_variants`、`n_retained`、`skipped_multiallelic` 和按坐标排序的 `retained_variants` 列表。
- **示例代码：**

```python
from dnakit.variants import prune_linkage_disequilibrium

result = prune_linkage_disequilibrium("sample.vcf", r2=0.2)
print(result["n_input_variants"], [row["pos"] for row in result["retained_variants"]])
```

- **示例输出：**

```text
3 [10]
```

## 7) `dnakit.variants.calculate_hwe_normalized_pca()`

- **作用：** 对样本基因型做 Hardy-Weinberg 标准化 PCA，返回特征值和样本得分。
- **处理方法：** 调用 `hl.hwe_normalized_pca`。Hail 会去掉全部为纯合参考或全部为纯合替代的位点。样本得分的符号由 SVD 决定，比较不同运行时应比较绝对值或样本间的相对关系。
- **API：** `dnakit.variants.calculate_hwe_normalized_pca(source, *, reference_genome="GRCh38", n_components=10, max_variants=100000, max_samples=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录；`n_components` 为正整数，不能超过样本数。
- **输出：** `eigenvalues`、`samples`（每项含 `sample` 和 `scores`）、`n_variants`、`n_samples` 和 `skipped_multiallelic`。
- **示例代码：**

```python
from dnakit.variants import calculate_hwe_normalized_pca

result = calculate_hwe_normalized_pca("sample.vcf", n_components=2)
print([round(value, 4) for value in result["eigenvalues"]])
first = result["samples"][0]
print(first["sample"], [round(value, 4) for value in first["scores"]])
```

- **示例输出：**

```text
[2.2761, 1.3333]
S1 [-0.4082, 0.5774]
```

## 8) `dnakit.variants.calculate_association_tests()`

- **作用：** 对每个双等位位点检验一个表型与等位基因剂量的关联，模型只含截距项，没有协变量。
- **处理方法：** 表型表在 Python 中读取，并通过 `hl.literal` 与样本对齐；数值表型用线性回归（`hl.linear_regression_rows`，t 检验），0/1 表型用逻辑回归（`hl.logistic_regression_rows`，Wald 检验）。缺失的剂量由 Hail 用均值填补。表型缺失的样本不参与计算。
- **API：** `dnakit.variants.calculate_association_tests(source, phenotype_path, *, phenotype_column, model="linear", sample_column="sample", reference_genome="GRCh38", max_samples=100000, max_variants=100000)`。
- **输入：** VCF 文件路径或 MatrixTable 目录；表型文件是制表符分隔的文本，必须包含 `sample_column` 和 `phenotype_column`；`model` 取 `"linear"` 或 `"logistic"`。
- **输出：** `results` 列表（每项含 `beta`、`standard_error`、`statistic`、`p_value`，逻辑回归还有 `converged`）、`n_samples_with_phenotype`、`statistic`（`t` 或 `z`）、`covariates`（固定为空列表）和 `skipped_multiallelic`。
- **示例代码：**

```python
from dnakit.variants import calculate_association_tests

result = calculate_association_tests("sample.vcf", "phenotype.tsv", phenotype_column="y")
for row in result["results"]:
    print(row["pos"], row["beta"], round(row["p_value"], 4))
```

- **示例输出：**

```text
10 -0.25 0.761
20 -0.25 0.761
30 1.0 0.0438
```

## 验证方式

测试文件是 `others/tests/unit/test_hail_backends.py`，共 14 项。在安装了 hail 的环境中全部运行；在没有 hail 的环境中，7 项参数校验与缺包测试运行，其余 7 项因缺少 hail 而跳过，并给出原因。

- 每个需要 Hail 的函数都有一项在真实 Hail 上运行的测试。
- 位点与样本的 call rate、等位基因频率、杂合与纯合计数，按样本基因型逐个数出；样本 S2 的转换与颠换按等位基因拷贝计数。
- Hardy-Weinberg 精确检验与用分数计算的 Levene–Haldane 中点 p 值一致，单侧检验也一致。
- 线性回归的 beta、标准误和 p 值与闭式公式一致（四个样本、两个参数，t 分布自由度为 2）。
- 逻辑回归的位点 20 和 30 与手工结果一致（beta 为 0，标准误为 sqrt(2)，p 为 1）。位点 10 的数据存在准完全分离，没有有限的极大似然解，因此返回 `converged: false` 和 `beta: None`。
- LD 剪枝结果中两两 r² 均小于阈值，被剔除的位点与保留位点的 r² 不小于阈值，并与直接调用 `hl.ld_prune` 的结果一致。
- PCA 的特征值和样本得分与直接调用 `hl.hwe_normalized_pca` 一致。
- 缺少 hail 时给出 `bio-hail` 安装提示；无效参数（布尔值、越界值、不支持的参考基因组、表型表中的重复样本、非数值表型、逻辑回归的非 0/1 表型、已存在的输出路径）在启动 Spark 之前被拒绝。

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-hail]"   # hail==0.2.139，需要 Java
```

项目地址：[hail](https://hail.is/)。

## 当前范围

- 已封装：VCF 导入与写出 MatrixTable、MatrixTable 汇总、位点与样本质控、Hardy-Weinberg 精确检验（双侧与杂合过剩单侧）、LD 剪枝、HWE 标准化 PCA、无协变量的线性与逻辑回归（Wald 检验）。
- 只处理双等位位点。多等位位点被跳过，并在结果的 `skipped_multiallelic` 中计数；汇总函数同时报告 `multiallelic_variants`。
- 未封装：带协变量的回归、亲缘关系与 GRM（如 `pc_relate`）、gVCF 与 VDS 合并、BGEN 和 PLINK 输入、插补流程，以及任意 Hail 表达式查询。这些需要单独的参数设计和大规模数据验证，本批次未覆盖。
- `max_variants` 与 `max_samples` 限制的是返回给 Python 的结果规模和表型行数，不限制 Spark 的计算量；Hail 仍会扫描完整 VCF。
- 与已有实现的区别：`calculate_hwe_normalized_pca()` 使用 Hail 的 HWE 标准化方法，不同于 `dnakit.genomics.genotype_pca()`（稠密剂量矩阵 PCA）和 `dnakit.variants.calculate_pca()`（scikit-allel 的 Patterson 标准化）。
- 首次调用时 Hail 以 1 GB driver 内存启动本地 Spark，日志写入系统临时目录中的 `dnakit-hail.log`。
