# 可选 VCF 后端与群体遗传统计

这一组函数是 `dnakit.variants` 中基于可选依赖的适配器。VCF 读取可以选择 `cyvcf2`、`pysam` 或 `vcfpy`，三者返回完全相同的记录结构；群体遗传统计使用 `scikit-allel`，输入就是这些记录。所有函数只处理显式给出的 VCF 或记录字典，不访问参考数据库，也不调用变异检测程序。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的示例都使用这个 VCF。它只有 5 个位点、4 个样本，只用于演示接口，统计值不代表真实群体。

**示例文件：** `sample.vcf`

```text
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	10	rs1	A	G	50	PASS	DP=20;AF=0.5;DB	GT	0/1	1/1	0/0	./.
chr1	20	.	C	T,A	.	q10	DP=7;AF=0.2,0.1	GT	0|1	0/0	0/2	1/2
chr1	30	rs3	G	C	99	.	DP=12	GT	0/0	0/1	1/1	0/1
chr1	40	rs4	T	A	60	PASS	DP=9;AF=0.4	GT	1/1	0/1	0/0	0/0
chr1	50	rs5	G	A	70	PASS	DP=11;AF=0.3	GT	0/1	0/1	./.	0/0
```

## 1) `dnakit.variants.read_vcf()`

- **作用：** 用指定后端读取 VCF，返回统一的变异记录。
- **处理方法：** 逐条解析位点，把基因型整理为 `"0/1"`、`"0|1"`、`"./."` 的形式。过滤器只保留未通过的过滤名，`PASS` 记为空列表。INFO 中只有一个值的字段记为标量。
- **API：** `dnakit.variants.read_vcf(source, *, backend="cyvcf2", max_records=100000)`。
- **输入：** VCF 文件路径；`backend` 取 `"cyvcf2"`、`"pysam"` 或 `"vcfpy"`；`max_records` 是显式的记录上限，超过时报错。
- **输出：** 含 `samples`、`variants` 和 `provenance` 的字典。每条变异包含 `chrom`、`pos`、`id`、`ref`、`alt`、`qual`、`filter`、`info` 和 `genotypes`。
- **示例代码：**

```python
from dnakit.variants import read_vcf

result = read_vcf("sample.vcf", backend="cyvcf2")
print(result["samples"])
print(result["variants"][0])
```

- **示例输出：**

```text
['S1', 'S2', 'S3', 'S4']
{'chrom': 'chr1', 'pos': 10, 'id': 'rs1', 'ref': 'A', 'alt': ['G'], 'qual': 50.0, 'filter': [], 'info': {'DP': 20, 'AF': 0.5, 'DB': True}, 'genotypes': ['0/1', '1/1', '0/0', './.']}
```

## 2) `dnakit.variants.fetch_vcf_region()`

- **作用：** 只读取指定区域内的变异，不需要扫描整个文件。
- **处理方法：** 调用 `pysam` 或 `cyvcf2` 的 tabix 区域查询。区域坐标是 1-based 闭区间，内部转换为 0-based 半开区间。
- **API：** `dnakit.variants.fetch_vcf_region(source, region, *, backend="pysam", max_records=100000)`。
- **输入：** 经 bgzip 压缩并建立 tabix 或 CSI 索引的 VCF；`region` 形如 `"chr1:1000-2000"`、`"chr1:1000"` 或 `"chr1"`。
- **输出：** 与 `read_vcf()` 相同的结构，并在 `region` 字段记录请求的区域。
- **示例代码：**

```python
import pysam
from dnakit.variants import fetch_vcf_region

pysam.tabix_compress("sample.vcf", "sample.vcf.gz", force=True)
pysam.tabix_index("sample.vcf.gz", preset="vcf", force=True)
result = fetch_vcf_region("sample.vcf.gz", "chr1:15-35", backend="pysam")
print([variant["pos"] for variant in result["variants"]])
```

- **示例输出：**

```text
[20, 30]
```

## 3) `dnakit.variants.calculate_observed_heterozygosity()`

- **作用：** 按群体计算平均观测杂合度。
- **处理方法：** 用 scikit-allel 的 `heterozygosity_observed` 计算每个位点的杂合比例，分母是有效基因型数。再对位点取平均，全部缺失的位点不计入。
- **API：** `dnakit.variants.calculate_observed_heterozygosity(data, *, populations=None)`。
- **输入：** `read_vcf()` 的输出或结构相同的字典；`populations` 是 `{群体名: [样本名, ...]}`，省略时所有样本作为一个群体。
- **输出：** 每个群体的样本数、有效位点数和平均杂合度。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_observed_heterozygosity

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_observed_heterozygosity(
    data, populations={"A": ["S1", "S2"], "B": ["S3", "S4"]}
)
print(result["populations"])
```

- **示例输出：**

```text
{'A': {'samples': 2, 'variants': 5, 'mean_observed_heterozygosity': 0.6}, 'B': {'samples': 2, 'variants': 5, 'mean_observed_heterozygosity': 0.3}}
```

## 4) `dnakit.variants.calculate_fst()`

- **作用：** 计算两个群体之间的全基因组 Fst。
- **处理方法：** 只使用双等位位点。`method="hudson"` 使用 Hudson 估计量（分子和分母分别求和后相除）；`method="weir_cockerham"` 使用 Weir–Cockerham 估计量。标准误由约 10 个块的 block jackknife 得到。Weir–Cockerham 要求每个群体至少有两个样本。
- **API：** `dnakit.variants.calculate_fst(data, *, populations, method="hudson")`。
- **输入：** 数据字典；`populations` 必须恰好包含两个群体；`method` 取 `"hudson"` 或 `"weir_cockerham"`。
- **输出：** Fst、标准误、位点数、块大小和被跳过的位点计数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_fst

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_fst(data, populations={"A": ["S1", "S2"], "B": ["S3", "S4"]})
print(result["fst"], result["standard_error"])
```

- **示例输出：**

```text
0.4920634920634921 0.12077896535034918
```

## 5) `dnakit.variants.calculate_pca()`

- **作用：** 对样本基因型做主成分分析。
- **处理方法：** 只使用双等位、二倍体且无缺失的位点，去掉单态位点后用 scikit-allel 的 PCA，并采用 Patterson 标准化。
- **API：** `dnakit.variants.calculate_pca(data, *, n_components=10, samples=None)`。
- **输入：** 数据字典；`n_components` 为正整数，不能超过样本数或多态位点数；`samples` 可选，用于选择样本子集。
- **输出：** 每个样本的坐标、各主成分的方差解释比例，以及参与计算的位点数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_pca

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_pca(data, n_components=2)
print(result["explained_variance_ratio"], result["polymorphic_variants"])
```

- **示例输出：**

```text
[0.9292512536048889, 0.07074874639511108] 2
```

## 6) `dnakit.variants.calculate_tajimas_d()`

- **作用：** 计算一条染色体上一组样本的 Tajima's D。
- **处理方法：** 只使用双等位位点，调用 scikit-allel 的 `tajima_d`。分离位点少于 3 个时，scikit-allel 返回 NaN，本函数返回 `None`。
- **API：** `dnakit.variants.calculate_tajimas_d(data, *, samples=None)`。
- **输入：** 数据字典，所有位点必须位于同一条染色体；`samples` 可选。
- **输出：** Tajima's D、样本数和位点数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_tajimas_d

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_tajimas_d(data)["tajimas_d"])
```

- **示例输出：**

```text
1.973246148487123
```

## 7) `dnakit.variants.calculate_watterson_theta()`

- **作用：** 计算每个碱基的 Watterson's theta。
- **处理方法：** 由 scikit-allel 的 `watterson_theta` 计算，对应 S / a_n / L。其中 S 是分离位点数，a_n 由染色体数决定，L 为 `sequence_length`。
- **API：** `dnakit.variants.calculate_watterson_theta(data, *, sequence_length, samples=None)`。
- **输入：** 数据字典；`sequence_length` 为正整数，且不能小于最大位点坐标；`samples` 可选。
- **输出：** 每碱基 theta、样本数和位点数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_watterson_theta

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_watterson_theta(data, sequence_length=100)["theta_per_base"])
```

- **示例输出：**

```text
0.015426997245179066
```

## 8) `dnakit.variants.calculate_nucleotide_diversity()`

- **作用：** 计算每个碱基的核苷酸多样性 π。
- **处理方法：** 对每个位点计算 ac × (an − ac) / C(an, 2)，求和后除以 `sequence_length`。计算由 scikit-allel 的 `sequence_diversity` 完成。
- **API：** `dnakit.variants.calculate_nucleotide_diversity(data, *, sequence_length, samples=None)`。
- **输入：** 与 Watterson's theta 相同。
- **输出：** 每碱基 π、样本数和位点数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_nucleotide_diversity

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_nucleotide_diversity(data, sequence_length=100)["pi_per_base"])
```

- **示例输出：**

```text
0.0224047619047619
```

## 9) `dnakit.variants.calculate_site_frequency_spectrum()`

- **作用：** 计算非折叠位点频谱（SFS）。
- **处理方法：** 统计双等位位点中 alt 等位基因的数量。在 `n` 个二倍体样本中，频谱长度为 `2n + 1`。
- **API：** `dnakit.variants.calculate_site_frequency_spectrum(data, *, samples=None)`。
- **输入：** 数据字典；`samples` 可选。
- **输出：** `spectrum` 列表，第 i 项是 alt 等位基因恰好出现 i 次的位点数；另有 `chromosomes` 字段表示染色体总数。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_site_frequency_spectrum

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_site_frequency_spectrum(data)["spectrum"])
```

- **示例输出：**

```text
[0, 0, 1, 2, 1, 0, 0, 0, 0]
```

## 10) `dnakit.variants.calculate_windowed_ld()`

- **作用：** 按固定窗口计算连锁不平衡 r² 的中位数。
- **处理方法：** 使用 scikit-allel 的 `windowed_r_squared`（Rogers 和 Huff 方法）。只使用双等位、二倍体且无缺失的位点。窗口从位置 1 开始，按 `window_size` 划分，最后一个窗口在最大保留位点处结束。
- **API：** `dnakit.variants.calculate_windowed_ld(data, *, window_size=10000, samples=None)`。
- **输入：** 单条染色体的数据字典；`window_size` 以 bp 为单位；`samples` 可选。
- **输出：** 每个窗口的起止坐标（1-based，闭区间）、位点数和中位 r²。
- **示例代码：**

```python
from dnakit.variants import read_vcf, calculate_windowed_ld

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_windowed_ld(data, window_size=100)["windows"])
```

- **示例输出：**

```text
[{'start': 1, 'end': 40, 'variants': 2, 'median_r_squared': 0.7272727489471436}]
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-cyvcf2]"          # cyvcf2，read_vcf / fetch_vcf_region
python -m pip install "dnakit[bio-reads]"           # pysam，read_vcf / fetch_vcf_region
python -m pip install "dnakit[bio-vcfpy]"           # vcfpy，read_vcf
python -m pip install "dnakit[bio-scikit-allel]"    # scikit-allel 与 scipy，calculate_* 群体统计
```

各包的项目地址：[cyvcf2](https://github.com/brentp/cyvcf2)、[pysam](https://github.com/pysam-developers/pysam)、[vcfpy](https://github.com/dahlo/vcfpy)、[scikit-allel](https://github.com/cggh/scikit-allel)。

## 验证方式

测试文件是 `others/tests/unit/test_variant_backends.py`，共 25 项，主要检查：

- 三个 VCF 后端的输出完全一致，过滤器、缺失基因型和相位信息都经过检查。
- Watterson's theta、π、杂合度、SFS 与手工按公式计算的结果一致。
- Tajima's D、Hudson Fst、Weir–Cockerham Fst 和 PCA 与直接调用 scikit-allel 的结果一致。
- 缺少可选依赖时给出安装提示；未知后端、超出 `max_records`、错误的基因型和未知样本名都会被拒绝。

## 当前范围

这一批只覆盖 VCF 读取和基于 scikit-allel 的群体遗传统计。其余 DNA 和基因功能会按领域分批补充，并同步更新文档。
