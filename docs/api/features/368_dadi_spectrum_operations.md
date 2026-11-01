# dadi 频谱操作、统计量、低覆盖度模型与扩散积分

`dnakit.population` 中封装 dadi 的 `Spectrum` 方法（屏蔽、投影、选择与重排群体、取对数、展开、重抽样、统计量、ms 与 sfs_code 解析）、`Numerics` 网格外推、`Misc.fragment_data_dict`、`Godambe` 导数、一到五个群体的扩散密度 `Integration`，以及 `dadi.LowPass` 的 GATK 多样本低覆盖度模型的函数。频谱为嵌套列表，被屏蔽项为 `None`。随机操作带 `seed` 并恢复 NumPy 的全局随机状态。需要 `dnakit[bio-dadi]`（测试版本 dadi 2.4.4）。不提供 `Inference.optimize_lbfgsb`：它对任何输入都抛出 `TypeError`（`fmin_l_bfgs_b` 不再接受 `iprint`）。

## 示例数据

**示例文件：** `toy.ms`

```text
ms 4 2 -t 2.0 -I 2 2 2
123 456 789

//
segsites: 2
positions: 0.1 0.5
01
10
01
10

//
segsites: 1
positions: 0.3
1
0
0
1
```

**示例文件：** `snps.txt`

```text
Ingroup	Outgroup	Allele1	A	B	Allele2	A	B	Gene	Position
TCT	TCT	C	4	3	T	2	2	g1	chr1_100
ACA	ACA	C	1	6	A	3	1	g1	chr1_200
GCG	GCG	C	1	6	G	3	1	g2	chr2_50
```

## 1) `dnakit.population.calculate_sfs_statistic()`

- **作用：** 站点频谱的群体遗传统计量。
- **处理方法：** 调用 `Spectrum` 中 `statistic` 指定的方法：`S`（分离位点数）、`Watterson_theta`、`theta_L`（Zeng 等 2006）、`Zengs_E`、`pi`（含 n / (n - 1) 校正的两两差异数）、`Tajima_D` 或 `Fst`（Weir 和 Cockerham，联合频谱）。除 `S` 和 `Fst` 外都需要单群体频谱。一个函数覆盖七条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.calculate_sfs_statistic(spectrum, statistic, *, folded=False, mask_corners=True)`.
- **输入：** `spectrum`（嵌套列表，被屏蔽项为 None）、`statistic`、`folded`、`mask_corners`。
- **输出：** `value` 和 `sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import calculate_sfs_statistic

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
for name in ("S", "Watterson_theta", "pi", "Tajima_D"):
    print(name, round(calculate_sfs_statistic(spectrum, name, verbose=False)["value"], 4))
```

- **示例输出：**

```text
S 28.29
Watterson_theta 10.0001
pi 10.0002
Tajima_D 0.0
```

## 2) `dnakit.population.mask_spectrum_corners()`

- **作用：** 屏蔽在所有样本中都不出现或都出现的项。
- **处理方法：** 调用 `Spectrum.mask_corners`。不提供 `Spectrum.unmask_all`，因为它在 NumPy 2 下对任何输入都抛出 `IndexError`。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.mask_spectrum_corners(spectrum, *, folded=False)`.
- **输入：** `spectrum` 和 `folded`。
- **输出：** 两个角点项为 None 的 `spectrum`、`sample_sizes`、`folded`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import mask_spectrum_corners

joint = [[(i + 1) * (j + 2) % 5 + 0.5 for j in range(5)] for i in range(5)]
result = mask_spectrum_corners(joint, verbose=False)
print(result["spectrum"][0][0], result["spectrum"][4][4], result["spectrum"][0][1])
```

- **示例输出：**

```text
None None 3.5
```

## 3) `dnakit.population.project_spectrum()`

- **作用：** 把频谱投影到更小的样本量。
- **处理方法：** 调用 `Spectrum.project`，按超几何分布重抽样。任何样本量都不能超过当前值。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.project_spectrum(spectrum, sample_sizes, *, folded=False, mask_corners=True)`.
- **输入：** `spectrum`、`sample_sizes`（每个群体一个）、`folded`、`mask_corners`。
- **输出：** 投影后的 `spectrum`、`projected_to`、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import project_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
result = project_spectrum(spectrum, [6], verbose=False)
print(result["sample_sizes"], [None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **示例输出：**

```text
[6] [None, 10.0, 5.0, 3.333, 2.5, 2.0, None]
```

## 4) `dnakit.population.subset_spectrum_populations()`

- **作用：** 保留联合频谱中的部分群体，对其余群体求和。
- **处理方法：** `keep` 调用 `Spectrum.filter_pops`；`marginalize_over` 调用 `Spectrum.marginalize`（dadi 的轴从 0 数起，DNAKit 两者都从 1 起给群体编号）。二者只能给一个。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.subset_spectrum_populations(spectrum, *, keep=None, marginalize_over=None, folded=False, mask_corners=True)`.
- **输入：** 联合 `spectrum` 以及 `keep` 或 `marginalize_over`（群体编号从 1 起）、`folded`、`mask_corners`。
- **输出：** `spectrum`、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import subset_spectrum_populations

joint = [[(i + 1) * (j + 2) % 5 + 0.5 for j in range(5)] for i in range(5)]
kept = subset_spectrum_populations(joint, keep=[2], verbose=False)
summed = subset_spectrum_populations(joint, marginalize_over=[1], verbose=False)
print(kept["sample_sizes"], kept["spectrum"] == summed["spectrum"])
print([None if v is None else round(v, 2) for v in kept["spectrum"]])
```

- **示例输出：**

```text
[4] True
[None, 12.5, 12.5, 2.5, None]
```

## 5) `dnakit.population.log_spectrum()`

- **作用：** 频谱的自然对数。
- **处理方法：** 调用 `Spectrum.log`；非正数的项被屏蔽。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.log_spectrum(spectrum, *, folded=False, mask_corners=True)`.
- **输入：** `spectrum`、`folded`、`mask_corners`。
- **输出：** 取对数后的 `spectrum`、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import log_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
result = log_spectrum(spectrum, verbose=False)
print([None if v is None else round(v, 3) for v in result["spectrum"][:4]])
```

- **示例输出：**

```text
[None, 2.303, 1.609, 1.204]
```

## 6) `dnakit.population.reorder_spectrum_populations()`

- **作用：** 重新排列联合频谱的群体顺序。
- **处理方法：** 调用 `Spectrum.reorder_pops`；`new_order` 按新顺序列出原群体编号（从 1 起）。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.reorder_spectrum_populations(spectrum, new_order, *, folded=False, mask_corners=True)`.
- **输入：** `spectrum`、`new_order`、`folded`、`mask_corners`。
- **输出：** 转置后的 `spectrum`、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import reorder_spectrum_populations

joint = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
result = reorder_spectrum_populations(joint, [2, 1], verbose=False)
print(result["sample_sizes"], result["spectrum"][1])
```

- **示例输出：**

```text
[2, 3] [2.0, 5.0, 8.0, 11.0]
```

## 7) `dnakit.population.unfold_spectrum()`

- **作用：** 展开折叠频谱。
- **处理方法：** 调用 `Spectrum.unfold`，假设每个等位基因同等可能为祖先型：结果为 `(spectrum + 反转的 spectrum) / 2`。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.unfold_spectrum(spectrum, *, mask_corners=True)`.
- **输入：** 折叠 `spectrum`（次要等位基因计数，被屏蔽项为 None）、`mask_corners`。
- **输出：** 展开后的 `spectrum`、`folded`（False）、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import unfold_spectrum

folded = [None, 4.0, 2.0, 1.0, None, None, None]
result = unfold_spectrum(folded, verbose=False)
print(result["folded"], result["spectrum"])
```

- **示例输出：**

```text
False [None, 2.0, 1.0, 1.0, 1.0, 2.0, None]
```

## 8) `dnakit.population.resample_spectrum()`

- **作用：** 对频谱做泊松重抽样或固定位点数重抽样。
- **处理方法：** `mode='poisson'` 调用 `Spectrum.sample`（每项变为以该项为均值的泊松抽样）；`mode='fixed_size'` 调用 `Spectrum.fixed_size_sample`（抽取 `num_sites` 个位点的多项抽样）。`seed` 使抽样可复现，并在结束后恢复 NumPy 的全局随机状态。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.resample_spectrum(spectrum, mode='poisson', *, num_sites=None, only_nonmasked=False, seed=None, folded=False, mask_corners=True)`.
- **输入：** `spectrum`、`mode`、`num_sites`、`only_nonmasked`、`seed`、`folded`、`mask_corners`。
- **输出：** 重抽样后的 `spectrum`、`mode`、`seed`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import resample_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
fixed = resample_spectrum(spectrum, "fixed_size", num_sites=100, seed=3, only_nonmasked=True, verbose=False)
print(sum(v for v in fixed["spectrum"] if v is not None))
again = resample_spectrum(spectrum, "fixed_size", num_sites=100, seed=3, only_nonmasked=True, verbose=False)
print(fixed["spectrum"] == again["spectrum"])
```

- **示例输出：**

```text
100.0
True
```

## 9) `dnakit.population.spectrum_from_simulation_file()`

- **作用：** 由 ms 或 sfs_code 模拟输出构建频谱。
- **处理方法：** `file_format='ms'` 调用 `Spectrum.from_ms_file`（群体来自 ms 命令行或 `pop_assignments`）；`file_format='sfs_code'` 调用 `Spectrum.from_sfscode_file`，`sites` 取 `all`、`syn` 或 `nonsyn`。`average` 为真时对重复运行取平均，否则求和。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.spectrum_from_simulation_file(path, *, file_format='ms', average=True, mask_corners=True, pop_assignments=None, sites='all')`.
- **输入：** `path`、`file_format`、`average`、`mask_corners`、`pop_assignments`、`sites`。
- **输出：** `spectrum`、`sample_sizes`，以及 ms 的 `command` 和 `seeds`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import spectrum_from_simulation_file

result = spectrum_from_simulation_file("toy.ms", verbose=False)
print(result["sample_sizes"], result["command"])
print(result["spectrum"])
```

- **示例输出：**

```text
[2, 2] ms 4 2 -t 2.0 -I 2 2 2
[[None, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, None]]
```

## 10) `dnakit.population.fragment_snp_data()`

- **作用：** 把 dadi SNP 文件切成基因组片段，用于分块自举。
- **处理方法：** 用 `Misc.make_data_dict` 读取 SNP 文件并调用 `Misc.fragment_data_dict`；每个片段包含一条染色体上长度为 `chunk_size` 碱基的窗口内的 SNP。SNP 的 ID 须形如 `染色体_位置`。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.fragment_snp_data(path, chunk_size, *, include_snp_ids=False)`.
- **输入：** `path`、以碱基数计的 `chunk_size`、`include_snp_ids`。
- **输出：** `num_snps`、`num_chunks` 和 `chunks`（编号、SNP 数、可选的 ID）。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import fragment_snp_data

result = fragment_snp_data("snps.txt", 1000, include_snp_ids=True, verbose=False)
print(result["num_snps"], result["num_chunks"])
print([chunk["snp_ids"] for chunk in result["chunks"]])
```

- **示例输出：**

```text
3 2
[['g1_chr1_100', 'g1_chr1_200'], ['g2_chr2_50']]
```

## 11) `dnakit.population.expected_sfs_extrapolated()`

- **作用：** 外推到无限细网格的期望频谱。
- **处理方法：** 在每个网格大小上计算人口历史模型并外推：`log_extrapolation=False` 用 `make_extrap_func`，`True` 用 `make_extrap_log_func`（对对数外推）。模型族（一、二、三个群体）由样本量个数决定。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.expected_sfs_extrapolated(model, params, sample_sizes, *, grid_points=(40, 50, 60), log_extrapolation=False)`.
- **输入：** `model`、`params`、`sample_sizes`、`grid_points`（三个或更多）、`log_extrapolation`。
- **输出：** `spectrum`、`model`、`grid_points`、`sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import expected_sfs_extrapolated

result = expected_sfs_extrapolated("two_epoch", [2.0, 0.1], [6], grid_points=[30, 40, 50], verbose=False)
print([None if v is None else round(v, 4) for v in result["spectrum"]])
```

- **示例输出：**

```text
[None, 1.1253, 0.506, 0.3331, 0.2497, 0.1999, None]
```

## 12) `dnakit.population.log_likelihood_derivatives()`

- **作用：** 模型对数似然关于参数的 Hessian 或梯度。
- **处理方法：** 对数据在外推后的人口历史模型下的对数似然（默认多项，`multinom=False` 时为泊松），按分数步长 `eps` 用有限差分调用 `Godambe.get_hess` 或 `Godambe.get_grad`。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.log_likelihood_derivatives(data, model, params, *, derivative='hessian', grid_points=(20, 30, 40), eps=0.01, multinom=True)`.
- **输入：** `data`、`model`、`params`、`derivative`（`hessian` 或 `gradient`）、`grid_points`、`eps`、`multinom`。
- **输出：** `value`：Hessian 为矩阵，梯度为列表。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import (
    expected_sfs_extrapolated,
    log_likelihood_derivatives,
)

truth = expected_sfs_extrapolated("two_epoch", [2.0, 0.3], [10], grid_points=[20, 30, 40], verbose=False)
data = [0.0 if v is None else 500 * v for v in truth["spectrum"]]
grad = log_likelihood_derivatives(data, "two_epoch", [2.0, 0.3], derivative="gradient", verbose=False)
hess = log_likelihood_derivatives(data, "two_epoch", [2.0, 0.3], verbose=False)
print([round(v, 2) for v in grad["value"]])
print([[round(v, 1) for v in row] for row in hess["value"]])
```

- **示例输出：**

```text
[0.0, 0.0]
[[-16.7, -30.1], [-30.1, -98.3]]
```

## 13) `dnakit.population.integrate_diffusion_phi()`

- **作用：** 对一到五个群体按时间积分扩散密度 phi。
- **处理方法：** 从祖先群体的平衡密度出发（`PhiManip.phi_1D`，X 连锁位点用 `phi_1D_X`），在时间零处分裂成所需数量的群体，再用 `Integration.one_pop`、`two_pops`、`three_pops`、`four_pops`、`five_pops` 或 `one_pop_X` 积分 `duration`（单位为 2N 代），人口大小、选择、显性度为常数并带迁移矩阵。给定 `sample_sizes` 时计算期望频谱（`Spectrum.from_phi`）。一个函数覆盖六条工作清单。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.integrate_diffusion_phi(populations, duration, *, nu=None, migration=None, gamma=None, h=None, theta0=1.0, grid_points=30, x_linked=False, breeding_ratio=1.0, mutation_ratio=1.0, sample_sizes=None, return_phi=False)`.
- **输入：** `populations`（1 到 5）、`duration`，以及可选的 `nu`、`migration`（矩阵，从 j 迁入 i 的速率）、`gamma`、`h`、`theta0`、`grid_points`、`x_linked`、`breeding_ratio`、`mutation_ratio`、`sample_sizes`、`return_phi`。
- **输出：** `phi_shape`、`phi_total`、`grid`，以及可选的 `phi` 和 `spectrum`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import integrate_diffusion_phi

one = integrate_diffusion_phi(1, 0.1, nu=[2.0], grid_points=20, sample_sizes=[6], verbose=False)
print(one["phi_shape"], [None if v is None else round(v, 3) for v in one["spectrum"]])
two = integrate_diffusion_phi(
    2, 0.1, nu=[2.0, 1.0], migration=[[0, 1.0], [0.5, 0]], grid_points=20, sample_sizes=[4, 4], verbose=False
)
print(two["phi_shape"], round(two["spectrum"][1][1], 3))
```

- **示例输出：**

```text
[20] [None, 1.253, 0.545, 0.349, 0.257, 0.203, None]
[20, 20] 0.234
```

## 14) `dnakit.population.simulate_low_coverage_reads()`

- **作用：** 在一种基因型划分上模拟读段计数。
- **处理方法：** 调用 `LowPass.simulate_reads`：每个个体的深度来自其群体的覆盖度分布，杂合子的读段在两个等位基因间按二项分布分配。`seed` 使抽样可复现。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.simulate_low_coverage_reads(coverage_distribution, genotype_partition, sequenced_haplotypes, number_simulations, *, seed=None)`.
- **输入：** `coverage_distribution`（群体到深度 0、1、2、... 概率的映射）、`genotype_partition`（每个个体 0、1 或 2）、`sequenced_haplotypes`、`number_simulations`、`seed`。
- **输出：** `reference_reads` 与 `alternate_reads`：行为模拟位点，列为个体。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import simulate_low_coverage_reads

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = simulate_low_coverage_reads(coverage, [0, 1, 2, 1], [8], 4, seed=2, verbose=False)
print(result["reference_reads"])
print(result["alternate_reads"])
```

- **示例输出：**

```text
[[3, 0, 0, 2], [3, 0, 0, 3], [2, 1, 0, 3], [1, 2, 0, 3]]
[[0, 0, 3, 1], [0, 2, 2, 1], [0, 1, 4, 0], [0, 1, 2, 1]]
```

## 15) `dnakit.population.subsample_genotypes()`

- **作用：** 把基因型调用子抽样到更少的单倍型。
- **处理方法：** 调用 `LowPass.subsample_genotypes_1D`：已调用个体少于 `subsample_haplotypes / 2` 的位点被丢弃，其余位点随机抽取该数量的已调用个体（99 表示缺失）。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.subsample_genotypes(genotype_calls, subsample_haplotypes, *, seed=None)`.
- **输入：** `genotype_calls`（行为位点，列为个体；0、1、2、99）、`subsample_haplotypes`、`seed`。
- **输出：** `genotypes`、`num_loci`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import subsample_genotypes

calls = [[0, 1, 2, 99], [1, 1, 0, 2], [99, 99, 99, 1]]
result = subsample_genotypes(calls, 4, seed=1, verbose=False)
print(result["num_loci"], result["genotypes"])
```

- **示例输出：**

```text
2 [[0, 1], [2, 0]]
```

## 16) `dnakit.population.simulate_gatk_calling_spectrum()`

- **作用：** 给定等位基因计数时，GATK 多样本调用产生的频谱。
- **处理方法：** 调用 `LowPass.simulate_GATK_multisample_calling`：按覆盖度分布模拟读段，像 GATK 一样调用（变异至少需要两条替代读段），再子抽样。第 0 项包含未被调用的位点。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.simulate_gatk_calling_spectrum(coverage_distribution, allele_counts, sequenced_haplotypes, subsample_haplotypes, number_simulations, *, inbreeding=None, seed=None)`.
- **输入：** `coverage_distribution`、`allele_counts`、`sequenced_haplotypes`、`subsample_haplotypes`（每个群体一个）、`number_simulations`、`inbreeding`、`seed`。
- **输出：** `probabilities`：每个观察到的衍生等位基因计数的概率。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import simulate_gatk_calling_spectrum

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = simulate_gatk_calling_spectrum(coverage, [3], [6], [4], 200, seed=1, verbose=False)
print(result["probabilities"], round(sum(result["probabilities"]), 6))
```

- **示例输出：**

```text
[0.095, 0.26, 0.45, 0.17, 0.025] 1.0
```

## 17) `dnakit.population.gatk_no_call_probability()`

- **作用：** GATK 多样本调用不给出基因型的概率。
- **处理方法：** 对每个衍生单倍型计数 0 到 n 调用 `LowPass.probability_of_no_call_1D_GATK_multisample`。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.gatk_no_call_probability(coverage_distribution, sequenced_haplotypes, *, inbreeding=0.0)`.
- **输入：** `coverage_distribution`（深度 0、1、... 的概率）、`sequenced_haplotypes`、`inbreeding`。
- **输出：** `probabilities`，第 k 项对应 k 个衍生单倍型。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import gatk_no_call_probability

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = gatk_no_call_probability(coverage["A"], 6, verbose=False)
print([round(v, 4) for v in result["probabilities"]])
```

- **示例输出：**

```text
[1.0, 0.5148, 0.1618, 0.0453, 0.0122, 0.0033, 0.0009]
```

## 18) `dnakit.population.probability_enough_individuals_covered()`

- **作用：** 有足够多个体被覆盖以得到子样本的概率。
- **处理方法：** 调用 `LowPass.probability_enough_individuals_covered`：对已被判为变异的位点，`sequenced_haplotypes / 2` 个个体中至少 `subsample_haplotypes / 2` 个有读段的概率。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.probability_enough_individuals_covered(coverage_distribution, sequenced_haplotypes, subsample_haplotypes)`.
- **输入：** `coverage_distribution`、`sequenced_haplotypes`、`subsample_haplotypes`。
- **输出：** `probability`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import probability_enough_individuals_covered

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = probability_enough_individuals_covered(coverage["A"], 6, 4, verbose=False)
print(round(result["probability"], 4))
```

- **示例输出：**

```text
0.9975
```

## 19) `dnakit.population.expected_sfs_low_pass()`

- **作用：** 受低覆盖度 GATK 多样本调用扭曲的人口历史模型期望频谱。
- **处理方法：** 用 `LowPass.make_low_pass_func_GATK_multisample` 包装外推后的模型并求值：测序单倍型数的频谱经不调用概率、基因型调用错误和投影变换到子样本大小，低于 `simulation_threshold` 的项用模拟处理。
- **API：** `dnakit.population.dadi_p22_spectrum_functions.expected_sfs_low_pass(model, params, coverage_distribution, sequenced_haplotypes, subsample_haplotypes, *, grid_points=(30, 40, 50), inbreeding=None, simulation_threshold=0.01, simulations_per_frequency=1000, seed=None)`.
- **输入：** `model`、`params`、`coverage_distribution`、`sequenced_haplotypes`、`subsample_haplotypes`、`grid_points`、`inbreeding`、`simulation_threshold`、`simulations_per_frequency`、`seed`。
- **输出：** 子样本大小的 `spectrum` 和 `sample_sizes`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_spectrum_functions import expected_sfs_low_pass

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = expected_sfs_low_pass(
    "two_epoch", [2.0, 0.1], coverage, [6], [4], grid_points=[20, 30, 40],
    simulations_per_frequency=20, seed=1, verbose=False,
)
print(result["sample_sizes"], [None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **示例输出：**

```text
[4] [1.023, 0.48, 0.474, 0.31, 0.115]
```
