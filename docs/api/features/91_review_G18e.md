# 群体遗传复审包装函数：dadi 与 tskit（第 91 组，G18e）

这一页记录第 91 组（G18e）对 dadi 与 tskit 中标记为 `other` 的函数的复审结果。dadi 中可离线运行、与 DNA 相关的函数已包装在 `dnakit.population` 中；其余函数保留了精确的状态与原因。tskit 的七个函数都不是 DNA 级操作，`dnakit.population` 不导出任何函数。

这些包装函数从 `dnakit.population` 导出，示例从该包导入。依赖都是可选的，只有调用对应函数时才会导入。输出文件从不覆盖已有文件；模型只能从固定名单中按名称选择，不执行任意代码；调用前后恢复全局随机状态。

## 示例数据

下面的文件在运行示例前会写入工作目录。其余示例数据在代码中由 dadi 模型加泊松噪声生成，随机种子固定。

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
chr1	100	.	A	G	.	PASS	AA=A	GT	0|0	0|1	1|1	0|0
chr1	200	.	C	T	.	PASS	AA=C	GT	0|1	1|1	0|0	0|0
chr1	300	.	G	A	.	PASS	AA=G	GT	0|0	0|0	0|1	1|1
chr1	400	.	T	C	.	PASS	AA=T	GT	1|0	0|1	0|0	1|1
```

**示例文件：** `popinfo.txt`

```text
s1 popA
s2 popA
s3 popB
s4 popB
```

**示例文件：** `slim_p1.txt`

```text
#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 2 m1 250 0.2
2 5 m2 400 0.0
Genomes:
p1:0 A 0 1
p1:1 A 1
p1:2 A
```

**示例文件：** `slim_p2.txt`

```text
#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 3 m1 300 0.2
Genomes:
p2:0 A 0
p2:1 A 1
p2:2 A 0 1
```

## 1) `dnakit.population.total_substitution_rate()`

- **作用：** 计算 DNA 马尔可夫替换模型的总瞬时替换速率。
- **处理方法：** 调用 dadi 的 `Misc.total_instantaneous_rate()`。对角线被忽略；结果为各状态频率乘以离开该状态的速率之和。
- **API：** `dnakit.population.total_substitution_rate(rate_matrix, stationary, *, verbose=True)`
- **输入：** `rate_matrix` 为 4、16 或 64 个状态的方阵（碱基、二核苷酸或三核苷酸，按 dadi 的 CGTA 顺序），`rate_matrix[i][j]` 为从状态 i 到 j 的速率；`stationary` 为各状态频率，和必须为 1。
- **输出：** 含 `num_states`、`total_rate` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import total_substitution_rate

jc = [[0.0 if i == j else 1 / 3 for j in range(4)] for i in range(4)]
result = total_substitution_rate(jc, [0.25] * 4)
print(result["num_states"], round(result["total_rate"], 6))
```

- **示例输出：**

```text
4 1.0
```

## 2) `dnakit.population.ms_command_line()`

- **作用：** 生成 `ms` 系合并模拟的命令行文本；只构造文本，不运行 `ms`。
- **处理方法：** 调用 dadi 的 `Misc.ms_command()`。设置重组时，重组位点默认为 `10 * theta`。
- **API：** `dnakit.population.ms_command_line(theta, sample_sizes, core, iterations, *, recombination=0.0, recombination_sites=None, seeds=None, verbose=True)`
- **输入：** `theta` 为正数；`sample_sizes` 为每个群体的单倍型样本量；`core` 为单行的群体历史命令片段；`iterations` 为模拟次数（最多 100000）；`seeds` 为三个整数或 None。
- **输出：** 含 `command`、`theta`、`sample_sizes`、`iterations` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import ms_command_line

result = ms_command_line(2.0, [3, 4], "-eN 0.5 2.0", 10, recombination=1.0, seeds=(1, 2, 3))
print(result["command"])
```

- **示例输出：**

```text
ms 7 10 -t 2.000000 -I 2 3 4 -eN 0.5 2.0 -r 1.000000 20 -seeds 1 2 3
```

## 3) `dnakit.population.ms_core_command_2d()`

- **作用：** 生成两群体模型的 `ms` 核心命令（群体历史部分）。
- **处理方法：** 调用 dadi 的 `Demographics2D.*_mscore()`。时间与速率采用 dadi 的尺度，而非合并单位。
- **API：** `dnakit.population.ms_core_command_2d(model, params, *, verbose=True)`
- **输入：** `model` 为 `split_mig`（参数 nu1, nu2, T, m）、`IM`（s, nu1, nu2, T, m12, m21）或 `IM_pre`（nuPre, TPre, s, nu1, nu2, T, m12, m21）；`s` 须在 0 与 1 之间。
- **输出：** 含 `model`、`core` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import ms_core_command_2d

print(ms_core_command_2d("split_mig", [1.0, 2.0, 0.5, 0.25])["core"])
```

- **示例输出：**

```text
-n 1 1.000000 -n 2 2.000000 -ma x 0.500000 0.500000 x -ej 0.250000 2 1 -en 0.250000 1 1
```

## 4) `dnakit.population.expected_neutral_sfs_2d()`

- **作用：** 两群体无分化时的中性期望频率谱。
- **处理方法：** 调用 dadi 的 `Demographics2D.snm()`。两个群体从不分化，各自是大小恒定的中性样本。
- **API：** `dnakit.population.expected_neutral_sfs_2d(sample_sizes, pts, *, verbose=True)`
- **输入：** `sample_sizes` 为两个样本量（各至少 2，最多 200）；`pts` 为积分网格点数（20 到 400）。
- **输出：** 含 `sample_sizes`、`pts`、`spectrum` 和 `provenance` 的字典；dadi 掩码的角点返回 None。
- **示例代码：**

```python
from dnakit.population import expected_neutral_sfs_2d

result = expected_neutral_sfs_2d((3, 3), 20)
print(result["sample_sizes"], result["spectrum"][0][0], round(result["spectrum"][0][1], 6))
```

- **示例输出：**

```text
[3, 3] None 0.552742
```

## 5) `dnakit.population.trinucleotide_misid_table()`

- **作用：** 写出祖先状态误判校正所用的 1 - fux 表。
- **处理方法：** 调用 dadi 的 `Misc.make_fux_table()`，内部使用矩阵指数。输入的速率矩阵先复制，调用方的数组不会被修改。输出文件须为新文件。
- **API：** `dnakit.population.trinucleotide_misid_table(path, *, expected_substitutions, rate_matrix, trinucleotide_frequencies, verbose=True)`
- **输入：** `path` 为新文件路径；`expected_substitutions` 为正数，表示内外群之间每位点的期望替换数；`rate_matrix` 为 64 x 64 三核苷酸速率矩阵（dadi 的 CGTA 顺序）；`trinucleotide_frequencies` 为 64 个三核苷酸（A、C、G、T）的频率。
- **输出：** 含 `path`、`rows`（固定为 256）、`expected_substitutions` 和 `provenance` 的字典。
- **示例代码：**

```python
from itertools import product
from pathlib import Path
from dnakit.population import trinucleotide_misid_table

keys = ["".join(p) for p in product("ACGT", repeat=3)]
rates = [[0.0 if i == j else 1 / 63 for j in range(64)] for i in range(64)]
Path("misid_table.txt").unlink(missing_ok=True)
result = trinucleotide_misid_table("misid_table.txt", expected_substitutions=0.1, rate_matrix=rates,
    trinucleotide_frequencies={key: 1 / 64 for key in keys})
print(result["rows"], Path(result["path"]).read_text().splitlines()[1])
```

- **示例输出：**

```text
256 CCC G 0.001626
```

## 6) `dnakit.population.slim_output_data_dict()`

- **作用：** 读取 SLiM 样本输出（每个群体一个文件），得到 dadi 的 SNP 字典。
- **处理方法：** 调用 dadi 的 `Misc.dd_from_SLiM_files()`。群体按文件顺序编号为 0、1、……；位点键为 `<前缀><位置>.<全局突变编号>`；结果按键排序。
- **API：** `dnakit.population.slim_output_data_dict(paths, *, mutation_types=None, prefix="SLIM_", verbose=True)`
- **输入：** `paths` 为每个群体一个 SLiM 输出文件（最多 10 个）；`mutation_types` 为要保留的突变类型列表，None 表示全部保留。
- **输出：** 含 `populations`、`sample_sizes`、`num_sites` 和 `sites` 的字典；每个位点的 `calls` 为各群体的 [祖先计数, 衍生计数]。
- **示例代码：**

```python
from dnakit.population import slim_output_data_dict

result = slim_output_data_dict(["slim_p1.txt", "slim_p2.txt"])
print(result["sample_sizes"], result["num_sites"])
print(result["sites"][0])
print(slim_output_data_dict(["slim_p1.txt", "slim_p2.txt"], mutation_types=["m2"])["num_sites"])
```

- **示例输出：**

```text
[3, 3] 4
{'key': 'SLiM_100.0', 'outgroup_allele': 0, 'calls': {'0': [2, 1], '1': [1, 2]}}
1
```

## 7) `dnakit.population.vcf_subsample_bootstrap_spectra()`

- **作用：** 从 VCF 中随机抽取个体子集并按染色体片段自助抽样，得到频率谱。
- **处理方法：** 调用 dadi 的 `Misc.bootstraps_subsample_vcf()`。`seed` 同时固定个体子集与片段抽样；调用前后保存并恢复 Python 与 NumPy 的全局随机状态。
- **API：** `dnakit.population.vcf_subsample_bootstrap_spectra(vcf_path, popinfo_path, *, subsample, pop_ids, num_boot, chunk_size, seed=0, filter_sites=True, polarized=True, verbose=True)`
- **输入：** `subsample` 为每个群体保留的二倍体个体数；`pop_ids` 为一到三个群体；`num_boot` 最多 200；`chunk_size` 为片段长度（碱基数）。
- **输出：** 含 `num_boot`、`seed`、`spectra`（每个自助谱为嵌套列表，角点为 None）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.population import vcf_subsample_bootstrap_spectra

result = vcf_subsample_bootstrap_spectra("variants.vcf", "popinfo.txt",
    subsample={"popA": 2, "popB": 2}, pop_ids=["popA", "popB"], num_boot=2, chunk_size=200, seed=3)
print(result["num_boot"], len(result["spectra"]), result["spectra"][0][0])
```

- **示例输出：**

```text
2 2 [None, 0.0, 0.0, 1.0, 0.0]
```

## 8) `dnakit.population.parameter_uncertainty_gim()`

- **作用：** 用 Godambe 信息矩阵（GIM）估计参数标准差。
- **处理方法：** 调用 dadi 的 `Godambe.GIM_uncert()`。模型按名称从固定名单中选取，并在 `grid_points` 上外推。`multinom=True` 时同时估计 theta，其不确定性位于最后一项。
- **API：** `dnakit.population.parameter_uncertainty_gim(data, bootstraps, *, model, p0, grid_points=(20, 30, 40), log=False, multinom=True, eps=0.01, verbose=True)`
- **输入：** `data` 为频率谱（嵌套列表，角点可为 None）；`bootstraps` 为与 `data` 形状相同的自助谱列表；`p0` 为正的起始参数。
- **输出：** 含 `params`、`uncertainties`、`num_boot` 和 `provenance` 的字典。
- **示例代码：**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import parameter_uncertainty_gim

result = parameter_uncertainty_gim(data, boots, model="two_epoch", p0=[1.5, 0.5])
print([round(value, 3) for value in result["uncertainties"]])
```

- **示例输出：**

```text
[0.112, 0.109, 6.849]
```

## 9) `dnakit.population.godambe_matrices()`

- **作用：** 得到 Godambe 信息矩阵、Hessian 矩阵和变异矩阵。
- **处理方法：** 调用 dadi 的 `Godambe.get_godambe()`，使用 dadi 的 Poisson 对数似然。模型按原样计算；若需估计 theta，应把它放进参数。Godambe 信息为 H J⁻¹ H。
- **API：** `dnakit.population.godambe_matrices(data, bootstraps, *, model, p0, grid_points=(20, 30, 40), log=False, eps=0.01, verbose=True)`
- **输入：** 与 `parameter_uncertainty_gim()` 相同。
- **输出：** 含 `godambe_information`、`hessian`、`variability` 和 `provenance` 的字典。
- **示例代码：**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import godambe_matrices

result = godambe_matrices(data, boots, model="two_epoch", p0=[1.5, 0.5])
print(len(result["godambe_information"]), round(result["godambe_information"][0][0], 3))
```

- **示例输出：**

```text
2 87951.897
```

## 10) `dnakit.population.nested_model_statistic()`

- **作用：** 带自助调整的嵌套模型检验：似然比检验（LRT）调整因子、调整后的 Wald 统计量与调整后的得分统计量。
- **处理方法：** 调用 dadi 的 `Godambe.LRT_adjust()`、`Wald_stat()` 与 `score_stat()`。`nested_indices` 为嵌套参数在复杂模型参数列表中的位置。嵌套值不可为零，否则矩阵奇异。
- **API：** `dnakit.population.nested_model_statistic(data, bootstraps, *, test, model, p0, nested_indices, grid_points=(20, 30, 40), full_params=None, multinom=True, eps=0.01, adjusted_and_original=False, verbose=True)`
- **输入：** `test` 为 `lrt`、`wald` 或 `score`；`p0` 为复杂模型的参数，嵌套参数取其嵌套值；Wald 检验还需 `full_params`，即复杂模型参数的完整列表。
- **输出：** `lrt` 返回 `adjustment`；`wald` 与 `score` 返回 `statistic`，`adjusted_and_original=True` 时另有 `unadjusted`。
- **示例代码：**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.three_epoch)
mean = np.ma.filled(model((2.0, 0.5, 0.2, 0.3), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(3)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import nested_model_statistic

result = nested_model_statistic(data, boots, test="score", model="three_epoch",
    p0=[1.0, 0.5, 0.2, 0.3], nested_indices=[0], adjusted_and_original=True)
print(round(result["statistic"], 6), round(result["unadjusted"], 6))
```

- **示例输出：**

```text
0.906362 2.542513
```

## 11) `dnakit.population.fit_demographic_model_to_residual()`

- **作用：** 拟合人口模型，使其残差谱逼近给定的目标残差谱。
- **处理方法：** 调用 dadi 的 `Inference.optimize_log_resid()`，对参数取对数后做 BFGS 优化；迭代次数不超过 `max_iterations`。
- **API：** `dnakit.population.fit_demographic_model_to_residual(data, target_residual, *, model, p0, grid_points=(20, 30, 40), lower_bound=None, upper_bound=None, multinom=True, max_iterations=50, verbose=True)`
- **输入：** `target_residual` 为目标残差谱，形状须与 `data` 相同（例如来自同义位点的拟合）；`p0` 为正数。
- **输出：** 含 `params`、`objective`、`converged`、`function_calls`、`gradient_calls` 和 `provenance` 的字典。
- **示例代码：**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import fit_demographic_model_to_residual

target = [value * 0.5 for value in data]
result = fit_demographic_model_to_residual(data, target, model="two_epoch", p0=[0.8, 0.3], max_iterations=2)
print([round(value, 1) for value in result["params"]], result["converged"])
```

- **示例输出：**

```text
[4060.9, 29915.2] False
```

## 12) tskit：没有包装函数

标记为 `other` 的七个 tskit 函数都不是 DNA 级操作，因此本模块不导出任何函数：

| 函数 | 状态 | 原因 |
|---|---|---|
| `tskit.naturalsize` | not-dna | 显示辅助函数，格式化字节数 |
| `tskit.format_number` | not-dna | 显示辅助函数，格式化数字 |
| `tskit.set_print_options` | not-dna | 修改全局打印选项，仅用于显示 |
| `tskit.jit.numba.jitwrap` | internal | 把树序列转换为 numba 可用的形式，不是 DNA 级操作 |
| `tskit.util.naturalsize`、`format_number`、`set_print_options` | not-dna | 与上面同名的函数，重复项 |

## 未包装的 dadi 函数

- **TwoLocus** 与 **Triallele**：在本环境中无法运行。Triallele 使用已被 NumPy 2 移除的 `np.float`（状态 `broken-upstream`）。TwoLocus 报错为 `dadi.tridiag` 不存在（dadi 2.4.4 把该模块命名为 `tridiag_cython`），这不是 NumPy 2 问题。
- **DFE 混合模型**：需要预先计算的 Cache1D/Cache2D 选择缓存（`needs-data`）；**Demes.Inference.optimize** 需要 dadi 的推断选项文件（`needs-data`）；**LowPass.compute_cov_dist** 需要数据字典中的逐位点等位基因深度（`needs-data`）。
- **LowPass** 其余函数为低覆盖校正的中间步骤，**Demes** 图编辑与输出函数依赖模块级全局状态，**Numerics.make_anc_state_misid_func** 与 **Inference.add_misid_param** 以可调用对象为参数（后者已弃用且调用时报错）。
