# dadi 选择效应分布积分谱与三等位基因频谱

`dnakit.population` 中封装 `dadi.DFE`（对选择效应分布积分的频谱）与 `dadi.Triallele`（中性三等位基因密度、抽样，以及三等位基因频谱的折叠、展开和屏蔽）的函数。频谱为嵌套列表，dadi 屏蔽的项为 `None`。需要 `dnakit[bio-dadi]`（测试版本 dadi 2.4.4）。dadi 其余三等位基因函数在上游对任何输入都会失败，因此不提供：人口历史模型调用了 `np.float`，向前积分返回 NaN 或 0，祖先错判、投影和最优缩放会抛出 `NameError` 或 `AttributeError`。

## 1) `dnakit.population.dfe_integrated_spectrum()`

- **作用：** 对选择效应分布（DFE）积分得到的期望频谱。
- **处理方法：** 为带选择的人口历史模型构建 dadi 的选择缓存（`Cache1D` 或 `Cache2D`），再按 `kind` 积分：`negative_1d` 对应 `Cache1D.integrate`，`point_pos_1d` 对应 `integrate_point_pos`，`bivariate_2d` 对应 `Cache2D.integrate`，`bivariate_point_pos_2d` 对应 `integrate_point_pos`，`bivariate_symmetric_point_pos_2d` 对应 `integrate_symmetric_point_pos`。一个函数覆盖五条工作清单，模型与密度函数按名称选择。
- **API：** `dnakit.population.dadi_p22_dfe_triallele_functions.dfe_integrated_spectrum(kind, model, demographic_params, sample_sizes, dfe_params, *, distribution=None, theta=1.0, n_positive=1, grid_points=(20, 25, 30), gamma_points=20, gamma_max=2000.0, exterior_int=True)`.
- **输入：** `kind`、`model`（`dadi.DFE.DemogSelModels` 中的名称）、人口历史参数、一到两个样本量、DFE 参数（点质量类的 `kind` 中最后几个值是正选择点质量），以及可选的 `distribution`（`gamma`、`lognormal`、`exponential`、`beta`、`biv_lognormal`、`biv_ind_gamma`）、`theta`、`grid_points`（三个网格大小）、`gamma_points`、`gamma_max`、`exterior_int`。
- **输出：** `spectrum`（嵌套列表，被屏蔽的项为 None）、`kind`、`model`、`distribution`、`sample_sizes`、缓存设置和 `provenance`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import dfe_integrated_spectrum

result = dfe_integrated_spectrum(
    "negative_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0],
    grid_points=[30, 40, 50], gamma_points=8, gamma_max=100.0, verbose=False,
)
print(result["kind"], result["distribution"], result["sample_sizes"])
print([None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **示例输出：**

```text
negative_1d gamma [6]
[None, 1.123, 0.42, 0.253, 0.179, 0.137, None]
```

## 2) `dnakit.population.triallelic_neutral_equilibrium()`

- **作用：** 中性平衡下的三等位基因密度与频谱。
- **处理方法：** 在网格上计算 Jenkins 等（2014）的精确中性密度（`Triallele.integration.equilibrium_neutral_exact`），再抽样到 `sample_size` 条染色体（`Triallele.numerics.sample`），可按主、次衍生等位基因折叠。这就是中性等位基因下 `Triallele.demographics.equilibrium` 的频谱；传入 `phi` 时直接对该密度抽样。dadi 的向前积分（`advance`，用于选择与分段历史模型）在当前环境对任何输入都返回 NaN 或 0，且模型调用了已被移除的 `np.float`，因此只提供中性平衡。一个函数覆盖三条工作清单。
- **API：** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_neutral_equilibrium(grid_points, sample_size=None, *, folded=False, phi=None, return_density=False)`.
- **输入：** `grid_points`（10 到 120）、可选的 `sample_size`、`folded`、`phi`（方阵嵌套列表）、`return_density`。
- **输出：** `spectrum`（方阵嵌套列表）、`sample_size`、`folded_major`、`folded_ancestral`、`density_total`，以及可选的 `phi`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import triallelic_neutral_equilibrium

result = triallelic_neutral_equilibrium(30, 8, verbose=False)
print(result["sample_size"], result["folded_major"], round(result["density_total"], 1))
print([None if v is None else round(v, 3) for v in result["spectrum"][1][:4]])
```

- **示例输出：**

```text
8 False 12675.3
[None, 0.758, 0.434, 0.291]
```

## 3) `dnakit.population.alternate_mutation_triallelic_spectrum()`

- **作用：** 两个突变同时发生在同一位点时的精确中性三等位基因频谱。
- **处理方法：** 调用 `Triallele.integration.alt_mut_mech_sample_spectrum`，即 Jenkins 等（2014）对另一种突变机制的解。群体大小变化不会改变该中性频谱。
- **API：** `dnakit.population.dadi_p22_dfe_triallele_functions.alternate_mutation_triallelic_spectrum(sample_size)`.
- **输入：** `sample_size`，抽样染色体数（3 到 100）。
- **输出：** `spectrum`（方阵嵌套列表，不可能的项为 None）和 `sample_size`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import alternate_mutation_triallelic_spectrum

result = alternate_mutation_triallelic_spectrum(6, verbose=False)
print(result["spectrum"][1])
```

- **示例输出：**

```text
[None, 0.5, 0.125, 0.05, 0.025, None, None]
```

## 4) `dnakit.population.triallelic_spectrum_transform()`

- **作用：** 对三等位基因频谱做屏蔽、展开、折叠或取对数。
- **处理方法：** 构建 `TriSpectrum` 并调用 `operation` 指定的方法：`mask_infeasible`、`unfold`、`fold_major`（按主、次衍生等位基因）、`fold_ancestral`（祖先状态未知，记录两个次要计数）或 `log`。一个函数覆盖五条工作清单。
- **API：** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_spectrum_transform(spectrum, operation, *, folded_major=False, folded_ancestral=False)`.
- **输入：** `spectrum`（方阵嵌套列表，被屏蔽项为 None）、`operation`，以及当前折叠状态 `folded_major`、`folded_ancestral`。
- **输出：** 变换后的 `spectrum`、`sample_size`、`folded_major`、`folded_ancestral`、`transform`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import (
    triallelic_neutral_equilibrium,
    triallelic_spectrum_transform,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
folded = triallelic_spectrum_transform(spectrum, "fold_major", verbose=False)
print(folded["folded_major"], folded["folded_ancestral"])
print([None if v is None else round(v, 3) for v in folded["spectrum"][2]])
back = triallelic_spectrum_transform(folded["spectrum"], "unfold", folded_major=True, verbose=False)
print(back["folded_major"])
```

- **示例输出：**

```text
True False
[None, 0.851, 0.247, None, None, None, None]
False
```

## 5) `dnakit.population.triallelic_segregating_sites()`

- **作用：** 三等位基因分离位点数。
- **处理方法：** 对频谱中未被屏蔽的项求和（`TriSpectrum.S`）。`apply_infeasible_mask=True` 时先屏蔽不可能出现的项，与 dadi 构建频谱时一致。
- **API：** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_segregating_sites(spectrum, *, folded_major=False, folded_ancestral=False, apply_infeasible_mask=False)`.
- **输入：** `spectrum`（方阵嵌套列表，被屏蔽项为 None）、`folded_major`、`folded_ancestral`、`apply_infeasible_mask`。
- **输出：** `segregating_sites` 和 `sample_size`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import (
    triallelic_neutral_equilibrium,
    triallelic_segregating_sites,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
result = triallelic_segregating_sites(spectrum, verbose=False)
print(result["sample_size"], round(result["segregating_sites"], 4))
```

- **示例输出：**

```text
6 3.1528
```
