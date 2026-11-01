# dadi：群体历史推断与频谱绘图

`dnakit.population` 中的这一组函数是对 [dadi](https://github.com/RyanGutenkunst/dadi)（用扩散近似做群体历史推断和适合度效应分布分析）的适配器，其中一个三维图使用 [MayaVi](https://github.com/enthought/mayavi)。demes 图通过 [demes](https://github.com/popsim-consortium/demes-python) 读取，图形用 [matplotlib](https://matplotlib.org) 绘制。这些函数把单群体、双群体、三群体位点频率谱，三等位基因频谱、双位点频谱和二元 DFE 画成 PNG 文件，计算 DFE 混合模型的期望频谱，把 demes 图的参数拟合到频谱，并计算测序深度分布。它们不读取 VCF，也不从序列数据构建频谱；频谱以嵌套列表传入，被屏蔽的项用 `None` 表示。

dadi 使用 0BSD 许可（来自已安装包的元数据）。MayaVi 的许可证在这里没有核查。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。示例从 `dnakit.population` 导入，并以 `verbose=False` 调用每个函数。绘图函数写出新的 PNG 文件，从不覆盖已有文件；示例写入临时目录，只打印可核对的信息（大小、字段、文件是否以 PNG 签名开头）。

在本环境安装的 matplotlib 3.11.2 和 dadi 2.4.4 下，dadi 的四个绘图或拟合函数会失败（`plot_trispectrum`、`plot_trispectrum_comparison`、`plot_two_locus_spectrum`、`fit_demes_graph_parameters`），`dfe_mixture_spectrum` 的 `mixture_point_pos` 变体也会失败。它们在本页都有可用的对应函数（`plot_triallelic_spectrum`、`plot_triallelic_spectrum_comparison`、`plot_two_locus_spectrum_3d`、`optimize_demes_graph_parameters` 和 `dfe_mixture_point_pos_spectrum`），四个失败函数的示例展示真实的错误。

## 示例数据

两个 demes 函数读取一个 demes 图和一个 dadi 推断选项文件。其余示例直接传入输入。

**示例文件：** `graph.yaml`

```text
description: split with two demes
time_units: generations
demes:
  - name: anc
    epochs:
      - start_size: 1000
        end_time: 500
  - name: A
    ancestors: [anc]
    epochs:
      - start_size: 800
        end_time: 0
  - name: B
    ancestors: [anc]
    epochs:
      - start_size: 600
        end_time: 0
```

**示例文件：** `options.yaml`

```text
parameters:
  - name: nuA
    values:
      - demes:
          A:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
  - name: nuB
    values:
      - demes:
          B:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
```

## 1) `dnakit.population.coverage_depth_distribution()`

- **作用：** 计算每个群体的测序深度分布。
- **处理方法：** 构造 dadi 的 LowPass 模块需要的逐位点覆盖度字典（每个位点保存每个群体的深度数组），再调用 `dadi.LowPass.LowPass.compute_cov_dist`。结果给出每个群体从 0 到最大深度的各个深度及其归一化频率。
- **API：** `coverage_depth_distribution(site_depths: 'Mapping[str, Sequence[Sequence[int]]]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。
- **输入：** `site_depths` 把每个群体 ID（唯一、非空字符串）映射到按位点排列的列表；每个位点是非空的非负整数深度列表，每个个体一个值。所有群体的位点数必须相同且不为 0（最多 1,000,000 个）。
- **输出：** `status`、`operation`、`sites`、`distribution`（每个群体的 `depth` 和 `frequency` 两个列表）和 `provenance`。
- **示例代码：**

```python
from dnakit.population import coverage_depth_distribution

site_depths = {"A": [[1, 2], [2, 2], [3, 1]], "B": [[1, 1], [1, 1], [2, 0]]}
result = coverage_depth_distribution(site_depths, verbose=False)
print(result["status"], result["sites"])
for population, row in result["distribution"].items():
    print(population, row["depth"], [round(f, 3) for f in row["frequency"]])
```

- **示例输出：**

```text
completed 3
A [0, 1, 2, 3] [0.0, 0.333, 0.5, 0.167]
B [0, 1, 2] [0.167, 0.667, 0.167]
```

## 2) `dnakit.population.dfe_mixture_spectrum()`

- **作用：** 计算 DFE（适合度效应分布）混合模型下的期望双群体频谱。
- **处理方法：** 为 `split_mig_sel_single_gamma` 构建 dadi 的 `Cache1D`，为 `split_mig_sel` 构建 `Cache2D`（gamma 范围 1e-4 到 `gamma_max`，`gamma_points` 个对数等距值，加上该变体需要的点质量 gamma），然后用对数正态与二元对数正态调用 `dadi.DFE.Cache2D_mod` 中对应的混合函数，或用 gamma 与二元独立 gamma 调用 `dadi.DFE.Vourlaki2022.Vourlaki_mixture`。频谱中的非有限值以 `None` 返回。
- **API：** `dfe_mixture_spectrum(variant: 'str', demographic_params: 'Sequence[float]', sample_sizes: 'Sequence[int]', dfe_params: 'Sequence[float]', theta: 'float' = 1.0, *, grid_points: 'Sequence[int]' = (20, 25, 30), gamma_points: 'int' = 20, gamma_max: 'float' = 2000.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。
- **输入：** `variant` 为 `mixture`（4 个 `dfe_params`）、`mixture_symmetric_point_pos`（6 个）、`mixture_point_pos`（8 个）或 `Vourlaki_mixture`（6 个）。`demographic_params` 为 nu1、nu2、T、m。`sample_sizes` 是 2 到 40 的两个整数，`grid_points` 是 8 到 200 的三个整数，`gamma_points` 是 5 到 60 的整数，`theta` 和 `gamma_max` 为正数。`demographic_params` 和 `dfe_params` 中每个数都不能小于 0，因此这里会拒绝负的对数正态 mu（`dfe_mixture_point_pos_spectrum` 接受负值）。变体 `mixture_point_pos` 在本函数中会抛出 `ConfigurationError`，因为 dadi 在需要数值的位置传入了 `None`；该变体请使用 `dfe_mixture_point_pos_spectrum`。`mixture_symmetric_point_pos` 在一次使用这些小参数的试验中返回了全为 `None` 的频谱，请检查结果中是否有 `None`。
- **输出：** `status`、`operation`、`variant`、`demographic_params`、`sample_sizes`、`dfe_params`、`theta`、`grid_points`、`gamma_points`、`gamma_max`、`spectrum`（嵌套列表，被屏蔽的项为 `None`）和 `provenance`。网格过粗时频谱可能数值不稳定，请使用接近默认值的网格。
- **示例代码：**

```python
import numpy as np

from dnakit.population import dfe_mixture_spectrum

result = dfe_mixture_spectrum(
    "mixture",
    [2.0, 2.0, 0.2, 0.1],
    [4, 4],
    [0.5, 1.0, 0.5, 0.5],
    grid_points=[20, 25, 30],
    gamma_points=10,
    gamma_max=100.0,
    verbose=False,
)
spectrum = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]])
print(result["variant"], result["sample_sizes"], spectrum.shape)
print(int(np.isnan(spectrum).sum()), round(float(np.nansum(spectrum)), 2))
```

- **示例输出：**

```text
mixture [4, 4] (5, 5)
2 3.05
```

## 3) `dnakit.population.fit_demes_graph_parameters()`

- **作用：** 把 demes 图中的自由参数拟合到观测频谱。
- **处理方法：** 调用 `dadi.Demes.Inference.optimize`，传入图文件、推断选项文件、频谱和积分网格，使用多项式似然（不设 `uL`）、所选优化器，并在 `log` 为真时使用对数参数化。
- **API：** `fit_demes_graph_parameters(graph: 'str | Path', inference_options: 'str | Path', sfs: 'Sequence[Any]', pop_ids: 'Sequence[str]', grid_points: 'Sequence[int]', *, method: 'str' = 'fmin', maxiter: 'int' = 1000, log: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。本环境安装的是 ruamel.yaml 0.19，dadi 的 `optimize` 用 `ruamel.yaml.load` 读取选项文件，而该版本已删除此函数，所以调用会抛出 `ConfigurationError`。需要仍提供 `load` 的 ruamel.yaml 版本；下面的 `optimize_demes_graph_parameters` 做同样的拟合，并且可在 ruamel.yaml 0.19 下运行。
- **输入：** `graph` 和 `inference_options` 是已存在的文件（demes YAML 图，以及带 `parameters` 列表的 dadi 选项 YAML）。`sfs` 是 1 维、2 维或 3 维的非负计数嵌套列表（被屏蔽的项为 `None`）；`pop_ids` 给各轴命名，必须与图中的 deme 对应；`grid_points` 是 8 到 200 的三个整数；`method` 为 `fmin`、`powell` 或 `lbfgsb`；`maxiter` 是 1 到 10000 的整数；`log` 为布尔值。
- **输出：** `status`、`operation`、`graph`、`inference_options`、`pop_ids`、`sample_sizes`、`grid_points`、`method`、`maxiter`、`log`、`parameter_names`、`params`（拟合值）、`fopt`（目标函数值，非有限时为 `None`）和 `provenance`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.population import fit_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
try:
    fit_demes_graph_parameters(
        "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12], maxiter=3, verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi Demes inference failed
```

示例展示了在本环境中的真实行为：拟合失败，错误以 `ConfigurationError` 报告。

- **安装依赖后的用法：**

```python
from dnakit.population import fit_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
result = fit_demes_graph_parameters(
    "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12], method="powell", maxiter=200
)
print(result["parameter_names"], result["params"], result["fopt"])
```

## 4) `dnakit.population.plot_dfe_bivariate()`

- **作用：** 把负选择系数的二元适合度效应分布画成热图。
- **处理方法：** 生成两个坐标轴相同的负 gamma 对数网格，用 `dadi.DFE.PDFs.biv_lognormal` 或 `biv_ind_gamma` 调用 `dadi.DFE.Plotting.plot_biv_dfe`，并把图保存为 PNG。
- **API：** `plot_dfe_bivariate(distribution: 'str', params: 'Sequence[float]', image_path: 'str | Path', *, gamma_min: 'float' = 0.001, gamma_max: 'float' = 100.0, grid_size: 'int' = 50, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `distribution` 为 `biv_lognormal`（`params` = mu、sigma、rho，mu 和 sigma 为对数尺度）或 `biv_ind_gamma`（`params` = alpha、beta 以及被忽略的第三个值）。`gamma_min` 和 `gamma_max` 为正数且 `gamma_min` 较小；`grid_size` 是 5 到 200 的整数；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`distribution`、`params` 和 `gamma_grid`（两个坐标轴使用的负 gamma 值）、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_dfe_bivariate

out = Path(tempfile.mkdtemp()) / "dfe.png"
result = plot_dfe_bivariate("biv_lognormal", [1.0, 1.0, 0.5], out, grid_size=20, verbose=False)
grid = result["gamma_grid"]
print(result["status"], result["distribution"], result["params"])
print(len(grid), round(grid[0], 6), round(grid[-1], 3))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed biv_lognormal [1.0, 1.0, 0.5]
20 -0.001 -100.0
True
```

## 5) `dnakit.population.plot_dfe_bivariate_point_pos()`

- **作用：** 把带正选择点质量的二元 DFE 画成热图。
- **处理方法：** 生成同样的负 gamma 网格，调用 `dadi.DFE.Plotting.plot_biv_point_pos_dfe`；它把 `params` 的最后四个值当作两个群体的正选择比例与正 gamma，前面的值作为二元分布的参数。
- **API：** `plot_dfe_bivariate_point_pos(distribution: 'str', params: 'Sequence[float]', image_path: 'str | Path', *, rho: 'float' = 0.0, gamma_min: 'float' = 0.001, gamma_max: 'float' = 100.0, grid_size: 'int' = 50, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `distribution` 为 `biv_lognormal` 或 `biv_ind_gamma`。`params` 至少 5 个数（分布参数之后接四个点质量值）。`rho`（不小于 -1）连接负选择和正选择部分；对 `biv_lognormal` 而言它与 `params` 内的 rho 是分开的。网格参数与 `plot_dfe_bivariate` 相同。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`distribution`、`params`、`rho` 和 `gamma_grid`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_dfe_bivariate_point_pos

out = Path(tempfile.mkdtemp()) / "dfe_pos.png"
params = [1.0, 1.0, 0.5, 0.1, 5.0, 0.1, 5.0]
result = plot_dfe_bivariate_point_pos(
    "biv_lognormal", params, out, rho=0.5, grid_size=20, colorbar=False, verbose=False
)
print(result["operation"], result["distribution"], result["rho"], len(result["params"]))
print(len(result["gamma_grid"]), round(result["gamma_grid"][-1], 3))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
plot_dfe_bivariate_point_pos biv_lognormal 0.5 7
20 -100.0
True
```

## 6) `dnakit.population.plot_sfs_1d()`

- **作用：** 绘制单群体位点频率谱。
- **处理方法：** 把列表包装成 `dadi.Spectrum`，调用 `dadi.Plotting.plot_1d_fs`，并把图保存为 PNG。
- **API：** `plot_sfs_1d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是非空的一维非负计数列表；`None` 表示被屏蔽的项。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`sample_sizes` 和 `spectrum`（与 dadi 持有的一致，被屏蔽的项为 `None`）、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_1d

out = Path(tempfile.mkdtemp()) / "sfs1d.png"
result = plot_sfs_1d([None, 10, 5, 3, 2, None], out, verbose=False)
print(result["status"], result["sample_sizes"], result["spectrum"])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed [5] [None, 10.0, 5.0, 3.0, 2.0, None]
True
```

## 7) `dnakit.population.plot_sfs_1d_comparison()`

- **作用：** 绘制单群体的模型频谱与数据频谱及其残差。
- **处理方法：** 按所选残差类型调用 `dadi.Plotting.plot_1d_comp_multinom` 或 `plot_1d_comp_Poisson`。
- **API：** `plot_sfs_1d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', residual: 'str' = 'Anscombe', plot_masked: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `model` 和 `data` 是一维非负计数列表。`likelihood` 为 `multinom` 或 `poisson`；`residual` 为 `Anscombe` 或 `linear`；`plot_masked` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`likelihood`、`residual`、`model` 和 `data`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_1d_comparison

model = [None, 9, 5, 3, 2, None]
data = [None, 10, 5, 3, 2, None]
folder = Path(tempfile.mkdtemp())
first = plot_sfs_1d_comparison(model, data, folder / "default.png", verbose=False)
second = plot_sfs_1d_comparison(
    model, data, folder / "poisson.png", likelihood="poisson", residual="linear", verbose=False
)
print(first["likelihood"], first["residual"], first["model"])
print(second["likelihood"], second["residual"], second["data"])
print((folder / "default.png").read_bytes()[:4] == b"\x89PNG", (folder / "poisson.png").exists())
```

- **示例输出：**

```text
multinom Anscombe [None, 9.0, 5.0, 3.0, 2.0, None]
poisson linear [None, 10.0, 5.0, 3.0, 2.0, None]
True True
```

## 8) `dnakit.population.plot_sfs_2d()`

- **作用：** 把双群体联合位点频率谱画成热图。
- **处理方法：** 把列表包装成 `dadi.Spectrum`，调用 `dadi.Plotting.plot_single_2d_sfs`。
- **API：** `plot_sfs_2d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是二维非负计数列表（被屏蔽的项为 `None`）。`vmin` 和 `vmax` 是可选的颜色范围且 `vmin` 较小；`pop_ids` 为两个互不相同的名称；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`sample_sizes` 和 `spectrum`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d

out = Path(tempfile.mkdtemp()) / "sfs2d.png"
spectrum = [[None, 3, 1, 0.5], [2, 1, 0.5, 0.2], [1, 0.5, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
result = plot_sfs_2d(spectrum, out, vmin=0.1, vmax=3.0, pop_ids=["A", "B"], verbose=False)
print(result["status"], result["sample_sizes"])
print(result["spectrum"][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed [3, 3]
[None, 3.0, 1.0, 0.5]
True
```

## 9) `dnakit.population.plot_sfs_2d_comparison()`

- **作用：** 并排绘制双群体频谱的模型、数据和残差。
- **处理方法：** 调用 `dadi.Plotting.plot_2d_comp_multinom` 或 `plot_2d_comp_Poisson`。
- **API：** `plot_sfs_2d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, residual: 'str' = 'Anscombe', adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `model` 和 `data` 是二维非负计数列表。`likelihood` 为 `multinom` 或 `poisson`；`residual` 为 `Anscombe` 或 `linear`；`vmin`/`vmax`（`vmin` 较小）、`resid_range`（正数）、`pop_ids` 和 `adjust`（布尔值）为可选。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`likelihood`、`residual`、`model` 和 `data`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_comparison

out = Path(tempfile.mkdtemp()) / "comparison2d.png"
model = [[None, 3, 1, 0.5], [2, 1, 0.5, 0.2], [1, 0.5, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
data = [[None, 2.5, 1, 0.5], [2, 1.2, 0.4, 0.2], [1, 0.6, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
result = plot_sfs_2d_comparison(model, data, out, pop_ids=["A", "B"], verbose=False)
print(result["likelihood"], result["residual"])
print(result["model"][1], result["data"][1])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
multinom Anscombe
[2.0, 1.0, 0.5, 0.2] [2.0, 1.2, 0.4, 0.2]
True
```

## 10) `dnakit.population.plot_sfs_2d_meta_residual()`

- **作用：** 绘制同义和非同义双群体残差频谱的元残差比较图。
- **处理方法：** 把两个残差分别包装成 `dadi.Spectrum`（允许负值），调用 `dadi.Plotting.plot_2d_meta_resid`。
- **API：** `plot_sfs_2d_meta_residual(synonymous_residual: 'Sequence[Any]', nonsynonymous_residual: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `synonymous_residual` 和 `nonsynonymous_residual` 是有限数的二维列表（被屏蔽的项为 `None`）。`resid_range` 若给出必须为正；`pop_ids` 为两个互不相同的名称；`adjust` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`synonymous_residual`、`nonsynonymous_residual` 和 `resid_range`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_meta_residual

out = Path(tempfile.mkdtemp()) / "meta.png"
synonymous = [[None, 0.5, -0.2, 0.1], [0.1, -0.3, 0.2, 0.0], [-0.1, 0.2, 0.0, -0.1], [0.2, 0.0, -0.1, None]]
nonsynonymous = [[None if v is None else -v for v in row] for row in synonymous]
result = plot_sfs_2d_meta_residual(
    synonymous, nonsynonymous, out, resid_range=1.0, pop_ids=["A", "B"], verbose=False
)
print(result["status"], result["resid_range"])
print(result["synonymous_residual"][0], result["nonsynonymous_residual"][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed 1.0
[None, 0.5, -0.2, 0.1] [None, -0.5, 0.2, -0.1]
True
```

## 11) `dnakit.population.plot_sfs_2d_residual()`

- **作用：** 把双群体残差频谱（带符号的值）画成热图。
- **处理方法：** 把残差包装成 `dadi.Spectrum`（允许负值），调用 `dadi.Plotting.plot_2d_resid`。
- **API：** `plot_sfs_2d_residual(residual: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `residual` 是有限数的二维列表（被屏蔽的项为 `None`）。`resid_range`（正数）把颜色标尺限制在正负该值之间；`pop_ids` 为两个互不相同的名称；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`residual` 和 `resid_range`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_residual

out = Path(tempfile.mkdtemp()) / "residual2d.png"
residual = [[None, 0.5, -0.2, 0.1], [0.1, -0.3, 0.2, 0.0], [-0.1, 0.2, 0.0, -0.1], [0.2, 0.0, -0.1, None]]
result = plot_sfs_2d_residual(residual, out, resid_range=1.0, pop_ids=["A", "B"], verbose=False)
print(result["status"], result["resid_range"])
print(result["residual"][1])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed 1.0
[0.1, -0.3, 0.2, 0.0]
True
```

## 12) `dnakit.population.plot_sfs_3d()`

- **作用：** 把三群体频谱画成按切片排列的对数热图。
- **处理方法：** 把列表包装成 `dadi.Spectrum`，调用 `dadi.Plotting.plot_3d_spectrum`。
- **API：** `plot_sfs_3d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是三维非负计数列表（被屏蔽的项为 `None`）。`vmin` 和 `vmax` 可选（`vmin` 较小）；`pop_ids` 为三个互不相同的名称。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`sample_sizes` 和 `spectrum`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "sfs3d.png"
result = plot_sfs_3d(cube(6.0), out, vmin=0.1, vmax=6.0, pop_ids=["A", "B", "C"], verbose=False)
print(result["status"], result["sample_sizes"])
print(result["spectrum"][0][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed [2, 2, 2]
[None, 3.0, 2.0]
True
```

## 13) `dnakit.population.plot_sfs_3d_comparison()`

- **作用：** 以切片面板绘制三群体频谱的模型、数据和残差。
- **处理方法：** 调用 `dadi.Plotting.plot_3d_comp_multinom` 或 `plot_3d_comp_Poisson`。
- **API：** `plot_sfs_3d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, residual: 'str' = 'Anscombe', adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `model` 和 `data` 是三维非负计数列表。其余参数与 `plot_sfs_2d_comparison` 相同。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`likelihood`、`residual`、`model` 和 `data`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d_comparison

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "comparison3d.png"
result = plot_sfs_3d_comparison(cube(6.0), cube(6.6), out, pop_ids=["A", "B", "C"], verbose=False)
print(result["likelihood"], result["residual"])
print(result["model"][0][0], [None if v is None else round(v, 2) for v in result["data"][0][0]])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
multinom Anscombe
[None, 3.0, 2.0] [None, 3.3, 2.2]
True
```

## 14) `dnakit.population.plot_sfs_3d_pairwise()`

- **作用：** 绘制三群体频谱的三个两两双群体边缘频谱。
- **处理方法：** 调用 `dadi.Plotting.plot_3d_pairwise`，它对每一对群体把频谱在第三个群体上做边缘化。
- **API：** `plot_sfs_3d_pairwise(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是三维非负计数列表（被屏蔽的项为 `None`）。`vmin`/`vmax`（`vmin` 较小）、`pop_ids`（三个名称）和 `adjust`（布尔值）为可选。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`sample_sizes` 和 `spectrum`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d_pairwise

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "pairwise3d.png"
result = plot_sfs_3d_pairwise(cube(6.0), out, pop_ids=["A", "B", "C"], verbose=False)
print(result["status"], result["sample_sizes"])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed [2, 2, 2]
True
```

## 15) `dnakit.population.plot_trispectrum()`

- **作用：** 用 dadi 的 Triallele 绘图模块把三等位基因频谱画成热图。
- **处理方法：** 把列表包装成 `dadi.Triallele` 的 `TriSpectrum`，调用 `dadi.Triallele.plotting.plot_single_trispectrum`。该上游函数把 `LogNorm` 与 `vmin`、`vmax` 同时传给 `pcolor`，而 matplotlib 3.11（本环境所装版本）会拒绝这种调用，所以在本环境中每次调用都会抛出 `ConfigurationError`。请改用 `plot_triallelic_spectrum`（本页第 23 节），它绘制同样的图，但不会触发这个调用。
- **API：** `plot_trispectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。已知在测试环境的 matplotlib 3.11.2 下会失败。
- **输入：** `spectrum` 是二维非负计数列表（被屏蔽的项为 `None`）；`folded`、`colorbar` 为布尔值；`vmin`/`vmax` 可选。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** 成功时返回 `status`、`operation`、`image_path`、`image_format`、`folded`、`spectrum` 和 `provenance`。在本环境中调用失败并抛出 `ConfigurationError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_trispectrum

n = 7
spectrum = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri.png"
try:
    plot_trispectrum(spectrum, out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_trispectrum
False
```

示例展示了在本环境中的真实行为：抛出错误，且不会写出文件。

## 16) `dnakit.population.plot_trispectrum_comparison()`

- **作用：** 用 dadi 的 Triallele 绘图模块并排绘制两个三等位基因频谱及其残差。
- **处理方法：** 把两个列表分别包装成 `TriSpectrum`，调用 `dadi.Triallele.plotting.plot_trispectrum_comp`。它的两个频谱面板使用同一个 `plot_single_trispectrum` 调用，所以在 matplotlib 3.11 下与 `plot_trispectrum` 一样失败。请改用 `plot_triallelic_spectrum_comparison`（第 24 节）。
- **API：** `plot_trispectrum_comparison(first: 'Sequence[Any]', second: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, title1: 'str' = 'sfs1', title2: 'str' = 'sfs2', title3: 'str' = 'Residual', colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。已知在测试环境的 matplotlib 3.11.2 下会失败。
- **输入：** `first` 和 `second` 是二维非负计数列表；`folded`、`colorbar` 为布尔值；`vmin`、`vmax`、`resid_range` 可选；`title1`、`title2`、`title3` 为字符串。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** 成功时返回 `status`、`operation`、`image_path`、`image_format`、`folded`、`first`、`second` 和 `provenance`。在本环境中调用失败并抛出 `ConfigurationError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_trispectrum_comparison

n = 7
spectrum = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri_comp.png"
try:
    plot_trispectrum_comparison(spectrum, spectrum, out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_trispectrum_comparison
False
```

示例展示了在本环境中的真实行为：抛出错误，且不会写出文件。

## 17) `dnakit.population.plot_two_locus_comparison()`

- **作用：** 以三维面板绘制双位点模型频谱与数据频谱及残差。
- **处理方法：** 把两个列表分别包装成 dadi 的 `TLSpectrum`，调用 `dadi.TwoLocus.plotting.plot_3d_comp`。
- **API：** `plot_two_locus_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float' = 1.0, max_sum: 'int | None' = None, max_index: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `model` 和 `data` 是三维非负计数列表。`resid_range` 为正数；`max_sum` 和 `max_index` 是用于切片的可选整数。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`resid_range`、`model` 和 `data`（保持 `TLSpectrum` 自带的屏蔽，因此保留的项比传入的少）、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_two_locus_comparison

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

def unmasked(cube3):
    return sum(v is not None for plane in cube3 for row in plane for v in row)


out = Path(tempfile.mkdtemp()) / "two_locus_comp.png"
result = plot_two_locus_comparison(cube(6.0), cube(6.6), out, resid_range=1.0, verbose=False)
print(result["status"], result["resid_range"])
print(unmasked(result["model"]), unmasked(result["data"]))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
completed 1.0
2 2
True
```

## 18) `dnakit.population.plot_two_locus_spectrum()`

- **作用：** 用 dadi 的 TwoLocus 绘图模块把双位点（三维）频谱画成三维散点图。
- **处理方法：** 把列表包装成 `TLSpectrum`，调用 `dadi.TwoLocus.plotting.plot_3d`。该上游函数总是以 `ax.colorbar()` 结尾，而三维坐标轴没有这个方法，所以每次调用都会抛出 `ConfigurationError`。请改用 `plot_two_locus_spectrum_3d`（第 25 节），它正确绘制体素并添加颜色条。
- **API：** `plot_two_locus_spectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, max_sum: 'int | None' = None, max_index: 'int | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。在 dadi 2.4.4 下对任何输入都会失败。
- **输入：** `spectrum` 是三维非负计数列表；`vmin`、`vmax`、`max_sum`、`max_index` 可选；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** 成功时返回 `status`、`operation`、`image_path`、`image_format`、`spectrum` 和 `provenance`。在本环境中调用失败并抛出 `ConfigurationError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_two_locus_spectrum

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "two_locus.png"
try:
    plot_two_locus_spectrum(cube(6.0), out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_two_locus_spectrum
False
```

示例展示了在本环境中的真实行为：抛出错误，且不会写出文件。

## 19) `dnakit.population.plot_3d_spectrum_mayavi()`

- **作用：** 把三群体频谱画成 MayaVi 三维柱状图并保存为 PNG。
- **处理方法：** 先检查输入，然后对每个未被屏蔽且不小于 `vmin` 的项用 `mayavi.mlab.barchart` 画一根柱子，颜色由计数的对数在 `vmin` 与 `vmax` 之间缩放得到，并用 `mlab.savefig` 保存。它遵循 dadi 的 `plot_3d_spectrum_mayavi`。
- **API：** `plot_3d_spectrum_mayavi(spectrum: 'Sequence[Any]', output_path: 'str', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 mayavi（额外依赖 `dnakit[bio-mayavi]`，同时会安装 dadi）以及显示器或离屏 VTK 构建。没有 mayavi 时抛出带安装提示的 `BackendUnavailableError`；测试环境没有安装 MayaVi。
- **输入：** `spectrum` 是三维非负计数列表，且至少有一个未被屏蔽的正值。`vmin` 和 `vmax` 为正数且 `vmin` 较小，默认取最小和最大的正计数。`pop_ids` 为三个互不相同的名称（默认 `pop0`、`pop1`、`pop2`）。`output_path` 必须以 `.png` 结尾且文件不能已存在；其父目录会被创建。
- **输出：** 返回包含 `path` 和 `data` 的字典；`data` 含 `pop_ids`、`sample_sizes`、`vmin`、`vmax`、`plotted_entries`、`spectrum` 和 `provenance`（本函数不返回 `status` 和 `operation`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.population import plot_3d_spectrum_mayavi

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "bars.png"
try:
    plot_3d_spectrum_mayavi(cube(6.0), out, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] mayavi.mlab could not be imported
```

示例展示了没有 MayaVi 时的真实行为。如果已安装 MayaVi，调用会绘制图形，且不会打印任何内容。

- **安装依赖后的用法：**

```python
from dnakit.population import plot_3d_spectrum_mayavi

corners = {(0, 0, 0), (2, 2, 2)}
spectrum = [
    [[None if (i, j, k) in corners else 6.0 / (1 + i + j + k) for k in range(3)] for j in range(3)]
    for i in range(3)
]
result = plot_3d_spectrum_mayavi(spectrum, "bars.png", vmin=0.5, pop_ids=["A", "B", "C"])
print(result["path"], result["data"]["plotted_entries"], result["data"]["vmin"], result["data"]["vmax"])
```

## 20) `dnakit.population.vourlaki_visual_validation_png()`

- **作用：** 把 dadi 的 Vourlaki（2022）DFE 验证图重绘为 PNG，并返回图中的频谱。
- **处理方法：** 为 dadi 的 visual-validation 模型（nu1 = nu2 = 2，T = 0.2，m = 0，gamma 范围 1e-4 到 200，gamma 为 10 处有一个点质量）构建 `Cache1D` 和 `Cache2D`，对四种情形（A 到 D 面板：ppos_wild、pchange、pchange_pos 分别为 0/0/0、1/0/0、0/1/0、0/1/0.9；alpha = 0.2，beta = 10）计算 `Vourlaki_mixture`，用 `plot_single_2d_sfs` 绘制每个双群体频谱，并在对数坐标面板中绘制野生与驯化群体的边缘频谱，最后保存为 10 乘 8 英寸、150 dpi 的一张图。
- **API：** `vourlaki_visual_validation_png(output_path: 'str', *, sample_sizes: 'Sequence[int]' = (10, 10), grid_points: 'Sequence[int]' = (30, 35, 40), gamma_points: 'int' = 30, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `output_path` 必须以 `.png` 结尾且文件不能已存在；其父目录会被创建。`sample_sizes` 是 2 到 40 的两个整数，`grid_points` 是 8 到 200 的三个整数，`gamma_points` 是 5 到 60 的整数。默认值与上游一致（(10, 10)、(30, 35, 40) 和 30），耗时比下面示例中的小参数更长。
- **输出：** 返回包含 `path` 和 `data` 的字典；`data` 含 `sample_sizes`、`grid_points`、`gamma_points`、`demographic_params`、`alpha`、`beta`、`gamma_pos`、`panels`（每个面板的 `panel`、`ppos_wild`、`pchange`、`pchange_pos`、`spectrum`、`wild_marginal`、`domesticate_marginal`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.population import vourlaki_visual_validation_png

out = Path(tempfile.mkdtemp()) / "figures" / "validation.png"
result = vourlaki_visual_validation_png(
    out, sample_sizes=[4, 4], grid_points=[20, 25, 30], gamma_points=10, verbose=False
)
data = result["data"]
print(sorted(result), out.read_bytes()[:4] == b"\x89PNG")
print(data["sample_sizes"], data["demographic_params"], data["alpha"], data["beta"], data["gamma_pos"])
for panel in data["panels"]:
    spectrum = np.array([[np.nan if v is None else v for v in row] for row in panel["spectrum"]])
    print(panel["panel"], panel["ppos_wild"], panel["pchange"], panel["pchange_pos"],
          spectrum.shape, round(float(np.nansum(spectrum)), 2), len(panel["wild_marginal"]))
```

- **示例输出：**

```text
['data', 'path'] True
[4, 4] [2, 2, 0.2, 0] 0.2 10.0 10.0
A 0.0 0.0 0.0 (5, 5) 3.43 5
B 1.0 0.0 0.0 (5, 5) 6.41 5
C 0.0 1.0 0.0 (5, 5) 4.93 5
D 0.0 1.0 0.9 (5, 5) 5.55 5
```

## 21) `dnakit.population.dfe_mixture_point_pos_spectrum()`

- **作用：** 计算 dadi 点质量 DFE 混合模型（`mixture_point_pos`）的期望双群体频谱。
- **处理方法：** 为 `split_mig_sel_single_gamma` 和 `split_mig_sel` 构建 `Cache1D` 与 `Cache2D`（两个正 gamma 作为额外缓存值），用 `integrate_point_pos` 计算一维部分，用二维缓存的 `integrate_point_pos` 计算二维部分，再按 (1 - 权重) 乘第一部分加权重乘第二部分合并。上游 `mixture_point_pos` 在该处传入 `None` 而失败，这里传入的是 `rho`。
- **API：** `dfe_mixture_point_pos_spectrum(demographic_params: 'Sequence[float]', sample_sizes: 'Sequence[int]', dfe_params: 'Sequence[float]', theta: 'float' = 1.0, *, grid_points: 'Sequence[int]' = (20, 25, 30), gamma_points: 'int' = 20, gamma_max: 'float' = 2000.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。
- **输入：** `demographic_params` 为 nu1、nu2、T、m（各不小于 0）。`dfe_params` 恰好是 8 个有限数：对数正态的 mu 和 sigma、rho、群体 1 的正选择比例和 gamma、群体 2 的正选择比例和 gamma，以及二维部分的权重（在 0 到 1 之间）；mu 可以为负。两个正 gamma 必须不同，否则缓存形状不同，dadi 会报错。其余参数与 `dfe_mixture_spectrum` 相同。
- **输出：** `status`、`operation`、`variant`（`mixture_point_pos`）、`demographic_params`、`sample_sizes`、`dfe_params`、`theta`、`grid_points`、`gamma_points`、`gamma_max`、`spectrum` 和 `provenance`。
- **示例代码：**

```python
import numpy as np

from dnakit.population import dfe_mixture_point_pos_spectrum

dfe_params = [-1.0, 0.6, 0.4, 0.2, 5.0, 0.3, 3.0, 0.6]
result = dfe_mixture_point_pos_spectrum(
    [2.0, 2.0, 0.2, 0.1],
    [4, 4],
    dfe_params,
    grid_points=[20, 25, 30],
    gamma_points=10,
    gamma_max=100.0,
    verbose=False,
)
spectrum = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]])
print(result["variant"], result["sample_sizes"], spectrum.shape)
print(int(np.isnan(spectrum).sum()), round(float(np.nansum(spectrum)), 2))
```

- **示例输出：**

```text
mixture_point_pos [4, 4] (5, 5)
2 5.22
```

## 22) `dnakit.population.optimize_demes_graph_parameters()`

- **作用：** 把 demes 图中的自由参数拟合到双群体或三群体频谱。
- **处理方法：** 自行执行 dadi `Demes.Inference.optimize` 的各个步骤（参数设置、边界、约束、可选扰动、对数参数化、`scipy.optimize.fmin`、`fmin_powell` 或 `fmin_l_bfgs_b`），但用 `ruamel.yaml.YAML(typ="safe")` 读取选项文件，因此可在 ruamel.yaml 0.19 下运行。若给出 `output_path`，则用 `demes.dump` 写出拟合后的图。
- **API：** `optimize_demes_graph_parameters(graph_path: 'str | Path', inference_options_path: 'str | Path', sfs: 'Sequence[Any]', pop_ids: 'Sequence[str]', grid_points: 'Sequence[int]', *, maxiter: 'int' = 1000, perturb: 'float' = 0, log: 'bool' = True, method: 'str' = 'fmin', uL: 'float | None' = None, output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。
- **输入：** `graph_path` 和 `inference_options_path` 是已存在的文件。`sfs` 是 1 维到 3 维的非负计数嵌套列表，`pop_ids` 给各轴命名。`grid_points` 是 8 到 200 的三个整数；`method` 为 `fmin`、`powell` 或 `lbfgsb`；`maxiter` 是 1 到 10000 的整数；`perturb` 为不小于 0 的数；`log` 为布尔值；`uL` 是可选的正数；`output_path` 是可选的新文件。
- **输出：** `status`、`operation`、`method`、`log`、`maxiter`、`parameter_names`、`parameters`（拟合值）、`objective`、`output_path`（拟合后 YAML 的路径或 `None`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import optimize_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
fitted = Path(tempfile.mkdtemp()) / "fitted.yaml"
result = optimize_demes_graph_parameters(
    "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12],
    maxiter=3, output_path=fitted, verbose=False,
)
print(result["method"], result["log"], result["parameter_names"])
print([round(v) for v in result["parameters"]], round(result["objective"], 2))
text = fitted.read_text()
print(text.splitlines()[0], "start_size: 868" in text)
```

- **示例输出：**

```text
fmin True ['nuA', 'nuB']
[800, 869] 11.73
description: split with two demes True
```

## 23) `dnakit.population.plot_triallelic_spectrum()`

- **作用：** 把三等位基因频谱画成对数尺度热图并写入新的 PNG。
- **处理方法：** 当 `folded` 为真时先折叠频谱（与 dadi Triallele 的 `_fold` 一致），然后只用 `LogNorm` 通过 matplotlib 的 `pcolor` 绘制。小于 `vmin` 的项不绘制。它取代了在 matplotlib 3.11 下会失败的 dadi `plot_single_trispectrum`。
- **API：** `plot_triallelic_spectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, colorbar: 'bool' = False, cmap: 'str' = 'cubehelix_r', verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是非负计数的二维列表（被屏蔽的项为 `None`；折叠时必须是方阵且至少 3 乘 3，并且 vmax 必须大于 vmin）。`vmin` 默认取最小的正计数，`vmax` 默认取最大值；`cmap` 是 matplotlib 颜色表名称（默认 `cubehelix_r`）；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`folded`、`vmin`、`vmax`、`cmap`、`spectrum`（实际绘制的项，其余为 `None`）、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_triallelic_spectrum

n = 7
first = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]

def drawn(grid):
    return sum(v is not None for row in grid for v in row)


folder = Path(tempfile.mkdtemp())
plain = plot_triallelic_spectrum(first, folder / "tri.png", verbose=False)
folded = plot_triallelic_spectrum(first, folder / "tri_folded.png", folded=True, colorbar=True, verbose=False)
print(plain["folded"], plain["vmin"], plain["vmax"], plain["cmap"], drawn(plain["spectrum"]))
print(folded["folded"], folded["vmin"], folded["vmax"], drawn(folded["spectrum"]))
print((folder / "tri.png").read_bytes()[:4] == b"\x89PNG", (folder / "tri_folded.png").exists())
```

- **示例输出：**

```text
False 1.0 22.0 cubehelix_r 36
True 14.0 36.0 6
True True
```

## 24) `dnakit.population.plot_triallelic_spectrum_comparison()`

- **作用：** 以三个纵向排列的面板绘制两个三等位基因频谱及其残差。
- **处理方法：** 用同一个对数颜色标尺绘制两个频谱，计算残差 (first - second) / sqrt(first)，并用以零为中心的发散色表 `RdBu` 绘制。它取代了在 matplotlib 3.11 下会失败的 dadi `plot_trispectrum_comp`。`folded=True` 目前会抛出 `ConfigurationError`（"the spectrum is already folded"），因为已折叠频谱的残差被再次折叠。
- **API：** `plot_triallelic_spectrum_comparison(first: 'Sequence[Any]', second: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, colorbar: 'bool' = False, titles: 'Sequence[str]' = ('sfs1', 'sfs2', 'Residual'), cmap: 'str' = 'cubehelix_r', verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `first` 和 `second` 是形状相同的非负计数方阵（二维列表）。`resid_range` 若给出必须为正（默认取最大残差绝对值）；`titles` 是三个字符串；`vmin`、`vmax`、`cmap` 和 `colorbar` 与 `plot_triallelic_spectrum` 相同。请保持 `folded=False`。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`folded`、`vmin`、`vmax`、`resid_range`、`first`、`second`、`residual`、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_triallelic_spectrum_comparison

n = 7
first = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
second = [[0.0 if i == 0 or j == 0 else max(1.0, 28.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri_comp.png"
result = plot_triallelic_spectrum_comparison(
    first, second, out, titles=["model", "data", "Residual"], verbose=False
)
print(result["vmin"], result["vmax"], round(result["resid_range"], 4))
print([round(v, 4) for v in result["residual"][1][1:4]])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **示例输出：**

```text
1.0 22.0 0.8165
[0.4264, 0.4714, 0.5345]
True
```

## 25) `dnakit.population.plot_two_locus_spectrum_3d()`

- **作用：** 把三维双位点频谱画成彩色体素并写入新的 PNG。
- **处理方法：** 隐藏小于 `vmin` 的项，应用 `max_sum` 和 `max_index`，用压缩后的 `cubehelix_r` 色表在对数尺度上给每个体素着色，只用 `Poly3DCollection` 绘制不与其他体素相邻的面，并用 `figure.colorbar` 添加颜色条。它取代了在三维坐标轴上会失败的 dadi `TwoLocus.plotting.plot_3d`。
- **API：** `plot_two_locus_spectrum_3d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, max_sum: 'int | None' = None, max_index: 'int | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 dadi（额外依赖 `dnakit[bio-dadi]`，同时会安装 demes；dadi 本身依赖 NumPy、SciPy 和 matplotlib）。图用 matplotlib 的 Agg 后端绘制，不需要显示器。
- **输入：** `spectrum` 是立方体形的三维非负计数列表（被屏蔽的项为 `None`）。`vmin` 默认取最小的正计数，`vmax` 默认取最大值；`max_sum` 隐藏下标之和超过它的体素；`max_index` 隐藏三个下标都超过它的体素（均为不小于 0 的可选整数）；`colorbar` 为布尔值。`image_path` 必须以 `.png` 结尾，文件不能已存在，且其父目录必须存在。
- **输出：** `status`、`operation`、`image_path`（新 PNG 的绝对路径）、`image_format`（`png`）、`vmin`、`vmax`、`max_sum`、`max_index`、`drawn_entries`（绘制的体素数）、`spectrum`（不小于 `vmin` 的项）、以及 `provenance`（dadi 版本和绘图约定）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_two_locus_spectrum_3d

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

folder = Path(tempfile.mkdtemp())
full = plot_two_locus_spectrum_3d(cube(6.0), folder / "voxels.png", colorbar=True, verbose=False)
cut = plot_two_locus_spectrum_3d(cube(6.0), folder / "voxels_cut.png", max_sum=2, verbose=False)
print(full["vmin"], full["vmax"], full["drawn_entries"])
print(cut["max_sum"], cut["drawn_entries"])
print((folder / "voxels.png").read_bytes()[:4] == b"\x89PNG", (folder / "voxels_cut.png").exists())
```

- **示例输出：**

```text
1.0 3.0 25
2 9
True True
```
