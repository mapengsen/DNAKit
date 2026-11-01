# dadi 双位点频谱与频谱文件

`dnakit.population` 中封装 `dadi.TwoLocus` 对已有频谱的操作（变换、边缘频谱、平均 r^2、各格 LD、对数似然、缩放、取条件、密度抽样）以及三等位基因和双位点频谱文本文件格式的函数。双位点频谱是 n 条染色体的单倍型计数立方体 `spectrum[AB][Ab][aB]`（ab 为 n 减去前三者），用嵌套列表表示，被屏蔽项为 `None`。需要 `dnakit[bio-dadi]`（测试版本 dadi 2.4.4）；导入 `dadi.TwoLocus` 时 dadi 会创建 `~/.dadi/TwoLocus_cache`。不提供 dadi 的双位点人口历史模型、向前积分及调用它们的三个优化器（它们抛出 `dadi.tridiag` 的 `AttributeError`，补上后返回 NaN，或使用了已移除的 `fmin_l_bfgs_b` 参数），也不提供为 Python 2 编写的基因型频谱函数。

## 1) `dnakit.population.two_locus_spectrum_transform()`

- **作用：** 对双位点频谱做屏蔽、折叠、展开、投影或祖先错判处理。
- **处理方法：** 按 `operation` 调用 dadi 例程：`mask_infeasible`、`unfold`、`fold`（`TLSpectrum` 方法），`fold_ancestral`、`fold_lr`、`project`、`misidentification`（`TwoLocus.numerics`），或 `array_to_spectrum`（对原始数值施加 dadi 的不可能项屏蔽）。不修改输入。一个函数覆盖八条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_spectrum_transform(spectrum, operation, *, folded=False, sample_size=None, p_misid=None, apply_infeasible_mask=False)`.
- **输入：** 边长为 n + 1 的双位点频谱 `spectrum[AB][Ab][aB]`（被屏蔽项为 None）、`operation`，以及按需提供的 `folded`、`sample_size`（投影）、`p_misid`（祖先错判）。
- **输出：** 新的 `spectrum`、`sample_size`、`folded`、`transform`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
folded = two_locus_spectrum_transform(spectrum, "fold", verbose=False)
projected = two_locus_spectrum_transform(spectrum, "project", sample_size=4, verbose=False)
print(folded["folded"], folded["sample_size"], projected["sample_size"])
print(round(projected["spectrum"][1][1][1], 4))
```

- **示例输出：**

```text
True 6 4
3.48
```

## 2) `dnakit.population.two_locus_marginal_spectra()`

- **作用：** 由双位点频谱得到位点 A 和位点 B 的单位点频谱。
- **处理方法：** `method='marginal'` 调用 `TLSpectrum.marginalA` 和 `marginalB`；`method='single_locus'` 调用 `TwoLocus.numerics.to_single_locus`，会跳过被屏蔽的项。一个函数覆盖三条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_marginal_spectra(spectrum, *, method='marginal', folded=False)`.
- **输入：** `spectrum`、`method`、`folded`。
- **输出：** `locus_a` 和 `locus_b`（第 k 项为携带 k 个衍生等位基因的位点数）、`sample_size`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_marginal_spectra,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
result = two_locus_marginal_spectra(spectrum, verbose=False)
print([None if v is None else round(v, 2) for v in result["locus_a"]])
```

- **示例输出：**

```text
[None, 13.5, 18.2, 20.3, 19.5, 15.5, None]
```

## 3) `dnakit.population.two_locus_mean_r2()`

- **作用：** 双位点频谱的平均 r^2。
- **处理方法：** `method='numerics'` 调用 `TwoLocus.numerics.mean_r2`（遍历未屏蔽项）；`method='spectrum_method'` 调用 `TLSpectrum.mean_r2`。二者都用频谱对 r^2 = D^2 / (pA (1 - pA) pB (1 - pB)) 加权。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_mean_r2(spectrum, *, method='numerics', folded=False)`.
- **输入：** `spectrum`、`method`、`folded`。
- **输出：** `mean_r2` 和 `sample_size`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_mean_r2,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
print(round(two_locus_mean_r2(spectrum, verbose=False)["mean_r2"], 4))
print(round(two_locus_mean_r2(spectrum, method="spectrum_method", verbose=False)["mean_r2"], 4))
```

- **示例输出：**

```text
0.3539
0.3539
```

## 4) `dnakit.population.two_locus_ld_per_bin()`

- **作用：** 每个单倍型计数格的连锁不平衡 D 与 r^2。
- **处理方法：** 调用 `TwoLocus.numerics.LD_per_bin(n)`；D = pAB - pA pB，r^2 = D^2 / (pA (1 - pA) pB (1 - pB))，其中 p 为单倍型计数除以 n。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_ld_per_bin(sample_size)`.
- **输入：** `sample_size`（2 到 40）。
- **输出：** `D` 和 `r2`，边长 n + 1 的立方体（未定义处为 None）。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import two_locus_ld_per_bin

result = two_locus_ld_per_bin(4, verbose=False)
print(result["D"][2][0][0], result["r2"][2][0][0])
```

- **示例输出：**

```text
0.25 1.0
```

## 5) `dnakit.population.two_locus_log_likelihood()`

- **作用：** 双位点数据频谱的泊松或多项对数似然。
- **处理方法：** 调用 `TwoLocus.inference.ll`（泊松）、`ll_multinom`（模型按最优标量缩放），`per_bin=True` 时调用 `ll_per_bin` 与 `ll_multinom_per_bin`。一个函数覆盖四条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_log_likelihood(model, data, *, multinom=False, per_bin=False)`.
- **输入：** 同样大小的双位点频谱 `model` 与 `data`、`multinom`、`per_bin`。
- **输出：** `log_likelihood`，`per_bin` 时为 `log_likelihood_per_bin`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_log_likelihood,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
data = [[[v * 2 if v is not None else None for v in row] for row in plane] for plane in spectrum]
print(round(two_locus_log_likelihood(spectrum, data, verbose=False)["log_likelihood"], 3))
print(round(two_locus_log_likelihood(spectrum, data, multinom=True, verbose=False)["log_likelihood"], 3))
```

- **示例输出：**

```text
-122.296
-88.688
```

## 6) `dnakit.population.two_locus_log_likelihood_over_rho_bins()`

- **作用：** 跨重组率分箱的双位点总对数似然。
- **处理方法：** 调用 `TwoLocus.inference.ll_over_rho_bins`，`multinom=True` 时调用 `ll_over_rho_bins_multinom`，对各分箱对数似然求和。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.two_locus_log_likelihood_over_rho_bins(models, data, *, multinom=False)`.
- **输入：** `models` 与 `data`：每个 rho 分箱一个双位点频谱的列表，以及 `multinom`。
- **输出：** `log_likelihood`、`num_bins`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_log_likelihood_over_rho_bins,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
total = two_locus_log_likelihood_over_rho_bins([spectrum, spectrum], [spectrum, spectrum], verbose=False)
print(total["num_bins"], round(total["log_likelihood"], 3))
```

- **示例输出：**

```text
2 -139.176
```

## 7) `dnakit.population.scale_two_locus_spectrum()`

- **作用：** 双位点模型与数据之间的最优标量，以及缩放后的模型。
- **处理方法：** 调用 `TwoLocus.inference.optimal_sfs_scaling`（两者都未屏蔽的格中，数据总和除以模型总和）和 `optimally_scaled_sfs`。一个函数覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.scale_two_locus_spectrum(model, data)`.
- **输入：** 同样大小的 `model` 与 `data` 频谱。
- **输出：** `scaling` 和 `scaled_model`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    scale_two_locus_spectrum,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
data = [[[v * 3 if v is not None else None for v in row] for row in plane] for plane in spectrum]
print(round(scale_two_locus_spectrum(spectrum, data, verbose=False)["scaling"], 3))
```

- **示例输出：**

```text
3.0
```

## 8) `dnakit.population.condition_two_locus_spectrum()`

- **作用：** 按位点 A 的衍生等位基因计数对双位点频谱取条件。
- **处理方法：** 调用 `TwoLocus.numerics.condition(F, i)`：第 `[a][b]` 项是 AB 单倍型 a 个、aB 单倍型 b 个且 AB + Ab = i 时的频谱值。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.condition_two_locus_spectrum(spectrum, derived_count)`.
- **输入：** `spectrum` 和 `derived_count`（1 到 n - 1）。
- **输出：** 形状为 (derived_count + 1) x (n - derived_count + 1) 的 `spectrum`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    condition_two_locus_spectrum,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
result = condition_two_locus_spectrum(spectrum, 2, verbose=False)
print(len(result["spectrum"]), len(result["spectrum"][0]))
print([round(v, 2) for v in result["spectrum"][1]])
```

- **示例输出：**

```text
3 5
[1.2, 1.3, 1.4, 1.5, 1.6]
```

## 9) `dnakit.population.sample_two_locus_density()`

- **作用：** 从网格上的密度抽样得到期望的双位点频谱。
- **处理方法：** 在均匀网格 `linspace(0, 1, g + 1)` 上调用 `TwoLocus.numerics.sample_cached(phi, n, x)`。不提供 dadi 自带的密度向前积分：它对任何输入都抛出 `AttributeError: module 'dadi' has no attribute 'tridiag'`，补上该名称后则返回 NaN。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.sample_two_locus_density(phi, sample_size)`.
- **输入：** `phi`：边长 g + 1（至多 41）的立方体，以及 `sample_size`。
- **输出：** `spectrum`（边长 n + 1 的立方体）、`grid_points`、`sample_size`。
- **示例代码：**

```python
from dnakit.population.dadi_p22_twolocus_functions import sample_two_locus_density

phi = [[[1.0] * 6 for _ in range(6)] for _ in range(6)]
result = sample_two_locus_density(phi, 4, verbose=False)
print(result["grid_points"], result["sample_size"], round(result["spectrum"][1][1][1], 4))
```

- **示例输出：**

```text
5 4 -2.7972
```

## 10) `dnakit.population.write_multilocus_spectrum_file()`

- **作用：** 把三等位基因或双位点频谱写入 dadi 文本文件。
- **处理方法：** 按 `kind` 调用 `TriSpectrum.to_file` 或 `TLSpectrum.to_file`。文件包含注释行、含大小、折叠和外推值的头部、数据行与屏蔽行（1 表示被屏蔽）。不会覆盖已有文件。与 `read_multilocus_spectrum_file` 共同覆盖四条工作清单。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.write_multilocus_spectrum_file(spectrum, destination, *, kind='triallelic', folded_major=False, folded_ancestral=False, extrap_x=None, extrap_t=None, precision=16, comment_lines=())`.
- **输入：** `spectrum`、新的 `destination` 路径、`kind`（`triallelic` 或 `two_locus`）、折叠标志、`extrap_x`、`extrap_t`、`precision`（1 到 17）、`comment_lines`。
- **输出：** `destination`、`bytes`、`kind`、`sample_size`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population.dadi_p22_dfe_triallele_functions import triallelic_neutral_equilibrium
from dnakit.population.dadi_p22_twolocus_functions import (
    read_multilocus_spectrum_file,
    write_multilocus_spectrum_file,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
path = Path(tempfile.mkdtemp()) / "tri.fs"
written = write_multilocus_spectrum_file(spectrum, str(path), comment_lines=["neutral"], verbose=False)
print(written["kind"], written["sample_size"])
print(path.read_text().splitlines()[:2])
back = read_multilocus_spectrum_file(str(path), verbose=False)
print(back["comments"], back["folded_major"], back["spectrum"] == spectrum)
```

- **示例输出：**

```text
triallelic 6
['# neutral', '6 unfolded_major unfolded_ancestral None None']
['neutral'] False False
```

## 11) `dnakit.population.read_multilocus_spectrum_file()`

- **作用：** 从 dadi 文本文件读取三等位基因或双位点频谱。
- **处理方法：** 按 `kind`（须与文件一致）调用 `TriSpectrum.from_file` 或 `TLSpectrum.from_file`。`mask_infeasible` 会屏蔽不可能出现的项。
- **API：** `dnakit.population.dadi_p22_twolocus_functions.read_multilocus_spectrum_file(path, *, kind='triallelic', mask_infeasible=True)`.
- **输入：** `path`、`kind`、`mask_infeasible`。
- **输出：** `spectrum`、`sample_size`、折叠标志、`extrap_x`、`extrap_t`、`comments`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population.dadi_p22_twolocus_functions import (
    read_multilocus_spectrum_file,
    two_locus_spectrum_transform,
    write_multilocus_spectrum_file,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(5)] for j in range(5)] for i in range(5)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
path = Path(tempfile.mkdtemp()) / "tl.fs"
write_multilocus_spectrum_file(spectrum, str(path), kind="two_locus", extrap_x=0.01, verbose=False)
back = read_multilocus_spectrum_file(str(path), kind="two_locus", verbose=False)
print(back["sample_size"], back["folded"], back["extrap_x"], back["extrap_t"])
print(round(back["spectrum"][1][1][1], 3))
```

- **示例输出：**

```text
4 False 0.01 None
1.3
```
