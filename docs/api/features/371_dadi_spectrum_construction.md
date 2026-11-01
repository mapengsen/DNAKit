# dadi 频谱构建

dadi 扩散流程中群体密度（phi）的各个步骤，输入输出为可 JSON 序列化的嵌套列表。phi 位于 dadi 默认频率网格上，每个坐标轴 `grid_points` 个点。需要 `dnakit[bio-dadi]`。

## 1) `dnakit.population.spectrum_from_phi()`

- **作用：** 从群体密度中抽样得到等位基因频率谱。
- **处理方法：** 调用 `Spectrum.from_phi`；给出 `inbreeding_coefficients` 与 `ploidies` 时改为调用 `Spectrum.from_phi_inbreeding`。覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.spectrum_from_phi(phi, *, grid_points, sample_sizes, inbreeding_coefficients=None, ploidies=None, mask_corners=True, pop_ids=None, admix_props=None, het_ascertained=None, force_direct=False)`。
- **输入：** `phi`、`grid_points`、`sample_sizes`（每群体一个）、可选近交系数与倍性、`mask_corners`、`pop_ids`、`admix_props`、`het_ascertained`、`force_direct`。
- **输出：** `spectrum`（嵌套列表，被遮蔽项为 None）、`sample_sizes`、`segregating_sites`、`mode`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, spectrum_from_phi

phi = phi_one_population(grid_points=30, model="neutral")["phi"]
plain = spectrum_from_phi(phi, grid_points=30, sample_sizes=[6])
inbred = spectrum_from_phi(phi, grid_points=30, sample_sizes=[6], inbreeding_coefficients=[0.2], ploidies=[2])
print(plain["mode"], [round(x, 3) for x in plain["spectrum"][1:-1]])
print(inbred["mode"], [round(x, 3) for x in inbred["spectrum"][1:-1]])
```

- **示例输出：**

```text
from_phi [1.039, 0.514, 0.339, 0.253, 0.201]
from_phi_inbreeding [0.881, 0.512, 0.327, 0.249, 0.194]
```

## 2) `dnakit.population.scramble_spectrum_pop_ids()`

- **作用：** 在频率谱的各群体之间随机重新分配所抽样的个体。
- **处理方法：** 在固定随机种子下调用 `Spectrum.scramble_pop_ids`（结束后恢复全局随机数状态）。
- **API：** `dnakit.population.dadi_p23_phi_functions.scramble_spectrum_pop_ids(spectrum, *, mask_corners=True, seed=0)`。
- **输入：** 至少两个群体的 `spectrum`、`mask_corners`、`seed`。
- **输出：** `spectrum`、`sample_sizes`、`segregating_sites`、`seed`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import scramble_spectrum_pop_ids

spectrum = [[None, 5.0, 3.0, 1.0], [4.0, 3.0, 2.0, 1.0], [2.0, 2.0, 1.0, 1.0], [1.0, 1.0, 1.0, None]]
first = scramble_spectrum_pop_ids(spectrum, seed=1)
second = scramble_spectrum_pop_ids(spectrum, seed=1)
print(first["sample_sizes"], round(first["segregating_sites"], 3), first["spectrum"] == second["spectrum"])
```

- **示例输出：**

```text
[3, 3] 28.0 True
```

## 3) `dnakit.population.spectrum_from_snp_data()`

**示例文件：** `toy.snp`

```text
Context Out Allele1 A B Allele2 A B ID
ACA ACA C 4 3 T 2 3 s1
CCG CCG C 5 2 T 1 4 s2
GCT GCT G 3 3 A 3 3 s3
TCA TCA C 6 1 T 0 5 s4
ACG ACG C 2 2 T 4 4 s5
```

- **作用：** 由 dadi SNP 文件构建等位基因频率谱，可选校正祖先等位基因误判。
- **处理方法：** 用 `Misc.make_data_dict` 读取文件后调用 `Spectrum.from_data_dict`；给出 `misidentification_path` 时调用 `Spectrum.from_data_dict_corrected`（Hernandez 等 2007）。覆盖两条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.spectrum_from_snp_data(path, *, pop_ids, projections, polarized=True, mask_corners=True, misidentification_path=None, force_positive=True)`。
- **输入：** SNP 文件路径、`pop_ids`、`projections`、`polarized`、`mask_corners`、可选误判概率文件、`force_positive`。
- **输出：** `spectrum`、`sample_sizes`、`segregating_sites`、`folded`、`mode`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import spectrum_from_snp_data

result = spectrum_from_snp_data("toy.snp", pop_ids=["A", "B"], projections=[4, 4])
print(result["mode"], result["sample_sizes"], round(result["segregating_sites"], 3))
folded = spectrum_from_snp_data("toy.snp", pop_ids=["A"], projections=[6], polarized=False)
print(folded["folded"], round(folded["segregating_sites"], 3))
```

- **示例输出：**

```text
from_data_dict [4, 4] 3.996
True 4.0
```
