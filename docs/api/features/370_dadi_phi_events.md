# dadi 群体密度：混合与群体选择

dadi 扩散流程中群体密度（phi）的各个步骤，输入输出为可 JSON 序列化的嵌套列表。phi 位于 dadi 默认频率网格上，每个坐标轴 `grid_points` 个点。需要 `dnakit[bio-dadi]`。

## 1) `dnakit.population.phi_admix_into_population()`

- **作用：** 用其他群体的混合替换密度中的某一个群体。
- **处理方法：** 按群体数与目标群体选择 dadi 例程：`phi_2D_admix_1_into_2` / `phi_2D_admix_2_into_1`、`phi_3D_admix_*_into_*`、`phi_4D_admix_into_1` 到 `_4`、`phi_5D_admix_into_1` 到 `_5`。比例按其余群体的升序给出。覆盖十四条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.phi_admix_into_population(phi, *, grid_points, target_population, proportions, output_path=None)`。
- **输入：** `phi`（2 到 5 个群体）、`grid_points`、从 1 开始的 `target_population`、`proportions`（n-1 个数）、可选 `output_path`。
- **输出：** 混合后的 `phi`（形状不变）、`target_population`、`proportions`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_admix_into_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
mixed = phi_admix_into_population(two["phi"], grid_points=20, target_population=2, proportions=[0.2])
print(mixed["shape"], mixed["target_population"], round(mixed["phi_sum"], 3))
```

- **示例输出：**

```text
[20, 20] 2 3554632.261
```

## 2) `dnakit.population.phi_select_populations()`

- **作用：** 删除、保留或重排密度中的群体。
- **处理方法：** `remove` 调用 `PhiManip.remove_pop`，`keep` 调用 `filter_pops`，`reorder` 调用 `reorder_pops`。覆盖三条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.phi_select_populations(phi, *, grid_points, operation, populations, output_path=None)`。
- **输入：** `phi`、`grid_points`、`operation`、从 1 开始的 `populations`、可选 `output_path`。
- **输出：** 新的 `phi`、`num_populations`、`shape`、`select`、`populations`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_select_populations

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
for op, pops in (("remove", [2]), ("keep", [2]), ("reorder", [2, 1])):
    out = phi_select_populations(two["phi"], grid_points=20, operation=op, populations=pops)
    print(op, out["shape"])
```

- **示例输出：**

```text
remove [20]
keep [20]
reorder [20, 20]
```
