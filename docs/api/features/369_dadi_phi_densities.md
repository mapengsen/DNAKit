# dadi 群体密度：构建与分裂

dadi 扩散流程中群体密度（phi）的各个步骤，输入输出为可 JSON 序列化的嵌套列表。phi 位于 dadi 默认频率网格上，每个坐标轴 `grid_points` 个点。需要 `dnakit[bio-dadi]`。

## 1) `dnakit.population.phi_one_population()`

- **作用：** 在 dadi 的频率网格上构建单个群体的平衡群体密度 phi。
- **处理方法：** 按 `model` 选择 dadi 例程并在 `Numerics.default_grid(grid_points)` 上调用：`dominance` 对应 `PhiManip.phi_1D`（数值积分），`genic` 对应 `phi_1D_genic`（解析级数），`neutral` 对应 `phi_1D_snm`，`x_linked` 对应 `phi_1D_X`。一个函数覆盖四条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.phi_one_population(*, grid_points, model='dominance', nu=1.0, theta0=1.0, gamma=0.0, h=0.5, theta=None, beta=1.0, alpha=1.0, output_path=None)`。
- **输入：** `grid_points`（8 到 400）、`model`、群体大小 `nu`、突变率 `theta0`、选择系数 `gamma`、显性度 `h`，以及可选的 `theta`、`beta`、`alpha` 和新的 `.npy` 文件 `output_path`。
- **输出：** `phi`（长度为 grid_points 的列表）、`phi_sum`、`shape`、`phi_path`、`provenance`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population

result = phi_one_population(grid_points=20, model="genic", gamma=1.0)
print(result["model"], result["shape"], round(result["phi_sum"], 3))
print([round(x, 3) for x in result["phi"][:3]])
```

- **示例输出：**

```text
genic [20] 5675.193
[2257.854, 2257.854, 681.0]
```

## 2) `dnakit.population.phi_split_population()`

- **作用：** 把密度中的一个群体分裂成两个群体。
- **处理方法：** 按 `event`（`1_to_2`、`2_to_3_split_2`、`2_to_3_split_1`）调用 `PhiManip.phi_1D_to_2D`、`phi_2D_to_3D_split_2` 或 `phi_2D_to_3D_split_1`。新群体起始频率与母群体相同。覆盖三条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.phi_split_population(phi, *, grid_points, event, output_path=None)`。
- **输入：** `phi`（嵌套列表或 `.npy` 路径）、`grid_points`、`event`，以及可选的 `output_path`。
- **输出：** 多一个坐标轴的新 `phi`、`num_populations`、`shape`、`event`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
three = phi_split_population(two["phi"], grid_points=20, event="2_to_3_split_1")
print(two["shape"], three["shape"], three["event"])
```

- **示例输出：**

```text
[20, 20] [20, 20, 20] 2_to_3_split_1
```

## 3) `dnakit.population.phi_create_admixed_population()`

- **作用：** 新增一个由现有所有群体混合形成的群体。
- **处理方法：** 对 2、3、4 个群体的密度分别调用 `PhiManip.phi_2D_to_3D_admix`、`phi_3D_to_4D`、`phi_4D_to_5D`。`proportions` 给出群体 1 到 n-1 的贡献比例，其余由最后一个群体提供。覆盖三条工作清单。
- **API：** `dnakit.population.dadi_p23_phi_functions.phi_create_admixed_population(phi, *, grid_points, proportions, output_path=None)`。
- **输入：** `phi`、`grid_points`、`proportions`（n-1 个位于 [0, 1] 且总和不超过 1 的数）、可选 `output_path`。
- **输出：** 新的 `phi`、`num_populations`、`shape`、`proportions`。
- **示例代码：**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_create_admixed_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
three = phi_create_admixed_population(two["phi"], grid_points=20, proportions=[0.3])
print(three["shape"], three["proportions"], round(three["phi_sum"], 3))
```

- **示例输出：**

```text
[20, 20, 20] [0.3] 4434397133.35
```
