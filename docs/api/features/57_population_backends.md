# 可选 tskit 与 dadi 群体遗传后端

这一组函数是 `dnakit.population` 中基于可选依赖的适配器。[tskit](https://tskit.dev/tskit) 负责由显式表格构建树序列、读写 `.trees` 文件，并计算位点模式的多样性、Tajima's D、Fst 和位点频率谱；[dadi](https://github.com/RyanGutenkunst/dadi) 负责恒定群体大小中性模型的期望频率谱、频率谱折叠，以及观测频率谱的 Poisson 对数似然。这里不做序列模拟，也不拟合人口历史模型。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

第 1 至 7 项的 tskit 示例读取下面的表格。它有 4 个样本、1 棵树和 3 个位点，只用于演示接口，统计值不代表真实群体。第 8 至 11 项的 dadi 示例直接使用列表，不需要文件。

**示例文件：** `tables.json`

```json
{
  "sequence_length": 100,
  "nodes": [{"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 1}, {"time": 2}, {"time": 3}],
  "edges": [{"left": 0, "right": 100, "parent": 4, "child": 0}, {"left": 0, "right": 100, "parent": 4, "child": 1}, {"left": 0, "right": 100, "parent": 5, "child": 2}, {"left": 0, "right": 100, "parent": 5, "child": 3}, {"left": 0, "right": 100, "parent": 6, "child": 4}, {"left": 0, "right": 100, "parent": 6, "child": 5}],
  "sites": [{"position": 10, "ancestral_state": "A"}, {"position": 20, "ancestral_state": "A"}, {"position": 30, "ancestral_state": "A"}],
  "mutations": [{"site": 0, "node": 0, "derived_state": "G"}, {"site": 1, "node": 4, "derived_state": "G"}, {"site": 2, "node": 5, "derived_state": "G"}]
}
```
## 1) `dnakit.population.build_tree_sequence()`

- **作用：** 由显式的节点、边、位点和突变行构建 tskit 树序列，并返回汇总信息。
- **处理方法：** 把各表的行写入 tskit 的 `TableCollection`，调用 `sort()` 排序，再用 `tree_sequence()` 生成树序列。区间是半开区间 `[left, right)`，坐标单位与 `sequence_length` 相同。tskit 检查不通过时报错。
- **API：** `dnakit.population.build_tree_sequence(tables, *, max_rows=100000)`。
- **输入：** `tables` 是映射，包含 `sequence_length`（大于 0）、`nodes`（每行有 `time`，可选 `sample`）、`edges`（`left`、`right`、`parent`、`child`）、`sites`（`position`、`ancestral_state`）和 `mutations`（`site`、`node`、`derived_state`，可选 `parent`、`time`）。`max_rows` 限制每张表的行数。
- **输出：** 含 `num_nodes`、`num_edges`、`num_samples`、`num_trees`、`num_sites`、`num_mutations`、`num_populations`、`sequence_length` 和 `provenance` 的字典。
- **示例代码：**

```python
import json
from dnakit.population import build_tree_sequence

tables = json.load(open("tables.json"))
result = build_tree_sequence(tables)
print(result["num_trees"], result["num_samples"], result["num_sites"])
```

- **示例输出：**

```text
1 4 3
```

## 2) `dnakit.population.write_tree_sequence()`

- **作用：** 由表格构建树序列，并写成 tskit 的 `.trees` 文件。
- **处理方法：** 与 `build_tree_sequence()` 使用相同的构建和校验规则，再调用 `TreeSequence.dump()`。已存在的文件不会被覆盖，父目录必须已经存在。
- **API：** `dnakit.population.write_tree_sequence(tables, destination, *, max_rows=100000)`。
- **输入：** 与 `build_tree_sequence()` 相同的 `tables`；`destination` 是输出路径。
- **输出：** 写入路径、文件字节数，以及与 `build_tree_sequence()` 相同的汇总字段。
- **示例代码：**

```python
import json
from pathlib import Path
from dnakit.population import write_tree_sequence

Path("sample.trees").unlink(missing_ok=True)  # the function never overwrites an existing file
tables = json.load(open("tables.json"))
result = write_tree_sequence(tables, "sample.trees")
print(result["destination"].endswith("sample.trees"), result["num_mutations"])
```

- **示例输出：**

```text
True 3
```

## 3) `dnakit.population.read_tree_sequence()`

- **作用：** 读取 tskit 的 `.trees` 文件，并返回汇总信息。
- **处理方法：** 先检查文件大小是否超过 `max_bytes`，再调用 `tskit.load()`。不是 tskit 格式的文件会报错。
- **API：** `dnakit.population.read_tree_sequence(source, *, max_bytes=100000000)`。
- **输入：** `.trees` 文件路径；`max_bytes` 是文件大小上限。
- **输出：** 含 `source`、`bytes` 以及与 `build_tree_sequence()` 相同汇总字段的字典。
- **示例代码：**

```python
from dnakit.population import read_tree_sequence

result = read_tree_sequence("sample.trees")
print(result["num_trees"], result["num_edges"])
```

- **示例输出：**

```text
1 6
```

## 4) `dnakit.population.calculate_diversity()`

- **作用：** 计算每个样本集的核苷酸多样性（每个碱基的 π）。
- **处理方法：** 调用 tskit 的 `diversity()`，使用位点模式，结果除以 `sequence_length`。未指定 `sample_sets` 时，全部样本作为名为 `"all"` 的一个集合。
- **API：** `dnakit.population.calculate_diversity(tables, *, sample_sets=None, max_rows=100000)`。
- **输入：** `tables`；`sample_sets` 是 `{名称: [样本节点编号, ...]}`。每个集合至少需要两个不重复的样本节点。
- **输出：** `diversity_per_base`（每个集合的值）、`sequence_length` 和 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population import calculate_diversity

tables = json.load(open("tables.json"))
result = calculate_diversity(tables, sample_sets={"A": [0, 1], "B": [2, 3]})
print({name: round(value, 6) for name, value in result["diversity_per_base"].items()})
```

- **示例输出：**

```text
{'A': 0.01, 'B': 0.0}
```

## 5) `dnakit.population.calculate_tajimas_d()`

- **作用：** 计算每个样本集的 Tajima's D 和分离位点数。
- **处理方法：** 调用 tskit 的 `Tajimas_D()` 和 `segregating_sites()`。tskit 返回 NaN 时（例如没有分离位点），结果记为 `None`。
- **API：** `dnakit.population.calculate_tajimas_d(tables, *, sample_sets=None, max_rows=100000)`。
- **输入：** 与 `calculate_diversity()` 相同。
- **输出：** `tajimas_d`、`segregating_sites` 和 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population import calculate_tajimas_d

tables = json.load(open("tables.json"))
result = calculate_tajimas_d(tables)
print(result["segregating_sites"]["all"], round(result["tajimas_d"]["all"], 4))
```

- **示例输出：**

```text
3 1.0898
```

## 6) `dnakit.population.calculate_fst()`

- **作用：** 计算两个样本集之间的 Fst。
- **处理方法：** 调用 tskit 的 `Fst()`，使用位点模式。tskit 的估计量为 `1 - 2(d_X + d_Y) / (d_X + 2 d_XY + d_Y)`，其中 `d_X`、`d_Y` 是组内多样性，`d_XY` 是组间分歧。分母为零时结果记为 `None`。
- **API：** `dnakit.population.calculate_fst(tables, *, sample_sets, max_rows=100000)`。
- **输入：** `tables`；`sample_sets` 必须恰好包含两个集合。
- **输出：** `fst`、两个群体的名称和样本数，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population import calculate_fst

tables = json.load(open("tables.json"))
result = calculate_fst(tables, sample_sets={"A": [0, 1], "B": [2, 3]})
print(round(result["fst"], 4))
```

- **示例输出：**

```text
0.6667
```

## 7) `dnakit.population.calculate_site_frequency_spectrum()`

- **作用：** 计算每个样本集的位点频率谱计数。
- **处理方法：** 调用 tskit 的 `allele_frequency_spectrum()`，输出计数而非频率。`spectrum[i]` 是衍生等位基因在该集合 `i` 个样本中出现的位点数。`folded=True` 时输出次等位基因频谱，把 `i` 与 `n - i` 的计数合并到 `min(i, n - i)`。
- **API：** `dnakit.population.calculate_site_frequency_spectrum(tables, *, sample_sets=None, folded=False, max_rows=100000)`。
- **输入：** `tables`；`sample_sets` 与 `calculate_diversity()` 相同；`folded` 取 `True` 或 `False`。
- **输出：** `spectra`（每个集合一个长度为 `n + 1` 的列表）、`folded` 和 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population import calculate_site_frequency_spectrum

tables = json.load(open("tables.json"))
result = calculate_site_frequency_spectrum(tables)
print(result["spectra"]["all"])
```

- **示例输出：**

```text
[0, 1, 2, 0, 0]
```

## 8) `dnakit.population.expected_neutral_sfs()`

- **作用：** 计算恒定群体大小的标准中性模型下，期望的位点频率谱。
- **处理方法：** 调用 dadi 的 `Demographics1D.snm()`，在 `grid_points` 个网格点上做数值积分，期望值近似 `theta / i`。结果乘以 `theta`。`folded=True` 时调用 dadi 的 `fold()`。第 0 类和第 n 类是单态位点，dadi 会将其掩蔽，输出为 `None`。网格越密，误差越小。
- **API：** `dnakit.population.expected_neutral_sfs(sample_size, *, theta=1.0, folded=False, grid_points=1000, max_sample_size=1000)`。
- **输入：** 样本染色体数 `sample_size`（至少为 2）；`theta` 为大于零的数；`grid_points` 取 50 至 10000 的整数。
- **输出：** `spectrum`（长度为 `sample_size + 1` 的列表，掩蔽项为 `None`）、`theta`、`folded` 和 `provenance`。
- **示例代码：**

```python
from dnakit.population import expected_neutral_sfs

result = expected_neutral_sfs(4, theta=1.0)
print([None if value is None else round(value, 3) for value in result["spectrum"]])
```

- **示例输出：**

```text
[None, 1.0, 0.5, 0.333, None]
```

## 9) `dnakit.population.fold_sfs()`

- **作用：** 把非折叠的位点频率谱折叠为次等位基因频率谱。
- **处理方法：** 用 dadi 的 `Spectrum(...).fold()` 折叠。第 0 类和第 n 类被掩蔽，不参与折叠。折叠后被掩蔽的位置输出为 `None`。
- **API：** `dnakit.population.fold_sfs(spectrum, *, max_sample_size=1000)`。
- **输入：** 长度为 `n + 1` 的非负计数列表，`n` 至少为 2。
- **输出：** `sample_size`、`spectrum`（折叠后的列表，掩蔽项为 `None`）和 `provenance`。
- **示例代码：**

```python
from dnakit.population import fold_sfs

result = fold_sfs([0, 2, 1, 1, 0])
print(result["spectrum"])
```

- **示例输出：**

```text
[None, 3.0, 1.0, None, None]
```

## 10) `dnakit.population.estimate_sfs_theta()`

- **作用：** 在中性模型下，估计观测频率谱的最大似然 `theta`。
- **处理方法：** 用 dadi 的 `Inference.optimal_sfs_scaling()`，只使用第 1 至 n-1 类。理论上它等于观测总数除以模型总数。
- **API：** `dnakit.population.estimate_sfs_theta(observed, *, grid_points=1000, max_sample_size=1000)`。
- **输入：** 观测计数 `observed`，长度为 `n + 1`；第 1 至 n-1 类之和必须大于零。
- **输出：** `sample_size`、`theta` 和 `provenance`。
- **示例代码：**

```python
from dnakit.population import estimate_sfs_theta

result = estimate_sfs_theta([0, 2, 1, 1, 0])
print(round(result["theta"], 3))
```

- **示例输出：**

```text
2.182
```

## 11) `dnakit.population.calculate_sfs_log_likelihood()`

- **作用：** 计算观测频率谱在中性模型下的 Poisson 对数似然。
- **处理方法：** 用 dadi 的 `Inference.ll()`，只使用第 1 至 n-1 类，包含 `log(x!)` 项。`theta` 省略时，先取最大似然 `theta`，与 dadi 的 `ll_multinom()` 结果一致。
- **API：** `dnakit.population.calculate_sfs_log_likelihood(observed, *, theta=None, grid_points=1000, max_sample_size=1000)`。
- **输入：** 与 `estimate_sfs_theta()` 相同的 `observed`；`theta` 可省略，也可以给出大于零的数。
- **输出：** `theta`、`theta_source`（`"maximum likelihood"` 或 `"given"`）、`log_likelihood` 和 `provenance`。
- **示例代码：**

```python
from dnakit.population import calculate_sfs_log_likelihood

result = calculate_sfs_log_likelihood([0, 2, 1, 1, 0])
print(result["theta_source"], round(result["log_likelihood"], 3))
```

- **示例输出：**

```text
maximum likelihood -3.364
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-tskit]"    # tskit，build/read/write_tree_sequence 与 calculate_* 群体统计
python -m pip install "dnakit[bio-dadi]"     # dadi，expected_neutral_sfs 等频率谱函数
```

项目地址：[tskit](https://tskit.dev/tskit)、[dadi](https://github.com/RyanGutenkunst/dadi)。

## 验证方式

测试文件是 `others/tests/unit/test_population_backends.py`，共 42 项，主要检查：

- tskit 的树序列汇总、多样性、Fst 和位点频率谱与直接调用 tskit 的结果一致；多样性、Tajima's D 和 Fst 的数值与手工按公式计算的结果一致。Tajima's D 另与按 Tajima (1989) 公式精确计算的结果比较。
- 位点行的顺序颠倒时，tskit 排序后的统计结果不变；`write_tree_sequence()` 写出的文件可由 `read_tree_sequence()` 读回。
- dadi 的中性期望频率谱与 `theta / i` 比较（相对误差小于 1e-4），并与直接调用 `Demographics1D.snm()` 一致。折叠结果与手工折叠一致；`theta` 的最大似然估计与 `24/11` 一致；对数似然与 Poisson 公式及 `Inference.ll_multinom()` 一致。
- 布尔值、越界值、未知键、非样本节点、超过 `max_rows` 或 `max_sample_size`，以及已存在的输出文件都会被拒绝。
- 缺少 `tskit` 或 `dadi` 时，报错信息包含对应的安装额外项名称。

## 当前范围

这一批只覆盖 tskit 的树序列读写、位点模式的多样性、Tajima's D、Fst（两个群体）和位点频率谱，以及 dadi 的恒定大小中性模型。树序列只能来自显式表格或已有的 `.trees` 文件，不做模拟（msprime 未纳入）。

未包装的 tskit 功能包括分支模式和窗口统计、divergence、连锁不平衡、基因型矩阵、树的遍历，以及 VCF 导出。未包装的 dadi 功能包括人口历史模型、频率谱数据字典读取、参数优化、Godambe 方差和二维频率谱。dadi 的期望频率谱是数值积分结果，误差随 `grid_points` 增大而减小。
