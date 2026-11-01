# tskit 样本集、相关度、性状与溯祖统计

这一组 `dnakit.population` 函数封装 [tskit](https://github.com/tskit-dev/tskit) 的统计量：连锁不平衡、Y1/Y2/Y3/f4、分离位点、通用统计量、遗传相关度、性状关联、平均后代数、谱系最近邻和配对溯祖。

tskit 是可选依赖（`dnakit[bio-tskit]`），只有调用函数时才导入；缺少时报错并给出安装命令。每个函数接受显式表行（与 `build_tree_sequence()` 相同的格式）或 `.trees` 文件路径作为树序列；编辑类函数以同样的行格式返回编辑后的表。区间为半开区间 `[left, right)`。不访问网络。

## 示例数据

所有示例使用同一棵树序列：4 个样本、100 bp；第一棵树覆盖 [0, 50)，样本 (0, 1) 与 (2, 3) 成对；第二棵树覆盖 [50, 100)，(0, 2) 与 (1, 3) 成对。共 4 个位点，每个位点 1 个突变。

**示例文件：** `tables.json`

```text
{
 "sequence_length": 100,
 "nodes": [
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 1
  },
  {
   "time": 2
  },
  {
   "time": 4
  },
  {
   "time": 1
  },
  {
   "time": 3
  },
  {
   "time": 5
  }
 ],
 "edges": [
  {
   "left": 0,
   "right": 50,
   "parent": 4,
   "child": 0
  },
  {
   "left": 0,
   "right": 50,
   "parent": 4,
   "child": 1
  },
  {
   "left": 0,
   "right": 50,
   "parent": 5,
   "child": 2
  },
  {
   "left": 0,
   "right": 50,
   "parent": 5,
   "child": 3
  },
  {
   "left": 0,
   "right": 50,
   "parent": 6,
   "child": 4
  },
  {
   "left": 0,
   "right": 50,
   "parent": 6,
   "child": 5
  },
  {
   "left": 50,
   "right": 100,
   "parent": 7,
   "child": 0
  },
  {
   "left": 50,
   "right": 100,
   "parent": 7,
   "child": 2
  },
  {
   "left": 50,
   "right": 100,
   "parent": 8,
   "child": 1
  },
  {
   "left": 50,
   "right": 100,
   "parent": 8,
   "child": 3
  },
  {
   "left": 50,
   "right": 100,
   "parent": 9,
   "child": 7
  },
  {
   "left": 50,
   "right": 100,
   "parent": 9,
   "child": 8
  }
 ],
 "sites": [
  {
   "position": 10,
   "ancestral_state": "A"
  },
  {
   "position": 30,
   "ancestral_state": "C"
  },
  {
   "position": 60,
   "ancestral_state": "G"
  },
  {
   "position": 80,
   "ancestral_state": "T"
  }
 ],
 "mutations": [
  {
   "site": 0,
   "node": 4,
   "derived_state": "T"
  },
  {
   "site": 1,
   "node": 5,
   "derived_state": "G"
  },
  {
   "site": 2,
   "node": 7,
   "derived_state": "A"
  },
  {
   "site": 3,
   "node": 8,
   "derived_state": "C"
  }
 ]
}
```

## 1) `dnakit.population.tskit_ld_r2()`

- **作用：** 位点间连锁不平衡 r2，可输出矩阵、位点对或数组。
- **处理方法：** 调用 tskit 的 `tskit.LdCalculator.r2_matrix` / `r2` / `r2_array`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_ld_r2(tables, *, output='matrix', site_a=None, site_b=None, focal_site=None, direction=1, max_mutations=None, max_distance=None, max_sites=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`output`、`site_a`、`site_b`、`focal_site`、`direction`、`max_mutations`、`max_distance`、`max_sites`。
- **输出：** `values`（矩阵、数值或数组）。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_ld_r2

tables = json.load(open("tables.json"))
print(tskit_ld_r2(tables, output="pair", site_a=0, site_b=1)["values"])
print(tskit_ld_r2(tables, output="array", focal_site=0)["values"])
print(tskit_ld_r2(tables)["values"][0])
```

- **示例输出：**

```text
1.0
[1.0, 0.0, 0.0]
[1.0, 1.0, 0.0, 0.0]
```

## 2) `dnakit.population.tskit_sample_set_statistic()`

- **作用：** 多样本集统计量 Y1、Y2、Y3 与 f4。
- **处理方法：** 调用 tskit 的 `TreeSequence.Y1` / `Y2` / `Y3` / `f4`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_sample_set_statistic(tables, *, statistic, sample_sets=None, indexes=None, windows=None, mode='site', span_normalise=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`statistic`、`sample_sets`、`indexes`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_sample_set_statistic

tables = json.load(open("tables.json"))
sets = [[0], [1], [2], [3]]
print(tskit_sample_set_statistic(tables, statistic="f4", sample_sets=sets, indexes=[[0, 1, 2, 3]])["values"])
print(tskit_sample_set_statistic(tables, statistic="Y3", sample_sets=sets, indexes=[[0, 1, 2]])["values"])
print(tskit_sample_set_statistic(tables, statistic="Y1", sample_sets=[[0, 1, 2]])["values"])
```

- **示例输出：**

```text
[0.02]
[0.0]
[0.013333333333333332]
```

## 3) `dnakit.population.tskit_segregating_sites()`

- **作用：** 每个窗口与样本集的分离位点密度（或计数）。
- **处理方法：** 调用 tskit 的 `TreeSequence.segregating_sites`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_segregating_sites(tables, *, sample_sets=None, windows=None, mode='site', span_normalise=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`sample_sets`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_segregating_sites

tables = json.load(open("tables.json"))
print(tskit_segregating_sites(tables, sample_sets=[[0, 1], [2, 3]], span_normalise=False)["values"])
print(tskit_segregating_sites(tables, windows=[0, 50, 100], span_normalise=False)["values"])
```

- **示例输出：**

```text
[2.0, 2.0]
[[2.0], [2.0]]
```

## 4) `dnakit.population.tskit_general_stat()`

- **作用：** 由样本权重或样本集计数得到的窗口统计量，汇总函数取自命名列表。
- **处理方法：** 调用 tskit 的 `TreeSequence.general_stat` / `sample_count_stat`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_general_stat(tables, *, weights=None, sample_sets=None, summary='identity', polarised=False, strict=True, windows=None, mode='site', span_normalise=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`weights` 或 `sample_sets`、`summary`（identity、square、product）、`polarised`、`strict`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_general_stat

tables = json.load(open("tables.json"))
w = [[1.0], [-1.0], [2.0], [-2.0]]
print(tskit_general_stat(tables, weights=w, summary="square", polarised=True)["values"])
print(tskit_general_stat(tables, sample_sets=[[0, 1], [2, 3]], summary="product", polarised=True, strict=False)["values"])
```

- **示例输出：**

```text
[0.18]
[0.02]
```

## 5) `dnakit.population.tskit_genetic_relatedness()`

- **作用：** 样本集之间的遗传相关度，输出成对值或完整矩阵。
- **处理方法：** 调用 tskit 的 `TreeSequence.genetic_relatedness` / `genetic_relatedness_matrix`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_genetic_relatedness(tables, *, sample_sets=None, indexes=None, windows=None, mode='site', span_normalise=True, centre=True, proportion=True, polarised=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`sample_sets`、`indexes`、`centre`、`proportion`、`polarised`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness

tables = json.load(open("tables.json"))
sets = [[0], [1], [2], [3]]
print(tskit_genetic_relatedness(tables, sample_sets=sets, indexes=[[0, 1], [0, 2]])["values"])
print(tskit_genetic_relatedness(tables, sample_sets=sets)["values"])
```

- **示例输出：**

```text
[0.0, 0.0]
[[0.01, -0.0, -0.0, -0.009999999999999998], [-0.0, 0.01, -0.009999999999999998, -0.0], [-0.0, -0.009999999999999998, 0.01, -0.0], [-0.009999999999999998, -0.0, -0.0, 0.01]]
```

## 6) `dnakit.population.tskit_genetic_relatedness_vector()`

- **作用：** 遗传相关矩阵与权重矩阵的乘积（枝长模式）。
- **处理方法：** 调用 tskit 的 `TreeSequence.genetic_relatedness_vector`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_genetic_relatedness_vector(tables, *, weights, windows=None, span_normalise=True, centre=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`weights`（每个样本一行）、`centre`、`windows`、`span_normalise`。
- **输出：** `values`，每个窗口形状与 `weights` 相同。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness_vector

tables = json.load(open("tables.json"))
r = tskit_genetic_relatedness_vector(tables, weights=[[1.0], [2.0], [0.5], [3.0]])
print(r["values"])
```

- **示例输出：**

```text
[[-4.03125], [2.59375], [-4.46875], [5.90625]]
```

## 7) `dnakit.population.tskit_genetic_relatedness_weighted()`

- **作用：** 权重列两两之间的遗传相关度。
- **处理方法：** 调用 tskit 的 `TreeSequence.genetic_relatedness_weighted`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_genetic_relatedness_weighted(tables, *, weights, indexes=None, windows=None, mode='site', span_normalise=True, polarised=False, centre=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`weights`、`indexes`、`polarised`、`centre`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness_weighted

tables = json.load(open("tables.json"))
r = tskit_genetic_relatedness_weighted(tables, weights=[[1.0, 0.5], [2.0, 1.0], [0.5, 3.0], [3.0, 2.0]], indexes=[[0, 1]])
print(r["values"])
```

- **示例输出：**

```text
[0.0]
```

## 8) `dnakit.population.tskit_trait_statistic()`

- **作用：** 基因型与性状的关联：相关系数、协方差或线性模型。
- **处理方法：** 调用 tskit 的 `TreeSequence.trait_correlation` / `trait_covariance` / `trait_linear_model`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_trait_statistic(tables, *, statistic, weights, covariates=None, windows=None, mode='site', span_normalise=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`statistic`、`weights`（样本 × 性状）、`covariates`、`windows`、`mode`、`span_normalise`。
- **输出：** `values`（窗口 × 集合或索引元组）以及 `mode`、`windows`、`span_normalise`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_trait_statistic

tables = json.load(open("tables.json"))
w = [[1.0], [2.0], [0.5], [3.0]]
for stat in ("correlation", "covariance", "linear_model"):
    print(stat, tskit_trait_statistic(tables, statistic=stat, weights=w)["values"])
```

- **示例输出：**

```text
correlation [0.01694915254237288]
covariance [0.006944444444444444]
linear_model [0.0625]
```

## 9) `dnakit.population.tskit_mean_descendants()`

- **作用：** 每个节点对每个样本集的平均后代样本数。
- **处理方法：** 调用 tskit 的 `TreeSequence.mean_descendants`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_mean_descendants(tables, *, sample_sets)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`sample_sets`。
- **输出：** `values`，每个节点一行、每个样本集一列。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_mean_descendants

tables = json.load(open("tables.json"))
r = tskit_mean_descendants(tables, sample_sets=[[0, 1], [2, 3]])
for row in r["values"]:
    print(row)
```

- **示例输出：**

```text
[1.0, 0.0]
[1.0, 0.0]
[0.0, 1.0]
[0.0, 1.0]
[2.0, 0.0]
[0.0, 2.0]
[2.0, 2.0]
[1.0, 1.0]
[1.0, 1.0]
[2.0, 2.0]
```

## 10) `dnakit.population.tskit_genealogical_nearest_neighbours()`

- **作用：** 每个焦点节点的谱系最近邻落在各样本集中的比例。
- **处理方法：** 调用 tskit 的 `TreeSequence.genealogical_nearest_neighbours`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_genealogical_nearest_neighbours(tables, *, focal, sample_sets)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`focal`、`sample_sets`（互不重叠）。
- **输出：** `values`，每个焦点节点一行。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genealogical_nearest_neighbours

tables = json.load(open("tables.json"))
r = tskit_genealogical_nearest_neighbours(tables, focal=[0, 3], sample_sets=[[1], [2]])
print(r["values"])
```

- **示例输出：**

```text
[[0.5, 0.5], [0.5, 0.5]]
```

## 11) `dnakit.population.tskit_pair_coalescence()`

- **作用：** 配对溯祖的计数、分位数或速率随时间的变化。
- **处理方法：** 调用 tskit 的 `TreeSequence.pair_coalescence_counts` / `pair_coalescence_quantiles` / `pair_coalescence_rates`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_pair_coalescence(tables, *, output='rates', time_windows='nodes', quantiles=None, sample_sets=None, indexes=None, windows=None, span_normalise=True, pair_normalise=False)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`output`、`time_windows`、`quantiles`、`sample_sets`、`indexes`、`windows`、`span_normalise`、`pair_normalise`。
- **输出：** `values` 以及所用的 `time_windows` 或 `quantiles`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_pair_coalescence

tables = json.load(open("tables.json"))
sets = [[0, 1, 2, 3]]
print(tskit_pair_coalescence(tables, output="counts", sample_sets=sets)["values"])
print(tskit_pair_coalescence(tables, output="rates", time_windows=[0, 2], sample_sets=sets)["values"])
print(tskit_pair_coalescence(tables, output="quantiles", quantiles=[0.5, 1.0], sample_sets=sets)["values"])
```

- **示例输出：**

```text
[0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 2.0, 0.5, 0.5, 2.0]
[0.09116077839697735, 0.4761904761904763]
[4.0, 5.0]
```
