# tskit 树序列编辑与突变表

这一组 `dnakit.population` 函数封装 [tskit](https://github.com/tskit-dev/tskit) 的编辑操作：表清理、按区间与时间裁剪、修剪、化简、取子集、合并与拼接、边拆分、单倍型延伸以及简约法突变映射。

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

## 1) `dnakit.population.tskit_link_ancestors()`

- **作用：** 把指定祖先节点的祖先片段链接到继承它们的样本。
- **处理方法：** 调用 tskit 的 `TableCollection.link_ancestors`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_link_ancestors(tables, *, samples, ancestors)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`samples`、`ancestors`（节点 ID）。
- **输出：** `num_edges` 与带 `left`、`right`、`parent`、`child` 的 `edges`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_link_ancestors

tables = json.load(open("tables.json"))
r = tskit_link_ancestors(tables, samples=[0, 1], ancestors=[4, 6])
print(r["num_edges"])
for e in r["edges"]:
    print(e)
```

- **示例输出：**

```text
3
{'left': 0.0, 'right': 50.0, 'parent': 4, 'child': 0}
{'left': 0.0, 'right': 50.0, 'parent': 4, 'child': 1}
{'left': 0.0, 'right': 50.0, 'parent': 6, 'child': 4}
```

## 2) `dnakit.population.tskit_compute_mutation_parents()`

- **作用：** 计算每个突变的父突变。
- **处理方法：** 调用 tskit 的 `TableCollection.compute_mutation_parents`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_compute_mutation_parents(tables)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）。
- **输出：** `mutation_parents`（-1 表示没有）和编辑后的表。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_parents

tables = json.load(open("tables.json"))
r = tskit_compute_mutation_parents(tables)
print(r["mutation_parents"])
```

- **示例输出：**

```text
[-1, -1, -1, -1]
```

## 3) `dnakit.population.tskit_compute_mutation_times()`

- **作用：** 由节点时间计算或补全突变时间（`method="tables"` 或 `"min"`）。
- **处理方法：** 调用 tskit 的 `TableCollection.compute_mutation_times` / `TreeSequence.impute_unknown_mutations_time`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_compute_mutation_times(tables, *, method='tables')`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`method`。
- **输出：** `mutation_times`；`method="tables"` 时还有编辑后的表。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_times

tables = json.load(open("tables.json"))
for method in ("tables", "min"):
    r = tskit_compute_mutation_times(tables, method=method)
    print(method, r["mutation_times"])
```

- **示例输出：**

```text
tables [2.5, 3.0, 3.0, 4.0]
min [1.0, 2.0, 1.0, 3.0]
```

## 4) `dnakit.population.tskit_deduplicate_sites()`

- **作用：** 合并位于同一位置的位点。
- **处理方法：** 调用 tskit 的 `TableCollection.deduplicate_sites`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_deduplicate_sites(tables)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）（按表读取，因为树序列不能含重复位点）。
- **输出：** `sites_before`、`sites_removed` 和编辑后的表。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_deduplicate_sites

tables = json.load(open("tables.json"))
dup = dict(tables, sites=tables["sites"] + [{"position": 10, "ancestral_state": "A"}], mutations=tables["mutations"] + [{"site": 4, "node": 1, "derived_state": "G"}])
r = tskit_deduplicate_sites(dup)
print(r["sites_before"], r["sites_removed"], r["result_summary"]["num_mutations"])
```

- **示例输出：**

```text
5 1 5
```

## 5) `dnakit.population.tskit_delete_sites()`

- **作用：** 按 ID 删除位点及其突变。
- **处理方法：** 调用 tskit 的 `TreeSequence.delete_sites`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_delete_sites(tables, *, site_ids)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`site_ids`。
- **输出：** `deleted` 和编辑后的表。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_delete_sites

tables = json.load(open("tables.json"))
r = tskit_delete_sites(tables, site_ids=[0, 3])
print(r["deleted"], r["result_summary"]["num_sites"], [s["position"] for s in r["tables"]["sites"]])
```

- **示例输出：**

```text
[0, 3] 2 [30.0, 60.0]
```

## 6) `dnakit.population.tskit_select_intervals()`

- **作用：** 删除或保留基因组区间（`action="delete"` 或 `"keep"`）。
- **处理方法：** 调用 tskit 的 `TreeSequence.delete_intervals` / `keep_intervals`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_select_intervals(tables, *, intervals, action='delete', simplify=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`intervals`（[left, right] 对）、`action`、`simplify`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_select_intervals

tables = json.load(open("tables.json"))
for action in ("delete", "keep"):
    r = tskit_select_intervals(tables, intervals=[[0, 50]], action=action)
    print(action, r["result_summary"]["num_trees"], r["result_summary"]["num_sites"])
```

- **示例输出：**

```text
delete 2 2
keep 2 2
```

## 7) `dnakit.population.tskit_remove_older_genealogy()`

- **作用：** 移除早于某时间的谱系，可选是否添加新根。
- **处理方法：** 调用 tskit 的 `TreeSequence.decapitate` / `TableCollection.delete_older`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_remove_older_genealogy(tables, *, time, add_roots=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`time`、`add_roots`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_remove_older_genealogy

tables = json.load(open("tables.json"))
r = tskit_remove_older_genealogy(tables, time=3.5)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_nodes"])
r = tskit_remove_older_genealogy(tables, time=3.5, add_roots=False)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_nodes"])
```

- **示例输出：**

```text
12 14
8 10
```

## 8) `dnakit.population.tskit_trim_genome()`

- **作用：** 裁剪基因组左侧、右侧或两侧的空白末端。
- **处理方法：** 调用 tskit 的 `TreeSequence.ltrim` / `rtrim` / `trim`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_trim_genome(tables, *, side='both')`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`side`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_trim_genome

tables = json.load(open("tables.json"))
from dnakit.population.tskit_p31_functions import tskit_shift_coordinates

shifted = tskit_shift_coordinates(tables, value=10)["tables"]
print(shifted["sequence_length"], shifted["edges"][0]["left"])
r = tskit_trim_genome(shifted, side="left")
print(r["tables"]["sequence_length"], r["tables"]["edges"][0]["left"])
```

- **示例输出：**

```text
110.0 10.0
100.0 0.0
```

## 9) `dnakit.population.tskit_shift_coordinates()`

- **作用：** 把所有基因组坐标整体平移一个常数。
- **处理方法：** 调用 tskit 的 `TreeSequence.shift`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_shift_coordinates(tables, *, value, sequence_length=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`value`、`sequence_length`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_shift_coordinates

tables = json.load(open("tables.json"))
r = tskit_shift_coordinates(tables, value=10)
print(r["tables"]["sequence_length"], r["tables"]["sites"][0]["position"])
```

- **示例输出：**

```text
110.0 20.0
```

## 10) `dnakit.population.tskit_simplify()`

- **作用：** 把谱系化简为所选样本的历史。
- **处理方法：** 调用 tskit 的 `TreeSequence.simplify`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_simplify(tables, *, samples=None, reduce_to_site_topology=False, keep_unary=False, keep_input_roots=False, filter_sites=True, filter_nodes=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`samples` 以及 tskit 的 simplify 开关。
- **输出：** 编辑后的表、`node_map` 和 `options`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_simplify

tables = json.load(open("tables.json"))
r = tskit_simplify(tables, samples=[0, 1])
print(r["result_summary"]["num_nodes"], r["node_map"])
```

- **示例输出：**

```text
4 [0, 1, -1, -1, 2, -1, -1, -1, -1, 3]
```

## 11) `dnakit.population.tskit_subset()`

- **作用：** 只保留所列节点及它们之间的内容。
- **处理方法：** 调用 tskit 的 `TreeSequence.subset`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_subset(tables, *, nodes, reorder_populations=True, remove_unreferenced=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`nodes`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_subset

tables = json.load(open("tables.json"))
r = tskit_subset(tables, nodes=[0, 1, 4])
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **示例输出：**

```text
3 2
```

## 12) `dnakit.population.tskit_union()`

- **作用：** 按节点对应关系并入第二棵树序列的非共享部分。
- **处理方法：** 调用 tskit 的 `TreeSequence.union`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_union(tables, other, *, node_mapping, all_edges=False, all_mutations=False, check_shared_equality=True)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`other`、`node_mapping`（`other` 每个节点对应 `tables` 的节点，-1 表示新节点）。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_union

tables = json.load(open("tables.json"))
r = tskit_union(tables, tables, node_mapping=list(range(10)))
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **示例输出：**

```text
10 12
```

## 13) `dnakit.population.tskit_concatenate()`

- **作用：** 沿基因组首尾连接多棵树序列。
- **处理方法：** 调用 tskit 的 `TreeSequence.concatenate`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_concatenate(tables, others)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`others`（树序列列表）。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_concatenate

tables = json.load(open("tables.json"))
r = tskit_concatenate(tables, [tables])
print(r["result_summary"]["sequence_length"], r["result_summary"]["num_trees"], r["result_summary"]["num_sites"])
```

- **示例输出：**

```text
200.0 4 8
```

## 14) `dnakit.population.tskit_extend_haplotypes()`

- **作用：** 让祖先单倍型片段延伸穿过单分支区域。
- **处理方法：** 调用 tskit 的 `TreeSequence.extend_haplotypes`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_extend_haplotypes(tables, *, max_iter=10)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）（突变时间需已知）、`max_iter`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_extend_haplotypes

tables = json.load(open("tables.json"))
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_times

timed = tskit_compute_mutation_times(tables)["tables"]
r = tskit_extend_haplotypes(timed)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_trees"])
```

- **示例输出：**

```text
12 2
```

## 15) `dnakit.population.tskit_split_edges()`

- **作用：** 在跨越给定时间的每条边上插入一个单分支节点。
- **处理方法：** 调用 tskit 的 `TreeSequence.split_edges`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_split_edges(tables, *, time)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`time`。
- **输出：** 编辑后树序列的 `result_summary` 和输入行格式的 `tables`，以及 `provenance`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_split_edges

tables = json.load(open("tables.json"))
r = tskit_split_edges(tables, time=1.5)
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **示例输出：**

```text
16 18
```

## 16) `dnakit.population.tskit_map_mutations()`

- **作用：** 在一棵树上用简约法放置能解释所观察基因型的突变。
- **处理方法：** 调用 tskit 的 `Tree.map_mutations`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_map_mutations(tables, *, genotypes, tree_index=0, alleles=None, ancestral_state=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`genotypes`（每个样本一个等位基因序号）、`tree_index`、`alleles`、`ancestral_state`。
- **输出：** `ancestral_state`、`alleles` 和 `mutations`（node、derived_state、parent、time）。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_map_mutations

tables = json.load(open("tables.json"))
r = tskit_map_mutations(tables, genotypes=[1, 1, 0, 0], tree_index=0)
print(r["ancestral_state"], r["mutations"])
```

- **示例输出：**

```text
0 [{'node': 4, 'derived_state': '1', 'parent': -1, 'time': None}]
```
