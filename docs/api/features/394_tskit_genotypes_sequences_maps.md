# tskit 列表、基因型、序列、IBD 与重组图

这一组 `dnakit.population` 函数封装 [tskit](https://github.com/tskit-dev/tskit) 的列表与导出操作：突变、位点、基因型矩阵、变异、比对、单倍型、断点、VCF 个体映射、IBD 片段和重组图查询。

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

**示例文件：** `map.txt`

```text
Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
chr1	0	1.0	0.0
chr1	1000000	2.0	1.0
chr1	2000000	0.0	3.0
```

## 1) `dnakit.population.tskit_list_mutations()`

- **作用：** 列出整个树序列或某一棵树的突变。
- **处理方法：** 调用 tskit 的 `TreeSequence.mutations` / `Tree.mutations`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_list_mutations(tables, *, tree_index=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`tree_index`。
- **输出：** `num_mutations` 与 `mutations` 行。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_list_mutations

tables = json.load(open("tables.json"))
r = tskit_list_mutations(tables, tree_index=1)
print(r["num_mutations"])
for m in r["mutations"]:
    print(m)
```

- **示例输出：**

```text
2
{'id': 2, 'site': 2, 'node': 7, 'derived_state': 'A', 'parent': -1, 'time': None, 'position': 60.0}
{'id': 3, 'site': 3, 'node': 8, 'derived_state': 'C', 'parent': -1, 'time': None, 'position': 80.0}
```

## 2) `dnakit.population.tskit_list_sites()`

- **作用：** 列出整个树序列或某一棵树的变异位点。
- **处理方法：** 调用 tskit 的 `TreeSequence.sites` / `Tree.sites`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_list_sites(tables, *, tree_index=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`tree_index`。
- **输出：** `num_sites` 与 `sites` 行。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_list_sites

tables = json.load(open("tables.json"))
r = tskit_list_sites(tables)
for s in r["sites"]:
    print(s)
```

- **示例输出：**

```text
{'id': 0, 'position': 10.0, 'ancestral_state': 'A', 'num_mutations': 1}
{'id': 1, 'position': 30.0, 'ancestral_state': 'C', 'num_mutations': 1}
{'id': 2, 'position': 60.0, 'ancestral_state': 'G', 'num_mutations': 1}
{'id': 3, 'position': 80.0, 'ancestral_state': 'T', 'num_mutations': 1}
```

## 3) `dnakit.population.tskit_genotype_matrix()`

- **作用：** 返回等位基因序号的基因型矩阵（位点 × 样本）。
- **处理方法：** 调用 tskit 的 `TreeSequence.genotype_matrix`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_genotype_matrix(tables, *, samples=None, isolated_as_missing=True, alleles=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`samples`、`isolated_as_missing`、`alleles`。
- **输出：** `matrix`（-1 为缺失）、`positions`、`num_sites`、`num_samples`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genotype_matrix

tables = json.load(open("tables.json"))
r = tskit_genotype_matrix(tables)
print(r["positions"])
for row in r["matrix"]:
    print(row)
```

- **示例输出：**

```text
[10.0, 30.0, 60.0, 80.0]
[1, 1, 0, 0]
[0, 0, 1, 1]
[1, 0, 1, 0]
[0, 1, 0, 1]
```

## 4) `dnakit.population.tskit_variants()`

- **作用：** 逐位点的基因型调用，含等位基因状态与计数。
- **处理方法：** 调用 tskit 的 `TreeSequence.variants`, `Variant.states`, `Variant.counts`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_variants(tables, *, samples=None, isolated_as_missing=True, alleles=None, left=None, right=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`samples`、`isolated_as_missing`、`alleles`、`left`、`right`。
- **输出：** `variants` 行，含 `site`、`position`、`alleles`、`genotypes`、`states`、`counts`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_variants

tables = json.load(open("tables.json"))
r = tskit_variants(tables)
for v in r["variants"][:2]:
    print(v["position"], v["alleles"], v["genotypes"], v["states"], v["counts"])
```

- **示例输出：**

```text
10.0 ['A', 'T'] [1, 1, 0, 0] [np.str_('T'), np.str_('T'), np.str_('A'), np.str_('A')] {'A': 2, 'T': 2}
30.0 ['C', 'G'] [0, 0, 1, 1] [np.str_('C'), np.str_('C'), np.str_('G'), np.str_('G')] {'C': 2, 'G': 2}
```

## 5) `dnakit.population.tskit_alignments()`

- **作用：** 基于参考序列的样本全长比对序列。
- **处理方法：** 调用 tskit 的 `TreeSequence.alignments`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_alignments(tables, *, reference_sequence, missing_data_character='N', isolated_as_missing=True, samples=None, left=None, right=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`reference_sequence`（长度等于序列长度）、`missing_data_character`、`samples`、`left`、`right`。
- **输出：** 带 `node` 和 `sequence` 的 `alignments` 行。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_alignments

tables = json.load(open("tables.json"))
r = tskit_alignments(tables, reference_sequence="A" * 100)
for a in r["alignments"]:
    print(a["node"], a["sequence"][8:12], a["sequence"][58:62])
```

- **示例输出：**

```text
0 AATA AAAA
1 AATA AAGA
2 AAAA AAAA
3 AAAA AAGA
```

## 6) `dnakit.population.tskit_haplotypes()`

- **作用：** 每个样本仅在变异位点上的单倍型字符串。
- **处理方法：** 调用 tskit 的 `TreeSequence.haplotypes`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_haplotypes(tables, *, missing_data_character='-', isolated_as_missing=True, samples=None, left=None, right=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`missing_data_character`、`isolated_as_missing`、`samples`、`left`、`right`。
- **输出：** 带 `node` 和 `haplotype` 的 `haplotypes` 行。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_haplotypes

tables = json.load(open("tables.json"))
r = tskit_haplotypes(tables)
for h in r["haplotypes"]:
    print(h["node"], h["haplotype"])
```

- **示例输出：**

```text
0 TCAT
1 TCGC
2 AGAT
3 AGGC
```

## 7) `dnakit.population.tskit_breakpoints()`

- **作用：** 相邻谱系之间的基因组断点。
- **处理方法：** 调用 tskit 的 `TreeSequence.breakpoints`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_breakpoints(tables)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）。
- **输出：** `breakpoints`（从 0 到序列长度）和 `num_trees`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_breakpoints

tables = json.load(open("tables.json"))
r = tskit_breakpoints(tables)
print(r["breakpoints"], r["num_trees"])
```

- **示例输出：**

```text
[0.0, 50.0, 100.0] 2
```

## 8) `dnakit.population.tskit_map_to_vcf_model()`

- **作用：** 把样本节点映射到 VCF 个体。
- **处理方法：** 调用 tskit 的 `TreeSequence.map_to_vcf_model`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_map_to_vcf_model(tables, *, ploidy=None, contig_id='1')`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`ploidy`、`contig_id`。
- **输出：** `individual_names`、`individual_nodes`、`contig_id`、`contig_length`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_map_to_vcf_model

tables = json.load(open("tables.json"))
r = tskit_map_to_vcf_model(tables, ploidy=2)
print(r["individual_names"], r["individual_nodes"], r["contig_length"])
```

- **示例输出：**

```text
['tsk_0', 'tsk_1'] [[0, 1], [2, 3]] 100
```

## 9) `dnakit.population.tskit_sample_nodes_by_ploidy()`

- **作用：** 把样本节点按倍性分组成行。
- **处理方法：** 调用 tskit 的 `TreeSequence.sample_nodes_by_ploidy`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_sample_nodes_by_ploidy(tables, *, ploidy)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`ploidy`。
- **输出：** `groups`，每个个体一个节点 ID 列表。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_sample_nodes_by_ploidy

tables = json.load(open("tables.json"))
r = tskit_sample_nodes_by_ploidy(tables, ploidy=2)
print(r["groups"])
```

- **示例输出：**

```text
[[0, 1], [2, 3]]
```

## 10) `dnakit.population.tskit_ibd_segments()`

- **作用：** 由谱系得到样本对之间的同源（IBD）片段。
- **处理方法：** 调用 tskit 的 `TreeSequence.ibd_segments`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_ibd_segments(tables, *, within=None, between=None, max_time=None, min_span=None)`。
- **输入：** `tables`（显式表行或 `.trees` 路径）、`within` 或 `between`、`max_time`、`min_span`。
- **输出：** `segments` 行，含 `node_a`、`node_b`、`left`、`right`、`mrca`。
- **示例代码：**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_ibd_segments

tables = json.load(open("tables.json"))
r = tskit_ibd_segments(tables, within=[0, 1, 2])
print(r["num_segments"])
for s in r["segments"]:
    print(s)
```

- **示例输出：**

```text
6
{'node_a': 0, 'node_b': 1, 'left': 0.0, 'right': 50.0, 'mrca': 4}
{'node_a': 0, 'node_b': 1, 'left': 50.0, 'right': 100.0, 'mrca': 9}
{'node_a': 0, 'node_b': 2, 'left': 50.0, 'right': 100.0, 'mrca': 7}
{'node_a': 0, 'node_b': 2, 'left': 0.0, 'right': 50.0, 'mrca': 6}
{'node_a': 1, 'node_b': 2, 'left': 0.0, 'right': 50.0, 'mrca': 6}
{'node_a': 1, 'node_b': 2, 'left': 50.0, 'right': 100.0, 'mrca': 9}
```

## 11) `dnakit.population.tskit_read_hapmap()`

- **作用：** 读取 HapMap 格式的重组图文件。
- **处理方法：** 调用 tskit 的 `RateMap.read_hapmap`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_read_hapmap(source, *, sequence_length=None, has_header=True, position_col=None, rate_col=None, map_col=None)`。
- **输入：** `source` 文件、`sequence_length`、`has_header`、`position_col`、`rate_col`、`map_col`。
- **输出：** 带 `left`、`right`、`rate` 的 `intervals`，以及 `total_mass`、`mean_rate`。
- **示例代码：**

```python
from dnakit.population.tskit_p31_functions import tskit_read_hapmap

r = tskit_read_hapmap("map.txt")
print(r["sequence_length"], r["total_mass"], r["mean_rate"])
for i in r["intervals"]:
    print(i)
```

- **示例输出：**

```text
2000000.0 0.03 1.5e-08
{'left': 0.0, 'right': 1000000.0, 'rate': 1e-08}
{'left': 1000000.0, 'right': 2000000.0, 'rate': 1.9999999999999997e-08}
```

## 12) `dnakit.population.tskit_recombination_map_lookup()`

- **作用：** 在给定位置查询重组率或累计遗传图距。
- **处理方法：** 调用 tskit 的 `RateMap.get_rate` / `get_cumulative_mass`，结果转换为普通 Python 值。
- **API：** `dnakit.population.tskit_recombination_map_lookup(positions, *, quantity='rate', hapmap_file=None, map_positions=None, map_rates=None, sequence_length=None, has_header=True)`。
- **输入：** `positions`、`quantity`，以及以 `hapmap_file` 或 `map_positions` 加 `map_rates` 给出的图。
- **输出：** `values`，每个位置一个值。
- **示例代码：**

```python
from dnakit.population.tskit_p31_functions import tskit_recombination_map_lookup

r = tskit_recombination_map_lookup([500000, 1500000], hapmap_file="map.txt")
print(r["values"])
r = tskit_recombination_map_lookup([1000000], quantity="cumulative_mass", hapmap_file="map.txt")
print(r["values"])
```

- **示例输出：**

```text
[1e-08, 1.9999999999999997e-08]
[0.01]
```
