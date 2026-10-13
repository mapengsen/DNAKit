# 可选系统发育后端：dendropy、toytree 与 pyvolve

这一组函数是 `dnakit.phylogeny` 中基于可选依赖的适配器。[dendropy](https://github.com/jeetsukumaran/DendroPy) 读取 Newick 或 Nexus 树，报告 tip 与分支长度，计算 tip 间的系统发育距离，并比较两棵树的分支长度；[toytree](https://github.com/eaton-lab/toytree) 计算布局坐标（不绘图）、quartet 距离并按外群定根；[pyvolve](https://github.com/sjspielman/pyvolve) 沿给定的树用固定随机种子模拟 DNA 序列。

这些函数不推断系统发育树，不绘制图形，也不访问网络。树的拓扑和分支长度必须由调用者提供。依赖都是可选的，只有调用对应函数时才会导入；缺少包时会报错并给出安装命令。

## 示例数据

下面的示例使用这几个小文件。树只有 4 个 tip，只用于演示接口。

**示例文件：** `tree.nwk`

```text
((a:1,b:2):0.5,(c:0.25,d:3):1.5);
```

**示例文件：** `trees.nex`

```text
#NEXUS
begin taxa;
  dimensions ntax=3;
  taxlabels X Y Z;
end;
begin trees;
  tree one = [&R] ((X:1,Y:2):3,Z:4);
  tree two = ((X:1,Z:2):3,Y:4);
end;
```

**示例文件：** `tree_b.nwk`

```text
((a:1,b:1):0.5,(c:0.25,d:2):1.5);
```

**示例文件：** `tree_c.nwk`

```text
((a:1,c:1):1,(b:1,d:1):1);
```

## 1) `dnakit.phylogeny.read_tree_summary()`

- **作用：** 读取 Newick 或 Nexus 文件中的一棵树，报告 tip 与分支长度。
- **处理方法：** 用 dendropy 的 `TreeList.get` 解析文件，`tree_index` 选择第几棵树。tip 按标签排序并附带各自的分支长度；`branch_lengths` 按前序列出除根边以外的每条边。标签中的下划线保持原样。缺失或非有限的长度记为 `None`。
- **API：** `dnakit.phylogeny.read_tree_summary(path, *, schema="newick", tree_index=0, max_tips=10000, max_file_bytes=10000000)`。
- **输入：** 已存在的树文件路径；`schema` 取 `"newick"` 或 `"nexus"`；`max_tips` 和 `max_file_bytes` 是显式上限，超过时报错。
- **输出：** 含 `tree_count`、`rooting`、`tips`、`branch_lengths`、`missing_branch_lengths` 和 `total_branch_length` 的字典。
- **示例代码：**

```python
from dnakit.phylogeny import read_tree_summary

result = read_tree_summary("tree.nwk")
print(result["result"]["tips"])
print(result["result"]["total_branch_length"])
nexus = read_tree_summary("trees.nex", schema="nexus", tree_index=1)
print(nexus["result"]["tree_count"], nexus["result"]["rooting"])
```

- **示例输出：**

```text
[{'label': 'a', 'branch_length': 1.0}, {'label': 'b', 'branch_length': 2.0}, {'label': 'c', 'branch_length': 0.25}, {'label': 'd', 'branch_length': 3.0}]
8.25
2 undefined
```

## 2) `dnakit.phylogeny.calculate_patristic_distances()`

- **作用：** 计算一棵树上所有 tip 两两之间的系统发育距离矩阵（patristic distance）。
- **处理方法：** 用 dendropy 的 `phylogenetic_distance_matrix` 计算每对 tip 之间路径上的分支长度之和。根边不在任何 tip 间的路径上。要求每条非根边都有分支长度，缺失时报错，而不是当作 0。
- **API：** `dnakit.phylogeny.calculate_patristic_distances(path, *, schema="newick", tree_index=0, max_taxa=200, max_file_bytes=10000000)`。
- **输入：** 树文件路径；`max_taxa` 限制 tip 数量。
- **输出：** 含 `taxa`（按标签排序）和 `matrix`（对称矩阵，对角线为 0）的字典。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_patristic_distances

result = calculate_patristic_distances("tree.nwk")
print(result["result"]["taxa"])
print(result["result"]["matrix"][0])
```

- **示例输出：**

```text
['a', 'b', 'c', 'd']
[0.0, 3.0, 3.25, 6.0]
```

## 3) `dnakit.phylogeny.calculate_branch_length_distance()`

- **作用：** 比较两棵 tip 标签相同的树的分支长度。
- **处理方法：** 用 dendropy 的 `euclidean_distance`，计算 Felsenstein (2004) 分支长度距离，即按无根二分（split）比较各分支长度差的平方和再开方。两棵树的根边会合并到无根的分支上，所以根的位置不影响结果。
- **API：** `dnakit.phylogeny.calculate_branch_length_distance(path_a, path_b, *, schema="newick", tree_index=0, max_taxa=10000, max_file_bytes=10000000)`。
- **输入：** 两个树文件路径，两棵树的 tip 标签必须完全相同；`tree_index` 同时作用于两个文件。
- **输出：** 含 `distance` 和 `tip_count` 的字典。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_branch_length_distance

result = calculate_branch_length_distance("tree.nwk", "tree_b.nwk")
print(result["result"]["distance"])
```

- **示例输出：**

```text
1.4142135623730951
```

## 4) `dnakit.phylogeny.calculate_node_coordinates()`

- **作用：** 计算一棵树在指定布局下每个节点的坐标，供外部绘图使用。
- **处理方法：** 用 toytree 的 `LinearLayout`、`CircularLayout` 或 `UnrootedLayout` 计算坐标。`r`、`l`、`u`、`d` 是直线布局，`c` 是圆形布局，`unrooted` 是无根布局。坐标在 toytree 的布局单位中，函数不绘制任何图形。
- **API：** `dnakit.phylogeny.calculate_node_coordinates(path, *, layout="r", use_edge_lengths=True, max_tips=10000, max_file_bytes=10000000)`。
- **输入：** Newick 树文件路径；`layout` 取 `"r"`、`"l"`、`"u"`、`"d"`、`"c"` 或 `"unrooted"`；`use_edge_lengths=False` 时所有边按单位长度布局。
- **输出：** 含 `layout`、`tip_count` 和 `nodes` 的字典。每个节点包含 `index`、`name`（内部节点为 `None`）、`is_tip`、`parent`、`x` 和 `y`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_node_coordinates

result = calculate_node_coordinates("tree.nwk", layout="r")
for node in result["result"]["nodes"][:4]:
    print(node["name"], node["x"], node["y"])
```

- **示例输出：**

```text
a -3.0 0.0
b -2.0 1.0
c -2.75 2.0
d 0.0 3.0
```

## 5) `dnakit.phylogeny.calculate_quartet_distance()`

- **作用：** 统计两棵 tip 标签相同的树在 4 个 tip 子集（quartet）上的一致与不一致数量。
- **处理方法：** 用 toytree 的 `distance.get_treedist_quartets` 计算。对每个 4 tip 子集，若两棵树诱导出的拓扑相同则计为 `shared`，不同则计为 `different`。`total_quartets` 等于 `C(n, 4)`。
- **API：** `dnakit.phylogeny.calculate_quartet_distance(path_a, path_b, *, max_taxa=40, max_file_bytes=10000000)`。
- **输入：** 两个 Newick 树文件路径，tip 标签必须相同，至少 4 个 tip。
- **输出：** 含 `total_quartets`、`shared_quartets`、`different_quartets` 和 `different_fraction` 的字典。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_quartet_distance

result = calculate_quartet_distance("tree.nwk", "tree_c.nwk")
print(result["result"]["total_quartets"], result["result"]["different_quartets"])
```

- **示例输出：**

```text
1 1
```

## 6) `dnakit.phylogeny.root_tree_on_outgroup()`

- **作用：** 按外群为无根或有根的 Newick 树定根，并写出结果文件。
- **处理方法：** 用 toytree 的 `root` 把新根放在外群 tip 的最近共同祖先（MRCA）上方的边上，根边的长度平分给两侧。外群必须是单系的，否则报错。输出文件必须是新文件，不覆盖已有文件。
- **API：** `dnakit.phylogeny.root_tree_on_outgroup(path, outgroup, output_path, *, max_tips=10000, max_file_bytes=10000000)`。
- **输入：** 树文件路径；`outgroup` 是一个 tip 标签或一组标签；`output_path` 是尚不存在的输出文件路径。
- **输出：** 含 `output_path`、`outgroup`、`tip_count` 和 `newick` 的字典。
- **示例代码：**

```python
from pathlib import Path
from dnakit.phylogeny import root_tree_on_outgroup

Path("rooted.nwk").unlink(missing_ok=True)
result = root_tree_on_outgroup("tree.nwk", "c", "rooted.nwk")
print(result["result"]["newick"].strip())
```

- **示例输出：**

```text
(c:0.125,(d:3,(a:1,b:2):2):0.125);
```

## 7) `dnakit.phylogeny.simulate_nucleotide_alignment()`

- **作用：** 沿给定的 Newick 树模拟 DNA 比对，并写成 FASTA 文件。
- **处理方法：** 用 pyvolve 的 `nucleotide` 模型。转换（A↔G、C↔T）的速率乘以 `kappa`，其余交换率为 1，即 HKY85 速率矩阵；碱基频率由 `state_freqs` 给定。分支长度的单位是每位点的期望替换数。`seed` 固定随机数，同样的输入得到同样的序列。
- **API：** `dnakit.phylogeny.simulate_nucleotide_alignment(tree_path, output_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`。
- **输入：** Newick 树文件路径，tip 名必须唯一且每条边都有长度；`output_path` 必须是尚不存在的文件；`length` 是位点数；`state_freqs` 为 A、C、G、T 的频率，和为 1。
- **输出：** 含 `output_path`、`taxa`、`sequence_count`、`site_count`、`model` 和 `seed` 的字典。
- **示例代码：**

```python
from pathlib import Path
from dnakit.phylogeny import simulate_nucleotide_alignment

Path("alignment.fasta").unlink(missing_ok=True)
result = simulate_nucleotide_alignment("tree.nwk", "alignment.fasta", length=12, seed=1)
print(result["result"]["taxa"], result["result"]["site_count"])
print(Path("alignment.fasta").read_text().splitlines()[:2])
```

- **示例输出：**

```text
['a', 'b', 'c', 'd'] 12
['>a', 'TTCCATCAGGTT']
```

## 8) `dnakit.phylogeny.simulate_nucleotide_sequences()`

- **作用：** 与上一个函数相同的模拟，但直接返回 tip 序列，不写文件。
- **处理方法：** 使用相同的 HKY85 速率矩阵和随机种子规则，从 pyvolve 的 `get_sequences` 取得 tip 序列。`tips × length` 超过 `max_bases` 时报错。
- **API：** `dnakit.phylogeny.simulate_nucleotide_sequences(tree_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`。
- **输入：** 与 `simulate_nucleotide_alignment()` 相同，但不需要输出路径。
- **输出：** 含 `taxa`、`site_count`、`sequences`（标签到序列的字典）、`model` 和 `seed` 的字典。
- **示例代码：**

```python
from dnakit.phylogeny import simulate_nucleotide_sequences

result = simulate_nucleotide_sequences("tree.nwk", length=12, kappa=2.0, seed=1)
print(result["result"]["sequences"]["a"])
print(result["result"]["sequences"]["d"])
```

- **示例输出：**

```text
TTCCCTCAGGTT
TCGAGTTTGCTA
```

## 9) `dnakit.phylogeny.simulate_substitution_counts()`

- **作用：** 用固定种子模拟序列，并统计每条分支上发生的替换次数。
- **处理方法：** 使用 pyvolve 的精确 Gillespie 算法（`algorithm=1`），而不是矩阵指数法，后者会低估替换数。分支名是 tip 名，内部分支为 pyvolve 生成的 `internalNode` 名；根边不计数。
- **API：** `dnakit.phylogeny.simulate_substitution_counts(tree_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`。
- **输入：** 与 `simulate_nucleotide_alignment()` 相同。
- **输出：** 含 `taxa`、`site_count`、`branch_substitution_counts`、`total_substitutions`、`model` 和 `seed` 的字典。
- **示例代码：**

```python
from dnakit.phylogeny import simulate_substitution_counts

result = simulate_substitution_counts("tree.nwk", length=1000, seed=1)
print(result["result"]["branch_substitution_counts"])
print(result["result"]["total_substitutions"])
```

- **示例输出：**

```text
{'a': 1038, 'b': 2025, 'c': 263, 'd': 2968, 'internalNode1': 498, 'internalNode2': 1538}
8330
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-dendropy]"   # dendropy，read_tree_summary 等树读取与距离函数
python -m pip install "dnakit[bio-toytree]"    # toytree，布局坐标、quartet 距离与定根
python -m pip install "dnakit[bio-pyvolve]"    # pyvolve，DNA 序列模拟
```

各包的项目地址：[dendropy](https://github.com/jeetsukumaran/DendroPy)、[toytree](https://github.com/eaton-lab/toytree)、[pyvolve](https://github.com/sjspielman/pyvolve)。

## 验证方式

测试文件是 `others/tests/unit/test_phylogeny_backends.py`，共 32 项，主要检查：

- dendropy 的 tip、分支长度、Nexus 树选择、patristic 矩阵与根边恒等式（两点路径长度 = 各自到根的距离之和减去两倍 MRCA 深度）一致；分支长度距离与直接调用 dendropy 的 `euclidean_distance` 一致，并用手工算出的 sqrt(2) 和根位置不变的 0 值校验。
- toytree 的直线布局坐标与手工计算的数值一致，并与直接调用 `LinearLayout` 一致；quartet 计数与 5 个 tip 的手工计数（共享 3、不同 2）一致；定根后的 Newick 与手工结果一致，且总分支长度保持不变。
- pyvolve 的 Jukes–Cantor 差异比例与解析公式一致；转换比例与 kappa/(kappa+2) 一致；每条分支的替换次数与期望值 `长度 × 位点数` 的偏差在随机误差范围内；碱基组成与给定频率一致；同一种子的输出可复现，文件与内存结果相同。
- 缺少任一可选依赖时给出安装提示；缺失分支长度、重复 tip 名、非法 kappa、频率和种子，以及输出文件已存在等情况都会被拒绝。

## 当前范围

这一批只包装树读取、距离、布局坐标、定根和 DNA 序列模拟。以下内容没有包装：

- dendropy 的树模拟（出生-死亡、coalescent）和字符模拟，因为超出了这一批的 DNA 序列范围。
- dendropy 与 toytree 的 Robinson–Foulds 距离。它们和 `dnakit.bioinformatics.phylogeny.robinson_foulds()` 计算同一个量，因此不重复包装。
- toytree 的绘图函数，以及 toytree 的树推断。
- pyvolve 的密码子、氨基酸模型，自定义速率矩阵，以及 PAML 或突变选择模型。
- pyvolve 的其他序列格式（Phylip、Nexus 等）。当前只写 FASTA。
