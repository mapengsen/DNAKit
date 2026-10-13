# 102 审阅补充（G18p）：toytree、dendropy 与 pyranges1

本页记录审阅组 G18p 新增的 42 个封装：toytree 的树节点查询、枚举、编辑、一致性支持与网络解析（26 个，[toytree](https://github.com/eaton-lab/toytree)）；dendropy 的分支划分比较与树形向量（6 个，[DendroPy](https://github.com/jeetsukumaran/DendroPy)）；pyranges1 的读取器、区间合并、ORF 与重叠统计（10 个，[pyranges1](https://pypi.org/project/pyranges1/)）。计算在本地完成，不调用网络或外部程序，也不绘图。树输入为 Newick 文本或单树文件；pyranges1 输入为本地文件。可选依赖：`bio-toytree`、`bio-dendropy`、`bio-pyranges1`，bigWig 与 ORF 另需 [pyBigWig](https://github.com/deeptools/pyBigWig) （`bio-pybigwig`）与 [pyfaidx](https://github.com/mdshw5/pyfaidx)（`bio-pyfaidx`）。模块为 `dnakit.phylogeny`、`dnakit.phylogeny` 与 `dnakit.annotation`，接入包导出由总负责人完成。

## 1) `dnakit.phylogeny.list_tree_node_path()`

- **作用：** 列出树上两个节点之间的路径，包括两端节点。
- **处理方法：** 调用 toytree.distance.get_node_path 查找路径，并给出每个节点的 idx、name 和是否为叶。
- **API：** `dnakit.phylogeny.list_tree_node_path(tree, start, end, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 文本或单树 Newick 文件；start 与 end 为节点名或整数节点索引。
- **输出：** nodes（按路径顺序排列的节点）、edge_count（边数）和 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import list_tree_node_path

result = list_tree_node_path("((A:1,B:2):1,(C:3,D:1):2);", "A", "C")["result"]
print([node["idx"] for node in result["nodes"]], result["edge_count"])
```

- **示例输出：**

```text
[0, 4, 6, 5, 2] 4
```

## 2) `dnakit.phylogeny.calculate_tree_node_distance()`

- **作用：** 计算树上两个节点之间的分支长度距离，或只按拓扑边数计算。
- **处理方法：** 调用 toytree.distance.get_node_distance；topology_only=True 时每条边计为 1。
- **API：** `dnakit.phylogeny.calculate_tree_node_distance(tree, start, end, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；start 与 end 为节点名或整数索引；topology_only 为布尔值。
- **输出：** distance（数值）、topology_only 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tree_node_distance

tree = "((A:1,B:2):1,(C:3,D:1):2);"
print(calculate_tree_node_distance(tree, "A", "C")["result"]["distance"])
print(calculate_tree_node_distance(tree, "A", "C", topology_only=True)["result"]["distance"])
```

- **示例输出：**

```text
7.0
4.0
```

## 3) `dnakit.phylogeny.calculate_all_node_distance_matrix()`

- **作用：** 计算树上所有节点（叶与内部节点）两两之间的距离矩阵。
- **处理方法：** 调用 toytree.distance.get_node_distance_matrix；行列按节点索引排列，未命名内部节点记为 node_<idx>。
- **API：** `dnakit.phylogeny.calculate_all_node_distance_matrix(tree, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；topology_only 为布尔值。
- **输出：** labels（节点标签）、matrix（对称矩阵）和 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_all_node_distance_matrix

result = calculate_all_node_distance_matrix("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print(result["labels"])
print(result["matrix"][0])
```

- **示例输出：**

```text
['A', 'B', 'C', 'D', 'node_4', 'node_5', 'node_6']
[0.0, 3.0, 7.0, 5.0, 1.0, 4.0, 2.0]
```

## 4) `dnakit.phylogeny.calculate_internal_node_distance_matrix()`

- **作用：** 计算内部节点（包括根）两两之间的距离矩阵。
- **处理方法：** 调用 toytree.distance.get_internal_node_distance_matrix，按节点索引排列。
- **API：** `dnakit.phylogeny.calculate_internal_node_distance_matrix(tree, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；topology_only 为布尔值。
- **输出：** labels、matrix 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_internal_node_distance_matrix

result = calculate_internal_node_distance_matrix("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print(result["labels"])
print(result["matrix"])
```

- **示例输出：**

```text
['node_4', 'node_5', 'node_6']
[[0.0, 3.0, 1.0], [3.0, 0.0, 2.0], [1.0, 2.0, 0.0]]
```

## 5) `dnakit.phylogeny.calculate_descendant_distances()`

- **作用：** 给出某节点（默认为根）到它每个后代节点的距离。
- **处理方法：** 调用 toytree.distance.get_descendant_dists；结果按节点索引排序。
- **API：** `dnakit.phylogeny.calculate_descendant_distances(tree, node=None, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；node 为节点名、整数索引或 None（根）。
- **输出：** node、descendants（每个后代的 idx、name、is_tip 与 distance）和 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_descendant_distances

result = calculate_descendant_distances("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print([(row["idx"], row["distance"]) for row in result["descendants"]])
```

- **示例输出：**

```text
[(0, 2.0), (1, 3.0), (2, 5.0), (3, 3.0), (4, 1.0), (5, 2.0), (6, 0.0)]
```

## 6) `dnakit.phylogeny.find_tree_farthest_node()`

- **作用：** 找出距指定节点最远的节点及其距离。
- **处理方法：** 调用 toytree.distance.get_farthest_node 和 get_farthest_node_distance。
- **API：** `dnakit.phylogeny.find_tree_farthest_node(tree, node, *, topology_only=False, descendants_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；node 为节点名或整数索引；topology_only 与 descendants_only 为布尔值。
- **输出：** node、farthest（最远节点）、distance 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import find_tree_farthest_node

result = find_tree_farthest_node("((A:1,B:2):1,(C:3,D:1):2);", "A")["result"]
print(result["farthest"]["name"], result["distance"])
```

- **示例输出：**

```text
C 7.0
```

## 7) `dnakit.phylogeny.calculate_farthest_node_distance()`

- **作用：** 计算指定节点（默认为根）到最远节点的距离。
- **处理方法：** 调用 toytree.distance.get_farthest_node_distance。
- **API：** `dnakit.phylogeny.calculate_farthest_node_distance(tree, node=None, *, topology_only=False, descendants_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树；node 为节点名、整数索引或 None；布尔选项同上。
- **输出：** node、distance、topology_only、descendants_only 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_farthest_node_distance

print(calculate_farthest_node_distance("((A:1,B:2):1,(C:3,D:1):2);", "A")["result"]["distance"])
```

- **示例输出：**

```text
7.0
```

## 8) `dnakit.phylogeny.list_tree_quadripartitions()`

- **作用：** 列出每条内部边两侧的四组叶集合（四分划分）。
- **处理方法：** 调用 toytree.enum.iter_quadripartitions，把每个四分划分展开为四组排序后的叶名列表。
- **API：** `dnakit.phylogeny.list_tree_quadripartitions(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树。
- **输出：** quadripartitions（每项为四组叶名）、count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import list_tree_quadripartitions

result = list_tree_quadripartitions("((A,B),(C,D),E);")["result"]
print(result["count"])
for groups in result["quadripartitions"]:
    print(groups)
```

- **示例输出：**

```text
2
[['A'], ['B'], ['C', 'D'], ['E']]
[['C'], ['D'], ['A', 'B'], ['E']]
```

## 9) `dnakit.phylogeny.list_tree_quartets()`

- **作用：** 列出由内部边诱导的四元组（边两侧各取两个叶的组合）。
- **处理方法：** 调用 toytree.enum.iter_quartets，并按 left、right 排序。
- **API：** `dnakit.phylogeny.list_tree_quartets(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树。
- **输出：** quartets（left 与 right 两组叶名）、count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import list_tree_quartets

result = list_tree_quartets("((A,B),(C,D),E);")["result"]
print(result["count"], result["quartets"][0])
```

- **示例输出：**

```text
5 {'left': ['A', 'B'], 'right': ['C', 'D']}
```

## 10) `dnakit.phylogeny.map_consensus_tree_features()`

- **作用：** 把一组树中的分支支持率（split support）和分支统计映射到目标树上。
- **处理方法：** 用 toytree.mtree 组成树列表，再调用 toytree.infer.get_consensus_features；support 为包含该分支的树所占比例。
- **API：** `dnakit.phylogeny.map_consensus_tree_features(tree, trees, *, features=('support',), max_trees=1000, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 目标 Newick 树、Newick 树列表（叶集须与目标树相同）和 features 名称列表。
- **输出：** newick、features、tree_count、nodes（每个节点的 tips 与 values）和 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import map_consensus_tree_features

trees = ["((A,B),(C,D));", "((A,B),(C,D));", "((A,C),(B,D));"]
result = map_consensus_tree_features("((A,B),(C,D));", trees)["result"]
print(result["newick"])
node = next(row for row in result["nodes"] if row["tips"] == ["C", "D"])
print(round(node["values"]["support"], 4))
```

- **示例输出：**

```text
(A:1,B:1,(C:1,D:1)0.666666666667:1);
0.6667
```

## 11) `dnakit.phylogeny.extend_tree_tips_to_align()`

- **作用：** 把每个叶的分支延长，使所有叶端点位于同一深度。
- **处理方法：** 调用 toytree.mod.edges_extend_tips_to_align，返回编辑后的树（Newick）。
- **API：** `dnakit.phylogeny.extend_tree_tips_to_align(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树。
- **输出：** newick、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import extend_tree_tips_to_align

print(extend_tree_tips_to_align("((A:1,B:2):1,(C:3,D:1):2);")["result"]["newick"])
```

- **示例输出：**

```text
((A:4,B:4):1,(C:3,D:3):2);
```

## 12) `dnakit.phylogeny.set_tree_node_heights()`

- **作用：** 把指定节点的高度（到根的距离）设为给定值。
- **处理方法：** 调用 toytree.mod.edges_set_node_heights；键为节点名或整数索引。
- **API：** `dnakit.phylogeny.set_tree_node_heights(tree, heights, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树和 {节点: 高度} 映射（高度为非负数）。
- **输出：** newick、heights（节点与高度列表）与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import set_tree_node_heights

print(set_tree_node_heights("((A:1,B:2):1,(C:3,D:1):2);", {"A": 2.0})["result"]["newick"])
```

- **示例输出：**

```text
((A:2,B:2):1,(C:3,D:1):2);
```

## 13) `dnakit.phylogeny.rotate_tree_node()`

- **作用：** 反转一个节点的子节点顺序（子节点左右互换）。
- **处理方法：** 调用 toytree.mod.rotate_node，返回编辑后的树（Newick）。
- **API：** `dnakit.phylogeny.rotate_tree_node(tree, node, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树和节点名或整数索引。
- **输出：** newick、node、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import rotate_tree_node

print(rotate_tree_node("((A:1,B:2):1,(C:3,D:1):2);", 4)["result"]["newick"])
```

- **示例输出：**

```text
((B:2,A:1):1,(C:3,D:1):2);
```

## 14) `dnakit.phylogeny.drop_tree_tips()`

- **作用：** 按叶名删除指定叶，并保留其余部分。
- **处理方法：** 调用 toytree.mod.drop_tips；删除后仅剩一个叶时报错。
- **API：** `dnakit.phylogeny.drop_tree_tips(tree, tips, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树和要删除的叶名列表（须存在于树中）。
- **输出：** newick、dropped（排序后的叶名）、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import drop_tree_tips

print(drop_tree_tips("((A:1,B:2):1,(C:3,D:1):2);", ["A"])["result"]["newick"])
```

- **示例输出：**

```text
((C:3,D:1):2,B:3);
```

## 15) `dnakit.phylogeny.bisect_tree()`

- **作用：** 在某节点上方的边处把树切成两部分：该节点一侧的子树和其余部分。
- **处理方法：** 调用 toytree.mod.bisect（reroot=False）；dist_partition 为正数。
- **API：** `dnakit.phylogeny.bisect_tree(tree, node, *, dist_partition=1.0, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点名或整数索引与 dist_partition。
- **输出：** node、subtree_newick 与 subtree_tips、remainder_newick 与 remainder_tips、provenance。
- **示例代码：**

```python
from dnakit.phylogeny import bisect_tree

result = bisect_tree("((A:1,B:2):1,(C:3,D:1):2);", 4)["result"]
print(result["subtree_newick"], result["remainder_newick"])
```

- **示例输出：**

```text
(A:1,B:2):1; ((C:3,D:1):2);
```

## 16) `dnakit.phylogeny.resolve_tree_polytomies()`

- **作用：** 以固定种子随机拆分多分支节点，得到二叉树。
- **处理方法：** 调用 toytree.mod.resolve_polytomies，并使用 numpy 的 default_rng(seed)；不改变全局随机状态。
- **API：** `dnakit.phylogeny.resolve_tree_polytomies(tree, *, node=None, seed, dist=1e-06, recursive=True, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、可选 node、必填整数 seed（0 到 2**32-1）、dist 与 recursive。
- **输出：** newick、node、seed、dist、recursive、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import resolve_tree_polytomies

print(resolve_tree_polytomies("(A,B,C,D);", seed=3)["result"]["newick"])
```

- **示例输出：**

```text
(D:1,(A:1,(B:1,C:1):1e-06):1e-06);
```

## 17) `dnakit.phylogeny.add_tree_internal_node()`

- **作用：** 在某节点上方的边中插入一个内部节点。
- **处理方法：** 调用 toytree.mod.add_internal_node；name 与 dist 可选。
- **API：** `dnakit.phylogeny.add_tree_internal_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点名或整数索引、可选 name 与 dist（非负）。
- **输出：** newick、node、name、dist、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import add_tree_internal_node

result = add_tree_internal_node("((A:1,B:2):1,(C:3,D:1):2);", "A", name="X", dist=0.5)["result"]
print(result["newick"], result["name"])
```

- **示例输出：**

```text
((B:2,(A:0.5):0.5):1,(C:3,D:1):2); X
```

## 18) `dnakit.phylogeny.add_tree_child_node()`

- **作用：** 为某节点添加一个新的叶子子节点。
- **处理方法：** 调用 toytree.mod.add_child_node。若目标是叶，该叶会变成只有新子节点的无名内部节点。
- **API：** `dnakit.phylogeny.add_tree_child_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点名或整数索引、可选 name 与 dist。
- **输出：** newick、node、name、dist、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import add_tree_child_node

print(add_tree_child_node("((A:1,B:2):1,(C:3,D:1):2);", 4, name="Z", dist=0.5)["result"]["newick"])
```

- **示例输出：**

```text
((A:1,B:2,Z:0.5):1,(C:3,D:1):2);
```

## 19) `dnakit.phylogeny.add_tree_sister_node()`

- **作用：** 在某节点旁边添加一个新的姐妹叶。
- **处理方法：** 调用 toytree.mod.add_sister_node。
- **API：** `dnakit.phylogeny.add_tree_sister_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点名或整数索引、可选 name 与 dist。
- **输出：** newick、node、name、dist、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import add_tree_sister_node

print(add_tree_sister_node("((A:1,B:2):1,(C:3,D:1):2);", "A", name="Z", dist=0.5)["result"]["newick"])
```

- **示例输出：**

```text
((A:1,B:2,Z:0.5):1,(C:3,D:1):2);
```

## 20) `dnakit.phylogeny.add_tree_internal_node_and_child()`

- **作用：** 在某节点上方的边中插入一个内部节点，并给它添加一个新叶子。
- **处理方法：** 调用 toytree.mod.add_internal_node_and_child；可选择命名父节点并设置父边长度。
- **API：** `dnakit.phylogeny.add_tree_internal_node_and_child(tree, node, *, name=None, dist=None, parent_name=None, parent_dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点名或整数索引、可选 name、dist、parent_name、parent_dist。
- **输出：** newick、node、name、dist、parent_name、parent_dist、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import add_tree_internal_node_and_child

result = add_tree_internal_node_and_child(
    "((A:1,B:2):1,(C:3,D:1):2);", "A", name="Z", dist=0.5, parent_name="P", parent_dist=0.5
)["result"]
print(result["newick"], result["tips"])
```

- **示例输出：**

```text
((B:2,(A:0.5,Z:0.5):0.5):1,(C:3,D:1):2); ['A', 'B', 'C', 'D', 'Z']
```

## 21) `dnakit.phylogeny.add_tree_internal_node_and_subtree()`

- **作用：** 把一棵 Newick 子树嫁接到某节点上方的边上，通过新的内部节点连接。
- **处理方法：** 调用 toytree.mod.add_internal_node_and_subtree；子树叶名不得与原树重复。
- **API：** `dnakit.phylogeny.add_tree_internal_node_and_subtree(tree, node, subtree, *, subtree_stem_dist=None, subtree_rescale=False, parent_name=None, parent_dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点、Newick 子树、可选 subtree_stem_dist、subtree_rescale、parent_name、parent_dist。
- **输出：** newick、node、subtree_tips、subtree_stem_dist、subtree_rescale、parent_name、parent_dist、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import add_tree_internal_node_and_subtree

result = add_tree_internal_node_and_subtree(
    "((A:1,B:2):1,(C:3,D:1):2);", "A", "(X:1,Y:1);", subtree_stem_dist=0.5, parent_name="P", parent_dist=0.5
)["result"]
print(result["newick"], result["subtree_tips"])
```

- **示例输出：**

```text
((B:2,(A:0.5,(X:1,Y:1):0.5):0.5):1,(C:3,D:1):2); ['X', 'Y']
```

## 22) `dnakit.phylogeny.remove_tree_nodes()`

- **作用：** 删除指定的一组节点，并保留它们的子节点。
- **处理方法：** 调用 toytree.mod.remove_nodes；preserve_dists=True 时把被删节点的分支长度加到子边上。删除叶后至少须保留两个叶。
- **API：** `dnakit.phylogeny.remove_tree_nodes(tree, nodes, *, preserve_dists=True, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树、节点列表（名称或索引）与 preserve_dists。
- **输出：** newick、removed、preserve_dists、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import remove_tree_nodes

print(remove_tree_nodes("((A:1,B:2):1,(C:3,D:1):2);", [4])["result"]["newick"])
```

- **示例输出：**

```text
((C:3,D:1):2,A:2,B:3);
```

## 23) `dnakit.phylogeny.merge_tree_nodes()`

- **作用：** 合并叶属性相同的分支：若某clade的所有叶共享同一节点属性值，则只保留一个代表。
- **处理方法：** 调用 toytree.mod.merge_nodes，merge_method 为节点属性名（须为标识符）。
- **API：** `dnakit.phylogeny.merge_tree_nodes(tree, feature, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树和一个标识符形式的节点属性名，例如 dist。
- **输出：** newick、feature、removed_node_count、tip_count、tips 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import merge_tree_nodes

result = merge_tree_nodes("((A:1,B:1):1,(C:1,D:1):1);", "dist")["result"]
print(result["newick"], result["removed_node_count"])
```

- **示例输出：**

```text
(A:2,C:2); 4
```

## 24) `dnakit.phylogeny.parse_tree_network()`

- **作用：** 解析带杂交节点（#H）的扩展 Newick 网络，得到主树和渐渗（admixture）边。
- **处理方法：** 调用 toytree.network.parse_network_to_tree_and_admix；输入可为网络文本或文件。
- **API：** `dnakit.phylogeny.parse_tree_network(network, *, max_file_bytes=10000000, verbose=True)`。
- **输入：** 网络 Newick 文本（含 ;）或文件路径。
- **输出：** newick、tips、admixture_edges（destination_tips、source_tips、gamma）、admixture_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import parse_tree_network

result = parse_tree_network("((A:1,(B:1,#H1:0::0.3):1):1,(#H1:0::0.7,C:2):1);")["result"]
print(result["tips"])
print(result["admixture_edges"])
```

- **示例输出：**

```text
['A', 'B', 'C']
[{'destination_tips': [], 'source_tips': ['B'], 'gamma': 0.3}]
```

## 25) `dnakit.phylogeny.calculate_tip_distance_matrix_from_vcv()`

- **作用：** 由叶的方差-协方差矩阵计算叶两两之间的距离矩阵。
- **处理方法：** 调用 toytree.pcm.get_distance_matrix_from_vcv_matrix；每个元素为 V[i,i]+V[j,j]-2V[i,j]，即平方欧氏距离（上游不开平方）。
- **API：** `dnakit.phylogeny.calculate_tip_distance_matrix_from_vcv(vcv, labels, *, verbose=True)`。
- **输入：** 对称方阵 vcv 与同长度的唯一叶名列表 labels。
- **输出：** labels、matrix 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tip_distance_matrix_from_vcv

result = calculate_tip_distance_matrix_from_vcv([[2, 1], [1, 2]], ["x", "y"])["result"]
print(result["labels"], result["matrix"])
```

- **示例输出：**

```text
['x', 'y'] [[0.0, 2.0], [2.0, 0.0]]
```

## 26) `dnakit.phylogeny.calculate_tip_diversification_rates()`

- **作用：** 计算每个叶的分化速率（Jetz 等 2012 的 tip-level 分化速率）。
- **处理方法：** 调用 toytree.pcm.get_tip_level_diversification；只保留叶的 DR 值，去掉 toytree 返回的无名根行。
- **API：** `dnakit.phylogeny.calculate_tip_diversification_rates(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树。
- **输出：** rates（每项为 tip 与 diversification_rate）、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tip_diversification_rates

for row in calculate_tip_diversification_rates("((A:1,B:2):1,(C:3,D:1):2);")["result"]["rates"]:
    print(row["tip"], row["diversification_rate"])
```

- **示例输出：**

```text
A 2.0
B 2.0
C 1.0
D 1.0
```

## 27) `dnakit.phylogeny.list_missing_bipartitions()`

- **作用：** 列出参考树中有、而比较树中没有的分支划分（split）。
- **处理方法：** 调用 dendropy.calculate.treecompare.find_missing_bipartitions；每个划分取两侧中较小的一侧（并列时取不含首个叶名的一侧）。
- **API：** `dnakit.phylogeny.list_missing_bipartitions(reference_tree, comparison_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵 Newick 树，叶名集合必须相同。
- **输出：** bipartitions（叶名列表）、missing_count、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import list_missing_bipartitions

result = list_missing_bipartitions("((A,B),(C,D),E);", "((A,C),(B,D),E);")["result"]
print(result["bipartitions"], result["missing_count"])
```

- **示例输出：**

```text
[['A', 'B'], ['C', 'D']] 2
```

## 28) `dnakit.phylogeny.calculate_n_bar_index()`

- **作用：** 计算树形统计量 N-bar：每个叶之上的平均节点数。
- **处理方法：** 调用 dendropy.calculate.treemeasure.N_bar。
- **API：** `dnakit.phylogeny.calculate_n_bar_index(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树。
- **输出：** n_bar、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_n_bar_index

print(calculate_n_bar_index("((A,B),C);")["result"]["n_bar"])
```

- **示例输出：**

```text
1.666666666667
```

## 29) `dnakit.phylogeny.calculate_node_ages()`

- **作用：** 计算各节点的年龄（逆向时间：节点到最深叶的距离）。
- **处理方法：** 调用 dendropy.calculate.treemeasure.node_ages；除根边外每条边都需要分支长度。
- **API：** `dnakit.phylogeny.calculate_node_ages(tree, *, internal_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树（带分支长度）与 internal_only 布尔值。
- **输出：** internal_only、values（升序数值列表）、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_node_ages

print(calculate_node_ages("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **示例输出：**

```text
[0.0, 2.0, 2.0, 3.0, 3.0, 4.0, 5.0]
```

## 30) `dnakit.phylogeny.calculate_node_depths()`

- **作用：** 计算各节点的深度（正向时间：根到节点的距离）。
- **处理方法：** 调用 dendropy.calculate.treemeasure.node_depths；需要分支长度。
- **API：** `dnakit.phylogeny.calculate_node_depths(tree, *, internal_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树（带分支长度）与 internal_only 布尔值。
- **输出：** internal_only、values（升序数值列表）、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_node_depths

print(calculate_node_depths("((A:1,B:2):1,(C:3,D:1):2);", internal_only=True)["result"]["values"])
```

- **示例输出：**

```text
[0.0, 1.0, 2.0]
```

## 31) `dnakit.phylogeny.calculate_coalescence_ages()`

- **作用：** 计算内部节点（合并事件）的逆向时间年龄。
- **处理方法：** 调用 dendropy.calculate.treemeasure.coalescence_ages；需要分支长度。
- **API：** `dnakit.phylogeny.calculate_coalescence_ages(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树（带分支长度）。
- **输出：** values（升序数值列表）、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_coalescence_ages

print(calculate_coalescence_ages("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **示例输出：**

```text
[3.0, 4.0, 5.0]
```

## 32) `dnakit.phylogeny.calculate_divergence_times()`

- **作用：** 计算内部节点（分化事件）的正向时间，从根开始计。
- **处理方法：** 调用 dendropy.calculate.treemeasure.divergence_times；需要分支长度。
- **API：** `dnakit.phylogeny.calculate_divergence_times(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 树（带分支长度）。
- **输出：** values（升序数值列表）、tip_count 与 provenance。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_divergence_times

print(calculate_divergence_times("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **示例输出：**

```text
[0.0, 1.0, 2.0]
```

## 33) `dnakit.annotation.read_pyranges1_gff()`

- **作用：** 用 pyranges1 读取 GTF 或 GFF 文件，每个特征输出一条记录（含属性列）。
- **处理方法：** 调用 pyranges1.read_gff，最多读取 max_rows+1 行；超出 max_rows 时报错。坐标为 0 起点、半开区间。
- **API：** `dnakit.annotation.read_pyranges1_gff(path, *, max_rows=100000, duplicate_attr=False, ignore_bad=False, verbose=True)`。
- **输入：** GTF 或 GFF 文件路径，可选 max_rows、duplicate_attr、ignore_bad。
- **输出：** records（chrom、start、end、strand 与属性列）、record_count 与 provenance。

**示例文件：** `a.gtf`

```text
chr1	src	gene	100	500	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	130	.	+	0	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	200	260	.	+	0	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	400	450	.	-	0	gene_id "g2"; transcript_id "t2";
```

- **示例代码：**

```python
from dnakit.annotation import read_pyranges1_gff

result = read_pyranges1_gff("a.gtf")
print(result["record_count"])
first = result["records"][0]
print(first["chrom"], first["start"], first["end"], first["Feature"], first["transcript_id"])
```

- **示例输出：**

```text
4
chr1 99 130 CDS t1
```

## 34) `dnakit.annotation.read_pyranges1_gff3()`

- **作用：** 用 pyranges1 读取 GFF3 文件，属性键成为记录字段。
- **处理方法：** 调用 pyranges1.read_gff3；起点由 1 起点转换为 0 起点。
- **API：** `dnakit.annotation.read_pyranges1_gff3(path, *, max_rows=100000, verbose=True)`。
- **输入：** GFF3 文件路径与可选 max_rows。
- **输出：** records、record_count 与 provenance。

**示例文件：** `a.gff3`

```text
##gff-version 3
chr1	src	gene	101	500	.	+	.	ID=g1;Name=G1
chr1	src	CDS	101	130	.	+	0	ID=c1;Parent=g1
```

- **示例代码：**

```python
from dnakit.annotation import read_pyranges1_gff3

for record in read_pyranges1_gff3("a.gff3")["records"]:
    print(record["Feature"], record["start"], record["end"], record["ID"])
```

- **示例输出：**

```text
CDS 100 130 c1
gene 100 500 g1
```

## 35) `dnakit.annotation.read_pyranges1_gtf_full()`

- **作用：** 用 pyranges1 的完整属性解析器读取 GTF 文件。
- **处理方法：** 调用 pyranges1.readers.read_gtf_full，最多读取 max_rows+1 行。
- **API：** `dnakit.annotation.read_pyranges1_gtf_full(path, *, max_rows=100000, duplicate_attr=False, ignore_bad=False, verbose=True)`。
- **输入：** GTF 文件路径，可选 max_rows、duplicate_attr、ignore_bad。
- **输出：** records、record_count 与 provenance。
- **示例代码：**

```python
from dnakit.annotation import read_pyranges1_gtf_full

result = read_pyranges1_gtf_full("a.gtf")
print(result["record_count"], sorted({record["transcript_id"] for record in result["records"]}))
```

- **示例输出：**

```text
4 ['t1', 't2']
```

## 36) `dnakit.annotation.read_pyranges1_bigwig()`

- **作用：** 读取 bigWig 信号轨道的全部区间，每个区间给出数值。
- **处理方法：** 调用 pyranges1.read_bigwig；先检查文件大小（max_file_bytes），再检查 max_rows。
- **API：** `dnakit.annotation.read_pyranges1_bigwig(path, *, max_file_bytes=50000000, max_rows=100000, verbose=True)`。
- **输入：** bigWig 文件路径，可选 max_file_bytes 与 max_rows。
- **输出：** records（chrom、start、end、value）、record_count 与 provenance。
- **示例代码：**

```python
import pyBigWig

from dnakit.annotation import read_pyranges1_bigwig

handle = pyBigWig.open("track.bw", "w")
handle.addHeader([("chr1", 1000), ("chr2", 500)])
handle.addEntries(["chr1", "chr1", "chr2"], [0, 200, 50], ends=[100, 260, 80], values=[1.5, 2.0, 0.25])
handle.close()
for record in read_pyranges1_bigwig("track.bw")["records"]:
    print(record["chrom"], record["start"], record["end"], record["value"])
```

- **示例输出：**

```text
chr1 0 100 1.5
chr1 200 260 2.0
chr2 50 80 0.25
```

## 37) `dnakit.annotation.concatenate_pyranges1_intervals()`

- **作用：** 把多个 BED 文件的区间合并为一个按坐标排序的列表。
- **处理方法：** 逐个调用 pyranges1.read_bed，再用 pyranges1.concat 合并；每个文件先检查行数。
- **API：** `dnakit.annotation.concatenate_pyranges1_intervals(paths, *, max_rows=100000, verbose=True)`。
- **输入：** BED 文件路径列表（最多 50 个）与可选 max_rows。
- **输出：** records、record_count、source_count 与 provenance。

**示例文件：** `x.bed`

```text
chr1	0	100	x1	0	+
chr1	200	300	x2	0	-
chr2	50	80	x3	0	+
```


**示例文件：** `y.bed`

```text
chr1	50	150	y1	0	+
chr2	60	70	y2	0	+
```

- **示例代码：**

```python
from dnakit.annotation import concatenate_pyranges1_intervals

result = concatenate_pyranges1_intervals(["x.bed", "y.bed"])
print(result["record_count"], result["source_count"])
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["Name"])
```

- **示例输出：**

```text
5 2
chr1 0 100 x1
chr1 50 150 y1
chr1 200 300 x2
chr2 50 80 x3
chr2 60 70 y2
```

## 38) `dnakit.annotation.calculate_pyranges1_cds_frame()`

- **作用：** 为 GTF 中的 CDS 特征分配阅读框（frame），按一个或多个属性分组。
- **处理方法：** 用 pyranges1.read_gff 读取，筛选 Feature 为 CDS，再调用 pyranges1.orfs.calculate_frame；frame 为同组之前 CDS 长度之和模 3。
- **API：** `dnakit.annotation.calculate_pyranges1_cds_frame(path, *, feature='CDS', group_by=('transcript_id',), max_rows=100000, verbose=True)`。
- **输入：** GTF/GFF 路径，可选 feature（默认 CDS）、group_by（默认 transcript_id）与 max_rows。
- **输出：** records（chrom、start、end、strand、分组列与 frame）、record_count 与 provenance。
- **示例代码：**

```python
from dnakit.annotation import calculate_pyranges1_cds_frame

for record in calculate_pyranges1_cds_frame("a.gtf")["records"]:
    print(record["transcript_id"], record["start"], record["end"], record["frame"])
```

- **示例输出：**

```text
t1 99 130 0
t1 199 260 1
t2 399 450 0
```

## 39) `dnakit.annotation.extend_pyranges1_orfs()`

- **作用：** 用 FASTA 序列把 BED6 格式的 CDS 区间延伸为完整的开放阅读框（ORF）。
- **处理方法：** 调用 pyranges1.orfs.extend_orfs（需要 pyfaidx）；向上游寻找 ATG 起始密码子，向下游寻找终止密码子 TAG、TGA、TAA。pyfaidx 会在 FASTA 旁写入 .fai 索引。
- **API：** `dnakit.annotation.extend_pyranges1_orfs(path, fasta_path, *, group_by='Name', direction=None, keep_off_bounds=False, chunk_size=900, max_rows=100000, verbose=True)`。
- **输入：** BED6 路径、FASTA 路径，可选 group_by（默认 Name）、direction（None、up 或 down）、keep_off_bounds、chunk_size（3 的倍数）。
- **输出：** records（延伸后的区间）、record_count 与 provenance。

**示例文件：** `cds.bed`

```text
chr1	9	21	g	0	+
```


**示例文件：** `g.fa`

```text
>chr1
GGGATGGCTGCAAAGTGATTTCCCATGCCCGCGTAAGGGCCC
```

- **示例代码：**

```python
from dnakit.annotation import extend_pyranges1_orfs

up = extend_pyranges1_orfs("cds.bed", "g.fa", direction="up")["records"][0]
down = extend_pyranges1_orfs("cds.bed", "g.fa", direction="down")["records"][0]
print(up["start"], up["end"], down["start"], down["end"])
```

- **示例输出：**

```text
3 21 9 36
```

## 40) `dnakit.annotation.calculate_pyranges1_forbes()`

- **作用：** 计算两组 BED 区间的 Forbes 系数，需要给出染色体长度。
- **处理方法：** 调用 pyranges1.stats.forbes；strand_behavior 可为 auto、same、opposite 或 ignore。
- **API：** `dnakit.annotation.calculate_pyranges1_forbes(source, other, *, chromsizes, strand_behavior='auto', max_rows=100000, verbose=True)`。
- **输入：** 两个 BED 路径、chromsizes（染色体到长度的映射）、可选 strand_behavior 与 max_rows。
- **输出：** forbes、strand_behavior 与 provenance。
- **示例代码：**

```python
from dnakit.annotation import calculate_pyranges1_forbes

result = calculate_pyranges1_forbes("x.bed", "y.bed", chromsizes={"chr1": 1000, "chr2": 500})
print(result["forbes"])
```

- **示例输出：**

```text
3.5573122529644268
```

## 41) `dnakit.annotation.calculate_pyranges1_mcc()`

- **作用：** 计算带标签的 BED 区间集之间碱基对重叠的 Matthews 相关系数（MCC）。
- **处理方法：** 调用 pyranges1.stats.mcc；每行把 truth_label 当作参考集、prediction_label 当作预测集，并给出 TP、FP、TN、FN（碱基对）。
- **API：** `dnakit.annotation.calculate_pyranges1_mcc(paths, labels, *, use_strand=False, max_rows=100000, verbose=True)`。
- **输入：** 2 到 10 个 BED 路径、与之等长且唯一的 labels、可选 use_strand 与 max_rows。
- **输出：** records（truth_label、prediction_label、tp、fp、tn、fn、mcc）、record_count 与 provenance。
- **示例代码：**

```python
from dnakit.annotation import calculate_pyranges1_mcc

result = calculate_pyranges1_mcc(["x.bed", "y.bed"], ["x", "y"])
for row in result["records"]:
    print(row["truth_label"], row["prediction_label"], row["tp"], row["fp"], row["tn"], row["fn"], round(row["mcc"], 4))
```

- **示例输出：**

```text
x x 230 0 150 0 1.0
x y 60 50 100 170 -0.0781
y x 60 170 100 50 -0.0781
y y 110 0 270 0 1.0
```

## 42) `dnakit.annotation.calculate_pyranges1_relative_distance()`

- **作用：** 计算两组 BED 区间之间的空间相关性（相对距离分箱）。
- **处理方法：** 调用 pyranges1.stats.relative_distance，并原样输出 reldist、count、total 与 fraction。
- **API：** `dnakit.annotation.calculate_pyranges1_relative_distance(source, other, *, max_rows=100000, verbose=True)`。
- **输入：** 两个 BED 路径与可选 max_rows。
- **输出：** records（reldist、count、total、fraction）、record_count 与 provenance。
- **示例代码：**

```python
from dnakit.annotation import calculate_pyranges1_relative_distance

result = calculate_pyranges1_relative_distance("x.bed", "y.bed")
for record in result["records"]:
    print(record["reldist"], record["count"], record["total"], record["fraction"])
```

- **示例输出：**

```text
0.0 1 1 1.0
```

## 验证方式

- 测试文件 `others/tests/unit/test_review_G18p.py`，共 31 个测试，全部离线运行，不使用网络。
- 树形与分支结果用小树的手算值检验（例如 `((A:1,B:2):1,(C:3,D:1):2);` 的节点距离 7、拓扑距离 4、节点深度 `[0, 1, 2, 2, 3, 3, 5]`；ORF 坐标由 42 bp 序列手算）；矩阵、网络解析、Forbes 系数、相对距离等与 toytree、dendropy、pyranges1 的直接调用逐项比较。
- 每个模块各有一个缺包测试（把上游模块置为 `None`，期望 `BackendUnavailableError` 且提示 extra 名称），另有输入校验测试（布尔值、节点查询、未知叶名、剩余叶数、VCV 对称性等）。

## 当前范围

- 已封装 42 个函数。toytree 的 `iter_node_path`、`iter_descendant_dists` 与 `io.parse_newick_string_custom` 分别与列表或读取封装重复；`mtree` 与 `io.mtree` 是构造器（internal）；`get_markov_model` 只返回模型对象（internal）；`simulate_discrete_data` 为任意状态数的离散性状模拟（not-dna，DNA 情形见 `simulate_nucleotide_sequences`）；`utils.show` 在浏览器中显示画布（plot）；`utils.network.parse_network` 为带调试输出的旧原型（deprecated）；`get_tree_from_vcv_matrix` 调用了 toytree 3.0.11 中不存在的 `infer_neighbor_joining_tree`（broken-upstream）。
- dendropy 的 `tree_source_iter` 与 `multi_tree_source_iter` 在 DendroPy 4 之后直接抛出 NotImplementedError（deprecated）；`mason_gamer_kellogg_score` 总是抛出 AttributeError（broken-upstream）；四个 `simulate_discrete_*` 为通用离散字符模拟（not-dna）。
- pyranges1 的 `read_bam` 需要未安装的 `bamread`，且缺包时上游直接 `sys.exit`，因此不封装（needs-binary:bamread）；`random` 为测试用随机区间生成器，`seqs.clear_kmer_memory` 为缓存函数（internal）；`version_info` 为版本报告（not-dna）；其余别名与上述函数相同（duplicate-of）。
- 上游 `pyranges1.orfs` 与 `extend_orfs` 通过 pyfaidx 写入 `.fai` 索引；toytree 的 `write()` 不写内部节点名，封装结果中的 `name` 字段是输入的回显。
