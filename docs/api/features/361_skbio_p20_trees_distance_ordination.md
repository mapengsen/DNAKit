# scikit-bio 系统发育树、距离矩阵与排序

本页介绍 P20 的树操作和距离矩阵统计接口。它们位于 `dnakit.phylogeny.skbio_p20_tree_functions` 与 `dnakit.microbiome.skbio_p20_stats_functions`，依赖可选项 `bio-skibio`（scikit-bio 0.7.4）。输入/输出均为 JSON 友好的 Newick 字符串、列表和字典；函数按需导入 scikit-bio。

## 1) `reroot_tree()`

- **作用：** 将树根移至指定节点，或其上方分支。
- **处理方法：** 使用 `TreeNode.root_at`、`unrooted_copy` 或 `unrooted_move`；节点可用名称或一组叶名称指定。
- **API：** `dnakit.phylogeny.skbio_p20_tree_functions.reroot_tree(newick, node, *, above=False, method='root_at', verbose=True)`。
- **输入：** Newick 字符串、节点名称或叶名列表。
- **输出：** 新 Newick 字符串、方法和叶数。
- **示例代码：**

```python
from dnakit.phylogeny.skbio_p20_tree_functions import reroot_tree, tree_diameter

tree = "((a:1,b:2)90:1,(c:3,d:1)0.5:2,e:1);"
print(reroot_tree(tree, "c", verbose=False)["newick"])
print(tree_diameter(tree, verbose=False)["distance"])
```

- **示例输出：**

```text
((d:1.0,((a:1.0,b:2.0)90:1.0,e:1.0):2.0)0.5:3.0)c;
8.0
```

## 2) 距离矩阵检验与排序

- **作用：** 对距离矩阵执行 ANOSIM、PERMDISP、Mantel 检验、PCoA 双标图投影和受约束排序。
- **处理方法：** `anosim_test`、`permdisp_test`、`mantel_test` 与 `pcoa_biplot_projection` 调用 scikit-bio 对应算法；`constrained_ordination` 支持 `rda` 和 `cca`，`plot_ordination` 保存图像并返回坐标数据。
- **API：** 例如 `dnakit.microbiome.skbio_p20_stats_functions.anosim_test(distance_matrix, grouping, *, ids=None, permutations=999, seed=None, verbose=True)`。
- **输入：** 对称距离矩阵、分组标签；不同排序函数还接受样本/特征表。
- **输出：** 检验统计量、p 值、排序坐标或图像路径与坐标数据。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p20_stats_functions import anosim_test

result = anosim_test(
    [[0, 1, 4, 5], [1, 0, 5, 4], [4, 5, 0, 1], [5, 4, 1, 0]],
    ["A", "A", "B", "B"],
    permutations=99,
    seed=7,
    verbose=False,
)
print(result["test_statistic_name"], result["number_of_groups"])
```

- **示例输出：**

```text
R 2
```
