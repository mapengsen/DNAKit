# Biopython GenePop 编辑、聚类与基序打分

本文档涵盖对 [Biopython](https://biopython.org) 1.88 的三组封装：`dnakit.population` 中的 GenePop 文件编辑（`Bio.PopGen.GenePop`）、同在 `dnakit.population` 的 `Bio.Cluster`（C Clustering Library）聚类例程，以及 `dnakit.regulation` 中的 `Bio.motifs` 打分。矩阵为行列表；距离名称为 `euclidean`、`cityblock`、`pearson`、`abs_pearson`、`uncentered`、`abs_uncentered`、`spearman` 和 `kendall`。k 均值、k 中心点和自组织映射的随机起点来自 C 库，无法设定种子；传 `initial_clusters` 可使其确定。用 `python -m pip install "dnakit[bio-formats]"` 安装 Biopython，只在函数运行时导入。WebLogo 示例需要联网。示例从模块 `dnakit.population.biopython_p19_cluster_functions` 与 `dnakit.regulation.biopython_p19_motif_functions` 导入。

## 示例数据

`pops.gen` 是一个小的 GenePop 文件，含三个群体和三个位点（两位数等位基因）。聚类和基序示例直接传入数据：五个条目、每个三个测量值，以及七个长度为 5 的 DNA 结合位点。

**示例文件：** `pops.gen`

```text
Test genepop
Loc1
Loc2
Loc3
Pop
Ind1, 0101 0202 0303
Ind2, 0102 0203 0000
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
Pop
Ind5, 0202 0101 0101
```

## 1) `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_population()`

- **作用：** 把 GenePop 文件的每个群体写成一个单独的 GenePop 文件。
- **处理方法：** 用 `Bio.PopGen.GenePop.read` 读入文件并调用 `Record.split_in_pops`。GenePop 文件中的群体名不可靠，因此默认按文件顺序命名为 `pop1.gen`、`pop2.gen` ……。每个输出保留所有位点。输出目录必须尚不存在。
- **API：** `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_population(genepop_path, output_dir, *, population_names=None)`。
- **输入：** GenePop 文件路径、新的 `output_dir` 和可选的 `population_names`（每个群体一个）。
- **输出：** `population_count`、`locus_count` 和 `files`（名称到路径）。
- **示例代码：**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import split_genepop_by_population

target = os.path.join(tempfile.mkdtemp(), "pops")
result = split_genepop_by_population("pops.gen", target)
print(result["population_count"], result["locus_count"], sorted(result["files"]))
print(open(result["files"]["pop2"]).read())
```

- **示例输出：**

```text
3 3 ['pop1', 'pop2', 'pop3']
Test genepop
Loc1
Loc2
Loc3
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
```

## 2) `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_locus()`

- **作用：** 把 GenePop 文件的每个位点写成一个单独的 GenePop 文件。
- **处理方法：** 调用 `Bio.PopGen.GenePop.Record.split_in_loci`。每个输出 `<位点>.gen` 含一个位点和全部群体。输出目录必须尚不存在。
- **API：** `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_locus(genepop_path, output_dir)`。
- **输入：** GenePop 文件路径和新的 `output_dir`。
- **输出：** `population_count`、`locus_count` 和 `files`（位点名到路径）。
- **示例代码：**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import split_genepop_by_locus

target = os.path.join(tempfile.mkdtemp(), "loci")
result = split_genepop_by_locus("pops.gen", target)
print(result["locus_count"], sorted(result["files"]))
print(open(result["files"]["Loc2"]).read())
```

- **示例输出：**

```text
3 ['Loc1', 'Loc2', 'Loc3']
Test genepop
Loc2
Pop
Ind1, 0202
Ind2, 0203
Pop
Ind3, 0202
Ind4, 0102
Pop
Ind5, 0101
```

## 3) `dnakit.population.biopython_p19_cluster_functions.remove_genepop_populations()`

- **作用：** 写出去掉部分群体的 GenePop 文件副本。
- **处理方法：** 对每个 0 起始的位置（从大到小）调用 `Bio.PopGen.GenePop.Record.remove_population` 并写出记录。不使用基于文件的 `FileRecord` 变体，因为它在 Biopython 1.88 中写出的等位基因补零有误。至少要保留一个群体；输出路径必须尚不存在。
- **API：** `dnakit.population.biopython_p19_cluster_functions.remove_genepop_populations(genepop_path, output_path, positions)`。
- **输入：** GenePop 文件路径、新的 `output_path` 和要去掉的群体 `positions`。
- **输出：** `removed_positions`、剩余的 `population_count` 和 `output` 路径。
- **示例代码：**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import remove_genepop_populations

target = os.path.join(tempfile.mkdtemp(), "kept.gen")
result = remove_genepop_populations("pops.gen", target, [0, 2])
print(result["removed_positions"], result["population_count"])
print(open(target).read())
```

- **示例输出：**

```text
[0, 2] 1
Test genepop
Loc1
Loc2
Loc3
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
```

## 4) `dnakit.population.biopython_p19_cluster_functions.remove_genepop_loci()`

- **作用：** 写出去掉部分位点的 GenePop 文件副本。
- **处理方法：** 调用 `Bio.PopGen.GenePop` 的 `Record.remove_locus_by_position` 或 `Record.remove_locus_by_name` 并写出记录。`positions`（0 起始）与 `names` 二选一；至少要保留一个位点。
- **API：** `dnakit.population.biopython_p19_cluster_functions.remove_genepop_loci(genepop_path, output_path, *, positions=None, names=None)`。
- **输入：** GenePop 文件路径、新的 `output_path`，以及 `positions` 或 `names`。
- **输出：** `removed_loci`、`remaining_loci` 和 `output` 路径。
- **示例代码：**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import remove_genepop_loci

target = os.path.join(tempfile.mkdtemp(), "kept.gen")
result = remove_genepop_loci("pops.gen", target, names=["Loc2"])
print(result["removed_loci"], result["remaining_loci"])
print(open(target).read())
```

- **示例输出：**

```text
['Loc2'] ['Loc1', 'Loc3']
Test genepop
Loc1
Loc3
Pop
Ind1, 0101 0303
Ind2, 0102 0000
Pop
Ind3, 0201 0303
Ind4, 0101 0103
Pop
Ind5, 0202 0101
```

## 5) `dnakit.population.biopython_p19_cluster_functions.kmeans_cluster()`

- **作用：** 用 k 均值或 k 中位数对条目聚类。
- **处理方法：** 调用 `Bio.Cluster.kcluster`。`method` 为 `mean`（k 均值）或 `median`（k 中位数）；`distance` 为 `euclidean`、`cityblock`、`pearson`、`abs_pearson`、`uncentered`、`abs_uncentered`、`spearman` 或 `kendall` 之一。EM 算法从随机起点运行 `n_passes` 次，或从 `initial_clusters` 出发只运行一次（完全确定）。
- **API：** `dnakit.population.biopython_p19_cluster_functions.kmeans_cluster(data, n_clusters, *, method='mean', distance='euclidean', n_passes=1, initial_clusters=None, mask=None, weights=None, transpose=False)`。
- **输入：** 矩阵（行列表，条目×测量值）、`n_clusters`、`method`、`distance`、`n_passes`、`initial_clusters`、`mask`（0 表示缺失）、`weights` 和 `transpose`（对列而非行聚类）。
- **输出：** `cluster_ids`、簇内距离之和 `error` 和 `times_found`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import kmeans_cluster

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
means = kmeans_cluster(data, 2, initial_clusters=[0, 0, 1, 1, 0])
print(means["cluster_ids"], round(means["error"], 4), means["times_found"])
medians = kmeans_cluster(data, 2, method="median", distance="cityblock", initial_clusters=[0, 0, 1, 1, 0])
print(medians["cluster_ids"], round(medians["error"], 4))
```

- **示例输出：**

```text
[0, 0, 1, 1, 0] 0.0256 1
[0, 0, 1, 1, 0] 0.3
```

## 6) `dnakit.population.biopython_p19_cluster_functions.kmedoids_cluster()`

- **作用：** 根据距离矩阵用 k 中心点法对条目聚类。
- **处理方法：** 调用 `Bio.Cluster.kmedoids`；只读取矩阵的下三角。每个条目的簇号是其中心点条目的编号。
- **API：** `dnakit.population.biopython_p19_cluster_functions.kmedoids_cluster(distance_matrix, n_clusters, *, n_passes=1, initial_clusters=None)`。
- **输入：** 方阵距离矩阵、`n_clusters`、`n_passes` 和 `initial_clusters`（每个条目一个小于 `n_clusters` 的起始簇号）。
- **输出：** `cluster_ids`（中心点条目编号）、`error` 和 `times_found`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import kmedoids_cluster, profile_distance_matrix

matrix = profile_distance_matrix([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]])["matrix"]
result = kmedoids_cluster(matrix, 2, initial_clusters=[0, 0, 1, 1, 0])
print(result["cluster_ids"], round(result["error"], 4), result["times_found"])
```

- **示例输出：**

```text
[4, 4, 2, 2, 4] 0.0467 1
```

## 7) `dnakit.population.biopython_p19_cluster_functions.hierarchical_cluster()`

- **作用：** 用单链、全链、平均链或质心链接做层次聚类。
- **处理方法：** 调用 `Bio.Cluster.treecluster`。输入数据矩阵（四种链接都可用）或方阵距离矩阵（只支持单链、全链和平均链）。树以合并节点列表表示：`i >= 0` 是叶子，`-k` 指第 `k` 步（从 1 起）产生的节点，与 Cluster 3.0 一致。
- **API：** `dnakit.population.biopython_p19_cluster_functions.hierarchical_cluster(data=None, *, distance_matrix=None, linkage='complete', distance='euclidean', mask=None, weights=None, transpose=False)`。
- **输入：** `data` 或 `distance_matrix`（二选一）、`linkage`、`distance`、`mask`、`weights` 和 `transpose`。
- **输出：** `item_count` 和 `nodes`（`left`、`right`、`distance`）。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for linkage in ("single", "complete", "average", "centroid"):
    nodes = hierarchical_cluster(data, linkage=linkage)["nodes"]
    print(linkage, [(n["left"], n["right"]) for n in nodes], round(nodes[-1]["distance"], 3))
```

- **示例输出：**

```text
single [(0, 4), (2, 3), (1, -1), (-2, -3)] 47.153
complete [(4, 0), (3, 2), (1, -1), (-2, -3)] 49.98
average [(4, 0), (3, 2), (1, -1), (-2, -3)] 48.866
centroid [(0, 4), (-1, 1), (3, 2), (-2, -3)] 48.855
```

## 8) `dnakit.population.biopython_p19_cluster_functions.cut_cluster_tree()`

- **作用：** 把层次聚类树切成指定数量的簇。
- **处理方法：** 由节点列表重建 `Bio.Cluster.Tree` 并调用 `Tree.cut`。
- **API：** `dnakit.population.biopython_p19_cluster_functions.cut_cluster_tree(nodes, n_clusters)`。
- **输入：** `hierarchical_cluster` 返回的 `nodes` 和 `n_clusters`。
- **输出：** 各条目的 `cluster_ids`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster, cut_cluster_tree

nodes = hierarchical_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], linkage="average")["nodes"]
for count in (2, 3):
    print(count, cut_cluster_tree(nodes, count)["cluster_ids"])
```

- **示例输出：**

```text
2 [1, 1, 0, 0, 1]
3 [2, 1, 0, 0, 2]
```

## 9) `dnakit.population.biopython_p19_cluster_functions.sort_cluster_tree()`

- **作用：** 为层次聚类树的叶子排序。
- **处理方法：** 调用 `Bio.Cluster.Tree.sort`：交换左右子节点，使叶子顺序尽量符合 `order`（默认为条目编号）。
- **API：** `dnakit.population.biopython_p19_cluster_functions.sort_cluster_tree(nodes, *, order=None)`。
- **输入：** 树的 `nodes` 和可选的 `order`（每个条目一个数）。
- **输出：** `leaf_order`（从左到右的条目编号）和排序后的 `nodes`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster, sort_cluster_tree

nodes = hierarchical_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], linkage="average")["nodes"]
print(sort_cluster_tree(nodes)["leaf_order"])
print(sort_cluster_tree(nodes, order=[4, 3, 2, 1, 0])["leaf_order"])
```

- **示例输出：**

```text
[3, 2, 1, 4, 0]
[3, 2, 4, 0, 1]
```

## 10) `dnakit.population.biopython_p19_cluster_functions.som_cluster()`

- **作用：** 在矩形网格上用自组织映射对条目聚类。
- **处理方法：** 调用 `Bio.Cluster.somcluster`。映射从随机向量开始，因此每次调用结果可能不同。
- **API：** `dnakit.population.biopython_p19_cluster_functions.som_cluster(data, *, nx_grid=2, ny_grid=1, init_tau=0.02, n_iterations=100, distance='euclidean', mask=None, weights=None, transpose=False)`。
- **输入：** 矩阵、网格大小 `nx_grid` 与 `ny_grid`、`init_tau`、`n_iterations`、`distance`、`mask`、`weights` 和 `transpose`。
- **输出：** `cell_ids`（每个条目所在网格的 x、y）和 `cell_data[x][y]`（每个网格的质心向量）。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import som_cluster

result = som_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], nx_grid=2, ny_grid=1, n_iterations=50)
print(len(result["cell_ids"]), len(result["cell_data"]), len(result["cell_data"][0][0]))
print(all(0 <= x < 2 and y == 0 for x, y in result["cell_ids"]))
```

- **示例输出：**

```text
5 2 3
True
```

## 11) `dnakit.population.biopython_p19_cluster_functions.cluster_distance()`

- **作用：** 计算两个条目簇之间的距离。
- **处理方法：** 调用 `Bio.Cluster.clusterdistance`。`method` 为 `mean` 或 `median`（中心之间的距离）、`single`（最小两两距离）、`complete`（最大）或 `average`（两两距离的平均）。
- **API：** `dnakit.population.biopython_p19_cluster_functions.cluster_distance(data, index1, index2, *, method='mean', distance='euclidean', mask=None, weights=None, transpose=False)`。
- **输入：** 矩阵、两个簇的成员编号 `index1` 与 `index2`、`method`、`distance`、`mask`、`weights` 和 `transpose`。
- **输出：** 簇间 `distance`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_distance

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for method in ("mean", "median", "single", "complete", "average"):
    print(method, round(cluster_distance(data, [0, 1], [2, 3], method=method)["distance"], 4))
```

- **示例输出：**

```text
mean 48.5417
median 48.5417
single 47.1533
complete 49.9533
average 48.5517
```

## 12) `dnakit.population.biopython_p19_cluster_functions.cluster_centroids()`

- **作用：** 计算每个簇的质心。
- **处理方法：** 调用 `Bio.Cluster.clustercentroids`：对每个簇的条目逐维取平均（`method="mean"`）或中位数。
- **API：** `dnakit.population.biopython_p19_cluster_functions.cluster_centroids(data, cluster_ids, *, method='mean', mask=None, transpose=False)`。
- **输入：** 矩阵、`cluster_ids`（每个条目一个非负簇号）、`method`、`mask` 和 `transpose`。
- **输出：** 每个簇一行质心（`centroids`）和对应的 0/1 `centroid_mask`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_centroids

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
print(cluster_centroids(data, [0, 0, 1, 1, 0])["centroids"])
print(cluster_centroids(data, [0, 0, 1, 1, 0], method="median")["centroids"])
```

- **示例输出：**

```text
[[1.0, 2.033333333333333, 3.1], [8.1, 9.05, 9.95]]
[[1.0, 2.0, 3.1], [8.1, 9.05, 9.95]]
```

## 13) `dnakit.population.biopython_p19_cluster_functions.profile_distance_matrix()`

- **作用：** 计算数据矩阵各谱（行或列）之间的距离矩阵。
- **处理方法：** 调用 `Bio.Cluster.distancematrix`。`distance` 为 `euclidean`、`cityblock`、`pearson`（1 - Pearson r）、`abs_pearson`（1 - |r|）、`uncentered`（1 - 未中心化相关）、`abs_uncentered`、`spearman`（1 - 秩相关）或 `kendall`（1 - Kendall tau）。缺失值通过 `mask` 排除。
- **API：** `dnakit.population.biopython_p19_cluster_functions.profile_distance_matrix(data, *, distance='euclidean', mask=None, weights=None, transpose=False)`。
- **输入：** 矩阵、`distance`、`mask`、`weights` 和 `transpose`。
- **输出：** `item_count` 和完整的对称 `matrix`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import profile_distance_matrix

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for distance in ("euclidean", "cityblock", "pearson", "abs_pearson", "uncentered", "abs_uncentered", "spearman", "kendall"):
    row = profile_distance_matrix(data, distance=distance)["matrix"][0]
    print(distance, [round(value, 4) for value in row])
```

- **示例输出：**

```text
euclidean [0.0, 0.02, 49.0, 49.9533, 0.0067]
cityblock [0.0, 0.1333, 7.0, 7.0667, 0.0667]
pearson [0.0, 0.0004, 0.0, 0.0006, -0.0]
abs_pearson [0.0, 0.0004, 0.0, 0.0006, -0.0]
uncentered [0.0, 0.0001, 0.0438, 0.048, 0.0006]
abs_uncentered [0.0, 0.0001, 0.0438, 0.048, 0.0006]
spearman [0.0, 0.0, 0.0, 0.0, 0.0]
kendall [0.0, 0.0, 0.0, 0.0, 0.0]
```

## 14) `dnakit.population.biopython_p19_cluster_functions.cluster_pca()`

- **作用：** 数据矩阵的主成分分析。
- **处理方法：** 调用 `Bio.Cluster.pca`（中心化奇异值分解）。结果按特征值从大到小排序；列均值加上 `coordinates @ components` 可还原数据。
- **API：** `dnakit.population.biopython_p19_cluster_functions.cluster_pca(data)`。
- **输入：** 至少两行两列的矩阵。
- **输出：** `column_means`、`components`、`coordinates` 和 `eigenvalues`。
- **示例代码：**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_pca

result = cluster_pca([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]])
print([round(value, 4) for value in result["column_means"]])
print([round(value, 4) for value in result["eigenvalues"]])
print([round(value, 4) for value in result["coordinates"][0]])
```

- **示例输出：**

```text
[3.84, 4.84, 5.84]
[13.2636, 0.181, 0.0368]
[-4.9185, -0.0734, 0.0035]
```

## 15) `dnakit.regulation.biopython_p19_motif_functions.motif_matrices()`

- **作用：** 由结合位点构建基序，并返回其计数、概率和对数几率矩阵。
- **处理方法：** `Bio.motifs.create` 生成计数矩阵，`Motif.pwm` 在加入 `pseudocounts` 后将其归一化为位置权重矩阵，`Motif.pssm` 给出相对 `background`（默认均匀）的 log2 几率矩阵。保留 Biopython 默认的零伪计数，因此未见过的碱基得分为负无穷，显示为 `null`。同时返回一致序列、简并一致序列和反一致序列。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.motif_matrices(instances=None, *, counts=None, pseudocounts=None, background=None)`。
- **输入：** 等长 DNA `instances` 或 `counts`（碱基到计数列表）、`pseudocounts`（数值或逐碱基映射）和 `background`（碱基概率）。
- **输出：** 形如 `{碱基: [每位置的值]}` 的 `counts`、`pwm`、`pssm`，以及 `background`、`consensus`、`degenerate_consensus` 和 `anticonsensus`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import motif_matrices

sites = ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]
plain = motif_matrices(sites)
print(plain["consensus"], plain["degenerate_consensus"], plain["anticonsensus"])
print(plain["counts"]["A"], [round(v, 3) for v in plain["pwm"]["A"]])
print(plain["pssm"]["C"])
smooth = motif_matrices(sites, pseudocounts=0.5)
print([round(v, 3) for v in smooth["pssm"]["C"]])
```

- **示例输出：**

```text
TACGC WACVC CCATG
[3.0, 7.0, 0.0, 2.0, 1.0] [0.429, 1.0, 0.0, 0.286, 0.143]
[None, None, 1.5145731728297582, 0.19264507794239583, 1.7776075786635521]
[-2.17, -2.17, 1.29, 0.152, 1.531]
```

## 16) `dnakit.regulation.biopython_p19_motif_functions.motif_anticonsensus()`

- **作用：** 返回基序的反一致序列：每个位置上最不常见的碱基。
- **处理方法：** 调用 `Bio.motifs.Motif.anticonsensus`。它是与基序匹配最差的序列；同时返回一致序列和简并一致序列以便比较。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.motif_anticonsensus(instances=None, *, counts=None)`。
- **输入：** 等长 DNA `instances` 或 `counts`。
- **输出：** `anticonsensus`、`consensus` 和 `degenerate_consensus`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import motif_anticonsensus

result = motif_anticonsensus(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])
print(result["anticonsensus"], result["consensus"], result["degenerate_consensus"])
```

- **示例输出：**

```text
CCATG TACGC WACVC
```

## 17) `dnakit.regulation.biopython_p19_motif_functions.pssm_search()`

- **作用：** 用基序扫描 DNA 序列，报告得分不低于阈值的命中。
- **处理方法：** 在对数几率矩阵上调用 `Bio.motifs.matrix.PositionSpecificScoringMatrix.search`（默认双链）。每个命中包含 Biopython 的 `position`（反链为负数，加上序列长度得到正向起点）、正向 0 起始半开区间的 `start` 与 `end`、`strand` 和 `score`。这里 `pseudocounts` 默认 0.5，使未见过的碱基不会得负无穷；传 0 即为 Biopython 默认。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.pssm_search(sequence, instances=None, *, counts=None, threshold=0.0, both_strands=True, pseudocounts=0.5, background=None, max_hits=100000)`。
- **输入：** DNA `sequence`（A、C、G、T、N）、基序 `instances` 或 `counts`、`threshold`、`both_strands`、`pseudocounts`、`background` 和 `max_hits`。
- **输出：** `hit_count`、`truncated` 和 `hits`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_search

result = pssm_search("TACACGTACGCAACGCATTTGTAAC", ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], threshold=3.0)
print(result["hit_count"])
for hit in result["hits"]:
    print(hit["position"], hit["start"], hit["end"], hit["strand"], round(hit["score"], 3))
```

- **示例输出：**

```text
5
0 0 5 + 5.709
6 6 11 + 6.194
11 11 16 + 5.832
-11 14 19 - 4.694
-7 18 23 - 3.594
```

## 18) `dnakit.regulation.biopython_p19_motif_functions.pssm_score_statistics()`

- **作用：** 计算基序 PSSM 得分的均值、标准差、最大值和最小值。
- **处理方法：** `mean` 与 `std` 是按 `background` 抽取的随机序列得分的期望和标准差（`PositionSpecificScoringMatrix.mean` 和 `.std`）；`max` 与 `min` 是一致序列和反一致序列的得分。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.pssm_score_statistics(instances=None, *, counts=None, pseudocounts=0.5, background=None)`。
- **输入：** 基序 `instances` 或 `counts`、`pseudocounts`（默认 0.5）和 `background`。
- **输出：** `mean`、`std`、`max` 和 `min`（无穷大为 `null`）。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_score_statistics

result = pssm_score_statistics(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])
print({key: round(result[key], 4) for key in ("mean", "std", "max", "min")})
skewed = pssm_score_statistics(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], background={"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4})
print(round(skewed["mean"], 4))
```

- **示例输出：**

```text
{'mean': 3.1603, 'std': 2.5129, 'max': 6.1944, 'min': -10.8496}
4.4366
```

## 19) `dnakit.regulation.biopython_p19_motif_functions.pssm_score_distribution()`

- **作用：** 近似计算 PSSM 在基序模型和背景模型下的得分分布。
- **处理方法：** 调用 `PositionSpecificScoringMatrix.distribution`（`Bio.motifs.thresholds.ScoreDistribution`），在 `precision * 宽度` 个得分区间上做动态规划。每个区间的得分为 `min_score + i * step`，并给出该得分在基序模型和背景下的概率。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.pssm_score_distribution(instances=None, *, counts=None, pseudocounts=0.5, background=None, precision=1000)`。
- **输入：** 基序 `instances` 或 `counts`、`pseudocounts`、`background` 和 `precision`。
- **输出：** `min_score`、`step`、`bin_count`、`scores`、`motif_density` 和 `background_density`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_score_distribution

result = pssm_score_distribution(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], precision=100)
print(result["bin_count"], round(result["min_score"], 3), round(result["step"], 4))
print(round(sum(result["motif_density"]), 6), round(sum(result["background_density"]), 6))
best = max(range(result["bin_count"]), key=lambda i: result["motif_density"][i])
print(round(result["scores"][best], 3))
```

- **示例输出：**

```text
500 -10.85 0.0342
1.0 1.0
5.716
```

## 20) `dnakit.regulation.biopython_p19_motif_functions.pssm_threshold()`

- **作用：** 为基序 PSSM 选择得分阈值。
- **处理方法：** 基于得分分布（`Bio.motifs.thresholds.ScoreDistribution`）：`fpr` 给出相对背景的假阳性率 `rate` 对应的阈值（`threshold_fpr`），`fnr` 给出相对基序的假阴性率对应的阈值（`threshold_fnr`），`balanced` 给出假阴性率为假阳性率 `rate_proportion` 倍时的阈值（`threshold_balanced`，同时返回率），`patser` 给出 Patser 程序使用的阈值（`threshold_patser`）。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.pssm_threshold(instances=None, *, counts=None, method='fpr', rate=None, rate_proportion=1.0, pseudocounts=0.5, background=None, precision=1000)`。
- **输入：** 基序 `instances` 或 `counts`、`method`、`rate`（`fpr` 与 `fnr` 必需）、`rate_proportion`、`pseudocounts`、`background` 和 `precision`。
- **输出：** `threshold`（`balanced` 还有 `rate`）。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_threshold

sites = ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]
print(round(pssm_threshold(sites, method="fpr", rate=0.01)["threshold"], 4))
print(round(pssm_threshold(sites, method="fnr", rate=0.1)["threshold"], 4))
balanced = pssm_threshold(sites, method="balanced")
print(round(balanced["threshold"], 4), round(balanced["rate"], 4))
print(round(pssm_threshold(sites, method="patser")["threshold"], 4))
```

- **示例输出：**

```text
4.2101
-0.4336
0.1698 0.126
0.4221
```

## 21) `dnakit.regulation.biopython_p19_motif_functions.jaspar_pseudocounts()`

- **作用：** 按 JASPAR 的方式计算伪计数。
- **处理方法：** 调用 `Bio.motifs.jaspar.calculate_pseudocounts`：位点数的平方根乘以每个碱基的背景概率。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.jaspar_pseudocounts(instances=None, *, counts=None, background=None)`。
- **输入：** 基序 `instances` 或 `counts` 和 `background`。
- **输出：** `pseudocounts`（每个碱基一个数）。
- **示例代码：**

```python
from dnakit.regulation.biopython_p19_motif_functions import jaspar_pseudocounts

print(jaspar_pseudocounts(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])["pseudocounts"])
print(jaspar_pseudocounts(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], background={"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3})["pseudocounts"])
```

- **示例输出：**

```text
{'A': 0.6614378277661477, 'C': 0.6614378277661477, 'G': 0.6614378277661477, 'T': 0.6614378277661477}
{'A': 0.7937253933193772, 'C': 0.5291502622129182, 'G': 0.5291502622129182, 'T': 0.7937253933193772}
```

## 22) `dnakit.regulation.biopython_p19_motif_functions.weblogo_service_logo()`

- **作用：** 用 WebLogo 网络服务绘制序列标识图。
- **处理方法：** 调用 `Bio.motifs.Motif.weblogo`，把基序发送到 https://weblogo.threeplusone.com（WebLogo 3.3），并把返回的图像保存为新文件。`options` 是 WebLogo 表单字段（未知名称被拒绝）。需要联网；不依赖服务的本地标识图用 `plot_pwm_logo`。
- **API：** `dnakit.regulation.biopython_p19_motif_functions.weblogo_service_logo(instances, output_path, *, counts=None, image_format='png', options=None)`。
- **输入：** 基序 `instances`（或 `counts`）、新的 `output_path`、`image_format`（`png`、`png_print`、`pdf`、`jpeg`、`eps`、`svg` 或 `txt`）和 `options`。
- **输出：** `output` 路径、`image_format`、`bytes` 和所用 `service`。
- **示例代码：**

```python
import os, tempfile
from dnakit.regulation.biopython_p19_motif_functions import weblogo_service_logo

target = os.path.join(tempfile.mkdtemp(), "logo.png")
result = weblogo_service_logo(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], target, options={"logo_title": "demo"})
print(result["status"], result["image_format"], result["bytes"] > 1000, open(target, "rb").read(4))
```

- **示例输出：**

```text
completed png True b'\x89PNG'
```
