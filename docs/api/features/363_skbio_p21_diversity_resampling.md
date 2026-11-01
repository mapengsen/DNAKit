# scikit-bio 多样性、重采样与协同进化检验

这一组函数位于 `dnakit.microbiome.skbio_p21_functions`，是对 [scikit-bio](https://scikit.bio/) 的适配器：17 个 alpha 多样性指标、UniFrac 距离、计数稀释、分组随机抽样、Hommola 共物种形成检验和梯度 ANOVA。scikit-bio 是可选依赖（`dnakit[bio-skibio]`），只有调用函数时才导入。所有输入输出都是 JSON 友好的列表和映射。

## 1) `dnakit.microbiome.skbio_p21_functions.calculate_alpha_diversity()`

- **作用：** 计算单个样本计数向量的一个 scikit-bio alpha 多样性指标。
- **处理方法：** `metric` 选择 `skbio.diversity.alpha` 中的函数：margalef、mcintosh_d、mcintosh_e、menhinick、michaelis_menten_fit、sobs（观察到的特征数）、osd、phydiv、pielou_e、renyi、robbins、shannon、simpson、simpson_e、singles、strong、tsallis。该函数合并了 17 个指标，各指标只使用相关的选项（`base`、`exp`、`order`、`finite`、`num_repeats`、`taxa`、`tree_newick`、`rooted`、`weight`）。michaelis_menten_fit 在 scikit-bio 内部随机抽样且没有种子，所以每次调用结果会有差异。
- **API：** `calculate_alpha_diversity(counts, *, metric='shannon', base=None, exp=False, order=2.0, finite=False, num_repeats=1, taxa=None, tree_newick=None, rooted=None, weight=False)`。
- **输入：** 整数计数列表；phydiv 还需要 `taxa`（与计数等长）和 Newick 树。
- **输出：** `value`（osd 返回 observed/singles/doubles 三项的字典）、`n_features`、`total_count` 和来源信息。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_alpha_diversity
counts = [1, 2, 0, 5, 1, 7]
for metric in ("shannon", "simpson", "pielou_e", "osd", "sobs"):
    print(metric, calculate_alpha_diversity(counts, metric=metric, verbose=False)["value"])
print(calculate_alpha_diversity(counts, metric="shannon", base=2, verbose=False)["value"])
```

- **示例输出：**

```text
shannon 1.331660286822432
simpson 0.6875
pielou_e 0.8274070571684498
osd {'observed': 5, 'singles': 2, 'doubles': 1}
sobs 5.0
1.9211796919474975
```

## 2) `dnakit.microbiome.skbio_p21_functions.calculate_unifrac_distance()`

- **作用：** 在系统发育树上计算样本两两之间的 UniFrac 距离矩阵。
- **处理方法：** `weighted=False` 为非加权 UniFrac（`beta.unweighted_unifrac`），`weighted=True` 为加权 UniFrac（`beta.weighted_unifrac`），`normalized` 选择归一化形式。通过 `beta_diversity` 对所有样本对计算。
- **API：** `calculate_unifrac_distance(counts, taxa, tree_newick, *, ids=None, weighted=False, normalized=False)`。
- **输入：** 样本×特征的计数矩阵（至少 2 个样本）、与列对应的树末端名称 `taxa`、Newick 字符串。
- **输出：** `ids` 与对称的 `distance_matrix`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_unifrac_distance
counts = [[1, 2, 0, 5, 1], [3, 0, 4, 1, 0], [0, 0, 1, 2, 9]]
tree = "(((a:1,b:1):1,c:2):1,(d:1,e:1):2);"
for weighted in (False, True):
    result = calculate_unifrac_distance(counts, list("abcde"), tree, weighted=weighted, normalized=True, verbose=False)
    print(result["metric"], [round(v, 4) for v in result["distance_matrix"][0]])
```

- **示例输出：**

```text
unweighted_unifrac [0.0, 0.4, 0.5]
weighted_unifrac [0.0, 0.6157, 0.4352]
```

## 3) `dnakit.microbiome.skbio_p21_functions.rarefy_count_vector()`

- **作用：** 把计数向量随机抽样（稀释）到指定深度。
- **处理方法：** 调用 `skbio.stats.subsample_counts`。默认不放回；`replace=True` 为放回抽样。给定 `seed` 可复现。
- **API：** `rarefy_count_vector(counts, depth, *, replace=False, seed=None)`。
- **输入：** 整数计数列表和目标深度（不放回时不能超过总数）。
- **输出：** 抽样后的 `counts`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import rarefy_count_vector
print(rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)["counts"])
print(rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)["counts"])
print(rarefy_count_vector([1, 1], 6, replace=True, seed=1, verbose=False)["counts"])
```

- **示例输出：**

```text
[2, 1, 0, 1]
[2, 1, 0, 1]
[3, 3]
```

## 4) `dnakit.microbiome.skbio_p21_functions.subsample_items_per_bin()`

- **作用：** 在每个分组中随机保留最多 N 个条目（流式随机抽样）。
- **处理方法：** 调用 `skbio.stats.isubsample`：每个分组不放回地保留至多 `maximum` 个条目，条目数少于 `minimum` 的分组被丢弃。`bins` 给出每个条目所属的分组，省略时所有条目属于同一分组。
- **API：** `subsample_items_per_bin(items, maximum, *, minimum=1, bins=None, seed=None)`。
- **输入：** 字符串条目列表，以及可选的等长分组标签列表。
- **输出：** `selected`，每项为 [分组, 条目]。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import subsample_items_per_bin
result = subsample_items_per_bin(list("abcdefgh"), 2, bins=list("xxxyyyzz"), seed=1, verbose=False)
print(result["selected"])
print(subsample_items_per_bin(list("abc"), 2, minimum=2, bins=["x", "y", "y"], seed=1, verbose=False)["selected"])
```

- **示例输出：**

```text
[['x', 'a'], ['x', 'b'], ['y', 'f'], ['y', 'd'], ['z', 'h'], ['z', 'g']]
[['y', 'c'], ['y', 'b']]
```

## 5) `dnakit.microbiome.skbio_p21_functions.run_hommola_cospeciation_test()`

- **作用：** 宿主-寄生物共物种形成的 Hommola 等（2009）置换检验。
- **处理方法：** 调用 `skbio.stats.evolve.hommola_cospeciation`：比较每条相互作用边上宿主距离与寄生物距离的相关性，通过置换相互作用连接获得 p 值。
- **API：** `run_hommola_cospeciation_test(host_distances, parasite_distances, interaction, *, host_ids=None, parasite_ids=None, permutations=999, seed=None)`。
- **输入：** 宿主距离矩阵、寄生物距离矩阵，以及 0/1 的 寄生物×宿主 相互作用矩阵。
- **输出：** `correlation`、`p_value`、`interaction_edges`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import run_hommola_cospeciation_test
host = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
parasite = [[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 3], [1, 2, 3, 0]]
link = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]]
result = run_hommola_cospeciation_test(host, parasite, link, seed=1, verbose=False)
print(round(result["correlation"], 4), result["p_value"], result["interaction_edges"])
```

- **示例输出：**

```text
0.2335 0.409 4
```

## 6) `dnakit.microbiome.skbio_p21_functions.calculate_gradient_anova()`

- **作用：** 按元数据类别计算排序坐标的梯度 ANOVA 轨迹。
- **处理方法：** `algorithm` 选择 scikit-bio 的类：`average`（AverageGradientANOVA，RMS 平均）、`trajectory`（TrajectoryGradientANOVA，RMS 轨迹）、`first_difference`（FirstDifferenceGradientANOVA）、`window_difference`（WindowDifferenceGradientANOVA，使用 `window_size`）。`axes=0` 使用全部坐标轴。
- **API：** `calculate_gradient_anova(coordinates, proportion_explained, metadata, *, algorithm='average', trajectory_categories=None, sort_category=None, axes=3, weighted=False, window_size=2)`。
- **输入：** 样本到坐标列表的映射、每个坐标轴解释的比例、样本到 {类别: 值} 的元数据映射。
- **输出：** `algorithm`、`weighted` 以及各类别的 `categories`（概率及每个分组的轨迹、均值和附加信息）。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_gradient_anova
coords = {f"s{i}": [0.1 * i, (i * 7 % 5) / 5, (i * 3 % 4) / 4] for i in range(1, 13)}
meta = {k: {"g": "A" if i % 2 else "B", "t": i // 2} for i, k in enumerate(coords)}
for algorithm in ("average", "trajectory", "first_difference", "window_difference"):
    result = calculate_gradient_anova(coords, [0.5, 0.3, 0.2], meta, algorithm=algorithm,
                                      sort_category="t", trajectory_categories=["g"], verbose=False)
    groups = result["categories"][0]["groups"]
    print(result["algorithm"], {k: round(v["mean"], 4) for k, v in groups.items()})
```

- **示例输出：**

```text
avg {'A': 0.1837, 'B': 0.1878}
trajectory {'A': 0.1786, 'B': 0.1786}
diff {'A': 0.0312, 'B': 0.0}
wdiff {'A': 0.0375, 'B': 0.0}
```
