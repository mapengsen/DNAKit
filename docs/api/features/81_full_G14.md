# DendroPy、toytree 与 pyvolve 系统发育函数（G14）

这一组函数是 `dnakit.phylogeny` 中基于可选依赖的包装，分为两个模块。`dnakit.phylogeny` 用 DendroPy 计算 DNA 比对的群体遗传统计（多态位点、平均两两差异、核苷酸多样性、Tajima's D、Watterson theta、衍生态与频率谱），以及两棵树的比较和树形统计。`dnakit.phylogeny` 用 toytree 完成定根、编辑、距离和类群矩阵，只返回 Newick 或数值，不绘图。

pyvolve 没有新增包装：它的 `read_tree` 已被现有的 `simulate_nucleotide_*` 函数使用，`print_tree` 只是文本渲染。

所有函数只处理显式给出的序列或树，不访问网络，也不调用外部程序。依赖是可选的，只有调用对应函数时才会导入；缺少包时会报错并给出 `dnakit[bio-dendropy]` 或 `dnakit[bio-toytree]` 的安装命令。这些函数目前需要从模块路径导入，尚未从 `dnakit.phylogeny` 统一导出。

## 示例数据

下面的示例大多把比对或树直接写在代码中。`tree_t1.nwk` 只用于演示文件路径输入。

**示例文件：** `tree_t1.nwk`

```text
((A:1,B:1):1,(C:2,D:2):0.5);
```

## 1) `dnakit.phylogeny.count_segregating_sites()`

- **作用：** 统计 DNA 比对中的多态位点数。
- **处理方法：** 调用 `dendropy.calculate.popgenstat.num_segregating_sites`。位点在任一序列上含歧义或缺失符号时，`ignore_uncertain=True` 按 dendropy 规则跳过该位点。
- **API：** `dnakit.phylogeny.count_segregating_sites(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 序列 ID 到序列的字典，或 `DNA` 批次；至少两条等长序列。`max_sequences` 与 `max_length` 是显式上限。
- **输出：** `result` 含 `segregating_sites`、`sequence_count`、`site_count` 和 `ignore_uncertain`；`provenance` 记录 dendropy 版本。
- **示例代码：**

```python
from dnakit.phylogeny import count_segregating_sites

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(count_segregating_sites(alignment)["result"])
```

- **示例输出：**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'segregating_sites': 3}
```

## 2) `dnakit.phylogeny.calculate_mean_pairwise_differences()`

- **作用：** 计算比对的平均两两差异数（Tajima 的 k）。
- **处理方法：** 调用 `dendropy.calculate.popgenstat.average_number_of_pairwise_differences`。对每一对序列计数不同位点，再对所有序列对取平均。
- **API：** `dnakit.phylogeny.calculate_mean_pairwise_differences(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `count_segregating_sites()` 相同的序列字典或 `DNA` 批次。
- **输出：** `result.mean_pairwise_differences` 为平均值；其余字段与 `count_segregating_sites()` 相同。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_mean_pairwise_differences

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_mean_pairwise_differences(alignment)["result"])
```

- **示例输出：**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'mean_pairwise_differences': 1.5}
```

## 3) `dnakit.phylogeny.calculate_nucleotide_diversity()`

- **作用：** 计算比对的核苷酸多样性 pi。
- **处理方法：** 调用 `dendropy.calculate.popgenstat.nucleotide_diversity`。数值由 dendropy 计算，本包装不再重算。
- **API：** `dnakit.phylogeny.calculate_nucleotide_diversity(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `count_segregating_sites()` 相同。
- **输出：** `result.nucleotide_diversity` 为 pi；其余字段相同。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_nucleotide_diversity

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_nucleotide_diversity(alignment)["result"])
```

- **示例输出：**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'nucleotide_diversity': 0.15}
```

## 4) `dnakit.phylogeny.calculate_tajimas_d()`

- **作用：** 计算比对的 Tajima's D。
- **处理方法：** 先用 dendropy 统计多态位点数 S；S 为 0 时 D 无定义，返回 `None`。否则调用 `dendropy.calculate.popgenstat.tajimas_d`。
- **API：** `dnakit.phylogeny.calculate_tajimas_d(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `count_segregating_sites()` 相同。
- **输出：** `result.tajimas_d` 为 D 或 `None`；`result.segregating_sites` 为 S。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tajimas_d

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_tajimas_d(alignment)["result"])
```

- **示例输出：**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'tajimas_d': -0.754451077653, 'segregating_sites': 3}
```

## 5) `dnakit.phylogeny.calculate_wattersons_theta()`

- **作用：** 计算比对的 Watterson theta。
- **处理方法：** 调用 `dendropy.calculate.popgenstat.wattersons_theta`，结果为 S / a_n，其中 a_n = 1 + 1/2 + ... + 1/(n-1)，不再除以位点数。
- **API：** `dnakit.phylogeny.calculate_wattersons_theta(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `count_segregating_sites()` 相同。
- **输出：** `result.wattersons_theta` 为 theta_W；`result.convention_note` 说明分母。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_wattersons_theta

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_wattersons_theta(alignment)["result"])
```

- **示例输出：**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'wattersons_theta': 1.636363636364, 'convention_note': 'segregating sites divided by a_n = sum(1/i) for i in 1..n-1'}
```

## 6) `dnakit.phylogeny.build_derived_state_matrix()`

- **作用：** 按祖先序列为每个类群位点标记 0（祖先态）或 1（衍生态）。
- **处理方法：** 把 `ancestral_id` 指定的序列作为祖先参考，并从类群中去掉它；调用 `dendropy.calculate.popgenstat.derived_state_matrix`，每行按 ID 排序。
- **API：** `dnakit.phylogeny.build_derived_state_matrix(sequences, *, ancestral_id, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `count_segregating_sites()` 相同，另加 `ancestral_id`，必须是输入中的某个 ID。
- **输出：** `result.rows` 为 `{id, states}` 列表，`states` 是 0/1 字符串；`result.ancestral_id` 回显祖先 ID。
- **示例代码：**

```python
from dnakit.phylogeny import build_derived_state_matrix

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(build_derived_state_matrix(alignment, ancestral_id="A")["result"]["rows"])
```

- **示例输出：**

```text
[{'id': 'B', 'states': '0000000100'}, {'id': 'C', 'states': '0000100100'}, {'id': 'D', 'states': '0001000100'}]
```

## 7) `dnakit.phylogeny.calculate_unfolded_site_frequency_spectrum()`

- **作用：** 计算类群相对祖先序列的未折叠位点频率谱。
- **处理方法：** 调用 `dendropy.calculate.popgenstat.unfolded_site_frequency_spectrum`（`pad=True`），输出长度为类群数加一的列表，第 k 项为恰有 k 个衍生等位基因的位点数。默认 `ignore_uncertain=False`。
- **API：** `dnakit.phylogeny.calculate_unfolded_site_frequency_spectrum(sequences, *, ancestral_id, ignore_uncertain=False, max_sequences=1000, max_length=100000, verbose=True)`。
- **输入：** 与 `build_derived_state_matrix()` 相同。
- **输出：** `result.site_frequency_spectrum` 为列表；`result.ingroup_count` 为类群序列数。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_unfolded_site_frequency_spectrum

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
result = calculate_unfolded_site_frequency_spectrum(alignment, ancestral_id="A")
print(result["result"]["site_frequency_spectrum"])
```

- **示例输出：**

```text
[7, 2, 0, 1]
```

## 8) `dnakit.phylogeny.calculate_weighted_robinson_foulds_distance()`

- **作用：** 计算两棵树（同一组类群）之间的加权 Robinson-Foulds 距离。
- **处理方法：** 两棵树共用一个 TaxonNamespace 读入；类群标签必须完全一致，且每条非根分支都要有长度。调用 `dendropy.calculate.treecompare.weighted_robinson_foulds_distance`。
- **API：** `dnakit.phylogeny.calculate_weighted_robinson_foulds_distance(tree_a, tree_b, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两个 Newick 字符串或单树 Newick 文件路径。`max_tips` 与 `max_file_bytes` 是显式上限。
- **输出：** `result.distance` 为距离，`result.tip_count` 为类群数。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_weighted_robinson_foulds_distance

print(calculate_weighted_robinson_foulds_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 5.0}
```

## 9) `dnakit.phylogeny.count_false_positive_and_negative_splits()`

- **作用：** 统计两棵树中只出现在一棵树里的二分（split）数。
- **处理方法：** 调用 `dendropy.calculate.treecompare.false_positives_and_negatives`。false_positives 是只在比较树中出现的二分，false_negatives 是只在参考树中出现的二分。
- **API：** `dnakit.phylogeny.count_false_positive_and_negative_splits(reference_tree, comparison_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 参考树与比较树，格式同 `calculate_weighted_robinson_foulds_distance()`。
- **输出：** `result.false_positives`、`result.false_negatives` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import count_false_positive_and_negative_splits

print(count_false_positive_and_negative_splits("((A,B),(C,D));", "((A,C),(B,D));")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'false_positives': 1, 'false_negatives': 1}
```

## 10) `dnakit.phylogeny.calculate_colless_tree_imbalance()`

- **作用：** 计算二叉树的 Colless 不平衡指数。
- **处理方法：** 调用 `dendropy.calculate.treemeasure.colless_tree_imbalance`。非二叉树会报 `ConfigurationError`。`normalize` 可取 `"max"`（默认）、`"yule"`、`"pda"` 或 `None`（不归一化）。
- **API：** `dnakit.phylogeny.calculate_colless_tree_imbalance(tree, *, normalize="max", max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.colless_imbalance`、`result.normalize` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_colless_tree_imbalance

print(calculate_colless_tree_imbalance("((A,B),(C,(D,E)));", normalize=None)["result"])
```

- **示例输出：**

```text
{'tip_count': 5, 'colless_imbalance': 2.0, 'normalize': 'none'}
```

## 11) `dnakit.phylogeny.calculate_sackin_index()`

- **作用：** 计算树的 Sackin 指数，即每个类群的祖先数之和。
- **处理方法：** 调用 `dendropy.calculate.treemeasure.sackin_index`。`normalize="tips"`（默认）除以类群数，等于 N-bar；`"yule"` 与 `"pda"` 按对应模型归一化；`None` 不归一化。
- **API：** `dnakit.phylogeny.calculate_sackin_index(tree, *, normalize="tips", max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.sackin_index`、`result.normalize` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_sackin_index

print(calculate_sackin_index("((A,B),(C,D));", normalize=None)["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'sackin_index': 8.0, 'normalize': 'none'}
```

## 12) `dnakit.phylogeny.calculate_b1_index()`

- **作用：** 计算 B1 树形统计量。
- **处理方法：** 调用 `dendropy.calculate.treemeasure.B1`。对每个非根内部节点，取 1/m 求和，其中 m 比其子节点中最大的 m 多 1，叶子的 m 为 0。
- **API：** `dnakit.phylogeny.calculate_b1_index(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.b1` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_b1_index

print(calculate_b1_index("((A,B),(C,(D,E)));")["result"])
```

- **示例输出：**

```text
{'tip_count': 5, 'b1': 2.5}
```

## 13) `dnakit.phylogeny.calculate_pybus_harvey_gamma()`

- **作用：** 计算 Pybus 和 Harvey 的 gamma 统计量，用于检验速率是否恒定。
- **处理方法：** 调用 `dendropy.calculate.treemeasure.pybus_harvey_gamma`。树必须是超度量的（所有叶子到根的距离在 `precision` 内相等），否则报 `ConfigurationError`。
- **API：** `dnakit.phylogeny.calculate_pybus_harvey_gamma(tree, *, precision=1e-05, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`precision` 取 (0, 1]。
- **输出：** `result.gamma` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_pybus_harvey_gamma

print(calculate_pybus_harvey_gamma("((A:1,B:1):1,(C:1,D:1):1);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'gamma': -0.816496580928}
```

## 14) `dnakit.phylogeny.calculate_treeness()`

- **作用：** 计算树的 treeness，即内部分支长度占总分支长度的比例。
- **处理方法：** 调用 `dendropy.calculate.treemeasure.treeness`。所有分支都必须有长度。
- **API：** `dnakit.phylogeny.calculate_treeness(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.treeness` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_treeness

print(calculate_treeness("((A:1,B:1):1,(C:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'treeness': 0.2}
```

## 15) `dnakit.phylogeny.calculate_parsimony_score()`

- **作用：** 计算 DNA 字符在树上的 Fitch 简约分值。
- **处理方法：** 把序列按 ID 对齐到树的类群标签上，调用 `dendropy.calculate.treescore.parsimony_score`。序列 ID 必须与树的叶子标签完全一致。
- **API：** `dnakit.phylogeny.calculate_parsimony_score(tree, sequences, *, gaps_as_missing=True, max_sequences=1000, max_length=100000, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径，以及与前述相同的序列字典或 `DNA` 批次。`gaps_as_missing` 决定缺口是否视为缺失数据。
- **输出：** `result.parsimony_score` 为整数，`result.tip_count` 与 `result.site_count` 为类群数与位点数。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_parsimony_score

chars = {"A": "AAAC", "B": "AAGC", "C": "TAGC", "D": "TCGC"}
print(calculate_parsimony_score("((A,B),(C,D));", chars)["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'site_count': 4, 'parsimony_score': 3, 'gaps_as_missing': True}
```

## 16) `dnakit.phylogeny.calculate_reconciliation_discordance()`

- **作用：** 计算基因树在物种树上调和（reconcile）时隐含的基因复制数。
- **处理方法：** 两棵树都按有根树读入（`rooting="force-rooted"`），共享类群集合必须相同，调用 `dendropy.model.reconcile.reconciliation_discordance`。
- **API：** `dnakit.phylogeny.calculate_reconciliation_discordance(gene_tree, species_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 基因树和物种树，均为 Newick 字符串或单树文件路径。
- **输出：** `result.duplications` 为复制数，`result.tip_count` 为类群数。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_reconciliation_discordance

print(calculate_reconciliation_discordance("((A,C),(B,(D,E)));", "((A,B),(C,(D,E)));")["result"])
```

- **示例输出：**

```text
{'tip_count': 5, 'duplications': 2}
```

## 17) `dnakit.phylogeny.calculate_monophyletic_partition_discordance()`

- **作用：** 计算若干类群组各自强制单系时需要的深度合并（deep coalescence）数。
- **处理方法：** `groups` 必须恰好划分树的全部类群（每个类群出现一次）。调用 `dendropy.model.reconcile.monophyletic_partition_discordance`。
- **API：** `dnakit.phylogeny.calculate_monophyletic_partition_discordance(tree, groups, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径，以及类群标签列表的列表。
- **输出：** `result.deep_coalescences`、`result.group_count` 和 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_monophyletic_partition_discordance

groups = [["A", "B"], ["C"], ["D", "E"]]
print(calculate_monophyletic_partition_discordance("((A,C),(B,(D,E)));", groups)["result"])
```

- **示例输出：**

```text
{'tip_count': 5, 'deep_coalescences': 1, 'group_count': 3}
```

## 18) `dnakit.phylogeny.root_tree_on_midpoint()`

- **作用：** 把树根放在最长两两路径的中点，返回 Newick。
- **处理方法：** 调用 `toytree.mod.root_on_midpoint`。先找最长路径，再把根放在其中点。
- **API：** `dnakit.phylogeny.root_tree_on_midpoint(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。`max_tips` 与 `max_file_bytes` 是显式上限。
- **输出：** `result.newick` 为有根 Newick，`result.tips` 为排序后的类群标签。
- **示例代码：**

```python
from dnakit.phylogeny import root_tree_on_midpoint

print(root_tree_on_midpoint("(A:1,B:1,(C:2,D:2):1.5);")["result"]["newick"])
```

- **示例输出：**

```text
((C:2,D:2):0.25,(A:1,B:1):1.25);
```

## 19) `dnakit.phylogeny.root_tree_on_balanced_midpoint()`

- **作用：** 按平衡中点法为树定根，使两侧深度大致相等。
- **处理方法：** 调用 `toytree.mod.root_on_balanced_midpoint`，`tolerance` 控制收敛精度。
- **API：** `dnakit.phylogeny.root_tree_on_balanced_midpoint(tree, *, tolerance=1e-06, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`tolerance` 为正数。
- **输出：** `result.newick`、`result.tolerance` 和类群信息。
- **示例代码：**

```python
from dnakit.phylogeny import root_tree_on_balanced_midpoint

print(root_tree_on_balanced_midpoint("(A:1,B:1,(C:2,D:2):1.5);", tolerance=1e-9)["result"]["newick"])
```

- **示例输出：**

```text
((C:2,D:2):0.249999999069,(A:1,B:1):1.25000000093);
```

## 20) `dnakit.phylogeny.root_tree_on_minimal_ancestor_deviation()`

- **作用：** 按最小祖先偏差（MAD）法为树定根。
- **处理方法：** 调用 `toytree.mod.root_on_minimal_ancestor_deviation`，选择使祖先偏差最小的边作为根的位置。
- **API：** `dnakit.phylogeny.root_tree_on_minimal_ancestor_deviation(tree, *, min_dist=1e-12, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`min_dist` 为正数。
- **输出：** `result.newick` 与 `result.min_dist`。
- **示例代码：**

```python
from dnakit.phylogeny import root_tree_on_minimal_ancestor_deviation

print(root_tree_on_minimal_ancestor_deviation("(A:1,B:1,(C:2,D:2):1.5);")["result"]["newick"])
```

- **示例输出：**

```text
((C:2,D:2):0.25,(A:1,B:1):1.25);
```

## 21) `dnakit.phylogeny.unroot_tree()`

- **作用：** 去掉根的分裂，返回无根树的 Newick。
- **处理方法：** 调用 `toytree.mod.unroot`，把根两侧的两条边合并为一条。
- **API：** `dnakit.phylogeny.unroot_tree(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.newick` 为无根 Newick。
- **示例代码：**

```python
from dnakit.phylogeny import unroot_tree

print(unroot_tree("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["newick"])
```

- **示例输出：**

```text
(A:1,B:1,(C:2,D:2):1.5);
```

## 22) `dnakit.phylogeny.ladderize_tree()`

- **作用：** 按子树大小排序各节点的子节点，返回 Newick。
- **处理方法：** 调用 `toytree.mod.ladderize`。`direction=True` 时反转排序方向。
- **API：** `dnakit.phylogeny.ladderize_tree(tree, *, direction=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`direction` 为布尔值。
- **输出：** `result.newick` 与 `result.direction`。
- **示例代码：**

```python
from dnakit.phylogeny import ladderize_tree

print(ladderize_tree("((A:1,B:1):1,(C:2,(D:2,E:2):1):0.5);", direction=True)["result"]["newick"])
```

- **示例输出：**

```text
(((D:2,E:2):1,C:2):0.5,(A:1,B:1):1);
```

## 23) `dnakit.phylogeny.extract_subtree()`

- **作用：** 提取由给定类群张成的子树（进化枝），返回 Newick。
- **处理方法：** 调用 `toytree.mod.extract_subtree`，子树根保留原有的茎长度。类群必须存在于树中。
- **API：** `dnakit.phylogeny.extract_subtree(tree, tips, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`tips` 为至少两个不重复的类群标签。
- **输出：** `result.newick`、`result.tips` 与 `result.requested_tips`。
- **示例代码：**

```python
from dnakit.phylogeny import extract_subtree

print(extract_subtree("((A:1,B:1):1,(C:2,D:2):0.5);", ["C", "D"])["result"]["newick"])
```

- **示例输出：**

```text
(C:2,D:2):0.5;
```

## 24) `dnakit.phylogeny.collapse_weak_nodes()`

- **作用：** 折叠支持度低于阈值的内部节点，返回 Newick。
- **处理方法：** 调用 `toytree.mod.collapse_nodes(min_support=...)`，支持度取自 Newick 中内部节点的标签。
- **API：** `dnakit.phylogeny.collapse_weak_nodes(tree, *, min_support, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`min_support` 为非负数。
- **输出：** `result.newick` 与 `result.min_support`。
- **示例代码：**

```python
from dnakit.phylogeny import collapse_weak_nodes

print(collapse_weak_nodes("((A:1,B:1)90:1,(C:2,D:2)30:0.5);", min_support=50)["result"]["newick"])
```

- **示例输出：**

```text
((A:1,B:1)90:1,C:2.5,D:2.5);
```

## 25) `dnakit.phylogeny.remove_unary_nodes()`

- **作用：** 去掉只有一个子节点的内部节点，并把分支长度合并。
- **处理方法：** 调用 `toytree.mod.remove_unary_nodes`。
- **API：** `dnakit.phylogeny.remove_unary_nodes(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.newick` 为去掉单子节点后的 Newick。
- **示例代码：**

```python
from dnakit.phylogeny import remove_unary_nodes

print(remove_unary_nodes("((A:1):1,(B:2,C:2):1);")["result"]["newick"])
```

- **示例输出：**

```text
((B:2,C:2):1,A:2);
```

## 26) `dnakit.phylogeny.scale_tree_to_root_height()`

- **作用：** 缩放所有分支长度，使根到叶的最大高度等于指定值。
- **处理方法：** 调用 `toytree.mod.edges_scale_to_root_height`。`include_stem=True` 时把根的茎边也计入缩放。
- **API：** `dnakit.phylogeny.scale_tree_to_root_height(tree, *, height=1.0, include_stem=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`height` 为正数；`include_stem` 为布尔值。
- **输出：** `result.newick`、`result.height` 与 `result.include_stem`。
- **示例代码：**

```python
from dnakit.phylogeny import scale_tree_to_root_height

print(scale_tree_to_root_height("((A:1,B:1):1,(C:2,D:2):0.5);", height=2.0)["result"]["newick"])
```

- **示例输出：**

```text
((A:0.8,B:0.8):0.8,(C:1.6,D:1.6):0.4);
```

## 27) `dnakit.phylogeny.calculate_tree_vcv_matrix()`

- **作用：** 计算类群的方差-协方差（VCV）矩阵。
- **处理方法：** 调用 `toytree.pcm.get_vcv_matrix_from_tree`。对角线是根到类群的距离，非对角线是两类群共享的从根出发路径长度。行和列按类群标签排序。
- **API：** `dnakit.phylogeny.calculate_tree_vcv_matrix(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 示例文件 `tree_t1.nwk` 或 Newick 字符串。
- **输出：** `result.tips` 与 `result.matrix`（嵌套列表）。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tree_vcv_matrix

print(calculate_tree_vcv_matrix("tree_t1.nwk")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'tips': ['A', 'B', 'C', 'D'], 'matrix': [[2.0, 1.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.5, 0.5], [0.0, 0.0, 0.5, 2.5]]}
```

## 28) `dnakit.phylogeny.calculate_tree_correlation_matrix()`

- **作用：** 计算类群的相关矩阵。
- **处理方法：** 调用 `toytree.pcm.get_corr_matrix_from_tree`，即 VCV 矩阵的每项除以两个类群方差的几何平均。
- **API：** `dnakit.phylogeny.calculate_tree_correlation_matrix(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.tips` 与 `result.matrix`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_tree_correlation_matrix

print(calculate_tree_correlation_matrix("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["matrix"])
```

- **示例输出：**

```text
[[1.0, 0.5, 0.0, 0.0], [0.5, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.2], [0.0, 0.0, 0.2, 1.0]]
```

## 29) `dnakit.phylogeny.calculate_information_rf_distance()`

- **作用：** 计算信息校正的 Robinson-Foulds 距离（rfi）。
- **处理方法：** 调用 `toytree.distance.get_treedist_rfi`。
- **API：** `dnakit.phylogeny.calculate_information_rf_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵树（Newick 字符串或文件路径），类群标签必须相同；`normalize` 默认为 `False`。
- **输出：** `result.distance`、`result.normalize` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_information_rf_distance

print(calculate_information_rf_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 3.169925001442, 'normalize': False}
```

## 30) `dnakit.phylogeny.calculate_matching_split_distance()`

- **作用：** 计算广义 Robinson-Foulds 的匹配二分距离（rfg_ms）。
- **处理方法：** 调用 `toytree.distance.get_treedist_rfg_ms`。toytree 对该度量不提供归一化。
- **API：** `dnakit.phylogeny.calculate_matching_split_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵树（Newick 字符串或文件路径），类群标签必须相同；`normalize` 默认为 `False`。
- **输出：** `result.distance`、`result.normalize` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_matching_split_distance

print(calculate_matching_split_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 31) `dnakit.phylogeny.calculate_matching_split_information_distance()`

- **作用：** 计算匹配二分的信息距离（rfg_msi）。
- **处理方法：** 调用 `toytree.distance.get_treedist_rfg_msi`。
- **API：** `dnakit.phylogeny.calculate_matching_split_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵树（Newick 字符串或文件路径），类群标签必须相同；`normalize` 默认为 `False`。
- **输出：** `result.distance`、`result.normalize` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_matching_split_information_distance

print(calculate_matching_split_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 32) `dnakit.phylogeny.calculate_shared_phylogenetic_information_distance()`

- **作用：** 计算广义 Robinson-Foulds 的共享系统发育信息距离（rfg_spi）。
- **处理方法：** 调用 `toytree.distance.get_treedist_rfg_spi`。两棵树含有完全相同的二分集合时，toytree 的相似度无定义，本包装报 `ConfigurationError`。
- **API：** `dnakit.phylogeny.calculate_shared_phylogenetic_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵树（Newick 字符串或文件路径），类群标签必须相同；`normalize` 默认为 `False`。
- **输出：** `result.distance`、`result.normalize` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_shared_phylogenetic_information_distance

print(calculate_shared_phylogenetic_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 3.169925001442, 'normalize': False}
```

## 33) `dnakit.phylogeny.calculate_mutual_clustering_information_distance()`

- **作用：** 计算广义 Robinson-Foulds 的互聚类信息距离（rfg_mci）。
- **处理方法：** 调用 `toytree.distance.get_treedist_rfg_mci`。
- **API：** `dnakit.phylogeny.calculate_mutual_clustering_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** 两棵树（Newick 字符串或文件路径），类群标签必须相同；`normalize` 默认为 `False`。
- **输出：** `result.distance`、`result.normalize` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_mutual_clustering_information_distance

print(calculate_mutual_clustering_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **示例输出：**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 34) `dnakit.phylogeny.list_tree_bipartitions()`

- **作用：** 列出树的二分（split），每个二分写成两侧排序后的类群列表。
- **处理方法：** 调用 `toytree.enum.iter_bipartitions`。每个二分中，首个类群较小的一侧排在前面；二分按字典序排序。`include_singletons=True` 时包含只有一个类群的二分。
- **API：** `dnakit.phylogeny.list_tree_bipartitions(tree, *, include_singletons=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径；`include_singletons` 为布尔值。
- **输出：** `result.splits`（二分列表）、`result.split_count` 与 `result.tip_count`。
- **示例代码：**

```python
from dnakit.phylogeny import list_tree_bipartitions

print(list_tree_bipartitions("((A:1,B:1):1,(C:2,D:2):0.5);", include_singletons=True)["result"]["splits"])
```

- **示例输出：**

```text
[[['A'], ['B', 'C', 'D']], [['A', 'B'], ['C', 'D']], [['A', 'B', 'C'], ['D']], [['A', 'B', 'D'], ['C']], [['A', 'C', 'D'], ['B']]]
```

## 35) `dnakit.phylogeny.calculate_relative_evolutionary_divergence()`

- **作用：** 计算每个节点的相对进化分歧（RED）。
- **处理方法：** 调用 `toytree.pcm.get_relative_evolutionary_divergence`，按节点索引排序输出。
- **API：** `dnakit.phylogeny.calculate_relative_evolutionary_divergence(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.nodes` 为列表，每项含 `node_index`、`name`（内部节点为 null）和 `red`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_relative_evolutionary_divergence

print(calculate_relative_evolutionary_divergence("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["nodes"])
```

- **示例输出：**

```text
[{'node_index': 0, 'name': 'A', 'red': 1.0}, {'node_index': 1, 'name': 'B', 'red': 1.0}, {'node_index': 2, 'name': 'C', 'red': 1.0}, {'node_index': 3, 'name': 'D', 'red': 1.0}, {'node_index': 5, 'name': None, 'red': 0.2}, {'node_index': 6, 'name': None, 'red': 0.0}]
```

## 36) `dnakit.phylogeny.calculate_equal_splits()`

- **作用：** 计算每个类群的等分（equal splits，ES）指标。
- **处理方法：** 调用 `toytree.pcm.get_equal_splits`，只保留类群（忽略根），按标签排序。
- **API：** `dnakit.phylogeny.calculate_equal_splits(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`。
- **输入：** Newick 字符串或单树文件路径。
- **输出：** `result.tips` 为列表，每项含 `tip` 与 `equal_splits`。
- **示例代码：**

```python
from dnakit.phylogeny import calculate_equal_splits

print(calculate_equal_splits("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["tips"])
```

- **示例输出：**

```text
[{'tip': 'A', 'equal_splits': 0.5}, {'tip': 'B', 'equal_splits': 0.5}, {'tip': 'C', 'equal_splits': 0.25}, {'tip': 'D', 'equal_splits': 0.25}]
```

## 当前范围

- pyvolve 没有新增包装。
- toytree 的 `calculate_shared_phylogenetic_information_distance()` 在两棵树的二分集合完全相同时报 `ConfigurationError`，因为 toytree 在该情形下的相似度无定义。
- 本文档只覆盖上述 36 个包装函数，不包括同一上游函数的别名、重复的现有功能和超出上限的函数。
