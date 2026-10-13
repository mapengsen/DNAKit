# 系统发育与分子进化

## 1) 调用 IQ-TREE 2 进行最大似然建树

- **作用：** 将已有的多序列比对交给调用者安装的 IQ-TREE 2，生成最大似然树和完整报告文件。
- **计算方法：** 固定调用 `iqtree2 -s ... -m ... -nt ... -pre ...`；可选 `-B` 超快速 bootstrap。DNAKit 不复制或安装 IQ-TREE 2 源码。
- **API：** `dnakit.phylogeny.infer_ml_tree(alignment, output_directory, model="MFP", bootstrap=0, threads=1, ...)`
- **输入：** FASTA、PHYLIP、NEXUS 等 IQ-TREE 2 支持的比对文件；新输出目录；调用者 PATH 中的 `iqtree2`。
- **输出：** JSON 可序列化阶段信封，包含 `.treefile`、`.iqtree`、`.log` 等产物路径，以及模型、bootstrap 和 provenance。
- **引用与许可：** IQ-TREE 2 应引用 Minh et al. (2020), DOI `10.1093/molbev/msaa015`；官方仓库记录 GPL-2.0。若使用 ModelFinder、UFBoot 或其他专门功能，还应按 IQ-TREE 2 文档补充对应论文。

```python
from dnakit.phylogeny import infer_ml_tree

result = infer_ml_tree(
    "alignment.fasta",
    "iqtree-result",
    model="MFP",
    bootstrap=1000,
)
print(result["result"]["treefile"])
```

## 2) 使用 FastTree 进行近似最大似然建树

- **作用：** 将已有比对交给调用者安装的 FastTree，捕获近似最大似然 Newick 树。
- **计算方法：** 固定调用 `FastTree [-nt] [-gtr] [-gamma] alignment`，标准输出写入 `tree.newick`，stderr 单独保留；不复制或安装 FastTree。
- **API：** `dnakit.phylogeny.infer_fast_tree(alignment, output_directory, nucleotide=True, gtr=False, gamma=False, ...)`
- **输入与输出：** 输入是本地比对文件和新输出目录；返回 JSON 可序列化信封，包含 `tree.newick`、stderr、固定参数和 provenance。
- **引用与许可：** 引用 Price, Dehal & Arkin (2010), [10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)；FastTree 2 官方仓库记录 GPL-3.0。

## 1) 从距离矩阵构建进化树 {#bio-6}

- **作用：** 使用 NJ 或 UPGMA 构建带枝长的树。
- **计算方法：** NJ 使用邻接法；UPGMA 按簇大小加权合并。
- **API：** `dnakit.bioinformatics.phylogeny.build_tree(distances, ids, method='nj')`
- **输入：** 对称、对角线为零的距离矩阵及对应 ID；UPGMA 假定分子钟。
- **输出：** 返回 `TreeNode`，示例展示 `tips`。
- **示例代码：**

```python
from dnakit.bioinformatics.phylogeny import build_tree

tree = build_tree([[0, 2, 4], [2, 0, 4], [4, 4, 0]], ["a", "b", "c"], method="upgma")
print(sorted(tree.tips))
```

- **示例结果：**

```text
['a', 'b', 'c']
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

## 2) 比较、合并和裁剪进化树 {#bio-8}

- **作用：** 计算树间 RF 距离、生成共识树或保留指定叶节点。
- **计算方法：** RF 比较分裂集合；共识树保留超过阈值的分支；剪枝合并单子节点。
- **API：** `dnakit.bioinformatics.phylogeny.robinson_foulds(left, right, rooted=False)`; `dnakit.bioinformatics.phylogeny.consensus_tree(trees, threshold=0.5)`; `dnakit.bioinformatics.phylogeny.prune_tree(tree, retain)`
- **输入：** RF/共识使用相同叶集合和一致的定根解释；`threshold=1` 表示严格共识。
- **输出：** `robinson_foulds` 返回距离数值；`consensus_tree` 和 `prune_tree` 返回进化树。
- **示例代码：**

```python
from dnakit.bioinformatics.phylogeny import build_tree, robinson_foulds, consensus_tree, prune_tree

tree = build_tree([[0, 2, 4], [2, 0, 4], [4, 4, 0]], ["a", "b", "c"])
print(robinson_foulds(tree, tree))
print(sorted(consensus_tree([tree, tree]).tips))
print(sorted(prune_tree(tree, ["a", "c"]).tips))
```

- **示例结果：**

```text
0.0
['a', 'b', 'c']
['a', 'c']
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

## 3) dN/dS 选择分析 {#12-dnds}

- **作用：** 对两条密码子对齐 CDS 计算 Nei–Gojobori dN、dS 和 dN/dS，并进行 Jukes–Cantor 修正。
- **处理方法：** 按 Nei–Gojobori 统计同义/非同义位点与差异，并作 Jukes–Cantor 修正。
- **API：** `dnakit.comparative.calculate_dn_ds(sequence_a[必须], sequence_b[必须])`
- **输入：** 等长、线性、canonical A/C/G/T、长度为 3 倍数的 `DNASequence`/`DNARecord`/`DNA`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.comparative import calculate_dn_ds

result = calculate_dn_ds(DNASequence("ATGGCTGAA"), DNASequence("ATGGCCGAG"))
print(result.data["dN"], result.data["dS"], result.data["dN_dS"])
```

- **示例结果：** `dN`、`dS`、`dN_dS`、位点/差异计数和解释；`dS=0` 时比值为 `None`。

[可选后端依赖与返回值说明](../../faq/tutorial_integrations.md#optional-bioinformatics)

<!-- BEGIN GENERATED API: selection -->

<span id="selection"></span>

## 4) 计算两条编码序列的 dN/dS {#dnakit.databases.sequence_dn_ds}

- **作用：** 计算两条编码序列的 dN/dS。
- **处理方法：** 对同框密码子比对使用 Nei–Gojobori 估计及 Jukes–Cantor 校正。
- **API：** `dnakit.databases.sequence_dn_ds(*, seq1=None, seq2=None, fasta1_path=None, fasta2_path=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `dN_dS`, `dN`, `dS`, `pN` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sequence_dn_ds(
    seq1='ATGGCCAAATTTGAAGATCTGAGCACCGTGCGTGGC',
    seq2='ATGGCCAAGTTCGAAAATCTGAGCACCGTACGTGGC',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "dN_dS": "<number>",
  "dN": "<number>",
  "dS": "<number>",
  "pN": "<number>",
  "pS": "<number>",
  "N_sites": "<number>",
  "S_sites": "<number>",
  "Nd": "<number>"
}
```

[完整参数与使用说明](../../faq/api/30_phylogeny.md#dnakit.databases.sequence_dn_ds-parameters)

<!-- END GENERATED API: selection -->
