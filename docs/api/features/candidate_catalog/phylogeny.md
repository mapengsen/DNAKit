# 系统发育推断、检验与绘图

本类收录 **8** 项有文档的功能；这 8 项都直接由输入（Newick 树、距离矩阵或 DNA 比对）用纯 Python 在本地计算出结果，不需要外部程序或可选库。只有 `infer` 可能返回 `status: "not_implemented"`，且仅当 `options.method` 要求最大似然或贝叶斯方法时（这类方法请使用会调用外部程序的 `infer_fast_tree` / `infer_ml_tree`）。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) annotate

- **作用：** 为 Newick 树的每个节点标注结构信息（边数深度、叶子数、到根的距离、所属分支的叶子集合），并可附加每个叶子的元数据。

- **处理方法：** 用纯 Python 解析 Newick 树，按前序为节点编号，并把没有名字的内部节点命名为 `N<前序编号>`。对每个节点给出父节点编号、子节点数、其下的叶子数、边数深度、自身分支长度和到根的距离（从根起分支长度之和；路径上有分支缺少长度时为 `None`）。叶子数不超过 `max_clade` 的内部节点还列出其叶子名称（已排序）。给出 `metadata` 时，每个节点列出其下叶子的不同元数据取值（转为字符串）以及这些取值是否一致（恰好只有一个取值）。不沿树做任何推断（不做祖先状态重建）。

- **API：** `dnakit.phylogeny.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：Newick 文本或 Newick 文件路径。`options`（可选）：`metadata`（叶子名到取值的映射，默认无；不是叶子的名称会报错）、`max_clade`（整数 1-1,000,000，默认 50；叶子数更多的内部节点不列出分支叶子）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`nodes`（按前序每个节点一行，含 `id`、`name`、`parent`、`is_leaf`、`child_count`、`leaf_count`、`depth_edges`、`branch_length`、`root_distance`，较小的内部节点含 `clade`，给出元数据时含 `metadata_values` 和 `metadata_homogeneous`）、`annotated_newick`（带生成标签的树）和 `internal_nodes_labelled`（生成的标签个数）。

- **示例代码：**

```python
import json

from dnakit.phylogeny import annotate

result = annotate(
    "((A:1,B:2):3,C:4);",
    options={"metadata": {"A": "x", "B": "x", "C": "y"}},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a96a137a9d8a4755cd277cb9193377aab112cd1f3719f4065ede30482e9bbd5b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "nodes": [
      {
        "id": 0,
        "name": "N0",
        "parent": null,
        "is_leaf": false,
        "child_count": 2,
        "leaf_count": 3,
        "depth_edges": 0,
        "branch_length": null,
        "root_distance": 0.0,
        "clade": ["A", "B", "C"],
        "metadata_values": ["x", "y"],
        "metadata_homogeneous": false
      },
      {
        "id": 1,
        "name": "N1",
        "parent": 0,
        "is_leaf": false,
        "child_count": 2,
        "leaf_count": 2,
        "depth_edges": 1,
        "branch_length": 3.0,
        "root_distance": 3.0,
        "clade": ["A", "B"],
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 2,
        "name": "A",
        "parent": 1,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 2,
        "branch_length": 1.0,
        "root_distance": 4.0,
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 3,
        "name": "B",
        "parent": 1,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 2,
        "branch_length": 2.0,
        "root_distance": 5.0,
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 4,
        "name": "C",
        "parent": 0,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 1,
        "branch_length": 4.0,
        "root_distance": 4.0,
        "metadata_values": ["y"],
        "metadata_homogeneous": true
      }
    ],
    "annotated_newick": "((A:1,B:2)N1:3,C:4)N0;",
    "internal_nodes_labelled": 2
  }
}
```

## 2) build

- **作用：** 根据成对距离矩阵用邻接法（neighbour-joining，Saitou 与 Nei 1987）或 UPGMA（平均连接）构建树。

- **处理方法：** 邻接法每次合并使 Q 判据 `(m-2)*d(i,j) - sum_k d(i,k) - sum_k d(j,k)` 最小的一对，由净分歧度计算两条分支长度，并把到新节点的距离更新为 `(d(i,k) + d(j,k) - d(i,j))/2`；最后三个节点在一个三分叉处合并。结果是无根树（以三分叉根的形式书写），分支长度可能为负，且不会被截断。UPGMA 每次合并最近的两个簇，新节点高度取二者距离的一半，到其他簇的距离按簇大小加权平均，因此得到有根的超度量树。距离相同时按输入顺序决定。矩阵必须是方阵、对称、对角线为零、元素为非负有限数，名称唯一，且最多 400 个类群。纯 Python 实现，不调用外部程序。

- **API：** `dnakit.phylogeny.build(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"names": [...], "matrix": [[...]]}`（缺少 `names` 时类群记为 t1、t2……）或嵌套映射 `{a: {b: 距离}}`。`options`（可选）：`method`（`"nj"` 为默认，或 `"upgma"`）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`newick`（树）、`method`、`rooted`（仅 UPGMA 为 true）、`leaf_count`、`total_branch_length`、`negative_branch_count` 以及说明文字 `note`。

- **示例代码：**

```python
import json

from dnakit.phylogeny import build

matrix = {
    "names": ["a", "b", "c", "d", "e"],
    "matrix": [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0],
    ],
}
result = build(matrix)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "build",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 137,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f8d804d5c145007c60bd2e00c2dd6c386da3e2bdccb29b0fececaf0e04ada64f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick": "(d:2,e:1,(c:4,(a:2,b:3):3):2);",
    "method": "nj",
    "rooted": false,
    "leaf_count": 5,
    "total_branch_length": 17.0,
    "negative_branch_count": 0,
    "note": "UPGMA gives a rooted ultrametric tree; neighbour-joining gives an unrooted tree drawn with a trifurcating root and may contain negative branch lengths"
  }
}
```

## 3) calculate

- **作用：** 从 Newick 树计算树的统计量或系统发育（patristic）距离矩阵，或从已比对的 DNA 序列计算成对距离矩阵。

- **处理方法：** `what="tree_statistics"` 给出叶子数、内部节点数和节点总数、叶子名称、根的度、是否为有根二叉树、多分支节点数（三分叉的根不计入）、叶子边数深度的最大值和平均值、Sackin 指数（叶子边数深度之和）、分支总长度、根到叶子距离的最小值、最大值和平均值、是否为超度量树（叶子距离在相对 1e-9 内相等）、treeness（内部分支长度 / 分支总长度），以及仅对有根二叉树给出的 Colless 指数（各内部节点两个子树叶子数之差的绝对值之和）。`what="patristic_matrix"` 返回经最近共同祖先的叶子间路径长度（所有分支都需要长度；最多 1,000 个叶子）。`what="distance_matrix"` 对每对序列计算 p 距离（错配数 / 参与比较的位点数）或 Jukes-Cantor 距离 `d = -3/4 ln(1 - 4p/3)`；只比较 A/C/G/T，因此空位和其他符号在成对比较时被跳过（U 视为 T，忽略大小写）。饱和的序列对（p >= 0.75）或没有可比较位点的序列对会报错，除非设置了 `saturation_distance`。

- **API：** `dnakit.phylogeny.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：计算树结果时为 Newick 树（文本或路径）；计算距离矩阵时为比对（FASTA 文本/路径、`{名称: 序列}` 或 `[{"name", "sequence"}]`，长度相等，最多 400 条序列）。`options`（可选）：`what`（默认 `"auto"`：Newick 字符串给出树统计量，比对给出距离矩阵；也可取 `"tree_statistics"`、`"patristic_matrix"`、`"distance_matrix"`）、`distance`（默认 `"jc69"`，或 `"p"`）、`saturation_distance`（为饱和或无法定义的序列对代替报错的正数，默认无）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：始终包含 `what`。`tree_statistics`：上述统计量（`leaf_count`、`internal_node_count`、`node_count`、`leaves`、`root_degree`、`rooted_binary`、`polytomy_count`、`max_depth_edges`、`mean_leaf_depth_edges`、`sackin_index`、`has_branch_lengths`、`all_branches_have_lengths`、`negative_branch_count`、`total_branch_length`；所有分支都有长度时还有 `max_root_to_tip`、`min_root_to_tip`、`mean_root_to_tip`、`ultrametric`、`treeness`；有根二叉树另有 `colless_index`）。`patristic_matrix`：`names`、`matrix`、`distance`（`"patristic"`）。`distance_matrix`：`names`、`matrix`、`distance`、`compared_sites`（参与比较的位点数矩阵）和 `alignment_length`。

- **示例代码：**

```python
import json

from dnakit.phylogeny import calculate

result = calculate("((A:1,B:2)X:3,C:4);")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 19,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "90a78a93da03a005ce42758202cb84a0e7b01040de38aa3be6ccff6b31fc785c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "what": "tree_statistics",
    "leaf_count": 3,
    "internal_node_count": 2,
    "node_count": 5,
    "leaves": ["A", "B", "C"],
    "root_degree": 2,
    "rooted_binary": true,
    "polytomy_count": 0,
    "max_depth_edges": 2,
    "mean_leaf_depth_edges": 1.6666666666666667,
    "sackin_index": 5,
    "has_branch_lengths": true,
    "all_branches_have_lengths": true,
    "negative_branch_count": 0,
    "total_branch_length": 10.0,
    "max_root_to_tip": 5.0,
    "min_root_to_tip": 4.0,
    "mean_root_to_tip": 4.333333333333333,
    "ultrametric": false,
    "treeness": 0.3,
    "colless_index": 1
  }
}
```

## 4) compare

- **作用：** 计算同一叶子集合上两棵树之间的 Robinson-Foulds 距离。

- **处理方法：** 收集每棵树的非平凡分割（把叶子分成两侧、且每侧至少两个叶子的二分；`rooted=True` 时改为至少含两个叶子的分支），RF 即两个分割集合对称差的大小（Robinson 与 Foulds 1981）。归一化值为 RF / (第一棵树分割数 + 第二棵树分割数)，而不是 RF / 最大 RF；二叉树的最大 RF（无根为 `2(n-3)`，有根为 `2(n-2)`）单独给出。分支长度和内部标签被忽略。叶子集合不同的两棵树会报错。

- **API：** `dnakit.phylogeny.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`[newick_a, newick_b]`、`{"a": ..., "b": ...}`（也可用 `left`/`right` 或 `tree1`/`tree2`），或单棵树加 `options["other"]`；树为 Newick 文本或路径。`options`（可选）：`rooted`（布尔，默认 False，即用无根分割）、`other`（第二棵树）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`robinson_foulds`、`normalized_rf`、`shared_splits`、`splits_only_first`、`splits_only_second`、`splits_first`、`splits_second`、`max_rf_binary`、`leaf_count`、`rooted` 以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.phylogeny import compare

result = compare(["((A,B),(C,D),(E,F));", "(((A,B),C),(D,(E,F)));"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 50,
    "sequence_bases": 12,
    "lengths": [6, 6],
    "sha256": "e6b7a61850243eb35ecb4e40acdf77a9dc6d014e0a4043ba5e351a3a04fcd760"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "robinson_foulds": 2,
    "normalized_rf": 0.3333333333333333,
    "shared_splits": 2,
    "splits_only_first": 1,
    "splits_only_second": 1,
    "splits_first": 3,
    "splits_second": 3,
    "max_rf_binary": 6,
    "leaf_count": 6,
    "rooted": false,
    "definition": "RF = |S1 symmetric-difference S2| over non-trivial splits (bipartitions); normalized = RF / (|S1| + |S2|)"
  }
}
```

## 5) export

- **作用：** 把 Newick 树序列化为规范化的 Newick、NEXUS、嵌套 JSON 或边列表。

- **处理方法：** 解析树后按所选格式重新输出。Newick 输出在需要时给标签重新加引号，分支长度保留 12 位有效数字；NEXUS 含 TAXA 块（`DIMENSIONS NTAX`、`TAXLABELS`）和带有该树的 TREES 块（标记为 `[&R]`）；`json` 为嵌套的 `{name, branch_length, children}` 结构；`edgelist` 用前序节点编号列出每条边的 `parent`、`child`、`child_name` 和 `length`。给出 `output` 时，文本同时写入该文件（JSON 以缩进文本写入）。纯 Python 实现。

- **API：** `dnakit.phylogeny.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：Newick 文本或 Newick 文件路径。`options`（可选）：`format`（默认 `"newick"`，或 `"nexus"`、`"json"`、`"edgelist"`）、`output`（要写入的文件路径，默认不写）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`format`、`content`（Newick 和 NEXUS 为字符串，`json` 为对象，`edgelist` 为列表）以及 `output_path`（写入的文件，未写入时为 null）。

- **示例代码：**

```python
import json

from dnakit.phylogeny import export

result = export("((A:1,B:2):3,'C d':4);", options={"format": "nexus"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 22,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "67603f444934d92f5102461f5902da759227011cabde49307f86bde1a39967c2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "nexus",
    "content": "#NEXUS\n\nBEGIN TAXA;\n  DIMENSIONS NTAX=3;\n  TAXLABELS A B 'C d';\nEND;\n\nBEGIN TREES;\n  TREE tree1 = [&R] ((A:1,B:2):3,'C d':4);\nEND;\n",
    "output_path": null
  }
}
```

## 6) infer

- **作用：** 从已比对的 DNA 序列推断基于距离的树（对 Jukes-Cantor 距离或 p 距离做邻接法或 UPGMA）。

- **处理方法：** 先按 `calculate` 的方式计算比对的成对距离矩阵（p 距离，或 Jukes-Cantor 距离 `d = -3/4 ln(1 - 4p/3)`，成对比较时只用 A/C/G/T 位点），再按 `build` 的方式用邻接法或 UPGMA 建树。这里不提供最大似然和贝叶斯推断：`options.method` 不是 `nj` / `upgma` 时，函数返回 `not_implemented`，并指向会调用外部程序（FastTree、IQ-TREE、RAxML、MrBayes）的 `infer_fast_tree` / `infer_ml_tree`。饱和的序列对（Jukes-Cantor 下 p >= 0.75）会报错，除非设置了 `saturation_distance`。

- **API：** `dnakit.phylogeny.infer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：FASTA 文本/路径、`{名称: 序列}` 或 `[{"name", "sequence"}]`；长度相等，最多 400 条序列。`options`（可选）：`method`（默认 `"nj"`，或 `"upgma"`）、`distance`（默认 `"jc69"`，或 `"p"`）、`saturation_distance`（为饱和或无法定义的序列对代替报错所用的距离，默认无）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`build` 的各字段（`newick`、`method`、`rooted`、`leaf_count`、`total_branch_length`、`negative_branch_count`、`note`）以及 `distance`、`distance_matrix` 和 `names`。`method` 不是 `nj` / `upgma` 时，返回的封装改为 `not_implemented`。

- **示例代码：**

```python
import json

from dnakit.phylogeny import infer

alignment = {
    "A": "ACGTACGTACGT",
    "B": "ACGTACGTACGA",
    "C": "TGCATGGTACGT",
    "D": "TGCATGGTACGA",
}
result = infer(alignment, options={"distance": "p"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "infer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 84,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "691762599c2228308bdd799aaa8d442e73fff378d08b81c0ea13b06176d43582"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick": "(C:0.0416666666667,D:0.0416666666667,(A:0.0416666666667,B:0.0416666666667):0.458333333333);",
    "method": "nj",
    "rooted": false,
    "leaf_count": 4,
    "total_branch_length": 0.6250000000000001,
    "negative_branch_count": 0,
    "note": "UPGMA gives a rooted ultrametric tree; neighbour-joining gives an unrooted tree drawn with a trifurcating root and may contain negative branch lengths",
    "distance": "p",
    "distance_matrix": [
      [0.0, 0.08333333333333333, 0.5, 0.5833333333333334],
      [0.08333333333333333, 0.0, 0.5833333333333334, 0.5],
      [0.5, 0.5833333333333334, 0.0, 0.08333333333333333],
      [0.5833333333333334, 0.5, 0.08333333333333333, 0.0]
    ],
    "names": ["A", "B", "C", "D"]
  }
}
```

## 7) plot

- **作用：** 把 Newick 树画成从左到右的矩形 SVG 系统发育图（phylogram）或分支图（cladogram）。

- **处理方法：** 叶子按树中的顺序放在等间距的行上，每个内部节点位于其子节点的中间，横向位置在 phylogram（`branch_lengths=True`）中为到根的距离，在 cladogram（`branch_lengths=False`）中为边数深度，并按绘图宽度缩放。树不会被重新排序。叶子名称以 SVG 文本写出。缺少分支长度的树若要画 phylogram 会报错。最多 2,000 个叶子；纯 Python 实现，不依赖绘图库。给出 `output` 时，SVG 同时写入该文件。

- **API：** `dnakit.phylogeny.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：Newick 文本或路径（最多 2,000 个叶子）。`options`（可选）：`branch_lengths`（布尔，默认 True）、`width`（整数像素 200-4000，默认 640）、`output`（SVG 的输出文件路径，默认不写）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`svg`（SVG 文本）、`leaf_count`、`branch_lengths`（实际使用的模式）和 `output_path`（写入的文件，未写入时为 null）。

- **示例代码：**

```python
import json

from dnakit.phylogeny import plot

result = plot("((A:1,B:2):3,C:4);", options={"width": 300})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a96a137a9d8a4755cd277cb9193377aab112cd1f3719f4065ede30482e9bbd5b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"300\" height=\"74\" viewBox=\"0 0 300 74\">\n<g stroke=\"#333\" stroke-width=\"1\" fill=\"none\">\n<line x1=\"10.0\" y1=\"28.0\" x2=\"10.0\" y2=\"55.0\"/>\n<line x1=\"172.0\" y1=\"19.0\" x2=\"172.0\" y2=\"37.0\"/>\n<line x1=...",
    "leaf_count": 3,
    "branch_lengths": true,
    "output_path": null
  }
}
```

（示例输出中的长数组和长文本已用 `...` 缩短。）

## 8) test

- **作用：** 为由比对构建的距离树的各个分支估计非参数自展（bootstrap）支持率。

- **处理方法：** 先用原始比对构建参考树（对 Jukes-Cantor 距离或 p 距离做邻接法或 UPGMA）。每次重复对比对的列做有放回重抽样（Felsenstein 1985），重新计算距离、重建树，并统计参考树的哪些分割再次出现（`nj` 用无根分割，`upgma` 用有根分支）。支持率是各次重复中出现的比例，取值 `[0, 1]`；返回的 Newick 把（四舍五入后的）百分比写作每个内部节点的标签。`seed` 固定时结果完全相同。重复中饱和或无法定义的序列对使用 `saturation_distance`，而不报错。计算量（序列对数 x 不同列模式数 x 重复次数）设有上限。这里只对距离树做自展，不对最大似然树做自展。

- **API：** `dnakit.phylogeny.test(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：FASTA 文本/路径、`{名称: 序列}` 或 `[{"name", "sequence"}]`；长度相等，最多 400 条序列。`options`（可选）：`replicates`（整数 1-10,000，默认 100）、`seed`（整数，默认 0）、`method`（默认 `"nj"`，或 `"upgma"`）、`distance`（默认 `"jc69"`，或 `"p"`）、`saturation_distance`（重复中饱和或无法定义的序列对所用的距离；`jc69` 默认 10.0，`p` 默认 1.0）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`newick_with_support`（带百分比标签的参考树）、`clades`（参考树的每个分割一项，含叶子集合 `clade` 和 `support`；无根树的 `clade` 取不含字母序最小叶子的那一侧）、`replicates`、`seed`、`method`、`distance` 和 `mean_support`。

- **示例代码：**

```python
import json

from dnakit.phylogeny import test

alignment = {
    "A": "AAAAAAGGGG",
    "B": "AAAAAATTTT",
    "C": "CCCCCCGGGG",
    "D": "CCCCCCTTTT",
}
result = test(alignment, options={"replicates": 100, "seed": 1, "distance": "p"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "phylogeny",
  "operation": "test",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "79dd2c6923081105920d7bb30423f9f3c0e984f0918d4115020eb9b7ca4fbf39"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick_with_support": "(C:0.2,D:0.2,(A:0.2,B:0.2)82:0.4);",
    "clades": [
      {
        "clade": ["C", "D"],
        "support": 0.82
      }
    ],
    "replicates": 100,
    "seed": 1,
    "method": "nj",
    "distance": "p",
    "mean_support": 0.82
  }
}
```
