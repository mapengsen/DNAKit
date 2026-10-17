# 三维基因组

本类收录 **9** 项有文档的候选功能，全部在本地原生计算出真实结果，处理对象是以稠密或稀疏文本、路径或内存矩阵给出的 Hi-C 接触矩阵（不读取二进制 `.hic`/`.cool` 文件）。`annotate`、`calculate`、`compare` 和 `normalize` 需要 numpy，未安装时返回 `status: "not_implemented"`；没有任何函数是无条件的 `not_implemented`。有两处结果是明确的规则而非统计调用程序：`annotate` 给出的 TAD 边界是绝缘分数极小值（基于窗口的规则，没有 p 值），区室标签是观测/期望相关矩阵首个特征向量的符号，除非提供参考轨道，否则其方向是任意的。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) aggregate

- **作用：** 通过对 `factor` x `factor` 的 bin 块求和来降低 Hi-C 接触矩阵的分辨率。

- **处理方法：** 结果的第 k 个 bin 在两个维度上都对原始 bin k*`factor` .. k*`factor`+`factor`-1 求和（n 不是 `factor` 的倍数时最后一个 bin 覆盖的原始 bin 较少），因此接触总数守恒。新的 bin 大小为 `bin_size` x `factor`。只是求和，不做任何归一化。

- **API：** `dnakit.genome3d.aggregate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`entity`（只支持 `"contact_matrix"`）；`factor`（必填，整数 >= 2）；`bin_size`（可选，bp）；`output`（可选，稠密 TSV 路径；已有文件需要 `overwrite=true`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `factor`、`n_bins`、`original_n_bins`、`bin_size`（仅在给出 `bin_size` 时）、`matrix`，写出文件时还有 `output_path`。超过 1000 个 bin 的矩阵不会放入 `result`（改为 `matrix_omitted: true`），此时需要 `output`。

- **示例代码：**

```python
import json
from dnakit.genome3d import aggregate

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = aggregate(matrix, options={"factor": 2, "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "aggregate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "factor": 2,
    "n_bins": 2,
    "original_n_bins": 4,
    "bin_size": 20000,
    "matrix": [
      [
        24.0,
        9.0
      ],
      [
        9.0,
        24.0
      ]
    ]
  }
}
```

## 2) annotate

- **作用：** 为接触矩阵的 bin 标注 TAD 边界（绝缘分数极小值）或 A/B 区室标签（`options["entity"]`）。

- **处理方法：** `tad`（默认）使用绝缘分数极小值规则，不是带统计检验的 TAD 调用程序。绝缘分数的算法同 `calculate`（窗口 `window` w，默认 5）：S(i) = M[a][b] 的均值，a 取 [i-w, i-1]、b 取 [i, i+w-1]，score = log2(S(i) / 有定义 bin 上 S 的均值)。bin i 是边界，当且仅当其得分低于 `threshold`（默认 0.0，即绝缘度低于整个矩阵的平均水平），并且在其两侧 `boundary_radius` 个 bin（默认 = `window`）内是最小值，并列时取最左侧的 bin。`tads` 是 bin 区间 [`start_bin`, `end_bin`)（右端不含）。设矩阵有 n 个 bin、边界为 b_1 < ... < b_m，则内部结构域为 [b_k, b_k+1)（`edge: false`）；`include_edge_domains`（默认 true）时还会报告两个边缘结构域，标记为 `edge: true`：[0, b_1)（b_1 = 0 时省略）和 [b_m, n)（最后一个边界开启最后一个结构域）。因此只有一个边界时最多得到两个边缘结构域，完全没有边界时整个矩阵 [0, n) 是一个边缘结构域。边缘结构域被矩阵边界截断，其大小与内部结构域不可比；`include_edge_domains: false` 则只返回内部结构域（边界少于两个时为空列表）。不提供 p 值或边界强度。`compartment`：对观测/期望接触（期望 = 非空 bin 上每条对角线的平均接触）的 Pearson 相关矩阵取首个特征向量，按符号标注 A/B；特征向量 >= 0 为 A，< 0 为 B。特征向量的符号本身是任意的：未给 `reference` 时取绝对值最大的元素为正，所以哪个区室叫 A 本身没有意义；给定逐 bin 的 `reference` 轨道（例如 GC 含量）时，会翻转符号使特征向量与之正相关。至少需要 3 个非空 bin。需要 numpy；缺少时函数返回 `status: "not_implemented"`，`required_dependency` 为 `pip install numpy`。

- **API：** `dnakit.genome3d.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options["entity"]` 取 `"tad"`（默认）或 `"compartment"`。`tad`：`window`（默认 5）、`boundary_radius`（默认等于 window）、`threshold`（默认 0.0）、`include_edge_domains`（默认 `true`）、`bin_size`（可选，bp）。`compartment`：`reference`（可选，逐 bin 的数值列表）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`tad`：`result` 含 `entity`、`window`、`boundary_radius`、`threshold`、`boundary_bins`（0 起始的 bin 索引）、`tads`（`start_bin`、`end_bin`（不含）、`edge`：结构域贴着矩阵起点或终点、只有一侧有边界或没有边界时为 `true`），给定 `bin_size` 时还有 `boundary_bp` 和 `bin_size`。`compartment`：`entity`、`eigenvector`（空 bin 为 `null`）、`labels`（`A`/`B`/`null`）、`eigenvalue`、`explained_variance`（首个特征值 / 相关矩阵特征值之和）和 `oriented_by_reference`。

- **示例代码：**

```python
import json
from dnakit.genome3d import annotate

n = 15
matrix = [[(10.0 / (1 + abs(i - j))) if i // 5 == j // 5 else 1.0 for j in range(n)] for i in range(n)]   # three 5-bin domains
result = annotate(matrix, options={"entity": "tad", "window": 2, "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 15,
    "text_length": 1440,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "78895d0d58ea88a971d7a5bed4fda84be9cc4baee71bf82b17abc664229d7c1b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "tad",
    "window": 2,
    "boundary_radius": 2,
    "threshold": 0.0,
    "boundary_bins": [
      5,
      10
    ],
    "tads": [
      {
        "start_bin": 0,
        "end_bin": 5,
        "edge": true
      },
      {
        "start_bin": 5,
        "end_bin": 10,
        "edge": false
      },
      {
        "start_bin": 10,
        "end_bin": 15,
        "edge": true
      }
    ],
    "boundary_bp": [
      50000,
      100000
    ],
    "bin_size": 10000
  }
}
```

## 3) calculate

- **作用：** 根据接触矩阵计算距离衰减曲线、绝缘分数或 A/B 区室特征向量（`options["entity"]`）。

- **处理方法：** `distance_decay`（默认）：P(s) = 对角线 s = |i - j| 上所有有效 bin 对的平均接触（行和为零的空 bin 被排除）；`slope` 是 `s_min` <= s <= `s_max` 时 log10 P(s) 对 log10 s 的最小二乘斜率。`insulation`：S(i) = M[a][b] 的均值，a 取 [i-w, i-1]、b 取 [i, i+w-1]（对角线上、角点位于 bin i 起点的方块；w <= i <= n-w 时有定义，否则为 `null`），score = log2(S(i) / 有定义 bin 上 S 的均值)；低分表示绝缘位置。`compartment`：对观测/期望矩阵（期望 = 每条对角线的平均值）的 Pearson 相关矩阵取首个特征向量（numpy `eigh`）；符号给出 A（>= 0）或 B（< 0），除非提供逐 bin 的 `reference`（如 GC），否则方向是任意的，此时取绝对值最大的元素为正。需要 numpy；缺少时函数返回 `status: "not_implemented"`，`required_dependency` 为 `pip install numpy`。

- **API：** `dnakit.genome3d.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options["entity"]` 取 `"distance_decay"`（默认）、`"insulation"` 或 `"compartment"`。`distance_decay`：`bin_size`（可选，bp，会加入 `distance_bp`）、`s_min`（默认 2）、`s_max`（默认 n/2）。`insulation`：`window`（默认 5）。`compartment`：`reference`（可选，逐 bin 的数值列表）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`distance_decay`：`result` 含 `entity`、`curve`（每条对角线有 `distance_bins`、`n_pairs`、`mean_contact`，给定 `bin_size` 时还有 `distance_bp`）和 `scaling`（`slope`、`points`）。`insulation`：`entity`、`window`、每个 bin 的 `mean_contact`（S(i)）和 `score`。`compartment`：`entity`、`eigenvector`、`labels`、`eigenvalue`、`explained_variance`、`oriented_by_reference`。

- **示例代码：**

```python
import json
from dnakit.genome3d import calculate

n = 15
matrix = [[(10.0 / (1 + abs(i - j))) if i // 5 == j // 5 else 1.0 for j in range(n)] for i in range(n)]   # three 5-bin domains
result = calculate(matrix, options={"entity": "distance_decay", "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "genome3d",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 15,
    "text_length": 1440,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "78895d0d58ea88a971d7a5bed4fda84be9cc4baee71bf82b17abc664229d7c1b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "distance_decay",
    "curve": [
      {
        "distance_bins": 0,
        "n_pairs": 15,
        "mean_contact": 10.0,
        "distance_bp": 0
      },
      {
        "distance_bins": 1,
        "n_pairs": 14,
        "mean_contact": 4.428571428571429,
        "distance_bp": 10000
      },
      {
        "distance_bins": 2,
        "n_pairs": 13,
        "mean_contact": 2.6153846153846154,
        "distance_bp": 20000
      },
      {
        "distance_bins": 3,
        "n_pairs": 12,
        "mean_contact": 1.75,
        "distance_bp": 30000
      },
      "..."
    ],
    "scaling": {
      "slope": -0.8199134580718668,
      "points": 6
    }
  }
}
```

## 4) compare

- **作用：** 通过上三角接触值的相关性比较两个接触矩阵。

- **处理方法：** 从两个矩阵中取出 `min_distance` <= j - i <= `max_distance` 的 bin 对的接触值（任一矩阵中含空 bin 的对被跳过），计算 Pearson 相关、Spearman 相关（对平均秩求 Pearson）以及平均绝对差。`distance_normalize` 为真时，先把每个元素除以其所在对角线的平均值（观测/期望变换）。这只是接触值的普通相关，不是按距离分层校正的（HiCRep），也不给出显著性。需要 numpy；缺少时函数返回 `status: "not_implemented"`，`required_dependency` 为 `pip install numpy`。

- **API：** `dnakit.genome3d.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{"a": 矩阵, "b": 矩阵}` 或两个同样大小的矩阵的列表，每个都是方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin。`options`：`entity`（只支持 `"contact_matrix"`）；`min_distance`（默认 1）；`max_distance`（可选）；`distance_normalize`（默认 `false`）。至少需要两个可比较的 bin 对。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `pairs`（被比较的 bin 对数）、`pearson`、`spearman`（任一侧方差为 0 时为 `null`）和 `mean_absolute_difference`。

- **示例代码：**

```python
import json
from dnakit.genome3d import compare

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
other = [[8, 3, 2, 1],
         [3, 9, 5, 2],
         [2, 5, 8, 3],
         [1, 2, 3, 7]]
result = compare({"a": matrix, "b": other}, options={"min_distance": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 126,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e1188dc81e4040ecb0c3c3932535a69005cf8552399b86faff6d9b4f1f159e6c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": 6,
    "pearson": 0.8443565629465051,
    "spearman": 0.9534625892455922,
    "mean_absolute_difference": 0.5
  }
}
```

## 5) export

- **作用：** 把接触矩阵序列化为稠密 TSV 或稀疏三元组文本并返回文本，可选写入文件。

- **处理方法：** `dense` 每个 bin 输出一行、以制表符分隔；`sparse` 对上三角（含对角线）的非零元素输出 `i<TAB>j<TAB>value`（0 起始 bin），可被 `read(format="sparse")` 读回。整数不带小数点。`write` 使用相同的序列化，但要求提供输出路径。

- **API：** `dnakit.genome3d.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`out_format`（`"dense"` 默认 | `"sparse"`）；`output`（可选文件路径）；`overwrite`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `format`、`text`、`output_path`（已写入的文件，未写入为 `null`）、`n_bins`、`total_contacts`（上三角含对角线之和）和 `nonzero_entries`（完整矩阵中的非零单元格数）。

- **示例代码：**

```python
import json
from dnakit.genome3d import export

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = export(matrix, options={"out_format": "sparse"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "sparse",
    "text": "0\t0\t8\n0\t1\t4\n0\t2\t2\n0\t3\t1\n1\t1\t8\n1\t2\t4\n1\t3\t2\n2\t2\t8\n2\t3\t4\n3\t3\t8\n",
    "output_path": null,
    "n_bins": 4,
    "total_contacts": 49.0,
    "nonzero_entries": 16
  }
}
```

## 6) normalize

- **作用：** 用迭代校正（ICE）或普通覆盖度归一化对 Hi-C 接触矩阵做归一化（平衡）。

- **处理方法：** `ice`（默认；Imakaev 等 2012）：非零元素少于 `min_nonzero` 的 bin 被屏蔽（并先把 `ignore_diagonals` 条近对角线清零）。每次迭代把有效 bin 的行和 s_i 除以其均值，并把 M_ij 除以 s_i s_j；所有 s_i 的乘积即偏差 b_i，结果为 原始值 / (b_i b_j)。当 max |s_i - 1| < `tol` 或达到 `max_iter` 次迭代时停止，最后缩放矩阵使每个有效行之和等于 `target_row_sum`（与 cooler 一致）。`coverage`：M_ij / (w_i w_j)，其中 w_i 是行和除以平均行和；不做迭代。`converged` 表示是否达到 ICE 容差，`max_row_sum_deviation` 是行和与 `target_row_sum` 的最大剩余偏差。需要 numpy；缺少时函数返回 `status: "not_implemented"`，`required_dependency` 为 `pip install numpy`。

- **API：** `dnakit.genome3d.normalize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`entity`（只支持 `"contact_matrix"`）；`method`（`"ice"` 默认 | `"coverage"`）；ICE：`max_iter`（默认 200）、`tol`（默认 1e-5）、`ignore_diagonals`（默认 0；1 表示只清主对角线）、`min_nonzero`（默认 1）、`target_row_sum`（默认 1.0）；`output`（可选，稠密 TSV 路径；已有文件需要 `overwrite=true`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `entity`、`method`、`matrix`（平衡后，对称）、`bias`（被屏蔽的 bin 为 `null`）、`masked_bins`，`ice` 还有 `iterations`、`converged` 和 `max_row_sum_deviation`；写出文件时另有 `output_path`。超过 1000 个 bin 的矩阵不会放入 `result`（改为 `matrix_omitted: true`），此时需要 `output`。

- **示例代码：**

```python
import json
from dnakit.genome3d import normalize

matrix = [[10, 4, 2],
          [4, 20, 6],
          [2, 6, 30]]
result = normalize(matrix, options={"method": "ice", "max_iter": 200, "tol": 1e-5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e5e798d23a5d2fead0ab177eab4c59bd3ac887cf29f21712f11738c75c293e99"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "contact_matrix",
    "method": "ice",
    "bias": [
      3.719048104659864,
      5.595833606989041,
      6.34110030203879
    ],
    "masked_bins": [],
    "iterations": 26,
    "converged": true,
    "max_row_sum_deviation": 9.753559970127057e-06,
    "matrix": [
      [
        0.7229968552636007,
        0.1922044344481455,
        0.08480736642436058
      ],
      [
        0.1922044344481455,
        0.6387051401499937,
        0.16909152282572365
      ],
      [
        0.08480736642436058,
        0.16909152282572365,
        0.7460913571899457
      ]
    ]
  }
}
```

## 7) plot

- **作用：** 把接触矩阵画成 SVG 热图（按 log10(1 + 接触数) 从白到红）。

- **处理方法：** 手写的 SVG（不依赖绘图库）。边长超过 `max_cells` 的矩阵先通过块求和降低分辨率（取使边长不超过 `max_cells` 的最小整数因子）。每个单元格是带 `data-i`、`data-j`、`data-value` 的 `<rect class="cell">`，颜色按 log10(1 + 值) 相对显示网格的最大值着色；根元素带有 `data-bins` 和 `data-factor`。

- **API：** `dnakit.genome3d.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`entity`（只支持 `"contact_matrix"`）；`max_cells`（默认 150，范围 2-400）；`output`（可选 `.svg` 路径）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `svg`（SVG 文本）和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.genome3d import plot

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = plot(matrix, options={"max_cells": 150})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "genome3d",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"40\" height=\"40\" data-bins=\"4\" data-factor=\"1\"><rect class=\"cell\" data-i=\"0\" data-j=...",
    "output_path": null
  }
}
```

## 8) read

- **作用：** 从稠密或稀疏文本、文本文件路径或内存中的矩阵读取 Hi-C 接触矩阵。

- **处理方法：** 以 `#` 开头的行和空行被跳过。稠密行按空白、逗号或制表符切分；非数值的首行和非数值的标签列被丢弃。稀疏文本中每个无序对列一次 `bin_i bin_j value`（重复的对会报错），并镜像为对称矩阵。数值必须有限且非负，矩阵必须是方阵且对称（容差 1e-9），否则抛出输入错误。不读取二进制 `.hic`、`.cool`、`.mcool`、`.scool` 文件：调用会抛出错误，并提示使用 `dnakit.epigenome` 中的 cooler 封装。

- **API：** `dnakit.genome3d.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`format`、`n_bins`（仅稀疏格式）、`output`（可选，稠密 TSV 路径；已有文件需要 `overwrite=true`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `n_bins`、`total_contacts`（上三角含对角线之和）、`nonzero_entries`（完整矩阵中的非零单元格数）和 `matrix`；写出文件时另有 `output_path`。超过 1000 个 bin 的矩阵不会放入 `result`（改为 `matrix_omitted: true`），此时需要 `output`。

- **示例代码：**

```python
import json
from dnakit.genome3d import read

sparse = "0\t0\t5\n0\t1\t2\n1\t1\t4\n1\t2\t1\n2\t2\t6\n"      # bin_i, bin_j, value (each unordered pair once)
result = read(sparse, options={"format": "sparse"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d1c502b64c132564d19debb6c077eb023535f38fa72064ccb2321743579d4340"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_bins": 3,
    "total_contacts": 18.0,
    "nonzero_entries": 7,
    "matrix": [
      [
        5.0,
        2.0,
        0.0
      ],
      [
        2.0,
        4.0,
        1.0
      ],
      [
        0.0,
        1.0,
        6.0
      ]
    ]
  }
}
```

## 9) write

- **作用：** 把接触矩阵写入稠密 TSV 或稀疏三元组文件。

- **处理方法：** 序列化方式与 `export` 相同（稠密：每个 bin 一行、制表符分隔；稀疏：上三角含对角线的非零元素 `i<TAB>j<TAB>value`，0 起始），但必须提供输出路径，且除非 `overwrite` 为真，否则已有文件会报错。

- **API：** `dnakit.genome3d.write(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：方形、对称、非负的接触矩阵：行列表、`{"matrix": ...}`、稠密或稀疏文本，或路径（`.hic`、`.cool`、`.mcool`、`.scool` 文件会被拒绝，并提示改用 cooler 封装）；数值须为有限的非负数，矩阵在 1e-9 内对称，至多 3000 个 bin；文本输入时 `options["format"]` 取 `"dense"`（默认：以空白/逗号/制表符分隔的行；跳过非数值的表头行和/或标签列）或 `"sparse"`（行 `bin_i bin_j value`，0 起始，每个无序对只列一次；可选 `n_bins`）。`options`：`output`（必填，文件路径）；`out_format`（`"dense"` 默认 | `"sparse"`）；`overwrite`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"genome3d"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256）、`warnings` 和 `result`。`result` 含 `format`、`output_path`、`bytes_written`、`n_bins`、`total_contacts` 和 `nonzero_entries`。文件作为副作用被写出。

- **示例代码：**

```python
import json
from dnakit.genome3d import write

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = write(matrix, options={"output": "contacts.tsv", "out_format": "dense"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "genome3d",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "dense",
    "output_path": "contacts.tsv",
    "bytes_written": 32,
    "n_bins": 4,
    "total_contacts": 49.0,
    "nonzero_entries": 16
  }
}
```
