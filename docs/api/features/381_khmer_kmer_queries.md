# khmer k-mer 查询：哈希、丰度、de Bruijn 图邻域与路径

本页对应工作清单 P27 包中查询 k-mer 表或读段 de Bruijn 图的 khmer 条目：k-mer 哈希值、逐 k-mer 丰度统计、按丰度修剪、谱错误位置、存在图并集、k-mer 的邻居与度、连通分量大小、分支点以及无分支路径。每个函数都从一种输入构建所需的 khmer 表：已保存的表（`graph_path`）、本地 FASTA/FASTQ 读段（`reads_path`）或 DNA 字符串（`sequences`）。仅参数不同的条目合并：两个哈希条目是带 `as_set` 的一个函数；最小值、最大值、中位数/均值/标准差、逐 k-mer 计数和中位数阈值条目是一个函数（`khmer_sequence_abundance`，同时对应 `count-median.py`）；邻居列表与度是一个函数；分量大小与半径计数是一个函数；遍历与组装线性路径是带 `output` 的一个函数。载入表已由 `query_saved_khmer_table` 和 `inspect_khmer_graph` 覆盖，保存已由 `khmer_build_countgraph` 和 `khmer_build_nodegraph` 覆盖，MurmurHash3 哈希已由 `hash_kmer_murmur3` 覆盖。

khmer 是可选依赖（`pip install "dnakit[bio-khmer]"`），只在函数运行时导入；已用 khmer 2.1.1 测试。默认值遵循 khmer：k = 32、4 张表、每张表 1e6 个桶；示例读段很短，所以使用 k = 15。k-mer 与其反向互补序列是图中的同一个节点。

## 示例文件

`reads.fa` 含三条读段：`a` 与 `c` 的前 30 个碱基相同，所以 de Bruijn 图在此分支；`b` 与它们无关。`other.fa` 含两条前 32 个碱基相同的读段。

- **示例文件：** `reads.fa`

```text
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

- **示例文件：** `other.fa`

```text
>x1
AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC
>x2
AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT
```

## 1) `dnakit.microbiome.khmer_kmer_hashes()`

- **作用：** 按顺序列出 DNA 序列每个 k-mer 的哈希值，或只列出不重复的值。
- **处理方法：** 创建指定 k 的空 khmer 计数图，调用 `Countgraph.get_kmer_hashes`（有序列表）或 `get_kmer_hashes_as_hashset`（`as_set=True`，排序后返回）。数值是 khmer 对规范 k-mer 的 2 位正向哈希，所以 k-mer 与其反向互补序列得到同一个值；`kmers` 给出这些值对应的 k-mer。输入转为大写。对应工作清单 G1219 与 G1220。
- **API：** `dnakit.microbiome.khmer_kmer_hashes(sequence, *, ksize=32, as_set=False)`。
- **输入：** `sequence`（由 A、C、G、T 组成的 DNA）、`ksize`（1 到 32，默认 32）、`as_set`。
- **输出：** `n_hashes`、`hashes`、`kmers`、`provenance`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_kmer_hashes

seq = "TCGCTGCTGTCGGACTCCTA"
r = khmer_kmer_hashes(seq, ksize=15, verbose=False)
print(r["n_hashes"], r["hashes"][:2])
print(r["kmers"][:2])
s = khmer_kmer_hashes(seq + seq, ksize=15, as_set=True, verbose=False)
print(s["n_hashes"])
```

- **示例输出：**

```text
6 [463369970, 229412411]
['TCGCTGCTGTCGGAC', 'AGTCCGACAGCAGCG']
20
```

## 2) `dnakit.microbiome.khmer_graph_summary()`

- **作用：** 载入已保存的 khmer 表，并估计其中不重复 k-mer 的数量。
- **处理方法：** 检查文件头，用 `khmer.load_countgraph`、`load_nodegraph`、`load_counttable` 或 `load_nodetable`（`graph_type`）载入表，并报告 k、各表大小和已占用的桶数。khmer 自带的 `n_unique_kmers()` 计数器不随文件保存，载入后读数为 0，所以不重复 k-mer 数用第一张哈希表上的线性计数估计：`-size * ln(1 - occupied / size)`。对应工作清单 G1215。
- **API：** `dnakit.microbiome.khmer_graph_summary(graph_path, *, graph_type='countgraph')`。
- **输入：** `graph_path`（已保存的表）、`graph_type`（`countgraph`、`nodegraph`、`counttable` 或 `nodetable`）。
- **输出：** `ksize`、`hashsizes`、`n_occupied`、`estimated_unique_kmers`（第一张表已满时为 None）、`provenance`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_gap_functions import khmer_build_countgraph
from dnakit.microbiome.khmer_p27_graph_functions import khmer_graph_summary

table = os.path.join(tempfile.mkdtemp(), "reads.countgraph")
khmer_build_countgraph("reads.fa", table, ksize=15, verbose=False)
r = khmer_graph_summary(table, verbose=False)
print(r["ksize"], r["n_occupied"], r["estimated_unique_kmers"])
```

- **示例输出：**

```text
15 122 122
```

## 3) `dnakit.microbiome.khmer_sequence_abundance()`

- **作用：** 在计数图中统计查询序列的丰度：最小值、最大值、中位数、均值、标准差、逐 k-mer 计数以及中位数阈值检验。
- **处理方法：** 构建或载入计数图（`graph_path`、`reads_path` 或 `sequences`），对每条查询调用 `get_min_count`、`get_max_count`、`get_median_count`（中位数、均值、标准差），可选调用 `get_kmer_counts`，给出 `median_threshold` 时调用 `median_at_least`。查询来自 `sequences_to_query` 和/或 FASTA/FASTQ `query_path`；与 khmer 的 `count-median.py` 一致，文件查询中的 N 按 A 处理，短于 k 的查询被跳过。`output_csv` 写出 `count-median.py` 的表（name、median、average、stddev、seqlen）。对应工作清单 G1221 至 G1225。
- **API：** `dnakit.microbiome.khmer_sequence_abundance(*, sequences_to_query=None, query_path=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, median_threshold=None, include_counts=False, output_csv=None, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequences_to_query` 和/或 `query_path`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`；`median_threshold`、`include_counts`、`output_csv`、`max_reads`、`max_bases`。
- **输出：** `n_queries`、`skipped_shorter_than_k`、`results`（`name`、`length`、`min`、`max`、`median`、`mean`、`sd`，可选 `median_at_least` 与 `counts`）、`output_csv`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_sequence_abundance

r = khmer_sequence_abundance(
    sequences_to_query=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA"],
    reads_path="reads.fa",
    ksize=15,
    median_threshold=2,
    include_counts=True,
    verbose=False,
)
row = r["results"][0]
print(row["min"], row["max"], row["median"], round(row["mean"], 3), row["median_at_least"])
print(row["counts"][:20])
```

- **示例输出：**

```text
1 2 1 1.348 False
[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1]
```

## 4) `dnakit.microbiome.khmer_trim_by_abundance()`

- **作用：** 在第一个 k-mer 丰度越过阈值的位置修剪序列。
- **处理方法：** 构建或载入计数图，调用 `Countgraph.trim_below_abundance`（`mode="trim_below_abundance"`，在第一个计数高于阈值的 k-mer 处截断）或 `trim_on_abundance`（`mode="trim_on_abundance"`，在第一个计数低于阈值的 k-mer 处截断）。返回修剪后的序列和 khmer 的修剪位置。对应工作清单 G1227。
- **API：** `dnakit.microbiome.khmer_trim_by_abundance(sequence, cutoff, *, mode='trim_below_abundance', graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`、`cutoff`、`mode`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`。
- **输出：** `trimmed_sequence`、`trim_position`、`input_length`、`mode`、`cutoff`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_trim_by_abundance

low = khmer_trim_by_abundance(
    "TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", 2, mode="trim_on_abundance", reads_path="reads.fa", ksize=15, verbose=False
)
print(low["trim_position"], low["trimmed_sequence"])
high = khmer_trim_by_abundance("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", 1, reads_path="reads.fa", ksize=15, verbose=False)
print(high["trim_position"], repr(high["trimmed_sequence"]))
```

- **示例输出：**

```text
30 TCGCTGCTGTCGGACTCCTAGTTACGTGGC
0 ''
```

## 5) `dnakit.microbiome.khmer_spectral_error_positions()`

- **作用：** 定位序列中落在低丰度 k-mer 内的位置，这是测序错误的特征。
- **处理方法：** 构建或载入计数图，调用 `Countgraph.find_spectral_error_positions(sequence, max_abundance)`。位置是查询序列中从 0 开始的偏移。对应工作清单 G1228。
- **API：** `dnakit.microbiome.khmer_spectral_error_positions(sequence, *, max_abundance=1, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`、`max_abundance`（默认 1）；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`。
- **输出：** `positions`、`max_abundance`、`graph`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_spectral_error_positions

truth = "TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"
error = truth[:30] + ("A" if truth[30] != "A" else "C") + truth[31:]
r = khmer_spectral_error_positions(error, sequences=[truth] * 4 + [error], ksize=15, verbose=False)
print(r["positions"])
```

- **示例输出：**

```text
[30]
```

## 6) `dnakit.microbiome.khmer_merge_nodegraphs()`

- **作用：** 两个已保存存在图的并集，写入新文件。
- **处理方法：** 载入两个存在图，调用 `Nodegraph.update`（把第二个并入第一个），再用 `Nodegraph.save` 保存结果。两个图的 k、表数量和表大小必须相同。输出文件必须不存在。对应工作清单 G1231。
- **API：** `dnakit.microbiome.khmer_merge_nodegraphs(first_graph_path, second_graph_path, output_path)`。
- **输入：** `first_graph_path`、`second_graph_path`、`output_path`。
- **输出：** `unique_kmers_first`、`unique_kmers_second`、`unique_kmers_merged`、`output_path`、`ksize`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_gap_functions import khmer_build_nodegraph
from dnakit.microbiome.khmer_p27_graph_functions import khmer_merge_nodegraphs

folder = tempfile.mkdtemp()
one = os.path.join(folder, "one")
two = os.path.join(folder, "two")
khmer_build_nodegraph("reads.fa", one, ksize=15, verbose=False)
khmer_build_nodegraph("other.fa", two, ksize=15, verbose=False)
r = khmer_merge_nodegraphs(one, two, os.path.join(folder, "union"), verbose=False)
print(r["unique_kmers_first"], r["unique_kmers_second"], r["unique_kmers_merged"])
```

- **示例输出：**

```text
0 0 0
```

## 7) `dnakit.microbiome.khmer_kmer_neighbors()`

- **作用：** 列出 de Bruijn 图中与某个 k-mer 相邻的 k-mer 及其度。
- **处理方法：** 构建或载入图（`graph_type` 为存在图或计数图），调用 `neighbors` 与 `kmer_degree`。邻居是图中存在、并与查询在任一侧重叠 k-1 个碱基的 k-mer。对应工作清单 G1232 与 G1234。
- **API：** `dnakit.microbiome.khmer_kmer_neighbors(kmer, *, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`。
- **输入：** `kmer`（长度为 k）；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`、`graph_type`。
- **输出：** `degree`、`neighbors`（k-mer 列表）、`kmer`、`graph`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_kmer_neighbors

fork = "TCCTAGTTACGTGGC"
r = khmer_kmer_neighbors(fork, reads_path="reads.fa", ksize=15, verbose=False)
print(r["degree"], r["neighbors"])
r = khmer_kmer_neighbors("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["degree"], r["neighbors"])
```

- **示例输出：**

```text
3 ['CTCCTAGTTACGTGG', 'TGCCACGTAACTAGG', 'CCTAGTTACGTGGCG']
1 ['CGTGCGAGTGTTGCG']
```

## 8) `dnakit.microbiome.khmer_connected_component_size()`

- **作用：** de Bruijn 图中与某个 k-mer 相连的 k-mer 数量，可选同时给出图半径内的数量。
- **处理方法：** 构建或载入图，调用 `calc_connected_graph_size(kmer, max_size)`；给出 `count_within_radius` 时还调用 `count_kmers_within_radius(kmer, radius, max_size)`。与 khmer 一致，`max_size` 为 0 表示不设上限。对应工作清单 G1233 与 G1235。
- **API：** `dnakit.microbiome.khmer_connected_component_size(kmer, *, max_size=0, count_within_radius=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`。
- **输入：** `kmer`（长度为 k）、`max_size`、`count_within_radius`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`、`graph_type`。
- **输出：** `component_size`、`kmers_within_radius`、`max_size`、`graph`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_connected_component_size

r = khmer_connected_component_size(
    "TCGCTGCTGTCGGAC", reads_path="reads.fa", ksize=15, count_within_radius=5, verbose=False
)
print(r["component_size"], r["kmers_within_radius"])
r = khmer_connected_component_size("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["component_size"])
```

- **示例输出：**

```text
76 6
46
```

## 9) `dnakit.microbiome.khmer_high_degree_kmers()`

- **作用：** 查找分支点：序列中在图里拥有两个以上邻居的 k-mer。
- **处理方法：** 构建或载入图，调用 `find_high_degree_nodes(sequence)`，并把哈希值转换回 k-mer。对应工作清单 G1236。
- **API：** `dnakit.microbiome.khmer_high_degree_kmers(sequence, *, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`、`graph_type`。
- **输出：** `n_high_degree`、`high_degree_kmers`（按哈希值排序）、`graph`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_high_degree_kmers

r = khmer_high_degree_kmers("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_high_degree"], r["high_degree_kmers"])
r = khmer_high_degree_kmers("CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA", reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_high_degree"])
```

- **示例输出：**

```text
1 ['TCCTAGTTACGTGGC']
0
```

## 10) `dnakit.microbiome.khmer_linear_path()`

- **作用：** 沿着从某个 k-mer 出发的无分支路径前进，返回组装后的单条序列或经过的 k-mer 列表。
- **处理方法：** 构建或载入图。`output="sequence"` 调用 `assemble_linear_path`，把路径作为 DNA 字符串返回。`output="kmers"` 调用 `traverse_linear_path`，返回经过的 k-mer、数量以及遇到的分支 k-mer；`stop_sequence` 提供终止遍历的分支点（该序列的高度数 k-mer），不给时会遍历整个连通分量。对应工作清单 G1237 与 G1238。
- **API：** `dnakit.microbiome.khmer_linear_path(kmer, *, output='sequence', stop_sequence=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`。
- **输入：** `kmer`（起点 k-mer）、`output`、`stop_sequence`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`、`graph_type`。
- **输出：** `path_sequence` 与 `path_length`，或 `path_size`、`path_kmers` 与 `branch_kmers`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_linear_path

r = khmer_linear_path("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["path_length"], r["path_sequence"] == "CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA")
k = khmer_linear_path(
    "TCGCTGCTGTCGGAC", output="kmers", stop_sequence="TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", reads_path="reads.fa", ksize=15, verbose=False
)
print(k["path_size"], k["branch_kmers"])
```

- **示例输出：**

```text
60 True
15 ['TCCTAGTTACGTGGC']
```
