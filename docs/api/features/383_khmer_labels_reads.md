# khmer 带标签图、HyperLogLog 合并、读段比对与读段文件工具

本页对应工作清单 P27 包中其余的 khmer 条目：带标签图（`GraphLabels`）、合并两个 HyperLogLog k-mer 计数器、用 khmer 的隐马尔可夫模型把读段比对到 k-mer 图，以及读段文件脚本 `extract-paired-reads.py`、`split-paired-reads.py` 与 `readstats.py`。带标签图是标签 k-mer 带有整数标签的存在图。文件中的读段以其在文件中的记录序号作标签（`consume_seqfile_and_tag_with_labels`）；字符串则使用显式标签。

仅参数不同的条目合并：从文件或字符串读入、载入与保存标签和标签 k-mer 都是 `khmer_label_reads` 的选项；扫描标签与扫描标签 k-mer 是带 `output` 的一个函数；双向延伸与仅向前延伸的比对是带 `forward_only` 的一个函数。解析 FASTA/FASTQ 读段和交错的配对读段（`ReadParser.iter_reads`、`iter_read_pairs`）已由 `read_paired_fastx` 覆盖。khmer 是可选依赖（`pip install "dnakit[bio-khmer]"`）；已用 khmer 2.1.1 测试。默认值遵循 khmer：k = 32、4 张表、每张表 1e6 个桶；示例使用 k = 15。khmer 的 `GraphLabels` 对象只是一层薄封装，它自己的图方法会崩溃，所以 DNAKit 只在它上面使用标签方法，图方法都用在它的 `graph` 属性上。

## 示例文件

`lab.fa` 含两条前 32 个碱基相同的读段。`mixed.fa` 含两对配对读段（`/1` 与 `/2`）和一条未配对读段。`paired.fa` 含两对交错排列的配对读段。

- **示例文件：** `lab.fa`

```text
>x1
AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC
>x2
AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT
```

- **示例文件：** `mixed.fa`

```text
>p1/1
ACGTACGTACGGTT
>p1/2
GGGGCCCCAATTGC
>s1
TTTTTTTTTTGGAA
>p2/1
AAAACCCCGGTTAC
>p2/2
CCCCAAAATTGGCA
```

- **示例文件：** `paired.fa`

```text
>p1/1
ACGTACGTACGGTT
>p1/2
GGGGCCCCAATTGC
>p2/1
AAAACCCCGGTTAC
>p2/2
CCCCAAAATTGGCA
```

## 1) `dnakit.microbiome.khmer_label_reads()`

- **作用：** 把读段插入带标签图并给其 k-mer 打上标签；可载入与保存标签和标签 k-mer。
- **处理方法：** 构建 `khmer.GraphLabels` 对象。`reads_path` 中的文件用 `consume_seqfile_and_tag_with_labels`（标签为文件中的记录序号）；`sequences` 用 `consume_sequence_and_tag_with_labels` 与 `labels`（默认取列表中的位置）。`load_labels_path` 先调用 `load_labels_and_tags`；`save_labels_path` 调用 `save_labels_and_tags`（文件必须不存在）。对应工作清单 G1265、G1272 与 G1273。
- **API：** `dnakit.microbiome.khmer_label_reads(*, reads_path=None, sequences=None, labels=None, load_labels_path=None, save_labels_path=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `reads_path`（文件或列表），或 `sequences` 及可选 `labels`；`load_labels_path`、`save_labels_path`；`ksize`、`table_size`、`n_tables`、`max_reads`、`max_bases`。
- **输出：** `n_labels`、`n_tags`、`tag_labels`（每个标签 k-mer 的标签）、`saved_labels`、`loaded_labels`、`graph`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_reads

r = khmer_label_reads(reads_path="lab.fa", ksize=15, verbose=False)
print(r["n_labels"], r["n_tags"], sorted({tuple(v) for v in r["tag_labels"].values()}))
saved = os.path.join(tempfile.mkdtemp(), "labels.bin")
s = khmer_label_reads(
    sequences=["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], labels=[7, 9], ksize=15, save_labels_path=saved, verbose=False
)
print(s["n_labels"], sorted({x for v in s["tag_labels"].values() for x in v}))
print(khmer_label_reads(load_labels_path=saved, ksize=15, verbose=False)["n_labels"])
```

- **示例输出：**

```text
2 4 [(0,), (1,)]
2 [7, 9]
2
```

## 2) `dnakit.microbiome.khmer_label_partitioned_reads()`

- **作用：** 把带分区标注的读段插入带标签图，以分区 ID 作为标签。
- **处理方法：** 对 `khmer_partition_reads` 写出的 `.part` 文件调用 `consume_partitioned_fasta_and_tag_with_labels`。对应工作清单 G1266。
- **API：** `dnakit.microbiome.khmer_label_partitioned_reads(partitioned_path, *, ksize=32, table_size=1000000, n_tables=4)`。
- **输入：** `partitioned_path`；`ksize`、`table_size`、`n_tables`。
- **输出：** `n_labels`、`n_tags`、`tag_labels`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_partition_reads
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_partitioned_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("lab.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_label_partitioned_reads(os.path.join(folder, "lab.fa.part"), ksize=15, verbose=False)
print(r["n_labels"], r["n_tags"])
```

- **示例输出：**

```text
1 4
```

## 3) `dnakit.microbiome.khmer_sweep_neighborhood()`

- **作用：** 在带标签图中查找某条序列邻域内的标签或标签 k-mer。
- **处理方法：** 按 `khmer_label_reads` 的方式构建带标签图，再调用 `sweep_label_neighborhood`（`output="labels"`，返回标签）或 `sweep_tag_neighborhood`（`output="tags"`，返回标签 k-mer）。`sweep_range` 是以 k-mer 计的遍历范围；`break_on_stop_tags` 与 `stop_big_traversals` 是 khmer 的选项。对应工作清单 G1267 与 G1268。
- **API：** `dnakit.microbiome.khmer_sweep_neighborhood(sequence, *, output='labels', sweep_range=0, break_on_stop_tags=False, stop_big_traversals=False, reads_path=None, sequences=None, labels=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`、`output`、`sweep_range`、`break_on_stop_tags`、`stop_big_traversals`；`reads_path`，或 `sequences` 及 `labels`；`ksize`、`table_size`、`n_tables`。
- **输出：** `labels` 或 `tags`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_sweep_neighborhood

query = "AGACTACTTATGAGGCCATACCAGGGTCGT"
r = khmer_sweep_neighborhood(query, reads_path="lab.fa", sweep_range=40, ksize=15, verbose=False)
print(r["labels"])
r = khmer_sweep_neighborhood(
    query, output="tags", reads_path="lab.fa", sweep_range=40, ksize=15, verbose=False
)
print(r["tags"])
```

- **示例输出：**

```text
[0, 1]
['AAAACGACCCTGGTA', 'TACCAGGGTCGTTTG']
```

## 4) `dnakit.microbiome.khmer_tag_labels()`

- **作用：** 查询某个标签 k-mer 所带的标签。
- **处理方法：** 按 `khmer_label_reads` 的方式构建带标签图并调用 `get_tag_labels(tag_kmer)`；该 k-mer 不是标签时结果为空。对应工作清单 G1269。
- **API：** `dnakit.microbiome.khmer_tag_labels(tag_kmer, *, reads_path=None, sequences=None, labels=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `tag_kmer`（长度为 k）；`reads_path`，或 `sequences` 及 `labels`；`ksize`、`table_size`、`n_tables`。
- **输出：** `labels`、`tag_kmer`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_reads, khmer_tag_labels

tags = khmer_label_reads(reads_path="lab.fa", ksize=15, verbose=False)["tag_labels"]
for tag in sorted(tags)[:3]:
    print(tag, khmer_tag_labels(tag, reads_path="lab.fa", ksize=15, verbose=False)["labels"])
print(khmer_tag_labels("A" * 15, reads_path="lab.fa", ksize=15, verbose=False)["labels"])
```

- **示例输出：**

```text
AAAACGACCCTGGTA [1]
ACGACTAGCCAAGCG [1]
CAGCCATAGTCCAAC [0]
[]
```

## 5) `dnakit.microbiome.khmer_label_across_high_degree_nodes()`

- **作用：** 跨越各条重叠群共有的分支点，为它们的 k-mer 打标签。
- **处理方法：** 把所有重叠群放入同一个带标签图，找出高度数 k-mer（`find_high_degree_nodes`），对每条重叠群调用 `label_across_high_degree_nodes(contig, hdn, label)`；标签默认按列表顺序取 1、2、……。返回分支 k-mer 及其获得的标签。对应工作清单 G1270。
- **API：** `dnakit.microbiome.khmer_label_across_high_degree_nodes(contigs, *, labels=None, ksize=32, table_size=1000000, n_tables=4)`。
- **输入：** `contigs`（每条长于 k）、`labels`；`ksize`、`table_size`、`n_tables`。
- **输出：** `n_high_degree`、`high_degree_labels`、`contig_labels`、`n_labels`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_across_high_degree_nodes

r = khmer_label_across_high_degree_nodes(["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], ksize=15, verbose=False)
print(r["n_high_degree"], r["contig_labels"], r["n_labels"])
print(r["high_degree_labels"])
```

- **示例输出：**

```text
1 [1, 2] 2
{'AAACGACCCTGGTAT': [1, 2]}
```

## 6) `dnakit.microbiome.khmer_assemble_labeled_path()`

- **作用：** 从起始 k-mer 出发，利用标签在分支间选择，组装穿过分支图的路径。
- **处理方法：** 按 `khmer_label_across_high_degree_nodes` 的方式给重叠群打标签，再调用 `assemble_labeled_path(start_kmer)`。每条一致的路径返回一条序列，已排序。对应工作清单 G1271。
- **API：** `dnakit.microbiome.khmer_assemble_labeled_path(start_kmer, contigs, *, labels=None, ksize=32, table_size=1000000, n_tables=4)`。
- **输入：** `start_kmer`（长度为 k）、`contigs`、`labels`；`ksize`、`table_size`、`n_tables`。
- **输出：** `n_paths`、`paths`、`contig_labels`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_assemble_labeled_path

r = khmer_assemble_labeled_path("AGACTACTTATGAGG", ["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], ksize=15, verbose=False)
print(r["n_paths"], [len(p) for p in r["paths"]])
print(all(p in ("AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT") for p in r["paths"]))
```

- **示例输出：**

```text
2 [60, 60]
True
```

## 7) `dnakit.microbiome.khmer_merge_hll_counters()`

- **作用：** 通过合并 HyperLogLog 计数器，估计两组读段及其并集的不重复 k-mer 数。
- **处理方法：** 为每个输入填充一个 `khmer.HLLCounter(error_rate, ksize)`（文件用 `consume_seqfile`，字符串用 `consume_string`），记录两个估计值（`estimate_cardinality`），把第二个计数器并入第一个（`merge`），再记录并集的估计值。默认值遵循 khmer 的 `unique-kmers.py`：误差率 0.01、k = 20。估计值带有所述的相对误差。对应工作清单 G1275。
- **API：** `dnakit.microbiome.khmer_merge_hll_counters(*, first_reads_path=None, first_sequences=None, second_reads_path=None, second_sequences=None, error_rate=0.01, ksize=20)`。
- **输入：** `first_reads_path` 或 `first_sequences`、`second_reads_path` 或 `second_sequences`、`error_rate`（0.0001 到 0.367）、`ksize`。
- **输出：** `first_estimate`、`second_estimate`、`merged_estimate`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_merge_hll_counters

r = khmer_merge_hll_counters(
    first_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"],
    second_sequences=["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC"],
    ksize=15,
    verbose=False,
)
print(r["first_estimate"], r["second_estimate"], r["merged_estimate"])
same = khmer_merge_hll_counters(
    first_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], second_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], ksize=15, verbose=False
)
print(same["first_estimate"], same["merged_estimate"])
```

- **示例输出：**

```text
45 45 90
45 45
```

## 8) `dnakit.microbiome.khmer_align_read()`

- **作用：** 用 khmer 的隐马尔可夫模型把读段比对到计数图的 k-mer 图上。
- **处理方法：** 构建或载入计数图，创建 `khmer.ReadAligner(graph, trusted_cov_cutoff, bits_theta)`。`forward_only=False` 调用 `align`，以读段为种子向两个方向延伸；`forward_only=True` 调用 `align_forward`，只向前延伸。`trusted_cov_cutoff`（默认 2）是 k-mer 被视为可信的丰度，`bits_theta`（默认 1）是得分阈值，均为 khmer sandbox 的取值。返回的图比对是图中的路径，因此能体现对测序错误的纠正；缺口用 `-` 表示。对应工作清单 G1276 与 G1277。
- **API：** `dnakit.microbiome.khmer_align_read(read, *, forward_only=False, trusted_cov_cutoff=2, bits_theta=1.0, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `read`、`forward_only`、`trusted_cov_cutoff`、`bits_theta`；`graph_path`、`reads_path` 或 `sequences`；`ksize`、`table_size`、`n_tables`。
- **输出：** `score`、`graph_alignment`、`read_alignment`、`truncated`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_align_read

genome = "TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"
read = genome[5:45]
error = read[:20] + ("A" if read[20] != "A" else "C") + read[21:]
for forward in (False, True):
    r = khmer_align_read(error, sequences=[genome] * 4, ksize=15, forward_only=forward, verbose=False)
    print(forward, round(r["score"], 1), r["graph_alignment"] == read, r["read_alignment"] == error, r["truncated"])
```

- **示例输出：**

```text
False 77.2 True True False
True 77.2 True True False
```

## 9) `dnakit.microbiome.khmer_extract_paired_reads()`

- **作用：** 把混合读段文件拆分为交错的配对读段和孤立读段。
- **处理方法：** 复现 `extract-paired-reads.py`：用 `khmer.ReadParser` 读取文件，用 `khmer.utils.broken_paired_reader` 按名称（`/1` 与 `/2`，或 Casava 的 `1:` 与 `2:`）把相邻读段配对，配对读段交错写入 `paired_output_path`，其余写入 `single_output_path`，格式与输入一致（FASTA 或 FASTQ）。输出文件必须不存在。与脚本一致，没有任何配对的文件会报错。已验证与脚本结果一致。对应工作清单 G1283。
- **API：** `dnakit.microbiome.khmer_extract_paired_reads(reads_path, paired_output_path, single_output_path)`。
- **输入：** `reads_path`、`paired_output_path`、`single_output_path`。
- **输出：** `n_reads`、`n_pairs`、`n_singletons`以及两个输出路径。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_extract_paired_reads

folder = tempfile.mkdtemp()
pairs, singles = os.path.join(folder, "pe.fa"), os.path.join(folder, "se.fa")
r = khmer_extract_paired_reads("mixed.fa", pairs, singles, verbose=False)
print(r["n_reads"], r["n_pairs"], r["n_singletons"])
print(open(singles).read().strip().replace("\n", " "))
```

- **示例输出：**

```text
5 2 1
>s1 TTTTTTTTTTGGAA
```

## 10) `dnakit.microbiome.khmer_split_paired_reads()`

- **作用：** 把交错的配对读段文件拆分为左读段文件和右读段文件。
- **处理方法：** 用 `khmer.utils.broken_paired_reader` 复现 `split-paired-reads.py`：每对读段的第一条写入 `first_output_path`，第二条写入 `second_output_path`。未给 `orphan_output_path` 时，未配对读段会报错（khmer 的 `require_paired`）并删除已写出的文件；给出时，孤立读段写入该文件。输出文件必须不存在。已验证与脚本结果一致。对应工作清单 G1284。
- **API：** `dnakit.microbiome.khmer_split_paired_reads(reads_path, first_output_path, second_output_path, *, orphan_output_path=None)`。
- **输入：** `reads_path`、`first_output_path`、`second_output_path`、`orphan_output_path`。
- **输出：** `n_reads`、`n_pairs`、`n_orphans`以及输出路径。
- **示例代码：**

```python
import os
import tempfile

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_split_paired_reads

folder = tempfile.mkdtemp()
left, right = os.path.join(folder, "left.fa"), os.path.join(folder, "right.fa")
r = khmer_split_paired_reads("paired.fa", left, right, verbose=False)
print(r["n_reads"], r["n_pairs"], r["n_orphans"])
print(open(left).read().strip().replace("\n", " "))
print(open(right).read().strip().replace("\n", " "))
try:
    khmer_split_paired_reads(
        "mixed.fa", os.path.join(folder, "l2.fa"), os.path.join(folder, "r2.fa"), verbose=False
    )
except ConfigurationError:
    print("unpaired read refused")
```

- **示例输出：**

```text
4 2 0
>p1/1 ACGTACGTACGGTT >p2/1 AAAACCCCGGTTAC
>p1/2 GGGGCCCCAATTGC >p2/2 CCCCAAAATTGGCA
unpaired read refused
```

## 11) `dnakit.microbiome.khmer_read_statistics()`

- **作用：** FASTA/FASTQ 文件的碱基数、序列数和平均序列长度。
- **处理方法：** 用 `khmer.ReadParser` 复现 `readstats.py`：逐文件及合计统计碱基数与序列数；平均长度按脚本的输出方式保留一位小数。khmer 无法打开或不含序列的文件，报告为序列数 0，并在 `error` 中给出解析器信息，且与脚本一致，不计入合计。`output_csv` 写出脚本 `--csv` 的表（bp、seqs、avg_len、filename）。已验证与脚本结果一致。对应工作清单 G1286。
- **API：** `dnakit.microbiome.khmer_read_statistics(reads_paths, *, output_csv=None)`。
- **输入：** `reads_paths`（文件或列表）、`output_csv`。
- **输出：** `files`（`filename`、`bp`、`seqs`、`avg_len`、`error`）、`total`、`output_csv`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_read_statistics

csv_path = os.path.join(tempfile.mkdtemp(), "stats.csv")
r = khmer_read_statistics(["lab.fa", "mixed.fa"], output_csv=csv_path, verbose=False)
for row in r["files"]:
    print(os.path.basename(row["filename"]), row["bp"], row["seqs"], row["avg_len"])
print(r["total"])
print(open(csv_path).read().splitlines()[0])
```

- **示例输出：**

```text
lab.fa 120 2 60.0
mixed.fa 70 5 14.0
{'bp': 190, 'seqs': 7, 'avg_len': 27.1}
bp,seqs,avg_len,filename
```
