# khmer 标签、分区与停止标签：把读段拆分为连通分量

本页对应工作清单 P27 包中关于标签、分区和停止标签的 khmer 条目。khmer 在读段的 de Bruijn 图中标记一组稀疏的标签 k-mer，把相连的标签归为分区，并为每条读段标注分区 ID（`do-partition.py` 流程）。高度连通的 k-mer（结）可以保存为停止标签，分区在这些位置被切开。读段是由 khmer 自己读取的本地 FASTA/FASTQ 文件，因此适用 khmer 的读段清洗（转大写、N 按 A 处理、短于 k 的读段被跳过）。

仅参数不同或只是脚本步骤不同的条目合并：载入、保存、打印和扩充标签集是 `khmer_tag_reads` 的选项；把标签分成子集、对子集分区（含按丰度分区）、合并子集映射、保存与载入分区映射、统计分区、子集大小分布、查询分区 ID 以及标注读段，都是 `khmer_partition_reads` 的选项，它对应 `do-partition.py`（一次调用完成 `load-graph.py`、`partition-graph.py`、`merge-partitions.py` 与 `annotate-partitions.py`）；寻找结、载入、保存、打印和添加停止标签是 `khmer_find_knots` 的选项（`make-initial-stoptags.py` 与 `find-knots.py`）。注意 khmer 的 `subset_count_partitions` 会破坏子集映射，所以函数在映射保存并合并之后才读取它。

khmer 是可选依赖（`pip install "dnakit[bio-khmer]"`）；已用 khmer 2.1.1 测试。默认值遵循 khmer：k = 32、4 张表、每张表 1e6 个桶、子集大小 1e5。示例使用 k = 15。

## 示例文件

`reads.fa` 含三条读段：`a` 与 `c` 的前 30 个碱基相同，所以落在同一个分区；`b` 与它们无关。

- **示例文件：** `reads.fa`

```text
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

## 1) `dnakit.microbiome.khmer_tag_reads()`

- **作用：** 把读段插入 khmer 图并标记标签 k-mer；可载入、扩充、保存并打印标签集。
- **处理方法：** 由 `reads_path` 和/或 `sequences` 构建图并打标签（`consume_seqfile_and_tag`、`consume_and_tag`）。可选的 khmer 调用：`load_tagset_path`（`load_tagset`）、`extra_tag_kmers`（`add_tag`）、`save_tagset_path`（`save_tagset`）、`tagset_text_path`（`print_tagset`，每行一个 k-mer）。输出文件必须不存在。对应工作清单 G1239、G1243、G1244、G1245 与 G1246。
- **API：** `dnakit.microbiome.khmer_tag_reads(*, reads_path=None, sequences=None, load_tagset_path=None, extra_tag_kmers=None, save_tagset_path=None, tagset_text_path=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `reads_path`（文件或列表）和/或 `sequences`；`load_tagset_path`、`extra_tag_kmers`、`save_tagset_path`、`tagset_text_path`；`graph_type`、`ksize`、`table_size`、`n_tables`、`max_reads`、`max_bases`。
- **输出：** `n_tags`、`tags`、`n_unique_kmers`、`extra_tags_added`、`saved_tagset`、`tagset_text`、`graph`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_tag_reads

folder = tempfile.mkdtemp()
binary = os.path.join(folder, "tags.bin")
text = os.path.join(folder, "tags.txt")
r = khmer_tag_reads(
    reads_path="reads.fa",
    ksize=15,
    extra_tag_kmers=["A" * 15],
    save_tagset_path=binary,
    tagset_text_path=text,
    verbose=False,
)
print(r["n_tags"], r["extra_tags_added"], "A" * 15 in r["tags"])
print(len(open(text).read().split()))
print(khmer_tag_reads(load_tagset_path=binary, ksize=15, verbose=False)["n_tags"])
```

- **示例输出：**

```text
7 1 True
7
7
```

## 2) `dnakit.microbiome.khmer_tag_positions()`

- **作用：** 序列上各标签 k-mer 的位置。
- **处理方法：** 由 `reads_path` 和/或 `sequences` 构建并标记图，然后调用 `get_tags_and_positions(sequence)`。位置是查询序列中从 0 开始的起始偏移。对应工作清单 G1240。
- **API：** `dnakit.microbiome.khmer_tag_positions(sequence, *, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`；`reads_path` 和/或 `sequences`；`graph_type`、`ksize`、`table_size`、`n_tables`。
- **输出：** `tags`：由 `position` 与 `kmer` 组成的列表。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_tag_positions

r = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)
for tag in r["tags"]:
    print(tag["position"], tag["kmer"])
```

- **示例输出：**

```text
19 AACGCCACGTAACTA
46 ACCACGACGGCAGGC
```

## 3) `dnakit.microbiome.khmer_reachable_tags()`

- **作用：** 查找与某个 k-mer 相连的标签，并可为它们指定同一个分区 ID。
- **处理方法：** 构建并标记图，然后调用 `find_all_tags_list(kmer)`。`assign_partition=True` 时还调用 `find_all_tags` 与 `assign_partition_id`：若其中某个标签已有分区则沿用其 ID，否则分配新 ID。对应工作清单 G1241 与 G1249。
- **API：** `dnakit.microbiome.khmer_reachable_tags(kmer, *, assign_partition=False, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `kmer`（长度为 k）、`assign_partition`；`reads_path` 和/或 `sequences`；`graph_type`、`ksize`、`table_size`、`n_tables`。
- **输出：** `n_tags`、`tags`、`partition_id`（未设 `assign_partition` 时为 None）。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_reachable_tags, khmer_tag_positions

tags = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["tags"]
start = tags[0]["kmer"]
r = khmer_reachable_tags(start, reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_tags"], r["partition_id"])
r = khmer_reachable_tags(
    start, reads_path="reads.fa", ksize=15, assign_partition=True, verbose=False
)
print(r["partition_id"])
```

- **示例输出：**

```text
2 None
2
```

## 4) `dnakit.microbiome.khmer_unique_paths()`

- **作用：** 序列中其 k-mer 大多不在图里的片段。
- **处理方法：** 由 `reads_path` 和/或 `sequences` 构建图（都没给时为空图），调用 `extract_unique_paths(sequence, min_length, min_unique_fraction)`。`min_length` 默认为 k，且不能小于 k。khmer 的文档把该方法标记为将移除；它就是工作清单 G1242 所指的 khmer 操作。
- **API：** `dnakit.microbiome.khmer_unique_paths(sequence, *, min_length=None, min_unique_fraction=1.0, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `sequence`、`min_length`、`min_unique_fraction`；`reads_path` 和/或 `sequences`；`graph_type`、`ksize`、`table_size`、`n_tables`。
- **输出：** `paths`（序列列表）、`min_length`、`min_unique_fraction`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_unique_paths

print(len(khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", ksize=15, verbose=False)["paths"][0]))
print(khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["paths"])
r = khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], ksize=15, verbose=False)
print([len(p) for p in r["paths"]])
```

- **示例输出：**

```text
60
[]
[44]
```

## 5) `dnakit.microbiome.khmer_partition_reads()`

- **作用：** 把读段划分为其 de Bruijn 图的连通分量，并给每条读段标注分区 ID。
- **处理方法：** 遵循 khmer 的 `do-partition.py`：`consume_seqfile_and_tag`、`divide_tags_into_subsets`（`subset_size`）、对每个子集 `do_subset_partition`、`merge_subset`，最后 `output_partitions`，写出 `<name>.part`，读段名后以制表符附上分区 ID。各选项都对应一次 khmer 调用：`stoptags_path`（`load_stop_tags`，分区在那里被切开）、`stop_big_traversals`、`abundance_range=[min, max]`（在计数图上 `do_subset_partition_with_abundance`）、`save_partitionmap_dir`（`save_subset_partitionmap`，文件 `<name>.subset.N.pmap`）、`partitionmap_paths`（跳过分区，用 `merge_subset_from_disk` 合并已保存的映射）、`query_kmers`（`get_partition_id`，对标签 k-mer 才有意义）。`annotate=False` 不写 `.part` 文件。分区 ID 0 表示未分配。对应工作清单 G1230、G1247、G1248、G1250、G1251、G1252、G1254、G1255、G1256、G1258 与 G1280。
- **API：** `dnakit.microbiome.khmer_partition_reads(reads_path, *, output_dir=None, subset_size=100000, stop_big_traversals=False, stoptags_path=None, abundance_range=None, partitionmap_paths=None, save_partitionmap_dir=None, query_kmers=None, annotate=True, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `reads_path`（文件或列表）、`output_dir`、`subset_size`、`stop_big_traversals`、`stoptags_path`、`abundance_range`、`partitionmap_paths`、`save_partitionmap_dir`、`query_kmers`、`annotate`、`ksize`、`table_size`、`n_tables`。
- **输出：** `partition_count`（`n_partitions`、`n_unassigned_tags`；合并已保存映射时为 None）、`n_subsets`、`per_subset`（各子集的大小分布）、`annotated`（输出文件及分区数）、`saved_partitionmaps`、`partition_ids`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_partition_reads

folder = tempfile.mkdtemp()
r = khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, subset_size=3, verbose=False)
print(r["partition_count"], r["n_subsets"])
print(r["annotated"][0]["n_partitions"])
for line in open(os.path.join(folder, "reads.fa.part")):
    if line.startswith(">"):
        print(line.rstrip().replace("\t", " -> partition "))
```

- **示例输出：**

```text
{'n_partitions': 2, 'n_unassigned_tags': 0} 2
2
>a -> partition 2
>b -> partition 3
>c -> partition 2
```

## 6) `dnakit.microbiome.khmer_join_partitions()`

- **作用：** 先对读段分区，再把两个标签 k-mer 所在的分区合并。
- **处理方法：** 按 `khmer_partition_reads` 的方式对读段分区，检查两个 k-mer 都是标签，用 `get_partition_id` 查询它们的分区 ID，用 `join_partitions` 合并，并报告合并前后的 ID 与分区数。对应工作清单 G1257。
- **API：** `dnakit.microbiome.khmer_join_partitions(reads_path, first_tag_kmer, second_tag_kmer, *, subset_size=100000, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** `reads_path`（文件或列表）、`first_tag_kmer`、`second_tag_kmer`；`subset_size`、`ksize`、`table_size`、`n_tables`。
- **输出：** `partition_ids_before`、`partition_ids_after`、`partition_count_before`、`partition_count_after`。
- **示例代码：**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_join_partitions, khmer_tag_positions

tag_a = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["tags"][0]["kmer"]
tag_b = khmer_tag_positions("CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA", reads_path="reads.fa", ksize=15, verbose=False)["tags"][0]["kmer"]
r = khmer_join_partitions("reads.fa", tag_a, tag_b, ksize=15, verbose=False)
print(r["partition_ids_before"], r["partition_count_before"])
print(r["partition_ids_after"], r["partition_count_after"])
```

- **示例输出：**

```text
[2, 3] [2, 0]
[2, 2] [1, 0]
```

## 7) `dnakit.microbiome.khmer_load_partitioned_reads()`

- **作用：** 把带分区标注的读段插入新图并重建其分区。
- **处理方法：** 对 `.part` 文件调用 `consume_partitioned_fasta`：每条读段被计数并打标签，其标签获得读段名中写的分区 ID；随后 `count_partitions`。返回各标签及其分区 ID；`query_kmers` 查询其他 k-mer（`get_partition_id`，除非该 k-mer 是重建后的标签，否则为 0）。对应工作清单 G1253。
- **API：** `dnakit.microbiome.khmer_load_partitioned_reads(partitioned_path, *, query_kmers=None, ksize=32, table_size=1000000, n_tables=4)`。
- **输入：** `partitioned_path`、`query_kmers`、`ksize`、`table_size`、`n_tables`。
- **输出：** `n_reads`、`n_kmers_consumed`、`n_partitions`、`n_unassigned_tags`、`n_tags`、`tag_partition_ids`、`partition_ids`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_load_partitioned_reads, khmer_partition_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_load_partitioned_reads(os.path.join(folder, "reads.fa.part"), ksize=15, verbose=False)
print(r["n_reads"], r["n_kmers_consumed"], r["n_partitions"], r["n_tags"])
print(sorted(set(r["tag_partition_ids"].values())))
```

- **示例输出：**

```text
3 138 2 2
[2, 3]
```

## 8) `dnakit.microbiome.khmer_find_knots()`

- **作用：** 通过对最大分区重新分区找出高度连通的 k-mer（结），并作为停止标签返回。
- **处理方法：** `stage="initial"` 遵循 `make-initial-stoptags.py`：对第一个子集（默认 10000 个标签）分区，并以计数阈值 5 调用 `repartition_largest_partition`。`stage="knots"` 遵循 `find-knots.py`：对所有子集（默认 100000 个标签）分区，逐个重新分区（计数阈值 2）、合并，再对合并后的映射重新分区。`excursion_distance`（40）与 `excursion_kmer_threshold`（200）是 khmer 的取值；游走用的计数图有 `counting_table_size` 个桶。停止标签选项：`load_stoptags_path`（`load_stop_tags`）、`extra_stop_kmers`（`add_stop_tag`）、`save_stoptags_path`（`save_stop_tags`）、`stoptags_text_path`（`print_stop_tags`）。已验证与两个脚本的结果一致。对应工作清单 G1259、G1260、G1261、G1262 与 G1264。
- **API：** `dnakit.microbiome.khmer_find_knots(reads_path, *, stage='knots', subset_size=None, excursion_distance=40, excursion_kmer_threshold=200, excursion_kmer_count_threshold=None, load_stoptags_path=None, extra_stop_kmers=None, save_stoptags_path=None, stoptags_text_path=None, ksize=32, table_size=1000000, n_tables=4, counting_table_size=1000000, max_reads=1000000, max_bases=100000000)`。
- **输入：** `reads_path`、`stage`、`subset_size`、`excursion_distance`、`excursion_kmer_threshold`、`excursion_kmer_count_threshold`、停止标签选项、`ksize`、`table_size`、`n_tables`、`counting_table_size`。
- **输出：** `n_stop_tags`、`stop_tags`、所用阈值、`saved_stoptags`、`stoptags_text`。
- **示例代码：**

```python
import os
import random
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_find_knots

rng = random.Random(11)


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


repeat = dna(100)
folder = tempfile.mkdtemp()
reads = os.path.join(folder, "knot.fa")
with open(reads, "w") as handle:
    for i in range(40):
        handle.write(f">k{i}\n{dna(120)}{repeat}{dna(120)}\n")
stoptags = os.path.join(folder, "knot.stoptags")
r = khmer_find_knots(reads, stage="initial", ksize=21, save_stoptags_path=stoptags, verbose=False)
print(r["stage"], r["n_stop_tags"], r["excursion_kmer_count_threshold"])
r = khmer_find_knots(reads, stage="knots", ksize=21, load_stoptags_path=stoptags, verbose=False)
print(r["stage"], r["n_stop_tags"], r["excursion_kmer_count_threshold"])
```

- **示例输出：**

```text
initial 1572 5
knots 1572 2
```

## 9) `dnakit.microbiome.khmer_filter_stoptags()`

- **作用：** 在停止标签处修剪读段，并写出保留的读段。
- **处理方法：** 遵循 `filter-stoptags.py`：载入停止标签（`stoptags_path` 用 `load_stop_tags`，和/或 `stop_kmers` 用 `add_stop_tag`），丢弃含 N 的读段，用 `trim_on_stoptags` 在第一个停止标签处截断其余读段，仅当剩余至少 k 个碱基时保留。输出为 FASTA（`>name` 与修剪后的序列）。含 A、C、G、T 以外符号的读段被丢弃。对应工作清单 G1263。
- **API：** `dnakit.microbiome.khmer_filter_stoptags(reads_path, output_path, *, stoptags_path=None, stop_kmers=None, ksize=32)`。
- **输入：** `reads_path`、`output_path`（必须不存在）、`stoptags_path` 和/或 `stop_kmers`、`ksize`。
- **输出：** `n_reads`、`n_kept`、`n_trimmed`、`n_dropped_with_n`、`output_path`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_filter_stoptags

out = os.path.join(tempfile.mkdtemp(), "filtered.fa")
r = khmer_filter_stoptags("reads.fa", out, stop_kmers=["GTTGCTCCACAGGTA"], ksize=15, verbose=False)
print(r["n_reads"], r["n_kept"], r["n_trimmed"])
for line in open(out):
    print(line.rstrip())
```

- **示例输出：**

```text
3 3 1
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGG
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

## 10) `dnakit.microbiome.khmer_extract_partitions()`

- **作用：** 按分区大小对带分区标注的读段分组，并写出分区大小分布。
- **处理方法：** 复现 `extract-partitions.py`。读取 `.part` 文件，统计每个分区的读段数，写出 `<prefix>.dist`，列为读段数、含该读段数的分区数、累计分区数和累计读段数。`write_groups=True`（即不带 `-n` 的脚本）时，读段数大于 `min_partition_size` 的分区按从小到大打包进约 `max_group_size` 条读段的 `<prefix>.groupNNNN.fa`（或 `.fq`）；`write_unassigned` 还会写出分区 0 的读段。已验证与脚本结果一致。对应工作清单 G1256（`-n` 选项）与 G1281。
- **API：** `dnakit.microbiome.khmer_extract_partitions(partitioned_paths, output_prefix, *, write_groups=True, min_partition_size=5, max_group_size=1000000, write_unassigned=False)`。
- **输入：** `partitioned_paths`（文件或列表）、`output_prefix`、`write_groups`、`min_partition_size`（5）、`max_group_size`（1000000）、`write_unassigned`。
- **输出：** `n_partitions`、`n_unassigned_reads`、`size_distribution`、`distribution_path`、`group_files`、`unassigned_path`。
- **示例代码：**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_extract_partitions, khmer_partition_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_extract_partitions(
    os.path.join(folder, "reads.fa.part"),
    os.path.join(folder, "groups"),
    min_partition_size=0,
    verbose=False,
)
print(r["n_partitions"], r["n_unassigned_reads"])
print(open(r["distribution_path"]).read().strip())
print(len(r["group_files"]), open(r["group_files"][0]).read().count(">"))
```

- **示例输出：**

```text
2 0
1 1 1 1
2 1 2 3
1 3
```
