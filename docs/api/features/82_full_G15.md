# 草图与 k-mer 辅助函数（sourmash 与 khmer）

这组函数是 `dnakit.microbiome` 中对 [sourmash](https://github.com/sourmash-bio/sourmash) 与 [khmer](https://github.com/dib-lab/khmer) 的 DNA 级包装。前八个函数处理 sourmash 的 scaled MinHash 草图：两两比较矩阵、突变率估计、草图大小误差概率、k-mer 哈希、降采样、交集、数据库搜索和谱系最低共同祖先。后八个处理 khmer 的 k-mer 哈希与 2 位编码、本地 reads 的低丰度修剪、已保存计数表的检查与查询、表大小估计以及假阳性率估计。这些函数只使用调用时给出的文件或草图，不下载数据库，也不做组装或物种分类。

khmer 对小写碱基会静默给出错误结果，因此这里的函数在调用 khmer 之前先把序列转为大写，并把 k 限制在 1 到 32。示例从定义函数的模块导入，即 `dnakit.microbiome` 与 `dnakit.microbiome`。依赖都是可选的，只有调用时才导入：sourmash 对应 `bio-sourmash`，khmer 对应 `bio-khmer`；缺少包时会报错并给出安装命令。

## 示例数据

下面的 sourmash 示例使用 `sample.fa`。两条基因组只在一个碱基上不同，第三条是短 read。khmer 示例使用 `reads.fa`：三条相同 read 构成高覆盖区域，`tailed` 在其后带有一段新序列，`novel` 的 k-mer 只出现一次，`nonacgt` 含非 ACGT 符号，`short` 比 k 短。两个文件只用于演示接口。

**示例文件：** `sample.fa`

```text
>genome_a
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGATCCATGGCA
>genome_b
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGTTCCATGGCA
>read_c
TTTTGGGGCCCCAAAAGGGGTTTTACGTACGTAAAA
```

**示例文件：** `reads.fa`

```text
>ref1
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>ref2
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>ref3
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>tailed
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGGGGGGGGGGGGGCCCCCCCCCCCC
>novel
CAGTCGATGCATGCAA
>nonacgt
ACGTNNNNACGTACGTACGT
>short
ACGT
```

## 1) `dnakit.microbiome.compare_minhash_matrix()`

- **作用：** 计算一组 scaled MinHash 草图两两之间的相似度矩阵。
- **处理方法：** `metric="jaccard"` 调用 `sourmash.compare.compare_serial`；`containment`、`max_containment` 和 `avg_containment` 分别调用 `compare_serial_containment`、`compare_serial_max_containment` 和 `compare_serial_avg_containment`。sourmash 把包含度矩阵存为转置形式，本函数转置之后使 `matrix[i][j]` 表示草图 i 中被草图 j 包含的比例。
- **API：** `dnakit.microbiome.compare_minhash_matrix(sketches, *, metric="jaccard")`。
- **输入：** 至少两个草图字典（含 `ksize`、`scaled`、`hashes`，可选 `id`），各草图的 `ksize` 与 `scaled` 必须相同；最多 500 个草图。
- **输出：** `ids`、`matrix`（按输入顺序排列的方阵）、`metric`、`ksize` 和 `scaled`。
- **示例代码：**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import compare_minhash_matrix

sketches = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][:2]
result = compare_minhash_matrix(sketches)
print(result["ids"], [round(value, 3) for value in result["matrix"][0]])
```

- **示例输出：**

```text
['genome_a', 'genome_b'] [1.0, 0.5]
```

## 2) `dnakit.microbiome.estimate_mutation_rate()`

- **作用：** 由草图的 Jaccard 指数或包含度估计每碱基突变率。
- **处理方法：** `measure="jaccard"` 调用 `sourmash.minhash.jaccard_to_distance`；`measure="containment"` 调用 `containment_to_distance` 并给出置信区间。需要给出序列长度（bp）或不同 k-mer 数之一。`ani` 等于 1 减去 `distance`。
- **API：** `dnakit.microbiome.estimate_mutation_rate(value, *, measure="jaccard", ksize=21, scaled=1000, n_unique_kmers=None, sequence_len_bp=None, confidence=0.95)`。
- **输入：** 开区间 (0, 1) 内的 `value`；`sequence_len_bp` 或 `n_unique_kmers` 至少给一个。这是近似估计，结果中的 `size_is_inaccurate` 提示草图可能太小。
- **输出：** `distance`、`ani`、包含度模式下的 `distance_low` 与 `distance_high`、`p_nothing_in_common` 和 `size_is_inaccurate`。
- **示例代码：**

```python
from dnakit.microbiome import estimate_mutation_rate

result = estimate_mutation_rate(0.8, ksize=21, scaled=10, n_unique_kmers=1000)
print(round(result["distance"], 4), round(result["ani"], 4), result["size_is_inaccurate"])
```

- **示例输出：**

```text
0.0056 0.9944 False
```

## 3) `dnakit.microbiome.set_size_error_probability()`

- **作用：** 计算 scaled 草图对 k-mer 集合大小的估计落在相对误差容许范围内的概率。
- **处理方法：** 调用 `sourmash.minhash.set_size_exact_prob`。草图大小服从参数为 `set_size` 与 `1 / scaled` 的二项分布。返回值是 P(|草图大小 × scaled − set_size| ≤ relative_error × set_size)。上游的文档字符串描述的是其补事件，但计算出的数值是这里所述的概率。
- **API：** `dnakit.microbiome.set_size_error_probability(set_size, scaled, *, relative_error=0.05)`。
- **输入：** `set_size` 为正整数；`scaled` 为 1 到 10^9 的整数；`relative_error` 在 0 与 1 之间（开区间）。
- **输出：** `expected_sketch_size`（即 `set_size / scaled`）和 `probability_within_tolerance`。
- **示例代码：**

```python
from dnakit.microbiome import set_size_error_probability

result = set_size_error_probability(1000, 10, relative_error=0.1)
print(result["expected_sketch_size"], round(result["probability_within_tolerance"], 4))
```

- **示例输出：**

```text
100.0 0.7318
```

## 4) `dnakit.microbiome.hash_canonical_kmer()`

- **作用：** 计算单个 DNA k-mer 的规范哈希，与 scaled=1 草图中存储的哈希一致。
- **处理方法：** 规范形式取 k-mer 与其反向互补序列中字符串较小者，再调用 `sourmash.minhash.hash_murmur(canonical, seed)`。小写输入先转为大写。
- **API：** `dnakit.microbiome.hash_canonical_kmer(kmer, *, seed=42)`。
- **输入：** 只含 A、C、G、T 的字符串，长度 1 到 100。
- **输出：** `canonical`、`strand`（forward、reverse_complement 或 self_reverse_complement）和 `hash`。
- **示例代码：**

```python
from dnakit.microbiome import hash_canonical_kmer

forward = hash_canonical_kmer("acgta")
reverse = hash_canonical_kmer("TACGT")
print(forward["canonical"], reverse["strand"], forward["hash"] == reverse["hash"])
```

- **示例输出：**

```text
ACGTA reverse_complement True
```

## 5) `dnakit.microbiome.downsample_minhash_sketch()`

- **作用：** 把 scaled 草图降采样到更粗的 scale。
- **处理方法：** 调用 `sourmash.minhash.flatten_and_downsample_scaled`。目标 `scaled` 不能小于草图自身的 scale，草图不能升采样。超过 sourmash 阈值的哈希会被丢弃；阈值由浮点运算得到，所以在 `2**64 / scaled` 附近几百个单位内的哈希，其归属以 sourmash 为准。
- **API：** `dnakit.microbiome.downsample_minhash_sketch(sketch, *, scaled)`。
- **输入：** 一个草图字典；`scaled` 必须大于等于草图的 scale。
- **输出：** 降采样后的 `scaled`、`hashes`、`input_hash_count` 和 `hash_count`。
- **示例代码：**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import downsample_minhash_sketch

sketch = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][0]
result = downsample_minhash_sketch(sketch, scaled=100)
print(sketch["hash_count"], result["hash_count"], result["scaled"])
```

- **示例输出：**

```text
5 1 100
```

## 6) `dnakit.microbiome.intersect_minhash_sketches()`

- **作用：** 返回两个草图共有的哈希。
- **处理方法：** 调用 `sourmash.minhash.flatten_and_intersect_scaled`：先把两者降采样到较粗的 scale，再取哈希交集。
- **API：** `dnakit.microbiome.intersect_minhash_sketches(left, right)`。
- **输入：** 两个草图字典，`ksize` 必须相同；两者的 scale 可以不同。
- **输出：** `hashes`、`hash_count` 和 `scaled`（较粗的那个 scale）。
- **示例代码：**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import intersect_minhash_sketches

first, second = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][:2]
result = intersect_minhash_sketches(first, second)
print(result["hash_count"], result["scaled"])
```

- **示例输出：**

```text
3 10
```

## 7) `dnakit.microbiome.search_minhash_database()`

- **作用：** 在草图列表中找出与查询草图共享至少 `threshold_bp` 个碱基的草图。
- **处理方法：** 用 `sourmash.index.LinearIndex` 装载数据库，调用 `sourmash.search.prefetch_database`。bp 重叠等于共享哈希数乘以较粗的 scale。结果按重叠 bp 降序排列，相同时按 id 排序。阈值超过查询草图自身的 bp 时，sourmash 会拒绝，本函数抛出 ConfigurationError。
- **API：** `dnakit.microbiome.search_minhash_database(query, database, *, threshold_bp=5000)`。
- **输入：** 一个查询草图和一个数据库草图列表；所有草图的 `ksize` 相同。空的草图会报错，需先过滤。
- **输出：** `matches` 列表，每项含 `id`、`md5`、`intersect_bp`、`jaccard`、`fraction_of_query_in_match`、`fraction_of_match_in_query` 和 `max_containment`。
- **示例代码：**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import search_minhash_database

sketches = [s for s in sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"]
            if s["hash_count"]]
result = search_minhash_database(sketches[0], sketches, threshold_bp=30)
print([(match["id"], match["intersect_bp"]) for match in result["matches"]])
```

- **示例输出：**

```text
[('genome_a', 50), ('genome_b', 30)]
```

## 8) `dnakit.microbiome.lowest_common_lineage()`

- **作用：** 求一组分类谱系的最低共同祖先（LCA）。
- **处理方法：** 用 `sourmash.lca.lca_utils.make_lineage` 解析，`build_tree` 建树，`find_lca` 求最低共同节点。`distinct_children_below_lca` 为 0 表示所有谱系在完整深度上一致；否则为下一层不同分类单元的数目。
- **API：** `dnakit.microbiome.lowest_common_lineage(lineages)`。
- **输入：** 以分号分隔、从界到株排列的谱系字符串列表，最多 8 层；不接受逗号分隔。
- **输出：** `lineage`（每层的 rank 与 name）、`lineage_string`、`lowest_rank` 和 `distinct_children_below_lca`。
- **示例代码：**

```python
from dnakit.microbiome import lowest_common_lineage

result = lowest_common_lineage(
    [
        "d__Bacteria;p__Firmicutes;c__Bacilli",
        "d__Bacteria;p__Firmicutes;c__Clostridia",
    ]
)
print(result["lineage_string"], result["distinct_children_below_lca"])
```

- **示例输出：**

```text
d__Bacteria;p__Firmicutes 2
```

## 9) `dnakit.microbiome.hash_kmer_murmur3()`

- **作用：** 计算 DNA k-mer 的 khmer MurmurHash3 哈希。
- **处理方法：** `canonical=True` 调用 `khmer.hash_murmur3`，正向与反向互补得到相同的值；`canonical=False` 调用 `khmer.hash_no_rc_murmur3`。khmer 对小写碱基会静默给出错误哈希，因此先转为大写。khmer 的哈希值与 sourmash 的不同。
- **API：** `dnakit.microbiome.hash_kmer_murmur3(kmer, *, canonical=True)`。
- **输入：** 只含 A、C、G、T 的字符串，长度 1 到 32。
- **输出：** `hash`（无符号 64 位整数）、`canonical` 和 `ksize`。
- **示例代码：**

```python
from dnakit.microbiome import hash_kmer_murmur3

print(hash_kmer_murmur3("acgta")["hash"] == hash_kmer_murmur3("TACGT")["hash"], hash_kmer_murmur3("ACGTA")["hash"])
```

- **示例输出：**

```text
True 7467422729906674866
```

## 10) `dnakit.microbiome.encode_kmer_2bit()`

- **作用：** 把 k-mer（k 最大 32）编码为精确的 2 位整数。
- **处理方法：** `canonical=True` 调用 `khmer.forward_hash`，取 k-mer 与其反向互补编码中较小者；`canonical=False` 调用 `khmer.forward_hash_no_rc`。khmer 的碱基顺序为 A=0、T=1、C=2、G=3，首碱基位于高位。
- **API：** `dnakit.microbiome.encode_kmer_2bit(kmer, *, canonical=True)`。
- **输入：** 只含 A、C、G、T 的字符串，长度 1 到 32。
- **输出：** `code` 和 `canonical`、`ksize`。
- **示例代码：**

```python
from dnakit.microbiome import encode_kmer_2bit

print(
    encode_kmer_2bit("ACGT")["code"],
    encode_kmer_2bit("GGGGA")["code"],
    encode_kmer_2bit("GGGGA", canonical=False)["code"],
)
```

- **示例输出：**

```text
45 426 1020
```

## 11) `dnakit.microbiome.decode_kmer_2bit()`

- **作用：** 由 `encode_kmer_2bit()` 的 2 位编码还原 k-mer。
- **处理方法：** 调用 `khmer.reverse_hash`。规范编码还原出的是规范 k-mer。取值必须在 0 到 4^k − 1 之间。
- **API：** `dnakit.microbiome.decode_kmer_2bit(value, *, ksize)`。
- **输入：** 整数 `value` 与 `ksize`（1 到 32）。
- **输出：** 大写的 `kmer`。
- **示例代码：**

```python
from dnakit.microbiome import decode_kmer_2bit

print(decode_kmer_2bit(45, ksize=4)["kmer"])
```

- **示例输出：**

```text
ACGT
```

## 12) `dnakit.microbiome.trim_low_abundance_reads()`

- **作用：** 对本地 reads 做低丰度修剪：在第一个计数低于 `cutoff` 的 k-mer 处截断。
- **处理方法：** 第一遍用 khmer `Countgraph` 统计所有可用 reads 的规范 k-mer；第二遍对每条 read 调用 `khmer.trimming.trim_record`，保留第一个低丰度 k-mer（从 5' 端数）之前的碱基。截断后短于 k 的 read 被丢弃。输出为 FASTA（`.fa`、`.fasta`、`.fna`）或 FASTQ（`.fq`、`.fastq`，且需要 FASTQ 输入）。输出碱基为大写。
- **API：** `dnakit.microbiome.trim_low_abundance_reads(source, output, *, ksize=21, cutoff=2, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** FASTA 或 FASTQ 文件；输出路径不得已存在。含 A、C、G、T 以外符号的 read 和短于 k 的 read 会被跳过并计数。
- **输出：** `reads`、`used_reads`、`kept_reads`、`trimmed_reads`、`dropped_short_reads` 等计数，以及输出文件路径。
- **示例代码：**

```python
from pathlib import Path
from dnakit.microbiome import trim_low_abundance_reads

Path("trimmed.fa").unlink(missing_ok=True)
result = trim_low_abundance_reads("reads.fa", "trimmed.fa", ksize=9, cutoff=2)
print(
    result["reads"],
    result["used_reads"],
    result["trimmed_reads"],
    result["kept_reads"],
    result["dropped_short_reads"],
)
```

- **示例输出：**

```text
7 5 1 4 1
```

## 13) `dnakit.microbiome.inspect_khmer_graph()`

- **作用：** 读取已保存 khmer 计数表或存在表的文件头，不加载整张表。
- **处理方法：** `kind="countgraph"` 调用 `khmer.extract_countgraph_info`；`kind="nodegraph"` 调用 `khmer.extract_nodegraph_info` 并映射到具名字段。khmer 不校验文件类型，因此本函数检查类型标记，并拒绝另一种文件。
- **API：** `dnakit.microbiome.inspect_khmer_graph(source, *, kind="countgraph")`。
- **输入：** khmer 保存的 `.cg` 或 `.ng` 文件。
- **输出：** `ksize`、`table_size`、`n_tables`、`n_occupied`、`format_version` 和 `table_type_flag` 等字段。
- **示例代码：**

```python
import khmer
from dnakit.microbiome import inspect_khmer_graph

table = khmer.Countgraph(7, 1009, 3)
table.consume("ACGTTGCAAGCTAGCTTTG")
table.save("counts.cg")
info = inspect_khmer_graph("counts.cg")
print(info["ksize"], info["n_tables"], info["table_size"], info["n_occupied"])
```

- **示例输出：**

```text
7 3 1000 11
```

## 14) `dnakit.microbiome.query_saved_countgraph()`

- **作用：** 查询已保存 khmer 计数表中指定 k-mer 的规范计数。
- **处理方法：** 调用 `khmer.load_countgraph` 加载，再对每个查询调用 `get()`。正向与反向互补共用一个计数；哈希冲突只会使计数偏大。查询 k-mer 的长度必须等于表的 ksize。
- **API：** `dnakit.microbiome.query_saved_countgraph(source, kmers)`。
- **输入：** khmer 计数表文件，以及长度相同的 k-mer 列表（最多 100000 个）。
- **输出：** 按 k-mer 排序的 `counts`，以及 `ksize` 和 `hashsizes`。
- **示例代码：**

```python
import khmer
from dnakit.microbiome import query_saved_countgraph

table = khmer.Countgraph(7, 1009, 3)
table.consume("ACGTTGCAAGCTAGCTTTG")
table.save("counts.cg")
result = query_saved_countgraph("counts.cg", ["ACGTTGC", "GCAACGT"])
print(result["counts"])
```

- **示例输出：**

```text
{'ACGTTGC': 1, 'GCAACGT': 1}
```

## 15) `dnakit.microbiome.estimate_kmer_table_size()`

- **作用：** 给定 k-mer 数目，按内存上限或假阳性率确定 khmer 表的布局。
- **处理方法：** 调用 `khmer.khmer_args.optimal_size`。`mem_cap` 与 `mem_use` 使用 khmer 的表单位（所有表的哈希位总数），khmer 不会换算为字节。
- **API：** `dnakit.microbiome.estimate_kmer_table_size(num_kmers, *, mem_cap=None, fp_rate=None)`。
- **输入：** `num_kmers` 为正整数；`mem_cap` 与 `fp_rate` 二选一。
- **输出：** `num_htables`、`htable_size`、`mem_use` 和 `expected_fp_rate`。
- **示例代码：**

```python
from dnakit.microbiome import estimate_kmer_table_size

result = estimate_kmer_table_size(1_000_000, fp_rate=0.01)
print(result["num_htables"], result["htable_size"], round(result["expected_fp_rate"], 4))
```

- **示例输出：**

```text
6 1602775 0.01
```

## 16) `dnakit.microbiome.expected_kmer_table_collisions()`

- **作用：** 用文件中的 reads 填充 khmer 存在表（Nodegraph），估计其假阳性率。
- **处理方法：** 调用 `khmer.calc_expected_collisions`。它用占用位数除以最小表大小，再取表数次幂，是快速近似。上游函数超过阈值时会打印错误并退出；本函数不会，而是返回 `within_threshold`。
- **API：** `dnakit.microbiome.expected_kmer_table_collisions(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000, max_false_positive=0.2)`。
- **输入：** FASTA 或 FASTQ 文件；`max_false_positive` 在 0（不含）与 1（含）之间。
- **输出：** `hashsizes`、`occupied_bins`、`expected_false_positive_rate`、`within_threshold` 以及 reads 计数。
- **示例代码：**

```python
from dnakit.microbiome import expected_kmer_table_collisions

result = expected_kmer_table_collisions("reads.fa", ksize=9, table_size=10007, n_tables=2)
print(result["within_threshold"], result["occupied_bins"], result["hashsizes"])
```

- **示例输出：**

```text
True 40 [9973, 9967]
```
