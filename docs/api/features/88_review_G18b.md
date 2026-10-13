# G18b 功能：sourmash、khmer 与 primer3

本页记录 G18b 组的包装函数。这些函数从 `dnakit.microbiome` 与 `dnakit.molecular_lab` 导出。例如 `from dnakit.microbiome import kmer_mutation_probability`。

- 本地函数：全部 9 个包装函数都可离线运行。输入为数值、本地 DNA 字符串、谱系字符串或本地 FASTA、FASTQ 文件。没有函数下载数据或访问网络，因此本组没有在线测试。
- khmer 的输入在交给 khmer 之前先转为大写。本组的 khmer 函数都不接受 k-mer 长度，因此 ksize 的 1 到 32 检查在此不适用。
- 示例中的输入文件以 `示例文件` 标签给出；其余示例直接运行。

## 已有包装与未包装的函数

- 本次审查前后 `other` 与 `cap` 行的台账计数（没有 `cap` 行）。
  - sourmash：121 行新增 6 个包装（`wrapped` 共 14 个）。其余 115 行为：内部辅助与常量 29 个；命令行或报告辅助 25 个；已包装函数的变体 22 个（签名加载与保存、JSON 字符串解析、搜索驱动、nodegraph 冲突与头部辅助）；没有序列含义的分类或 picklist 辅助 20 个；与新包装重复的 2 个；需要预先构建文件的 17 个（10 个 LCA 数据库、3 个 SBT 索引、2 个 SQLite 索引、2 组 gather CSV 文件）。
  - khmer：4 行新增 2 个包装（`read_paired_fastx` 与 `clean_khmer_read_sequences`）。其余 2 个为 `get_n_primes_near_x` 与 `is_prime`，只用于确定哈希表大小。
  - primer3：4 行新增 1 个包装（`validate_acgt_sequence`）。`ensure_acgt_uppercase_b` 是其 bytes 形式（变体）。`includes` 与 `get_libprimer3_version` 是打包和溯源辅助。
- 若干同类函数是现有包装的变体。`search_minhash_database` 覆盖包含度搜索；搜索驱动（`search_databases_with_flat_query`、`search_databases_with_abund_query`）与 SBT 树搜索（`search_sbt_index_private`、`sbtmh.search_sbt_index`）执行同类搜索，只是索引或评分不同。`read_sourmash_signatures` 读取与签名加载函数相同的 `.sig` 文件，`expected_kmer_table_collisions` 覆盖 nodegraph 的冲突估计。
- 因缺少文件而未包装：LCA 数据库函数（`summarize`、`gather_assignments`、`classify_signature`、`make_lca_counts` 及其加载函数）、SBT 与 SQLite 索引加载函数，以及 gather CSV 加载函数。这些函数都需要预先构建的文件，包装函数不会生成这些文件。
- 每个清单键的完整状态见 G18b 台账。

## 1) `kmer_mutation_probability`

- **作用：** 给定每个碱基的突变率，计算一个 k-mer 至少含一个突变碱基的概率。
- **处理方法：** 调用 sourmash 的 `distance_utils.r1_to_q`，计算 `1 - (1 - rate) ** ksize`。假设各碱基独立突变。
- **API：** `kmer_mutation_probability(mutation_rate[必须], ksize[必须], verbose[可选])`
- **输入：** `mutation_rate` 为 [0, 1] 内的数；`ksize` 为 1 到 100 的整数。
- **输出：** 字典，含 `ksize`、`mutation_rate`、`kmer_mutation_probability` 与 `provenance`。

**示例代码：**

```python
from dnakit.microbiome import kmer_mutation_probability

for k in (1, 2, 21):
    result = kmer_mutation_probability(0.01, k, verbose=False)
    print(k, round(result["kmer_mutation_probability"], 6))
```

**示例输出：**

```text
1 0.01
2 0.0199
21 0.190272
```

## 2) `expected_mutated_kmers`

- **作用：** 给定每碱基突变率，求一条序列中会被改变的不同 k-mer 数目的期望值。
- **处理方法：** 对不同 k-mer 的数目调用 sourmash 的 `distance_utils.exp_n_mutated`，结果为该数目乘以每个 k-mer 的突变概率。长度为 L 且无重复的序列，不同 k-mer 数为 `L - ksize + 1`；重复的 k-mer 不合并。
- **API：** `expected_mutated_kmers(ksize[必须], mutation_rate[必须], sequence_len_bp[可选], n_unique_kmers[可选], verbose[可选])`
- **输入：** `sequence_len_bp`（不小于 `ksize`）与 `n_unique_kmers`（正整数）二者必须恰给一个。`mutation_rate` 在 [0, 1] 内。
- **输出：** 字典，含 `ksize`、`mutation_rate`、`n_unique_kmers`、`kmer_mutation_probability`、`expected_mutated_kmers`、`expected_unmutated_kmers` 与 `provenance`。输出值是期望而非计数。

**示例代码：**

```python
from dnakit.microbiome import expected_mutated_kmers

result = expected_mutated_kmers(21, 0.01, sequence_len_bp=5000, verbose=False)
print(result["n_unique_kmers"], round(result["expected_mutated_kmers"], 2))
```

**示例输出：**

```text
4980 947.56
```

## 3) `variance_mutated_kmers`

- **作用：** 在同一突变模型下，求被突变的不同 k-mer 数目的方差。
- **处理方法：** 调用 sourmash 的 `distance_utils.var_n_mutated`，即突变 k-mer 计数的闭式方差。突变率为 0 时结果为 0。对于极小的突变率，闭式在浮点运算中可能变为负值；sourmash 对此会抛出异常，本包装将其报告为配置错误。
- **API：** `variance_mutated_kmers(ksize[必须], mutation_rate[必须], sequence_len_bp[可选], n_unique_kmers[可选], verbose[可选])`
- **输入：** 与 `expected_mutated_kmers` 相同。
- **输出：** 字典，含 `ksize`、`mutation_rate`、`n_unique_kmers`、`variance_mutated_kmers`、`standard_deviation_mutated_kmers` 与 `provenance`。

**示例代码：**

```python
from dnakit.microbiome import variance_mutated_kmers

result = variance_mutated_kmers(2, 0.1, sequence_len_bp=11, verbose=False)
print(round(result["variance_mutated_kmers"], 4))
```

**示例输出：**

```text
2.8512
```

## 4) `probability_no_shared_kmers`

- **作用：** 求一条序列的 scaled 草图与其突变副本的草图不共享任何哈希的概率。
- **处理方法：** 调用 sourmash 的 `distance_utils.get_exp_probability_nothing_common`，计算 `exp((n - E[突变数]) * ln(1 - 1/scaled))`。突变率为 0 时结果为 0，突变率为 1 时结果为 1。这是期望值，较大的结果提示两个草图可能毫无共同之处。
- **API：** `probability_no_shared_kmers(mutation_rate[必须], ksize[必须], scaled[必须], sequence_len_bp[可选], n_unique_kmers[可选], verbose[可选])`
- **输入：** `scaled` 为不小于 2 的整数；scale 为 1 时保留全部 k-mer，对数无定义。`sequence_len_bp` 与 `n_unique_kmers` 二者必须恰给一个。
- **输出：** 字典，含 `mutation_rate`、`ksize`、`scaled`、`n_unique_kmers`、`probability_nothing_in_common` 与 `provenance`。

**示例代码：**

```python
from dnakit.microbiome import probability_no_shared_kmers

for rate in (0.0, 0.01, 0.1):
    result = probability_no_shared_kmers(rate, 21, 1000, sequence_len_bp=5000, verbose=False)
    print(rate, f"{result['probability_nothing_in_common']:.3e}")
```

**示例输出：**

```text
0.0 0.000e+00
0.01 1.770e-02
0.1 5.797e-01
```

## 5) `bp_threshold_to_hashes`

- **作用：** 把以 bp 计的重叠阈值换算为哈希数，以及相对查询草图的包含度。
- **处理方法：** 调用 sourmash 的 `search.calc_threshold_from_bp`。哈希数为 `threshold_bp / scaled`，包含度为该哈希数除以 `query_size`。阈值超过查询本身的 bp 数时永远无法达到，sourmash 会抛出异常，本包装将其报告为配置错误。
- **API：** `bp_threshold_to_hashes(threshold_bp[必须], scaled[必须], query_size[必须], verbose[可选])`
- **输入：** `threshold_bp` 为非负整数（bp）；`scaled` 为正整数；`query_size` 为查询草图中的哈希数。
- **输出：** 字典，含 `threshold_bp`、`scaled`、`query_size`、`threshold_hashes`（可能为小数）、`containment_threshold`（在 [0, 1] 内）与 `provenance`。

**示例代码：**

```python
from dnakit.microbiome import bp_threshold_to_hashes

result = bp_threshold_to_hashes(50_000, 1000, 200, verbose=False)
print(result["threshold_hashes"], result["containment_threshold"])
```

**示例输出：**

```text
50.0 0.25
```

## 6) `count_kmer_lca_assignments`

- **作用：** 对每个 k-mer 哈希，求其分类赋值的最近公共祖先，再统计每个祖先被多少个哈希共享。
- **处理方法：** 调用 sourmash 的 `lca.lca_utils.count_lca_for_assignments`。每个哈希贡献其最近公共祖先，未给权重时计为 1。`lowest_common_lineage` 处理单个谱系列表；本包装对一个样本中的每个哈希执行同样的操作。
- **API：** `count_kmer_lca_assignments(assignments[必须], hashval_counts[可选], verbose[可选])`
- **输入：** `assignments` 把 64 位哈希映射到非空的 `;` 分隔谱系字符串列表，从域（superkingdom）往下，至多 8 个等级，且不含逗号。`hashval_counts` 若给出，则须把完全相同的哈希映射到正整数权重。
- **输出：** 字典，含 `hash_count`、`weighted`、`total_count` 与 `counts`。`counts` 的每项含 `lineage`（空字符串表示根）、`depth`、`lowest_rank` 与 `count`。各项先按计数从大到小排序，再按谱系排序。

**示例代码：**

```python
from dnakit.microbiome import count_kmer_lca_assignments

assignments = {
    101: ["d__Bacteria;p__Firmicutes;c__Bacilli", "d__Bacteria;p__Firmicutes;c__Clostridia"],
    102: ["d__Bacteria;p__Bacteroidetes"],
    103: ["d__Archaea", "d__Bacteria"],
}
result = count_kmer_lca_assignments(assignments, verbose=False)
for item in result["counts"]:
    print(repr(item["lineage"]), item["lowest_rank"], item["count"])
```

**示例输出：**

```text
'' None 1
'd__Bacteria;p__Bacteroidetes' phylum 1
'd__Bacteria;p__Firmicutes' phylum 1
```

## 7) `clean_khmer_read_sequences`

- **作用：** 与 khmer 自身的清理一致地把 reads 转为大写并将 N 替换为 A，同时报告替换了多少个 N。
- **处理方法：** 调用 `khmer.utils.clean_input_reads`，它把 `cleaned_seq` 设为大写且 N 换成 A 的 read。这种替换有损，清理后的 read 不再是原始序列。每条 read 报告其 `replaced_n_count`。dnakit 自身的 khmer 计数与修剪包装遇到含 N 的 reads 时会跳过它们，而不是替换。
- **API：** `clean_khmer_read_sequences(sequences[必须], verbose[可选])`
- **输入：** 从 read ID 到 DNA 字符串的映射，允许大小写的 A、C、G、T 与 N。至多 10000 条 reads，每条至多 100000 个碱基；其他符号会报错。
- **输出：** 字典，含 `read_count`、`reads_with_n`、`replaced_n_total` 与按 ID 排序的 `reads`。每条 read 含 `id`、`sequence`（大写）、`cleaned_sequence` 与 `replaced_n_count`。

**示例代码：**

```python
from dnakit.microbiome import clean_khmer_read_sequences

result = clean_khmer_read_sequences({"read_b": "ACGT", "read_a": "acgNNt"}, verbose=False)
for item in result["reads"]:
    print(item["id"], item["sequence"], item["cleaned_sequence"], item["replaced_n_count"])
print(result["replaced_n_total"])
```

**示例输出：**

```text
read_a ACGNNT ACGAAT 2
read_b ACGT ACGT 0
2
```

## 8) `read_paired_fastx`

- **作用：** 配对本地 FASTA 或 FASTQ 文件中的记录，未配对的记录作为单条保留。
- **处理方法：** 用 `khmer.ReadParser` 读取文件，再把记录交给 `khmer.utils.broken_paired_reader`。两条相邻记录的名称若作为配对 mate 匹配（如 `name/1` 与 `name/2`，或 `1:`/`2:` 与 `seq/1`/`seq/2` 形式），则组成一对。序列转为大写。
- **API：** `read_paired_fastx(source[必须], min_length[可选], require_paired[可选], force_single[可选], max_reads[可选], verbose[可选])`
- **输入：** 本地 FASTA 或 FASTQ 文件。`min_length` 丢弃短于它的 reads；`require_paired=True` 把任何孤立 read 变为错误；`force_single=True` 把每条记录都作为单条返回；`max_reads`（默认 100000）限制记录数。`require_paired` 与 `force_single` 不能同时为真。
- **输出：** 字典，含 `unit_count`、`pair_count`、`singleton_count` 与按文件顺序排列的 `units`。每个单元含 `index`（按 khmer 计记录数的方式计数，一对前进 2）、`is_pair`、`read1` 与 `read2`（单条时为 `None`）。每条 read 含 `name` 与 `sequence`。

**示例文件：** `pairs.fq`

```text
@r1/1
ACGTACGT
+
IIIIIIII
@r1/2
TTGGCCAA
+
IIIIIIII
@r2/1
GGGGAAAA
+
IIIIIIII
@r3/1
CCCCTTTT
+
IIIIIIII
@r3/2
AAAACCCC
+
IIIIIIII
```

**示例代码：**

```python
from dnakit.microbiome import read_paired_fastx

result = read_paired_fastx("pairs.fq", verbose=False)
for unit in result["units"]:
    names = [unit["read1"]["name"]] + ([unit["read2"]["name"]] if unit["read2"] else [])
    print(unit["index"], unit["is_pair"], " ".join(names))
print(result["pair_count"], result["singleton_count"])
```

**示例输出：**

```text
0 True r1/1 r1/2
2 False r2/1
3 True r3/1 r3/2
2 1
```

## 9) `validate_acgt_sequence`

- **作用：** 把 DNA 序列转为大写，并要求其只含 A、C、G、T，与引物设计的输入要求一致。
- **处理方法：** 调用 primer3-py 的 `p3helpers.ensure_acgt_uppercase`。与 dnakit 的规范化字母表检查不同，它不把 U 转为 T，也不去除空白，因此遇到这些符号会拒绝，而不是改写。
- **API：** `validate_acgt_sequence(sequence[必须], verbose[可选])`
- **输入：** 至多 100000 个字符的非空字符串。接受小写 a、c、g、t。
- **输出：** 字典，含 `length`、`sequence`（大写）与 `provenance`，其中记录 primer3-py 版本与 GPL-2.0 许可证。

**示例代码：**

```python
from dnakit.molecular_lab import validate_acgt_sequence

print(validate_acgt_sequence("acgtacgt", verbose=False)["sequence"])
try:
    validate_acgt_sequence("ACGU", verbose=False)
except Exception as exc:
    print(type(exc).__name__)
```

**示例输出：**

```text
ACGTACGT
ConfigurationError
```
