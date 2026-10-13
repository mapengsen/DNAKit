# 可选 MinHash 与 k-mer 后端

这一组函数是 `dnakit.microbiome` 中基于可选依赖的 DNA 序列适配器。[sourmash](https://github.com/sourmash-bio/sourmash) 用于计算 scaled MinHash 草图、Jaccard 相似度和包含度，并读写 sourmash 签名文件；[khmer](https://github.com/dib-lab/khmer) 用于规范 k-mer 计数、唯一 k-mer 估计、丰度直方图和数字归一化（digital normalization）。这些函数只读取显式给出的 FASTA 或 FASTQ 文件，不下载数据库，也不做物种分类或组装；读入的 reads 数和碱基数都有显式上限。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。sourmash 对应 `bio-sourmash`，khmer 对应 `bio-khmer`。

## 示例数据

下面的示例都使用这个 FASTA。两条基因组只在一个碱基上不同，第三条是短 read。它只用于演示接口，统计值不代表真实微生物群落。

**示例文件：** `sample.fa`

```text
>genome_a
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGATCCATGGCA
>genome_b
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGTTCCATGGCA
>read_c
TTTTGGGGCCCCAAAAGGGGTTTTACGTACGTAAAA
```

## 1) `dnakit.microbiome.sketch_fasta_minhash()`

- **作用：** 把 FASTA 中每条记录切成 scaled MinHash 草图。
- **处理方法：** 调用 `sourmash.MinHash(n=0, ksize, scaled)` 的 `add_sequence(..., force=True)`。k-mer 取规范形式，正向与反向互补共用一个哈希；含 A、C、G、T 以外符号的 k-mer 由 sourmash 跳过。FASTA 解析是纯 Python 实现，并有显式上限。
- **API：** `dnakit.microbiome.sketch_fasta_minhash(source, *, ksize=31, scaled=1000, max_records=1000, max_bases=10000000)`。
- **输入：** FASTA 文件；`ksize` 取 1 到 100；`scaled` 为 1 到 10^9 的整数，越小保留的哈希越多；超过 `max_records` 或 `max_bases` 时报错。
- **输出：** 含 `records` 和 `sketches` 的字典。每条草图包含 `id`、`length`、`ksize`、`scaled`、`hash_count` 和排序后的 `hashes`（无符号 64 位整数）。
- **示例代码：**

```python
from dnakit.microbiome import sketch_fasta_minhash

result = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)
print(result["records"], [(s["id"], s["hash_count"]) for s in result["sketches"]])
```

- **示例输出：**

```text
3 [('genome_a', 5), ('genome_b', 4), ('read_c', 0)]
```

## 2) `dnakit.microbiome.compare_minhash_sketches()`

- **作用：** 比较两条草图的 Jaccard 相似度和包含度。
- **处理方法：** 用 sourmash 的 `MinHash.similarity(ignore_abundance=True)` 计算 Jaccard，用 `contained_by()` 计算两个方向的包含度。两条草图的 `ksize` 与 `scaled` 必须相同。包含度不对称：`left_contained_by_right` 等于 |L∩R| / |L|。
- **API：** `dnakit.microbiome.compare_minhash_sketches(left, right)`。
- **输入：** 两个草图字典，可以来自 `sketch_fasta_minhash()` 或 `read_sourmash_signatures()` 的 `sketches` 元素。哈希不能为空，且必须落在 `scaled` 的范围内。
- **输出：** Jaccard 相似度、两个方向的包含度、共享哈希数以及两边的哈希数。
- **示例代码：**

```python
from dnakit.microbiome import compare_minhash_sketches, sketch_fasta_minhash

sketches = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"]
result = compare_minhash_sketches(sketches[0], sketches[1])
print(round(result["jaccard"], 3), round(result["left_contained_by_right"], 3))
```

- **示例输出：**

```text
0.5 0.603
```

## 3) `dnakit.microbiome.write_sourmash_signatures()`

- **作用：** 把草图写成一个 sourmash JSON 签名文件。
- **处理方法：** 用 `sourmash.signature.SourmashSignature` 包装每条草图，再用 `sourmash.save_signatures_to_json()` 写入。文件名以 `.gz` 结尾时自动压缩。已存在的文件不会被覆盖。
- **API：** `dnakit.microbiome.write_sourmash_signatures(sketches, output)`。
- **输入：** 草图列表，每条草图至少含一个哈希；`output` 必须以 `.sig` 或 `.sig.gz` 结尾，所在目录必须已存在。
- **输出：** 写入的文件路径和签名数量。
- **示例代码：**

```python
from pathlib import Path
from dnakit.microbiome import sketch_fasta_minhash, write_sourmash_signatures

Path("sketches.sig.gz").unlink(missing_ok=True)
sketches = [s for s in sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"] if s["hash_count"]]
result = write_sourmash_signatures(sketches, "sketches.sig.gz")
print(result["signatures"], Path(result["output"]).name)
```

- **示例输出：**

```text
2 sketches.sig.gz
```

## 4) `dnakit.microbiome.read_sourmash_signatures()`

- **作用：** 读取 sourmash 签名文件，返回与 `sketch_fasta_minhash()` 相同结构的草图。
- **处理方法：** 用 `sourmash.load_file_as_signatures()` 加载。只接受 DNA 的 scaled 签名，蛋白质签名和 num 型草图会报错。丢弃 abundance 值，并记录每个签名的 `md5sum`。
- **API：** `dnakit.microbiome.read_sourmash_signatures(source, *, max_signatures=1000)`。
- **输入：** `.sig` 或 `.sig.gz` 文件；签名数量超过 `max_signatures` 时报错。
- **输出：** `sketches` 列表。`length` 为 `None`，因为签名不保存序列长度。
- **示例代码：**

```python
from dnakit.microbiome import read_sourmash_signatures

result = read_sourmash_signatures("sketches.sig.gz")
print([(s["id"], s["hash_count"]) for s in result["sketches"]])
```

示例读取上一节写出的 `sketches.sig.gz`，因此两段示例要按顺序运行。

- **示例输出：**

```text
[('genome_a', 5), ('genome_b', 4)]
```

## 5) `dnakit.microbiome.query_kmer_counts()`

- **作用：** 统计指定 k-mer 在 reads 中出现的次数。
- **处理方法：** 用 khmer `Countgraph` 计数，k-mer 取规范形式，正向与反向互补共用一个计数。读入时统一转为大写；含 A、C、G、T 以外符号的 reads 和短于 k 的 reads 会被跳过并计数。计数表是近似的，哈希冲突只会使计数偏大。
- **API：** `dnakit.microbiome.query_kmer_counts(source, kmers, *, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** FASTA 或 FASTQ 文件；`kmers` 为等长的大写 ACGT 字符串列表，长度 1 到 32。
- **输出：** 每个查询 k-mer 的计数，以及读入、使用、跳过的 reads 数。
- **示例代码：**

```python
from dnakit.microbiome import query_kmer_counts

result = query_kmer_counts("sample.fa", ["ACGTTG", "TTTTGG"])
print(result["counts"], result["used_reads"])
```

- **示例输出：**

```text
{'ACGTTG': 2, 'TTTTGG': 2} 3
```

## 6) `dnakit.microbiome.estimate_unique_kmers()`

- **作用：** 估计文件中不同的规范 k-mer 数量。
- **处理方法：** 用 khmer `Nodegraph`（存在性表）记录 k-mer，唯一数量由 `n_unique_kmers()` 根据已占用的槽位估计。结果是近似值，表越大越接近精确数量。
- **API：** `dnakit.microbiome.estimate_unique_kmers(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** FASTA 或 FASTQ 文件；`ksize` 取 1 到 32。
- **输出：** k-mer 位置总数、唯一数量估计、已占用槽位数和哈希表大小。
- **示例代码：**

```python
from dnakit.microbiome import estimate_unique_kmers

result = estimate_unique_kmers("sample.fa", ksize=21, table_size=1_000_000)
print(result["kmer_positions"], result["unique_kmers_estimate"])
```

- **示例输出：**

```text
56 46
```

## 7) `dnakit.microbiome.kmer_abundance_histogram()`

- **作用：** 统计每个出现次数对应多少个不同的规范 k-mer。
- **处理方法：** 先把所有 reads 计入 `Countgraph`，再调用 khmer 的 `abundance_distribution()` 扫描同一个文件。khmer 会重新读取文件，因此每条 read 都必须是大写 A、C、G、T，且长度不小于 `ksize`；否则报错，而不是跳过。
- **API：** `dnakit.microbiome.kmer_abundance_histogram(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** FASTA 或 FASTQ 文件；`ksize` 取 1 到 32。
- **输出：** `distinct_kmers` 和 `histogram` 列表，每项包含 `count` 与 `kmers`。
- **示例代码：**

```python
from dnakit.microbiome import kmer_abundance_histogram

result = kmer_abundance_histogram("sample.fa", ksize=11)
print(result["distinct_kmers"], result["histogram"])
```

- **示例输出：**

```text
60 [{'count': 1, 'kmers': 34}, {'count': 2, 'kmers': 26}]
```

## 8) `dnakit.microbiome.normalize_reads_by_median()`

- **作用：** 数字归一化：丢弃覆盖度过高的重复 reads。
- **处理方法：** 按文件顺序处理每条 read。先用 `get_median_count()` 取规范 k-mer 计数的中位数，中位数小于 `cutoff` 时保留并计入计数表，否则丢弃。输出为 FASTA，只保留标题和序列。
- **API：** `dnakit.microbiome.normalize_reads_by_median(source, output, *, ksize=20, cutoff=20, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`。
- **输入：** FASTA 或 FASTQ 文件，以及输出路径（以 `.fa`、`.fasta` 或 `.fna` 结尾，不覆盖已存在的文件）。
- **输出：** 读入、保留、因中位数过高而丢弃、非 ACGT 跳过和过短跳过的 reads 数。
- **示例代码：**

```python
from pathlib import Path
from dnakit.microbiome import normalize_reads_by_median

Path("normalized.fa").unlink(missing_ok=True)
result = normalize_reads_by_median("sample.fa", "normalized.fa", ksize=11, cutoff=1)
print(result["kept_reads"], result["discarded_high_median_reads"])
```

- **示例输出：**

```text
2 1
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-sourmash]"    # sourmash：草图、比较、签名文件
python -m pip install "dnakit[bio-khmer]"       # khmer：k-mer 计数、估计、直方图、归一化
```

项目地址：[sourmash](https://github.com/sourmash-bio/sourmash)、[khmer](https://github.com/dib-lab/khmer)。

## 验证方式

测试文件是 `others/tests/unit/test_microbiome_backends.py`，共 19 项，主要检查：

- sourmash 草图与直接调用 `sourmash.MinHash` 的结果一致；规范 k-mer 哈希数用手算校验（`ACGTACGT` 的 3-mer 得到 2 个哈希）。
- Jaccard 与包含度用手算集合校验，并与 sourmash 直接计算的结果一致；签名文件写入后可被 sourmash 原样读回。
- khmer 的 k-mer 计数与 Python 精确统计的规范 k-mer 计数一致；丰度直方图与精确统计一致；数字归一化的手算例子保留 2 条、丢弃 1 条。
- 唯一 k-mer 估计与 khmer 直接调用的结果一致，并与精确数量相差不超过 1%。
- 缺少 sourmash 或 khmer 时给出安装提示；ksize、max_reads、非 ACGT 输入、不合法的草图等都会被拒绝。
- 蛋白质签名的拒绝规则使用替身加载器测试，因为 sourmash 在本环境中无法重新加载极短的蛋白质签名。

## 当前范围

这一批覆盖 sourmash 的 scaled MinHash 草图、比较、签名文件读写，以及 khmer 的 k-mer 计数、唯一 k-mer 估计、丰度直方图和数字归一化。

- 未包装 sourmash 的 `search`、`gather`、`taxonomy` 和 LCA/SBT 索引：它们需要预先构建的数据库，超出文件级接口。
- 未包装 num 型草图、蛋白质草图和丰度加权比较（`track_abundance`）：DNA 工作流只使用 scaled 草图，本函数只返回集合关系。
- 未包装 khmer 的组装器（`LinearAssembler`、`JunctionCountAssembler`）、stop tag 修剪和图文件保存加载：它们生成 contig 或修改图文件，不属于计数接口。
- 未包装 FASTQ 质量裁剪和双端 reads 处理。
- `dnakit.integrations.compare_sourmash_sequences()` 已比较两条字符串；本批的 `compare_minhash_sketches()` 比较文件中的草图，二者用途不同。
