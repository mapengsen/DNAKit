# sourmash：谱系数据库、分类学报告与签名工具

这一组位于 `dnakit.microbiome` 的函数封装了 [sourmash](https://github.com/sourmash-bio/sourmash)，它是用于比较和分类序列数据的 FracMinHash 工具集。这些函数处理你提供的本地文件：LCA 谱系数据库、签名文件、SBT 与 SQLite 索引、`sourmash compare` 矩阵，以及 `sourmash gather` 与谱系 CSV 文件。它们可以载入并描述谱系数据库，在谱系树上对签名分类和汇总，比较谱系表，写出 LCA 索引构建的报告，把 gather 结果转换为 Krona、CSV、文本、Kraken 风格、CAMI 和 LIN 谱系组报告，并处理签名与 SBT 索引（按清单做草图、grep、inflate、合并、迁移、转换）。它们不下载任何内容，也没有函数会根据测序读段构建数据库。

sourmash 使用 BSD-3-Clause 许可（已在已安装包的元数据中核对）。该依赖是可选的，只有调用函数时才会导入；缺少时函数会抛出 `BackendUnavailableError` 并给出安装命令。示例从 `dnakit.microbiome` 导入。

## 示例数据

示例使用微型的合成数据而不是真实基因组，因此每个数字都可以手工核对：三个基因组，k=31、scaled=1000，哈希是很小的整数（哈希 6 至 10 由两个 *Escherichia* 基因组共有）。需要比较矩阵、SBT 索引或 SQLite 索引的示例，会在临时目录中用 NumPy 和 sourmash 的接口创建它们，所有输出文件都写入临时目录。

含三个基因组的 LCA 数据库：GCF_A（大肠杆菌，哈希 1 至 10）、GCF_B（Escherichia fergusonii，哈希 6 至 15）和 GCF_C（枯草芽孢杆菌，哈希 21 至 30）。

**示例文件：** `db.lca.json`

```text
{"version": "2.1", "type": "sourmash_lca", "license": "CC0", "ksize": 31, "scaled": 1000, "moltype": "DNA", "lid_to_lineage": {"0": [["superkingdom", "Bacteria"], ["phylum", "Pseudomonadota"], ["class", "Gammaproteobacteria"], ["order", "Enterobacterales"], ["family", "Enterobacteriaceae"], ["genus", "Escherichia"], ["species", "Escherichia coli"]], "1": [["superkingdom", "Bacteria"], ["phylum", "Pseudomonadota"], ["class", "Gammaproteobacteria"], ["order", "Enterobacterales"], ["family", "Enterobacteriaceae"], ["genus", "Escherichia"], ["species", "Escherichia fergusonii"]], "2": [["superkingdom", "Bacteria"], ["phylum", "Bacillota"], ["class", "Bacilli"], ["order", "Bacillales"], ["family", "Bacillaceae"], ["genus", "Bacillus"], ["species", "Bacillus subtilis"]]}, "hashval_to_idx": {"1": [0], "2": [0], "3": [0], "4": [0], "5": [0], "6": [0, 1], "7": [0, 1], "8": [0, 1], "9": [0, 1], "10": [0, 1], "11": [1], "12": [1], "13": [1], "14": [1], "15": [1], "21": [2], "22": [2], "23": [2], "24": [2], "25": [2], "26": [2], "27": [2], "28": [2], "29": [2], "30": [2]}, "ident_to_name": {"GCF_A": "GCF_A Escherichia coli", "GCF_B": "GCF_B Escherichia fergusonii", "GCF_C": "GCF_C Bacillus subtilis"}, "ident_to_idx": {"GCF_A": 0, "GCF_B": 1, "GCF_C": 2}, "idx_to_lid": {"0": 0, "1": 1, "2": 2}}
```

查询签名 `query1`，哈希为 1 至 8、22 和 23。

**示例文件：** `query.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"query1","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,4,5,6,7,8,22,23],"md5sum":"80370ff3293372ea9e60f59626212d6a","molecule":"DNA"}],"version":0.4}]
```

查询签名 `query2`，哈希为 11 至 15。

**示例文件：** `query2.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"query2","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[11,12,13,14,15],"md5sum":"556d99975b863d98753e9796ed9cd6a5","molecule":"DNA"}],"version":0.4}]
```

数据库中三个基因组的签名，合并为一个签名文件。

**示例文件：** `genomes.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_A Escherichia coli","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,4,5,6,7,8,9,10],"md5sum":"41276f1cde0088f537f94eb1456e7dc7","molecule":"DNA"}],"version":0.4},{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_B Escherichia fergusonii","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[6,7,8,9,10,11,12,13,14,15],"md5sum":"28cd92c421ea26a1c81c441fcbd134b4","molecule":"DNA"}],"version":0.4},{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_C Bacillus subtilis","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[21,22,23,24,25,26,27,28,29,30],"md5sum":"44c0573131a0d93f2277857df1d1b155","molecule":"DNA"}],"version":0.4}]
```

记录丰度的签名 `reads_abund`。

**示例文件：** `abund.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"reads_abund","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,21,30],"md5sum":"e3fb529d70ad3c826c56d757251058fd","abundances":[5,3,1,2,7],"molecule":"DNA"}],"version":0.4}]
```

`query1` 的 gather CSV，含两个匹配。

**示例文件：** `gather.csv`

```text
query_name,name,f_unique_weighted,f_unique_to_query,unique_intersect_bp,remaining_bp,query_md5,query_filename,query_bp,ksize,scaled,query_n_hashes,sum_weighted_found,total_weighted_hashes
query1,GCF_A Escherichia coli,0.8,0.8,8000,2000,80370ff3,query.sig,10000,31,1000,10,8,10
query1,GCF_C Bacillus subtilis,0.2,0.2,2000,0,80370ff3,query.sig,10000,31,1000,10,10,10
```

三个基因组的谱系表。

**示例文件：** `tax.csv`

```text
ident,superkingdom,phylum,class,order,family,genus,species
GCF_A,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
GCF_B,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii
GCF_C,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
```

同样的谱系，另加 NCBI `taxpath` 列。

**示例文件：** `taxid.csv`

```text
ident,superkingdom,phylum,class,order,family,genus,species,taxpath
GCF_A,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli,2|1224|1236|91347|543|561|562
GCF_B,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii,2|1224|1236|91347|543|561|564
GCF_C,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis,2|1239|91061|1385|186817|1386|1423
```

三个基因组的 LIN 谱系。

**示例文件：** `lin.csv`

```text
ident,lin
GCF_A,0;0;0
GCF_B,0;0;1
GCF_C,1;0;0
```

两个 LIN 谱系组。

**示例文件：** `lingroup.csv`

```text
name,lin
group_ecoli_like,0
group_bacillus_like,1
```

用于比较的 classify 风格表格。

**示例文件：** `classify.csv`

```text
ID,status,superkingdom,phylum,class,order,family,genus,species
q1,found,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q2,found,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,
q3,found,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
q4,found,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
```

用于比较的自定义谱系表。

**示例文件：** `custom.csv`

```text
ID,superkingdom,phylum,class,order,family,genus,species
q1,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q2,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii
q3,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q5,Archaea,Euryarchaeota,Methanobacteria,Methanobacteriales,Methanobacteriaceae,Methanobacterium,Methanobacterium formicicum
```

供 `sketch_fromfile_signatures()` 使用的清单。

**示例文件：** `genomes.csv`

```text
name,genome_filename,protein_filename
genome1,genome1.fa,
```

一个 93 个碱基的基因组。

**示例文件：** `genome1.fa`

```text
>genome1
ATGGCGTACGTTAGCCTAGGATCCGATTACAGGCTTAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC
```

## 1) `dnakit.microbiome.check_and_load_gather_csvs_summary()`

- **作用：** 把多个 `sourmash gather` CSV 与谱系表一起载入，并为每个查询返回一行摘要；遇到空文件或缺失标识符时尽早报错。
- **处理方法：** 先用 sourmash 的 `MultiLineageDB.load()` 读取谱系表，再把 gather 文件交给 `sourmash.tax.tax_utils.check_and_load_gather_csvs()`。该上游函数会拒绝空的 gather 文件；设置 `fail_on_missing_taxonomy=True` 时，还会拒绝没有谱系的匹配标识符。DNAKit 随后把每个已载入的查询压缩成一行摘要。`force=True` 会同时传给谱系表加载器和 gather 加载器：空的或无法读取的 gather 文件会被跳过，而不是抛出 `ValueError`；但同一个查询出现在两个 gather 文件中时仍会抛出该错误。
- **API：** `check_and_load_gather_csvs_summary(gather_csvs: 'Sequence[str]', taxonomy_files: 'Sequence[str]', *, fail_on_missing_taxonomy: 'bool' = False, force: 'bool' = False, keep_full_identifiers: 'bool' = False, keep_identifier_versions: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `gather_csvs` 与 `taxonomy_files` 是非空的已存在文件路径列表（直接传字符串会被拒绝）。谱系表是 sourmash 谱系格式的 CSV 或 sqlite 文件（一个 `ident` 列加每个分类等级一列）。gather CSV 需要 sourmash 4.4 及以上版本的列：`query_name`、`name`、`f_unique_weighted`、`f_unique_to_query`、`unique_intersect_bp`、`remaining_bp`、`query_md5`、`query_filename`、`query_bp`、`ksize` 和 `scaled`。`keep_full_identifiers` 与 `keep_identifier_versions` 控制如何把匹配名称转换为谱系标识符。
- **输出：** 返回字典包含 `status`、`backend`、解析后的 `gather_csvs` 路径，以及 `queries`：每个查询一行，含 `query_name`、`query_md5`、`query_filename`、`ksize`、`scaled`、`match_count`、`missed_identifiers` 和 `skipped_identifiers`。
- **示例代码：**

```python
from dnakit.microbiome import check_and_load_gather_csvs_summary

result = check_and_load_gather_csvs_summary(["gather.csv"], ["tax.csv"], verbose=False)
print(result["status"], len(result["gather_csvs"]))
for query in result["queries"]:
    print(query["query_name"], query["match_count"], query["missed_identifiers"], query["scaled"])
```

- **示例输出：**

```text
completed 1
query1 2 0 1000
```

## 2) `dnakit.microbiome.classify_signature_lineage()`

- **作用：** 把一个文件中的每个签名归入其 k-mer 的最低公共祖先谱系，功能等同于 `sourmash lca classify`。
- **处理方法：** 用 `sourmash.lca.lca_utils.load_databases()` 载入 LCA 数据库，按数据库的 k-mer 大小读取文件中的签名；若签名比数据库更精细，则降采样到数据库的 `scaled` 值，然后对每个签名调用 `sourmash.lca.command_classify.classify_signature()`。k-mer 数少于 `threshold` 的谱系会被忽略。`majority=True` 时使用出现最多的谱系，而不是最低公共祖先。若查询比数据库更粗，会抛出 `ConfigurationError`。
- **API：** `classify_signature_lineage(signature: 'str', lca_databases: 'Sequence[str]', *, threshold: 'int' = 5, majority: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `signature` 是已存在的签名文件（可含多个签名）。`lca_databases` 是非空的已存在 LCA 数据库文件列表，k-mer 大小必须一致。`threshold` 是非负整数（sourmash 默认值为 5）。`majority` 是布尔值。
- **输出：** 返回字典包含 `status`、`backend`、`ksize`、`scaled`、`threshold`、`majority` 和 `classifications`：每个签名一行，含 `name`、`status`（`found`、`disagree` 或 `nomatch`）和 `lineage`；`lineage` 是 `{"rank", "name"}` 映射的列表，从最一般的等级排到最具体的等级。
- **示例代码：**

```python
from dnakit.microbiome import classify_signature_lineage

for threshold in (5, 2):
    result = classify_signature_lineage("query.sig", ["db.lca.json"], threshold=threshold, verbose=False)
    row = result["classifications"][0]
    print(threshold, row["name"], row["status"], row["lineage"][-1]["rank"], row["lineage"][-1]["name"])
```

- **示例输出：**

```text
5 query1 found species Escherichia coli
2 query1 disagree superkingdom Bacteria
```

## 3) `dnakit.microbiome.gather_lineage_assignments()`

- **作用：** 针对一个签名的每个 k-mer 哈希，列出 LCA 数据库为它分配的谱系。
- **处理方法：** 用 `load_databases()` 载入数据库；要求文件中恰有一个符合数据库 k-mer 大小的签名，将其降采样到数据库的 `scaled` 值，再对其哈希调用 `sourmash.lca.lca_utils.gather_assignments()`。哈希值转换为字符串，每个哈希的谱系经过排序，因此结果稳定。
- **API：** `gather_lineage_assignments(signature: 'str', lca_databases: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `signature` 是已存在的签名文件，其中恰有一个符合数据库 k-mer 大小的签名（否则抛出 `ConfigurationError`）。`lca_databases` 是非空的已存在 LCA 数据库文件列表。
- **输出：** 返回字典包含 `status`、`backend`、`signature`、`ksize`、`scaled`、`assigned_hash_count` 和 `assignments`：每个在数据库中找到的哈希一行，含 `hashval`（字符串）和 `lineages`（谱系列表）。没有任何数据库包含的哈希不会列出。
- **示例代码：**

```python
from dnakit.microbiome import gather_lineage_assignments

result = gather_lineage_assignments("query.sig", ["db.lca.json"], verbose=False)
print(result["ksize"], result["scaled"], result["assigned_hash_count"])
for row in result["assignments"]:
    species = [lineage[-1]["name"] for lineage in row["lineages"]]
    print(row["hashval"], species)
```

- **示例输出：**

```text
31 1000 10
1 ['Escherichia coli']
2 ['Escherichia coli']
3 ['Escherichia coli']
4 ['Escherichia coli']
5 ['Escherichia coli']
6 ['Escherichia coli', 'Escherichia fergusonii']
7 ['Escherichia coli', 'Escherichia fergusonii']
8 ['Escherichia coli', 'Escherichia fergusonii']
22 ['Bacillus subtilis']
23 ['Bacillus subtilis']
```

## 4) `dnakit.microbiome.lca_rank_counts()`

- **作用：** 统计每个谱系作为最低公共祖先所对应的数据库 k-mer 哈希数。
- **处理方法：** 载入数据库并调用 `sourmash.lca.command_rankinfo.make_lca_counts()`，它遍历数据库中所有哈希并记录每个哈希的最低公共祖先。DNAKit 将谱系序列化，并按计数从大到小、再按谱系排序。
- **API：** `lca_rank_counts(lca_databases: 'Sequence[str]', *, min_num: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `lca_databases` 是非空的已存在 LCA 数据库文件列表。`min_num` 是非负整数：哈希数更少的谱系会被丢弃（默认 0）。数据库按其自身的 `scaled` 值使用。
- **输出：** 返回字典包含 `status`、`backend`、`ksize`、`scaled`、`min_num` 和 `counts`：每行含 `lineage`（`{"rank", "name"}` 列表）和 `count`。
- **示例代码：**

```python
from dnakit.microbiome import lca_rank_counts

result = lca_rank_counts(["db.lca.json"], verbose=False)
for row in result["counts"]:
    print(row["lineage"][-1]["rank"], row["lineage"][-1]["name"], row["count"])
```

- **示例输出：**

```text
species Bacillus subtilis 10
genus Escherichia 5
species Escherichia coli 5
species Escherichia fergusonii 5
```

## 5) `dnakit.microbiome.load_comparison_matrix()`

- **作用：** 读取 `sourmash compare -o BASE` 写出的相似度矩阵与标签。
- **处理方法：** 检查 `basefile` 与 `<basefile>.labels.txt` 是否存在，然后调用 `sourmash.fig.load_matrix_and_labels()`。矩阵被转换为浮点数，并以嵌套列表返回。
- **API：** `load_comparison_matrix(basefile: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）和 numpy（extra `dnakit[bio-analysis]`）。
- **输入：** `basefile` 是传给 `sourmash compare -o` 的路径：以该名称原样保存的 NumPy `.npy` 数组（不会自动添加扩展名），标签文件位于旁边，名为 `<basefile>.labels.txt`。
- **输出：** 返回字典包含 `status`、`backend`、`basefile`、`size`（矩阵边长）、`labels` 和 `matrix`（按行排列的列表）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy

from dnakit.microbiome import load_comparison_matrix

work = Path(tempfile.mkdtemp())
base = work / "cmp"
with open(base, "wb") as handle:
    numpy.save(handle, numpy.array([[1.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 1.0]]))
(work / "cmp.labels.txt").write_text("genome_A\ngenome_B\ngenome_C\n")
result = load_comparison_matrix(str(base), verbose=False)
print(result["size"], result["labels"])
print(result["matrix"][1])
```

- **示例输出：**

```text
3 ['genome_A', 'genome_B', 'genome_C']
[0.5, 1.0, 0.25]
```

## 6) `dnakit.microbiome.load_gather_results_summary()`

- **作用：** 读取一个 `sourmash gather` CSV，从谱系表挂接谱系，并为每个查询生成摘要。
- **处理方法：** 用 `MultiLineageDB.load()` 载入谱系表，再调用 `sourmash.tax.tax_utils.load_gather_results()`，把每个查询压缩成一行摘要。设置 `fail_on_missing_taxonomy=True` 时，没有谱系的匹配会报错；否则计入 `missed_identifiers`。
- **API：** `load_gather_results_summary(gather_csv: 'str', taxonomy_files: 'Sequence[str]', *, fail_on_missing_taxonomy: 'bool' = False, keep_full_identifiers: 'bool' = False, keep_identifier_versions: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `gather_csv` 是已存在的 gather CSV，需包含 sourmash 4.4 及以上版本的列（见函数 1）。`taxonomy_files` 是非空的 CSV 或 sqlite 谱系表列表。`keep_full_identifiers` 与 `keep_identifier_versions` 控制如何把匹配名称简化为标识符。
- **输出：** 返回字典包含 `status`、`backend`、`gather_csv`、`columns`（gather CSV 的表头）和 `queries`，每行字段与函数 1 相同。
- **示例代码：**

```python
from dnakit.microbiome import load_gather_results_summary

result = load_gather_results_summary("gather.csv", ["tax.csv"], verbose=False)
print(len(result["columns"]), result["columns"][:3])
for query in result["queries"]:
    print(query["query_name"], query["query_md5"], query["match_count"], query["ksize"])
```

- **示例输出：**

```text
14 ['query_name', 'name', 'f_unique_weighted']
query1 80370ff3 2 31
```

## 7) `dnakit.microbiome.load_lca_database()`

- **作用：** 载入一个 LCA 数据库，并报告其 k-mer 大小、scaled 值和规模。
- **处理方法：** 对该文件调用 `sourmash.lca.lca_db.load_single_database()`，并统计其中的签名数和哈希数。
- **API：** `load_lca_database(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `path` 是已存在的 LCA 数据库文件（JSON，或 sqlite 形式的 LCA 数据库），由 `sourmash lca index` 生成。
- **输出：** 返回字典包含 `status`、`backend`、`path`、`ksize`、`scaled`、`signature_count` 和 `hash_count`。
- **示例代码：**

```python
from dnakit.microbiome import load_lca_database

result = load_lca_database("db.lca.json", verbose=False)
print(result["ksize"], result["scaled"], result["signature_count"], result["hash_count"])
```

- **示例输出：**

```text
31 1000 3 25
```

## 8) `dnakit.microbiome.load_lca_databases()`

- **作用：** 一次载入多个 LCA 数据库，并描述每个数据库。
- **处理方法：** 调用 `sourmash.lca.lca_utils.load_databases()`，它要求所有数据库的 k-mer 大小一致；返回共同的 k-mer 大小与 `scaled` 值，以及每个数据库的摘要。指定 `scaled` 时，数据库按这个更粗的分辨率载入。
- **API：** `load_lca_databases(lca_databases: 'Sequence[str]', *, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `lca_databases` 是非空的已存在 LCA 数据库文件列表。`scaled` 是可选的正整数。比数据库自身更大的值会对数据库降采样；更小或相等的值则保持数据库自身的 `scaled` 不变。
- **输出：** 返回字典包含 `status`、`backend`、`ksize`、`scaled` 和 `databases`：每个文件一行，含 `path`、`ksize`、`scaled`、`signature_count` 和 `hash_count`。
- **示例代码：**

```python
from dnakit.microbiome import load_lca_databases

result = load_lca_databases(["db.lca.json"], verbose=False)
print(result["ksize"], result["scaled"], result["databases"][0]["hash_count"])
coarse = load_lca_databases(["db.lca.json"], scaled=2000, verbose=False)
print(coarse["scaled"], coarse["databases"][0]["scaled"])
```

- **示例输出：**

```text
31 1000 25
2000 2000
```

## 9) `dnakit.microbiome.load_sbt_index_summary()`

- **作用：** 打开 sourmash 的 SBT 索引（`.sbt.zip`），列出其中的签名。
- **处理方法：** 调用 `sourmash.sbtmh.load_sbt_index()`（不显示版本警告），遍历索引中的签名，收集名称、k-mer 大小和 `scaled` 值。
- **API：** `load_sbt_index_summary(path: 'str', *, cache_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `path` 是已存在的 SBT 索引文件。`cache_size` 是可选的缓存大小，直接传给 sourmash。
- **输出：** 返回字典包含 `status`、`backend`、`path`、`index_type`（固定为 `SBT`）、`signature_count`、`signature_names`、`ksizes` 和 `scaled_values`（后两者为排序后的列表）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import load_sbt_index_summary

path = Path(tempfile.mkdtemp()) / "genomes.sbt.zip"
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(path))
result = load_sbt_index_summary(str(path), verbose=False)
print(result["index_type"], result["signature_count"], result["ksizes"], result["scaled_values"])
print(sorted(result["signature_names"]))
```

- **示例输出：**

```text
SBT 3 [31] [1000]
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 10) `dnakit.microbiome.load_sqlite_index_summary()`

- **作用：** 打开 sourmash 的 SQLite 签名集合，列出其中的签名。
- **处理方法：** 调用 `sourmash.index.sqlite_index.load_sqlite_index()`，遍历签名并收集名称、k-mer 大小和 `scaled` 值。`index_type` 是 sourmash 返回对象的类名。
- **API：** `load_sqlite_index_summary(path: 'str', *, request_manifest: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `path` 是已存在的 SQLite 签名集合。`request_manifest` 直接传给 sourmash。
- **输出：** 返回字典包含 `status`、`backend`、`path`、`index_type`、`signature_count`、`signature_names`、`ksizes` 和 `scaled_values`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.index.sqlite_index import SqliteIndex

from dnakit.microbiome import load_sqlite_index_summary

path = Path(tempfile.mkdtemp()) / "genomes.sqldb"
index = SqliteIndex.create(str(path))
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    index.insert(signature)
index.commit()
index.close()
result = load_sqlite_index_summary(str(path), verbose=False)
print(result["index_type"], result["signature_count"], result["ksizes"], result["scaled_values"])
print(sorted(result["signature_names"]))
```

- **示例输出：**

```text
SqliteIndex 3 [31] [1000]
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 11) `dnakit.microbiome.plot_comparison_matrix_composite()`

- **作用：** 为 `sourmash compare` 矩阵绘制带树状图的重排热图，并保存为图片。
- **处理方法：** 用 `sourmash.fig.load_matrix_and_labels()` 载入矩阵，再用 `sourmash.fig.plot_composite_matrix()`（单连接聚类，YlGnBu 配色）绘图，以 150 dpi 保存后关闭图形。函数返回树状图的叶子顺序和重排后的矩阵。矩阵超出 [0, 1] 范围会报错，除非设置 `force=True` 让 sourmash 重新缩放。
- **API：** `plot_comparison_matrix_composite(basefile: 'str', output_image: 'str', *, show_labels: 'bool' = True, vmin: 'float' = 0.0, vmax: 'float' = 1.0, force: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）、matplotlib（extra `dnakit[bio-plots]`）和 scipy（extra `dnakit[bio-analysis]`）。
- **输入：** `basefile` 同函数 5。`output_image` 是已存在目录中的新文件名，扩展名为 `.png`、`.pdf` 或 `.svg`。`show_labels` 控制是否显示标签。`vmin` 必须小于 `vmax`（默认 0.0 和 1.0）。
- **输出：** 返回字典包含 `status`、`backend`、`image`、`labels`（按树状图顺序）、`matrix`（重排后的行）、`vmin` 和 `vmax`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy

from dnakit.microbiome import plot_comparison_matrix_composite

work = Path(tempfile.mkdtemp())
base = work / "cmp"
with open(base, "wb") as handle:
    numpy.save(handle, numpy.array([[1.0, 0.1, 0.8], [0.1, 1.0, 0.2], [0.8, 0.2, 1.0]]))
(work / "cmp.labels.txt").write_text("genome_A\ngenome_B\ngenome_C\n")
image = work / "cmp.png"
result = plot_comparison_matrix_composite(str(base), str(image), verbose=False)
print(result["labels"])
print([[round(value, 2) for value in row] for row in result["matrix"]])
print(image.read_bytes()[:4] == b"\x89PNG", image.stat().st_size > 1000)
```

- **示例输出：**

```text
['genome_B', 'genome_A', 'genome_C']
[[1.0, 0.1, 0.2], [0.1, 1.0, 0.8], [0.2, 0.8, 1.0]]
True True
```

## 12) `dnakit.microbiome.summarize_lineage_counts()`

- **作用：** 在谱系树上统计一个签名的 k-mer：每个谱系及其所有祖先都会得到分配到其下方的 k-mer 数。
- **处理方法：** 载入数据库，要求文件中恰有一个符合其 k-mer 大小的签名，并降采样到数据库的 `scaled` 值，然后调用 `sourmash.lca.command_summarize.summarize()`。k-mer 数少于 `threshold` 的谱系不计入。结果按计数从大到小、再按谱系排序。
- **API：** `summarize_lineage_counts(signature: 'str', lca_databases: 'Sequence[str]', *, threshold: 'int' = 5, ignore_abundance: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `signature` 是已存在的文件，其中恰有一个符合数据库 k-mer 大小的签名。`lca_databases` 是非空的 LCA 数据库文件列表。`threshold` 是非负整数（sourmash 默认值为 5）。`ignore_abundance=True` 时，即使签名记录了丰度，每个哈希也只计一次。
- **输出：** 返回字典包含 `status`、`backend`、`signature`、`ksize`、`scaled`、`threshold`、`ignore_abundance` 和 `counts`：每行含 `lineage` 与 `count`。
- **示例代码：**

```python
from dnakit.microbiome import summarize_lineage_counts

result = summarize_lineage_counts("query.sig", ["db.lca.json"], threshold=2, verbose=False)
for row in result["counts"]:
    last = row["lineage"][-1]
    print(last["rank"], last["name"], row["count"])
```

- **示例输出：**

```text
superkingdom Bacteria 10
phylum Pseudomonadota 8
class Gammaproteobacteria 8
order Enterobacterales 8
family Enterobacteriaceae 8
genus Escherichia 8
species Escherichia coli 5
phylum Bacillota 2
class Bacilli 2
order Bacillales 2
family Bacillaceae 2
genus Bacillus 2
species Bacillus subtilis 2
```

## 13) `dnakit.microbiome.classify_lca_signatures()`

- **作用：** 用 LCA 数据库对多个文件中的每个签名分类，并把结果写成 CSV 文件。
- **处理方法：** 执行 `sourmash lca classify` 的计算：载入数据库（可用更粗的 `scaled`），按数据库的 k-mer 大小读取每个查询文件，必要时降采样，并对每个签名调用 `classify_signature()`。CSV 逐行写出，列为 `ID`、`status` 和每个分类等级一列。
- **API：** `classify_lca_signatures(query_files: 'Sequence[str]', lca_databases: 'Sequence[str]', output_csv: 'str', *, threshold: 'int' = 5, majority: 'bool' = False, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `query_files` 与 `lca_databases` 是非空的已存在文件列表。`output_csv` 是已存在目录中的新文件。`threshold` 是非负整数（默认 5）。`majority` 表示使用出现最多的谱系。`scaled` 是可选的正整数，用于以更粗的分辨率载入数据库。
- **输出：** 返回字典包含 `status`、`backend`、`output_csv`、`ksize`、`scaled`、`threshold`、`majority`、`classification_count` 和 `classifications`：与 CSV 内容相同的行，每行含 `id`、`status` 和 `lineage`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import classify_lca_signatures

output = Path(tempfile.mkdtemp()) / "classified.csv"
result = classify_lca_signatures(
    ["query.sig", "query2.sig"], ["db.lca.json"], str(output), verbose=False
)
print(result["classification_count"], result["ksize"], result["scaled"])
print(output.read_text().splitlines()[0])
for row in result["classifications"]:
    print(row["id"], row["status"], row["lineage"][-1]["name"])
```

- **示例输出：**

```text
2 31 1000
ID,status,superkingdom,phylum,class,order,family,genus,species,strain
query1 found Escherichia coli
query2 found Escherichia fergusonii
```

## 14) `dnakit.microbiome.compare_lineage_assignment_csvs()`

- **作用：** 比较两张表为每个标识符分配的谱系，功能等同于 `sourmash lca compare-csv`。
- **处理方法：** 用 `sourmash.lca.command_index.load_taxonomy_assignments()` 读取第一张表（classify 风格的 CSV，从第 3 列开始），再按指定的分隔符与起始列读取第二张表。对每个同时出现在两表中且谱系不同的标识符，构建谱系树并调用 `find_lca()`：一方是另一方的祖先时判为 `compatible`，否则判为 `incompatible`。不兼容的谱系按公共祖先的等级计数（没有公共祖先时按第一个等级计）。
- **API：** `compare_lineage_assignment_csvs(classify_csv: 'str', custom_csv: 'str', *, start_column: 'int' = 2, tabs: 'bool' = False, no_headers: 'bool' = False, force: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `classify_csv` 与 `custom_csv` 是已存在的文件。`start_column` 是第二张表中谱系起始列（从 1 计数，至少为 2，默认 2）。`tabs` 表示第二张表以制表符分隔。`no_headers` 表示第二张表没有表头行。`force` 直接传给 sourmash，遇到问题行时继续处理。
- **输出：** 返回字典包含 `status`、`backend`、两个路径、`classify_rows`、`custom_rows`、`classify_lineage_count`、`custom_lineage_count`、`missing_in_custom`、`missing_in_classify`、`common_count`、`different_count`、`compatible_count`、`incompatible_count`、`incompatible_by_rank`（每个等级一项）和 `differences`（每行含 `identifier`、`verdict` 和 `lca`）。
- **示例代码：**

```python
from dnakit.microbiome import compare_lineage_assignment_csvs

result = compare_lineage_assignment_csvs("classify.csv", "custom.csv", verbose=False)
print(result["classify_rows"], result["custom_rows"], result["common_count"])
print(result["missing_in_custom"], result["missing_in_classify"])
print(result["compatible_count"], result["incompatible_count"])
for row in result["differences"]:
    print(row["identifier"], row["verdict"], [pair["name"] for pair in row["lca"]][-1])
print([item for item in result["incompatible_by_rank"] if item["count"]])
```

- **示例输出：**

```text
4 4 3
['q4'] ['q5']
1 1
q2 compatible Escherichia fergusonii
q3 incompatible Bacteria
[{'rank': 'superkingdom', 'count': 1}]
```

## 15) `dnakit.microbiome.lca_rank_info()`

- **作用：** 统计最低公共祖先落在各分类等级上的数据库 k-mer 哈希数，功能等同于 `sourmash lca rankinfo`。
- **处理方法：** 载入数据库（可用更粗的 `scaled`），用 `make_lca_counts()` 收集每个谱系的计数，按每个谱系最后一项的等级汇总，并按分类顺序给出每个等级的计数与百分比。百分比不做四舍五入。
- **API：** `lca_rank_info(lca_databases: 'Sequence[str]', *, min_num: 'int' = 0, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `lca_databases` 是非空的已存在 LCA 数据库文件列表。`min_num` 是非负整数（默认 0）。`scaled` 是可选的正整数。
- **输出：** 返回字典包含 `status`、`backend`、`ksize`、`scaled`、`min_num`、`total`、`note`（为 `None`；没有任何哈希带谱系时为提示信息）和 `ranks`：每个等级一行，含 `rank`、`count` 和 `percent`。总数为 0 时 `ranks` 为空。
- **示例代码：**

```python
from dnakit.microbiome import lca_rank_info

result = lca_rank_info(["db.lca.json"], verbose=False)
print(result["total"], result["note"])
for row in result["ranks"]:
    if row["count"]:
        print(row["rank"], row["count"], round(row["percent"], 1))
```

- **示例输出：**

```text
25 None
genus 5 20.0
species 20 80.0
```

## 16) `dnakit.microbiome.summarize_lca_signatures()`

- **作用：** 在谱系树上汇总多个文件中每个签名的 k-mer，可选写出 CSV 文件，功能等同于 `sourmash lca summarize`。
- **处理方法：** 载入数据库，用 `load_singletons_and_count()` 逐个读取查询文件中的签名，用 `summarize()` 分配哈希，并让每个谱系及其所有祖先得到相应的哈希数。给出 `output_csv` 时，使用 sourmash 自己的 `output_csv()` 写出各行（列为 `count`、每个等级一列、`filename`、`sig_name`、`sig_md5`、`total_counts`）。
- **API：** `summarize_lca_signatures(query_files: 'Sequence[str]', lca_databases: 'Sequence[str]', *, output_csv: 'str | None' = None, threshold: 'int' = 5, ignore_abundance: 'bool' = False, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** `query_files` 与 `lca_databases` 是非空的已存在文件列表。`output_csv` 是可选的、位于已存在目录中的新文件。`threshold` 是非负整数（默认 5）。`ignore_abundance=True` 时每个哈希只计一次。`scaled` 是数据库可选的正整数。
- **输出：** 返回字典包含 `status`、`backend`、`output_csv`（路径或 `None`）、`ksize`、`scaled`、`threshold`、`ignore_abundance`、`signature_count` 和 `signatures`：每个签名一行，含 `filename`、`name`、`md5`、`total`（丰度之和；使用 `ignore_abundance` 时为哈希数）和 `counts`（每行含 `lineage` 与 `count`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import summarize_lca_signatures

output = Path(tempfile.mkdtemp()) / "summary.csv"
result = summarize_lca_signatures(
    ["query.sig", "query2.sig"], ["db.lca.json"], output_csv=str(output), threshold=2, verbose=False
)
print(result["signature_count"], output.read_text().splitlines()[0])
for signature in result["signatures"]:
    top = max(signature["counts"], key=lambda row: (len(row["lineage"]), row["count"]))
    print(signature["name"], signature["total"], len(signature["counts"]), top["lineage"][-1]["name"], top["count"])
```

- **示例输出：**

```text
2 count,superkingdom,phylum,class,order,family,genus,species,strain,filename,sig_name,sig_md5,total_counts
query1 10.0 13 Escherichia coli 5
query2 5.0 7 Escherichia fergusonii 5
```

## 17) `dnakit.microbiome.write_lca_index_report()`

- **作用：** 把 LCA 索引构建时发现的异常情况写成文本报告。
- **处理方法：** 把未使用的谱系转换为 sourmash 的 `LineagePair` 元组，然后调用 `sourmash.lca.command_index.generate_report()`，即 `sourmash lca index` 的 `--report` 文件所用的函数。它不构建索引，也不读取任何签名。
- **API：** `write_lca_index_report(duplicate_signatures: 'Sequence[str]', identifiers_without_lineage: 'Sequence[str]', identifiers_without_signature: 'Sequence[str]', unused_identifiers: 'Sequence[str]', unused_lineages: 'Sequence[Sequence[Mapping[str, str]]]', output_report: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。
- **输入：** 前四个参数是标识符列表（重复的签名、没有谱系的标识符、没有签名的标识符、未使用的标识符）。`unused_lineages` 是谱系列表，每个谱系是 `{"rank", "name"}` 映射的列表，从最一般的等级排到最具体的等级。`output_report` 是已存在目录中的新文件。
- **输出：** 返回字典包含 `status`、`backend`、`report`（路径）以及计数 `duplicate_count`、`no_lineage_count`、`no_signature_count`、`unused_identifier_count` 和 `unused_lineage_count`。报告文件按组列出数量与条目；未使用的谱系以 `;` 分隔的名称写出。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import write_lca_index_report

lineage = [{"rank": "superkingdom", "name": "Bacteria"}, {"rank": "phylum", "name": "Bacillota"}]
report = Path(tempfile.mkdtemp()) / "report.txt"
result = write_lca_index_report(
    ["GCF_A"], ["GCF_X"], ["GCF_Y"], ["GCF_Z"], [lineage], str(report), verbose=False
)
print({key: value for key, value in result.items() if key.endswith("_count")})
print(report.read_text().rstrip())
```

- **示例输出：**

```text
{'duplicate_count': 1, 'no_lineage_count': 1, 'no_signature_count': 1, 'unused_identifier_count': 1, 'unused_lineage_count': 1}
Duplicate signatures: 1
GCF_A
----
Unused identifiers: 1
GCF_Z
----
No lineage provided for these identifiers: 1
GCF_X
----
No signatures found for these identifiers: 1
GCF_Y
----
Unused lineages: 1
Bacteria;Bacillota;;;;;;
```

## 18) `dnakit.microbiome.sketch_fromfile_signatures()`

- **作用：** 对 CSV 文件中列出的每个基因组做草图，并把签名保存到一个新的集合中。
- **处理方法：** 执行 `sourmash sketch fromfile CSV... --output-signatures OUT --license LICENSE`，每个参数字符串对应一个 `-p`；需要时加上 `--already-done` 和 `--check-sequence`。命令在当前工作目录中运行，因此 CSV 里的相对基因组路径从当前工作目录解析。集合的格式由 `output_signatures` 的后缀决定（例如 `.zip`）。
- **API：** `sketch_fromfile_signatures(csv_paths: 'Sequence[str]', output_signatures: 'str', *, param_strings: 'Sequence[str]' = ('dna,k=31,scaled=1000',), already_done: 'Sequence[str]' = (), license: 'str' = 'CC0', check_sequence: 'bool' = False, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `csv_paths` 是非空的 CSV 文件列表，列为 `name`、`genome_filename` 和 `protein_filename`（蛋白路径可为空）。`output_signatures` 是已存在目录中的新文件。`param_strings` 是非空的 sourmash 参数字符串列表，例如 `dna,k=31,scaled=1000`（默认值）。`already_done` 列出已存在的集合，其中已有的签名不会重新计算。`license` 是签名许可（默认 `CC0`）。`check_sequence` 让 sourmash 检查序列。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import sketch_fromfile_signatures

output = Path(tempfile.mkdtemp()) / "sketched.zip"
result = sketch_fromfile_signatures(
    ["genomes.csv"], str(output), param_strings=["dna,k=21,scaled=1"], verbose=False
)
print(result["status"], result["return_code"], len(result["outputs"]))
for signature in sourmash.load_file_as_signatures(str(output)):
    print(signature.name, signature.minhash.ksize, signature.minhash.scaled, len(signature.minhash))
```

- **示例输出：**

```text
completed 0 1
genome1 21 1 73
```

## 19) `dnakit.microbiome.grep_signatures()`

- **作用：** 保留名称匹配正则表达式的签名，并保存到新文件。
- **处理方法：** 执行 `sourmash sig grep PATTERN SIGNATURE... -o OUT`，按已设置的选项加上 `-i`、`-v`、`-k`、分子类型标志和 `--csv`。
- **API：** `grep_signatures(pattern: 'str', signatures: 'Sequence[str]', output: 'str', *, ignore_case: 'bool' = False, invert_match: 'bool' = False, ksize: 'int | None' = None, moltype: 'str | None' = None, csv_manifest: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `pattern` 是非空的正则表达式。`signatures` 是非空的已存在签名文件列表。`output` 是已存在目录中的新文件。`ignore_case` 表示忽略大小写。`invert_match` 表示保留不匹配的签名。`ksize` 是正整数。`moltype` 取 `DNA`、`RNA`、`protein`、`dayhoff` 或 `hp` 之一。`csv_manifest` 是接收所保留签名 CSV 清单的新文件。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import grep_signatures

work = Path(tempfile.mkdtemp())
result = grep_signatures("escherichia", ["genomes.sig"], str(work / "kept.sig"), ignore_case=True, verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print([s.name for s in sourmash.load_file_as_signatures(str(work / "kept.sig"))])
result = grep_signatures("escherichia", ["genomes.sig"], str(work / "rest.sig"), ignore_case=True, invert_match=True, verbose=False)
print([s.name for s in sourmash.load_file_as_signatures(str(work / "rest.sig"))])
```

- **示例输出：**

```text
completed ['kept.sig']
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii']
['GCF_C Bacillus subtilis']
```

## 20) `dnakit.microbiome.inflate_signatures()`

- **作用：** 让其他签名获得记录丰度的参考签名中的丰度。
- **处理方法：** 执行 `sourmash sig inflate REFERENCE OTHER... -o OUT`，带 `-k KSIZE`；未给出 `ksize` 时改用 `--dna`。上游对每个其他签名调用 `MinHash.inflate()`，因此只有同时出现在两个签名中的哈希才保留丰度，输出签名会记录丰度。
- **API：** `inflate_signatures(signature_from: 'str', other_signatures: 'Sequence[str]', output: 'str', *, ksize: 'int | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `signature_from` 是带丰度的已存在签名文件。`other_signatures` 是非空的已存在签名文件列表。`output` 是已存在目录中的新文件。`ksize` 是正整数。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import inflate_signatures

output = Path(tempfile.mkdtemp()) / "inflated.sig"
result = inflate_signatures("abund.sig", ["genomes.sig"], str(output), ksize=31, verbose=False)
print(result["status"], result["return_code"])
for signature in sourmash.load_file_as_signatures(str(output)):
    print(signature.name, signature.minhash.track_abundance, sorted(signature.minhash.hashes.items()))
```

- **示例输出：**

```text
completed 0
GCF_A Escherichia coli True [(1, 5), (2, 3), (3, 1)]
GCF_B Escherichia fergusonii True []
GCF_C Bacillus subtilis True [(21, 2), (30, 7)]
```

## 21) `dnakit.microbiome.combine_sbt_indexes()`

- **作用：** 把多个 SBT 索引合并成一个新的 SBT 索引。
- **处理方法：** 执行 `sourmash sbt_combine OUT INPUT... -x BF_SIZE`。
- **API：** `combine_sbt_indexes(output_sbt: 'str', input_sbts: 'Sequence[str]', *, bf_size: 'float' = 100000.0, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `output_sbt` 是已存在目录中的新文件，通常以 `.sbt.zip` 结尾。`input_sbts` 是非空的已存在 SBT 索引文件列表。`bf_size` 是保存索引时使用的正数 Bloom 过滤器大小（默认 100000）。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import combine_sbt_indexes

work = Path(tempfile.mkdtemp())
signatures = list(sourmash.load_file_as_signatures("genomes.sig"))
inputs = []
for number, chunk in enumerate((signatures[:1], signatures[1:])):
    tree = create_sbt_index(1000)
    for signature in chunk:
        tree.insert(signature)
    inputs.append(str(work / f"part{number}.sbt.zip"))
    tree.save(inputs[-1])
result = combine_sbt_indexes(str(work / "all.sbt.zip"), inputs, verbose=False)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
merged = sourmash.load_file_as_index(str(work / "all.sbt.zip"))
print(sorted(signature.name for signature in merged.signatures()))
```

- **示例输出：**

```text
completed 0 ['all.sbt.zip']
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 22) `dnakit.microbiome.migrate_sbt_index()`

- **作用：** 把 SBT 索引改写为最新格式，原文件保持不变。
- **处理方法：** 必要时创建 `output_dir`，以相同文件名把索引复制进去（`.sbt.zip` 索引的存储标签取决于文件名，改名后的副本会丢失签名），再对副本运行 `sourmash migrate`。如果索引已是最新格式，则只是重写这份副本。
- **API：** `migrate_sbt_index(sbt_path: 'str', output_dir: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `sbt_path` 是已存在的自包含 `.sbt.zip` 索引。`output_dir` 是存放副本的目录；其中同名文件必须尚不存在。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import migrate_sbt_index

work = Path(tempfile.mkdtemp())
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(work / "genomes.sbt.zip"))
result = migrate_sbt_index(str(work / "genomes.sbt.zip"), str(work / "migrated"), verbose=False)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).parent.name)
migrated = sourmash.load_file_as_index(result["outputs"][0])
print(len(list(migrated.signatures())))
```

- **示例输出：**

```text
completed 0 ['genomes.sbt.zip']
migrated
3
```

## 23) `dnakit.microbiome.convert_sbt_storage()`

- **作用：** 把 SBT 索引的存储转换为另一种后端，在副本上进行。
- **处理方法：** 把索引复制到 `output_sbt`，再对副本运行 `sourmash storage convert OUT -b BACKEND(PATH)`。目标路径会被解析为绝对路径，且必须尚不存在。
- **API：** `convert_sbt_storage(sbt_path: 'str', output_sbt: 'str', backend: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `sbt_path` 是已存在的 SBT 索引文件。`output_sbt` 是已存在目录中的新文件。`backend` 为 `zip(PATH)` 或 `fs(PATH)`，其中 PATH 是接收签名数据的新文件或目录；其他形式会抛出 `ConfigurationError`。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。 转换后的索引和新的存储目标都属于预期输出。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import convert_sbt_storage

work = Path(tempfile.mkdtemp())
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(work / "genomes.sbt.zip"))
result = convert_sbt_storage(
    str(work / "genomes.sbt.zip"), str(work / "converted.sbt.zip"), f"zip({work / 'signatures.zip'})", verbose=False
)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
```

- **示例输出：**

```text
completed 0 ['converted.sbt.zip', 'signatures.zip']
```

## 24) `dnakit.microbiome.krona_taxonomy_profile()`

- **作用：** 写出 Krona 表，给出每个查询在某一等级上分配到各谱系的比例。
- **处理方法：** 执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F krona -o BASE --rank RANK`。输出文件为 `BASE.krona.tsv`；谱系在该等级上的汇总由 sourmash 完成。
- **API：** `krona_taxonomy_profile(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str', sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。 `rank` 必填，且必须是 `superkingdom`、`phylum`、`class`、`order`、`family`、`genus`、`species` 或 `strain` 之一。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import krona_taxonomy_profile

base = Path(tempfile.mkdtemp()) / "profile"
result = krona_taxonomy_profile(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **示例输出：**

```text
completed ['profile.krona.tsv']
fraction	superkingdom	phylum	class	order	family	genus
0.8	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia
0.2	Bacteria	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus
```

## 25) `dnakit.microbiome.csv_taxonomy_summary()`

- **作用：** 把 gather 结果的分等级汇总写成 CSV 文件。
- **处理方法：** 执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F csv_summary -o BASE`。输出文件为 `BASE.summarized.csv`，列包括 `query_name`、`rank`、`fraction`、`lineage`、查询标识和匹配统计。
- **API：** `csv_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import csv_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = csv_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
lines = Path(result["outputs"][0]).read_text().splitlines()
print(lines[0])
for line in lines[1:]:
    fields = line.split(",")
    if fields[1] == "genus":
        print(fields[0], fields[1], fields[2], fields[3].split(";")[-1])
```

- **示例输出：**

```text
completed ['profile.summarized.csv']
query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank,total_weighted_hashes
query1 genus 0.8 Escherichia
query1 genus 0.2 Bacillus
```

## 26) `dnakit.microbiome.human_taxonomy_summary()`

- **作用：** 写出一张人类可读的表，列出每个查询中发现的谱系。
- **处理方法：** 执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F human -o BASE`，给出等级时加上 `--rank RANK`。输出文件为 `BASE.human.txt`。不指定等级时，sourmash 报告它找到的谱系。
- **API：** `human_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。 `rank` 可选，给出时必须是函数 24 中列出的八个等级之一。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import human_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = human_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **示例输出：**

```text
completed ['profile.human.txt']
sample name    proportion   cANI   lineage
-----------    ----------   ----   -------
query1            80.0%     99.3%  Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia
query1            20.0%     94.9%  Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus
```

## 27) `dnakit.microbiome.lineage_taxonomy_summary()`

- **作用：** 把每个查询在某一等级上分配到各谱系的比例写成表格，每个查询一列。
- **处理方法：** 执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F lineage_summary -o BASE --rank RANK`。输出文件为 `BASE.lineage_summary.tsv`。
- **API：** `lineage_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str', sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。 `rank` 必填，且是函数 24 中列出的八个等级之一。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import lineage_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = lineage_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **示例输出：**

```text
completed ['profile.lineage_summary.tsv']
lineage	query1
Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus	0.2
Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia	0.8
```

## 28) `dnakit.microbiome.bioboxes_taxonomy_profile()`

- **作用：** 写出 gather 结果的 CAMI bioboxes 分类学谱。
- **处理方法：** 先检查每个谱系 CSV 的表头含有 `taxpath` 列，然后执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F bioboxes -o BASE`。输出文件为 `BASE.bioboxes.profile`。`taxpath` 列保存谱系的 NCBI taxid，以 `|` 分隔，每个等级一个。在 sourmash 4.9.4 中，如果被报告的谱系没有 taxid，命令会因 `TypeError` 失败，此时函数返回状态 `failed` 和日志末尾。
- **API：** `bioboxes_taxonomy_profile(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。 每个谱系 CSV 都需要 `taxpath` 列，否则抛出 `ConfigurationError`。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import bioboxes_taxonomy_profile

base = Path(tempfile.mkdtemp()) / "profile"
result = bioboxes_taxonomy_profile(["gather.csv"], ["taxid.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
lines = Path(result["outputs"][0]).read_text().splitlines()
print(lines[1:5])
for line in lines[6:]:
    fields = line.split("\t")
    if fields[1] == "species":
        print(fields[0], fields[1], fields[2], fields[4])
```

- **示例输出：**

```text
completed ['profile.bioboxes.profile']
['@SampleID:query1', '@Version:0.10.0', '@Ranks:superkingdom|phylum|class|order|family|genus|species|strain', '@__program__:sourmash']
562 species 2|1224|1236|91347|543|561|562 80.00
1423 species 2|1239|91061|1385|186817|1386|1423 20.00
```

## 29) `dnakit.microbiome.kreport_taxonomy_report()`

- **作用：** 写出 gather 结果的 Kraken 风格报告。
- **处理方法：** 执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F kreport -o BASE`。输出文件为 `BASE.kreport.txt`，列为百分比、包含的碱基数、直接归属的碱基数、等级代码、taxid 和名称。taxid 列为空，查询中未分类的部分以 `U` 报告。
- **API：** `kreport_taxonomy_report(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 与 `taxonomy_csvs` 是非空的已存在文件列表：`sourmash gather` 生成的 gather CSV，以及谱系 CSV（一个 `ident` 列加每个等级一列）。`output_base` 是位于已存在目录中的路径前缀；输出文件为 `output_base` 加上上述后缀，且必须尚不存在。 gather CSV 必须带有非零的 `total_weighted_hashes` 列（sourmash 4.5 及以后版本写出）；否则 sourmash 以 `ValueError` 终止，状态为 `failed`，并留下一个空的输出文件。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import kreport_taxonomy_report

base = Path(tempfile.mkdtemp()) / "profile"
result = kreport_taxonomy_report(["gather.csv"], ["tax.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
for line in Path(result["outputs"][0]).read_text().splitlines():
    if line.split("\t")[3] in ("D", "S"):
        print(line)
```

- **示例输出：**

```text
completed ['profile.kreport.txt']
100.00	10000	0	D		Bacteria
80.00	8000	8000	S		Escherichia coli
20.00	2000	2000	S		Bacillus subtilis
```

## 30) `dnakit.microbiome.lingroup_taxonomy_report()`

- **作用：** 按 LIN 谱系组汇总 gather 结果。
- **处理方法：** 先检查每个谱系 CSV 的表头含有 `lin` 列，然后执行 `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F lingroup -o BASE --lins --lingroup LINGROUP_CSV`。输出文件为 `BASE.lingroup.tsv`，列为 `name`、`lin`、`percent_containment` 和 `num_bp_contained`。
- **API：** `lingroup_taxonomy_report(gather_csvs: 'Sequence[str]', lin_taxonomy_csvs: 'Sequence[str]', lingroup_csv: 'str', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 sourmash（extra `dnakit[bio-sourmash]`）。函数会运行 `sourmash` 命令，该命令随包一起安装；查找顺序为 PATH、Python 解释器所在目录，或使用 `sourmash_path`。找不到时抛出 `BackendUnavailableError`。
- **输入：** `gather_csvs` 是非空的 gather CSV 列表；`lin_taxonomy_csvs` 是非空的谱系 CSV 列表，需包含 `ident` 与 `lin` 列（否则抛出 `ConfigurationError`）。`lingroup_csv` 是已存在的 CSV，列为 `name` 和 `lin`，给出每个组的 LIN 前缀。`output_base` 是位于已存在目录中的路径前缀；`BASE.lingroup.tsv` 必须尚不存在。与函数 29 相同，gather CSV 需要非零的 `total_weighted_hashes` 列。
- **输出：** 返回字典包含 `status`（命令退出码为 0 且所有预期输出都存在时为 `completed`，否则为 `failed`）、`backend`、`command`（命令参数，程序名显示为 `sourmash`）、`return_code`、`outputs`（实际存在的预期输出路径）和 `log_tail`（命令输出的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import lingroup_taxonomy_report

base = Path(tempfile.mkdtemp()) / "profile"
result = lingroup_taxonomy_report(["gather.csv"], ["lin.csv"], "lingroup.csv", str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
header, *rows = Path(result["outputs"][0]).read_text().splitlines()
print(header)
for row in sorted(rows):
    print(row)
```

- **示例输出：**

```text
completed ['profile.lingroup.tsv']
name	lin	percent_containment	num_bp_contained
group_bacillus_like	1	20.00	2000
group_ecoli_like	0	80.00	8000
```
