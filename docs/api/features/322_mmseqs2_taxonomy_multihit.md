# MMseqs2 分类学归类与多命中搜索

本组函数位于 `dnakit.alignment`，是对 [MMseqs2](https://github.com/soedinglab/MMseqs2) 命令行程序的封装：每个函数以固定参数列表（不经 shell）运行一个 MMseqs2 模块。需要用户自行安装 MMseqs2（conda 或发布二进制）并放到 PATH，或传入 `mmseqs_path`；测试所用版本为提交 76de5d8ae114750adc59288516d947b6a72a4b05。数据库用路径前缀表示；未设置的选项不传给 MMseqs2，沿用其默认值；新输出必须不存在。示例通过模块路径导入。

## 示例数据

**示例文件：** `proteins.fa`

```text
>fam0_m0
MMFPCDVENWCTHCDQQDIDVQCWEIWCWWPCICVFLQFVEWLVGEWWHNEVDWCYHSVQ
>fam0_m1
MMFPCDVLNWCTHCDQQDDDVQCWEIWCWWPLICVFLQRVEWLVGEWWHNEVDWCYHSVF
>fam0_m2
MMFPCQVENWCTHCDQQDIDVQCWEIWCWWPCIHVFLQFVEWLVEEWWHNEVDWCYHSKQ
>fam1_m0
MGNIVVTMIYHIPIHTSNAAKSKHYNRNNDIEISHMHSYYASNDEPHSGQMDPRPDGGFA
>fam1_m1
MGNIVVTMIYEIPIHTSNAAKSKHYNRNNDIEISRMHSYYASTDEPHHGQMDPRPDGGFA
>fam1_m2
MGNIVVTMIYHIPIHTSNAAKSKHYNSNNDIEISHNVSYYASNDEPHKGQMDPRPDGGFC
>fam2_m0
MGQDKADKDYIDKERAMVQKYFCTIEGKCGHLLTHLRTGKNAKCAATVHTSIREQSVPTL
>fam2_m1
MGQDKADKFKIPKERAMVQKYFCTIEGKPKHLLKHLRTGKTAKCAATVHTSIREQSVPTL
>fam2_m2
MGWDKADCDYIDKERAMDQKYFSTKEGKCGHLLTHFRTGKNAKCHATVHTSIREQRTPTL
>fam3_m0
MHHDWDFTKNFYTKENISSPAGASRPLFQNPMEMAMMPEHALKNDPPWDNQKCKECLFIK
>fam3_m1
MHHDWPFTKNQYTKEVISSPAGASRPLFQNPMEMAMMPMHALKNQPPWDNQKCKECLIIK
>fam3_m2
MHHDWDFTKNTYTKENISSPAGASFPLFQNPMEMTMMPEVALKNDPPWDNQSCKICLFMK
```

**示例文件：** `contigs.fa`

```text
>contig0
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAATAAAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCAGAATGAGGTTCATCATTAGAAGCATAATAAGAATGCATATGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATATGATAAATCATAGTAACAACAATATTACCCATTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>contig1
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAATTATTTAATAAACAGACATTCTTTACATTTTTGATTATCCCAAGGAGGATCATTTTTCAGAGCATGTTCAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATATTTTCTTTAGTATAAAAATTTTTAGTAAAATCCCAATCATGATGCATGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
>contig2
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGAATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTTTAACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGATTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCATGATGAGGTTCATCAGTAGAAGCATAATAAGAATGCATACGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATTTCATAAATCATAGTAACAACAATATTACCCATATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>contig3
TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTAATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAATTATTTAATAATCAGACATTCTTTACATTTTTGATTATCCCAAGGAGGTTGATTTTTCAGAGCATGCATAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATAACTTCTTTAGTATATTGATTTTTAGTAAAAGGCCAATCATGATGCATATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGAT
```

**示例文件：** `taxmap.tsv`

```text
fam0_m0	561
fam0_m1	561
fam0_m2	561
fam1_m0	561
fam1_m1	561
fam1_m2	561
fam2_m0	9605
fam2_m1	9605
fam2_m2	9605
fam3_m0	9605
fam3_m1	9605
fam3_m2	9605
```

**示例文件：** `nodes.dmp`

```text
1	|	1	|	no rank	|		|
2	|	1	|	superkingdom	|		|
2759	|	1	|	superkingdom	|		|
561	|	2	|	genus	|		|
9605	|	2759	|	genus	|		|
```

**示例文件：** `names.dmp`

```text
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
2759	|	Eukaryota	|		|	scientific name	|
561	|	Escherichia	|		|	scientific name	|
9605	|	Homo	|		|	scientific name	|
```

**示例文件：** `merged.dmp`

```text
9999	|	561	|
```

**示例文件：** `set_a.fa`

```text
>contig0
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAATAAAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCAGAATGAGGTTCATCATTAGAAGCATAATAAGAATGCATATGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATATGATAAATCATAGTAACAACAATATTACCCATTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>contig1
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAATTATTTAATAAACAGACATTCTTTACATTTTTGATTATCCCAAGGAGGATCATTTTTCAGAGCATGTTCAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATATTTTCTTTAGTATAAAAATTTTTAGTAAAATCCCAATCATGATGCATGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
```

**示例文件：** `set_b.fa`

```text
>contig2
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGAATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTTTAACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGATTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCATGATGAGGTTCATCAGTAGAAGCATAATAAGAATGCATACGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATTTCATAAATCATAGTAACAACAATATTACCCATATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>contig3
TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTAATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAATTATTTAATAATCAGACATTCTTTACATTTTTGATTATCCCAAGGAGGTTGATTTTTCAGAGCATGCATAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATAACTTCTTTAGTATATTGATTTTTAGTAAAAGGCCAATCATGATGCATATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGAT
```

## 1) `dnakit.alignment.mmseqs_createtaxdb()`

- **作用：** 为序列数据库加入分类学标签（原地写入分类文件）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_createtaxdb(sequence_db, *, ncbi_tax_dump=None, tax_mapping_file=None, tax_mapping_mode=None, tax_db_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
wf.mmseqs_createdb("proteins.fa", "work/ref", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/ctg", verbose=False)
r = tx.mmseqs_createtaxdb("work/ref", ncbi_tax_dump=".", tax_mapping_file="taxmap.tsv", verbose=False)
al.mmseqs_extractorfs("work/ctg", "work/orfs", translate=True, min_length=30, verbose=False)
al.mmseqs_prefilter("work/orfs", "work/ref", "work/orfpf", verbose=False)
al.mmseqs_rescorediagonal("work/orfs", "work/ref", "work/orfpf", "work/orfaln", rescore_mode=2, add_backtrace=True, verbose=False)
print([Path(p).name for p in r["taxonomy_files"]])
```

- **示例输出：**

```text
['ref_mapping', 'ref_taxonomy']
```

## 2) `dnakit.alignment.mmseqs_taxonomy()`

- **作用：** 对序列做分类学归类（FASTA 或数据库两种输入）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_taxonomy(query, target_db, output, *, input_kind='fasta', lca_mode=None, majority=None, vote_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query`, `target_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_taxonomy("contigs.fa", "work/ref", "work/tax", verbose=False)
print([(row["query"], row["name"]) for row in r["records"]])
d = tx.mmseqs_taxonomy("work/orfs", "work/ref", "work/taxdb", input_kind="database", tax_output_mode=2, verbose=False)
print(d["entry_counts"])
```

- **示例输出：**

```text
[('contig0', 'Escherichia'), ('contig1', 'Homo'), ('contig2', 'Escherichia'), ('contig3', 'unclassified')]
{'taxonomy_db': 30}
```

## 3) `dnakit.alignment.mmseqs_taxonomyreport()`

- **作用：** 由分类结果生成 Kraken 或 Krona 报告。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs taxonomyreport`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_taxonomyreport(sequence_db, result_db, output, *, report_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `result_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_taxonomyreport("work/ref", "work/taxdb", "work/report.txt", verbose=False)
for row in r["rows"]:
    print(row[3:])
```

- **示例输出：**

```text
['no rank', '0', 'unclassified']
['no rank', '1', 'root']
['superkingdom', '2', '  Bacteria']
['genus', '561', '    Escherichia']
['superkingdom', '2759', '  Eukaryota']
['genus', '9605', '    Homo']
```

## 4) `dnakit.alignment.mmseqs_addtaxonomy()`

- **作用：** 给结果数据库加上分类学列。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs addtaxonomy`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_addtaxonomy(target_db, result_db, annotated_db, *, tax_lineage=None, lca_ranks=None, pick_id_from=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `target_db`, `result_db`, `annotated_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_addtaxonomy("work/ref", "work/orfaln", "work/annotated", tax_lineage=1, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'annotated_db': 8}
```

## 5) `dnakit.alignment.mmseqs_filtertaxdb()`

- **作用：** 按分类单元表达式过滤分类结果。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs filtertaxdb`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_filtertaxdb(target_db, taxonomy_db, filtered_db, *, taxon_list=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `target_db`, `taxonomy_db`, `filtered_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_filtertaxdb("work/ref", "work/taxdb", "work/only_561", taxon_list="561", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'filtered_db': 30}
```

## 6) `dnakit.alignment.mmseqs_filtertaxseqdb()`

- **作用：** 按分类单元表达式过滤序列数据库。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs filtertaxseqdb`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_filtertaxseqdb(sequence_db, filtered_db, *, taxon_list=None, subdb_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `filtered_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_filtertaxseqdb("work/ref", "work/only_homo", taxon_list="9605", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'filtered_db': 6}
```

## 7) `dnakit.alignment.mmseqs_majoritylca()`

- **作用：** 按多数投票计算最近共同祖先。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs majoritylca`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_majoritylca(target_db, result_db, taxonomy_db, *, majority=None, vote_mode=None, lca_ranks=None, blacklist=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `target_db`, `result_db`, `taxonomy_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_majoritylca("work/ref", "work/orfaln", "work/mlca", majority=0.6, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'taxonomy_db': 30}
```

## 8) `dnakit.alignment.mmseqs_lcaalign()`

- **作用：** 为 LCA 计算做有空位比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs lcaalign`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_lcaalign(query_db, target_db, result_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `result_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_lcaalign("work/orfs", "work/ref", "work/orfpf", "work/lcaaln", add_backtrace=True, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 30}
```

## 9) `dnakit.alignment.mmseqs_aggregate_taxonomy()`

- **作用：** 把同一集合（如同一重叠群的 ORF）的分类标签聚合为一个标签。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_aggregate_taxonomy(sequence_db, set_to_sequence_db, taxonomy_per_sequence_db, output, *, alignment_per_sequence_db=None, swap_headers=False, majority=None, vote_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `set_to_sequence_db`, `taxonomy_per_sequence_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = tx.mmseqs_aggregate_taxonomy("work/ref", "work/orfs_h", "work/taxdb", "work/per_contig", swap_headers=True, verbose=False)
w = tx.mmseqs_aggregate_taxonomy("work/ref", "work/orfs_h", "work/taxdb", "work/per_contig_w", alignment_per_sequence_db="work/taxdb_aln", swap_headers=True, majority=0.6, verbose=False)
t = wf.mmseqs_createtsv("work/ctg", "work/per_contig", "work/per_contig.tsv", target_column=0, verbose=False)
print(r["entry_counts"], w["entry_counts"])
for row in t["rows"]:
    print(row[0], row[2], row[4])
```

- **示例输出：**

```text
{'taxonomy_db': 4} {'taxonomy_db': 4}
contig0 561 Escherichia
contig1 9605 Homo
contig2 561 Escherichia
contig3 0 unclassified
```

## 10) `dnakit.alignment.mmseqs_multihitdb()`

- **作用：** 为多命中搜索建立分组（集合）数据库，每个输入文件为一个集合。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs multihitdb`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_multihitdb(fasta_files, set_db, *, dbtype=None, shuffle=None, createdb_mode=None, min_length=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `fasta_files`, `set_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

a = tx.mmseqs_multihitdb(["set_a.fa", "set_b.fa"], "work/sets1", min_length=30, verbose=False)
tx.mmseqs_multihitdb(["set_b.fa", "set_a.fa"], "work/sets2", min_length=30, verbose=False)
print(a["entry_counts"])
```

- **示例输出：**

```text
{'set_db': 30}
```

## 11) `dnakit.alignment.mmseqs_multihitsearch()`

- **作用：** 分组查询集合对分组目标集合的多命中搜索。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs multihitsearch`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_multihitsearch(query_set_db, target_set_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_set_db`, `target_set_db`, `result_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_multihitsearch("work/sets1", "work/sets2", "work/hits", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'result_db': 2}
```

## 12) `dnakit.alignment.mmseqs_besthitperset()`

- **作用：** 每个集合保留最佳命中并更新 p 值。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs besthitperset`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_besthitperset(query_set_db, target_set_db, result_db, best_hit_db, *, simple_best_hit=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_set_db`, `target_set_db`, `result_db`, `best_hit_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

al.mmseqs_prefilter("work/sets1", "work/sets2", "work/spf", verbose=False)
al.mmseqs_rescorediagonal("work/sets1", "work/sets2", "work/spf", "work/saln", rescore_mode=2, verbose=False)
r = tx.mmseqs_besthitperset("work/sets1", "work/sets2", "work/saln", "work/best", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'best_hit_db': 30}
```

## 13) `dnakit.alignment.mmseqs_combinepvalperset()`

- **作用：** 合并每对集合的命中 p 值。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs combinepvalperset`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_combinepvalperset(query_set_db, target_set_db, result_db, pvalue_db, *, alpha=None, aggregation_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_set_db`, `target_set_db`, `result_db`, `pvalue_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_combinepvalperset("work/sets1", "work/sets2", "work/hits", "work/pvalues", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'pvalue_db': 2}
```

## 14) `dnakit.alignment.mmseqs_mergeresultsbyset()`

- **作用：** 把各 ORF 的结果合并回其所属集合。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs mergeresultsbyset`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_mergeresultsbyset(set_to_member_db, result_db, merged_db, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `set_to_member_db`, `result_db`, `merged_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_mergeresultsbyset("work/sets1_set_to_member", "work/best", "work/merged", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'merged_db': 2}
```
