# MMseqs2 数据库、搜索、聚类与导出流程

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

**示例文件：** `features.gff3`

```text
##gff-version 3
contig0	toy	CDS	41	223	.	+	0	ID=cds0
contig0	toy	CDS	264	446	.	-	0	ID=cds1
contig1	toy	CDS	41	223	.	+	0	ID=cds2
contig1	toy	CDS	264	446	.	-	0	ID=cds3
contig2	toy	CDS	41	223	.	+	0	ID=cds4
contig2	toy	CDS	264	446	.	-	0	ID=cds5
contig3	toy	CDS	41	223	.	+	0	ID=cds6
contig3	toy	CDS	264	446	.	-	0	ID=cds7
```

## 1) `dnakit.alignment.mmseqs_createdb()`

- **作用：** 把 FASTA/FASTQ 文件转换为 MMseqs2 序列数据库。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs createdb`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_createdb(fasta_files, database, *, dbtype=None, shuffle=None, createdb_mode=None, write_lookup=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `fasta_files`, `database`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
result = wf.mmseqs_createdb("proteins.fa", "work/prot", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/nt", dbtype=2, verbose=False)
al.mmseqs_prefilter("work/prot", "work/prot", "work/pf", verbose=False)
al.mmseqs_rescorediagonal("work/prot", "work/prot", "work/pf", "work/aln", rescore_mode=2, add_backtrace=True, verbose=False)
print(result["status"], result["entry_counts"])
```

- **示例输出：**

```text
completed {'database': 12}
```

## 2) `dnakit.alignment.mmseqs_map()`

- **作用：** 把近乎相同的序列映射到目标数据库。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs map`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_map(query_db, target_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_map("work/prot", "work/prot", "work/map", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 12}
```

## 3) `dnakit.alignment.mmseqs_convertalis()`

- **作用：** 把比对数据库导出为 BLAST-tab、SAM 或自定义列表格。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs convertalis`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_convertalis(query_db, target_db, alignment_db, output, *, format_mode=None, format_output=None, search_type=None, translation_table=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `alignment_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_convertalis("work/prot", "work/prot", "work/aln", "work/aln.tsv", format_output="query,target,fident,bits", verbose=False)
print(r["line_counts"])
print(r["records"][1])
```

- **示例输出：**

```text
{'output': 36}
{'query': 'fam0_m0', 'target': 'fam0_m2', 'fident': '0.933', 'bits': '173'}
```

## 4) `dnakit.alignment.mmseqs_clust()`

- **作用：** 对比对或预过滤结果图做聚类。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs clust`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_clust(sequence_db, result_db, cluster_db, *, cluster_mode=None, max_iterations=None, similarity_type=None, set_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `result_db`, `cluster_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_clust("work/prot", "work/aln", "work/clu", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'cluster_db': 4}
```

## 5) `dnakit.alignment.mmseqs_createtsv()`

- **作用：** 把结果数据库导出为制表符分隔文件。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs createtsv`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_createtsv(query_db, result_db, output, *, target_db=None, first_seq_as_repr=None, target_column=None, full_header=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `result_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_createtsv("work/prot", "work/clu", "work/clu.tsv", target_db="work/prot", verbose=False)
print(r["rows"][:3])
```

- **示例输出：**

```text
[['fam0_m2', 'fam0_m2'], ['fam0_m2', 'fam0_m0'], ['fam0_m2', 'fam0_m1']]
```

## 6) `dnakit.alignment.mmseqs_convert2fasta()`

- **作用：** 把序列数据库写成 FASTA 文件。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs convert2fasta`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_convert2fasta(sequence_db, output, *, use_header_file=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_convert2fasta("work/prot", "work/prot.fa", verbose=False)
print(r["sequence_count"])
```

- **示例输出：**

```text
12
```

## 7) `dnakit.alignment.mmseqs_result2repseq()`

- **作用：** 提取每个簇的代表序列。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs result2repseq`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_result2repseq(sequence_db, result_db, representative_db, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `result_db`, `representative_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_result2repseq("work/prot", "work/clu", "work/rep", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'representative_db': 4}
```

## 8) `dnakit.alignment.mmseqs_pickconsensusrep()`

- **作用：** 按簇共识重新挑选每簇代表序列。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs pickconsensusrep`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_pickconsensusrep(sequence_db, cluster_db, new_cluster_db, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `cluster_db`, `new_cluster_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_pickconsensusrep("work/prot", "work/clu", "work/clu_consensus", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'new_cluster_db': 4}
```

## 9) `dnakit.alignment.mmseqs_clusthash()`

- **作用：** 对等长序列做哈希聚类。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs clusthash`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_clusthash(sequence_db, alignment_db, *, min_seq_id=None, max_seq_len=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_clusthash("work/prot", "work/hash", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 12}
```

## 10) `dnakit.alignment.mmseqs_clusterupdate()`

- **作用：** 用新序列增量更新已有聚类（旧簇数据库会被 MMseqs2 删除）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs clusterupdate`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_clusterupdate(old_sequence_db, new_sequence_db, old_cluster_db, mapped_sequence_db, new_cluster_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `old_sequence_db`, `new_sequence_db`, `old_cluster_db`, `mapped_sequence_db`, `new_cluster_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_clust("work/prot", "work/aln", "work/old_clu", verbose=False)
r = wf.mmseqs_clusterupdate("work/prot", "work/prot", "work/old_clu", "work/mapped", "work/new_clu", verbose=False)
print(r["entry_counts"], r["old_cluster_db_removed"])
```

- **示例输出：**

```text
{'mapped_sequence_db': 12, 'new_cluster_db': 4} True
```

## 11) `dnakit.alignment.mmseqs_summarizealis()`

- **作用：** 把每条查询序列的比对汇总为一行。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs summarizealis`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_summarizealis(alignment_db, summary_db, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `alignment_db`, `summary_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_summarizealis("work/aln", "work/summary", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'summary_db': 12}
```

## 12) `dnakit.alignment.mmseqs_filterresult()`

- **作用：** 按冗余度、一致性、得分与覆盖度过滤比对结果。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs filterresult`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_filterresult(query_db, target_db, result_db, filtered_db, *, max_seq_id=None, qid=None, qsc=None, cov=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `result_db`, `filtered_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_filterresult("work/prot", "work/prot", "work/aln", "work/filtered", qid="0.9", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'filtered_db': 12}
```

## 13) `dnakit.alignment.mmseqs_extractalignedregion()`

- **作用：** 截取比对覆盖的查询或目标区域。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs extractalignedregion`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_extractalignedregion(query_db, target_db, alignment_db, region_db, *, extract_mode=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `alignment_db`, `region_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_extractalignedregion("work/prot", "work/prot", "work/aln", "work/region", extract_mode=1, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'region_db': 36}
```

## 14) `dnakit.alignment.mmseqs_result2msa()`

- **作用：** 由结果数据库生成每条查询的多序列比对（DNA 用 dna=True）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_result2msa(query_db, target_db, result_db, msa_db, *, dna=False, msa_format_mode=None, filter_msa=None, filter_min_enable=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `result_db`, `msa_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_result2msa("work/prot", "work/prot", "work/aln", "work/msa", verbose=False)
al.mmseqs_prefilter("work/nt", "work/nt", "work/ntpf", verbose=False)
al.mmseqs_rescorediagonal("work/nt", "work/nt", "work/ntpf", "work/ntaln", rescore_mode=2, add_backtrace=True, verbose=False)
d = wf.mmseqs_result2msa("work/nt", "work/nt", "work/ntaln", "work/dnamsa", dna=True, verbose=False)
print(r["entry_counts"], d["entry_counts"])
```

- **示例输出：**

```text
{'msa_db': 12} {'msa_db': 4}
```

## 15) `dnakit.alignment.mmseqs_gff2db()`

- **作用：** 按 GFF3 特征提取序列区域。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs gff2db`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_gff2db(gff_files, sequence_db, feature_db, *, gff_type=None, id_offset=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `gff_files`, `sequence_db`, `feature_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_gff2db(["features.gff3"], "work/nt", "work/features", gff_type="CDS", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'feature_db': 8}
```

## 16) `dnakit.alignment.mmseqs_maskbygff()`

- **作用：** 按 GFF3 特征屏蔽序列数据库中的区域。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_maskbygff(gff_file, sequence_db, masked_db, *, gff_type=None, id_offset=None, names_as_ids=True, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `gff_file`, `sequence_db`, `masked_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_maskbygff("features.gff3", "work/nt", "work/masked", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'masked_db': 4}
```

## 17) `dnakit.alignment.mmseqs_createlinindex()`

- **作用：** 为序列数据库建立 k-mer 线性索引。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_createlinindex(sequence_db, *, kmer_size=None, kmer_per_seq=None, min_seq_id=None, search_type=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_createlinindex("work/nt", search_type=3, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'kmer_index_db': 18}
```

## 18) `dnakit.alignment.mmseqs_linsearch()`

- **作用：** 线性时间、较低灵敏度的序列搜索。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_linsearch(query, target, output, *, input_kind='fasta', search_type=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query`, `target`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_linsearch("contigs.fa", "contigs.fa", "work/lin.tsv", search_type=3, verbose=False)
print(sorted({row["query"] for row in r["records"]}))
d = wf.mmseqs_linsearch("work/nt", "work/nt", "work/lin_db", input_kind="database", search_type=3, verbose=False)
print(d["entry_counts"])
```

- **示例输出：**

```text
['contig0', 'contig1', 'contig2', 'contig3']
{'alignment_db': 4}
```

## 19) `dnakit.alignment.mmseqs_linclust()`

- **作用：** 线性时间聚类。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_linclust(sequences, output, *, input_kind='fasta', min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequences`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_linclust("proteins.fa", "work/lc", verbose=False)
print(r["cluster_count"], r["sequence_count"])
```

- **示例输出：**

```text
7 12
```

## 20) `dnakit.alignment.mmseqs_rbh()`

- **作用：** 互为最佳匹配（reciprocal best hit）搜索。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：``）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_rbh(query, target, output, *, input_kind='fasta', sensitivity=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query`, `target`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_rbh("proteins.fa", "proteins.fa", "work/rbh.tsv", verbose=False)
print(r["pair_count"], r["records"][0]["query"])
```

- **示例输出：**

```text
12 fam0_m0
```
