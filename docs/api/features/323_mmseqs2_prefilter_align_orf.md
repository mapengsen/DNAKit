# MMseqs2 预过滤、比对模块与 ORF 工具

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

**示例文件：** `cds.fa`

```text
>fam0_m0
ATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAA
>fam0_m1
ATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTT
>fam0_m2
ATGATGTTTCCTTGTCAAGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTCATGTTTTTCTGCAATTTGTTGAATGGCTGGTTGAAGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTAAACAA
>fam1_m0
ATGGGTAATATTGTTGTTACTATGATTTATCATATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATCGTAATAATGATATTGAAATTTCTCATATGCATTCTTATTATGCTTCTAATGATGAACCTCATTCTGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTGCT
>fam1_m1
ATGGGTAATATTGTTGTTACTATGATTTATGAAATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATCGTAATAATGATATTGAAATTTCTCGTATGCATTCTTATTATGCTTCTACTGATGAACCTCATCATGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTGCT
>fam1_m2
ATGGGTAATATTGTTGTTACTATGATTTATCATATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATTCTAATAATGATATTGAAATTTCTCATAATGTTTCTTATTATGCTTCTAATGATGAACCTCATAAAGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTTGT
>fam2_m0
ATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTG
>fam2_m1
ATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTG
>fam2_m2
ATGGGTTGGGATAAAGCTGATTGTGATTATATTGATAAAGAACGTGCTATGGATCAAAAATATTTTTCTACTAAAGAAGGTAAATGTGGTCATCTGCTGACTCATTTTCGTACTGGTAAAAATGCTAAATGTCATGCTACTGTTCATACTTCTATTCGTGAACAACGTACTCCTACTCTG
>fam3_m0
ATGCATCATGATTGGGATTTTACTAAAAATTTTTATACTAAAGAAAATATTTCTTCTCCTGCTGGTGCTTCTCGTCCTCTGTTTCAAAATCCTATGGAAATGGCTATGATGCCTGAACATGCTCTGAAAAATGATCCTCCTTGGGATAATCAAAAATGTAAAGAATGTCTGTTTATTAAA
>fam3_m1
ATGCATCATGATTGGCCTTTTACTAAAAATCAATATACTAAAGAAGTTATTTCTTCTCCTGCTGGTGCTTCTCGTCCTCTGTTTCAAAATCCTATGGAAATGGCTATGATGCCTATGCATGCTCTGAAAAATCAACCTCCTTGGGATAATCAAAAATGTAAAGAATGTCTGATTATTAAA
>fam3_m2
ATGCATCATGATTGGGATTTTACTAAAAATACTTATACTAAAGAAAATATTTCTTCTCCTGCTGGTGCTTCTTTTCCTCTGTTTCAAAATCCTATGGAAATGACTATGATGCCTGAAGTTGCTCTGAAAAATGATCCTCCTTGGGATAATCAATCTTGTAAAATTTGTCTGTTTATGAAA
```

## 1) `dnakit.alignment.mmseqs_prefilter()`

- **作用：** k-mer 双连续对角线预过滤。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs prefilter`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_prefilter(query_db, target_db, prefilter_db, *, sensitivity=None, kmer_size=None, max_seqs=None, min_ungapped_score=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `prefilter_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
wf.mmseqs_createdb("proteins.fa", "work/prot", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/nt", dbtype=2, verbose=False)
wf.mmseqs_createdb("proteins.fa", "work/prot0", shuffle=False, verbose=False)
wf.mmseqs_createdb("cds.fa", "work/cds", shuffle=False, verbose=False)
r = al.mmseqs_prefilter("work/prot", "work/prot", "work/pf", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'prefilter_db': 12}
```

## 2) `dnakit.alignment.mmseqs_ungappedprefilter()`

- **作用：** 最优无空位对角线预过滤。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs ungappedprefilter`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_ungappedprefilter(query_db, target_db, prefilter_db, *, max_seqs=None, min_ungapped_score=None, comp_bias_corr=None, coverage=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `prefilter_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_ungappedprefilter("work/prot", "work/prot", "work/upf", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'prefilter_db': 12}
```

## 3) `dnakit.alignment.mmseqs_gappedprefilter()`

- **作用：** Smith-Waterman 有空位预过滤（慢）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs gappedprefilter`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_gappedprefilter(query_db, target_db, prefilter_db, *, max_seqs=None, min_ungapped_score=None, comp_bias_corr=None, coverage=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `prefilter_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_gappedprefilter("work/prot", "work/prot", "work/gpf", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'prefilter_db': 12}
```

## 4) `dnakit.alignment.mmseqs_kmermatcher()`

- **作用：** 在同一数据库内查找 bottom-m 哈希 k-mer 匹配。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs kmermatcher`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_kmermatcher(sequence_db, prefilter_db, *, kmer_size=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `prefilter_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_kmermatcher("work/prot", "work/kmm", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'prefilter_db': 12}
```

## 5) `dnakit.alignment.mmseqs_kmersearch()`

- **作用：** 在查询库与目标 k-mer 索引之间查找 k-mer 匹配。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs kmersearch`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_kmersearch(query_db, kmer_index_db, prefilter_db, *, coverage=None, cov_mode=None, kmer_per_seq=None, pick_n_sim_kmer=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `kmer_index_db`, `prefilter_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_createlinindex("work/nt", search_type=3, verbose=False)
r = al.mmseqs_kmersearch("work/nt", "work/nt.linidx", "work/kms", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'prefilter_db': 4}
```

## 6) `dnakit.alignment.mmseqs_rescorediagonal()`

- **作用：** 按对角线计算序列一致性（汉明距离或局部/全局比对）。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs rescorediagonal`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_rescorediagonal(query_db, target_db, prefilter_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `prefilter_db`, `result_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_rescorediagonal("work/prot", "work/prot", "work/pf", "work/aln", rescore_mode=2, add_backtrace=True, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'result_db': 12}
```

## 7) `dnakit.alignment.mmseqs_alignall()`

- **作用：** 对每个结果条目内部的成员做全对全局部比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs alignall`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_alignall(sequence_db, result_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `result_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_clust("work/prot", "work/aln", "work/clu", verbose=False)
r = al.mmseqs_alignall("work/prot", "work/clu", "work/all", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 4}
```

## 8) `dnakit.alignment.mmseqs_transitivealign()`

- **作用：** 通过传递性转移比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs transitivealign`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_transitivealign(sequence_db, alignment_db, transitive_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `alignment_db`, `transitive_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_transitivealign("work/prot", "work/aln", "work/trans", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'transitive_db': 12}
```

## 9) `dnakit.alignment.mmseqs_fwbw()`

- **作用：** 前向-后向（最大准确度）比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs fwbw`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_fwbw(query_db, target_db, alignment_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `alignment_db`, `result_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_fwbw("work/prot", "work/prot", "work/aln", "work/fwbw", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'result_db': 12}
```

## 10) `dnakit.alignment.mmseqs_alignbykmer()`

- **作用：** 以 k-mer 为种子的启发式有空位局部比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs alignbykmer`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_alignbykmer(query_db, target_db, result_db, alignment_db, *, kmer_size=None, spaced_kmer_mode=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `target_db`, `result_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_alignbykmer("work/prot", "work/prot", "work/pf", "work/abk", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 12}
```

## 11) `dnakit.alignment.mmseqs_extractorfs()`

- **作用：** 从核苷酸序列提取六框开放阅读框。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs extractorfs`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_extractorfs(sequence_db, orf_db, *, min_length=None, max_length=None, max_gaps=None, contig_start_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `orf_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_extractorfs("work/nt", "work/orfs", translate=True, min_length=30, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'orf_db': 30}
```

## 12) `dnakit.alignment.mmseqs_extractframes()`

- **作用：** 提取核苷酸序列的阅读框。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs extractframes`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_extractframes(sequence_db, frame_db, *, forward_frames=None, reverse_frames=None, translation_table=None, translate=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `frame_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_extractframes("work/nt", "work/frames", forward_frames="1,2", reverse_frames="1", translate=True, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'frame_db': 12}
```

## 13) `dnakit.alignment.mmseqs_orftocontig()`

- **作用：** 以比对格式写出 ORF 在重叠群上的位置。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs orftocontig`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_orftocontig(contig_db, orf_db, alignment_db, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `contig_db`, `orf_db`, `alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_orftocontig("work/nt", "work/orfs", "work/orf2contig", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'alignment_db': 30}
```

## 14) `dnakit.alignment.mmseqs_recoverlongestorf()`

- **作用：** 为分类学恢复每个重叠群的最长 ORF。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs recoverlongestorf`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_recoverlongestorf(orf_db, result_db, output, *, , threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `orf_db`, `result_db`, `output`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

al.mmseqs_prefilter("work/orfs", "work/prot", "work/orfpf", verbose=False)
al.mmseqs_rescorediagonal("work/orfs", "work/prot", "work/orfpf", "work/orfaln", rescore_mode=2, verbose=False)
r = al.mmseqs_recoverlongestorf("work/orfs", "work/orfaln", "work/longest.tsv", verbose=False)
print(r["line_counts"])
```

- **示例输出：**

```text
{'output': 0}
```

## 15) `dnakit.alignment.mmseqs_splitsequence()`

- **作用：** 把长序列切成带重叠的片段。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs splitsequence`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_splitsequence(sequence_db, split_db, *, sequence_overlap=None, max_seq_len=None, sequence_split_mode=None, headers_split_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `split_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_splitsequence("work/nt", "work/split", max_seq_len=300, sequence_overlap=30, verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'split_db': 8}
```

## 16) `dnakit.alignment.mmseqs_masksequence()`

- **作用：** 用 tantan 软屏蔽低复杂度区域。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs masksequence`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_masksequence(sequence_db, masked_db, *, mask=None, mask_prob=None, mask_lower_case=None, mask_n_repeat=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `sequence_db`, `masked_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_masksequence("work/prot", "work/masked", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'masked_db': 12}
```

## 17) `dnakit.alignment.mmseqs_offsetalignment()`

- **作用：** 按 ORF 起点偏移比对坐标。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs offsetalignment`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_offsetalignment(query_db, query_orf_db, target_db, target_orf_db, alignment_db, offset_db, *, search_type=None, chain_alignments=None, merge_query=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `query_db`, `query_orf_db`, `target_db`, `target_orf_db`, `alignment_db`, `offset_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_offsetalignment("work/nt", "work/orfs", "work/prot", "work/prot", "work/orfaln", "work/offset", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'offset_db': 4}
```

## 18) `dnakit.alignment.mmseqs_proteinaln2nucl()`

- **作用：** 把蛋白比对映射为核苷酸比对。
- **处理方法：** 以固定参数列表运行 MMseqs2 对应模块（见下方英文版的 MMseqs2 命令名：`mmseqs proteinaln2nucl`）；未设置的选项沿用 MMseqs2 默认值，输出前缀/文件必须不存在。
- **API：** `dnakit.alignment.mmseqs_proteinaln2nucl(nucleotide_query_db, nucleotide_target_db, protein_query_db, protein_target_db, alignment_db, nucleotide_alignment_db, *, gap_open=None, gap_extend=None, threads=1, mmseqs_path=None, verbose=True)`。
- **输入：** 位置参数 `nucleotide_query_db`, `nucleotide_target_db`, `protein_query_db`, `protein_target_db`, `alignment_db`, `nucleotide_alignment_db`（数据库为路径前缀，FASTA/GFF3 为文件路径）；关键字参数为可选的 MMseqs2 选项，取值范围见函数文档字符串。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（新文件或数据库前缀）、`entry_counts` 或 `line_counts`、`log_tail` 与 `provenance`；个别函数另有 `records`、`rows`、`clusters` 等字段。
- **示例代码：**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

al.mmseqs_prefilter("work/prot0", "work/prot0", "work/pf0", verbose=False)
al.mmseqs_rescorediagonal("work/prot0", "work/prot0", "work/pf0", "work/aln0", rescore_mode=2, add_backtrace=True, verbose=False)
r = al.mmseqs_proteinaln2nucl("work/cds", "work/cds", "work/prot0", "work/prot0", "work/aln0", "work/nucl_aln", verbose=False)
print(r["entry_counts"])
```

- **示例输出：**

```text
{'nucleotide_alignment_db': 12}
```
