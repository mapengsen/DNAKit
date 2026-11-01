# VSEARCH：搜索、聚类、嵌合体、去冗余、读段过滤与合并

`dnakit.fasta_fastq` 中的这组函数是对 [VSEARCH](https://github.com/torognes/vsearch)（用 2.32.0 测试）的带类型适配器。每个函数以固定的、带类型的选项列表不经 shell 运行一条 VSEARCH 命令；只因参数不同而区分的命令（嵌合体检测方法、去冗余变体、读段过滤、UDB 报告）合并为带该参数的一个函数。它们涵盖全局比对与搜索、UNOISE3 去噪、嵌合体检测、去冗余与再复制、排序、FASTQ 统计、读段过滤、双端合并/连接/同步、记录与子序列提取、定向、限制性酶切、打乱、SFF 转换、SINTAX 分类和 UDB 索引。它们没有暴露 VSEARCH 的全部选项，也不做聚类（聚类请用 `dnakit.fasta_fastq.cluster_vsearch()`）。

VSEARCH 为双许可（BSD-2-Clause 或 GPL-3.0）。DNAKit 把它作为独立进程运行，不复制其代码。

这些函数位于模块 `dnakit.fasta_fastq.vsearch_p09_cli_functions`。VSEARCH 不是 Python 包：用 `conda install -c bioconda vsearch`、`apt install vsearch` 或发布的二进制安装，并放入 PATH 或传入 `vsearch_path`。缺失时函数在创建任何输出目录之前抛出带安装提示的 `BackendUnavailableError`。每个函数都写入新的 `output_directory`（该目录不能已存在），使用文档说明的文件名，校验输入（`ConfigurationError`），VSEARCH 失败时抛出带日志末尾的 `BackendExecutionError`。可选参数默认为 `None`，此时不传给程序，因此沿用 VSEARCH 的默认值。VSEARCH 会拒绝其命令不接受的选项。

## 示例数据

所有示例都从工作目录读取下面这些小文件。`toy.fa` 含两条不同的 90 碱基序列（带丰度注释）、第一条的一个变体、两条序列的一个嵌合体和第一条的完全重复。`site.fa` 有一条含 EcoRI 位点（`GAATTC`）的记录。`tax.fa` 是含两条记录的物种分类数据库。`reads.fq` 有四条末端碱基质量较低的读段，`r1.fq`/`r2.fq` 是三对有 30 碱基重叠的双端读段（`r1_orphan.fq` 多一条正向读段）。每个示例都把输出写入新的临时目录，并假设 `vsearch` 在 PATH 中。

**示例文件：** `toy.fa`

```text
>a;size=50
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
>b;size=30
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
>a2;size=5
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAT
>chim;size=3
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATATTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
>dup1;size=2
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
```

**示例文件：** `site.fa`

```text
>withsite
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGGAATTCCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAA
>nosite
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCA
```

**示例文件：** `tax.fa`

```text
>t1;tax=d:Bacteria,p:Firmicutes,c:Bacilli
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
>t2;tax=d:Bacteria,p:Proteobacteria,c:Gamma
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
```

**示例文件：** `reads.fq`

```text
@read0
TGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((
@read1
GCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((
@read2
AGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((((((((((((
@read3
ACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((((((((((((((((((((((
```

**示例文件：** `r1.fq`

```text
@pair0/1
TTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/1
GGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/1
TCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**示例文件：** `r2.fq`

```text
@pair0/2
AAGGTTCGACAAACTTAGTTCGGGTGGGATTTATGCTCCGTTAGGAATAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/2
AGATTATACCCCTTAGTAACGTTTATGAATGAGCCCGTCTTACCTGGAGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/2
GCGCCCCTAGGCGCACAAGACTGCTGTAATTGGATTGCAAGTCTACTAAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**示例文件：** `r1_orphan.fq`

```text
@pair0/1
TTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/1
GGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/1
TCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@lonely/1
CCCAAAGGTAAACGAACCGTTGCGGTCAATCTTGTCGCGGCTGATGAATT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_allpairs_global()`

- **作用：** 对 FASTA 文件中每一对序列做全局比对，并报告每一对的相似度。
- **处理方法：** 运行 `vsearch --allpairs_global`，输出 `--alnout` 和 `--blast6out`。`accept_all=True` 加 `--acceptall`，报告所有序列对；为 False 时只报告达到 `min_identity`（`--id`）的序列对。`identity_definition` 为 `--iddef`，`match_score`、`mismatch_score`、`gap_open`、`gap_extend` 为打分选项。类 BLAST 表格解析到 `records`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_allpairs_global(input_fasta: 'str | Path', output_directory: 'str | Path', *, accept_all: 'bool' = True, min_identity: 'float | None' = None, identity_definition: 'int | None' = None, match_score: 'int | None' = None, mismatch_score: 'int | None' = None, gap_open: 'str | None' = None, gap_extend: 'str | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 一个 FASTA 文件、一个新的输出目录和上述可选参数。
- **输出：** 字典，含 `status`、`version`、`command`、`outputs`（`alignments` 为 `alignments.aln`，`hits` 为 `hits.tsv`）、`records`（查询、目标、相似度、比对长度、错配数、gap 数、位置、e 值、bit 分数）、`record_count` 和 `log_tail`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_allpairs_global

result = vsearch_allpairs_global("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
best = max(result["records"], key=lambda row: row["identity"])
print(result["status"], result["record_count"])
print(best["query"], best["target"], best["identity"])
```

- **示例输出：**

```text
completed 10
a;size=50 dup1;size=2 100.0
```

## 2) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_search_global()`

- **作用：** 用全局比对在 FASTA 或 UDB 数据库中搜索查询序列。
- **处理方法：** 运行 `vsearch --usearch_global`，使用 `--db`、必需的 `--id`（`min_identity`）并输出类 BLAST 表格。搜索设置对应 `--strand`、`--iddef`、`--maxaccepts`、`--maxrejects`、`--minwordmatches`、`--wordlength`、`--query_cov`、`--target_cov`、`--top_hits_only`、`--output_no_hits`、`--qmask`、`--dbmask` 和 `--sizein`。`write_alignments`、`write_uc`、`write_otu_table`、`write_sam`、`write_matched_queries` 会增加输出文件。VSEARCH 库 API 中的单序列和批量搜索调用就是这个操作。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_search_global(query_fasta: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, min_identity: 'float', strand: "Literal['plus', 'both'] | None" = None, identity_definition: 'int | None' = None, max_accepts: 'int | None' = None, max_rejects: 'int | None' = None, min_word_matches: 'int | None' = None, word_length: 'int | None' = None, query_coverage: 'float | None' = None, target_coverage: 'float | None' = None, top_hits_only: 'bool' = False, output_no_hits: 'bool' = False, query_mask: "Literal['none', 'dust', 'soft'] | None" = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, abundance_in: 'bool' = False, write_alignments: 'bool' = False, write_uc: 'bool' = False, write_otu_table: 'bool' = False, write_sam: 'bool' = False, write_matched_queries: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 查询 FASTA、数据库（FASTA 或 UDB）、新的输出目录，以及 0 到 1 之间的 `min_identity`。
- **输出：** `status`、`outputs`（`hits` 及可选文件）、`records`（类 BLAST 行）、`record_count` 和通用字段。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_search_global

result = vsearch_search_global("toy.fa", "toy.fa", tempfile.mkdtemp() + "/run", min_identity=0.9, abundance_in=True, write_otu_table=True, verbose=False)
for row in sorted(result["records"], key=lambda row: row["query"]):
    print(row["query"], row["target"], row["identity"])
print(sorted(result["outputs"]))
```

- **示例输出：**

```text
a2;size=5 a2;size=5 100.0
a;size=50 a;size=50 100.0
b;size=30 b;size=30 100.0
chim;size=3 chim;size=3 100.0
dup1;size=2 a;size=50 100.0
['hits', 'otu_table']
```

## 3) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cluster_unoise()`

- **作用：** 用 UNOISE3 算法对扩增子序列去噪，每条真实生物序列保留一个中心序列。
- **处理方法：** 运行 `vsearch --cluster_unoise`，输出 `--centroids` 和 `--uc`。`min_abundance` 为 `--minsize`（VSEARCH 默认 8），`unoise_alpha` 为 `--unoise_alpha`（2.0），`identity` 为 `--id`；`abundance_in`/`abundance_out` 处理 `;size=` 注释，`relabel` 重命名中心序列。`write_otu_table` 增加 `--otutabout`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cluster_unoise(input_fasta: 'str | Path', output_directory: 'str | Path', *, min_abundance: 'int | None' = None, unoise_alpha: 'float | None' = None, identity: 'float | None' = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, relabel: 'str | None' = None, write_otu_table: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 标签带 `;size=` 丰度的 FASTA 文件和新的输出目录。
- **输出：** `outputs`（`centroids` 为 `centroids.fasta`，`clusters` 为 `clusters.uc`，可选 `otu_table`）、`centroid_count` 和通用字段。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_cluster_unoise

result = vsearch_cluster_unoise("toy.fa", tempfile.mkdtemp() + "/run", min_abundance=2, abundance_in=True, abundance_out=True, verbose=False)
print(result["status"], result["centroid_count"])
print([line for line in open(result["outputs"]["centroids"]) if line.startswith(">")])
```

- **示例输出：**

```text
completed 3
['>a;size=57\n', '>b;size=30\n', '>chim;size=3\n']
```

## 4) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_detect_chimeras()`

- **作用：** 用 UCHIME、UCHIME2/3、基于参考库的 UCHIME 或新的从头算法检测嵌合序列。
- **处理方法：** `method` 选择 VSEARCH 命令：`uchime_denovo`、`uchime2_denovo`、`uchime3_denovo`、`uchime_ref`（需要 `database`）或 `chimeras_denovo`（长的精确序列，选项 `chimeras_diff_pct`、`chimeras_length_min`、`chimeras_parents_max`、`chimeras_parts`）。`abundance_skew` 为 `--abskew`；`min_score`、`min_divergence`、`min_diffs`、`no_vote_pseudocount`、`no_vote_weight` 分别为 `--minh`、`--mindiv`、`--mindiffs`、`--dn`、`--xn`。属于其他方法的参数会在运行程序前被拒绝。18 列报告解析到 `records`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_detect_chimeras(input_fasta: 'str | Path', output_directory: 'str | Path', *, method: "Literal['uchime_denovo', 'uchime2_denovo', 'uchime3_denovo', 'uchime_ref', 'chimeras_denovo']" = 'uchime_denovo', database: 'str | Path | None' = None, abundance_skew: 'float | None' = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_score: 'float | None' = None, min_divergence: 'float | None' = None, min_diffs: 'int | None' = None, no_vote_pseudocount: 'float | None' = None, no_vote_weight: 'float | None' = None, exclude_same_label: 'bool' = False, exclude_same_sequence: 'bool' = False, chimeras_diff_pct: 'float | None' = None, chimeras_length_min: 'int | None' = None, chimeras_parents_max: 'int | None' = None, chimeras_parts: 'int | None' = None, relabel: 'str | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 带 `;size=` 丰度的 FASTA 文件（`abundance_in=True`）、新的输出目录，`uchime_ref` 还需要参考 FASTA。
- **输出：** `outputs`（`chimeras`、`nonchimeras`、`report`、`alignments`，UCHIME 方法还有 `borderline`）、`records`（分数、查询、亲本、相似度、投票、差异度、`chimeric` 为 Y/N/?）、`chimera_count`、`nonchimera_count` 和 `method`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_detect_chimeras

result = vsearch_detect_chimeras("toy.fa", tempfile.mkdtemp() + "/run", method="uchime_denovo", abundance_in=True, verbose=False)
print(result["chimera_count"], result["nonchimera_count"])
print([row["query"] for row in result["records"] if row["chimeric"] == "Y"])
```

- **示例输出：**

```text
1 4
['chim;size=3']
```

## 5) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cut_restriction_sites()`

- **作用：** 在限制性内切酶位点处切割序列。
- **处理方法：** 运行 `vsearch --cut`，`--cut_pattern` 用 `^` 和 `_` 标出正链和反链的切点（例如 `G^AATT_C`）。正链和反链的片段分别输出；不含该位点的序列写入 discarded 文件。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cut_restriction_sites(input_fasta: 'str | Path', output_directory: 'str | Path', *, cut_pattern: 'str', vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 文件、新的输出目录和 `cut_pattern`（核苷酸字母，必须同时含 `^` 和 `_`）。
- **输出：** `outputs`（`fragments_forward`、`fragments_reverse`、`discarded_forward`、`discarded_reverse`）、`forward_fragment_count`、`reverse_fragment_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_cut_restriction_sites

result = vsearch_cut_restriction_sites("site.fa", tempfile.mkdtemp() + "/run", cut_pattern="G^AATT_C", verbose=False)
print(result["forward_fragment_count"], result["reverse_fragment_count"])
print([line.strip() for line in open(result["outputs"]["fragments_forward"]) if line.startswith(">")])
```

- **示例输出：**

```text
2 2
['>withsite', '>withsite']
```

## 6) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_fulllength()`

- **作用：** 合并完全相同的序列并计数（全长去冗余）。
- **处理方法：** `method='fulllength'` 运行 `--derep_fulllength`，`'smallmem'` 运行 `--derep_smallmem`（省内存，无 UC 文件），`'fastx_uniques'` 对 FASTA 或 FASTQ 输入运行 `--fastx_uniques`（`fastq_output=True` 保留质量值）。`strand`、`abundance_in`、`abundance_out`、`min_unique_size`、`max_unique_size`、`top_n`、`relabel`、`relabel_mode` 对应 `--strand`、`--sizein`、`--sizeout`、`--minuniquesize`、`--maxuniquesize`、`--topn`、`--relabel` 和摘要重命名开关。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_fulllength(input_fasta: 'str | Path', output_directory: 'str | Path', *, method: "Literal['fulllength', 'smallmem', 'fastx_uniques']" = 'fulllength', fastq_output: 'bool' = False, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA（`fastx_uniques` 可为 FASTQ）文件和新的输出目录。
- **输出：** `outputs`（`unique`，除 `smallmem` 外还有 `clusters`）和 `unique_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_fulllength

result = vsearch_dereplicate_fulllength("toy.fa", tempfile.mkdtemp() + "/run", abundance_in=True, abundance_out=True, verbose=False)
print(result["unique_count"])
print([line.strip() for line in open(result["outputs"]["unique"]) if line.startswith(">")])
```

- **示例输出：**

```text
4
['>a;size=52', '>b;size=30', '>a2;size=5', '>chim;size=3']
```

## 7) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_by_id()`

- **作用：** 按标识符和序列同时去冗余。
- **处理方法：** 运行 `vsearch --derep_id`。只合并标签和序列都相同的记录，序列相同但标签不同的记录保持分开。选项与上一个函数的去冗余选项相同。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_by_id(input_fasta: 'str | Path', output_directory: 'str | Path', *, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 文件和新的输出目录。
- **输出：** `outputs`（`unique`、`clusters`）和 `unique_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_by_id

result = vsearch_dereplicate_by_id("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print(result["unique_count"])
```

- **示例输出：**

```text
5
```

## 8) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_prefix()`

- **作用：** 把作为更长序列前缀的序列合并到长序列中。
- **处理方法：** 运行 `vsearch --derep_prefix`。与更长序列开头完全相同的序列计入长序列。`threads` 作为 `--threads` 传入。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_prefix(input_fasta: 'str | Path', output_directory: 'str | Path', *, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 文件和新的输出目录。
- **输出：** `outputs`（`unique`、`clusters`）和 `unique_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_prefix

prefix = tempfile.mkdtemp() + "/prefixes.fa"
lines = open("toy.fa").read().splitlines()
open(prefix, "w").write(">long\n" + lines[1] + "\n>short\n" + lines[1][:60] + "\n>other\n" + lines[3] + "\n")
result = vsearch_dereplicate_prefix(prefix, tempfile.mkdtemp() + "/run", abundance_out=True, verbose=False)
print(result["unique_count"])
print([line.strip() for line in open(result["outputs"]["unique"]) if line.startswith(">")])
```

- **示例输出：**

```text
2
['>long;size=2', '>other;size=1']
```

## 9) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_rereplicate()`

- **作用：** 按每条记录 `;size=` 丰度指示的次数展开成相应数量的副本。
- **处理方法：** 运行 `vsearch --rereplicate`。`strip_abundance=True` 加 `--xsize`，从标签中去掉丰度注释。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_rereplicate(input_fasta: 'str | Path', output_directory: 'str | Path', *, strip_abundance: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 带 `;size=` 注释的 FASTA 文件和新的输出目录。
- **输出：** `outputs`（`rereplicated`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_rereplicate

result = vsearch_rereplicate("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
90
```

## 10) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sort_by_size()`

- **作用：** 按丰度从高到低对 FASTA 记录排序。
- **处理方法：** 运行 `vsearch --sortbysize`。`min_size` 和 `max_size` 去掉丰度范围之外的记录，`top_n` 保留前 n 条，`abundance_out` 写出 `;size=`，`relabel`/`relabel_mode` 重命名记录。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sort_by_size(input_fasta: 'str | Path', output_directory: 'str | Path', *, min_size: 'int | None' = None, max_size: 'int | None' = None, top_n: 'int | None' = None, abundance_out: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 带 `;size=` 注释的 FASTA 文件和新的输出目录。
- **输出：** `outputs`（`sorted`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sort_by_size

result = vsearch_sort_by_size("toy.fa", tempfile.mkdtemp() + "/run", min_size=3, verbose=False)
print(result["record_count"])
print([line.strip() for line in open(result["outputs"]["sorted"]) if line.startswith(">")])
```

- **示例输出：**

```text
4
['>a;size=50', '>b;size=30', '>a2;size=5', '>chim;size=3']
```

## 11) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats()`

- **作用：** FASTQ 文件按位置统计的质量值和期望错误数。
- **处理方法：** 运行 `vsearch --fastq_eestats` 并解析表格。`fastq_ascii`、`fastq_qmin`、`fastq_qmax` 设置质量值编码。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats(input_fastq: 'str | Path', output_directory: 'str | Path', *, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTQ 文件和新的输出目录。
- **输出：** `outputs`（`statistics` 为 `eestats.tsv`）、`records`（每个位置一行，含 `pos`、`recs`，以及 Q、Pe、EE 的最小/低/中位/平均/高/最大值）和 `position_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_eestats

result = vsearch_fastq_eestats("reads.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["position_count"])
print(result["records"][0]["pos"], result["records"][0]["recs"], result["records"][-1]["recs"])
```

- **示例输出：**

```text
70
1 4 1
```

## 12) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats2()`

- **作用：** 在各个长度阈值和最大期望错误数下能保留多少条读段。
- **处理方法：** 运行 `vsearch --fastq_eestats2`。`ee_cutoffs` 为 MaxEE 值（默认 0.5、1、2），`length_cutoffs` 为 `min,max,increment`（默认 `50,*,50`）。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats2(input_fastq: 'str | Path', output_directory: 'str | Path', *, ee_cutoffs: 'Sequence[float] | None' = None, length_cutoffs: 'str | None' = None, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTQ 文件、新的输出目录和可选阈值。
- **输出：** `outputs`（`statistics`）、`records`（每个长度阈值一项，含 `retained` 列表，其中有 `max_ee`、`reads`、`percent`）以及摘要 `reads`、`max_length`、`mean_length`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_eestats2

result = vsearch_fastq_eestats2("reads.fq", tempfile.mkdtemp() + "/run", ee_cutoffs=[1, 3], length_cutoffs="30,70,20", verbose=False)
print(result["reads"], result["max_length"])
for row in result["records"]:
    print(row["length"], [item["reads"] for item in row["retained"]])
```

- **示例输出：**

```text
4 70
30 [4, 4]
50 [0, 0]
70 [0, 0]
```

## 13) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_stats()`

- **作用：** FASTQ 文件的详细统计：读段长度、质量值和截短情况。
- **处理方法：** 运行 `vsearch --fastq_stats`，报告通过 `--log` 写入 `stats.txt`，并从中解析读长和质量值分布表。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_stats(input_fastq: 'str | Path', output_directory: 'str | Path', *, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTQ 文件和新的输出目录。
- **输出：** `outputs`（`report`）、`read_length_distribution`、`quality_distribution`、`read_count`、`mean_length` 和 `bases`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_stats

result = vsearch_fastq_stats("reads.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["read_count"], result["mean_length"])
print(result["read_length_distribution"][0])
```

- **示例输出：**

```text
4 55.0
{'length': 70, 'reads': 1, 'percent': 25.0, 'cumulative_percent': 25.0}
```

## 14) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_filter_reads()`

- **作用：** 按质量值、期望错误数、长度、N 含量和丰度对 FASTA 或 FASTQ 记录做截短和过滤。
- **处理方法：** 运行 `vsearch --fastx_filter`。一个函数涵盖这些条件：`min_quality`（`--fastq_minqual`）、`truncate_quality`（`--fastq_truncqual`）、`truncate_ee` 和 `truncate_ee_rate`（`--fastq_truncee`、`--fastq_truncee_rate`）、`max_ee`、`max_ee_rate`、`min_length`、`max_length`、`max_ns`、`strip_left`、`strip_right`、`truncate_length`、`truncate_length_keep`，以及丰度过滤 `min_size` 和 `max_size`。`reverse_path` 对双端读段一起过滤，`write_discarded` 保留被拒绝的记录。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_filter_reads(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, reverse_path: 'str | Path | None' = None, min_quality: 'int | None' = None, truncate_quality: 'int | None' = None, truncate_ee: 'float | None' = None, truncate_ee_rate: 'float | None' = None, max_ee: 'float | None' = None, max_ee_rate: 'float | None' = None, min_length: 'int | None' = None, max_length: 'int | None' = None, max_ns: 'int | None' = None, strip_left: 'int | None' = None, strip_right: 'int | None' = None, truncate_length: 'int | None' = None, truncate_length_keep: 'int | None' = None, min_size: 'int | None' = None, max_size: 'int | None' = None, write_discarded: 'bool' = False, include_expected_errors: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, abundance_out: 'bool' = False, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTQ 或 FASTA 文件、新的输出目录和过滤条件。除非给出 `output_format`，输出格式与输入一致。
- **输出：** `outputs`（`filtered`，可选 `filtered_reverse`、`discarded`、`discarded_reverse`）、`passed_count`，使用 `write_discarded` 时还有 `discarded_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_filter_reads

result = vsearch_filter_reads("reads.fq", tempfile.mkdtemp() + "/run", truncate_quality=20, min_length=30, write_discarded=True, verbose=False)
print(result["passed_count"], result["discarded_count"])
lengths = [len(line.strip()) for line in open(result["outputs"]["filtered"]).readlines()[1::4]]
print(lengths)
```

- **示例输出：**

```text
4 0
[30, 30, 30, 30]
```

## 15) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sync_pairs()`

- **作用：** 让双端读段的两个文件顺序一致，并把没有配对的读段单独存放。
- **处理方法：** 运行 `vsearch --fastx_syncpairs`，输出已同步读段和孤儿读段共四个文件。`read_separators` 为配对编号前的字符（默认 `/`）。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sync_pairs(forward_path: 'str | Path', reverse_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, read_separators: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 正向和反向 FASTQ（或 FASTA）文件和新的输出目录。
- **输出：** `outputs`（`synced_forward`、`synced_reverse`、`orphans_forward`、`orphans_reverse`）、`synced_pairs`、`forward_orphans`、`reverse_orphans`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sync_pairs

result = vsearch_sync_pairs("r1_orphan.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["synced_pairs"], result["forward_orphans"], result["reverse_orphans"])
```

- **示例输出：**

```text
3 1 0
```

## 16) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_merge_pairs()`

- **作用：** 把有重叠的双端读段合并成一条序列。
- **处理方法：** 运行 `vsearch --fastq_mergepairs`。参数有 `allow_stagger`、`max_diffs`、`max_diff_pct`、`max_ee`、`max_ns`、`min_length`、`min_overlap`、`min_merge_length`、`max_merge_length`、`truncate_quality`、`label_suffix`、`include_expected_errors`、质量值范围选项和 `write_error_table`（`--eetabbedout`）。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_merge_pairs(forward_fastq: 'str | Path', reverse_fastq: 'str | Path', output_directory: 'str | Path', *, allow_stagger: 'bool' = False, max_diffs: 'int | None' = None, max_diff_pct: 'float | None' = None, max_ee: 'float | None' = None, max_ns: 'int | None' = None, min_length: 'int | None' = None, min_overlap: 'int | None' = None, min_merge_length: 'int | None' = None, max_merge_length: 'int | None' = None, truncate_quality: 'int | None' = None, label_suffix: 'str | None' = None, include_expected_errors: 'bool' = False, write_error_table: 'bool' = False, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, fastq_qmin_out: 'int | None' = None, fastq_qmax_out: 'int | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 正向和反向 FASTQ 文件和新的输出目录。
- **输出：** `outputs`（`merged`、`notmerged_forward`、`notmerged_reverse`，可选 `errors`）、`merged_count`、`not_merged_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_merge_pairs

result = vsearch_merge_pairs("r1.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["merged_count"], result["not_merged_count"])
print([len(line.strip()) for line in open(result["outputs"]["merged"]).readlines()[1::4]])
```

- **示例输出：**

```text
3 0
[70, 70, 70]
```

## 17) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_join_pairs()`

- **作用：** 把双端读段用一段间隔序列连接成一条序列。
- **处理方法：** 运行 `vsearch --fastq_join`。反向读段取反向互补后接在 `pad_sequence`（`--join_padgap`，默认 `NNNNNNNN`）之后，间隔的质量值为 `pad_quality`（`--join_padgapq`，默认 `IIIIIIII`）。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_join_pairs(forward_fastq: 'str | Path', reverse_fastq: 'str | Path', output_directory: 'str | Path', *, pad_sequence: 'str | None' = None, pad_quality: 'str | None' = None, write_fasta: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 正向和反向 FASTQ 文件和新的输出目录。
- **输出：** `outputs`（`joined`，可选 `joined_fasta`）和 `joined_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_join_pairs

result = vsearch_join_pairs("r1.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
sequence = open(result["outputs"]["joined"]).readlines()[1].strip()
print(result["joined_count"], len(sequence), sequence[50:58])
```

- **示例输出：**

```text
3 108 NNNNNNNN
```

## 18) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_records()`

- **作用：** 按标题标签从 FASTA 或 FASTQ 文件中提取记录。
- **处理方法：** 运行 `vsearch --fastx_getseqs`；只给一个 `label` 即单条记录的 `--fastx_getseq` 情形。选择器有 `label`、`labels_file`、`label_word`、`label_words_file`（配合 `label_field`）；`substring_match` 允许标题任意位置匹配；`write_unmatched` 同时输出其余记录。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_records(input_path: 'str | Path', output_directory: 'str | Path', *, label: 'str | None' = None, labels_file: 'str | Path | None' = None, label_word: 'str | None' = None, label_words_file: 'str | Path | None' = None, label_field: 'str | None' = None, substring_match: 'bool' = False, output_format: "Literal['fasta', 'fastq'] | None" = None, write_unmatched: 'bool' = False, abundance_in: 'bool' = False, abundance_out: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 或 FASTQ 文件、新的输出目录和至少一个选择器。
- **输出：** `outputs`（`extracted`，可选 `unmatched`）和 `extracted_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_extract_records

result = vsearch_extract_records("toy.fa", tempfile.mkdtemp() + "/run", label="b;size=30", write_unmatched=True, verbose=False)
print(result["extracted_count"])
print([line.strip() for line in open(result["outputs"]["unmatched"]) if line.startswith(">")])
```

- **示例输出：**

```text
1
['>a;size=50', '>a2;size=5', '>chim;size=3', '>dup1;size=2']
```

## 19) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_subsequence()`

- **作用：** 按标签和位置提取单条记录的子序列。
- **处理方法：** 运行 `vsearch --fastx_getsubseq`，使用 `--label`、`--subseq_start` 和 `--subseq_end`（从 1 开始，含两端）。`label_suffix` 追加到标签后。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_subsequence(input_path: 'str | Path', output_directory: 'str | Path', *, label: 'str', start: 'int', end: 'int', substring_match: 'bool' = False, output_format: "Literal['fasta', 'fastq'] | None" = None, label_suffix: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 或 FASTQ 文件、新的输出目录、`label`、`start` 和 `end`。
- **输出：** `outputs`（`subsequence`）和 `extracted_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_extract_subsequence

result = vsearch_extract_subsequence("toy.fa", tempfile.mkdtemp() + "/run", label="a;size=50", start=1, end=20, verbose=False)
print(open(result["outputs"]["subsequence"]).read().split()[1])
```

- **示例输出：**

```text
TTTCCTCATGCAATTCAAAA
```

## 20) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_orient_sequences()`

- **作用：** 按参考数据库的链方向对序列定向。
- **处理方法：** 运行 `vsearch --orient`，以参考序列作为 `--db`；与参考反向匹配的序列会取反向互补。`word_length`、`query_mask`、`database_mask` 对应 `--wordlength`、`--qmask`、`--dbmask`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_orient_sequences(input_path: 'str | Path', reference_fasta: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, word_length: 'int | None' = None, query_mask: "Literal['none', 'dust', 'soft'] | None" = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 或 FASTQ 文件、参考 FASTA 和新的输出目录。
- **输出：** `outputs`（`oriented`、`undetermined`、`orientation`）、`records`（查询、链 `+`/`-`/`?`、正向和反向匹配数）、`record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_orient_sequences

comp = str.maketrans("ACGT", "TGCA")
forward = open("toy.fa").read().splitlines()[1]
flipped = tempfile.mkdtemp() + "/flipped.fa"
open(flipped, "w").write(">flipped\n" + forward.translate(comp)[::-1] + "\n>straight\n" + forward + "\n")
result = vsearch_orient_sequences(flipped, "toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print([(row["query"], row["strand"]) for row in result["records"]])
```

- **示例输出：**

```text
[('flipped', '-'), ('straight', '+')]
```

## 21) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_scramble_sequences()`

- **作用：** 随机打乱每条序列内部的核苷酸顺序。
- **处理方法：** 运行 `vsearch --scramble`。`kmer_size`（`--scramble_kmer`）保持不超过该长度的词频，`random_seed`（`--randseed`）使结果可重复。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_scramble_sequences(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, kmer_size: 'int | None' = None, random_seed: 'int | None' = None, label_suffix: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 或 FASTQ 文件和新的输出目录。
- **输出：** `outputs`（`scrambled`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_scramble_sequences

result = vsearch_scramble_sequences("toy.fa", tempfile.mkdtemp() + "/run", random_seed=3, verbose=False)
original = open("toy.fa").read().splitlines()[1]
scrambled = "".join(open(result["outputs"]["scrambled"]).read().split(">")[1].split("\n")[1:])
print(result["record_count"], sorted(original) == sorted(scrambled), original != scrambled)
```

- **示例输出：**

```text
5 True True
```

## 22) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_convert_sff()`

- **作用：** 把 SFF 流图文件（454、Ion Torrent）转换为 FASTQ。
- **处理方法：** 运行 `vsearch --sff_convert`。`clip_ends=True` 应用文件中保存的裁剪位置（`--sff_clip`）；`fastq_ascii_out`、`fastq_qmin_out`、`fastq_qmax_out` 设置输出质量值编码。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_convert_sff(input_sff: 'str | Path', output_directory: 'str | Path', *, clip_ends: 'bool' = False, fastq_ascii_out: 'int | None' = None, fastq_qmin_out: 'int | None' = None, fastq_qmax_out: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** SFF 文件和新的输出目录。
- **输出：** `outputs`（`converted`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_convert_sff

import struct

flows, key = "TACG" * 3, "TCAG"
head_len = 31 + len(flows) + len(key)
head_len += -head_len % 8
blob = struct.pack(">IIQIIHHHB", 0x2E736666, 1, 0, 0, 1, head_len, len(key), len(flows), 1) + flows.encode() + key.encode()
blob += bytes(head_len - len(blob))
bases = "ACGTACGTAC"
read_head = struct.pack(">HHIHHHH", 24, 2, len(bases), 1, len(bases), 0, 0) + b"r1"
read_head += bytes(24 - len(read_head))
data = struct.pack(">12H", *[100] * 12) + bytes([1] * 10) + bases.encode() + bytes([30] * 10)
data += bytes(-len(data) % 8)
sff = tempfile.mkdtemp() + "/toy.sff"
open(sff, "wb").write(blob + read_head + data)
result = vsearch_convert_sff(sff, tempfile.mkdtemp() + "/run", verbose=False)
print(open(result["outputs"]["converted"]).read().split())
```

- **示例输出：**

```text
['@r1', 'ACGTACGTAC', '+', '??????????']
```

## 23) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sintax()`

- **作用：** 用 SINTAX 算法对序列做物种分类。
- **处理方法：** 运行 `vsearch --sintax`，数据库为标签带 `;tax=d:...,p:...` 注释的 FASTA 或 UDB。`cutoff` 为 `--sintax_cutoff`，`random_seed` 为 `--randseed`，`random_ties` 为 `--sintax_random`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sintax(input_path: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, cutoff: 'float | None' = None, random_seed: 'int | None' = None, random_ties: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 查询 FASTA 或 FASTQ、物种分类数据库和新的输出目录。
- **输出：** `outputs`（`classification`）、`records`（查询、各层级带置信度的命中、链、超过阈值的分类）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sintax

result = vsearch_sintax("toy.fa", "tax.fa", tempfile.mkdtemp() + "/run", cutoff=0.5, random_seed=1, verbose=False)
for row in result["records"][:2]:
    print(row["query"], row["taxonomy"])
```

- **示例输出：**

```text
a;size=50 d:Bacteria,p:Firmicutes,c:Bacilli
b;size=30 d:Bacteria,p:Proteobacteria,c:Gamma
```

## 24) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_make_udb()`

- **作用：** 由 FASTA 文件构建 UDB 词索引数据库。
- **处理方法：** 运行 `vsearch --makeudb_usearch`。`word_length`（`--wordlength`，默认 8）和 `database_mask`（`--dbmask`，默认 dust）保存在索引中；`hardmask` 用 N 屏蔽。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_make_udb(input_fasta: 'str | Path', output_directory: 'str | Path', *, word_length: 'int | None' = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, hardmask: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** FASTA 文件和新的输出目录。
- **输出：** `outputs`（`database` 为 `database.udb`）和 `database_bytes`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb

result = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/run", word_length=6, verbose=False)
print(result["status"], result["database_bytes"] > 0)
```

- **示例输出：**

```text
completed True
```

## 25) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_to_fasta()`

- **作用：** 把 UDB 文件中的序列导出为 FASTA。
- **处理方法：** 运行 `vsearch --udb2fasta`。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_to_fasta(input_udb: 'str | Path', output_directory: 'str | Path', *, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** UDB 文件和新的输出目录。
- **输出：** `outputs`（`sequences`）和 `record_count`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb, vsearch_udb_to_fasta

built = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/udb", verbose=False)
result = vsearch_udb_to_fasta(built["outputs"]["database"], tempfile.mkdtemp() + "/run", verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
5
```

## 26) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_report()`

- **作用：** 描述 UDB 文件：摘要信息或 k-mer 词统计。
- **处理方法：** `report='info'` 运行 `vsearch --udbinfo`（序列数、字母表、词宽、索引设置，解析到 `fields`）；`report='stats'` 运行 `vsearch --udbstats`（词大小分布和最常见的词）。报告通过 `--log` 读取。
- **API：** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_report(input_udb: 'str | Path', output_directory: 'str | Path', *, report: "Literal['info', 'stats']" = 'info', vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** UDB 文件、新的输出目录和 `report`。
- **输出：** `outputs`（`report`）、`report_text`，`info` 时还有 `fields`。
- **示例代码：**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb, vsearch_udb_report

built = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/udb", word_length=6, verbose=False)
info = vsearch_udb_report(built["outputs"]["database"], tempfile.mkdtemp() + "/info", verbose=False)
print(info["fields"]["Seqs"], info["fields"]["Word width"])
stats = vsearch_udb_report(built["outputs"]["database"], tempfile.mkdtemp() + "/stats", report="stats", verbose=False)
print(stats["report_text"].splitlines()[0])
```

- **示例输出：**

```text
5 6
450 nt in 5 seqs, min 90, max 90, avg 90
```
