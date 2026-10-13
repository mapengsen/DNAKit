# G18f 功能：polars-bio 与 gtars

本页记录 G18f 组的包装函数。8 个 polars-bio 函数和 6 个 gtars 函数都从 `dnakit.annotation` 导出。例如 `from dnakit.annotation import lola_check_universe`。

- 本地函数：全部 14 个包装函数只读写本地文件，不需要网络，测试也不访问网络。
- 示例输入文件以 `示例文件` 标签给出。二进制输入（CRAM、BAM、bigBed、GTOK）由示例代码生成，因此每个示例都可直接运行。
- mygene 与 myvariant 在本组没有 DNA 级函数，两个模块都不导出函数，见下方未包装部分。

## 已有包装与未包装的函数

- polars-bio：68 个 `other` 行。8 个已包装（本页的函数）。35 个是上述函数的族变体，如惰性的 `scan_*`、`sink_*`，`write_bam` / `write_cram` / `write_sam`，以及 `query_genomic_tables` 的 16 个 `register_*` 形式。6 个属于会话或全局状态辅助函数（`get_option`、`set_option`、`set_loglevel`、`set_source_metadata`、`register_view`、`from_polars`）。15 个需要未安装的测试数据写入器：BGEN、PGEN 与 VCF Zarr 的读取器和描述函数（`needs-data:bgen-fixture`、`needs-data:pgen-fixture`、`needs-data:vcf-zarr-store`），以及不接受参考序列参数的 `register_cram`（`needs-data:reference-fasta`）。2 个是通用表格读取器（`read_table`、`scan_table`；`not-dna`）。`read_bigwig` 与 `scan_bigwig` 与已有的 `read_pybigwig_intervals` 重复（`duplicate-of`）。
- gtars：11 个 `other` 行。6 个已包装。`read_tokens_from_gtok_as_strings` 是 `read_region_tokens` 的族变体。`hgvs_to_vrs_id` 需要转录本存储（`needs-data:cdot-transcript-store`）。`median_abs_distance` 是通用中位数（`not-dna`）。`force_unlock_store` 与 `store_lock_status` 管理 RefgetStore 的锁（`internal`）。
- mygene：10 个 `other` 行，全部是客户端类或工厂、缓存与端点设置（`internal`），或元数据列表（`get_fields`、`metadata`；`not-dna`）。没有任何函数返回尚未包装的基因记录。
- myvariant：10 个 `other` 行，分类与 mygene 相同。没有任何函数返回尚未包装的变异记录。

## 1) `read_file_header`

- **作用：** 按 polars-bio 的解析结果返回 VCF、BCF、SAM 或 BAM 文件的头部信息，不读取记录。
- **处理方法：** 用对应的惰性 `scan_*` 读取器打开文件，再从数据框读取 `get_metadata`。内部字段 `_datafusion_table_name` 会被删除。
- **API：** `read_file_header(path[必须], file_format[可选], verbose[可选])`
- **输入：** 本地 VCF、BCF、SAM 或 BAM 文件。`file_format` 取 `vcf`（默认）、`bcf`、`sam` 或 `bam`。
- **输出：** 字典，含 `file_format`、`coordinate_system_zero_based`、`header`（按格式给出的字段，如 `sample_names`、`contigs` 或 `sort_order`）与 `provenance`。

**示例文件:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**示例代码：**

```python
from dnakit.annotation import read_file_header

result = read_file_header("sample.vcf", file_format="vcf", verbose=False)
print(result["header"]["sample_names"])
print(result["header"]["contigs"])
```

**示例输出：**

```text
['s1', 's2']
[{'id': 'chr1', 'length': 1000}]
```

## 2) `describe_genomic_file`

- **作用：** 描述 VCF、BCF、SAM、BAM、CRAM、cool 或 Stockholm（STO）文件的结构或头部条目。
- **处理方法：** 调用 polars-bio 对应格式的 `describe_*` 函数。SAM、BAM 与 CRAM 的标签类型来自一批样本记录。CRAM 需要 `reference_path`；其他格式传入该参数会报错。
- **API：** `describe_genomic_file(path[必须], file_format[可选], reference_path[可选], verbose[可选])`
- **输入：** 本地文件。`file_format` 取 `vcf`、`bcf`、`sam`、`bam`、`cram`、`cool` 或 `sto` 之一。
- **输出：** 字典，含 `file_format`、`entry_count`、`entries`（每个条目一个字典，键随格式而变）与 `provenance`。

**示例文件:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**示例代码：**

```python
from dnakit.annotation import describe_genomic_file

result = describe_genomic_file("sample.vcf", file_format="vcf", verbose=False)
for entry in result["entries"]:
    print(entry["name"], entry["field_type"], entry["data_type"])
```

**示例输出：**

```text
DP INFO Integer
genotypes FORMAT Struct
```

## 3) `read_cram_records`

- **作用：** 读取 CRAM 比对记录，每条读段一条记录，并用参考 FASTA 解码。
- **处理方法：** 调用 polars-bio 的 `read_cram`，传入参考序列，坐标以 0 起始、半开区间报告。参考序列需要 `.fai` 索引，可由 `pysam.faidx` 生成。
- **API：** `read_cram_records(path[必须], reference_path[必须], tag_fields[可选], max_rows[可选], verbose[可选])`
- **输入：** 一个 CRAM 文件及其参考 FASTA。示例中用 pysam 由 SAM 文件生成 CRAM，因为 CRAM 是二进制格式。
- **输出：** 字典，含 `record_count`、`alignments`（name、chrom、start、end、flags、cigar 与序列字段）与 `provenance`。示例中 SAM 的 POS 5 读为起点 4，`5M2I3M` 的 POS 40 读为终点 47。

**示例文件:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:73
@RG	ID:rg1	SM:s1
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	16	chr1	20	40	8M	*	0	0	GCCGATTA	IIIIIIII	RG:Z:rg1
r3	0	chr1	40	30	5M2I3M	*	0	0	ACGTAGCTAG	IIIIIIIIII	NM:i:2
```

**示例文件:** `ref.fa`

```text
>chr1
ACGTACGTTAGCCGATTAGCATGCAATGCCGTAGCTAGCTAGGATCGATCGTAGCTAGCTAGCTACGATCGAT
```

**示例代码：**

```python
import pathlib

import pysam
from dnakit.annotation import read_cram_records

pysam.faidx("ref.fa")
pathlib.Path("reads.cram").unlink(missing_ok=True)
pysam.view("-C", "-T", "ref.fa", "-o", "reads.cram", "reads.sam", catch_stdout=False)
result = read_cram_records("reads.cram", "ref.fa", verbose=False)
for row in result["alignments"]:
    print(row["name"], row["start"], row["end"], row["cigar"])
```

**示例输出：**

```text
r1 4 14 10M
r2 19 27 8M
r3 39 47 5M2I3M
```

## 4) `read_bigbed_intervals`

- **作用：** 读取 bigBed 文件的每个区间及其存储的字段。
- **处理方法：** 调用 polars-bio 的 `read_bigbed`，使用 `schema='rest'` 与 0 起始坐标。`end` 之后以制表符分隔的文本被拆分为 `fields`。
- **API：** `read_bigbed_intervals(path[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 bigBed 文件。示例用 pybigtools 生成，因为 bigBed 是二进制格式。
- **输出：** 字典，含 `interval_count`、`intervals`（chrom、start、end 与 `fields`）与 `provenance`。

**示例代码：**

```python
import pathlib

import pybigtools
from dnakit.annotation import read_bigbed_intervals

pathlib.Path("features.bb").unlink(missing_ok=True)
writer = pybigtools.open("features.bb", "w")
writer.write(
    {"chr1": 100, "chr2": 100},
    [
        ("chr1", 0, 10, "featA\t500\t+"),
        ("chr1", 20, 35, "featB\t300\t-"),
        ("chr2", 5, 15, "featC\t100\t+"),
    ],
)
result = read_bigbed_intervals("features.bb", verbose=False)
for row in result["intervals"]:
    print(row["chrom"], row["start"], row["end"], row["fields"][0])
```

**示例输出：**

```text
chr1 0 10 featA
chr1 20 35 featB
chr2 5 15 featC
```

## 5) `read_alignment_sequences`

- **作用：** 读取 A2M 或 A3M 多序列比对，每条比对序列一条记录。
- **处理方法：** 调用 polars-bio 的 `read_a2m` 或 `read_a3m`。空位字符与 A3M 小写插入按原样保留。只有当所有序列长度相同时才给出 `alignment_columns`。
- **API：** `read_alignment_sequences(path[必须], file_format[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 A2M 或 A3M 文件，带 `>` 标题行。`file_format` 取 `a2m`（默认）或 `a3m`。
- **输出：** 字典，含 `file_format`、`sequence_count`、`alignment_columns`、`sequences`（name、description、sequence）与 `provenance`。

**示例文件:** `align.a2m`

```text
>seq1
ACGT-ACGT
>seq2
AC-TACGTT
```

**示例代码：**

```python
from dnakit.annotation import read_alignment_sequences

result = read_alignment_sequences("align.a2m", file_format="a2m", verbose=False)
print(result["sequence_count"], result["alignment_columns"])
for row in result["sequences"]:
    print(row["name"], row["sequence"])
```

**示例输出：**

```text
2 9
seq1 ACGT-ACGT
seq2 AC-TACGTT
```

## 6) `write_fastq_reads`

- **作用：** 把带名称、DNA 序列与质量字符串的读段写入新的 FASTQ 文件。
- **处理方法：** 在调用 polars-bio 的 `write_fastq` 之前，检查每个名称（不含空白）、每条序列（DNA 字母，转为大写）与 `quality_scores`（Phred+33 字符，每个碱基一个）。已存在的文件不会被覆盖。
- **API：** `write_fastq_reads(reads[必须], destination[必须], verbose[可选])`
- **输入：** 字典列表，含 `name`、`sequence` 与 `quality_scores`，可选 `description`。`read_fastq_reads` 返回的记录可直接传入。
- **输出：** 字典，含 `records_written`、`destination` 与 `provenance`。

**示例代码：**

```python
import pathlib

from dnakit.annotation import write_fastq_reads

reads = [
    {"name": "r1", "description": "first", "sequence": "acgtn", "quality_scores": "IIII#"},
    {"name": "r2", "description": None, "sequence": "GG", "quality_scores": "II"},
]
pathlib.Path("out.fastq").unlink(missing_ok=True)
result = write_fastq_reads(reads, "out.fastq", verbose=False)
print(result["records_written"])
print(pathlib.Path("out.fastq").read_text(encoding="utf-8"))
```

**示例输出：**

```text
2
@r1 first
ACGTN
+
IIII#
@r2
GG
+
II
```

## 7) `write_alignment_file`

- **作用：** 在 SAM、BAM 与 CRAM 之间转换比对记录，并保持记录顺序。
- **处理方法：** 用 polars-bio 的 `read_sam`、`read_bam` 或 `read_cram` 读取源文件，再用 `write_sam`、`write_bam` 或 `write_cram` 写出。CRAM 作为输入或输出时需要 `reference_path`。`@HD SO` 由写出器自行设置。
- **API：** `write_alignment_file(source[必须], destination[必须], source_format[可选], output_format[可选], reference_path[可选], verbose[可选])`
- **输入：** 本地 SAM、BAM 或 CRAM 文件。`source_format` 与 `output_format` 取 `sam`（源默认）、`bam`（输出默认）或 `cram`。
- **输出：** 字典，含 `source_format`、`output_format`、`records_written`、`destination` 与 `provenance`。

**示例文件:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:73
@RG	ID:rg1	SM:s1
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	16	chr1	20	40	8M	*	0	0	GCCGATTA	IIIIIIII	RG:Z:rg1
r3	0	chr1	40	30	5M2I3M	*	0	0	ACGTAGCTAG	IIIIIIIIII	NM:i:2
```

**示例代码：**

```python
import pathlib

import pysam
from dnakit.annotation import write_alignment_file

pathlib.Path("reads_out.bam").unlink(missing_ok=True)
result = write_alignment_file(
    "reads.sam",
    "reads_out.bam",
    source_format="sam",
    output_format="bam",
    verbose=False,
)
print(result["output_format"], result["records_written"])
with pysam.AlignmentFile("reads_out.bam", "rb", check_sq=False) as handle:
    print([record.reference_start + 1 for record in handle.fetch(until_eof=True)])
```

**示例输出：**

```text
bam 3
[5, 20, 40]
```

## 8) `query_genomic_tables`

- **作用：** 对本地基因组文件运行一条只读的 SELECT 或 WITH 查询。
- **处理方法：** 每个文件仅在本次调用中以小写表名注册到 polars-bio 会话，然后执行查询。即使查询失败，表也会在 `finally` 块中被删除。其他 SQL 语句以及已注册的表名都会被拒绝。若 SQL 含 `ORDER BY`，结果保持查询顺序，否则排序。
- **API：** `query_genomic_tables(tables[必须], query[必须], max_rows[可选], verbose[可选])`；`tables` 将表名映射到 `(path, file_format)`。格式：`bed`、`gtf`、`gff`、`vcf`、`bcf`、`sam`、`bam`、`fasta`、`fastq`、`a2m`、`a3m`、`sto`、`pairs`、`cool`、`bigwig`、`bigbed`。
- **输入：** 本地文件与一条 SQL 语句。查询使用 DataFusion SQL；`end` 是保留字，因此加了双引号。
- **输出：** 字典，含 `table_names`、`columns`、`row_count`、`rows` 与 `provenance`。BED 与 VCF 的起点为 0 起始，因此 VCF 的 POS 10 显示为 9。

**示例文件:** `iv.bed`

```text
chr1	0	10	featA	500	+
chr1	20	35	featB	300	-
chr2	5	15	featC	100	+
```

**示例文件:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**示例代码：**

```python
from dnakit.annotation import query_genomic_tables

result = query_genomic_tables(
    {"peaks": ("iv.bed", "bed"), "vars": ("sample.vcf", "vcf")},
    'SELECT chrom, start, "end", name FROM peaks ORDER BY start',
    verbose=False,
)
for row in result["rows"]:
    print(row["chrom"], row["start"], row["end"], row["name"])
```

**示例输出：**

```text
chr1 0 10 featA
chr2 5 15 featC
chr1 20 35 featB
```

## 9) `lola_restricted_universe`

- **作用：** 从多个用户区域集合构建受限的 LOLA 背景集。
- **处理方法：** 对 BED3 用户集合调用 gtars 的 `lola.build_restricted_universe`。重叠区域被切成互不重叠的片段，因此结果可作为 `lola_enrichment` 的背景集。
- **API：** `lola_restricted_universe(user_regions[必须], max_rows[可选], verbose[可选])`
- **输入：** BED 文件列表，每个文件是一个用户集合。每行至少包含 chrom、start 与 end，且 0 <= start < end。
- **输出：** 字典，含 `user_set_count`、`region_count`、`regions`（chrom、start、end）与 `provenance`。示例中 chr1 [0, 20) 在 5 与 10 处被切开。

**示例文件:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**示例文件:** `setB.bed`

```text
chr1	0	10
chr2	1	4
```

**示例代码：**

```python
from dnakit.annotation import lola_restricted_universe

result = lola_restricted_universe(["setA.bed", "setB.bed"], verbose=False)
for region in result["regions"]:
    print(region["chrom"], region["start"], region["end"])
```

**示例输出：**

```text
chr1 0 5
chr1 5 10
chr1 10 20
chr2 1 4
```

## 10) `lola_check_universe`

- **作用：** 检查背景 BED 文件对每个用户区域集合的覆盖情况，并给出 LOLA 的警告。
- **处理方法：** 调用 gtars 的 `lola.check_universe`。`coverage` 是与背景集重叠的用户区域所占比例。`many_to_many` 统计与多个背景区域重叠的用户区域数，`lola_redefine_user_sets` 可消除这种情况。
- **API：** `lola_check_universe(user_regions[必须], universe[必须], max_rows[可选], verbose[可选])`
- **输入：** 用户集合 BED 文件列表与一个背景 BED 文件。
- **输出：** 字典，含 `user_set_count`、`user_sets`（user_set、total_regions、regions_in_universe、coverage、many_to_many）、`warnings` 与 `provenance`。

**示例文件:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**示例文件:** `uni.bed`

```text
chr1	0	10
chr1	5	20
chr2	0	50
chr1	100	110
```

**示例代码：**

```python
from dnakit.annotation import lola_check_universe

result = lola_check_universe(["setA.bed"], "uni.bed", verbose=False)
for item in result["user_sets"]:
    print(
        item["user_set"],
        item["total_regions"],
        item["regions_in_universe"],
        item["coverage"],
        item["many_to_many"],
    )
print(len(result["warnings"]))
```

**示例输出：**

```text
0 2 2 1.0 2
1
```

## 11) `lola_redefine_user_sets`

- **作用：** 把每个用户区域集合重新定义为其重叠的背景区域。
- **处理方法：** 调用 gtars 的 `lola.redefine_user_sets`。每个用户集合被其所触及的背景区域替换，并按 chrom、start、end 排序。
- **API：** `lola_redefine_user_sets(user_regions[必须], universe[必须], max_rows[可选], verbose[可选])`
- **输入：** 用户集合 BED 文件列表与一个背景 BED 文件。
- **输出：** 字典，含 `user_set_count`、`user_sets`（user_set 与 regions）与 `provenance`。

**示例文件:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**示例文件:** `setB.bed`

```text
chr1	0	10
chr2	1	4
```

**示例文件:** `uni.bed`

```text
chr1	0	10
chr1	5	20
chr2	0	50
chr1	100	110
```

**示例代码：**

```python
from dnakit.annotation import lola_redefine_user_sets

result = lola_redefine_user_sets(["setA.bed", "setB.bed"], "uni.bed", verbose=False)
for item in result["user_sets"]:
    spans = [(region["chrom"], region["start"], region["end"]) for region in item["regions"]]
    print(item["user_set"], spans)
```

**示例输出：**

```text
0 [('chr1', 0, 10), ('chr1', 5, 20)]
1 [('chr1', 0, 10), ('chr1', 5, 20), ('chr2', 0, 50)]
```

## 12) `write_region_tokens`

- **作用：** 把区域 token 编号写入新的 GTOK 二进制文件。
- **处理方法：** 检查每个编号都是 0 到 2**32 - 1 之间的整数，然后调用 gtars 的 `utils.write_tokens_to_gtok`。编号按给定顺序存储。
- **API：** `write_region_tokens(tokens[必须], destination[必须], verbose[可选])`
- **输入：** 非空的非负整数列表，例如 `tokenize_fragments` 得到的背景集 token 编号。
- **输出：** 字典，含 `token_count`、`destination` 与 `provenance`。

**示例代码：**

```python
import pathlib

from dnakit.annotation import write_region_tokens

pathlib.Path("token_ids.gtok").unlink(missing_ok=True)
result = write_region_tokens([10, 0, 7], "token_ids.gtok", verbose=False)
print(result["token_count"], pathlib.Path(result["destination"]).name)
```

**示例输出：**

```text
3 token_ids.gtok
```

## 13) `read_region_tokens`

- **作用：** 从 GTOK 二进制文件读取区域 token 编号。
- **处理方法：** 调用 gtars 的 `utils.read_tokens_from_gtok`。编号按存储时的顺序返回。gtars 提供的字符串形式（`read_tokens_from_gtok_as_strings`）返回相同编号的文本，因此不单独包装。
- **API：** `read_region_tokens(path[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 GTOK 文件，例如由 `write_region_tokens` 写出的文件。
- **输出：** 字典，含 `token_count`、`tokens`（按存储顺序的整数）与 `provenance`。

**示例代码：**

```python
import pathlib

from dnakit.annotation import read_region_tokens, write_region_tokens

pathlib.Path("tokens.gtok").unlink(missing_ok=True)
write_region_tokens([3, 1, 2], "tokens.gtok", verbose=False)
result = read_region_tokens("tokens.gtok", verbose=False)
print(result["token_count"], result["tokens"])
```

**示例输出：**

```text
3 [3, 1, 2]
```

## 14) `vrs_location_digest`

- **作用：** 返回序列上某个区间的 GA4GH VRS SequenceLocation 摘要。
- **处理方法：** 以序列摘要（GA4GH sha512t24u）、0 起始的起点与开放的终点调用 gtars 的 `vrs.location_digest`。除非 `prefixed=False`，标识符为 `ga4gh:SL.<digest>`。
- **API：** `vrs_location_digest(sequence_digest_value[必须], start[必须], end[必须], prefixed[可选], verbose[可选])`
- **输入：** 32 个字符的 sha512t24u 序列摘要（例如来自 `sequence_digest`），以及满足 0 <= start < end 的整数。
- **输出：** 字典，含 `identifier`、`digest`、`prefixed`、`sequence_digest`、`start`、`end` 与 `provenance`。

**示例代码：**

```python
from gtars.refget import sha512t24u_digest

from dnakit.annotation import vrs_location_digest

digest = str(sha512t24u_digest("ACGTACGTAC"))
result = vrs_location_digest(digest, 0, 4, verbose=False)
print(result["identifier"])
bare = vrs_location_digest(digest, 0, 4, prefixed=False, verbose=False)
print(bare["digest"] == result["digest"])
```

**示例输出：**

```text
ga4gh:SL.xT-rAZC83CmO1DXxq6GGIm327TB1xI52
True
```
