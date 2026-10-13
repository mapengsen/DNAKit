# 70 全量 G4：polars-bio、gtars、ncls 与 genomicranges 注释函数

这一组函数把 polars-bio、gtars、ncls 与 genomicranges 中适合 DNA 与基因层面、可在本地离线运行的公共函数包装为 DNAKit 接口。每个函数只读取显式给出的本地文件，结果为可直接序列化的字典，并包含 `status`、`operation`、结果字段与 `provenance`。坐标默认为 BED 约定（0 起点、半开区间），例外会在对应小节说明。

## 使用说明

**导入路径：** 函数位于 `dnakit.annotation`、`dnakit.annotation` 与 `dnakit.annotation`。示例直接从这些模块导入。
**依赖：** 缺少对应包时调用会报错并给出安装提示；polars-bio、gtars、genomicranges 与 ncls 只在调用时导入。`read_cooler_matrix` 的示例还需要 `cooler` 与 `pandas`。
**gtars：** 上游函数定义在 `gtars.genomic_distributions`、`gtars.lola`、`gtars.refget`、`gtars.tokenizers`、`gtars.vrs` 与 `gtars.vrs.hgvs` 子模块中，顶层命名空间不直接导出函数，因此普通的顶层扫描得到的清单为空。
**ncls：** 清单中的两个索引构造函数分别是整数索引（已由 `dnakit.annotation.backends` 包装）与浮点索引（坐标不是碱基位置，不包装），本组不新增函数。
**genomicranges：** 上游 GTF 读取器保存的终止位置为 GTF 终止位置减一，`read_gtf_ranges` 会恢复 GTF 原值；`compute_up_down` 返回的是宽度而非终点，`upstream_downstream_regions` 已完成换算。

## 示例数据

下面的示例共用这些文件，示例代码运行时生成的文件（如 `copy.vcf`、`contacts.cool`）会在运行前清除旧文件。

**示例文件：** `a.bed`

```text
chr1	100	200	a
chr1	150	300	b
chr1	400	500	c
chr2	0	50	d
```

**示例文件：** `b.bed`

```text
chr1	180	220	x
chr1	450	460	y
chr2	10	20	z
```

**示例文件：** `genes.gtf`

```text
chr1	ens	gene	1001	2000	.	+	.	gene_id "G1"; gene_name "ABC"; gene_type "protein_coding";
chr1	ens	exon	1001	1200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";
chr1	ens	exon	1801	2000	.	+	.	gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";
chr1	ens	gene	5001	6000	.	-	.	gene_id "G2"; gene_name "DEF"; gene_type "protein_coding";
chr1	ens	exon	5001	6000	.	-	.	gene_id "G2"; transcript_id "T2"; gene_type "protein_coding";
```

**示例文件：** `genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	gene	800	900	.	-	.	ID=gene2;Name=GENE2
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**示例文件：** `genome.fa`

```text
>chr1 first
ACGTACGTAC
GGTT
>chr2
TTTTAAAACCCC
```

**示例文件：** `reads.fq`

```text
@r1
ACGTN
+
IIII#
@r2
GGCCAA
+
IIIIII
```

**示例文件：** `aln.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	11	60	5M	*	0	0	ACGTA	IIIII	NM:i:0
r2	0	chr1	13	60	4M	*	0	0	GTAC	IIII	NM:i:1
```

**示例文件：** `contacts.pairs`

```text
## pairs format v1.0
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
r1	chr1	100	chr1	500	+	-
r2	chr1	200	chr2	40	+	+
```

**示例文件：** `align.sto`

```text
# STOCKHOLM 1.0
#=GF ID test_family
seq1  ACGT-ACGT
seq2  ACGTTACGA
//
```

**示例文件：** `gc.bed`

```text
chr1	0	10
chr1	10	14
chr2	0	4
chr2	8	12
```

**示例文件：** `regions.bed`

```text
chr1	950	1010	p1
chr1	1500	1600	p2
chr1	5950	6100	p3
chr1	9000	9100	p4
```

**示例文件：** `signal.tsv`

```text
region	condA	condB
chr1_0_100	1.0	2.0
chr2_0_100	3.0	4.0
```

**示例文件：** `frags.tsv`

```text
chr1	10	20	AAA	1
chr1	150	160	AAA	2
chr2	5	9	BBB	3
```

**示例文件：** `universe.bed`

```text
chr1	0	100
chr1	100	200
chr2	0	50
```

**示例文件：** `user.bed`

```text
chr1	0	10
chr1	20	30
```

**示例文件：** `db1.bed`

```text
chr1	0	10
```

**示例文件：** `db2.bed`

```text
chr1	20	30
chr2	0	50
```

**示例文件：** `set1.bed`

```text
chr1	100	200
chr1	300	400
```

**示例文件：** `set2.bed`

```text
chr1	100	200
chr1	500	600
```

## polars-bio：区间、基因注释、变异、序列与比对

### 1) `dnakit.annotation.overlap_intervals()`

- **作用：** 列出两个 BED 文件中所有重叠的区间对。
- **处理方法：** 用 polars-bio 读取两个文件（0 起点、半开区间），调用 `overlap`，结果按坐标排序。
- **API：** `overlap_intervals(left, right, *, max_rows=100000, verbose=True)`
- **输入：** 两个 BED 文件路径；`max_rows` 限制输出行数。
- **输出：** `pairs` 列表，每项包含两侧的 `chrom_*`、`start_*`、`end_*`、`name_*`；另有 `pair_count`。
- **示例代码：**

```python
from dnakit.annotation import overlap_intervals

result = overlap_intervals("a.bed", "b.bed")
print(result["pair_count"])
print([(row["start_1"], row["end_1"], row["start_2"], row["end_2"]) for row in result["pairs"]])
```

- **示例输出：**

```text
4
[(100, 200, 180, 220), (150, 300, 180, 220), (400, 500, 450, 460), (0, 50, 10, 20)]
```

### 2) `dnakit.annotation.count_interval_overlaps()`

- **作用：** 统计 `other` 中与 `source` 每个区间重叠的区间数。
- **处理方法：** 调用 polars-bio 的 `count_overlaps`，每个源区间保留一行，没有重叠的区间计数为 0。
- **API：** `count_interval_overlaps(source, other, *, max_rows=100000, verbose=True)`
- **输入：** 源 BED 文件与比较 BED 文件；`max_rows` 限制输出行数。
- **输出：** `intervals` 列表（含 `count`），以及 `total_overlaps` 与 `interval_count`。
- **示例代码：**

```python
from dnakit.annotation import count_interval_overlaps

result = count_interval_overlaps("a.bed", "b.bed")
print(result["total_overlaps"])
print([(row["chrom"], row["start"], row["count"]) for row in result["intervals"]])
```

- **示例输出：**

```text
4
[('chr1', 100, 1), ('chr1', 150, 1), ('chr1', 400, 1), ('chr2', 0, 1)]
```

### 3) `dnakit.annotation.cluster_intervals()`

- **作用：** 把重叠或距离不超过 `min_distance` 的区间归为一个簇。
- **处理方法：** 调用 polars-bio 的 `cluster`，再按坐标顺序重新编号簇，使结果与上游处理顺序无关。
- **API：** `cluster_intervals(source, *, min_distance=0, max_rows=100000, verbose=True)`
- **输入：** BED 文件；`min_distance` 为非负整数（碱基数）。
- **输出：** `intervals`（含 `cluster`、`cluster_start`、`cluster_end`）与 `cluster_count`。
- **示例代码：**

```python
from dnakit.annotation import cluster_intervals

result = cluster_intervals("a.bed", min_distance=0)
print(result["cluster_count"])
print([(row["name"], row["cluster"]) for row in result["intervals"]])
```

- **示例输出：**

```text
3
[('a', 0), ('b', 0), ('c', 1), ('d', 2)]
```

### 4) `dnakit.annotation.complement_intervals()`

- **作用：** 返回每条染色体上未被区间覆盖的间隙。
- **处理方法：** 用 `chrom_sizes` 构造视图区间 [0, 长度)，调用 polars-bio 的 `complement` 求差集。
- **API：** `complement_intervals(source, *, chrom_sizes, max_rows=100000, verbose=True)`
- **输入：** BED 文件；`chrom_sizes` 为 {染色体: 长度} 映射。
- **输出：** `gaps` 列表（`chrom`、`start`、`end`）与 `gap_count`。
- **示例代码：**

```python
from dnakit.annotation import complement_intervals

result = complement_intervals("a.bed", chrom_sizes={"chr1": 1000, "chr2": 100})
print(result["gap_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

- **示例输出：**

```text
4
[('chr1', 0, 100), ('chr1', 300, 400), ('chr1', 500, 1000), ('chr2', 50, 100)]
```

### 5) `dnakit.annotation.subtract_intervals()`

- **作用：** 从 `source` 的区间中去掉 `other` 覆盖的碱基。
- **处理方法：** 调用 polars-bio 的 `subtract`；被部分覆盖的区间会被切成剩余片段。
- **API：** `subtract_intervals(source, other, *, max_rows=100000, verbose=True)`
- **输入：** 源 BED 文件与被减去的 BED 文件；`max_rows` 限制输出行数。
- **输出：** `pieces` 列表（保留源文件的列）与 `piece_count`。
- **示例代码：**

```python
from dnakit.annotation import subtract_intervals

result = subtract_intervals("a.bed", "b.bed")
print(result["piece_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["pieces"]][:3])
```

- **示例输出：**

```text
7
[('chr1', 100, 180), ('chr1', 150, 180), ('chr1', 220, 300)]
```

### 6) `dnakit.annotation.read_gene_annotation()`

- **作用：** 读取 GTF 或 GFF3 基因注释，每个特征一条记录。
- **处理方法：** 按 `file_format` 调用 polars-bio 的 `read_gtf` 或 `read_gff`，属性字段由 `attributes` 指定。
- **API：** `read_gene_annotation(path, *, file_format="gtf", attributes=None, max_rows=100000, verbose=True)`
- **输入：** GTF 或 GFF3 文件；`file_format` 为 `gtf` 或 `gff`。
- **输出：** `features` 列表（`chrom`、`start`、`end`、`type`、`strand` 与属性列）与 `feature_count`。
- **示例代码：**

```python
from dnakit.annotation import read_gene_annotation

result = read_gene_annotation("genes.gtf", attributes=["gene_id", "gene_name"])
print(result["feature_count"])
print([(row["type"], row["gene_id"], row["start"], row["end"]) for row in result["features"]])
gff = read_gene_annotation("genes.gff3", file_format="gff")
print(gff["feature_count"], sorted({row["type"] for row in gff["features"]}))
```

- **示例输出：**

```text
5
[('exon', 'G1', 1000, 1200), ('gene', 'G1', 1000, 2000), ('exon', 'G1', 1800, 2000), ('exon', 'G2', 5000, 6000), ('gene', 'G2', 5000, 6000)]
4 ['exon', 'gene', 'mRNA']
```

### 7) `dnakit.annotation.read_variant_calls()`

- **作用：** 读取 VCF 或 BCF 变异文件，每个变异一条记录，包含 INFO 与 FORMAT 字段。
- **处理方法：** 按 `file_format` 调用 polars-bio 的 `read_vcf` 或 `read_bcf`，坐标转为 0 起点半开区间。
- **API：** `read_variant_calls(path, *, file_format="vcf", info_fields=None, format_fields=None, max_rows=100000, verbose=True)`
- **输入：** VCF 或 BCF 文件；`info_fields` 与 `format_fields` 指定要读取的字段。
- **输出：** `variants` 列表（`chrom`、`start`、`end`、`id`、`ref`、`alt`、INFO 字段与 `genotypes`）与 `variant_count`。
- **示例代码：**

```python
from dnakit.annotation import read_variant_calls

result = read_variant_calls("variants.vcf", info_fields=["DP"], format_fields=["GT"])
print(result["variant_count"])
print([(row["start"], row["end"], row["id"], row["DP"]) for row in result["variants"]])
```

- **示例输出：**

```text
2
[(9, 10, 'rs1', 12), (24, 25, '', 8)]
```

### 8) `dnakit.annotation.write_vcf_variants()`

- **作用：** 读取一个 VCF 文件并把记录写成新的 VCF 文件。
- **处理方法：** 用 polars-bio 读取源文件，再调用 `write_vcf`；已存在的输出文件不会被覆盖。
- **API：** `write_vcf_variants(source, destination, *, info_fields=None, format_fields=None, verbose=True)`
- **输入：** 源 VCF 路径与输出路径；`info_fields` 与 `format_fields` 指定要保留的字段。
- **输出：** `records_written`、`destination` 与 `status`。
- **示例代码：**

```python
import pathlib

from dnakit.annotation import write_vcf_variants

pathlib.Path("copy.vcf").unlink(missing_ok=True)
result = write_vcf_variants("variants.vcf", "copy.vcf", info_fields=["DP"], format_fields=["GT"])
print(result["records_written"], pathlib.Path(result["destination"]).name)
```

- **示例输出：**

```text
2 copy.vcf
```

### 9) `dnakit.annotation.read_fasta_sequences()`

- **作用：** 读取 FASTA 文件中的每条序列及其长度。
- **处理方法：** 调用 polars-bio 的 `read_fasta`，序列保留原始大小写，结果按名称排序。
- **API：** `read_fasta_sequences(path, *, max_rows=100000, verbose=True)`
- **输入：** FASTA 文件路径。
- **输出：** `sequences` 列表（`name`、`description`、`sequence`、`length`）与 `sequence_count`。
- **示例代码：**

```python
from dnakit.annotation import read_fasta_sequences

result = read_fasta_sequences("genome.fa")
print(result["sequence_count"], [(row["name"], row["length"]) for row in result["sequences"]])
```

- **示例输出：**

```text
2 [('chr1', 14), ('chr2', 12)]
```

### 10) `dnakit.annotation.read_fastq_reads()`

- **作用：** 读取 FASTQ 文件中的每条测序读段及其质量字符串。
- **处理方法：** 调用 polars-bio 的 `read_fastq`，质量字符串按 Phred+33 原样保留。
- **API：** `read_fastq_reads(path, *, max_rows=100000, verbose=True)`
- **输入：** FASTQ 文件路径。
- **输出：** `reads` 列表（`name`、`sequence`、`quality_scores`、`length`）与 `read_count`。
- **示例代码：**

```python
from dnakit.annotation import read_fastq_reads

result = read_fastq_reads("reads.fq")
print(result["read_count"], [(row["name"], row["sequence"], row["quality_scores"]) for row in result["reads"]])
```

- **示例输出：**

```text
2 [('r1', 'ACGTN', 'IIII#'), ('r2', 'GGCCAA', 'IIIIII')]
```

### 11) `dnakit.annotation.write_fasta_sequences()`

- **作用：** 把 DNA 序列（`DNA` 记录或 ID 到序列的映射）写成 FASTA 文件。
- **处理方法：** 序列经 DNA 字母表检查后转为大写，用 polars-bio 的 `write_fasta` 写出；已存在的输出文件不会被覆盖。
- **API：** `write_fasta_sequences(sequences, destination, *, verbose=True)`
- **输入：** 序列映射或 `DNA` 对象；输出路径。
- **输出：** `records_written`、`destination` 与 `status`。
- **示例代码：**

```python
import pathlib

from dnakit.annotation import write_fasta_sequences

pathlib.Path("written.fa").unlink(missing_ok=True)
result = write_fasta_sequences({"s1": "acgtn", "s2": "GGCC"}, "written.fa")
print(result["records_written"])
print(pathlib.Path("written.fa").read_text(encoding="utf-8"))
```

- **示例输出：**

```text
2
>s1
ACGTN
>s2
GGCC
```

### 12) `dnakit.annotation.read_alignment_records()`

- **作用：** 读取 SAM 或 BAM 比对记录，每条读段一行，包含 CIGAR 与序列。
- **处理方法：** 按 `file_format` 调用 polars-bio 的 `read_sam` 或 `read_bam`，坐标为 0 起点半开区间。
- **API：** `read_alignment_records(path, *, file_format="sam", tag_fields=None, max_rows=100000, verbose=True)`
- **输入：** SAM 或 BAM 文件；`tag_fields` 指定要读取的 SAM 标签。
- **输出：** `alignments` 列表（`name`、`chrom`、`start`、`end`、`cigar`、`sequence` 等）与 `record_count`。
- **示例代码：**

```python
from dnakit.annotation import read_alignment_records

result = read_alignment_records("aln.sam")
print(result["record_count"])
print([(row["name"], row["start"], row["end"], row["cigar"]) for row in result["alignments"]])
```

- **示例输出：**

```text
2
[('r1', 10, 15, '5M'), ('r2', 12, 16, '4M')]
```

### 13) `dnakit.annotation.alignment_depth()`

- **作用：** 返回比对文件中读段深度保持不变的区段。
- **处理方法：** 调用 polars-bio 的 `depth`，按 `filter_flag` 与 `min_mapping_quality` 过滤读段（默认过滤未比对、重复和未通过质控的读段）。
- **API：** `alignment_depth(path, *, min_mapping_quality=0, filter_flag=1796, max_rows=100000, verbose=True)`
- **输入：** SAM 或 BAM 文件；过滤参数。
- **输出：** `runs` 列表（`chrom`、`start`、`end`、`depth`）与 `run_count`。
- **示例代码：**

```python
from dnakit.annotation import alignment_depth

result = alignment_depth("aln.sam")
print([(row["start"], row["end"], row["depth"]) for row in result["runs"]])
```

- **示例输出：**

```text
[(10, 12, 1), (12, 15, 2), (15, 16, 1)]
```

### 14) `dnakit.annotation.fastqc_summary()`

- **作用：** 对 FASTQ 文件做一次流式扫描，给出 FastQC 风格的基础统计与模块状态。
- **处理方法：** 调用 polars-bio 的 `fastqc`，汇总 `basic_stats`、模块状态与每条序列的 GC 分布；GC 百分比不计 N 碱基。
- **API：** `fastqc_summary(path, *, modules=None, verbose=True)`
- **输入：** FASTQ 文件；`modules` 可选，限定计算的模块名。
- **输出：** `basic_statistics`、`module_status` 与 `per_sequence_gc`。
- **示例代码：**

```python
from dnakit.annotation import fastqc_summary

result = fastqc_summary("reads.fq")
print(result["basic_statistics"])
print(result["module_status"]["basic_stats"])
```

- **示例输出：**

```text
{'gc_pct': 60.0, 'max_len': 6, 'min_len': 5, 'n_seq': 2, 'total_bases': 11}
PASS
```

### 15) `dnakit.annotation.read_contact_pairs()`

- **作用：** 读取 4DN Pairs 格式的 Hi-C 读段对。
- **处理方法：** 调用 polars-bio 的 `read_pairs`，位置转为 0 起点坐标，列名沿用文件表头。
- **API：** `read_contact_pairs(path, *, max_rows=100000, verbose=True)`
- **输入：** Pairs 文件路径。
- **输出：** `pairs` 列表与 `pair_count`。
- **示例代码：**

```python
from dnakit.annotation import read_contact_pairs

result = read_contact_pairs("contacts.pairs")
print([(row["readID"], row["pos1"], row["pos2"]) for row in result["pairs"]])
```

- **示例输出：**

```text
[('r1', 99, 499), ('r2', 199, 39)]
```

### 16) `dnakit.annotation.read_bed_intervals()`

- **作用：** 读取 BED 文件，每个区间一条记录，保留名称、分数等列。
- **处理方法：** 调用 polars-bio 的 `read_bed`，坐标为 0 起点半开区间，结果按坐标排序。
- **API：** `read_bed_intervals(path, *, max_rows=100000, verbose=True)`
- **输入：** BED 文件路径。
- **输出：** `intervals` 列表与 `interval_count`。
- **示例代码：**

```python
from dnakit.annotation import read_bed_intervals

result = read_bed_intervals("a.bed")
print(result["interval_count"], [(row["name"], row["start"], row["end"]) for row in result["intervals"]][:2])
```

- **示例输出：**

```text
4 [('a', 100, 200), ('b', 150, 300)]
```

### 17) `dnakit.annotation.read_cooler_matrix()`

- **作用：** 读取 Cooler（.cool 或 .mcool）格式的 Hi-C 接触矩阵，每个非零的 bin 对一条记录。
- **处理方法：** 调用 polars-bio 的 `read_cool`；`.mcool` 文件需通过 `resolution` 指定分辨率。
- **API：** `read_cooler_matrix(path, *, resolution=None, include_weights=False, max_rows=100000, verbose=True)`
- **输入：** Cooler 文件；`resolution` 与 `include_weights` 可选。
- **输出：** `contacts` 列表（`chrom1`、`start1`、`chrom2`、`start2`、`count`）、`contact_count` 与 `total_counts`。
- **示例代码：**

```python
import pathlib

import cooler
import pandas as pd

from dnakit.annotation import read_cooler_matrix

bins = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 100], "end": [100, 200]})
pixels = pd.DataFrame({"bin1_id": [0, 0, 1], "bin2_id": [0, 1, 1], "count": [3, 2, 5]})
pathlib.Path("contacts.cool").unlink(missing_ok=True)
cooler.create_cooler("contacts.cool", bins, pixels, ordered=False)
result = read_cooler_matrix("contacts.cool")
print(result["contact_count"], result["total_counts"])
print([(row["start1"], row["start2"], row["count"]) for row in result["contacts"]])
```

- **示例输出：**

```text
3 10
[(0, 0, 3), (0, 100, 2), (100, 100, 5)]
```

### 18) `dnakit.annotation.read_stockholm_alignments()`

- **作用：** 读取 Stockholm 多序列比对，每条比对序列一条记录。
- **处理方法：** 调用 polars-bio 的 `read_sto`，保留比对中的间隙字符。
- **API：** `read_stockholm_alignments(path, *, max_rows=100000, verbose=True)`
- **输入：** Stockholm 文件路径（`.sto` 或 `.stk`）。
- **输出：** `alignments` 列表（`alignment_id`、`name`、`sequence`）、`sequence_count` 与 `alignment_columns`。
- **示例代码：**

```python
from dnakit.annotation import read_stockholm_alignments

result = read_stockholm_alignments("align.sto")
print(result["sequence_count"], result["alignment_columns"])
print([(row["name"], row["sequence"]) for row in result["alignments"]])
```

- **示例输出：**

```text
2 9
[('seq1', 'ACGT-ACGT'), ('seq2', 'ACGTTACGA')]
```

## gtars：基因组序列、区域集、分区、LOLA、token 与 VRS

### 19) `dnakit.annotation.region_gc_content()`

- **作用：** 计算每个 BED 区域在参考基因组中的 GC 比例。
- **处理方法：** 用 gtars 的 `GenomeAssembly` 读取 FASTA，调用 `calc_gc_content`；N 碱基不计入 GC。
- **API：** `region_gc_content(regions, genome_fasta, *, ignore_unknown_chromosomes=False, max_rows=100000, verbose=True)`
- **输入：** BED 区域文件与 FASTA 基因组；`ignore_unknown_chromosomes` 决定是否跳过 FASTA 中不存在的染色体。
- **输出：** `regions`（`chrom`、`start`、`end`、`gc_fraction`）与 `mean_gc_fraction`。
- **示例代码：**

```python
from dnakit.annotation import region_gc_content

result = region_gc_content("gc.bed", "genome.fa")
print(result["mean_gc_fraction"])
print([(row["chrom"], row["start"], row["end"], row["gc_fraction"]) for row in result["regions"]])
```

- **示例输出：**

```text
0.5
[('chr1', 0, 10, 0.5), ('chr1', 10, 14, 0.5), ('chr2', 0, 4, 0.0), ('chr2', 8, 12, 1.0)]
```

### 20) `dnakit.annotation.region_dinucleotide_frequencies()`

- **作用：** 统计每个 BED 区域内 16 种二核苷酸的出现次数或百分比。
- **处理方法：** 调用 gtars 的 `calc_dinucl_freq`；二核苷酸名称转为大写，并给出所有区域的合计。
- **API：** `region_dinucleotide_frequencies(regions, genome_fasta, *, raw_counts=False, ignore_unknown_chromosomes=False, max_rows=100000, verbose=True)`
- **输入：** BED 区域文件与 FASTA 基因组；`raw_counts=True` 返回计数，否则返回百分比。
- **输出：** `regions`（每个区域的 `frequencies`）、`pooled` 合计、`dinucleotides` 与 `value_unit`。
- **示例代码：**

```python
from dnakit.annotation import region_dinucleotide_frequencies

result = region_dinucleotide_frequencies("gc.bed", "genome.fa", raw_counts=True)
print(result["value_unit"], result["regions"][0]["label"])
print(result["regions"][0]["frequencies"]["AC"], result["regions"][0]["frequencies"]["CG"])
```

- **示例输出：**

```text
count chr1_0_10
3.0 2.0
```

### 21) `dnakit.annotation.partition_region_counts()`

- **作用：** 把 BED 区域归入核心启动子、近端启动子、外显子、内含子或基因间区，并统计数目。
- **处理方法：** 用 gtars 的 `PartitionList.from_gtf` 建立分区，再调用 `calc_partitions`；给出 `chrom_sizes` 时改用 `calc_expected_partitions`，并返回卡方检验 p 值。
- **API：** `partition_region_counts(regions, gene_annotation, *, core_promoter_bp=100, proximal_promoter_bp=2000, protein_coding_only=True, chrom_sizes=None, bp_proportion=False, verbose=True)`
- **输入：** BED 区域文件与 GTF 基因注释；启动子长度与是否只用蛋白编码基因可调。
- **输出：** `partitions`（各分区的观察值，若有 `chrom_sizes` 则含期望值与 p 值）与 `total`。
- **示例代码：**

```python
from dnakit.annotation import partition_region_counts

observed = partition_region_counts("regions.bed", "genes.gtf")
print(observed["total"], [(row["partition"], row["observed"]) for row in observed["partitions"]])
expected = partition_region_counts("regions.bed", "genes.gtf", chrom_sizes={"chr1": 10000})
print(expected["mode"], [round(row["expected"], 2) for row in expected["partitions"]])
```

- **示例输出：**

```text
4 [('promoterCore', 2), ('promoterProx', 0), ('exon', 0), ('intron', 1), ('intergenic', 1)]
expected [0.08, 1.12, 0.56, 0.24, 2.0]
```

### 22) `dnakit.annotation.signal_matrix_summary()`

- **作用：** 取出与查询区域重叠的标记区域信号值，并给出每个条件的汇总统计。
- **处理方法：** 用 gtars 的 `SignalMatrix.from_tsv` 读取信号表，调用 `calc_summary_signal`。信号表第一列为 `chrom_start_end` 标签，其余列是各条件的值。
- **API：** `signal_matrix_summary(regions, signal_table, *, verbose=True)`
- **输入：** BED 区域文件与信号表（TSV）。
- **输出：** `regions`（每个区域的 `values`）、`condition_names` 与 `matrix_stats`。
- **示例代码：**

```python
from dnakit.annotation import signal_matrix_summary

result = signal_matrix_summary("gc.bed", "signal.tsv")
print(result["condition_names"])
print([(row["label"], row["values"]) for row in result["regions"]])
```

- **示例输出：**

```text
['condA', 'condB']
[('chr1_0_10', {'condA': 1.0, 'condB': 2.0}), ('chr1_10_14', {'condA': 1.0, 'condB': 2.0}), ('chr2_0_4', {'condA': 3.0, 'condB': 4.0}), ('chr2_8_12', {'condA': 3.0, 'condB': 4.0})]
```

### 23) `dnakit.annotation.consensus_regions()`

- **作用：** 统计每个不同区域被多少个输入 BED 文件包含。
- **处理方法：** 把每个文件读为 gtars `RegionSet`，调用 `consensus` 汇总，结果按坐标排序。
- **API：** `consensus_regions(region_files, *, max_rows=100000, verbose=True)`
- **输入：** 至少两个 BED 文件路径。
- **输出：** `regions`（`chrom`、`start`、`end`、`count`）、`file_count` 与 `region_count`。
- **示例代码：**

```python
from dnakit.annotation import consensus_regions

result = consensus_regions(["set1.bed", "set2.bed"])
print([(row["chrom"], row["start"], row["end"], row["count"]) for row in result["regions"]])
```

- **示例输出：**

```text
[('chr1', 100, 200, 2), ('chr1', 300, 400, 1), ('chr1', 500, 600, 1)]
```

### 24) `dnakit.annotation.lola_enrichment()`

- **作用：** 检验用户区域与每个数据库区域集的重叠是否超过随机预期（LOLA 风格的富集分析）。
- **处理方法：** 用 gtars 的 `RegionDB.from_bed_files` 建库，再调用 `run_lola`；用户区域与数据库区域都必须落在给定的背景区域内。
- **API：** `lola_enrichment(user_regions, universe, database_files, *, min_overlap=1, direction="enrichment", max_rows=100000, verbose=True)`
- **输入：** 用户 BED 文件、背景 BED3 文件、数据库 BED 文件列表（文件名需唯一）。
- **输出：** `results`（每个数据库区域集的 `p_value_log10`、`odds_ratio`、`support`、列联表 b、c、d 等）。
- **示例代码：**

```python
from dnakit.annotation import lola_enrichment

result = lola_enrichment("user.bed", "universe.bed", ["db1.bed", "db2.bed"])
print([(row["db_name"], row["support"], row["b"], row["c"], row["d"]) for row in result["results"]])
```

- **示例输出：**

```text
[('db1.bed', 1, 0, 1, 1), ('db2.bed', 1, 1, 1, 0)]
```

### 25) `dnakit.annotation.tokenize_fragments()`

- **作用：** 把片段文件中的记录映射到通用区域（BED3）上，每个条形码得到一组区域 token。
- **处理方法：** 用 gtars 的 `Tokenizer.from_bed` 构建分词器，再调用 `tokenize_fragment_file`；token 编号即通用区域的行号。
- **API：** `tokenize_fragments(fragments, universe, *, max_rows=100000, verbose=True)`
- **输入：** 片段文件（至少含 chrom、start、end、barcode 列）与 BED3 通用区域文件。
- **输出：** `barcodes`（`barcode`、`token_ids`、`token_count`）与 `universe_size`。
- **示例代码：**

```python
from dnakit.annotation import tokenize_fragments

result = tokenize_fragments("frags.tsv", "universe.bed")
print(result["universe_size"], [(row["barcode"], row["token_ids"]) for row in result["barcodes"]])
```

- **示例输出：**

```text
3 [('AAA', [0, 1]), ('BBB', [2])]
```

### 26) `dnakit.annotation.digest_fasta_sequences()`

- **作用：** 为 FASTA 文件中的每条序列计算 GA4GH sha512t24u 与 MD5 摘要，并给出集合摘要。
- **处理方法：** 调用 gtars 的 `refget.digest_fasta`；摘要按 GA4GH refget 规范计算。
- **API：** `digest_fasta_sequences(fasta, *, max_rows=100000, verbose=True)`
- **输入：** FASTA 文件路径（未压缩）。
- **输出：** `sequences`（`name`、`length`、`sha512t24u`、`md5`）、`collection_digest` 与 `sequence_count`。
- **示例代码：**

```python
from dnakit.annotation import digest_fasta_sequences

result = digest_fasta_sequences("genome.fa")
print(result["collection_digest"])
print([(row["name"], row["length"], row["sha512t24u"]) for row in result["sequences"]])
```

- **示例输出：**

```text
rM-CtSsxdCWLtxc1sQpCDYBRFxBRXGP2
[('chr1', 14, 'bxzNj2Cq11mJJ1GuetSIcq2lO4sKlIds'), ('chr2', 12, 'yPlv_G-_BuLFIFYisu1dBipO_I3_bLNz')]
```

### 27) `dnakit.annotation.fasta_index_records()`

- **作用：** 返回 FASTA 文件的索引记录（序列名与长度），与 samtools faidx 的名称和长度一致。
- **处理方法：** 调用 gtars 的 `refget.compute_fai`；仅支持未压缩的 FASTA。
- **API：** `fasta_index_records(fasta, *, max_rows=100000, verbose=True)`
- **输入：** 未压缩的 FASTA 文件路径。
- **输出：** `records`（`name`、`length`）与 `sequence_count`。
- **示例代码：**

```python
from dnakit.annotation import fasta_index_records

result = fasta_index_records("genome.fa")
print([(row["name"], row["length"]) for row in result["records"]])
```

- **示例输出：**

```text
[('chr1', 14), ('chr2', 12)]
```

### 28) `dnakit.annotation.sequence_digest()`

- **作用：** 计算一条 DNA 序列的 GA4GH sha512t24u 或 MD5 摘要。
- **处理方法：** 序列先经 DNA 字母表检查并转为大写，再调用 gtars 的 `sha512t24u_digest` 或 `md5_digest`。
- **API：** `sequence_digest(sequence, *, algorithm="sha512t24u", verbose=True)`
- **输入：** DNA 序列字符串；`algorithm` 为 `sha512t24u` 或 `md5`。
- **输出：** `digest`、`algorithm` 与 `length`。
- **示例代码：**

```python
from dnakit.annotation import sequence_digest

print(sequence_digest("acgt")["digest"])
print(sequence_digest("ACGT", algorithm="md5")["digest"])
```

- **示例输出：**

```text
aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2
f1f8f4bf413b16ad135722aa4591043e
```

### 29) `dnakit.annotation.normalize_variant_allele()`

- **作用：** 去掉等位基因两端共有的碱基，并按 VRS 规则把变异扩展到相邻重复区域，得到规范表示。
- **处理方法：** 先检查 `ref_allele` 与参考序列在 `start` 处一致，再调用 gtars 的 `vrs.normalize_allele`；例如在 TTTT 中删除一个 T 会表示为覆盖整段的 TTT。
- **API：** `normalize_variant_allele(sequence, start, ref_allele, alt_allele, *, verbose=True)`
- **输入：** 参考序列、0 起点位置、参考等位基因与变异等位基因。
- **输出：** `start`、`end`（0 起点，`end` 为开区间）与 `allele`。
- **示例代码：**

```python
from dnakit.annotation import normalize_variant_allele

snv = normalize_variant_allele("ACGTACGTAC", 2, "G", "T")
print(snv["start"], snv["end"], snv["allele"])
deletion = normalize_variant_allele("ACGTTTTA", 2, "GT", "G")
print(deletion["start"], deletion["end"], repr(deletion["allele"]))
```

- **示例输出：**

```text
2 3 T
3 7 'TTT'
```

### 30) `dnakit.annotation.vrs_allele_id()`

- **作用：** 计算一个等位基因在序列上的 GA4GH VRS 标识符（`ga4gh:VA.<摘要>`）。
- **处理方法：** 调用 gtars 的 `vrs.vrs_id`；`prefixed=False` 时调用 `vrs.vrs_digest`，只返回不带前缀的摘要。
- **API：** `vrs_allele_id(sequence_digest_value, start, end, alt, *, prefixed=True, verbose=True)`
- **输入：** 序列的 sha512t24u 摘要、0 起点的 `start`、开区间 `end` 与变异等位基因 `alt`。
- **输出：** `identifier`、`prefixed` 与输入的摘要和坐标。
- **示例代码：**

```python
from dnakit.annotation import sequence_digest, vrs_allele_id

digest = sequence_digest("ACGTACGTACGGTT")["digest"]
print(vrs_allele_id(digest, 1, 2, "G")["identifier"])
print(vrs_allele_id(digest, 1, 2, "G", prefixed=False)["identifier"])
```

- **示例输出：**

```text
ga4gh:VA.BwLTTHJnaEI7FfDFxp8JjUw08nmooAjb
BwLTTHJnaEI7FfDFxp8JjUw08nmooAjb
```

### 31) `dnakit.annotation.parse_hgvs_variant()`

- **作用：** 解析一条 HGVS 表达式，例如 `NC_000001.11:g.12345A>G`，并给出登录号、参考类型与编辑内容。
- **处理方法：** 调用 gtars 的 `vrs.hgvs.parse_hgvs`，把结果转为 JSON 兼容的字典。
- **API：** `parse_hgvs_variant(text, *, verbose=True)`
- **输入：** HGVS 字符串。
- **输出：** `accession`、`reference_type` 与完整的 `variant` 字典。
- **示例代码：**

```python
from dnakit.annotation import parse_hgvs_variant

result = parse_hgvs_variant("NC_000001.11:g.12345A>G")
print(result["accession"], result["reference_type"])
print(result["variant"]["pos_edit"]["edit"])
```

- **示例输出：**

```text
NC_000001.11 g
{'kind': 'substitution', 'ref': 'A', 'alt': 'G'}
```

## genomicranges：GTF 区间、启动子窗口与分箱

### 32) `dnakit.annotation.read_gtf_ranges()`

- **作用：** 用 genomicranges 读取 GTF 文件，每个特征一条记录，并保留全部属性列。
- **处理方法：** 调用 genomicranges 的 `read_gtf`，把 1 起点的起始位置转为 0 起点；上游保存的终止位置为 GTF 终止位置减一，本函数将其恢复为 GTF 原值。
- **API：** `read_gtf_ranges(path, *, max_rows=100000, verbose=True)`
- **输入：** GTF 文件路径。
- **输出：** `features` 列表（`chrom`、`start`、`end`、`strand` 与属性列）与 `feature_count`。
- **示例代码：**

```python
from dnakit.annotation import read_gtf_ranges

result = read_gtf_ranges("genes.gtf")
print(result["feature_count"])
print([(row["feature"], row["gene_id"], row["start"], row["end"], row["strand"]) for row in result["features"]])
```

- **示例输出：**

```text
5
[('exon', 'G1', 1000, 1200, '+'), ('gene', 'G1', 1000, 2000, '+'), ('exon', 'G1', 1800, 2000, '+'), ('exon', 'G2', 5000, 6000, '-'), ('gene', 'G2', 5000, 6000, '-')]
```

### 33) `dnakit.annotation.upstream_downstream_regions()`

- **作用：** 在每个 BED 区间的 TSS（转录起始位点）或 TES（转录终止位点）周围生成启动子或终止子窗口。
- **处理方法：** 调用 genomicranges 的 `compute_up_down`，输入转为 1 起点的闭区间与整数链方向（-1 为负链）；上游返回的是宽度而非终点，本函数将其换算为窗口终点，并把越界的起点截断为 0。
- **API：** `upstream_downstream_regions(source, *, upstream=1000, downstream=0, site="TSS", max_rows=100000, verbose=True)`
- **输入：** BED 文件；`upstream` 与 `downstream` 为碱基数；`site` 为 `TSS` 或 `TES`。
- **输出：** `windows`（`chrom`、`start`、`end`、`strand` 与原区间坐标）、`window_count` 与 `site`。
- **示例代码：**

```python
from dnakit.annotation import upstream_downstream_regions

tss = upstream_downstream_regions("regions.bed", upstream=100, downstream=50)
print([(row["source_start"], row["start"], row["end"]) for row in tss["windows"]])
tes = upstream_downstream_regions("regions.bed", upstream=100, downstream=50, site="TES")
print([(row["source_start"], row["start"], row["end"]) for row in tes["windows"]][:2])
```

- **示例输出：**

```text
[(950, 850, 1000), (1500, 1400, 1550), (5950, 5850, 6000), (9000, 8900, 9050)]
[(950, 909, 1059), (1500, 1499, 1649)]
```

### 34) `dnakit.annotation.split_interval_bins()`

- **作用：** 把半开区间 [start, end) 切成相邻的等宽 bin，最后一个 bin 可能较短。
- **处理方法：** 调用 genomicranges.utils 的 `split_intervals`（传入闭区间终点 end-1），并把上游的偏移量换算为半开区间；检查 bin 合计长度与原区间一致。
- **API：** `split_interval_bins(start, end, bin_size, *, verbose=True)`
- **输入：** 整数 `start`、`end`（开区间）与 `bin_size`。
- **输出：** `bins`（`index`、`start`、`end`）与 `bin_count`。
- **示例代码：**

```python
from dnakit.annotation import split_interval_bins

result = split_interval_bins(0, 100, 30)
print([(row["start"], row["end"]) for row in result["bins"]])
```

- **示例输出：**

```text
[(0, 30), (30, 60), (60, 90), (90, 100)]
```

### 35) `dnakit.annotation.slide_interval_windows()`

- **作用：** 在半开区间 [start, end) 上以 `step` 为步长滑动，生成长度为 `width` 的窗口。
- **处理方法：** 调用 genomicranges.utils 的 `slide_intervals`，并把上游的闭区间结果换算为半开区间；区间短于 `width` 时只返回一个窗口。
- **API：** `slide_interval_windows(start, end, width, step, *, verbose=True)`
- **输入：** 整数 `start`、`end`（开区间）、`width` 与 `step`。
- **输出：** `windows`（`index`、`start`、`end`）与 `window_count`。
- **示例代码：**

```python
from dnakit.annotation import slide_interval_windows

result = slide_interval_windows(0, 100, 30, 10)
print(result["window_count"], [(row["start"], row["end"]) for row in result["windows"]][:3])
```

- **示例输出：**

```text
8 [(0, 30), (10, 40), (20, 50)]
```

## 覆盖与跳过

下表汇总各包清单中每个公共函数的状态（完整记录见批次台账）。未包装的函数均给出了原因。

| 包 | 清单函数数 | 状态统计 |
| --- | --- | --- |
| polars_bio | 211 | duplicate-of 3, family-variant 17, internal 94, not-dna 10, other 68, plot 1, wrapped 18 |
| ncls | 5 | already-wrapped 1, internal 3, not-dna 1 |
| genomicranges | 22 | duplicate-of 3, family-variant 1, internal 8, needs-network 4, not-dna 2, wrapped 4 |
| gtars | 29 | family-variant 5, other 11, wrapped 13 |
