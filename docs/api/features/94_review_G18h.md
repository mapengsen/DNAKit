# 审阅组 G18h：pysam 与 hgvs 变异区辅助函数

本页记录 pysam（samtools 与 bcftools 子命令）和 hgvs 中未被覆盖的 DNA 级函数的离线封装。所有函数只读取本地输入，输出写入新路径。bed-reader、pandas-plink 和 bgen-reader 剩余的行都是内部辅助或测试入口，没有公开的基因型函数，因此不在本页封装。hgvs 的参考类函数需要 seqrepo 与 PostgreSQL，标记为 needs-data，不在本页。

模块：`dnakit.variants`、`dnakit.variants`。

## 示例文件

以下文件在运行示例前写入当前目录；每个示例都在新的临时目录中写出结果。

**示例文件：**
`reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:80
@RG	ID:rg1	SM:sampleA
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	0	chr1	5	60	10M	*	0	0	ACGTTTAGCC	IIIIIIIIII	RG:Z:rg1
r3	16	chr1	20	30	8M2I	*	0	0	GATCGATCGG	IIIIIIIIII	RG:Z:rg1
r4	0	chr2	3	60	10M	*	0	0	GGCCTTAAGG	IIIIIIIIII	RG:Z:rg1
r5	4	*	0	0	*	*	0	0	ACGTACGTAC	IIIIIIIIII	RG:Z:rg1
```

**示例文件：**
`ref.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
```

**示例文件：**
`primers.bed`

```text
chr1	0	10	amp1_LEFT	60	+
chr1	100	110	amp1_RIGHT	60	-
```

**示例文件：**
`intervals.bed`

```text
chr1	0	10
chr1	0	120
```

**示例文件：**
`genes.gff3`

```text
##gff-version 3
chr1	test	gene	1	120	.	+	.	ID=gene1;Name=G1
chr1	test	mRNA	1	120	.	+	.	ID=tx1;Parent=gene1;Name=T1
chr1	test	CDS	1	120	.	+	0	ID=cds1;Parent=tx1
```

**示例文件：**
`names.txt`

```text
A1
A2
```

**示例文件：**
`header.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:80
@RG	ID:new	SM:sampleB
```

**示例文件：**
`reads.fq`

```text
@r1
ACGTACGTAC
+
IIIIIIIIII
@r2
GGGGCCCCAA
+
IIIIIIIIII
```

**示例文件：**
`calls.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0/1:20	1/1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0/0:20	0/1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0/1:20	0/0:18
```

**示例文件：**
`calls2.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S3	S4
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0/1:20	1/1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0/0:20	0/1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0/1:20	0/0:18
```

**示例文件：**
`phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0|1:20	1|1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0|0:20	0|1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0|1:20	0|0:18
```

## samtools 比对文件函数

## 1) 设置读段组并标记每条读段

- **作用：** 为 BAM 写入 @RG 头行，并给每条读段加上 RG 标签。
- **计算方法：** 调用 `samtools addreplacerg`；`mode` 可为 overwrite_all 或 orphan_only。
- **API：** `dnakit.variants.add_bam_read_group(source, output_path, read_group, *, mode='overwrite_all', force=False, verbose=True)`
- **输入：** 比对文件 (SAM/BAM/CRAM)、@RG 头行、模式；输出路径必须是新文件。
- **输出：** 返回 read_group_id、mode、records 和输出路径。

**示例代码：**

```python
from dnakit.variants import add_bam_read_group

import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = add_bam_read_group(bam, work / "rg.bam", "@RG\tID:rgX\tSM:sampleB")
print(result["read_group_id"], result["records"])
```

**示例输出：**

```text
rgX 5
```

## 2) 软剪切扩增子引物

- **作用：** 去掉落在引物区间内的比对碱基，保留读段数。
- **计算方法：** 调用 `samtools ampliconclip -b BED`，区间为 0-based BED 坐标。
- **API：** `dnakit.variants.clip_amplicon_primers(source, primers_bed, output_path, *, force=False, verbose=True)`
- **输入：** 比对文件和引物 BED；输出为新文件。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import clip_amplicon_primers, select_alignment_records

work = Path(tempfile.mkdtemp())
shutil.copy("primers.bed", work / "primers.bed")
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = clip_amplicon_primers(bam, work / "primers.bed", work / "clipped.bam")
print(result["records"])
```

**示例输出：**

```text
5
```

## 3) 扩增子布局与统计

- **作用：** 读取扩增子位置和每条参考序列的读段统计。
- **计算方法：** 调用 `samtools ampliconstats`，解析 AMPLICON 与 FSS 行。
- **API：** `dnakit.variants.summarize_amplicon_stats(source, primers_bed, *, verbose=True)`
- **输入：** 比对文件和引物 BED。
- **输出：** 返回 amplicons、summary 和 count。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records, summarize_amplicon_stats

work = Path(tempfile.mkdtemp())
shutil.copy("primers.bed", work / "primers.bed")
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = summarize_amplicon_stats(bam, work / "primers.bed")
print(result["count"], result["amplicons"][0]["left"], result["amplicons"][0]["right"])
```

**示例输出：**

```text
1 1-10 101-110
```

## 4) 从 BAM 导出 FASTQ

- **作用：** 按读段角色把 BAM 写成四个 FASTQ 文件；序列和质量保持不变。
- **计算方法：** 调用 `samtools fastq -0 -1 -2 -s`，四个输出都显式指定，不会写到标准输出。
- **API：** `dnakit.variants.extract_fastq_from_bam(source, output_directory, *, verbose=True)`
- **输入：** 比对文件和新目录。
- **输出：** 返回各角色的记录数 records、文件 files 和目录 directory。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import extract_fastq_from_bam

result = extract_fastq_from_bam(bam, work / "fastq")
print(result["records"])
```

**示例输出：**

```text
{'other': 5, 'read1': 0, 'read2': 0, 'singletons': 0}
```

## 5) 从 BAM 导出 FASTA

- **作用：** 按读段角色把 BAM 序列写成四个 FASTA 文件。
- **计算方法：** 调用 `samtools fasta -0 -1 -2 -s`，输出路径全部显式指定。
- **API：** `dnakit.variants.extract_fasta_from_bam(source, output_directory, *, verbose=True)`
- **输入：** 比对文件和新目录。
- **输出：** 返回各角色的记录数、文件和目录。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import extract_fasta_from_bam

result = extract_fasta_from_bam(bam, work / "fasta")
print(result["records"])
```

**示例输出：**

```text
{'other': 5, 'read1': 0, 'read2': 0, 'singletons': 0}
```

## 6) 按读段名分组

- **作用：** 把同名读段排在一起，适合成对读段的后续处理。
- **计算方法：** 调用 `samtools collate -o`。
- **API：** `dnakit.variants.collate_bam_reads(source, output_path, *, force=False, verbose=True)`
- **输入：** 比对文件和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import collate_bam_reads

result = collate_bam_reads(bam, work / "collated.bam")
print(result["records"])
```

**示例输出：**

```text
5
```

## 7) 拼接 BAM 文件

- **作用：** 按输入顺序拼接多个具有相同头部的 BAM。
- **计算方法：** 调用 `samtools cat -o`。
- **API：** `dnakit.variants.concatenate_bam_files(sources, output_path, *, force=False, verbose=True)`
- **输入：** 两个或更多比对文件，输出为新文件。
- **输出：** 返回 records、inputs 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import concatenate_bam_files

result = concatenate_bam_files([bam, bam], work / "cat.bam")
print(result["records"])
```

**示例输出：**

```text
10
```

## 8) 统计 BED 区间内的覆盖碱基数

- **作用：** 对每个 BED 区间，统计所有读段在区间内对齐的碱基数之和。
- **计算方法：** 调用 `samtools bedcov`，需要 BAM 索引。同一碱基被多条读段覆盖时会重复计数。
- **API：** `dnakit.variants.count_bam_reads_per_interval(intervals_bed, source, *, verbose=True)`
- **输入：** BED 区间文件和比对文件。
- **输出：** 返回每个区间的 chrom、start、end 和 covered_bases。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.variants import count_bam_reads_per_interval, index_alignment_file

index_alignment_file(bam)
shutil.copy("intervals.bed", work / "intervals.bed")
result = count_bam_reads_per_interval(work / "intervals.bed", bam)
print([row["covered_bases"] for row in result["intervals"]])
```

**示例输出：**

```text
[16, 28]
```

## 9) 重算 NM 与 MD 标签

- **作用：** 用参考 FASTA 重新计算每条比对的 NM 与 MD 标签。
- **计算方法：** 调用 `samtools calmd -b`；参考需先有 .fai 索引。
- **API：** `dnakit.variants.recompute_nm_md_tags(source, reference, output_path, *, force=False, verbose=True)`
- **输入：** 比对文件、参考 FASTA 和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import recompute_nm_md_tags

shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
result = recompute_nm_md_tags(bam, work / "ref.fa", work / "calmd.bam")
print(result["records"])
```

**示例输出：**

```text
5
```

## 10) 按序列统计覆盖度

- **作用：** 每条序列的读段数、覆盖碱基数和平均深度。
- **计算方法：** 调用 `samtools coverage`，解析表格行。
- **API：** `dnakit.variants.summarize_bam_coverage(source, *, verbose=True)`
- **输入：** 比对文件。
- **输出：** 返回 contigs 列表与 count。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import summarize_bam_coverage

result = summarize_bam_coverage(bam)
print({row["rname"]: row["covbases"] for row in result["contigs"]})
```

**示例输出：**

```text
{'chr1': 22, 'chr2': 10}
```

## 11) 逐碱基读段深度

- **作用：** 输出每个位置的读段深度；可指定区域。
- **计算方法：** 调用 `samtools depth -a`；区域查询需要索引。
- **API：** `dnakit.variants.calculate_bam_read_depth(source, *, region=None, all_positions=True, max_positions=100000, verbose=True)`
- **输入：** 比对文件、可选区域与最大行数。
- **输出：** 返回 positions（contig、position、depth）和 count。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import calculate_bam_read_depth, index_alignment_file

index_alignment_file(bam)
result = calculate_bam_read_depth(bam, region="chr1:1-10")
print([row["depth"] for row in result["positions"]])
```

**示例输出：**

```text
[1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
```

## 12) 写出 FASTA 序列字典

- **作用：** 写出参考的 .dict 文件，并列出每条序列的名称、长度和 MD5。
- **计算方法：** 调用 `samtools dict`。
- **API：** `dnakit.variants.write_fasta_sequence_dictionary(reference, output_path, *, force=False, verbose=True)`
- **输入：** 参考 FASTA 和以 .dict 结尾的新输出路径。
- **输出：** 返回 contigs（name、length、md5）和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import write_fasta_sequence_dictionary

work = Path(tempfile.mkdtemp())
shutil.copy("ref.fa", work / "ref.fa")
result = write_fasta_sequence_dictionary(work / "ref.fa", work / "ref.dict")
print([(item["name"], item["length"]) for item in result["contigs"]])
```

**示例输出：**

```text
[('chr1', 120), ('chr2', 80)]
```

## 13) 解码 SAM FLAG

- **作用：** 把 FLAG 整数转换为位名称，如 99 为 PAIRED、PROPER_PAIR 等。
- **计算方法：** 调用 `samtools flags`，结果按位名称排序。
- **API：** `dnakit.variants.decode_sam_flags(flags, *, verbose=True)`
- **输入：** 0 到 65535 之间的整数列表。
- **输出：** 返回每个 FLAG 的 flag、hex 和 names。

**示例代码：**

```python
from dnakit.variants import decode_sam_flags

result = decode_sam_flags([99, 4])
print([(item["flag"], item["names"]) for item in result["flags"]])
```

**示例输出：**

```text
[(99, ['MREVERSE', 'PAIRED', 'PROPER_PAIR', 'READ1']), (4, ['UNMAP'])]
```

## 14) 读取头部与前若干条记录

- **作用：** 返回比对文件的头部行和前 N 条记录。
- **计算方法：** 调用 `samtools head -n`。
- **API：** `dnakit.variants.read_bam_header_and_records(source, *, records=5, verbose=True)`
- **输入：** 比对文件和记录数（0 到 1000）。
- **输出：** 返回 header 与 records（SAM 列）。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import read_bam_header_and_records

result = read_bam_header_and_records(bam, records=2)
print([row[0] for row in result["records"]], result["header"][0])
```

**示例输出：**

```text
['r1', 'r2'] @HD	VN:1.6	SO:coordinate
```

## 15) 为 BAM 或 CRAM 建索引

- **作用：** 写出 BAI（或 CSI）索引，供区域查询使用。
- **计算方法：** 调用 `samtools index`，索引写在输入文件旁边。
- **API：** `dnakit.variants.index_alignment_file(source, *, csi=False, force=False, verbose=True)`
- **输入：** 坐标排序的比对文件。
- **输出：** 返回索引路径和格式。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file

print(index_alignment_file(bam)["format"], index_alignment_file(bam, csi=True)["format"])
```

**示例输出：**

```text
BAI CSI
```

## 16) 按序列统计已比对与未比对读段

- **作用：** 从索引读取每条序列的已比对和未比对读段数。
- **计算方法：** 调用 `samtools idxstats`；需要先有索引。
- **API：** `dnakit.variants.summarize_bam_index_stats(source, *, verbose=True)`
- **输入：** 已建索引的比对文件。
- **输出：** 返回 contigs、unplaced_unmapped 和 provenance。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file, summarize_bam_index_stats

index_alignment_file(bam)
result = summarize_bam_index_stats(bam)
print([(row["contig"], row["mapped"]) for row in result["contigs"]], result["unplaced_unmapped"])
```

**示例输出：**

```text
[('chr1', 3), ('chr2', 1)] 1
```

## 17) 标记重复读段

- **作用：** 给重复读段加 0x400 标志，不删除读段。
- **计算方法：** 调用 `samtools markdup`，输入需按坐标排序。
- **API：** `dnakit.variants.mark_bam_duplicates(source, output_path, *, force=False, verbose=True)`
- **输入：** 坐标排序的比对文件和新输出路径。
- **输出：** 返回 records、duplicates 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import mark_bam_duplicates

result = mark_bam_duplicates(bam, work / "marked.bam")
print(result["records"], result["duplicates"])
```

**示例输出：**

```text
5 0
```

## 18) 合并 BAM 文件

- **作用：** 合并多个坐标排序的 BAM，输出仍按坐标排序。
- **计算方法：** 调用 `samtools merge`。
- **API：** `dnakit.variants.merge_bam_files(sources, output_path, *, force=False, verbose=True)`
- **输入：** 两个或更多比对文件和新输出路径。
- **输出：** 返回 records、inputs 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import merge_bam_files

result = merge_bam_files([bam, bam], work / "merged.bam")
print(result["records"])
```

**示例输出：**

```text
10
```

## 19) 逐位置碱基堆叠文本

- **作用：** 输出每个被覆盖位置的碱基堆叠文本，用于后续变异检测。
- **计算方法：** 调用 `samtools mpileup`；给出参考时需要 .fai 索引。
- **API：** `dnakit.variants.pileup_alignment_positions(source, output_path, *, reference=None, force=False, verbose=True)`
- **输入：** 比对文件、可选参考 FASTA 和新输出路径。
- **输出：** 返回 positions（被覆盖的位置数）和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import select_alignment_records
from dnakit.variants import pileup_alignment_positions

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
result = pileup_alignment_positions(bam, work / "pileup.txt", reference=work / "ref.fa")
print(result["positions"])
```

**示例输出：**

```text
32
```

## 20) 补全成对读段的配对字段

- **作用：** 为按名称排序的成对读段补写配对位置、插入长度等标签。
- **计算方法：** 调用 `samtools fixmate -m`；输入必须按读段名排序，否则报错。
- **API：** `dnakit.variants.fix_bam_mate_fields(source, output_path, *, force=False, verbose=True)`
- **输入：** 按名称排序的比对文件和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import fix_bam_mate_fields, sort_alignment_file

named = sort_alignment_file(bam, work / "byname.bam", by_name=True)["output"]
print(fix_bam_mate_fields(named, work / "fixed.bam")["records"])
```

**示例输出：**

```text
5
```

## 21) 检查比对文件是否有效

- **作用：** 逐个检查 BAM/CRAM 能否打开、头部是否有效。
- **计算方法：** 调用 `samtools quickcheck`，每个文件单独判断。
- **API：** `dnakit.variants.check_alignment_file_integrity(sources, *, verbose=True)`
- **输入：** 本地文件列表。
- **输出：** 返回每个文件的 valid 状态及有效和无效数量。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from pathlib import Path
from dnakit.variants import check_alignment_file_integrity

(work / "junk.bam").write_text("not a bam", encoding="utf-8")
result = check_alignment_file_integrity([bam, work / "junk.bam"])
print([(Path(item["path"]).name, item["valid"]) for item in result["files"]])
```

**示例输出：**

```text
[('reads.bam', True), ('junk.bam', False)]
```

## 22) 替换 BAM 头部

- **作用：** 只替换头部，读段记录原样复制。
- **计算方法：** 调用 `samtools reheader`，结果写入新文件。
- **API：** `dnakit.variants.replace_bam_header(header_sam, source, output_path, *, force=False, verbose=True)`
- **输入：** 以 @ 开头的 SAM 头部文件、比对文件和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.variants import replace_bam_header

shutil.copy("header.sam", work / "header.sam")
print(replace_bam_header(work / "header.sam", bam, work / "rh.bam")["records"])
```

**示例输出：**

```text
5
```

## 23) 清除比对信息

- **作用：** 清除位置、CIGAR 与比对标志，使读段变为未比对；序列、质量和读段组标签保留。
- **计算方法：** 调用 `samtools reset -o`。
- **API：** `dnakit.variants.reset_alignments_to_unmapped(source, output_path, *, force=False, verbose=True)`
- **输入：** 比对文件和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import reset_alignments_to_unmapped

print(reset_alignments_to_unmapped(bam, work / "reset.bam")["records"])
```

**示例输出：**

```text
5
```

## 24) 列出读段组样本名

- **作用：** 列出每个比对文件 @RG 头行中的 SM 样本名。
- **计算方法：** 调用 `samtools samples`。
- **API：** `dnakit.variants.list_read_group_samples(sources, *, verbose=True)`
- **输入：** 一个或多个比对文件。
- **输出：** 返回 samples（sample 与 source）列表。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import list_read_group_samples

print([item["sample"] for item in list_read_group_samples([bam])["samples"]])
```

**示例输出：**

```text
['sampleA']
```

## 25) 排序比对文件

- **作用：** 按坐标排序，或用 by_name=True 按读段名排序。
- **计算方法：** 调用 `samtools sort -o`，可选 `-n`。
- **API：** `dnakit.variants.sort_alignment_file(source, output_path, *, by_name=False, force=False, verbose=True)`
- **输入：** 比对文件和新输出路径。
- **输出：** 返回 order、records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import sort_alignment_file

print(sort_alignment_file(bam, work / "coord.bam")["order"], sort_alignment_file(bam, work / "name.bam", by_name=True)["order"])
```

**示例输出：**

```text
coordinate queryname
```

## 26) 按读段组拆分 BAM

- **作用：** 每个读段组写一个 BAM；所有读段都必须带 RG 标签。
- **计算方法：** 调用 `samtools split -f`，文件名为读段组 ID。
- **API：** `dnakit.variants.split_bam_by_read_group(source, output_directory, *, verbose=True)`
- **输入：** 比对文件和新目录。
- **输出：** 返回每个文件的名称与记录数。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import split_bam_by_read_group

result = split_bam_by_read_group(bam, work / "split")
print([(item["file"], item["records"]) for item in result["files"]])
```

**示例输出：**

```text
[('rg1.bam', 5)]
```

## 27) 比对统计摘要

- **作用：** 从 samtools stats 的 SN 行取出读段数、比对数和长度等汇总。
- **计算方法：** 调用 `samtools stats` 并解析 SN 行。
- **API：** `dnakit.variants.summarize_alignment_stats(source, *, verbose=True)`
- **输入：** 比对文件。
- **输出：** 返回 summary 字典。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import summarize_alignment_stats

print(summarize_alignment_stats(bam)["summary"]["raw total sequences"])
```

**示例输出：**

```text
5
```

## 28) 筛选比对记录

- **作用：** 按区域、FLAG 掩码和最低比对质量保留读段。
- **计算方法：** 调用 `samtools view`；区域需要索引。输出以 .sam 结尾时写文本，否则写 BAM。
- **API：** `dnakit.variants.select_alignment_records(source, output_path, *, region=None, required_flags=0, excluded_flags=0, min_mapping_quality=0, force=False, verbose=True)`
- **输入：** 比对文件、可选区域、FLAG 掩码和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file, select_alignment_records

index_alignment_file(bam)
print(select_alignment_records(bam, work / "sel.bam", region="chr1:1-10", min_mapping_quality=30)["records"])
```

**示例输出：**

```text
2
```

## 29) 由比对重建参考碱基

- **作用：** 用带 MD 标签的比对重建参考序列；无读段覆盖的位置为 N。
- **计算方法：** 调用 `samtools reference`；没有 MD 标签时结果几乎全为 N。
- **API：** `dnakit.variants.reconstruct_reference_from_alignments(source, output_path, *, force=False, verbose=True)`
- **输入：** 带 MD 标签的比对文件和新 FASTA 路径。
- **输出：** 返回每条序列的长度与已确定碱基数。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import recompute_nm_md_tags, reconstruct_reference_from_alignments, select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
tagged = recompute_nm_md_tags(bam, work / "ref.fa", work / "calmd.bam")["output"]
result = reconstruct_reference_from_alignments(tagged, work / "rebuilt.fa")
print({item["name"]: item["determined_bases"] for item in result["contigs"]})
```

**示例输出：**

```text
{'chr1': 22, 'chr2': 10}
```

## 30) 按单倍型拆分读段

- **作用：** 给读段加 HP 标签，并按单倍型写入不同文件。
- **计算方法：** 调用 `samtools phase -k -b`；输出为 prefix.0.bam、prefix.1.bam 和 prefix.chimera.bam。
- **API：** `dnakit.variants.phase_bam_haplotypes(source, output_prefix, *, block_length=13, verbose=True)`
- **输入：** 比对文件和输出前缀；同一前缀下已有文件时拒绝覆盖。
- **输出：** 返回每个输出文件的名称与记录数。

**示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import phase_bam_haplotypes

result = phase_bam_haplotypes(bam, work / "phased")
print(sorted(item["file"].rsplit("/", 1)[-1] for item in result["files"]))
```

**示例输出：**

```text
['phased.0.bam', 'phased.1.bam', 'phased.chimera.bam']
```

## 31) 为 FASTQ 建索引

- **作用：** 写出每条记录一行的 .fai 索引。
- **计算方法：** 调用 `samtools fqidx`，索引写在 FASTQ 旁边。
- **API：** `dnakit.variants.index_fastq_file(source, *, force=False, verbose=True)`
- **输入：** FASTQ 文件。
- **输出：** 返回索引路径和条目数。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_fastq_file

work = Path(tempfile.mkdtemp())
shutil.copy("reads.fq", work / "reads.fq")
print(index_fastq_file(work / "reads.fq")["entries"])
```

**示例输出：**

```text
2
```

## 32) FASTQ 导入为未比对 BAM

- **作用：** 把一个（单端）或两个（双端）FASTQ 转为未比对 BAM。
- **计算方法：** 调用 `samtools import -o`。
- **API：** `dnakit.variants.import_fastq_to_bam(sources, output_path, *, force=False, verbose=True)`
- **输入：** 一个或两个 FASTQ 文件和新输出路径。
- **输出：** 返回 records、layout 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import import_fastq_to_bam

work = Path(tempfile.mkdtemp())
shutil.copy("reads.fq", work / "reads.fq")
result = import_fastq_to_bam([work / "reads.fq"], work / "imported.bam")
print(result["records"], result["layout"])
```

**示例输出：**

```text
2 single
```

## bcftools 变异文件函数

## 33) 按表达式筛选 VCF 位点

- **作用：** 保留（include）或去除（exclude）满足 bcftools 表达式的位点。
- **计算方法：** 调用 `bcftools filter -i` 或 `-e`；表达式作为单个参数传递，不经过 shell。
- **API：** `dnakit.variants.filter_vcf_by_expression(source, expression, output_path, *, mode='include', force=False, verbose=True)`
- **输入：** bgzip 并建索引的 VCF、表达式和新输出路径。
- **输出：** 返回 records（保留的位点数）和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import filter_vcf_by_expression

print(filter_vcf_by_expression(gz, "QUAL>45", work / "kept.vcf")["records"])
```

**示例输出：**

```text
2
```

## 34) 注释 VCF 的蛋白质后果

- **作用：** 依据 GFF3 基因模型和参考 FASTA，为每个位点写入 BCSQ 标签。
- **计算方法：** 调用 `bcftools csq -g -f`。输入需已分相，否则 bcftools 会在未分相杂合位点上停止。
- **API：** `dnakit.variants.annotate_vcf_consequence(source, gene_annotation, reference, output_path, *, force=False, verbose=True)`
- **输入：** 已分相的 VCF、GFF3、参考 FASTA 和新输出路径。
- **输出：** 返回 records、annotated_records 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import annotate_vcf_consequence

work = Path(tempfile.mkdtemp())
for name in ("phased.vcf", "genes.gff3", "ref.fa"):
    shutil.copy(name, work / name)
result = annotate_vcf_consequence(work / "phased.vcf", work / "genes.gff3", work / "ref.fa", work / "csq.vcf")
print(result["annotated_records"])
```

**示例输出：**

```text
3
```

## 35) 检查样本间基因型一致性

- **作用：** 比较每对样本的基因型不一致率；默认与自身比较。
- **计算方法：** 调用 `bcftools gtcheck -g`。
- **API：** `dnakit.variants.check_vcf_sample_concordance(source, *, genotype_source=None, verbose=True)`
- **输入：** VCF，可选第二个基因型来源。
- **输出：** 返回 pairs（query、genotyped、discordance、sites、matching）。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import check_vcf_sample_concordance

result = check_vcf_sample_concordance(gz)
print([(item["query"], item["genotyped"], item["discordance"]) for item in result["pairs"]])
```

**示例输出：**

```text
[('S1', 'S1', 0.0), ('S1', 'S2', 27.63102), ('S2', 'S1', 27.63102), ('S2', 'S2', 0.0)]
```

## 36) 读取 VCF 头部

- **作用：** 返回 VCF 或 BCF 的头部行和样本名。
- **计算方法：** 调用 `bcftools head`。
- **API：** `dnakit.variants.read_vcf_header(source, *, verbose=True)`
- **输入：** VCF 或 BCF 文件。
- **输出：** 返回 header_lines 和 samples。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import read_vcf_header

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
print(read_vcf_header(work / "calls.vcf")["samples"])
```

**示例输出：**

```text
['S1', 'S2']
```

## 37) 为 VCF 建 CSI 索引

- **作用：** 为 bgzip 压缩的 VCF 或 BCF 写出 CSI 索引。
- **计算方法：** 调用 `bcftools index -c`。
- **API：** `dnakit.variants.index_vcf_csi(source, *, force=False, verbose=True)`
- **输入：** 以 .vcf.gz 或 .bcf 结尾的文件。
- **输出：** 返回索引路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import index_vcf_csi

print(Path(index_vcf_csi(gz)["index"]).name)
```

**示例输出：**

```text
calls.vcf.gz.csi
```

## 38) 合并已索引 VCF 的样本

- **作用：** 把多个来源不同样本的 VCF 合并为一个多样本 VCF。
- **计算方法：** 调用 `bcftools merge`；每个输入都需要 .tbi 或 .csi 索引。
- **API：** `dnakit.variants.merge_indexed_vcfs(sources, output_path, *, force=False, verbose=True)`
- **输入：** 两个或更多 bgzip VCF 和新输出路径。
- **输出：** 返回 samples、records 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file
from dnakit.variants import merge_indexed_vcfs

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
shutil.copy("calls2.vcf", work / "calls2.vcf")
first = index_tabix_file(work / "calls.vcf")["compressed"]
second = index_tabix_file(work / "calls2.vcf")["compressed"]
print(merge_indexed_vcfs([first, second], work / "merged.vcf")["samples"])
```

**示例输出：**

```text
['S1', 'S2', 'S3', 'S4']
```

## 39) BAM 生成 VCF 似然

- **作用：** 输出 BAM 在每个位点的基因型似然（PL），调用前的中间结果。
- **计算方法：** 调用 `bcftools mpileup -f -Ov`；参考需要 .fai 索引。
- **API：** `dnakit.variants.pileup_bam_to_vcf(source, reference, output_path, *, force=False, verbose=True)`
- **输入：** 比对文件、参考 FASTA 和新输出路径。
- **输出：** 返回 records 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import pileup_bam_to_vcf, select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
print(pileup_bam_to_vcf(bam, work / "ref.fa", work / "pileup.vcf")["records"])
```

**示例输出：**

```text
32
```

## 40) 重命名 VCF 样本

- **作用：** 用文件中的新名称替换样本名，其他内容不变。
- **计算方法：** 调用 `bcftools reheader -s`。
- **API：** `dnakit.variants.rename_vcf_samples(source, sample_names, output_path, *, force=False, verbose=True)`
- **输入：** VCF、每行一个名称的文件（名称不得重复）和新输出路径。
- **输出：** 返回 samples 和输出路径。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import rename_vcf_samples

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
shutil.copy("names.txt", work / "names.txt")
print(rename_vcf_samples(work / "calls.vcf", work / "names.txt", work / "renamed.vcf")["samples"])
```

**示例输出：**

```text
['A1', 'A2']
```

## 41) 用 bcftools 检测纯合区段

- **作用：** 用 bcftools 的 HMM 检测一个样本的纯合区段（ROH）。
- **计算方法：** 调用 `bcftools roh -G`，只需要 GT 字段。小数据集可能没有区段。
- **API：** `dnakit.variants.detect_bcftools_roh(source, *, sample, genotype_quality=30.0, verbose=True)`
- **输入：** bgzip VCF、样本名与 PL 替代值（0 到 1000）。
- **输出：** 返回 runs（每段的 chromosome、start、end、length、markers、quality）和 count。

**示例代码：**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file
from dnakit.variants import detect_bcftools_roh

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
result = detect_bcftools_roh(gz, sample="S1")
print(result["sample"], result["count"])
```

**示例输出：**

```text
S1 0
```

## hgvs 辅助函数

## 42) 去除等位基因的公共后缀

- **作用：** 去掉等位基因字符串共有的后缀，并报告去掉的长度。
- **计算方法：** 调用 hgvs 的 `trim_common_suffixes`。
- **API：** `dnakit.variants.trim_common_suffix_strings(strings, *, min_length=0, verbose=True)`
- **输入：** 至少两个由 A、C、G、T、N 组成的等位基因字符串。
- **输出：** 返回 trimmed_length 和 strings。

**示例代码：**

```python
from dnakit.variants import trim_common_suffix_strings

result = trim_common_suffix_strings(["GGCTA", "AGCTA"])
print(result["trimmed_length"], result["strings"])
```

**示例输出：**

```text
4 ['G', 'A']
```

## 43) HGVS 变异的 0 基半开区间

- **作用：** 解析一条 HGVS 字符串（不需要参考数据），返回 0 基半开（interbase）起止坐标。
- **计算方法：** 用 hgvs 解析器读取位置，再调用 `get_start_end_interbase`。
- **API：** `dnakit.variants.hgvs_interbase_span(variant, *, outer_confidence=True, verbose=True)`
- **输入：** 一条 HGVS 变异字符串，如 NC_000001.11:g.100_102del。
- **输出：** 返回 start 与 end（0 基、半开）。

**示例代码：**

```python
from dnakit.variants import hgvs_interbase_span

result = hgvs_interbase_span("NC_000001.11:g.100_102del")
print(result["start"], result["end"])
```

**示例输出：**

```text
99 102
```

## 44) 等位基因的重复序列记号

- **作用：** 根据重复块，为一条等位基因写出 HGVS 风格的重复记号，如 AC[3]。
- **计算方法：** 调用 hgvs 的 `get_repeat_str`；没有记号时结果为 None。
- **API：** `dnakit.variants.repeat_notation_for_alleles(sequence, other_sequence, blocks, other_blocks=(), *, reverse=False, verbose=True)`
- **输入：** 两条序列及其重复块（detect_repeat_blocks 的输出）。
- **输出：** 返回 formatted 字符串或 None。

**示例代码：**

```python
from dnakit.variants import detect_repeat_blocks
from dnakit.variants import repeat_notation_for_alleles

blocks = detect_repeat_blocks("ACACAC")["blocks"]
print(repeat_notation_for_alleles("ACACAC", "ACAC", blocks)["formatted"])
```

**示例输出：**

```text
AC[3]
```
