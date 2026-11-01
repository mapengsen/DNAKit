# 读段工具：pysam、pyfaidx、cutadapt 与 UMI-tools

这一组函数封装了四个处理测序读段和 FASTA 文件的包：[pysam](https://github.com/pysam-developers/pysam)（连同它自带的 samtools 与 bcftools 代码）、[pyfaidx](https://github.com/mdshw5/pyfaidx)、[cutadapt](https://github.com/marcelm/cutadapt) 和 [UMI-tools](https://github.com/CGATOxford/UMI-tools)。pysam 部分的函数可转换质量字符串、读取 tabix 风格文件、运行 samtools 与 bcftools 命令（去填充、导出 FASTQ、整理读段顺序、CRAM 大小、拷贝数推断），并把比对渲染成文本视图。pyfaidx 部分的函数解析区域与 BED 行，并提取、屏蔽、制表和拆分 FASTA 序列。cutadapt 部分的函数对 FASTQ 读段修剪接头、低质量末端、固定数量的碱基和 N 碱基，并提供多种报告格式。UMI-tools 部分的函数把 UMI 移到读段名中、估计细胞白名单、按 UMI 对 BAM 分组与去重、按基因统计 UMI，并为 RSEM 准备 BAM。从 `dnakit.alignment`、`dnakit.sequence`、`dnakit.fasta_fastq` 和 `dnakit.sequencing_qc` 导入，例如 `from dnakit.alignment import qualitystring_to_array`。

已安装各包的许可证（取自其包元数据）：pysam、cutadapt 和 UMI-tools 为 MIT，pyfaidx 为 BSD-3-Clause；均不是 GPL。

依赖都是可选的，只有调用对应函数时才会导入。缺少包或程序时会报错并给出安装提示。示例对每个函数都传入 `verbose=False`，把输出写入新的临时目录，输入则来自下面声明的示例文件。

## 示例数据

示例使用下面声明的小文件：一个 SAM 文件（`reads.sam`）及其参考序列（`ref.fa`）、一个 BED 文件、一个填充比对及其填充参考序列、带 `BAF` 与 `LRR` 值的 VCF、一个 FASTA 文件（`genes.fa`）、两个 FASTQ 文件以及 UMI 用的 SAM 文件。需要 BAM 或 CRAM 的示例会在代码里用 `pysam.view()` 由 SAM 文件生成。

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:40
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTG	IIIIIHHHHH
r2	16	chr1	12	60	8M	*	0	0	TTGACCGA	FFFFFFFF
r3	0	chr1	20	60	6M	*	0	0	GGGCCC	IIIIII
```

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTACGTTGACCGAGGGCCCAAATTTGGGCCCAAA
```

**示例文件：** `regions.bed`

```text
g1	2	10
g2	0	8
```

**示例文件：** `padded.fa`

```text
>chr1
ACGT**ACGTAC
```

**示例文件：** `pad.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:12
r1	0	chr1	1	60	4M2D6M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	2M2D4M	*	0	0	GTCGTA	IIIIII
```

**示例文件：** `cnv.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele frequency">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	2000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	3000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	4000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	5000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	6000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	7000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	8000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	9000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	10000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	11000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	12000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	13000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	14000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	15000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	16000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	17000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	18000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	19000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	20000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	21000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	22000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	23000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	24000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	25000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	26000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	27000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	28000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	29000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	30000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	31000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	32000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	33000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	34000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	35000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	36000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
```

**示例文件：** `genes.fa`

```text
>g1
ACGTACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
>g3
GATTACA
```

**示例文件：** `reads.fq`

```text
@r1
ACGTTGCAAGATCGGAAGAGCACACG
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@r2
TTGACCGATCAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIII
@r3
GGGCCCAAATTTGGGAAA
+
IIIIIIIIIIIIIIIIII
@r4
NNACGTACGTTGAACGNN
+
##IIIIIIIIIIIIII##
```

**示例文件：** `umi_reads.fq`

```text
@r1
ACGTAAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r2
TTGCCAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r3
ACGTAAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r4
ACGTTTGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
```

**示例文件：** `umi.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:100
a1_AAAA	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a2_AAAA	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a3_AAAT	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a4_CCCC	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a5_GGGG	0	chr1	30	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneB
a6_GGGG	0	chr1	30	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneB
```

**示例文件：** `pair.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:tx1	LN:100
p1_AAAA	99	tx1	5	60	10M	=	40	45	ACGTACGTAC	IIIIIIIIII	RX:Z:AAAA
p1_AAAA	147	tx1	40	60	10M	=	5	-45	ACGTACGTAC	IIIIIIIIII
p2_CCCC	99	tx1	10	60	10M	=	50	45	ACGTACGTAC	IIIIIIIIII	RX:Z:CCCC
p2_CCCC	147	tx1	50	60	10M	=	10	-45	ACGTACGTAC	IIIIIIIIII
```

## 1) `dnakit.alignment.array_to_qualitystring()`

- **作用：** 把整数 Phred 质量值转换为质量字符串。
- **处理方法：** 把整数放进无符号字节数组，调用 `pysam.array_to_qualitystring()`，它给每个值加上 `offset` 后返回对应字符。
- **API：** `array_to_qualitystring(qualities: 'Iterable[int]', *, offset: 'int' = 33, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `qualities` 是整数可迭代对象（逐个用 `int()` 转换并存为无符号字节，范围 0 到 255）。`offset` 对 Sanger/Illumina 1.8 为 33，对旧版 Illumina 编码为 64。
- **输出：** 字符串，每个质量值对应一个字符。
- **示例代码：**

```python
from dnakit.alignment import array_to_qualitystring

print(array_to_qualitystring([40, 40, 20, 2], verbose=False))
print(array_to_qualitystring([40, 40, 20, 2], offset=64, verbose=False))
```

- **示例输出：**

```text
II5#
hhTB
```

## 2) `dnakit.alignment.bam_to_fastq()`

- **作用：** 用 samtools `fastq` 命令把比对记录转换为 FASTQ 文本。
- **处理方法：** 用你给的参数再加上 `-0 <output_path>` 调用 `pysam.bam2fq()`。pysam 会把输出写到进程的标准输出，所以函数让 samtools 把结果写入 `output_path`，再把该文件读回。
- **API：** `bam_to_fastq(arguments: 'Sequence[str]', output_path: 'str', *, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是传给 `samtools fastq` 的字符串列表（传入字符串或含非字符串项会抛出 `ConfigurationError`），至少要给出输入的 BAM、SAM 或 CRAM 路径。`output_path` 必须尚不存在。在 samtools 中，`-0` 接收既没有 READ1 也没有 READ2 标志的读段，例如示例中的单端读段。
- **输出：** FASTQ 文本字符串。同样的文本也写入新文件 `output_path`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import bam_to_fastq

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
text = bam_to_fastq(["-n", str(bam)], str(work / "reads.fq"), verbose=False)
print(text)
```

- **示例输出：**

```text
@r1
ACGTACGTTG
+
IIIIIHHHHH
@r2
TCGGTCAA
+
FFFFFFFF
@r3
GGGCCC
+
IIIIII
```

## 3) `dnakit.alignment.cram_container_sizes()`

- **作用：** 报告 CRAM 文件的容器与数据块大小。
- **处理方法：** 用你给的参数调用 `pysam.cram_size()`（即 samtools 的 `cram_size` 命令）并返回其文本。
- **API：** `cram_container_sizes(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是传给 `samtools cram_size` 的字符串列表，通常是 CRAM 路径，可加 `-v` 等选项。非 CRAM 文件会抛出 pysam 的 `SamtoolsError`（"Input is not a CRAM file"）。
- **输出：** 文本报告：先是数据块表（内容 id、未压缩与压缩大小、比率、方法、数据序列），然后是 `Number of containers`、`Number of slices`、`Number of sequences`、`Number of bases`、`Total file size` 和 `Format overhead size` 几行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import cram_container_sizes

work = Path(tempfile.mkdtemp())
cram = work / "reads.cram"
pysam.view("-C", "-T", "ref.fa", "-o", str(cram), "reads.sam", catch_stdout=False)
report = cram_container_sizes([str(cram)], verbose=False)
for line in report.splitlines():
    if line.startswith("Number of"):
        print(line)
```

- **示例输出：**

```text
Number of containers                   1
Number of slices                       1
Number of sequences                    3
Number of bases                       24
```

## 4) `dnakit.alignment.pad_to_unpadded()`

- **作用：** 用旧版 samtools 命令名 `pad2unpad` 和你自己的参数把填充（padded）比对转换为非填充坐标。
- **处理方法：** 用你的参数调用 `pysam.pad2unpad()`。这个 pysam 包装不捕获标准输出：不带 `-o` 时比对结果直接写到进程标准输出，函数返回空字符串；带 `-o <文件>` 时 samtools 写出该文件，函数同样返回空字符串。
- **API：** `pad_to_unpadded(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是传给 `samtools pad2unpad` 的字符串列表：用 `-T <fasta>` 给出填充后的参考序列，`-s` 要求输出 SAM，`-o <文件>` 指定输出，最后是输入比对文件。读段不能含 P 或 I 的 CIGAR 操作。
- **输出：** 通常是空字符串，结果要从 `-o` 文件读取。若想直接得到数据或结构化结果，请用 `remove_padding()` 或 `pad_to_unpad_alignments()`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pad_to_unpadded

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
out = work / "unpadded.sam"
returned = pad_to_unpadded(["-T", "padded.fa", "-s", "-o", str(out), str(bam)], verbose=False)
print(repr(returned))
for line in out.read_text().splitlines():
    if not line.startswith("@"):
        print(line)
```

- **示例输出：**

```text
''
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	6M	*	0	0	GTCGTA	IIIIII
```

## 5) `dnakit.alignment.qualities_to_qualitystring()`

- **作用：** 把质量分数列表转换为质量字符串。
- **处理方法：** 把整数放进无符号字节数组并调用 `pysam.qualities_to_qualitystring()`。对整数输入，结果与 `array_to_qualitystring()` 相同。
- **API：** `qualities_to_qualitystring(qualities: 'Iterable[int]', *, offset: 'int' = 33, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `qualities` 是整数可迭代对象（列表、`range` 等都可以，逐个用 `int()` 转换并存为无符号字节）。`offset` 默认为 33。
- **输出：** 字符串，每个质量值对应一个字符。
- **示例代码：**

```python
from dnakit.alignment import qualities_to_qualitystring

print(qualities_to_qualitystring([40, 40, 20, 2], verbose=False))
print(qualities_to_qualitystring(range(0, 41, 10), offset=64, verbose=False))
```

- **示例输出：**

```text
II5#
@JT^h
```

## 6) `dnakit.alignment.qualitystring_to_array()`

- **作用：** 把 Phred 质量字符串转换为整数质量值列表。
- **处理方法：** 调用 `pysam.qualitystring_to_array()`，它把每个字符的编码减去 `offset`，再把结果转换为 Python 整数列表。
- **API：** `qualitystring_to_array(quality: 'str', *, offset: 'int' = 33, verbose: bool = True) -> 'list[int]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `quality` 必须是字符串（否则抛出 `ConfigurationError`）。`offset` 默认为 33；旧版 Illumina 编码用 64。
- **输出：** 整数列表，每个字符对应一个整数。
- **示例代码：**

```python
from dnakit.alignment import qualitystring_to_array

print(qualitystring_to_array("II5#", verbose=False))
print(qualitystring_to_array("hhTB", offset=64, verbose=False))
```

- **示例输出：**

```text
[40, 40, 20, 2]
[40, 40, 20, 2]
```

## 7) `dnakit.alignment.remove_padding()`

- **作用：** 用你自己的 samtools 参数去掉填充比对中的填充（padding）。
- **处理方法：** 用你的参数调用 `pysam.depad()`（即 samtools 的 `depad` 命令），返回 pysam 捕获的标准输出。pysam 自己会加一个临时的 `-o` 选项，所以返回的是完整 SAM 文本，其中的 `@PG` 头行含有临时文件名。
- **API：** `remove_padding(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是传给 `samtools depad` 的字符串列表：`-T <fasta>` 给出填充后的参考序列（填充列写作 `*`），`-s` 要求输出 SAM 文本，最后是输入比对文件。不要自己加 `-o`，捕获的文本就是结果。读段不能含 P 或 I 的 CIGAR 操作。
- **输出：** 加 `-s` 时是字符串形式的 SAM 文本。不加 `-s` 时 samtools 输出 BAM，不是文本。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import remove_padding

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
text = remove_padding(["-T", "padded.fa", "-s", str(bam)], verbose=False)
for line in text.splitlines():
    if not line.startswith("@"):
        print(line)
```

- **示例输出：**

```text
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	6M	*	0	0	GTCGTA	IIIIII
```

## 8) `dnakit.alignment.shuffle_bam()`

- **作用：** 用你自己的参数运行 samtools 的 `bamshuf` 命令。
- **处理方法：** 用你的参数调用 `pysam.bamshuf()`，返回 pysam 捕获的标准输出。在随包的 samtools 1.21 中 `bamshuf` 就是 `collate` 命令（它打印的用法文本是 `samtools collate` 的），作用是把同名读段排到一起。
- **API：** `shuffle_bam(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是字符串列表：选项如 `-o <文件>`（输出文件）、`-n <整数>`（临时文件数）、`-T <前缀>`（临时文件前缀），然后是输入 BAM；不用 `-o` 时再加一个输出前缀。请给 `-o` 或前缀，让结果写入文件；`-O` 会把 BAM 写到进程标准输出，不会被返回。
- **输出：** 捕获到的标准输出；结果写入文件时为空字符串。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import shuffle_bam

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
out = work / "shuffled.bam"
returned = shuffle_bam(
    ["-o", str(out), "-n", "2", "-T", str(work / "tmp"), str(bam)], verbose=False
)
print(repr(returned))
print(sorted(read.query_name for read in pysam.AlignmentFile(str(out))))
```

- **示例输出：**

```text
''
['r1', 'r2', 'r3']
```

## 9) `dnakit.alignment.tabix_records()`

- **作用：** 用 pysam 的 tabix 解析器读取文本文件的全部记录。
- **处理方法：** 以二进制模式打开文件，返回 `list(pysam.tabix_iterator(handle, parser))`。
- **API：** `tabix_records(path: 'str', parser: 'Any', *, verbose: bool = True) -> 'list[Any]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `path` 是已存在的文件（直接打开；示例使用普通 BED 文件，没有索引）。`parser` 是 pysam 解析器对象，如 `pysam.asBed()`、`pysam.asGTF()`、`pysam.asVCF()` 或 `pysam.asTuple()`。
- **输出：** 解析后的记录列表。使用 `pysam.asBed()` 时，每条记录有 `contig`、`start`（0 起始）、`end` 和 `name`。
- **示例代码：**

```python
import pysam
from dnakit.alignment import tabix_records

records = tabix_records("regions.bed", pysam.asBed(), verbose=False)
print(len(records))
for record in records:
    print(record.contig, record.start, record.end)
```

- **示例输出：**

```text
2
g1 2 10
g2 0 8
```

## 10) `dnakit.alignment.target_cut()`

- **作用：** 用你自己的参数运行 samtools 的 `targetcut` 命令。
- **处理方法：** 用你的参数调用 `pysam.targetcut()`，返回 pysam 捕获的标准输出。该命令读取 BAM 文件并打印 SAM 行；在示例中，它为被覆盖的区域输出一条由各读段的碱基与质量合成的伪读段（没有覆盖的位置是 `N`，质量为 `!`）。
- **API：** `target_cut(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `arguments` 是传给 `samtools targetcut` 的字符串列表：选项如 `-Q <整数>`（最低碱基质量），后面跟 BAM 路径。若参数中没有输入文件，samtools 会等待标准输入，所以一定要给出 BAM。
- **输出：** 字符串形式的 SAM 文本，每个切出的区域一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import target_cut

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
text = target_cut([str(bam)], verbose=False)
for line in text.splitlines():
    fields = line.split("\t")
    print(fields[0], fields[1], fields[3], fields[5], fields[9])
```

- **示例输出：**

```text
chr1:2-25 0 2 24M NNNACGTACGTTGACCGAGGGCCC
```

## 11) `dnakit.alignment.call_copy_number_states()`

- **作用：** 根据 VCF 中的 B 等位基因频率（BAF）和 log R ratio（LRR）沿染色体推断拷贝数状态。
- **处理方法：** 通过 `pysam.bcftools.cnv` 运行 `bcftools cnv`，参数为 `-s <query_sample> -o <output_dir>`（给出时再加 `-c`、`-f`、`-r`），并设置 `catch_stdout=False`。然后读取 `cn.<sample>.tab`（每个位点一行）和 `summary.<sample>.tab`（其中的 `RG` 行是片段）。
- **API：** `call_copy_number_states(vcf: 'str', output_dir: 'str', *, query_sample: 'str', control_sample: 'str | None' = None, allele_frequency_file: 'str | None' = None, regions: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `vcf` 是已存在的 VCF，且查询样本带有 FORMAT 字段 `BAF` 和 `LRR`。`output_dir` 必须尚不存在且其上级目录必须存在，函数会创建它。`query_sample` 是不以 `-` 开头的非空名称。`control_sample`（`-c`）、`allele_frequency_file`（`-f`，已存在的文件）和 `regions`（`-r`，不以 `-` 开头的区域字符串）都是可选的。
- **输出：** 返回字典：`status`、`backend`、`command`、`sample`、`outputs`（存在的 `cn`、`dat`、`summary`、`plot_script` 文件路径）、`segments`（每项含 `chromosome`、`start`、`end`、`copy_number`、`quality`、`n_sites`、`n_hets`）和 `sites`（每项含 `chromosome`、`position`、`copy_number` 以及四个 `probabilities`）。片段的 `quality` 可能是 `nan`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.alignment import call_copy_number_states

out_dir = Path(tempfile.mkdtemp()) / "cnv"
result = call_copy_number_states("cnv.vcf", str(out_dir), query_sample="S1", verbose=False)
print(result["status"], result["sample"], sorted(result["outputs"]))
for seg in result["segments"]:
    print(seg["chromosome"], seg["start"], seg["end"], seg["copy_number"], seg["n_sites"])
print(len(result["sites"]), result["sites"][0]["position"], result["sites"][0]["copy_number"])
```

- **示例输出：**

```text
completed S1 ['cn', 'dat', 'plot_script', 'summary']
chr1 1000 13999 2 14
chr1 14000 25999 3 12
chr1 26000 36000 2 10
36 1000 2
```

## 12) `dnakit.alignment.depad_alignments()`

- **作用：** 把填充（padded）比对转换为非填充坐标并返回记录。
- **处理方法：** 用 `-T <padded_reference> -o <output_alignment> <alignment>`（要输出 SAM 时再加 `-s`）和 `catch_stdout=False` 运行 `pysam.depad()`，再用 `pysam.AlignmentFile` 读回输出文件。
- **API：** `depad_alignments(alignment: 'str', padded_reference: 'str', output_alignment: 'str', *, output_sam: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** `alignment` 和 `padded_reference`（填充列为 `*` 的 FASTA）必须是已存在的文件。`output_alignment` 必须尚不存在，其所在目录必须存在。`output_sam=True` 写 SAM，否则写 BAM。读段不能含 P 或 I 的 CIGAR 操作。
- **输出：** 返回字典：`status`、`backend`、`command`、`output_format`（`BAM` 或 `SAM`）、`outputs`、`record_count` 和 `records`（最多前 200 条，每条含 `query_name`、`flag`、`reference_name`、1 起始的 `reference_start`、`cigar`、`sequence`）。samtools 失败时抛出 `BackendExecutionError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import depad_alignments

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
result = depad_alignments(str(bam), "padded.fa", str(work / "unpadded.bam"), verbose=False)
print(result["status"], result["output_format"], result["record_count"])
for read in result["records"]:
    print(read["query_name"], read["reference_start"], read["cigar"], read["sequence"])
```

- **示例输出：**

```text
completed BAM 2
r1 1 10M ACGTACGTAC
r2 3 6M GTCGTA
```

## 13) `dnakit.alignment.pad_to_unpad_alignments()`

- **作用：** 通过旧版 samtools 命令名 `pad2unpad` 把填充（padded）比对转换为非填充坐标并返回记录。
- **处理方法：** 代码路径与 `depad_alignments()` 相同，只是调用 `pysam.pad2unpad()`。用随包的 samtools，这两个命令在本页数据上给出相同的记录。
- **API：** `pad_to_unpad_alignments(alignment: 'str', padded_reference: 'str', output_alignment: 'str', *, output_sam: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pysam（extra `dnakit[bio-reads]`）。samtools 代码已随 pysam 打包，不需要外部程序。
- **输入：** 与 `depad_alignments()` 相同：已存在的 `alignment` 和 `padded_reference`，新的 `output_alignment`，可选的 `output_sam`。
- **输出：** 与 `depad_alignments()` 相同的字典；`command` 中显示的是 `pad2unpad`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pad_to_unpad_alignments

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
result = pad_to_unpad_alignments(
    str(bam), "padded.fa", str(work / "unpadded.sam"), output_sam=True, verbose=False
)
print(result["status"], result["output_format"], result["record_count"])
print(result["command"][:2])
for read in result["records"]:
    print(read["query_name"], read["reference_start"], read["cigar"])
```

- **示例输出：**

```text
completed SAM 2
['samtools', 'pad2unpad']
r1 1 10M
r2 3 6M
```

## 14) `dnakit.alignment.pysam_tview_text()`

- **作用：** 用 `samtools tview` 把 BAM 或 CRAM 比对的一屏内容渲染为文本。
- **处理方法：** 运行 PATH 中（或 `samtools=` 指定）的 samtools：`samtools tview -d T -p <chrom:pos> -w <width> <alignment> [reference]`，把标准输出写入新文件 `output`，并返回各行。pysam 自带的 samtools 没有 `tview`，所以函数不使用 `pysam.tview`。版本号来自 `samtools --version`。
- **API：** `pysam_tview_text(alignment: 'str | Path', position: 'str', *, output: 'str | Path', reference: 'str | Path | None' = None, width: 'int' = 80, samtools: 'str | Path' = 'samtools', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 samtools（源码中的安装提示说测试用的是 1.24），放在 PATH 中或用 `samtools=` 指定；pysam 自带的那份不能用。缺少时函数抛出带安装提示的 `BackendUnavailableError`。本环境的 PATH 中没有 samtools，所以示例展示的是这个错误；真正运行的调用写在示例之后。
- **输入：** `alignment` 是已存在的 BAM 或 CRAM，`position` 是 `chrom:pos`，位置 1 起始。`output` 必须是新文件且其目录存在，`<output>.stderr` 也不能存在。`width` 是 20 到 1000 的整数。`reference` 是可选的 FASTA。这些检查先于查找 samtools，所以比对文件不存在时先得到 `ConfigurationError`。
- **输出：** 返回字典：`status`、`backend`、`program`、`version`、`alignment`、`reference`、`position`、`width`、`output`、`line_count`、`lines`（文本视图）和 `log`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pysam_tview_text
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
try:
    pysam_tview_text(
        bam, "chr1:12", output=work / "view.txt", samtools="samtools_not_installed", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "samtools executable 'samtools_not_installed' was not found" in str(error))
print((work / "view.txt").exists())
```

- **示例输出：**

```text
BackendUnavailableError True
False
```

安装依赖后的用法（不执行）：

```python
from dnakit.alignment import pysam_tview_text

result = pysam_tview_text("reads.bam", "chr1:12", output="view.txt", reference="ref.fa", width=60)
print(result["version"], result["line_count"])
print("\n".join(result["lines"][:5]))
```

## 15) `dnakit.sequence.faidx_check_default_seq()`

- **作用：** 检查用于缺失位置和屏蔽的单字符默认碱基。
- **处理方法：** `None` 返回 `None`；恰好为一个字符（不是 NUL）时返回该值本身；其余情况抛出 `ConfigurationError`。规则与 pyfaidx 的 `pyfaidx.cli.check_seq_length`（`--default-seq`）一致；函数本身不调用 pyfaidx。
- **API：** `faidx_check_default_seq(value: 'str | None', *, verbose: bool = True) -> 'str | None'`。
- **运行要求：** 安装 DNAKit 后即可使用（不导入任何可选包）。
- **输入：** `value` 是 `None` 或字符串。短横线 `-` 可作为值。
- **输出：** `None` 或这个单字符字符串。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import faidx_check_default_seq

print(faidx_check_default_seq(None, verbose=False))
print(faidx_check_default_seq("-", verbose=False))
try:
    faidx_check_default_seq("NN", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "default_seq must be a single character" in str(error))
```

- **示例输出：**

```text
None
-
ConfigurationError True
```

## 16) `dnakit.sequence.faidx_fetch_sequences()`

- **作用：** 把指定区域（或全部序列）提取到新的 FASTA 文件。
- **处理方法：** 运行 `faidx`，旗标由参数构建（`-c` 取互补、`-r` 反向、`-y` 自动链方向、`-n` 无名称、`-f` 长名称、`-t` 无坐标、`-l` 惰性、`-s` 默认碱基、`-d` 分隔符、`--no-rebuild`，BED 文件用 `-b`），并把其标准输出写入 `output`。函数不使用 `faidx -o`，因为 pyfaidx 0.8.2 用该方式写较大文件时可能丢失末尾。
- **API：** `faidx_fetch_sequences(fasta: 'str | Path', output: 'str | Path', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, complement: 'bool' = False, reverse: 'bool' = False, auto_strand: 'bool' = False, no_names: 'bool' = False, long_names: 'bool' = False, no_coords: 'bool' = False, lazy: 'bool' = False, default_seq: 'str | None' = None, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），它会安装 `faidx` 程序。DNAKit 依次在 `faidx_path`、PATH 和当前 Python 所在目录中查找 `faidx`。 没有索引时，faidx 会在 FASTA 旁写出 `.fai` 索引。
- **输入：** `fasta` 是已存在的文件，`output` 是已存在目录中的新文件。`regions`（`name:start-end`，1 起始且含两端）与 `bed` 二选一，不能同时给；两者都不给则写出全部序列。`no_names` 与 `long_names` 不能同时为真。`default_seq` 是单个字符。
- **输出：** 返回字典：`status`（退出码为 0 时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（新文件）和 `log_tail`（faidx 的标准错误）。 FASTA 中不存在的区域会记录在 `log_tail` 中，且不产生记录。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_fetch_sequences

out = Path(tempfile.mkdtemp()) / "regions.fa"
result = faidx_fetch_sequences(
    "genes.fa", str(out), regions=["g1:13-20", "g2:9-16"], complement=True, reverse=True, verbose=False
)
print(result["status"], result["backend"], result["return_code"])
print(out.read_text(), end="")
```

- **示例输出：**

```text
completed faidx 0
>g1:20-13 (complement)
CTCGGTCA
>g2:16-9 (complement)
TTTTGGGG
```

## 17) `dnakit.sequence.faidx_mask_sequences()`

- **作用：** 在 FASTA 的新副本中，把指定区域的碱基替换为屏蔽字符。
- **处理方法：** 先把 FASTA 复制到 `output`，再对副本针对这些区域运行 `faidx -m -s <default_seq>`。输入文件不会被修改。不提供小写屏蔽（`-M`），因为 pyfaidx 0.8.2 在该选项上会失败。
- **API：** `faidx_mask_sequences(fasta: 'str | Path', output: 'str | Path', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, default_seq: 'str' = 'N', delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），它会安装 `faidx` 程序。DNAKit 依次在 `faidx_path`、PATH 和当前 Python 所在目录中查找 `faidx`。
- **输入：** `fasta` 是已存在的文件，`output` 是新文件。至少需要一个区域（`name:start-end`，1 起始且含两端）或一个 `bed` 文件（否则抛出 `ConfigurationError`）。`default_seq` 是替换用的单个字符（默认 `N`）。
- **输出：** 返回字典：`status`、`backend`、`version`、`command`、`return_code`、`outputs`（屏蔽后的副本）和 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_mask_sequences

out = Path(tempfile.mkdtemp()) / "masked.fa"
result = faidx_mask_sequences("genes.fa", str(out), regions=["g1:1-4", "g3:5-7"], verbose=False)
print(result["status"])
print(out.read_text(), end="")
print(Path("genes.fa").read_text().splitlines()[1])
```

- **示例输出：**

```text
completed
>g1
NNNNACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
>g3
GATTNNN
ACGTACGTACGTTGACCGAGGGCCC
```

## 18) `dnakit.sequence.faidx_parse_bed_entry()`

- **作用：** 把一行 BED 解析为染色体名及 0 起始的起止位置。
- **处理方法：** 在 Python 进程中调用 `pyfaidx.bed_split()`。注释行（以 `#` 开头）三个字段都是 `None`。格式错误的行抛出 `ConfigurationError`。
- **API：** `faidx_parse_bed_entry(line: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），不运行外部程序。
- **输入：** `line` 是字符串（制表符分隔，使用前三列）。
- **输出：** 含 `name`、`start`、`end` 的字典（整数，BED 坐标，0 起始、半开区间）；注释行三项均为 `None`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import faidx_parse_bed_entry

print(faidx_parse_bed_entry("g1\t2\t10", verbose=False))
print(faidx_parse_bed_entry("# a comment", verbose=False))
try:
    faidx_parse_bed_entry("g1\tx\t10", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Malformed BED entry" in str(error))
```

- **示例输出：**

```text
{'name': 'g1', 'start': 2, 'end': 10}
{'name': None, 'start': None, 'end': None}
ConfigurationError True
```

## 19) `dnakit.sequence.faidx_parse_ucsc_region()`

- **作用：** 把 UCSC 区域字符串解析为名称及 0 起始的起点和不含终点的终点。
- **处理方法：** 调用 `pyfaidx.ucsc_split()`。输入是 1 起始且含两端的坐标（`g1:3-10`），结果是 0 起始、终点不含（start 为 2，end 为 10）。只有名称时 start 和 end 为 `None`，表示整条序列。
- **API：** `faidx_parse_ucsc_region(region: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），不运行外部程序。
- **输入：** `region` 是非空且不以 `-` 开头的字符串。
- **输出：** 含 `name`、`start`、`end` 的字典。
- **示例代码：**

```python
from dnakit.sequence import faidx_parse_ucsc_region

print(faidx_parse_ucsc_region("g1:3-10", verbose=False))
print(faidx_parse_ucsc_region("g1", verbose=False))
```

- **示例输出：**

```text
{'name': 'g1', 'start': 2, 'end': 10}
{'name': 'g1', 'start': None, 'end': None}
```

## 20) `dnakit.sequence.faidx_split_regions()`

- **作用：** 把区域字符串或 BED 文件转换为含名称、起点、终点的记录。
- **处理方法：** 区域字符串逐个调用 `faidx_parse_ucsc_region()`；BED 文件逐行调用 `faidx_parse_bed_entry()` 并丢弃注释行。保持输入顺序。
- **API：** `faidx_split_regions(regions: 'Sequence[str]' = (), *, bed: 'str | Path | None' = None, verbose: bool = True) -> 'list[dict[str, Any]]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），不运行外部程序。
- **输入：** `regions` 是 `name:start-end` 字符串序列；`bed` 是已存在的 BED 文件。两者同时给出会抛出 `ConfigurationError`。两者都不给则返回空列表。
- **输出：** 含 `name`、`start`、`end` 的字典列表。
- **示例代码：**

```python
from dnakit.sequence import faidx_split_regions

print(faidx_split_regions(["g1:3-10", "g2"], verbose=False))
print(faidx_split_regions(bed="regions.bed", verbose=False))
print(faidx_split_regions(verbose=False))
```

- **示例输出：**

```text
[{'name': 'g1', 'start': 2, 'end': 10}, {'name': 'g2', 'start': None, 'end': None}]
[{'name': 'g1', 'start': 2, 'end': 10}, {'name': 'g2', 'start': 0, 'end': 8}]
[]
```

## 21) `dnakit.sequence.faidx_transform_sequences()`

- **作用：** 把 FASTA 的序列或区域描述为表格，而不是 FASTA。
- **处理方法：** 运行 `faidx -i <transform>`，并把其标准输出写入 `output`。`bed` 写出名称、0 起始的起点和终点；`chromsizes` 写出名称和长度；`nucleotide` 写出一行表头以及 A、T、C、G、N 和其他字符的计数；`transposed` 写出名称、起点、终点和序列。
- **API：** `faidx_transform_sequences(fasta: 'str | Path', output: 'str | Path', transform: 'str', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, complement: 'bool' = False, reverse: 'bool' = False, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），它会安装 `faidx` 程序。DNAKit 依次在 `faidx_path`、PATH 和当前 Python 所在目录中查找 `faidx`。
- **输入：** `transform` 取 `bed`、`chromsizes`、`nucleotide`、`transposed` 之一。`fasta` 是已存在的文件，`output` 是新文件。`regions`（或 `bed`）、`complement`、`reverse`、`delimiter`、`no_rebuild` 可选。
- **输出：** 返回字典：`status`（退出码为 0 时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（新文件）和 `log_tail`（faidx 的标准错误）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_transform_sequences

out = Path(tempfile.mkdtemp()) / "composition.tsv"
result = faidx_transform_sequences("genes.fa", str(out), "nucleotide", verbose=False)
print(result["status"])
for line in out.read_text().splitlines():
    print(line.rstrip("	"))
```

- **示例输出：**

```text
completed
name	start	end	A	T	C	G	N	others
g1	1	25	5	4	8	8	0
g2	1	26	6	6	6	6	2
g3	1	7	3	2	1	1	0
```

## 22) `dnakit.sequence.faidx_wrap_sequence()`

- **作用：** 把序列按固定宽度折行。
- **处理方法：** 调用 `pyfaidx.wrap_sequence()`，并以列表返回各行。每一行（包括最后一行）都以换行符结尾。
- **API：** `faidx_wrap_sequence(sequence: 'str', width: 'int', *, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），不运行外部程序。
- **输入：** `sequence` 是字符串；`width` 是不小于 1 的整数。
- **输出：** 字符串列表，每行最多 `width` 个字符外加换行符。
- **示例代码：**

```python
from dnakit.sequence import faidx_wrap_sequence

print(faidx_wrap_sequence("ACGTACGTAC", 4, verbose=False))
```

- **示例输出：**

```text
['ACGT\n', 'ACGT\n', 'AC\n']
```

## 23) `dnakit.sequence.faidx_write_sequences()`

- **作用：** 按名称正则表达式和长度范围选出 FASTA 序列，写入一个文件或每条序列一个文件。
- **处理方法：** 运行 `faidx -g <regex>`（`-v` 反选，`-a low,high` 为长度范围，`-x` 为拆分文件，输出旗标与 `faidx_fetch_sequences()` 相同）。给 `output` 时，标准输出写入新文件；给 `split_directory` 时，faidx 在新目录中为每条选中的序列写一个文件，以序列名和坐标命名，函数会列出这些文件。
- **API：** `faidx_write_sequences(fasta: 'str | Path', *, output: 'str | Path | None' = None, split_directory: 'str | Path | None' = None, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, regex: 'str' = '.*', invert_match: 'bool' = False, size_range: 'Sequence[int] | None' = None, complement: 'bool' = False, reverse: 'bool' = False, auto_strand: 'bool' = False, no_names: 'bool' = False, long_names: 'bool' = False, no_coords: 'bool' = False, lazy: 'bool' = False, default_seq: 'str | None' = None, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyfaidx（extra `dnakit[bio-pyfaidx]`），它会安装 `faidx` 程序。DNAKit 依次在 `faidx_path`、PATH 和当前 Python 所在目录中查找 `faidx`。
- **输入：** `output`（新文件）和 `split_directory`（新目录）必须恰好给一个，否则抛出 `ConfigurationError`。`regex` 用于选序列名且不锚定，要精确匹配请写 `^...$`。`size_range=(low, high)` 含两端，保留长度（或区域长度）在范围内的序列。`regions`/`bed` 限定区域；其他旗标同 `faidx_fetch_sequences()`。
- **输出：** 返回字典：`status`（退出码为 0 时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（新文件）和 `log_tail`（faidx 的标准错误）。 使用 `split_directory` 时，`outputs` 列出该目录中写出的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_write_sequences

work = Path(tempfile.mkdtemp())
result = faidx_write_sequences(
    "genes.fa", output=str(work / "long.fa"), size_range=(10, 100), verbose=False
)
print(result["status"], [Path(path).name for path in result["outputs"]])
print((work / "long.fa").read_text(), end="")
split = faidx_write_sequences("genes.fa", split_directory=str(work / "split"), regex="^g3$", verbose=False)
print(split["status"], [Path(path).name for path in split["outputs"]])
```

- **示例输出：**

```text
completed ['long.fa']
>g1
ACGTACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
completed ['g3.fa']
```

## 24) `dnakit.fasta_fastq.cutadapt_full_report()`

- **作用：** 修剪 FASTQ 文件的接头和低质量末端，并返回 cutadapt 的完整文本报告。
- **处理方法：** 与 `cutadapt_trim()` 使用同一段代码，报告类型为 `--report full`：不经 shell 运行 `cutadapt -e <error_rate> -j <cores> [-a/-g/-b 接头] [-q] [-m] [-M] [--discard-untrimmed] -o <output> <input>`。报告包含摘要、读段去向统计、每个接头一节以及长度直方图。
- **API：** `cutadapt_full_report(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `input_path` 是已存在的 FASTQ，`output_path` 是新文件。接头（`adapter_3p` 对应 `-a`，`adapter_5p` 对应 `-g`，`anywhere_adapter` 对应 `-b`）是 IUPAC DNA 字母，会转为大写。`error_rate` 在 0 到 1 之间；`quality_cutoff` 形如 `20` 或 `10,20`；`minimum_length` 与 `maximum_length` 非负（最大值不小于最小值）；`cores` 至少为 1。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（cutadapt 的标准输出）和 `log_tail`（标准错误的末尾）。 报告中含命令行，其中有文件路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_full_report

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_full_report("reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", verbose=False)
print(result["status"], result["return_code"])
for line in result["report"].splitlines():
    if line.startswith(("Total reads", "Reads with adapters", "Reads written", "Sequence:")):
        print(line)
```

- **示例输出：**

```text
completed 0
Total reads processed:                       4
Reads with adapters:                         2 (50.0%)
Reads written (passing filters):             4 (100.0%)
Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 2 times
```

## 25) `dnakit.fasta_fastq.cutadapt_json_report()`

- **作用：** 修剪 FASTQ 文件的接头和低质量末端，并同时写出 cutadapt 的 JSON 报告。
- **处理方法：** 与 `cutadapt_trim()` 使用同一段代码，参数为 `--report full --json <json_path>`。运行后解析 JSON 文件，并放在返回值的 `json` 中。
- **API：** `cutadapt_json_report(input_path: 'str', output_path: 'str', json_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** 同 `cutadapt_full_report()`，另有 `json_path`，是新文件路径（必须尚不存在）。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（cutadapt 的标准输出）和 `log_tail`（标准错误的末尾）。 另外有 `json`，即解析后的 JSON 报告（键如 `read_counts`、`basepair_counts`、`adapters_read1`、`command_line_arguments`）；此时 `outputs` 含修剪后的 FASTQ 和 JSON 文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_json_report

work = Path(tempfile.mkdtemp())
result = cutadapt_json_report(
    "reads.fq", str(work / "trimmed.fq"), str(work / "report.json"),
    adapter_3p="AGATCGGAAGAGC", verbose=False,
)
print(result["status"], len(result["outputs"]))
counts = result["json"]["read_counts"]
print(counts["input"], counts["output"], counts["read1_with_adapter"])
```

- **示例输出：**

```text
completed 2
4 4 2
```

## 26) `dnakit.fasta_fastq.cutadapt_minimal_report()`

- **作用：** 修剪 FASTQ 文件的接头和低质量末端，并返回 cutadapt 的两行制表符分隔的精简报告。
- **处理方法：** 与 `cutadapt_trim()` 使用同一段代码，报告类型为 `--report minimal`。
- **API：** `cutadapt_minimal_report(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** 同 `cutadapt_full_report()`。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（cutadapt 的标准输出）和 `log_tail`（标准错误的末尾）。 `report` 的第一行是列名（`status`、`in_reads`、`in_bp`、`too_short`、`too_long`、`too_many_n`、`out_reads`、`w/adapters`、`qualtrim_bp`、`out_bp`），第二行是对应的值。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_minimal_report

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_minimal_report("reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", verbose=False)
header, values = result["report"].strip().splitlines()
for name, value in zip(header.split("\t"), values.split("\t")):
    print(name, value)
```

- **示例输出：**

```text
status OK
in_reads 4
in_bp 85
too_short 0
too_long 0
too_many_n 0
out_reads 4
w/adapters 2
qualtrim_bp 0
out_bp 54
```

## 27) `dnakit.fasta_fastq.cutadapt_trim()`

- **作用：** 修剪单个 FASTQ 文件的 3' 和 5' 接头以及低质量末端。
- **处理方法：** 运行 `cutadapt --report <full|minimal> -e <error_rate> -j <cores>`，并按需加上 `-a`、`-g`、`-b`、`-q`、`-m`、`-M` 和 `--discard-untrimmed`，把读段写入 `output_path`。其余三个报告函数就是固定了报告类型的本函数。
- **API：** `cutadapt_trim(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, report: 'str' = 'full', cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** 同 `cutadapt_full_report()`，另有 `report`，取 `"full"`（默认）或 `"minimal"`。无效值会在 cutadapt 启动前抛出 `ConfigurationError`。`discard_untrimmed=True` 会丢弃没有匹配到接头的读段。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（cutadapt 的标准输出）和 `log_tail`（标准错误的末尾）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_trim

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_trim(
    "reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", quality_cutoff="10",
    minimum_length=10, report="minimal", verbose=False,
)
print(result["status"], result["return_code"])
lines = out.read_text().splitlines()
for name, sequence in zip(lines[0::4], lines[1::4]):
    print(name, sequence)
```

- **示例输出：**

```text
completed 0
@r2 TTGACCGATC
@r3 GGGCCCAAATTTGGGAAA
@r4 NNACGTACGTTGAACG
```

## 28) `dnakit.fasta_fastq.cutadapt_adapter_action()`

- **作用：** 在每条读段中查找接头，并对每个匹配执行一种操作：trim（切除）、retain（保留）、mask（用 N 屏蔽）、lowercase（改小写）、crop（裁剪）或 none（不处理）。
- **处理方法：** 运行 `cutadapt --report minimal -e <error_rate> --action <action> -n <times> [--rc]`，接头由 `-a`、`-g`、`-b` 给出。加 `--rc` 时，cutadapt 还会搜索每条读段的反向互补序列，并把这样的读段以反向互补形式写出，名称后加 ` rc`。
- **API：** `cutadapt_adapter_action(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, action: 'str' = 'trim', times: 'int' = 1, reverse_complement: 'bool' = False, error_rate: 'float' = 0.1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `input_path` 是已存在的 FASTQ，`output_path` 是新文件。至少需要一个接头。`action` 取 `trim`（默认）、`retain`、`mask`、`lowercase`、`crop` 或 `none`。`times`（至少 1）是在一条读段中搜索接头的次数；`reverse_complement` 设置 `--rc`；`error_rate` 在 0 到 1 之间。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（cutadapt 的标准输出）和 `log_tail`（标准错误的末尾）。 这里的 `report` 是单行精简报告。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_adapter_action

work = Path(tempfile.mkdtemp())
result = cutadapt_adapter_action(
    "reads.fq", str(work / "masked.fq"), adapter_3p="AGATCGGAAGAGC", action="mask", verbose=False
)
print(result["status"])
lines = (work / "masked.fq").read_text().splitlines()
print(lines[1::4])
rc = cutadapt_adapter_action(
    "reads.fq", str(work / "rc.fq"), adapter_5p="TTTCCCAAATTTGGG",
    reverse_complement=True, verbose=False,
)
lines = (work / "rc.fq").read_text().splitlines()
print(lines[8], lines[9])
```

- **示例输出：**

```text
completed
['ACGTTGCANNNNNNNNNNNNNNNNNN', 'TTGACCGATCNNNNNNNNNNNNN', 'GGGCCCAAATTTGGGAAA', 'NNACGTACGTTGAACGNN']
@r3 rc CCC
```

## 29) `dnakit.fasta_fastq.cutadapt_shorten()`

- **作用：** 把每条读段缩短到固定长度。
- **处理方法：** 运行 `cutadapt --report minimal -l <length>`。正的长度保留每条读段的前面若干碱基，负的长度保留最后若干碱基。已经更短的读段不变。
- **API：** `cutadapt_shorten(input_path: 'str', output_path: 'str', *, length: 'int', cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `input_path` 是已存在的 FASTQ，`output_path` 是新文件。`length` 是非零整数（为 0 或非整数时抛出 `ConfigurationError`）。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（单行精简报告）和 `log_tail`（标准错误的末尾）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_shorten

out = Path(tempfile.mkdtemp()) / "short.fq"
result = cutadapt_shorten("reads.fq", str(out), length=-6, verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4])
```

- **示例输出：**

```text
completed
['CACACG', 'AAGAGC', 'GGGAAA', 'AACGNN']
```

## 30) `dnakit.fasta_fastq.cutadapt_trim_n_ends()`

- **作用：** 去掉每条读段两端的 N 碱基。
- **处理方法：** 运行 `cutadapt --report minimal --trim-n`。
- **API：** `cutadapt_trim_n_ends(input_path: 'str', output_path: 'str', *, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `input_path` 是已存在的 FASTQ，`output_path` 是新文件。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（单行精简报告）和 `log_tail`（标准错误的末尾）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_trim_n_ends

out = Path(tempfile.mkdtemp()) / "no_n.fq"
result = cutadapt_trim_n_ends("reads.fq", str(out), verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4][-1])
```

- **示例输出：**

```text
completed
ACGTACGTTGAACG
```

## 31) `dnakit.fasta_fastq.cutadapt_unconditional_cut()`

- **作用：** 从每条读段的 5' 端和/或 3' 端去掉固定数量的碱基。
- **处理方法：** 5' 端运行 `cutadapt --report minimal -u <cut_5p>`，3' 端运行 `-u -<cut_3p>`，只处理非零的一端或两端。
- **API：** `cutadapt_unconditional_cut(input_path: 'str', output_path: 'str', *, cut_5p: 'int' = 0, cut_3p: 'int' = 0, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 cutadapt（extra `dnakit[bio-cutadapt]`），它会安装 `cutadapt` 程序。DNAKit 依次在 `cutadapt_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `input_path` 是已存在的 FASTQ，`output_path` 是新文件。`cut_5p` 和 `cut_3p` 是非负整数，至少有一个为正数（否则抛出 `ConfigurationError`）。
- **输出：** 返回字典：`status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`report`（单行精简报告）和 `log_tail`（标准错误的末尾）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_unconditional_cut

out = Path(tempfile.mkdtemp()) / "cut.fq"
result = cutadapt_unconditional_cut("reads.fq", str(out), cut_5p=2, cut_3p=3, verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4])
```

- **示例输出：**

```text
completed
['GTTGCAAGATCGGAAGAGCAC', 'GACCGATCAGATCGGAAG', 'GCCCAAATTTGGG', 'ACGTACGTTGAAC']
```

## 32) `dnakit.sequencing_qc.umi_count_per_gene()`

- **作用：** 统计 BAM 文件中每个基因的 UMI 组数，并写出基因计数表。
- **处理方法：** 运行 `umi_tools count`，带 UMI 选项（`--extract-umi-method`、`--umi-separator` 或 `--umi-tag`、`--method`、`--edit-distance-threshold`、`--paired`），`--gene-tag <tag>` 或 `--per-contig`，以及 `--mapping-quality`。读段按 UMI 和位置分组（同 `umi_group_bam()`），再统计每个基因的组数。
- **API：** `umi_count_per_gene(bam: 'str', output: 'str', *, gene_tag: 'str | None' = None, per_contig: 'bool' = False, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, mapping_quality: 'int' = 0, paired: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `bam` 是已存在、按坐标排序且有索引的 BAM，`output` 是新文件。`gene_tag`（用 BAM 标签标明每条读段的基因）与 `per_contig=True` 必须二选一，不能同时给也不能都不给（否则 `ConfigurationError`）。`extract_umi_method` 取 `read_id`（UMI 位于读段名中 `umi_separator` 之后）、`tag`（`umi_tag`，默认 `RX`）或 `umis`。`method` 取 `adjacency`、`directional`、`percentile`、`unique` 或 `cluster`。`edit_distance_threshold` 与 `mapping_quality` 非负。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 输出 TSV 的列为 `gene` 和 `count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_count_per_gene

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_count_per_gene(bam, str(work / "counts.tsv"), gene_tag="XT", verbose=False)
print(result["status"])
print((work / "counts.tsv").read_text(), end="")
```

- **示例输出：**

```text
completed
gene	count
geneA	2
geneB	1
```

## 33) `dnakit.sequencing_qc.umi_count_table()`

- **作用：** 从由读段 id 和基因组成的两列表中统计每个基因的 UMI 组数。
- **处理方法：** 运行 `umi_tools count_tab`，带 `--barcode-separator`、`--method`、`--edit-distance-threshold`，`per_cell` 时再加 `--per-cell`。对每个基因的 UMI 聚类后统计组数。
- **API：** `umi_count_table(table: 'str', output: 'str', *, per_cell: 'bool' = False, barcode_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `table` 是已存在的制表符分隔文件，第 1 列是读段 id，第 2 列是基因；UMI 是读段 id 中 `barcode_separator`（默认 `_`）之后的最后一段。`per_cell=True` 时读段 id 必须以 `cell_barcode<sep>UMI` 结尾：代码把倒数第二段当作细胞（程序自带的帮助文字说的正好相反）。`output` 是新文件。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 输出 TSV 的列为 `gene` 和 `count`，使用 `per_cell` 时为 `cell`、`gene` 和 `count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_count_table

work = Path(tempfile.mkdtemp())
table = work / "read_gene.tsv"
table.write_text("r1_CELL1_AAAA\tgeneA\nr2_CELL1_AAAA\tgeneA\nr3_CELL2_AAAA\tgeneA\nr4_CELL2_CCCC\tgeneB\n")
result = umi_count_table(str(table), str(work / "counts.tsv"), per_cell=True, verbose=False)
print(result["status"])
print((work / "counts.tsv").read_text(), end="")
```

- **示例输出：**

```text
completed
cell	gene	count
CELL1	geneA	1
CELL2	geneA	1
CELL2	geneB	1
```

## 34) `dnakit.sequencing_qc.umi_dedup_bam()`

- **作用：** 去除按坐标排序的 BAM 文件中的 PCR 重复，每个 UMI 组在每个位置只保留一条读段。
- **处理方法：** 运行 `umi_tools dedup`，带 UMI 选项（同 `umi_group_bam()`），并按需加上 `--ignore-umi`、`--per-contig` 或 `--per-gene --gene-tag <tag>`，把读段写入 `output_bam`。
- **API：** `umi_dedup_bam(bam: 'str', output_bam: 'str', *, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, paired: 'bool' = False, ignore_umi: 'bool' = False, per_contig: 'bool' = False, gene_tag: 'str | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `bam` 是已存在、按坐标排序且有索引的 BAM，`output_bam` 是新文件。`extract_umi_method`、`umi_tag`、`umi_separator`、`method`、`edit_distance_threshold`、`paired` 同 `umi_count_per_gene()`。`ignore_umi=True` 只按位置去重；`per_contig=True` 或 `gene_tag` 在各染色体或各基因内去重。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 输出是 BAM，每个位置上每个 UMI 组只有一条读段。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_dedup_bam

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_dedup_bam(bam, str(work / "dedup.bam"), verbose=False)
print(result["status"])
print(sorted(read.query_name for read in pysam.AlignmentFile(str(work / "dedup.bam"))))
```

- **示例输出：**

```text
completed
['a2_AAAA', 'a4_CCCC', 'a6_GGGG']
```

## 35) `dnakit.sequencing_qc.umi_extract_to_names()`

- **作用：** 把每条读段开头的 UMI 和细胞条形码移到读段名中。
- **处理方法：** 运行 `umi_tools extract --extract-method=<string|regex> --bc-pattern=<pattern> --umi-separator=<sep>`（再按需加 `--3prime`，双端读段还有 `--read2-in`、`--read2-out`、`--bc-pattern2`），以 FASTQ 作为标准输入，新 FASTQ 作为标准输出。条形码碱基从读段中移除，并追加到读段名：细胞条形码和 UMI，以 `_` 分隔。
- **API：** `umi_extract_to_names(fastq: 'str', output: 'str', *, bc_pattern: 'str', extract_method: 'str' = 'string', umi_separator: 'str' = '_', three_prime: 'bool' = False, read2_in: 'str | None' = None, read2_out: 'str | None' = None, bc_pattern2: 'str | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `fastq` 是已存在的文件，`output` 是新文件。`extract_method="string"` 时，`bc_pattern` 由 `C`（细胞条形码）、`N`（UMI）和 `X`（丢弃）组成，且必须含 `N`；`"regex"` 时是带命名分组的正则表达式。双端读段需给出 `read2_in` 和 `read2_out`（read 2 有条形码时再给 `bc_pattern2`）；没有 `read2_in` 却给了 `read2_out` 或 `bc_pattern2` 会抛出 `ConfigurationError`。`three_prime=True` 表示条形码在 3' 端。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 输出 FASTQ 的读段名中带有条形码。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_extract_to_names

out = Path(tempfile.mkdtemp()) / "extracted.fq"
result = umi_extract_to_names("umi_reads.fq", str(out), bc_pattern="CCCCCNNN", verbose=False)
print(result["status"], result["return_code"])
lines = out.read_text().splitlines()
for name, sequence in zip(lines[0::4], lines[1::4]):
    print(name, sequence)
```

- **示例输出：**

```text
completed 0
@r1_ACGTA_AGG GCCCAAATTTGGG
@r2_TTGCC_AGG GCCCAAATTTGGG
@r3_ACGTA_AGG GCCCAAATTTGGG
@r4_ACGTT_TGG GCCCAAATTTGGG
```

## 36) `dnakit.sequencing_qc.umi_group_bam()`

- **作用：** 按 UMI 和比对位置对 BAM 文件的读段分组。
- **处理方法：** 运行 `umi_tools group --output-bam`，带 UMI 选项、`--umi-group-tag=<tag>`，需要时再加 `--group-out=<group_tsv>`。输出 BAM 的每条读段把所属组写入标签 `umi_group_tag`（默认 `BX`）。
- **API：** `umi_group_bam(bam: 'str', output_bam: 'str', *, group_tsv: 'str | None' = None, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', umi_group_tag: 'str' = 'BX', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, paired: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `bam` 是已存在、按坐标排序且有索引的 BAM，`output_bam` 是新文件；`group_tsv` 是可选的新文件。`extract_umi_method`、`umi_tag`、`umi_separator`、`method`、`edit_distance_threshold`、`paired` 同 `umi_count_per_gene()`。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 `group_tsv` 的列为 `read_id`、`contig`、`position`、`gene`、`umi`、`umi_count`、`final_umi`、`final_umi_count` 和 `unique_id`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_group_bam

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_group_bam(bam, str(work / "grouped.bam"), group_tsv=str(work / "groups.tsv"), verbose=False)
print(result["status"], len(result["outputs"]))
for read in pysam.AlignmentFile(str(work / "grouped.bam")):
    print(read.query_name, read.get_tag("BX"))
rows = (work / "groups.tsv").read_text().splitlines()
print(rows[0].split("\t")[-4:])
```

- **示例输出：**

```text
completed 2
a1_AAAA AAAA
a2_AAAA AAAA
a3_AAAT AAAA
a4_CCCC CCCC
a5_GGGG GGGG
a6_GGGG GGGG
['umi_count', 'final_umi', 'final_umi_count', 'unique_id']
```

## 37) `dnakit.sequencing_qc.umi_prepare_bam_for_rsem()`

- **作用：** 为 RSEM 准备双端、比对到转录组的 BAM。
- **处理方法：** 运行 `umi_tools prepare-for-rsem [--tags=<tags>] [--sam]`。该程序把每个读段名的比对配对起来，把逗号分隔的标签从 read 1 复制到 read 2，并按 RSEM 所需的形式写出 BAM。
- **API：** `umi_prepare_bam_for_rsem(bam: 'str', output_bam: 'str', *, tags: 'str | None' = None, sam: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `bam` 是已存在的双端文件（示例用的是按名称排序的），`output_bam` 是新文件。`tags` 是逗号分隔的标签名列表（如 `RX`）。`sam=True` 时读写 SAM 文本。作者检查时单端输入得到的 BAM 里没有读段，所以请使用双端数据。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_prepare_bam_for_rsem

work = Path(tempfile.mkdtemp())
bam = str(work / "pair.bam")
pysam.view("-b", "-o", bam, "pair.sam", catch_stdout=False)
result = umi_prepare_bam_for_rsem(bam, str(work / "rsem.bam"), tags="RX", verbose=False)
print(result["status"])
for read in pysam.AlignmentFile(str(work / "rsem.bam")):
    print(read.query_name, read.flag, read.get_tag("RX"))
```

- **示例输出：**

```text
completed
p1_AAAA 99 AAAA
p1_AAAA 147 AAAA
p2_CCCC 99 CCCC
p2_CCCC 147 CCCC
```

## 38) `dnakit.sequencing_qc.umi_whitelist_cells()`

- **作用：** 估计 FASTQ 文件中存在哪些细胞条形码，并写成白名单。
- **处理方法：** 运行 `umi_tools whitelist --extract-method=<...> --bc-pattern=<...> --method=<reads|umis> --knee-method=<distance|density>`（给出时再加 `--expect-cells` 和 `--set-cell-number`）。不使用 umi_tools 的绘图选项。
- **API：** `umi_whitelist_cells(fastq: 'str', output: 'str', *, bc_pattern: 'str', extract_method: 'str' = 'string', method: 'str' = 'reads', knee_method: 'str' = 'distance', expect_cells: 'int | None' = None, set_cell_number: 'int | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 umi_tools（extra `dnakit[bio-umi-tools]`），它会安装 `umi_tools` 程序。DNAKit 依次在 `umi_tools_path`、PATH 和当前 Python 所在目录中查找。
- **输入：** `fastq` 是已存在的文件，`output` 是新文件。`bc_pattern` 与 `extract_method` 同 `umi_extract_to_names()`，但 `string` 方式下模式必须含 `C`。`method` 按条形码统计读段数（`reads`）或唯一 UMI 数（`umis`）。`knee_method="density"` 需要 `expect_cells`。`expect_cells` 与 `set_cell_number` 是不小于 1 的整数。
- **输出：** 返回字典：`status`（退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）和 `log_tail`（umi_tools 日志的末尾）。结果在输出文件中。 白名单文件每行一个条形码，有四个制表符分隔的列：条形码、纠错后归并到它的条形码、它的计数、以及这些条形码的计数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_whitelist_cells

out = Path(tempfile.mkdtemp()) / "whitelist.txt"
result = umi_whitelist_cells("umi_reads.fq", str(out), bc_pattern="CCCCCNNN", verbose=False)
print(result["status"])
for line in out.read_text().splitlines():
    print(line.split("\t"))
```

- **示例输出：**

```text
completed
['ACGTA', '', '2', '']
['ACGTT', '', '1', '']
['TTGCC', '', '1', '']
```
