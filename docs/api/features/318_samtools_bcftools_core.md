# samtools 与 bcftools 命令行封装：FASTQ 区域、覆盖度直方图、校验和、convert、index、norm、polysomy

这一组函数位于 `dnakit.alignment`（模块 `samtools_p06_cli_functions`）和 `dnakit.variants`（模块 `bcftools_p06_cli_functions`），封装外部程序 [samtools](https://github.com/samtools/samtools) 和 [bcftools](https://github.com/samtools/bcftools) 的命令（测试版本为 1.24）。每个函数校验参数，用固定的参数列表且不经过 shell 运行一条命令，分别捕获标准输出和标准错误，并返回可 JSON 序列化的字典。输出写入明确指定的新路径，已存在的文件不会被覆盖，运行失败时不会留下不完整的输出。

samtools 与 bcftools 是外部程序：请自行安装（例如 `conda install -c bioconda samtools bcftools`，或从 https://www.htslib.org/download/ 下载发行版），并放到 PATH 中，或传入 `samtools_path` / `bcftools_path`。找不到程序时抛出带安装提示的 `BackendUnavailableError`；程序报错退出时抛出 `BackendExecutionError`（代码 `HTSLIB_EXECUTION_FAILED`），并附带日志尾部。两个程序采用 MIT 许可，不随 DNAKit 分发。`bcftools polysomy` 仅存在于用 GSL 编译的版本中。

所有结果都包含 `status`（`completed`）、`backend`、`operation`、`version`（读取自程序）、`command`、`outputs`（写出的文件）、`log_tail`（标准错误的尾部）和 `provenance`（工具、版本、测试版本、可执行文件、引用）；下面各节只列出额外字段。`samtools consensus`、`flagstat`、`stats`、`ampliconstats`、`phase`、`reference -e` 以及 `bcftools call -c` 的 pysam 内置封装是 `dnakit.bioinformatics.pysam_gap_functions` 中已有的函数，插件在后面两份文档中。

## 示例数据

示例从工作目录读取这些文件，并把输出写入临时目录。运行时生成的索引位于输入文件旁或临时目录中。

**示例文件：** `reads.fastq`

```text
@read1
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG
+
55F5D5?D55I5::5DDIDD5FI:I??5??5D5::I555D
@read2
TCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCT
+
F:D:I:?:F:I::IIF:IDDF5D5555F?:IID?DFD?F:
@read3
ACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCA
+
5DF?I5III:??5ID5D5DF5FIDDF5F555I5?DI?DII
@read4
TTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGA
+
I5?:5DD?:DDIFDI?IFD?55?F5II?F??IF5I:DD:5
@read5
GCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAAC
+
:D5?F??I5?F55D?FD?D5:5D5:F?I:D:FIFIIDDFI
@read6
ACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGAC
+
DF5?:I55:::5I??I?:FD5DIF5FF?I:IFD5DIDFF:
```

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:400
r1	0	chr1	1	60	50M	*	0	0	GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCA	??II:I:I::II::5F:I?5D:IF5?:5?DID:??:DFFDI?DI?I5?FI
r2	0	chr1	41	60	50M	*	0	0	TGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACG	?IDFFI:DIIDF:?::F???D?F?D:FDF:F::I?:5DIIDD5ID5::DI
r3	0	chr1	81	60	50M	*	0	0	GATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTT	DI?:IDI??:II?FII:IF:ID:F5FIF5?FIF?:?::IF5::5DIFFFI
r4	0	chr1	121	60	50M	*	0	0	TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTC	5D?I5F5:?5D??FFI5DF?55??:DFFFD:I:?5DFI:I:?DI:?DF:?
r5	0	chr1	161	60	50M	*	0	0	TTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCT	:F?D:I?D5?DD555F5I?:DDFD?5?FD::55?::??FD?5?:5?DDFF
r6	0	chr1	201	60	50M	*	0	0	GACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCG	D5?F?:5?:FI:5:I:5D?5??:5I:I:FD:FDF5F55F:F?FD?:F??5
r7	0	chr1	241	60	50M	*	0	0	CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGC	?55:I55FI?5?ID:FIII?555?DDDD?::D?5F5F???::?D:DD?:D
r8	0	chr1	281	60	50M	*	0	0	ATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTT	FFD:DID:5F5DF::55F:?5DFFI5:IIIID5F5I:?F:I5IFD?:D5F
r9	0	chr2	1	60	50M	*	0	0	TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTG	5:?FF5FD?I:F5FID:FI55?:D?:I?I:I?F5?I::I:F55FI:FIF5
r10	0	chr2	61	60	50M	*	0	0	TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTC	D?DFFDFFD5DF??DI5:I?DF5DD?5IDDDDFI55?:5:F??F:DF?5:
r11	0	chr2	121	60	50M	*	0	0	GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGC	I5I5I:5FDI:?:5:I:55:DF:IIDI5DF:55?F5FFI:5ID?D?:??F
r12	0	chr2	181	60	50M	*	0	0	AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACAT	IFD??II:5I5:55?::F?55DDF55IFFD:ID:D??55ID5I5?D5D:I
```

**示例文件：** `genome.fa`

```text
>chr1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC
TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC
GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA
CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT
ACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAA
GACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATTTGCCGCCTGACAAGTCAATG
CGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATT
AACTGATAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAA
CGGCGAGCTTTACATTTGCTGTGAGAGGTACAGGGATTAGTGAGAAGCCGTGCGTATCAA
>chr2
TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGC
TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTA
GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCA
AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAG
CCGGTGACTCCTAATGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAA
ATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTA
GCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGAT
```

**示例文件：** `raw.tsv`

```text
rs1	chr1	60	GT
rs2	chr1	120	CC
rs3	chr1	180	--
```

**示例文件：** `sample.g.vcf`

```text
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=600>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	G	<*>	.	.	END=5	GT	0/0
chr1	6	.	A	C,<*>	50	.	.	GT	1/1
```

**示例文件：** `multi.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	10	.	A	C,G	.	.	.	GT	1/2
chr1	20	.	CA	GC	.	.	.	GT	0/1
chr1	30	.	A	C	.	.	.	GT	0/1
chr1	30	.	A	C	.	.	.	GT	0/1
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	60	var1	G	T	40	PASS	DP=50	GT	0/1	0/0	1/1	0/1	0/0	0/1
chr1	120	var2	C	T	41	PASS	DP=55	GT	0/0	0/1	0/1	1/1	0/0	0/0
chr1	180	var3	C	A	42	PASS	DP=60	GT	1/1	0/1	0/0	0/1	0/1	./.
chr1	240	var4	ACT	A	43	PASS	DP=65	GT	0/1	0/1	0/0	0/0	1/1	0/1
chr1	300	var5	T	C	44	PASS	DP=70	GT	0/0	0/0	0/1	0/1	0/1	1/1
chr1	360	var6	A	T	45	PASS	DP=75	GT	0/1	1/1	0/1	0/0	0/0	0/1
chr1	420	var7	G	GCG	46	PASS	DP=80	GT	0/0	0/1	0/0	0/1	1/1	0/1
chr1	480	var8	T	C	47	PASS	DP=85	GT	0/1	0/0	1/1	1/1	0/1	0/0
chr1	520	var9	G	C	48	PASS	DP=90	GT	0/0	0/0	0/1	0/0	0/1	0/1
chr1	570	var10	A	C	49	PASS	DP=95	GT	1/1	0/1	0/1	0/0	0/0	0/1
```

## 1) `dnakit.alignment.samtools_extract_fastq_regions()`

- **作用：** 通过索引从 FASTQ 文件中提取整条记录或子区域（对应工作清单中的 `samtools faidx -f REGION` 与 `samtools fqidx REGION`，两者运行的是同一段代码）。
- **处理方法：** 用固定参数列表运行 `samtools fqidx`。FASTQ 索引（`.fai`）会建在输入文件旁边，或建在 `index_path`（文件夹只读时使用）。区域写作 `name` 或 `name:start-end`（1 起始，闭区间）；`reverse_complement` 对序列取反向互补并反转质量串，`mark_strand` 取 `rc`、`no` 或 `sign`，`continue_on_missing` 让程序在某个区域取不到时继续。输出文件必须不存在。
- **API：** `dnakit.alignment.samtools_extract_fastq_regions(fastq_path, output_path, *, regions=(), region_file=None, line_length=None, reverse_complement=False, mark_strand=None, continue_on_missing=False, index_path=None, samtools_path=None)`.
- **输入：** 一个 FASTQ 文件、一个新的输出路径，以及 `regions` 和/或每行一个区域的 `region_file`。
- **输出：** `status`、`backend`、`version`、`command`、`outputs`（`output`）、`log_tail`、`provenance`，以及 `records`：写出的每条 FASTQ 记录的名称和序列长度。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_extract_fastq_regions

output = Path(tempfile.mkdtemp()) / "regions.fq"
result = samtools_extract_fastq_regions(
    "reads.fastq", output, regions=["read2:5-14", "read3"], reverse_complement=True, mark_strand="sign"
)
print(result["status"], [(row["name"], row["length"]) for row in result["records"]])
print(output.read_text().splitlines()[1])
```

- **示例输出：**

```text
completed [('read2:5-14(-)', 10), ('read3(-)', 40)]
CTGAATAAGG
```

## 2) `dnakit.alignment.samtools_coverage_histogram()`

- **作用：** 把每条参考序列的覆盖度画成文本直方图（`samtools coverage -m`）。
- **处理方法：** 运行 `samtools coverage -m` 并保存各个面板。`bins` 是直方图的列数（`-w`）；samtools 自己使用终端宽度减 40，在没有终端时取决于输出方式，所以这个封装把它固定为 40（80 列终端的值）。`ascii_only` 把方块字符换成 ASCII（`-A`）。`region`、读长与质量过滤以及标志掩码（`--rf`、`--ff`，数字或 `UNMAP,SECONDARY` 这样的标志名）原样传递。
- **API：** `dnakit.alignment.samtools_coverage_histogram(alignment_paths, output_path, *, region=None, bins=40, ascii_only=False, min_read_length=0, min_mapping_quality=0, min_base_quality=0, required_flags=None, excluded_flags=None, max_depth=1000000, min_depth=1, reference_fasta=None, samtools_path=None)`.
- **输入：** 一个或多个已建索引的 BAM/CRAM 文件（使用 region 需要索引）和一个新的输出路径。
- **输出：** 通用字段和 `contigs`：每条参考序列的名称、长度文本和 `statistics`（读段数、覆盖碱基数、覆盖百分比、平均覆盖度、平均碱基与比对质量、直方图分箱宽度和最大分箱）。
- **示例代码：**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_coverage_histogram

work = Path(tempfile.mkdtemp())
subprocess.run(["samtools", "view", "-b", "-o", str(work / "reads.bam"), "reads.sam"], check=True)
result = samtools_coverage_histogram([work / "reads.bam"], work / "histogram.txt", ascii_only=True, bins=20)
for contig in result["contigs"]:
    stats = contig["statistics"]
    print(contig["name"], contig["length"], stats["Number of reads"], stats["Percent covered"])
print((work / "histogram.txt").read_text().splitlines()[1])
```

- **示例输出：**

```text
chr1 600bp 8 55
chr2 400bp 4 50
>  90.00% |:::::::::::         | Number of reads: 8
```

## 3) `dnakit.alignment.samtools_checksum()`

- **作用：** 计算比对文件中序列、质量、名称和标签内容的与顺序无关的校验和，并合并已有的校验和报告（对应工作清单中的 `samtools checksum` 与 `samtools checksum -m`）。
- **处理方法：** 运行 `samtools checksum`。除非设置 `in_order`（`-O`），校验和与记录顺序无关；`check_position`、`check_cigar`、`check_mate` 分别加入 CHR/POS、MAPQ/CIGAR 和 PNEXT/RNEXT/TLEN（`-P -C -M`）。`exclude_flags`、`require_flags`、`flag_mask` 对应 `-F`、`-f`、`-b`；`tags` 对应 `-t`；`max_records` 对应 `-N`。`merge=True`（`-m`）时，输入是本函数先前写出的校验和报告，输出合并后的报告；报告选项必须与生成输入时一致。
- **API：** `dnakit.alignment.samtools_checksum(input_paths, output_path, *, merge=False, exclude_flags=None, require_flags=None, flag_mask=None, no_reverse_complement=False, tags=None, in_order=False, check_position=False, check_cigar=False, check_mate=False, max_records=0, show_qc=False, tab_separated=False, bamseqchksum_format=False, samtools_path=None)`.
- **输入：** SAM/BAM/CRAM 文件（`merge=True` 时为报告文件）和一个新的输出路径。
- **输出：** 通用字段、`merged`，以及 `groups`：每个分组一行，含 `Group`、`QC`、`count` 和校验和 `flag+seq`、`+name`、`+qual`、`+aux`、`combined`。
- **示例代码：**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_checksum

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
by_name = work / "by_name.bam"
subprocess.run(["samtools", "view", "-b", "-o", str(bam), "reads.sam"], check=True)
subprocess.run(["samtools", "sort", "-n", "-o", str(by_name), str(bam)], check=True)
first = samtools_checksum([bam], work / "a.chk")
second = samtools_checksum([by_name], work / "b.chk")
print(first["groups"][0]["count"], first["groups"][0]["combined"] == second["groups"][0]["combined"])
merged = samtools_checksum([work / "a.chk", work / "b.chk"], work / "merged.chk", merge=True)
print(merged["merged"], merged["groups"][0]["count"])
```

- **示例输出：**

```text
12 True
True 24
```

## 4) `dnakit.variants.bcftools_convert_to_vcf()`

- **作用：** 把基因型文本格式转换为 VCF/BCF：gen/sample、hap/sample、hap/legend/sample、23andMe 这类 TSV 基因型表，以及 gVCF 参考块（`bcftools convert -G`、`--hapsample2vcf`、`-H`、`--tsv2vcf`、`--gvcf2vcf`）。
- **处理方法：** `input_format` 决定数据来源，`input_files` 必须按下列顺序给出文件：`gensample`（gen、sample）、`hapsample`（hap、sample）、`haplegendsample`（hap、legend、sample）、`tsv`（一个表；需要 `reference_fasta` 和 `sample_names`；`columns` 同 `-c`，默认 `ID,CHROM,POS,AA`）、`gvcf`（一个 gVCF，其参考块被展开为每个位置一条记录；需要 `reference_fasta`）。`haploid_to_diploid` 把 hap 的单倍型基因型变成二倍体纯合子。`output_type` 取 `v`、`z`、`b` 或 `u`（默认由文件名决定）。
- **API：** `dnakit.variants.bcftools_convert_to_vcf(input_format, input_files, output_path, *, reference_fasta=None, columns=None, sample_names=None, haploid_to_diploid=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** `input_format`、输入文件、一个新的输出路径，以及 `tsv` 和 `gvcf` 所需的参考 FASTA。
- **输出：** 通用字段、`input_format`、`output_type`、`records`（仅纯文本 VCF 输出）和 `summary`：bcftools 打印的计数（TSV 输入时为总行数、跳过行数、写出位点数、缺失基因型数、hom RR、het RA、hom AA、het AA）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_convert_to_vcf

output = Path(tempfile.mkdtemp()) / "sample.vcf"
result = bcftools_convert_to_vcf(
    "tsv", ["raw.tsv"], output, reference_fasta="genome.fa", sample_names=["SAMP1"]
)
print(result["records"], result["summary"]["Missing GTs"])
print([line.split("\t")[9] for line in output.read_text().splitlines() if not line.startswith("#")])
```

- **示例输出：**

```text
3 1
['0/1', '0/0', './.']
```

- **示例代码：**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_convert_to_vcf

work = Path(tempfile.mkdtemp())
subprocess.run(["bcftools", "convert", "-g", str(work / "panel"), "variants.vcf"], check=True, capture_output=True)
result = bcftools_convert_to_vcf(
    "gensample", [work / "panel.gen.gz", work / "panel.samples"], work / "back.vcf"
)
print(result["records"], result["summary"])
expanded = bcftools_convert_to_vcf("gvcf", ["sample.g.vcf"], work / "expanded.vcf", reference_fasta="genome.fa")
print(expanded["records"])
```

- **示例输出：**

```text
10 {'Number of processed rows': 10}
6
```

## 5) `dnakit.variants.bcftools_index_records()`

- **作用：** 从 VCF/BCF 的索引中读取每个 contig 的记录数或总记录数（`bcftools index -s` 与 `bcftools index -n`）。
- **处理方法：** 运行 `bcftools index -s`（`per_contig=True`，`include_empty_contigs` 时加 `-a`）或 `bcftools index -n`。文件需要已有 CSI 或 TBI 索引，不会重写任何文件。
- **API：** `dnakit.variants.bcftools_index_records(vcf_path, *, per_contig=True, include_empty_contigs=False, bcftools_path=None)`.
- **输入：** 带索引的 bgzip 压缩 VCF 或 BCF。
- **输出：** 通用字段和 `total_records`；`per_contig` 时还有 `contigs`，每个 contig 一行，含 `contig`、`length`（头部没有长度时为 None）和 `records`。
- **示例代码：**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_index_records

work = Path(tempfile.mkdtemp())
compressed = work / "variants.vcf.gz"
subprocess.run(["bcftools", "view", "-Oz", "-o", str(compressed), "variants.vcf"], check=True)
subprocess.run(["bcftools", "index", str(compressed)], check=True)
result = bcftools_index_records(compressed)
print(result["contigs"], result["total_records"])
print(bcftools_index_records(compressed, per_contig=False)["total_records"])
```

- **示例输出：**

```text
[{'contig': 'chr1', 'length': 600, 'records': 10}] 10
10
```

## 6) `dnakit.variants.bcftools_norm_variants()`

- **作用：** 规范化 VCF/BCF：左对齐 indel、拆分或合并多等位位点、拆解复杂变异、去除重复记录并检查 REF 等位基因（`bcftools norm`）。它补全 `bcftools_normalize_vcf`（后者只提供 `-f` 与 `-m`）。
- **处理方法：** 一个函数覆盖 `bcftools norm` 的所有模式。`reference_fasta`（`-f`）启用左对齐与 REF 检查；`check_ref` 取 `e`、`w`、`x` 或 `s`（`-c`）；`multiallelics` 取 `-snps`、`-indels`、`-both`、`-any`（拆分）或 `+snps`、`+indels`、`+both`、`+any`（合并）（`-m`）；`atomize` 即 `-a`，`atom_overlaps` 取 `*` 或 `.`；`remove_duplicates` 取 `snps`、`indels`、`both`、`all`、`exact`（`-d`）；`no_realign` 即 `-N`；`gff_annotation` 在转录本中右对齐（`-g`）。至少要请求一个动作。
- **API：** `dnakit.variants.bcftools_norm_variants(vcf_path, output_path, *, reference_fasta=None, multiallelics=None, atomize=False, atom_overlaps='*', check_ref=None, remove_duplicates=None, no_realign=False, gff_annotation=None, old_record_tag=None, strict_filter=False, site_window=1000, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 一个 VCF/BCF（只有 `regions` 才需要 bgzip 和索引）、一个新的输出路径和上述选项。
- **输出：** 通用字段、`output_type`、`records`（仅纯文本 VCF 输出）和 `summary`（标准错误上打印的计数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_norm_variants

work = Path(tempfile.mkdtemp())
result = bcftools_norm_variants("multi.vcf", work / "split.vcf", multiallelics="-any")
print(result["records"])
atoms = bcftools_norm_variants("multi.vcf", work / "atoms.vcf", atomize=True, remove_duplicates="exact")
print(atoms["records"])
for line in (work / "atoms.vcf").read_text().splitlines():
    if not line.startswith("#"):
        print(*line.split("\t")[1:5])
```

- **示例输出：**

```text
5
5
10 . A C,*
10 . A G,*
20 . C G
21 . A C
30 . A C
```

## 7) `dnakit.variants.bcftools_polysomy()`

- **作用：** 根据 B 等位基因频率检测整条染色体的拷贝数异常（`bcftools polysomy`）。
- **处理方法：** 运行 `bcftools polysomy`。VCF/BCF 需要 `FORMAT/BAF`。对每条染色体分别拟合 2、3、4 个拷贝的 BAF 分布峰，拟合最好的拷贝数即估计值。`peak_size`、`cn_penalty`、`fit_threshold`、`min_fraction`、`peak_symmetry`、`include_aa` 对应 `-b -c -f -m -p -i`。该命令只存在于带 GSL 编译的 bcftools（`--enable-libgsl`，conda 版本带有）中，否则抛出 `BackendUnavailableError`。绘图脚本 `dist.py` 会被写出但不会运行。
- **API：** `dnakit.variants.bcftools_polysomy(vcf_path, output_dir, *, sample=None, regions=None, targets=None, peak_size=0.1, cn_penalty=0.7, fit_threshold=3.3, min_fraction=0.1, peak_symmetry=0.5, include_aa=False, bcftools_path=None)`.
- **输入：** 带 `FORMAT/BAF` 的 VCF/BCF 和一个新的输出文件夹。
- **输出：** 通用字段、`files`（含 BAF 分布、拟合函数和 `CN` 行的 `dist.dat`，以及绘图脚本 `dist.py`）和 `copy_numbers`：染色体、估计拷贝数和绝对拟合偏差。
- **示例代码：**

```python
import random
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_polysomy

rng = random.Random(3)
lines = [
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000000>", "##contig=<ID=2,length=1000000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]
for chrom, centers in (("1", [0.0, 0.5, 1.0]), ("2", [0.0, 1 / 3, 2 / 3, 1.0])):  # two and three copies
    for index in range(2500):
        center = rng.choice(centers)
        value = abs(rng.gauss(0, 0.02)) if center == 0 else 1 - abs(rng.gauss(0, 0.02)) if center == 1 else rng.gauss(center, 0.025)
        genotype = "0/0" if center == 0 else "1/1" if center == 1 else "0/1"
        lines.append(f"{chrom}\t{1000 + index * 100}\t.\tA\tG\t.\t.\t.\tGT:BAF\t{genotype}:{min(1.0, max(0.0, value)):.4f}")
work = Path(tempfile.mkdtemp())
(work / "baf.vcf").write_text("\n".join(lines) + "\n")
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "baf.vcf.gz"), str(work / "baf.vcf")], check=True)
subprocess.run(["bcftools", "index", str(work / "baf.vcf.gz")], check=True)
result = bcftools_polysomy(work / "baf.vcf.gz", work / "polysomy", sample="S1")
print([(row["chromosome"], row["copy_number"]) for row in result["copy_numbers"]], sorted(result["files"]))
```

- **示例输出：**

```text
[('1', 2.0), ('2', 3.0)] ['dist.dat', 'dist.py']
```
