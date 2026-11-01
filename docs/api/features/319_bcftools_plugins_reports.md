# bcftools 插件封装：统计、报告与注释

这一组函数位于 `dnakit.variants`（模块 `bcftools_plugins_p06_cli_functions`），封装外部程序 [bcftools](https://github.com/samtools/bcftools) 的插件（`bcftools +plugin`，测试版本为 1.24）：等位基因与基因型统计、注释、基因型比较和子集提取。每个函数校验参数，用固定的参数列表且不经过 shell 运行一个插件，分别捕获标准输出和标准错误，并返回可 JSON 序列化的字典。打印报告的插件会把报告保存到 `output_path` 并解析成 `report`（注释行、列名和数值行）。修改记录的插件把 VCF/BCF 写到 `output_path`（`output_type` 取 `v`、`z`、`b` 或 `u`，默认由文件名决定），纯文本 VCF 输出时返回记录数。已存在的输出不会被覆盖，运行失败时不会留下任何文件。

bcftools 是带插件的外部程序：请自行安装（例如 `conda install -c bioconda bcftools`，或从 https://www.htslib.org/download/ 下载发行版），并放到 PATH 中，或传入 `bcftools_path`。插件通过环境变量 `BCFTOOLS_PLUGINS` 查找；未设置时使用安装目录旁的 `libexec/bcftools` 文件夹。找不到程序时抛出带安装提示的 `BackendUnavailableError`；程序报错退出时抛出 `BackendExecutionError`（代码 `HTSLIB_EXECUTION_FAILED`），并附带日志尾部。

常用的 bcftools 选项是 `regions`（`-r`，需要索引）、`include`（`-i`）和 `exclude`（`-e`）。不提供的内容：`bcftools +frameshifts`（1.24 中它对每个 indel 都写 `OOF=-1`，因为 htslib 现在把 indel 标记为 `VCF_INS`/`VCF_DEL`，插件的类型判断 `var[i].type != VCF_INDEL` 永远不成立），以及 indel-stats 的 `--nvaf` 选项。所有结果都包含 `status`、`backend`、`operation`、`version`、`command`、`outputs`、`log_tail` 和 `provenance`；下面各节只列出额外字段。

## 示例数据

示例从工作目录读取这些文件，并把输出写入临时目录。`variants.vcf` 有六个样本、十个位点（其中有一个缺失和一个插入）；`trio.vcf` 有孩子 C、父亲 F、母亲 M 和无关样本 U，含 `GT`、`PL`、`AD`、`DP`、`GQ`、`QS`、`SP`；`trio_indel.vcf` 与之相同，但每隔两个位点有一个 indel。

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

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	60	var1	G	T	40	PASS	DP=50;AF=0.417	GT	0/1	0/0	1/1	0/1	0/0	0/1
chr1	120	var2	C	T	41	PASS	DP=55;AF=0.333	GT	0/0	0/1	0/1	1/1	0/0	0/0
chr1	180	var3	C	A	42	PASS	DP=60;AF=0.5	GT	1/1	0/1	0/0	0/1	0/1	./.
chr1	240	var4	ACT	A	43	PASS	DP=65;AF=0.417	GT	0/1	0/1	0/0	0/0	1/1	0/1
chr1	300	var5	T	C	44	PASS	DP=70;AF=0.417	GT	0/0	0/0	0/1	0/1	0/1	1/1
chr1	360	var6	A	T	45	PASS	DP=75;AF=0.417	GT	0/1	1/1	0/1	0/0	0/0	0/1
chr1	420	var7	G	GCG	46	PASS	DP=80;AF=0.417	GT	0/0	0/1	0/0	0/1	1/1	0/1
chr1	480	var8	T	C	47	PASS	DP=85;AF=0.5	GT	0/1	0/0	1/1	1/1	0/1	0/0
chr1	520	var9	G	C	48	PASS	DP=90;AF=0.25	GT	0/0	0/0	0/1	0/0	0/1	0/1
chr1	570	var10	A	C	49	PASS	DP=95;AF=0.417	GT	1/1	0/1	0/1	0/0	0/0	0/1
```

**示例文件：** `trio.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-scaled quality sum">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Strand bias">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	F	M	U
chr1	40	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
```

**示例文件：** `trio_indel.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-scaled quality sum">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Strand bias">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	F	M	U
chr1	40	.	T	TCG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	TTT	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	ACG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	TTG	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	GCG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
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

**示例文件：** `sites.txt`

```text
chr1	60	G	T
chr1	120	C	G
```

## 1) `dnakit.variants.bcftools_ad_bias()`

- **作用：** 用 Fisher 精确检验找出某个样本与其对照之间 ALT 等位基因比例（`FORMAT/AD`）不同的位点（`bcftools +ad-bias`）。
- **处理方法：** 把 `[样本, 对照]` 对写入临时配对文件并运行插件。p 值低于 `threshold`（`-t`）的位点输出一行 `FT`，含两个样本的等位基因深度。`min_alt_depth` 与 `min_depth` 对应 `-a` 与 `-d`，`variant_type` 取 `snp` 或 `indel`（`-v`），`extra_columns` 是追加到每行的 `bcftools query` 格式（`-f`）。`clean_vcf=True`（`-c`）时输出的是去掉未通过阈值的位点和 ALT 等位基因的 VCF。
- **API：** `dnakit.variants.bcftools_ad_bias(vcf_path, output_path, *, sample_pairs, threshold=0.001, min_alt_depth=1, min_depth=0, variant_type=None, extra_columns=None, clean_vcf=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/AD` 的 VCF/BCF、一个新的输出路径和样本对。
- **输出：** 通用字段和 `report`（注释行、列名、数值 `rows`）；`clean_vcf` 模式下改为 `records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_ad_bias

output = Path(tempfile.mkdtemp()) / "adbias.txt"
result = bcftools_ad_bias("trio.vcf", output, sample_pairs=[["F", "M"]])
for row in result["report"]["rows"][:3]:
    print(row[:7], f"{row[-1]:.2e}")
```

- **示例输出：**

```text
['FT', 'F', 'M', 'chr1', 40, 'T', 'A'] 5.83e-06
['FT', 'F', 'M', 'chr1', 120, 'C', 'G'] 1.69e-17
['FT', 'F', 'M', 'chr1', 280, 'A', 'C'] 5.83e-06
```

## 2) `dnakit.variants.bcftools_af_distribution()`

- **作用：** 统计哈迪-温伯格平衡下的等位基因频率偏差分布和基因型概率分布（`bcftools +af-dist`）。
- **处理方法：** 只使用频率已知位点上的非参考基因型；杂合的概率为 `2*AF*(1-AF)`，纯合替代为 `AF**2`。`deviation_bins` 与 `probability_bins` 是分箱边界（`-d`、`-p`，默认 0、0.1、…、1），`af_tag` 是存放频率的 INFO 标签（`-t`），`per_sample` 追加每个样本的几何平均概率（`-s`）。`regions`、`include`、`exclude` 即 bcftools 的位点筛选选项 `-r`、`-i`、`-e`（`regions` 需要已建索引的文件）。
- **API：** `dnakit.variants.bcftools_af_distribution(vcf_path, output_path, *, deviation_bins=None, probability_bins=None, per_sample=False, af_tag='AF', regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** 带 INFO 等位基因频率标签和 `FORMAT/GT` 的 VCF/BCF，以及一个新的输出路径。
- **输出：** 通用字段和 `report`；`PROB_DIST` 与 `DEV_DIST` 行依次为分箱起点、分箱终点和计数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_af_distribution

output = Path(tempfile.mkdtemp()) / "afdist.txt"
result = bcftools_af_distribution("variants.vcf", output, deviation_bins=[0, 0.25, 0.5, 1], probability_bins=[0, 0.5, 1])
for row in result["report"]["rows"]:
    print(row)
```

- **示例输出：**

```text
['PROB_DIST', 0.0, 0.5, 33]
['PROB_DIST', 0.5, 1.0, 5]
['DEV_DIST', 0.0, 0.25, 10]
['DEV_DIST', 0.25, 0.5, 0]
['DEV_DIST', 0.5, 1.0, 0]
```

## 3) `dnakit.variants.bcftools_allele_length_counts()`

- **作用：** 统计 REF 和 ALT 等位基因各种长度出现的次数（`bcftools +allele-length`）。
- **处理方法：** 运行插件并解析其表格。`regions`、`include`、`exclude` 即 bcftools 的位点筛选选项 `-r`、`-i`、`-e`。
- **API：** `dnakit.variants.bcftools_allele_length_counts(vcf_path, output_path, *, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段和 `report`，列为 `LENGTH`、`REF`、`ALT`、`REF+ALT`、`REF+ALT WITH NON-BASE NUCLEOTIDES`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_allele_length_counts

output = Path(tempfile.mkdtemp()) / "lengths.txt"
result = bcftools_allele_length_counts("variants.vcf", output)
print(result["report"]["columns"][:3])
print([row[:3] for row in result["report"]["rows"] if any(row[1:])])
```

- **示例输出：**

```text
['LENGTH', 'REF', 'ALT']
[[1, 9, 9], [2, 0, 0], [3, 1, 1], [4, 0, 0], ['', '', '']]
```

## 4) `dnakit.variants.bcftools_check_ploidy()`

- **作用：** 报告每个样本在各条染色体上的倍性，以及它在各位点之间是否一致（`bcftools +check-ploidy`）。
- **处理方法：** 运行插件。除非设置 `use_missing`（`-m`），缺失和半缺失基因型会被忽略。
- **API：** `dnakit.variants.bcftools_check_ploidy(vcf_path, output_path, *, use_missing=False, regions=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段和 `report`：每个样本、每个区域一行，含样本、染色体、区域起点、区域终点和倍性。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_check_ploidy

output = Path(tempfile.mkdtemp()) / "ploidy.txt"
result = bcftools_check_ploidy("variants.vcf", output)
print(result["report"]["columns"])
print(result["report"]["rows"][:2], len(result["report"]["rows"]))
```

- **示例输出：**

```text
['Sample', 'Chromosome', 'Region Start', 'Region End', 'Ploidy']
[['S1', 'chr1', 60, 570, 2], ['S2', 'chr1', 60, 570, 2]] 6
```

## 5) `dnakit.variants.bcftools_check_sparsity()`

- **作用：** 列出在某个区域或染色体上基因型数量少于最小值的样本（`bcftools +check-sparsity`）。
- **处理方法：** 用 `-n min_markers` 运行插件；没有 `regions` 或 `regions_file`（`-r`、`-R`）时检查每条染色体。报告为空表示每个样本都有足够的标记。
- **API：** `dnakit.variants.bcftools_check_sparsity(vcf_path, output_path, *, min_markers=1, regions=None, regions_file=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段，以及保存在 `output_path` 的样本列表。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_check_sparsity

work = Path(tempfile.mkdtemp())
bcftools_check_sparsity("variants.vcf", work / "all.txt", min_markers=1)
bcftools_check_sparsity("variants.vcf", work / "dense.txt", min_markers=10)
print(repr((work / "all.txt").read_text()))
print((work / "dense.txt").read_text().split())
```

- **示例输出：**

```text
''
['chr1', 'S6']
```

## 6) `dnakit.variants.bcftools_genotype_dosage()`

- **作用：** 按位点和样本输出基因型剂量（0 到 2，缺失为 -1），来源为 PL、GL 或 GT（`bcftools +dosage`）。
- **处理方法：** 运行插件；`tags` 是剂量所取 FORMAT 标签的优先顺序（`-t`）。
- **API：** `dnakit.variants.bcftools_genotype_dosage(vcf_path, output_path, *, tags=('PL', 'GL', 'GT'), regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段和 `report`：每个位点一行，含 CHROM、POS、REF、ALT 和每个样本的剂量。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_genotype_dosage

output = Path(tempfile.mkdtemp()) / "dosage.txt"
result = bcftools_genotype_dosage("variants.vcf", output, tags=["GT"])
for row in result["report"]["rows"][:3]:
    print(row)
```

- **示例输出：**

```text
['chr1', 60, 'G', 'T', 1.0, 0.0, 2.0, 1.0, 0.0, 1.0]
['chr1', 120, 'C', 'T', 0.0, 1.0, 1.0, 2.0, 0.0, 0.0]
['chr1', 180, 'C', 'A', 2.0, 1.0, 0.0, 1.0, 1.0, -1.0]
```

## 7) `dnakit.variants.bcftools_genotype_intersections()`

- **作用：** 统计样本所有子集上的基因型交集数量（`bcftools +GTisec`）。
- **处理方法：** 运行插件。`include_missing` 追加每个样本的缺失基因型数（`-m`），`annotate_subsets` 在每行标出该子集的样本名（`-v`），`human_readable` 按样本排序输出（`-H`，隐含 `annotate_subsets`）。
- **API：** `dnakit.variants.bcftools_genotype_intersections(vcf_path, output_path, *, include_missing=False, annotate_subsets=False, human_readable=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** 多样本 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段和 `report`：每个样本子集一行，给出恰好这些样本共享一个基因型的位点数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_genotype_intersections

output = Path(tempfile.mkdtemp()) / "gtisec.txt"
result = bcftools_genotype_intersections("trio.vcf", output, annotate_subsets=True)
for row in result["report"]["rows"][:5]:
    print(row)
```

- **示例输出：**

```text
['@SMPS U M F C']
[4, 'U']
[2, 'M']
[2, 'F']
[10, 'C']
```

## 8) `dnakit.variants.bcftools_indel_stats()`

- **作用：** 按样本统计 indel，或统计 PED 文件中三联体的新生（de novo）indel（`bcftools +indel-stats`）。
- **处理方法：** 运行插件。`max_length` 是统计的最长 indel（`--max-len`），`csq_tag` 是存放 VEP/BCSQ 后果的 INFO 标签，用于区分移码与框内 indel（`-c`）。给出 `ped_path`（`-p`）时只统计新生 indel；`alt_to_ref_dnm` 还把 0/1 + 1/1 -> 0/0 这类基因型错误计为新生。`include` 可以使用 `GQ>{10,20}` 这种花括号范围同时扫描多个阈值。插件需要 `FORMAT/AD`。bcftools 1.24 文档中的 `--nvaf` 选项会被插件拒绝（它按小数而不是分箱数解析），所以不提供。
- **API：** `dnakit.variants.bcftools_indel_stats(vcf_path, output_path, *, ped_path=None, csq_tag='CSQ', max_length=20, alt_to_ref_dnm=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 和 `FORMAT/AD` 的 VCF/BCF，以及一个新的输出路径。
- **输出：** 通用字段和 `report`，含 `SN*`、`FRM`、`VAF`、`NFRAC` 等行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_indel_stats

output = Path(tempfile.mkdtemp()) / "indels.txt"
result = bcftools_indel_stats("trio_indel.vcf", output, include="GQ>{10,20}")
for row in result["report"]["rows"]:
    if str(row[0]).startswith("SN"):
        print(row[:6])
```

- **示例输出：**

```text
['SN0', 4, 5, 5, 10, 3]
['SN1', 4, 5, 5, 10, 3]
```

## 9) `dnakit.variants.bcftools_vcf_to_table()`

- **作用：** 把变异打印成文本表格：变异表、INFO 表、基因型表和后果表（`bcftools +vcf2table`）。
- **处理方法：** 运行插件并保存其文本。`hide` 列出要省略的表：`HOM_REF`/`RR`、`HET`/`AR`、`NO_CALL`/`MISSING`、`CSQ`/`VEP`、`SPLICEAI`、`ANN`/`SNPEFF`、`LOF`、`VC`、`INFO`、`GT`、`GTTYPES`、`URL`。
- **API：** `dnakit.variants.bcftools_vcf_to_table(vcf_path, output_path, *, hide=None, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`report`（`variants`：变异块的数量）和 `text_preview`（文本的前 2000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_vcf_to_table

output = Path(tempfile.mkdtemp()) / "table.txt"
result = bcftools_vcf_to_table("variants.vcf", output, hide=["HOM_REF", "NO_CALL", "INFO"], include="POS<=120")
print(result["report"]["variants"])
print(output.read_text().splitlines()[0])
print(output.read_text().splitlines()[3])
```

- **示例输出：**

```text
2
<<< chr1:60:G (n. 1)
+--------+-------+
```

## 10) `dnakit.variants.bcftools_variant_read_frequency_score()`

- **作用：** 根据许多未患病样本的 BAM/CRAM 文件评估每个位点的噪声程度（变异读段频率评分，`bcftools +vrfs`）。
- **处理方法：** 评分模式：`alignment_list` 是每行一个 BAM/CRAM 路径的文本文件（`-a`），`reference_fasta` 是参考序列（`-f`），`sites_file` 是由染色体、位置、ref、alt 组成的制表符分隔文件（`-s`）；`min_depth` 与 `vaf_bins` 对应 `-d` 与 `-n`；`batch` 取 `I/N`，只运行 N 批中的第 I 批（`-b`）。合并模式：`merge_batch_files` 列出要合并的批输出（`-M`），此时不能给出评分参数。
- **API：** `dnakit.variants.bcftools_variant_read_frequency_score(output_path, *, alignment_list=None, reference_fasta=None, sites_file=None, min_depth=10, vaf_bins=20, batch=None, merge_batch_files=None, bcftools_path=None)`.
- **输入：** 列在文本文件中的比对文件、已建索引的参考 FASTA 和位点文件；或要合并的批输出。
- **输出：** 通用字段和 `report`，评分表保存在 `output_path`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variant_read_frequency_score

import subprocess

work = Path(tempfile.mkdtemp())
subprocess.run(["samtools", "view", "-b", "-o", str(work / "reads.bam"), "reads.sam"], check=True)
(work / "bams.txt").write_text(f"{work / 'reads.bam'}\n")
result = bcftools_variant_read_frequency_score(
    work / "scores.txt", alignment_list=work / "bams.txt", reference_fasta="genome.fa", sites_file="sites.txt", min_depth=1
)
print(result["status"], len(result["report"]["rows"]))
print([row for row in result["report"]["rows"] if row and row[0] == "SC"][:1])
```

- **示例输出：**

```text
completed 4
[]
```

## 11) `dnakit.variants.bcftools_add_variantkey()`

- **作用：** 添加 INFO 字段 `VKX`（十六进制 64 位 VariantKey）和 `RSX`（十六进制 rsID 编号）（`bcftools +add-variantkey`）。
- **处理方法：** 运行插件并写出带注释的 VCF/BCF；`output_type` 取 `v`、`z`、`b` 或 `u`（默认由文件名决定）。
- **API：** `dnakit.variants.bcftools_add_variantkey(vcf_path, output_path, *, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_add_variantkey

output = Path(tempfile.mkdtemp()) / "keys.vcf"
result = bcftools_add_variantkey("variants.vcf", output)
print(result["records"])
print(output.read_text().splitlines()[-1].split("\t")[7])
```

- **示例输出：**

```text
10
DP=95;AF=0.417;VKX=0800011c88880000;RSX=00000000
```

## 12) `dnakit.variants.bcftools_variantkey_tables()`

- **作用：** 生成未排序的十六进制 VariantKey 查找表（`bcftools +variantkey-hex`）。
- **处理方法：** 在新文件夹 `output_dir` 中运行插件，它写出 `vkrs.unsorted.hex`（VariantKey 到 rsID）、`rsvk.unsorted.hex`（rsID 到 VariantKey）和 `nrvk.unsorted.tsv`（不可逆的 VariantKey）。
- **API：** `dnakit.variants.bcftools_variantkey_tables(vcf_path, output_dir, *, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出文件夹。
- **输出：** 通用字段、`files` 和 `counters`（bcftools 打印的数字：VariantKeys 与 Non-reversible VariantKeys）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variantkey_tables

result = bcftools_variantkey_tables("variants.vcf", Path(tempfile.mkdtemp()) / "keys")
print(result["files"], result["counters"])
```

- **示例输出：**

```text
['nrvk.unsorted.tsv', 'rsvk.unsorted.hex', 'vkrs.unsorted.hex'] {'VariantKeys': 10, 'Non-reversible VariantKeys': 0}
```

## 13) `dnakit.variants.bcftools_contrast_groups()`

- **作用：** 在两组样本之间做病例/对照关联检验，并检测新的等位基因或基因型（`bcftools +contrast`）。
- **处理方法：** `control_samples` 是背景组（`-0`），`case_samples` 是预期出现新等位基因或基因型的组（`-1`）。`annotations` 选择写入的 INFO 标签：`PASSOC`（Fisher 检验，REF 对非 REF）、`FASSOC`（对照与病例中非 REF 的比例）、`NASSOC`（等位基因计数）、`NOVELAL`、`NOVELGT`（默认 PASSOC、FASSOC、NOVELAL）。`max_allele_freq`（`-f`）检验稀有等位基因的富集：(0, 1) 之间的值是 ALT 频率，整数是 ALT 计数。
- **API：** `dnakit.variants.bcftools_contrast_groups(vcf_path, output_path, *, control_samples, case_samples, annotations=None, max_allele_freq=None, force_samples=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 多样本 VCF/BCF、两组样本列表和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_contrast_groups

output = Path(tempfile.mkdtemp()) / "contrast.vcf"
result = bcftools_contrast_groups("trio.vcf", output, control_samples=["F", "M"], case_samples=["C"])
print(result["records"])
for line in output.read_text().splitlines()[-40:-37]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **示例输出：**

```text
14
```

## 14) `dnakit.variants.bcftools_fill_from_fasta()`

- **作用：** 用 FASTA 文件填充 REF 列或某个 INFO 标签（`bcftools +fill-from-fasta`）。
- **处理方法：** `column` 取 `REF`（修复为 N 或错误的 REF 等位基因）或 `AA` 这类 INFO 标签（祖先等位基因）。FASTA 必须已建索引（`samtools faidx`）。`header_lines_path` 存放新标签的 `##INFO` 行（`-h`）。`include` / `exclude`（插件选项 `-i`、`-e`）限制被注释的记录，`replace_non_acgtn` 把其他字符替换为 N（`-N`）。
- **API：** `dnakit.variants.bcftools_fill_from_fasta(vcf_path, output_path, *, fasta_path, column, header_lines_path=None, include=None, exclude=None, replace_non_acgtn=False, output_type=None, bcftools_path=None)`.
- **输入：** VCF/BCF、已建索引的 FASTA、列名和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fill_from_fasta

import subprocess

subprocess.run(["samtools", "faidx", "genome.fa"], check=True)
Path("aa.hdr").write_text('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n')
output = Path(tempfile.mkdtemp()) / "aa.vcf"
result = bcftools_fill_from_fasta(
    "variants.vcf", output, fasta_path="genome.fa", column="AA", header_lines_path="aa.hdr", include='TYPE="snp"'
)
print(result["records"])
for line in output.read_text().splitlines()[-10:-7]:
    print(line.split("\t")[1], line.split("\t")[3], line.split("\t")[7])
```

- **示例输出：**

```text
10
60 G DP=50;AF=0.417;AA=G
120 C DP=55;AF=0.333;AA=C
180 C DP=60;AF=0.5;AA=C
```

## 15) `dnakit.variants.bcftools_fill_tags()`

- **作用：** 根据基因型计算 INFO 与 FORMAT 标签（AF、AC、AN、MAF、HWE、ExcHet、NS、F_MISSING、VAF 以及自定义表达式）（`bcftools +fill-tags`）。
- **处理方法：** `tags`（`-t`）列出要填充的内容：`AF`、`AC`、`AC_Hom`、`AC_Het`、`AC_Hemi`、`AN`、`ExcHet`、`HWE`、`MAF`、`NS`、`F_MISSING`、`END`、`TYPE`、`VAF`、`VAF1`、`all`，或 `DP:1=int(sum(FORMAT/DP))` 这样的自定义表达式。`populations_file`（`-S`）列出样本及逗号分隔的群体名称，用于按群体计算标签。`drop_missing` 不把半缺失基因型计为半合子（`-d`）。
- **API：** `dnakit.variants.bcftools_fill_tags(vcf_path, output_path, *, tags, populations_file=None, drop_missing=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** VCF/BCF、标签列表和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fill_tags

output = Path(tempfile.mkdtemp()) / "tags.vcf"
result = bcftools_fill_tags("variants.vcf", output, tags=["AN", "AC", "NS", "F_MISSING"])
print(result["records"])
for line in output.read_text().splitlines()[-10:-8]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **示例输出：**

```text
10
60 DP=50;AF=0.417;F_MISSING=0;NS=6;AN=12;AC=5
120 DP=55;AF=0.333;F_MISSING=0;NS=6;AN=12;AC=4
```

## 16) `dnakit.variants.bcftools_exclusive_genotype_sites()`

- **作用：** 只输出所列样本都共享同一个基因型、且其他样本都没有该基因型的位点（`bcftools +GTsubset`）。
- **处理方法：** 用逗号连接的 `samples`（`-s`）运行插件。
- **API：** `dnakit.variants.bcftools_exclusive_genotype_sites(vcf_path, output_path, *, samples, regions=None, output_type=None, bcftools_path=None)`.
- **输入：** 多样本 VCF/BCF、样本列表和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_exclusive_genotype_sites

output = Path(tempfile.mkdtemp()) / "exclusive.vcf"
result = bcftools_exclusive_genotype_sites("variants.vcf", output, samples=["S1", "S2"])
print(result["records"])
for line in output.read_text().splitlines():
    if not line.startswith("#"):
        print(line.split("\t")[1], line.split("\t")[9:15])
```

- **示例输出：**

```text
1
300 ['0/0', '0/0', '0/1', '0/1', '0/1', '1/1']
```

## 17) `dnakit.variants.bcftools_isec_genotypes()`

- **作用：** 比较两个文件，把第一个文件中与第二个文件不同的基因型设为缺失（`bcftools +isecGT`）。
- **处理方法：** 对两个文件运行插件；两者都必须是 bgzip 压缩的 VCF 或 BCF 并已建索引，样本和位点按名称和位置匹配。
- **API：** `dnakit.variants.bcftools_isec_genotypes(vcf_path, other_vcf_path, output_path, *, regions=None, output_type=None, bcftools_path=None)`.
- **输入：** 两个已建索引的 VCF/BCF 文件和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_isec_genotypes

import subprocess

work = Path(tempfile.mkdtemp())
text = Path("trio.vcf").read_text().replace("0/1:60,0,60:15,15", "1/1:60,0,60:15,15", 3)
(work / "other.vcf").write_text(text)
for name, source in (("first", "trio.vcf"), ("other", str(work / "other.vcf"))):
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / f"{name}.vcf.gz"), source], check=True)
    subprocess.run(["bcftools", "index", str(work / f"{name}.vcf.gz")], check=True)
output = work / "isec.vcf"
result = bcftools_isec_genotypes(work / "first.vcf.gz", work / "other.vcf.gz", output)
rows = [line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")]
print(result["records"], sum(call.startswith("./.") for row in rows for call in row[9:]))
```

- **示例输出：**

```text
14 3
```
