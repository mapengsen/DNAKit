# bcftools 插件封装：VCF 编辑、拆分与三联体分析

这一组函数位于 `dnakit.variants`（模块 `bcftools_plugins_p06_cli_functions` 与 `bcftools_trio_p06_cli_functions`），封装外部程序 [bcftools](https://github.com/samtools/bcftools) 中修改或拆分 VCF/BCF 的插件（`bcftools +plugin`，测试版本为 1.24），以及三联体插件：孟德尔一致性、新生突变筛查、传递统计、相位切换率和拷贝数变异的亲本来源。每个函数校验参数，用固定的参数列表且不经过 shell 运行一个插件，分别捕获标准输出和标准错误，并返回可 JSON 序列化的字典。修改记录的插件把 VCF/BCF 写到 `output_path`（`output_type` 取 `v`、`z`、`b` 或 `u`，默认由文件名决定）；写出多个文件的插件使用新文件夹；报告类插件保存报告并解析为 `report`。已存在的输出不会被覆盖，运行失败时不会留下任何文件。

bcftools 是带插件的外部程序：请自行安装（例如 `conda install -c bioconda bcftools`，或从 https://www.htslib.org/download/ 下载发行版），并放到 PATH 中，或传入 `bcftools_path`。插件通过 `BCFTOOLS_PLUGINS` 或安装目录旁的 `libexec/bcftools` 查找。找不到程序时抛出带安装提示的 `BackendUnavailableError`；程序报错退出时抛出 `BackendExecutionError`（代码 `HTSLIB_EXECUTION_FAILED`），并附带日志尾部。三联体用 PED 文件（`ped_path`）或样本名列表（`trio`）给出。`bcftools +trio-dnm3` 的 `dmm` 模型需要 bcftools 1.24 或更新版本。bcftools 支持的地方可使用通用选项 `regions`（`-r`）、`include`（`-i`）、`exclude`（`-e`）。所有结果都包含 `status`、`backend`、`operation`、`version`、`command`、`outputs`、`log_tail` 和 `provenance`；下面各节只列出额外字段。

## 示例数据

示例从工作目录读取这些文件，并把输出写入临时目录。`trio.vcf` 有孩子 C、父亲 F、母亲 M 和无关样本 U（十四个位点，循环七种基因型模式，其中两种违反孟德尔遗传）；`trio_phased.vcf` 相同但基因型已定相；`trio.ped` 给出三联体；`csq.vcf` 带有 VEP 风格的 `CSQ` 标签。

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

**示例文件：** `trio_phased.vcf`

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
chr1	40	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0
```

**示例文件：** `csq.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Feature">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	G	.	.	CSQ=G|missense_variant|MODERATE|GENE1|T1,G|synonymous_variant|LOW|GENE1|T2
chr1	20	.	C	T	.	.	CSQ=T|stop_gained|HIGH|GENE2|T3
chr1	30	.	G	A	.	.	.
```

**示例文件：** `trio.ped`

```text
fam1	C	F	M	1	0
```

## 1) `dnakit.variants.bcftools_prune_sites()`

- **作用：** 按连锁不平衡（r2、Lewontin D、Ragsdale D）或窗口内位点密度修剪或注释位点（`bcftools +prune`）。
- **处理方法：** `statistic` 取 `r2`、`LD`、`RD` 或 `count`；`max_value` 删除高于它的位点（`-m`），给出 `set_filter`（`-f`）时改为用该 FILTER 标记。`annotate` 添加最大的 `r2`/`LD`/`RD` 值及其位置，或 `count`（`-a`）。`window` 取 `100kb`、`1000`（位点数）等（`-w`）。`sites_per_window` 每个窗口最多保留 N 个位点（`-n`），按 `sites_per_window_mode`（`maxAF`、`1st`、`rand`）挑选。`af_tag`、`keep_sites`、`random_seed`、`randomize_missing` 为对应选项。`max_value`、`sites_per_window`、`annotate` 至少给一个。
- **API：** `dnakit.variants.bcftools_prune_sites(vcf_path, output_path, *, statistic='r2', max_value=None, annotate=None, window='100kb', sites_per_window=None, sites_per_window_mode='maxAF', af_tag=None, set_filter=None, keep_sites=False, random_seed=None, randomize_missing=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_prune_sites

work = Path(tempfile.mkdtemp())
pruned = bcftools_prune_sites("variants.vcf", work / "pruned.vcf", max_value=0.3, window="1000")
thinned = bcftools_prune_sites("variants.vcf", work / "thin.vcf", sites_per_window=1, window="150bp", sites_per_window_mode="1st")
annotated = bcftools_prune_sites("variants.vcf", work / "r2.vcf", annotate=["r2"], window="1", set_filter=".")
print(pruned["records"], thinned["records"], annotated["records"])
print(annotated["outputs"]["output"].endswith("r2.vcf"), "R2=" in (work / "r2.vcf").read_text())
```

- **示例输出：**

```text
4 4 10
True True
```

## 2) `dnakit.variants.bcftools_remove_overlaps()`

- **作用：** 删除、列出或标记重叠和重复的变异（`bcftools +remove-overlaps`）。
- **处理方法：** `mark` 选择位点（`-m`）：`overlap`、`dup`（所有重复位点）或 `min(QUAL)`（QUAL 最低的位点，直到重叠消除；`missing_value` 取 `0` 或 `DP`）。被选中的位点会删除；给出 `mark_tag` 时则保留并标记 `INFO/TAG`（`-M`）。`reverse` 使用相反逻辑。`output_type` 为 `t` 或 `tz` 时写出染色体与位置的纯文本（压缩）列表，而不是 VCF。
- **API：** `dnakit.variants.bcftools_remove_overlaps(vcf_path, output_path, *, mark='overlap', mark_tag=None, missing_value='0', reverse=False, output_type=None, regions=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type`、`records`（仅纯文本 VCF 输出）和 `messages`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_remove_overlaps

work = Path(tempfile.mkdtemp())
dup = Path("variants.vcf").read_text().splitlines()
row = dup[-1].split("\t")
row[2] = "copy"
(work / "dup.vcf").write_text("\n".join(dup + ["\t".join(row)]) + "\n")
removed = bcftools_remove_overlaps(work / "dup.vcf", work / "removed.vcf", mark="dup")
marked = bcftools_remove_overlaps(work / "dup.vcf", work / "marked.vcf", mark="dup", mark_tag="DUP")
listed = bcftools_remove_overlaps(work / "dup.vcf", work / "sites.txt", mark="dup", output_type="t")
print(removed["records"], marked["records"], (work / "sites.txt").read_text().split())
```

- **示例输出：**

```text
9 11 ['chr1', '60', 'chr1', '120', 'chr1', '180', 'chr1', '240', 'chr1', '300', 'chr1', '360', 'chr1', '420', 'chr1', '480', 'chr1', '520']
```

## 3) `dnakit.variants.bcftools_scatter_vcf()`

- **作用：** 按每块 N 个位点或按区域窗口把 VCF/BCF 拆成多个文件（`bcftools +scatter`）。
- **处理方法：** `sites_per_chunk`（`-n`）、`regions_to_scatter`（`-s`，每个区域一个输出）、`scatter_file`（`-S`，区域加可选的第二列用于命名输出）三者必须恰好给一个。`extra_name` 把不在所列区域内的记录写到单独文件（`-x`），`prefix` 加在输出文件名前（`-p`）。文件写入新文件夹 `output_dir`。
- **API：** `dnakit.variants.bcftools_scatter_vcf(vcf_path, output_dir, *, sites_per_chunk=None, regions_to_scatter=None, scatter_file=None, extra_name=None, prefix=None, output_type='v', bcftools_path=None)`.
- **输入：** VCF/BCF（按区域拆分需要索引）和一个新的输出文件夹。
- **输出：** 通用字段、`output_type` 和 `files`（写出的文件名）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_scatter_vcf

result = bcftools_scatter_vcf("variants.vcf", Path(tempfile.mkdtemp()) / "chunks", sites_per_chunk=4, prefix="part-")
print(result["files"])
```

- **示例输出：**

```text
['part-0.vcf', 'part-1.vcf', 'part-2.vcf']
```

## 4) `dnakit.variants.bcftools_set_genotypes()`

- **作用：** 按规则设置基因型：缺失改为参考、主要或次要等位基因、定相或去相、自定义基因型（`bcftools +setGT` 与 `bcftools +missing2ref`）。
- **处理方法：** `target_genotypes`（`-t`）：`./.`、`./x`、`.`、`a`、`q`（`genotype_include` / `genotype_exclude` 为真的基因型，即插件的 `-i` / `-e`）、`b`（二项检验未通过的杂合子，如 `b:AD<1e-3`）或 `r:FLOAT`（随机比例）。`new_genotype`（`-n`）：`.`、`0`、`c:GT`、`i`、`m`、`M`、`X`、`p`、`u`，可写成 `0p`/`Mp`。`+missing2ref` 即 `target_genotypes='.'` 配合 `new_genotype` 为 `0`（其 `-p` 为 `0p`）或 `M`（`-m -p` 为 `Mp`），输出完全相同。
- **API：** `dnakit.variants.bcftools_set_genotypes(vcf_path, output_path, *, target_genotypes='.', new_genotype, genotype_include=None, genotype_exclude=None, random_seed=0, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_set_genotypes

work = Path(tempfile.mkdtemp())
bcftools_set_genotypes("variants.vcf", work / "ref.vcf", target_genotypes=".", new_genotype="0p")
bcftools_set_genotypes("trio.vcf", work / "low_depth.vcf", target_genotypes="q", new_genotype=".", genotype_include="FMT/DP<31")
print([line.split("\t")[11] for line in (work / "ref.vcf").read_text().splitlines()[-8:-6]])
print(sum("./." in line for line in (work / "low_depth.vcf").read_text().splitlines() if not line.startswith("#")))
```

- **示例输出：**

```text
['0/0', '0/0']
14
```

## 5) `dnakit.variants.bcftools_split_vcf_by_sample()`

- **作用：** 按样本或样本组拆分 VCF/BCF（`bcftools +split`）。
- **处理方法：** 默认每个样本写出一个以样本名命名的单样本文件。`samples_file`（`-S`）列出每个输出文件的样本（可选第二列重命名，可选第三列指定文件名）；`groups_file`（`-G`）把每个样本分配到一个或多个输出。`keep_tags` 限制保留的标签，如 `['INFO', 'FMT/GT']`（`-k`）。文件写入新文件夹 `output_dir`。
- **API：** `dnakit.variants.bcftools_split_vcf_by_sample(vcf_path, output_dir, *, samples_file=None, groups_file=None, keep_tags=None, regions=None, include=None, exclude=None, output_type='v', bcftools_path=None)`.
- **输入：** 多样本 VCF/BCF 和一个新的输出文件夹。
- **输出：** 通用字段、`output_type` 和 `files`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_split_vcf_by_sample

work = Path(tempfile.mkdtemp())
(work / "groups.txt").write_text("S1,S2\t-\tpair\nS3\n")
print(bcftools_split_vcf_by_sample("variants.vcf", work / "each", keep_tags=["INFO", "FMT/GT"])["files"])
print(sorted(bcftools_split_vcf_by_sample("variants.vcf", work / "grouped", samples_file=work / "groups.txt")["files"]))
```

- **示例输出：**

```text
['S1.vcf', 'S2.vcf', 'S3.vcf', 'S4.vcf', 'S5.vcf', 'S6.vcf']
['S3.vcf', 'pair.vcf']
```

## 6) `dnakit.variants.bcftools_split_vep()`

- **作用：** 查询或提取结构化的后果注释：VEP 的 `CSQ`、bcftools 的 `BCSQ` 和 SnpEff 的 `ANN`（`bcftools +split-vep`）。
- **处理方法：** `annotation` 指定 INFO 标签（可自动检测）。`list_fields` 列出注释字段（`-l`）。`columns` 按名称或 0 起始的序号把字段提取为新的 INFO 标签（`Consequence,IMPACT`、`1-3`、`-` 表示全部；`:Integer`/`:Float` 指定类型）。`select` 为 `TR:CSQ[:PRN]`，如 `worst` 或 `primary:missense+`。给出 `text_format`（`bcftools query` 格式，`-f`）时输出文本表而不是 VCF。`annotation_prefix`、`all_fields_delimiter`、`severity_file`、`gene_list`、`duplicate`、`drop_sites`、`allow_undefined_tags`、`print_header` 对应 `-p -A -S -g -d -x/-X -u -H`。
- **API：** `dnakit.variants.bcftools_split_vep(vcf_path, output_path, *, annotation=None, columns=None, select=None, text_format=None, list_fields=False, annotation_prefix=None, all_fields_delimiter=None, severity_file=None, gene_list=None, duplicate=False, drop_sites=None, allow_undefined_tags=False, print_header=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 带后果标签的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段；VCF 输出时有 `records`，文本输出时有 `report`，`list_fields` 时还有 `fields`（序号和名称）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_split_vep

work = Path(tempfile.mkdtemp())
fields = bcftools_split_vep("csq.vcf", work / "fields.txt", list_fields=True)
print([field["name"] for field in fields["fields"]])
table = bcftools_split_vep("csq.vcf", work / "table.txt", select="worst", text_format="%CHROM\t%POS\t%Consequence\t%IMPACT\n")
print(table["report"]["rows"])
vcf = bcftools_split_vep("csq.vcf", work / "vep.vcf", columns="Consequence,IMPACT", select="worst")
print(vcf["records"], (work / "vep.vcf").read_text().splitlines()[-3].split("\t")[7])
```

- **示例输出：**

```text
['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Feature']
[['chr1', 10, 'missense_variant', 'MODERATE'], ['chr1', 20, 'stop_gained', 'HIGH']]
3 CSQ=G|missense_variant|MODERATE|GENE1|T1,G|synonymous_variant|LOW|GENE1|T2;Consequence=missense_variant;IMPACT=MODERATE
```

## 7) `dnakit.variants.bcftools_convert_format_tags()`

- **作用：** 在等价的 FORMAT 标签之间转换：GL、PL、GP 与 GT，QR/QA 到 QS，以及局部化等位基因（`bcftools +tag2tag`）。
- **处理方法：** `conversion` 取 `GL-to-PL`、`GL-to-GP`、`GL-to-GT`、`PL-to-GL`、`PL-to-GP`、`PL-to-GT`、`GP-to-GL`、`GP-to-PL`、`GP-to-GT`、`QR-QA-to-QS`、`LXX-to-XX`、`LPL-to-PL` 或 `LAD-to-AD`。`replace_source` 删除源标签（`-r`），`gp_threshold` 是 GP 转 GT 硬判的最大误差（`-t`），`skip_nalt` 使局部化转换中 ALT 等位基因更多的位点保持不变（`-s`），`defaults` 填充缺失的局部化字段（`-d`）。相反方向的 `--XX-to-LXX` 在 bcftools 1.24 帮助中列出但运行时只输出 "todo"，所以不提供。
- **API：** `dnakit.variants.bcftools_convert_format_tags(vcf_path, output_path, *, conversion, replace_source=False, gp_threshold=0.1, skip_nalt=0, defaults=None, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 含源标签的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_convert_format_tags

output = Path(tempfile.mkdtemp()) / "gl.vcf"
result = bcftools_convert_format_tags("trio.vcf", output, conversion="PL-to-GL")
row = output.read_text().splitlines()[-1].split("\t")
print(result["records"], row[8])
print(row[9])
```

- **示例输出：**

```text
14 GT:PL:AD:DP:GQ:QS:SP:GL
1/1:255,60,0:0,30:30:60:0,1140:0:-25.5,-6,-0
```

## 8) `dnakit.variants.bcftools_variant_distance()`

- **作用：** 为位点标注到最近变异的距离（`bcftools +variant-distance`）。
- **处理方法：** `direction` 取 `nearest`（默认）、`fwd`、`rev` 或 `both`（两个值：到前一个和后一个变异的距离）；`tag_name` 是写入的 INFO 标签（`-n`）。
- **API：** `dnakit.variants.bcftools_variant_distance(vcf_path, output_path, *, direction='nearest', tag_name='DIST', regions=None, output_type=None, bcftools_path=None)`.
- **输入：** VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variant_distance

output = Path(tempfile.mkdtemp()) / "distance.vcf"
bcftools_variant_distance("variants.vcf", output, direction="both")
for line in output.read_text().splitlines()[-9:-6]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **示例输出：**

```text
120 DP=55;AF=0.333;DIST=60,60
180 DP=60;AF=0.5;DIST=60,60
240 DP=65;AF=0.417;DIST=60,60
```

## 9) `dnakit.variants.bcftools_fixref()`

- **作用：** 对照参考序列报告或修复 REF/ALT 的链方向与互换（`bcftools +fixref`；对应工作清单中 `-m flip|flip-all|id|ref-alt|swap|top` 与 `-m stats` 两类）。
- **处理方法：** `mode="stats"`（默认）保存链方向与不匹配统计（`SC`、`ST`、`NS` 行，解析为 `sections`）。转换模式写出 VCF/BCF：`flip`、`flip-all`（也处理 A/T 和 C/G 位点）、`id`（用 ID 列对照 `dbsnp_vcf`）、`ref-alt`、`swap`、`top`（Illumina TOP 链转正向链）。`discard_unresolved` 丢弃无法修复的位点（`-d`）；`tag_name` 重命名 `INFO/FIXREF`（`-t`）。不了解数据的链约定时请勿转换。
- **API：** `dnakit.variants.bcftools_fixref(vcf_path, output_path, *, reference_fasta, mode='stats', dbsnp_vcf=None, discard_unresolved=False, tag_name=None, regions=None, output_type=None, bcftools_path=None)`.
- **输入：** VCF/BCF、已建索引的参考 FASTA 和一个新的输出路径。
- **输出：** 通用字段和 `mode`；`stats` 模式有 `report` 与 `sections`，转换模式有 `records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fixref

import subprocess

subprocess.run(["samtools", "faidx", "genome.fa"], check=True)
work = Path(tempfile.mkdtemp())
stats = bcftools_fixref("variants.vcf", work / "stats.txt", reference_fasta="genome.fa")
print({str(row[0]).strip(): row[1] for row in stats["sections"]["NS"]})
fixed = bcftools_fixref("variants.vcf", work / "flip.vcf", reference_fasta="genome.fa", mode="flip")
print(fixed["mode"], fixed["records"])
```

- **示例输出：**

```text
{'total': 10, 'ref match': 8, 'ref mismatch': 0, 'errors': 0, 'skipped': 2, 'non-ACGT': 0, 'non-SNP': 2, 'non-biallelic': 0}
flip 10
```

## 10) `dnakit.variants.bcftools_color_chromosomes()`

- **作用：** 根据已定相的基因型，找出三联体中或两个无关样本之间共享的染色体片段（`bcftools +color-chrs`）。
- **处理方法：** `trio` 为 `[母亲, 父亲, 孩子]`（`-t`），`unrelated` 为两个样本名（`-u`）；两者必须恰好给一个。片段写入新文件夹 `output_dir` 中的 `<prefix>.dat`，这是 bcftools 的 `color-chrs.pl` 作图所用的格式。
- **API：** `dnakit.variants.bcftools_color_chromosomes(vcf_path, output_dir, *, trio=None, unrelated=None, prefix='color-chrs', regions=None, bcftools_path=None)`.
- **输入：** 带定相基因型的 VCF/BCF 和一个新的输出文件夹。
- **输出：** 通用字段、`files`、`segments`（`SG` 行：染色体、起点、终点、单倍型标签）和 `switches`（`SW` 行）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_color_chromosomes

import subprocess

work = Path(tempfile.mkdtemp())
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "phased.vcf.gz"), "trio_phased.vcf"], check=True)
subprocess.run(["bcftools", "index", str(work / "phased.vcf.gz")], check=True)
result = bcftools_color_chromosomes(work / "phased.vcf.gz", work / "colors", trio=["M", "F", "C"])
print(result["files"], result["segments"], result["switches"])
```

- **示例输出：**

```text
['color-chrs.dat'] [['chr1', 1, 399, 'M:2', 'F:2'], ['chr1', 400, 439, 'M:1', 'F:2'], ['chr1', 440, 479, 'M:1', 'F:1']] [['M', 'chr1', 4, 1, 0.333333], ['F', 'chr1', 6, 1, 0.2]]
```

## 11) `dnakit.variants.bcftools_parental_origin()`

- **作用：** 判断三联体中拷贝数变异的亲本来源（`bcftools +parental-origin`）。
- **处理方法：** `trio` 为 `[先证者, 父亲, 母亲]`（`-p`），`cnv_type` 取 `del` 或 `dup`（`-t`），`region` 为 `chr:start-end` 形式的 CNV（`-r`）。ALT 比例偏斜的亲本杂合位点在低于 `min_binomial_probability` 时被排除（`-b`）；`greedy` 也使用有歧义的位点（`-g`）；`list_informative_sites` 列出这些位点（`-d`）。VCF 需要索引、`FORMAT/PL` 和 `FORMAT/AD`。
- **API：** `dnakit.variants.bcftools_parental_origin(vcf_path, output_path, *, trio, cnv_type, region, min_binomial_probability=0.01, greedy=False, list_informative_sites=False, include=None, exclude=None, bcftools_path=None)`.
- **输入：** 已建索引的 VCF/BCF 和一个新的输出路径。
- **输出：** 通用字段、`report`、`predicted_origin`（`paternal`、`maternal` 或 `uncertain`）、`quality` 和 `markers`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_parental_origin

import subprocess

work = Path(tempfile.mkdtemp())
lines = [
    "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\tF\tM",
]
for index in range(30):  # proband carries an extra ALT-bearing copy that the father's heterozygous sites explain
    lines.append(f"chr1\t{1000 + 100 * index}\t.\tA\tG\t50\tPASS\t.\tGT:AD:PL\t0/1:10,20:60,0,60\t0/1:15,15:60,0,60\t0/0:30,0:0,60,255")
(work / "cnv.vcf").write_text("\n".join(lines) + "\n")
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "cnv.vcf.gz"), str(work / "cnv.vcf")], check=True)
subprocess.run(["bcftools", "index", str(work / "cnv.vcf.gz")], check=True)
result = bcftools_parental_origin(work / "cnv.vcf.gz", work / "origin.txt", trio=["P", "F", "M"], cnv_type="dup", region="chr1:1-5000")
print(result["predicted_origin"], result["markers"])
```

- **示例输出：**

```text
paternal 30
```

## 12) `dnakit.variants.bcftools_trio_stats()`

- **作用：** 统计三联体中孩子的传递率和孟德尔统计量（`bcftools +trio-stats`）。
- **处理方法：** 三联体来自 `ped_path`（`-p`）或 `trio=[先证者, 父亲, 母亲]`（`-P`）。`include` / `exclude` 可以使用 `GQ>{10,20}` 这种花括号范围同时扫描多个阈值，每个过滤器和每个三联体各一行 `FLT`。`max_alt_trios` 把传递率限制在至多有这么多个含替代等位基因三联体的位点（`-a`）；`debug` 取自 `mendel-errors` 和 `transmitted`（`-d`）。
- **API：** `dnakit.variants.bcftools_trio_stats(vcf_path, output_path, *, ped_path=None, trio=None, max_alt_trios=0, debug=None, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **输入：** VCF/BCF、PED 文件或三联体名称，以及一个新的输出路径。
- **输出：** 通用字段、`report`、`trio_columns` 和 `trios`：每个过滤器和三联体一个字典，含有效基因型数、非参考基因型数、DNM 或孟德尔错误数、单例计数、ts、tv、ts/tv、纯合与重复 DNM。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_stats

result = bcftools_trio_stats("trio.vcf", Path(tempfile.mkdtemp()) / "stats.txt", ped_path="trio.ped", include="GQ>{10,20}")
for row in result["trios"]:
    print(row["filter"], row["child"], row["valid_gts"], row["dnm_or_mendel_errors"], row["transmitted_singletons"])
```

- **示例输出：**

```text
FLT0 C 14 4 2
FLT1 C 14 4 2
```

## 13) `dnakit.variants.bcftools_trio_switch_rate()`

- **作用：** 根据已定相基因型计算三联体孩子的相位切换率（`bcftools +trio-switch-rate`）。
- **处理方法：** `ped_path` 列出三联体；可选的第 7 列按群体对结果分组。
- **API：** `dnakit.variants.bcftools_trio_switch_rate(vcf_path, output_path, *, ped_path, regions=None, bcftools_path=None)`.
- **输入：** 带定相基因型的 VCF/BCF、PED 文件和一个新的输出路径。
- **输出：** 通用字段、`report` 和 `trios`：父亲、母亲、孩子、检验的基因型数、孟德尔错误数、切换数和切换百分比。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_switch_rate

result = bcftools_trio_switch_rate("trio_phased.vcf", Path(tempfile.mkdtemp()) / "switch.txt", ped_path="trio.ped")
print(result["trios"])
```

- **示例输出：**

```text
[{'father': 'F', 'mother': 'M', 'child': 'C', 'tested': 4, 'mendelian_errors': 2, 'switches': 0, 'switch_percent': 0.0}]
```

## 14) `dnakit.variants.bcftools_mendelian_consistency()`

- **作用：** 按孟德尔遗传一致性对三联体的基因型计数、注释或过滤（`bcftools +mendelian2`；对应工作清单中 `-m c` 与 `-m a|d|e|E|g|m|M|S` 两类）。
- **处理方法：** 三联体来自 `ped_path`（`-P`）或 `trio=[孩子, 父亲, 母亲]`（`-p`；`1X:`/`2X:` 前缀选择 X 染色体模式）。`mode="c"`（默认）保存文本摘要，含每个三联体的良好、错误和缺失基因型数。其他情况下 `mode` 是由 VCF 输出字母组成的字符串：`a` 添加 `INFO/MERR` 等计数，`d` 删除不一致三联体的基因型，`e` 只保留有错误的位点，`E` 丢弃它们，`g` 保留有良好三联体的位点，`m` 保留有缺失基因型的位点，`M` 丢弃它们，`S` 丢弃被跳过的位点。`rules` 是预定义的遗传规则集，`rules_file` 是自定义规则文件。
- **API：** `dnakit.variants.bcftools_mendelian_consistency(vcf_path, output_path, *, ped_path=None, trio=None, mode='c', rules='GRCh37', rules_file=None, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 带 `FORMAT/GT` 的 VCF/BCF、三联体和一个新的输出路径。
- **输出：** 通用字段和 `mode`；`c` 模式有 `report`，VCF 模式有 `records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_mendelian_consistency

work = Path(tempfile.mkdtemp())
counts = bcftools_mendelian_consistency("trio.vcf", work / "counts.txt", ped_path="trio.ped")
print({row[0]: row[1] for row in counts["report"]["rows"] if row[0] in ("sites_merr", "sites_good")})
for mode in ("E", "e"):
    print(mode, bcftools_mendelian_consistency("trio.vcf", work / f"{mode}.vcf", ped_path="trio.ped", mode=mode)["records"])
```

- **示例输出：**

```text
{'sites_merr': 4, 'sites_good': 10}
E 10
e 4
```

## 15) `dnakit.variants.bcftools_trio_denovo()`

- **作用：** 在三联体中筛查可能的新生突变并写入 `FORMAT/DNM`（`bcftools +trio-dnm3`；对应工作清单中 `--use-DMM`、`--use-ALM`、`--use-DNG`、`--use-NAIVE` 四个模型）。
- **处理方法：** `model` 取 `dmm`（考虑位点噪声的 Dirichlet 多项模型，bcftools 1.24 的默认）、`alm`（等位基因似然模型，1.23 及以前的默认）、`dng`（原始 DeNovoGear）或 `naive`（只用基因型，`DNM:flag`）。三联体来自 `ped_path`（`-P`）或 `trio=[孩子, 父亲, 母亲]`（`-p`）。`dnm_tag` 是输出标签和分数类型（默认 `DNM:log`；可选 `log`、`flag`、`phred`、`prob`），`min_score` 省略低于它的分数（`-m`），`chrx_regions` 是 X 染色体区域或预设（`-X`），`mutation_rate` 是先验（`--mrate`），`strictly_novel` 即 `-n`。模型选项：`noise_prior`、`phi`（dmm、alm），`min_vaf`、`strand_bias`（dmm、alm；`strand_bias` 也适用于 dng），`max_qm`（dmm；负数表示不使用 `FORMAT/QM`），`allele_dropout`（alm），`use_parental_ad`、`use_child_ad`（dmm），`force_ad`。不属于所选模型的选项会被拒绝。
- **API：** `dnakit.variants.bcftools_trio_denovo(vcf_path, output_path, *, ped_path=None, trio=None, model='dmm', min_score=None, dnm_tag=None, chrx_regions='GRCh37', mutation_rate=1e-08, strictly_novel=False, noise_prior=None, phi=None, min_vaf=None, max_qm=None, strand_bias=None, allele_dropout=None, use_parental_ad=False, use_child_ad=False, force_ad=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **输入：** 含 `GT`、`PL`、`AD`（及 `QS`）的 VCF/BCF、三联体和一个新的输出路径。
- **输出：** 通用字段、`model`、`output_type` 和 `records`（仅纯文本 VCF 输出）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_denovo

work = Path(tempfile.mkdtemp())
naive = bcftools_trio_denovo("trio.vcf", work / "naive.vcf", trio=["C", "F", "M"], model="naive")
flagged = [
    line.split("\t")[1]
    for line in (work / "naive.vcf").read_text().splitlines()
    if not line.startswith("#") and "DNM" in line.split("\t")[8].split(":")
]
alm = bcftools_trio_denovo("trio.vcf", work / "alm.vcf", trio=["C", "F", "M"], model="alm")
scores = {}
for line in (work / "alm.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        scores[int(fields[1])] = float(fields[9].split(":")[fields[8].split(":").index("DNM")])
print(naive["model"], flagged)
print(alm["records"], max(scores, key=scores.get), round(max(scores.values()), 2))
```

- **示例输出：**

```text
naive ['200', '280', '480', '560']
14 280 -3.93
```
