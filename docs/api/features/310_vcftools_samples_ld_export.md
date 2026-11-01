# VCFtools：个体统计、连锁不平衡、PCA、格式导出与比较

这一组 `dnakit.variants` 中的函数延续上一篇文档（`309_vcftools_filters_site_reports`）里带类型的 `vcftools` 封装；通用行为、结果字典、`filters` 映射和 `output_prefix` 规则见该文档。它们计算每个个体的杂合度与近交系数、变异负荷、长纯合片段和亲缘关系；六种形式的成对连锁不平衡；每个 BED 区域内的唯一单倍型数；主成分与 SNP 载荷；把基因型导出为 0/1/2、IMPUTE、LDhat、LDhelmet、BEAGLE、PLINK PED/MAP 和 TPED/TFAM 文件；以及比较两个 VCF 文件（个体、基因型不一致、相位切换错误）。

VCFtools 使用 LGPL-3.0 许可；DNAKit 把 `vcftools` 作为独立进程运行，不复制其任何代码。请自行安装（`conda install -c bioconda vcftools`、`apt install vcftools` 或从源码编译），放到 PATH 上或传入 `vcftools_path`；`vcftools_pca()` 需要用 `./configure --enable-pca`（LAPACK）编译的版本，否则抛出带该提示的 `BackendUnavailableError`。

VCFtools 0.1.16 的三个上游行为被原样保留。`--hapcount` 把 BED 文件的第一行当作表头跳过，并且当 VCF 的最后一条记录落在最后一条染色体的最后一个区域内时会读取已释放的内存，可能使最后一行变成 0；`vcftools_haplotype_count()` 对这两种情况都在 `warnings` 中提示。长纯合片段、LDhat、LDhelmet 和 BEAGLE 输出需要指定一条染色体。

示例从 `dnakit.variants` 导入，并写入新的临时目录。

## 示例数据

`toy.vcf` 与上一篇文档中的 6 样本 VCF 相同：`chr1`（14 个）和 `chr2`（6 个）共 20 个已定相位点，3 个 indel、2 个 `LowQual` 位点、1 个缺失基因型，FORMAT 字段为 GT、DP、GQ 和 PL。`win.bed` 有一行表头和三个区域；`pos.txt` 列出两个位点。

**示例文件：** `toy.vcf`

```text
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##contig=<ID=chr1,length=20000>
##contig=<ID=chr2,length=20000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	100	v1	G	T	38	PASS	DP=91;AF=0.1	GT:DP:GQ:PL	1|0:9:45:30,0,30	0|1:31:60:30,0,30	0|1:9:30:30,0,30	0|1:20:99:30,0,30	1|0:8:12:30,0,30	0|0:14:99:0,30,300
chr1	800	v2	G	A	55	PASS	DP=120;AF=0.25	GT:DP:GQ:PL	1|0:25:99:30,0,30	1|1:21:60:300,30,0	0|0:17:45:0,30,300	0|0:21:60:0,30,300	0|0:11:45:0,30,300	0|0:25:99:0,30,300
chr1	1500	v3	A	G	67	PASS	DP=99;AF=0.5	GT:DP:GQ:PL	1|1:10:45:300,30,0	0|1:11:99:30,0,30	0|1:23:60:30,0,30	0|0:36:45:0,30,300	0|1:10:60:30,0,30	1|1:9:45:300,30,0
chr1	2200	v4	A	AAC	38	PASS	DP=166;AF=0.5	GT:DP:GQ:PL	1|0:33:30:30,0,30	1|1:23:12:300,30,0	1|0:32:30:30,0,30	1|0:28:30:30,0,30	1|1:30:30:300,30,0	1|1:20:30:300,30,0
chr1	2900	v5	A	G	94	PASS	DP=172;AF=0.25	GT:DP:GQ:PL	0|1:31:60:30,0,30	0|0:31:12:0,30,300	0|0:31:30:0,30,300	0|0:31:12:0,30,300	0|0:12:30:0,30,300	1|0:36:60:30,0,30
chr1	3600	v6	C	T	22	LowQual	DP=180;AF=0.1	GT:DP:GQ:PL	0|0:37:30:0,30,300	0|0:35:12:0,30,300	0|1:36:45:30,0,30	0|0:36:45:0,30,300	0|0:25:60:0,30,300	0|0:11:30:0,30,300
chr1	4300	v7	A	G	81	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	0|0:27:30:0,30,300	1|1:20:99:300,30,0	0|0:18:60:0,30,300	0|0:21:45:0,30,300	0|0:31:12:0,30,300	0|0:20:12:0,30,300
chr1	5000	v8	G	A	81	PASS	DP=124;AF=0.25	GT:DP:GQ:PL	0|0:6:30:0,30,300	0|0:36:60:0,30,300	0|0:28:30:0,30,300	0|1:11:60:30,0,30	1|1:13:45:300,30,0	1|0:30:12:30,0,30
chr1	5700	v9	TGT	T	55	PASS	DP=93;AF=0.5	GT:DP:GQ:PL	1|0:15:30:30,0,30	0|1:14:60:30,0,30	1|1:7:30:300,30,0	0|0:6:30:0,30,300	1|0:12:12:30,0,30	0|0:39:45:0,30,300
chr1	6400	v10	C	T	38	PASS	DP=173;AF=0.25	GT:DP:GQ:PL	0|0:40:30:0,30,300	0|0:15:99:0,30,300	1|0:39:12:30,0,30	0|0:38:30:0,30,300	./.:.:.:.	0|0:34:30:0,30,300
chr1	7100	v11	C	T	22	LowQual	DP=177;AF=0.1	GT:DP:GQ:PL	0|0:34:60:0,30,300	0|0:26:99:0,30,300	0|0:38:99:0,30,300	0|0:38:60:0,30,300	0|0:18:99:0,30,300	0|0:23:30:0,30,300
chr1	7800	v12	G	GAC	94	PASS	DP=135;AF=0.25	GT:DP:GQ:PL	0|1:29:12:30,0,30	1|0:15:60:30,0,30	1|1:22:60:300,30,0	1|0:14:30:30,0,30	0|0:35:30:0,30,300	1|0:20:30:30,0,30
chr1	8500	v13	T	G	38	PASS	DP=93;AF=0.4	GT:DP:GQ:PL	0|1:12:45:30,0,30	1|1:11:30:300,30,0	1|0:22:60:30,0,30	0|1:23:45:30,0,30	0|1:8:60:30,0,30	0|1:17:30:30,0,30
chr1	9200	v14	T	A	55	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	1|0:32:12:30,0,30	0|1:23:30:30,0,30	0|0:14:45:0,30,300	0|1:8:12:30,0,30	0|1:39:30:30,0,30	0|0:21:30:0,30,300
chr2	300	v15	G	C	94	PASS	DP=204;AF=0.4	GT:DP:GQ:PL	0|0:33:30:0,30,300	1|0:37:30:30,0,30	0|1:40:45:30,0,30	1|0:31:30:30,0,30	0|0:38:30:0,30,300	1|0:25:60:30,0,30
chr2	1200	v16	G	C	38	PASS	DP=141;AF=0.1	GT:DP:GQ:PL	1|0:23:99:30,0,30	0|0:34:45:0,30,300	1|0:6:30:30,0,30	0|0:22:12:0,30,300	0|0:29:45:0,30,300	0|0:27:30:0,30,300
chr2	2100	v17	G	A	38	PASS	DP=158;AF=0.4	GT:DP:GQ:PL	1|0:11:30:30,0,30	0|1:39:45:30,0,30	0|1:15:99:30,0,30	1|1:30:30:300,30,0	1|1:26:12:300,30,0	1|1:37:99:300,30,0
chr2	3000	v18	T	A	94	PASS	DP=133;AF=0.1	GT:DP:GQ:PL	0|0:7:60:0,30,300	0|0:40:12:0,30,300	0|0:21:99:0,30,300	1|0:37:99:30,0,30	0|0:22:12:0,30,300	0|0:6:99:0,30,300
chr2	3900	v19	A	T	94	PASS	DP=117;AF=0.4	GT:DP:GQ:PL	0|0:18:99:0,30,300	1|0:10:99:30,0,30	1|0:15:30:30,0,30	0|0:27:12:0,30,300	1|0:22:60:30,0,30	1|1:25:12:300,30,0
chr2	4800	v20	T	C	67	PASS	DP=154;AF=0.1	GT:DP:GQ:PL	0|0:10:30:0,30,300	0|0:38:12:0,30,300	0|0:34:99:0,30,300	0|0:23:12:0,30,300	0|0:30:30:0,30,300	0|1:19:99:30,0,30
```

**示例文件：** `win.bed`

```text
#CHROM	START	END
chr1	0	5000
chr1	5000	12000
chr2	0	4000
```

**示例文件：** `pos.txt`

```text
chr1	100
chr1	800
```

## 1) `dnakit.variants.vcftools_heterozygosity()`

- **作用：** 计算每个个体的纯合度和近交系数 F。
- **处理方法：** 运行 `vcftools --het`；表为 `<output_prefix>.het`。F 为矩估计（观察到的纯合位点数减去期望值，再除以位点数减期望值），与 PLINK `--het` 相同。
- **API：** `dnakit.variants.vcftools_heterozygosity(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `INDV, O(HOM), E(HOM), N_SITES, F`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_heterozygosity

result = vcftools_heterozygosity("toy.vcf", str(Path(tempfile.mkdtemp()) / "het"))
print(result["table"]["columns"])
for row in result["table"]["rows"][:3]:
    print(row)
```

- **示例输出：**

```text
['INDV', 'O(HOM)', 'E(HOM)', 'N_SITES', 'F']
['S1', 9, 11.5, 19, -0.33549]
['S2', 11, 11.5, 19, -0.06839]
['S3', 9, 11.5, 19, -0.33549]
```

## 2) `dnakit.variants.vcftools_individual_burden()`

- **作用：** 按基因型类别或等位基因频率类别统计每个个体的变异负荷。
- **处理方法：** 运行 `vcftools --indv-burden`（`burden="genotype_classes"`，`.iburden`：纯合参考、杂合、纯合变异和缺失基因型）、`--indv-freq-burden`（`"frequency_classes"`，`.ifreqburden`：按位点等位基因计数对非参考等位基因分类，纯合基因型的两个等位基因都计数）或 `--indv-freq-burden2`（`"frequency_classes_hom_once"`，文件名相同，纯合基因型只计一次）。
- **API：** `dnakit.variants.vcftools_individual_burden(vcf_path: 'str', output_prefix: 'str', *, burden: 'str' = 'genotype_classes', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`burden`、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典（`INDV, N_HOM_REF, N_HET, N_HOM_ALT, N_MISS`，或 `INDV` 加每个等位基因计数一列）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_individual_burden

out = Path(tempfile.mkdtemp())
by_genotype = vcftools_individual_burden("toy.vcf", str(out / "g"))
print(by_genotype["table"]["columns"])
print(by_genotype["table"]["rows"][:2])
both = vcftools_individual_burden("toy.vcf", str(out / "a"), burden="frequency_classes")
once = vcftools_individual_burden("toy.vcf", str(out / "b"), burden="frequency_classes_hom_once")
print(sum(sum(r[1:]) for r in both["table"]["rows"]), sum(sum(r[1:]) for r in once["table"]["rows"]))
```

- **示例输出：**

```text
['INDV', 'N_HOM_REF', 'N_HET', 'N_HOM_ALT', 'N_MISS']
[['S1', 9, 10, 1, 0], ['S2', 8, 8, 4, 0]]
78 63
```

## 3) `dnakit.variants.vcftools_long_runs_of_homozygosity()`

- **作用：** 在一条染色体上寻找长纯合片段。
- **处理方法：** 运行 `vcftools --LROH --chr <chromosome>`；表为 `<output_prefix>.LROH`。VCFtools 要求指定染色体。
- **API：** `dnakit.variants.vcftools_long_runs_of_homozygosity(vcf_path: 'str', output_prefix: 'str', chromosome: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`chromosome`（第三个位置参数；`filters` 中不能再有 `chromosomes`）、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典（VCFtools 报告片段的起点、终点和长度；没有片段时表为空）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_long_runs_of_homozygosity

result = vcftools_long_runs_of_homozygosity("toy.vcf", str(Path(tempfile.mkdtemp()) / "roh"), "chr1")
print(result["status"], result["sites_kept"], Path(result["outputs"][0]).suffix)
print(result["table"]["columns"], result["table"]["row_count"])
```

- **示例输出：**

```text
completed 14 .LROH
['CHROM', 'AUTO_START', 'AUTO_END', 'MIN_START', 'MAX_END', 'N_VARIANTS_BETWEEN_MAX_BOUNDARIES', 'N_MISMATCHES', 'INDV'] 0
```

## 4) `dnakit.variants.vcftools_relatedness()`

- **作用：** 估计个体两两之间的亲缘关系。
- **处理方法：** 运行 `vcftools --relatedness`（`method="yang"`，Yang 等 2010，`.relatedness`，统计量 `RELATEDNESS_AJK`）或 `--relatedness2`（`"manichaikul"`，Manichaikul 等 2010，`.relatedness2`，统计量 `RELATEDNESS_PHI`）。
- **API：** `dnakit.variants.vcftools_relatedness(vcf_path: 'str', output_prefix: 'str', *, method: 'str' = 'yang', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`method`、`filters`、`max_rows`。
- **输出：** 带 `table`（包括自身配对在内的所有个体对）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_relatedness

out = Path(tempfile.mkdtemp())
yang = vcftools_relatedness("toy.vcf", str(out / "y"), max_rows=2)
print(yang["table"]["columns"], yang["table"]["row_count"], yang["table"]["rows"])
phi = vcftools_relatedness("toy.vcf", str(out / "m"), method="manichaikul", max_rows=2)
print(phi["table"]["columns"])
print(phi["table"]["rows"])
```

- **示例输出：**

```text
['INDV1', 'INDV2', 'RELATEDNESS_AJK'] 21 [['S1', 'S1', 0.632232], ['S1', 'S2', -0.0434305]]
['INDV1', 'INDV2', 'N_AaAa', 'N_AAaa', 'N1_Aa', 'N2_Aa', 'RELATEDNESS_PHI']
[['S1', 'S1', 10, 0, 10, 10, 0.5], ['S1', 'S2', 5, 1, 10, 8, 0.166667]]
```

## 5) `dnakit.variants.vcftools_linkage_disequilibrium()`

- **作用：** 计算位点之间的成对连锁不平衡。
- **处理方法：** `statistic="geno_r2"` 运行 `--geno-r2`（由 0/1/2 基因型相关系数得到 r2，`.geno.ld`），`"geno_chisq"` 运行 `--geno-chisq`（基因型独立性卡方检验，`.geno.chisq`），`"hap_r2"` 运行 `--hap-r2`（由定相单倍型得到 r2，`.hap.ld`）。给出 `positions_file` 时只把列出的位点与所有其他位点比较（`--geno-r2-positions` `.list.geno.ld`，`--hap-r2-positions` `.list.hap.ld`）。`interchromosomal=True` 时比较不同染色体上的位点（`--interchrom-geno-r2` `.interchrom.geno.ld`，`--interchrom-hap-r2` `.interchrom.hap.ld`）。`ld_window_bp`、`ld_window_snps`、`ld_window_bp_min`、`ld_window_snps_min` 和 `min_r2` 对应 `--ld-window-bp`、`--ld-window`、`--ld-window-bp-min`、`--ld-window-min` 和 `--min-r2`。
- **API：** `dnakit.variants.vcftools_linkage_disequilibrium(vcf_path: 'str', output_prefix: 'str', *, statistic: 'str' = 'geno_r2', positions_file: 'str | None' = None, interchromosomal: 'bool' = False, ld_window_bp: 'int | None' = None, ld_window_snps: 'int | None' = None, ld_window_bp_min: 'int | None' = None, ld_window_snps_min: 'int | None' = None, min_r2: 'float | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`statistic`、`positions_file`（制表符分隔的 CHROM、POS）、`interchromosomal`、窗口选项、`min_r2`（0 到 1）、`filters`、`max_rows`。卡方检验没有列表或跨染色体形式；`positions_file` 和 `interchromosomal` 不能同时使用。
- **输出：** 带 `table` 的通用结果字典（`CHR, POS1, POS2, N_INDV, R^2`；单倍型为 `N_CHR, R^2, D, Dprime`；`geno_chisq` 为卡方统计列）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_linkage_disequilibrium

out = Path(tempfile.mkdtemp())
pairs = vcftools_linkage_disequilibrium("toy.vcf", str(out / "g"), max_rows=2)
print(pairs["table"]["columns"], pairs["table"]["row_count"])
print(pairs["table"]["rows"])
near = vcftools_linkage_disequilibrium("toy.vcf", str(out / "w"), ld_window_bp=1500, min_r2=0.1)
print(near["table"]["row_count"])
hap = vcftools_linkage_disequilibrium("toy.vcf", str(out / "h"), statistic="hap_r2", max_rows=1)
print(hap["table"]["columns"], hap["table"]["rows"])
listed = vcftools_linkage_disequilibrium("toy.vcf", str(out / "l"), positions_file="pos.txt")
print(Path(listed["outputs"][0]).name, listed["table"]["row_count"])
inter = vcftools_linkage_disequilibrium("toy.vcf", str(out / "i"), interchromosomal=True)
print(Path(inter["outputs"][0]).name, inter["table"]["row_count"])
```

- **示例输出：**

```text
['CHR', 'POS1', 'POS2', 'N_INDV', 'R^2'] 106
[['chr1', 100, 800, 6, 0.0857143], ['chr1', 100, 1500, 6, 0.294118]]
16
['CHR', 'POS1', 'POS2', 'N_CHR', 'R^2', 'D', 'Dprime'] [['chr1', 100, 800, 12, 0.0857143, 0.0625, 0.428571]]
l.list.geno.ld 13
i.interchrom.geno.ld 84
```

## 6) `dnakit.variants.vcftools_haplotype_count()`

- **作用：** 统计已定相 VCF 中每个 BED 区域内的唯一单倍型数。
- **处理方法：** 运行 `vcftools --hapcount <regions_bed>`；表为 `<output_prefix>.hapcount`。区域不能重叠。VCFtools 会跳过 BED 的第一行（预期是表头），并且当最后一条 VCF 记录落在最后一个区域内时可能把最后一行置零；`warnings` 会对此提示。
- **API：** `dnakit.variants.vcftools_haplotype_count(vcf_path: 'str', output_prefix: 'str', regions_bed: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`（已定相）、`output_prefix`、`regions_bed`、`filters`、`max_rows`。
- **输出：** 带 `table`（`CHROM, BIN_START, BIN_END, N_SNP, N_UNIQ_HAPS, N_GROUPS` 以及 multiplicity:frequency 对）和 `warnings` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_haplotype_count

result = vcftools_haplotype_count("toy.vcf", str(Path(tempfile.mkdtemp()) / "hap"), "win.bed")
for row in result["table"]["rows"]:
    print(row[:5])
print(len(result["warnings"]))
```

- **示例输出：**

```text
['chr1', 0, 5000, 8, 11]
['chr1', 5000, 12000, 5, 9]
['chr2', 0, 4000, 5, 9]
1
```

## 7) `dnakit.variants.vcftools_pca()`

- **作用：** 对基因型做主成分分析，或写出 SNP 载荷。
- **处理方法：** 运行 `vcftools --pca`（`output="components"`，`.pca`；`normalise=False` 时运行 `--pca-no-norm`）或 `vcftools --pca-snp-loadings <n>`（`output="snp_loadings"`，`.pca.loadings`，只能与标准化一起使用），方法为 Patterson、Price 和 Reich 2006。VCFtools 需要双等位二倍体位点、位点数多于个体数，并且编译时带有 PCA 支持（`--enable-pca`）。
- **API：** `dnakit.variants.vcftools_pca(vcf_path: 'str', output_prefix: 'str', *, output: 'str' = 'components', normalise: 'bool' = True, snp_loadings: 'int | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`output`、`normalise`、`snp_loadings`（成分数，载荷必填）、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典：`INDV, EIG_0, ...`（第一行 `EIGENVALUE` 是特征值，之后每个个体一行），或 `CHROM, POS, GAMMA_0, ...`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_pca

out = Path(tempfile.mkdtemp())
pcs = vcftools_pca("toy.vcf", str(out / "pca"))
print(pcs["status"], pcs["table"]["columns"][:3], pcs["table"]["row_count"])
print(pcs["table"]["rows"][0][:3])
loadings = vcftools_pca("toy.vcf", str(out / "load"), output="snp_loadings", snp_loadings=2)
print(loadings["table"]["columns"], loadings["table"]["row_count"], loadings["table"]["rows"][0])
```

- **示例输出：**

```text
completed ['INDV', 'EIG_0', 'EIG_1'] 7
['EIGENVALUE', 3.03327, 2.69104]
['CHROM', 'POS', 'GAMMA_0', 'GAMMA_1'] 19 ['chr1', 100, 0.493707, -0.370979]
```

## 8) `dnakit.variants.vcftools_export_genotypes()`

- **作用：** 把 VCF 的基因型导出为其他程序的输入格式。
- **处理方法：** `export_format` 选择 VCFtools 选项：`"012"`（`--012`；`.012`、`.012.indv`、`.012.pos`）、`"impute"`（`--IMPUTE`；`.impute.hap`、`.impute.hap.indv`、`.impute.legend`；完全定相的双等位 SNP）、`"ldhat"`（`--ldhat`；`.ldhat.sites`、`.ldhat.locs`）、`"ldhelmet"`（`--ldhelmet`；`.ldhelmet.snps`、`.ldhelmet.pos`）、`"ldhat_geno"`（`--ldhat-geno`；未定相的 LDhat）、`"beagle_gl"` 和 `"beagle_pl"`（`--BEAGLE-GL`、`--BEAGLE-PL`；`.BEAGLE.GL`、`.BEAGLE.PL`；需要 GL 或 PL FORMAT 字段）、`"plink"`（`--plink`；`.ped`、`.map`）和 `"plink_tped"`（`--plink-tped`；`.tped`、`.tfam`）。LDhat、LDhelmet 和 BEAGLE 需要指定一条 `chromosome`。
- **API：** `dnakit.variants.vcftools_export_genotypes(vcf_path: 'str', output_prefix: 'str', *, export_format: 'str' = '012', chromosome: 'str | None' = None, filters: 'Mapping[str, object] | None' = None, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`export_format`、`chromosome`、`filters`。
- **输出：** 通用结果字典；`outputs` 列出所有写出的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_export_genotypes

out = Path(tempfile.mkdtemp())
for fmt in ("012", "impute", "plink", "plink_tped"):
    result = vcftools_export_genotypes("toy.vcf", str(out / fmt), export_format=fmt)
    print(fmt, result["status"], sorted(Path(p).name.split(".", 1)[1] for p in result["outputs"]))
beagle = vcftools_export_genotypes("toy.vcf", str(out / "bgl"), export_format="beagle_pl", chromosome="chr1")
print(beagle["status"], beagle["sites_kept"])
print(Path(out / "plink.ped").read_text().splitlines()[0][:60])
print(Path(out / "012.012").read_text().splitlines()[0])
```

- **示例输出：**

```text
012 completed ['012', '012.indv', '012.pos']
impute completed ['impute.hap', 'impute.hap.indv', 'impute.legend']
plink completed ['map', 'ped']
plink_tped completed ['tfam', 'tped']
completed 14
S1	S1	0	0	0	0	T	G	A	G	G	G	AAC	A	A	G	C	C	A	A	G	G	T	TGT	C	C	C	
0	1	1	2	1	1	0	0	0	1	0	0	1	1	1	0	1	1	0	0	0
```

## 9) `dnakit.variants.vcftools_compare_vcfs()`

- **作用：** 比较两个 VCF 文件：共有个体、每个位点的基因型不一致，或相位切换错误。
- **处理方法：** 运行 `vcftools --vcf <第一个> --diff <第二个>`，并加 `--diff-indv`（`mode="individuals"`，`.diff.indv_in_files`，FILES 为 `1`、`2` 或 `B`（两者都有））、`--diff-site-discordance`（`"site_discordance"`，`.diff.sites`）或 `--diff-switch-error`（`"switch_error"`，`.diff.switch` 和 `.diff.indv.switch`）。第二个文件可以是 VCF、VCF.gz 或 BCF。`individual_map_file` 重命名第二个文件的个体（`--diff-indv-map`）。`filters` 作用于第一个文件。
- **API：** `dnakit.variants.vcftools_compare_vcfs(vcf_path: 'str', second_vcf_path: 'str', output_prefix: 'str', *, mode: 'str' = 'individuals', individual_map_file: 'str | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`second_vcf_path`、`output_prefix`、`mode`、`individual_map_file`、`filters`、`max_rows`。
- **输出：** 带第一个输出文件 `table` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_compare_vcfs

lines = Path("toy.vcf").read_text().splitlines()
changed = [x.replace("1|1:", "0|0:", 1) if x.startswith("chr1\t800") else x for x in lines]
Path("toy_changed.vcf").write_text("\n".join(changed) + "\n")
out = Path(tempfile.mkdtemp())
people = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "i"))
print(people["table"]["columns"], people["table"]["rows"][:2])
sites = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "s"), mode="site_discordance")
print(sites["table"]["columns"])
print([row for row in sites["table"]["rows"] if row[5] > 0])
switch = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "w"), mode="switch_error")
print(switch["status"], sorted(Path(p).name.split(".", 1)[1] for p in switch["outputs"]))
```

- **示例输出：**

```text
['INDV', 'FILES'] [['S1', 'B'], ['S2', 'B']]
['CHROM', 'POS', 'FILES', 'MATCHING_ALLELES', 'N_COMMON_CALLED', 'N_DISCORD', 'DISCORDANCE']
[['chr1', 800, 'B', 1, 6, 1, 0.166667]]
completed ['diff.indv.switch', 'diff.switch']
```
