# VCFtools：过滤与位点、深度、质量报告

这一组 `dnakit.variants` 中的函数是 [VCFtools](https://github.com/vcftools/vcftools) 的 `vcftools` 程序的带类型命令行封装（测试版本 0.1.16）。每个函数用带类型的参数列表运行一条固定的 `vcftools` 命令，不经过 shell。它们过滤位点、个体和基因型，列出通过过滤的位点，并写出深度、质量、缺失率、singleton、indel 长度、SNP 密度、Ts/Tv、位点多样性以及 INFO/FORMAT 字段表。只暴露它们所覆盖的 VCFtools 操作的选项。

VCFtools 使用 LGPL-3.0 许可。DNAKit 把 `vcftools` 作为独立进程运行，不复制其任何代码。该程序不是 Python 依赖：请自行安装（`conda install -c bioconda vcftools`、`apt install vcftools` 或从源码编译），并放到 PATH 上或传入 `vcftools_path`。缺少程序时函数抛出带安装提示的 `BackendUnavailableError`。

通用行为。`vcf_path` 可以是 `.vcf`、`.vcf.gz`（用 `--gzvcf` 读取）或 `.bcf`。`output_prefix` 是文件名前缀，与 VCFtools 的 `--out` 完全一致：结果文件为 `<output_prefix>.<后缀>`，后缀在各小节中给出。前缀所在文件夹不存在会抛出 `ConfigurationError`，已存在的输出文件不会被覆盖。报告类函数带 `filters` 参数，这是一个映射，键名与 `vcftools_filter_vcf()` 的关键字参数相同（例如 `{"min_maf": 0.05, "remove_filtered_all": True}`），因此每个报告都可以在过滤后的子集上计算。默认值就是 VCFtools 的默认值：只有给出的选项才会传给程序。程序报错退出时不抛异常，`status` 为 `failed`，`log_tail` 里是错误信息。每个结果都是可 JSON 序列化的字典，包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的文件）、`log_tail`、`output_prefix`，以及 VCFtools 报告时的 `sites_kept`、`sites_total`、`individuals_kept`、`individuals_total`。报告类函数还返回 `table`：`columns`、前 `max_rows` 行（`rows`，数值已转换，`nan` 为 `None`）、`row_count` 和 `truncated`。

示例从 `dnakit.variants` 导入，使用默认的详细输出，并写入新的临时目录。

## 示例数据

唯一需要声明的输入是 `toy.vcf`：6 个已定相的样本（S1 到 S6），`chr1` 和 `chr2` 上共 20 个位点，其中 3 个 indel、2 个 `LowQual` 位点、4 个 singleton 和 1 个 doubleton，1 个缺失基因型，FORMAT 字段为 GT、DP、GQ 和 PL。`chr1` 有 14 个位点，`chr2` 有 6 个。

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

## 1) `dnakit.variants.vcftools_filter_vcf()`

- **作用：** 用 VCFtools 过滤 VCF 或 BCF 的位点、个体和基因型，并把结果写成 VCF 或 BCF。
- **处理方法：** 根据给出的选项构造一条 `vcftools` 命令，并加上 `--recode`（`output_format="vcf"`）或 `--recode-bcf`（`"bcf"`）。所有给出的过滤条件与 VCFtools 一样同时生效。位点过滤：`chromosomes`、`exclude_chromosomes`（`--chr`、`--not-chr`）、`from_bp`、`to_bp`、`positions_file`、`exclude_positions_file`、`positions_overlap_file`、`exclude_positions_overlap_file`、`bed_file`、`exclude_bed_file`、`thin`（最小间距，bp）、`mask_file`、`invert_mask_file`、`mask_min`、`snp_ids`、`snp_ids_file`、`exclude_ids_file`、`only_indels`、`remove_indels`、`remove_filtered_all`、`keep_filtered`、`remove_filtered`、`keep_info`、`remove_info`、`min_maf`、`max_maf`、`min_mac`、`max_mac`、`min_alleles`、`max_alleles`、`min_mean_depth`、`max_mean_depth`、`hwe_pvalue`、`phased_only`、`min_quality`、`max_missing`。个体过滤：`keep_individuals`、`keep_individuals_files`、`remove_individuals`、`remove_individuals_files`、`max_individuals`（随机抽稀）。基因型过滤：`remove_filtered_genotypes`、`remove_filtered_genotype_flags`、`min_genotype_quality`、`min_genotype_depth`、`max_genotype_depth`。与 VCFtools 一致，输出会丢弃 INFO，除非 `recode_info` 为 `True`（保留全部 INFO）或 INFO 键列表。
- **API：** `dnakit.variants.vcftools_filter_vcf(vcf_path: 'str', output_prefix: 'str', *, output_format: 'str' = 'vcf', recode_info: 'bool | Sequence[str]' = False, chromosomes: 'Sequence[str] | None' = None, exclude_chromosomes: 'Sequence[str] | None' = None, from_bp: 'int | None' = None, to_bp: 'int | None' = None, positions_file: 'str | None' = None, exclude_positions_file: 'str | None' = None, positions_overlap_file: 'str | None' = None, exclude_positions_overlap_file: 'str | None' = None, bed_file: 'str | None' = None, exclude_bed_file: 'str | None' = None, thin: 'int | None' = None, mask_file: 'str | None' = None, invert_mask_file: 'str | None' = None, mask_min: 'int | None' = None, snp_ids: 'Sequence[str] | None' = None, snp_ids_file: 'str | None' = None, exclude_ids_file: 'str | None' = None, only_indels: 'bool' = False, remove_indels: 'bool' = False, remove_filtered_all: 'bool' = False, keep_filtered: 'Sequence[str] | None' = None, remove_filtered: 'Sequence[str] | None' = None, keep_info: 'Sequence[str] | None' = None, remove_info: 'Sequence[str] | None' = None, min_maf: 'float | None' = None, max_maf: 'float | None' = None, min_mac: 'int | None' = None, max_mac: 'int | None' = None, min_alleles: 'int | None' = None, max_alleles: 'int | None' = None, min_mean_depth: 'float | None' = None, max_mean_depth: 'float | None' = None, hwe_pvalue: 'float | None' = None, phased_only: 'bool' = False, min_quality: 'float | None' = None, max_missing: 'float | None' = None, keep_individuals: 'Sequence[str] | None' = None, keep_individuals_files: 'Sequence[str] | None' = None, remove_individuals: 'Sequence[str] | None' = None, remove_individuals_files: 'Sequence[str] | None' = None, max_individuals: 'int | None' = None, remove_filtered_genotypes: 'bool' = False, remove_filtered_genotype_flags: 'Sequence[str] | None' = None, min_genotype_quality: 'float | None' = None, min_genotype_depth: 'int | None' = None, max_genotype_depth: 'int | None' = None, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`（VCF、VCF.gz 或 BCF）、`output_prefix`、`output_format`（`"vcf"` 或 `"bcf"`）、`recode_info`，以及上面的过滤选项。列表为 `Sequence[str]`；文件选项必须是已存在的文件（位置列表为制表符分隔的 CHROM、POS；BED 为 CHROM、START、END）。数字在程序启动前检查。
- **输出：** 通用结果字典。`outputs` 为 `<output_prefix>.recode.vcf` 或 `<output_prefix>.recode.bcf`；`sites_kept`/`sites_total` 与 `individuals_kept`/`individuals_total` 来自 VCFtools 日志。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_filter_vcf

out = Path(tempfile.mkdtemp())
kept = vcftools_filter_vcf(
    "toy.vcf", str(out / "pass"), min_maf=0.2, remove_filtered_all=True, min_genotype_quality=20
)
print(kept["status"], kept["sites_kept"], kept["sites_total"], kept["individuals_kept"])
subset = vcftools_filter_vcf(
    "toy.vcf", str(out / "sub"), keep_individuals=["S1", "S2"], remove_indels=True, recode_info=["DP"]
)
print(subset["sites_kept"], subset["individuals_kept"], Path(subset["outputs"][0]).name)
print(Path(subset["outputs"][0]).read_text().splitlines()[-1])
bcf = vcftools_filter_vcf("toy.vcf", str(out / "b"), output_format="bcf", chromosomes=["chr2"])
print(Path(bcf["outputs"][0]).name, bcf["sites_kept"])
```

- **示例输出：**

```text
completed 14 20 6
17 2 sub.recode.vcf
chr2	4800	v20	T	C	67	PASS	DP=154	GT:DP:GQ:PL	0|0:10:30:0,30,300	0|0:38:12:0,30,300
b.recode.bcf 6
```

## 2) `dnakit.variants.vcftools_list_sites()`

- **作用：** 列出通过或未通过一组 VCFtools 过滤条件的位点。
- **处理方法：** 用 `filters` 中的过滤条件加 `--kept-sites`（`which="kept"`，文件 `<output_prefix>.kept.sites`）或 `--removed-sites`（`which="removed"`，文件 `<output_prefix>.removed.sites`）运行 `vcftools`，再读取两列表。
- **API：** `dnakit.variants.vcftools_list_sites(vcf_path: 'str', output_prefix: 'str', *, which: 'str' = 'kept', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`which`（`"kept"` 或 `"removed"`）、`filters`（`vcftools_filter_vcf()` 的键）、`max_rows`。
- **输出：** 带 `table`（列 `CHROM`、`POS`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_list_sites

out = Path(tempfile.mkdtemp())
flt = {"min_maf": 0.3, "remove_filtered_all": True}
kept = vcftools_list_sites("toy.vcf", str(out / "k"), filters=flt)
removed = vcftools_list_sites("toy.vcf", str(out / "r"), which="removed", filters=flt)
print(kept["table"]["columns"], kept["table"]["row_count"], removed["table"]["row_count"])
print(kept["table"]["rows"][:3])
print(Path(removed["outputs"][0]).name)
```

- **示例输出：**

```text
['CHROM', 'POS'] 9 11
[['chr1', 100], ['chr1', 1500], ['chr1', 5000]]
r.removed.sites
```

## 3) `dnakit.variants.vcftools_depth()`

- **作用：** 根据基因型 DP 字段报告每个个体、每个位点或每个基因型的测序深度。
- **处理方法：** 用 `--depth`（`level="individual_mean"`，`.idepth`）、`--site-depth`（`"site_sum"`，`.ldepth`）、`--site-mean-depth`（`"site_mean"`，`.ldepth.mean`）或 `--geno-depth`（`"genotype"`，`.gdepth`）运行 `vcftools`。
- **API：** `dnakit.variants.vcftools_depth(vcf_path: 'str', output_prefix: 'str', *, level: 'str' = 'individual_mean', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`level`、`filters`、`max_rows`。VCF 必须带有 FORMAT/DP。
- **输出：** 带 `table` 的通用结果字典；列随 `level` 而变（`INDV, N_SITES, MEAN_DEPTH`；`CHROM, POS, SUM_DEPTH, SUMSQ_DEPTH`；`CHROM, POS, MEAN_DEPTH, VAR_DEPTH`；或每个个体一列）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_depth

out = Path(tempfile.mkdtemp())
per_sample = vcftools_depth("toy.vcf", str(out / "i"))
print(per_sample["table"]["columns"], per_sample["table"]["rows"][:2])
per_site = vcftools_depth("toy.vcf", str(out / "m"), level="site_mean", max_rows=2)
print(per_site["table"]["columns"], per_site["table"]["rows"], per_site["table"]["truncated"])
genotype = vcftools_depth("toy.vcf", str(out / "g"), level="genotype", max_rows=1)
print(genotype["table"]["columns"], genotype["table"]["rows"][0])
```

- **示例输出：**

```text
['INDV', 'N_SITES', 'MEAN_DEPTH'] [['S1', 20, 22.1], ['S2', 20, 25.5]]
['CHROM', 'POS', 'MEAN_DEPTH', 'VAR_DEPTH'] [['chr1', 100, 15.1667, 80.5667], ['chr1', 800, 20, 28.4]] True
['CHROM', 'POS', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'] ['chr1', 100, 9, 31, 9, 20, 8, 14]
```

## 4) `dnakit.variants.vcftools_site_quality()`

- **作用：** 提取每个位点的 QUAL 值。
- **处理方法：** 运行 `vcftools --site-quality`；表为 `<output_prefix>.lqual`。
- **API：** `dnakit.variants.vcftools_site_quality(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `CHROM`、`POS`、`QUAL`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_site_quality

result = vcftools_site_quality("toy.vcf", str(Path(tempfile.mkdtemp()) / "q"), max_rows=3)
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **示例输出：**

```text
['CHROM', 'POS', 'QUAL']
[['chr1', 100, 38], ['chr1', 800, 55], ['chr1', 1500, 67]]
```

## 5) `dnakit.variants.vcftools_missingness()`

- **作用：** 报告每个个体或每个位点的基因型缺失比例。
- **处理方法：** 运行 `vcftools --missing-indv`（`by="individual"`，`.imiss`）或 `--missing-site`（`by="site"`，`.lmiss`）。
- **API：** `dnakit.variants.vcftools_missingness(vcf_path: 'str', output_prefix: 'str', *, by: 'str' = 'individual', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`by`（`"individual"` 或 `"site"`）、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典（`INDV, N_DATA, N_GENOTYPES_FILTERED, N_MISS, F_MISS` 或 `CHR, POS, N_DATA, N_GENOTYPE_FILTERED, N_MISS, F_MISS`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_missingness

out = Path(tempfile.mkdtemp())
sites = vcftools_missingness("toy.vcf", str(out / "s"), by="site")
print([row for row in sites["table"]["rows"] if row[4] > 0])
people = vcftools_missingness("toy.vcf", str(out / "i"))
print(people["table"]["columns"])
print(people["table"]["rows"][4])
```

- **示例输出：**

```text
[['chr1', 6400, 12, 0, 2, 0.166667]]
['INDV', 'N_DATA', 'N_GENOTYPES_FILTERED', 'N_MISS', 'F_MISS']
['S5', 20, 0, 1, 0.05]
```

## 6) `dnakit.variants.vcftools_singletons()`

- **作用：** 定位 singleton 和 doubleton 等位基因。
- **处理方法：** 运行 `vcftools --singletons`；表为 `<output_prefix>.singletons`。singleton 等位基因在样本中只出现一次，doubleton 在同一个体中出现两次。
- **API：** `dnakit.variants.vcftools_singletons(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `CHROM, POS, SINGLETON/DOUBLETON, ALLELE, INDV`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_singletons

result = vcftools_singletons("toy.vcf", str(Path(tempfile.mkdtemp()) / "s"))
for row in result["table"]["rows"]:
    print(row)
```

- **示例输出：**

```text
['chr1', 3600, 'S', 'T', 'S3']
['chr1', 4300, 'D', 'G', 'S2']
['chr1', 6400, 'S', 'T', 'S3']
['chr2', 3000, 'S', 'A', 'S4']
['chr2', 4800, 'S', 'C', 'S6']
```

## 7) `dnakit.variants.vcftools_indel_length_histogram()`

- **作用：** 统计 VCF 中 indel 的长度分布。
- **处理方法：** 运行 `vcftools --hist-indel-len`；表为 `<output_prefix>.indel.hist`。负数长度为缺失，正数长度为插入，0 为 SNP。
- **API：** `dnakit.variants.vcftools_indel_length_histogram(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `LENGTH, COUNT, PRCT`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_indel_length_histogram

result = vcftools_indel_length_histogram("toy.vcf", str(Path(tempfile.mkdtemp()) / "h"))
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **示例输出：**

```text
['LENGTH', 'COUNT', 'PRCT']
[[-2, 1, 5], [0, 17, 85], [2, 2, 10]]
```

## 8) `dnakit.variants.vcftools_snp_density()`

- **作用：** 沿每条染色体统计不重叠窗口中的 SNP 数。
- **处理方法：** 运行 `vcftools --SNPdensity <bin_size>`；表为 `<output_prefix>.snpden`。
- **API：** `dnakit.variants.vcftools_snp_density(vcf_path: 'str', output_prefix: 'str', bin_size: 'int', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`bin_size`（bp，至少为 1）、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `CHROM, BIN_START, SNP_COUNT, VARIANTS/KB`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_snp_density

result = vcftools_snp_density("toy.vcf", str(Path(tempfile.mkdtemp()) / "d"), 5000)
print(result["table"]["columns"])
for row in result["table"]["rows"]:
    print(row)
```

- **示例输出：**

```text
['CHROM', 'BIN_START', 'SNP_COUNT', 'VARIANTS/KB']
['chr1', 0, 7, 1.4]
['chr1', 5000, 7, 1.4]
['chr2', 0, 6, 1.2]
```

## 9) `dnakit.variants.vcftools_transition_transversion()`

- **作用：** 统计双等位 SNP 的转换（Ts）和颠换（Tv），并按四种分组报告 Ts/Tv 比值。
- **处理方法：** 运行 `vcftools --TsTv <bin_size>`（`mode="distance_bins"`，`.TsTv`）、`--TsTv-by-count`（`"alt_allele_count"`，`.TsTv.count`）、`--TsTv-by-qual`（`"qual_threshold"`，`.TsTv.qual`）或 `--FILTER-summary`（`"filter_category"`，`.FILTER.summary`）。VCFtools 0.1.16 在 `--TsTv-by-qual` 最高 QUAL 阈值那一行（阈值之上没有位点）的三个 `GT` 列里输出未初始化内存，这些单元格没有意义。
- **API：** `dnakit.variants.vcftools_transition_transversion(vcf_path: 'str', output_prefix: 'str', *, mode: 'str' = 'distance_bins', bin_size: 'int | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`mode`、`bin_size`（`distance_bins` 必填）、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典；列与 VCFtools 文件一致（`CHROM, BinStart, SNP_count, Ts/Tv`；`ALT_ALLELE_COUNT, N_Ts, N_Tv, Ts/Tv`；`QUAL_THRESHOLD, N_Ts_LT_QUAL_THRESHOLD, ...`；或 `FILTER, N_VARIANTS, N_Ts, N_Tv, Ts/Tv`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_transition_transversion

out = Path(tempfile.mkdtemp())
bins = vcftools_transition_transversion("toy.vcf", str(out / "b"), bin_size=5000)
print(bins["table"]["columns"], bins["table"]["rows"])
by_filter = vcftools_transition_transversion("toy.vcf", str(out / "f"), mode="filter_category")
print(by_filter["table"]["rows"])
```

- **示例输出：**

```text
['CHROM', 'BinStart', 'SNP_count', 'Ts/Tv'] [['chr1', 0, 6, 5], ['chr1', 5000, 5, 1.5], ['chr2', 0, 6, 0.5]]
[['PASS', 18, 8, 7, 1.14286], ['LowQual', 2, 2, 0, None]]
```

## 10) `dnakit.variants.vcftools_site_nucleotide_diversity()`

- **作用：** 计算每个位点的核苷酸多样性（pi）。
- **处理方法：** 运行 `vcftools --site-pi`；表为 `<output_prefix>.sites.pi`。
- **API：** `dnakit.variants.vcftools_site_nucleotide_diversity(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`filters`、`max_rows`。
- **输出：** 带 `table`（列 `CHROM, POS, PI`）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_site_nucleotide_diversity

result = vcftools_site_nucleotide_diversity("toy.vcf", str(Path(tempfile.mkdtemp()) / "pi"), max_rows=3)
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **示例输出：**

```text
['CHROM', 'POS', 'PI']
[['chr1', 100, 0.530303], ['chr1', 800, 0.409091], ['chr1', 1500, 0.530303]]
```

## 11) `dnakit.variants.vcftools_extract_fields()`

- **作用：** 把每个位点的 INFO 键，或每个基因型的一个 FORMAT 键提取成表。
- **处理方法：** 对每个键运行一次 `vcftools --get-INFO <key>`（`source="info"`，`.INFO`），或运行 `--extract-FORMAT-info <key>`（`source="format"`，`.<key>.FORMAT`）。VCFtools 每次只提取一个 FORMAT 键，所以 `source="format"` 时 `fields` 只有一项。
- **API：** `dnakit.variants.vcftools_extract_fields(vcf_path: 'str', output_prefix: 'str', fields: 'Sequence[str]', *, source: 'str' = 'info', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要 `vcftools` 程序（VCFtools 0.1.16 或更新；`conda install -c bioconda vcftools`）。DNAKit 依次在 `vcftools_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_prefix`、`fields`（由字母、数字、`_` 和 `.` 组成的名称）、`source`、`filters`、`max_rows`。
- **输出：** 带 `table` 的通用结果字典：`CHROM, POS, REF, ALT` 加 INFO 键，或 `CHROM, POS` 加每个个体一列。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_extract_fields

out = Path(tempfile.mkdtemp())
info = vcftools_extract_fields("toy.vcf", str(out / "i"), ["DP", "AF"], max_rows=2)
print(info["table"]["columns"], info["table"]["rows"])
fmt = vcftools_extract_fields("toy.vcf", str(out / "f"), ["GQ"], source="format", max_rows=1)
print(Path(fmt["outputs"][0]).name, fmt["table"]["columns"], fmt["table"]["rows"])
```

- **示例输出：**

```text
['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AF'] [['chr1', 100, 'G', 'T', 91, 0.1], ['chr1', 800, 'G', 'A', 120, 0.25]]
f.GQ.FORMAT ['CHROM', 'POS', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'] [['chr1', 100, 45, 60, 30, 99, 12, 99]]
```
