# VCFtools Perl 脚本：vcf-annotate、vcf-subset、vcf-compare 及其他 vcf-* 工具

这一组 `dnakit.variants` 中的函数是随 [VCFtools](https://github.com/vcftools/vcftools) 一起发布的 Perl 脚本的带类型封装（测试版本 0.1.16）：`vcf-annotate`、`vcf-subset`、`vcf-convert`、`vcf-to-tab`、`vcf-fix-ploidy`、`vcf-shuffle-cols`、`vcf-contrast`、`vcf-compare`、`vcf-indel-stats`、`vcf-phased-join`、`fill-aa`、`fill-fs` 和 `fill-ref-md5`。每个函数用由带类型参数构造的固定参数列表运行一个脚本，不经过 shell；标准输入和输出连接到文件，因此每个函数读一个文件并写一个新文件。同包的 `vcftools` 程序由前两篇文档覆盖（`309_vcftools_filters_site_reports`、`310_vcftools_samples_ld_export`）。

VCFtools 使用 LGPL-3.0 许可；DNAKit 把脚本作为独立进程运行，不复制其任何代码。请自行安装 VCFtools（`conda install -c bioconda vcftools`、`apt install vcftools` 或从源码编译），把脚本放到 PATH 上或传入 `script_path`（单个脚本的完整路径）。脚本需要它们的 `Vcf.pm` 库：当库尚未被找到时，DNAKit 会把 `PERL5LIB` 设为安装目录旁的文件夹（`../share/perl*`）。读取压缩文件的脚本需要 PATH 上有 `tabix`（htslib），`fill-aa`、`fill-fs` 和 `fill-ref-md5` 需要 `samtools`。缺少脚本时函数抛出带安装提示的 `BackendUnavailableError`。

通用行为。`output_path` 必须尚不存在且其文件夹必须存在，运行失败不会留下残缺文件（`status` 为 `failed`，`log_tail` 里是 Perl 的消息）。结果是可 JSON 序列化的字典，含 `status`、`backend`、`version`（`Vcf.pm` 的版本）、`command`、`return_code`、`outputs`、`log_tail`、`output_path`，VCF 输出还有 `records`（数据行数）。若干脚本遍历 Perl 哈希，因此 INFO 键、头部行和每个样本行的顺序在不同运行之间不稳定。

示例从 `dnakit.variants` 导入，用 pysam 对需要的文件做 bgzip 压缩和 tabix 索引，并写入新的临时目录。

## 示例数据

唯一需要声明的输入是 `toy.vcf`：6 个已定相的样本（S1 到 S6），`chr1`（14 个）和 `chr2`（6 个）共 20 个位点，3 个 indel、2 个 `LowQual` 位点、4 个 singleton、1 个 doubleton、1 个缺失基因型，FORMAT 字段为 GT、DP、GQ 和 PL。需要参考 FASTA 的示例会写出一个碱基与 `toy.vcf` 的 REF 等位基因一致的文件。

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

## 1) `dnakit.variants.vcf_annotate()`

- **作用：** 用 `vcf-annotate` 添加 FILTER 标记、填充 AC、AN、HWE 和 TYPE、规范化等位基因、删除标签或硬过滤位点。
- **处理方法：** 用 `-f <filters>`（内置位点过滤，例如 `+/Qual=40/MinDP=5`：`+` 以默认值启用所有过滤，`-X` 排除过滤 X；不支持用户自定义过滤文件）、`-H`（`hard_filter`：删除 FILTER 不是 PASS 或 "." 的位点）、`-n`（`normalize_alleles`：压缩 REF 和 ALT，例如 TA,TAA 变为 T,TA）、`-r <tags>`（`remove_tags`：`ID`、`QUAL`、`FILTER`、`INFO/DP`、`FORMAT/DP`、`FILTER/MinDP` 或单独的标签）、`--fill-AC-AN`（`fill_ac_an`，等同于 `fill-an-ac` 脚本）、`--fill-HWE`（`fill_hwe`：同时加入 HWE p 值）和 `--fill-type`（`fill_type`：INFO/TYPE 为 snp、del、ins、mnp 或 complex，即 `Vcf::event_type` 的分类）运行 `vcf-annotate`。至少要选一个操作。
- **API：** `dnakit.variants.vcf_annotate(vcf_path: 'str', output_path: 'str', *, filters: 'str | None' = None, hard_filter: 'bool' = False, normalize_alleles: 'bool' = False, remove_tags: 'Sequence[str] | None' = None, fill_ac_an: 'bool' = False, fill_hwe: 'bool' = False, fill_type: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`（VCF 或 VCF.gz）、`output_path`、`filters`、`hard_filter`、`normalize_alleles`、`remove_tags`、`fill_ac_an`、`fill_hwe`、`fill_type`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典，注释后的 VCF 位于 `output_path`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_annotate

out = Path(tempfile.mkdtemp())
flagged = vcf_annotate("toy.vcf", str(out / "flag.vcf"), filters="+/Qual=40")
rows = [line.split("\t") for line in Path(flagged["output_path"]).read_text().splitlines() if line[0] != "#"]
print(flagged["status"], [row[6] for row in rows][:6])
hard = vcf_annotate("toy.vcf", str(out / "hard.vcf"), filters="+/Qual=40", hard_filter=True)
print(hard["records"])
filled = vcf_annotate("toy.vcf", str(out / "fill.vcf"), fill_ac_an=True, fill_type=True)
print(sorted(Path(filled["output_path"]).read_text().splitlines()[-1].split("\t")[7].split(";")))
bare = vcf_annotate("toy.vcf", str(out / "bare.vcf"), remove_tags=["INFO/AF", "QUAL"])
print(Path(bare["output_path"]).read_text().splitlines()[-1].split("\t")[5:8])
```

- **示例输出：**

```text
completed ['Qual', 'PASS', 'PASS', 'Qual', 'PASS', 'LowQual;Qual']
12
['AC=1', 'AF=0.1', 'AN=12', 'DP=154', 'TYPE=snp']
['.', 'PASS', 'DP=154']
```

## 2) `dnakit.variants.vcf_subset()`

- **作用：** 用 `vcf-subset` 选择 VCF 的样本和变异类别。
- **处理方法：** 用 `-c`（`samples`，或每行一个名称的 `samples_file`）、`-e`（`exclude_ref_only`：删除所选样本中没有 ALT 等位基因的位点）、`-p`（`private_only`：保留只有所选样本带有变异等位基因的位点）、`-a`（`trim_alt_alleles`：删除所选样本未使用的 ALT 等位基因）、`-t`（`variant_types`，取自 `ref`、`SNPs`、`indels`、`MNPs`、`other`）、`-r`（`replace_with_ref`：排除的类型变为参考等位基因而不是 "."）、`-u`（`keep_uncalled`）和 `-f`（`force`：样本不存在时继续）运行 `vcf-subset`。
- **API：** `dnakit.variants.vcf_subset(vcf_path: 'str', output_path: 'str', *, samples: 'Sequence[str] | None' = None, samples_file: 'str | None' = None, exclude_ref_only: 'bool' = False, private_only: 'bool' = False, trim_alt_alleles: 'bool' = False, variant_types: 'Sequence[str] | None' = None, replace_with_ref: 'bool' = False, keep_uncalled: 'bool' = False, force: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`samples` 或 `samples_file`（二选一）、`exclude_ref_only`、`private_only`、`trim_alt_alleles`、`variant_types`、`replace_with_ref`、`keep_uncalled`、`force`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_subset

out = Path(tempfile.mkdtemp())
two = vcf_subset("toy.vcf", str(out / "two.vcf"), samples=["S1", "S2"])
header = [x for x in Path(two["output_path"]).read_text().splitlines() if x.startswith("#CHROM")][0]
print(two["records"], header.split("\t")[9:])
indels = vcf_subset("toy.vcf", str(out / "indel.vcf"), samples=["S1", "S2"], variant_types=["indels"], exclude_ref_only=True)
print(indels["records"])
private = vcf_subset("toy.vcf", str(out / "private.vcf"), samples=["S3"], private_only=True)
print(private["records"])
```

- **示例输出：**

```text
20 ['S1', 'S2']
3
2
```

## 3) `dnakit.variants.vcf_convert()`

- **作用：** 用 `vcf-convert` 在 VCF 格式版本 4.0、4.1 和 4.2 之间转换。
- **处理方法：** 运行 `vcf-convert -v <target_version>`，VCF 作为标准输入，给出时加 `-r <reference_fasta>`。只有文件含 indel 时才需要参考序列（用 `samtools faidx` 建立索引的 FASTA）。脚本把降级转换标为实验性功能。
- **API：** `dnakit.variants.vcf_convert(vcf_path: 'str', output_path: 'str', *, target_version: 'str' = '4.2', reference_fasta: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`target_version`（`"4.0"`、`"4.1"` 或 `"4.2"`，默认 `"4.2"`）、`reference_fasta`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典；含 indel 又没有参考序列时 `status` 为 `failed`，脚本消息在 `log_tail` 中。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_convert

lines = Path("toy.vcf").read_text().splitlines()
snps = [x for x in lines if x[0] == "#" or (len(x.split("\t")[3]) == 1 and len(x.split("\t")[4]) == 1)]
Path("snps.vcf").write_text("\n".join(snps) + "\n")
out = Path(tempfile.mkdtemp())
converted = vcf_convert("snps.vcf", str(out / "v40.vcf"), target_version="4.0")
print(converted["status"], converted["records"], Path(converted["output_path"]).read_text().splitlines()[0])
refused = vcf_convert("toy.vcf", str(out / "indels.vcf"), target_version="4.0")
print(refused["status"], "Indels present" in refused["log_tail"])
```

- **示例输出：**

```text
completed 17 ##fileformat=VCFv4.0
failed True
```

## 4) `dnakit.variants.vcf_to_tab()`

- **作用：** 用 `vcf-to-tab` 把 VCF 转换成以碱基表示基因型的制表符分隔表。
- **处理方法：** 运行 `vcf-to-tab`，VCF 作为标准输入；`iupac=True` 时加 `-i`（杂合基因型用单字母 IUPAC 码表示；遇到 indel 脚本会停止）。列为 CHROM、POS、REF 以及每个样本一列碱基。脚本自己说明将来应改用 `bcftools query`。
- **API：** `dnakit.variants.vcf_to_tab(vcf_path: 'str', output_path: 'str', *, iupac: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`iupac`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records`（数据行数）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_to_tab

out = Path(tempfile.mkdtemp())
table = vcf_to_tab("toy.vcf", str(out / "t.tab"))
lines = Path(table["output_path"]).read_text().splitlines()
print(table["status"], table["records"])
print(lines[0])
print(lines[1])
iupac = vcf_to_tab("toy.vcf", str(out / "t2.tab"), iupac=True)
print(iupac["status"], "IUPAC" in iupac["log_tail"])
```

- **示例输出：**

```text
completed 20
#CHROM	POS	REF	S1	S2	S3	S4	S5	S6
chr1	100	G	T/G	G/T	G/T	G/T	T/G	G/G
failed True
```

## 5) `dnakit.variants.vcf_fix_ploidy()`

- **作用：** 用 `vcf-fix-ploidy` 按样本性别设置 X、Y 和 MT 染色体的基因型倍性。
- **处理方法：** 运行 `vcf-fix-ploidy -s <samples_file>`，VCF 作为标准输入，另加 `-a <M|F>`（`assumed_sex`，用于列表中缺失的样本）、`-l`（`fix_likelihoods`：增加或删除杂合似然）和 `-p <ploidy_file>`（替换默认人类定义的倍性定义）。默认定义让男性在拟常染色体区之外的 X 为单倍体，女性没有 Y、男性 Y 为单倍体，MT 为单倍体。
- **API：** `dnakit.variants.vcf_fix_ploidy(vcf_path: 'str', output_path: 'str', samples_file: 'str', *, assumed_sex: 'str | None' = None, fix_likelihoods: 'bool' = False, ploidy_file: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`samples_file`（每行 `样本名 M|F`）、`assumed_sex`、`fix_likelihoods`、`ploidy_file`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_fix_ploidy

lines = Path("toy.vcf").read_text().splitlines()
x_lines = [x.replace("ID=chr1", "ID=X") if x[0] == "#" else x.replace("chr1", "X", 1) for x in lines]
Path("x.vcf").write_text("\n".join(x_lines) + "\n")
Path("sex.txt").write_text("S1 M\nS2 F\nS3 M\nS4 F\nS5 M\nS6 F\n")
fixed = vcf_fix_ploidy("x.vcf", str(Path(tempfile.mkdtemp()) / "fixed.vcf"), "sex.txt")
row = Path(fixed["output_path"]).read_text().splitlines()[-20].split("\t")
print(fixed["status"], fixed["records"])
print([cell.split(":")[0] for cell in row[9:15]])
```

- **示例输出：**

```text
completed 20
['1', '0|1', '0', '0|1', '1', '0|0']
```

## 6) `dnakit.variants.vcf_shuffle_columns()`

- **作用：** 用 `vcf-shuffle-cols` 把 VCF 的样本列重排成模板 VCF 的顺序。
- **处理方法：** 运行 `vcf-shuffle-cols -t <template_vcf_path> <vcf_path>`。两个文件含相同的样本，通常经 bgzip 压缩并做了 tabix 索引。
- **API：** `dnakit.variants.vcf_shuffle_columns(vcf_path: 'str', template_vcf_path: 'str', output_path: 'str', *, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`template_vcf_path`、`output_path`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典。
- **示例代码：**

```python
import shutil
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_shuffle_columns

lines = Path("toy.vcf").read_text().splitlines()
columns = [x for x in lines if x.startswith("#CHROM")][0].split("\t")
order = [*range(9), *range(len(columns) - 1, 8, -1)]
template = [
    "\t".join(x.split("\t")[i] for i in order) if x[:2] != "##" else x for x in lines
]
Path("template.vcf").write_text("\n".join(template) + "\n")
shutil.copy("toy.vcf", "toy_c.vcf")
for name in ("template.vcf", "toy_c.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_shuffle_columns("toy_c.vcf.gz", "template.vcf.gz", str(Path(tempfile.mkdtemp()) / "shuffled.vcf"))
header = [x for x in Path(result["output_path"]).read_text().splitlines() if x.startswith("#CHROM")][0]
print(result["status"], result["records"], header.split("\t")[9:])
```

- **示例输出：**

```text
completed 20 ['S6', 'S5', 'S4', 'S3', 'S2', 'S1']
```

## 7) `dnakit.variants.vcf_contrast()`

- **作用：** 用 `vcf-contrast` 寻找只在部分样本中出现的基因型和等位基因。
- **处理方法：** 运行 `vcf-contrast +<变异样本> -<背景样本> [-d <min_depth>] [-f] [-n] <vcf>`。当变异样本的基因型或等位基因在所有背景样本中都不存在时，在 INFO 中加入 NOVELGT、NOVELAL 和 NOVELTY。`min_depth` 为背景样本的最小深度（`-d`），`apply_filters` 跳过 FILTER 不是 PASS 或 "." 的位点（`-f`），`novel_sites_only` 只输出含新基因型的记录（`-n`）。单倍体基因型按纯合二倍体处理。
- **API：** `dnakit.variants.vcf_contrast(vcf_path: 'str', output_path: 'str', variant_samples: 'Sequence[str]', background_samples: 'Sequence[str]', *, min_depth: 'int | None' = None, apply_filters: 'bool' = False, novel_sites_only: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`（VCF 或 bgzip 压缩的 VCF）、`output_path`、`variant_samples`、`background_samples`、`min_depth`、`apply_filters`、`novel_sites_only`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_contrast

result = vcf_contrast(
    "toy.vcf", str(Path(tempfile.mkdtemp()) / "novel.vcf"), ["S1"], ["S2", "S3"], novel_sites_only=True
)
print(result["status"], result["records"])
for line in Path(result["output_path"]).read_text().splitlines():
    if line[0] != "#":
        fields = line.split("\t")
        print(fields[0], fields[1], sorted(fields[7].split(";")))
```

- **示例输出：**

```text
completed 5
chr1 800 ['AF=0.25', 'DP=120', 'NOVELGT=S1', 'NOVELTY=30']
chr1 1500 ['AF=0.5', 'DP=99', 'NOVELGT=S1', 'NOVELTY=30']
chr1 2900 ['AF=0.25', 'DP=172', 'NOVELAL=S1', 'NOVELTY=30']
chr2 300 ['AF=0.4', 'DP=204', 'NOVELGT=S1', 'NOVELTY=30']
chr2 3900 ['AF=0.4', 'DP=117', 'NOVELGT=S1', 'NOVELTY=30']
```

## 8) `dnakit.variants.vcf_compare()`

- **作用：** 用 `vcf-compare` 比较两个或多个 VCF 文件的位点和基因型。
- **处理方法：** 用 `-g`（`compare_genotypes`）、`--ignore-indels`（`ignore_indels`）、`-a`（`apply_filters`：忽略 FILTER 不是 PASS 或 "." 的位点）、`-r`（`regions`）、`-s`（`samples`）、`-w`（`window`：相距不超过该 bp 数的记录视为匹配）和 `-R`（`reference_fasta`：比较序列而不只是位置）运行 `vcf-compare`。文件必须经 bgzip 压缩并有 tabix 索引。报告写入 `output_path`，并解析为 `venn`（每种文件组合的位点数）、`summary`（`SN` 行）和 `genotype_classes`（`GS` 行）。
- **API：** `dnakit.variants.vcf_compare(vcf_paths: 'Sequence[str]', output_path: 'str', *, compare_genotypes: 'bool' = False, apply_filters: 'bool' = False, ignore_indels: 'bool' = False, regions: 'Sequence[str] | None' = None, samples: 'Sequence[str] | None' = None, window: 'int | None' = None, reference_fasta: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_paths`（至少两个）、`output_path`、`compare_genotypes`、`apply_filters`、`ignore_indels`、`regions`、`samples`、`window`、`reference_fasta`、`script_path`、`timeout_seconds`。
- **输出：** 带 `venn`、`summary` 和 `genotype_classes` 的通用结果字典；完整报告在 `output_path`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_compare

lines = Path("toy.vcf").read_text().splitlines()
head = [x for x in lines if x[0] == "#"]
body = [x for x in lines if x[0] != "#"]
Path("all.vcf").write_text("\n".join(head + body) + "\n")
Path("half.vcf").write_text("\n".join(head + body[::2]) + "\n")
for name in ("all.vcf", "half.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_compare(["all.vcf.gz", "half.vcf.gz"], str(Path(tempfile.mkdtemp()) / "cmp.txt"), compare_genotypes=True)
print(result["status"])
for row in sorted(result["venn"], key=lambda r: len(r["files"])):
    print(row["sites"], [Path(name.split(" ")[0]).name for name in row["files"]])
print(result["summary"]["Number of REF matches"], result["summary"]["Number of samples in GT comparison"])
```

- **示例输出：**

```text
completed
10 ['all.vcf.gz']
10 ['all.vcf.gz', 'half.vcf.gz']
10 6
```

## 9) `dnakit.variants.vcf_indel_stats()`

- **作用：** 用 `vcf-indel-stats` 根据外显子坐标统计 VCF 中的框内和移码 indel。
- **处理方法：** 运行 `vcf-indel-stats -e <exons_file>`，VCF 作为标准输入。脚本把落在外显子内的 indel 按长度是否为 3 的倍数分类，并报告框内比例。
- **API：** `dnakit.variants.vcf_indel_stats(vcf_path: 'str', output_path: 'str', exons_file: 'str', *, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`exons_file`（制表符分隔的 CHROM、FROM、TO；1 起始，含两端）、`script_path`、`timeout_seconds`。
- **输出：** 带 `statistics`（报告各行组成的字典）的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_indel_stats

Path("exons.txt").write_text("chr1\t2000\t3000\nchr1\t5600\t6000\n")
result = vcf_indel_stats("toy.vcf", str(Path(tempfile.mkdtemp()) / "stats.txt"), "exons.txt")
print(result["status"], result["statistics"])
```

- **示例输出：**

```text
completed {'total': 3, 'in-frame': 0, 'frameshift': 2, 'ratio': 1.0}
```

## 10) `dnakit.variants.vcf_phased_join()`

- **作用：** 用 `vcf-phased-join` 把相互重叠、已预先定相的 VCF 片段拼接成一个 VCF。
- **处理方法：** 运行 `vcf-phased-join -o <output_path> [-j <min_join_quality>] [-q <min_phase_quality>] <chunks>`。重叠区中的杂合调用决定每个片段的相位。`min_join_quality` 是拼接相位块的质量阈值（脚本默认 10），`min_phase_quality` 会在 PQ 值更低处断开预先定相的片段（脚本默认 0.6）。片段按基因组顺序给出，通常经 bgzip 压缩并有 tabix 索引。
- **API：** `dnakit.variants.vcf_phased_join(vcf_paths: 'Sequence[str]', output_path: 'str', *, min_join_quality: 'float | None' = None, min_phase_quality: 'float | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_paths`（至少两个片段）、`output_path`、`min_join_quality`、`min_phase_quality`、`script_path`、`timeout_seconds`。
- **输出：** 带拼接后 VCF 的 `records` 的通用结果字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_phased_join

lines = Path("toy.vcf").read_text().splitlines()
head = [x for x in lines if x[0] == "#"]
body = [x for x in lines if x[0] != "#"]
Path("chunk1.vcf").write_text("\n".join(head + body[:12]) + "\n")
Path("chunk2.vcf").write_text("\n".join(head + body[7:14]) + "\n")
for name in ("chunk1.vcf", "chunk2.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_phased_join(["chunk1.vcf.gz", "chunk2.vcf.gz"], str(Path(tempfile.mkdtemp()) / "joined.vcf"))
print(result["status"], result["records"])
```

- **示例输出：**

```text
completed 14
```

## 11) `dnakit.variants.vcf_fill_ancestral_allele()`

- **作用：** 用 `fill-aa` 从按染色体分开的祖先序列 FASTA 文件中添加祖先等位基因（INFO/AA）。
- **处理方法：** 运行 `fill-aa -a <ancestral_prefix> -t <variant_types>`，VCF 作为标准输入。文件名为 `<ancestral_prefix><染色体>.fa.gz`，必须压缩并用 `samtools faidx` 建立索引。`variant_types` 取自 `all`、`indel`、`ref` 和 `snp`。
- **API：** `dnakit.variants.vcf_fill_ancestral_allele(vcf_path: 'str', output_path: 'str', ancestral_prefix: 'str', *, variant_types: 'Sequence[str]' = ('all',), script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`ancestral_prefix`、`variant_types`（默认 `("all",)`）、`script_path`、`timeout_seconds`。脚本需要 PATH 上有 `samtools`。
- **输出：** 带 `records` 的通用结果字典；`log_tail` 里是脚本对已填充位点的汇总。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_ancestral_allele

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
for name, bases in sequences.items():
    Path(f"anc_{name}.fa").write_text(f">{name}\n" + "".join(bases) + "\n")
    pysam.tabix_compress(f"anc_{name}.fa", f"anc_{name}.fa.gz", force=True)
    pysam.faidx(f"anc_{name}.fa.gz")
result = vcf_fill_ancestral_allele("toy.vcf", str(Path(tempfile.mkdtemp()) / "aa.vcf"), "anc_")
print(result["status"], result["records"])
print(Path(result["output_path"]).read_text().splitlines()[-1].split("\t")[7])
```

- **示例输出：**

```text
completed 20
DP=154;AF=0.1;AA=T
```

## 12) `dnakit.variants.vcf_fill_flanking_sequence()`

- **作用：** 用 `fill-fs` 在每个变异周围添加侧翼参考序列（INFO/FS）。
- **处理方法：** 运行 `fill-fs -r <reference_fasta> -l <length> [-m <mask_character> -b <bed 掩蔽>] [-m <mask_character> -v <vcf 掩蔽>] [-m <mask_character> -c <cluster_distance>] <vcf>`。参考序列是用 `samtools faidx` 建立索引的 FASTA。两侧各取 `length` 个碱基（脚本默认 100）；`bed_mask_files` 和 `vcf_mask_files`（带 tabix 索引）用 `mask_character`（单个字符或表示小写的 `lc`，默认 `N`）掩蔽侧翼中的区域和已知变异，`cluster_distance` 对成簇变异做自掩蔽。脚本要求 `-m` 位于每个掩蔽选项之前，所以 DNAKit 在每个掩蔽前重复它；没有掩蔽时它不起作用。
- **API：** `dnakit.variants.vcf_fill_flanking_sequence(vcf_path: 'str', output_path: 'str', reference_fasta: 'str', *, length: 'int' = 100, mask_character: 'str' = 'N', bed_mask_files: 'Sequence[str] | None' = None, vcf_mask_files: 'Sequence[str] | None' = None, cluster_distance: 'int | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`、`output_path`、`reference_fasta`、`length`、`mask_character`、`bed_mask_files`、`vcf_mask_files`、`cluster_distance`、`script_path`、`timeout_seconds`。
- **输出：** 带 `records` 的通用结果字典；INFO/FS 的形式为 `左侧[REF/ALT]右侧`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_flanking_sequence

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
Path("ref.fa").write_text("".join(f">{n}\n{''.join(b)}\n" for n, b in sequences.items()))
pysam.faidx("ref.fa")
result = vcf_fill_flanking_sequence("toy.vcf", str(Path(tempfile.mkdtemp()) / "fs.vcf"), "ref.fa", length=5)
rows = [x.split("\t") for x in Path(result["output_path"]).read_text().splitlines() if x[0] != "#"]
print(result["status"], result["records"])
print(rows[0][7])
print(rows[1][7])
```

- **示例输出：**

```text
completed 20
DP=91;AF=0.1;FS=AAAAA[G/T]AAAAA
DP=120;AF=0.25;FS=AAAAA[G/A]AAAAA
```

## 13) `dnakit.variants.vcf_fill_reference_md5()`

- **作用：** 用 `fill-ref-md5` 把参考序列各 contig 的 MD5 写入 VCF 头部。
- **处理方法：** 运行 `fill-ref-md5 -r <reference_fasta> [-d <dictionary_path>] [-i AS:..,SP:..,TX:..] <vcf.gz> <out.vcf.gz>`。脚本按 VCF 4.1 的推荐加入 `reference` 和 `contig` 头部行，并用 `tabix -r` 替换头部。输入必须经 bgzip 压缩并有 tabix 索引。`dictionary_path` 缓存已计算的 MD5，以追加模式打开，所以文件可以已存在。脚本会忽略 `tabix` 的失败，因此 DNAKit 把空输出报告为 `failed`。
- **API：** `dnakit.variants.vcf_fill_reference_md5(vcf_path: 'str', output_path: 'str', reference_fasta: 'str', *, dictionary_path: 'str | None' = None, assembly: 'str | None' = None, species: 'str | None' = None, taxonomy: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖：** 需要本节的 VCFtools Perl 脚本（VCFtools 0.1.16 或更新）、perl，以及文中说明的 `tabix`/`samtools`。DNAKit 依次在 `script_path`、PATH 和当前 Python 解释器旁查找；缺少时在写入任何文件之前抛出 `BackendUnavailableError`。
- **输入：** `vcf_path`（`.vcf.gz`）、`output_path`（新的 `.vcf.gz`）、`reference_fasta`、`dictionary_path`、`assembly`、`species`、`taxonomy`（AS、SP 和 TX 头部值；不含逗号）、`script_path`、`timeout_seconds`。需要 PATH 上有 `tabix` 和 `samtools`。
- **输出：** 通用结果字典；新的压缩 VCF 位于 `output_path`。
- **示例代码：**

```python
import gzip
import shutil
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_reference_md5

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
Path("ref.fa").write_text("".join(f">{n}\n{''.join(b)}\n" for n, b in sequences.items()))
pysam.faidx("ref.fa")
shutil.copy("toy.vcf", "toy_md5.vcf")
pysam.tabix_index("toy_md5.vcf", preset="vcf", force=True)
result = vcf_fill_reference_md5("toy_md5.vcf.gz", str(Path(tempfile.mkdtemp()) / "md5.vcf.gz"), "ref.fa", assembly="toy")
print(result["status"])
for line in gzip.open(result["output_path"], "rt").read().splitlines():
    if line.startswith("##contig"):
        print(line)
```

- **示例输出：**

```text
completed
##contig=<ID=chr1,assembly=toy,length=20000,md5=93e08f94707a244b4a8ccfd891707c3d>
##contig=<ID=chr2,assembly=toy,length=20000,md5=6a1cea00ce86defc019f8e80d5492f86>
```
