# GATK 的 VCF 工具：头部、字典、liftover、拆分与转换

`dnakit.variants.gatk_p02_cli_functions` 中的这些函数通过固定参数白名单（不经 shell）调用 GATK 命令行工具（测试版本 4.7.0.0）。需要安装 GATK 与 Java 17（`conda install -c bioconda gatk4` 或发布包），并把 `gatk` 放入 PATH 或传入 `gatk_path`。输出写入新路径（不覆盖已有文件），GATK 生成的索引文件列于 `index_files`。结果为字典，含 `status`（仅当退出码为 0 且输出存在时为 `completed`）、`version`、`command`、`outputs`、`log_tail` 及各工具特有的值。未设置的参数沿用 GATK 默认值。找不到 `gatk` 时抛出带安装提示的 `BackendUnavailableError`。

## 示例数据

示例在同一目录中按顺序运行，后面的示例会使用前面生成的文件。需要 `gatk` 在 PATH 中。

**示例文件：** `calls.vcf`

```text
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##contig=<ID=chr1,length=600>
##contig=<ID=chr2,length=300>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	50	.	G	A	50	.	QD=1.5;DP=40	GT:DP	0/1:20	0/0:20
chr1	120	.	T	A	50	.	QD=18.0;DP=40	GT:DP	1/1:20	0/1:20
chr1	200	.	CTA	C	50	.	QD=12.0;DP=40	GT:DP	0/1:20	0/1:20
chr1	310	.	A	C	50	.	QD=25.0;DP=40	GT:DP	0/1:20	1/1:20
chr1	400	.	C	CAC	50	.	QD=9.0;DP=40	GT:DP	0/1:20	0/0:20
chr2	30	.	A	C	50	.	QD=21.0;DP=40	GT:DP	1/1:20	0/0:20
chr2	150	.	A	C	50	.	QD=3.0;DP=40	GT:DP	0/1:20	0/1:20
```

**示例文件：** `ref.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:300
```

**示例文件：** `shift.chain`

```text
chain 1000 chr1 600 + 5 600 chr1 595 + 0 595 1
595

chain 900 chr2 300 + 0 300 chr2 300 + 0 300 2
300
```

**示例文件：** `new.fa`

```text
>chr1
TTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATTTGGTC
TTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCGAGGTT
TCGATTTCCGCTATTTATCCTGATTTCCTGACCTGATTGTCGTAATTTTGTAGTCTAAAACAAGCAAGAT
ATTAAATCGTCAGCATCATCAGTGCAATAGTACATCATTATAACGTATTTAAAGTTTTAACTGTTGCAAC
TCATTTCTTTTACAAAATCGGATAACGTATTGTGTCTAAGTATCTGGTTATACAATCGTTGATGCTAATT
GCTTAAAAAAGGTTTAAATTGGACCACAGCCGGATATCGGCAAGCTGGCTTTCCCGTTCTATCCTACTTA
CTAGTATTGGTCCAGACTCTCTGTAAGCATCCGCACTTCTTTGATAACAATATATTTCTTTTTAACGGAA
ACGTAGGATTTAATCAATTGTTATTAGAGGACCCTCTACTTAAGTTATCATGACAAAAGCCCTGCCTGTA
ATCTTTATATATATGTGAATACGCTCGTATATTAT
>chr2
GTGCTGGCTTACTTGCGCATTGAGGGACAACGGTAGTACCCCCGGGGATCGTGGCATTACGGGCCGTTAG
GCTTCCTTGCTGTCATATACCGTGCTGACCTTCCAATCAGGTACAGATCTTCTAGAGGCTTTCATCACGG
CGACAACGGAGCTGAGGGGTAGCGGGCATTAACCGCGCGTCCACGGTAAGTTACAGGTTATGCGAGGATT
ATCGCTGAGTAGATTCAGAGACCCCGAAGTGAACTGAGGTCAGATAGGACTGACGCCGTATCTTGAACAA
ATGCTGGCGTCCGATTGACA
```

**示例文件：** `new.fa.fai`

```text
chr1	595	6	70	71
chr2	300	616	70	71
```

**示例文件：** `new.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:595
@SQ	SN:chr2	LN:300
```

## 1) `fix_vcf_header_gatk()`

- **作用：** 替换 VCF 头部，或补全记录用到但头部缺失的 INFO/FORMAT/FILTER 定义。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk FixVcfHeader`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.fix_vcf_header_gatk(vcf: 'str', output_path: 'str', *, header: 'str | None' = None, enforce_same_samples: 'bool' = True, check_first_n_records: 'int | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF，可选用于替换头部的 VCF。
- **输出：** 修复后的 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.fix_vcf_header_gatk("calls.vcf", "fixed.vcf")
print(r["status"], r["summary"]["records"], r["summary"]["samples"])
```

- **示例输出：**

```text
completed 7 ['S1', 'S2']
```

## 2) `update_vcf_sequence_dictionary_gatk()`

- **作用：** 用字典、VCF、BAM 或参考序列更新 VCF 的 contig 行。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk UpdateVCFSequenceDictionary`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.update_vcf_sequence_dictionary_gatk(vcf: 'str', output_path: 'str', *, source_dictionary: 'str | None' = None, reference: 'str | None' = None, replace: 'bool' = False, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF 与字典来源（或参考序列）。
- **输出：** 更新后的 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.update_vcf_sequence_dictionary_gatk("calls.vcf", "updated.vcf", source_dictionary="ref.dict", replace=True)
print(r["status"], [x for x in Path("updated.vcf").read_text().splitlines() if x.startswith("##contig")])
```

- **示例输出：**

```text
completed ['##contig=<ID=chr1,length=600>', '##contig=<ID=chr2,length=300>']
```

## 3) `liftover_vcf_gatk()`

- **作用：** 用 chain 文件把 VCF 转换到另一参考基因组版本。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk LiftoverVcf`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.liftover_vcf_gatk(vcf: 'str', chain_file: 'str', reference: 'str', output_path: 'str', reject_path: 'str', *, min_match: 'float | None' = None, recover_swapped_ref_alt: 'bool' = False, write_original_position: 'bool' = False, write_original_alleles: 'bool' = False, tags_to_drop: 'list[str] | None' = None, tags_to_reverse: 'list[str] | None' = None, allow_missing_fields_in_header: 'bool' = False, disable_sort: 'bool' = False, warn_on_missing_contig: 'bool' = False, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF、chain 文件与目标参考 FASTA（含 .fai 与 .dict）。
- **输出：** 已转换与被拒绝的 VCF 及其摘要。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.liftover_vcf_gatk("calls.vcf", "shift.chain", "new.fa", "lifted.vcf", "rejected.vcf")
rows = [x.split("\t")[:2] for x in Path("lifted.vcf").read_text().splitlines() if not x.startswith("#")]
print(r["status"], r["summary"]["records"], r["reject_summary"]["records"], rows[:3])
```

- **示例输出：**

```text
completed 7 0 [['chr1', '45'], ['chr1', '115'], ['chr1', '195']]
```

## 4) `make_sites_only_vcf_gatk()`

- **作用：** 去除基因型列但保留位点级信息。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk MakeSitesOnlyVcf`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.make_sites_only_vcf_gatk(vcf: 'str', output_path: 'str', *, samples: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF/BCF，可选使用的样本。
- **输出：** 仅含位点的 VCF 与 `summary`（样本列表为空）。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.make_sites_only_vcf_gatk("calls.vcf", "sites.vcf")
print(r["status"], r["summary"]["records"], r["summary"]["samples"])
```

- **示例输出：**

```text
completed 7 []
```

## 5) `split_vcf_snps_indels_gatk()`

- **作用：** 把 VCF 拆分为头部一致的 SNP 与 indel 两个文件。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk SplitVcfs`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.split_vcf_snps_indels_gatk(vcf: 'str', snp_output_path: 'str', indel_output_path: 'str', *, strict: 'bool' = True, sequence_dictionary: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF/BCF 以及遇到其他类型是否报错（`strict`）。
- **输出：** SNP 与 indel VCF 及其摘要。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.split_vcf_snps_indels_gatk("calls.vcf", "snps.vcf", "indels.vcf")
print(r["status"], r["snp_summary"]["records"], r["indel_summary"]["records"])
```

- **示例输出：**

```text
completed 5 2
```

## 6) `convert_vcf_format_gatk()`

- **作用：** VCF 与 BCF 互转（格式由扩展名决定）。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk VcfFormatConverter`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.convert_vcf_format_gatk(input_path: 'str', output_path: 'str', *, require_index: 'bool' = True, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** 扩展名为 .vcf、.vcf.gz 或 .bcf 的输入与输出路径。
- **输出：** 转换后的文件（VCF 输出附 `summary`）。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

a = g.convert_vcf_format_gatk("calls.vcf", "calls.bcf", require_index=False)
b = g.convert_vcf_format_gatk("calls.bcf", "back.vcf", require_index=False)
print(a["status"], b["status"], b["summary"]["records"])
```

- **示例输出：**

```text
completed completed 7
```

## 7) `vcf_to_interval_list_gatk()`

- **作用：** 由 VCF 记录生成 Picard 区间列表。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk VcfToIntervalList`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.vcf_to_interval_list_gatk(vcf: 'str', output_path: 'str', *, include_filtered: 'bool' = False, variant_id_method: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** VCF/BCF、是否包含已过滤记录以及区间命名方式。
- **输出：** 区间列表文件及 `intervals` 行 `[contig, start, end, strand, name]`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.vcf_to_interval_list_gatk("calls.vcf", "calls.interval_list")
print(r["status"], len(r["intervals"]), r["intervals"][:2])
```

- **示例输出：**

```text
completed 7 [['chr1', 50, 50, '+', 'interval-1'], ['chr1', 120, 120, '+', 'interval-2']]
```

## 8) `print_variants_spark()`

- **作用：** 通过 GATK Spark 引擎把 VCF 变异重写为单个文件（4.7.0.0 需仅含位点的输入）。
- **处理方法：** 用带类型的参数构造固定参数列表并运行 `gatk PrintVariantsSpark`，只传递你设置的参数。
- **API：** `dnakit.variants.gatk_p02_cli_functions.print_variants_spark(vcf: 'str', output_path: 'str', *, reference: 'str | None' = None, variant_shard_size: 'int | None' = None, use_nio: 'bool' = False, spark_master: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **输入：** 仅含位点的 VCF，可选分片大小与 Spark master。
- **输出：** 输出 VCF 与 `summary`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.print_variants_spark("sites.vcf", "printed.vcf", spark_master="local[1]")
print(r["status"], r["summary"]["records"])
```

- **示例输出：**

```text
completed 7
```
