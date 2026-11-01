# Ensembl VEP：变异后果、过滤、单倍型与缓存工具

`dnakit.variants` 中的这组函数是对 [Ensembl VEP](https://github.com/Ensembl/ensembl-vep)（用 116.3 版测试）各程序的带类型适配器：`vep`、`filter_vep`、`haplo`、`convert_cache.pl` 和输入格式校验脚本。`vep_annotate_variants()` 把 `vep` 的注释选项作为参数（VCF、tab、JSON 输出格式也可由较早的 `annotate_vcf_ensembl_vep()` 生成）。`variant_recoder` 没有封装：它需要实时连接 Ensembl MySQL 服务器，而这一点无法验证。VEP 使用 Apache-2.0 许可；DNAKit 把它作为独立进程运行，不复制其代码。

这些函数位于模块 `dnakit.variants.ensembl_vep_p09_cli_functions`。用 `conda install -c bioconda ensembl-vep htslib`（或从 GitHub 仓库运行 `perl INSTALL.pl`）安装 VEP 并把程序放入 PATH，或传入 `vep_path`、`filter_vep_path`、`haplo_path`、`convert_cache_path`。程序缺失时在创建任何输出目录之前抛出带安装提示的 `BackendUnavailableError`。缓存、FASTA 和插件由调用者负责，不会下载任何东西。每个函数都写入新的 `output_directory`，校验输入（`ConfigurationError`），VEP 失败时抛出带日志末尾的 `BackendExecutionError`。可选参数为 `None` 或 `False` 时不传，因此沿用 VEP 的默认值。

## 示例数据

所有示例都从工作目录读取这几个小文件：`genome.fa` 是 3000 碱基的 contig，`genes.gff3` 是其上一个含两个外显子的编码基因，`variants.vcf` 有四个 SNV（两个在 CDS，一个在内含子，一个在下游），`phased.vcf` 有一个样本的三个定相 SNV，`input.txt` 是 VEP 默认输入格式的三行文件，其中有两处错误。示例用 pysam 把 GFF3 转成 bgzip 压缩并建立 tabix 索引的文件，并假设 `vep`、`filter_vep` 和 `haplo` 在 PATH 中。

**示例文件：** `genome.fa`

```text
>chr1
GTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTC
AATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGG
TGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGC
TTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGG
TTAACGGGGATGGCAAGTACATTTTTTCGTAGATGTGCCTTGCTAACGAAAGTATTAAAC
ACGTCCCTCACAATAGAATCATAGTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAA
TACTCAATCCTGCGGGTTCGGTGACCTAAAACCCATTGATTGTGTTACCCAGTTCGAGCG
CATAGGGAATTCAGGTCCACACATGGCTGGATCCCCATGATATTCAAGAACTATACATTA
AGTTGAACCTCCAGAACACAATGGCTCAAGTTCCGTCAGTATTTTCTGTTGGCTTGGACA
ACTCTCTCTATTTTATTGGCAATTTAAGGACTCCTCTGTGGATTCGAGAGACGGCAATTA
AAGTGATCCGGCACAGTTTCTGGTATCGAAGGTACTCCACAAGAGCTGCTGTGAGTTTTA
CTGACGGAAAGTGAGCAAAGGCTAACGTTATTCCGTGAGCACGGGACATCCATTCTTCGT
GAGCTACAGCTCGAGAATCAGCTTCTAACCAAGCGATGCAGAACCGGCTACTTTAAGCAT
TGATGAATGCGTCGTAAGTGATACTCGACGATTCTCATGCAACGAAGTTAACCTATAGTA
ACTTACATTTTACGCGCTAGCTTCGCTGGAACTAATATCCATGTCTCAGAACTAGCGGCC
GAGAATGGGTTCCGAATCCTAAACTCCGACATGAGTTAAGGTTGCATACTAGGTCTGATA
CTAAAAGCGGGGTCAGGAGTCCGTCCAGAATATAATATAGCTAAGAAGAAAATGGAAAAT
CTATATCAGGTGGAGGCAGCCAATAGAGTGTGGTGGAATCAGTTTTGTAATATTACGACA
TAACCCTCGGAGGTTGAGCCAAACACCCCATTCGGCGGCCAAAGCGGGTGTCCCCACACA
TCCGTCGTAATGTCTAGGGAACAGTTGGAGCTGTGCGCAATCGTGTTGGACTATTGACAT
ACCCCTATTCGACCAGTGCAGTAGACTATACCACTTTTGAATCATCAGCAGACCTAATAG
CTTCATCCCTTCTAGTCGACTTTCCGGACCGATGCGCACTAATGATCGAAGTGTGTCTTT
ACTGAATCAGAAGTCGGAGAAAATTCTGCTGTACGAGATGTACGTAGCGATGATGATACG
GGTGATTCTTACATGAGGTACATCGAAAAAGGGTCATTGCGTTTACGTGAATGGGATTTG
CCTGGCCTATGGCCTTAGTACCTCTAAGAGGGCAACTAGTCACGGCGTAAACAGACAAGG
GTCGGATCCTAGTCACTAGTAGATCAACGGCTAAGTGGGCGTCACCAGACCCTCGCCCAT
CTGGACTAGTAGACCGTGTGCCTCAATCGACGCGTGAGGACACGTTGATTTCTTATACGG
CGTCTCAGCTTTCTTTTCCGCCAATTACATCAAATTCAAGTGCCTTGAGATACCGGGCCT
ATTCGTGCGTAGACCGGCGCGATGAGACAATGCCGGCATTTAAGTAAAGGTCCGGGATTA
TTGTAGACTAACTCTGCCAAATTATGCGTGCCTTGAGAGTATACACTCCGCATCCTTCTC
GTATGATGGCAGGGAGGGCCTTGGTTCCAGAACTTATCGCTAATGCGCTGGGTGCCATGC
ATGGAGCGTCTCTTGACATCAGGAAGTACCGGCTCCTGTTGGGGAGTTCGAAAATAGCTT
ATGGGACATCAGGCCATACTACTAGAACCTTGACCTGTACTGCATATCGTTTTCGGTTAA
ACTCGGTAGGTAAAAACACTTGTCGCCAGGATTCGCCACTCAAGGTTTATAGCAAAGAAT
ACTACTGAAACTGGCGTCCCCCATTAGTCACATCCTGCGAGGCGGCTCTCGACGGCATAA
CGGGATCTCTCCGGCTAAGGTCGAATGAGTTTAGTCTGAAGGTGGCAGCGCACACGACGG
TATGTTTCACCGCCGTCTGTTCGGAGCTATAGCATCATTAGGCCTAGCTATGGCGCCCCC
GTTTCCTACCGCAGCTACTGCATGCACTCGTTACGACAAAGTTCCTTTATTGGCTACAGG
ATTCTATCGTTCCGGCCTAGTGCGGGGGCATCTGCAGGCCGTCGTCCGTCCTACGTAAGG
CTCATGCGTACGGTAACAAAAGTTGGTTAGACCTGACCACGCGTGATCCGCGCCTCCTTC
CCCGCGTCACGTTAGGTGAAACCCTGCCCCGCTTAGCACGTCTCATGCTAACCCTCTTTC
CTGAGGGTAAGGGAACGTGCAATCCCGAAAAAGGTGAATGAAAGGCTGGGACGCTCCTGG
TTTGGGGCGCTCTAGGTGTTGTGGCTTGAACGTTAGCACCTCGGCTTCCGAATTGATACT
AGGCATCCTATCGAACACTATTAGCCGCCCTTCTTCATATTACTGTGGGTACAGACAAAT
AGAGGAGGTGTACGGTGACTAGCGGGTATCTATCATTAGGAAAGGGTACCGGCATCAGGA
AGAAAGCACCGTCCAATGATCACCACGGCTCCTTGCTTGCCACCGAAGCCATAGACATAT
GAATAGCCCGCTGTGGCTCTCCGACGAAGTGCGCGGCAGATCAGTTTTCTCCCCCTGAGA
TCACCAGACGACGGAACAGAACTGCGAGGCTGATCCACTGTTTTGCGCCTCTGAAATAGA
TTAAAGGGAATATCCCAGTGCAGGCCGAGGGGGTCTGATAATATACAGTCTAAAGAGTTA
GAATATCATCGTAACAGCCCATAGAACAATTCCCGTATTTAAGCATAGCATACGGATAAT
```

**示例文件：** `genes.gff3`

```text
##gff-version 3
chr1	test	gene	501	1150	.	+	.	ID=gene1;Name=GENE1;biotype=protein_coding
chr1	test	mRNA	501	1150	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201;biotype=protein_coding
chr1	test	exon	501	650	.	+	.	ID=exon1;Parent=tx1
chr1	test	CDS	501	650	.	+	0	ID=cds1;Parent=tx1
chr1	test	exon	1001	1150	.	+	.	ID=exon2;Parent=tx1
chr1	test	CDS	1001	1150	.	+	0	ID=cds2;Parent=tx1
```

**示例文件：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=3000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	521	.	A	C	.	.	.
chr1	530	.	T	A	.	.	.
chr1	700	.	C	A	.	.	.
chr1	2000	.	T	A	.	.	.
```

**示例文件：** `phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=3000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	522	.	T	A	.	.	.	GT	0|1
chr1	531	.	G	A	.	.	.	GT	1|1
chr1	540	.	A	C	.	.	.	GT	0|1
```

**示例文件：** `input.txt`

```text
5 140532 140532 T/C +
5 140540 140540 A/G f
3 100 100 T +
```

## 1) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_annotate_variants()`

- **作用：** 用 Ensembl VEP 预测变异后果，注释选项以带类型的参数给出。
- **处理方法：** 用固定参数列表运行 `vep`。一个函数涵盖 `summary`、`variant_class`、`hgvs`/`hgvsg`/`spdi`、3' 端移位（`shift_hgvs`、`shift_3prime`、`shift_genomic`）、`check_ref`/`lookup_ref`、按样本的基因型（`individual`、`individual_zyg`）、已知变异（`check_existing`）及 `af`、`af_1kg`、`af_gnomade`、`af_gnomadg`、`max_af`、频率过滤（`check_frequency`、`freq_pop`、`freq_freq`、`freq_gt_lt`、`freq_filter`、`filter_common`）、ClinVar（`clin_sig_allele`、`clinvar_somatic_classification`）、`gene_phenotype`、`mirna`、`regulatory`（`cell_type`）、UCSC 保守性分数（`phylop`、`phastcons`、`ucsc_assembly`、`ucsc_data_root`）、自定义轨道（`custom_annotations`、`gff`、`gtf`、`bigwig`）、用 `bam` 校正转录本、`domains`、`numbers`、`nearest` 和 `stats` 报告。`output_format` 为 `vcf`、`tab` 或 `json`。转录本来自缓存（`cache_directory`）或经 bgzip 压缩并建立 tabix 索引的 `gff`/`gtf` 加参考 `fasta`。需要缓存、数据库或网络的选项（频率、ClinVar、UCSC 分数；不提供 Condel 标志，因为当前 VEP 忽略它）原样传给 VEP，只有相应资源存在时才有数据；这里测试了它们的命令行映射，但没有用真实缓存运行。VEP 会在 FASTA 文件旁写索引。
- **API：** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_annotate_variants(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['vcf', 'tab', 'json']" = 'vcf', input_format: "Literal['ensembl', 'vcf', 'hgvs', 'id', 'region', 'vep', 'guess'] | None" = None, species: 'str' = 'homo_sapiens', assembly: 'str | None' = None, offline: 'bool | None' = None, cache_directory: 'str | Path | None' = None, cache_version: 'int | None' = None, fasta: 'str | Path | None' = None, gff: 'str | Path | None' = None, gtf: 'str | Path | None' = None, bigwig: 'str | Path | None' = None, bam: 'str | Path | None' = None, custom_annotations: 'Sequence[Mapping[str, object]] | None' = None, summary: 'bool' = False, variant_class: 'bool' = False, hgvs: 'bool' = False, hgvsg: 'bool' = False, spdi: 'bool' = False, shift_hgvs: 'int | None' = None, shift_3prime: 'int | None' = None, shift_genomic: 'int | None' = None, check_ref: 'bool' = False, lookup_ref: 'bool' = False, individual: 'str | Sequence[str] | None' = None, individual_zyg: 'str | Sequence[str] | None' = None, check_existing: 'bool' = False, af: 'bool' = False, af_1kg: 'bool' = False, af_gnomade: 'bool' = False, af_gnomadg: 'bool' = False, max_af: 'bool' = False, check_frequency: 'bool' = False, freq_pop: 'str | None' = None, freq_freq: 'float | None' = None, freq_gt_lt: "Literal['gt', 'lt'] | None" = None, freq_filter: "Literal['exclude', 'include'] | None" = None, filter_common: 'bool' = False, clin_sig_allele: 'int | None' = None, clinvar_somatic_classification: 'bool' = False, gene_phenotype: 'bool' = False, mirna: 'bool' = False, regulatory: 'bool' = False, cell_type: 'str | None' = None, phylop: 'Sequence[str] | None' = None, phastcons: 'Sequence[str] | None' = None, ucsc_assembly: 'str | None' = None, ucsc_data_root: 'str | None' = None, domains: 'bool' = False, numbers: 'bool' = False, nearest: "Literal['transcript', 'gene', 'symbol'] | None" = None, stats: "Literal['html', 'text'] | None" = None, fork: 'int | None' = None, vep_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 变异文件（默认 VCF）、新的输出目录和选项。既没有 `cache_directory` 也没有 `gff`/`gtf` 时，VEP 离线运行并需要默认目录中的缓存；`offline=False` 读取 Ensembl 数据库。
- **输出：** `outputs`（`output` 为 `annotations.vcf`/`.tsv`/`.json`，可选 `stats`、`warnings`）、`columns`、`records`（前几行；VCF 行带解析后的 `CSQ` 条目）和 `record_count`，以及 `version` 和 `command`。
- **示例代码：**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_annotate_variants

work = tempfile.mkdtemp()
result = vep_annotate_variants("variants.vcf", work + "/vcf", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
for row in result["records"]:
    print(row["pos"], row["annotations"][0]["Consequence"])
table = vep_annotate_variants("variants.vcf", work + "/tab", output_format="tab", gff="genes.gff3.gz", fasta="genome.fa", variant_class=True, hgvs=True, hgvsg=True, spdi=True, numbers=True, verbose=False)
first = table["records"][0]
print(first["VARIANT_CLASS"], first["EXON"], first["HGVSc"], first["HGVSg"], first["SPDI"])
```

- **示例输出：**

```text
521 synonymous_variant
530 synonymous_variant
700 intron_variant
2000 downstream_gene_variant
SNV 1/2 tx1:c.21A>C chr1:g.521A>C chr1:520:A:C
```

## 2) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_filter_results()`

- **作用：** 按表达式过滤 VEP 注释后的变异记录。
- **处理方法：** 运行 `filter_vep`。`filters` 为形如 `Consequence is missense_variant and AF < 0.01` 的表达式，多个表达式是“或”关系。`list_fields` 运行 `--list`，`only_matched` 把 VCF 的 CSQ 字段改写为只含匹配的条目，`soft_filter` 只标记不丢弃，`count_only` 返回匹配数，`input_format`、`vcf_info_field`、`start`、`limit`、`gzipped` 对应 `--format`、`--vcf_info_field`、`--start`、`--limit`、`--gz`。不提供本体选项（需要数据库）。
- **API：** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_filter_results(input_path: 'str | Path', output_directory: 'str | Path', *, filters: 'Sequence[str] | None' = None, list_fields: 'bool' = False, only_matched: 'bool' = False, soft_filter: 'bool' = False, count_only: 'bool' = False, input_format: "Literal['vcf', 'tab'] | None" = None, vcf_info_field: 'str | None' = None, start: 'int | None' = None, limit: 'int | None' = None, gzipped: 'bool' = False, filter_vep_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** VEP 结果文件（VCF 或 tab）、新的输出目录和表达式（只允许字母、数字和 `_ .,:<>=!()"'+*/%-`）。
- **输出：** `outputs`、`record_count` 和 `records`（按列拆分的数据行），`list_fields` 时为 `fields`，`count_only` 时为 `count`。
- **示例代码：**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_annotate_variants, vep_filter_results

work = tempfile.mkdtemp()
annotated = vep_annotate_variants("variants.vcf", work + "/vcf", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
path = annotated["outputs"]["output"]
print("IMPACT" in vep_filter_results(path, work + "/fields", list_fields=True, verbose=False)["fields"])
kept = vep_filter_results(path, work + "/kept", filters=["Consequence is intron_variant"], verbose=False)
print(kept["record_count"], kept["records"][0][1])
print(vep_filter_results(path, work + "/count", filters=["IMPACT is LOW"], count_only=True, verbose=False)["count"])
```

- **示例输出：**

```text
True
1 700
2
```

## 3) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_predict_haplotypes()`

- **作用：** 对定相基因型中完整转录本单倍型的后果做预测。
- **处理方法：** 运行 `haplo --json`：重建每个样本单倍型的转录本序列，并报告整个单倍型的氨基酸变化。转录本来自 `cache_directory` 或 `gff`/`gtf` 加 `fasta`；`haplotype_frequencies`、`phased`、`chromosomes` 对应 `--haplotype_frequencies`、`--phased`、`--chr`。
- **API：** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_predict_haplotypes(input_path: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['vcf', 'ensembl', 'region', 'vep', 'guess'] | None" = None, species: 'str' = 'homo_sapiens', assembly: 'str | None' = None, offline: 'bool | None' = None, cache_directory: 'str | Path | None' = None, cache_version: 'int | None' = None, fasta: 'str | Path | None' = None, gff: 'str | Path | None' = None, gtf: 'str | Path | None' = None, haplotype_frequencies: 'str | Path | None' = None, phased: 'bool' = False, chromosomes: 'str | None' = None, haplo_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 带基因型的 VCF 和新的输出目录。
- **输出：** `outputs`（`output` 为 `haplotypes.json`）、`records`（每个转录本一个对象，含 `protein_haplotypes` 和 `cds_haplotypes`）和 `record_count`。
- **示例代码：**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_predict_haplotypes

result = vep_predict_haplotypes("phased.vcf", tempfile.mkdtemp() + "/run", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
transcript = result["records"][0]
print(result["record_count"], transcript["transcript_id"])
for haplotype in transcript["protein_haplotypes"]:
    print(haplotype["name"], haplotype["frequency"])
```

- **示例输出：**

```text
1 tx1
tx1.p:8F>I,11G>S,14N>H 0.5
tx1.p:11G>S 0.5
```

## 4) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_validate_input_format()`

- **作用：** 检查 VEP 默认输入格式的文件。
- **处理方法：** 用当前 Python 解释器运行 Ensembl 的 `vep_input_format_validator.py`。它检查列数、等位基因字符与 `REF/ALT` 形式、链、坐标与等位基因长度是否一致、插入以及变异的排序。`validator_path` 是 ensembl-vep 仓库 `validator/` 目录中的脚本（conda 包不含它）。不写任何文件。
- **API：** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_validate_input_format(input_path: 'str | Path', *, validator_path: 'str | Path | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 输入文件和 `validator_path`。
- **输出：** `valid`、`error_count` 和 `errors`（行号、消息、行内容）。
- **示例代码：**

```python
import os
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_validate_input_format

script = os.environ.get("VEP_INPUT_VALIDATOR", "vep_input_format_validator.py")
result = vep_validate_input_format("input.txt", validator_path=script, verbose=False)
print(result["valid"], result["error_count"])
for error in result["errors"]:
    print(error["line"], error["message"])
```

- **示例输出：**

```text
False 2
2 Strand is not in the expected format (+ or -)
3 Allele 'T' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants
```

## 5) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_convert_cache()`

- **作用：** 转换下载的 VEP 缓存，使其通过 tabix 读取。
- **处理方法：** 对 `<cache_directory>/<species>/<version>` 运行 `convert_cache.pl`（conda 包中名为 `vep_convert_cache`）；`bgzip` 和 `tabix` 必须在 PATH 中。它会原地修改缓存；`remove_original` 为 `--remove`，`use_sereal` 为 `--sereal`，`compress_command` 为 `--compress`，`fork` 为 `--fork`。这里只检查了空缓存目录的运行，因为真实缓存有数 GB。
- **API：** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_convert_cache(cache_directory: 'str | Path', *, species: 'str', version: 'str', compress_command: 'str | None' = None, use_sereal: 'bool' = False, remove_original: 'bool' = False, fork: 'int | None' = None, convert_cache_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 包含 `<species>/<version>` 的目录、`species` 和 `version`（名称、逗号分隔的列表或 `all`）。
- **输出：** `status`、`command`、`cache_directory` 和程序的 `log_tail`。
- **示例代码：**

```python
import os
import tempfile
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_convert_cache

cache = tempfile.mkdtemp()
os.makedirs(cache + "/homo_sapiens/116_GRCh38/1")
result = vep_convert_cache(cache, species="homo_sapiens", version="116_GRCh38", verbose=False)
print(result["status"], "All done" in result["log_tail"])
```

- **示例输出：**

```text
completed True
```
