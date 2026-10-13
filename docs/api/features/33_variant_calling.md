# 变异检测与统计

## 1) 调用 FreeBayes 进行单倍型变异检测

- **作用：** 将调用者准备好的 BAM/CRAM 和参考 FASTA 交给 FreeBayes，输出 VCF。
- **计算方法：** 固定调用 `freebayes -f reference -p ploidy ...`，可选最小替代等位基因比例和计数；不复制 FreeBayes 源码。
- **API：** `dnakit.variants.call_freebayes(reference, alignments, output_directory, ploidy=2, ...)`
- **输入/输出：** 已准备的比对文件、参考 FASTA 和新输出目录；返回 JSON 可序列化阶段信封及 `variants.vcf`。
- **引用与许可：** Garrison & Marth (2012), [arXiv:1207.3907](https://arxiv.org/abs/1207.3907)；官方仓库记录 MIT。FreeBayes 仍需由调用者单独安装。

## 2) 调用 VCFtools 过滤并重新编码 VCF

- **作用：** 对本地 VCF 应用固定的去 indel、去 SNP、缺失率和最小等位基因频率过滤。
- **API：** `dnakit.variants.filter_vcf_vcftools(source, output_directory, remove_indels=False, remove_snps=False, max_missing=None, min_allele_frequency=None, ...)`
- **输出：** VCFtools 的 `filtered.recode.vcf` 和完整诊断产物；输入文件不会覆盖。
- **引用与许可：** Danecek et al. (2011), DOI `10.1093/bioinformatics/btr330`；官方仓库记录 LGPL-3.0。完整 VCFtools 参数面没有全部暴露。

## 3) 使用 vcfanno 注释 VCF

- **作用：** 按调用者提供的 TOML 配置，从 VCF、BED 或其他索引注释文件向查询 VCF 的 INFO 字段写入注释。
- **API：** `dnakit.variants.annotate_vcf_vcfanno(source, configuration, output_directory, threads=1, ...)`。
- **输出：** `annotated.vcf`、stderr 和结构化运行信封；`configuration` 与注释数据由调用者准备，DNAKit 不下载数据。
- **引用与许可：** Pedersen et al. (2016), [10.1186/s13059-016-0973-5](https://doi.org/10.1186/s13059-016-0973-5)；官方仓库 MIT。完整 TOML/Lua 参数面没有全部重导出。

## 4) 使用 SURVIVOR 合并结构变异 VCF

- **作用：** 用 SURVIVOR 的 `merge` 模式，按断点距离、支持样本数、类型、链方向和最小 SV 长度合并至少两个 SV VCF。
- **API：** `dnakit.variants.merge_structural_variants_survivor(sources, output_directory, max_distance=100, min_support=1, min_sv_size=30, ...)`。
- **输出：** `merged.vcf`、输入清单和诊断日志；输出目录必须是新目录。
- **引用与许可：** Jeffares et al. (2017), [10.1038/NCOMMS14061](https://doi.org/10.1038/NCOMMS14061)；官方仓库 MIT。SURVIVOR 的模拟、格式转换和其他过滤子命令不在此接口中。

## 5) 使用 GEMINI 建库并查询变异

- **作用：** 将已经准备好的 VCF 装载为 GEMINI 本地变异数据库，并对该数据库执行有界 SQL 查询。
- **计算方法：** 分别调用 GEMINI 的 `load` 和 `query` 子命令；注释数据库、注释资源和下载由调用者准备，DNAKit 不隐式下载或伪造注释结果。
- **API：** `dnakit.variants.load_gemini_database(source_vcf, output_directory, annotation_type=None, cores=1, ...)`；`dnakit.variants.query_gemini_database(database, query, output_directory, header=True, result_format=None, ...)`。
- **输入/输出：** `load` 需要本地 VCF 和新输出目录，返回包含 `variants.db` 的 JSON 可序列化信封；`query` 需要本地 GEMINI 数据库和 SQL 文本，返回 `query.tsv`、查询参数和诊断日志。
- **示例代码：**

```python
from dnakit.variants import load_gemini_database, query_gemini_database

loaded = load_gemini_database("prepared.vcf", "gemini_database")
queried = query_gemini_database(
    loaded["result"]["database"],
    "SELECT chrom, start, ref, alt FROM variants LIMIT 10",
    "gemini_query",
)
print(queried["result"]["output"])
```

- **引用与许可：** GEMINI 官方仓库记录 MIT；请引用 Paila et al. (2013), [10.1371/journal.pcbi.1003153](https://doi.org/10.1371/journal.pcbi.1003153)。GEMINI、其注释资源和数据库构建环境仍需由调用者单独安装和准备。

## 6) 使用 Ensembl VEP 注释 VCF

- **作用：** 对调用者准备的 VCF 运行 Ensembl Variant Effect Predictor，生成带功能后果注释的 VCF、表格或 JSON。
- **计算方法：** 使用 VEP 的 `--input_file`/`--output_file` 边界，可选物种、组装版本、离线缓存、参考 FASTA、输出格式和受限并行；缓存、转录本资源、插件与数据库连接由调用者准备。
- **API：** `dnakit.variants.annotate_vcf_ensembl_vep(source_vcf, output_directory, species="homo_sapiens", assembly=None, cache_directory=None, fasta=None, offline=True, output_format="vcf", fork=1, ...)`。
- **输入/输出：** 需要本地 VCF 和新输出目录；默认返回包含 `annotated.vcf` 的 JSON 可序列化信封，也可选择 `annotated.tsv` 或 `annotated.json`。
- **示例代码：**

```python
from dnakit.variants import annotate_vcf_ensembl_vep

result = annotate_vcf_ensembl_vep(
    "variants.vcf",
    "vep_output",
    assembly="GRCh38",
    cache_directory="/data/vep_cache",
    offline=True,
)
print(result["result"]["output"])
```

- **引用与许可：** Ensembl VEP 官方仓库记录 Apache-2.0；请引用 McLaren et al. (2016), [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)。DNAKit 不下载缓存、参考序列或插件。

## 1) 从 BAM 检测变异 {#bio-17}

- **作用：** 输出规范化、可索引查询的 VCF。
- **计算方法：** bcftools 依次运行 mpileup、call 和 norm。
- **API：** `dnakit.bioinformatics.ngs.call_variants(bam, reference, output_directory, ploidy=2, config=None)`
- **输入：** 已排序并索引的 BAM、带 `.fai` 的参考 FASTA、明确的倍性；bcftools 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **Agent/MCP API：** `dnakit.ngs.call_variants(bam, reference, output_directory, ploidy=2, ...)`；返回 JSON 可序列化的 VCF、索引和诊断信息。
- **引用与许可：** BCFtools 使用 Danecek et al. (2021), DOI `10.1093/gigascience/giab008`（MIT）；后端由调用者安装。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import call_variants

result = call_variants("aligned.bam", "reference.fa", "new_variants", ploidy=2)
print(result.artifacts["variants.vcf.gz"])
```

- **示例结果：**

返回 `ToolResult`；`artifacts['variants.vcf.gz']` 是压缩 VCF 路径，同时生成索引。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

<!-- BEGIN GENERATED API: vcf -->

<span id="vcf"></span>

## 2) 按过滤条件统计 VCF 变异数量 {#dnakit.databases.vcf_count_variants}

- **作用：** 按过滤条件统计 VCF 变异数量。
- **处理方法：** 先用 bcftools view 过滤，再用 bcftools stats 统计变异记录。
- **API：** `dnakit.databases.vcf_count_variants(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`
- **输入：** 必须提供 `vcf_path`。 需本地 VCF/BCF 文件，且 PATH 中有 `bcftools`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 返回字典，包含 `status`, `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.vcf_count_variants(
    operation='count_variants',
    vcf_path='/tmp/sample.vcf.gz',
    min_qual=30,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "success",
  "data": "<variant counts and filters_applied>",
  "metadata": "<dict>"
}
```

[完整参数与使用说明](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_count_variants-parameters)

## 3) 规范化 VCF 变异表示与多等位位点 {#dnakit.databases.vcf_normalize}

- **作用：** 规范化 VCF 变异表示与多等位位点。
- **处理方法：** 用 bcftools norm 拆分或合并多等位位点；提供参考时左对齐，并比较前后统计。
- **API：** `dnakit.databases.vcf_normalize(*, operation=None, vcf_path, multiallelics=None, reference_fasta=None, verbose=True)`
- **输入：** 必须提供 `vcf_path`。 需本地 VCF/BCF 文件，且 PATH 中有 `bcftools`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 返回字典，包含 `status`, `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.vcf_normalize(
    operation='normalize',
    vcf_path='/tmp/sample.vcf.gz',
    multiallelics='split',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "success",
  "data": {
    "before": "<variant counts>",
    "after": "<variant counts>",
    "left_aligned": "<bool>",
    "multiallelics": "<str>"
  },
  "metadata": "<dict>"
}
```

[完整参数与使用说明](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_normalize-parameters)

## 4) 汇总 VCF 变异类型与样本基因型统计 {#dnakit.databases.vcf_summary_stats}

- **作用：** 汇总 VCF 变异类型与样本基因型统计。
- **处理方法：** 对筛选后的 VCF 调用 bcftools stats，解析变异类型和样本统计。
- **API：** `dnakit.databases.vcf_summary_stats(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`
- **输入：** 必须提供 `vcf_path`。 需本地 VCF/BCF 文件，且 PATH 中有 `bcftools`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 返回字典，包含 `status`, `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.vcf_summary_stats(
    operation='summary_stats',
    vcf_path='/tmp/sample.vcf.gz',
    pass_only=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "success",
  "data": "<variant and sample statistics; filters_applied>",
  "metadata": "<dict>"
}
```

[完整参数与使用说明](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_summary_stats-parameters)

<!-- END GENERATED API: vcf -->

<!-- BEGIN GENERATED API: variant-statistics -->

<span id="variant-statistics"></span>

## 5) 计算满足注释与频率条件的编码变异比例 {#dnakit.databases.coding_variant_fraction}

- **作用：** 计算满足注释与频率条件的编码变异比例。
- **处理方法：** 以编码区变异为分母，统计满足注释及 VAF 阈值的比例。
- **API：** `dnakit.databases.coding_variant_fraction(*, file, vaf_threshold=None, annotation=None, header_rows=None, vaf_col=None, so_col=None, verbose=True)`
- **输入：** 必须提供 `file`。
- **输出：** 返回字典，包含 `total_variants`, `vaf_filtered`, `coding_variants_below_vaf`, `matching_annotation` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.coding_variant_fraction(
    file='/tmp/tu_variants_demo.csv',
    vaf_threshold=0.3,
    annotation='synonymous_variant',
    header_rows=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_variants": "<integer>",
  "vaf_filtered": "<integer>",
  "coding_variants_below_vaf": "<integer>",
  "matching_annotation": "<integer>",
  "fraction": "<number>",
  "naive_fraction": "<number>",
  "vaf_column": "<string>",
  "so_column": "<string>"
}
```

[完整参数与使用说明](../../faq/api/33_variant_calling.md#dnakit.databases.coding_variant_fraction-parameters)

<!-- END GENERATED API: variant-statistics -->
