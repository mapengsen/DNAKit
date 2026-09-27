# 变异检测与统计

## 1) 从 BAM 检测变异 {#bio-17}

- **作用：** 输出规范化、可索引查询的 VCF。
- **计算方法：** bcftools 依次运行 mpileup、call 和 norm。
- **API：** `dnakit.bioinformatics.ngs.call_variants(bam, reference, output_directory, ploidy=2, config=None)`
- **输入：** 已排序并索引的 BAM、带 `.fai` 的参考 FASTA、明确的倍性；bcftools 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
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
