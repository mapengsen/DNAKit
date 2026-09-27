# 变异记录与群体频率查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_variant_frequency.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询 dbSNP 等位基因频率 {#dnakit.databases.dbsnp_get_frequencies-parameters}

`dnakit.databases.dbsnp_get_frequencies(*, rsid, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | 必填 | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 按 rsID 查询 dbSNP 变异坐标与等位基因 {#dnakit.databases.dbsnp_get_variant_by_rsid-parameters}

`dnakit.databases.dbsnp_get_variant_by_rsid(*, rsid, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | 必填 | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 gnomAD 基因变异约束指标 {#dnakit.databases.gnomad_get_gene_constraints-parameters}

`dnakit.databases.gnomad_get_gene_constraints(*, gene_symbol, reference_genome=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填 | Gene symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;) |
| `reference_genome` | string | None | Reference genome. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 gnomAD 变异基本信息 {#dnakit.databases.gnomad_get_variant-parameters}

`dnakit.databases.gnomad_get_variant(*, variant_id, dataset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | Variant ID (e.g., &#x27;19-44908822-C-T&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 gnomAD 变异的分群等位基因频率 {#dnakit.databases.gnomad_get_variant_populations-parameters}

`dnakit.databases.gnomad_get_variant_populations(*, variant_id, dataset=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | Variant ID (format &#x27;chrom-pos-ref-alt&#x27;, e.g., &#x27;1-55051215-G-GA&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r4&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 gnomAD 变异编号 {#dnakit.databases.gnomad_search_variants-parameters}

`dnakit.databases.gnomad_search_variants(*, query=None, dataset=None, variant=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Variant search query (e.g., &#x27;rs7412&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `variant` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | 运行提示开关。 |
