# 队列表型关联与精细定位查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_phewas.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询日本生物样本库中的变异表型关联 {#dnakit.databases.bio_bank_japan_phewas_by_variant-parameters}

`dnakit.databases.bio_bank_japan_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh37 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh37, e.g. &#x27;10:114758349:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 FinnGen 变异精细定位结果 {#dnakit.databases.finn_gen_get_variant_finemapping-parameters}

`dnakit.databases.finn_gen_get_variant_finemapping(*, variant, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | Genomic variant in chr:pos:ref:alt format (GRCh38). Examples: &#x27;19:44908684:T:C&#x27; (APOE rs429358), &#x27;10:112998590:C:T&#x27; (TCF7L2 rs7903146). Also accepts chr-pos-ref-alt format. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因稀有变异负荷与表型的关联 {#dnakit.databases.genebass_gene_burden_phewas-parameters}

`dnakit.databases.genebass_gene_burden_phewas(*, gene, burden_set=None, limit=None, max_pval=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填 | Ensembl gene ID (e.g. &#x27;ENSG00000148737&#x27;) or gene symbol (e.g. &#x27;TCF7L2&#x27;). Symbols are auto-resolved via Ensembl. |
| `burden_set` | string | None | Variant category aggregated in the burden test. One of: &#x27;pLoF&#x27; (default), &#x27;missense&#124;LC&#x27;, &#x27;synonymous&#x27;. Aliases accepted: lof/plof -&gt; pLoF, missense -&gt; missense&#124;LC. {&quot;default&quot;: &quot;pLoF&quot;} |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending burden p-value). {&quot;default&quot;: 25} |
| `max_pval` | number &#124; null | None | Optional: only return associations with burden p-value &lt;= this threshold (e.g. 2.5e-6 for exome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询台湾精准医疗队列的变异表型关联 {#dnakit.databases.tpmi_phewas_by_variant-parameters}

`dnakit.databases.tpmi_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 UKB-TOPMed 的变异表型关联 {#dnakit.databases.ukbtopmed_phewas_by_variant-parameters}

`dnakit.databases.ukbtopmed_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | 运行提示开关。 |
