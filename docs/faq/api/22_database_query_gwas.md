# GWAS 研究与关联查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_gwas.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询指定 SNP 的 GWAS 关联 {#dnakit.databases.gwas_get_associations_for_snp-parameters}

`dnakit.databases.gwas_get_associations_for_snp(*, rs_id, sort=None, direction=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | 必填 | dbSNP rs identifier |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;) |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;) |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |


## 查询指定 GWAS 研究的关联结果 {#dnakit.databases.gwas_get_associations_for_study-parameters}

`dnakit.databases.gwas_get_associations_for_study(*, accession_id, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession_id` | string | 必填 | Study accession identifier |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |


## 查询指定性状的 GWAS 关联结果 {#dnakit.databases.gwas_get_associations_for_trait-parameters}

`dnakit.databases.gwas_get_associations_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;hypertension&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

以下参数组必须且只能提供一组：`disease_trait`; `efo_uri`; `efo_id`; `efo_trait`。


## 按 rsID 查询 GWAS SNP 详情 {#dnakit.databases.gwas_get_snp_by_id-parameters}

`dnakit.databases.gwas_get_snp_by_id(*, rs_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | 必填 | dbSNP rs identifier |
| `verbose` | bool | True | 运行提示开关。 |


## 查询映射到指定基因的 GWAS SNP {#dnakit.databases.gwas_get_snps_for_gene-parameters}

`dnakit.databases.gwas_get_snps_for_gene(*, size=None, page=None, gene_symbol=None, gene=None, mapped_gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `gene_symbol` | string | None | Gene name or symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;GBA&#x27;) |
| `gene` | string | None | Alias for gene_symbol |
| `mapped_gene` | string | None | Alias for gene_symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 查询指定性状的 GWAS 研究 {#dnakit.databases.gwas_get_studies_for_trait-parameters}

`dnakit.databases.gwas_get_studies_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;alzheimer disease&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | None | Cohort name (e.g., &#x27;UKB&#x27; for UK Biobank) |
| `gxe` | boolean | None | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | None | Filter for studies with full summary statistics |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |

以下参数组必须且只能提供一组：`disease_trait`; `efo_uri`; `efo_id`; `efo_trait`。


## 按研究编号查询 GWAS 详情 {#dnakit.databases.gwas_get_study_by_id-parameters}

`dnakit.databases.gwas_get_study_by_id(*, study_id=None, accession_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `study_id` | string | None | GWAS study identifier (GCST accession, e.g. &#x27;GCST000392&#x27;) |
| `accession_id` | string | None | Alias for study_id -- the same GCST accession identifier, accepted so results from other gwas_* tools (which use the field name &#x27;accession_id&#x27;) can be passed straight in. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索疾病或性状相关的 GWAS 变异 {#dnakit.databases.gwas_get_variants_for_trait-parameters}

`dnakit.databases.gwas_get_variants_for_trait(*, disease_trait=None, trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, limit=None, page=None, sort=None, direction=None, p_value_threshold=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `trait` | string | None | Alias for disease_trait. Disease or trait name (e.g., &#x27;type 2 diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | None | Number of results to return per page |
| `limit` | integer | None | Alias for size. Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;). Defaults to &#x27;p_value&#x27; so the most significant associations appear first. |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;). Defaults to &#x27;asc&#x27;. |
| `p_value_threshold` | number | None | 对当前返回页应用的 P 值过滤阈值；不代表已遍历全库。 {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |

以下参数组必须且只能提供一组：`disease_trait`; `trait`; `efo_uri`; `efo_id`; `efo_trait`。


## 按关键词检索 GWAS 关联 {#dnakit.databases.gwas_search_associations-parameters}

`dnakit.databases.gwas_search_associations(*, disease_trait=None, query=None, efo_uri=None, efo_id=None, efo_trait=None, rs_id=None, accession_id=None, p_value=None, p_value_threshold=None, sort=None, direction=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;coronary artery disease&#x27;, &#x27;breast cancer&#x27;). NOTE: GWAS associations are indexed by trait/disease, not by gene or drug name — passing a gene symbol or drug name here will fail. Use efo_id for precise filtering. |
| `query` | string | None | Alias for disease_trait. Disease or trait name (e.g., &#x27;type 2 diabetes&#x27;). NOTE: gene symbols and drug names are not valid GWAS trait queries. |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label (e.g., &#x27;serum alanine aminotransferase amount&#x27;). Use when you know the canonical trait string. |
| `rs_id` | string | None | dbSNP rs identifier |
| `accession_id` | string | None | Study accession identifier |
| `p_value` | number | None | Maximum p-value threshold for filtering (e.g., 5e-8 for genome-wide significance). IMPORTANT: the GWAS Catalog REST API v2 has NO server-side p-value filter (its OpenAPI document lists only pubmed_id, rs_id, full_pvalue_set, accession_id, efo_trait, efo_id, show_child_trait, mapped_gene, extended_geneset, sort, direction, page and size on GET /v2/associations), so this threshold is applied CLIENT-SIDE to the page of associations fetched. To make it meaningful the page is fetched p_value-ascending (most-significant-first) whenever sort is &#x27;p_value&#x27; (the default), and the surviving rows are then re-ordered into the direction you asked for. metadata.p_value_filter_scope, metadata.filtered_total and metadata.p_value_filter_note report exactly what was filtered. |
| `p_value_threshold` | number | None | Alias for p_value: maximum p-value threshold (also applied client-side; see p_value). |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;). Defaults to &#x27;p_value&#x27;. When a p_value threshold is supplied and sort is &#x27;p_value&#x27;, the upstream page is always fetched ascending so the threshold sees the significant end; metadata.sort_override records this and the returned rows are ordered per your &#x27;direction&#x27;. |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;). Defaults to &#x27;asc&#x27;. Combined with a p_value threshold this orders the RETURNED rows; see metadata.sort_override. |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因检索 GWAS 风险等位基因与显著性 {#dnakit.databases.gwas_search_associations_by_gene-parameters}

`dnakit.databases.gwas_search_associations_by_gene(*, gene_name, size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_name` | string | 必填 | Gene symbol (e.g., BRCA1). |
| `size` | integer | None | Max associations to return (default 100; results are sorted by p-value, most significant first). {&quot;default&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 rsID 或基因检索 GWAS SNP {#dnakit.databases.gwas_search_snps-parameters}

`dnakit.databases.gwas_search_snps(*, rs_id=None, mapped_gene=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `rs_id` | string | None | dbSNP rs identifier |
| `mapped_gene` | string | None | Gene name or symbol |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |


## 按性状与队列条件检索 GWAS 研究 {#dnakit.databases.gwas_search_studies-parameters}

`dnakit.databases.gwas_search_studies(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;breast cancer&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `cohort` | string | None | Cohort name (e.g., &#x27;UKB&#x27; for UK Biobank) |
| `gxe` | boolean | None | Filter for Gene-by-Environment interaction studies |
| `full_pvalue_set` | boolean | None | Filter for studies with full summary statistics |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | 运行提示开关。 |
