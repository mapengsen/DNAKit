# GWAS Studies and Associations: parameters and notes

See the [tutorial](../../api/features/22_database_query_gwas.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve GWAS associations for a SNP {#dnakit.databases.gwas_get_associations_for_snp-parameters}

`dnakit.databases.gwas_get_associations_for_snp(*, rs_id, sort=None, direction=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rs_id` | string | Required | dbSNP rs identifier |
| `sort` | string | None | Sort field (e.g., &#x27;p_value&#x27;, &#x27;or_value&#x27;) |
| `direction` | string | None | Sort direction (&#x27;asc&#x27; or &#x27;desc&#x27;) |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | Progress flag. |


## Retrieve associations from a GWAS study {#dnakit.databases.gwas_get_associations_for_study-parameters}

`dnakit.databases.gwas_get_associations_for_study(*, accession_id, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession_id` | string | Required | Study accession identifier |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS associations for a trait {#dnakit.databases.gwas_get_associations_for_trait-parameters}

`dnakit.databases.gwas_get_associations_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `disease_trait` | string | None | Disease or trait name for text-based search (e.g., &#x27;diabetes&#x27;, &#x27;hypertension&#x27;) |
| `efo_uri` | string | None | Full EFO ontology URI (e.g., &#x27;http://www.ebi.ac.uk/efo/EFO_0001645&#x27;) |
| `efo_id` | string | None | EFO/OBA term ID (e.g., &#x27;EFO_0001645&#x27;, &#x27;OBA_2050062&#x27;). Recommended for reliable trait filtering. |
| `efo_trait` | string | None | Exact EFO trait label. Use when you know the canonical trait string. |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | Progress flag. |

Provide exactly one of these parameter groups: `disease_trait`; `efo_uri`; `efo_id`; `efo_trait`.


## Retrieve GWAS SNP details by rsID {#dnakit.databases.gwas_get_snp_by_id-parameters}

`dnakit.databases.gwas_get_snp_by_id(*, rs_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rs_id` | string | Required | dbSNP rs identifier |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS SNPs mapped to a gene {#dnakit.databases.gwas_get_snps_for_gene-parameters}

`dnakit.databases.gwas_get_snps_for_gene(*, size=None, page=None, gene_symbol=None, gene=None, mapped_gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `size` | integer | None | Number of results to return per page |
| `page` | integer | None | Page number for pagination |
| `gene_symbol` | string | None | Gene name or symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;GBA&#x27;) |
| `gene` | string | None | Alias for gene_symbol |
| `mapped_gene` | string | None | Alias for gene_symbol |
| `verbose` | bool | True | Progress flag. |


## Retrieve GWAS studies for a trait {#dnakit.databases.gwas_get_studies_for_trait-parameters}

`dnakit.databases.gwas_get_studies_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
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
| `verbose` | bool | True | Progress flag. |

Provide exactly one of these parameter groups: `disease_trait`; `efo_uri`; `efo_id`; `efo_trait`.


## Retrieve a GWAS study by accession {#dnakit.databases.gwas_get_study_by_id-parameters}

`dnakit.databases.gwas_get_study_by_id(*, study_id=None, accession_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `study_id` | string | None | GWAS study identifier (GCST accession, e.g. &#x27;GCST000392&#x27;) |
| `accession_id` | string | None | Alias for study_id -- the same GCST accession identifier, accepted so results from other gwas_* tools (which use the field name &#x27;accession_id&#x27;) can be passed straight in. |
| `verbose` | bool | True | Progress flag. |


## Find GWAS variants associated with a disease or trait {#dnakit.databases.gwas_get_variants_for_trait-parameters}

`dnakit.databases.gwas_get_variants_for_trait(*, disease_trait=None, trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, limit=None, page=None, sort=None, direction=None, p_value_threshold=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
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
| `p_value_threshold` | number | None | P-value filter applied to the current page, not an exhaustive database scan. {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `verbose` | bool | True | Progress flag. |

Provide exactly one of these parameter groups: `disease_trait`; `trait`; `efo_uri`; `efo_id`; `efo_trait`.


## Search GWAS associations by keyword {#dnakit.databases.gwas_search_associations-parameters}

`dnakit.databases.gwas_search_associations(*, disease_trait=None, query=None, efo_uri=None, efo_id=None, efo_trait=None, rs_id=None, accession_id=None, p_value=None, p_value_threshold=None, sort=None, direction=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
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
| `verbose` | bool | True | Progress flag. |


## Search gene GWAS associations and risk alleles {#dnakit.databases.gwas_search_associations_by_gene-parameters}

`dnakit.databases.gwas_search_associations_by_gene(*, gene_name, size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_name` | string | Required | Gene symbol (e.g., BRCA1). |
| `size` | integer | None | Max associations to return (default 100; results are sorted by p-value, most significant first). {&quot;default&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## Search GWAS SNPs by rsID or gene {#dnakit.databases.gwas_search_snps-parameters}

`dnakit.databases.gwas_search_snps(*, rs_id=None, mapped_gene=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rs_id` | string | None | dbSNP rs identifier |
| `mapped_gene` | string | None | Gene name or symbol |
| `size` | integer | None | Number of results to return |
| `page` | integer | None | Page number for pagination |
| `verbose` | bool | True | Progress flag. |


## Search GWAS studies by trait and cohort {#dnakit.databases.gwas_search_studies-parameters}

`dnakit.databases.gwas_search_studies(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
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
| `verbose` | bool | True | Progress flag. |
