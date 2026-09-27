# ClinVar Clinical Significance: parameters and notes

See the [tutorial](../../api/features/22_database_query_clinical_variants.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve ClinVar clinical significance {#dnakit.databases.clin_var_get_clinical_significance-parameters}

`dnakit.databases.clin_var_get_clinical_significance(*, variant_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | ClinVar variant ID (e.g., &#x27;12345&#x27;, &#x27;123456&#x27;) |
| `verbose` | bool | True | Progress flag. |


## Retrieve ClinVar variant summaries and review status {#dnakit.databases.clin_var_get_variant_details-parameters}

`dnakit.databases.clin_var_get_variant_details(*, variant_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | ClinVar variant ID (e.g., &#x27;12345&#x27;, &#x27;123456&#x27;) |
| `verbose` | bool | True | Progress flag. |


## Search ClinVar variants by gene, disease or identifier {#dnakit.databases.clin_var_search_variants-parameters}

`dnakit.databases.clin_var_search_variants(*, gene=None, condition=None, variant_id=None, variant_name=None, max_results=None, limit=None, clinical_significance=None, gene_symbol=None, significance=None, query=None, variant=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Gene name or symbol (e.g., &#x27;BRCA1&#x27;, &#x27;BRCA2&#x27;) At least one of gene, condition, or variant_id must be provided. |
| `condition` | string | None | Disease or condition name (e.g., &#x27;breast cancer&#x27;, &#x27;diabetes&#x27;) At least one of gene, condition, or variant_id must be provided. Matched against ClinVar&#x27;s disease index, which carries MedGen&#x27;s own condition names exactly -- a synonym, a locus name or a misspelling (&#x27;Lebers congenital amaurosis&#x27; for the indexed &#x27;Leber congenital amaurosis&#x27;) matches nothing there, so the search falls back to an unrestricted free-text match. When that happens the response says so: &#x27;condition_filter_applied&#x27;: false plus a &#x27;condition_filter_warning&#x27;, meaning total_count is NOT a disease-filtered count. Use MedGen_search_conditions to get the indexed name. |
| `variant_id` | string | None | ClinVar variant ID (e.g., &#x27;12345&#x27;) At least one of gene, condition, or variant_id must be provided. |
| `variant_name` | string &#124; array &#124; null | None | Protein change or HGVS notation to search for (e.g., &#x27;Glu6Val&#x27;, &#x27;V600E&#x27;, &#x27;c.20A&gt;T&#x27;), matched against ClinVar&#x27;s own &#91;Variant name&#93; index -- combine with &#x27;gene&#x27; for an exact-match lookup of a specific known variant instead of browsing all of a gene&#x27;s ClinVar entries. A gene-level browse can silently miss well-known but old/low-ID records once a gene has more entries than max_results (e.g. HBB&#x27;s canonical sickle-cell record, NM_000518.5(HBB):c.20A&gt;T (p.Glu7Val), isn&#x27;t within the first 100-425 gene-level rows). Accepts a list to try multiple candidate names in one call (e.g. every allele at a multi-allelic site) -- combined with OR. |
| `max_results` | integer | None | Maximum number of results to return (default: 20). Alias: limit. {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `limit` | integer | None | Alias for max_results: maximum number of results to return. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `clinical_significance` | string | None | Filter by clinical significance. Must be one of the classes ClinVar&#x27;s index actually recognizes: &#x27;Pathogenic&#x27;, &#x27;Likely pathogenic&#x27;, &#x27;Uncertain significance&#x27;, &#x27;VUS&#x27;, &#x27;Likely benign&#x27;, &#x27;Benign&#x27;, &#x27;Established risk allele&#x27;, &#x27;Likely risk allele&#x27;, &#x27;Uncertain risk allele&#x27;, &#x27;drug response&#x27;, &#x27;risk factor&#x27;, &#x27;not provided&#x27;, &#x27;other&#x27;. Case-insensitive; underscores/hyphens are treated as spaces. Combine two classes with &#x27;/&#x27; to search their union (e.g. &#x27;Pathogenic/Likely pathogenic&#x27;). Any other value is rejected with an error rather than silently filtering the results to nothing. |
| `gene_symbol` | string &#124; null | None | Alias for gene. HGNC gene symbol (e.g., &quot;DPYD&quot;, &quot;CYP2C19&quot;). |
| `significance` | string &#124; null | None | Alias for clinical_significance (e.g., &quot;pathogenic&quot;, &quot;benign&quot;, &quot;uncertain_significance&quot;). Same accepted value set and same validation as clinical_significance. |
| `query` | string &#124; null | None | Alias for condition. Free-text search mapped to condition/disease field. |
| `variant` | string &#124; array &#124; null | None | DNAKit optional alias for variant_name |
| `verbose` | bool | True | Progress flag. |
