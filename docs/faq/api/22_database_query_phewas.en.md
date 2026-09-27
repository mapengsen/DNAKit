# Cohort Phenotype Associations and Fine-mapping: parameters and notes

See the [tutorial](../../api/features/22_database_query_phewas.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Find variant–phenotype associations in BioBank Japan {#dnakit.databases.bio_bank_japan_phewas_by_variant-parameters}

`dnakit.databases.bio_bank_japan_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh37 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh37, e.g. &#x27;10:114758349:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | Progress flag. |


## Retrieve FinnGen fine-mapping results for a variant {#dnakit.databases.finn_gen_get_variant_finemapping-parameters}

`dnakit.databases.finn_gen_get_variant_finemapping(*, variant, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant` | string | Required | Genomic variant in chr:pos:ref:alt format (GRCh38). Examples: &#x27;19:44908684:T:C&#x27; (APOE rs429358), &#x27;10:112998590:C:T&#x27; (TCF7L2 rs7903146). Also accepts chr-pos-ref-alt format. |
| `verbose` | bool | True | Progress flag. |


## Query gene burden associations across phenotypes {#dnakit.databases.genebass_gene_burden_phewas-parameters}

`dnakit.databases.genebass_gene_burden_phewas(*, gene, burden_set=None, limit=None, max_pval=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Ensembl gene ID (e.g. &#x27;ENSG00000148737&#x27;) or gene symbol (e.g. &#x27;TCF7L2&#x27;). Symbols are auto-resolved via Ensembl. |
| `burden_set` | string | None | Variant category aggregated in the burden test. One of: &#x27;pLoF&#x27; (default), &#x27;missense&#124;LC&#x27;, &#x27;synonymous&#x27;. Aliases accepted: lof/plof -&gt; pLoF, missense -&gt; missense&#124;LC. {&quot;default&quot;: &quot;pLoF&quot;} |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending burden p-value). {&quot;default&quot;: 25} |
| `max_pval` | number &#124; null | None | Optional: only return associations with burden p-value &lt;= this threshold (e.g. 2.5e-6 for exome-wide significance). |
| `verbose` | bool | True | Progress flag. |


## Query variant–phenotype associations in TPMI {#dnakit.databases.tpmi_phewas_by_variant-parameters}

`dnakit.databases.tpmi_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | Progress flag. |


## Query variant–phenotype associations in UKB-TOPMed {#dnakit.databases.ukbtopmed_phewas_by_variant-parameters}

`dnakit.databases.ukbtopmed_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rsid` | string | None | dbSNP rsID, e.g. &#x27;rs7903146&#x27;. Auto-resolved to GRCh38 coordinates via Ensembl. Provide this OR &#x27;variant&#x27;. |
| `variant` | string | None | Variant in chr:pos:ref:alt or chr-pos-ref-alt format, in GRCh38, e.g. &#x27;10:112998590:C:T&#x27;. Provide this OR &#x27;rsid&#x27;. A bare rsID is also accepted here. |
| `limit` | integer | None | Maximum number of phenotype associations to return (sorted by ascending p-value). {&quot;default&quot;: 20} |
| `max_pval` | number &#124; null | None | Optional: only return associations with p-value &lt;= this threshold (e.g. 5e-8 for genome-wide significance). |
| `verbose` | bool | True | Progress flag. |
