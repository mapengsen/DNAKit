# Cancer Variants and Treatment Evidence: parameters and notes

See the [tutorial](../../api/features/22_database_query_cancer_variants.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve a CIViC molecular profile {#dnakit.databases.civic_get_molecular_profile-parameters}

`dnakit.databases.civic_get_molecular_profile(*, molecular_profile_id, evidence_limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `molecular_profile_id` | integer | Required | CIViC molecular profile ID (e.g., 12 for BRAF V600E) |
| `evidence_limit` | integer | None | Maximum number of evidence items to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |


## Retrieve a CIViC cancer variant {#dnakit.databases.civic_get_variant-parameters}

`dnakit.databases.civic_get_variant(*, variant_id, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | integer | Required | CIViC variant ID (e.g., 4170) |
| `verbose` | bool | True | Progress flag. |


## Find CIViC variants associated with a gene {#dnakit.databases.civic_get_variants_by_gene-parameters}

`dnakit.databases.civic_get_variants_by_gene(*, gene_id=None, gene_name=None, gene_symbol=None, gene=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_id` | integer | None | CIViC gene ID (e.g., 19 for EGFR, 12 for BRAF). Find gene IDs using civic_search_genes. |
| `gene_name` | string | None | Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;BRAF&#x27;, &#x27;TP53&#x27;). Will be looked up automatically. Aliases: gene, gene_symbol, query. |
| `gene_symbol` | string | None | Alias for gene_name. Standard gene symbol (e.g., &#x27;KRAS&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;). |
| `gene` | string | None | Alias for gene_name. Gene symbol (e.g., &#x27;EGFR&#x27;, &#x27;KRAS&#x27;). |
| `limit` | integer | None | Maximum number of variants to return (default: 500, uses cursor pagination to bypass CIViC&#x27;s 100/page server cap) {&quot;default&quot;: 500} |
| `verbose` | bool | True | Progress flag. |


## Search CIViC clinical assertions {#dnakit.databases.civic_search_assertions-parameters}

`dnakit.databases.civic_search_assertions(*, limit=None, therapy=None, therapy_name=None, disease=None, disease_name=None, variant_name=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of assertions to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `therapy` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;ponatinib&#x27;). Alias: therapy_name. |
| `therapy_name` | string | None | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | None | Filter by disease name (e.g., &#x27;leukemia&#x27;, &#x27;melanoma&#x27;). Alias: disease_name. |
| `disease_name` | string | None | Alias for disease. Filter by disease name. |
| `variant_name` | string | None | Filter by variant name (e.g., &#x27;V600E&#x27;, &#x27;T315I&#x27;). |
| `verbose` | bool | True | Progress flag. |


## Search CIViC clinical evidence items {#dnakit.databases.civic_search_evidence_items-parameters}

`dnakit.databases.civic_search_evidence_items(*, limit=None, status=None, therapy=None, therapy_name=None, disease=None, disease_name=None, molecular_profile=None, evidence_type=None, significance=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of evidence items to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `status` | string &#124; null | None | Filter by curation status. Default: ACCEPTED (peer-reviewed). Options: ACCEPTED, SUBMITTED, REJECTED, ALL (returns all statuses combined). |
| `therapy` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;pembrolizumab&#x27;). Alias: therapy_name. Note: matches any evidence item where the therapy appears, including combination regimens — results may include multi-drug combinations. |
| `therapy_name` | string | None | Alias for therapy. Filter by therapy/drug name. |
| `disease` | string | None | Filter by disease name (e.g., &#x27;leukemia&#x27;, &#x27;melanoma&#x27;, &#x27;lung cancer&#x27;). Alias: disease_name. Note: CIViC uses specific disease names (e.g., &#x27;Lung Non-small Cell Carcinoma&#x27;, not &#x27;NSCLC&#x27;); try partial names or multiple searches if results are empty. |
| `disease_name` | string | None | Alias for disease. Filter by disease name. |
| `molecular_profile` | string | None | Filter by molecular profile name (e.g., &#x27;BRAF V600E&#x27;, &#x27;EGFR T790M&#x27;, &#x27;KRAS G12C&#x27;). Uses substring matching — &#x27;FLT3 ITD&#x27; will also match &#x27;FLT3 ITD AND FLT3 D835Y&#x27;. For gene fusions, CIViC uses double-colon notation: &#x27;GENE::PARTNER Fusion&#x27; (e.g., &#x27;FGFR2::BICC1 Fusion&#x27;, &#x27;ALK::EML4 Fusion&#x27;). Use civic_search_molecular_profiles to discover exact profile names. Gene-prefixed clinical spellings (&#x27;EGFRvIII&#x27;, &#x27;BRAFV600E&#x27;) are also handled: CIViC names the profile &#x27;&lt;gene&gt; &lt;designation&gt;&#x27; (&#x27;EGFR VIII&#x27;), so that separated form is queried too and both result sets are merged, reported in normalization_note. |
| `evidence_type` | string &#124; null | None | Filter by evidence type. Values: PREDICTIVE (drug response), DIAGNOSTIC (disease diagnosis), PROGNOSTIC (patient outcomes), PREDISPOSING (disease risk), ONCOGENIC (variant pathogenicity), FUNCTIONAL (molecular function). |
| `significance` | string &#124; null | None | Filter by clinical significance. Values depend on evidence_type, e.g. SENSITIVITYRESPONSE, RESISTANCE (predictive); BETTER_OUTCOME, POOR_OUTCOME (prognostic); POSITIVE, NEGATIVE (diagnostic); ONCOGENIC, PROTECTIVE. Use to separate, e.g., resistance from sensitivity evidence for a variant. |
| `verbose` | bool | True | Progress flag. |


## Search genes in CIViC {#dnakit.databases.civic_search_genes-parameters}

`dnakit.databases.civic_search_genes(*, limit=None, name=None, query=None, gene_name=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of genes to return (default: 10, recommended max: 100) {&quot;default&quot;: 10} |
| `name` | string | None | Gene symbol to search for (e.g., &quot;EGFR&quot;, &quot;BRAF&quot;, &quot;BRCA1&quot;). Alias: use &#x27;query&#x27; or &#x27;gene_name&#x27; instead. |
| `query` | string | None | Gene symbol to search for (e.g., &quot;FLT3&quot;, &quot;KRAS&quot;, &quot;TP53&quot;). Alias for &#x27;name&#x27;. |
| `gene_name` | string | None | Gene symbol to search for. Alias for &#x27;name&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Search therapies in CIViC {#dnakit.databases.civic_search_therapies-parameters}

`dnakit.databases.civic_search_therapies(*, limit=None, name=None, query=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of therapies to return (default: 20, recommended max: 100) {&quot;default&quot;: 20} |
| `name` | string | None | Filter by therapy/drug name (e.g., &#x27;imatinib&#x27;, &#x27;ponatinib&#x27;, &#x27;pembrolizumab&#x27;). Alias: query. |
| `query` | string | None | Alias for name. Filter by therapy/drug name. |
| `verbose` | bool | True | Progress flag. |


## Retrieve COSMIC somatic mutations for a gene {#dnakit.databases.cosmic_get_mutations_by_gene-parameters}

`dnakit.databases.cosmic_get_mutations_by_gene(*, operation=None, gene=None, gene_name=None, max_results=None, genome_build=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_by_gene) {&quot;const&quot;: &quot;get_by_gene&quot;, &quot;default&quot;: &quot;get_by_gene&quot;} |
| `gene` | string | None | Gene symbol (e.g., BRAF, TP53, EGFR, KRAS, PIK3CA). Alias: gene_name also accepted. |
| `gene_name` | string | None | Alias for gene parameter. Gene symbol (e.g., FLT3, BRAF, TP53). |
| `max_results` | integer | None | Maximum number of mutations to return (default: 100, max: 500) {&quot;default&quot;: 100, &quot;minimum&quot;: 1, &quot;maximum&quot;: 500} |
| `genome_build` | integer | None | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 {&quot;default&quot;: 37, &quot;enum&quot;: &#91;37, 38&#93;} |
| `verbose` | bool | True | Progress flag. |


## Search COSMIC cancer mutations {#dnakit.databases.cosmic_search_mutations-parameters}

`dnakit.databases.cosmic_search_mutations(*, operation=None, terms=None, query=None, max_results=None, genome_build=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `terms` | string | None | Search query - gene name (e.g., BRAF), mutation (e.g., V600E), or mutation ID (e.g., COSM476). Aliases: query or gene also accepted. |
| `query` | string | None | Alias for terms. Search query - gene name, mutation, or COSMIC ID. |
| `max_results` | integer | None | Maximum number of results to return (default: 20, max: 500) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 500} |
| `genome_build` | integer | None | Genome build version: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 37 {&quot;default&quot;: 37, &quot;enum&quot;: &#91;37, 38&#93;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve OncoKB oncogenicity and treatment evidence {#dnakit.databases.onco_kb_annotate_variant-parameters}

`dnakit.databases.onco_kb_annotate_variant(*, operation=None, gene=None, gene_symbol=None, variant=None, tumor_type=None, alteration=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: annotate_variant) {&quot;const&quot;: &quot;annotate_variant&quot;, &quot;default&quot;: &quot;annotate_variant&quot;} |
| `gene` | string | None | Gene symbol (e.g., BRAF, EGFR, TP53, KRAS) |
| `gene_symbol` | string | None | Alias for gene. Gene symbol (e.g., BRAF, EGFR, TP53, KRAS). |
| `variant` | string | None | Variant notation - protein change (e.g., V600E, T790M, G12D) |
| `tumor_type` | string | None | Optional OncoTree tumor type code (e.g., MEL for melanoma, LUAD for lung adenocarcinoma) |
| `alteration` | string | None | Alias for variant. Protein change notation (e.g., &quot;V600E&quot;, &quot;T790M&quot;, &quot;G12D&quot;). |
| `verbose` | bool | True | Progress flag. |


## Search cancer samples with regional copy number changes {#dnakit.databases.progenetix_cnv_search-parameters}

`dnakit.databases.progenetix_cnv_search(*, reference_name, start, end, variant_type=None, filters=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `reference_name` | string | Required | RefSeq chromosome accession. Examples: &#x27;refseq:NC_000007.14&#x27; (chr7/GRCh38), &#x27;refseq:NC_000017.11&#x27; (chr17/GRCh38), &#x27;refseq:NC_000001.11&#x27; (chr1/GRCh38). |
| `start` | integer | Required | Start position (1-based, GRCh38). Example: 55019017 for EGFR start. |
| `end` | integer | Required | End position (1-based, GRCh38). Example: 55211628 for EGFR end. |
| `variant_type` | string | None | CNV type: &#x27;DUP&#x27; for amplification/duplication, &#x27;DEL&#x27; for deletion. Leave empty for any CNV type. {&quot;default&quot;: &quot;&quot;} |
| `filters` | string | None | Optional NCIt ontology code to filter by cancer type. Example: &#x27;NCIT:C4017&#x27; for breast cancer. {&quot;default&quot;: &quot;&quot;} |
| `limit` | integer | None | Maximum number of biosamples to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |


## List Progenetix cancer CNV cohorts {#dnakit.databases.progenetix_list_cohorts-parameters}

`dnakit.databases.progenetix_list_cohorts(*, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `limit` | integer | None | Maximum number of cohorts to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `verbose` | bool | True | Progress flag. |


## List Progenetix filtering terms {#dnakit.databases.progenetix_list_filtering_terms-parameters}

`dnakit.databases.progenetix_list_filtering_terms(*, prefixes=None, limit=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `prefixes` | string | None | Ontology prefix to filter terms by. Use &#x27;NCIT&#x27; for NCI Thesaurus cancer types (default), &#x27;EDAM&#x27; for operation types, &#x27;HP&#x27; for HPO phenotypes. {&quot;default&quot;: &quot;NCIT&quot;} |
| `limit` | integer | None | Maximum number of terms to return (default: 25). {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `verbose` | bool | True | Progress flag. |


## Search cancer biosamples by disease ontology {#dnakit.databases.progenetix_search_biosamples-parameters}

`dnakit.databases.progenetix_search_biosamples(*, filters, limit=None, skip=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `filters` | string | Required | NCIt ontology code to filter biosamples by cancer type. Examples: &#x27;NCIT:C3058&#x27; (Glioblastoma), &#x27;NCIT:C4017&#x27; (Breast Carcinoma), &#x27;NCIT:C9145&#x27; (Acute Myeloid Leukemia), &#x27;NCIT:C3224&#x27; (Melanoma), &#x27;NCIT:C4872&#x27; (Lung Adenocarcinoma), &#x27;NCIT:C3167&#x27; (Neuroblastoma). |
| `limit` | integer | None | Maximum number of biosamples to return (default: 10). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `skip` | integer | None | Number of results to skip for pagination (default: 0). {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |
