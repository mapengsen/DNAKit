# Gene–Disease Evidence, Dosage and Actionability: parameters and notes

See the [tutorial](../../api/features/22_database_query_genetic_evidence.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve ClinGen gene dosage sensitivity curations {#dnakit.databases.clin_gen_dosage_by_gene-parameters}

`dnakit.databases.clin_gen_dosage_by_gene(*, gene=None, gene_symbol=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol to search for (e.g., &#x27;BRCA1&#x27;, &#x27;MECP2&#x27;, &#x27;TP53&#x27;). Case-insensitive. |
| `gene_symbol` | string | None | DNAKit optional alias for gene |
| `verbose` | bool | True | Progress flag. |


## Search ClinGen dosage sensitivity by genomic interval {#dnakit.databases.clin_gen_dosage_region_search-parameters}

`dnakit.databases.clin_gen_dosage_region_search(*, chromosome, start, end, assembly=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chromosome` | string | Required | Chromosome (e.g., &#x27;17&#x27;, &#x27;X&#x27;). |
| `start` | integer | Required | Start position. |
| `end` | integer | Required | End position. |
| `assembly` | string | None | Genome assembly version. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve adult clinical actionability curations {#dnakit.databases.clin_gen_get_actionability_adult-parameters}

`dnakit.databases.clin_gen_get_actionability_adult(*, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | Progress flag. |


## Retrieve pediatric clinical actionability curations {#dnakit.databases.clin_gen_get_actionability_pediatric-parameters}

`dnakit.databases.clin_gen_get_actionability_pediatric(*, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | Progress flag. |


## Retrieve ClinGen dosage sensitivity records {#dnakit.databases.clin_gen_get_dosage_sensitivity-parameters}

`dnakit.databases.clin_gen_get_dosage_sensitivity(*, gene=None, include_regions=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `include_regions` | boolean | None | Include recurrent CNV regions (default: false) {&quot;default&quot;: false} |
| `verbose` | bool | True | Progress flag. |


## Retrieve ClinGen gene–disease validity curations {#dnakit.databases.clin_gen_get_gene_validity-parameters}

`dnakit.databases.clin_gen_get_gene_validity(*, gene=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | Progress flag. |


## Retrieve ClinGen variant classifications {#dnakit.databases.clin_gen_get_variant_classifications-parameters}

`dnakit.databases.clin_gen_get_variant_classifications(*, gene=None, variant=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `variant` | string | None | Optional: filter to one variant. Accepts a ClinGen canonical allele ID (&#x27;CA16020993&#x27;, with or without a &#x27;CAR:&#x27; prefix), a ClinVar VariationID (&#x27;370074&#x27;), full HGVS (&#x27;NM_000277.3:c.1315+1G&gt;T&#x27;), or a protein change (&#x27;p.Arg408Trp&#x27;). The form supplied selects the matching server-side filter. HGVS and protein-change strings match on substring, so a short fragment can return many rows -- combine with &#x27;gene&#x27; to scope them. |
| `verbose` | bool | True | Progress flag. |


## Search gene actionability across adult and pediatric contexts {#dnakit.databases.clin_gen_search_actionability-parameters}

`dnakit.databases.clin_gen_search_actionability(*, gene, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol to search (e.g., BRCA1, MLH1, LDLR) |
| `verbose` | bool | True | Progress flag. |


## Search dosage sensitivity scores by gene symbol {#dnakit.databases.clin_gen_search_dosage_sensitivity-parameters}

`dnakit.databases.clin_gen_search_dosage_sensitivity(*, gene, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol to search (e.g., MECP2, PMP22, RAI1) |
| `verbose` | bool | True | Progress flag. |


## Search gene–disease validity by gene symbol {#dnakit.databases.clin_gen_search_gene_validity-parameters}

`dnakit.databases.clin_gen_search_gene_validity(*, gene, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene` | string | Required | Gene symbol to search (e.g., BRCA1, TP53, CFTR) |
| `verbose` | bool | True | Progress flag. |
