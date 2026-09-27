# Structural Variant Queries: parameters and notes

See the [tutorial](../../api/features/22_database_query_structural_variants.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve Ensembl structural variants in an interval {#dnakit.databases.ensembl_get_structural_variants-parameters}

`dnakit.databases.ensembl_get_structural_variants(*, species, region=None, feature=None, chrom=None, start=None, end=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `species` | string | Required | Species name (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;). |
| `region` | string | None | Genomic region in format &#x27;chr:start-end&#x27; (e.g., &#x27;17:43044295-43125370&#x27;). Maximum 5Mb. |
| `feature` | string | None | Feature type to retrieve. Fixed to &#x27;structural_variation&#x27; for this tool. {&quot;default&quot;: &quot;structural_variation&quot;} |
| `chrom` | string | None | Chromosome identifier; provide with start and end. |
| `start` | integer | None | Region start; this parameter branch uses Ensembl 1-based inclusive coordinates. {&quot;minimum&quot;: 1} |
| `end` | integer | None | Region end; this parameter branch uses Ensembl 1-based inclusive coordinates. {&quot;minimum&quot;: 1} |
| `verbose` | bool | True | Progress flag. |


## Retrieve Ensembl structural variant details {#dnakit.databases.ensembl_get_sv_detail-parameters}

`dnakit.databases.ensembl_get_sv_detail(*, species, id=None, sv_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `species` | string | Required | Species name (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;). |
| `id` | string | None | Structural variant accession from Ensembl/DGVa (e.g., &#x27;nsv2769779&#x27;, &#x27;esv3647486&#x27;). Use ensembl_get_structural_variants to find IDs. |
| `sv_id` | string | None | DNAKit optional alias for id |
| `verbose` | bool | True | Progress flag. |


## Find gnomAD structural variants by gene {#dnakit.databases.gnomad_get_sv_by_gene-parameters}

`dnakit.databases.gnomad_get_sv_by_gene(*, gene_symbol, reference_genome=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_symbol` | string | Required | Gene symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;, &#x27;MECP2&#x27;) |
| `reference_genome` | string | None | Reference genome assembly. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | Progress flag. |


## Find gnomAD structural variants by genomic interval {#dnakit.databases.gnomad_get_sv_by_region-parameters}

`dnakit.databases.gnomad_get_sv_by_region(*, chrom, start, stop=None, reference_genome=None, end=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chrom` | string | Required | Chromosome (e.g., &#x27;17&#x27;, &#x27;X&#x27;). |
| `start` | integer | Required | Start position (1-based, GRCh38 by default). |
| `stop` | integer | None | Stop position. |
| `reference_genome` | string | None | Reference genome assembly. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `end` | integer | None | DNAKit optional alias for stop |
| `verbose` | bool | True | Progress flag. |


## Retrieve gnomAD structural variant details {#dnakit.databases.gnomad_get_sv_detail-parameters}

`dnakit.databases.gnomad_get_sv_detail(*, variant_id=None, dataset=None, sv_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | None | gnomAD structural variant ID. GRCh38/v4 IDs look like &#x27;DEL_chr17_24e4872b&#x27;; GRCh37/v2.1 IDs look like &#x27;DEL_1_1937&#x27;. |
| `dataset` | string | None | Optional SV callset override. Defaults to the callset auto-detected from the variant ID format (gnomad_sv_r2_1 for GRCh37-style IDs, gnomad_sv_r4 otherwise). Set explicitly only to query a v2.1 subset such as gnomad_sv_r2_1_controls. {&quot;enum&quot;: &#91;&quot;gnomad_sv_r2_1&quot;, &quot;gnomad_sv_r2_1_controls&quot;, &quot;gnomad_sv_r2_1_non_neuro&quot;, &quot;gnomad_sv_r4&quot;&#93;} |
| `sv_id` | string | None | DNAKit optional alias for variant_id |
| `verbose` | bool | True | Progress flag. |
