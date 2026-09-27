# Variant Records and Population Frequencies: parameters and notes

See the [tutorial](../../api/features/22_database_query_variant_frequency.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve dbSNP allele frequencies {#dnakit.databases.dbsnp_get_frequencies-parameters}

`dnakit.databases.dbsnp_get_frequencies(*, rsid, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rsid` | string | Required | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | Progress flag. |


## Retrieve dbSNP coordinates and alleles by rsID {#dnakit.databases.dbsnp_get_variant_by_rsid-parameters}

`dnakit.databases.dbsnp_get_variant_by_rsid(*, rsid, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `rsid` | string | Required | dbSNP rsID (e.g., &#x27;rs12345&#x27;, &#x27;12345&#x27;) |
| `verbose` | bool | True | Progress flag. |


## Retrieve gnomAD gene constraint metrics {#dnakit.databases.gnomad_get_gene_constraints-parameters}

`dnakit.databases.gnomad_get_gene_constraints(*, gene_symbol, reference_genome=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `gene_symbol` | string | Required | Gene symbol (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;) |
| `reference_genome` | string | None | Reference genome. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve gnomAD variant metadata {#dnakit.databases.gnomad_get_variant-parameters}

`dnakit.databases.gnomad_get_variant(*, variant_id, dataset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | Variant ID (e.g., &#x27;19-44908822-C-T&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `verbose` | bool | True | Progress flag. |


## Retrieve ancestry-specific gnomAD allele frequencies {#dnakit.databases.gnomad_get_variant_populations-parameters}

`dnakit.databases.gnomad_get_variant_populations(*, variant_id, dataset=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | Variant ID (format &#x27;chrom-pos-ref-alt&#x27;, e.g., &#x27;1-55051215-G-GA&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r4&quot;} |
| `verbose` | bool | True | Progress flag. |


## Search gnomAD variant identifiers {#dnakit.databases.gnomad_search_variants-parameters}

`dnakit.databases.gnomad_search_variants(*, query=None, dataset=None, variant=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Variant search query (e.g., &#x27;rs7412&#x27;). |
| `dataset` | string | None | gnomAD dataset ID. Allowed values: gnomad_r4, gnomad_r4_non_ukb, gnomad_r3, gnomad_r3_controls_and_biobanks, gnomad_r3_non_cancer, gnomad_r3_non_neuro, gnomad_r3_non_topmed, gnomad_r3_non_v2, gnomad_r2_1, gnomad_r2_1_controls, gnomad_r2_1_non_neuro, gnomad_r2_1_non_cancer, gnomad_r2_1_non_topmed, exac. Coverage: gnomad_r4* = GRCh38 with genomes + exomes; gnomad_r3* = GRCh38 with genomes ONLY (no exome component, so `exome` is always null); gnomad_r2_1* = GRCh37 with genomes + exomes; exac = GRCh37 with exomes only. {&quot;default&quot;: &quot;gnomad_r3&quot;} |
| `variant` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | Progress flag. |
