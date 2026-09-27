# Variant Records and Population Frequencies

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

## 1) dbSNP variant information {#7-dbsnp-variant-information}

- **Purpose:** Retrieve alleles and GRCh37/GRCh38 coordinates for one rsID.
- **Method:** Query dbSNP by rsID and retain alleles and assembly coordinates.
- **API:** `dnakit.annotation.get_dbsnp_variant(rsid[required])`
- **Input:** An rsID with or without the `rs` prefix.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_dbsnp_variant

result = get_dbsnp_variant("rs28934578")
print(result.data)
```

- **Example output:** dbSNP coordinates, assembly, and allele records.

## 2) dbSNP allele frequencies {#8-dbsnp-allele-frequencies}

- **Purpose:** Retrieve population allele-frequency records for one rsID.
- **Method:** Read dbSNP allele counts and frequencies by study or population.
- **API:** `dnakit.annotation.get_dbsnp_frequencies(rsid[required])`
- **Input:** One rsID.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_dbsnp_frequencies

result = get_dbsnp_frequencies("rs28934578")
print(result.data)
```

- **Example output:** Study/population allele-count and frequency records.

## 3) gnomAD variant search {#9-gnomad-variant-search}

- **Purpose:** Resolve an rsID or variant query to canonical `chrom-pos-ref-alt` identifiers.
- **Method:** Search the chosen gnomAD dataset for matching variant IDs.
- **API:** `dnakit.annotation.search_gnomad_variants(query[required], dataset[optional])`
- **Input:** Query text; DNAKit consistently defaults to `gnomad_r4`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import search_gnomad_variants

result = search_gnomad_variants("rs28934578")
print(result.data)
```

- **Example output:** Matching gnomAD variant IDs and basic metadata.

## 4) gnomAD aggregate variant data {#10-gnomad-aggregate-variant-data}

- **Purpose:** Retrieve aggregate frequency data and metadata for one gnomAD variant.
- **Method:** Retrieve gnomAD AC, AN and AF by variant ID.
- **API:** `dnakit.annotation.get_gnomad_variant(variant_id[required], dataset[optional])`
- **Input:** A `chrom-pos-ref-alt` identifier; dataset defaults to `gnomad_r4`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_gnomad_variant

result = get_gnomad_variant("17-7675088-C-G")
print(result.data)
```

- **Example output:** Aggregate genome/exome allele counts, allele numbers, and frequencies.

## 5) gnomAD population frequencies {#11-gnomad-population-frequencies}

- **Purpose:** Retrieve ancestry- and sex-stratified allele frequencies.
- **Method:** Extract ancestry/sex-stratified frequencies from the gnomAD variant record.
- **API:** `dnakit.annotation.get_gnomad_population_frequencies(variant_id[required], dataset[optional])`
- **Input:** A gnomAD variant ID; dataset defaults to `gnomad_r4`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_gnomad_population_frequencies

result = get_gnomad_population_frequencies("17-7675088-C-G")
print(result.data)
```

- **Example output:** AC, AN, and AF values stratified by ancestry and sex.

<!-- BEGIN GENERATED API: variant-records -->

<span id="variant-records"></span>

## 6) Retrieve dbSNP allele frequencies {#dnakit.databases.dbsnp_get_frequencies}

- **Purpose:** Retrieve dbSNP allele frequencies.
- **Method:** Use `dbSNPGetFrequencies` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dbsnp_get_frequencies(*, rsid, verbose=True)`
- **Input:** Provide `rsid`.
- **Output:** Returns a mapping with `status`, `data`, `rsid`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dbsnp_get_frequencies(
    rsid='rs429358',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "frequencies": [
      "<object>"
    ],
    "total_studies": "<integer>"
  },
  "rsid": "<string>",
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.dbsnp_get_frequencies-parameters)

## 7) Retrieve dbSNP coordinates and alleles by rsID {#dnakit.databases.dbsnp_get_variant_by_rsid}

- **Purpose:** Retrieve dbSNP coordinates and alleles by rsID.
- **Method:** Use `dbSNPGetVariantByRsID` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dbsnp_get_variant_by_rsid(*, rsid, verbose=True)`
- **Input:** Provide `rsid`.
- **Output:** Returns a mapping with `status`, `data`, `rsid`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dbsnp_get_variant_by_rsid(
    rsid='rs12345',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "refsnp_id": "<string>",
    "snp_id": "<integer>",
    "chromosome": "<string>",
    "position": "<string>",
    "assembly": "<string>",
    "position_grch38": "<string>",
    "position_grch37": "<string>",
    "allele": "<string>"
  },
  "rsid": "<string>",
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.dbsnp_get_variant_by_rsid-parameters)

## 8) Retrieve gnomAD gene constraint metrics {#dnakit.databases.gnomad_get_gene_constraints}

- **Purpose:** Retrieve gnomAD gene constraint metrics.
- **Method:** Use `gnomADGetGeneConstraints` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_gene_constraints(*, gene_symbol, reference_genome=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `gene`, `note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_gene_constraints(
    gene_symbol='BRCA1',
    reference_genome='GRCh38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene": {
    "symbol": "<string>",
    "gene_id": "<string>",
    "exac_constraint": "<object>",
    "gnomad_constraint": {
      "exp_lof": "<number>",
      "obs_lof": "<integer>",
      "oe_lof": "<number>",
      "pLI": "<number>",
      "exp_mis": "<number>",
      "obs_mis": "<integer>",
      "oe_mis": "<number>",
      "exp_syn": "<number>"
    },
    "mitochondrial_constraint": "<object>"
  },
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_gene_constraints-parameters)

## 9) Retrieve gnomAD variant metadata {#dnakit.databases.gnomad_get_variant}

- **Purpose:** Retrieve gnomAD variant metadata.
- **Method:** Use `gnomADGraphQLQuery` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_variant(*, variant_id, dataset=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `variant`, `dataset`, `reference_genome`, `dataset_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_variant(
    variant_id='19-44908822-C-T',
    dataset='gnomad_r3',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": {
    "variant_id": "<string>",
    "chrom": "<string>",
    "pos": "<integer>",
    "ref": "<string>",
    "alt": "<string>",
    "rsid": "<string>",
    "genome": {
      "af": "<number>",
      "an": "<integer>",
      "ac": "<integer>"
    },
    "exome": {
      "af": "<number>",
      "an": "<integer>",
      "ac": "<integer>"
    }
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_variant-parameters)

## 10) Retrieve ancestry-specific gnomAD allele frequencies {#dnakit.databases.gnomad_get_variant_populations}

- **Purpose:** Retrieve ancestry-specific gnomAD allele frequencies.
- **Method:** Use `gnomADGetVariantPopulations` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_variant_populations(*, variant_id, dataset=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `variant_id`, `chrom`, `pos`, `ref` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_variant_populations(
    variant_id='1-55051215-G-GA',
    dataset='gnomad_r4',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant_id": "<string>",
  "chrom": "<string>",
  "pos": "<integer>",
  "ref": "<string>",
  "alt": "<string>",
  "rsid": "<string>",
  "dataset": "<string>",
  "genome": {
    "ac": "<integer>",
    "an": "<integer>",
    "af": "<number>",
    "populations": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_variant_populations-parameters)

## 11) Search gnomAD variant identifiers {#dnakit.databases.gnomad_search_variants}

- **Purpose:** Search gnomAD variant identifiers.
- **Method:** Use `gnomADGraphQLQuery` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_search_variants(*, query=None, dataset=None, variant=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `variant_search`, `dataset`, `reference_genome`, `dataset_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_search_variants(
    query='rs7412',
    dataset='gnomad_r3',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant_search": [
    {
      "variant_id": "<string>"
    }
  ],
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_search_variants-parameters)

<!-- END GENERATED API: variant-records -->
