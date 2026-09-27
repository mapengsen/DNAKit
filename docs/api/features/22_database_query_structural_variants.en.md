# Structural Variant Queries

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: structural-variants -->

<span id="structural-variants"></span>

## 1) Retrieve Ensembl structural variants in an interval {#dnakit.databases.ensembl_get_structural_variants}

- **Purpose:** Retrieve Ensembl structural variants in an interval.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_structural_variants(*, species, region=None, feature=None, chrom=None, start=None, end=None, verbose=True)`
- **Input:** Provide `species`.
- **Output:** Returns a list of records; each record includes `id`, `seq_region_name`, `start`, `end` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_structural_variants(
    species='human',
    region='22:27000000-27050000',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "id": "<string>",
    "seq_region_name": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "strand": "<integer>",
    "feature_type": "<string>",
    "assembly_name": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.ensembl_get_structural_variants-parameters)

## 2) Retrieve Ensembl structural variant details {#dnakit.databases.ensembl_get_sv_detail}

- **Purpose:** Retrieve Ensembl structural variant details.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_sv_detail(*, species, id=None, sv_id=None, verbose=True)`
- **Input:** Provide `species`.
- **Output:** Returns a mapping with `name`, `var_class`, `source`, `clinical_significance` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_sv_detail(
    species='human',
    id='esv3647486',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "name": "<string>",
  "var_class": "<string>",
  "source": "<string>",
  "clinical_significance": [
    "<string>"
  ],
  "alias": "<string>",
  "supporting_evidence": [
    {
      "name": "<string>",
      "var_class": "<string>",
      "source": "<string>",
      "clinical_significance": "<array>",
      "copy_number": "<integer>",
      "study": "<object>",
      "mappings": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.ensembl_get_sv_detail-parameters)

## 3) Find gnomAD structural variants by gene {#dnakit.databases.gnomad_get_sv_by_gene}

- **Purpose:** Find gnomAD structural variants by gene.
- **Method:** Use `gnomADGraphQLQuery` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_sv_by_gene(*, gene_symbol, reference_genome=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `gene`, `dataset`, `reference_genome`, `dataset_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_by_gene(
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
    "chrom": "<string>",
    "start": "<integer>",
    "stop": "<integer>",
    "structural_variants": [
      "<object>"
    ]
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_by_gene-parameters)

## 4) Find gnomAD structural variants by genomic interval {#dnakit.databases.gnomad_get_sv_by_region}

- **Purpose:** Find gnomAD structural variants by genomic interval.
- **Method:** Use `gnomADGraphQLQuery` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_sv_by_region(*, chrom, start, stop=None, reference_genome=None, end=None, verbose=True)`
- **Input:** Provide `chrom`, `start`.
- **Output:** Returns a mapping with `region`, `dataset`, `reference_genome`, `dataset_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_by_region(
    chrom='17',
    start=43044295,
    stop=43070295,
    reference_genome='GRCh38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "region": {
    "chrom": "<string>",
    "start": "<integer>",
    "stop": "<integer>",
    "genes": [
      "<object>"
    ],
    "structural_variants": [
      "<object>"
    ]
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_by_region-parameters)

## 5) Retrieve gnomAD structural variant details {#dnakit.databases.gnomad_get_sv_detail}

- **Purpose:** Retrieve gnomAD structural variant details.
- **Method:** Use `gnomADGraphQLQuery` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gnomad_get_sv_detail(*, variant_id=None, dataset=None, sv_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `structural_variant`, `dataset`, `reference_genome`, `dataset_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_detail(
    variant_id='DEL_chr17_24e4872b',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "structural_variant": {
    "variant_id": "<string>",
    "chrom": "<string>",
    "pos": "<integer>",
    "end": "<integer>",
    "type": "<string>",
    "length": "<integer>",
    "ac": "<integer>",
    "an": "<integer>"
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_detail-parameters)

<!-- END GENERATED API: structural-variants -->
