# Cohort Phenotype Associations and Fine-mapping

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: phewas -->

<span id="phewas"></span>

## 1) Find variant–phenotype associations in BioBank Japan {#dnakit.databases.bio_bank_japan_phewas_by_variant}

- **Purpose:** Find variant–phenotype associations in BioBank Japan.
- **Method:** Use `PheWebPheWAS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.bio_bank_japan_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `variant`, `rsids`, `nearest_genes`, `associations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.bio_bank_japan_phewas_by_variant(
    rsid='rs7903146',
    max_pval=5e-08,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phewas.md#dnakit.databases.bio_bank_japan_phewas_by_variant-parameters)

## 2) Retrieve FinnGen fine-mapping results for a variant {#dnakit.databases.finn_gen_get_variant_finemapping}

- **Purpose:** Retrieve FinnGen fine-mapping results for a variant.
- **Method:** Use `FinnGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.finn_gen_get_variant_finemapping(*, variant, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a list of records; each record includes `phenocode`, `chromosome`, `start`, `end` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.finn_gen_get_variant_finemapping(
    variant='19:44908684:T:C',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "phenocode": "<string>",
    "chromosome": "<integer>",
    "start": "<integer>",
    "end": "<integer>",
    "type": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phewas.md#dnakit.databases.finn_gen_get_variant_finemapping-parameters)

## 3) Query gene burden associations across phenotypes {#dnakit.databases.genebass_gene_burden_phewas}

- **Purpose:** Query gene burden associations across phenotypes.
- **Method:** Use `Genebass` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.genebass_gene_burden_phewas(*, gene, burden_set=None, limit=None, max_pval=None, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a mapping with `gene_id`, `gene_symbol`, `burden_set`, `associations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.genebass_gene_burden_phewas(
    gene='TCF7L2',
    burden_set='pLoF',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_id": "<string>",
  "gene_symbol": "<string>",
  "burden_set": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "description": "<string>",
      "trait_type": "<string>",
      "pval": "<number>",
      "pval_skat": "<number>",
      "beta": "<number>",
      "total_variants": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phewas.md#dnakit.databases.genebass_gene_burden_phewas-parameters)

## 4) Query variant–phenotype associations in TPMI {#dnakit.databases.tpmi_phewas_by_variant}

- **Purpose:** Query variant–phenotype associations in TPMI.
- **Method:** Use `PheWebPheWAS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.tpmi_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `variant`, `rsids`, `nearest_genes`, `associations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.tpmi_phewas_by_variant(
    rsid='rs7903146',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phewas.md#dnakit.databases.tpmi_phewas_by_variant-parameters)

## 5) Query variant–phenotype associations in UKB-TOPMed {#dnakit.databases.ukbtopmed_phewas_by_variant}

- **Purpose:** Query variant–phenotype associations in UKB-TOPMed.
- **Method:** Use `PheWebPheWAS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ukbtopmed_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `variant`, `rsids`, `nearest_genes`, `associations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ukbtopmed_phewas_by_variant(
    rsid='rs7903146',
    max_pval=5e-08,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phewas.md#dnakit.databases.ukbtopmed_phewas_by_variant-parameters)

<!-- END GENERATED API: phewas -->
