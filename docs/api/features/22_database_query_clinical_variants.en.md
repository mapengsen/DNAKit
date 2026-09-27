# ClinVar Clinical Significance

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

## 1) ClinVar variant search {#4-clinvar-variant-search}

- **Purpose:** Search ClinVar by gene, condition, ClinVar ID, HGVS, or protein change.
- **Method:** Search ClinVar with the supplied filters and result limit.
- **API:** `dnakit.annotation.search_clinvar_variants(gene[optional], condition[optional], variant_id[optional], variant_name[optional], clinical_significance[optional], limit[optional])`
- **Input:** At least one search field; `limit` must be 1–100.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import search_clinvar_variants

result = search_clinvar_variants(gene="TP53", variant_name="c.215C>G", limit=10)
print(result.data)
```

- **Example output:** ClinVar identifiers and summary records.

## 2) ClinVar variant details {#5-clinvar-variant-details}

- **Purpose:** Retrieve accession, gene, coordinates, and review status for one ClinVar variation.
- **Method:** Retrieve the ClinVar detail record by variation ID.
- **API:** `dnakit.annotation.get_clinvar_variant(variant_id[required])`
- **Input:** A ClinVar variation ID such as `12345`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_clinvar_variant

result = get_clinvar_variant("12345")
print(result.data)
```

- **Example output:** ClinVar variant details.

## 3) ClinVar clinical significance {#6-clinvar-clinical-significance}

- **Purpose:** Retrieve clinical-significance assertions for one ClinVar variation.
- **Method:** Read clinical significance and interpretation from the selected ClinVar record.
- **API:** `dnakit.annotation.get_clinvar_significance(variant_id[required])`
- **Input:** A ClinVar variation ID.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import get_clinvar_significance

result = get_clinvar_significance("12345")
print(result.data)
```

- **Example output:** Pathogenicity classifications and clinical interpretations.

<!-- BEGIN GENERATED API: clinical-variants -->

<span id="clinical-variants"></span>

## 4) Retrieve ClinVar clinical significance {#dnakit.databases.clin_var_get_clinical_significance}

- **Purpose:** Retrieve ClinVar clinical significance.
- **Method:** Use `ClinVarGetClinicalSignificance` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_var_get_clinical_significance(*, variant_id, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `status`, `data`, `variant_id`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_var_get_clinical_significance(
    variant_id='12345',
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
    "variant_id": "<string>",
    "germline_classification": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>",
      "fda_recognized": "<string>",
      "traits": "<array>"
    },
    "clinical_impact": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>"
    },
    "oncogenicity": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>"
    },
    "raw_data": "<object>"
  },
  "variant_id": "<string>",
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_get_clinical_significance-parameters)

## 5) Retrieve ClinVar variant summaries and review status {#dnakit.databases.clin_var_get_variant_details}

- **Purpose:** Retrieve ClinVar variant summaries and review status.
- **Method:** Use `ClinVarGetVariantDetails` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_var_get_variant_details(*, variant_id, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `status`, `data`, `variant_id`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_var_get_variant_details(
    variant_id='12345',
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
    "variant_id": "<string>",
    "accession": "<string>",
    "obj_type": "<string>",
    "chromosome": "<string>",
    "location": "<string>",
    "variation_name": "<string>",
    "title": "<string>",
    "genes": [
      "<string>"
    ]
  },
  "variant_id": "<string>",
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_get_variant_details-parameters)

## 6) Search ClinVar variants by gene, disease or identifier {#dnakit.databases.clin_var_search_variants}

- **Purpose:** Search ClinVar variants by gene, disease or identifier.
- **Method:** Use `ClinVarSearchVariants` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_var_search_variants(*, gene=None, condition=None, variant_id=None, variant_name=None, max_results=None, limit=None, clinical_significance=None, gene_symbol=None, significance=None, query=None, variant=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data`, `search_params`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_var_search_variants(
    gene='BRCA1',
    max_results=10,
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
    "esearchresult": {
      "count": "<string>",
      "idlist": "<array>",
      "querytranslation": "<string>"
    }
  },
  "search_params": {
    "gene": "<string>",
    "condition": "<string>",
    "variant_id": "<string>"
  },
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_search_variants-parameters)

<!-- END GENERATED API: clinical-variants -->
