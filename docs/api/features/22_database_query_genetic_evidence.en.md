# Gene–Disease Evidence, Dosage and Actionability

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: genetic-evidence -->

<span id="genetic-evidence"></span>

## 1) Retrieve ClinGen gene dosage sensitivity curations {#dnakit.databases.clin_gen_dosage_by_gene}

- **Purpose:** Retrieve ClinGen gene dosage sensitivity curations.
- **Method:** Use `ClinGenDosage` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_dosage_by_gene(*, gene=None, gene_symbol=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `symbol`, `hgnc_id`, `location`, `grch37` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_dosage_by_gene(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "symbol": "<string>",
    "hgnc_id": "<string>",
    "location": "<string>",
    "grch37": "<string>",
    "grch38": "<string>",
    "haplo_assertion": "<integer>",
    "triplo_assertion": "<integer>",
    "pli": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_dosage_by_gene-parameters)

## 2) Search ClinGen dosage sensitivity by genomic interval {#dnakit.databases.clin_gen_dosage_region_search}

- **Purpose:** Search ClinGen dosage sensitivity by genomic interval.
- **Method:** Use `ClinGenDosage` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_dosage_region_search(*, chromosome, start, end, assembly=None, verbose=True)`
- **Input:** Provide `chromosome`, `start`, `end`.
- **Output:** Returns a list of records; each record includes `symbol`, `hgnc_id`, `type`, `location` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_dosage_region_search(
    chromosome='17',
    start=43044295,
    end=43125370,
    assembly='GRCh38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "symbol": "<string>",
    "hgnc_id": "<string>",
    "type": "<integer>",
    "location": "<string>",
    "grch37": "<string>",
    "grch38": "<string>",
    "haplo_assertion": "<integer>",
    "triplo_assertion": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_dosage_region_search-parameters)

## 3) Retrieve adult clinical actionability curations {#dnakit.databases.clin_gen_get_actionability_adult}

- **Purpose:** Retrieve adult clinical actionability curations.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_get_actionability_adult(*, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `docId`, `topicIri`, `curationType`, `latestSearchDate` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_get_actionability_adult(
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "docId": "<string>",
    "topicIri": "<string>",
    "curationType": "<string>",
    "latestSearchDate": "<string>",
    "lastUpdated": "<string>",
    "lastAuthor": "<string>",
    "context": "<string>",
    "contextIri": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_actionability_adult-parameters)

## 4) Retrieve pediatric clinical actionability curations {#dnakit.databases.clin_gen_get_actionability_pediatric}

- **Purpose:** Retrieve pediatric clinical actionability curations.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_get_actionability_pediatric(*, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `docId`, `topicIri`, `curationType`, `latestSearchDate` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_get_actionability_pediatric(
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "docId": "<string>",
    "topicIri": "<string>",
    "curationType": "<string>",
    "latestSearchDate": "<string>",
    "lastUpdated": "<string>",
    "lastAuthor": "<string>",
    "context": "<string>",
    "contextIri": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_actionability_pediatric-parameters)

## 5) Retrieve ClinGen dosage sensitivity records {#dnakit.databases.clin_gen_get_dosage_sensitivity}

- **Purpose:** Retrieve ClinGen dosage sensitivity records.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_get_dosage_sensitivity(*, gene=None, include_regions=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “Retrieve ClinGen dosage sensitivity records”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_get_dosage_sensitivity(
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<object>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_dosage_sensitivity-parameters)

## 6) Retrieve ClinGen gene–disease validity curations {#dnakit.databases.clin_gen_get_gene_validity}

- **Purpose:** Retrieve ClinGen gene–disease validity curations.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_get_gene_validity(*, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “Retrieve ClinGen gene–disease validity curations”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_get_gene_validity(
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<object>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_gene_validity-parameters)

## 7) Retrieve ClinGen variant classifications {#dnakit.databases.clin_gen_get_variant_classifications}

- **Purpose:** Retrieve ClinGen variant classifications.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_get_variant_classifications(*, gene=None, variant=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `Variation`, `ClinVar Variation Id`, `HGNC Gene Symbol`, `Disease` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_get_variant_classifications(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "Variation": "<string>",
    "ClinVar Variation Id": "<string>",
    "HGNC Gene Symbol": "<string>",
    "Disease": "<string>",
    "Mondo Id": "<string>",
    "Assertion": "<string>",
    "Expert Panel": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_variant_classifications-parameters)

## 8) Search gene actionability across adult and pediatric contexts {#dnakit.databases.clin_gen_search_actionability}

- **Purpose:** Search gene actionability across adult and pediatric contexts.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_search_actionability(*, gene, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a mapping with `data`, `gene_searched`, `adult_count`, `pediatric_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_search_actionability(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "Adult": [
      "<value>"
    ],
    "Pediatric": [
      "<value>"
    ]
  },
  "gene_searched": "<string>",
  "adult_count": "<integer>",
  "pediatric_count": "<integer>",
  "failed_contexts": "<object>",
  "note": "<string>",
  "source": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_actionability-parameters)

## 9) Search dosage sensitivity scores by gene symbol {#dnakit.databases.clin_gen_search_dosage_sensitivity}

- **Purpose:** Search dosage sensitivity scores by gene symbol.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_search_dosage_sensitivity(*, gene, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a list of records; each record includes `GENE SYMBOL`, `HGNC ID`, `HAPLOINSUFFICIENCY`, `TRIPLOSENSITIVITY` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_search_dosage_sensitivity(
    gene='MECP2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "GENE SYMBOL": "<string>",
    "HGNC ID": "<string>",
    "HAPLOINSUFFICIENCY": "<string>",
    "TRIPLOSENSITIVITY": "<string>",
    "ONLINE REPORT": "<string>",
    "DATE": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_dosage_sensitivity-parameters)

## 10) Search gene–disease validity by gene symbol {#dnakit.databases.clin_gen_search_gene_validity}

- **Purpose:** Search gene–disease validity by gene symbol.
- **Method:** Use `ClinGen` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clin_gen_search_gene_validity(*, gene, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a list of records; each record includes `GENE SYMBOL`, `GENE ID (HGNC)`, `DISEASE LABEL`, `DISEASE ID (MONDO)` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clin_gen_search_gene_validity(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "GENE SYMBOL": "<string>",
    "GENE ID (HGNC)": "<string>",
    "DISEASE LABEL": "<string>",
    "DISEASE ID (MONDO)": "<string>",
    "CLASSIFICATION": "<string>",
    "MOI": "<string>",
    "SOP": "<string>",
    "GCEP": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_gene_validity-parameters)

<!-- END GENERATED API: genetic-evidence -->
