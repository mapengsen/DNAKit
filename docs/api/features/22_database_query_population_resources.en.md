# Populations, Samples and Linkage Disequilibrium

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: population-resources -->

<span id="population-resources"></span>

## 1) List IGSR and 1000 Genomes data collections {#dnakit.databases.igsr_list_data_collections}

- **Purpose:** List IGSR and 1000 Genomes data collections.
- **Method:** Use `IGSR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.igsr_list_data_collections(*, operation=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `returned`, `collections` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.igsr_list_data_collections(
    operation='list_data_collections',
    limit=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "collections": [
    {
      "code": "<string>",
      "title": "<string>",
      "short_title": "<string>",
      "sample_count": "<integer>",
      "population_count": "<integer>",
      "data_types": "<array>",
      "website": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_list_data_collections-parameters)

## 2) Search 1000 Genomes populations {#dnakit.databases.igsr_search_populations}

- **Purpose:** Search 1000 Genomes populations.
- **Method:** Use `IGSR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.igsr_search_populations(*, operation=None, superpopulation=None, query=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `returned`, `populations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.igsr_search_populations(
    operation='search_populations',
    superpopulation='EAS',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "populations": [
    {
      "code": "<string>",
      "name": "<string>",
      "description": "<string>",
      "sample_count": "<integer>",
      "superpopulation_code": "<string>",
      "superpopulation_name": "<string>",
      "latitude": "<number>",
      "longitude": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_search_populations-parameters)

## 3) Search 1000 Genomes samples by population or collection {#dnakit.databases.igsr_search_samples}

- **Purpose:** Search 1000 Genomes samples by population or collection.
- **Method:** Use `IGSR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.igsr_search_samples(*, operation=None, population=None, data_collection=None, sample_name=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `returned`, `samples` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.igsr_search_samples(
    operation='search_samples',
    population='GBR',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "samples": [
    {
      "name": "<string>",
      "sex": "<string>",
      "biosample_id": "<string>",
      "populations": "<array>",
      "data_collections": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_search_samples-parameters)

## 4) Find population-specific LD proxy variants {#dnakit.databases.ldlink_get_proxies}

- **Purpose:** Find population-specific LD proxy variants.
- **Method:** Use `LDlink` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ldlink_get_proxies(*, variant, population=None, r2_threshold=None, genome_build=None, limit=None, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a list of records; each record includes `rsid`, `r2`, `coord`, `alleles` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ldlink_get_proxies(
    variant='rs7903146',
    population='CEU',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "rsid": "<string>",
    "r2": "<number>",
    "coord": "<string>",
    "alleles": "<string>",
    "maf": "<string>",
    "distance": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_population_resources.md#dnakit.databases.ldlink_get_proxies-parameters)

<!-- END GENERATED API: population-resources -->
