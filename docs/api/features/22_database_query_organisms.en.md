# Organism, taxonomy and microbiome queries

Query model organisms, taxonomy and microbiome data.

## 1) Query taxonomy and TaxIDs {#search-taxonomy}

- **Purpose:** Resolve a species name, TaxID or alias into taxonomy records for subsequent queries and downloads.
- **Method:** Query NCBI Taxonomy and normalize its reports into QueryResult.
- **API:** `dnakit.search.taxonomy(query[required], report[optional], children[optional], include_lineage[optional], ranks[optional], page_size[optional], page_token[optional], config[optional], verbose[optional])`
- **Input:** A species name or TaxID; optional report type, children, lineage, ranks, pagination and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import SearchConfig, taxonomy

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = taxonomy("human", page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
taxonomy NCBI Datasets <record count returned>
```

<!-- BEGIN GENERATED API: organisms -->

<span id="organisms"></span>

## 2) Retrieve Drosophila gene details {#dnakit.databases.fly_base_get_gene}

- **Purpose:** Retrieve Drosophila gene details.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `id`, `symbol`, `name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene(
    gene_id='FB:FBgn0003996',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene-parameters)

## 3) Retrieve Drosophila gene alleles and phenotypes {#dnakit.databases.fly_base_get_gene_alleles}

- **Purpose:** Retrieve Drosophila gene alleles and phenotypes.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene_alleles(*, gene_id, limit=None, page=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `id`, `symbol`, `symbol_text`, `category` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_alleles(
    gene_id='FB:FBgn0000490',
    limit=5,
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
    "symbol": "<string>",
    "symbol_text": "<string>",
    "category": "<string>",
    "has_disease": "<boolean>",
    "has_phenotype": "<boolean>",
    "variants": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_alleles-parameters)

## 4) Retrieve Drosophila gene disease models {#dnakit.databases.fly_base_get_gene_disease_models}

- **Purpose:** Retrieve Drosophila gene disease models.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene_disease_models(*, gene_id, limit=None, page=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `model_id`, `model_name`, `gene_symbol`, `gene_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_disease_models(
    gene_id='FB:FBgn0000490',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "model_id": "<string>",
    "model_name": "<string>",
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "data_provider": "<string>",
    "diseases": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_disease_models-parameters)

## 5) Retrieve Drosophila developmental and anatomical expression {#dnakit.databases.fly_base_get_gene_expression}

- **Purpose:** Retrieve Drosophila developmental and anatomical expression.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene_expression(*, gene_id, limit=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `total_annotations`, `returned`, `annotations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_expression(
    gene_id='FB:FBgn0000490',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_annotations": "<integer>",
  "returned": "<integer>",
  "annotations": [
    {
      "stage": "<string>",
      "location": "<string>",
      "data_provider": "<string>",
      "reference": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_expression-parameters)

## 6) Retrieve Drosophila molecular and genetic interactions {#dnakit.databases.fly_base_get_gene_interactions}

- **Purpose:** Retrieve Drosophila molecular and genetic interactions.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene_interactions(*, gene_id, interaction_type=None, limit=None, page=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `subject_gene_id`, `subject_symbol`, `interactor_gene_id`, `interactor_symbol` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_interactions(
    gene_id='FB:FBgn0000490',
    interaction_type='genetic',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "subject_gene_id": "<string>",
    "subject_symbol": "<string>",
    "interactor_gene_id": "<string>",
    "interactor_symbol": "<string>",
    "interactor_species": "<string>",
    "interaction_type": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_interactions-parameters)

## 7) Retrieve cross-species orthologs of a Drosophila gene {#dnakit.databases.fly_base_get_gene_orthologs}

- **Purpose:** Retrieve cross-species orthologs of a Drosophila gene.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_base_get_gene_orthologs(*, gene_id, limit=None, page=None, stringency=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `subject_gene_id`, `subject_symbol`, `subject_species`, `ortholog_gene_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_orthologs(
    gene_id='FB:FBgn0004647',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "subject_gene_id": "<string>",
    "subject_symbol": "<string>",
    "subject_species": "<string>",
    "ortholog_gene_id": "<string>",
    "ortholog_symbol": "<string>",
    "ortholog_species": "<string>",
    "is_best_score": "<string>",
    "is_best_score_reverse": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_orthologs-parameters)

## 8) Search integrated Drosophila genomics records {#dnakit.databases.fly_mine_search}

- **Purpose:** Search integrated Drosophila genomics records.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fly_mine_search(*, q, size=None, format=None, facet_Category=None, verbose=True)`
- **Input:** Provide `q`.
- **Output:** Returns a mapping with `totalHits`, `facets`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fly_mine_search(
    q='Notch',
    size=5,
    format='json',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "totalHits": "<integer>",
  "facets": "<object>",
  "results": [
    {
      "id": "<integer>",
      "type": "<string>",
      "relevance": "<number>",
      "fields": "<object>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_mine_search-parameters)

## 9) Retrieve GBIF species details by taxon key {#dnakit.databases.gbif_get_species}

- **Purpose:** Retrieve GBIF species details by taxon key.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gbif_get_species(*, speciesKey, verbose=True)`
- **Input:** Provide `speciesKey`.
- **Output:** Returns a mapping with `key`, `scientificName`, `canonicalName`, `vernacularName` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gbif_get_species(
    speciesKey=2435099,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "key": "<integer>",
  "scientificName": "<string>",
  "canonicalName": "<string>",
  "vernacularName": "<string>",
  "authorship": "<string>",
  "rank": "<string>",
  "taxonomicStatus": "<string>",
  "kingdom": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gbif_get_species-parameters)

## 10) Search GBIF species names and taxonomy {#dnakit.databases.gbif_search_species}

- **Purpose:** Search GBIF species names and taxonomy.
- **Method:** Use `GBIF` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gbif_search_species(*, query, limit=None, offset=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** The data field contains the GBIF species search response and pagination information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gbif_search_species(
    query='Arabidopsis thaliana',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The data field contains the GBIF species search response and pagination information.

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gbif_search_species-parameters)

## 11) Retrieve GTDB species clusters and their genomes {#dnakit.databases.gtdb_get_species}

- **Purpose:** Retrieve GTDB species clusters and their genomes.
- **Method:** Use `GTDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtdb_get_species(*, operation=None, species, max_genomes=None, verbose=True)`
- **Input:** Provide `species`.
- **Output:** Returns a mapping with `species_name`, `total_genomes`, `genomes_shown`, `genomes` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtdb_get_species(
    operation='get_species',
    species='Escherichia coli',
    max_genomes=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "species_name": "<string>",
  "total_genomes": "<integer>",
  "genomes_shown": "<integer>",
  "genomes": [
    {
      "accession": "<string>",
      "ncbi_org_name": "<string>",
      "ncbi_tax": "<string>",
      "gtdb_species_rep": "<boolean>",
      "ncbi_type_material": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_get_species-parameters)

## 12) Retrieve GTDB taxonomy and release history {#dnakit.databases.gtdb_get_taxon_info}

- **Purpose:** Retrieve GTDB taxonomy and release history.
- **Method:** Use `GTDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtdb_get_taxon_info(*, operation=None, taxon, verbose=True)`
- **Input:** Provide `taxon`.
- **Output:** Returns a mapping with `taxon`, `rank`, `n_genomes`, `higher_ranks` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtdb_get_taxon_info(
    operation='get_taxon_info',
    taxon='f__Lachnospiraceae',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "taxon": "<string>",
  "rank": "<string>",
  "n_genomes": "<integer>",
  "higher_ranks": [
    "<string>"
  ],
  "in_releases": [
    "<string>"
  ],
  "lineage": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_get_taxon_info-parameters)

## 13) Search prokaryotic genomes by organism name {#dnakit.databases.gtdb_search_genomes}

- **Purpose:** Search prokaryotic genomes by organism name.
- **Method:** Use `GTDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtdb_search_genomes(*, operation=None, query, page=None, items_per_page=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `query`, `match_type`, `total_results`, `total_results_scope` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtdb_search_genomes(
    operation='search_genomes',
    query='Escherichia',
    items_per_page=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": "<string>",
  "match_type": "<string>",
  "total_results": "<integer>",
  "total_results_scope": "<string>",
  "note": "<string>",
  "page": "<integer>",
  "results": [
    {
      "gid": "<string>",
      "accession": "<string>",
      "ncbiOrgName": "<string>",
      "ncbiTaxonomy": "<string>",
      "gtdbTaxonomy": "<string>",
      "isGtdbSpeciesRep": "<boolean>",
      "isNcbiTypeMaterial": "<boolean>"
    }
  ],
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_search_genomes-parameters)

## 14) Search prokaryotic taxa by partial name {#dnakit.databases.gtdb_search_taxon}

- **Purpose:** Search prokaryotic taxa by partial name.
- **Method:** Use `GTDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtdb_search_taxon(*, operation=None, query, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `query`, `matches`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtdb_search_taxon(
    operation='search_taxon',
    query='Lachnospiraceae',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": "<string>",
  "matches": [
    "<string>"
  ],
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_search_taxon-parameters)

## 15) Retrieve mouse gene details by MGI identifier {#dnakit.databases.mgi_get_gene}

- **Purpose:** Retrieve mouse gene details by MGI identifier.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgi_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `id`, `symbol`, `name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgi_get_gene(
    gene_id='MGI:98834',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_get_gene-parameters)

## 16) Retrieve mouse gene mutant phenotype annotations {#dnakit.databases.mgi_get_phenotypes}

- **Purpose:** Retrieve mouse gene mutant phenotype annotations.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgi_get_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `gene_symbol`, `gene_id`, `phenotype_statement` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgi_get_phenotypes(
    gene_id='MGI:98834',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "phenotype_statement": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_get_phenotypes-parameters)

## 17) Search mouse genes by name or keyword {#dnakit.databases.mgi_search_genes}

- **Purpose:** Search mouse genes by name or keyword.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgi_search_genes(*, query, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `symbol`, `name`, `gene_id`, `category` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgi_search_genes(
    query='Trp53',
    limit=5,
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
    "name": "<string>",
    "gene_id": "<string>",
    "category": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_search_genes-parameters)

## 18) Retrieve microbial genome quality and taxonomy {#dnakit.databases.mgnify_get_genome}

- **Purpose:** Retrieve microbial genome quality and taxonomy.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_get_genome(*, genome_id, verbose=True)`
- **Input:** Provide `genome_id`.
- **Output:** Returns a mapping with `genome_id`, `type`, `taxonomy`, `completeness` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_get_genome(
    genome_id='MGYG000000001',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "genome_id": "<string>",
  "type": "<string>",
  "taxonomy": "<string>",
  "completeness": "<number>",
  "contamination": "<number>",
  "length": "<integer>",
  "num_proteins": "<integer>",
  "gc_content": "<number>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_genome-parameters)

## 19) Retrieve GO annotations from a metagenomics analysis {#dnakit.databases.mgnify_get_go_terms}

- **Purpose:** Retrieve GO annotations from a metagenomics analysis.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_get_go_terms(*, analysis_id, page_size=None, verbose=True)`
- **Input:** Provide `analysis_id`.
- **Output:** Returns a list of records; each record includes `go_id`, `description`, `count`, `category` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_get_go_terms(
    analysis_id='MGYA00585482',
    page_size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "go_id": "<string>",
    "description": "<string>",
    "count": "<integer>",
    "category": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_go_terms-parameters)

## 20) Retrieve metagenomic protein domain annotations {#dnakit.databases.mgnify_get_interpro}

- **Purpose:** Retrieve metagenomic protein domain annotations.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_get_interpro(*, analysis_id, page_size=None, verbose=True)`
- **Input:** Provide `analysis_id`.
- **Output:** Returns a list of records; each record includes `interpro_id`, `description`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_get_interpro(
    analysis_id='MGYA00585482',
    page_size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "interpro_id": "<string>",
    "description": "<string>",
    "count": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_interpro-parameters)

## 21) Retrieve MGnify metagenomics study details {#dnakit.databases.mgnify_get_study_detail}

- **Purpose:** Retrieve MGnify metagenomics study details.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_get_study_detail(*, study_accession, verbose=True)`
- **Input:** Provide `study_accession`.
- **Output:** Returns a mapping with `study_id`, `study_name`, `study_abstract`, `bioproject` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_get_study_detail(
    study_accession='MGYS00002008',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "study_id": "<string>",
  "study_name": "<string>",
  "study_abstract": "<string>",
  "bioproject": "<string>",
  "centre_name": "<string>",
  "biomes": [
    "<string>"
  ],
  "is_public": "<boolean>",
  "samples_count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_study_detail-parameters)

## 22) Retrieve metagenomic taxonomic composition and abundance {#dnakit.databases.mgnify_get_taxonomy}

- **Purpose:** Retrieve metagenomic taxonomic composition and abundance.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_get_taxonomy(*, analysis_id, rna_type=None, page_size=None, verbose=True)`
- **Input:** Provide `analysis_id`.
- **Output:** Returns a list of records; each record includes `organism_id`, `name`, `rank`, `domain` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_get_taxonomy(
    analysis_id='MGYA00585482',
    page_size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "organism_id": "<string>",
    "name": "<string>",
    "rank": "<string>",
    "domain": "<string>",
    "count": "<integer>",
    "lineage": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_taxonomy-parameters)

## 23) List analyses belonging to a metagenomics study {#dnakit.databases.mgnify_list_analyses}

- **Purpose:** List analyses belonging to a metagenomics study.
- **Method:** Use `MGnifyAnalyses` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_list_analyses(*, study_accession, size=None, verbose=True)`
- **Input:** Provide `study_accession`.
- **Output:** data contains MGnify analysis records for the study and the service's pagination information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_list_analyses(
    study_accession='MGYS00000596',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data contains MGnify analysis records for the study and the service's pagination information.

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_list_analyses-parameters)

## 24) Browse the metagenomics biome hierarchy {#dnakit.databases.mgnify_list_biomes}

- **Purpose:** Browse the metagenomics biome hierarchy.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_list_biomes(*, depth=None, page=None, page_size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `biome_id`, `biome_name`, `samples_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_list_biomes(
    page_size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "biome_id": "<string>",
    "biome_name": "<string>",
    "samples_count": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_list_biomes-parameters)

## 25) Search metagenome-assembled and isolate genomes {#dnakit.databases.mgnify_search_genomes}

- **Purpose:** Search metagenome-assembled and isolate genomes.
- **Method:** Use `MGnifyExpanded` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_search_genomes(*, taxonomy=None, page=None, page_size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `genome_id`, `type`, `taxonomy`, `completeness` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_search_genomes(
    page=1,
    page_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "genome_id": "<string>",
    "type": "<string>",
    "taxonomy": "<string>",
    "completeness": "<number>",
    "contamination": "<number>",
    "length": "<integer>",
    "num_proteins": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_search_genomes-parameters)

## 26) Search metagenomics studies by biome and keyword {#dnakit.databases.mgnify_search_studies}

- **Purpose:** Search metagenomics studies by biome and keyword.
- **Method:** Use `MGnifyStudies` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mgnify_search_studies(*, biome=None, search=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data contains MGnify study search results and the service's pagination information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mgnify_search_studies(
    search='human gut',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data contains MGnify study search results and the service's pagination information.

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_search_studies-parameters)

## 27) Search global plant names and taxonomy {#dnakit.databases.powo_search_plants}

- **Purpose:** Search global plant names and taxonomy.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.powo_search_plants(*, q, f=None, perPage=None, cursor=None, verbose=True)`
- **Input:** Provide `q`.
- **Output:** Returns a mapping with `totalResults`, `page`, `totalPages`, `perPage` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.powo_search_plants(
    q='Rosa canina',
    f='species_f:true',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "totalResults": "<integer>",
  "page": "<integer>",
  "totalPages": "<integer>",
  "perPage": "<integer>",
  "cursor": "<string>",
  "results": [
    {
      "accepted": "<boolean>",
      "author": "<string>",
      "kingdom": "<string>",
      "family": "<string>",
      "name": "<string>",
      "rank": "<string>",
      "snippet": "<string>",
      "url": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.powo_search_plants-parameters)

## 28) Retrieve rat gene disease, phenotype and functional annotations {#dnakit.databases.rgd_get_annotations}

- **Purpose:** Retrieve rat gene disease, phenotype and functional annotations.
- **Method:** Use `RGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_annotations(*, rgd_id, aspect=None, verbose=True)`
- **Input:** Provide `rgd_id`.
- **Output:** Returns a list for “Retrieve rat gene disease, phenotype and functional annotations”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_annotations(
    rgd_id='2218',
    aspect='D',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_annotations-parameters)

## 29) Retrieve rat gene details by RGD identifier {#dnakit.databases.rgd_get_gene}

- **Purpose:** Retrieve rat gene details by RGD identifier.
- **Method:** Use `RGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_gene(*, rgd_id, verbose=True)`
- **Input:** Provide `rgd_id`.
- **Output:** Returns rat gene names, functional description, gene type and Ensembl cross-references.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_gene(
    rgd_id='2218',
    verbose=False,
)
print(result)
```

- **Example output:** Returns rat gene names, functional description, gene type and Ensembl cross-references.

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_gene-parameters)

## 30) Retrieve cross-species orthologs of a rat gene {#dnakit.databases.rgd_get_orthologs}

- **Purpose:** Retrieve cross-species orthologs of a rat gene.
- **Method:** Use `RGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_orthologs(*, rgd_id, verbose=True)`
- **Input:** Provide `rgd_id`.
- **Output:** Returns a list for “Retrieve cross-species orthologs of a rat gene”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_orthologs(
    rgd_id='2218',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_orthologs-parameters)

## 31) Retrieve rat QTLs overlapping a genomic interval {#dnakit.databases.rgd_get_qtls_in_region}

- **Purpose:** Retrieve rat QTLs overlapping a genomic interval.
- **Method:** Use `RGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_qtls_in_region(*, chromosome, start, stop, map_key, verbose=True)`
- **Input:** Provide `chromosome`, `start`, `stop`, `map_key`.
- **Output:** Returns a list of records; each record includes `rgd_id`, `symbol`, `name`, `chromosome` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_qtls_in_region(
    chromosome='10',
    start=1,
    stop=50000000,
    map_key=360,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "rgd_id": "<integer>",
    "symbol": "<string>",
    "name": "<string>",
    "chromosome": "<string>",
    "lod": "<number>",
    "p_value": "<number>",
    "variance": "<number>",
    "inheritance_type": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_qtls_in_region-parameters)

## 32) Retrieve rat strain details by strain symbol {#dnakit.databases.rgd_get_strain}

- **Purpose:** Retrieve rat strain details by strain symbol.
- **Method:** Use `RGDStrain` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_strain(*, symbol, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns a mapping with `rgd_id`, `symbol`, `name`, `full_strain_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_strain(
    symbol='SHR',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "rgd_id": "<integer>",
  "symbol": "<string>",
  "name": "<string>",
  "full_strain_name": "<string>",
  "substrain": "<string>",
  "strain_type": "<string>",
  "genetics": "<string>",
  "origin": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_strain-parameters)

## 33) Retrieve rat strain disease and phenotype annotations {#dnakit.databases.rgd_get_strain_annotations}

- **Purpose:** Retrieve rat strain disease and phenotype annotations.
- **Method:** Use `RGDStrain` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_get_strain_annotations(*, symbol, category=None, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns a list of records; each record includes `term`, `term_acc`, `category`, `aspect` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_get_strain_annotations(
    symbol='SHR',
    category='disease',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "term": "<string>",
    "term_acc": "<string>",
    "category": "<string>",
    "aspect": "<string>",
    "qualifier": "<string>",
    "evidence": "<string>",
    "data_source": "<string>",
    "ref_rgd_id": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_strain_annotations-parameters)

## 34) Search RGD genes by name or keyword {#dnakit.databases.rgd_search_genes}

- **Purpose:** Search RGD genes by name or keyword.
- **Method:** Use `RGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_search_genes(*, query, species=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list for “Search RGD genes by name or keyword”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_search_genes(
    query='Brca1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_search_genes-parameters)

## 35) Search rat strains by keyword and type {#dnakit.databases.rgd_search_strains}

- **Purpose:** Search rat strains by keyword and type.
- **Method:** Use `RGDStrain` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rgd_search_strains(*, query=None, strain_type=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `rgd_id`, `symbol`, `name`, `full_strain_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rgd_search_strains(
    query='hypertensive',
    strain_type='inbred',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "rgd_id": "<integer>",
    "symbol": "<string>",
    "name": "<string>",
    "full_strain_name": "<string>",
    "substrain": "<string>",
    "strain_type": "<string>",
    "genetics": "<string>",
    "origin": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_search_strains-parameters)

## 36) Retrieve yeast gene and locus overviews {#dnakit.databases.sgd_get_gene}

- **Purpose:** Retrieve yeast gene and locus overviews.
- **Method:** Use `SGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sgd_get_gene(*, sgd_id, verbose=True)`
- **Input:** Provide `sgd_id`.
- **Output:** Returns a mapping with `sgd_id`, `display_name`, `gene_name`, `systematic_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sgd_get_gene(
    sgd_id='S000000259',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "sgd_id": "<string>",
  "display_name": "<string>",
  "gene_name": "<string>",
  "systematic_name": "<string>",
  "locus_type": "<string>",
  "description": "<string>",
  "uniprot_id": "<string>",
  "aliases": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_gene-parameters)

## 37) Retrieve yeast gene GO annotations and evidence {#dnakit.databases.sgd_get_go_annotations}

- **Purpose:** Retrieve yeast gene GO annotations and evidence.
- **Method:** Use `SGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sgd_get_go_annotations(*, sgd_id, verbose=True)`
- **Input:** Provide `sgd_id`.
- **Output:** Returns a list of records; each record includes `go_id`, `go_term`, `go_aspect`, `evidence_code` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sgd_get_go_annotations(
    sgd_id='S000003219',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "go_id": "<string>",
    "go_term": "<string>",
    "go_aspect": "<string>",
    "evidence_code": "<string>",
    "qualifier": "<string>",
    "annotation_type": "<string>",
    "reference": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_go_annotations-parameters)

## 38) Retrieve yeast genetic and physical interactions {#dnakit.databases.sgd_get_interactions}

- **Purpose:** Retrieve yeast genetic and physical interactions.
- **Method:** Use `SGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sgd_get_interactions(*, sgd_id, verbose=True)`
- **Input:** Provide `sgd_id`.
- **Output:** Returns a list of records; each record includes `interaction_type`, `experiment_type`, `bait_gene`, `hit_gene` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sgd_get_interactions(
    sgd_id='S000003219',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "interaction_type": "<string>",
    "experiment_type": "<string>",
    "bait_gene": "<string>",
    "hit_gene": "<string>",
    "annotation_type": "<string>",
    "reference": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_interactions-parameters)

## 39) Retrieve yeast mutant phenotypes and experimental conditions {#dnakit.databases.sgd_get_phenotypes}

- **Purpose:** Retrieve yeast mutant phenotypes and experimental conditions.
- **Method:** Use `SGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sgd_get_phenotypes(*, sgd_id, verbose=True)`
- **Input:** Provide `sgd_id`.
- **Output:** Returns a list of records; each record includes `phenotype`, `mutant_type`, `experiment_type`, `strain_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sgd_get_phenotypes(
    sgd_id='S000003219',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "phenotype": "<string>",
    "mutant_type": "<string>",
    "experiment_type": "<string>",
    "strain_name": "<string>",
    "allele": "<string>",
    "chemical": "<string>",
    "reference": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_phenotypes-parameters)

## 40) Search the yeast database across data types {#dnakit.databases.sgd_search}

- **Purpose:** Search the yeast database across data types.
- **Method:** Use `SGD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sgd_search(*, query, category=None, limit=None, offset=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `name`, `category`, `description`, `href` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sgd_search(
    query='GAL4',
    category='locus',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "name": "<string>",
    "category": "<string>",
    "description": "<string>",
    "href": "<string>",
    "aliases": [
      "<value>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_search-parameters)

## 41) Retrieve plant morphology, growth and reproduction characteristics {#dnakit.databases.usda_plants_get_characteristics}

- **Purpose:** Retrieve plant morphology, growth and reproduction characteristics.
- **Method:** Use `USDAPlantsCharacteristics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.usda_plants_get_characteristics(*, symbol, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns a list for “Retrieve plant morphology, growth and reproduction characteristics”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.usda_plants_get_characteristics(
    symbol='ABBA',
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

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.usda_plants_get_characteristics-parameters)

## 42) Retrieve plant names, taxonomy and growth habits {#dnakit.databases.usda_plants_get_profile}

- **Purpose:** Retrieve plant names, taxonomy and growth habits.
- **Method:** Use `USDAPlantsProfile` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.usda_plants_get_profile(*, symbol, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns scientific/common names, taxonomy, growth habits and regional native status; fine-grained distributions are not always available.
- **Example code:**

```python
import dnakit

result = dnakit.databases.usda_plants_get_profile(
    symbol='ABBA',
    verbose=False,
)
print(result)
```

- **Example output:** Returns scientific/common names, taxonomy, growth habits and regional native status; fine-grained distributions are not always available.

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.usda_plants_get_profile-parameters)

## 43) Retrieve nematode tissue and developmental expression {#dnakit.databases.worm_base_get_expression}

- **Purpose:** Retrieve nematode tissue and developmental expression.
- **Method:** Use `WormBase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.worm_base_get_expression(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `wormbase_id`, `gene_name`, `expressed_in_count`, `expressed_in` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.worm_base_get_expression(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "expressed_in_count": "<integer>",
  "expressed_in": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "expressed_during": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "subcellular_localization": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "expression_clusters_count": "<integer>",
  "expression_clusters": [
    {
      "cluster_id": "<string>",
      "cluster_label": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_expression-parameters)

## 44) Retrieve C. elegans gene details {#dnakit.databases.worm_base_get_gene}

- **Purpose:** Retrieve C. elegans gene details.
- **Method:** Use `WormBase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.worm_base_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `wormbase_id`, `gene_name`, `sequence_name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.worm_base_get_gene(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "sequence_name": "<string>",
  "species": "<string>",
  "description": "<string>",
  "gene_type": "<string>",
  "status": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_gene-parameters)

## 45) Retrieve experimentally observed nematode gene phenotypes {#dnakit.databases.worm_base_get_phenotypes}

- **Purpose:** Retrieve experimentally observed nematode gene phenotypes.
- **Method:** Use `WormBase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.worm_base_get_phenotypes(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `wormbase_id`, `gene_name`, `phenotype_count`, `phenotypes` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.worm_base_get_phenotypes(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "phenotype_count": "<integer>",
  "phenotypes": [
    {
      "phenotype_id": "<string>",
      "phenotype_name": "<string>",
      "evidence_type": "<string>"
    }
  ],
  "not_observed_count": "<integer>",
  "phenotypes_not_observed": [
    {
      "phenotype_id": "<string>",
      "phenotype_name": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_phenotypes-parameters)

## 46) Retrieve Xenopus gene details and cross-references {#dnakit.databases.xenbase_get_gene}

- **Purpose:** Retrieve Xenopus gene details and cross-references.
- **Method:** Use `Xenbase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.xenbase_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `gene_id`, `symbol`, `name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.xenbase_get_gene(
    gene_id='XB-GENE-484286',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "synonyms": [
    "<string>"
  ],
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "so_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.xenbase_get_gene-parameters)

## 47) Search genes in Xenopus tropicalis and Xenopus laevis {#dnakit.databases.xenbase_search_genes}

- **Purpose:** Search genes in Xenopus tropicalis and Xenopus laevis.
- **Method:** Use `Xenbase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.xenbase_search_genes(*, query, species=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `gene_id`, `symbol`, `name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.xenbase_search_genes(
    query='tp53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gene_id": "<string>",
    "symbol": "<string>",
    "name": "<string>",
    "species": "<string>",
    "synonyms": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.xenbase_search_genes-parameters)

## 48) Retrieve zebrafish gene details {#dnakit.databases.zfin_get_gene}

- **Purpose:** Retrieve zebrafish gene details.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.zfin_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `id`, `symbol`, `name`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.zfin_get_gene(
    gene_id='ZFIN:ZDB-GENE-990415-8',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene-parameters)

## 49) Retrieve zebrafish developmental and anatomical expression {#dnakit.databases.zfin_get_gene_expression}

- **Purpose:** Retrieve zebrafish developmental and anatomical expression.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.zfin_get_gene_expression(*, gene_id, limit=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `total_annotations`, `returned`, `annotations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.zfin_get_gene_expression(
    gene_id='ZFIN:ZDB-GENE-980526-166',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_annotations": "<integer>",
  "returned": "<integer>",
  "annotations": [
    {
      "stage": "<string>",
      "location": "<string>",
      "data_provider": "<string>",
      "reference": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene_expression-parameters)

## 50) Retrieve zebrafish gene mutant phenotypes {#dnakit.databases.zfin_get_gene_phenotypes}

- **Purpose:** Retrieve zebrafish gene mutant phenotypes.
- **Method:** Use `AllianceGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.zfin_get_gene_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `gene_symbol`, `gene_id`, `phenotype_statement` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.zfin_get_gene_phenotypes(
    gene_id='ZFIN:ZDB-GENE-990415-8',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "phenotype_statement": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene_phenotypes-parameters)

<!-- END GENERATED API: organisms -->
