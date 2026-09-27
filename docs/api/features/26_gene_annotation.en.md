# Gene names and intervals from genomic DNA

Map genomic DNA to a species reference and return overlapping genes; [FAQ: details](../../faq/tutorial_integrations.md#gene-annotation).

## 1) Annotate genes {#annotate-genes}

- **Purpose:** Return gene names, IDs, strands and intervals in the query sequence.
- **Method:** Align genomic DNA with minimap2 and project existing GFF3 gene annotations onto the query.
- **API:** `dnakit.annotation.annotate_genes(sequence, *, species, assembly=None, sequence_id="query", config=None, ...)`
- **Input:** One genomic DNA sequence and a species name or TaxID; requires external `minimap2`. The first call may download a reference. Supports reference annotation transfer, not de novo prediction or spliced transcripts.
- **Output:** Returns a `GeneAnnotationResult`; the example shows `reference`, `gene_records`, `segments`.
- **Example code:**

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",  # NCBI TaxIDs such as 9606 are also accepted
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(cache_dir=".dnakit/gene_annotation", threads=4),
)
print(result.reference["assembly"])
for gene in result.gene_records:
    print(gene.to_dict())
for segment in result.segments:
    print(segment.to_dict())
```

- **Example output:** The selected assembly, gene records and `gene` / `intergenic` / `unknown` segments; values depend on `query.fasta` and the reference. Missing names remain `null`. Python intervals are 0-based half-open; `to_dict()` and CLI JSON use 1-based closed coordinates.

## 2) Query genes by symbol or GeneID {#search-gene}

- **Purpose:** Retrieve NCBI Gene identifiers, summaries or product reports to resolve gene queries.
- **Method:** Query genes through NCBI Datasets and normalize the selected report into QueryResult.
- **API:** `dnakit.search.gene(query[required], taxon[optional], by[optional], report[optional], gene_types[optional], page_size[optional], page_token[optional], config[optional], verbose[optional])`
- **Input:** A GeneID or gene symbol; symbol queries must include taxon. Optional lookup mode, report type, gene types, pagination and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import SearchConfig, gene

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = gene("BRCA2", taxon="human", report="product", page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
gene NCBI Datasets <record count returned>
```

## 3) Query transcript models for a gene {#search-transcripts}

- **Purpose:** Retrieve Ensembl transcript models and related annotations by gene name or identifier.
- **Method:** Resolve the Ensembl query identifier and retrieve transcript details according to the expand option.
- **API:** `dnakit.search.transcripts(query[required], species[optional], by[optional], expand[optional], config[optional], verbose[optional])`
- **Input:** A gene name or identifier; optional species, lookup mode, expansion and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import transcripts

result = transcripts("BRCA2", species="human", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
transcripts Ensembl <record count returned>
```

## 4) Query annotation features in a genomic interval {#search-annotation}

- **Purpose:** Retrieve genes and other annotation features overlapping a region for interval annotation.
- **Method:** Translate DNAKit interval coordinates, query the Ensembl overlap endpoint and retain coordinate metadata.
- **API:** `dnakit.search.annotation(species[required], region[required], features[optional], strand[optional], biotype[optional], config[optional], verbose[optional])`
- **Input:** A species and a **0-based half-open** interval; optional feature types, strand, biotype and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import annotation

result = annotation("human", "13:32315000-32320000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
annotation Ensembl <record count returned>
```

## 5) Query UCSC track hits in an interval {#search-ucsc-track-data}

- **Purpose:** Retrieve records from a UCSC track that overlap a genomic interval, such as CpG islands.
- **Method:** Query the UCSC API by genome, track and region and normalize track hits into QueryResult.
- **API:** `dnakit.search.ucsc_track_data(genome[required], track[required], region[required], max_items[optional], config[optional], verbose[optional])`
- **Input:** A UCSC genome name, track name and **0-based half-open** interval; optional max_items limit and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import ucsc_track_data

result = ucsc_track_data("hg38", "cpgIslandExt", "chr13:32315000-32320000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
annotation UCSC <record count returned>
```

## 6) Map interval coordinates between reference assemblies {#search-map-coordinates}

- **Purpose:** Map a genomic interval from one reference assembly to another.
- **Method:** Call the Ensembl mapping endpoint and retain source assembly, target assembly and coordinate metadata.
- **API:** `dnakit.search.map_coordinates(species[required], region[required], source_assembly[required], target_assembly[required], strand[optional], config[optional], verbose[optional])`
- **Input:** A species, a **0-based half-open** interval, source assembly and target assembly; optional strand and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import map_coordinates

result = map_coordinates(
    "human",
    "1:100000-101000",
    source_assembly="GRCh38",
    target_assembly="GRCh37",
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
coordinate_mapping Ensembl <record count returned>
```

<!-- BEGIN GENERATED API: gene-annotation -->

<span id="gene-annotation"></span>

## 7) Retrieve a cross-species gene family tree {#dnakit.databases.ensembl_compara_get_gene_tree}

- **Purpose:** Retrieve a cross-species gene family tree.
- **Method:** Use `EnsemblCompara` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_compara_get_gene_tree(*, gene, species=None, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a mapping with `tree_id`, `newick`, `rooted`, `members` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_gene_tree(
    gene='BRCA1',
    species='human',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "tree_id": "<string>",
  "newick": "<string>",
  "rooted": "<boolean>",
  "members": [
    "<object>"
  ],
  "total_members": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_gene_tree-parameters)

## 8) Retrieve orthologous genes across species {#dnakit.databases.ensembl_compara_get_orthologues}

- **Purpose:** Retrieve orthologous genes across species.
- **Method:** Use `EnsemblCompara` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_compara_get_orthologues(*, gene, species=None, target_species=None, target_taxon=None, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a list of records; each record includes `source_gene`, `target_gene`, `target_protein`, `target_species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_orthologues(
    gene='BRCA1',
    target_species='mouse',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "source_gene": "<string>",
    "target_gene": "<string>",
    "target_protein": "<string>",
    "target_species": "<string>",
    "homology_type": "<string>",
    "taxonomy_level": "<string>",
    "method": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_orthologues-parameters)

## 9) Retrieve paralogous genes {#dnakit.databases.ensembl_compara_get_paralogues}

- **Purpose:** Retrieve paralogous genes.
- **Method:** Use `EnsemblCompara` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_compara_get_paralogues(*, gene, species=None, verbose=True)`
- **Input:** Provide `gene`.
- **Output:** Returns a list of records; each record includes `source_gene`, `paralogue_gene`, `paralogue_protein`, `species` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_paralogues(
    gene='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "source_gene": "<string>",
    "paralogue_gene": "<string>",
    "paralogue_protein": "<string>",
    "species": "<string>",
    "paralogy_type": "<string>",
    "taxonomy_level": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_paralogues-parameters)

## 10) Retrieve gene homology by symbol {#dnakit.databases.ensembl_get_homology}

- **Purpose:** Retrieve gene homology by symbol.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_homology(*, species, symbol, target_species=None, target_taxon=None, type=None, sequence=None, aligned=None, verbose=True)`
- **Input:** Provide `species`, `symbol`.
- **Output:** Returns a list of records; each record includes `homologies` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_homology(
    species='human',
    symbol='BRCA2',
    target_species='mouse',
    type='orthologues',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "homologies": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_homology-parameters)

## 11) Retrieve annotations overlapping a genomic interval {#dnakit.databases.ensembl_get_overlap_features}

- **Purpose:** Retrieve annotations overlapping a genomic interval.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_overlap_features(*, species, region, feature=None, biotype=None, verbose=True)`
- **Input:** Provide `species`, `region`.
- **Output:** Returns a list of records; each record includes `id`, `external_name`, `biotype`, `start` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_overlap_features(
    species='human',
    region='7:140424943-140624564',
    feature='gene',
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
    "external_name": "<string>",
    "biotype": "<string>",
    "start": "<integer>",
    "end": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_overlap_features-parameters)

## 12) Retrieve regulatory features in a genomic interval {#dnakit.databases.ensembl_get_regulatory_features}

- **Purpose:** Retrieve regulatory features in a genomic interval.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_regulatory_features(*, region, species=None, feature, verbose=True)`
- **Input:** Provide `region`, `feature`.
- **Output:** Returns a list of records; each record includes `id`, `feature_type`, `start`, `end` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_regulatory_features(
    region='1:1000000-2000000',
    species='human',
    feature='regulatory',
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
    "feature_type": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "strand": "<integer>",
    "seq_region_name": "<string>",
    "description": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_regulatory_features-parameters)

## 13) Look up Ensembl genes by identifier or symbol {#dnakit.databases.ensembl_lookup_gene}

- **Purpose:** Look up Ensembl genes by identifier or symbol.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_lookup_gene(*, gene_id=None, species=None, id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `id`, `display_name`, `description`, `biotype` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_lookup_gene(
    gene_id='ENSG00000012048',
    species='homo_sapiens',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "id": "<string>",
  "display_name": "<string>",
  "description": "<string>",
  "biotype": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "strand": "<integer>",
  "seq_region_name": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_lookup_gene-parameters)

## 14) Retrieve HGNC annotations by official gene symbol {#dnakit.databases.hgnc_fetch_gene_by_symbol}

- **Purpose:** Retrieve HGNC annotations by official gene symbol.
- **Method:** Use `HGNC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hgnc_fetch_gene_by_symbol(*, symbol, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns a mapping with `hgnc_id`, `symbol`, `name`, `locus_type` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hgnc_fetch_gene_by_symbol(
    symbol='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "hgnc_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "locus_type": "<string>",
  "locus_group": "<string>",
  "status": "<string>",
  "location": "<string>",
  "entrez_id": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.hgnc_fetch_gene_by_symbol-parameters)

## 15) Batch-query gene annotations and identifier mappings {#dnakit.databases.my_gene_batch_query}

- **Purpose:** Batch-query gene annotations and identifier mappings.
- **Method:** Use `MyGene` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_gene_batch_query(*, gene_ids, species=None, fields=None, verbose=True)`
- **Input:** Provide `gene_ids`.
- **Output:** Returns a mapping with `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_gene_batch_query(
    gene_ids=['CDK2', 'TP53', 'BRCA1'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "results": [
    {
      "_id": "<string>",
      "query": "<string>",
      "symbol": "<string>",
      "name": "<string>",
      "entrezgene": "<integer>",
      "notfound": "<boolean>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_batch_query-parameters)

## 16) Retrieve comprehensive annotations for a gene identifier {#dnakit.databases.my_gene_get_gene_annotation}

- **Purpose:** Retrieve comprehensive annotations for a gene identifier.
- **Method:** Use `MyGene` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_gene_get_gene_annotation(*, gene_id, fields=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `_id`, `symbol`, `name`, `entrezgene` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_gene_get_gene_annotation(
    gene_id='1017',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "entrezgene": "<integer>",
  "ensembl": {
    "gene": "<string>",
    "transcript": [
      "<string>"
    ],
    "protein": [
      "<string>"
    ]
  },
  "summary": "<string>",
  "go": {
    "BP": [
      "<value>"
    ],
    "CC": [
      "<value>"
    ],
    "MF": [
      "<value>"
    ]
  },
  "pathway": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_get_gene_annotation-parameters)

## 17) Search genes by keyword or identifier {#dnakit.databases.my_gene_query_genes}

- **Purpose:** Search genes by keyword or identifier.
- **Method:** Use `MyGene` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_gene_query_genes(*, query=None, species=None, fields=None, size=None, q=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `took`, `total`, `max_score`, `hits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_gene_query_genes(
    query='CDK2',
    species='human',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "took": "<integer>",
  "total": "<integer>",
  "max_score": "<number>",
  "hits": [
    {
      "_id": "<string>",
      "_score": "<number>",
      "symbol": "<string>",
      "name": "<string>",
      "entrezgene": "<integer>",
      "ensembl": "<object>",
      "summary": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_query_genes-parameters)

## 18) Retrieve gene details by NCBI Gene ID {#dnakit.databases.ncbidatasets_get_gene}

- **Purpose:** Retrieve gene details by NCBI Gene ID.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_gene(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `gene_id`, `symbol`, `description`, `tax_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_gene(
    gene_id='7157',
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
  "description": "<string>",
  "tax_id": "<string>",
  "taxname": "<string>",
  "common_name": "<string>",
  "type": "<string>",
  "chromosomes": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_gene-parameters)

## 19) Look up NCBI genes by symbol and organism {#dnakit.databases.ncbidatasets_get_gene_by_symbol}

- **Purpose:** Look up NCBI genes by symbol and organism.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_gene_by_symbol(*, symbol, taxon=None, organism=None, species=None, verbose=True)`
- **Input:** Provide `symbol`.
- **Output:** Returns a list of records; each record includes `gene_id`, `symbol`, `description`, `tax_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_gene_by_symbol(
    symbol='BRCA1',
    taxon='human',
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
    "description": "<string>",
    "tax_id": "<string>",
    "taxname": "<string>",
    "type": "<string>",
    "chromosomes": [
      "<string>"
    ],
    "swiss_prot_accessions": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_gene_by_symbol-parameters)

## 20) Retrieve orthologs for an NCBI Gene ID {#dnakit.databases.ncbidatasets_get_orthologs}

- **Purpose:** Retrieve orthologs for an NCBI Gene ID.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_orthologs(*, gene_id, page_size=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `gene_id`, `symbol`, `description`, `tax_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_orthologs(
    gene_id='7157',
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
    "gene_id": "<string>",
    "symbol": "<string>",
    "description": "<string>",
    "tax_id": "<string>",
    "taxname": "<string>",
    "common_name": "<string>",
    "type": "<string>",
    "chromosomes": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_orthologs-parameters)

## 21) Retrieve NCBI gene summaries and functional information {#dnakit.databases.ncbigene_get_summary}

- **Purpose:** Retrieve NCBI gene summaries and functional information.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbigene_get_summary(*, id, verbose=True)`
- **Input:** Provide `id`.
- **Output:** Returns a mapping with `result` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbigene_get_summary(
    id='7157',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "result": {
    "uids": [
      "<string>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbigene_get_summary-parameters)

## 22) Search NCBI Gene IDs by organism and keyword {#dnakit.databases.ncbigene_search}

- **Purpose:** Search NCBI Gene IDs by organism and keyword.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbigene_search(*, term, retmax=None, verbose=True)`
- **Input:** Provide `term`.
- **Output:** Returns a mapping with `header`, `esearchresult` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbigene_search(
    term='TP53[Symbol] AND Homo sapiens[Organism]',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "header": {
    "type": "<string>",
    "version": "<string>"
  },
  "esearchresult": {
    "count": "<string>",
    "retmax": "<string>",
    "retstart": "<string>",
    "idlist": [
      "<string>"
    ],
    "querytranslation": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbigene_search-parameters)

## 23) Retrieve MANE transcripts and RefSeq mappings {#dnakit.databases.tark_get_mane_transcripts}

- **Purpose:** Retrieve MANE transcripts and RefSeq mappings.
- **Method:** Use `TarkManeTranscripts` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.tark_get_mane_transcripts(*, gene=None, ensembl_id=None, refseq_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `gene`, `mane_type`, `ensembl_transcript`, `refseq_transcript` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.tark_get_mane_transcripts(
    gene='BRCA2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gene": "<string>",
    "mane_type": "<string>",
    "ensembl_transcript": "<string>",
    "refseq_transcript": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.tark_get_mane_transcripts-parameters)

## 24) Retrieve archived Ensembl transcript records {#dnakit.databases.tark_get_transcript}

- **Purpose:** Retrieve archived Ensembl transcript records.
- **Method:** Use `TarkTranscript` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.tark_get_transcript(*, stable_id, verbose=True)`
- **Input:** Provide `stable_id`.
- **Output:** Returns a list for “Retrieve archived Ensembl transcript records”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.tark_get_transcript(
    stable_id='ENST00000380152',
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

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.tark_get_transcript-parameters)

## 25) Retrieve gene transcripts and MANE annotations {#dnakit.databases.variant_validator_gene2transcripts}

- **Purpose:** Retrieve gene transcripts and MANE annotations.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.variant_validator_gene2transcripts(*, gene_symbol, gene=None, gene_name=None, transcript_set=None, genome_build=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a list of records; each record includes `current_name`, `current_symbol`, `hgnc`, `previous_symbol` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.variant_validator_gene2transcripts(
    gene_symbol='TP53',
    transcript_set='mane',
    genome_build='GRCh38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "current_name": "<string>",
    "current_symbol": "<string>",
    "hgnc": "<string>",
    "previous_symbol": "<string>",
    "requested_symbol": "<string>",
    "transcripts": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/26_gene_annotation.md#dnakit.databases.variant_validator_gene2transcripts-parameters)

<!-- END GENERATED API: gene-annotation -->
