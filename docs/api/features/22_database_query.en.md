# Data query

Query public databases such as NCBI, Ensembl, ENA, ENCODE, and UCSC using a unified interface to obtain DNA sequences and related metadata.

<span id="1"></span>**Providers and Entrances**

- **Function:** Aggregate public databases accessible to DNAKit, and unify timeouts, retries, caches, rate limits, and return objects so that query results from different sources have consistent provenance and error handling.
- **API:** `dnakit.search.SearchConfig(timeout[optional], max_response_bytes[optional], max_records[optional], api_key[optional], email[optional], tool[optional])`, `dnakit.search.QueryResult`.
- **Input:** Optional timeout, maximum number of response bytes, maximum number of records, and NCBI contact information; specific query conditions are provided by the subsequent query API.

| Data source | Main entrance | Purpose |
| ------------------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| NCBI Datasets | `taxonomy()`, `assembly()`, `gene()`, `virus()` | Classification, assembly, gene and virus queries |
| NCBI Entrez | `entrez()` and project/sample/variant/expression/literature packages | accession, BioProject, BioSample, ClinVar, GEO, PubMed, version |
| NCBI BLAST | `submit_blast()`, `identify()`, `novelty()` | Asynchronous similarity, provenance clues, and remote novelty checks |
| Ensembl REST | `sequence()`, `transcripts()`, `annotation()`, `variant()`, `homology()` | Coordinate sequences, transcripts, annotations, VEPs, homologous and comparative genomes |
| ENA Portal | `reads()`, `ena_project()`, `ena_sample()` | Study/Sample/Experiment/Run/Analysis and public documentation |
| ENCODE Portal | `encode_search()` | Epigenomic experiments and documentation |
| UCSC REST/Downloads | `ucsc_*()` | Chromosomes, sequences, CpG/repeat/regulatory/conserved tracks and file directories |

These entries are implemented based on [NCBI Datasets REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/), [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/), [Ensembl REST](https://rest.ensembl.org/), [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api), [ENCODE REST API](https://www.encodeproject.org/help/rest-api/) and [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html).

## 1) Basic query <span id="2"></span>

- **Function:** Query NCBI Taxonomy, Genome and Gene by species, TaxID, assembly accession, GeneID or gene symbol, and return standardized identifiers, abstracts and source links for parsing subsequent download targets.
- **API:** `dnakit.search.taxonomy(query[required], report[optional], children[optional], include_lineage[optional], ranks[optional], page_size[optional], page_token[optional], config[optional])`, `dnakit.search.assembly(query[required], by[optional], reference_only[optional], assembly_source[optional], current_only[optional], page_size[optional], page_token[optional], config[optional])`, `dnakit.search.gene(query[required], taxon[optional], by[optional], report[optional], gene_types[optional], page_size[optional], page_token[optional], config[optional])`.
- **Input:** Required species name, TaxID, assembly accession, GeneID, or gene name; optional report type, species range, paging parameters, and unification `SearchConfig`. Gene symbol queries must also provide `taxon`.
- **Sample code:**

```python
from dnakit.search import SearchConfig, assembly, gene, taxonomy

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)

taxon = taxonomy("human", page_size=5, config=config)
genomes = assembly("human", reference_only=True, page_size=5, config=config)
brca2 = gene("BRCA2", taxon="human", report="product", page_size=5, config=config)

for result in (taxon, genomes, brca2):
    print(result.query_type, result.provider, len(result.records))
```

- **Example results:** The following only shows the stable return structure; the number of real-time records is determined by the database at the time of query.

```text
taxonomy NCBI Datasets <record count returned>
assembly NCBI Datasets <record count returned>
gene NCBI Datasets <record count returned>
```

## 2) Coordinates, transcripts and area annotations <span id="3"></span>

- **Function:** Query regional sequences, transcripts, gene annotations and UCSC track hits by assembly and genome coordinates, and return structured records of bound coordinate systems for local region annotation.
- **API:** `dnakit.search.sequence(species[required], region[required], strand[optional], upstream[optional], downstream[optional], mask[optional], progress[optional], config[optional])`, `dnakit.search.transcripts(query[required], species[optional], by[optional], expand[optional], config[optional])`, `dnakit.search.annotation(species[required], region[required], features[optional], strand[optional], biotype[optional], config[optional])`, `dnakit.search.ucsc_track_data(genome[required], track[required], region[required], max_items[optional], config[optional])`.
- **Input:** Required species or UCSC genome name, `chromosome:start-end` region, and query target; optional chain direction, upstream and downstream extensions, annotation type, track name, progress callback, and resource configuration. DNAKit input coordinates are unified into **0-based half-open intervals**.
- **Sample code:**

```python
from dnakit.search import annotation, sequence, transcripts, ucsc_track_data

region = sequence(
    "human",
    ("1:0-100", "X:1000-1100"),
    strand=-1,
    upstream=20,
)

gene_model = transcripts("BRCA2", species="human")
features = annotation("human", "13:32315000-32320000")
cpg = ucsc_track_data("hg38", "cpgIslandExt", "chr13:32315000-32320000")

for result in (region, gene_model, features, cpg):
    print(result.query_type, result.provider, len(result.records))
```

- **Example results:** The following only shows the stable return structure; the number of real-time records is determined by the database at the time of query.

```text
sequence Ensembl <record count returned>
transcripts Ensembl <record count returned>
annotation Ensembl <record count returned>
annotation UCSC <record count returned>
```

Batch Ensembl sequence query supports progress callback. The following uses the existing Rich dependency of the project to display the progress bar:

```python
from rich.progress import Progress

from dnakit.search import sequence

regions = ("1:0-100", "2:0-100", "X:0-100")
with Progress() as bar:
    task = bar.add_task("Query Ensembl", total=len(regions))

    def update(event):
        bar.update(task, completed=event.completed, description=event.item)

    result = sequence("human", regions, progress=update)
```

## 3) Sequencing, expression, regulation and comparison of genomes <span id="4"></span>

- **Function:** Query ENA sequencing runs, GEO expression data, ENCODE regulatory records, cross-assembly coordinate mapping and comparative genome alignment, and return filterable metadata and file links.
- **API:** `dnakit.search.reads(query[required], fields[optional], result[optional], limit[optional], offset[optional], config[optional])`, `dnakit.search.expression(term[required], retmax[optional], config[optional])`, `dnakit.search.encode_search(object_type[optional], search_term[optional], filters[optional], limit[optional], config[optional])`, `dnakit.search.map_coordinates(species[required], region[required], source_assembly[required], target_assembly[required], strand[optional], config[optional])`, `dnakit.search.comparative_alignment(species[required], region[required], strand[optional], method[optional], species_set_group[optional], display_species[optional], aligned[optional], mask[optional], config[optional])`.
- **Input:** The query expression or genomic region corresponding to the provider is required; the coordinate map must also provide the source assembly and the target assembly. Optional return fields, filters, record caps, comparison methods, and unification `SearchConfig`.
- **Sample code:**

```python
from dnakit.search import (
    comparative_alignment,
    encode_search,
    expression,
    map_coordinates,
    reads,
)

runs = reads(
    'library_strategy="RNA-Seq" AND instrument_platform="ILLUMINA"',
    limit=20,
)
geo = expression("single cell RNA sequencing[All Fields]", retmax=20)
encode_files = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=20,
)
mapped = map_coordinates(
    "human",
    "1:100000-101000",
    source_assembly="GRCh38",
    target_assembly="GRCh37",
)
alignment = comparative_alignment("human", "1:100000-101000")

for result in (runs, geo, encode_files, mapped, alignment):
    print(result.query_type, result.provider, len(result.records))
```

- **Example results:** The following only shows the stable return structure; the number of real-time records is determined by the database at the time of query.

```text
reads ENA <record count returned>
expression NCBI Entrez <record count returned>
regulation ENCODE <record count returned>
coordinate_mapping Ensembl <record count returned>
comparative_alignment Ensembl <record count returned>
```
