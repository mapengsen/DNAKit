# Sequence, assembly and sequencing queries

Query sequences, reference assemblies and sequencing projects.

## 1) Query reference assemblies and assembly accessions {#search-assembly}

- **Purpose:** Retrieve assembly records by species or assembly accession to select a reference genome.
- **Method:** Query NCBI Genome with filters for reference assemblies, assembly source and current versions.
- **API:** `dnakit.search.assembly(query[required], by[optional], reference_only[optional], assembly_source[optional], current_only[optional], page_size[optional], page_token[optional], config[optional], verbose[optional])`
- **Input:** A species name, TaxID or assembly accession; optional lookup mode, reference filters, pagination and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import SearchConfig, assembly

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = assembly("human", reference_only=True, page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
assembly NCBI Datasets <record count returned>
```

## 2) Query sequences by genomic interval {#search-sequence}

- **Purpose:** Retrieve DNA from genomic intervals with optional strand selection and flanking sequence.
- **Method:** Translate 0-based half-open intervals into Ensembl coordinates, retrieve sequences and retain coordinate metadata.
- **API:** `dnakit.search.sequence(species[required], region[required], strand[optional], upstream[optional], downstream[optional], mask[optional], progress[optional], config[optional], verbose[optional])`
- **Input:** A species and one or more chromosome:start-end intervals in **0-based half-open coordinates**; optional strand, flanks, masking, progress callback and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import sequence

result = sequence(
    "human",
    ("1:0-100", "X:1000-1100"),
    strand=-1,
    upstream=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
sequence Ensembl <record count returned>
```

## 3) Query ENA sequencing runs and file metadata {#search-reads}

- **Purpose:** Find ENA sequencing runs by library strategy, platform and other conditions, returning metadata and file information.
- **Method:** Submit an expression to ENA Portal and return QueryResult with the requested fields, offset and record limit.
- **API:** `dnakit.search.reads(query[required], fields[optional], result[optional], limit[optional], offset[optional], config[optional], verbose[optional])`
- **Input:** An ENA query expression; optional fields, result type, limit, offset and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import reads

result = reads(
    'library_strategy="RNA-Seq" AND instrument_platform="ILLUMINA"',
    limit=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
reads ENA <record count returned>
```

## 4) Query comparative genomic alignment blocks {#search-comparative-alignment}

- **Purpose:** Retrieve cross-species alignment blocks for a genomic region to compare genome sequences.
- **Method:** Retrieve regional alignments from Ensembl Compara using the selected alignment method and species set.
- **API:** `dnakit.search.comparative_alignment(species[required], region[required], strand[optional], method[optional], species_set_group[optional], display_species[optional], aligned[optional], mask[optional], config[optional], verbose[optional])`
- **Input:** A species and a **0-based half-open** genomic interval; optional strand, method, species set, displayed species, alignment and masking options. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import comparative_alignment

result = comparative_alignment("human", "1:100000-101000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
comparative_alignment Ensembl <record count returned>
```

<!-- BEGIN GENERATED API: sequences -->

<span id="sequences"></span>

## 5) Search ENA biological samples by taxonomy or keyword {#dnakit.databases.enaportal_search_samples}

- **Purpose:** Search ENA biological samples by taxonomy or keyword.
- **Method:** Use `ENAPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.enaportal_search_samples(*, query, limit=None, fields=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `sample_accession`, `sample_alias`, `description`, `tax_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.enaportal_search_samples(
    query='tax_tree(562)',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "sample_accession": "<string>",
    "sample_alias": "<string>",
    "description": "<string>",
    "tax_id": "<string>",
    "scientific_name": "<string>",
    "first_public": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.enaportal_search_samples-parameters)

## 6) Search ENA sequencing studies {#dnakit.databases.enaportal_search_studies}

- **Purpose:** Search ENA sequencing studies.
- **Method:** Use `ENAPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.enaportal_search_studies(*, query, limit=None, fields=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `study_accession`, `study_title`, `center_name`, `first_public` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.enaportal_search_studies(
    query='description="SARS-CoV-2"',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "study_accession": "<string>",
    "study_title": "<string>",
    "center_name": "<string>",
    "first_public": "<string>",
    "description": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.enaportal_search_studies-parameters)

## 7) Retrieve gene, transcript or genomic interval sequences {#dnakit.databases.ensembl_get_sequence}

- **Purpose:** Retrieve gene, transcript or genomic interval sequences.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_get_sequence(*, id, type=None, species=None, multiple_sequences=None, verbose=True)`
- **Input:** Provide `id`.
- **Output:** Returns a list of records; each record includes `id`, `seq`, `desc`, `molecule` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_get_sequence(
    id='ENSG00000012048',
    type='genomic',
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
    "seq": "<string>",
    "desc": "<string>",
    "molecule": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_get_sequence-parameters)

## 8) Retrieve a sequence by Ensembl identifier {#dnakit.databases.ensembl_seq_get_id_sequence}

- **Purpose:** Retrieve a sequence by Ensembl identifier.
- **Method:** Use `EnsemblSequence` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_seq_get_id_sequence(*, ensembl_id, type=None, verbose=True)`
- **Input:** Provide `ensembl_id`.
- **Output:** Returns a mapping with `ensembl_id`, `molecule`, `description`, `sequence_length` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_seq_get_id_sequence(
    ensembl_id='ENSP00000269305',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "ensembl_id": "<string>",
  "molecule": "<string>",
  "description": "<string>",
  "sequence_length": "<integer>",
  "sequence": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_seq_get_id_sequence-parameters)

## 9) Retrieve Ensembl reference sequence by genomic coordinates {#dnakit.databases.ensembl_seq_get_region_sequence}

- **Purpose:** Retrieve Ensembl reference sequence by genomic coordinates.
- **Method:** Use `EnsemblSequence` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_seq_get_region_sequence(*, region, species=None, verbose=True)`
- **Input:** Provide `region`.
- **Output:** Returns a mapping with `species`, `region`, `molecule`, `description` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_seq_get_region_sequence(
    region='17:7668421-7668520',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "species": "<string>",
  "region": "<string>",
  "molecule": "<string>",
  "description": "<string>",
  "sequence_length": "<integer>",
  "sequence": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_seq_get_region_sequence-parameters)

## 10) Convert GenBank UIDs to sequence accessions {#dnakit.databases.ncbi_fetch_accessions}

- **Purpose:** Convert GenBank UIDs to sequence accessions.
- **Method:** Use `NCBINucleotideSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbi_fetch_accessions(*, operation=None, uids, verbose=True)`
- **Input:** Provide `uids`.
- **Output:** Returns a list for “Convert GenBank UIDs to sequence accessions”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbi_fetch_accessions(
    uids=['545778205'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<string>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_fetch_accessions-parameters)

## 11) Retrieve NCBI nucleotide records by accession {#dnakit.databases.ncbi_get_sequence}

- **Purpose:** Retrieve NCBI nucleotide records by accession.
- **Method:** Use `NCBINucleotideSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbi_get_sequence(*, operation=None, accession, format=None, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a string.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbi_get_sequence(
    accession='U00096',
    format='fasta',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
"<string>"
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_get_sequence-parameters)

## 12) Search nucleotide records by organism, gene or strain {#dnakit.databases.ncbi_search_nucleotide}

- **Purpose:** Search nucleotide records by organism, gene or strain.
- **Method:** Use `NCBINucleotideSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbi_search_nucleotide(*, operation=None, organism=None, gene=None, strain=None, keywords=None, seq_type=None, query=None, limit=None, sort=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `uids`, `count`, `returned`, `search_term` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbi_search_nucleotide(
    organism='Escherichia coli',
    strain='K-12',
    keywords='complete genome',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "uids": [
    "<string>"
  ],
  "count": "<integer>",
  "returned": "<integer>",
  "search_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_search_nucleotide-parameters)

## 13) Search SRA runs by study and experimental conditions {#dnakit.databases.ncbi_sra_search_runs}

- **Purpose:** Search SRA runs by study and experimental conditions.
- **Method:** Use `NCBISRA` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbi_sra_search_runs(*, operation=None, study=None, organism=None, strategy=None, platform=None, source=None, query=None, limit=None, sort=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `uids`, `count`, `returned`, `search_term` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbi_sra_search_runs(
    operation='search',
    organism='Homo sapiens',
    strategy='RNA-Seq',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "uids": [
    "<string>"
  ],
  "count": "<integer>",
  "returned": "<integer>",
  "search_term": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_sra_search_runs-parameters)

## 14) Retrieve genome assembly versions and metadata {#dnakit.databases.ncbidatasets_get_genome_assembly}

- **Purpose:** Retrieve genome assembly versions and metadata.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_genome_assembly(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `accession`, `paired_accession`, `source_database`, `organism_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_genome_assembly(
    accession='GCF_000005845.2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "accession": "<string>",
  "paired_accession": "<string>",
  "source_database": "<string>",
  "organism_name": "<string>",
  "tax_id": "<integer>",
  "strain": "<string>",
  "assembly_name": "<string>",
  "assembly_level": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_genome_assembly-parameters)

## 15) Retrieve chromosome and scaffold sequence reports {#dnakit.databases.ncbidatasets_get_sequence_reports}

- **Purpose:** Retrieve chromosome and scaffold sequence reports.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_sequence_reports(*, accession, page_size=None, page_token=None, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a list of records; each record includes `chr_name`, `sequence_name`, `role`, `location_type` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_sequence_reports(
    accession='GCF_000005845.2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "chr_name": "<string>",
    "sequence_name": "<string>",
    "role": "<string>",
    "location_type": "<string>",
    "refseq_accession": "<string>",
    "genbank_accession": "<string>",
    "length": "<integer>",
    "gc_percent": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_sequence_reports-parameters)

## 16) Retrieve taxonomy lineages and counts by taxon ID {#dnakit.databases.ncbidatasets_get_taxonomy}

- **Purpose:** Retrieve taxonomy lineages and counts by taxon ID.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_get_taxonomy(*, tax_id, verbose=True)`
- **Input:** Provide `tax_id`.
- **Output:** Returns a mapping with `tax_id`, `organism_name`, `genbank_common_name`, `rank` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_taxonomy(
    tax_id='9606',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "tax_id": "<integer>",
  "organism_name": "<string>",
  "genbank_common_name": "<string>",
  "rank": "<string>",
  "blast_name": "<string>",
  "lineage": [
    "<integer>"
  ],
  "lineage_names": [
    {
      "tax_id": "<integer>",
      "organism_name": "<string>",
      "rank": "<string>"
    }
  ],
  "children": [
    "<integer>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_taxonomy-parameters)

## 17) List genome assemblies for a taxon {#dnakit.databases.ncbidatasets_list_genomes_by_taxon}

- **Purpose:** List genome assemblies for a taxon.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_list_genomes_by_taxon(*, taxon, limit=None, reference_only=None, verbose=True)`
- **Input:** Provide `taxon`.
- **Output:** Returns a list of records; each record includes `accession`, `paired_accession`, `source_database`, `organism_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_list_genomes_by_taxon(
    taxon='Saccharomyces cerevisiae',
    limit=5,
    reference_only=True,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "accession": "<string>",
    "paired_accession": "<string>",
    "source_database": "<string>",
    "organism_name": "<string>",
    "tax_id": "<integer>",
    "strain": "<string>",
    "assembly_name": "<string>",
    "assembly_level": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_list_genomes_by_taxon-parameters)

## 18) Find taxonomy suggestions from a partial organism name {#dnakit.databases.ncbidatasets_suggest_taxonomy}

- **Purpose:** Find taxonomy suggestions from a partial organism name.
- **Method:** Use `NCBIDatasets` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbidatasets_suggest_taxonomy(*, query, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `scientific_name`, `tax_id`, `common_name`, `rank` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbidatasets_suggest_taxonomy(
    query='drosophila',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "scientific_name": "<string>",
    "tax_id": "<string>",
    "common_name": "<string>",
    "rank": "<string>",
    "group_name": "<string>",
    "matched_term": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_suggest_taxonomy-parameters)

## 19) Retrieve SRA experiment metadata by accession {#dnakit.databases.sra_get_experiment}

- **Purpose:** Retrieve SRA experiment metadata by accession.
- **Method:** Use `SRA` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sra_get_experiment(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sra_get_experiment(
    accession='SRX515811',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "uid": "<string>",
    "title": "<string>",
    "organism": "<string>",
    "taxid": "<string>",
    "platform": "<string>",
    "library_strategy": "<string>",
    "library_source": "<string>",
    "library_selection": "<string>"
  },
  "metadata": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.sra_get_experiment-parameters)

## 20) Search SRA sequencing experiments {#dnakit.databases.sra_search_experiments}

- **Purpose:** Search SRA sequencing experiments.
- **Method:** Use `SRA` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.sra_search_experiments(*, query=None, organism=None, library_strategy=None, platform=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sra_search_experiments(
    query='BRCA1',
    organism='Homo sapiens',
    library_strategy='RNA-Seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "total": "<integer>",
    "returned": "<integer>",
    "query_used": "<string>",
    "experiments": [
      "<object>"
    ]
  },
  "metadata": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.databases.sra_search_experiments-parameters)

## 21) Resolve sequence UIDs to versioned accessions in batches {#dnakit.search.ncbi.accessions}

- **Purpose:** Resolve sequence UIDs to versioned accessions in batches.
- **Method:** Resolve UIDs to versioned accessions using NCBI ESummary.
- **API：** `dnakit.search.ncbi.accessions(uids, *, database='nuccore', config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`
- **Input:** Provide `uids`.
- **Output:** tuple[str, ...]; empty input gives an empty tuple, missing UIDs or versions raise QueryError.
- **Example code:**

```python
import dnakit

result = dnakit.search.ncbi.accessions(
    uids=['556503834'],
    verbose=False,
)
print(result)
```

- **Example output:** tuple[str, ...]; empty input gives an empty tuple, missing UIDs or versions raise QueryError.

[Full parameters and usage notes](../../faq/api/22_database_query_sequences.md#dnakit.search.ncbi.accessions-parameters)

<!-- END GENERATED API: sequences -->
