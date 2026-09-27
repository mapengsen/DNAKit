# Universal search

Perform exact, subsequence, approximate, reverse complement and nearest neighbor searches in single or multiple DNA sequences.

## 1) Exact search {#1-sim-001-exact-search}

- **Purpose:** Search the record set for records that are identical to the query sequence character by character, and return the target index, start and end coordinates, and chain direction, used for member checking, exact copy positioning, and small-scale reference queries.
- **Method:** Compare the complete query and target symbol strings character by character.
- **API:** `dnakit.similarity.exact_search(query[required], targets[required], reverse_complement[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.
- **Input:** Required query sequence and target sequence/`DNASet`; optional reverse complementary search, chain merging and resource limit.
- **Output:** Returns `SearchResult` with positions and matched records in `matches`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import exact_search

result = exact_search(
    DNASequence("AC"),
    [DNASequence("AC"), DNASequence("TACG")],
)
print([(hit.target_index, hit.start, hit.end) for hit in result.matches])
```

- **Example output:**

```text
[(0, 0, 2)]
```


## 2) Subsequence search {#2-sim-002-subsequence-search}

- **Purpose:** Find all precise occurrence positions of short query sequences in the target long sequence, and return the start and end coordinates and matching directions, which can be used to locate motifs, primers or known fragments.
- **Method:** Find exact substrings, optionally retaining overlapping matches.
- **API:** `dnakit.similarity.subsequence_search(query[required], target[required], strand[optional], overlapping[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.
- **Input:** Required query and target; optional `strand`, `overlapping`, chain merging and resource limit.
- **Output:** Returns `SearchResult` with positions and matched records in `matches`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import subsequence_search

result = subsequence_search(DNASequence("ANA", alphabet="iupac"),
                            DNASequence("ANANA", alphabet="iupac"))
print([(hit.start, hit.end) for hit in result.matches])
```

- **Example output:**

```text
[(0, 3), (2, 5)]
```


## 3) Approximate matching {#3-sim-004-approximate-matching}

- **Purpose:** Find approximate hits in the target sequence that allow a specified number of substitutions, insertions, or deletions, and return coordinates, distances, and matching fragments for searches that tolerate mutations or sequencing errors.
- **Method:** Select matches within the configured weighted substitution/insertion/deletion distance.
- **API:** `dnakit.similarity.approximate_search(query[required], targets[required], max_distance[required], substitution_cost[optional], insertion_cost[optional], deletion_cost[optional], reverse_complement[optional], max_targets[optional], max_matches[optional], max_cells[optional])`.
- **Input:** Required query, target and `max_distance`; optional three types of editing costs, reverse complementation and resource limit.
- **Output:** Returns an `ApproximateSearchResult`; the example shows `matches`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import approximate_search

result = approximate_search(DNASequence("ACG"), DNASequence("TTACGATG"),
                            max_distance=1)
print(any(hit.start == 2 and hit.distance == 0 for hit in result.matches))
```

- **Example output:**

```text
True
```

## 4) Reverse complementary search {#4-sim-005-reverse-complementary-search}

- **Purpose:** Search using the query sequence and its reverse complement simultaneously, returning the coordinates and strand direction of each hit to avoid missing the same fragment located in the reverse strand.
- **Method:** Search the query and its reverse complement and report strand.
- **API:** `dnakit.similarity.reverse_complement_search(query[required], target[required], overlapping[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.
- **Input:** Required query and target; optional whether to allow overlap and whether to merge palindromic duplicate hits.
- **Output:** Returns `SearchResult` with positions and matched records in `matches`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import reverse_complement_search

result = reverse_complement_search(DNASequence("ATG"), DNASequence("GGCATCC"))
print([(hit.start, hit.end, hit.strand.value) for hit in result.matches])
```

- **Example output:**

```text
[(2, 5, 'reverse')]
```

## 5) Nearest neighbor search {#5-sim-014-nearest-neighbor-search}

- **Purpose:** Use k-mer Sketch of the query sequence to filter the Top-k records with the highest similarity in the existing index, and return the record ID, ranking and approximate similarity, which is suitable for fast candidate retrieval.
- **Method:** Compare k-mer sketches and rank the top-k approximate similarities.
- **API:** `dnakit.similarity.build_sketch_index(records[required], k[optional], num_hashes[optional], canonical[optional], seed[optional], max_records[optional])`, `dnakit.similarity.nearest_neighbors(query[required], index[required], top_k[optional], min_similarity[optional])`.
- **Input:** Required query and `SketchIndex`; optional `top_k`, minimum similarity and sketch parameters when building the index.
- **Output:** Index construction returns `SketchIndex`; querying returns nearest-neighbor hits ranked by similarity.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.similarity import build_sketch_index, nearest_neighbors

records = [
    DNARecord(DNASequence("ACGTACGT"), "a"),
    DNARecord(DNASequence("ACGTTCGT"), "b"),
]
index = build_sketch_index(records, k=3, num_hashes=100)
result = nearest_neighbors(records[0], index, top_k=2)
print([hit.record_id for hit in result.hits])
```

- **Example output:**

```text
['a', 'b']
```


## 6) Database Index {#6-sim-015-database-index}

- **Purpose:** Create a reusable Sketch index with parameters and record IDs for the reference sequence collection, and support saving, verification and reloading to avoid recalculating the reference summary for each Top-k query.
- **Method:** Compute parameterized sketches and save the index with integrity information as JSON.
- **API:** `dnakit.similarity.build_sketch_index(records[required], k[optional], num_hashes[optional], canonical[optional], seed[optional], max_records[optional])`, `dnakit.similarity.save_sketch_index(index[required], path[required], overwrite[optional])`, `dnakit.similarity.load_sketch_index(path[required])`.
- **Input:** Required `DNARecord` collection; optional k, hash number, canonical, seed, save path and record limit.
- **Output:** Returns `SketchIndex` with indexed records and parameters for later similarity queries.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNARecord, DNASequence
from dnakit.similarity import build_sketch_index, load_sketch_index, save_sketch_index

records = [DNARecord(DNASequence("ACGT"), "r1")]
with TemporaryDirectory() as directory:
    path = Path(directory) / "sketch-index.json"
    index = build_sketch_index(records, k=2, num_hashes=20)
    digest = save_sketch_index(index, path)
    print(load_sketch_index(path) == index, len(digest))
```

- **Example output:**

```text
True 64
```

## 7) Build a local BLAST database and read hits {#bio-4}

- **Purpose:** Search your own nucleotide sequence database.
- **Method:** Run `makeblastdb` and `blastn`, then parse XML with SearchIO.
- **API:** `dnakit.bioinformatics.alignment.build_blast_database(fasta, output_directory)`; `dnakit.bioinformatics.alignment.local_blast(query, database, output_directory)`; `dnakit.bioinformatics.alignment.parse_search_results(path, format='blast-xml')`
- **Input:** Existing `reference.fa` and `query.fa`; BLAST+ on PATH and `bio-formats`. Output directories must be new.
- **Output:** Database construction and search return `ToolResult` with the XML path in `artifacts`; parsing yields queries and their hits.
- **Example code:**

```python
from dnakit.bioinformatics.alignment import build_blast_database, local_blast, parse_search_results

database = build_blast_database("reference.fa", "new_blast_db")
result = local_blast("query.fa", database.directory / "database", "new_blast_search")
for query in parse_search_results(result.artifacts["results.xml"]):
    print(query.id, len(query))
```

- **Example output:**

Returns `ToolResult` with XML at `artifacts['results.xml']`; the loop prints each query ID and hit count. Actual hits depend on the input; this example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

## 8) Build a profile HMM and search DNA {#bio-5}

- **Purpose:** Build a model from an alignment and locate related sequence regions.
- **Method:** Build with `hmmbuild`, search with `nhmmer`, and parse nucleotide hits; use `hmmsearch` for proteins.
- **API:** `dnakit.bioinformatics.alignment.build_profile(alignment, output_directory, alphabet='dna')`; `dnakit.bioinformatics.alignment.search_profile(profile, sequences, output_directory, method='nhmmer')`; `dnakit.bioinformatics.alignment.parse_nhmmer_table(path)`
- **Input:** Existing `aligned.sto` and `targets.fa`; HMMER on PATH; new output directories.
- **Output:** Model building and search return `ToolResult` with the hit-table path in `artifacts`; parsing yields hit mappings with positions, strand, and score.
- **Example code:**

```python
from dnakit.bioinformatics.alignment import build_profile, search_profile, parse_nhmmer_table

model = build_profile("aligned.sto", "new_profile", alphabet="dna")
result = search_profile(model.artifacts["profile.hmm"], "targets.fa", "new_hmm_search")
for hit in parse_nhmmer_table(result.artifacts["hits.tbl"]):
    print(hit["target"], hit["start"], hit["end"], hit["strand"])
```

- **Example output:**

Each hit is a dictionary with `target`, `query`, `start/end`, `strand`, `evalue` and `score`; coordinates are zero-based, half-open. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

<!-- BEGIN GENERATED API: blast -->

<span id="blast"></span>

## 9) Search similar protein sequences with remote BLAST {#dnakit.databases.blast_protein_search}

- **Purpose:** Search similar protein sequences with remote BLAST.
- **Method:** Use `NCBIBlast` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.blast_protein_search(*, sequence, database=None, expect=None, hitlist_size=None, verbose=True)`
- **Input:** Provide `sequence`.
- **Output:** Returns a mapping with `query_id`, `query_length`, `database`, `algorithm` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.blast_protein_search(
    sequence='MKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL',
    database='swissprot',
    hitlist_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query_id": "<string>",
  "query_length": "<integer>",
  "database": "<string>",
  "algorithm": "<string>",
  "alignments": [
    {
      "hit_id": "<string>",
      "hit_def": "<string>",
      "hit_length": "<integer>",
      "hsps": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/09_search.md#dnakit.databases.blast_protein_search-parameters)

<!-- END GENERATED API: blast -->
