# Universal search

Perform exact, subsequence, approximate, reverse complement and nearest neighbor searches in single or multiple DNA sequences.

## 1) SIM-001 · Exact search

**Function:** Search the record set for records that are identical to the query sequence character by character, and return the target index, start and end coordinates, and chain direction, used for member checking, exact copy positioning, and small-scale reference queries.

**API:** `dnakit.similarity.exact_search(query[required], targets[required], reverse_complement[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.

**Input:** Required query sequence and target sequence/`DNASet`; optional reverse complementary search, chain merging and resource limit.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import exact_search

result = exact_search(
    DNASequence("AC"),
    [DNASequence("AC"), DNASequence("TACG")],
)
print([(hit.target_index, hit.start, hit.end) for hit in result.matches])
```

**Example results:**

```text
[(0, 0, 2)]
```


## 2) SIM-002 · Subsequence search

**Function:** Find all precise occurrence positions of short query sequences in the target long sequence, and return the start and end coordinates and matching directions, which can be used to locate motifs, primers or known fragments.

**API:** `dnakit.similarity.subsequence_search(query[required], target[required], strand[optional], overlapping[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.

**Input:** Required query and target; optional `strand`, `overlapping`, chain merging and resource limit.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import subsequence_search

result = subsequence_search(DNASequence("ANA", alphabet="iupac"),
                            DNASequence("ANANA", alphabet="iupac"))
print([(hit.start, hit.end) for hit in result.matches])
```

**Example results:**

```text
[(0, 3), (2, 5)]
```


## 3) SIM-004 · Approximate matching

**Function:** Find approximate hits in the target sequence that allow a specified number of substitutions, insertions, or deletions, and return coordinates, distances, and matching fragments for searches that tolerate mutations or sequencing errors.

**API:** `dnakit.similarity.approximate_search(query[required], targets[required], max_distance[required], substitution_cost[optional], insertion_cost[optional], deletion_cost[optional], reverse_complement[optional], max_targets[optional], max_matches[optional], max_cells[optional])`.

**Input:** Required query, target and `max_distance`; optional three types of editing costs, reverse complementation and resource limit.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import approximate_search

result = approximate_search(DNASequence("ACG"), DNASequence("TTACGATG"),
                            max_distance=1)
print(any(hit.start == 2 and hit.distance == 0 for hit in result.matches))
```

**Example results:**

```text
True
```

## 4) SIM-005 · Reverse complementary search

**Function:** Search using the query sequence and its reverse complement simultaneously, returning the coordinates and strand direction of each hit to avoid missing the same fragment located in the reverse strand.

**API:** `dnakit.similarity.reverse_complement_search(query[required], target[required], overlapping[optional], merge_strands[optional], max_targets[optional], max_matches[optional])`.

**Input:** Required query and target; optional whether to allow overlap and whether to merge palindromic duplicate hits.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import reverse_complement_search

result = reverse_complement_search(DNASequence("ATG"), DNASequence("GGCATCC"))
print([(hit.start, hit.end, hit.strand.value) for hit in result.matches])
```

**Example results:**

```text
[(2, 5, 'reverse')]
```

## 5) SIM-014 · Nearest neighbor search

**Function:** Use k-mer Sketch of the query sequence to filter the Top-k records with the highest similarity in the existing index, and return the record ID, ranking and approximate similarity, which is suitable for fast candidate retrieval.

**API:** `dnakit.similarity.build_sketch_index(records[required], k[optional], num_hashes[optional], canonical[optional], seed[optional], max_records[optional])`, `dnakit.similarity.nearest_neighbors(query[required], index[required], top_k[optional], min_similarity[optional])`.

**Input:** Required query and `SketchIndex`; optional `top_k`, minimum similarity and sketch parameters when building the index.

**Sample code:**

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

**Example results:**

```text
['a', 'b']
```


## 6) SIM-015 · Database Index

**Function:** Create a reusable Sketch index with parameters and record IDs for the reference sequence collection, and support saving, verification and reloading to avoid recalculating the reference summary for each Top-k query.

**API:** `dnakit.similarity.build_sketch_index(records[required], k[optional], num_hashes[optional], canonical[optional], seed[optional], max_records[optional])`, `dnakit.similarity.save_sketch_index(index[required], path[required], overwrite[optional])`, `dnakit.similarity.load_sketch_index(path[required])`.

**Input:** Required `DNARecord` collection; optional k, hash number, canonical, seed, save path and record limit.

**Sample code:**

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

**Example results:**

```text
True 64
```
