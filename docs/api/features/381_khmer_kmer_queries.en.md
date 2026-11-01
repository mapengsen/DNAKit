# khmer k-mer queries: hashes, abundances, de Bruijn graph neighbourhoods and paths

This page covers the khmer rows of work-list pack P27 that query a k-mer table or the de Bruijn graph of a set of reads: k-mer hash values, per-k-mer abundance statistics, abundance trimming, spectral error positions, nodegraph union, neighbours and degree of a k-mer, connected-component size, branch points and unbranched paths. Every function builds the khmer table it needs from one input: a saved table (`graph_path`), local FASTA/FASTQ reads (`reads_path`) or DNA strings (`sequences`). Rows that differ only by a parameter are merged: the two hash rows are one function with `as_set`; the minimum, maximum, median/mean/SD, per-k-mer count and median-threshold rows are one function (`khmer_sequence_abundance`, which also covers `count-median.py`); neighbour list and degree are one function; component size and radius count are one function; traversing and assembling a linear path are one function with `output`. Loading a table is already covered by `query_saved_khmer_table` and `inspect_khmer_graph`, saving by `khmer_build_countgraph` and `khmer_build_nodegraph`, and MurmurHash3 hashing by `hash_kmer_murmur3`.

khmer is an optional dependency (`pip install "dnakit[bio-khmer]"`) and is imported only when a function runs; tested with khmer 2.1.1. Defaults follow khmer: k = 32, 4 tables, 1e6 buckets per table; the examples use k = 15 because the reads are short. A k-mer and its reverse complement are one node of the graph.

## Example files

`reads.fa` holds three reads: `a` and `c` share their first 30 bases, so the de Bruijn graph branches there, and `b` is unrelated. `other.fa` holds two reads that share their first 32 bases.

- **Example file:** `reads.fa`

```text
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

- **Example file:** `other.fa`

```text
>x1
AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC
>x2
AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT
```

## 1) `dnakit.microbiome.khmer_kmer_hashes()`

- **Purpose:** List the hash value of every k-mer of a DNA sequence in order, or only the distinct values.
- **Method:** Creates an empty khmer countgraph of the requested k and calls `Countgraph.get_kmer_hashes` (ordered list) or `get_kmer_hashes_as_hashset` (`as_set=True`, returned sorted). The values are khmer's 2-bit forward hash of the canonical k-mer, so a k-mer and its reverse complement give one value; `kmers` lists the k-mers the values stand for. Input is upper-cased. Covers work-list rows G1219 and G1220.
- **API:** `dnakit.microbiome.khmer_kmer_hashes(sequence, *, ksize=32, as_set=False)`.
- **Input:** `sequence` (DNA of A, C, G, T), `ksize` (1 to 32, default 32), `as_set`.
- **Output:** `n_hashes`, `hashes`, `kmers`, `provenance`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_kmer_hashes

seq = "TCGCTGCTGTCGGACTCCTA"
r = khmer_kmer_hashes(seq, ksize=15, verbose=False)
print(r["n_hashes"], r["hashes"][:2])
print(r["kmers"][:2])
s = khmer_kmer_hashes(seq + seq, ksize=15, as_set=True, verbose=False)
print(s["n_hashes"])
```

- **Example output:**

```text
6 [463369970, 229412411]
['TCGCTGCTGTCGGAC', 'AGTCCGACAGCAGCG']
20
```

## 2) `dnakit.microbiome.khmer_graph_summary()`

- **Purpose:** Load a saved khmer table and estimate how many distinct k-mers it holds.
- **Method:** Checks the file header, loads the table with `khmer.load_countgraph`, `load_nodegraph`, `load_counttable` or `load_nodetable` (`graph_type`) and reports k, the table sizes and the occupied buckets. khmer's own `n_unique_kmers()` counter is not stored in the file and reads 0 after loading, so the distinct k-mer number is estimated by linear counting on the first hash table: `-size * ln(1 - occupied / size)`. Covers work-list row G1215.
- **API:** `dnakit.microbiome.khmer_graph_summary(graph_path, *, graph_type='countgraph')`.
- **Input:** `graph_path` (saved table), `graph_type` (`countgraph`, `nodegraph`, `counttable` or `nodetable`).
- **Output:** `ksize`, `hashsizes`, `n_occupied`, `estimated_unique_kmers` (None when the first table is full), `provenance`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_gap_functions import khmer_build_countgraph
from dnakit.microbiome.khmer_p27_graph_functions import khmer_graph_summary

table = os.path.join(tempfile.mkdtemp(), "reads.countgraph")
khmer_build_countgraph("reads.fa", table, ksize=15, verbose=False)
r = khmer_graph_summary(table, verbose=False)
print(r["ksize"], r["n_occupied"], r["estimated_unique_kmers"])
```

- **Example output:**

```text
15 122 122
```

## 3) `dnakit.microbiome.khmer_sequence_abundance()`

- **Purpose:** Abundance statistics of query sequences against a countgraph: minimum, maximum, median, mean, SD, per-k-mer counts and a median threshold test.
- **Method:** Builds or loads a countgraph (`graph_path`, `reads_path` or `sequences`) and, for each query, calls `get_min_count`, `get_max_count`, `get_median_count` (median, mean, SD), optionally `get_kmer_counts` and, with `median_threshold`, `median_at_least`. Queries come from `sequences_to_query` and/or a FASTA/FASTQ `query_path`; as in khmer's `count-median.py`, N in a file query is read as A and queries shorter than k are skipped. `output_csv` writes the `count-median.py` table (name, median, average, stddev, seqlen). Covers work-list rows G1221 to G1225.
- **API:** `dnakit.microbiome.khmer_sequence_abundance(*, sequences_to_query=None, query_path=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, median_threshold=None, include_counts=False, output_csv=None, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequences_to_query` and/or `query_path`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`; `median_threshold`, `include_counts`, `output_csv`, `max_reads`, `max_bases`.
- **Output:** `n_queries`, `skipped_shorter_than_k`, `results` (`name`, `length`, `min`, `max`, `median`, `mean`, `sd`, optional `median_at_least` and `counts`), `output_csv`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_sequence_abundance

r = khmer_sequence_abundance(
    sequences_to_query=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA"],
    reads_path="reads.fa",
    ksize=15,
    median_threshold=2,
    include_counts=True,
    verbose=False,
)
row = r["results"][0]
print(row["min"], row["max"], row["median"], round(row["mean"], 3), row["median_at_least"])
print(row["counts"][:20])
```

- **Example output:**

```text
1 2 1 1.348 False
[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1]
```

## 4) `dnakit.microbiome.khmer_trim_by_abundance()`

- **Purpose:** Trim a sequence at the first k-mer whose countgraph abundance crosses a cutoff.
- **Method:** Builds or loads a countgraph and calls `Countgraph.trim_below_abundance` (`mode="trim_below_abundance"`, cut at the first k-mer with a count above the cutoff) or `trim_on_abundance` (`mode="trim_on_abundance"`, cut at the first k-mer with a count below the cutoff). Returns the trimmed sequence and khmer's trim position. Covers work-list row G1227.
- **API:** `dnakit.microbiome.khmer_trim_by_abundance(sequence, cutoff, *, mode='trim_below_abundance', graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`, `cutoff`, `mode`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`.
- **Output:** `trimmed_sequence`, `trim_position`, `input_length`, `mode`, `cutoff`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_trim_by_abundance

low = khmer_trim_by_abundance(
    "TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", 2, mode="trim_on_abundance", reads_path="reads.fa", ksize=15, verbose=False
)
print(low["trim_position"], low["trimmed_sequence"])
high = khmer_trim_by_abundance("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", 1, reads_path="reads.fa", ksize=15, verbose=False)
print(high["trim_position"], repr(high["trimmed_sequence"]))
```

- **Example output:**

```text
30 TCGCTGCTGTCGGACTCCTAGTTACGTGGC
0 ''
```

## 5) `dnakit.microbiome.khmer_spectral_error_positions()`

- **Purpose:** Locate positions of a sequence that lie in low-abundance k-mers, the signature of a sequencing error.
- **Method:** Builds or loads a countgraph and calls `Countgraph.find_spectral_error_positions(sequence, max_abundance)`. Positions are 0-based offsets in the query. Covers work-list row G1228.
- **API:** `dnakit.microbiome.khmer_spectral_error_positions(sequence, *, max_abundance=1, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`, `max_abundance` (default 1); `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`.
- **Output:** `positions`, `max_abundance`, `graph`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_spectral_error_positions

truth = "TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"
error = truth[:30] + ("A" if truth[30] != "A" else "C") + truth[31:]
r = khmer_spectral_error_positions(error, sequences=[truth] * 4 + [error], ksize=15, verbose=False)
print(r["positions"])
```

- **Example output:**

```text
[30]
```

## 6) `dnakit.microbiome.khmer_merge_nodegraphs()`

- **Purpose:** Union of two saved nodegraphs, written to a new file.
- **Method:** Loads both nodegraphs, calls `Nodegraph.update` (the second is merged into the first) and saves the result with `Nodegraph.save`. The graphs must have the same k, number of tables and table sizes. The output must not exist. Covers work-list row G1231.
- **API:** `dnakit.microbiome.khmer_merge_nodegraphs(first_graph_path, second_graph_path, output_path)`.
- **Input:** `first_graph_path`, `second_graph_path`, `output_path`.
- **Output:** `unique_kmers_first`, `unique_kmers_second`, `unique_kmers_merged`, `output_path`, `ksize`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_gap_functions import khmer_build_nodegraph
from dnakit.microbiome.khmer_p27_graph_functions import khmer_merge_nodegraphs

folder = tempfile.mkdtemp()
one = os.path.join(folder, "one")
two = os.path.join(folder, "two")
khmer_build_nodegraph("reads.fa", one, ksize=15, verbose=False)
khmer_build_nodegraph("other.fa", two, ksize=15, verbose=False)
r = khmer_merge_nodegraphs(one, two, os.path.join(folder, "union"), verbose=False)
print(r["unique_kmers_first"], r["unique_kmers_second"], r["unique_kmers_merged"])
```

- **Example output:**

```text
0 0 0
```

## 7) `dnakit.microbiome.khmer_kmer_neighbors()`

- **Purpose:** List the k-mers adjacent to a k-mer in the de Bruijn graph, with its degree.
- **Method:** Builds or loads a graph (`graph_type` nodegraph or countgraph) and calls `neighbors` and `kmer_degree`. A neighbour is a k-mer present in the graph that overlaps the query by k-1 bases on either side. Covers work-list rows G1232 and G1234.
- **API:** `dnakit.microbiome.khmer_kmer_neighbors(kmer, *, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`.
- **Input:** `kmer` (length k); `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`, `graph_type`.
- **Output:** `degree`, `neighbors` (k-mers), `kmer`, `graph`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_kmer_neighbors

fork = "TCCTAGTTACGTGGC"
r = khmer_kmer_neighbors(fork, reads_path="reads.fa", ksize=15, verbose=False)
print(r["degree"], r["neighbors"])
r = khmer_kmer_neighbors("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["degree"], r["neighbors"])
```

- **Example output:**

```text
3 ['CTCCTAGTTACGTGG', 'TGCCACGTAACTAGG', 'CCTAGTTACGTGGCG']
1 ['CGTGCGAGTGTTGCG']
```

## 8) `dnakit.microbiome.khmer_connected_component_size()`

- **Purpose:** Number of k-mers connected to a k-mer in the de Bruijn graph, optionally also the number within a graph radius.
- **Method:** Builds or loads a graph and calls `calc_connected_graph_size(kmer, max_size)`; with `count_within_radius` it also calls `count_kmers_within_radius(kmer, radius, max_size)`. `max_size` 0 means no limit, as in khmer. Covers work-list rows G1233 and G1235.
- **API:** `dnakit.microbiome.khmer_connected_component_size(kmer, *, max_size=0, count_within_radius=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`.
- **Input:** `kmer` (length k), `max_size`, `count_within_radius`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`, `graph_type`.
- **Output:** `component_size`, `kmers_within_radius`, `max_size`, `graph`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_connected_component_size

r = khmer_connected_component_size(
    "TCGCTGCTGTCGGAC", reads_path="reads.fa", ksize=15, count_within_radius=5, verbose=False
)
print(r["component_size"], r["kmers_within_radius"])
r = khmer_connected_component_size("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["component_size"])
```

- **Example output:**

```text
76 6
46
```

## 9) `dnakit.microbiome.khmer_high_degree_kmers()`

- **Purpose:** Find branch points: the k-mers of a sequence that have more than two neighbours in the graph.
- **Method:** Builds or loads a graph and calls `find_high_degree_nodes(sequence)`; the hash values are converted back to k-mers. Covers work-list row G1236.
- **API:** `dnakit.microbiome.khmer_high_degree_kmers(sequence, *, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`, `graph_type`.
- **Output:** `n_high_degree`, `high_degree_kmers` (sorted by hash value), `graph`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_high_degree_kmers

r = khmer_high_degree_kmers("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_high_degree"], r["high_degree_kmers"])
r = khmer_high_degree_kmers("CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA", reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_high_degree"])
```

- **Example output:**

```text
1 ['TCCTAGTTACGTGGC']
0
```

## 10) `dnakit.microbiome.khmer_linear_path()`

- **Purpose:** Follow the unbranched path that starts at a k-mer, as one assembled sequence or as the list of k-mers it passes.
- **Method:** Builds or loads a graph. `output="sequence"` calls `assemble_linear_path` and returns the path as a DNA string. `output="kmers"` calls `traverse_linear_path` and returns the traversed k-mers, their count and the branch k-mers met; `stop_sequence` supplies the branch points that end the traversal (the high-degree k-mers of that sequence); without it the whole connected component is followed. Covers work-list rows G1237 and G1238.
- **API:** `dnakit.microbiome.khmer_linear_path(kmer, *, output='sequence', stop_sequence=None, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, graph_type='nodegraph', max_reads=1000000, max_bases=100000000)`.
- **Input:** `kmer` (the start k-mer), `output`, `stop_sequence`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`, `graph_type`.
- **Output:** `path_sequence` and `path_length`, or `path_size`, `path_kmers` and `branch_kmers`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_graph_functions import khmer_linear_path

r = khmer_linear_path("CCGCAACACTCGCAC", reads_path="reads.fa", ksize=15, verbose=False)
print(r["path_length"], r["path_sequence"] == "CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA")
k = khmer_linear_path(
    "TCGCTGCTGTCGGAC", output="kmers", stop_sequence="TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", reads_path="reads.fa", ksize=15, verbose=False
)
print(k["path_size"], k["branch_kmers"])
```

- **Example output:**

```text
60 True
15 ['TCCTAGTTACGTGGC']
```
