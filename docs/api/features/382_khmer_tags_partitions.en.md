# khmer tags, partitions and stop tags: splitting reads into connected components

This page covers the khmer rows of work-list pack P27 about tags, partitions and stop tags. khmer marks a sparse set of tag k-mers in the de Bruijn graph of the reads, groups connected tags into partitions and annotates every read with its partition ID (the `do-partition.py` workflow). Highly connected k-mers (knots) can be saved as stop tags, at which partitions are cut. Reads are local FASTA/FASTQ files consumed by khmer itself, so its read cleaning applies (upper-casing, N read as A, reads shorter than k skipped).

Rows that differ only by a parameter or by a script step are merged: loading, saving, printing and extending the tag set are options of `khmer_tag_reads`; dividing tags into subsets, partitioning a subset (with or without abundance), merging subset maps, saving and loading partition maps, counting partitions, the subset size distribution, looking up a partition ID and annotating the reads are options of `khmer_partition_reads`, which is `do-partition.py` (`load-graph.py`, `partition-graph.py`, `merge-partitions.py` and `annotate-partitions.py` in one call); finding knots, loading, saving, printing and adding stop tags are options of `khmer_find_knots` (`make-initial-stoptags.py` and `find-knots.py`). Note that khmer's `subset_count_partitions` damages a subset map, so the function reads it only after the map was saved and merged.

khmer is an optional dependency (`pip install "dnakit[bio-khmer]"`); tested with khmer 2.1.1. Defaults follow khmer: k = 32, 4 tables, 1e6 buckets per table, subset size 1e5. The examples use k = 15.

## Example files

`reads.fa` holds three reads: `a` and `c` share their first 30 bases, so they fall in one partition, and `b` is unrelated.

- **Example file:** `reads.fa`

```text
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

## 1) `dnakit.microbiome.khmer_tag_reads()`

- **Purpose:** Insert reads into a khmer graph and mark tag k-mers; load, extend, save and print the tag set.
- **Method:** Builds a graph from `reads_path` and/or `sequences` and tags it (`consume_seqfile_and_tag`, `consume_and_tag`). Optional khmer calls: `load_tagset_path` (`load_tagset`), `extra_tag_kmers` (`add_tag`), `save_tagset_path` (`save_tagset`), `tagset_text_path` (`print_tagset`, one k-mer per line). Output files must not exist. Covers work-list rows G1239, G1243, G1244, G1245 and G1246.
- **API:** `dnakit.microbiome.khmer_tag_reads(*, reads_path=None, sequences=None, load_tagset_path=None, extra_tag_kmers=None, save_tagset_path=None, tagset_text_path=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `reads_path` (file or list) and/or `sequences`; `load_tagset_path`, `extra_tag_kmers`, `save_tagset_path`, `tagset_text_path`; `graph_type`, `ksize`, `table_size`, `n_tables`, `max_reads`, `max_bases`.
- **Output:** `n_tags`, `tags`, `n_unique_kmers`, `extra_tags_added`, `saved_tagset`, `tagset_text`, `graph`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_tag_reads

folder = tempfile.mkdtemp()
binary = os.path.join(folder, "tags.bin")
text = os.path.join(folder, "tags.txt")
r = khmer_tag_reads(
    reads_path="reads.fa",
    ksize=15,
    extra_tag_kmers=["A" * 15],
    save_tagset_path=binary,
    tagset_text_path=text,
    verbose=False,
)
print(r["n_tags"], r["extra_tags_added"], "A" * 15 in r["tags"])
print(len(open(text).read().split()))
print(khmer_tag_reads(load_tagset_path=binary, ksize=15, verbose=False)["n_tags"])
```

- **Example output:**

```text
7 1 True
7
7
```

## 2) `dnakit.microbiome.khmer_tag_positions()`

- **Purpose:** Positions of the tag k-mers along a sequence.
- **Method:** Builds and tags a graph from `reads_path` and/or `sequences`, then calls `get_tags_and_positions(sequence)`. Positions are 0-based start offsets in the query. Covers work-list row G1240.
- **API:** `dnakit.microbiome.khmer_tag_positions(sequence, *, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`; `reads_path` and/or `sequences`; `graph_type`, `ksize`, `table_size`, `n_tables`.
- **Output:** `tags`: list of `position` and `kmer`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_tag_positions

r = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)
for tag in r["tags"]:
    print(tag["position"], tag["kmer"])
```

- **Example output:**

```text
19 AACGCCACGTAACTA
46 ACCACGACGGCAGGC
```

## 3) `dnakit.microbiome.khmer_reachable_tags()`

- **Purpose:** Find the tags connected to a k-mer, and optionally assign them one partition ID.
- **Method:** Builds and tags a graph, then calls `find_all_tags_list(kmer)`. With `assign_partition=True` it also calls `find_all_tags` and `assign_partition_id`, which gives the tags the ID of an existing partition if one of them has one, and a new ID otherwise. Covers work-list rows G1241 and G1249.
- **API:** `dnakit.microbiome.khmer_reachable_tags(kmer, *, assign_partition=False, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `kmer` (length k), `assign_partition`; `reads_path` and/or `sequences`; `graph_type`, `ksize`, `table_size`, `n_tables`.
- **Output:** `n_tags`, `tags`, `partition_id` (None unless `assign_partition`).
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_reachable_tags, khmer_tag_positions

tags = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["tags"]
start = tags[0]["kmer"]
r = khmer_reachable_tags(start, reads_path="reads.fa", ksize=15, verbose=False)
print(r["n_tags"], r["partition_id"])
r = khmer_reachable_tags(
    start, reads_path="reads.fa", ksize=15, assign_partition=True, verbose=False
)
print(r["partition_id"])
```

- **Example output:**

```text
2 None
2
```

## 4) `dnakit.microbiome.khmer_unique_paths()`

- **Purpose:** Pieces of a sequence whose k-mers are mostly absent from a graph.
- **Method:** Builds a graph from `reads_path` and/or `sequences` (an empty graph if none is given) and calls `extract_unique_paths(sequence, min_length, min_unique_fraction)`. `min_length` defaults to k and may not be smaller than k. khmer's documentation marks this method for removal; it is the khmer operation named in work-list row G1242.
- **API:** `dnakit.microbiome.khmer_unique_paths(sequence, *, min_length=None, min_unique_fraction=1.0, reads_path=None, sequences=None, graph_type='nodegraph', ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`, `min_length`, `min_unique_fraction`; `reads_path` and/or `sequences`; `graph_type`, `ksize`, `table_size`, `n_tables`.
- **Output:** `paths` (list of sequences), `min_length`, `min_unique_fraction`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_unique_paths

print(len(khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", ksize=15, verbose=False)["paths"][0]))
print(khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["paths"])
r = khmer_unique_paths("TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA", sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], ksize=15, verbose=False)
print([len(p) for p in r["paths"]])
```

- **Example output:**

```text
60
[]
[44]
```

## 5) `dnakit.microbiome.khmer_partition_reads()`

- **Purpose:** Partition reads into connected components of their de Bruijn graph and annotate every read with its partition ID.
- **Method:** Follows khmer's `do-partition.py`: `consume_seqfile_and_tag`, `divide_tags_into_subsets` (`subset_size`), `do_subset_partition` for every subset, `merge_subset`, then `output_partitions`, which writes `<name>.part` with the partition ID after a tab in each read name. Options, each a khmer call: `stoptags_path` (`load_stop_tags`; partitions are cut there), `stop_big_traversals`, `abundance_range=[min, max]` (`do_subset_partition_with_abundance` on a countgraph), `save_partitionmap_dir` (`save_subset_partitionmap`, files `<name>.subset.N.pmap`), `partitionmap_paths` (skip partitioning and merge saved maps with `merge_subset_from_disk`), `query_kmers` (`get_partition_id`; meaningful for tag k-mers). `annotate=False` skips the `.part` files. Partition ID 0 means unassigned. Covers work-list rows G1230, G1247, G1248, G1250, G1251, G1252, G1254, G1255, G1256, G1258 and G1280.
- **API:** `dnakit.microbiome.khmer_partition_reads(reads_path, *, output_dir=None, subset_size=100000, stop_big_traversals=False, stoptags_path=None, abundance_range=None, partitionmap_paths=None, save_partitionmap_dir=None, query_kmers=None, annotate=True, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `reads_path` (file or list), `output_dir`, `subset_size`, `stop_big_traversals`, `stoptags_path`, `abundance_range`, `partitionmap_paths`, `save_partitionmap_dir`, `query_kmers`, `annotate`, `ksize`, `table_size`, `n_tables`.
- **Output:** `partition_count` (`n_partitions`, `n_unassigned_tags`; None when saved maps were merged), `n_subsets`, `per_subset` (size distribution of each subset), `annotated` (output file and number of partitions), `saved_partitionmaps`, `partition_ids`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_partition_reads

folder = tempfile.mkdtemp()
r = khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, subset_size=3, verbose=False)
print(r["partition_count"], r["n_subsets"])
print(r["annotated"][0]["n_partitions"])
for line in open(os.path.join(folder, "reads.fa.part")):
    if line.startswith(">"):
        print(line.rstrip().replace("\t", " -> partition "))
```

- **Example output:**

```text
{'n_partitions': 2, 'n_unassigned_tags': 0} 2
2
>a -> partition 2
>b -> partition 3
>c -> partition 2
```

## 6) `dnakit.microbiome.khmer_join_partitions()`

- **Purpose:** Partition reads, then join the partitions of two tag k-mers.
- **Method:** Partitions the reads as `khmer_partition_reads` does, checks that both k-mers are tags, looks up their partition IDs (`get_partition_id`), joins them with `join_partitions` and reports the IDs and partition counts before and after. Covers work-list row G1257.
- **API:** `dnakit.microbiome.khmer_join_partitions(reads_path, first_tag_kmer, second_tag_kmer, *, subset_size=100000, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `reads_path` (file or list), `first_tag_kmer`, `second_tag_kmer`; `subset_size`, `ksize`, `table_size`, `n_tables`.
- **Output:** `partition_ids_before`, `partition_ids_after`, `partition_count_before`, `partition_count_after`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_partition_functions import khmer_join_partitions, khmer_tag_positions

tag_a = khmer_tag_positions("TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT", reads_path="reads.fa", ksize=15, verbose=False)["tags"][0]["kmer"]
tag_b = khmer_tag_positions("CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA", reads_path="reads.fa", ksize=15, verbose=False)["tags"][0]["kmer"]
r = khmer_join_partitions("reads.fa", tag_a, tag_b, ksize=15, verbose=False)
print(r["partition_ids_before"], r["partition_count_before"])
print(r["partition_ids_after"], r["partition_count_after"])
```

- **Example output:**

```text
[2, 3] [2, 0]
[2, 2] [1, 0]
```

## 7) `dnakit.microbiome.khmer_load_partitioned_reads()`

- **Purpose:** Insert partition-annotated reads into a new graph and rebuild its partitions.
- **Method:** Calls `consume_partitioned_fasta` on a `.part` file: each read is counted and tagged and its tag receives the partition ID written in the read name; then `count_partitions`. Returns the tags with their partition IDs; `query_kmers` looks up other k-mers (`get_partition_id`; 0 unless the k-mer is one of the rebuilt tags). Covers work-list row G1253.
- **API:** `dnakit.microbiome.khmer_load_partitioned_reads(partitioned_path, *, query_kmers=None, ksize=32, table_size=1000000, n_tables=4)`.
- **Input:** `partitioned_path`, `query_kmers`, `ksize`, `table_size`, `n_tables`.
- **Output:** `n_reads`, `n_kmers_consumed`, `n_partitions`, `n_unassigned_tags`, `n_tags`, `tag_partition_ids`, `partition_ids`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_load_partitioned_reads, khmer_partition_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_load_partitioned_reads(os.path.join(folder, "reads.fa.part"), ksize=15, verbose=False)
print(r["n_reads"], r["n_kmers_consumed"], r["n_partitions"], r["n_tags"])
print(sorted(set(r["tag_partition_ids"].values())))
```

- **Example output:**

```text
3 138 2 2
[2, 3]
```

## 8) `dnakit.microbiome.khmer_find_knots()`

- **Purpose:** Find highly connected k-mers (knots) by repartitioning the largest partition and return them as stop tags.
- **Method:** `stage="initial"` follows `make-initial-stoptags.py`: partition the first subset (default 10000 tags) and call `repartition_largest_partition` with count threshold 5. `stage="knots"` follows `find-knots.py`: partition all subsets (default 100000 tags), repartition each (count threshold 2), merge, and repartition the merged map. `excursion_distance` (40) and `excursion_kmer_threshold` (200) are khmer's values; the excursion countgraph has `counting_table_size` buckets. Stop-tag options: `load_stoptags_path` (`load_stop_tags`), `extra_stop_kmers` (`add_stop_tag`), `save_stoptags_path` (`save_stop_tags`), `stoptags_text_path` (`print_stop_tags`). Verified identical to the two scripts. Covers work-list rows G1259, G1260, G1261, G1262 and G1264.
- **API:** `dnakit.microbiome.khmer_find_knots(reads_path, *, stage='knots', subset_size=None, excursion_distance=40, excursion_kmer_threshold=200, excursion_kmer_count_threshold=None, load_stoptags_path=None, extra_stop_kmers=None, save_stoptags_path=None, stoptags_text_path=None, ksize=32, table_size=1000000, n_tables=4, counting_table_size=1000000, max_reads=1000000, max_bases=100000000)`.
- **Input:** `reads_path`, `stage`, `subset_size`, `excursion_distance`, `excursion_kmer_threshold`, `excursion_kmer_count_threshold`, stop-tag options, `ksize`, `table_size`, `n_tables`, `counting_table_size`.
- **Output:** `n_stop_tags`, `stop_tags`, thresholds used, `saved_stoptags`, `stoptags_text`.
- **Example code:**

```python
import os
import random
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_find_knots

rng = random.Random(11)


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


repeat = dna(100)
folder = tempfile.mkdtemp()
reads = os.path.join(folder, "knot.fa")
with open(reads, "w") as handle:
    for i in range(40):
        handle.write(f">k{i}\n{dna(120)}{repeat}{dna(120)}\n")
stoptags = os.path.join(folder, "knot.stoptags")
r = khmer_find_knots(reads, stage="initial", ksize=21, save_stoptags_path=stoptags, verbose=False)
print(r["stage"], r["n_stop_tags"], r["excursion_kmer_count_threshold"])
r = khmer_find_knots(reads, stage="knots", ksize=21, load_stoptags_path=stoptags, verbose=False)
print(r["stage"], r["n_stop_tags"], r["excursion_kmer_count_threshold"])
```

- **Example output:**

```text
initial 1572 5
knots 1572 2
```

## 9) `dnakit.microbiome.khmer_filter_stoptags()`

- **Purpose:** Trim reads at stop tags and write the kept reads.
- **Method:** Follows `filter-stoptags.py`: loads stop tags (`stoptags_path` with `load_stop_tags`, and/or `stop_kmers` with `add_stop_tag`), drops reads containing N, cuts every other read at the first stop tag with `trim_on_stoptags` and keeps it only if at least k bases remain. Output is FASTA (`>name` and the trimmed sequence). Reads with symbols other than A, C, G, T are dropped. Covers work-list row G1263.
- **API:** `dnakit.microbiome.khmer_filter_stoptags(reads_path, output_path, *, stoptags_path=None, stop_kmers=None, ksize=32)`.
- **Input:** `reads_path`, `output_path` (must not exist), `stoptags_path` and/or `stop_kmers`, `ksize`.
- **Output:** `n_reads`, `n_kept`, `n_trimmed`, `n_dropped_with_n`, `output_path`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_filter_stoptags

out = os.path.join(tempfile.mkdtemp(), "filtered.fa")
r = khmer_filter_stoptags("reads.fa", out, stop_kmers=["GTTGCTCCACAGGTA"], ksize=15, verbose=False)
print(r["n_reads"], r["n_kept"], r["n_trimmed"])
for line in open(out):
    print(line.rstrip())
```

- **Example output:**

```text
3 3 1
>a
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGG
>b
CCGCAACACTCGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGA
>c
TCGCTGCTGTCGGACTCCTAGTTACGTGGCAGACAACCCTCTGGTTCTTTCCCGTCCGTA
```

## 10) `dnakit.microbiome.khmer_extract_partitions()`

- **Purpose:** Group partition-annotated reads by partition size and write the partition size distribution.
- **Method:** Reproduces `extract-partitions.py`. It reads `.part` files, counts reads per partition and writes `<prefix>.dist` with the columns reads, partitions with that many reads, cumulative partitions and cumulative reads. With `write_groups=True` (the script without `-n`) partitions with more than `min_partition_size` reads are packed, smallest first, into `<prefix>.groupNNNN.fa` (or `.fq`) of about `max_group_size` reads; `write_unassigned` also writes the reads of partition 0. Verified identical to the script. Covers work-list rows G1256 (the `-n` option) and G1281.
- **API:** `dnakit.microbiome.khmer_extract_partitions(partitioned_paths, output_prefix, *, write_groups=True, min_partition_size=5, max_group_size=1000000, write_unassigned=False)`.
- **Input:** `partitioned_paths` (file or list), `output_prefix`, `write_groups`, `min_partition_size` (5), `max_group_size` (1000000), `write_unassigned`.
- **Output:** `n_partitions`, `n_unassigned_reads`, `size_distribution`, `distribution_path`, `group_files`, `unassigned_path`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_extract_partitions, khmer_partition_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("reads.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_extract_partitions(
    os.path.join(folder, "reads.fa.part"),
    os.path.join(folder, "groups"),
    min_partition_size=0,
    verbose=False,
)
print(r["n_partitions"], r["n_unassigned_reads"])
print(open(r["distribution_path"]).read().strip())
print(len(r["group_files"]), open(r["group_files"][0]).read().count(">"))
```

- **Example output:**

```text
2 0
1 1 1 1
2 1 2 3
1 3
```
