# khmer labelled graphs, HyperLogLog merging, read alignment and read-file tools

This page covers the remaining khmer rows of work-list pack P27: labelled graphs (`GraphLabels`), merging two HyperLogLog k-mer counters, aligning a read to a k-mer graph with khmer's hidden Markov model, and the read-file scripts `extract-paired-reads.py`, `split-paired-reads.py` and `readstats.py`. A labelled graph is a nodegraph whose tag k-mers carry integer labels. Reads in files get their record index in the file as label (`consume_seqfile_and_tag_with_labels`); strings get explicit labels.

Rows that differ only by a parameter are merged: reading reads from a file or from strings, loading and saving the labels and tags are options of `khmer_label_reads`; sweeping labels and sweeping tags are one function with `output`; the seeded-both-directions and forward-only alignments are one function with `forward_only`. Parsing FASTA/FASTQ reads and interleaved pairs (`ReadParser.iter_reads`, `iter_read_pairs`) is already covered by `read_paired_fastx`. khmer is an optional dependency (`pip install "dnakit[bio-khmer]"`); tested with khmer 2.1.1. Defaults follow khmer: k = 32, 4 tables, 1e6 buckets per table; the examples use k = 15. khmer's `GraphLabels` object is a thin wrapper whose own graph methods crash, so DNAKit only uses the label methods on it and the graph methods on its `graph` attribute.

## Example files

`lab.fa` holds two reads that share their first 32 bases. `mixed.fa` holds two read pairs (`/1` and `/2`) and one unpaired read. `paired.fa` holds two interleaved pairs.

- **Example file:** `lab.fa`

```text
>x1
AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC
>x2
AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT
```

- **Example file:** `mixed.fa`

```text
>p1/1
ACGTACGTACGGTT
>p1/2
GGGGCCCCAATTGC
>s1
TTTTTTTTTTGGAA
>p2/1
AAAACCCCGGTTAC
>p2/2
CCCCAAAATTGGCA
```

- **Example file:** `paired.fa`

```text
>p1/1
ACGTACGTACGGTT
>p1/2
GGGGCCCCAATTGC
>p2/1
AAAACCCCGGTTAC
>p2/2
CCCCAAAATTGGCA
```

## 1) `dnakit.microbiome.khmer_label_reads()`

- **Purpose:** Insert reads into a labelled graph and tag their k-mers with labels; load and save the labels and tags.
- **Method:** Builds a `khmer.GraphLabels` object. Files in `reads_path` use `consume_seqfile_and_tag_with_labels` (label = record index in the file); `sequences` use `consume_sequence_and_tag_with_labels` with `labels` (default: position in the list). `load_labels_path` first calls `load_labels_and_tags`; `save_labels_path` calls `save_labels_and_tags` (the file must not exist). Covers work-list rows G1265, G1272 and G1273.
- **API:** `dnakit.microbiome.khmer_label_reads(*, reads_path=None, sequences=None, labels=None, load_labels_path=None, save_labels_path=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `reads_path` (file or list) or `sequences` with optional `labels`; `load_labels_path`, `save_labels_path`; `ksize`, `table_size`, `n_tables`, `max_reads`, `max_bases`.
- **Output:** `n_labels`, `n_tags`, `tag_labels` (labels of every tag k-mer), `saved_labels`, `loaded_labels`, `graph`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_reads

r = khmer_label_reads(reads_path="lab.fa", ksize=15, verbose=False)
print(r["n_labels"], r["n_tags"], sorted({tuple(v) for v in r["tag_labels"].values()}))
saved = os.path.join(tempfile.mkdtemp(), "labels.bin")
s = khmer_label_reads(
    sequences=["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], labels=[7, 9], ksize=15, save_labels_path=saved, verbose=False
)
print(s["n_labels"], sorted({x for v in s["tag_labels"].values() for x in v}))
print(khmer_label_reads(load_labels_path=saved, ksize=15, verbose=False)["n_labels"])
```

- **Example output:**

```text
2 4 [(0,), (1,)]
2 [7, 9]
2
```

## 2) `dnakit.microbiome.khmer_label_partitioned_reads()`

- **Purpose:** Insert partition-annotated reads into a labelled graph, using the partition ID as label.
- **Method:** Calls `consume_partitioned_fasta_and_tag_with_labels` on a `.part` file written by `khmer_partition_reads`. Covers work-list row G1266.
- **API:** `dnakit.microbiome.khmer_label_partitioned_reads(partitioned_path, *, ksize=32, table_size=1000000, n_tables=4)`.
- **Input:** `partitioned_path`; `ksize`, `table_size`, `n_tables`.
- **Output:** `n_labels`, `n_tags`, `tag_labels`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_partition_functions import khmer_partition_reads
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_partitioned_reads

folder = tempfile.mkdtemp()
khmer_partition_reads("lab.fa", output_dir=folder, ksize=15, verbose=False)
r = khmer_label_partitioned_reads(os.path.join(folder, "lab.fa.part"), ksize=15, verbose=False)
print(r["n_labels"], r["n_tags"])
```

- **Example output:**

```text
1 4
```

## 3) `dnakit.microbiome.khmer_sweep_neighborhood()`

- **Purpose:** Find the labels or the tags in the neighbourhood of a sequence in a labelled graph.
- **Method:** Builds a labelled graph as `khmer_label_reads` does, then calls `sweep_label_neighborhood` (`output="labels"`, returns the labels) or `sweep_tag_neighborhood` (`output="tags"`, returns the tag k-mers). `sweep_range` is the traversal range in k-mers; `break_on_stop_tags` and `stop_big_traversals` are khmer's options. Covers work-list rows G1267 and G1268.
- **API:** `dnakit.microbiome.khmer_sweep_neighborhood(sequence, *, output='labels', sweep_range=0, break_on_stop_tags=False, stop_big_traversals=False, reads_path=None, sequences=None, labels=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `sequence`, `output`, `sweep_range`, `break_on_stop_tags`, `stop_big_traversals`; `reads_path` or `sequences` with `labels`; `ksize`, `table_size`, `n_tables`.
- **Output:** `labels` or `tags`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_sweep_neighborhood

query = "AGACTACTTATGAGGCCATACCAGGGTCGT"
r = khmer_sweep_neighborhood(query, reads_path="lab.fa", sweep_range=40, ksize=15, verbose=False)
print(r["labels"])
r = khmer_sweep_neighborhood(
    query, output="tags", reads_path="lab.fa", sweep_range=40, ksize=15, verbose=False
)
print(r["tags"])
```

- **Example output:**

```text
[0, 1]
['AAAACGACCCTGGTA', 'TACCAGGGTCGTTTG']
```

## 4) `dnakit.microbiome.khmer_tag_labels()`

- **Purpose:** Look up the labels carried by one tag k-mer.
- **Method:** Builds a labelled graph as `khmer_label_reads` does and calls `get_tag_labels(tag_kmer)`; the result is empty if the k-mer is not a tag. Covers work-list row G1269.
- **API:** `dnakit.microbiome.khmer_tag_labels(tag_kmer, *, reads_path=None, sequences=None, labels=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `tag_kmer` (length k); `reads_path` or `sequences` with `labels`; `ksize`, `table_size`, `n_tables`.
- **Output:** `labels`, `tag_kmer`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_reads, khmer_tag_labels

tags = khmer_label_reads(reads_path="lab.fa", ksize=15, verbose=False)["tag_labels"]
for tag in sorted(tags)[:3]:
    print(tag, khmer_tag_labels(tag, reads_path="lab.fa", ksize=15, verbose=False)["labels"])
print(khmer_tag_labels("A" * 15, reads_path="lab.fa", ksize=15, verbose=False)["labels"])
```

- **Example output:**

```text
AAAACGACCCTGGTA [1]
ACGACTAGCCAAGCG [1]
CAGCCATAGTCCAAC [0]
[]
```

## 5) `dnakit.microbiome.khmer_label_across_high_degree_nodes()`

- **Purpose:** Label the k-mers of contigs across the branch points they share.
- **Method:** Inserts all contigs into one labelled graph, finds the high-degree k-mers (`find_high_degree_nodes`) and calls `label_across_high_degree_nodes(contig, hdn, label)` for each contig; labels default to 1, 2, ... in list order. Returns the branch k-mers and the labels each received. Covers work-list row G1270.
- **API:** `dnakit.microbiome.khmer_label_across_high_degree_nodes(contigs, *, labels=None, ksize=32, table_size=1000000, n_tables=4)`.
- **Input:** `contigs` (each longer than k), `labels`; `ksize`, `table_size`, `n_tables`.
- **Output:** `n_high_degree`, `high_degree_labels`, `contig_labels`, `n_labels`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_label_across_high_degree_nodes

r = khmer_label_across_high_degree_nodes(["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], ksize=15, verbose=False)
print(r["n_high_degree"], r["contig_labels"], r["n_labels"])
print(r["high_degree_labels"])
```

- **Example output:**

```text
1 [1, 2] 2
{'AAACGACCCTGGTAT': [1, 2]}
```

## 6) `dnakit.microbiome.khmer_assemble_labeled_path()`

- **Purpose:** Assemble the paths through a branching graph from a start k-mer, using labels to choose between branches.
- **Method:** Labels the contigs as `khmer_label_across_high_degree_nodes` does, then calls `assemble_labeled_path(start_kmer)`. Returns one sequence per consistent path, sorted. Covers work-list row G1271.
- **API:** `dnakit.microbiome.khmer_assemble_labeled_path(start_kmer, contigs, *, labels=None, ksize=32, table_size=1000000, n_tables=4)`.
- **Input:** `start_kmer` (length k), `contigs`, `labels`; `ksize`, `table_size`, `n_tables`.
- **Output:** `n_paths`, `paths`, `contig_labels`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_assemble_labeled_path

r = khmer_assemble_labeled_path("AGACTACTTATGAGG", ["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT"], ksize=15, verbose=False)
print(r["n_paths"], [len(p) for p in r["paths"]])
print(all(p in ("AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC", "AGACTACTTATGAGGCCATACCAGGGTCGTTTTCCGGGTATTGGCCGCTTGGCTAGTCGT") for p in r["paths"]))
```

- **Example output:**

```text
2 [60, 60]
True
```

## 7) `dnakit.microbiome.khmer_merge_hll_counters()`

- **Purpose:** Estimate the distinct k-mers of two read sets and of their union by merging HyperLogLog counters.
- **Method:** Fills one `khmer.HLLCounter(error_rate, ksize)` per input (`consume_seqfile` for a file, `consume_string` for strings), records both estimates (`estimate_cardinality`), merges the second counter into the first (`merge`) and records the estimate of the union. Defaults follow khmer's `unique-kmers.py`: error rate 0.01, k = 20. Estimates carry the stated relative error. Covers work-list row G1275.
- **API:** `dnakit.microbiome.khmer_merge_hll_counters(*, first_reads_path=None, first_sequences=None, second_reads_path=None, second_sequences=None, error_rate=0.01, ksize=20)`.
- **Input:** `first_reads_path` or `first_sequences`, `second_reads_path` or `second_sequences`, `error_rate` (0.0001 to 0.367), `ksize`.
- **Output:** `first_estimate`, `second_estimate`, `merged_estimate`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_merge_hll_counters

r = khmer_merge_hll_counters(
    first_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"],
    second_sequences=["AGACTACTTATGAGGCCATACCAGGGTCGTTTGCAAAGTCAATAGCAGCCATAGTCCAAC"],
    ksize=15,
    verbose=False,
)
print(r["first_estimate"], r["second_estimate"], r["merged_estimate"])
same = khmer_merge_hll_counters(
    first_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], second_sequences=["TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"], ksize=15, verbose=False
)
print(same["first_estimate"], same["merged_estimate"])
```

- **Example output:**

```text
45 45 90
45 45
```

## 8) `dnakit.microbiome.khmer_align_read()`

- **Purpose:** Align a read to the k-mer graph of a countgraph with khmer's hidden Markov model.
- **Method:** Builds or loads a countgraph and creates `khmer.ReadAligner(graph, trusted_cov_cutoff, bits_theta)`. `forward_only=False` calls `align`, which seeds on the read and extends in both directions; `forward_only=True` calls `align_forward`, which extends forward only. `trusted_cov_cutoff` (default 2) is the abundance from which a k-mer is trusted and `bits_theta` (default 1) the score threshold, both khmer's sandbox values. The returned graph alignment is the path in the graph, so it shows the correction of a sequencing error; gaps are `-`. Covers work-list rows G1276 and G1277.
- **API:** `dnakit.microbiome.khmer_align_read(read, *, forward_only=False, trusted_cov_cutoff=2, bits_theta=1.0, graph_path=None, reads_path=None, sequences=None, ksize=32, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** `read`, `forward_only`, `trusted_cov_cutoff`, `bits_theta`; `graph_path`, `reads_path` or `sequences`; `ksize`, `table_size`, `n_tables`.
- **Output:** `score`, `graph_alignment`, `read_alignment`, `truncated`.
- **Example code:**

```python
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_align_read

genome = "TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGT"
read = genome[5:45]
error = read[:20] + ("A" if read[20] != "A" else "C") + read[21:]
for forward in (False, True):
    r = khmer_align_read(error, sequences=[genome] * 4, ksize=15, forward_only=forward, verbose=False)
    print(forward, round(r["score"], 1), r["graph_alignment"] == read, r["read_alignment"] == error, r["truncated"])
```

- **Example output:**

```text
False 77.2 True True False
True 77.2 True True False
```

## 9) `dnakit.microbiome.khmer_extract_paired_reads()`

- **Purpose:** Separate a mixed read file into interleaved pairs and orphan reads.
- **Method:** Reproduces `extract-paired-reads.py`: reads the file with `khmer.ReadParser`, pairs consecutive reads by name (`/1` and `/2`, or Casava `1:` and `2:`) with `khmer.utils.broken_paired_reader`, writes pairs interleaved to `paired_output_path` and the rest to `single_output_path`, in FASTA or FASTQ as the input is. Output files must not exist. As in the script, a file without any pair is an error. Verified identical to the script. Covers work-list row G1283.
- **API:** `dnakit.microbiome.khmer_extract_paired_reads(reads_path, paired_output_path, single_output_path)`.
- **Input:** `reads_path`, `paired_output_path`, `single_output_path`.
- **Output:** `n_reads`, `n_pairs`, `n_singletons`, the two output paths.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_extract_paired_reads

folder = tempfile.mkdtemp()
pairs, singles = os.path.join(folder, "pe.fa"), os.path.join(folder, "se.fa")
r = khmer_extract_paired_reads("mixed.fa", pairs, singles, verbose=False)
print(r["n_reads"], r["n_pairs"], r["n_singletons"])
print(open(singles).read().strip().replace("\n", " "))
```

- **Example output:**

```text
5 2 1
>s1 TTTTTTTTTTGGAA
```

## 10) `dnakit.microbiome.khmer_split_paired_reads()`

- **Purpose:** Split an interleaved paired-read file into a left-read file and a right-read file.
- **Method:** Reproduces `split-paired-reads.py` with `khmer.utils.broken_paired_reader`: the first read of every pair goes to `first_output_path`, the second to `second_output_path`. Without `orphan_output_path` an unpaired read is an error (khmer's `require_paired`) and the outputs are removed; with it, orphans are written there. Output files must not exist. Verified identical to the script. Covers work-list row G1284.
- **API:** `dnakit.microbiome.khmer_split_paired_reads(reads_path, first_output_path, second_output_path, *, orphan_output_path=None)`.
- **Input:** `reads_path`, `first_output_path`, `second_output_path`, `orphan_output_path`.
- **Output:** `n_reads`, `n_pairs`, `n_orphans`, the output paths.
- **Example code:**

```python
import os
import tempfile

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome.khmer_p27_label_read_functions import khmer_split_paired_reads

folder = tempfile.mkdtemp()
left, right = os.path.join(folder, "left.fa"), os.path.join(folder, "right.fa")
r = khmer_split_paired_reads("paired.fa", left, right, verbose=False)
print(r["n_reads"], r["n_pairs"], r["n_orphans"])
print(open(left).read().strip().replace("\n", " "))
print(open(right).read().strip().replace("\n", " "))
try:
    khmer_split_paired_reads(
        "mixed.fa", os.path.join(folder, "l2.fa"), os.path.join(folder, "r2.fa"), verbose=False
    )
except ConfigurationError:
    print("unpaired read refused")
```

- **Example output:**

```text
4 2 0
>p1/1 ACGTACGTACGGTT >p2/1 AAAACCCCGGTTAC
>p1/2 GGGGCCCCAATTGC >p2/2 CCCCAAAATTGGCA
unpaired read refused
```

## 11) `dnakit.microbiome.khmer_read_statistics()`

- **Purpose:** Number of bases, number of sequences and mean sequence length of FASTA/FASTQ files.
- **Method:** Reproduces `readstats.py` with `khmer.ReadParser`: counts bases and sequences per file and in total; the average length is rounded to one decimal as the script prints it. A file that khmer cannot open or that holds no sequence is reported with zero sequences and the parser message in `error`, and is left out of the total, as in the script. `output_csv` writes the script's `--csv` table (bp, seqs, avg_len, filename). Verified identical to the script. Covers work-list row G1286.
- **API:** `dnakit.microbiome.khmer_read_statistics(reads_paths, *, output_csv=None)`.
- **Input:** `reads_paths` (file or list), `output_csv`.
- **Output:** `files` (`filename`, `bp`, `seqs`, `avg_len`, `error`), `total`, `output_csv`.
- **Example code:**

```python
import os
import tempfile

from dnakit.microbiome.khmer_p27_label_read_functions import khmer_read_statistics

csv_path = os.path.join(tempfile.mkdtemp(), "stats.csv")
r = khmer_read_statistics(["lab.fa", "mixed.fa"], output_csv=csv_path, verbose=False)
for row in r["files"]:
    print(os.path.basename(row["filename"]), row["bp"], row["seqs"], row["avg_len"])
print(r["total"])
print(open(csv_path).read().splitlines()[0])
```

- **Example output:**

```text
lab.fa 120 2 60.0
mixed.fa 70 5 14.0
{'bp': 190, 'seqs': 7, 'avg_len': 27.1}
bp,seqs,avg_len,filename
```
