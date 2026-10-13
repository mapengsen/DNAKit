# Sourmash and khmer k-mer helpers

This group in `dnakit.microbiome` wraps [sourmash](https://github.com/sourmash-bio/sourmash) and [khmer](https://github.com/dib-lab/khmer) at the DNA level. The first eight functions work on sourmash scaled MinHash sketches: pairwise comparison matrices, mutation-rate estimates, set-size error probabilities, k-mer hashes, downsampling, intersection, database search and lowest common lineage. The last eight work on khmer: k-mer hashes and 2-bit codes, low-abundance trimming of local reads, inspection and lookup of saved count tables, table sizing and an expected false-positive rate. They use only the file or sketches passed in. They download no database and do not assemble genomes or classify taxa.

khmer silently returns wrong results for lowercase bases, so these functions uppercase each sequence before khmer sees it, and limit k to 1 to 32. The examples import from the module that defines each function, `dnakit.microbiome` or `dnakit.microbiome`. The dependencies are optional and are imported only when a function runs: sourmash corresponds to `bio-sourmash` and khmer to `bio-khmer`. A missing package raises an error with the install command.

## Example data

The sourmash examples use `sample.fa`. Two genomes differ at one base, and the third record is a short read. The khmer examples use `reads.fa`: three identical reads form a high-coverage region, `tailed` has a novel tail after it, `novel` has k-mers seen once, `nonacgt` contains a non-ACGT symbol, and `short` is shorter than k. These files only demonstrate the interface.

**Example file:** `sample.fa`

```text
>genome_a
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGATCCATGGCA
>genome_b
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGTTCCATGGCA
>read_c
TTTTGGGGCCCCAAAAGGGGTTTTACGTACGTAAAA
```

**Example file:** `reads.fa`

```text
>ref1
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>ref2
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>ref3
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG
>tailed
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGGGGGGGGGGGGGCCCCCCCCCCCC
>novel
CAGTCGATGCATGCAA
>nonacgt
ACGTNNNNACGTACGTACGT
>short
ACGT
```

## 1) `dnakit.microbiome.compare_minhash_matrix()`

- **Purpose:** Compare every pair of scaled MinHash sketches and return a similarity matrix.
- **Method:** `metric="jaccard"` calls `sourmash.compare.compare_serial`. The containment metrics call `compare_serial_containment`, `compare_serial_max_containment` or `compare_serial_avg_containment`. sourmash stores its containment matrix transposed, so this function transposes it: `matrix[i][j]` is the share of sketch i found in sketch j.
- **API:** `dnakit.microbiome.compare_minhash_matrix(sketches, *, metric="jaccard")`.
- **Input:** At least two sketch mappings with `ksize`, `scaled` and `hashes`, and an optional `id`. All sketches share `ksize` and `scaled`. At most 500 sketches.
- **Output:** `ids`, `matrix` (a square list in input order), `metric`, `ksize` and `scaled`.
- **Example code:**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import compare_minhash_matrix

sketches = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][:2]
result = compare_minhash_matrix(sketches)
print(result["ids"], [round(value, 3) for value in result["matrix"][0]])
```

- **Example output:**

```text
['genome_a', 'genome_b'] [1.0, 0.5]
```

## 2) `dnakit.microbiome.estimate_mutation_rate()`

- **Purpose:** Estimate the per-base mutation rate from a sketch Jaccard index or containment.
- **Method:** `measure="jaccard"` calls `sourmash.minhash.jaccard_to_distance`. `measure="containment"` calls `containment_to_distance` and returns an interval. Give the sequence length in bp or the number of distinct k-mers. `ani` is 1 minus `distance`.
- **API:** `dnakit.microbiome.estimate_mutation_rate(value, *, measure="jaccard", ksize=21, scaled=1000, n_unique_kmers=None, sequence_len_bp=None, confidence=0.95)`.
- **Input:** `value` strictly between 0 and 1; at least one of `sequence_len_bp` and `n_unique_kmers`. The estimate is approximate; `size_is_inaccurate` flags a sketch that may be too small.
- **Output:** `distance`, `ani`, `distance_low` and `distance_high` (containment only), `p_nothing_in_common` and `size_is_inaccurate`.
- **Example code:**

```python
from dnakit.microbiome import estimate_mutation_rate

result = estimate_mutation_rate(0.8, ksize=21, scaled=10, n_unique_kmers=1000)
print(round(result["distance"], 4), round(result["ani"], 4), result["size_is_inaccurate"])
```

- **Example output:**

```text
0.0056 0.9944 False
```

## 3) `dnakit.microbiome.set_size_error_probability()`

- **Purpose:** Probability that a scaled sketch estimates a k-mer set size within a relative tolerance.
- **Method:** Calls `sourmash.minhash.set_size_exact_prob`. The sketch size is binomial with parameters `set_size` and `1 / scaled`. The returned value is P(|sketch_size × scaled − set_size| ≤ relative_error × set_size). The upstream docstring describes the complement, but the computed value is the one stated here.
- **API:** `dnakit.microbiome.set_size_error_probability(set_size, scaled, *, relative_error=0.05)`.
- **Input:** `set_size` a positive integer; `scaled` an integer from 1 to 10^9; `relative_error` strictly between 0 and 1.
- **Output:** `expected_sketch_size` (`set_size / scaled`) and `probability_within_tolerance`.
- **Example code:**

```python
from dnakit.microbiome import set_size_error_probability

result = set_size_error_probability(1000, 10, relative_error=0.1)
print(result["expected_sketch_size"], round(result["probability_within_tolerance"], 4))
```

- **Example output:**

```text
100.0 0.7318
```

## 4) `dnakit.microbiome.hash_canonical_kmer()`

- **Purpose:** Canonical MurmurHash of one DNA k-mer, matching the hash a scaled=1 sketch stores.
- **Method:** The canonical form is the smaller of the k-mer and its reverse complement as strings. The hash is `sourmash.minhash.hash_murmur(canonical, seed)`. Lowercase input is uppercased first.
- **API:** `dnakit.microbiome.hash_canonical_kmer(kmer, *, seed=42)`.
- **Input:** A string of A, C, G and T only, 1 to 100 bases long.
- **Output:** `canonical`, `strand` (forward, reverse_complement or self_reverse_complement) and `hash`.
- **Example code:**

```python
from dnakit.microbiome import hash_canonical_kmer

forward = hash_canonical_kmer("acgta")
reverse = hash_canonical_kmer("TACGT")
print(forward["canonical"], reverse["strand"], forward["hash"] == reverse["hash"])
```

- **Example output:**

```text
ACGTA reverse_complement True
```

## 5) `dnakit.microbiome.downsample_minhash_sketch()`

- **Purpose:** Downsample a scaled sketch to a coarser scale.
- **Method:** Calls `sourmash.minhash.flatten_and_downsample_scaled`. The target `scaled` must be at least the sketch's own scale, because sketches cannot be upsampled. Hashes above sourmash's cut are dropped. The cut is computed in floating point, so hashes within a few hundred units of `2**64 / scaled` are decided by sourmash.
- **API:** `dnakit.microbiome.downsample_minhash_sketch(sketch, *, scaled)`.
- **Input:** One sketch mapping. `scaled` must be at least the sketch's scale.
- **Output:** The new `scaled`, `hashes`, `input_hash_count` and `hash_count`.
- **Example code:**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import downsample_minhash_sketch

sketch = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][0]
result = downsample_minhash_sketch(sketch, scaled=100)
print(sketch["hash_count"], result["hash_count"], result["scaled"])
```

- **Example output:**

```text
5 1 100
```

## 6) `dnakit.microbiome.intersect_minhash_sketches()`

- **Purpose:** Hashes shared by two sketches.
- **Method:** Calls `sourmash.minhash.flatten_and_intersect_scaled`: both sketches are downsampled to the coarser scale, then their hashes are intersected.
- **API:** `dnakit.microbiome.intersect_minhash_sketches(left, right)`.
- **Input:** Two sketch mappings with the same `ksize`. Their scales may differ.
- **Output:** `hashes`, `hash_count` and `scaled` (the coarser of the two).
- **Example code:**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import intersect_minhash_sketches

first, second = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"][:2]
result = intersect_minhash_sketches(first, second)
print(result["hash_count"], result["scaled"])
```

- **Example output:**

```text
3 10
```

## 7) `dnakit.microbiome.search_minhash_database()`

- **Purpose:** Find database sketches that share at least `threshold_bp` of k-mers with a query sketch.
- **Method:** Loads the database into a `sourmash.index.LinearIndex` and calls `sourmash.search.prefetch_database`. Overlap in bp is shared hashes times the coarser scale. Matches are sorted by overlap, largest first, then by id. A threshold above the query's own bp is unattainable, and sourmash raises an error, which this function reports as ConfigurationError.
- **API:** `dnakit.microbiome.search_minhash_database(query, database, *, threshold_bp=5000)`.
- **Input:** A query sketch and a list of database sketches with the same `ksize`. Sketches with no hashes are rejected, so filter them first.
- **Output:** A `matches` list. Each item has `id`, `md5`, `intersect_bp`, `jaccard`, `fraction_of_query_in_match`, `fraction_of_match_in_query` and `max_containment`.
- **Example code:**

```python
from dnakit.microbiome import sketch_fasta_minhash
from dnakit.microbiome import search_minhash_database

sketches = [s for s in sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"]
            if s["hash_count"]]
result = search_minhash_database(sketches[0], sketches, threshold_bp=30)
print([(match["id"], match["intersect_bp"]) for match in result["matches"]])
```

- **Example output:**

```text
[('genome_a', 50), ('genome_b', 30)]
```

## 8) `dnakit.microbiome.lowest_common_lineage()`

- **Purpose:** Lowest common ancestor of a set of taxonomic lineages.
- **Method:** Parses with `sourmash.lca.lca_utils.make_lineage`, builds a tree with `build_tree` and finds the lowest shared node with `find_lca`. `distinct_children_below_lca` is 0 when every lineage agrees to full depth; otherwise it counts the different taxa at the next rank.
- **API:** `dnakit.microbiome.lowest_common_lineage(lineages)`.
- **Input:** A list of `;`-separated lineage strings from superkingdom down, at most 8 ranks. Commas are rejected.
- **Output:** `lineage` (rank and name per level), `lineage_string`, `lowest_rank` and `distinct_children_below_lca`.
- **Example code:**

```python
from dnakit.microbiome import lowest_common_lineage

result = lowest_common_lineage(
    [
        "d__Bacteria;p__Firmicutes;c__Bacilli",
        "d__Bacteria;p__Firmicutes;c__Clostridia",
    ]
)
print(result["lineage_string"], result["distinct_children_below_lca"])
```

- **Example output:**

```text
d__Bacteria;p__Firmicutes 2
```

## 9) `dnakit.microbiome.hash_kmer_murmur3()`

- **Purpose:** MurmurHash3 of a DNA k-mer, computed by khmer.
- **Method:** `canonical=True` calls `khmer.hash_murmur3`, which gives the same value for a k-mer and its reverse complement. `canonical=False` calls `khmer.hash_no_rc_murmur3`. khmer hashes lowercase bases silently wrong, so the input is uppercased first. These values differ from sourmash hashes.
- **API:** `dnakit.microbiome.hash_kmer_murmur3(kmer, *, canonical=True)`.
- **Input:** A string of A, C, G and T only, 1 to 32 bases long.
- **Output:** `hash` (unsigned 64-bit integer), `canonical` and `ksize`.
- **Example code:**

```python
from dnakit.microbiome import hash_kmer_murmur3

print(hash_kmer_murmur3("acgta")["hash"] == hash_kmer_murmur3("TACGT")["hash"], hash_kmer_murmur3("ACGTA")["hash"])
```

- **Example output:**

```text
True 7467422729906674866
```

## 10) `dnakit.microbiome.encode_kmer_2bit()`

- **Purpose:** Exact 2-bit integer code of a DNA k-mer (k up to 32).
- **Method:** `canonical=True` calls `khmer.forward_hash` and keeps the smaller code of the k-mer and its reverse complement. `canonical=False` calls `khmer.forward_hash_no_rc`. khmer's base order is A=0, T=1, C=2, G=3, with the first base in the high bits.
- **API:** `dnakit.microbiome.encode_kmer_2bit(kmer, *, canonical=True)`.
- **Input:** A string of A, C, G and T only, 1 to 32 bases long.
- **Output:** `code`, `canonical` and `ksize`.
- **Example code:**

```python
from dnakit.microbiome import encode_kmer_2bit

print(
    encode_kmer_2bit("ACGT")["code"],
    encode_kmer_2bit("GGGGA")["code"],
    encode_kmer_2bit("GGGGA", canonical=False)["code"],
)
```

- **Example output:**

```text
45 426 1020
```

## 11) `dnakit.microbiome.decode_kmer_2bit()`

- **Purpose:** Recover the k-mer from a 2-bit code produced by `encode_kmer_2bit()`.
- **Method:** Calls `khmer.reverse_hash`. A canonical code decodes to the canonical k-mer. The value must lie between 0 and 4^k − 1.
- **API:** `dnakit.microbiome.decode_kmer_2bit(value, *, ksize)`.
- **Input:** Integer `value` and `ksize` (1 to 32).
- **Output:** The uppercase `kmer`.
- **Example code:**

```python
from dnakit.microbiome import decode_kmer_2bit

print(decode_kmer_2bit(45, ksize=4)["kmer"])
```

- **Example output:**

```text
ACGT
```

## 12) `dnakit.microbiome.trim_low_abundance_reads()`

- **Purpose:** Low-abundance trimming of local reads: cut each read at the first k-mer whose count is below `cutoff`.
- **Method:** Pass one counts canonical k-mers from all usable reads with a khmer `Countgraph`. Pass two calls `khmer.trimming.trim_record` on each read, keeping the bases before the first low-count k-mer (counted from the 5' end). A read left shorter than k is dropped. Output is FASTA (`.fa`, `.fasta`, `.fna`) or FASTQ (`.fq`, `.fastq`, which needs FASTQ input). Output bases are uppercase.
- **API:** `dnakit.microbiome.trim_low_abundance_reads(source, output, *, ksize=21, cutoff=2, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** A FASTA or FASTQ file. The output path must not exist. Reads with symbols other than A, C, G and T, and reads shorter than k, are skipped and counted.
- **Output:** Counts `reads`, `used_reads`, `kept_reads`, `trimmed_reads`, `dropped_short_reads` and others, plus the output path.
- **Example code:**

```python
from pathlib import Path
from dnakit.microbiome import trim_low_abundance_reads

Path("trimmed.fa").unlink(missing_ok=True)
result = trim_low_abundance_reads("reads.fa", "trimmed.fa", ksize=9, cutoff=2)
print(
    result["reads"],
    result["used_reads"],
    result["trimmed_reads"],
    result["kept_reads"],
    result["dropped_short_reads"],
)
```

- **Example output:**

```text
7 5 1 4 1
```

## 13) `dnakit.microbiome.inspect_khmer_graph()`

- **Purpose:** Read the header of a saved khmer count table or presence table without loading it.
- **Method:** `kind="countgraph"` calls `khmer.extract_countgraph_info`; `kind="nodegraph"` calls `khmer.extract_nodegraph_info` and maps the tuple to named fields. khmer does not check the file type, so this function checks the type flag and rejects the other kind.
- **API:** `dnakit.microbiome.inspect_khmer_graph(source, *, kind="countgraph")`.
- **Input:** A `.cg` or `.ng` file saved by khmer.
- **Output:** `ksize`, `table_size`, `n_tables`, `n_occupied`, `format_version` and `table_type_flag`.
- **Example code:**

```python
import khmer
from dnakit.microbiome import inspect_khmer_graph

table = khmer.Countgraph(7, 1009, 3)
table.consume("ACGTTGCAAGCTAGCTTTG")
table.save("counts.cg")
info = inspect_khmer_graph("counts.cg")
print(info["ksize"], info["n_tables"], info["table_size"], info["n_occupied"])
```

- **Example output:**

```text
7 3 1000 11
```

## 14) `dnakit.microbiome.query_saved_countgraph()`

- **Purpose:** Look up canonical k-mer counts in a saved khmer count table.
- **Method:** Loads the table with `khmer.load_countgraph` and calls `get()` for each query. A k-mer and its reverse complement share one count, and a hash collision can only raise a count. Query k-mers must have the table's ksize.
- **API:** `dnakit.microbiome.query_saved_countgraph(source, kmers)`.
- **Input:** A khmer count-table file and a list of equal-length k-mers (at most 100000).
- **Output:** `counts` sorted by k-mer, plus `ksize` and `hashsizes`.
- **Example code:**

```python
import khmer
from dnakit.microbiome import query_saved_countgraph

table = khmer.Countgraph(7, 1009, 3)
table.consume("ACGTTGCAAGCTAGCTTTG")
table.save("counts.cg")
result = query_saved_countgraph("counts.cg", ["ACGTTGC", "GCAACGT"])
print(result["counts"])
```

- **Example output:**

```text
{'ACGTTGC': 1, 'GCAACGT': 1}
```

## 15) `dnakit.microbiome.estimate_kmer_table_size()`

- **Purpose:** Table layout for a number of distinct k-mers, from a memory cap or a false-positive rate.
- **Method:** Calls `khmer.khmer_args.optimal_size`. `mem_cap` and `mem_use` are in khmer's table units (total hash bins over all tables); khmer does not convert them to bytes.
- **API:** `dnakit.microbiome.estimate_kmer_table_size(num_kmers, *, mem_cap=None, fp_rate=None)`.
- **Input:** `num_kmers` a positive integer; give exactly one of `mem_cap` and `fp_rate`.
- **Output:** `num_htables`, `htable_size`, `mem_use` and `expected_fp_rate`.
- **Example code:**

```python
from dnakit.microbiome import estimate_kmer_table_size

result = estimate_kmer_table_size(1_000_000, fp_rate=0.01)
print(result["num_htables"], result["htable_size"], round(result["expected_fp_rate"], 4))
```

- **Example output:**

```text
6 1602775 0.01
```

## 16) `dnakit.microbiome.expected_kmer_table_collisions()`

- **Purpose:** Expected false-positive rate of a khmer presence table filled with the reads of a file.
- **Method:** Calls `khmer.calc_expected_collisions`, which divides occupied bins by the smallest table and raises the result to the number of tables. It is a quick approximation. The upstream function prints an error and exits above its threshold; this wrapper does neither and reports `within_threshold`.
- **API:** `dnakit.microbiome.expected_kmer_table_collisions(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000, max_false_positive=0.2)`.
- **Input:** A FASTA or FASTQ file. `max_false_positive` is greater than 0 and at most 1.
- **Output:** `hashsizes`, `occupied_bins`, `expected_false_positive_rate`, `within_threshold` and read counts.
- **Example code:**

```python
from dnakit.microbiome import expected_kmer_table_collisions

result = expected_kmer_table_collisions("reads.fa", ksize=9, table_size=10007, n_tables=2)
print(result["within_threshold"], result["occupied_bins"], result["hashsizes"])
```

- **Example output:**

```text
True 40 [9973, 9967]
```
