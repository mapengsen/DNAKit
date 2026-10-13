# G18b functions: sourmash, khmer and primer3

This page documents the G18b wrappers. They are exported from `dnakit.microbiome` and `dnakit.molecular_lab`. For example, `from dnakit.microbiome import kmer_mutation_probability`.

- Local functions: all 9 wrappers run offline. Their inputs are numbers, local DNA strings, lineage strings or a local FASTA or FASTQ file. No function downloads data or uses the network, so this group has no live test.
- khmer inputs are uppercased before khmer sees them. None of the khmer functions in this group takes a k-mer length, so the ksize check (1 to 32) does not apply here.
- Input files are declared with the `Example file` label; the other examples run as written.

## Existing and unwrapped functions

- Ledger counts for `other` and `cap` rows, before and after this review. There were no `cap` rows.
  - sourmash: 121 rows become 6 new wrappers (14 `wrapped` in total). The other 115 rows are 29 internal helpers and constants; 25 command-line or report helpers; 22 family variants (signature loaders and savers, JSON-string parsers, search drivers, nodegraph collision and header helpers); 20 taxonomy or picklist helpers with no sequence meaning; 2 duplicates of the new wrappers; and 17 that need a prebuilt file (10 LCA databases, 3 SBT indexes, 2 SQLite indexes and 2 sets of gather CSV files).
  - khmer: 4 rows become 2 wrappers (`read_paired_fastx` and `clean_khmer_read_sequences`). The other 2 are `get_n_primes_near_x` and `is_prime`, which only size hash tables.
  - primer3: 4 rows become 1 wrapper (`validate_acgt_sequence`). `ensure_acgt_uppercase_b` is its bytes form (family variant). `includes` and `get_libprimer3_version` are packaging and provenance helpers.
- Several sibling functions are variants of existing wrappers. `search_minhash_database` covers the containment search; the search drivers (`search_databases_with_flat_query`, `search_databases_with_abund_query`) and the SBT tree searches (`search_sbt_index_private`, `sbtmh.search_sbt_index`) run the same kind of search with another index or score. `read_sourmash_signatures` reads the same `.sig` files as the signature loaders, and `expected_kmer_table_collisions` covers the nodegraph collision estimate.
- Not wrapped for lack of a file: the LCA database functions (`summarize`, `gather_assignments`, `classify_signature`, `make_lca_counts` and the loaders), the SBT and SQLite index loaders, and the gather CSV loaders. Each needs a prebuilt file that the wrapper does not create.
- The full status for every inventory key is in the G18b ledger.

## 1) `kmer_mutation_probability`

- **Purpose:** Probability that a k-mer contains at least one mutated base, given a per-base mutation rate.
- **Method:** Calls sourmash's `distance_utils.r1_to_q`, which computes `1 - (1 - rate) ** ksize`. The bases are assumed to mutate independently.
- **API:** `kmer_mutation_probability(mutation_rate[required], ksize[required], verbose[optional])`
- **Input:** `mutation_rate` is a number in [0, 1]; `ksize` is an integer from 1 to 100.
- **Output:** A dictionary with `ksize`, `mutation_rate`, `kmer_mutation_probability` and `provenance`.

**Example code:**

```python
from dnakit.microbiome import kmer_mutation_probability

for k in (1, 2, 21):
    result = kmer_mutation_probability(0.01, k, verbose=False)
    print(k, round(result["kmer_mutation_probability"], 6))
```

**Example output:**

```text
1 0.01
2 0.0199
21 0.190272
```

## 2) `expected_mutated_kmers`

- **Purpose:** Expected number of distinct k-mers of a sequence that a per-base mutation rate changes.
- **Method:** Calls sourmash's `distance_utils.exp_n_mutated` on the number of distinct k-mers. The result is that count times the per-k-mer mutation probability. For a sequence of length L without repeats, the number of distinct k-mers is `L - ksize + 1`; repeats are not collapsed.
- **API:** `expected_mutated_kmers(ksize[required], mutation_rate[required], sequence_len_bp[optional], n_unique_kmers[optional], verbose[optional])`
- **Input:** Give exactly one of `sequence_len_bp` (at least `ksize`) and `n_unique_kmers` (a positive integer). `mutation_rate` is in [0, 1].
- **Output:** A dictionary with `ksize`, `mutation_rate`, `n_unique_kmers`, `kmer_mutation_probability`, `expected_mutated_kmers`, `expected_unmutated_kmers` and `provenance`. The values are expectations, not counts.

**Example code:**

```python
from dnakit.microbiome import expected_mutated_kmers

result = expected_mutated_kmers(21, 0.01, sequence_len_bp=5000, verbose=False)
print(result["n_unique_kmers"], round(result["expected_mutated_kmers"], 2))
```

**Example output:**

```text
4980 947.56
```

## 3) `variance_mutated_kmers`

- **Purpose:** Variance of the number of mutated distinct k-mers, from the same mutation model.
- **Method:** Calls sourmash's `distance_utils.var_n_mutated`, the closed-form variance of the mutated k-mer count. A rate of 0 gives 0. For very small rates the closed form can go negative in floating point; sourmash raises for that case, and this wrapper reports it as a configuration error.
- **API:** `variance_mutated_kmers(ksize[required], mutation_rate[required], sequence_len_bp[optional], n_unique_kmers[optional], verbose[optional])`
- **Input:** The same as `expected_mutated_kmers`.
- **Output:** A dictionary with `ksize`, `mutation_rate`, `n_unique_kmers`, `variance_mutated_kmers`, `standard_deviation_mutated_kmers` and `provenance`.

**Example code:**

```python
from dnakit.microbiome import variance_mutated_kmers

result = variance_mutated_kmers(2, 0.1, sequence_len_bp=11, verbose=False)
print(round(result["variance_mutated_kmers"], 4))
```

**Example output:**

```text
2.8512
```

## 4) `probability_no_shared_kmers`

- **Purpose:** Probability that a scaled sketch of a sequence and a sketch of its mutated copy share no hash.
- **Method:** Calls sourmash's `distance_utils.get_exp_probability_nothing_common`, which evaluates `exp((n - E[mutated]) * ln(1 - 1/scaled))`. A rate of 0 gives 0 and a rate of 1 gives 1. This is an expectation, so a large value suggests the two sketches may have nothing in common.
- **API:** `probability_no_shared_kmers(mutation_rate[required], ksize[required], scaled[required], sequence_len_bp[optional], n_unique_kmers[optional], verbose[optional])`
- **Input:** `scaled` is an integer of at least 2; a scale of 1 keeps every k-mer and the logarithm is undefined. Give exactly one of `sequence_len_bp` and `n_unique_kmers`.
- **Output:** A dictionary with `mutation_rate`, `ksize`, `scaled`, `n_unique_kmers`, `probability_nothing_in_common` and `provenance`.

**Example code:**

```python
from dnakit.microbiome import probability_no_shared_kmers

for rate in (0.0, 0.01, 0.1):
    result = probability_no_shared_kmers(rate, 21, 1000, sequence_len_bp=5000, verbose=False)
    print(rate, f"{result['probability_nothing_in_common']:.3e}")
```

**Example output:**

```text
0.0 0.000e+00
0.01 1.770e-02
0.1 5.797e-01
```

## 5) `bp_threshold_to_hashes`

- **Purpose:** Convert a bp overlap threshold into a hash count and a containment fraction of the query sketch.
- **Method:** Calls sourmash's `search.calc_threshold_from_bp`. The hash count is `threshold_bp / scaled`, and the containment fraction is that count divided by `query_size`. A threshold above the query's own bp can never be reached, and sourmash raises for it; this wrapper reports that as a configuration error.
- **API:** `bp_threshold_to_hashes(threshold_bp[required], scaled[required], query_size[required], verbose[optional])`
- **Input:** `threshold_bp` is a non-negative integer (bp); `scaled` is a positive integer; `query_size` is the number of hashes in the query sketch.
- **Output:** A dictionary with `threshold_bp`, `scaled`, `query_size`, `threshold_hashes` (may be fractional), `containment_threshold` (in [0, 1]) and `provenance`.

**Example code:**

```python
from dnakit.microbiome import bp_threshold_to_hashes

result = bp_threshold_to_hashes(50_000, 1000, 200, verbose=False)
print(result["threshold_hashes"], result["containment_threshold"])
```

**Example output:**

```text
50.0 0.25
```

## 6) `count_kmer_lca_assignments`

- **Purpose:** For each k-mer hash, find the lowest common ancestor of its taxonomic assignments, then count how many hashes share each ancestor.
- **Method:** Calls sourmash's `lca.lca_utils.count_lca_for_assignments`. Each hash contributes its lowest common ancestor, or one unit of weight when no weights are given. `lowest_common_lineage` covers a single list of lineages; this wrapper does the same for every hash of a sample.
- **API:** `count_kmer_lca_assignments(assignments[required], hashval_counts[optional], verbose[optional])`
- **Input:** `assignments` maps a 64-bit hash to a non-empty list of `;`-separated lineage strings, from superkingdom down, with at most 8 ranks and no commas. `hashval_counts`, if given, maps exactly the same hashes to positive integer weights.
- **Output:** A dictionary with `hash_count`, `weighted`, `total_count` and `counts`. Each entry of `counts` has `lineage` (an empty string is the root), `depth`, `lowest_rank` and `count`. The entries are sorted by count (largest first), then by lineage.

**Example code:**

```python
from dnakit.microbiome import count_kmer_lca_assignments

assignments = {
    101: ["d__Bacteria;p__Firmicutes;c__Bacilli", "d__Bacteria;p__Firmicutes;c__Clostridia"],
    102: ["d__Bacteria;p__Bacteroidetes"],
    103: ["d__Archaea", "d__Bacteria"],
}
result = count_kmer_lca_assignments(assignments, verbose=False)
for item in result["counts"]:
    print(repr(item["lineage"]), item["lowest_rank"], item["count"])
```

**Example output:**

```text
'' None 1
'd__Bacteria;p__Bacteroidetes' phylum 1
'd__Bacteria;p__Firmicutes' phylum 1
```

## 7) `clean_khmer_read_sequences`

- **Purpose:** Uppercase reads and replace N with A, as khmer's own cleaning does, and report how many Ns were replaced.
- **Method:** Calls `khmer.utils.clean_input_reads`, which sets `cleaned_seq` to the uppercase read with N replaced by A. The replacement is lossy, so the cleaned read is not the original sequence. Each read reports its `replaced_n_count`. dnakit's own khmer counting and trimming wrappers skip reads that contain N instead of replacing it.
- **API:** `clean_khmer_read_sequences(sequences[required], verbose[optional])`
- **Input:** A mapping of read ID to DNA string, using A, C, G, T and N in either case. At most 10000 reads, each at most 100000 bases; other symbols raise an error.
- **Output:** A dictionary with `read_count`, `reads_with_n`, `replaced_n_total` and `reads`, sorted by ID. Each read has `id`, `sequence` (uppercase), `cleaned_sequence` and `replaced_n_count`.

**Example code:**

```python
from dnakit.microbiome import clean_khmer_read_sequences

result = clean_khmer_read_sequences({"read_b": "ACGT", "read_a": "acgNNt"}, verbose=False)
for item in result["reads"]:
    print(item["id"], item["sequence"], item["cleaned_sequence"], item["replaced_n_count"])
print(result["replaced_n_total"])
```

**Example output:**

```text
read_a ACGNNT ACGAAT 2
read_b ACGT ACGT 0
2
```

## 8) `read_paired_fastx`

- **Purpose:** Pair the records of a local FASTA or FASTQ file, keeping unpaired records as singletons.
- **Method:** Reads the file with `khmer.ReadParser` and passes the records to `khmer.utils.broken_paired_reader`. Two consecutive records form a pair when their names match as mates, such as `name/1` and `name/2`, or in the `1:`/`2:` and `seq/1`/`seq/2` forms. Sequences are uppercased.
- **API:** `read_paired_fastx(source[required], min_length[optional], require_paired[optional], force_single[optional], max_reads[optional], verbose[optional])`
- **Input:** A local FASTA or FASTQ file. `min_length` drops reads shorter than it; `require_paired=True` turns any orphan into an error; `force_single=True` returns every record as a singleton; `max_reads` (default 100000) bounds the number of records. `require_paired` and `force_single` cannot both be true.
- **Output:** A dictionary with `unit_count`, `pair_count`, `singleton_count` and `units`, in file order. Each unit has `index` (counted as khmer counts records, advancing by 2 for a pair), `is_pair`, `read1` and `read2` (`None` for a singleton). Each read has `name` and `sequence`.

**Example file:** `pairs.fq`

```text
@r1/1
ACGTACGT
+
IIIIIIII
@r1/2
TTGGCCAA
+
IIIIIIII
@r2/1
GGGGAAAA
+
IIIIIIII
@r3/1
CCCCTTTT
+
IIIIIIII
@r3/2
AAAACCCC
+
IIIIIIII
```

**Example code:**

```python
from dnakit.microbiome import read_paired_fastx

result = read_paired_fastx("pairs.fq", verbose=False)
for unit in result["units"]:
    names = [unit["read1"]["name"]] + ([unit["read2"]["name"]] if unit["read2"] else [])
    print(unit["index"], unit["is_pair"], " ".join(names))
print(result["pair_count"], result["singleton_count"])
```

**Example output:**

```text
0 True r1/1 r1/2
2 False r2/1
3 True r3/1 r3/2
2 1
```

## 9) `validate_acgt_sequence`

- **Purpose:** Uppercase a DNA sequence and require it to contain only A, C, G and T, as primer design input does.
- **Method:** Calls primer3-py's `p3helpers.ensure_acgt_uppercase`. Unlike dnakit's normalised alphabet check, it does not convert U to T or strip whitespace, so those symbols are refused rather than changed.
- **API:** `validate_acgt_sequence(sequence[required], verbose[optional])`
- **Input:** A non-empty string of at most 100000 characters. Lowercase a, c, g and t are accepted.
- **Output:** A dictionary with `length`, `sequence` (uppercase) and `provenance`, which records the primer3-py version and its GPL-2.0 license.

**Example code:**

```python
from dnakit.molecular_lab import validate_acgt_sequence

print(validate_acgt_sequence("acgtacgt", verbose=False)["sequence"])
try:
    validate_acgt_sequence("ACGU", verbose=False)
except Exception as exc:
    print(type(exc).__name__)
```

**Example output:**

```text
ACGTACGT
ConfigurationError
```
