# Optional MinHash and k-mer backends

This group of functions in `dnakit.microbiome` is a set of adapters over optional packages for DNA sequences. [sourmash](https://github.com/sourmash-bio/sourmash) computes scaled MinHash sketches, Jaccard similarity and containment, and reads or writes sourmash signature files. [khmer](https://github.com/dib-lab/khmer) provides canonical k-mer counting, unique k-mer estimates, abundance histograms and digital normalization. These functions read only the explicit FASTA or FASTQ file they are given. They download no database and do not classify taxa or assemble genomes. Read counts and base counts have explicit upper limits.

The dependencies are optional and are imported only when the matching function is called. A missing package raises an error with the install command. sourmash corresponds to `bio-sourmash`, and khmer corresponds to `bio-khmer`.

## Example data

All examples below use this FASTA. Two genomes differ at one base, and the third record is a short read. It only demonstrates the interface, and its statistics do not describe a real microbial community.

**Example file:** `sample.fa`

```text
>genome_a
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGATCCATGGCA
>genome_b
ACGTTGCAAGCTAGCTTTGCAAGGCCTAGGTTCCATGGCA
>read_c
TTTTGGGGCCCCAAAAGGGGTTTTACGTACGTAAAA
```

## 1) `dnakit.microbiome.sketch_fasta_minhash()`

- **Purpose:** Sketch each FASTA record as a scaled MinHash.
- **Method:** Call `sourmash.MinHash(n=0, ksize, scaled)` and `add_sequence(..., force=True)`. k-mers are canonical, so a k-mer and its reverse complement share one hash. sourmash skips k-mers that contain symbols other than A, C, G and T. The FASTA parser is pure Python and has explicit limits.
- **API:** `dnakit.microbiome.sketch_fasta_minhash(source, *, ksize=31, scaled=1000, max_records=1000, max_bases=10000000)`.
- **Input:** A FASTA file. `ksize` is 1 to 100. `scaled` is an integer from 1 to 10^9; smaller values keep more hashes. Exceeding `max_records` or `max_bases` raises an error.
- **Output:** A dictionary with `records` and `sketches`. Each sketch has `id`, `length`, `ksize`, `scaled`, `hash_count`, and sorted `hashes` (unsigned 64-bit integers).
- **Example code:**

```python
from dnakit.microbiome import sketch_fasta_minhash

result = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)
print(result["records"], [(s["id"], s["hash_count"]) for s in result["sketches"]])
```

- **Example output:**

```text
3 [('genome_a', 5), ('genome_b', 4), ('read_c', 0)]
```

## 2) `dnakit.microbiome.compare_minhash_sketches()`

- **Purpose:** Compare two sketches by Jaccard similarity and containment.
- **Method:** Use sourmash's `MinHash.similarity(ignore_abundance=True)` for Jaccard and `contained_by()` for containment in both directions. Both sketches must share `ksize` and `scaled`. Containment is directional: `left_contained_by_right` equals |L∩R| / |L|.
- **API:** `dnakit.microbiome.compare_minhash_sketches(left, right)`.
- **Input:** Two sketch dictionaries, either from `sketch_fasta_minhash()` or as items of the `sketches` list from `read_sourmash_signatures()`. Hashes cannot be empty and must lie within the `scaled` range.
- **Output:** Jaccard similarity, containment in both directions, shared hash count, and hash counts.
- **Example code:**

```python
from dnakit.microbiome import compare_minhash_sketches, sketch_fasta_minhash

sketches = sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"]
result = compare_minhash_sketches(sketches[0], sketches[1])
print(round(result["jaccard"], 3), round(result["left_contained_by_right"], 3))
```

- **Example output:**

```text
0.5 0.603
```

## 3) `dnakit.microbiome.write_sourmash_signatures()`

- **Purpose:** Write sketches as one sourmash JSON signature file.
- **Method:** Wrap each sketch in `sourmash.signature.SourmashSignature`, then write with `sourmash.save_signatures_to_json()`. A name ending in `.gz` is compressed. An existing file is never overwritten.
- **API:** `dnakit.microbiome.write_sourmash_signatures(sketches, output)`.
- **Input:** A list of sketches, each with at least one hash. `output` must end with `.sig` or `.sig.gz`, and its directory must exist.
- **Output:** The written file path and the number of signatures.
- **Example code:**

```python
from pathlib import Path
from dnakit.microbiome import sketch_fasta_minhash, write_sourmash_signatures

Path("sketches.sig.gz").unlink(missing_ok=True)
sketches = [s for s in sketch_fasta_minhash("sample.fa", ksize=21, scaled=10)["sketches"] if s["hash_count"]]
result = write_sourmash_signatures(sketches, "sketches.sig.gz")
print(result["signatures"], Path(result["output"]).name)
```

- **Example output:**

```text
2 sketches.sig.gz
```

## 4) `dnakit.microbiome.read_sourmash_signatures()`

- **Purpose:** Read a sourmash signature file and return sketches with the same structure as `sketch_fasta_minhash()`.
- **Method:** Load with `sourmash.load_file_as_signatures()`. Only DNA scaled signatures are accepted; protein signatures and num-based sketches raise an error. Abundance values are dropped, and each signature's `md5sum` is recorded.
- **API:** `dnakit.microbiome.read_sourmash_signatures(source, *, max_signatures=1000)`.
- **Input:** A `.sig` or `.sig.gz` file. More signatures than `max_signatures` raise an error.
- **Output:** A `sketches` list. `length` is `None` because signatures do not store sequence length.
- **Example code:**

```python
from dnakit.microbiome import read_sourmash_signatures

result = read_sourmash_signatures("sketches.sig.gz")
print([(s["id"], s["hash_count"]) for s in result["sketches"]])
```

The example reads the `sketches.sig.gz` written in the previous section, so run the two examples in order.

- **Example output:**

```text
[('genome_a', 5), ('genome_b', 4)]
```

## 5) `dnakit.microbiome.query_kmer_counts()`

- **Purpose:** Count how often each query k-mer occurs in the reads.
- **Method:** Count with a khmer `Countgraph`. k-mers are canonical, so a k-mer and its reverse complement share one count. Reads are uppercased first. Reads containing symbols other than A, C, G and T, and reads shorter than k, are skipped and counted. The count table is approximate; hash collisions can only raise counts.
- **API:** `dnakit.microbiome.query_kmer_counts(source, kmers, *, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** A FASTA or FASTQ file. `kmers` is a list of uppercase ACGT strings of one length from 1 to 32.
- **Output:** The count for each query k-mer, plus the numbers of reads read, used and skipped.
- **Example code:**

```python
from dnakit.microbiome import query_kmer_counts

result = query_kmer_counts("sample.fa", ["ACGTTG", "TTTTGG"])
print(result["counts"], result["used_reads"])
```

- **Example output:**

```text
{'ACGTTG': 2, 'TTTTGG': 2} 3
```

## 6) `dnakit.microbiome.estimate_unique_kmers()`

- **Purpose:** Estimate the number of distinct canonical k-mers in a file.
- **Method:** Record k-mers in a khmer `Nodegraph` (a presence table). The distinct count comes from `n_unique_kmers()`, which uses the number of occupied slots. The result is approximate, and a larger table brings it closer to the exact count.
- **API:** `dnakit.microbiome.estimate_unique_kmers(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** A FASTA or FASTQ file. `ksize` is 1 to 32.
- **Output:** Total k-mer positions, the distinct-count estimate, occupied slots, and hash table sizes.
- **Example code:**

```python
from dnakit.microbiome import estimate_unique_kmers

result = estimate_unique_kmers("sample.fa", ksize=21, table_size=1_000_000)
print(result["kmer_positions"], result["unique_kmers_estimate"])
```

- **Example output:**

```text
56 46
```

## 7) `dnakit.microbiome.kmer_abundance_histogram()`

- **Purpose:** Count how many distinct canonical k-mers occur each number of times.
- **Method:** All reads are first counted into a `Countgraph`. Then khmer's `abundance_distribution()` scans the same file. khmer reads the file again, so every read must be uppercase A, C, G and T, and at least `ksize` long. Otherwise the function raises an error instead of skipping.
- **API:** `dnakit.microbiome.kmer_abundance_histogram(source, *, ksize=21, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** A FASTA or FASTQ file. `ksize` is 1 to 32.
- **Output:** `distinct_kmers` and a `histogram` list whose items have `count` and `kmers`.
- **Example code:**

```python
from dnakit.microbiome import kmer_abundance_histogram

result = kmer_abundance_histogram("sample.fa", ksize=11)
print(result["distinct_kmers"], result["histogram"])
```

- **Example output:**

```text
60 [{'count': 1, 'kmers': 34}, {'count': 2, 'kmers': 26}]
```

## 8) `dnakit.microbiome.normalize_reads_by_median()`

- **Purpose:** Digital normalization: discard reads whose coverage is already high.
- **Method:** Visit each read in file order. Take the median canonical k-mer count with `get_median_count()`. Keep the read and count its k-mers when the median is below `cutoff`; otherwise discard it. The output is FASTA with only the header and sequence.
- **API:** `dnakit.microbiome.normalize_reads_by_median(source, output, *, ksize=20, cutoff=20, table_size=1000000, n_tables=4, max_reads=1000000, max_bases=100000000)`.
- **Input:** A FASTA or FASTQ file, and an output path ending in `.fa`, `.fasta` or `.fna`. An existing file is never overwritten.
- **Output:** Counts of reads read, kept, discarded for high median, skipped as non-ACGT, and skipped as short.
- **Example code:**

```python
from pathlib import Path
from dnakit.microbiome import normalize_reads_by_median

Path("normalized.fa").unlink(missing_ok=True)
result = normalize_reads_by_median("sample.fa", "normalized.fa", ksize=11, cutoff=1)
print(result["kept_reads"], result["discarded_high_median_reads"])
```

- **Example output:**

```text
2 1
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-sourmash]"    # sourmash: sketches, comparison, signature files
python -m pip install "dnakit[bio-khmer]"       # khmer: k-mer counts, estimates, histograms, normalization
```

Project pages: [sourmash](https://github.com/sourmash-bio/sourmash), [khmer](https://github.com/dib-lab/khmer).

## How these functions are verified

The test file is `others/tests/unit/test_microbiome_backends.py`, with 19 tests. They check the following:

- sourmash sketches match direct calls to `sourmash.MinHash`, and the canonical k-mer hash count is checked by hand (`ACGTACGT` gives 2 hashes for 3-mers).
- Jaccard similarity and containment match hand-counted sets and direct sourmash calls; written signature files are read back by sourmash unchanged.
- khmer k-mer counts match exact canonical counts computed in Python; the abundance histogram matches exact counts; the hand-worked normalization example keeps 2 reads and discards 1.
- The unique k-mer estimate matches a direct khmer call and lies within 1% of the exact count.
- A missing sourmash or khmer raises an install hint. Invalid ksize, read limits, non-ACGT input and malformed sketches are rejected.
- The rule that rejects protein signatures is tested with a stand-in loader, because sourmash could not reload very short protein signatures in this environment.

## Current scope

This batch covers sourmash scaled MinHash sketches, comparison and signature file I/O, and khmer k-mer counts, unique k-mer estimates, abundance histograms and digital normalization.

- Not wrapped: sourmash `search`, `gather`, `taxonomy` and LCA/SBT indexes. They need prebuilt databases, which is beyond a file-level interface.
- Not wrapped: num-based and protein sketches, and abundance-weighted comparison (`track_abundance`). The DNA workflow uses scaled sketches, and these functions return set relations only.
- Not wrapped: khmer assemblers (`LinearAssembler`, `JunctionCountAssembler`), stop-tag trimming, and graph save/load. They produce contigs or modify graph files, which is outside a counting interface.
- Not wrapped: FASTQ quality trimming and paired-read handling.
- `dnakit.integrations.compare_sourmash_sequences()` already compares two strings. `compare_minhash_sketches()` compares sketches read from files, so the two have different uses.
