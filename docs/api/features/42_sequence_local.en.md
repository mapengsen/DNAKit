# Dependency-free sequence primitives

These functions operate on in-memory DNA/RNA strings. They do not parse FASTA/FASTQ, align reads, or infer sequencing quality, so they can be reported as deterministic preprocessing or descriptive-statistics steps in a methods section.

## 1) Calculate sequence metrics

- **Purpose:** Calculate sequence length, GC content, N count, IUPAC ambiguity count, entropy, and composition.
- **Method:** Remove whitespace, uppercase the sequence, validate it against the IUPAC DNA alphabet, and exclude ambiguity symbols from GC counts.
- **API:** `dnakit.sequence.calculate_sequence_metrics(sequence[required])`.
- **Input:** One DNA string, including optional whitespace and IUPAC ambiguity symbols.
- **Output:** A serializable mapping with `length`, `gc_fraction`, `gc_percent`, `n_count`, `ambiguous_count`, `entropy_bits`, and `composition`.
- **Example code:**

```python
from dnakit.sequence import calculate_sequence_metrics

print(calculate_sequence_metrics("ACGTN", verbose=False))
```

- **Example output:**

```json
{
  "length": 5,
  "gc_fraction": 0.4,
  "gc_percent": 40.0,
  "n_count": 1,
  "ambiguous_count": 1,
  "entropy_bits": 2.321928094887362,
  "composition": {"A": 1, "C": 1, "G": 1, "N": 1, "T": 1}
}
```

## 2) Count overlapping k-mers

- **Purpose:** Count every overlapping DNA k-mer.
- **Method:** Enumerate left-to-right windows and retain each k-mer count without canonicalization, correction, or deduplication.
- **API:** `dnakit.sequence.count_kmers(sequence[required], k[optional])`.
- **Input:** A DNA string and a positive integer `k`.
- **Output:** A dictionary sorted by k-mer text.
- **Example code:**

```python
from dnakit.sequence import count_kmers

print(count_kmers("ACACA", k=2, verbose=False))
```

- **Example output:**

```json
{"AC": 2, "CA": 2}
```

## 3) Deduplicate sequences

- **Purpose:** Remove duplicate DNA sequences while preserving first occurrence order.
- **Method:** Normalize whitespace and case; with `canonical=True`, use the lexicographically smaller of each sequence and its reverse complement as the deduplication key.
- **API:** `dnakit.sequence.deduplicate_sequences(sequences[required], canonical[optional])`.
- **Input:** A sequence of DNA strings and an optional Boolean `canonical`.
- **Output:** A list of normalized unique sequences.
- **Example code:**

```python
from dnakit.sequence import deduplicate_sequences

print(deduplicate_sequences(["ACG", " ACG ", "CGT"], canonical=True, verbose=False))
```

- **Example output:**

```json
["ACG"]
```

## 4) Reverse-complement DNA

- **Purpose:** Return the reverse complement of DNA, including IUPAC ambiguity symbols.
- **Method:** Validate the IUPAC DNA alphabet, reverse the symbols, and apply the DNA complement table.
- **API:** `dnakit.sequence.reverse_complement_dna(sequence[required])`.
- **Input:** A DNA string with optional IUPAC ambiguity symbols.
- **Output:** An uppercase DNA string.
- **Example code:**

```python
from dnakit.sequence import reverse_complement_dna

print(reverse_complement_dna("ACGTRN", verbose=False))
```

- **Example output:**

```text
NYACGT
```

## 5) Summarize sequence lengths and GC content

- **Purpose:** Summarize count, length distribution, and mean GC fraction for a collection of DNA/RNA strings in one pass.
- **Method:** Ignore whitespace, uppercase the input, convert U to T, and count ambiguity symbols toward length but not GC. The iterable is consumed once.
- **API:** `dnakit.sequence.sequence_length_gc_summary(sequences[required])`.
- **Input:** An iterable of DNA/RNA strings.
- **Output:** A serializable mapping with sequence count, per-sequence lengths, minimum/maximum/mean/median length, and mean GC fraction.
- **Example code:**

```python
from dnakit.sequence import sequence_length_gc_summary

print(sequence_length_gc_summary(["ACGT", "GGNN"], verbose=False))
```

- **Example output:**

```json
{
  "count": 2,
  "lengths": [4, 4],
  "min_length": 4,
  "max_length": 4,
  "mean_length": 4,
  "median_length": 4.0,
  "mean_gc_fraction": 0.5
}
```
