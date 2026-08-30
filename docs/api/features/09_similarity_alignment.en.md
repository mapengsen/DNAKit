# Similarity calculation

Calculate k-mer, fingerprint, and sketch similarities between DNA sequences and generate batch similarity matrices and reference library similarity results.

## 1) SIM-010 · k-mer similarity

**Function:** Convert two sequences into k-mer sets or count vectors, calculate Jaccard, Containment or Cosine similarity, which is used to quantify the degree of overlap of local components.

**API:** `dnakit.similarity.kmer_similarity(left[required], right[required], k[required], metric[optional], mode[optional], canonical[optional], overlapping[optional], zero_vector_policy[optional])`, `dnakit.similarity.kmer_vector_similarity(left[required], right[required], metric[optional], mode[optional], zero_vector_policy[optional])`.

**Input:** Two sequences and `k` are required; metric, set/count, canonical, overlap and zero vector strategies are optional.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import kmer_similarity

result = kmer_similarity(DNASequence("AAA"), DNASequence("AAC"), k=1)
print(result.value, result.components["shared_weight"])
```

**Example results:**

```text
0.5 1.0
```

## 2) SIM-011 · Fingerprint similarity

**Function:** Compare two DNA fingerprints or numerical vectors with consistent schema, calculate Tanimoto, Jaccard, Cosine or numerical distance, and return scores that can be used for sorting, clustering and threshold judgment.

**API:** `dnakit.similarity.fingerprint_similarity(left[required], right[required], metric[optional], weights[optional], zero_vector_policy[optional])`.

**Input:** Two isomorphic fingerprints/vectors are required; metric, feature weight and zero vector strategy are optional.

**Sample code:**

```python
from dnakit.similarity import fingerprint_similarity

result = fingerprint_similarity((1.0, 1.0), (1.0, 0.0), metric="tanimoto")
print(result.value)
```

**Example results:**

```text
0.5
```

## 3) SIM-012 · Sketch similarity

**Function:** Compare the hash values retained in two MinHash/FracMinHash Sketch, approximately estimate Jaccard or Containment, and report compatibility parameters, suitable for fast comparison of large sequences.

**API:** `dnakit.fingerprints.minhash(value[required], k[required], num_hashes[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, `dnakit.fingerprints.fracminhash(value[required], k[required], scaled[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, and `dnakit.similarity.sketch_similarity(left[required], right[required], metric[optional], min_shared_hashes[optional])`.

**Input:** Two compatible sketches are required; optional metric and minimum number of shared hashes.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash
from dnakit.similarity import sketch_similarity

left = minhash(DNASequence("ACGTAC"), k=2, num_hashes=100)
right = minhash(DNASequence("ACGTTC"), k=2, num_hashes=100)
result = sketch_similarity(left, right, min_shared_hashes=1)
print(result.value, result.shared_hash_count)
```

**Example results:**

```text
0.4 2
```

## 4) SIM-013 · Dashing sketch similarity

**Function:** Call external Dashing to build high-performance Sketch for multiple sequences and calculate approximate Jaccard matrices or Top-k neighbors, which is suitable for data scale that is unbearable for native exhaustive calculations.

**API:** `dnakit.similarity.DashingAdapter(executable_path[required])`, `dnakit.similarity.DashingAdapter.matrix(inputs[required], k[optional], mode[optional], sketch_size_log2[optional], canonical[optional], threads[optional], temp_dir[optional], output_path[optional], overwrite[optional], timeout_seconds[optional], max_items[optional], max_input_bytes[optional], max_output_bytes[optional], max_capture_bytes[optional], max_sketch_memory_bytes[optional])`, `dnakit.similarity.DashingAdapter.top_k(inputs[required], top_k[required], k[optional], mode[optional], sketch_size_log2[optional], canonical[optional], threads[optional], temp_dir[optional], output_path[optional], overwrite[optional], timeout_seconds[optional], max_items[optional], max_input_bytes[optional], max_output_bytes[optional], max_capture_bytes[optional], max_sketch_memory_bytes[optional])`.

**Input:** Required executable and sequence/FASTA/FASTQ explicitly configured by the caller; optional k, mode, sketch size, canonical, threads, temporary directory and output path.

**Sample code:**

```python
import os
from pathlib import Path

from dnakit import DNASequence
from dnakit.similarity import DashingAdapter

configured = os.environ.get("DNAKIT_DASHING_EXECUTABLE")
if configured is None:
    print("Skipped: configure DNAKIT_DASHING_EXECUTABLE explicitly first")
else:
    result = DashingAdapter(Path(configured)).matrix(
        (DNASequence("AACCGG"), DNASequence("AACCTT")), k=2, mode="exact"
    )
    print(result.values)
```

**Example results:**

```text
Skipped: configure DNAKIT_DASHING_EXECUTABLE explicitly first
```

## 6) `EVAL-011` Reference similarity

- **Function:** Find the most similar record in the versioned reference library for each query sequence, and return the reference ID, similarity, coverage and ranking for reference attribution and novelty analysis.
- **API**: `dnakit.evaluation.evaluate_reference_similarity(queries[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.ReferenceSearchConfig`.
- **Input**: query and `ReferenceLibrary`; optional similarity and coverage configuration.
- **Example code**: First construct `queries` and `reference` using the [reference library example](12_evaluation.md#eval-008-novelty), and then run the following code.

```python
from dnakit.evaluation import evaluate_reference_similarity

report = evaluate_reference_similarity(queries, reference)
print(report.metrics["mean_nearest_similarity"])
```

- **Example results:**

```text
0.5
```
