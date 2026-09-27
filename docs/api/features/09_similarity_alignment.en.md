# Similarity calculation

Calculate k-mer, fingerprint, and sketch similarities between DNA sequences and generate batch similarity matrices and reference library similarity results.

## 1) k-mer similarity {#1-sim-010-k-mer-similarity}

- **Purpose:** Convert sequences into k-mer sets or count vectors to calculate Jaccard, Containment, or Cosine similarity; collections can also be compared by their mean similarity.
- **Method:** For one sequence on each side, compute set Jaccard/containment or cosine similarity of k-mer count vectors. For collections, compare each query with every reference; `value` is the arithmetic mean of all pair scores: `sum(sim(query, reference)) / (query count × reference count)`. This is not the mean nearest-neighbor score.

- **API:** `dnakit.similarity.kmer_similarity(left[required], right[required], k[required], metric[optional], mode[optional], canonical[optional], overlapping[optional], zero_vector_policy[optional], max_pairwise_comparisons[optional])`, `dnakit.similarity.kmer_vector_similarity(left[required], right[required], metric[optional], mode[optional], zero_vector_policy[optional])`.

- **Input:** `left` holds queries and `right` holds references. Either side accepts a single `DNASequence`/`DNARecord`, `DNA([...])`, `DNASet`, or a list of `DNASequence`/`DNARecord` objects. Collections must be nonempty. `k` is required; metric, set/count, canonical, overlap, and zero-vector policy are optional. At most 1,000,000 pairs are compared by default; adjust `max_pairwise_comparisons` if needed.
- **Output:** Returns `SimilarityResult`. A single pair retains fields such as `components["shared_weight"]`; for collections, `value` is the overall mean and `components["query_means"]` gives each query's mean against all references in input order, together with `query_count`, `reference_count`, and `pair_count`.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import kmer_similarity

result = kmer_similarity(DNASequence("AAA"), DNASequence("AAC"), k=1)
print(result.value, result.components["shared_weight"])
```

- **Example output:**

```text
0.5 1.0
```

- **Collection example:**

```python
from dnakit import DNA
from dnakit.similarity import kmer_similarity

queries = DNA(["AAAA", "GGGG"])
references = DNA(["AAAA", "AAAT"])
result = kmer_similarity(queries, references, k=2)
print(result.value, result.components["query_means"])
```

```text
0.375 (0.75, 0.0)
```

## 2) Fingerprint similarity {#2-sim-011-fingerprint-similarity}

- **Purpose:** Compare two DNA fingerprints or numerical vectors with consistent schema, calculate Tanimoto, Jaccard, Cosine or numerical distance, and return scores that can be used for sorting, clustering and threshold judgment.
- **Method:** Validate matching schemas, then compare bit sets or numeric vectors with the chosen metric.

- **API:** `dnakit.similarity.fingerprint_similarity(left[required], right[required], metric[optional], weights[optional], zero_vector_policy[optional])`.

- **Input:** Two isomorphic fingerprints/vectors are required; metric, feature weight and zero vector strategy are optional.
- **Output:** Returns `SimilarityResult` with similarity `value` and comparison `method`.

- **Example code:**

```python
from dnakit.similarity import fingerprint_similarity

result = fingerprint_similarity((1.0, 1.0), (1.0, 0.0), metric="tanimoto")
print(result.value)
```

- **Example output:**

```text
0.5
```

## 3) Sketch similarity {#3-sim-012-sketch-similarity}

- **Purpose:** Compare the hash values retained in two MinHash/FracMinHash Sketch, approximately estimate Jaccard or Containment, and report compatibility parameters, suitable for fast comparison of large sequences.
- **Method:** Check k, seed and sampling compatibility, then estimate overlap from shared hashes.

- **API:** `dnakit.fingerprints.minhash(value[required], k[required], num_hashes[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, `dnakit.fingerprints.fracminhash(value[required], k[required], scaled[optional], canonical[optional], seed[optional], max_hashes[optional], max_unique_hashes[optional])`, and `dnakit.similarity.sketch_similarity(left[required], right[required], metric[optional], min_shared_hashes[optional])`.

- **Input:** Two compatible sketches are required; optional metric and minimum number of shared hashes.
- **Output:** `minhash`/`fracminhash` return `SketchResult`; `sketch_similarity` returns `SketchSimilarityResult` with similarity `value` and shared-hash count.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash
from dnakit.similarity import sketch_similarity

left = minhash(DNASequence("ACGTAC"), k=2, num_hashes=100)
right = minhash(DNASequence("ACGTTC"), k=2, num_hashes=100)
result = sketch_similarity(left, right, min_shared_hashes=1)
print(result.value, result.shared_hash_count)
```

- **Example output:**

```text
0.4 2
```

## 4) Dashing sketch similarity {#4-sim-013-dashing-sketch-similarity}

- **Purpose:** Call external Dashing to build high-performance Sketch for multiple sequences and calculate approximate Jaccard matrices or Top-k neighbors, which is suitable for data scale that is unbearable for native exhaustive calculations.
- **Method:** Run the explicitly configured Dashing executable and parse the sketch matrix or Top-k results.

- **API:** `dnakit.similarity.DashingAdapter(executable_path[required])`, `dnakit.similarity.DashingAdapter.matrix(inputs[required], k[optional], mode[optional], sketch_size_log2[optional], canonical[optional], threads[optional], temp_dir[optional], output_path[optional], overwrite[optional], timeout_seconds[optional], max_items[optional], max_input_bytes[optional], max_output_bytes[optional], max_capture_bytes[optional], max_sketch_memory_bytes[optional])`, `dnakit.similarity.DashingAdapter.top_k(inputs[required], top_k[required], k[optional], mode[optional], sketch_size_log2[optional], canonical[optional], threads[optional], temp_dir[optional], output_path[optional], overwrite[optional], timeout_seconds[optional], max_items[optional], max_input_bytes[optional], max_output_bytes[optional], max_capture_bytes[optional], max_sketch_memory_bytes[optional])`.

- **Input:** Required executable and sequence/FASTA/FASTQ explicitly configured by the caller; optional k, mode, sketch size, canonical, threads, temporary directory and output path.
- **Output:** Returns similarity results from the Dashing sketch comparison.

- **Example code:**

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

- **Example output:**

```text
Skipped: configure DNAKIT_DASHING_EXECUTABLE explicitly first
```

## 5) Reference similarity {#6-eval-011-reference-similarity}

- **Purpose:** Find the most similar record in the versioned reference library for each query sequence, and return the reference ID, similarity, coverage and ranking for reference attribution and novelty analysis.
- **Method:** Search the reference library with the chosen metric, then filter and rank by score and coverage.
- **API:**  `dnakit.evaluation.evaluate_reference_similarity(queries[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.ReferenceSearchConfig`.
- **Input:**  query and `ReferenceLibrary`; optional similarity and coverage configuration.
- **Output:** Returns `EvaluationReport` with aggregate `metrics` and per-record `entries`.
- **Example code:**  First construct `queries` and `reference` using the [reference library example](12_evaluation.md#eval-008-novelty), and then run the following code.

```python
from dnakit.evaluation import evaluate_reference_similarity

report = evaluate_reference_similarity(queries, reference)
print(report.metrics["mean_nearest_similarity"])
```

- **Example output:**

```text
0.5
```
