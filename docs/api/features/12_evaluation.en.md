# Commonly used evaluation indicators

Comprehensive evaluation of DNA sequences or data sets in terms of legality, uniqueness, diversity, novelty, Fréchet representation distribution distance, fragment distribution, nearest neighbor similarity, ambiguity and redundancy.

`EVAL-011` and `EVAL-012` are listed separately in [similarity calculation](09_similarity_alignment.md) under data evaluation. All "reference-dependent" results must be bound to the local `ReferenceLibrary` provided by the caller; novelty is not experimental conclusions or mission model predictions.

The following API brackets list all direct call parameters; `[required]` and `[optional]` are parameter descriptions, not Python syntax.

## 1) `EVAL-001` Validity

- **Function:** Check DNA sequences and records according to rules such as alphabet, length, gap, metadata, etc., and return valid proportions and item-by-item questions, which are used to confirm whether the input is legal before evaluation or modeling.
- **API**: `dnakit.evaluation.evaluate_validity(value[required], config[optional], limits[optional])`; `config` uses `dnakit.ValidationConfig`, `limits` uses `dnakit.evaluation.EvaluationLimits`.
- **Input**: `DNASequence`, `DNARecord`, `DNASet` or corresponding iterator; optionally `ValidationConfig` and resource cap.
- **Sample Code**:

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_validity

report = evaluate_validity(DNASequence("ACGT"))
print(report.metrics["valid_fraction"])  # 1.0
```

- **Example results:**

```text
1.0
```

## 2) `EVAL-005` Uniqueness

- **Function:** Group sequences according to exact, reverse complement, IUPAC or approximate equivalent rules, return the number of unique sequences, proportions and repeat groups, used to measure the proportion of independent records in the set.
- **API**: `dnakit.evaluation.evaluate_uniqueness(value[required], config[optional])`; `config` uses `dnakit.evaluation.UniquenessEvaluationConfig`.
- **Input**: `DNASet`; Optional equivalence rules, methods and thresholds.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_uniqueness

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAA"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
])
report = evaluate_uniqueness(records)
print(report.metrics["uniqueness_score"])  # 0.666666...
print(report.metrics["duplicate_groups"])  # (("a", "b"),)
```

- **Example results:**

```text
0.6666666666666666
(('a', 'b'),)
```

## 3) `EVAL-006` Diversity

- **Function:** By default, calculate normalized similarity distances, nearest-neighbor distances, and threshold clusters. A second method follows the paper's mean raw pairwise Levenshtein-distance formula.
- **API**: `dnakit.evaluation.evaluate_diversity(value[required], config[optional])`; `config` uses `dnakit.evaluation.DiversityEvaluationConfig`.
- **Input**: `DNASet`; the default is `calculation="similarity"`, while the paper method uses `calculation="levenshtein"`; set `show_progress=True` to display pairwise progress.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_diversity

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
])
report = evaluate_diversity(records)
print(report.metrics["mean_pair_distance"])
print(report.metrics["cluster_count"])
```

- **Example results:**

```text
1.0
3
```

- **Second method (mean pairwise Levenshtein distance):**

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import DiversityEvaluationConfig, evaluate_diversity

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
])
report = evaluate_diversity(
    records,
    config=DiversityEvaluationConfig(calculation="levenshtein"),
)
print(report.metrics["mean_pairwise_levenshtein_distance"])
```

- **Example results:**

```text
3.0
```

The value is an unnormalized count of edit operations; larger values mean greater within-set sequence differences. The formula is undefined for fewer than two records and returns `None`. Avoid directly comparing datasets with different sequence-length distributions.

<a id="eval-008-novelty"></a>

## 4) `EVAL-008` Novelty

**Reference library example preparation:**

This project uses a versioned local reference library. The minimal construction method is as follows:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import create_reference_library

reference_records = DNASet([
    DNARecord(DNASequence("AAAA"), "ref-a"),
    DNARecord(DNASequence("CCCC"), "ref-c"),
])
reference = create_reference_library(
    reference_records,
    name="training",
    version="1",
    source="local:example",
)
queries = DNASet([
    DNARecord(DNASequence("AAAA"), "copy"),
    DNARecord(DNASequence("GGGG"), "query-new"),
])
```

- **Function:** By default, evaluate novelty as `1 - nearest-reference similarity`. A second method computes each query's minimum raw Levenshtein distance to the reference library and averages those distances.
- **API**: `dnakit.evaluation.evaluate_novelty(queries[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.ReferenceSearchConfig`.
- **Input**: query collection and a non-empty `ReferenceLibrary`; the default is `novelty_calculation="similarity"`, while the paper method uses `novelty_calculation="levenshtein"`; `show_progress=True` is optional.
- **Example Code**: Run after "Reference Library Example Preparation".

```python
from dnakit.evaluation import ReferenceSearchConfig, evaluate_novelty

report = evaluate_novelty(
    queries,
    reference,
    config=ReferenceSearchConfig(method="identity", copy_threshold=0.9),
)
print(report.metrics["novel_fraction"])  # 0.5
```

- **Example results:**

```text
0.5
```

- **Second method (mean nearest-reference Levenshtein distance):**

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import (
    ReferenceSearchConfig,
    create_reference_library,
    evaluate_novelty,
)

reference = create_reference_library(
    DNASet([
        DNARecord(DNASequence("AAAA"), "ref-a"),
        DNARecord(DNASequence("CCCC"), "ref-c"),
    ]),
    name="training",
    version="1",
    source="local:example",
)
queries = DNASet([
    DNARecord(DNASequence("AAAA"), "copy"),
    DNARecord(DNASequence("GGGG"), "query-new"),
])
report = evaluate_novelty(
    queries,
    reference,
    config=ReferenceSearchConfig(novelty_calculation="levenshtein"),
)
print(report.metrics["mean_nearest_levenshtein_distance"])
```

- **Example results:**

```text
2.0
```

The value is an unnormalized count of edit operations; larger values mean the queries are farther from the reference set. This mode does not use similarity, coverage, or `copy_threshold`; per-query `is_novel` only means that the nearest edit distance is greater than `0`. DNAKit does not pad or truncate sequences automatically.

## 5) `EVAL-016` Fréchet DNA distance

- **Function:** Use the same basic DNA model to represent two sequence sets respectively, approximate the two sets of vectors into multivariate Gaussian distributions, and then calculate the Fréchet distance of the mean and covariance. The smaller the value, the closer the distribution of the two sets in the representation space is, and the exact same representation distribution approaches `0`.
- **API**: `dnakit.evaluation.evaluate_frechet_distance(left[required], right[required], config[optional], backend[optional])`; `config` uses `dnakit.evaluation.FrechetDistanceConfig`, where the representation configuration uses `dnakit.representations.RepresentationConfig`.
- **Input**: Two sets each containing at least 2 `DNARecord`; default reuse `DATA-027`'s `lucaone`, mean pooling and L2 normalization. Additional registered models, checkpoints, pooling, devices, dtypes, and batch sizes can be selected.
- **Formula**: `||μ_left - μ_right||² + Tr(Σ_left + Σ_right - 2(Σ_left^(1/2) Σ_right Σ_left^(1/2))^(1/2))`. The implementation uses the mathematically equivalent sample space cross-Gram kernel norm without explicitly creating a high-dimensional dense covariance matrix.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import FrechetDistanceConfig, evaluate_frechet_distance
from dnakit.representations import RepresentationConfig

reference = DNASet([
    DNARecord(DNASequence("ACGTACGT"), "ref-1"),
    DNARecord(DNASequence("AACCGGTT"), "ref-2"),
])
generated = DNASet([
    DNARecord(DNASequence("ACGTTCGT"), "gen-1"),
    DNARecord(DNASequence("AACCAGTT"), "gen-2"),
])
report = evaluate_frechet_distance(
    generated,
    reference,
    config=FrechetDistanceConfig(
        representation=RepresentationConfig(allow_remote_code=True),
    ),
)
print(report.metrics["frechet_distance"])
```

- **Example results:**

```text
Non-negative floating-point number; the exact value depends on the model, checkpoint, pooling, normalization, and input sets.
```

- **Progress:** Checkpoint download and sequence-by-sequence characterization extraction show progress bars by default; can be turned off with `RepresentationConfig(show_progress=False)`.
- **Based on**: [Preuer et al.’s original paper on FCD](https://doi.org/10.1021/acs.jcim.8b00234); [LucaOne’s original paper](https://doi.org/10.1038/s42256-025-01044-4).

## 6) `EVAL-017` Frag

- **Function:** Refer to the Frag in the molecule generation evaluation, convert the generated set and the reference set into fragment occurrence number vectors, and then calculate the cosine similarity. The value range is `[0,1]`, and the higher it is, the closer the local fragment distributions of the two sets are.
- **API**: `dnakit.evaluation.evaluate_fragment_similarity(generated[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.FragmentSimilarityConfig`.
- **DNA Adapter:** Molecular fragments use BRICS fragments; DNA has no corresponding chemical bond breaking rules, so this implementation explicitly uses overlapping fixed-length k-mers. Default `k=3`, canonical k-mer, ignore windows containing IUPAC ambiguous bases.
- **Formula**: `Σ_f c_generated(f)c_reference(f) / sqrt(Σ_f c_generated(f)² × Σ_f c_reference(f)²)`.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import FragmentSimilarityConfig, evaluate_fragment_similarity

generated = DNASet([
    DNARecord(DNASequence("AAAA"), "gen-a"),
    DNARecord(DNASequence("CCCC"), "gen-c"),
])
reference = DNASet([
    DNARecord(DNASequence("AAAA"), "ref-a"),
    DNARecord(DNASequence("GGGG"), "ref-g"),
])
report = evaluate_fragment_similarity(
    generated,
    reference,
    config=FragmentSimilarityConfig(k=2, canonical=False, show_progress=False),
)
print(report.metrics["frag"])  # 0.5
```

- **Example results:**

```text
0.5
```

- **Progress:** Display k-mer statistics progress for the generated and reference sets by default; can be turned off with `FragmentSimilarityConfig(show_progress=False)`.
- **Based on**: [MOSES's definition of Frag](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full).

## 7) `EVAL-018` SNN

- **Function:** For each generated DNA, find the sequence with the highest binary fingerprint Tanimoto similarity in the reference set, and then take the arithmetic average of these nearest neighbor similarities. The value range is `[0,1]`, and the higher it is, the closer the generated sample is to the representation space covered by the reference set.
- **API**: `dnakit.evaluation.evaluate_snn(generated[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.SNNConfig`.
- **DNA adaptation:** Molecular SNN uses a 1024-bit, radius 2 Morgan fingerprint; this implementation defaults to a 1024-bit DNA binary fingerprint obtained by SHA-256 mapping of the canonical 7-mer, and calculates the Tanimoto similarity equivalent to the binary Jaccard.
- **Formula**: `mean_g max_r Tanimoto(fp(g), fp(r))`; This indicator uses the generated set as a query, so it is usually asymmetric.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import SNNConfig, evaluate_snn

generated = DNASet([
    DNARecord(DNASequence("AAAAAAA"), "copy"),
    DNARecord(DNASequence("ATATATA"), "far"),
])
reference = DNASet([
    DNARecord(DNASequence("AAAAAAA"), "ref-a"),
    DNARecord(DNASequence("CCCCCCC"), "ref-c"),
])
report = evaluate_snn(
    generated,
    reference,
    config=SNNConfig(show_progress=False),
)
print(report.metrics["snn"])  # 0.5
print(report.entries[0].metrics["nearest_reference_id"])  # ref-a
```

- **Example results:**

```text
0.5
ref-a
```

- **Progress:** Displays two sets of fingerprint construction and nearest neighbor scan progress by default; can be turned off with `SNNConfig(show_progress=False)`.
- **Based on**: [MOSES's definition of SNN](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full).

## 8) `EVAL-002` Ambiguity

- **Function:** Count the number, position and proportion of IUPAC ambiguous bases in each sequence, and determine whether it exceeds the limit according to the configuration, which is used to quantify the uncertainty base burden of the data.
- **API**: `dnakit.evaluation.evaluate_ambiguity(value[required], config[optional])`; `config` uses `dnakit.evaluation.AmbiguityEvaluationConfig`.
- **Input**: A line or set of DNA; configurable maximum scale, sign weight, and Gap denominator strategies.
- **Sample Code**:

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_ambiguity

report = evaluate_ambiguity(DNASequence("ACNT", alphabet="iupac"))
entry = report.entries[0]
print(entry.metrics["ambiguity_count"])     # 1
print(entry.metrics["ambiguity_fraction"])  # 0.25
```

- **Example results:**

```text
1
0.25
```

## 9) `EVAL-007` Redundancy

- **Function:** Count the proportion of completely repeated and approximately repeated sequences in the sequence set, and report the approximate sequence pairs and cluster compression ratio, which are used to determine whether there is excessive repetition in the data.
- **API**: `dnakit.evaluation.evaluate_redundancy(value[required], config[optional])`; `config` uses `dnakit.evaluation.DiversityEvaluationConfig`.
- **Input**: `DNASet`; Optional similarity method and threshold.
- **Sample Code**:

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_redundancy

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAA"), "b"),
    DNARecord(DNASequence("AAAT"), "c"),
])
report = evaluate_redundancy(records)
print(report.metrics["score"])
print(report.metrics["exact_duplicate_fraction"])
```

- **Example results:**

```text
0.3333333333333333
0.3333333333333333
```
