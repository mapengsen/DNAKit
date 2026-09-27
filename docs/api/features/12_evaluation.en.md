# Commonly used evaluation indicators

Comprehensive evaluation of DNA sequences or data sets in terms of legality, uniqueness, diversity, novelty, Fréchet representation distribution distance, fragment distribution, nearest neighbor similarity, ambiguity and redundancy.

Reference similarity and distribution similarity are listed separately in [similarity calculation](09_similarity_alignment.md) under data evaluation. Reference-library similarity and novelty use a caller-supplied `ReferenceLibrary`; reference-genome consistency uses a local FASTA or PAF. Novelty is not an experimental conclusion or task-model prediction.

The following API brackets list all direct call parameters; `[required]` and `[optional]` are parameter descriptions, not Python syntax.

For DNA sequence evaluation, use `DNA(["ACGT", "ACGA", ...])` as the standard batch input. To retain custom IDs or metadata, pass record mappings to `DNA`. A bare `list[str]` is not a shared evaluation input type. Per-sequence results appear in `report.entries`, and aggregate metrics in `report.metrics`. Collection-level measures such as uniqueness, diversity, fragment similarity, and distribution distance are not averages of per-sequence scores. Reference comparisons still require a reference library or a second `DNA` collection; numeric matrices, FASTQ paths, and scorecards keep their specialized input types.

| Input shape | Evaluation functions | Result shape |
| --- | --- | --- |
| One `DNA` collection | `evaluate_validity`, `evaluate_ambiguity`, `evaluate_complexity`, `evaluate_quality`, `evaluate_synthesis_risk` | Per-sequence `entries` and summary `metrics` |
| One `DNA` collection | `evaluate_uniqueness`, `evaluate_diversity`, `evaluate_kmer_diversity`, `evaluate_redundancy` | Collection-level `metrics` |
| `DNA` + `ReferenceLibrary` | `nearest_reference`, `evaluate_novelty`, `evaluate_memorization`, `evaluate_reference_similarity` | Per-query nearest-reference results and summary metrics |
| Two `DNA` collections | `evaluate_fragment_similarity`, `evaluate_snn`, `evaluate_distribution_similarity`, `evaluate_frechet_distance` | SNN has per-query entries; the others are mainly collection-level |
| `DNA` + reference FASTA/PAF | `evaluate_reference_consistency` | Per-query rates and overall mean |

`evaluate_scorecard` takes a mapping from metric names to values or reports, rather than DNA sequences.
For other operations that work on a single sequence, pass `DNA([...])` to [`run_batch`](14_engineering.en.md#3-eng-006-batch-calculation) and apply the operation per record; those operations do not implicitly return a collection mean.

## 1) Validity {#1-eval-001-validity}

Validity has two submodules: syntax validity and reference-genome consistency. They return separate reports; passing syntax checks does not establish consistency with a target reference genome.

### 1.1 Syntax validity {#1-1-syntax-validity}

- **Purpose:** Check DNA sequences and records according to rules such as alphabet, length, gap, metadata, etc., and return valid proportions and item-by-item questions, which are used to confirm whether the input is legal before evaluation or modeling.
- **Method:** Apply validation rules per record and divide the valid record count by the total count.
- **API:**  `dnakit.evaluation.evaluate_validity(value[required], sequence_length[optional], config[optional], limits[optional])`; `config` uses `dnakit.ValidationConfig`, `limits` uses `dnakit.evaluation.EvaluationLimits`.
- **Length requirement:** `sequence_length` is a non-negative integer requiring each valid sequence to contain exactly that many nucleotide symbols, excluding `Gap` spans. For nonempty input, `valid_fraction` is `1.0` only when every sequence meets the length requirement and all other validation rules; a mismatch produces `STD_SEQUENCE_LENGTH_MISMATCH`. The default `None` preserves the config; an explicit value overrides `config.sequence_length`, preserving other rules without modifying the original config.
- **Input:**  `DNA`, `DNASequence`, `DNARecord`, `DNASet`, or a record iterator; optionally `ValidationConfig` and resource cap.
- **Output:** Returns `EvaluationReport`; `metrics["valid_fraction"]` is the valid fraction, and `entries` holds per-record validity and issues.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_validity

sequences = DNA(["ACGT", "ACGA", "ACG"])
report = evaluate_validity(sequences, sequence_length=4)
print(report.metrics["valid_fraction"])  # 2 / 3
print([(entry.subject_id, entry.metrics["is_valid"]) for entry in report.entries])
```

- **Example output:**

```text
0.6666666666666666
[('sequence_1', True), ('sequence_2', True), ('sequence_3', False)]
```

### 1.2 Reference-genome consistency {#1-2-reference-consistency}

- **Purpose:** Measure how closely a group of DNA sequences matches a specified reference genome. Return a consistency rate for each query and the mean over all input queries.
- **Method:** Choose the primary alignment with the highest matched-base fraction for each query. Per-query consistency rate = `PAF matching bases / query length`; a query without a primary alignment contributes `0`. Overall `consistency_rate` is the mean across all input queries, or `None` for empty input. Unaligned bases remain in the query-length denominator.
- **API:** `dnakit.evaluation.evaluate_reference_consistency(queries[required], paf_path[optional], *, reference_fasta[optional], reference_dir[optional], reference_name[optional], reference_version[optional], reference_sha256[optional], config[optional], verbose[optional])`; `config` uses `dnakit.evaluation.ReferenceConsistencyConfig`.
- **Input:** `DNA([...])` or a DNA record collection. By default, DNAKit downloads and reuses the hg38 reference FASTA. You can supply a local `reference_fasta`; both paths require minimap2. For an existing `paf_path`, also provide `reference_name` and `reference_version`, and ensure PAF query IDs and lengths match the input.
- **Output:** Returns an `EvaluationReport`. `report.metrics["consistency_rate"]` is the mean over all input queries; each item's `metrics["consistency_rate"]` in `report.entries` is that query's rate.
- **Example code:** This illustrative PAF shows how the mean is calculated.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNA
from dnakit.evaluation import evaluate_reference_consistency

sequences = DNA(["ACGTACGT", "TGCATGCA", "CCCCCCCC"])
with TemporaryDirectory() as folder:
    paf = Path(folder) / "example.paf"
    paf.write_text(
        "sequence_1\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P\tcg:Z:8M\n"
        "sequence_2\t8\t0\t8\t-\tchr2\t100\t30\t38\t7\t8\t30\ttp:A:P\tcg:Z:8M\n",
        encoding="utf-8",
    )
    report = evaluate_reference_consistency(
        sequences, paf, reference_name="example", reference_version="demo-1",
        verbose=False,
    )

print(report.metrics["consistency_rate"])
print([(entry.subject_id, entry.metrics["consistency_rate"]) for entry in report.entries])
```

- **Example output:**

```text
0.625
[('sequence_1', 1.0), ('sequence_2', 0.875), ('sequence_3', 0.0)]
```

## 2) Uniqueness {#2-eval-005-uniqueness}

- **Purpose:** Group sequences according to exact, reverse complement, IUPAC or approximate equivalent rules, return the number of unique sequences, proportions and repeat groups, used to measure the proportion of independent records in the set.
- **Method:** Deduplicate using the equivalence rule and divide the unique count by the total count.
- **API:**  `dnakit.evaluation.evaluate_uniqueness(value[required], config[optional], plot[optional])`; `config` uses `dnakit.evaluation.UniquenessEvaluationConfig`.
- **Input:**  `DNA([...])` or another DNA record collection; optional equivalence rules, methods and thresholds.
- **Output:** Returns `EvaluationReport`; `metrics` contains the uniqueness score and duplicate groups.
- **Optional plot:** `plot=True` saves `uniqueness_ratio.png` in the current working directory (default `False`; an existing file with the same name is overwritten). The default `exact` equivalence uses the title “Exact uniqueness”; other equivalences are named in the title. Plotting requires `dnakit[bio-plots]`.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_uniqueness

records = DNA(["AAAA", "AAAA", "CCCC"])
report = evaluate_uniqueness(records)
print(report.metrics["uniqueness_score"])  # 0.666666...
print(report.metrics["duplicate_groups"])  # (("sequence_1", "sequence_2"),)
```

- **Example output:**

```text
0.6666666666666666
(('sequence_1', 'sequence_2'),)
```

## 3) Diversity {#3-eval-006-diversity}

- **Purpose:** Calculate the mean raw Levenshtein distance over all sequence pairs using the paper's formula.
- **Method:** `Σ(i≠j) Levenshtein(xᵢ,xⱼ) / [n(n−1)]`, equivalent to the mean raw edit distance over all unordered pairs; returns `None` for fewer than two sequences.
- **API:**  `dnakit.evaluation.evaluate_diversity(value[required], config[optional], plot[optional])`; `config` uses `dnakit.evaluation.DiversityEvaluationConfig`.
- **Input:**  `DNA([...])` or another DNA record collection; set `show_progress=True` to display pairwise progress. The default uses up to 16 processes. Set `DiversityEvaluationConfig(workers=8)` to adjust this, or `workers=1` for serial execution; the metric is unchanged.
- **Output:** Returns `EvaluationReport`; `metrics` contains mean pairwise Levenshtein distance.
- **Optional plots:** `plot=True` saves `diversity_mean_bar.png` and `diversity_pairwise_distance_histogram.png` in the current working directory (default `False`; existing files with the same names are overwritten). The histogram includes every unordered sequence pair; plotting requires `dnakit[bio-plots]`.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_diversity

records = DNA(["AAAA", "AAAT", "CCCC"])
report = evaluate_diversity(records)
print(report.metrics["mean_pairwise_levenshtein_distance"])
```

Multicore execution is the default. To change the process count, import `DiversityEvaluationConfig` from `dnakit.evaluation` and call `evaluate_diversity(records, config=DiversityEvaluationConfig(workers=8))`. This uses processes because the pure-Python edit-distance loop is limited by the GIL.

- **Example output:**

```text
3
```

For canonical k-mer Jaccard distance, use `evaluate_kmer_diversity(records, k=6)`. Its `mean_pair_distance` is the mean `1 − Jaccard` over unordered pairs; `mean_nearest_neighbor_distance` averages each sequence's shortest distance to another member. Both are `None` for fewer than two records. Reverse-complement canonicalization is enabled by default, and `EvaluationLimits` bounds pairwise work. This is a different metric from raw Levenshtein distance.

```python
from dnakit.evaluation import evaluate_kmer_diversity

report = evaluate_kmer_diversity(records, k=6)
print(report.metrics["mean_nearest_neighbor_distance"])
```

[FAQ: details](../../faq/tutorial_basics.md#levenshtein-3)

## 4) Novelty {#4-eval-008-novelty}

- **Purpose:** Report continuous novelty and full-record exact-match rate side by side. By default, each query has `novelty = 1 - nearest-reference similarity`, and `score` is their mean. `exact_match_rate` is the fraction of queries literally equal to any full reference record. A second method reports each query's minimum raw Levenshtein distance and its mean, alongside the exact-match rate.
- **Method:** Check literal equality against every reference record independently of nearest-hit ranking. Novelty no longer uses a `0.9` cutoff to produce `is_novel` or `novel_fraction`. `copy_threshold` applies only to memorization/near-copy evaluation and has no effect on these novelty results.
- **API:**  `dnakit.evaluation.evaluate_novelty(queries[required], reference[required], config[optional], plot[optional])`; `config` uses `dnakit.evaluation.ReferenceSearchConfig`.
- **Input:**  a `DNA([...])` query collection and a non-empty `ReferenceLibrary`; the default is `novelty_calculation="similarity"`, while the paper method uses `novelty_calculation="levenshtein"`; `show_progress=True` is optional.
- **Output:** Returns `EvaluationReport`; `metrics` contains mean novelty score and exact-match rate.
- **Optional plots:** `plot=True` saves `novelty_batch_summary.png` (mean novelty and exact-match rate), `novelty_distribution.png` (per-query novelty), and `novelty_nearest_levenshtein.png` (raw nearest-reference edit distances) in the current working directory. The default is `False`; existing files with these names are overwritten. Nonempty text in query-record `metadata["batch"]` groups the summary bars by batch; otherwise they show one evaluated set. The third plot performs an additional Levenshtein scan in similarity mode. Plotting requires `dnakit[bio-plots]`. In Levenshtein mode, per-query novelty already is the nearest raw edit distance, so the two histograms show the same distances.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import ReferenceSearchConfig, create_reference_library, evaluate_novelty

reference_records = DNA(["AAAA", "CCCC"])
reference = create_reference_library(
    reference_records,
    name="training",
    version="1",
    source="local:example",
)
queries = DNA(["AAAA", "GGGG"])

report = evaluate_novelty(
    queries,
    reference,
    config=ReferenceSearchConfig(method="identity"),
)
print(report.metrics["score"])             # mean novelty: 0.5
print(report.metrics["exact_match_rate"])  # full-record exact-match rate: 0.5
```

- **Example output:**

```text
0.5
0.5
```

[FAQ: details](../../faq/tutorial_basics.md#levenshtein-4)

## 5) Fréchet DNA distance {#5-eval-016-frechet-dna-distance}

- **Purpose:** Use the same basic DNA model to represent two sequence sets respectively, approximate the two sets of vectors into multivariate Gaussian distributions, and then calculate the Fréchet distance of the mean and covariance. The smaller the value, the closer the distribution of the two sets in the representation space is, and the exact same representation distribution approaches `0`.
- **Method:** Combine squared mean differences with the covariance Fréchet term using an equivalent sample-space cross-Gram nuclear norm.
- **API:** `dnakit.evaluation.evaluate_frechet_distance(left[required], right[required], config[optional], backend[optional], plot_distance[optional], plot_pca[optional], plot_output_dir[optional], plot_overwrite[optional])`; both plot flags default to `False`. `config` uses `dnakit.evaluation.FrechetDistanceConfig`, whose representation configuration uses `dnakit.representations.RepresentationConfig`.
- **Input:** Two `DNA([...])` collections with at least two sequences each, or other record collections; this metric defaults to `lucaone`, mean pooling, and L2 normalization. Other registered models, checkpoints, pooling methods, devices, dtypes, and batch sizes remain selectable.
- **Output:** Returns `EvaluationReport`; `metrics["frechet_distance"]` is the Fréchet distance. Enabled plots are saved as `frechet_distance.png` and/or `frechet_pca.png` under `plot_output_dir` (the current directory by default). Absolute paths are in `parameters["plot_files"]`. Existing files are preserved unless `plot_overwrite=True`.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import FrechetDistanceConfig, evaluate_frechet_distance
from dnakit.representations import RepresentationConfig

reference = DNA(["ACGTACGT", "AACCGGTT"])
generated = DNA(["ACGTTCGT", "AACCAGTT"])
report = evaluate_frechet_distance(
    generated,
    reference,
    config=FrechetDistanceConfig(
        representation=RepresentationConfig(allow_remote_code=True),
    ),
    plot_distance=True,
    plot_pca=True,
    plot_output_dir="frechet_plots",
)
print(report.metrics["frechet_distance"])
print(report.parameters["plot_files"])
```

- **Example output:**

```text
Non-negative floating-point number; the exact value depends on the model, checkpoint, pooling, normalization, and input sets. The example also saves two PNG figures.
```

- **Based on**: [Preuer et al.’s original paper on FCD](https://doi.org/10.1021/acs.jcim.8b00234).


[FAQ: details](../../faq/tutorial_basics.md#evaluation-background-5)

## 6) Frag {#6-eval-017-frag}

- **Purpose:** Refer to the Frag in the molecule generation evaluation, convert the generated set and the reference set into fragment occurrence number vectors, and then calculate the cosine similarity. The value range is `[0,1]`, and the higher it is, the closer the local fragment distributions of the two sets are.
- **Method:** Count overlapping k-mers in each collection and calculate cosine similarity between count vectors.
- **API:**  `dnakit.evaluation.evaluate_fragment_similarity(generated[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.FragmentSimilarityConfig`.
- **Input:** Generated and reference collections can both be `DNA([...])`; default k=3 and canonical k-mers, skipping ambiguous windows.
- **Output:** Returns `EvaluationReport`; `metrics["frag"]` is the fragment-distribution score.
- **DNA Adapter:** Molecular fragments use BRICS fragments; DNA has no corresponding chemical bond breaking rules, so this implementation explicitly uses overlapping fixed-length k-mers. Default `k=3`, canonical k-mer, ignore windows containing IUPAC ambiguous bases.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import FragmentSimilarityConfig, evaluate_fragment_similarity

generated = DNA(["AAAA", "CCCC"])
reference = DNA(["AAAA", "GGGG"])
report = evaluate_fragment_similarity(
    generated,
    reference,
    config=FragmentSimilarityConfig(k=2, canonical=False, show_progress=False),
)
print(report.metrics["frag"])  # 0.5
```

- **Example output:**

```text
0.5
```

- **Based on**: [MOSES's definition of Frag](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full).


[FAQ: details](../../faq/tutorial_basics.md#evaluation-background-6)

## 7) SNN {#7-eval-018-snn}

- **Purpose:** For each generated DNA, find the sequence with the highest binary fingerprint Tanimoto similarity in the reference set, and then take the arithmetic average of these nearest neighbor similarities. The value range is `[0,1]`, and the higher it is, the closer the generated sample is to the representation space covered by the reference set.
- **Method:** Take the highest reference Tanimoto fingerprint similarity per query and average across queries; the metric is generally asymmetric.
- **API:**  `dnakit.evaluation.evaluate_snn(generated[required], reference[required], config[optional])`; `config` uses `dnakit.evaluation.SNNConfig`.
- **Input:** Generated and reference collections can both be `DNA([...])`; default canonical 7-mer hashes in a 1024-bit fingerprint.
- **Output:** Returns `EvaluationReport`; `metrics["snn"]` is the nearest-neighbor similarity measure.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import SNNConfig, evaluate_snn

generated = DNA(["AAAAAAA", "ATATATA"])
reference = DNA(["AAAAAAA", "CCCCCCC"])
report = evaluate_snn(
    generated,
    reference,
    config=SNNConfig(show_progress=False),
)
print(report.metrics["snn"])  # 0.5
print(report.entries[0].metrics["nearest_reference_id"])  # sequence_1
```

- **Example output:**

```text
0.5
sequence_1
```

- **Based on**: [MOSES's definition of SNN](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full).


[FAQ: details](../../faq/tutorial_basics.md#evaluation-background-7)

## 8) Ambiguity {#8-eval-002-ambiguity}

- **Purpose:** Count the number, position and proportion of IUPAC ambiguous bases in each sequence, and determine whether it exceeds the limit according to the configuration, which is used to quantify the uncertainty base burden of the data.
- **Method:** Count ambiguous IUPAC symbols, calculate their fraction with the configured gap denominator policy and check thresholds.
- **API:**  `dnakit.evaluation.evaluate_ambiguity(value[required], config[optional])`; `config` uses `dnakit.evaluation.AmbiguityEvaluationConfig`.
- **Input:**  One or more DNA sequences; for a batch, use `DNA([...], alphabet="iupac")`. Threshold, symbol weights, and Gap denominator policy are configurable.
- **Output:** Returns `EvaluationReport`; per-record `entries` contain ambiguous-base counts and fractions, while `metrics` holds the summary.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_ambiguity

report = evaluate_ambiguity(DNA(["ACNT", "ACGT"], alphabet="iupac"))
entry = report.entries[0]
print(entry.metrics["ambiguity_count"])     # 1
print(entry.metrics["ambiguity_fraction"])  # 0.25
print(report.metrics["ambiguity_fraction_summary"]["mean"])  # 0.125
```

- **Example output:**

```text
1
0.25
0.125
```

## 9) Redundancy {#9-eval-007-redundancy}

- **Purpose:** Count the proportion of completely repeated and approximately repeated sequences in the sequence set, and report the approximate sequence pairs and cluster compression ratio, which are used to determine whether there is excessive repetition in the data.
- **Method:** Summarize exact duplicate groups and threshold-based approximate groups into redundancy and compression fractions.
- **API:**  `dnakit.evaluation.evaluate_redundancy(value[required], config[optional])`; `config` uses `dnakit.evaluation.DiversityEvaluationConfig`.
- **Input:**  `DNA([...])` or another DNA record collection; optional similarity method and threshold.
- **Output:** Returns `EvaluationReport`; `metrics` contains redundancy score and exact-duplicate fraction.
- **Example code:**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_redundancy

records = DNA(["AAAA", "AAAA", "AAAT"])
report = evaluate_redundancy(records)
print(report.metrics["score"])
print(report.metrics["exact_duplicate_fraction"])
```

- **Example output:**

```text
0.3333333333333333
0.3333333333333333
```

<span id="eval-008-novelty"></span>

## 10) Generate FastQC quality reports {#bio-19}

- **Purpose:** Inspect read-quality distributions and sequencing QC reports.
- **Method:** Run FastQC on the supplied FASTQ files.
- **API:** `dnakit.bioinformatics.ngs.quality_control(reads, output_directory, config=None)`
- **Input:** A list of FASTQ paths, FastQC on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import quality_control

result = quality_control(["R1.fq.gz", "R2.fq.gz"], "new_qc")
print(sorted(result.artifacts))
```

- **Example output:**

Returns `ToolResult` with HTML/ZIP reports for the input files in `artifacts`. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

## 11) Calculate properties for a sequence batch {#bio-26}

- **Purpose:** Build a property matrix from custom prediction functions.
- **Method:** Apply named scalar predictors to each sequence while preserving ID and column order.
- **API:** `dnakit.bioinformatics.metrics.sequence_properties(sequences, predictors)`
- **Input:** `DNA([...])` or an ID-to-sequence mapping, plus named Python scalar functions; requires `bio-analysis`. The result preserves `DNA` record IDs and order.
- **Output:** Returns a mapping with sequence IDs, property names, and a per-sequence property matrix.
- **Example code:**

```python
from dnakit import DNA
from dnakit.bioinformatics.metrics import sequence_properties

result = sequence_properties(DNA(["ACGT", "GGGG"]),
                             {"gc": lambda s: sum(x in "GC" for x in s) / len(s)})
print(result["ids"], result["names"], result["values"].tolist())
```

- **Example output:**

```text
('sequence_1', 'sequence_2') ('gc',) [[0.5], [1.0]]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 12) Calculate MMD² or KID distribution distance {#bio-27}

- **Purpose:** Compare generated and reference embedding distributions.
- **Method:** MMD uses a Gaussian kernel; KID uses a cubic polynomial kernel, with unbiased estimation by default.
- **API:** `dnakit.bioinformatics.metrics.kernel_distance(generated, reference, metric='mmd', sigma=1.0, unbiased=True)`
- **Input:** Two embedding matrices with matching feature counts; unbiased estimates can be negative; requires `bio-analysis`.
- **Output:** Returns the MMD² or KID value as a float.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import kernel_distance

values = [[0, 0], [1, 1], [2, 0]]
print(kernel_distance(values, values, unbiased=False))
```

- **Example output:**

```text
0.0
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 13) Calculate PRDC precision, recall, density and coverage {#bio-28}

- **Purpose:** Measure local overlap between generated and reference distributions.
- **Method:** Use k-neighbor radii and strict open balls to measure coverage between the sets.
- **API:** `dnakit.bioinformatics.metrics.manifold_metrics(generated, reference, neighbors=5)`
- **Input:** Embeddings with matching dimensions; each set must contain more than `neighbors` samples. Fix the embedding model and scaling; requires `bio-analysis`.
- **Output:** Returns a mapping with `precision`, `recall`, `density`, and `coverage`.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import manifold_metrics

values = [[0], [1], [3], [6]]
print(manifold_metrics(values, values, neighbors=1))
```

- **Example output:**

```text
{'precision': 1.0, 'recall': 1.0, 'density': 1.0, 'coverage': 1.0}
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 14) Calculate clipped density and coverage {#bio-29}

- **Purpose:** Measure coverage with clipping and finite-sample calibration.
- **Method:** Use median-clipped radii, leave-one-out density and beta-binomial coverage calibration.
- **API:** `dnakit.bioinformatics.metrics.clipped_manifold_metrics(generated, reference, neighbors=5)`
- **Input:** Two embedding sets with equal size and dimension, each larger than `neighbors`; requires `bio-analysis`.
- **Output:** Returns a mapping with `clipped_density` and `clipped_coverage`.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import clipped_manifold_metrics

values = [[0], [1], [3], [6]]
print(clipped_manifold_metrics(values, values, neighbors=1))
```

- **Example output:**

```text
{'clipped_density': 1.0, 'clipped_coverage': 1.0}
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 15) Measure authenticity against the training set {#bio-30}

- **Purpose:** Estimate the fraction of generated samples that are not excessively close to training samples.
- **Method:** Compare each generated point's nearest-training distance with that training point's nearest-other-training distance.
- **API:** `dnakit.bioinformatics.metrics.authenticity(generated, training)`
- **Input:** Generated embeddings and the generator's actual training embeddings; at least two training samples; requires `bio-analysis`.
- **Output:** Returns an authenticity score as a `float`.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import authenticity

print(authenticity([[0], [20]], [[0], [1], [2]]))
```

- **Example output:**

```text
0.5
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 16) Calculate Pareto ranks and hypervolume {#bio-31}

- **Purpose:** Compare candidates optimized for multiple objectives.
- **Method:** Assign nondominated fronts by per-column direction and compute exact hypervolume against a worse reference point.
- **API:** `dnakit.bioinformatics.metrics.pareto_ranks(values, maximize=...)`; `dnakit.bioinformatics.metrics.hypervolume(values, reference_point, maximize=...)`
- **Input:** A samples×objectives matrix and per-column maximize/minimize flags; hypervolume supports up to 5 dimensions and 1000 points; requires `bio-analysis`.
- **Output:** `pareto_ranks` returns sample ranks; `hypervolume` returns a hypervolume value.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import pareto_ranks, hypervolume

values = [[1, 1], [2, 2], [3, 1]]
print(pareto_ranks(values, maximize=[True, True]).tolist())
print(hypervolume(values, [0, 0], maximize=[True, True]))
```

- **Example output:**

```text
[1, 0, 0]
5.0
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 17) Calculate property hit rates {#bio-32}

- **Purpose:** Measure how many samples satisfy all requested property ranges.
- **Method:** Test closed intervals per column and return joint/per-property rates and a mask.
- **API:** `dnakit.bioinformatics.metrics.property_hit_rate(values, lower=..., upper=...)`
- **Input:** A property matrix and per-column bounds; requires `bio-analysis`.
- **Output:** Returns a mapping with the joint hit rate and hit rates for each property.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import property_hit_rate

result = property_hit_rate([[0.5, 10], [0.7, 15]], lower=[0.4, 5], upper=[0.6, 20])
print(result["joint"])
print(sorted(result))
```

- **Example output:**

```text
0.5
['joint', 'mask', 'per_property']
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 18) Calculate KL divergence of property distributions {#bio-33}

- **Purpose:** Compare generated and reference distributions for each property.
- **Method:** Use shared histogram edges and a pseudocount to compute natural-log KL(generated∥reference).
- **API:** `dnakit.bioinformatics.metrics.property_kl(generated, reference, bins=20, pseudocount=1e-8)`
- **Input:** Two property matrices with matching dimensions; fix bins and pseudocount before comparison; requires `bio-analysis`.
- **Output:** Returns KL divergence values for each property.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import property_kl

values = [[0], [1], [2], [3]]
print(property_kl(values, values, bins=2).tolist())
```

- **Example output:**

```text
[0.0]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).

## 19) Calculate a distribution conformity score {#bio-34}

- **Purpose:** Score generated samples using the reference density distribution.
- **Method:** Use fixed-bandwidth Gaussian KDE, cross-fit reference data, then score against validation-density empirical distributions.
- **API:** `dnakit.bioinformatics.metrics.conformity_score(generated, reference, folds=5, bandwidth=1.0, seed=0)`
- **Input:** Generated and independent reference embeddings, with enough samples for the folds; requires `bio-analysis`. The score is not a class probability.
- **Output:** Returns a mapping with the consistency score, standard error, and cross-validation scores.
- **Example code:**

```python
from dnakit.bioinformatics.metrics import conformity_score

result = conformity_score([[0], [1], [2]], [[0], [1], [2], [3], [4], [5]], folds=3)
print(sorted(result))
print(round(result["value"], 6))
```

- **Example output:**

```text
['fold_scores', 'standard_error', 'value']
0.62963
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#metrics).
