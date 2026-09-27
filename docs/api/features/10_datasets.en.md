# Data partitioning

Partition DNA datasets by random, stratified, similarity, or metadata constraints to create isolated data subsets for training, validation, and testing.

All stochastic or heuristic processes log seeds, configurations, groupings, and resource caps.

These record-collection APIs accept `DNA(["ACGT", "TTAA", ...])` directly: `DNA` yields its `DNARecord` values in order. For stable IDs or per-record metadata, use `DNA([{"id": "r1", "sequence": "ACGT", "metadata": {...}}, ...])`, or continue to pass a `DNARecord` iterator.

## 1) Random and stable hash partitioning {#1-data-012-random-and-stable-hash-partitioning}

- **Purpose:** Partition records into named subsets such as `train / valid / test` by ratios, or generate a fixed K-fold partition and select one fold as `valid` and the rest as `train`. Return the subsets and assignment manifest.
- **Method:** Use seeded random ordering, stable ID hashing, or fixed K-fold partitions with a selected validation fold.
- **API:** `dnakit.datasets.split(records[required], config[required])`; use `dnakit.datasets.SplitConfig` with `method="random"`, `method="hash"`, or `method="kfold"`.
- **Input:** Records are required. Ratio splitting uses fractions summing to 1, defaulting to `train=0.8, valid=0.1, test=0.1`. K-fold requires scikit-learn and uses `n_splits` for the fold count and `fold_id` for the validation fold. Defaults are `seed=42` and `shuffle=True`. Hash splitting requires stable, unique record IDs.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.

The default seed for `SplitConfig` and `dnakit split` changed from `0` to `42`. Pass `seed=0` explicitly to reproduce earlier ratio splits that omitted the seed.

- **Example code:**

```python
from dnakit import DNA
from dnakit.datasets import SplitConfig, split

records = DNA([{"id": f"r{i}", "sequence": "AC"} for i in range(10)])
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

- **Example output:**

```text
{'train': 8, 'test': 2}
```

<span id="1-hash"></span>Hash and K-fold examples, installation and seed compatibility: [FAQ: details](../../faq/tutorial_integrations.md#dataset-splits).

For paired strands sharing a `pair_id`, store it in each record's metadata and use `SplitConfig(method="group", metadata_key="pair_id")` for a group-aware ratio split. For row files, `export_grouped_split` assigns each `pair_id` once and writes the assignment manifest and split CSVs without building a `DNARecord` for every row. `method="kfold"` assigns input records individually; to keep paired chunks in the same fold, make one representative record per `pair_id` and map its fold back to both strands. Fig.4 task A05 uses the file API for hg38 chunks, while `--n-splits` controls the two benchmarks' official training folds.

### File-based grouped splits and K-fold exports

`dnakit.datasets.export_grouped_split(source, output_dir, group_field=..., config=...)` performs a random ratio split by group field. The default `SplitConfig()` means an 8:1:1 split with seed 42. It writes `assignment.csv` and one CSV per split, keeping every row in a group together and preserving input row order. Optional `sample_sizes` selects reservoir samples per split, `sample_filter` chooses eligible rows, and `sample_digest_field` adds a `<field>_sha256` column. Samples are written under `analysis/`. The return value contains `seed`, `group_counts`, and `analysis_counts`. This file API currently accepts `SplitConfig(method="random")`. Input may be any supported row file, such as CSV or Parquet; output is CSV.

```python
from dnakit.datasets import SplitConfig, export_grouped_split, export_kfold_table

summary = export_grouped_split(
    "chunks.csv", "split_output", group_field="pair_id",
    config=SplitConfig(seed=42),
    sample_sizes={"train": 800, "valid": 100, "test": 100},
    sample_filter=lambda row: row["strand"] == "+",
    sample_digest_field="sequence",
    overwrite=True,
)
fold_counts = export_kfold_table(
    "benchmark_train.csv", "benchmark_train_with_folds.csv",
    config=SplitConfig(method="kfold", n_splits=10, seed=42),
    label_field="label", overwrite=True,
)
```

`export_kfold_table` preserves input row order and appends zero-based `fold` values. With `label_field`, it returns per-fold training and validation label counts, which can be saved with `dnakit.io.save_json`. It defaults to 10 folds and seed 42 and requires scikit-learn. It exports every fold, so `SplitConfig.fold_id` does not select one fold here.

## 2) Stratified random division {#2-data-013-stratified-random-division}

- **Purpose:** Try to maintain the original proportion of the specified category label in train, valid, and test, return the stratified subset and label statistics, and reduce the category distribution shift.
- **Method:** Group records by label and allocate each group to the requested split ratios.
- **API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="stratified"`.
- **Input:** Required record, scale, and `metadata_key`; optional seed, missing field policy, and order policy.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), f"p{i}", metadata={"label": "p"})
    for i in range(4)
] + [
    DNARecord(DNASequence("GT"), f"n{i}", metadata={"label": "n"})
    for i in range(4)
]
result = split(
    records,
    config=SplitConfig(
        method="stratified", ratios={"train": 0.5, "test": 0.5},
        metadata_key="label", seed=2,
    ),
)
print(dict(result.counts))
```

- **Example output:**

```text
{'train': 4, 'test': 4}
```

## 3) Similarity division {#3-data-014-similarity-division}

- **Purpose:** First establish a sequence group based on the similarity threshold, and then assign train, valid, and test in units of the entire group to prevent data leakage of approximate sequences across sets.
- **Method:** Group by the k-mer similarity threshold and assign complete groups to splits.
- **API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="similarity"`.
- **Input:** Required record and scale; optional k, threshold, IUPAC/Gap policy, seed, and size cap.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
    DNARecord(DNASequence("CCCG"), "d"),
]
result = split(
    records,
    config=SplitConfig(
        method="similarity", ratios={"train": 0.5, "test": 0.5},
        similarity_k=2, similarity_threshold=0.5, seed=5,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'a': 'train', 'b': 'train', 'c': 'test', 'd': 'test'}
```

## 4) Cluster split {#4-data-015-cluster-split}

- **Purpose:** Use the existing cluster label as an inseparable group to divide the data, ensuring that all records of the same cluster only enter a subset, and reporting the actual proportion deviation.
- **Method:** Group records by cluster labels and keep each complete cluster in one split.
- **API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="group", metadata_key="cluster"`.
- **Input:** Required records, cluster metadata and proportion for each record; optional seed, missing value policy.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "a1", metadata={"cluster": "A"}),
    DNARecord(DNASequence("AG"), "a2", metadata={"cluster": "A"}),
    DNARecord(DNASequence("GT"), "b1", metadata={"cluster": "B"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group", ratios={"train": 0.5, "test": 0.5},
        metadata_key="cluster", seed=11,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'a1': 'train', 'a2': 'train', 'b1': 'test'}
```

## 5) Species classification {#5-data-016-species-classification}

- **Purpose:** Group records according to species metadata and then divide them as a whole to ensure that the same species does not appear in the training set and evaluation set at the same time, which is used to test cross-species generalization.
- **Method:** Group by species metadata and keep each species in one split.
- **API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.
- **Input:** Required species metadata key and proportion; optional seed, missing value strategy.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "h1", metadata={"species": "human"}),
    DNARecord(DNASequence("AG"), "h2", metadata={"species": "human"}),
    DNARecord(DNASequence("GT"), "m1", metadata={"species": "mouse"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="species"
))
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'h1': 'train', 'h2': 'train', 'm1': 'test'}
```

## 6) Chromosome Division {#6-data-017-chromosome-division}

- **Purpose:** Group records according to chromosome metadata and then divide them as a whole to avoid leakage caused by highly correlated regions of the same chromosome across train, valid, and test.
- **Method:** Group by chromosome metadata and keep each chromosome in one split.
- **API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.
- **Input:** Required chromosome metadata key and proportion; optional seed, missing value strategy.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "c1-a", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("AG"), "c1-b", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("GT"), "c2-a", metadata={"chromosome": "chr2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="chromosome"
))
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'c1-a': 'train', 'c1-b': 'train', 'c2-a': 'test'}
```

## 7) Individual division {#7-data-018-individual-division}

- **Purpose:** Group samples by individuals or donor metadata and divide them as a whole to ensure that the same individual does not span data subsets to avoid the leakage of individual-specific information.
- **Method:** Group by individual metadata and keep each individual in one split.
- **API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.
- **Input:** Required individual metadata key and proportion; optional seed, missing value strategy.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "p1-a", metadata={"individual": "p1"}),
    DNARecord(DNASequence("AG"), "p1-b", metadata={"individual": "p1"}),
    DNARecord(DNASequence("GT"), "p2-a", metadata={"individual": "p2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="individual"
))
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'p1-a': 'train', 'p1-b': 'train', 'p2-a': 'test'}
```

## 8) Divide by custom label

- **Purpose:** Group records according to any metadata label (such as family, locus, batch) specified by the caller, ensuring that the same label value only enters a subset, and returns the group allocation result.
- **Method:** Group by the chosen metadata label and assign complete groups to splits.
- **API:** `dnakit.datasets.split(records[required], config[required])`; use `dnakit.datasets.SplitConfig(method="group", metadata_key=label)`.
- **Input:** Required records, division ratios, and customization `label`. `label` is the metadata field name in each record, such as `donor`, `family`, `locus`, `batch` or other caller-provided fields.
- **Output:** Returns `SplitResult` with named `subsets`, their `counts`, and per-record `assignments`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

label = "family"  # Can also be donor, locus, or another metadata field name
records = [
    DNARecord(DNASequence("AC"), "f1-a", metadata={label: "f1"}),
    DNARecord(DNASequence("AG"), "f1-b", metadata={label: "f1"}),
    DNARecord(DNASequence("GT"), "f2-a", metadata={label: "f2"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group",
        ratios={"train": 0.5, "test": 0.5},
        metadata_key=label,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **Example output:**

```text
{'f1-a': 'train', 'f1-b': 'train', 'f2-a': 'test'}
```

## 9) Leak detection {#9-data-023-leak-detection}

- **Purpose:** Perform cross-set exact or approximate comparisons on train, valid, test and other sets, and return leaked sequence pairs, similarities, involved subsets and summary proportions, which are used to verify partition independence.
- **Method:** Compare sequences or grouping information across splits and report leakage matches.
- **API:** `dnakit.datasets.detect_leakage(splits[required], config[optional])`; `config` uses `dnakit.datasets.LeakageConfig`.
- **Input:** Required Split name to `DNASet` mapping; optional identity/edit/k-mer/fingerprint, threshold and resource cap.
- **Output:** Returns `LeakageReport` with a leakage flag, event counts, and cross-split hits.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.datasets import LeakageConfig, detect_leakage

splits = {
    "train": DNASet([DNARecord(DNASequence("AAAA"), "a")]),
    "test": DNASet([DNARecord(DNASequence("AAAA"), "same")]),
}
report = detect_leakage(
    splits, config=LeakageConfig(method="identity", threshold=0.9)
)
print(report.has_leakage, report.exact_event_count)
```

- **Example output:**

```text
True 1
```

<!-- BEGIN GENERATED API: encode-collection -->

<span id="encode-collection"></span>

## 10) Compare ENCODE experiments for joint analysis compatibility {#dnakit.encode.compare_experiments}

- **Purpose:** Compare ENCODE experiments for joint analysis compatibility.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.compare_experiments(accession1, accession2, *, verbose=True)`
- **Input:** Provide `accession1`, `accession2`. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON compatibility report with verdict, issues, warnings, and recommendations.
- **Example code:**

```python
import dnakit

result = dnakit.encode.compare_experiments(
    accession1='ENCSR133RZO',
    accession2='ENCSR000AKS',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON compatibility report with verdict, issues, warnings, and recommendations.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.compare_experiments-parameters)

## 11) Export tracked ENCODE experiments as a table {#dnakit.encode.export_data}

- **Purpose:** Export tracked ENCODE experiments as a table.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.export_data(format='csv', assay_title=None, organism=None, organ=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  Formatted table data in the requested format.
- **Example code:**

```python
import dnakit

result = dnakit.encode.export_data(
    format='csv',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  Formatted table data in the requested format.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.export_data-parameters)

## 12) Retrieve publications and citations for tracked experiments {#dnakit.encode.get_citations}

- **Purpose:** Retrieve publications and citations for tracked experiments.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.get_citations(accession=None, export_format='json', *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  Publications in the requested format.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_citations(
    export_format='bibtex',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  Publications in the requested format.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.get_citations-parameters)

## 13) Retrieve derived-file provenance {#dnakit.encode.get_provenance}

- **Purpose:** Retrieve derived-file provenance.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.get_provenance(file_path=None, source_accession=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON provenance chain or list of derived files.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_provenance(
    file_path='analysis/consensus_peaks.bed',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON provenance chain or list of derived files.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.get_provenance-parameters)

## 14) Retrieve external references linked to experiments {#dnakit.encode.get_references}

- **Purpose:** Retrieve external references linked to experiments.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.get_references(experiment_accession=None, reference_type=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON with linked external references.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_references(
    reference_type='doi',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with linked external references.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.get_references-parameters)

## 15) Link an external reference to a tracked experiment {#dnakit.encode.link_reference}

- **Purpose:** Link an external reference to a tracked experiment.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.link_reference(experiment_accession, reference_type, reference_id, description='', *, verbose=True)`
- **Input:** Provide `experiment_accession`, `reference_type`, `reference_id`. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON with linking result.
- **Example code:**

```python
import dnakit

result = dnakit.encode.link_reference(
    experiment_accession='ENCSR133RZO',
    reference_type='doi',
    reference_id='10.1038/nature11247',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with linking result.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.link_reference-parameters)

## 16) List locally tracked ENCODE experiments {#dnakit.encode.list_tracked}

- **Purpose:** List locally tracked ENCODE experiments.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.list_tracked(assay_title=None, organism=None, organ=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON with tracked experiments metadata table and tracker stats.
- **Example code:**

```python
import dnakit

result = dnakit.encode.list_tracked(
    assay_title='ATAC-seq',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with tracked experiments metadata table and tracker stats.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.list_tracked-parameters)

## 17) Record a derived file and its provenance {#dnakit.encode.log_derived_file}

- **Purpose:** Record a derived file and its provenance.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.log_derived_file(file_path, source_accessions, description='', file_type='', tool_used='', parameters='', *, verbose=True)`
- **Input:** Provide `file_path`, `source_accessions`. Retrieve a real file accession with the file query API and replace `ENCFF_REPLACE_ME`. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON with the provenance record ID.
- **Example code:**

```python
import dnakit

result = dnakit.encode.log_derived_file(
    file_path='analysis/consensus_peaks.bed',
    source_accessions=['ENCFF_REPLACE_ME'],
    description='Consensus peak set',
    file_type='bed',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with the provenance record ID.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.log_derived_file-parameters)

## 18) Summarize tracked experiment collections by field {#dnakit.encode.summarize_collection}

- **Purpose:** Summarize tracked experiment collections by field.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.summarize_collection(assay_title=None, organism=None, organ=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON summary with experiment counts grouped by multiple dimensions.
- **Example code:**

```python
import dnakit

result = dnakit.encode.summarize_collection(
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON summary with experiment counts grouped by multiple dimensions.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.summarize_collection-parameters)

## 19) Track an ENCODE experiment and its provenance locally {#dnakit.encode.track_experiment}

- **Purpose:** Track an ENCODE experiment and its provenance locally.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.track_experiment(accession, fetch_publications=True, fetch_pipelines=True, notes='', *, verbose=True)`
- **Input:** Provide `accession`. Collection operations use local tracking records. Track experiments before comparison, reference linking or export, and register derived files before provenance queries.
- **Output:** Python return type: str.  JSON with tracking result including publications and pipeline info found.
- **Example code:**

```python
import dnakit

result = dnakit.encode.track_experiment(
    accession='ENCSR133RZO',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with tracking result including publications and pipeline info found.

[Full parameters and usage notes](../../faq/api/10_datasets.md#dnakit.encode.track_experiment-parameters)

<!-- END GENERATED API: encode-collection -->
