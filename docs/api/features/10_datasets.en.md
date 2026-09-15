# Data partitioning

Public processing APIs share [file and text input rules](00_processing_inputs.md).

Partition DNA datasets by random, stratified, similarity, or metadata constraints to create isolated data subsets for training, validation, and testing.

All stochastic or heuristic processes log seeds, configurations, groupings, and resource caps.

## 1) DATA-012 · Random and stable hash partitioning

**Function:** Partition records into named subsets such as `train / valid / test` by ratios, or generate a fixed K-fold partition and select one fold as `valid` and the rest as `train`. Return the subsets and assignment manifest.


**API:** `dnakit.datasets.split(records[required], config[required])`; use `dnakit.datasets.SplitConfig` with `method="random"`, `method="hash"`, or `method="kfold"`.

**Input:** Records are required. Ratio splitting uses fractions summing to 1, defaulting to `train=0.8, valid=0.1, test=0.1`. K-fold uses `n_splits` for the fold count and `fold_id` for the validation fold. Defaults are `seed=42` and `shuffle=True`. Hash splitting requires stable, unique record IDs.

The default seed for `SplitConfig` and `dnakit split` changed from `0` to `42`. Pass `seed=0` explicitly to reproduce earlier ratio splits that omitted the seed.

**Sample code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

**Example results:**

```text
{'train': 8, 'test': 2}
```

<span id="1-hash"></span>**Order-independent `hash` partitioning**

The `hash` pattern uses versioned SHA-256 to calculate the stable sort key for `seed + record.id` when `shuffle=True` and the stable endianness for `record.id` when `shuffle=False` and allocates the records proportionally. Therefore, even if the input order of the same batch of records changes, the classification will still be exactly the same when viewed by `record.id`; it does not use Python's built-in process randomization `hash()`.

**Sample code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
config = SplitConfig(
    method="hash",
    ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
    seed=17,
    preserve_order=False,
)
first = split(records, config=config)
second = split(list(reversed(records)), config=config)
first_by_id = {item.record_id: item.split for item in first.assignments}
second_by_id = {item.record_id: item.split for item in second.assignments}
print(first_by_id == second_by_id)
print(first.get("train").ids == second.get("train").ids)
```

**Example results:**

```text
True
True
```

### Fixed K-fold partitioning: train / valid

Set `method="kfold"` to use [scikit-learn KFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html) semantics:

```python
KFold(n_splits=10, shuffle=True, random_state=42)
```

**One random seed → one fixed K-fold partition → K validation folds selected in turn.** Keep the records and their input order, `seed`, `shuffle`, and `n_splits` unchanged; change only `fold_id`. For 10 folds, use `fold_id=0..9`: each sample serves as validation exactly once. Do not use `seed + fold_id`. Each call reconstructs the same partition and returns the selected pair.

| `SplitConfig` argument | Default | Meaning |
| --- | --- | --- |
| `method` | `"random"` | Select `"kfold"` for K-fold; `"random"` / `"hash"` retain ratio splitting |
| `n_splits` | `10` in K-fold mode | Fold count, such as 5 or 10; requires `2 <= n_splits <= sample count` |
| `fold_id` | `0` in K-fold mode | Validation fold index, from `0` through `n_splits - 1` |
| `seed` | `42` | Fixed seed for the whole partition; K-fold requires `0 <= seed < 2**32` |
| `shuffle` | `True` | Shuffle before partitioning; `False` gives consecutive validation blocks and seed has no effect |
| `ratios` | `0.8 / 0.1 / 0.1` in ratio mode | Omit for K-fold; target train/valid ratios are derived from the fold count |

K-fold always returns subsets named `train` and `valid`. Records within each subset retain input order, matching scikit-learn, even with `preserve_order=False`. If the sample count is not divisible by K, the first `sample count % K` validation folds contain one extra sample; fold sizes differ by at most 1. Results and exported manifests record `seed`, `shuffle`, `n_splits`, `fold_id`, and each record's assignment.

This mode requires scikit-learn: install `python -m pip install -e ".[bio-analysis]"`. Environments with `.[ml]` or `.[training]` already include it. Ratio splitting does not require this dependency.

**Example: hold out a fixed test set, then apply 10-fold splitting to train.**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("ACGT"), f"r{i}") for i in range(100)]
seed = 42
n_splits = 10  # Change to 5 for five folds.

# Split off the fixed test set once.
outer = split(
    records,
    config=SplitConfig(method="random", ratios={"train": 0.8, "test": 0.2}, seed=seed),
)
training_pool = outer.get("train")
test = outer.get("test")

# Only fold_id changes; never pass the test set into K-fold.
for fold_id in range(n_splits):
    fold = split(
        training_pool,
        config=SplitConfig(
            method="kfold", n_splits=n_splits, fold_id=fold_id,
            shuffle=True, seed=seed,
        ),
    )
    train = fold.get("train")
    valid = fold.get("valid")
    print(f"fold={fold_id}: train={len(train)}, valid={len(valid)}, test={len(test)}")
    # Create and train a fresh model here, then evaluate on valid.
```

Folds 0 through 9 each contain `train=72, valid=8, test=20`; validation samples rotate while the 20 test samples stay fixed. If you already have training data, use it directly as `training_pool`. `n_splits` sets the total fold count, and `fold_id` selects one fold. Iterate all folds for a complete K-fold procedure; `split()` partitions data and does not train models or aggregate metrics.

To use `train_custom_property_model()`, write the current train/valid subsets and fixed test set using the training API's supported split labels, and map that column with `DataMapping.split_field`. Values such as `fold_0` are not training split labels. See [custom model training](23_deep_learning_property_prediction.en.md).

**CLI:** Pass the existing training pool to export the selected train/valid pair and manifest:

```bash
dnakit split train.fasta runs/fold_0 --method kfold --n-splits 10 --fold-id 0 --seed 42
```

Keep the input and seed fixed for subsequent runs, changing `--fold-id` and the output directory. Use `--no-shuffle` to disable shuffling. Agent/MCP calls to `dnakit.datasets.split` accept the same `config`. JSON/YAML workflow `split.params` also accepts `method`, `n_splits`, and `fold_id`; its seed comes from the workflow's global `seed`, which should be set to `42` for this example.

## 2) DATA-013 · Stratified random division

**Function:** Try to maintain the original proportion of the specified category label in train, valid, and test, return the stratified subset and label statistics, and reduce the category distribution shift.


**API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="stratified"`.

**Input:** Required record, scale, and `metadata_key`; optional seed, missing field policy, and order policy.

**Sample code:**

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

**Example results:**

```text
{'train': 4, 'test': 4}
```

## 3) DATA-014 · Similarity division

**Function:** First establish a sequence group based on the similarity threshold, and then assign train, valid, and test in units of the entire group to prevent data leakage of approximate sequences across sets.


**API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="similarity"`.

**Input:** Required record and scale; optional k, threshold, IUPAC/Gap policy, seed, and size cap.

**Sample code:**

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

**Example results:**

```text
{'a': 'train', 'b': 'train', 'c': 'test', 'd': 'test'}
```

## 4) DATA-015 · Cluster split

**Function:** Use the existing cluster label as an inseparable group to divide the data, ensuring that all records of the same cluster only enter a subset, and reporting the actual proportion deviation.


**API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="group", metadata_key="cluster"`.

**Input:** Required records, cluster metadata and proportion for each record; optional seed, missing value policy.

**Sample code:**

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

**Example results:**

```text
{'a1': 'train', 'a2': 'train', 'b1': 'test'}
```

## 5) DATA-016 · Species classification

**Function:** Group records according to species metadata and then divide them as a whole to ensure that the same species does not appear in the training set and evaluation set at the same time, which is used to test cross-species generalization.


**API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.

**Input:** Required species metadata key and proportion; optional seed, missing value strategy.

**Sample code:**

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

**Example results:**

```text
{'h1': 'train', 'h2': 'train', 'm1': 'test'}
```

## 6) DATA-017 · Chromosome Division

**Function:** Group records according to chromosome metadata and then divide them as a whole to avoid leakage caused by highly correlated regions of the same chromosome across train, valid, and test.


**API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.

**Input:** Required chromosome metadata key and proportion; optional seed, missing value strategy.

**Sample code:**

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

**Example results:**

```text
{'c1-a': 'train', 'c1-b': 'train', 'c2-a': 'test'}
```

## 7) DATA-018 · Individual division

**Function:** Group samples by individuals or donor metadata and divide them as a whole to ensure that the same individual does not span data subsets to avoid the leakage of individual-specific information.


**API:** `dnakit.datasets.split(records[required], config[required])`; uses `group` mode for `dnakit.datasets.SplitConfig`.

**Input:** Required individual metadata key and proportion; optional seed, missing value strategy.

**Sample code:**

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

**Example results:**

```text
{'p1-a': 'train', 'p1-b': 'train', 'p2-a': 'test'}
```

## 8) Divide by custom label

**Function:** Group records according to any metadata label (such as family, locus, batch) specified by the caller, ensuring that the same label value only enters a subset, and returns the group allocation result.

**API:** `dnakit.datasets.split(records[required], config[required])`; use `dnakit.datasets.SplitConfig(method="group", metadata_key=label)`.

**Input:** Required records, division ratios, and customization `label`. `label` is the metadata field name in each record, such as `donor`, `family`, `locus`, `batch` or other caller-provided fields.

**Sample code:**

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

**Example results:**

```text
{'f1-a': 'train', 'f1-b': 'train', 'f2-a': 'test'}
```

## 9) DATA-023 · Leak detection

**Function:** Perform cross-set exact or approximate comparisons on train, valid, test and other sets, and return leaked sequence pairs, similarities, involved subsets and summary proportions, which are used to verify partition independence.


**API:** `dnakit.datasets.detect_leakage(splits[required], config[optional])`; `config` uses `dnakit.datasets.LeakageConfig`.

**Input:** Required Split name to `DNASet` mapping; optional identity/edit/k-mer/fingerprint, threshold and resource cap.

**Sample code:**

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

**Example results:**

```text
True 1
```
