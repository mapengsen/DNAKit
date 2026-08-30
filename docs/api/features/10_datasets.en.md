# Data partitioning

Partition DNA datasets by random, stratified, similarity, or metadata constraints to create isolated data subsets for training, validation, and testing.

All stochastic or heuristic processes log seeds, configurations, groupings, and resource caps.

## 1) DATA-012 · Random and stable hash partitioning

**Function:** Randomly allocate records to train, valid, and test according to the target proportion, and combine stable ID and seed to generate reproducible results, and return each subset and allocation list.


**API:** `dnakit.datasets.split(records[required], config[required])`; `config` uses `dnakit.datasets.SplitConfig`, this item is set to `method="random"` or `method="hash"`.

**Input:** Required records and proportion that sum to 1; optional seed, shuffle and whether to maintain the original order within the subset. The `hash` pattern requires a stable and unique `record.id` for each record.

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
