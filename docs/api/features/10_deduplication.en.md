# Remove duplicates

Deduplication is performed on the DNA sequence to obtain the DNA sequence after deduplication.

## 1) Standard deduplication {#1-data-001-standard-deduplication}

- **Purpose:** Group the sequence by character-by-character equality, retain only one representative record in each group, and return the merged record and its representative mapping for lossless deletion of exact copies.
- **Method:** Group identical sequences, select representatives and resolve metadata conflicts by configuration.

- **API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; `config` uses `dnakit.datasets.DeduplicationConfig`.

- **Input:** Required `DNARecord` collection; optional represents policy, conflicting fields/policies, and metadata merging.
- **Output:** Returns `DeduplicationResult` with retained `records`, equivalence `groups`, and deduplication counts.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("ACGT"), "a"),
    DNARecord(DNASequence("TTAA"), "b"),
    DNARecord(DNASequence("ACGT"), "c"),
]
result = deduplicate(records, equivalence="exact")
print(result.records.ids, result.groups[0].member_ids)
```

- **Example output:**

```text
('a', 'b') ('a', 'c')
```

## 2) Reverse complementary deduplication {#2-data-002-reverse-complementary-deduplication}

- **Purpose:** Compare the forward and reverse complementary forms of sequences, and group records with only different strand directions but the same content into one group to avoid repeated counting caused by double-stranded directions.
- **Method:** Group sequences with their reverse complements and select one representative.

- **API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; this setting is `equivalence="reverse_complement"`.

- **Input:** Required `DNARecord` collection; optional represents selection and metadata conflict configuration.
- **Output:** Returns `DeduplicationResult` with retained `records`, equivalence `groups`, and deduplication counts.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("AAGC"), "forward"),
    DNARecord(DNASequence("GCTT"), "reverse"),
]
result = deduplicate(records, equivalence="reverse_complement")
print(result.groups[0].member_ids, result.groups[0].orientations)
```

- **Example output:**

```text
('forward', 'reverse') ('forward', 'reverse_complement')
```

## 3) Ring equivalent deduplication {#3-data-003-ring-equivalent-deduplication}

- **Purpose:** Perform rotational equivalence comparison on circular sequences, merge records with different storage starting points but the same base sequence on the loop, and retain representative records and membership relationships.
- **Method:** Normalize circular rotations to an equivalence key, optionally including reverse complements.

- **API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; use `equivalence="circular"` or `equivalence="circular_reverse_complement"` for this item.

- **Input:** Required Records explicitly marked `topology="circular"`; optional whether to also consider reverse complement.
- **Output:** Returns `DeduplicationResult` with retained `records`, equivalence `groups`, and deduplication counts.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("AACG", topology="circular"), "origin-0"),
    DNARecord(DNASequence("CGAA", topology="circular"), "origin-2"),
]
result = deduplicate(records, equivalence="circular")
print(result.records.ids, result.groups[0].rotation_offsets)
```

- **Example output:**

```text
('origin-0',) (0, 2)
```

## 4) IUPAC-aware deduplication {#4-data-004-iupac-aware-deduplication}

- **Purpose:** Interpret IUPAC ambiguity symbols as a set of bases, determine bit by bit whether two sequences are exactly the same, may be compatible, or clearly conflict, and are used for equivalence checking of data containing ambiguous bases.
- **Method:** Compare the base sets represented by IUPAC symbols position by position and apply equivalence/compatibility rules.

- **API:** `dnakit.datasets.deduplicate_iupac(records[required], config[optional])`; `config` uses `dnakit.datasets.IUPACDeduplicationConfig`.

- **Input:** Required record that may contain IUPAC characters; optional represents strategy and pairwise comparison upper limit.
- **Output:** Returns `IUPACDeduplicationResult` with equivalence groups and relations within each group.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate_iupac

records = [
    DNARecord(DNASequence("A", alphabet="iupac"), "a"),
    DNARecord(DNASequence("N", alphabet="iupac"), "n"),
    DNARecord(DNASequence("G", alphabet="iupac"), "g"),
]
result = deduplicate_iupac(records)
print([(group.member_ids, group.relation) for group in result.groups])
```

- **Example output:**

```text
[(('a', 'n'), 'compatible'), (('g',), 'singleton')]
```

## 5) Approximate deduplication {#5-data-005-approximate-deduplication}

- **Purpose:** Use Identity, Edit distance or k-mer similarity to compare sequences, merge approximate copies that exceed the threshold into groups, and return representative records and grouping basis.
- **Method:** Compute pairwise similarities, form threshold clusters and retain a representative per cluster.

- **API:** `dnakit.datasets.deduplicate_approximate(records[required], config[required])`; `config` uses `dnakit.datasets.ClusterConfig`.

- **Input:** Required record set and similarity method/threshold; optional k, canonical, representation policy and resource cap.
- **Output:** Returns `ClusteringResult` with near-duplicate groups and representative records.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, deduplicate_approximate

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = deduplicate_approximate(
    records, config=ClusterConfig(method="identity", threshold=0.7)
)
print(result.labels, result.representatives.ids)
```

- **Example output:**

```text
(0, 0, 1) ('a', 'c')
```
