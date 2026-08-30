# Remove duplicates

Deduplication is performed on the DNA sequence to obtain the DNA sequence after deduplication.

## 1) DATA-001 · Standard deduplication

**Function:** Group the sequence by character-by-character equality, retain only one representative record in each group, and return the merged record and its representative mapping for lossless deletion of exact copies.

**API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; `config` uses `dnakit.datasets.DeduplicationConfig`.

**Input:** Required `DNARecord` collection; optional represents policy, conflicting fields/policies, and metadata merging.

**Sample code:**

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

**Example results:**

```text
('a', 'b') ('a', 'c')
```

## 2) DATA-002 · Reverse complementary deduplication

**Function:** Compare the forward and reverse complementary forms of sequences, and group records with only different strand directions but the same content into one group to avoid repeated counting caused by double-stranded directions.

**API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; this setting is `equivalence="reverse_complement"`.

**Input:** Required `DNARecord` collection; optional represents selection and metadata conflict configuration.

**Sample code:**

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

**Example results:**

```text
('forward', 'reverse') ('forward', 'reverse_complement')
```

## 3) DATA-003 · Ring equivalent deduplication

**Function:** Perform rotational equivalence comparison on circular sequences, merge records with different storage starting points but the same base sequence on the loop, and retain representative records and membership relationships.

**API:** `dnakit.datasets.deduplicate(records[required], equivalence[optional], config[optional])`; use `equivalence="circular"` or `equivalence="circular_reverse_complement"` for this item.

**Input:** Required Records explicitly marked `topology="circular"`; optional whether to also consider reverse complement.

**Sample code:**

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

**Example results:**

```text
('origin-0',) (0, 2)
```

## 4) DATA-004 · IUPAC-aware deduplication

**Function:** Interpret IUPAC fuzzy characters as a set of bases, determine bit by bit whether two sequences are exactly the same, may be compatible, or clearly conflict, and are used for equivalence checking of data containing fuzzy bases.

**API:** `dnakit.datasets.deduplicate_iupac(records[required], config[optional])`; `config` uses `dnakit.datasets.IUPACDeduplicationConfig`.

**Input:** Required record that may contain IUPAC characters; optional represents strategy and pairwise comparison upper limit.

**Sample code:**

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

**Example results:**

```text
[(('a', 'n'), 'compatible'), (('g',), 'singleton')]
```

## 5) DATA-005 · Approximate deduplication

**Function:** Use Identity, Edit distance or k-mer similarity to compare sequences, merge approximate copies that exceed the threshold into groups, and return representative records and grouping basis.

**API:** `dnakit.datasets.deduplicate_approximate(records[required], config[required])`; `config` uses `dnakit.datasets.ClusterConfig`.

**Input:** Required record set and similarity method/threshold; optional k, canonical, representation policy and resource cap.

**Sample code:**

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

**Example results:**

```text
(0, 0, 1) ('a', 'c')
```
