# Sequence distance and alignment

Compare two DNA sequences: sequence distance "how much they differ", pairwise sequence alignment "how the bases correspond".

## 1. Sequence distance

### 1.1 SIM-006 · Hamming distance

**Function:** Compare two equal-length sequences position by position and return the number and specific positions of different bases, which can be used to count point mutations or quickly determine the differences between equal-length sequences.

**API:** `dnakit.similarity.hamming_distance(left[required], right[required], max_distance[optional])`.

**Input:** Two equal-length sequences are required; optional `max_distance` audit threshold.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import hamming_distance

result = hamming_distance(DNASequence("ACGT"), DNASequence("ACCT"))
print(result.distance, [item.position for item in result.mismatches])
```

**Example results:**

```text
1.0 [2]
```

### 1.2 SIM-007 · Edit distance

**Function:** Calculate the minimum replacement, insertion and deletion costs required to change one sequence into another sequence, and return the optimal editing path to quantify the difference between sequences of different lengths.

**API:** `dnakit.similarity.edit_distance(left[required], right[required], substitution_cost[optional], insertion_cost[optional], deletion_cost[optional], max_distance[optional], return_path[optional], max_cells[optional])`.

**Input:** Two sequences are required; optional operation cost, maximum distance, editing path and DP upper limit.

**Sample code:**

```python
from dnakit import DNASequence
from dnakit.similarity import edit_distance

result = edit_distance(DNASequence("ACGT"), DNASequence("AGT"), return_path=True)
print(result.distance, [step.operation for step in result.edit_path or ()])
```

**Example results:**

```text
1.0 ['match', 'delete', 'match', 'match']
```

## 2. Pairwise sequence alignment (SIM-008)

Pairwise sequence alignment adds gaps at appropriate positions, generates two sequences that correspond to each other, and calculates score, identity, and coverage.

**API:** `dnakit.alignment.align_pairwise(query[required], target[required], config[optional])`.

**Configuration:** Use `dnakit.alignment.AlignmentConfig` to set the pattern, match/mismatch score, linear or affine gap score and `max_cells`, return `dnakit.alignment.AlignmentResult`.

### 2.1 Global comparison

**Function:** Align two complete sequences from beginning to end, and return the corresponding relationship, score, identity and coverage after adding Gap. It is suitable for comparing sequences with similar lengths and overall correlation.

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("AGT"),
    config=AlignmentConfig(mode="global", gap_score=-1),
)
print(result.aligned_query, result.aligned_target, result.score, result.identity)
```

```text
ACGT A-GT 2.0 0.75
```

### 2.2 Local comparison

**Function:** Find the local alignment area with the highest score in the two sequences, and return the coordinates, Identity and Coverage of the fragment in the two original sequences, which is suitable for finding common fragments.

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("TTACGTAA"),
    DNASequence("GGACGTCC"),
    config=AlignmentConfig(mode="local", mismatch_score=-2, gap_score=-2),
)
print(
    result.aligned_query,
    result.aligned_target,
    result.identity,
    result.query_coverage,
    result.target_coverage,
)
```

```text
ACGT ACGT 1.0 0.5 0.5
```

### 2.3 Semi-global comparison

**Function:** No Gap points will be deducted for unaligned regions at both ends, and the core alignment and its original sequence coordinates will be returned. It is suitable for comparing primers, amplicons or short sequences with target regions in long sequences.

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("TTACGTAA"),
    config=AlignmentConfig(mode="semi_global", mismatch_score=-2, gap_score=-2),
)
print(
    result.aligned_query,
    result.aligned_target,
    result.query_coverage,
    result.target_coverage,
    (result.target_start, result.target_end),
)
```

```text
ACGT ACGT 1.0 0.5 (2, 6)
```
