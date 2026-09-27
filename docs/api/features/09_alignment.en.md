# Sequence distance and alignment

Compare two DNA sequences: sequence distance "how much they differ", pairwise sequence alignment "how the bases correspond".

<span id="1-sequence-distance"></span>

## 1) Hamming distance {#11-sim-006-hamming-distance}

- **Purpose:** Compare two equal-length sequences position by position and return the number and specific positions of different bases, which can be used to count point mutations or quickly determine the differences between equal-length sequences.
- **Method:** Count differing characters at corresponding positions; sequences must have equal lengths.

- **API:** `dnakit.similarity.hamming_distance(left[required], right[required], max_distance[optional])`.

- **Input:** Two equal-length sequences are required; optional `max_distance` audit threshold.
- **Output:** Returns a `DistanceResult`; the example shows `distance`, `mismatches`.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import hamming_distance

result = hamming_distance(DNASequence("ACGT"), DNASequence("ACCT"))
print(result.distance, [item.position for item in result.mismatches])
```

- **Example output:**

```text
1.0 [2]
```

## 2) Edit distance {#12-sim-007-edit-distance}

- **Purpose:** Calculate the minimum replacement, insertion and deletion costs required to change one sequence into another sequence, and return the optimal editing path to quantify the difference between sequences of different lengths.
- **Method:** Use dynamic programming to minimize substitution, insertion and deletion costs, optionally tracing the edit path.

- **API:** `dnakit.similarity.edit_distance(left[required], right[required], substitution_cost[optional], insertion_cost[optional], deletion_cost[optional], max_distance[optional], return_path[optional], max_cells[optional])`.

- **Input:** Two sequences are required; optional operation cost, maximum distance, editing path and DP upper limit.
- **Output:** Returns a `DistanceResult`; the example shows `distance`, `edit_path`.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.similarity import edit_distance

result = edit_distance(DNASequence("ACGT"), DNASequence("AGT"), return_path=True)
print(result.distance, [step.operation for step in result.edit_path or ()])
```

- **Example output:**

```text
1.0 ['match', 'delete', 'match', 'match']
```

<span id="2-pairwise-sequence-alignment-sim-008"></span>

## 3) Global comparison {#21-global-comparison}

- **Purpose:** Align two complete sequences from beginning to end, and return the corresponding relationship, score, identity and coverage after adding Gap. It is suitable for comparing sequences with similar lengths and overall correlation.
- **Method:** Use dynamic programming to optimize match, mismatch and gap scores over both complete sequences.

- **API:** `dnakit.alignment.align_pairwise(query[required], target[required], config[optional])`.
- **Input:** Two DNA sequences; optional `AlignmentConfig(mode="global")` sets match/mismatch scores, linear or affine gap scores and `max_cells`. Returns `AlignmentResult` with aligned sequences, score, identity and coverage.
- **Output:** Returns `AlignmentResult` with aligned sequences, score, and coordinates.
- **Example code:**

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

- **Example output:**

```text
ACGT A-GT 2.0 0.75
```

## 4) Local comparison {#22-local-comparison}

- **Purpose:** Find the local alignment area with the highest score in the two sequences, and return the coordinates, Identity and Coverage of the fragment in the two original sequences, which is suitable for finding common fragments.
- **Method:** Allow dynamic-programming scores to restart at zero and trace the highest-scoring local segment.

- **API:** `dnakit.alignment.align_pairwise(query[required], target[required], config[optional])`.
- **Input:** Two DNA sequences; optional `AlignmentConfig(mode="local")` sets match/mismatch scores, linear or affine gap scores and `max_cells`. Returns `AlignmentResult` with aligned sequences, score, identity and coverage.
- **Output:** Returns `AlignmentResult` with aligned sequences, score, and coordinates.
- **Example code:**

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

- **Example output:**

```text
ACGT ACGT 1.0 0.5 0.5
```

## 5) Semi-global comparison {#23-semi-global-comparison}

- **Purpose:** No Gap points will be deducted for unaligned regions at both ends, and the core alignment and its original sequence coordinates will be returned. It is suitable for comparing primers, amplicons or short sequences with target regions in long sequences.
- **Method:** Optimize alignment with free terminal gaps and trace the aligned core.

- **API:** `dnakit.alignment.align_pairwise(query[required], target[required], config[optional])`.
- **Input:** Two DNA sequences; optional `AlignmentConfig(mode="semi_global")` sets match/mismatch scores, linear or affine gap scores and `max_cells`. Returns `AlignmentResult` with aligned sequences, score, identity and coverage.
- **Output:** Returns `AlignmentResult` with aligned sequences, score, and coordinates.
- **Example code:**

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

- **Example output:**

```text
ACGT ACGT 1.0 0.5 (2, 6)
```

## 6) Calculate distances from a multiple alignment {#bio-1}

- **Purpose:** Convert aligned sequences into a pairwise distance matrix.
- **Method:** Compute pairwise mismatch fractions, excluding gap positions by default.
- **API:** `dnakit.bioinformatics.alignment.MultipleAlignment.from_sequences(sequences)`; `MultipleAlignment.distances(gap_policy='pairwise_delete')`
- **Input:** Equal-width sequences with unique IDs; gaps (`-`) are allowed.
- **Output:** Returns a pairwise distance matrix for the multiple alignment.
- **Example code:**

```python
from dnakit.bioinformatics.alignment import MultipleAlignment

msa = MultipleAlignment.from_sequences({"a": "ACGT", "b": "ACTT", "c": "TCGT"})
print(msa.ids)
print(msa.distances().tolist())
```

- **Example output:**

```text
('a', 'b', 'c')
[[0.0, 0.25, 0.25], [0.25, 0.0, 0.5], [0.25, 0.5, 0.0]]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

## 7) Align sequences with MAFFT, MUSCLE or Clustal Omega {#bio-3}

- **Purpose:** Create a multiple alignment from unaligned DNA sequences.
- **Method:** Run the selected executable, then validate returned IDs and sequence content.
- **API:** `dnakit.bioinformatics.alignment.align_sequences(sequences, output_directory, method='mafft', config=None)`
- **Input:** `DNA([...])` or an ID-to-sequence mapping and a new output directory; requires MAFFT, MUSCLE 5 or Clustal Omega on PATH and `bio-formats`.
- **Output:** Returns the `MultipleAlignment` generated by the external tool, containing aligned sequences.
- **Example code:**

```python
from dnakit import DNA
from dnakit.bioinformatics.alignment import align_sequences

print(msa.ids, msa.sequences)
```

- **Example output:**

Returns `MultipleAlignment`: `ids` preserves input IDs and `sequences` contains equal-width aligned sequences. Gap placement depends on the external program; this example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

## 8) Map reads to a reference genome {#bio-15}

- **Purpose:** Map single-end or paired-end reads and produce SAM.
- **Method:** Build a private index in a new directory, then run BWA MEM or Bowtie2.
- **API:** `dnakit.bioinformatics.ngs.map_reads(reference, first, output_directory, second=None, method='bwa', config=None)`
- **Input:** Reference FASTA, R1/R2 FASTQ, the chosen aligner on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Example code:**

```python
from dnakit.bioinformatics import ExecutionConfig
from dnakit.bioinformatics.ngs import map_reads

result = map_reads("reference.fa", "R1.fq.gz", "new_mapping",
                   second="R2.fq.gz", config=ExecutionConfig(threads=4))
print(result.artifacts["aligned.sam"])
```

- **Example output:**

Returns `ToolResult`; `artifacts['aligned.sam']` is the SAM path. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).
