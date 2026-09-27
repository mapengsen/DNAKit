# Sequencing read processing

## 1) Pair reads by ID {#bio-10}

- **Purpose:** Match R1/R2 records by name and identify orphans.
- **Method:** Match after removing trailing `/1` or `/2`; a missing mate is `None`.
- **API:** `dnakit.bioinformatics.reads.pair_reads(first, second, max_pending=100000, include_orphans=True)`
- **Input:** Two iterators of `Read` objects with unique normalized IDs in each mate stream.
- **Output:** Returns an iterator of `ReadPair` entries with a pair ID and both reads.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Read, pair_reads

pairs = pair_reads([Read("r1/1", "ACGT")], [Read("r1/2", "TGCA")])
print([(pair.id, pair.first.sequence, pair.second.sequence) for pair in pairs])
```

- **Example output:**

```text
[('r1', 'ACGT', 'TGCA')]
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 2) Recover valid records from damaged FASTQ {#bio-11}

- **Purpose:** Recover recognizable four-line FASTQ records from a damaged file.
- **Method:** Scan for headers, sequences, separators and matching quality lengths, counting discarded lines.
- **API:** `dnakit.bioinformatics.reads.repair_fastq(source, destination, phred_offset=33)`
- **Input:** A damaged FASTQ and a new output path; use strict `read_fastx` for multiline FASTQ.
- **Output:** Returns `RepairResult` with the repaired FASTQ path, recovered-read count, and discarded-line count.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import repair_fastq

with TemporaryDirectory() as folder:
    source = Path(folder) / "broken.fq"
    source.write_text("invalid\n@r1\nACGT\n+\nIIII\n")
    result = repair_fastq(source, Path(folder) / "repaired.fq")
    print(result.recovered_reads, result.discarded_lines)
```

- **Example output:**

```text
1 1
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 3) Sort, intersect, rename and fingerprint reads {#bio-12}

- **Purpose:** Organize records, extract shared sequences and compute a sequence digest.
- **Method:** Sort or shuffle with a seed, intersect by ID/sequence, rename with regex, and hash length-delimited sequences.
- **API:** `dnakit.bioinformatics.reads.arrange_reads(reads, operation='sort')`; `dnakit.bioinformatics.reads.common_reads(collections)`; `dnakit.bioinformatics.reads.rename_reads(reads, pattern, replacement)`; `dnakit.bioinformatics.reads.sequence_digest(reads)`
- **Input:** Collections of `Read` objects; the digest depends on record order, not IDs.
- **Output:** Returns sorted, intersected, or renamed `Read` records; the check function returns a comparison result.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Read, arrange_reads, common_reads, rename_reads, sequence_digest

reads = [Read("b", "TGCA"), Read("a", "ACGT")]
ordered = arrange_reads(reads)
shared = common_reads([ordered, [Read("a", "ACGT")]])
renamed = list(rename_reads(shared, "^a$", "sample_a"))
print([read.id for read in ordered])
print([read.id for read in renamed])
print(sequence_digest(shared) == sequence_digest(renamed))
```

- **Example output:**

```text
['a', 'b']
['sample_a']
True
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 4) Prepare, mark duplicates and index BAM {#bio-16}

- **Purpose:** Prepare SAM/BAM for variant calling.
- **Method:** Run samtools name-sort, fixmate, coordinate sort, markdup, index and stats.
- **API:** `dnakit.bioinformatics.ngs.prepare_bam(source, output_directory, mark_duplicates=True, config=None)`
- **Input:** Existing `aligned.sam` or BAM, samtools on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import prepare_bam

result = prepare_bam("aligned.sam", "new_bam")
print(result.artifacts["aligned.bam"])
```

- **Example output:**

Returns `ToolResult` containing `aligned.bam`, its index and statistics files. This external-program example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).
