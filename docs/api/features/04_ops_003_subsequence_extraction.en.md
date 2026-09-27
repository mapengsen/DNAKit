# subsequence extraction

Extracts a specified fragment from a DNA sequence.

## 1) Extract subsequences

- **Purpose:** Intercept DNA fragments according to the start and end coordinates and retain the source range. It also supports interception of circular sequences across origins for extracting genes, functional regions or specified windows.
- **Method:** Slice the 0-based half-open interval [start,end); concatenate the tail and head for circular intervals crossing the origin.
- **API:** `dnakit.ops.subsequence(sequence[required], start[required], end[required], allow_gaps[optional])`, `dnakit.ops.circular_subsequence(sequence[required], start[required], end[required], allow_gaps[optional])`.
- **Input:** Required `DNASequence`, `start`, `end`; Optional `allow_gaps=True` Keep or truncate known gaps.
- **Output:** Returns the processed `DNASequence`; read its symbols from `symbols`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import subsequence

seq = DNASequence("AACCGG")
selected = subsequence(seq, 1, 5)
print(selected.symbols)
```

- **Example output:**

```text
ACCG
```

## 2) Read regions and fetch reference sequences {#bio-13}

- **Purpose:** Fetch reference FASTA segments from BED/GTF regions, including reverse complements.
- **Method:** Normalize intervals to zero-based, half-open coordinates and fetch through a pysam FASTA index.
- **API:** `dnakit.bioinformatics.reads.read_regions(path, format='bed', feature='exon')`; `dnakit.bioinformatics.reads.index_reference(path)`; `dnakit.bioinformatics.reads.IndexedReference(path).fetch(region)`
- **Input:** Reference FASTA plus BED/GTF or `Region`; requires `bio-reads`. `index_reference` creates `.fai` beside the reference.
- **Output:** `read_regions` yields `Region` objects; `index_reference` creates a `.fai` index; `IndexedReference.fetch` returns a `Read` containing the requested sequence.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import read_regions, index_reference, IndexedReference

with TemporaryDirectory() as folder:
    fasta, bed = Path(folder) / "ref.fa", Path(folder) / "regions.bed"
    fasta.write_text(">chr1\nAACCGGTT\n")
    bed.write_text("chr1\t0\t4\tregion1\t0\t-\n")
    index_reference(fasta)
    region = next(read_regions(bed))
    with IndexedReference(fasta) as reference:
        print(region.start, region.end, reference.fetch(region).sequence)
```

- **Example output:**

```text
0 4 GGTT
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 3) Calculate interval union, intersection or difference {#bio-60}

- **Purpose:** Combine region sets to extract shared or unique spans.
- **Method:** Group by reference and strand, scan endpoints and merge adjacent spans.
- **API:** `dnakit.bioinformatics.workflows.query_regions(left, right, operation='intersection')`
- **Input:** Two collections of `Region` objects using zero-based, half-open coordinates.
- **Output:** Returns a tuple of intervals, each with a reference name, start, and end.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Region
from dnakit.bioinformatics.workflows import query_regions

result = query_regions([Region("chr1", 0, 10)], [Region("chr1", 5, 15)])
print([(region.reference, region.start, region.end) for region in result])
```

- **Example output:**

```text
[('chr1', 5, 10)]
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#workflows).
