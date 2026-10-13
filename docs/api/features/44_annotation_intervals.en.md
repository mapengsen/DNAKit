# Annotation coordinates and gene models

These functions operate only on explicit annotation intervals. They do not read a reference genome or query an annotation database. Each output records its coordinate convention so the operation can be reproduced in a methods section.

For interval-processing methodology, cite the [pybedtools paper](https://doi.org/10.1093/bioinformatics/btr539). The implementation here is dependency-free deterministic interval arithmetic; report coordinate conventions, versions, and input data separately.

## 1) `convert_annotation_coordinates`

- **Purpose:** Convert explicit intervals between `1-based-closed` and `0-based-half-open` coordinates.
- **Method:** Apply the deterministic start/end conversion to every record while preserving other fields.
- **API:** `dnakit.annotation.convert_annotation_coordinates(input, *, options=None, verbose=True)`
- **Input:** `{"records": [{"chrom": "chr1", "start": 2, "end": 5}]}`; set `options.coordinate_system` and `target_coordinate_system` explicitly.
- **Output:** A serializable completion envelope with converted `records`, both coordinate conventions, and provenance.
- **Example code:**

```python
from dnakit.annotation import convert_annotation_coordinates

result = convert_annotation_coordinates(
    {"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
    options={
        "coordinate_system": "1-based-closed",
        "target_coordinate_system": "0-based-half-open",
    },
    verbose=False,
)
print(result["result"]["records"])
```

- **Example output:**

```json
[{"chrom": "chr1", "start": 1, "end": 5, "id": "g1", "coordinate_system": "0-based-half-open"}]
```

## 2) `make_intergenic_regions`

- **Purpose:** Generate intergenic gaps within declared contig boundaries.
- **Method:** Sort and merge overlapping or adjacent gene intervals, then emit the gaps before, between, and after them.
- **API:** `dnakit.annotation.make_intergenic_regions(input, *, options=None, verbose=True)`
- **Input:** `input.genes` contains `chrom/start/end` intervals and `input.contig_lengths` declares contig lengths; the default is `1-based-closed`.
- **Output:** `result.regions` plus counts, contig lengths, and the merge rule.
- **Example code:**

```python
from dnakit.annotation import make_intergenic_regions

result = make_intergenic_regions(
    {"genes": [{"chrom": "chr1", "start": 2, "end": 4}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed"},
    verbose=False,
)
print(result["result"]["regions"])
```

- **Example output:**

```json
[
  {"id": "intergenic_chr1_1", "chrom": "chr1", "start": 1, "end": 1, "feature_type": "intergenic", "coordinate_system": "1-based-closed"},
  {"id": "intergenic_chr1_2", "chrom": "chr1", "start": 5, "end": 10, "feature_type": "intergenic", "coordinate_system": "1-based-closed"}
]
```

## 3) `make_promoter_regions`

- **Purpose:** Generate strand-aware promoter windows around explicit genes.
- **Method:** Use `start` as the plus-strand TSS and `end` as the minus-strand TSS, then clip to declared contig boundaries.
- **API:** `dnakit.annotation.make_promoter_regions(input, *, options=None, verbose=True)`
- **Input:** `input.genes` contains `chrom/start/end/strand`; `options.upstream` and `downstream` are non-negative window lengths.
- **Output:** `result.promoters` with gene IDs, strand, window lengths, and coordinate convention.
- **Example code:**

```python
from dnakit.annotation import make_promoter_regions

result = make_promoter_regions(
    {"genes": [{"id": "g1", "chrom": "chr1", "start": 4, "end": 7, "strand": "-"}], "contig_lengths": {"chr1": 12}},
    options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1},
    verbose=False,
)
print(result["result"]["promoters"])
```

- **Example output:**

```json
[{"id": "promoter_g1", "gene_id": "g1", "chrom": "chr1", "start": 7, "end": 9, "strand": "-", "feature_type": "promoter", "coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1}]
```

## 4) `merge_gene_models`

- **Purpose:** Merge explicit exon or segment intervals by gene, chromosome, and strand.
- **Method:** Sort intervals for each gene and merge overlaps; set `merge_adjacent=True` to merge touching intervals too.
- **API:** `dnakit.annotation.merge_gene_models(input, *, options=None, verbose=True)`
- **Input:** `input.segments` contains `gene_id` (or `gene`/`id`) and `chrom/start/end`, with optional `strand`.
- **Output:** `result.models` with model bounds, merged `segments`, and segment counts.
- **Example code:**

```python
from dnakit.annotation import merge_gene_models

result = merge_gene_models(
    {"segments": [{"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3}, {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5}]},
    options={"coordinate_system": "1-based-closed", "merge_adjacent": True},
    verbose=False,
)
print(result["result"]["models"])
```

- **Example output:**

```json
[{"gene_id": "g1", "gene_name": "g1", "chrom": "chr1", "strand": "+", "start": 1, "end": 5, "segments": [{"start": 1, "end": 5}], "segment_count": 1, "coordinate_system": "1-based-closed"}]
```
