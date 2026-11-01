# Genomic intervals and assembly

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.annotation.genomicranges_p40_functions`

| Function | Operation |
| --- | --- |
| `genomicranges_binned_average` | Average the BED score (column 5) of ``source`` over each interval of ``bins`` (GenomicRanges.binned_average). |
| `genomicranges_disjoint_bins` | Assign each interval to a bin so that intervals within one bin do not overlap (GenomicRanges.disjoint_bins). |
| `genomicranges_interval_coverage` | Per-position coverage depth of the intervals on each chromosome (GenomicRanges.coverage). |
| `genomicranges_is_disjoint` | Test whether no two intervals overlap on the same chromosome (and strand) (GenomicRanges.is_disjoint). |
| `genomicranges_match_intervals` | Find, for every query interval, the identical interval in ``source`` (GenomicRanges.match). |
| `genomicranges_narrow_intervals` | Narrow every interval by relative ``start``, ``width`` or ``end`` (GenomicRanges.narrow). |
| `genomicranges_order_intervals` | Indices that sort the intervals by chromosome, strand and position (GenomicRanges.order). |
| `genomicranges_out_of_bound_intervals` | Find intervals that extend past the chromosome length (GenomicRanges.get_out_of_bound_index). |
| `genomicranges_rank_intervals` | Rank of each interval in coordinate order (GenomicRanges.rank). |
| `genomicranges_restrict_intervals` | Restrict intervals to a coordinate window (genomicranges GenomicRanges.restrict). |
| `genomicranges_span_intervals` | Span (minimum start to maximum end) of the intervals per chromosome and strand (GenomicRanges.range). |
| `genomicranges_tile_intervals` | Tile every interval into ``n`` equal sub-intervals or into pieces of ``width`` (GenomicRanges.tile). |

`dnakit.annotation.polars_bio_p40_functions`

| Function | Operation |
| --- | --- |
| `read_bigwig_intervals` | Read the intervals of a BigWig signal track with polars-bio (read_bigwig). |
| `read_msa_records` | Read an A2M or A3M alignment file with polars-bio (read_a2m, read_a3m). |

`dnakit.bioinformatics.bioframe_p40_functions`

| Function | Operation |
| --- | --- |
| `bioframe_fetch_ucsc_table` | Fetch a UCSC assembly table through bioframe's UCSCClient. |
| `bioframe_is_complete_ucsc_string` | Test whether a UCSC-style region string names a chromosome with a known end (bioframe.is_complete_ucsc_string). |
| `bioframe_load_fasta` | Open an indexed FASTA lazily with bioframe.load_fasta and report its records. |

`dnakit.assembly.spades_p40_cli_functions`

| Function | Operation |
| --- | --- |
| `spades_assemble` | Assemble reads with SPAdes in one of its modes (``spades.py`` and the ``*spades.py`` wrappers). |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from dnakit.bioinformatics.bioframe_p40_functions import bioframe_is_complete_ucsc_string
print(bioframe_is_complete_ucsc_string('chr1:10-20', verbose=False)['is_complete'])
```

**Example output:**

```text
True
```
