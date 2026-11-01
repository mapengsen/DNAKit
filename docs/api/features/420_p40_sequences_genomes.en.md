# Sequence files and genome annotation

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.annotation.pyrodigal_p40_functions`

| Function | Operation |
| --- | --- |
| `pyrodigal_candidate_nodes` | Extract (and optionally score) the candidate start/stop nodes of a sequence. |
| `pyrodigal_codon_probability` | Start- or stop-codon probability expected from the GC content of a sequence. |
| `pyrodigal_find_genes` | Predict genes with pyrodigal and write the chosen annotation files. |
| `pyrodigal_gc_frame_plot` | Maximum-GC frame plot of one sequence (pyrodigal Sequence.max_gc_frame_plot). |
| `pyrodigal_shine_dalgarno` | Score the Shine-Dalgarno motif upstream of a start codon (pyrodigal Sequence.shine_dalgarno). |
| `pyrodigal_train_gene_finder` | Train the Prodigal gene finder on a genome FASTA (pyrodigal GeneFinder.train, ``pyrodigal -t``). |

`dnakit.sequence.pyfaidx_p40_functions`

| Function | Operation |
| --- | --- |
| `faidx_overwrite_region` | Overwrite a FASTA region with a sequence of the same length (pyfaidx Faidx.to_file). |
| `faidx_read_index` | Read an existing FASTA ``.fai`` index (pyfaidx Faidx.read_fai) or BGZF ``.gzi`` index (Faidx.read_gzi). |
| `faidx_write_index` | Build a FASTA index and write it to a new file. |

`dnakit.sequence.pyfastx_p40_functions`

| Function | Operation |
| --- | --- |
| `pyfastx_count_sequences` | Count the FASTA sequences that are at least ``min_length`` bp long (pyfastx Fasta.count). |
| `pyfastx_extract_records` | Extract FASTA sequences or FASTQ reads by name into a new file (``pyfastx extract``). |
| `pyfastx_list_sequence_names` | List the sequence names of a FASTA file (pyfastx Fasta.keys). |
| `pyfastx_nl_statistic` | Compute the N50/L50-style statistic of a FASTA file (pyfastx Fasta.nl). |
| `pyfastx_split_fastq` | Split a FASTQ file into several files (``pyfastx split``). |

`dnakit.sequence.dnaio_p40_functions`

| Function | Operation |
| --- | --- |
| `dnaio_read_parallel_files` | Read three or more parallel FASTA/FASTQ files record by record (dnaio.open with several files). |
| `dnaio_reverse_complement_reads` | Reverse-complement every read of a FASTA/FASTQ file (dnaio SequenceRecord.reverse_complement). |
| `dnaio_split_into_chunks` | Split a FASTA/FASTQ file into chunks of complete records, without parsing them (dnaio.read_chunks). |
| `dnaio_write_reads` | Write parallel read lists with dnaio: paired-end files, parallel FASTA/FASTQ files or interleaved reads. |

`dnakit.databases.genomepy_p40_functions`

| Function | Operation |
| --- | --- |
| `genomepy_list_installed` | List genomes installed under an explicit local directory, with no provider lookup. |
| `genomepy_import_local_annotation` | Copy a user-supplied local GTF/BED annotation to a new path and parse it with genomepy. |
| `genomepy_annotation_attributes` | List GTF attribute names using Annotation.attributes (upstream feature sampling). |
| `genomepy_annotation_genes` | List unique GTF gene_name values or BED names with Annotation.genes. |
| `genomepy_annotation_lengths` | Gene/transcript exon lengths with Annotation.lengths, preserving upstream coordinate arithmetic. |
| `genomepy_annotation_dictionary` | Make a key-value dictionary from GTF columns/attributes with Annotation.gtf_dict. |
| `genomepy_sanitize_annotation` | Match GTF contigs against full FASTA headers, then drop unmatched rows using genomepy sanitize. |
| `genomepy_download_blacklist` | Download and decompress the Boyle Lab blacklist for an assembly supported by genomepy (network). |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.sequence.pyfastx_p40_functions import pyfastx_count_sequences
with TemporaryDirectory() as folder:
    fasta = Path(folder) / 'example.fa'
    fasta.write_text('>a\nACGT\n>b\nN\n')
    print(pyfastx_count_sequences(fasta, min_length=2, verbose=False)['sequence_count'])
```

**Example output:**

```text
1
```

`pyrodigal_gc_frame_plot` preserves upstream frame `-1` for bases with no assigned frame, counts those separately in `unknown_frame_count`, and reports omitted runs with `segments_truncated`. `dnaio_split_into_chunks` recompresses chunks when the output filename has a compression suffix.
