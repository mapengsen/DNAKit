# DNA sketches and HMMs

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.alignment.pyhmmer_p39_functions`

| Function | Operation |
| --- | --- |
| `build_hmm_from_sequence` | Build a DNA profile HMM from one sequence and write it as a HMMER3 file (pyhmmer.plan7.Builder.build). |
| `build_sequence_ssi_index` | Index the records of a FASTA file by name in a new SSI file next to it (SSIWriter.add_file, add_key, add_alias). |
| `compute_alignment_weights` | Compute sequence weights of a DNA alignment with Easel (MSA.compute_weights, esl-weight). |
| `emit_hmm_alignment` | Sample a multiple alignment of sequences emitted from the core model of a DNA profile HMM (HMM.emit_alignment). |
| `emit_hmm_sequences` | Sample sequences from a DNA profile HMM and write them to a new FASTA file (HMM.emit_sequence, Profile.emit_sequence). |
| `fetch_sequence_by_ssi` | Fetch one record of a FASTA file by name or alias through its SSI index (SSIReader.find_name). |
| `filter_alignment_by_identity` | Remove alignment sequences whose pairwise identity reaches ``max_identity`` and write the rest to a new file (DigitalMSA.identity_filter, esl-weight -f). |
| `guess_file_alphabet` | Detect whether a sequence file or an alignment file holds DNA, RNA or amino acids (SequenceFile.guess_alphabet, MSAFile.guess_alphabet). |
| `pyhmmer_hmmscan` | Scan each DNA sequence of a FASTA file against a database of DNA profile HMMs with the HMMER 3 hmmscan pipeline (pyhmmer.hmmscan). |
| `pyhmmer_hmmsearch` | Search a DNA sequence file with one or more DNA profile HMMs using the HMMER 3 hmmsearch pipeline (pyhmmer.hmmsearch). |
| `reverse_complement_alignment` | Reverse-complement a DNA multiple alignment, including per-column annotation, and write it to a new file (DigitalMSA.reverse_complement). |

`dnakit.microbiome.sourmash_p39_functions`

| Function | Operation |
| --- | --- |
| `build_signature_index` | Build an on-disk searchable signature index of a chosen type (sourmash index -F). |
| `ingest_mash_json_signatures` | Convert JSON sketch dumps written by ``mash info -d`` into one sourmash signature file (sourmash sig ingest). |
| `murmur_hash_string` | MurmurHash64 of a string exactly as given, with sourmash's seed 42 by default (sourmash.minhash.hash_murmur). |
| `sequence_kmer_hashes` | Hash every k-mer of a DNA sequence without making a sketch (MinHash.seq_to_hashes, kmers_and_hashes). |
| `sketch_intersection_union_size` | Number of shared hashes and the size of the union of two sketches (MinHash.count_common and intersection_and_union_size). |
| `sketch_size_accuracy` | Tell whether a scaled sketch is large enough to estimate its k-mer set size (MinHash.size_is_accurate). |
| `sketch_weighted_containment` | Abundance-weighted fraction of the query sketch found in the target sketch (MinHash.contained_by_weighted). |
| `tax_annotate_gather_results` | Add a lineage column to gather CSV files (sourmash tax annotate); writes ``<name>.with-lineages.csv``. |
| `tax_genome_classification` | Assign one lineage to each query genome from gather results (sourmash tax genome). |
| `tax_grep_lineages` | Select the taxonomy lineages whose rank names match a string (sourmash tax grep). |
| `tax_prepare_taxonomy_database` | Combine taxonomy spreadsheets into one new taxonomy database, SQLite (``sql``) or CSV (sourmash tax prepare). |
| `tax_summarize_lineages` | Summarize taxonomy databases or lineage spreadsheets: number of distinct lineages per rank (sourmash tax summarize). |
| `watch_reads_against_index` | Stream FASTA/FASTQ reads into a growing sketch and stop at the first index match (sourmash watch). |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from dnakit.microbiome.sourmash_p39_functions import murmur_hash_string
print(murmur_hash_string('ACGT', verbose=False)['hash'])
```

**Example output:**

```text
2597925387403686983
```
