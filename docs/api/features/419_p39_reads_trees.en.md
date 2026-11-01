# Reads and phylogenetic trees

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.fasta_fastq.cutadapt_p39_functions`

| Function | Operation |
| --- | --- |
| `cutadapt_filter_reads` | Filter reads by N content, expected errors, average error rate, adapter presence or CASAVA chastity. |
| `cutadapt_rewrite_read_names` | Rewrite read headers using a rename template, length tag or stripped suffix. |
| `cutadapt_clamp_negative_qualities` | Clamp negative Phred qualities to zero with cutadapt --zero-cap. |
| `cutadapt_interleaved_reads` | Read and write interleaved paired-end FASTQ, optionally trimming each mate's 3' adapter. |
| `cutadapt_demultiplex_reads` | Single-index or combinatorial dual-index demultiplexing into a new directory. |

`dnakit.phylogeny.toytree_p39_functions`

| Function | Operation |
| --- | --- |
| `toytree_rf_distance` | Robinson-Foulds bipartition distance between two Newick trees with the same tips. |
| `toytree_consensus` | Majority-rule consensus; min_frequency=0.5 keeps splits present in at least half the trees. |
| `toytree_neighbor_joining` | Infer a neighbor-joining tree from labelled pairwise distances. |
| `toytree_upgma` | Infer a UPGMA tree from labelled pairwise distances. |
| `toytree_mrca` | Return the index and descendant tip labels of the most recent common ancestor. |
| `toytree_monophyly` | Test whether selected tips form a monophyletic group in a rooted Newick tree. |
| `toytree_unique_topologies` | Count distinct rooted or unrooted topologies in a Newick tree sample. |
| `toytree_local_genealogy` | Extract one local ToyTree genealogy from a local tskit file, by tree index or genomic site. |
| `toytree_k80_matrix` | K80 nucleotide transition matrix in toytree's upstream base order. |
| `toytree_tn93_matrix` | TN93 nucleotide transition matrix using upstream parameter and base ordering. |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from dnakit.phylogeny.toytree_p39_functions import toytree_monophyly
print(toytree_monophyly('((a,b),c);', ['a', 'b'], verbose=False)['monophyletic'])
```

**Example output:**

```text
True
```
