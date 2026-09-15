# Shared inputs and outputs for data processing

Public standardization, sequence-operation, deduplication, clustering, and splitting APIs accept
DNA text, file paths, and existing DNA objects. File suffixes identify FASTA/FNA, TXT, CSV/TSV,
and other DNAKit I/O formats, including gzip. Existing object calls retain their return types,
coordinate conventions, random parameters, and algorithmic limits.

```python
from dnakit.ops import reverse_complement, trim, transcribe, translate
from dnakit.standardize import normalize

if __name__ == "__main__":
    print(reverse_complement("AATCGC"))  # GCGATT
    normalize("input.txt")
    reverse_complement("hg38.fna")
    trim("input.csv", left=2, right=2, sequence_column="sequence")
    transcribe("input.fna", output_path="rna.fa")
    translate("input.csv", output_path="proteins.csv")
```

## Input rules

- DNA text is treated as one sequence. Reverse complementation returns a string; other operations retain their native result types.
- A `str` containing a dot or path separator is generally treated as a path. Use `Path` for extensionless files or names that could be mistaken for DNA. Standardization still accepts punctuation in raw text and switches to file mode only for recognizable sequence-file suffixes.
- TXT defaults to one sequence per non-empty line. Use `ReadConfig(text_mode="sequence")` when the entire file is one wrapped sequence.
- CSV defaults to the `sequence` column; change it with `sequence_column="dna"`. Labels and extra columns enter record metadata and are written back with their original names and order. CSV cells retain string semantics. Conflicting values in metadata and extra columns with the same name cause an error.
- Uncompressed inputs use their actual file size to set default reading limits, so hg38 needs no manual size configuration. Gzip defaults to a 20 GB decompressed reading limit. Explicit `read_config=ReadConfig(...)` values take precedence.
- Circular operations require `topology="circular"` or circular topology declared in `ReadConfig` or the DNA object.

The underlying `read_inputs()` API supports custom workflows and closes its files:

```python
from dnakit.io import read_inputs

with read_inputs(["first.fna", "second.csv", "AATCGC"]) as records:
    for record in records:
        print(record.id, record.sequence.symbol_length)
```

## File output rules

Per-record operations default to `original_operation.original_suffix` in the input directory;
set `output_path` to override it. They return `WriteResult`: the output path is in
`result.target_artifact.relative_path` and the record count in `result.record_count`.
Records are read incrementally, with start, record-progress, and completion messages enabled by
default. Set `show_progress=False` to disable them.

The input cannot also be the output, and existing results are not overwritten by default.
Use `write_config=WriteConfig(overwrite=True)` to regenerate them. Incomplete output is not
installed after a write failure. Coordinates, random seeds, N/IUPAC support, and optional
dependencies continue to follow the underlying algorithm.

Transcription and translation write RNA and protein text, respectively, supporting FASTA, TXT,
CSV/TSV, and JSONL. FASTA sources default to `.fa` outputs. Automatically derived CSV outputs
remove DNA structural columns that do not apply to RNA/proteins. Complex features or per-base
annotations that cannot be mapped reliably cause explicit errors; record-aware editing,
reverse complementation, circular transforms, and chunking use their corresponding mappings.

## Covered public APIs

| Category | Entry points and processing units |
| --- | --- |
| Standardization | `standardize.normalize`, `normalize_gaps`; the former can first clean invalid raw characters from FASTA/FNA, TXT, and CSV/TSV |
| Validation | `standardize.validate` writes per-record JSONL reports; `validate_set` writes a dataset report |
| Direction | `ops.reverse`, `complement`, `reverse_complement` |
| Transcription and translation | `ops.transcribe`, `translate` |
| Extraction and editing | `ops.subsequence`, `circular_subsequence`, `insert`, `delete`, `substitute`, `mask`, `trim` |
| Generation | `ops.mutate`, `evolution_generate`, `indel_generate`, `rearrange_generate`, `kmer_shuffle` |
| Concatenation | `ops.concat`, `concat_overlap` follow file and record order; `crossover` pairs records from two inputs in order and rejects unequal counts |
| Circular operations | `ops.rotate`, `canonical_origin` |
| Record operations | Editing, reverse-complement, and circular `*_record` entry points exported by `ops` |
| Chunking | New `chunking.chunk_sequences`; existing `iter_sequence_chunks` and `iter_fasta_chunks` retain lazy iteration |
| Deduplication | `datasets.deduplicate`, `deduplicate_iupac`, `deduplicate_approximate` |
| Clustering | `datasets.cluster_sequences`, `hierarchical_cluster`, `neural_cluster_sequences`, `select_representatives` |
| Preparation and splitting | `datasets.exclude_by_metadata`, `exclude_chromosomes`, `exclude_species`, `split`, `temporal_split`, `joint_split` |
| Dataset evaluation | `datasets.detect_leakage` accepts a mapping of subset names to input files; `evaluate_split_quality` accepts files and its existing assignments argument |

Structured auxiliary inputs such as AGP assemblies and BED regions retain their dedicated APIs
and are not guessed to be DNA text.

## Deduplication, clustering, and splitting

Collection operations accept a file or a list of files. Original record IDs are retained; when
combining files, ensure IDs meet the algorithm's uniqueness requirements. These algorithms operate
on datasets and retain their existing memory and record limits. Except for random file splitting,
file mode returns the native result object and writes reports/data to `original_operation/`, or
the directory specified by `output_path`. `result.json` stores the result. CSV subsets and
representatives retain original table columns; other inputs use annotation-preserving JSONL.
The directory is installed only after all writes succeed. `overwrite=True` replaces the result
directory, so use a dedicated output directory.

```python
from dnakit.datasets import SplitConfig, deduplicate, split

if __name__ == "__main__":
    dedup = deduplicate("sequences.csv")
    result = split(
        "sequences.csv",
        config=SplitConfig(ratios={"train": 0.8, "valid": 0.1, "test": 0.1}, seed=42),
        output_path="splits",
    )
    # splits/train.csv, valid.csv, test.csv, result.json
```

Write sequence chunks directly without a manual loop:

```python
from dnakit.chunking import ChunkingConfig, chunk_sequences

result = chunk_sequences(
    "hg38.fna", config=ChunkingConfig(length=200), output_path="windows.fna"
)
```

## Multicore processing

Shared processing APIs accept `workers=-1` (default: all available logical cores), `workers=1`
(serial), or a positive integer limiting worker processes. CPU affinity is respected; concurrency
is also limited by task count and available memory.

```python
from dnakit.ops import reverse_complement

if __name__ == "__main__":
    result = reverse_complement("hg38.fna", workers=8)
```

Place processing calls inside `if __name__ == "__main__":` in ordinary Python scripts to avoid
re-execution when child processes start. Use `workers=1` in interactive environments or with
non-serializable custom callbacks. Paper examples `02_reverse_complement.py` and
`03_split_dataset.py` set `workers` inside the script. The CLI uses a global option:

```bash
dnakit --workers 8 deduplicate input.fna output.fna
```

Place the global CLI `--workers` option before the subcommand. It controls supported CPU steps.

| Computational step | Execution strategy |
| --- | --- |
| Sequence operations, standardization, and validation on large, multi-record uncompressed FASTA | With default reading settings, partition complete record regions; workers read/write directly and the parent merges in order. Custom reading, compression, or output-size limits use the general path |
| Reversal, complementation, and reverse complementation of contiguous in-memory sequences | C string operations avoid interprocess copying; large FASTA reverse complementation also has a base-region processing path |
| Exact/reverse-complement deduplication | Compute lightweight equivalence keys locally to avoid copying sequences between processes |
| Circular deduplication | Compute equivalence keys in parallel and merge groups sequentially |
| IUPAC deduplication | Compute positional compatibility in parallel while retaining stable grouping rules |
| Threshold/hierarchical clustering, approximate deduplication, representative selection | Compute each sequence's k-mer set or binary fingerprint once; other methods still compute pairwise similarities in parallel, followed by sequential grouping/merging |
| Similarity-based and joint splitting | Parallel fingerprints or pairwise similarities, followed by sequential subset assignment |
| Crossover between two files | Pair rows in parallel while preserving output order |
| Concatenation, lazy windows over in-memory objects, simple metadata filtering, random/temporal assignment itself | Retain streaming or sequential processing for lightweight or order-dependent steps |

Neural clustering and inference retain backend-managed CPU/GPU execution; this setting does not
replicate GPU models automatically. More cores do not guarantee faster execution. Short tasks run
serially. Pending work is bounded, so very large records do not accumulate without limit in process
queues. Multiprocess output is merged in input order; random seeds follow the original algorithm,
and scheduling order does not enter random-number generation. Explicit mutable `rng` objects retain
serial invocation and state advancement; custom chunk-progress callbacks remain in the calling
process. Single-worker file operations use direct streaming instead of unnecessary shard merging.
Large FASTA direction operations retain their specialized region path; such files up to 512 MiB
use at most eight workers to reduce startup costs and disk contention.

## Memory usage for random file splitting

`split(file_or_files, config=SplitConfig(method="random", ...))` preferentially indexes
uncompressed FASTA with default reader settings, retaining IDs, offsets, and assignments. It reads
records on demand and writes subsets in parallel without an intermediate sequence copy. Default
JSONL export also groups records so that even the train subset can use multiple processes.
Regular FASTA layouts allow JSONL record positions to be computed in advance; workers write
non-overlapping regions directly, avoiding another complete output copy. Irregular layouts fall
back to general writing. Other formats, custom reader configurations, and layouts unsuitable for
indexing retain disk staging.

The entire genome is not held as sequence objects at once; each worker still needs room for its
largest record and conversion buffers. The final directory is installed only after every subset
is written successfully. Temporary data is cleaned after success or failure. Other splitting
methods and native object calls retain their algorithmic and memory semantics.

Random file splitting returns `FileSplitResult`. `counts`, `assignments`, `seed`, and `to_dict()`
retain the original split audit information. `get("train")` and `subsets[i].records` return
re-iterable `FileRecords` that read the written files on demand. Use `DNASet(result.get("train"))`
only when the entire subset must be loaded. Ratios, record order, and assignments for `seed=42`
match the in-memory version; `preserve_order=False` also retains the original ordering semantics.

Random file splitting displays percentages, data volume, speed, elapsed time, and ETA for reading
and each subset write by default. Reading measures source-file bytes; writing measures the data
volume represented by completed records. ETA estimates the current stage, not the entire task.
It is unknown at startup before speed samples exist and during file verification.

Existing download/extraction, sequence-processing, prediction, representation, evaluation, and
workflow progress bars display `elapsed` and `ETA`. Known totals use observed completion speed;
streaming with unknown totals and external commands neither pre-read data nor guess total duration,
so ETA remains unknown. Speed samples reset when downloads switch to extraction or another file.
Existing `verbose=False`, `show_progress=False`, and custom callbacks still apply; scripts need
no repeated time-display parameters.

### Corrupt compressed inputs

The shared reader converts EOFError from truncated gzip into
`InputFormatError(code="INPUT_READ_FAILED")`, retaining the source and underlying cause.
Standardization's raw-text reader likewise handles compression, decoding, and CSV parsing errors.
File operations using these readers retain staged output installation, so failure does not replace
an existing complete output.
