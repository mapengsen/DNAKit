# sequence segmentation

All chunk coordinates are 0-based and half-open. Results are yielded lazily.

## 1) Fixed length, no overlap <span id="2"></span> {#fixed-non-overlapping-chunks}

- **Purpose:** Generate nonoverlapping fixed-length fragments.
- **Method:** Advance by length from the start; discard incomplete trailing chunks unless include_partial=True.
- **API:** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **Input:** A DNASequence or DNARecord and the lengths/step/sampling settings shown below.
- **Output:** Returns an iterator of `SequenceChunk` entries with sequence and start/end coordinates.
- **Example code:**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(length=4)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **Example output:**

```text
[(0, 4, 'AACC'), (4, 8, 'GGTT')]
```

## 2) Fixed length, overlapping sliding window <span id="3"></span> {#overlapping-sliding-chunks}

- **Purpose:** Generate overlapping windows at a fixed step.
- **Method:** Advance by step with 1 ≤ step ≤ length.
- **API:** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **Input:** A DNASequence or DNARecord and the lengths/step/sampling settings shown below.
- **Output:** Returns an iterator of `SequenceChunk` entries with overlapping window sequences and coordinates.
- **Example code:**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="sliding", length=4, step=2)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **Example output:**

```text
[(0, 4, 'AACC'), (2, 6, 'CCGG'), (4, 8, 'GGTT')]
```

## 3) Random interval <span id="4"></span> {#random-interval-chunks}

- **Purpose:** Sample reproducible random fragments.
- **Method:** Seed the random source and sample num_samples intervals in the requested length range; repeated intervals are possible.
- **API:** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **Input:** A DNASequence or DNARecord and the lengths/step/sampling settings shown below.
- **Output:** Returns an iterator of `SequenceChunk` entries with sampled sequence and coordinates.
- **Example code:**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="random", min_length=2, max_length=4, num_samples=3, seed=19)
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **Example output:**

```text
[(0, 4, 'AACC'), (0, 4, 'AACC'), (1, 5, 'ACCG')]
```

## 4) Multi-scale <span id="5"></span> {#multiscale-chunks}

- **Purpose:** Split the same sequence at multiple lengths.
- **Method:** Generate windows in lengths order; windows are nonoverlapping within each scale by default.
- **API:** `dnakit.iter_sequence_chunks(value, config=ChunkingConfig(...))`
- **Input:** A DNASequence or DNARecord and the lengths/step/sampling settings shown below.
- **Output:** Returns an iterator of `SequenceChunk` entries with window length, sequence, and coordinates.
- **Example code:**

```python
from dnakit import DNASequence, ChunkingConfig, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
config = ChunkingConfig(strategy="multiscale", lengths=(2, 4))
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **Example output:**

```text
[(0, 2, 'AA'), (2, 4, 'CC'), (4, 6, 'GG'), (6, 8, 'TT'), (0, 4, 'AACC'), (4, 8, 'GGTT')]
```

## 5) Short first and then long <span id="6"></span> {#length-curriculum-chunks}

- **Purpose:** Create a short-to-long sequence-length curriculum.
- **Method:** Pair lengths with stage_steps and call to_config(); the plan does not start training.
- **API:** `dnakit.make_length_curriculum(lengths, stage_steps=...)` → `curriculum.to_config()` → `iter_sequence_chunks(...)`
- **Input:** A DNASequence or DNARecord and the lengths/step/sampling settings shown below.
- **Output:** Returns `LengthCurriculum` with selected chunk lengths by stage.
- **Example code:**

```python
from dnakit import DNASequence, make_length_curriculum, iter_sequence_chunks

sequence = DNASequence("AACCGGTT")
curriculum = make_length_curriculum((2, 4), stage_steps=(10, 20))
config = curriculum.to_config()
chunks = iter_sequence_chunks(sequence, config=config)
print([(c.source_start, c.source_end, c.sequence.symbols) for c in chunks])
```

- **Example output:**

```text
[(0, 2, 'AA'), (2, 4, 'CC'), (4, 6, 'GG'), (6, 8, 'TT'), (0, 4, 'AACC'), (4, 8, 'GGTT')]
```

<span id="1"></span><span id="7-fastabed"></span>[FAQ: details](../../faq/tutorial_basics.md#chunking-details)

## Paired strands and multi-format export {#paired-chunk-export}

`dnakit.chunking.iter_paired_chunks(source, reverse_source=None, config=None, read_config=None, source_ids=None)` reads forward FASTA and yields adjacent forward/reverse-complement row dictionaries. Fields are `id, sequence, pair_id, source_id, source_start, source_end, strand`. Both strands share a coordinate-based pair ID and use zero-based, half-open coordinates on the forward source.

An optional `reverse_source` must contain whole-record reverse complements with matching IDs, order and lengths. These structural conditions are checked; callers are responsible for supplying actual reverse complements. A forward interval `[start:end)` corresponds to `[L-end:L-start)` in that file. When omitted, reverse complements are computed from each chunk. `source_ids` filters FASTA IDs. For file paths without an explicit `read_config`, the read limits are sized to the larger input file; streams may need an explicit `ReadConfig`. The reader holds one chromosome pair at a time, without collecting all output chunks.

`export_paired_chunks(source, target, ..., overwrite=False)` writes the same iterator to one or several paths and returns `pair_count`, `record_count`, the configured `length`, and `coordinates`. CSV, TSV, JSON, JSONL, Parquet and FASTA suffixes are supported, including empty outputs. Tables retain every field; FASTA retains IDs and sequences.

```python
from dnakit.chunking import ChunkingConfig, export_paired_chunks
from dnakit.io import save_json

result = export_paired_chunks(
    "normalized.fasta",
    ["output/chunks.csv", "output/chunks.fasta"],
    reverse_source="reverse_complement.fasta",
    config=ChunkingConfig(length=1000),
    overwrite=True,
)
save_json(result, "output/result.json", overwrite=True)
```

Each pair consists of two separate 1000 bp records, totaling 2000 bp. File-path read limits adapt to hg38 automatically; callers may still set `read_config` for streams or special inputs. Normalization, reverse complementation and chunking remain separate operations, while parsing and serialization are handled by the library.
