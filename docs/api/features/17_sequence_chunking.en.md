# OPS-010 sequence segmentation

**Function:** Cut long DNA into fragments with source coordinates according to a fixed window, overlapping sliding window, random interval or multi-scale strategy, and return them one by one with an iterator for model input, partition analysis and large sequence flow processing.

Split DNA sequences according to fixed length, sliding window, random interval, or multi-scale strategies to generate sequence fragments with coordinates.

All splits use 0-based, half-open interval coordinates, and lazily return `SequenceChunk`. The entire result will not be read into memory at once.

<span id="1"></span>**Sequence Segmentation Overview**

The API entrance is as follows:

- `dnakit.iter_sequence_chunks(value[required], config[optional], start[optional], end[optional], source_id[optional], split[optional], region_index[optional], progress[optional])`: split into single `DNASequence` or `DNARecord`;
- `dnakit.iter_fasta_chunks(source[required], config[optional], bed[optional], read_config[optional], progress[optional])`: Read and split FASTA one by one;
- `dnakit.ChunkingConfig(strategy[optional], length[optional], step[optional], min_length[optional], max_length[optional], num_samples[optional], lengths[optional], steps[optional], stage_steps[optional], include_partial[optional], seed[optional], split[optional], allow_gaps[optional])`: Define the segmentation strategy;
- `dnakit.make_length_curriculum(lengths[required], stage_steps[optional])`: Create a training plan of short length first and then long length;
- `dnakit.LengthCurriculum.to_config(window_step[optional], include_partial[optional], seed[optional], split[optional], allow_gaps[optional])`: Convert length plan to split configuration;
- `dnakit.SequenceChunk.to_record()`: Convert the segmentation result to `DNARecord`, no calling parameters.

All strategies return a lazy `SequenceChunk` and retain provenance information such as source ID, start and end coordinates, split, strategy and scale numbers.

## 1) Fixed length, no overlap <span id="2"></span> {#fixed-non-overlapping-chunks}

**Function:** Generate fragments of the same length and non-overlapping from the starting point of the sequence, and return the source coordinates of each fragment, suitable for establishing training samples or analysis units without repeated coverage.

This is the default policy for `OPS-010 Sequence Chunking`. Each window length is specified by `length`; windows do not overlap, and an incomplete trailing window is discarded by default.

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(length=1024)
for chunk in iter_fasta_chunks("genome.fa", config=config):
    print(chunk.id, chunk.source_start, chunk.source_end)
```

If you need to keep a window with insufficient length at the end, set `include_partial=True`.

## 2) Fixed length, overlapping sliding window <span id="3"></span> {#overlapping-sliding-chunks}

**Function:** Generate overlapping fragments according to the specified window length and step size, so that bases near the boundaries appear in multiple windows, suitable for local predictions that require continuous context coverage.

Use `strategy="sliding"` and `step` to set the distance between the start points of the windows. Windows overlap when `step` is smaller than `length`.

```python
from dnakit import ChunkingConfig, iter_sequence_chunks

config = ChunkingConfig(strategy="sliding", length=1024, step=512)
for chunk in iter_sequence_chunks(sequence, config=config):
    print(chunk.source_start, chunk.source_end)
```

`step` cannot be larger than the window length, otherwise uncovered gaps will occur.

## 3) Random interval <span id="4"></span> {#random-interval-chunks}

**Function:** Extract a specified number of random fragments within the allowed length range and coordinate range, and ensure reproduction through seed. It is suitable for random training sampling or sampling inspection of long sequences.

The random strategy generates a specified number of random intervals within the range `min_length` and `max_length`. Setting `seed` gives reproducible results.

```python
from dnakit import ChunkingConfig, iter_sequence_chunks

config = ChunkingConfig(
    strategy="random",
    min_length=512,
    max_length=2048,
    num_samples=100,
    seed=19,
)
chunks = iter_sequence_chunks(sequence, config=config)
```

## 4) Multi-scale <span id="5"></span> {#multiscale-chunks}

**Function:** Use multiple window lengths to generate short, medium, and long scale fragments and mark scale numbers, which is used to simultaneously model local motifs and longer range contexts.

The multi-scale strategy generates windows for each length in `lengths`, suitable for extracting both local and longer-range context.

```python
from dnakit import ChunkingConfig, iter_fasta_chunks

config = ChunkingConfig(
    strategy="multiscale",
    lengths=(1024, 4096, 16384),
)
chunks = iter_fasta_chunks("genome.fa", config=config, bed="splits.bed")
for chunk in chunks:
    print(chunk.id, chunk.level_index, chunk.requested_length)
```

The scales do not overlap by default. You can also specify the step size for each scale through `steps`.

## 5) Short first and then long <span id="6"></span> {#length-curriculum-chunks}

**Function:** Define the sequence length and number of continuous steps used in the training phase, and then convert it into the corresponding segmentation configuration, which is used to implement the curriculum data plan that gradually increases from short segments to long segments.

Use `make_length_curriculum()` to create a training plan that gradually increases sequence length, and then use `to_config()` to convert to a split configuration.

```python
from dnakit import iter_fasta_chunks, make_length_curriculum

curriculum = make_length_curriculum(
    (1024, 4096, 16384),
    stage_steps=(1000, 2000, 3000),
)
config = curriculum.to_config()
chunks = iter_fasta_chunks("genome.fa", config=config)
```

This plan only describes the segmentation phases and does not automatically start training.

<span id="7-fastabed"></span>**FASTA, BED and result metadata**

`iter_sequence_chunks()` is used to split a single piece `DNASequence` or `DNARecord`; `iter_fasta_chunks()` is used to read and split FASTA one by one.

When `bed` is not passed, all FASTA records use the `train` tag by default; when BED is passed in, the first column of the BED is required to be exactly the same as the FASTA ID, and the fourth column can use `train`, `valid`, `test` and other tags.

Each `SequenceChunk` contains `source_id`, `source_start`, `source_end`, `split`, `strategy`, and `level_index`. After calling `to_record()`, the source ID, start and end coordinates, split and scale numbers are written to the record metadata.

Pass in `progress=callback` to receive `started`, `yielded`, `completed` events. The current FASTA reader retains one FASTA record at a time; when handling large chromosome-level records, adjust the input limit for `ReadConfig` based on actual memory.
