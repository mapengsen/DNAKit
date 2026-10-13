# G18d review: kipoi, genomic-benchmarks and gReLU

This page records the G18d review. The ledgers have 121 rows with status `other` or `cap`: 86 in kipoi, 3 in genomic-benchmarks and 32 in gReLU (31 `other` and 1 `cap`). After the review every row has a final status and a one-line reason. This group adds 20 gReLU wrappers in `src/dnakit/predictions/grelu_review_functions.py`. kipoi and genomic-benchmarks have no offline DNA-level function; their modules `kipoi_review_functions.py` and `genomic_benchmarks_review_functions.py` export nothing.

- The 86 kipoi rows are plugin registries, YAML and file helpers, working-directory context managers, CLI argument parsers and writers, none of them DNA-level. `kipoi.list_sources` and `kipoi.config.list_sources` clone the 202 MB model repository on first use, so they are `needs-network`. Tests replace them with mocks that fail if called, and they are never run live. The existing `list_kipoi_models` and `describe_kipoi_model` are not wrapped again.
- The 3 genomic-benchmarks rows: a label-order helper for the cached test split (no sequence content, filesystem-order dependent) is `not-dna`; two cache helpers that delete a generated dataset folder are `internal`.
- Of the 32 gReLU rows, 20 offline DNA-level functions are wrapped. `filter_cells` is `not-dna`; `update_ckpt`, `check_bigwig`, `check_fasta`, `filter_intervals`, `filter_obs` and two pairs of resource path helpers are `internal`; two rows are `duplicate-of` wrapped upstream objects.

## Review results

| Package | Other/cap rows before | New status | New wrapped |
| --- | --- | --- | --- |
| kipoi | 86 | internal 58, not-dna 24, needs-network 2, cli 2 | 0 |
| genomic_benchmarks | 3 | internal 2, not-dna 1 | 0 |
| grelu | 32 (31 other + 1 cap) | wrapped 20 (including the cap row), duplicate-of 2, internal 9, not-dna 1 | 20 |

## Upstream issues and how the wrappers handle them

- `grelu.io.bed.read_narrowpeak`: upstream treats the first line as a header, so the first peak of a headerless file is lost. The wrapper defaults to `has_header=False`.
- `grelu.variant.predict_variant_effects`: its `check_reference` flag shadows a function of the same name, so upstream's check cannot run. The wrapper checks reference bases against the local FASTA itself; `pos` is 1-based.
- `grelu.interpret.simulate.marginalize_patterns`: raw output and `n_shuffles=0` both fail with a concat error on the installed torch. The wrapper requires `compare_func` to be subtract or divide and `n_shuffles` to be at least 1.
- `grelu.design.ledidi`: upstream passes a 0-dimensional target, which ledidi 2.2.0 rejects with an IndexError. The wrapper calls `ledidi.Ledidi` directly with a one-element target and the same loss.
- `grelu.data.preprocess.filter_coverage`: percentile mode is missing its quantile argument and cannot run; negative sampling depends on a global seed. Neither is exposed.
- `grelu.design.evolve`: the model must carry task names in `data_params['tasks']['name']`; the wrapper rejects a model without them.
- `grelu.interpret.modisco.run_modisco`: random sequences usually give no patterns, and then nothing is written; the result has `patterns_found` False.

## Remaining unwrapped rows

- kipoi: the two `list_sources` rows are `needs-network`; the rest are `not-dna`, `internal` or `cli` (see above).
- genomic-benchmarks: `labels_in_order` is `not-dna`; the two `remove_dataset_from_disk` rows are `internal`.
- gReLU: `filter_cells` is `not-dna`; `filter_intervals`, `filter_obs`, `update_ckpt`, `check_bigwig`, `check_fasta` and the two resource path pairs are `internal`; `data.dataset.check_intervals` and `data.preprocess.get_chromosomes` are `duplicate-of`.

## Wrapped functions

The test file is `others/tests/unit/test_review_G18d.py`. Each wrapped function is checked against the real upstream call or a hand calculation. The example outputs below come from running each example.

## 1) `read_narrowpeak_intervals`

- **Purpose:** Read a narrowPeak (BED6+4) peak file in file order. Each peak keeps its chrom, start, end, summit, signal, p-value and q-value.
- **Method:** Calls `grelu.io.bed.read_narrowpeak`. Upstream treats the first line as a header by default and silently drops the first peak of a headerless file, so the wrapper defaults to `has_header=False` and passes `header=None` to pandas. Files with fewer than ten columns are rejected.
- **API:** `read_narrowpeak_intervals(path[required], skiprows[optional], has_header[optional], verbose[optional])`
- **Input:** Path to a local narrowPeak file.
- **Output:** A dictionary with `count`, `peaks` (one dictionary per peak) and `provenance`.

**Example file:** `peaks.narrowPeak`

```text
chr1	100	200	p1	500	.	5.5	3.2	1.1	50
chr2	300	420	p2	80	+	2.0	-1	-1	60
```

**Example code:**

```python
from dnakit.predictions import read_narrowpeak_intervals

result = read_narrowpeak_intervals("peaks.narrowPeak")
print(result["count"])
for peak in result["peaks"]:
    print(peak["chrom"], peak["start"], peak["end"], peak["summit"], peak["signal"])
```

**Example output:**

```text
2
chr1 100 200 50 5.5
chr2 300 420 60 2.0
```

## 2) `expand_chromosome_names`

- **Purpose:** Expand a chromosome shortcut (`autosomes`, `autosomesX`, `autosomesXY`) into GRCh38-style names, or check an explicit list.
- **Method:** Calls `grelu.data.utils.get_chromosomes`. Upstream returns an unknown shortcut unchanged as a string, so the wrapper rejects unknown shortcuts instead.
- **API:** `expand_chromosome_names(chroms[required], verbose[optional])`
- **Input:** A shortcut string, or a list of chromosome name strings.
- **Output:** A dictionary with `count`, `chromosomes` and `provenance`.

**Example code:**

```python
from dnakit.predictions import expand_chromosome_names

print(expand_chromosome_names("autosomesXY")["chromosomes"][-3:])
print(expand_chromosome_names(["chrM", "chr1"])["count"])
```

**Example output:**

```text
['chr22', 'chrX', 'chrY']
2
```

## 3) `check_genomic_interval_table`

- **Purpose:** Check that an interval table has gReLU's genomic-interval layout: the first three columns are `chrom`, `start`, `end`, and `start` and `end` are integers.
- **Method:** Calls `grelu.sequence.format.check_intervals`. Rows are built in the given order, so a wrong column order or a float coordinate makes `valid` False.
- **API:** `check_genomic_interval_table(intervals[required], verbose[optional])`
- **Input:** A list of mappings, each with `chrom`, `start` and `end`.
- **Output:** A dictionary with `valid`, `columns`, `dtypes`, `row_count` and `provenance`.

**Example code:**

```python
from dnakit.predictions import check_genomic_interval_table

good = [{"chrom": "chr1", "start": 1, "end": 5}]
bad = [{"chrom": "chr1", "start": 1.5, "end": 5}]
print(check_genomic_interval_table(good)["valid"])
print(check_genomic_interval_table(bad)["dtypes"]["start"])
print(check_genomic_interval_table(bad)["valid"])
```

**Example output:**

```text
True
float64
False
```

## 4) `filter_interval_coverage`

- **Purpose:** Keep genomic intervals by coverage across tasks: `cutoff` keeps intervals at or above a threshold, `top` keeps the n intervals with the highest coverage.
- **Method:** Calls `grelu.data.preprocess.filter_coverage`. The per-interval value is the aggregate over tasks (mean, sum, max or min). Upstream percentile mode is broken and negative sampling needs a seed, so neither is exposed.
- **API:** `filter_interval_coverage(coverage[required], intervals[required], tasks[optional], aggfunc[optional], method[optional], cutoff[optional], verbose[optional])`
- **Input:** coverage has one row per task and one column per interval; intervals is a list of interval mappings.
- **Output:** A dictionary with `count`, `intervals` (with original indices) and `provenance`.

**Example code:**

```python
from dnakit.predictions import filter_interval_coverage

intervals = [
    {"chrom": "chr1", "start": 0, "end": 5},
    {"chrom": "chr1", "start": 5, "end": 9},
    {"chrom": "chr2", "start": 1, "end": 4},
]
coverage = [[0.5, 3.0, 1.0]]
cutoff = filter_interval_coverage(coverage, intervals, cutoff=1.0)
print(cutoff["count"], [item["index"] for item in cutoff["intervals"]])
top = filter_interval_coverage(coverage, intervals, method="top", cutoff=1)
print([item["chrom"] + ":" + str(item["start"]) for item in top["intervals"]])
```

**Example output:**

```text
2 [1, 2]
['chr1:5']
```

## 5) `add_negative_intervals`

- **Purpose:** Append negative control intervals after the positive ones. Negative columns have zero coverage, and the `class` field labels each interval as positive or negative.
- **Method:** Calls `grelu.data.preprocess.add_negatives` with `negative_labels=0`, the only value upstream handles.
- **API:** `add_negative_intervals(coverage[required], intervals[required], negatives[required], tasks[optional], verbose[optional])`
- **Input:** coverage is a task-by-positive-interval matrix; intervals are the positives; negatives are the negative intervals.
- **Output:** A dictionary with `count`, `positive_count`, `negative_count`, `intervals`, `coverage` and `provenance`.

**Example code:**

```python
from dnakit.predictions import add_negative_intervals

positives = [{"chrom": "chr1", "start": 0, "end": 5}, {"chrom": "chr1", "start": 5, "end": 9}]
negatives = [{"chrom": "chr3", "start": 0, "end": 7}]
result = add_negative_intervals([[1.0, 2.0]], positives, negatives)
print(result["count"], result["negative_count"])
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["class"])
print(result["coverage"])
```

**Example output:**

```text
3 1
chr1 0 5 positive
chr1 5 9 positive
chr3 0 7 negative
[[1.0, 2.0, 0.0]]
```

## 6) `read_tomtom_matches`

- **Purpose:** Read every `tomtom.tsv` file in a TOMTOM output folder and keep motif matches with a q-value up to the threshold.
- **Method:** Calls `grelu.io.read_tomtom`, which drops rows with missing values as upstream does. Matches are sorted by query, target, orientation and offset. A folder without a `tomtom.tsv` file is rejected.
- **API:** `read_tomtom_matches(tomtom_dir[required], qthresh[optional], verbose[optional])`
- **Input:** A local folder that contains a `tomtom.tsv` file.
- **Output:** A dictionary with `count`, `matches` (one dictionary per row) and `provenance`.

**Example file:** `tomtom.tsv`

```text
Target_ID	Optimal_offset	p-value	E-value	q-value	Overlap	Query_consensus	Target_consensus	Orientation
MotifA	0	0.001	0.01	0.02	8	ACGTACGT	ACGTACGA	+
MotifB	1	0.2	0.5	0.4	7	ACGTACGT	TTTTACGT	-
```

**Example code:**

```python
from dnakit.predictions import read_tomtom_matches

result = read_tomtom_matches(".", qthresh=0.05)
print(result["count"])
print(result["matches"][0]["Target_ID"], result["matches"][0]["q-value"])
```

**Example output:**

```text
1
MotifA 0.02
```

## 7) `read_modisco_motifs`

- **Purpose:** Read TF-MoDISco motifs from an H5 report as trimmed position probability matrices (PPMs).
- **Method:** Calls `grelu.io.motifs.read_modisco_report`. `group` selects pos or neg, `names` selects by name and overrides `group`, and `trim_threshold` defaults to 0.3.
- **API:** `read_modisco_motifs(h5_path[required], group[optional], names[optional], trim_threshold[optional], verbose[optional])`
- **Input:** Path to a local TF-MoDISco H5 file.
- **Output:** A dictionary with `count`, `motifs` (each with `name`, `length` and `ppm`) and `provenance`.

**Example code:**

```python
import h5py
import numpy as np

with h5py.File("modisco_report.h5", "w") as handle:
    group = handle.create_group("pos_patterns").create_group("pattern_0")
    group.create_dataset("sequence", data=np.tile([4.0, 0.0, 0.0, 0.0], (6, 1)))
    group.create_dataset("contrib_scores", data=np.ones((6, 4)))
from dnakit.predictions import read_modisco_motifs

result = read_modisco_motifs("modisco_report.h5")
print(result["count"], result["motifs"][0]["name"], result["motifs"][0]["length"])
```

**Example output:**

```text
1 pos_pattern_0 6
```

## 8) `score_motif_sites`

- **Purpose:** Score each motif hit by the mean attribution over the hit interval (site attribution score).
- **Method:** Calls `grelu.interpret.motifs.score_sites`. Attributions and sequences are matched by list order or by key, and each hit must lie inside its sequence.
- **API:** `score_motif_sites(hits[required], sequences[required], attributions[required], verbose[optional])`
- **Input:** hits is the hit list from `scan_motifs_with_pwms`; sequences are the scanned sequences; attributions hold one 4 x L matrix per sequence.
- **Output:** A dictionary with `count`, `hits` (each with `site_attr_score`) and `provenance`.

**Example code:**

```python
from dnakit.predictions import score_motif_sites

sequences = {"s0": "ACGTAC"}
attributions = {"s0": [[1, 1, 1, 1, 0, 0], [0] * 6, [0] * 6, [0] * 6]}
hits = [{"motif": "m1", "sequence": "s0", "start": 0, "end": 4, "strand": "+"}]
result = score_motif_sites(hits, sequences, attributions)
print(result["hits"][0]["site_attr_score"])
```

**Example output:**

```text
0.25
```

## 9) `score_motif_pwm_products`

- **Purpose:** Score each motif hit by the sum of the element-wise product of its PWM and the attributions. The score is high when the motif shape fits the attribution profile.
- **Method:** Calls `grelu.interpret.motifs.score_motifs`. Minus-strand hits use the reverse-complement PWM, and the PWM length must equal the hit length.
- **API:** `score_motif_pwm_products(hits[required], sequences[required], attributions[required], motifs[required], verbose[optional])`
- **Input:** As for `score_motif_sites`, plus motifs: a mapping from motif name to a 4 x L probability matrix.
- **Output:** A dictionary with `count`, `hits` (each with `motif_attr_score`) and `provenance`.

**Example code:**

```python
from dnakit.predictions import score_motif_pwm_products

pwm = [[1, 1], [0, 0], [0, 0], [0, 0]]
attributions = [[[2, 3, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]
hits = [{"motif": "m1", "sequence": "seq0", "start": 0, "end": 2, "strand": "+"}]
result = score_motif_pwm_products(hits, ["ACGT"], attributions, {"m1": pwm})
print(result["hits"][0]["motif_attr_score"])
```

**Example output:**

```text
2.5
```

## 10) `run_modisco_motifs`

- **Purpose:** Discover sequence patterns from a model's attributions with TF-MoDISco. If patterns are found, the report `modisco_report.h5` and its logos are written to a new output folder.
- **Method:** Calls `grelu.interpret.modisco.run_modisco`. Attributions for one task are computed first, then TF-MoDISco runs. `out_dir` must be absent or empty. Random sequences usually give no patterns, and then nothing is written.
- **API:** `run_modisco_motifs(model[required], sequences[required], out_dir[required], task_index[optional], window[optional], method[optional], n_shuffles[optional], seed[optional], verbose[optional])`
- **Input:** A LightningModel object and a list of equal-length DNA sequences; out_dir is the output folder.
- **Output:** A dictionary with `patterns_found`, `motifs`, `files` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
from dnakit.predictions import run_modisco_motifs

rng = random.Random(1)
sequences = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(40)]
result = run_modisco_motifs(model, sequences, out_dir="modisco_out", window=60, n_shuffles=2, seed=0)
print(result["patterns_found"], result["files"])
```

**Example output:**

```text
False []
```

## 11) `get_sequence_attributions`

- **Purpose:** Compute per-base attributions of one task for equal-length DNA sequences, shaped (B, 4, L).
- **Method:** Calls `grelu.interpret.score.get_attributions` with deepshap, saliency, inputxgradient or integratedgradients. DeepSHAP dinucleotide shuffles use the global random state, so the wrapper restores NumPy, torch and Python random state afterwards.
- **API:** `get_sequence_attributions(model[required], sequences[required], task_index[optional], method[optional], n_shuffles[optional], seed[optional], correct_grad[optional], verbose[optional])`
- **Input:** A LightningModel object and equal-length DNA sequences.
- **Output:** A dictionary with `count`, `length`, `attributions` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import get_sequence_attributions

result = get_sequence_attributions(model, [SEQ], method="saliency")
print(result["count"], result["length"])
print([round(value, 4) for value in result["attributions"][0][0][:5]])
```

**Example output:**

```text
1 60
[0.0029, 0.0025, -0.0017, -0.0013, -0.0002]
```

## 12) `in_silico_mutagenesis_scores`

- **Purpose:** In silico mutagenesis: for each position and each alternative base, the change of one task's prediction.
- **Method:** Calls `grelu.interpret.score.ISM_predict`. `compare_func` is subtract, divide or log2FC (mutant against reference). Under subtract the reference base scores zero.
- **API:** `in_silico_mutagenesis_scores(model[required], sequence[required], task_index[optional], compare_func[optional], start[optional], end[optional], batch_size[optional], verbose[optional])`
- **Input:** A LightningModel object and one DNA sequence; start and end bound the positions (end is exclusive).
- **Output:** A dictionary with `positions`, `bases` and `scores` (one row per position, columns A, C, G, T).

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import in_silico_mutagenesis_scores

result = in_silico_mutagenesis_scores(model, SEQ, compare_func="subtract", start=0, end=3)
print(result["positions"], result["bases"])
for row in result["scores"]:
    print([round(value, 4) for value in row])
```

**Example output:**

```text
3 ['A', 'C', 'G', 'T']
[0.0, -0.0029, -0.0032, -0.0023]
[0.0032, 0.0, 0.0023, 0.0009]
[0.0009, -0.0019, 0.0, -0.0026]
```

## 13) `get_attention_scores`

- **Purpose:** Return a transformer model's attention maps for one sequence.
- **Method:** Calls `grelu.interpret.score.get_attention_scores`. The model needs `n_transformers` of at least 1. With `block_idx` unset, all layers are returned.
- **API:** `get_attention_scores(model[required], sequence[required], block_idx[optional], verbose[optional])`
- **Input:** A LightningModel with transformer blocks and one DNA sequence.
- **Output:** A dictionary with `shape`, `attention` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
attention_model = LightningModel(
    model_params={
        "model_type": "ConvTransformerModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "n_transformers": 1,
        "n_heads": 2,
        "key_len": 4,
    },
    train_params={"task": "regression", "loss": "mse"},
)
attention_model.eval()
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import get_attention_scores

result = get_attention_scores(attention_model, SEQ, block_idx=0)
print(result["shape"])
print(round(sum(result["attention"][0][0]), 4))
```

**Example output:**

```text
[2, 60, 60]
1.0
```

## 14) `marginalize_pattern_insertions`

- **Purpose:** Insert DNA patterns into background sequences (optionally dinucleotide-shuffled) and compare predictions before and after insertion.
- **Method:** Calls `grelu.interpret.simulate.marginalize_patterns`. `compare_func` is subtract or divide, and `n_shuffles` must be at least 1. Upstream raw output and `n_shuffles=0` both fail with a concat error on the installed torch, so they are not offered.
- **API:** `marginalize_pattern_insertions(model[required], sequences[required], patterns[required], task_index[optional], compare_func[optional], n_shuffles[optional], seed[optional], rc[optional], batch_size[optional], verbose[optional])`
- **Input:** A LightningModel, a list of background sequences and a list of patterns; each pattern must be shorter than every background.
- **Output:** A dictionary with `shape`, `effects` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import marginalize_pattern_insertions

result = marginalize_pattern_insertions(model, [SEQ], ["GATTACA"], compare_func="subtract", n_shuffles=2, seed=0)
print(result["shape"])
values = [value for shuffle in result["effects"][0] for pattern in shuffle for task in pattern for value in task]
print([round(value, 4) for value in values])
```

**Example output:**

```text
[1, 2, 1, 1, 1]
[-0.0051, -0.0008]
```

## 15) `marginalize_pattern_spacing`

- **Purpose:** Study how the spacing between two patterns affects prediction: the fixed pattern sits at the centre and the moving pattern is placed every `stride` bases.
- **Method:** Calls `grelu.interpret.simulate.marginalize_pattern_spacing`, returning the spacing list and the prediction change at each spacing.
- **API:** `marginalize_pattern_spacing(model[required], sequence[required], fixed_pattern[required], moving_pattern[required], task_index[optional], stride[optional], n_shuffles[optional], seed[optional], rc[optional], compare_func[optional], batch_size[optional], verbose[optional])`
- **Input:** A LightningModel, one background sequence and two patterns.
- **Output:** A dictionary with `spacings`, `shape`, `effects` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import marginalize_pattern_spacing

result = marginalize_pattern_spacing(model, SEQ, "GATTACA", "CCC", stride=20, n_shuffles=1, compare_func="subtract")
print(result["spacings"], result["shape"])
```

**Example output:**

```text
[-26, -6, 14] [1, 1, 3, 1, 1]
```

## 16) `shuffle_tile_effects`

- **Purpose:** Cut equal-length sequences into tiles, dinucleotide-shuffle each tile in turn and measure the prediction change, to locate regulatory elements.
- **Method:** Calls `grelu.interpret.simulate.shuffle_tiles`. Tile starts are every `stride` bases (default `tile_len`). `protect_center` leaves a central window unshuffled.
- **API:** `shuffle_tile_effects(model[required], sequences[required], tile_len[required], task_index[optional], stride[optional], protect_center[optional], n_shuffles[optional], seed[optional], compare_func[optional], batch_size[optional], verbose[optional])`
- **Input:** A LightningModel and a list of equal-length DNA sequences.
- **Output:** A dictionary with `tile_count`, `tiles` (each with start and end), `effects` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import shuffle_tile_effects

result = shuffle_tile_effects(model, [SEQ], 20, stride=20, n_shuffles=1, compare_func="subtract")
print(result["tile_count"], [(tile["start"], tile["end"]) for tile in result["tiles"]])
```

**Example output:**

```text
3 [(0, 20), (20, 40), (40, 60)]
```

## 17) `grelu_predict_variant_effects`

- **Purpose:** Predict the effect of each SNV or indel on one task, using a local reference genome.
- **Method:** Calls `grelu.variant.predict_variant_effects`. The genome is a local FASTA with `.fai` and `.sizes` files. Upstream's reference-allele check cannot run because a flag shadows a function of the same name, so the wrapper checks reference bases against the FASTA itself (`pos` is 1-based).
- **API:** `grelu_predict_variant_effects(model[required], variants[required], genome_fasta[required], task_index[optional], seq_len[optional], compare_func[optional], rc[optional], max_seq_shift[optional], check_reference[optional], batch_size[optional], verbose[optional])`
- **Input:** A LightningModel, variant mappings (chrom, pos, ref, alt) and the path to a local FASTA.
- **Output:** A dictionary with `count`, `variants`, `effects` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
rng = random.Random(9)
sequence = "".join(rng.choice("ACGT") for _ in range(300))
with open("genome.fa", "w") as handle:
    handle.write(">chr1\n" + sequence + "\n")
with open("genome.fa.fai", "w") as handle:
    handle.write("chr1\t300\t6\t300\t301\n")
with open("genome.fa.sizes", "w") as handle:
    handle.write("chr1\t300\n")
from dnakit.predictions import grelu_predict_variant_effects

ref = sequence[99]
alt = "A" if ref != "A" else "C"
result = grelu_predict_variant_effects(model, [{"chrom": "chr1", "pos": 100, "ref": ref, "alt": alt}], "genome.fa", seq_len=60)
print(result["count"], round(result["effects"][0][0], 4))
```

**Example output:**

```text
1 0.9891
```

## 18) `grelu_marginalize_variants`

- **Purpose:** Run a marginalization test for each variant: compare the reference and alternate windows against shuffled backgrounds.
- **Method:** Calls `grelu.variant.marginalize_variants`. The genome and variant rules match `grelu_predict_variant_effects`. The output is upstream's result table.
- **API:** `grelu_marginalize_variants(model[required], variants[required], genome_fasta[required], task_index[optional], seq_len[optional], n_shuffles[optional], seed[optional], compare_func[optional], rc[optional], max_seq_shift[optional], verbose[optional])`
- **Input:** A LightningModel, variant mappings and the path to a local FASTA.
- **Output:** A dictionary with `count`, `result` (upstream's result table) and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
import random
rng = random.Random(9)
sequence = "".join(rng.choice("ACGT") for _ in range(300))
with open("genome.fa", "w") as handle:
    handle.write(">chr1\n" + sequence + "\n")
with open("genome.fa.fai", "w") as handle:
    handle.write("chr1\t300\t6\t300\t301\n")
with open("genome.fa.sizes", "w") as handle:
    handle.write("chr1\t300\n")
from dnakit.predictions import grelu_marginalize_variants

ref = sequence[99]
alt = "A" if ref != "A" else "C"
result = grelu_marginalize_variants(model, [{"chrom": "chr1", "pos": 100, "ref": ref, "alt": alt}], "genome.fa", seq_len=60, n_shuffles=2, seed=0)
print(result["count"], sorted(result["result"]))
```

**Example output:**

```text
1 ['effect_size', 'mean', 'pvalue', 'sd']
```

## 19) `design_sequences_evolve`

- **Purpose:** Greedy directed evolution of each sequence: try single-base substitutions (ism) or pattern insertions (pattern) and keep candidates that raise one task's prediction.
- **Method:** Calls `grelu.design.evolve`. The model must list task names in `data_params['tasks']['name']`, which gReLU sets when a model is trained or loaded. The wrapper rejects a model without them.
- **API:** `design_sequences_evolve(sequences[required], model[required], task_index[optional], method[optional], patterns[optional], max_iter[optional], positions[optional], batch_size[optional], return_seqs[optional], verbose[optional])`
- **Input:** A list of DNA sequences and a LightningModel; positions restrict which bases may change.
- **Output:** A dictionary with `count`, `records` (one per iteration and candidate) and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import design_sequences_evolve

sequence = SEQ[:30]
result = design_sequences_evolve([sequence], model, max_iter=2, positions=[0, 1, 2])
changed = sorted({index for record in result["records"] for index, (a, b) in enumerate(zip(sequence, record["seq"])) if a != b})
print(result["count"], changed)
```

**Example output:**

```text
19 [0, 1, 2]
```

## 20) `design_sequence_ledidi`

- **Purpose:** Design edits of one sequence with Ledidi to raise one task's prediction. Only the listed positions may change.
- **Method:** Calls `ledidi.Ledidi` directly with the loss of `grelu.design.ledidi` (the negative mean prediction). Upstream passes a 0-dimensional target that ledidi 2.2.0 rejects with an IndexError on `y_bar.shape[0]`, so the wrapper passes a one-element target. `edit_penalty` is the L1 edit penalty (upstream default 0.1). Global torch random state is restored afterwards.
- **API:** `design_sequence_ledidi(sequence[required], model[required], task_index[optional], max_iter[optional], batch_size[optional], positions[optional], lr[optional], edit_penalty[optional], seed[optional], verbose[optional])`
- **Input:** One DNA sequence and a LightningModel; positions are the bases allowed to change.
- **Output:** A dictionary with `designs`, `edit_counts`, `count` and `provenance`.

**Example code:**

```python
import torch
from grelu.lightning import LightningModel

torch.manual_seed(0)
model = LightningModel(
    model_params={
        "model_type": "ConvModel",
        "n_tasks": 2,
        "stem_channels": 8,
        "stem_kernel_size": 5,
        "n_conv": 1,
        "channel_init": 8,
        "kernel_size": 3,
        "final_pool_func": "avg",
    },
    train_params={"task": "regression", "loss": "mse"},
)
model.eval()
model.data_params["tasks"] = {"name": ["t0", "t1"]}
model.data_params["train"] = {"seq_len": 60}
SEQ = "ACGTTAGCAGGTCATTAGCCATGAGGCTAACGTTCAAGTACCGTAGGATCCAGTTGACAT"
from dnakit.predictions import design_sequence_ledidi

result = design_sequence_ledidi(SEQ, model, max_iter=30, positions=[0, 1, 2, 3, 4], edit_penalty=0.0)
print(result["count"], sorted(set(result["edit_counts"])))
changed = sorted({index for design in result["designs"] for index, (a, b) in enumerate(zip(SEQ, design)) if a != b})
print(changed)
```

**Example output:**

```text
16 [3]
[1, 3, 4]
```
