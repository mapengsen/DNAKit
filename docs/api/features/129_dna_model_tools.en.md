# DNA model toolkits: AlphaGenome, Enformer, Selene, tangermeme and genomic-benchmarks

This group of functions is a set of adapters over five optional packages for sequence-based deep-learning models: [AlphaGenome](https://github.com/google-deepmind/alphagenome) (Google DeepMind's SDK for the hosted AlphaGenome model and its plotting helpers), [enformer-pytorch](https://github.com/lucidrains/enformer-pytorch) (the Enformer model), [Selene](https://github.com/FunctionLab/selene) (selene-sdk, a framework for training and applying sequence models), [tangermeme](https://github.com/jmschrei/tangermeme) (attribution and motif plots for sequence models) and [genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks) (benchmark datasets of genomic sequences). They connect to the AlphaGenome service, select and pad predicted tracks, draw tracks, logos, contact maps, transcripts and ROC curves, extract Enformer embeddings, run Selene pipelines and model predictions, compute attributions, and load benchmark datasets. The AlphaGenome, Enformer, Selene and genomic-benchmarks functions are in `dnakit.predictions`, and the tangermeme functions are in `dnakit.regulation`. The examples import them with `from dnakit.predictions import ...` and `from dnakit.regulation import ...`.

The dependencies are optional and are imported only when a function is called. If a package is missing, the function raises `BackendUnavailableError` with an install command. The package metadata gives these licences for the Python code: Apache-2.0 for alphagenome and genomic-benchmarks, MIT for enformer-pytorch and tangermeme, BSD for selene-sdk; none of them is GPL. These licences cover the packages only. The hosted AlphaGenome service has its own terms of use and needs an API key, and the pretrained Enformer weights and the benchmark datasets are downloaded from third-party servers under their own terms.

`predict_selene_sequences()` unpickles a whole model file, and the Selene pipeline functions import the model Python file and the classes named in the configuration. Use only files that you trust.

## Example data

The Selene pipeline examples (sections 15, 16 and 17) share three declared files: a tiny model definition, a YAML configuration that scores the sequences of a FASTA file with the `analyze` operation, and the two sequences. All other examples build their inputs in the code.

**Example file:** `model.py`

```text
import torch
import torch.nn as nn


class TinyNet(nn.Module):
    def __init__(self, sequence_length, n_targets):
        super().__init__()
        self.linear = nn.Linear(4 * sequence_length, n_targets)

    def forward(self, x):
        return torch.sigmoid(self.linear(x.flatten(1)))


def criterion():
    return nn.BCELoss()


def get_optimizer(lr):
    return torch.optim.SGD, {"lr": lr}
```

**Example file:** `seqs.fa`

```text
>s1
ACGTACGT
>s2
TTTTGGGG
```

**Example file:** `config.yml`

```text
ops: [analyze]
model: {
    path: model.py,
    class: TinyNet,
    class_args: {sequence_length: 8, n_targets: 2}
}
analyze_sequences: !obj:selene_sdk.predict.AnalyzeSequences {
    trained_model_path: weights.pth,
    sequence_length: 8,
    features: [f1, f2],
    batch_size: 2,
    use_cuda: False
}
prediction: {
    input_path: seqs.fa,
    output_dir: preds
}
random_seed: 123
```

## 1) `dnakit.predictions.create_alphagenome_atlas_client()`

- **Purpose:** Connect to the AlphaGenome Atlas service with an API key and describe the connection.
- **Method:** Call `alphagenome.atlas.atlas.create()` with the key, the timeout and the address. alphagenome opens a secure gRPC channel and waits until it is ready. The function returns the client class name, the address and the timeout. The client object itself is not returned (it keeps the channel open) and the key is not part of the result.
- **API:** `create_alphagenome_atlas_client(api_key: 'str', *, timeout: 'float | None' = None, address: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package alphagenome (extra `dnakit[bio-alphagenome]`), internet access and an AlphaGenome API key.
- **Input:** `api_key` is a non-empty string. `timeout` is a positive number of seconds, or `None` to wait without a limit. `address` is an optional non-empty server address; the default is `dns:///gdmscience.googleapis.com:443`.
- **Output:** `service` (`alphagenome.atlas`), `source` (the address used) and `result` with `client_type`, `address` and `timeout_seconds`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import create_alphagenome_atlas_client

try:
    create_alphagenome_atlas_client("   ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **Example output:**

```text
ConfigurationError api_key must be a non-empty string.
```

The example does not open a connection: the key is checked before the network is used. Usage with the dependency installed:

```python
from dnakit.predictions import create_alphagenome_atlas_client

info = create_alphagenome_atlas_client("YOUR_API_KEY", timeout=30.0)
print(info["result"]["client_type"])
```

## 2) `dnakit.predictions.create_alphagenome_dna_client()`

- **Purpose:** Connect to the hosted AlphaGenome DNA model service with an API key and describe the client.
- **Method:** Map `model_version` to alphagenome's `ModelVersion` enum and call `alphagenome.models.dna_client.create()`. alphagenome opens a secure gRPC channel and waits until it is ready. The client object is not returned and the key is not part of the result.
- **API:** `create_alphagenome_dna_client(api_key: 'str', *, model_version: 'str | None' = None, timeout: 'float | None' = None, address: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package alphagenome (extra `dnakit[bio-alphagenome]`), internet access and an AlphaGenome API key.
- **Input:** `api_key` is a non-empty string. `model_version` is `None` (alphagenome's default, reported as `all_folds`) or one of `all_folds`, `fold_0`, `fold_1`, `fold_2`, `fold_3`. `timeout` and `address` are as in `create_alphagenome_atlas_client()`.
- **Output:** `service` (`alphagenome.models.dna_client`), `source` (the address used) and `result` with `client_type`, `address`, `model_version`, `model_version_source` (`requested` or `upstream default`) and `timeout_seconds`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import create_alphagenome_dna_client

try:
    create_alphagenome_dna_client("key", model_version="fold_9", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **Example output:**

```text
ConfigurationError model_version must be one of: all_folds, fold_0, fold_1, fold_2, fold_3.
```

The example does not open a connection: the model version is checked before the network is used. Usage with the dependency installed:

```python
from dnakit.predictions import create_alphagenome_dna_client

info = create_alphagenome_dna_client("YOUR_API_KEY", model_version="fold_0", timeout=30.0)
print(info["result"]["model_version"])
```

## 3) `dnakit.predictions.get_alphagenome_junctions_to_plot()`

- **Purpose:** Select the splice junctions of one predicted track and strand with alphagenome's rule.
- **Method:** Build an alphagenome `JunctionData` object from the junction list, the value matrix and the track names, then call `alphagenome.data.junction_data.get_junctions_to_plot()`. It keeps the junctions on the requested strand whose value for the named track is at or above the threshold, and sets `k` to that value rounded to two decimals. The function only selects data; `plot_alphagenome_sashimi()` draws.
- **API:** `get_alphagenome_junctions_to_plot(junctions: 'list[dict[str, Any]]', values: 'list[list[float]]', track_names: 'list[str]', name: 'str', strand: 'str', *, k_threshold: 'float | None' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`) and numpy (extra `dnakit[bio-analysis]`). Works offline.
- **Input:** `junctions` is a list of objects with `chromosome`, `start`, `end` (greater than `start`) and `strand` (`+` or `-`). `values` has one row per junction and one column per track. `track_names` has unique names, one per column. `name` is the track to select and `strand` is `+` or `-`. `k_threshold` is a number, or `None` for alphagenome's rule of 5 % of the maximum value.
- **Output:** `track`, `strand`, `k_threshold`, `count`, and `junctions` (objects with `chromosome`, `start`, `end`, `strand`, `k`).
- **Example code:**

```python
from dnakit.predictions import get_alphagenome_junctions_to_plot

junctions = [
    {"chromosome": "chr1", "start": 100, "end": 200, "strand": "+"},
    {"chromosome": "chr1", "start": 300, "end": 400, "strand": "+"},
    {"chromosome": "chr1", "start": 500, "end": 600, "strand": "-"},
]
values = [[0.5, 2.0], [3.0, 0.1], [1.0, 4.0]]
result = get_alphagenome_junctions_to_plot(
    junctions, values, ["liver", "brain"], "liver", "+", k_threshold=0.4, verbose=False
)
print(result["count"])
for junction in result["junctions"]:
    print(junction["start"], junction["end"], junction["k"])
```

- **Example output:**

```text
2
100 200 0.5
300 400 3.0
```

## 4) `dnakit.predictions.pad_alphagenome_track()`

- **Purpose:** Centre a positions-by-channels track inside a longer window.
- **Method:** Call alphagenome's `pad_track()`, which puts equal padding on both sides when the length difference is even and the extra position at the end when it is odd. The function reports the index of the first original row as `offset`.
- **API:** `pad_alphagenome_track(track: 'list[list[float]]', new_len: 'int', *, value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`) and numpy (extra `dnakit[bio-analysis]`). Works offline.
- **Input:** `track` is a non-empty rectangular list of rows of finite numbers. `new_len` is an integer at least as large as the number of rows. `value` is the finite fill number.
- **Output:** `length`, `new_length`, `offset`, `value` and `padded` (the padded rows).
- **Example code:**

```python
from dnakit.predictions import pad_alphagenome_track

result = pad_alphagenome_track([[1.0, 2.0], [3.0, 4.0]], 5, value=-1.0, verbose=False)
print(result["length"], result["new_length"], result["offset"])
print(result["padded"])
```

- **Example output:**

```text
2 5 1
[[-1.0, -1.0], [1.0, 2.0], [3.0, 4.0], [-1.0, -1.0], [-1.0, -1.0]]
```

## 5) `dnakit.predictions.plot_alphagenome_components()`

- **Purpose:** Draw several model outputs as stacked components over one genomic interval and save them as a PNG.
- **Method:** Convert each component object to an alphagenome component (`Tracks`, `SeqLogo` or `EmptyComponent`) and call `alphagenome.visualization.plot_components.plot()` for the interval. The figure is saved at 150 dpi and closed.
- **API:** `plot_alphagenome_components(components: 'list[dict[str, Any]]', output_path: 'str', *, chromosome: 'str', start: 'int', end: 'int', fig_width: 'float' = 20.0, fig_height_scale: 'float' = 1.0, title: 'str | None' = None, despine: 'bool' = True, despine_keep_bottom: 'bool' = False, hspace: 'float' = 0.3, xlabel: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `components` is a non-empty list of objects. A `tracks` object has `names`, `values` (one row per base, one column per track), and optional `strands` (`+`, `-` or `.`), `track_height` and `filled`. A `seq_logo` object has `scores` (one row per base, one column per letter) and optional `alphabet`, `fig_height` and `ylabel`. An `empty` object has an optional `fig_height`. The number of rows must equal `end - start`. Transcript and contact-map components are not supported here. `output_path` must end with `.png`; missing folders are created and an existing file is overwritten.
- **Output:** `path` (the PNG) and `data` with the `interval`, the described `components` and the `title`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_components

components = [
    {"type": "tracks", "names": ["a", "b"], "strands": ["+", "-"],
     "values": [[1, 2], [2, 1], [3, 0], [0, 1]]},
    {"type": "seq_logo", "scores": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]},
    {"type": "empty"},
]
result = plot_alphagenome_components(
    components, "components.png", chromosome="chr1", start=100, end=104, verbose=False
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print([item["type"] for item in result["data"]["components"]], result["data"]["interval"])
```

- **Example output:**

```text
True
['tracks', 'seq_logo', 'empty'] {'chromosome': 'chr1', 'start': 100, 'end': 104}
```

## 6) `dnakit.predictions.plot_alphagenome_contact_map()`

- **Purpose:** Draw a square contact map on a logarithmic colour scale and save it as a PNG.
- **Method:** Call alphagenome's `plot_contact_map()` on the matrix. The colour limits are the given `vmin` and `vmax`, or the smallest positive value and the largest value of the matrix; the scale that was used is reported.
- **API:** `plot_alphagenome_contact_map(contact_map: 'list[list[float]]', output_path: 'str', *, vmin: 'float | None' = None, vmax: 'float | None' = None, square: 'bool' = True, cbar_shrink: 'float' = 0.4, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `contact_map` is a square matrix of finite numbers with at least one positive value (or an explicit `vmin`). `vmin` must be positive and below `vmax`. `square` is a boolean and `cbar_shrink` is a positive number. `output_path` must end with `.png`.
- **Output:** `path` and `data` with the matrix, `vmin`, `vmax`, `square` and `cbar_shrink`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_contact_map

matrix = [[1.0, 0.5, 0.1], [0.5, 1.0, 0.4], [0.1, 0.4, 1.0]]
result = plot_alphagenome_contact_map(matrix, "contacts.png", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["vmin"], result["data"]["vmax"])
```

- **Example output:**

```text
True
0.1 1.0
```

## 7) `dnakit.predictions.plot_alphagenome_sashimi()`

- **Purpose:** Draw splice junctions as a Sashimi plot and save it as a PNG.
- **Method:** Build alphagenome `Junction` objects with their read counts and call `alphagenome.visualization.plot.sashimi_plot()` with a random generator seeded by `seed`, so that the arc jitter and the PNG are reproducible. The result also lists the junctions at or above `filter_threshold`, which are the ones drawn.
- **API:** `plot_alphagenome_sashimi(junctions: 'list[dict[str, Any]]', output_path: 'str', *, interval: 'dict[str, Any] | None' = None, filter_threshold: 'float' = 0.01, annotate_counts: 'bool' = True, seed: 'int | None' = 0, color: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `junctions` is a list of objects with `chromosome`, `start`, `end`, `strand` (`+` or `-`) and `k` (a non-negative read count); at least one `k` must be above zero. `interval` is an optional object with `chromosome`, `start` and `end` that limits the count labels. `filter_threshold` is a non-negative number. `seed` is an integer or `None`. `color` is an optional colour name.
- **Output:** `path` and `data` with `junctions`, `drawn`, `drawn_count`, `total_k`, `filter_threshold`, `annotate_counts`, `seed` and `ylim`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_sashimi

junctions = [
    {"chromosome": "chr1", "start": 100, "end": 200, "strand": "+", "k": 5},
    {"chromosome": "chr1", "start": 300, "end": 400, "strand": "+", "k": 0.001},
    {"chromosome": "chr1", "start": 150, "end": 380, "strand": "+", "k": 2},
]
result = plot_alphagenome_sashimi(
    junctions, "sashimi.png",
    interval={"chromosome": "chr1", "start": 50, "end": 450}, verbose=False,
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["drawn_count"], result["data"]["total_k"], result["data"]["seed"])
```

- **Example output:**

```text
True
2 7.001 0
```

## 8) `dnakit.predictions.plot_alphagenome_seqlogo()`

- **Purpose:** Draw a sequence logo from per-position letter heights and save it as a PNG.
- **Method:** Call alphagenome's `seqlogo()` on a matplotlib axis. Positive heights stack above the axis and negative heights below it.
- **API:** `plot_alphagenome_seqlogo(letter_heights: 'list[list[float]]', output_path: 'str', *, alphabet: 'str' = 'ACGT', one_based: 'bool' = True, start: 'int' = 0, letter_colors: 'str' = 'default', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `letter_heights` has one row per position and one column per letter of `alphabet`. `one_based` and `start` set the numbering of the x axis. `letter_colors` is `default` or `factorbook`. `output_path` must end with `.png`.
- **Output:** `path` and `data` with `letter_heights`, `alphabet`, `one_based`, `start`, `letter_colors` and `positions`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_seqlogo

heights = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, -1]]
result = plot_alphagenome_seqlogo(heights, "seqlogo.png", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["positions"], result["data"]["alphabet"], result["data"]["letter_colors"])
```

- **Example output:**

```text
True
3 ACGT default
```

## 9) `dnakit.predictions.plot_alphagenome_track()`

- **Purpose:** Draw one track and save it as a PNG, choosing the drawing from the shape of the track.
- **Method:** Call alphagenome's `plot_track()`. A 1D list, or a 2D list with one column, is drawn as a line (filled when `filled` is true); a 2D list with four columns is drawn as a sequence logo; a 2D list with two columns is drawn as two lines (positive and negative strand). The chosen drawing is reported as `rendering`.
- **API:** `plot_alphagenome_track(track: 'list[Any]', output_path: 'str', *, x: 'list[float] | None' = None, legend: 'bool' = False, ylim: 'list[float] | None' = None, color: 'str | list[str] | None' = None, filled: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `track` is a 1D list or a 2D list with 1, 2 or 4 columns. `x` is an optional list with one position per row. `ylim` is a `[low, high]` pair with `high` above `low`. `color` is a colour or a list of colours. `legend` and `filled` are booleans.
- **Output:** `path` and `data` with `track`, `length`, `rendering` (`line`, `sequence_logo` or `two_strand_lines`), `x`, `ylim` and `filled`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_track

result = plot_alphagenome_track([0.1, 0.5, 0.9, 0.4], "track.png", filled=True, verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["rendering"], result["data"]["length"], result["data"]["filled"])
```

- **Example output:**

```text
True
line 4 True
```

## 10) `dnakit.predictions.plot_alphagenome_tracks()`

- **Purpose:** Draw several tracks as stacked panels in one PNG.
- **Method:** Call alphagenome's `plot_tracks()` on the named tracks. Each track becomes one panel with its own y scale (`ylim="auto"`), one shared scale (`"same"`) or the given limits.
- **API:** `plot_alphagenome_tracks(tracks: 'dict[str, list[Any]]', output_path: 'str', *, x: 'list[float] | None' = None, title: 'str | None' = None, legend: 'bool' = False, fig_width: 'float' = 20.0, fig_track_height: 'float | dict[str, float]' = 1.5, ylim: 'str | list[float] | None' = 'auto', yticks_min_max_only: 'bool' = False, ylab: 'bool' = True, color: 'str | dict[str, str] | None' = None, horizontal_ylab: 'bool' = True, filled_tracks: 'list[str] | None' = None, despine: 'bool' = True, despine_keep_bottom: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages alphagenome (extra `dnakit[bio-alphagenome]`), numpy (extra `dnakit[bio-analysis]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `tracks` maps each track name to a 1D list or a 2D list; all tracks must have the same length. `ylim` is `auto`, `same`, a `[low, high]` pair or `None`. `fig_track_height` is one number or a dict with one height per track. `color` is one colour or a dict that names only existing tracks. `filled_tracks` is a list of track names. The other arguments are booleans, numbers and a title.
- **Output:** `path` and `data` with `tracks`, `track_names`, `length`, `x`, `ylim`, `title` and `filled_tracks`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_tracks

tracks = {"liver": [0, 1, 2, 1], "brain": [2, 1, 0, 1]}
result = plot_alphagenome_tracks(
    tracks, "tracks.png", title="Two tracks", filled_tracks=["brain"], verbose=False
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["track_names"], result["data"]["length"], result["data"]["filled_tracks"])
```

- **Example output:**

```text
True
['liver', 'brain'] 4 ['brain']
```

## 11) `dnakit.predictions.plot_alphagenome_transcripts()`

- **Purpose:** Draw transcript models over one genomic interval and save them as a PNG.
- **Method:** Build alphagenome `Transcript` objects from the exon and coding-region intervals and call `alphagenome.visualization.plot_transcripts.plot_transcripts()` on a matplotlib axis. The figure height grows with the number of transcripts.
- **API:** `plot_alphagenome_transcripts(transcripts: 'list[dict[str, Any]]', output_path: 'str', *, chromosome: 'str', start: 'int', end: 'int', zero_origin: 'bool' = False, label_name: 'str | None' = None, plot_labels_once: 'bool' = False, fig_width: 'float' = 12.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package alphagenome (extra `dnakit[bio-alphagenome]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `transcripts` is a non-empty list of objects with `transcript_id`, `strand` (`+` or `-`), `exons` (a list of `[start, end]` pairs), and optional `cds` (pairs, or null for a non-coding transcript) and `info` (an object of simple values). `chromosome`, `start` and `end` give the interval. `label_name` names the `info` key shown next to each transcript and every transcript must have it. `zero_origin` and `plot_labels_once` are booleans.
- **Output:** `path` and `data` with `interval`, the described `transcripts`, `zero_origin`, `label_name` and `plot_labels_once`.
- **Example code:**

```python
from dnakit.predictions import plot_alphagenome_transcripts

transcripts = [
    {"transcript_id": "T1", "strand": "+", "exons": [[100, 200], [300, 400]],
     "cds": [[150, 200], [300, 350]], "info": {"gene_name": "GENE1"}},
    {"transcript_id": "T2", "strand": "-", "exons": [[120, 180]], "info": {"gene_name": "GENE2"}},
]
result = plot_alphagenome_transcripts(
    transcripts, "transcripts.png", chromosome="chr1", start=50, end=450,
    label_name="gene_name", verbose=False,
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print([item["transcript_id"] for item in result["data"]["transcripts"]], result["data"]["label_name"])
```

- **Example output:**

```text
True
['T1', 'T2'] gene_name
```

## 12) `dnakit.predictions.download_pretrained_enformer()`

- **Purpose:** Download pretrained Enformer weights from the Hugging Face hub and save them to local files.
- **Method:** Call `enformer_pytorch.modeling_enformer.from_pretrained()` (passing `use_tf_gamma` only when it is not `None`), save `model.state_dict()` with `torch.save()` and write the model configuration as sorted JSON next to it.
- **API:** `download_pretrained_enformer(name: 'str', output_state_dict: 'str', *, use_tf_gamma: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages enformer-pytorch and torch (extra `dnakit[neural-enformer]`) and internet access the first time (later calls use the Hugging Face cache). The install hint printed by the function names an extra called `enformer`; the extra in `pyproject.toml` is `neural-enformer`.
- **Input:** `name` is a non-empty Hugging Face model id such as `EleutherAI/enformer-official-rough`. `output_state_dict` is a new file path in an existing folder; `<output_state_dict>.config.json` must not exist either. `use_tf_gamma` is `None`, true or false.
- **Output:** `status`, `backend`, `model`, `state_dict` and `config` (the two file paths), `parameter_count` and `output_heads` (sorted head names).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import download_pretrained_enformer

try:
    download_pretrained_enformer("", "enformer.pt", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **Example output:**

```text
ConfigurationError name must be a Hugging Face model id such as EleutherAI/enformer-official-rough.
```

The example does not download anything: the model name is checked before the network is used. Usage with the dependency installed:

```python
from dnakit.predictions import download_pretrained_enformer

result = download_pretrained_enformer(
    "EleutherAI/enformer-official-rough", "enformer_official.pt"
)
print(result["parameter_count"], result["output_heads"])
```

## 13) `dnakit.predictions.enformer_embeddings()`

- **Purpose:** Compute Enformer embeddings for DNA sequences from a local state dict and write them to a `.npy` file.
- **Method:** Build `Enformer(EnformerConfig(**config))`, load the state dict (`torch.load(..., weights_only=True)` on the CPU), put the model in eval mode, encode the sequences with `enformer_pytorch.data.str_to_seq_indices()` and call `enformer_pytorch.finetune.get_enformer_embeddings()`. The array is saved with `numpy.save()`.
- **API:** `enformer_embeddings(state_dict: 'str', sequences: 'Sequence[str]', output_npy: 'str', *, config: 'Mapping[str, Any]', freeze: 'bool' = False, train_layernorms_only: 'bool' = False, train_last_n_layers_only: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages enformer-pytorch and torch (extra `dnakit[neural-enformer]`) and numpy (extra `dnakit[bio-analysis]`). Works offline once the weights are on disk. The install hint printed by the function names an extra called `enformer`; the extra in `pyproject.toml` is `neural-enformer`.
- **Input:** `state_dict` is an existing weights file and `config` the non-empty mapping of `EnformerConfig` arguments the weights were trained with (the JSON written by `download_pretrained_enformer()` can be loaded and passed). `sequences` is a non-empty list of strings of A, C, G, T and N with one common length. `output_npy` is a new file path in an existing folder. `freeze` cannot be combined with `train_layernorms_only` or `train_last_n_layers_only`, which must not be negative.
- **Output:** `status`, `backend`, `output` (the file path), `shape`, `sequence_count`, `sequence_length`, and the three training options.
- **Example code:**

```python
import numpy as np
import torch
from enformer_pytorch import Enformer, EnformerConfig

from dnakit.predictions import enformer_embeddings

config = {"dim": 48, "depth": 1, "heads": 4, "target_length": -1, "num_downsamples": 3,
          "dim_divisible_by": 4, "output_heads": {"human": 4}}
torch.save(Enformer(EnformerConfig(**config)).state_dict(), "tiny_enformer.pt")
result = enformer_embeddings(
    "tiny_enformer.pt", ["ACGT" * 32, "ACGN" * 32], "embeddings.npy", config=config, verbose=False
)
print(result["shape"], result["sequence_count"], result["sequence_length"])
print(np.load("embeddings.npy").shape)
```

- **Example output:**

```text
[2, 16, 96] 2 128
(2, 16, 96)
```

The example uses a tiny randomly initialised Enformer, so only the shapes are printed. With a real checkpoint, pass the configuration saved by `download_pretrained_enformer()`.

## 14) `dnakit.predictions.predict_selene_sequences()`

- **Purpose:** Run a whole trained Selene-style torch model on DNA sequences on the CPU.
- **Method:** Load the file with `torch.load(..., weights_only=False)`, require a `torch.nn.Module`, put it in eval mode, one-hot encode the sequences with selene's `Genome.sequence_to_encoding()` (N becomes a uniform 0.25 row) and call `selene_sdk.predict.model_predict.predict()` on a (batch, length, 4) array without CUDA; selene transposes it to (batch, 4, length) before it calls the model.
- **API:** `predict_selene_sequences(model_path: 'str', sequences: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`, which installs torch, numpy and matplotlib as dependencies). Works offline.
- **Input:** `model_path` is an existing file written with `torch.save(model)` for a whole module; a state dict is rejected. The file is unpickled, so use only trusted files. `sequences` is one DNA string or a list of at most 256 strings of one common length of at most 10,000 bases.
- **Output:** `status`, `operation`, `model_path`, `sequence_count`, `length`, `shape` (batch, features), `predictions` (nested lists) and `provenance`.
- **Example code:**

```python
import torch

from dnakit.predictions import predict_selene_sequences


class CountA(torch.nn.Module):
    def forward(self, x):
        return x[:, 0, :].sum(dim=1, keepdim=True)


torch.save(CountA(), "count_a.pt")
result = predict_selene_sequences("count_a.pt", ["AACG", "TTTT", "ANNA"], verbose=False)
print(result["shape"], result["predictions"])
```

- **Example output:**

```text
[3, 1] [[2.0], [0.0], [2.5]]
```

The model receives a (batch, 4, length) tensor and counts the A positions of each sequence; the `N` bases of the third sequence count as 0.25 each, so `ANNA` gives 2.5.

## 15) `dnakit.predictions.run_selene_config_file()`

- **Purpose:** Run the Selene pipeline described by a YAML configuration file and write its outputs under an output folder.
- **Method:** Load the YAML with selene's `load_path()` (without instantiating the objects), replace its `output_dir` key by the given folder and call `selene_sdk.utils.config_utils.parse_configs_and_run()`. This is what selene's `load_and_parse_configs_and_run()` does after reading the file, except that the copy of the configuration in the run folder is a re-dumped `configuration.yaml` instead of a copy of the original text. selene prints progress messages to standard output.
- **API:** `run_selene_config_file(config_path: 'str', output_dir: 'str', *, create_subdirectory: 'bool' = True, lr: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`). Works offline.
- **Input:** `config_path` is an existing YAML file with a mapping at the top level, an `ops` list of up to 16 names from `train`, `evaluate`, `analyze`, and `model: {path: ...}` naming an existing model file. `output_dir` is a folder (created when missing). `create_subdirectory` puts the run into a new timestamped subfolder. `lr` is a positive number, needed only when the pipeline trains.
- **Output:** `status`, `operation`, `source` (the configuration path), `operations`, `output_dir`, `run_directory`, `files` (the sorted names of the files in the run folder) and `provenance`.
- **Example code:**

```python
import contextlib
import io
import os

import torch

import model
from dnakit.predictions import run_selene_config_file

net = model.TinyNet(8, 2)
with torch.no_grad():
    net.linear.weight.zero_()
    net.linear.bias.copy_(torch.tensor([0.0, 1.0986123]))
torch.save(net.state_dict(), "weights.pth")

with contextlib.redirect_stdout(io.StringIO()):
    result = run_selene_config_file("config.yml", "run_file", verbose=False)
print(result["operations"], result["files"])
with open(os.path.join("preds", "seqs_predictions.tsv")) as handle:
    print(handle.read())
```

- **Example output:**

```text
['analyze'] ['configuration.yaml', 'model.py']
index	name	f1	f2
0	s1	5.00e-01	7.50e-01
1	s2	5.00e-01	7.50e-01
```

The weights are fixed (zero weights and a bias of 0 and ln 3), so the two features of every sequence are sigmoid(0) = 0.5 and sigmoid(ln 3) = 0.75.

## 16) `dnakit.predictions.run_selene_configs()`

- **Purpose:** Run a Selene pipeline from a parsed configuration dictionary and write its outputs under an output folder.
- **Method:** Copy the dictionary, set its `output_dir` and call `selene_sdk.utils.config_utils.parse_configs_and_run()`. The dictionary has the structure that selene's `load_path()` returns, so that objects given with `!obj:` tags are still unresolved when the pipeline starts.
- **API:** `run_selene_configs(configs: 'Mapping[str, Any]', output_dir: 'str', *, create_subdirectory: 'bool' = True, lr: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`). Works offline.
- **Input:** `configs` is a mapping with string keys, an `ops` list and `model: {path: ...}` as in `run_selene_config_file()`; the caller's object is not changed. `output_dir`, `create_subdirectory` and `lr` are as in `run_selene_config_file()`. When `ops` contains `train`, either `configs` has `lr` or the `lr` argument is given.
- **Output:** The same keys as `run_selene_config_file()`, with `source` set to `configs dict`.
- **Example code:**

```python
import contextlib
import io
import os

import torch
from selene_sdk.utils.config_utils import load_path

import model
from dnakit.predictions import run_selene_configs

net = model.TinyNet(8, 2)
with torch.no_grad():
    net.linear.weight.zero_()
    net.linear.bias.copy_(torch.tensor([0.0, 1.0986123]))
torch.save(net.state_dict(), "weights.pth")

configs = load_path("config.yml", instantiate=False)
configs["prediction"]["output_dir"] = "preds_dict"
with contextlib.redirect_stdout(io.StringIO()):
    result = run_selene_configs(configs, "run_dict", verbose=False)
print(result["source"], result["operations"], result["files"])
print(sorted(os.listdir("preds_dict")))
```

- **Example output:**

```text
configs dict ['analyze'] ['configuration.yaml', 'model.py']
['seqs_predictions.tsv']
```

## 17) `dnakit.predictions.run_selene_operations()`

- **Purpose:** Run Selene's configuration-defined operations (`train`, `evaluate`, `analyze`) in the given order.
- **Method:** Call selene's `execute()` with the operation names, the configuration and the output folder. Unlike `run_selene_configs()` it does not create a run folder, copy the model or write a configuration file; each operation reads the configuration sections it needs.
- **API:** `run_selene_operations(operations: 'Sequence[str]', configs: 'Mapping[str, Any]', output_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`). Works offline.
- **Input:** `operations` is a list of up to 16 names from `train`, `evaluate`, `analyze`. `configs` is a mapping with the sections those operations read (`model` for all of them, `analyze_sequences` and `prediction` for the example below, and `lr` when `train` is listed). `output_dir` is a folder (created when missing).
- **Output:** `status`, `operation`, `operations`, `output_dir` and `provenance`. The files written by the operations are reported by the operations themselves, in the folders named in the configuration.
- **Example code:**

```python
import contextlib
import io
import os

import torch
from selene_sdk.utils.config_utils import load_path

import model
from dnakit.predictions import run_selene_operations

net = model.TinyNet(8, 2)
with torch.no_grad():
    net.linear.weight.zero_()
    net.linear.bias.copy_(torch.tensor([0.0, 1.0986123]))
torch.save(net.state_dict(), "weights.pth")

configs = load_path("config.yml", instantiate=False)
configs["prediction"]["output_dir"] = "preds_ops"
with contextlib.redirect_stdout(io.StringIO()):
    result = run_selene_operations(["analyze"], configs, "run_ops", verbose=False)
print(result["status"], result["operations"])
print(sorted(os.listdir("preds_ops")))
```

- **Example output:**

```text
completed ['analyze']
['seqs_predictions.tsv']
```

## 18) `dnakit.predictions.selene_precision_recall_curves()`

- **Purpose:** Write one precision-recall curve per model feature to an SVG file.
- **Method:** Select the features by selene's rule, compute the average precision of each with scikit-learn, and call `selene_sdk.utils.performance_metrics.visualize_precision_recall_curves()`, which writes `precision_recall_curves.svg` into the output folder. The matplotlib style is applied only for this call and open figures are closed.
- **API:** `selene_precision_recall_curves(prediction: 'Sequence[Sequence[float]]', target: 'Sequence[Sequence[int]]', output_dir: 'str', *, report_gt_feature_n_positives: 'int' = 50, style: 'str' = 'seaborn-v0_8-colorblind', fig_title: 'str' = 'Feature precision-recall curves', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`, which installs scikit-learn and matplotlib as dependencies). Works offline.
- **Input:** `prediction` and `target` are (samples, features) tables of the same shape with at least two samples; `target` holds only 0 and 1. A feature is drawn only when both labels occur and its number of positives is larger than `report_gt_feature_n_positives`. `output_dir` is a folder (created when missing) in which `precision_recall_curves.svg` must not exist yet. `style` is a matplotlib style name.
- **Output:** `status`, `operation`, `path` (the SVG), `image_format`, `data` (`sample_count`, `feature_count`, `positives_per_feature`, `plotted_features`, `metric` = `average_precision`, `metric_per_plotted_feature`, `style`) and `provenance`.
- **Example code:**

```python
import os

from dnakit.predictions import selene_precision_recall_curves

prediction = [[0.9, 0.1], [0.4, 0.3], [0.5, 0.6], [0.1, 0.7], [0.7, 0.2], [0.3, 0.9]]
target = [[1, 0], [1, 0], [0, 1], [0, 1], [1, 0], [0, 1]]
result = selene_precision_recall_curves(
    prediction, target, "pr_curves", report_gt_feature_n_positives=1, verbose=False
)
print(os.path.basename(result["path"]), result["data"]["plotted_features"])
print([round(value, 3) for value in result["data"]["metric_per_plotted_feature"]])
```

- **Example output:**

```text
precision_recall_curves.svg [0, 1]
[0.917, 1.0]
```

## 19) `dnakit.predictions.selene_roc_curves()`

- **Purpose:** Write one ROC curve per model feature to an SVG file.
- **Method:** Select the features by selene's rule, compute the ROC AUC of each with scikit-learn, and call `selene_sdk.utils.performance_metrics.visualize_roc_curves()`, which writes `roc_curves.svg` into the output folder. The default style is `seaborn-v0_8-colorblind` because selene's own default style name no longer exists in matplotlib.
- **API:** `selene_roc_curves(prediction: 'Sequence[Sequence[float]]', target: 'Sequence[Sequence[int]]', output_dir: 'str', *, report_gt_feature_n_positives: 'int' = 50, style: 'str' = 'seaborn-v0_8-colorblind', fig_title: 'str' = 'Feature ROC curves', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`, which installs scikit-learn and matplotlib as dependencies). Works offline.
- **Input:** The same as `selene_precision_recall_curves()`; the file is `roc_curves.svg`. With the default `report_gt_feature_n_positives=50`, a feature needs more than 50 positives to be drawn, so small tables draw nothing.
- **Output:** The same keys as `selene_precision_recall_curves()`, with `metric` = `roc_auc`.
- **Example code:**

```python
import os

from dnakit.predictions import selene_roc_curves

prediction = [[0.9, 0.1], [0.4, 0.3], [0.5, 0.6], [0.1, 0.7], [0.7, 0.2], [0.3, 0.9]]
target = [[1, 0], [1, 0], [0, 1], [0, 1], [1, 0], [0, 1]]
result = selene_roc_curves(
    prediction, target, "roc_curves", report_gt_feature_n_positives=1, verbose=False
)
print(os.path.basename(result["path"]), result["data"]["positives_per_feature"])
print(result["data"]["plotted_features"], [round(v, 3) for v in result["data"]["metric_per_plotted_feature"]])
```

- **Example output:**

```text
roc_curves.svg [3, 3]
[0, 1] [0.889, 1.0]
```

## 20) `dnakit.predictions.selene_score_heatmap()`

- **Purpose:** Draw a positions-by-bases in silico mutagenesis score matrix as a heatmap PNG.
- **Method:** Call selene's `selene_sdk.interpret.vis.heatmap()` on a matplotlib axis with selene's `Genome` as the sequence type. The rows of the figure are the bases T, G, C, A from top to bottom, as in selene.
- **API:** `selene_score_heatmap(score_matrix: 'Sequence[Sequence[float]]', output_path: 'str', *, mask: 'Sequence[Sequence[bool]] | None' = None, cmap: 'str' = 'Blues', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`). Works offline.
- **Input:** `score_matrix` has at most 5,000 rows (positions) of four numbers for A, C, G, T. `mask` is an optional grid of the same shape with True for the cells to hide. `cmap` is a matplotlib colour-map name. `output_path` must be a new `.png` path.
- **Output:** `status`, `operation`, `path`, `image_format` and `data` (`position_count`, `cmap`, `row_labels_top_to_bottom`, `score_matrix`, `mask`) and `provenance`.
- **Example code:**

```python
import os

from dnakit.predictions import selene_score_heatmap

scores = [[0.1, 0.2, 0.3, 0.4], [0.4, 0.1, 0.0, 0.2], [0.0, 0.5, 0.1, 0.1]]
mask = [[False, False, False, False], [False, True, False, False], [False, False, False, False]]
result = selene_score_heatmap(scores, "ism_heatmap.png", mask=mask, verbose=False)
print(os.path.basename(result["path"]), open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["position_count"], result["data"]["row_labels_top_to_bottom"])
```

- **Example output:**

```text
ism_heatmap.png True
3 ['T', 'G', 'C', 'A']
```

## 21) `dnakit.predictions.plot_sequence_logo_scores()`

- **Purpose:** Draw a DNA sequence logo from a positions-by-bases score matrix, for example attribution scores, and save it as a PNG.
- **Method:** Reproduce selene's `sequence_logo()` with ordinary matplotlib patches, because selene's own call fails with the installed matplotlib (its glyph path effect uses a removed method). The bases of each position are stacked by value (`order="value"`, negatives below zero) or in alphabet order, and every bar is filled with selene's glyph outline of its base (`selene_sdk.interpret.vis._SVG_PATHS`).
- **API:** `plot_sequence_logo_scores(score_matrix: 'Sequence[Sequence[float]]', image_path: 'str | Path', *, order: 'str' = 'value', width: 'float' = 1.0, color_scheme: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`). Works offline.
- **Input:** `score_matrix` has at most 5,000 rows of four finite numbers for A, C, G, T. `order` is `value` or `alpha`. `width` is a non-negative bar width. `color_scheme` is a list of at least four colours for A, C, G, T; the default is the first four colours of seaborn's Set1 palette. `image_path` must be a new `.png` path in an existing folder.
- **Output:** `status`, `operation`, `image_path`, `image_format`, `order`, `bases` and `scores` (both with one row per stack level and one column per position, in drawing order) and `provenance`.
- **Example code:**

```python
from dnakit.predictions import plot_sequence_logo_scores

scores = [[0.2, 0.1, -0.1, 0.5], [0.0, 0.3, 0.1, 0.1], [-0.2, 0.4, 0.0, 0.2]]
result = plot_sequence_logo_scores(scores, "logo_scores.png", verbose=False)
print(open(result["image_path"], "rb").read(4) == b"\x89PNG", result["order"])
for bases, values in zip(result["bases"], result["scores"]):
    print("".join(bases), values)
```

- **Example output:**

```text
True value
GAA [-0.1, 0.0, -0.2]
CGG [0.1, 0.1, 0.0]
ATT [0.2, 0.1, 0.2]
TCC [0.5, 0.3, 0.4]
```

Each printed row is one stack level and each column one position. In the first position the stack runs from G (-0.1, below zero) through C (0.1) and A (0.2) to T (0.5).

## 22) `dnakit.predictions.plot_variant_diffs_scatter()`

- **Purpose:** Plot each variant's largest absolute difference score as a point and write the plot to a new HTML file.
- **Method:** Order the variants by chromosome (selene's `sort_standard_chrs`) and position, take the maximum over the chosen feature columns, optionally keep the variants at or above a percentile, and draw a plotly scatter written with `plotly.offline.plot()`. Selene's own function fails with plotly 6 and mis-pairs labels and scores when a feature filter is used, so this function writes the layout with plotly 6 title objects and takes rows and columns in label order.
- **API:** `plot_variant_diffs_scatter(scores: 'Sequence[Sequence[float]]', labels: 'Sequence[Sequence[Any]]', features: 'Sequence[str]', output_path: 'str | Path', *, feature_indices: 'Sequence[int] | None' = None, nth_percentile: 'float | None' = None, threshold_line: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package selene-sdk (extra `dnakit[bio-selene]`, which installs plotly and numpy as dependencies). Works offline.
- **Input:** `scores` is a variants-by-features matrix of finite numbers (at most 5,000 rows) and `features` names its columns. `labels` has one `(chrom, pos, name, ref, alt)` entry per row. `feature_indices` keeps only those columns. `nth_percentile` is between 0 and 100. `threshold_line` is a number. `output_path` must be a new `.html` path in an existing folder.
- **Output:** `status`, `operation`, `output_path`, `output_format`, `features` (the names used), `variants` (objects with `chrom`, `pos`, `ref`, `alt`, `max_abs_diff` in plotted order) and `provenance`.
- **Example code:**

```python
import os

from dnakit.predictions import plot_variant_diffs_scatter

scores = [[0.1, 0.5], [0.9, 0.2], [0.3, 0.3]]
labels = [("chr2", 50, "v1", "A", "G"), ("chr1", 200, "v2", "C", "T"), ("chr1", 100, "v3", "G", "A")]
result = plot_variant_diffs_scatter(
    scores, labels, ["f1", "f2"], "variant_scatter.html",
    nth_percentile=50, threshold_line=0.4, verbose=False,
)
print(os.path.basename(result["output_path"]), os.path.getsize("variant_scatter.html") > 1000)
for variant in result["variants"]:
    print(variant["chrom"], variant["pos"], variant["max_abs_diff"])
```

- **Example output:**

```text
variant_scatter.html True
chr1 200 0.9
chr2 50 0.5
```

## 23) `dnakit.regulation.plot_attribution_logo()`

- **Purpose:** Draw a sequence logo of a positions-by-bases attribution matrix and save it as a PNG.
- **Method:** Call tangermeme's `plot_logo()` on the transposed matrix (tangermeme uses a (4, length) layout), with the optional window, annotation table and symmetric y limit. The limit is applied to the axis after drawing, because tangermeme ignores it when there are no annotations.
- **API:** `plot_attribution_logo(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, start: 'int | None' = None, end: 'int | None' = None, ylim: 'float | None' = None, annotations: 'Sequence[Mapping[str, Any]] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages tangermeme and torch (extra `dnakit[bio-tangermeme]`) and matplotlib and pandas. Works offline.
- **Input:** `attributions` has at most 10,000 rows of four finite numbers for A, C, G, T. `start` and `end` select a window (`start` below `end`, both within the sequence). `ylim` is a positive bound; the axis runs from `-ylim` to `ylim`. `annotations` is a list of at most 1,000 records with `name`, `start`, `end` and `score`; tangermeme draws only annotations that start after the window start and end before its end. `output_path` must be a new `.png` path.
- **Output:** `status`, `operation`, `path`, `image_format`, `data` (`position_count`, `start`, `end`, `ylim` as `[-ylim, ylim]` or `None`, `alphabet`, `attributions`, `annotations`) and `provenance`.
- **Example code:**

```python
from dnakit.regulation import plot_attribution_logo

attributions = [
    [0.5, 0.0, 0.1, -0.1], [0.0, 0.7, 0.0, 0.1], [0.1, 0.0, -0.4, 0.0],
    [0.0, 0.0, 0.0, 0.9], [0.2, 0.2, 0.0, 0.0],
]
result = plot_attribution_logo(
    attributions, "attribution_logo.png", ylim=1.0,
    annotations=[{"name": "motif", "start": 1, "end": 3, "score": 0.9}], verbose=False,
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["position_count"], result["data"]["ylim"], result["data"]["alphabet"])
```

- **Example output:**

```text
True
5 [-1.0, 1.0] ['A', 'C', 'G', 'T']
```

## 24) `dnakit.regulation.plot_categorical_attribution_scatter()`

- **Purpose:** Plot the summed attribution at every position, coloured by the dominant base, and save it as a PNG.
- **Method:** Call tangermeme's `plot_categorical_scatter()`. The y value is the sum over the four bases and the colour is that of the base with the largest absolute value (A green, C blue, G orange, T red). The colours are fixed because tangermeme's `colors` argument fails. The function also reports the dominant base and the sum of every position. It is meant for sequences that are too long for a logo.
- **API:** `plot_categorical_attribution_scatter(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, size: 'float' = 20.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages tangermeme and torch (extra `dnakit[bio-tangermeme]`) and matplotlib. Works offline.
- **Input:** `attributions` has at most 10,000 rows of four finite numbers for A, C, G, T. `size` is a positive marker size. `output_path` must be a new `.png` path.
- **Output:** `status`, `operation`, `path`, `image_format`, `data` (`position_count`, `alphabet`, `dominant_base_index_by_position` with 0 to 3 for A, C, G, T, `summed_attribution_by_position`, `attributions`) and `provenance`.
- **Example code:**

```python
from dnakit.regulation import plot_categorical_attribution_scatter

attributions = [
    [0.5, 0.0, 0.1, -0.1], [0.0, 0.7, 0.0, 0.1], [0.1, 0.0, -0.4, 0.0],
    [0.0, 0.0, 0.0, 0.9], [0.2, 0.2, 0.0, 0.0],
]
result = plot_categorical_attribution_scatter(attributions, "attribution_scatter.png", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["dominant_base_index_by_position"])
print([round(value, 2) for value in result["data"]["summed_attribution_by_position"]])
```

- **Example output:**

```text
True
[0, 1, 2, 3, 0]
[0.5, 0.8, -0.3, 0.9, 0.4]
```

## 25) `dnakit.regulation.plot_pwm_logo()`

- **Purpose:** Draw a position weight matrix as an information-content logo and save it as a PNG.
- **Method:** Call tangermeme's `plot_pwm()` on the transposed matrix with the optional name and the logarithm guard `eps`.
- **API:** `plot_pwm_logo(pwm: 'Sequence[Sequence[float]]', output_path: 'str', *, name: 'str | None' = None, eps: 'float' = 1e-07, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages tangermeme and torch (extra `dnakit[bio-tangermeme]`), numpy and matplotlib. Works offline.
- **Input:** `pwm` has at most 10,000 rows of four probabilities for A, C, G, T; each row is non-negative and sums to 1 within 1e-4. `name` is an optional title of at most 200 characters. `eps` is a number between 0 and 1. `output_path` must be a new `.png` path.
- **Output:** `status`, `operation`, `path`, `image_format`, `data` (`position_count`, `alphabet`, `name`, `pwm`) and `provenance`.
- **Example code:**

```python
from dnakit.regulation import plot_pwm_logo

pwm = [[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25], [0.0, 1.0, 0.0, 0.0]]
result = plot_pwm_logo(pwm, "pwm_logo.png", name="demo motif", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["position_count"], result["data"]["name"])
```

- **Example output:**

```text
True
3 demo motif
```

## 26) `dnakit.regulation.plot_sequence_attributions()`

- **Purpose:** Compute DeepLIFT-SHAP attributions for DNA sequences with a torch model and draw one logo per sequence.
- **Method:** One-hot encode the sequences, then call tangermeme's `plot_attributions()` with `deep_lift_shap` as the attribution method and dinucleotide-shuffled references (`dinucleotide_shuffle`, seeded with `random_state`). The logos are stacked in one PNG. The model is not loaded by this function.
- **API:** `plot_sequence_attributions(model: 'Any', sequences: 'Sequence[str]', output_path: 'str', *, target: 'int' = 0, n_shuffles: 'int' = 20, random_state: 'int' = 0, hypothetical: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages tangermeme and torch (extra `dnakit[bio-tangermeme]`) and matplotlib. Works offline.
- **Input:** `model` is a `torch.nn.Module` already in eval mode. `sequences` is a list of 1 to 32 DNA strings of one common length (at most 10,000 bases). `target` is the non-negative index of the model output to explain. `n_shuffles` is between 1 and 200, `random_state` a non-negative integer, `hypothetical` a boolean. `output_path` must be a new `.png` path.
- **Output:** `status`, `operation`, `path`, `image_format`, `data` (`sequence_count`, `length`, `target`, `n_shuffles`, `random_state`, `hypothetical`, `sequences`, `attributions` as (sequence, position, base) lists) and `provenance`.
- **Example code:**

```python
import torch

from dnakit.regulation import plot_sequence_attributions


class TinyConv(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = torch.nn.Conv1d(4, 1, 3)
        with torch.no_grad():
            self.conv.weight.copy_(torch.tensor([[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0], [0, 0, 0]]]))
            self.conv.bias.zero_()

    def forward(self, x):
        return self.conv(x).sum(dim=-1)


result = plot_sequence_attributions(
    TinyConv().eval(), ["ACGTA", "TTGCA"], "sequence_attributions.png", n_shuffles=5, verbose=False
)
data = result["data"]
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(data["sequence_count"], data["length"], len(data["attributions"]), len(data["attributions"][0]), len(data["attributions"][0][0]))
```

- **Example output:**

```text
True
2 5 2 5 4
```

## 27) `dnakit.regulation.render_interactive_logo_png()`

- **Purpose:** Draw tangermeme's interactive attribution logo to a static PNG (the hover tooltips are not captured).
- **Method:** Call tangermeme's `interactive_logo()` on a matplotlib axis and save the figure. tangermeme's interactive logo needs mpld3; when it cannot be imported, the function raises `BackendUnavailableError` and writes nothing.
- **API:** `render_interactive_logo_png(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, start: 'int | None' = None, end: 'int | None' = None, ylim: 'Sequence[float] | None' = None, annotations: 'Sequence[Mapping[str, Any]] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages tangermeme, torch and mpld3, and matplotlib. `pyproject.toml` lists mpld3 in the `dnakit[bio-tangermeme]` extra, but the function's docstring and error hint say the extra does not include it; if the error appears, install it with `python -m pip install mpld3`. mpld3 is not installed in the environment where this page was checked, so the example shows the error.
- **Input:** The same as `plot_attribution_logo()`, except that `ylim` is a `(low, high)` pair with `low` below `high`.
- **Output:** `status`, `operation`, `path`, `image_format`, `data` (`position_count`, `start`, `end`, `ylim`, `alphabet`, `attributions`) and `provenance`.
- **Example code:**

```python
from dnakit.exceptions import BackendUnavailableError
from dnakit.regulation import render_interactive_logo_png

attributions = [[0.5, 0.0, 0.1, -0.1], [0.0, 0.7, 0.0, 0.1], [0.1, 0.0, -0.4, 0.0]]
try:
    result = render_interactive_logo_png(attributions, "interactive_logo.png", verbose=False)
    print("drawn", result["data"]["position_count"])
except BackendUnavailableError as error:
    print(type(error).__name__, error.message.split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError mpld3 could not be imported
```

The output above is from an environment without mpld3. With mpld3 installed the example prints `drawn 3`. Usage with the dependency installed:

```python
from dnakit.regulation import render_interactive_logo_png

attributions = [[0.5, 0.0, 0.1, -0.1], [0.0, 0.7, 0.0, 0.1], [0.1, 0.0, -0.4, 0.0]]
result = render_interactive_logo_png(attributions, "interactive_logo.png", ylim=[-1.0, 1.0])
print(result["data"]["position_count"])
```

## 28) `dnakit.predictions.download_genomic_benchmark_cloud_cache()`

- **Purpose:** Download a pre-built full-sequence benchmark dataset from the Zenodo cloud cache and describe it.
- **Method:** Check that the (name, version) pair is in the package's `CLOUD_CACHE`, then call `genomic_benchmarks.loc2seq.cloud_caching.download_from_cloud_cache()`, which downloads a zip archive and unpacks it into the parent of `dest_path`. Names that are not in the cache are rejected before any connection is made. With `force_download` an existing `dest_path` is removed first.
- **API:** `download_genomic_benchmark_cloud_cache(dataset_name: 'str', version: 'int', dest_path: 'str', *, force_download: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and internet access.
- **Input:** `dataset_name` and `version` form a pair that is in the cloud cache (for example `demo_human_or_worm` with version 0). `dest_path` is the dataset folder to create; its parent folder is created when missing. `force_download` is a boolean.
- **Output:** `service`, `source` (the Zenodo record URL) and `result` with `dataset`, `local_path` and `summary` (`record_count`, `split_counts`, `label_counts`, `first_record` with `split`, `label`, `length` and `sequence`).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import download_genomic_benchmark_cloud_cache

try:
    download_genomic_benchmark_cloud_cache("not_a_dataset", 0, "full_data", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **Example output:**

```text
ConfigurationError not_a_dataset version 0 is not in the genomic_benchmarks cloud cache.
```

The example does not download anything: names outside the cloud cache are rejected first. Usage with the dependency installed:

```python
from dnakit.predictions import download_genomic_benchmark_cloud_cache

result = download_genomic_benchmark_cloud_cache("demo_human_or_worm", 0, "data/demo_human_or_worm")
print(result["result"]["summary"]["record_count"])
```

## 29) `dnakit.predictions.download_genomic_benchmark_sequences()`

- **Purpose:** Build the full-sequence version of an interval-list benchmark dataset and describe it.
- **Method:** Call `genomic_benchmarks.loc2seq.download_dataset()`, which reads the dataset's `metadata.yaml` and interval tables, takes the reference genomes from `cache_path` (downloading them when missing), cuts the sequences and writes `<dest_path>/<name>/<split>/<label>/<id>.txt`. With the cloud cache enabled, a pre-built copy is used when one exists for the name and version. The function then counts the files by split and label.
- **API:** `download_genomic_benchmark_sequences(interval_list_dataset: 'str', *, version: 'int | None' = None, dest_path: 'str | None' = None, cache_path: 'str | None' = None, force_download: 'bool' = False, use_cloud_cache: 'bool' = True, local_repo: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`). Internet access is needed for a benchmark name from the package's list and for reference genomes that are not yet in `cache_path`; the example below uses a local dataset folder and a local reference, so it runs offline.
- **Input:** `interval_list_dataset` is a benchmark name or, with `local_repo=True`, the path of a local dataset folder with `metadata.yaml` and `<split>/<class>.csv.gz` files. `version` must match the metadata when given. `dest_path` and `cache_path` default to folders under `~/.genomic_benchmarks`. `force_download` forces new reference downloads. `use_cloud_cache` and `local_repo` are booleans.
- **Output:** `service`, `source` (the dataset name) and `result` with `dataset`, `local_path` and `summary` (`record_count`, `split_counts`, `label_counts` by class folder name, `first_record`).
- **Example code:**

```python
import contextlib
import gzip
import io
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import download_genomic_benchmark_sequences

Path("refs").mkdir()
with gzip.open("refs/mini.fa.gz", "wt") as handle:
    handle.write(">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTT\n")
dataset = Path("mini_ds")
for split in ("train", "test"):
    (dataset / split).mkdir(parents=True)
    with gzip.open(dataset / split / "pos.csv.gz", "wt") as handle:
        handle.write("id,region,start,end,strand\n0,chr1,0,8,+\n1,chr2,4,12,-\n")
(dataset / "metadata.yaml").write_text(
    "version: 0\nclasses:\n  pos:\n    url: http://example.invalid/mini.fa.gz\n    type: fa.gz\n"
)
with contextlib.redirect_stdout(io.StringIO()):
    result = download_genomic_benchmark_sequences(
        "mini_ds", version=0, dest_path="full", cache_path="refs",
        use_cloud_cache=False, local_repo=True, verbose=False,
    )
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["split_counts"], summary["label_counts"])
print(summary["first_record"]["length"], summary["first_record"]["sequence"])
```

- **Example output:**

```text
mini_ds 4 {'train': 2, 'test': 2} {'pos': 4}
8 ACGTACGT
```

The example writes a two-sequence reference and a two-row interval table; the second interval is on the minus strand, so its sequence is reverse-complemented.

## 30) `dnakit.predictions.get_genomic_benchmark_dataset()`

- **Purpose:** Load a named genomic-benchmarks classification dataset for one split and describe it.
- **Method:** Call `genomic_benchmarks.dataset_getters.pytorch_datasets.get_dataset()`. If `~/.genomic_benchmarks/<dataset_name>` is missing (or `force_download` is set), the package downloads the dataset first, from the Zenodo cloud cache when `use_cloud_cache` is true and the name is listed. The function counts the records and labels and returns the first record; the PyTorch dataset object is not returned. The label is the index of the class folder in the order the file system lists them, so the numbering is not guaranteed to be sorted.
- **API:** `get_genomic_benchmark_dataset(dataset_name: 'str', split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example below places a tiny dataset in a folder that stands in for the home directory, so it runs offline.
- **Input:** `dataset_name` is a non-empty benchmark name. `split` is `train` or `test`. `version` is a non-negative integer or `None`. `force_download` and `use_cloud_cache` are booleans.
- **Output:** `service`, `source` (the cloud-cache URL, or the GitHub dataset URL for names outside the cache) and `result` with `dataset`, `split`, `local_path` and `summary` (`record_count`, `split_counts`, `label_counts`, `first_record`).
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import get_genomic_benchmark_dataset

for split in ("train", "test"):
    for label, sequence in (("alpha", "ACGTACGT"), ("beta", "TTTTCCCC")):
        folder = Path(".genomic_benchmarks/tiny_set") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text(sequence)
result = get_genomic_benchmark_dataset("tiny_set", "train", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], result["result"]["split"])
print(summary["record_count"], summary["split_counts"], summary["label_counts"], summary["first_record"]["length"])
```

- **Example output:**

```text
tiny_set train
4 {'train': 4} {'0': 2, '1': 2} 8
```

## 31) `dnakit.predictions.load_demo_coding_vs_intergenomic_seqs()`

- **Purpose:** Load the demo dataset of coding versus intergenomic sequences for one split and describe it.
- **Method:** Call the package's `DemoCodingVsIntergenomicSeqs` getter (dataset `demo_coding_vs_intergenomic_seqs`) and summarise it as in `get_genomic_benchmark_dataset()`. The dataset is downloaded to `~/.genomic_benchmarks` when it is not on disk.
- **API:** `load_demo_coding_vs_intergenomic_seqs(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `demo_coding_vs_intergenomic_seqs`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_demo_coding_vs_intergenomic_seqs

for split in ("train", "test"):
    for label in ("coding", "intergenomic"):
        folder = Path(".genomic_benchmarks/demo_coding_vs_intergenomic_seqs") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_demo_coding_vs_intergenomic_seqs("test", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], result["source"])
print(summary["record_count"], summary["split_counts"], summary["label_counts"])
```

- **Example output:**

```text
demo_coding_vs_intergenomic_seqs https://zenodo.org/records/16605299
4 {'test': 4} {'0': 2, '1': 2}
```

## 32) `dnakit.predictions.load_demo_human_or_worm()`

- **Purpose:** Load the demo dataset that separates human from worm sequences for one split and describe it.
- **Method:** Call the package's `DemoHumanOrWorm` getter (dataset `demo_human_or_worm`) and summarise it as in `get_genomic_benchmark_dataset()`.
- **API:** `load_demo_human_or_worm(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `demo_human_or_worm`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_demo_human_or_worm

for split in ("train", "test"):
    for label in ("human", "worm"):
        folder = Path(".genomic_benchmarks/demo_human_or_worm") / split / label
        folder.mkdir(parents=True)
        for index in range(3):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_demo_human_or_worm("train", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
print(summary["first_record"]["length"])
```

- **Example output:**

```text
demo_human_or_worm 6 {'0': 3, '1': 3}
8
```

## 33) `dnakit.predictions.load_demo_mouse_enhancers()`

- **Purpose:** Load the demo mouse enhancer dataset for one split and describe it.
- **Method:** Call the package's `DemoMouseEnhancers` getter. The dataset is stored upstream as `dummy_mouse_enhancers_ensembl`, and that is the folder name used under `~/.genomic_benchmarks`. The summary is built as in `get_genomic_benchmark_dataset()`.
- **API:** `load_demo_mouse_enhancers(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `dummy_mouse_enhancers_ensembl`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_demo_mouse_enhancers

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/dummy_mouse_enhancers_ensembl") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_demo_mouse_enhancers("test", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["split_counts"])
```

- **Example output:**

```text
dummy_mouse_enhancers_ensembl 4 {'test': 4}
```

## 34) `dnakit.predictions.load_drosophila_enhancers_stark()`

- **Purpose:** Load the Drosophila enhancer dataset of Stark et al. for one split and describe it.
- **Method:** Call the package's `DrosophilaEnhancersStark` getter (dataset `drosophila_enhancers_stark`). That upstream getter has no `use_cloud_cache` argument and always uses the cloud cache when it can, so this function has none either. The summary is built as in `get_genomic_benchmark_dataset()`.
- **API:** `load_drosophila_enhancers_stark(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download` and `version` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `drosophila_enhancers_stark`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_drosophila_enhancers_stark

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/drosophila_enhancers_stark") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_drosophila_enhancers_stark("train", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
```

- **Example output:**

```text
drosophila_enhancers_stark 4 {'0': 2, '1': 2}
```

## 35) `dnakit.predictions.load_human_enhancers_cohn()`

- **Purpose:** Load the human enhancer dataset of Cohn et al. for one split and describe it.
- **Method:** Call the package's `HumanEnhancersCohn` getter (dataset `human_enhancers_cohn`) and summarise it as in `get_genomic_benchmark_dataset()`.
- **API:** `load_human_enhancers_cohn(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `human_enhancers_cohn`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_human_enhancers_cohn

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/human_enhancers_cohn") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_human_enhancers_cohn("test", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
```

- **Example output:**

```text
human_enhancers_cohn 4 {'0': 2, '1': 2}
```

## 36) `dnakit.predictions.load_human_enhancers_ensembl()`

- **Purpose:** Load the human enhancer dataset built from Ensembl regulatory regions for one split and describe it.
- **Method:** Call the package's `HumanEnhancersEnsembl` getter (dataset `human_enhancers_ensembl`) and summarise it as in `get_genomic_benchmark_dataset()`.
- **API:** `load_human_enhancers_ensembl(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `human_enhancers_ensembl`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_human_enhancers_ensembl

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/human_enhancers_ensembl") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_human_enhancers_ensembl("train", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
```

- **Example output:**

```text
human_enhancers_ensembl 4 {'0': 2, '1': 2}
```

## 37) `dnakit.predictions.load_human_nontata_promoters()`

- **Purpose:** Load the human TATA-less promoter dataset for one split and describe it.
- **Method:** Call the package's `HumanNontataPromoters` getter (dataset `human_nontata_promoters`) and summarise it as in `get_genomic_benchmark_dataset()`.
- **API:** `load_human_nontata_promoters(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `human_nontata_promoters`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_human_nontata_promoters

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/human_nontata_promoters") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_human_nontata_promoters("test", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
```

- **Example output:**

```text
human_nontata_promoters 4 {'0': 2, '1': 2}
```

## 38) `dnakit.predictions.load_human_ocr_ensembl()`

- **Purpose:** Load the human open-chromatin region dataset built from Ensembl for one split and describe it.
- **Method:** Call the package's `HumanOcrEnsembl` getter (dataset `human_ocr_ensembl`) and summarise it as in `get_genomic_benchmark_dataset()`.
- **API:** `load_human_ocr_ensembl(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package genomic-benchmarks (extra `dnakit[bio-genomic-benchmarks]`) and torch, which the extra does not install. Needs internet access unless the dataset is already on disk; the example places a tiny stand-in dataset in a folder used as the home directory, so it runs offline.
- **Input:** `split` is `train` or `test`. `force_download`, `version` and `use_cloud_cache` are as in `get_genomic_benchmark_dataset()`.
- **Output:** The same structure as `get_genomic_benchmark_dataset()`, with `dataset` set to `human_ocr_ensembl`.
- **Example code:**

```python
import os
from pathlib import Path

os.environ["HOME"] = os.getcwd()
from dnakit.predictions import load_human_ocr_ensembl

for split in ("train", "test"):
    for label in ("negative", "positive"):
        folder = Path(".genomic_benchmarks/human_ocr_ensembl") / split / label
        folder.mkdir(parents=True)
        for index in range(2):
            (folder / f"{index}.txt").write_text("ACGTACGT")
result = load_human_ocr_ensembl("train", verbose=False)
summary = result["result"]["summary"]
print(result["result"]["dataset"], summary["record_count"], summary["label_counts"])
```

- **Example output:**

```text
human_ocr_ensembl 4 {'0': 2, '1': 2}
```
