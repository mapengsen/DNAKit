# SnapATAC2: datasets and single-cell chromatin plots

This group of functions in `dnakit.epigenome` is a set of adapters over [SnapATAC2](https://github.com/scverse/SnapATAC2), the single-cell chromatin accessibility toolkit (MIT licence in the package metadata). Nine functions download example datasets, motif collections and genome files through the snapatac2 registry; fifteen functions turn data or Plotly-style figures into PNG plots (heatmap, density contour, scatter, 2D/3D embedding, network scores, quality-control curves, peak accessibility, motif enrichment and coverage tracks). They only fetch or draw data: nothing here runs peak calling, clustering or other analysis.

The functions are imported with `from dnakit.epigenome import ...`. The dependencies are optional and are imported only when the matching function is called; if one is missing, the function raises `BackendUnavailableError` with the install command. The plot functions draw the PNG with Matplotlib, because Plotly's own PNG export needs a Chrome installation. The download functions need internet access, so their examples show the error raised for an invalid argument or for a missing dependency; the real calls are given as unexecuted usage blocks.

## Example data

No example file is needed. Every example passes small inputs directly or writes a tiny fragment file to a temporary directory.

## 1) `dnakit.epigenome.load_snapatac2_dataset()`

- **Purpose:** Fetch one file of a snapatac2 registry dataset into the local cache.
- **Method:** Checks that `dataset` (and `file`, if given) are non-empty strings, then calls `snapatac2.datasets.load(dataset[, file=...])`, which downloads through pooch. If the result is a path or a list of paths they are returned in `paths`; otherwise the object type and, for a list, its length are returned in `data`. An unknown dataset or file name (`KeyError`) and a rejected request (`TypeError`) become `ConfigurationError`; `OSError`, `ValueError` and `RuntimeError` become `DownloadError` with a hint about `SNAP_DATA_DIR`.
- **API:** `dnakit.epigenome.load_snapatac2_dataset(dataset: 'str', *, file: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `dataset` is a registry name: `HEA_cCRE`, `Meuleman_2020`, `atac_pbmc_500`, `atac_pbmc_5k`, `cisBP_human`, `colon_transverse`, `gencode_GRCh37`, `gencode_GRCh38`, `gencode_GRCm38`, `gencode_GRCm39` or `pbmc_10k_multiome`. `file` optionally names one file of that dataset. A `ConfigurationError` is raised for invalid arguments before any download; unknown registry names also raise `ConfigurationError`, and network or file failures raise `DownloadError`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `parameters` is `{"file": file}`; `data` is empty unless snapatac2 returned a non-path object.
- **Example code:**

```python
from dnakit.epigenome import load_snapatac2_dataset
from dnakit.exceptions import ConfigurationError

try:
    load_snapatac2_dataset("no_such_dataset", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error.args[0]).split(":")[0])
```

- **Example output:**

```text
ConfigurationError snapatac2 does not know the dataset or file 'no_such_dataset'
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import load_snapatac2_dataset

result = load_snapatac2_dataset("atac_pbmc_500", file="atac_pbmc_500_downsample.tsv.gz")
print(result["paths"])
```

## 2) `dnakit.epigenome.fetch_snapatac2_pbmc500()`

- **Purpose:** Download the 10x 500-cell PBMC scATAC-seq example as fragment, BAM or FASTQ files.
- **Method:** Validates the arguments, then calls `snapatac2.datasets.pbmc500(type=data_type, downsample=downsample)` and returns the local paths under the dataset name `atac_pbmc_500`.
- **API:** `dnakit.epigenome.fetch_snapatac2_pbmc500(*, data_type: 'str' = 'fragment', downsample: bool = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `data_type` is `fragment`, `bam` or `fastq`. `downsample` must be a boolean and may be true only with `data_type="fragment"`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `parameters` holds `data_type` and `downsample`; `data` is empty.
- **Example code:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc500
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc500(data_type="bam", downsample=True, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **Example output:**

```text
ConfigurationError downsample is only available with data_type='fragment'.
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc500

result = fetch_snapatac2_pbmc500(data_type="fragment", downsample=True)
print(result["dataset"], result["paths"])
```

## 3) `dnakit.epigenome.fetch_snapatac2_pbmc5k()`

- **Purpose:** Download the 10x 5k PBMC scATAC-seq example as a fragment, h5ad or annotated h5ad file.
- **Method:** Validates `data_type`, then calls `snapatac2.datasets.pbmc5k(type=data_type)` and returns the local paths under the dataset name `atac_pbmc_5k`.
- **API:** `dnakit.epigenome.fetch_snapatac2_pbmc5k(*, data_type: 'str' = 'fragment', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `data_type` is `fragment`, `h5ad` or `annotated_h5ad`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `parameters` holds `data_type`; `data` is empty.
- **Example code:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc5k
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc5k(data_type="csv", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **Example output:**

```text
ConfigurationError data_type must be one of: fragment, h5ad, annotated_h5ad.
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc5k

result = fetch_snapatac2_pbmc5k(data_type="h5ad")
print(result["dataset"], result["paths"])
```

## 4) `dnakit.epigenome.fetch_snapatac2_pbmc10k_multiome()`

- **Purpose:** Download the RNA or ATAC part of the 10x 10k PBMC multiome example.
- **Method:** Validates both choices, then calls `snapatac2.datasets.pbmc10k_multiome(modality=..., type=...)` and returns the local paths under the dataset name `pbmc_10k_multiome`.
- **API:** `dnakit.epigenome.fetch_snapatac2_pbmc10k_multiome(*, modality: 'str' = 'RNA', data_type: 'str' = 'h5ad', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `modality` is `RNA` or `ATAC`; `data_type` is `fragment` or `h5ad`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `parameters` holds `modality` and `data_type`; `data` is empty.
- **Example code:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc10k_multiome
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc10k_multiome(modality="ADT", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **Example output:**

```text
ConfigurationError modality must be one of: RNA, ATAC.
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc10k_multiome

result = fetch_snapatac2_pbmc10k_multiome(modality="ATAC", data_type="h5ad")
print(result["dataset"], result["paths"])
```

## 5) `dnakit.epigenome.fetch_snapatac2_colon_transverse()`

- **Purpose:** Download and extract the transverse colon scATAC-seq fragment files, one per sample.
- **Method:** Calls `snapatac2.datasets.colon()`, which yields (sample name, fragment file) pairs, and lists them.
- **API:** `dnakit.epigenome.fetch_snapatac2_colon_transverse(*, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example simulates a missing `snapatac2.datasets` module (a `None` entry in `sys.modules`) so that it works offline and shows the error raised when the dependency is absent. To run it for real, install the dependency and have internet access.
- **Input:** No data arguments. A missing snapatac2 raises `BackendUnavailableError` with the install command; download failures raise `DownloadError`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `paths` lists the fragment files and `data["samples"]` lists `{"sample", "fragment_file"}` rows. The dataset name is `colon_transverse`.
- **Example code:**

```python
import sys

# Simulate a missing snapatac2.datasets: a None entry makes the import fail.
sys.modules["snapatac2.datasets"] = None

from dnakit.epigenome import fetch_snapatac2_colon_transverse
from dnakit.exceptions import BackendUnavailableError

try:
    fetch_snapatac2_colon_transverse(verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.args[0].split(":")[0])
    print(error.hint)
```

- **Example output:**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_colon_transverse

result = fetch_snapatac2_colon_transverse()
for row in result["data"]["samples"]:
    print(row["sample"], row["fragment_file"])
```

## 6) `dnakit.epigenome.fetch_snapatac2_hea_ccre()`

- **Purpose:** Download the curated human colon cis-regulatory element BED file.
- **Method:** Calls `snapatac2.datasets.cre_HEA()` and returns the local path under the dataset name `HEA_cCRE`.
- **API:** `dnakit.epigenome.fetch_snapatac2_hea_ccre(*, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example simulates a missing `snapatac2.datasets` module (a `None` entry in `sys.modules`) so that it works offline and shows the error raised when the dependency is absent. To run it for real, install the dependency and have internet access.
- **Input:** No data arguments. A missing snapatac2 raises `BackendUnavailableError` with the install command; download failures raise `DownloadError`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. `paths` holds the BED file; `data` is empty.
- **Example code:**

```python
import sys

# Simulate a missing snapatac2.datasets: a None entry makes the import fail.
sys.modules["snapatac2.datasets"] = None

from dnakit.epigenome import fetch_snapatac2_hea_ccre
from dnakit.exceptions import BackendUnavailableError

try:
    fetch_snapatac2_hea_ccre(verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.args[0].split(":")[0])
    print(error.hint)
```

- **Example output:**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_hea_ccre

result = fetch_snapatac2_hea_ccre()
print(result["dataset"], result["paths"])
```

## 7) `dnakit.epigenome.fetch_snapatac2_cis_bp_motifs()`

- **Purpose:** Download the CIS-BP transcription factor motifs and return their names and identifiers.
- **Method:** Calls `snapatac2.datasets.cis_bp(unique=unique)` and summarises each returned motif by its `id` and `name`. The motif objects themselves are not returned.
- **API:** `dnakit.epigenome.fetch_snapatac2_cis_bp_motifs(*, unique: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `unique` must be a boolean (snapatac2 then keeps one motif per transcription factor).
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. The dataset name is `cisBP_human`; `paths` is empty; `data` holds `registry_file` (`cisBP_human.meme`), `motif_count` and `motifs` (a list of `{"id", "name"}`).
- **Example code:**

```python
from dnakit.epigenome import fetch_snapatac2_cis_bp_motifs
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_cis_bp_motifs(unique="yes", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **Example output:**

```text
ConfigurationError unique must be a boolean.
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_cis_bp_motifs

result = fetch_snapatac2_cis_bp_motifs(unique=True)
print(result["data"]["motif_count"], result["data"]["motifs"][0])
```

## 8) `dnakit.epigenome.fetch_snapatac2_meuleman_2020_motifs()`

- **Purpose:** Download the Meuleman 2020 grouped motifs and return their names, identifiers and families.
- **Method:** Calls `snapatac2.datasets.Meuleman_2020()` and summarises each motif by `id`, `name` and, when the motif has one, `family`.
- **API:** `dnakit.epigenome.fetch_snapatac2_meuleman_2020_motifs(*, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example simulates a missing `snapatac2.datasets` module (a `None` entry in `sys.modules`) so that it works offline and shows the error raised when the dependency is absent. To run it for real, install the dependency and have internet access.
- **Input:** No data arguments. A missing snapatac2 raises `BackendUnavailableError` with the install command; download failures raise `DownloadError`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. The dataset name is `Meuleman_2020`; `paths` is empty; `data` holds `registry_file` (`Meuleman_2020.meme`), `motif_count` and `motifs`.
- **Example code:**

```python
import sys

# Simulate a missing snapatac2.datasets: a None entry makes the import fail.
sys.modules["snapatac2.datasets"] = None

from dnakit.epigenome import fetch_snapatac2_meuleman_2020_motifs
from dnakit.exceptions import BackendUnavailableError

try:
    fetch_snapatac2_meuleman_2020_motifs(verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.args[0].split(":")[0])
    print(error.hint)
```

- **Example output:**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_meuleman_2020_motifs

result = fetch_snapatac2_meuleman_2020_motifs()
print(result["data"]["motif_count"], result["data"]["motifs"][0])
```

## 9) `dnakit.epigenome.fetch_snapatac2_genome_files()`

- **Purpose:** Download the FASTA and/or gene annotation of a genome preset.
- **Method:** Maps the genome name or alias to a GENCODE registry dataset, then calls `snapatac2.genome.load(dataset, file=...)` for the FASTA (gzip-decompressed with `pooch.Decompress`) and/or the GFF3 annotation (kept as downloaded). The genome and `kind` are checked before snapatac2 is imported.
- **API:** `dnakit.epigenome.fetch_snapatac2_genome_files(genome: 'str', *, kind: 'str' = 'both', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`; pooch comes with it) and internet access for the first download. snapatac2 keeps the files in a cache folder that `SNAP_DATA_DIR` can change. No key is needed. The example needs no download: it shows the error raised for an invalid argument, which is checked before snapatac2 is called. To run it for real, install the dependency and have internet access.
- **Input:** `genome` is `GRCh37`, `GRCh38`, `GRCm39`, `GRCm38` or an alias (`hg19`, `hg38`, `mm39`, `mm10`). `kind` is `fasta`, `annotation` or `both`.
- **Output:** A dictionary with `status` (`"completed"`), `service`, `dataset`, `parameters`, `paths` (local file paths), `retrieved_utc` (UTC time of the call) and `data`. The dataset name is `gencode_<genome>`; `parameters` holds `genome` (canonical name), `requested_name` and `kind`; `data` maps `fasta` and/or `annotation` to the local path.
- **Example code:**

```python
from dnakit.epigenome import fetch_snapatac2_genome_files
from dnakit.exceptions import ConfigurationError

for arguments in ({"genome": "hg00"}, {"genome": "hg38", "kind": "bed"}):
    try:
        fetch_snapatac2_genome_files(verbose=False, **arguments)
    except ConfigurationError as error:
        print(type(error).__name__, error.args[0])
```

- **Example output:**

```text
ConfigurationError genome must be one of: GRCh37, GRCh38, GRCm38, GRCm39, hg19, hg38, mm10, mm39.
ConfigurationError kind must be one of: fasta, annotation, both.
```

- **Usage with the dependency installed:**

```python
from dnakit.epigenome import fetch_snapatac2_genome_files

result = fetch_snapatac2_genome_files("hg38", kind="annotation")
print(result["parameters"]["genome"], result["data"]["annotation"])
```

## 10) `dnakit.epigenome.render_snapatac2_figure_png()`

- **Purpose:** Render a Plotly figure dictionary with snapatac2's layout rules and write it to a new PNG file.
- **Method:** Builds a `plotly.graph_objects.Figure` from `figure["data"]` and `figure["layout"]`, passes it to `snapatac2.plotting.render_plot(width=..., height=..., show=False)` and draws the result with Matplotlib. Only scatter, bar, violin, heatmap, contour and scatter3d traces can be drawn; other trace types raise `ConfigurationError`.
- **API:** `dnakit.epigenome.render_snapatac2_figure_png(figure: 'Mapping[str, Any]', output_png: 'str', *, width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `figure` is a mapping with a `data` list and an optional `layout`. `width` and `height` are positive integers. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `width` and `height`; `operation` is `render_plot`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import render_snapatac2_figure_png

figure = {
    "data": [{"type": "bar", "x": ["a", "b", "c"], "y": [1, 3, 2]}],
    "layout": {"title": {"text": "Demo"}, "xaxis": {"title": {"text": "Group"}}},
}
out = Path(tempfile.mkdtemp())
result = render_snapatac2_figure_png(figure, out / "figure.png", width=500, height=300, verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print(result["data"]["title"], result["data"]["xaxis_title"])
print(result["data"]["traces"])
```

- **Example output:**

```text
completed render_plot {'width': 500, 'height': 300}
True True
Demo Group
[{'type': 'bar', 'name': None, 'x': ['a', 'b', 'c'], 'y': [1, 3, 2]}]
```

## 11) `dnakit.epigenome.plot_snapatac2_heatmap_png()`

- **Purpose:** Draw a numeric matrix as a heatmap with optional hierarchical clustering dendrograms.
- **Method:** Calls `snapatac2.plotting.heatmap(...)` (which clusters with SciPy) and snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG. When rows or columns are not clustered, the given names are set as tick labels.
- **API:** `dnakit.epigenome.plot_snapatac2_heatmap_png(matrix: 'Sequence[Sequence[float]]', output_png: 'str', *, row_names: 'Sequence[str] | None' = None, column_names: 'Sequence[str] | None' = None, cluster_rows: 'bool' = True, cluster_columns: 'bool' = True, colorscale: 'str' = 'Blues', linkage: 'str' = 'ward', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `matrix` is a non-empty two-dimensional matrix of finite numbers. `row_names` and `column_names` need one name per row or column. Clustering needs at least two rows or columns. `linkage` is one of `single`, `complete`, `average`, `weighted`, `centroid`, `median`, `ward`. `colorscale` is a Plotly or Matplotlib colour-map name (an unknown name is drawn with viridis). `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `rows`, `columns`, `cluster_rows`, `cluster_columns`, `colorscale` and `linkage`. With clustering, the tick texts follow the clustered order and the dendrograms are extra scatter traces.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_heatmap_png

out = Path(tempfile.mkdtemp())
result = plot_snapatac2_heatmap_png(
    [[1, 2, 3], [4, 5, 6], [7, 8, 10]],
    out / "heatmap.png",
    row_names=["r1", "r2", "r3"],
    column_names=["c1", "c2", "c3"],
    verbose=False,
)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print([trace["type"] for trace in result["data"]["traces"]])
print(result["data"]["xaxis_ticktext"], result["data"]["yaxis_ticktext"])
```

- **Example output:**

```text
completed heatmap {'rows': 3, 'columns': 3, 'cluster_rows': True, 'cluster_columns': True, 'colorscale': 'Blues', 'linkage': 'ward'}
True True
['scatter', 'scatter', 'scatter', 'scatter', 'heatmap']
['c3', 'c1', 'c2'] ['r3', 'r1', 'r2']
```

## 12) `dnakit.epigenome.plot_snapatac2_kde2d_png()`

- **Purpose:** Draw a two-dimensional kernel density contour of paired observations.
- **Method:** Calls `snapatac2.plotting.kde2d(x, y, log_x=..., log_y=...)`, which returns a density contour, and snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_kde2d_png(x: 'Sequence[float]', y: 'Sequence[float]', output_png: 'str', *, log_x: 'bool' = False, log_y: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `x` and `y` are equally long sequences of at least three finite numbers; with `log_x` or `log_y` the matching sequence must be positive. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `points`, `log_x` and `log_y`; the trace is a `contour`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_kde2d_png

out = Path(tempfile.mkdtemp())
x = [1, 2, 3, 4, 5, 6, 7, 8]
y = [2, 1, 4, 3, 6, 5, 8, 9]
result = plot_snapatac2_kde2d_png(x, y, out / "kde2d.png", log_x=True, verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print([trace["type"] for trace in result["data"]["traces"]], result["data"]["xaxis_type"])
```

- **Example output:**

```text
completed kde2d {'points': 8, 'log_x': True, 'log_y': False}
True True
['contour'] log
```

## 13) `dnakit.epigenome.plot_snapatac2_scatter_png()`

- **Purpose:** Draw a two-dimensional scatter plot with optional categorical or numeric colours.
- **Method:** Calls `snapatac2.plotting.scatter(...)` with `show=False`. Text colours become one trace per category; all-numeric colours become a colour bar. snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_scatter_png(x: 'Sequence[float]', y: 'Sequence[float]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, x_label: 'str' = 'X', y_label: 'str' = 'Y', color_label: 'str' = 'Color', marker_size: 'float' = 2.0, marker_opacity: 'float' = 0.5, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `x` and `y` are equally long sequences of finite numbers. `color` has one value per point (all numbers or text). `x_label`, `y_label` and `color_label` are non-empty strings. `marker_size` must be positive and `marker_opacity` between 0 and 1. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `points`, `marker_size` and `marker_opacity`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_scatter_png

out = Path(tempfile.mkdtemp())
result = plot_snapatac2_scatter_png(
    [1, 2, 3, 4],
    [4, 3, 2, 1],
    out / "scatter.png",
    color=["a", "a", "b", "b"],
    x_label="PC1",
    y_label="PC2",
    verbose=False,
)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print([(trace["name"], trace["x"], trace["y"]) for trace in result["data"]["traces"]])
print(result["data"]["xaxis_title"], result["data"]["yaxis_title"])
```

- **Example output:**

```text
completed scatter {'points': 4, 'marker_size': 2.0, 'marker_opacity': 0.5}
True True
[('a', [1.0, 2.0], [4.0, 3.0]), ('b', [3.0, 4.0], [2.0, 1.0])]
PC1 PC2
```

## 14) `dnakit.epigenome.plot_snapatac2_scatter3d_png()`

- **Purpose:** Draw a three-dimensional scatter plot with optional categorical or numeric colours.
- **Method:** Calls `snapatac2.plotting.scatter3d(...)` with `show=False` and draws the plotted points on Matplotlib 3D axes. The 3D drawing uses one colour per trace, so numeric colours are kept in the returned figure data but not as a colour bar in the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_scatter3d_png(x: 'Sequence[float]', y: 'Sequence[float]', z: 'Sequence[float]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, x_label: 'str' = 'X', y_label: 'str' = 'Y', z_label: 'str' = 'Z', color_label: 'str' = 'Color', marker_size: 'float' = 2.0, marker_opacity: 'float' = 0.5, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `x`, `y` and `z` are equally long sequences of finite numbers. `color`, labels, `marker_size` and `marker_opacity` follow the same rules as in the 2D scatter plot. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `points`, `marker_size` and `marker_opacity`; the trace type is `scatter3d`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_scatter3d_png

out = Path(tempfile.mkdtemp())
result = plot_snapatac2_scatter3d_png(
    [1, 2, 3, 4],
    [4, 3, 2, 1],
    [1, 1, 2, 2],
    out / "scatter3d.png",
    color=[0.1, 0.2, 0.3, 0.4],
    verbose=False,
)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
trace = result["data"]["traces"][0]
print(trace["type"], trace["x"], trace["y"], trace["z"])
```

- **Example output:**

```text
completed scatter3d {'points': 4, 'marker_size': 2.0, 'marker_opacity': 0.5}
True True
scatter3d [1.0, 2.0, 3.0, 4.0] [4.0, 3.0, 2.0, 1.0] [1.0, 1.0, 2.0, 2.0]
```

## 15) `dnakit.epigenome.plot_snapatac2_network_scores_png()`

- **Purpose:** Plot the mean edge score of a regulatory network in distance-to-TSS bins.
- **Method:** Builds a `rustworkx.PyDiGraph` from the node types and edges, calls `snapatac2.plotting.network_scores(graph, score_name=...)`, and draws the bar chart with Matplotlib.
- **API:** `dnakit.epigenome.plot_snapatac2_network_scores_png(node_types: 'Sequence[str]', edges: 'Sequence[Mapping[str, Any]]', output_png: 'str', *, score_name: 'str', width: 'int' = 800, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed. rustworkx comes with snapatac2.
- **Input:** `node_types` is a non-empty list of strings. `edges` is a non-empty list of mappings with integer `source` and `target` node indexes and optional numeric `distance`, `cor_score` and `regr_score`. Every edge needs a `distance`, and at least one edge needs the chosen score. `score_name` is `cor_score` or `regr_score`. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `score_name` and `edges`; the trace is a bar chart with distance-bin labels as x and mean scores as y.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_network_scores_png

types = ["gene", "peak", "peak", "gene"]
edges = [
    {"source": 1, "target": 0, "distance": 1000, "cor_score": 0.5, "regr_score": 0.2},
    {"source": 2, "target": 0, "distance": 5000, "cor_score": 0.3, "regr_score": 0.1},
    {"source": 2, "target": 3, "distance": 20000, "cor_score": -0.2, "regr_score": 0.4},
    {"source": 1, "target": 3, "distance": 30000, "cor_score": 0.1, "regr_score": 0.3},
]
out = Path(tempfile.mkdtemp())
result = plot_snapatac2_network_scores_png(types, edges, out / "scores.png", score_name="cor_score", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
trace = result["data"]["traces"][0]
print(trace["type"], trace["x"], trace["y"])
print(result["data"]["xaxis_title"], "|", result["data"]["yaxis_title"])
```

- **Example output:**

```text
completed network_scores {'score_name': 'cor_score', 'edges': 4}
True True
bar ['500 - 2K', '2K - 20K', '20K - 50K'] [0.5, 0.3, -0.05]
Distance to TSS (bp) | Average score
```

## 16) `dnakit.epigenome.plot_snapatac2_network_edge_stat_png()`

- **Purpose:** Draw correlation-score violins for each source and target node type pair.
- **Method:** Builds a `rustworkx.PyDiGraph` as in the network score plot, calls `snapatac2.plotting.network_edge_stat(graph)` and draws one violin per `source -> target` type pair with Matplotlib. A group with a single value or no spread is drawn as a point.
- **API:** `dnakit.epigenome.plot_snapatac2_network_edge_stat_png(node_types: 'Sequence[str]', edges: 'Sequence[Mapping[str, Any]]', output_png: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed. rustworkx comes with snapatac2.
- **Input:** `node_types` and `edges` follow the same rules as in the network score plot, except that `distance` is not required; at least one edge needs a `cor_score`. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `edges`; each violin trace is named `<source type> -> <target type>` and its `y` holds the correlation scores.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_network_edge_stat_png

types = ["gene", "peak", "peak", "gene"]
edges = [
    {"source": 1, "target": 0, "distance": 1000, "cor_score": 0.5, "regr_score": 0.2},
    {"source": 2, "target": 0, "distance": 5000, "cor_score": 0.3, "regr_score": 0.1},
    {"source": 2, "target": 3, "distance": 20000, "cor_score": -0.2, "regr_score": 0.4},
    {"source": 1, "target": 3, "distance": 30000, "cor_score": 0.1, "regr_score": 0.3},
]
out = Path(tempfile.mkdtemp())
result = plot_snapatac2_network_edge_stat_png(types, edges, out / "edge_stat.png", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
for trace in result["data"]["traces"]:
    print(trace["type"], trace["name"], trace["y"])
```

- **Example output:**

```text
completed network_edge_stat {'edges': 4}
True True
violin peak -> gene [0.5, 0.3, -0.2, 0.1]
```

## 17) `dnakit.epigenome.plot_snapatac2_valid_cells_png()`

- **Purpose:** Plot ranked barcode counts on log-log axes to choose a cell-calling cutoff.
- **Method:** Calls `snapatac2.plotting.valid_cells(counts, ...)`, which ranks the counts from high to low, and snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_valid_cells_png(counts: 'Sequence[float]', output_png: 'str', *, width: 'int' = 500, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `counts` has at least two positive finite numbers, one per barcode, in any order. `width` and `height` are positive integers. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `barcodes`; the trace has the rank as x and the sorted counts as y.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_valid_cells_png

out = Path(tempfile.mkdtemp())
counts = [10, 1000, 50, 800, 500, 100]
result = plot_snapatac2_valid_cells_png(counts, out / "valid_cells.png", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
trace = result["data"]["traces"][0]
print(trace["x"], trace["y"])
print(result["data"]["xaxis_type"], result["data"]["yaxis_type"])
```

- **Example output:**

```text
completed valid_cells {'barcodes': 6}
True True
[1.0, 2.0, 3.0, 4.0, 5.0, 6.0] [1000.0, 800.0, 500.0, 100.0, 50.0, 10.0]
log log
```

## 18) `dnakit.epigenome.plot_snapatac2_tsse_png()`

- **Purpose:** Plot TSS enrichment against the number of unique fragments per cell.
- **Method:** Puts the two sequences into the `obs` table of an in-memory `AnnData`, calls `snapatac2.plotting.tsse(adata, min_fragment=...)` (a density contour) and snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_tsse_png(n_fragment: 'Sequence[float]', tsse: 'Sequence[float]', output_png: 'str', *, min_fragment: 'int' = 500, width: 'int' = 500, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `n_fragment` and `tsse` are equally long sequences of at least two finite numbers, one per cell. `min_fragment` is a non-negative integer. `width` and `height` are positive integers. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `cells` and `min_fragment`; `data["cells_plotted"]` counts the cells with `n_fragment >= min_fragment`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_tsse_png

out = Path(tempfile.mkdtemp())
n_fragment = [1000, 2000, 3000, 600, 400]
tsse = [5, 8, 10, 6, 2]
result = plot_snapatac2_tsse_png(n_fragment, tsse, out / "tsse.png", min_fragment=500, verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print(result["data"]["cells_plotted"], [trace["type"] for trace in result["data"]["traces"]])
print(result["data"]["xaxis_title"], "|", result["data"]["yaxis_title"])
```

- **Example output:**

```text
completed tsse {'cells': 5, 'min_fragment': 500}
True True
4 ['contour']
Number of unique fragments | TSS enrichment score
```

## 19) `dnakit.epigenome.plot_snapatac2_frag_size_distr_png()`

- **Purpose:** Plot the fragment size distribution from counts per fragment length.
- **Method:** Calls `snapatac2.plotting.frag_size_distr(counts, max_recorded_size=...)`, a line plot of count against fragment size, and snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG. In the example below, the value at index 0 is not plotted and the sizes run from 1 to `max_recorded_size`.
- **API:** `dnakit.epigenome.plot_snapatac2_frag_size_distr_png(distribution: 'Sequence[float]', output_png: 'str', *, max_recorded_size: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `distribution` has at least two non-negative finite counts, where the index is the fragment length. `max_recorded_size` is a positive integer. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `max_recorded_size` and `sizes` (the number of counts given).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_frag_size_distr_png

out = Path(tempfile.mkdtemp())
counts = [0, 0, 5, 10, 30, 20, 10, 5, 2, 1]
result = plot_snapatac2_frag_size_distr_png(counts, out / "fragment_sizes.png", max_recorded_size=9, verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
trace = result["data"]["traces"][0]
print(trace["x"])
print(trace["y"])
```

- **Example output:**

```text
completed frag_size_distr {'max_recorded_size': 9, 'sizes': 10}
True True
[1, 2, 3, 4, 5, 6, 7, 8, 9]
[0.0, 5.0, 10.0, 30.0, 20.0, 10.0, 5.0, 2.0, 1.0]
```

## 20) `dnakit.epigenome.plot_snapatac2_spectral_eigenvalues_png()`

- **Purpose:** Plot spectral embedding eigenvalues and the elbow that snapatac2 picks.
- **Method:** Stores the values in `adata.uns["spectral_eigenvalue"]` of an empty `AnnData`, calls `snapatac2.plotting.spectral_eigenvalues(adata)`, and reads the elbow snapatac2 writes to `adata.uns["num_eigen"]`. snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_spectral_eigenvalues_png(eigenvalues: 'Sequence[float]', output_png: 'str', *, width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `eigenvalues` has at least two finite numbers. `width` and `height` are positive integers. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `eigenvalues` (the count); `data["num_eigen"]` is the elbow position, or `None` when snapatac2 finds no elbow (for example, a smoothly decreasing series).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_spectral_eigenvalues_png

out = Path(tempfile.mkdtemp())
eigenvalues = [10, 9.5, 9, 8.5, 2, 1.5, 1.2, 1.0, 0.9, 0.85, 0.8, 0.75]
result = plot_snapatac2_spectral_eigenvalues_png(eigenvalues, out / "eigenvalues.png", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print("elbow:", result["data"]["num_eigen"])
print(result["data"]["xaxis_title"], "|", result["data"]["yaxis_title"])
```

- **Example output:**

```text
completed spectral_eigenvalues {'eigenvalues': 12}
True True
elbow: 9
Component | Eigenvalue
```

## 21) `dnakit.epigenome.plot_snapatac2_regions_png()`

- **Purpose:** Draw the grouped, RPKM-normalised accessibility of selected peaks.
- **Method:** Builds an in-memory `AnnData` from the cell-by-peak matrix (peak names as variables, group labels in `obs[groupby]`) and calls `snapatac2.plotting.regions(adata, groupby=..., peaks=...)`, which returns a heatmap of groups against the selected peaks. snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_regions_png(matrix: 'Sequence[Sequence[float]]', peak_names: 'Sequence[str]', cell_groups: 'Sequence[str]', peaks: 'Mapping[str, Sequence[str]]', output_png: 'str', *, groupby: 'str' = 'group', width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `matrix` is a cells-by-peaks matrix of finite numbers. `peak_names` has one unique name per column in the form `chrom:start-end`. `cell_groups` has one label per row. `peaks` maps each group name to a non-empty list of names taken from `peak_names`. `groupby` names the `obs` column. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `cells`, `peaks` and `groupby`; `data["selected_peaks"]` counts the peaks listed in `peaks`. The heatmap trace has the groups as x, the selected peaks as y and the normalised values as z.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_regions_png

out = Path(tempfile.mkdtemp())
matrix = [[5, 0, 1], [6, 1, 0], [0, 7, 2], [1, 8, 1]]
peak_names = ["chr1:100-200", "chr1:300-400", "chr2:10-90"]
cell_groups = ["A", "A", "B", "B"]
peaks = {"A": ["chr1:100-200"], "B": ["chr1:300-400", "chr2:10-90"]}
result = plot_snapatac2_regions_png(matrix, peak_names, cell_groups, peaks, out / "regions.png", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
trace = result["data"]["traces"][0]
print(trace["type"], trace["x"], trace["y"])
print(len(trace["z"]), len(trace["z"][0]), result["data"]["selected_peaks"])
```

- **Example output:**

```text
completed regions {'cells': 4, 'peaks': 3, 'groupby': 'group'}
True True
heatmap ['A', 'B'] ['chr1:100-200', 'chr1:300-400', 'chr2:10-90']
3 2 3
```

## 22) `dnakit.epigenome.plot_snapatac2_umap_png()`

- **Purpose:** Plot a two- or three-dimensional embedding such as UMAP, coloured by group or value.
- **Method:** Calls `snapatac2.plotting.umap(embedding, color=...)` with `show=False`. In the source, a three-column embedding without `color` is given the single label `all cells` for every point, because snapatac2 2.10 fails for an uncoloured 3D embedding. snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_umap_png(embedding: 'Sequence[Sequence[float]]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, marker_size: 'float | None' = None, marker_opacity: 'float' = 1.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `embedding` is a matrix of finite numbers with two or three columns. `color` has one value per row. `marker_size`, if given, must be positive; `marker_opacity` is between 0 and 1. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `points` and `dimensions`. A 2D embedding has the axis titles `UMAP-1` and `UMAP-2`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_umap_png

out = Path(tempfile.mkdtemp())
embedding = [[0, 0], [1, 1], [2, 0], [3, 1]]
result = plot_snapatac2_umap_png(embedding, out / "umap.png", color=["a", "a", "b", "b"], verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print([(trace["name"], trace["x"], trace["y"]) for trace in result["data"]["traces"]])
print(result["data"]["xaxis_title"], result["data"]["yaxis_title"])
```

- **Example output:**

```text
completed umap {'points': 4, 'dimensions': 2}
True True
[('a', [0.0, 1.0], [0.0, 1.0]), ('b', [2.0, 3.0], [0.0, 1.0])]
UMAP-1 UMAP-2
```

## 23) `dnakit.epigenome.plot_snapatac2_motif_enrichment_png()`

- **Purpose:** Draw motif enrichment scores of several groups as a clustered heatmap.
- **Method:** Turns each group's columns into a Polars `DataFrame` and calls `snapatac2.plotting.motif_enrichment(tables, min_log_fc=..., max_fdr=...)`, which keeps the motifs whose absolute log2 fold change reaches `min_log_fc` in at least one group and whose adjusted p-value is at most `max_fdr` in at least one group (not necessarily the same group), and clusters them. snapatac2 builds a Plotly figure; DNAKit draws the plotted figure again with Matplotlib (Plotly's own PNG export needs Chrome) and writes the PNG.
- **API:** `dnakit.epigenome.plot_snapatac2_motif_enrichment_png(enrichment: 'Mapping[str, Mapping[str, Sequence[Any]]]', output_png: 'str', *, min_log_fc: 'float' = 1.0, max_fdr: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed. Polars comes with snapatac2.
- **Input:** `enrichment` maps at least two group names to column mappings with the columns `id`, `log2(fold change)`, `adjusted p-value` and `p-value` of equal length; every group must list the same motif ids in the same order. `min_log_fc` is not negative and `max_fdr` is in (0, 1]. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `image_bytes`, `parameters` and `data`. `data` holds the plotted `traces` (type, name and the x/y/z values; an array above 20000 values is replaced by `{"omitted": true, "length": n}`), the figure `title`, the axis titles and types (`xaxis_title`, `xaxis_type`, `yaxis_title`, `yaxis_type`) and, when the figure has them, the tick texts. `parameters` holds `groups`, `motifs`, `min_log_fc` and `max_fdr`. The tick texts are the group names (x) and the motifs that were kept (y).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_snapatac2_motif_enrichment_png

out = Path(tempfile.mkdtemp())
enrichment = {
    "A": {
        "id": ["M1", "M2", "M3"],
        "log2(fold change)": [2.0, 0.1, 1.5],
        "adjusted p-value": [0.001, 0.5, 0.0001],
        "p-value": [0.0001, 0.4, 0.00001],
    },
    "B": {
        "id": ["M1", "M2", "M3"],
        "log2(fold change)": [-1.0, 0.2, 2.5],
        "adjusted p-value": [0.001, 0.5, 0.0001],
        "p-value": [0.0001, 0.4, 0.00001],
    },
}
result = plot_snapatac2_motif_enrichment_png(enrichment, out / "motifs.png", verbose=False)
print(result["status"], result["operation"], result["parameters"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print(result["data"]["xaxis_ticktext"], result["data"]["yaxis_ticktext"])
```

- **Example output:**

```text
completed motif_enrichment {'groups': 2, 'motifs': 3, 'min_log_fc': 1.0, 'max_fdr': 0.01}
True True
['A', 'B'] ['M1', 'M3']
```

## 24) `dnakit.epigenome.plot_snapatac2_coverage_png()`

- **Purpose:** Draw fragment coverage tracks of grouped cells across one region of an h5ad file.
- **Method:** Opens the file read-only with `snapatac2.read(..., backed="r")`, takes the labels from `obs[groupby]`, computes the per-group coverage with snapatac2's internal `_snapatac2.get_coverage(adata, region, groups)` and draws one filled track per group with Matplotlib (y axis label `RPM`). This function does not call the public snapatac2 plotting function, and an internal function may change between snapatac2 versions.
- **API:** `dnakit.epigenome.plot_snapatac2_coverage_png(input_h5ad: 'str', output_png: 'str', *, region: 'str', groupby: 'str', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package snapatac2 (extra `dnakit[bio-snapatac2]`), plus numpy (extra `dnakit[bio-analysis]`) and Matplotlib (extra `dnakit[bio-plots]`). No internet access or Chrome is needed.
- **Input:** `input_h5ad` is an existing h5ad file that holds fragments (for example from `snapatac2.pp.import_fragments`). `region` is `chrom:start-end` with start below end. `groupby` is an `obs` column. `output_png` must end in `.png`, must not exist yet, and its folder must exist; otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `operation` (`coverage`), `image_path`, `image_bytes`, `parameters` (`input_h5ad` as an absolute path, `region`, `groupby`) and `data`: `region`, `groupby`, `positions` (the number of positions in the region) and `tracks`, a list sorted by group with `group`, `max` and `signal` (the coverage values, or an `omitted` marker above 20000 values).
- **Example code:**

```python
import tempfile
from pathlib import Path

import snapatac2 as snap

from dnakit.epigenome import plot_snapatac2_coverage_png

work = Path(tempfile.mkdtemp())
lines = []
for i in range(40):
    start = 1000 + (i * 37) % 800
    cell = "cellA" if i % 2 == 0 else "cellB"
    lines.append(f"chr1\t{start}\t{start + 100}\t{cell}\t1")
(work / "fragments.tsv").write_text("\n".join(lines) + "\n")
adata = snap.pp.import_fragments(
    str(work / "fragments.tsv"),
    chrom_sizes={"chr1": 5000},
    file=str(work / "cells.h5ad"),
    sorted_by_barcode=False,
    min_num_fragments=1,
)
adata.obs["cluster"] = ["c1" if name == "cellA" else "c2" for name in adata.obs_names]
adata.close()

result = plot_snapatac2_coverage_png(
    str(work / "cells.h5ad"), str(work / "coverage.png"), region="chr1:900-2000", groupby="cluster", verbose=False
)
print(result["status"], result["operation"])
print(Path(result["image_path"]).read_bytes()[:4] == b"\x89PNG", result["image_bytes"] > 0)
print(result["data"]["region"], result["data"]["groupby"], result["data"]["positions"])
for track in result["data"]["tracks"]:
    print(track["group"], track["max"], len(track["signal"]))
```

- **Example output:**

```text
completed coverage
True True
chr1:900-2000 cluster 1100
c1 2000.0 1100
c2 2000.0 1100
```
