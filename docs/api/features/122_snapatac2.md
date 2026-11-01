# SnapATAC2：数据集与单细胞染色质绘图

这一组函数是 `dnakit.epigenome` 中基于 [SnapATAC2](https://github.com/scverse/SnapATAC2)（单细胞染色质可及性工具包，包元数据中的许可为 MIT）的适配器。其中九个函数通过 snapatac2 注册表下载示例数据集、motif 集合和基因组文件；十五个函数把数据或 Plotly 风格的图转成 PNG（热图、密度等高线、散点图、二维/三维嵌入、网络得分、质控曲线、peak 可及性、motif 富集和覆盖度轨迹）。它们只负责获取或绘制数据：这里没有任何函数做 peak calling、聚类等分析。

导入方式为 `from dnakit.epigenome import ...`。依赖都是可选的，只有调用对应函数时才会导入；缺少依赖时函数会抛出带安装命令的 `BackendUnavailableError`。绘图函数用 Matplotlib 生成 PNG，因为 Plotly 自带的 PNG 导出需要安装 Chrome。下载函数需要联网，所以它们的示例展示无效参数或缺少依赖时抛出的错误；真实调用以不执行的用法代码块给出。

## 示例数据

不需要示例文件。每个示例直接传入小型输入，或把一个很小的 fragment 文件写入临时目录。

## 1) `dnakit.epigenome.load_snapatac2_dataset()`

- **作用：** 把 snapatac2 注册表数据集中的一个文件下载到本地缓存。
- **处理方法：** 检查 `dataset`（以及给出的 `file`）是否为非空字符串，然后调用 `snapatac2.datasets.load(dataset[, file=...])`，由 pooch 下载。结果若是路径或路径列表，则放入 `paths`；否则在 `data` 中返回对象类型，列表还返回其长度。未知的数据集或文件名（`KeyError`）和被拒绝的请求（`TypeError`）转为 `ConfigurationError`；`OSError`、`ValueError`、`RuntimeError` 转为 `DownloadError`，并附带关于 `SNAP_DATA_DIR` 的提示。
- **API：** `dnakit.epigenome.load_snapatac2_dataset(dataset: 'str', *, file: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `dataset` 是注册表名称：`HEA_cCRE`、`Meuleman_2020`、`atac_pbmc_500`、`atac_pbmc_5k`、`cisBP_human`、`colon_transverse`、`gencode_GRCh37`、`gencode_GRCh38`、`gencode_GRCm38`、`gencode_GRCm39` 或 `pbmc_10k_multiome`。`file` 可选，指定该数据集中的某个文件。参数无效时在下载前抛出 `ConfigurationError`；注册表中不存在的名称同样抛出 `ConfigurationError`，网络或文件失败抛出 `DownloadError`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`parameters` 为 `{"file": file}`；除非 snapatac2 返回了非路径对象，否则 `data` 为空。
- **示例代码：**

```python
from dnakit.epigenome import load_snapatac2_dataset
from dnakit.exceptions import ConfigurationError

try:
    load_snapatac2_dataset("no_such_dataset", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error.args[0]).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError snapatac2 does not know the dataset or file 'no_such_dataset'
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import load_snapatac2_dataset

result = load_snapatac2_dataset("atac_pbmc_500", file="atac_pbmc_500_downsample.tsv.gz")
print(result["paths"])
```

## 2) `dnakit.epigenome.fetch_snapatac2_pbmc500()`

- **作用：** 下载 10x 500 细胞 PBMC scATAC-seq 示例，格式可选 fragment、BAM 或 FASTQ。
- **处理方法：** 校验参数，然后调用 `snapatac2.datasets.pbmc500(type=data_type, downsample=downsample)`，以数据集名 `atac_pbmc_500` 返回本地路径。
- **API：** `dnakit.epigenome.fetch_snapatac2_pbmc500(*, data_type: 'str' = 'fragment', downsample: bool = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `data_type` 为 `fragment`、`bam` 或 `fastq`。`downsample` 必须是布尔值，且只有在 `data_type="fragment"` 时才能为真。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`parameters` 含 `data_type` 和 `downsample`；`data` 为空。
- **示例代码：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc500
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc500(data_type="bam", downsample=True, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **示例输出：**

```text
ConfigurationError downsample is only available with data_type='fragment'.
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc500

result = fetch_snapatac2_pbmc500(data_type="fragment", downsample=True)
print(result["dataset"], result["paths"])
```

## 3) `dnakit.epigenome.fetch_snapatac2_pbmc5k()`

- **作用：** 下载 10x 5k PBMC scATAC-seq 示例，格式可选 fragment、h5ad 或带注释的 h5ad。
- **处理方法：** 校验 `data_type`，然后调用 `snapatac2.datasets.pbmc5k(type=data_type)`，以数据集名 `atac_pbmc_5k` 返回本地路径。
- **API：** `dnakit.epigenome.fetch_snapatac2_pbmc5k(*, data_type: 'str' = 'fragment', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `data_type` 为 `fragment`、`h5ad` 或 `annotated_h5ad`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`parameters` 含 `data_type`；`data` 为空。
- **示例代码：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc5k
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc5k(data_type="csv", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **示例输出：**

```text
ConfigurationError data_type must be one of: fragment, h5ad, annotated_h5ad.
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc5k

result = fetch_snapatac2_pbmc5k(data_type="h5ad")
print(result["dataset"], result["paths"])
```

## 4) `dnakit.epigenome.fetch_snapatac2_pbmc10k_multiome()`

- **作用：** 下载 10x 10k PBMC multiome 示例中的 RNA 或 ATAC 部分。
- **处理方法：** 校验两个选项，然后调用 `snapatac2.datasets.pbmc10k_multiome(modality=..., type=...)`，以数据集名 `pbmc_10k_multiome` 返回本地路径。
- **API：** `dnakit.epigenome.fetch_snapatac2_pbmc10k_multiome(*, modality: 'str' = 'RNA', data_type: 'str' = 'h5ad', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `modality` 为 `RNA` 或 `ATAC`；`data_type` 为 `fragment` 或 `h5ad`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`parameters` 含 `modality` 和 `data_type`；`data` 为空。
- **示例代码：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc10k_multiome
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_pbmc10k_multiome(modality="ADT", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **示例输出：**

```text
ConfigurationError modality must be one of: RNA, ATAC.
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_pbmc10k_multiome

result = fetch_snapatac2_pbmc10k_multiome(modality="ATAC", data_type="h5ad")
print(result["dataset"], result["paths"])
```

## 5) `dnakit.epigenome.fetch_snapatac2_colon_transverse()`

- **作用：** 下载并解压横结肠 scATAC-seq fragment 文件，每个样本一个。
- **处理方法：** 调用 `snapatac2.datasets.colon()`，它给出（样本名，fragment 文件）对，并将其列出。
- **API：** `dnakit.epigenome.fetch_snapatac2_colon_transverse(*, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例通过在 `sys.modules` 中放入 `None` 来模拟缺少 `snapatac2.datasets` 模块，因此可离线运行，并展示缺少依赖时抛出的错误。要真正运行，请安装依赖并联网。
- **输入：** 没有数据参数。缺少 snapatac2 时抛出带安装命令的 `BackendUnavailableError`；下载失败抛出 `DownloadError`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`paths` 列出 fragment 文件，`data["samples"]` 列出 `{"sample", "fragment_file"}` 行。数据集名为 `colon_transverse`。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_colon_transverse

result = fetch_snapatac2_colon_transverse()
for row in result["data"]["samples"]:
    print(row["sample"], row["fragment_file"])
```

## 6) `dnakit.epigenome.fetch_snapatac2_hea_ccre()`

- **作用：** 下载整理好的人类结肠顺式调控元件 BED 文件。
- **处理方法：** 调用 `snapatac2.datasets.cre_HEA()`，以数据集名 `HEA_cCRE` 返回本地路径。
- **API：** `dnakit.epigenome.fetch_snapatac2_hea_ccre(*, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例通过在 `sys.modules` 中放入 `None` 来模拟缺少 `snapatac2.datasets` 模块，因此可离线运行，并展示缺少依赖时抛出的错误。要真正运行，请安装依赖并联网。
- **输入：** 没有数据参数。缺少 snapatac2 时抛出带安装命令的 `BackendUnavailableError`；下载失败抛出 `DownloadError`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。`paths` 为 BED 文件；`data` 为空。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_hea_ccre

result = fetch_snapatac2_hea_ccre()
print(result["dataset"], result["paths"])
```

## 7) `dnakit.epigenome.fetch_snapatac2_cis_bp_motifs()`

- **作用：** 下载 CIS-BP 转录因子 motif，并返回它们的名称和标识符。
- **处理方法：** 调用 `snapatac2.datasets.cis_bp(unique=unique)`，用 `id` 和 `name` 概括每个返回的 motif。不返回 motif 对象本身。
- **API：** `dnakit.epigenome.fetch_snapatac2_cis_bp_motifs(*, unique: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `unique` 必须是布尔值（此时 snapatac2 每个转录因子只保留一个 motif）。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。数据集名为 `cisBP_human`；`paths` 为空；`data` 含 `registry_file`（`cisBP_human.meme`）、`motif_count` 和 `motifs`（`{"id", "name"}` 列表）。
- **示例代码：**

```python
from dnakit.epigenome import fetch_snapatac2_cis_bp_motifs
from dnakit.exceptions import ConfigurationError

try:
    fetch_snapatac2_cis_bp_motifs(unique="yes", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.args[0])
```

- **示例输出：**

```text
ConfigurationError unique must be a boolean.
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_cis_bp_motifs

result = fetch_snapatac2_cis_bp_motifs(unique=True)
print(result["data"]["motif_count"], result["data"]["motifs"][0])
```

## 8) `dnakit.epigenome.fetch_snapatac2_meuleman_2020_motifs()`

- **作用：** 下载 Meuleman 2020 分组 motif，并返回它们的名称、标识符和家族。
- **处理方法：** 调用 `snapatac2.datasets.Meuleman_2020()`，用 `id`、`name` 以及（若有）`family` 概括每个 motif。
- **API：** `dnakit.epigenome.fetch_snapatac2_meuleman_2020_motifs(*, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例通过在 `sys.modules` 中放入 `None` 来模拟缺少 `snapatac2.datasets` 模块，因此可离线运行，并展示缺少依赖时抛出的错误。要真正运行，请安装依赖并联网。
- **输入：** 没有数据参数。缺少 snapatac2 时抛出带安装命令的 `BackendUnavailableError`；下载失败抛出 `DownloadError`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。数据集名为 `Meuleman_2020`；`paths` 为空；`data` 含 `registry_file`（`Meuleman_2020.meme`）、`motif_count` 和 `motifs`。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError snapatac2.datasets could not be imported
Install the optional dependency with: python -m pip install "dnakit[bio-snapatac2]"
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_meuleman_2020_motifs

result = fetch_snapatac2_meuleman_2020_motifs()
print(result["data"]["motif_count"], result["data"]["motifs"][0])
```

## 9) `dnakit.epigenome.fetch_snapatac2_genome_files()`

- **作用：** 下载某个基因组预设的 FASTA 和/或基因注释。
- **处理方法：** 把基因组名称或别名映射到 GENCODE 注册表数据集，然后调用 `snapatac2.genome.load(dataset, file=...)` 获取 FASTA（用 `pooch.Decompress` 做 gzip 解压）和/或 GFF3 注释（保持下载时的样子）。基因组和 `kind` 在导入 snapatac2 之前检查。
- **API：** `dnakit.epigenome.fetch_snapatac2_genome_files(genome: 'str', *, kind: 'str' = 'both', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`，pooch 随它一起安装），首次下载需要联网。snapatac2 把文件保存在缓存目录，可用环境变量 `SNAP_DATA_DIR` 修改。不需要密钥。 示例不需要下载：它展示无效参数触发的错误，该检查发生在调用 snapatac2 之前。要真正运行，请安装依赖并联网。
- **输入：** `genome` 为 `GRCh37`、`GRCh38`、`GRCm39`、`GRCm38` 或别名（`hg19`、`hg38`、`mm39`、`mm10`）。`kind` 为 `fasta`、`annotation` 或 `both`。
- **输出：** 字典，含 `status`（`"completed"`）、`service`、`dataset`、`parameters`、`paths`（本地文件路径）、`retrieved_utc`（调用时的 UTC 时间）和 `data`。数据集名为 `gencode_<genome>`；`parameters` 含 `genome`（规范名称）、`requested_name` 和 `kind`；`data` 把 `fasta` 和/或 `annotation` 映射到本地路径。
- **示例代码：**

```python
from dnakit.epigenome import fetch_snapatac2_genome_files
from dnakit.exceptions import ConfigurationError

for arguments in ({"genome": "hg00"}, {"genome": "hg38", "kind": "bed"}):
    try:
        fetch_snapatac2_genome_files(verbose=False, **arguments)
    except ConfigurationError as error:
        print(type(error).__name__, error.args[0])
```

- **示例输出：**

```text
ConfigurationError genome must be one of: GRCh37, GRCh38, GRCm38, GRCm39, hg19, hg38, mm10, mm39.
ConfigurationError kind must be one of: fasta, annotation, both.
```

- **安装依赖后的用法：**

```python
from dnakit.epigenome import fetch_snapatac2_genome_files

result = fetch_snapatac2_genome_files("hg38", kind="annotation")
print(result["parameters"]["genome"], result["data"]["annotation"])
```

## 10) `dnakit.epigenome.render_snapatac2_figure_png()`

- **作用：** 用 snapatac2 的版式规则渲染 Plotly 图字典，并写入一个新的 PNG 文件。
- **处理方法：** 用 `figure["data"]` 和 `figure["layout"]` 构建 `plotly.graph_objects.Figure`，交给 `snapatac2.plotting.render_plot(width=..., height=..., show=False)`，再用 Matplotlib 绘制结果。只能绘制 scatter、bar、violin、heatmap、contour 和 scatter3d 轨迹，其他类型抛出 `ConfigurationError`。
- **API：** `dnakit.epigenome.render_snapatac2_figure_png(figure: 'Mapping[str, Any]', output_png: 'str', *, width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `figure` 是含 `data` 列表和可选 `layout` 的映射。`width` 和 `height` 为正整数。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `width` 和 `height`；`operation` 为 `render_plot`。
- **示例代码：**

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

- **示例输出：**

```text
completed render_plot {'width': 500, 'height': 300}
True True
Demo Group
[{'type': 'bar', 'name': None, 'x': ['a', 'b', 'c'], 'y': [1, 3, 2]}]
```

## 11) `dnakit.epigenome.plot_snapatac2_heatmap_png()`

- **作用：** 把数值矩阵画成热图，可选层次聚类树状图。
- **处理方法：** 调用 `snapatac2.plotting.heatmap(...)`（用 SciPy 聚类），snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。行或列不聚类时，把给定名称设为刻度标签。
- **API：** `dnakit.epigenome.plot_snapatac2_heatmap_png(matrix: 'Sequence[Sequence[float]]', output_png: 'str', *, row_names: 'Sequence[str] | None' = None, column_names: 'Sequence[str] | None' = None, cluster_rows: 'bool' = True, cluster_columns: 'bool' = True, colorscale: 'str' = 'Blues', linkage: 'str' = 'ward', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `matrix` 是由有限数值组成的非空二维矩阵。`row_names` 和 `column_names` 需要每行或每列一个名称。聚类至少需要两行或两列。`linkage` 为 `single`、`complete`、`average`、`weighted`、`centroid`、`median`、`ward` 之一。`colorscale` 是 Plotly 或 Matplotlib 的色图名称（未知名称按 viridis 绘制）。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `rows`、`columns`、`cluster_rows`、`cluster_columns`、`colorscale` 和 `linkage`。聚类时刻度文字按聚类后的顺序排列，树状图是额外的 scatter 轨迹。
- **示例代码：**

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

- **示例输出：**

```text
completed heatmap {'rows': 3, 'columns': 3, 'cluster_rows': True, 'cluster_columns': True, 'colorscale': 'Blues', 'linkage': 'ward'}
True True
['scatter', 'scatter', 'scatter', 'scatter', 'heatmap']
['c3', 'c1', 'c2'] ['r3', 'r1', 'r2']
```

## 12) `dnakit.epigenome.plot_snapatac2_kde2d_png()`

- **作用：** 为成对观测画二维核密度等高线图。
- **处理方法：** 调用 `snapatac2.plotting.kde2d(x, y, log_x=..., log_y=...)`，它返回密度等高线图，snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_kde2d_png(x: 'Sequence[float]', y: 'Sequence[float]', output_png: 'str', *, log_x: 'bool' = False, log_y: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `x` 和 `y` 是等长、至少三个有限数值的序列；使用 `log_x` 或 `log_y` 时，对应序列必须为正。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `points`、`log_x` 和 `log_y`；轨迹为 `contour`。
- **示例代码：**

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

- **示例输出：**

```text
completed kde2d {'points': 8, 'log_x': True, 'log_y': False}
True True
['contour'] log
```

## 13) `dnakit.epigenome.plot_snapatac2_scatter_png()`

- **作用：** 画二维散点图，可按类别或数值着色。
- **处理方法：** 以 `show=False` 调用 `snapatac2.plotting.scatter(...)`。文本颜色会按类别各生成一条轨迹；全为数字的颜色会生成颜色条。snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_scatter_png(x: 'Sequence[float]', y: 'Sequence[float]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, x_label: 'str' = 'X', y_label: 'str' = 'Y', color_label: 'str' = 'Color', marker_size: 'float' = 2.0, marker_opacity: 'float' = 0.5, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `x` 和 `y` 是等长的有限数值序列。`color` 每个点一个值（全为数字或全为文本）。`x_label`、`y_label` 和 `color_label` 为非空字符串。`marker_size` 必须为正，`marker_opacity` 在 0 到 1 之间。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `points`、`marker_size` 和 `marker_opacity`。
- **示例代码：**

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

- **示例输出：**

```text
completed scatter {'points': 4, 'marker_size': 2.0, 'marker_opacity': 0.5}
True True
[('a', [1.0, 2.0], [4.0, 3.0]), ('b', [3.0, 4.0], [2.0, 1.0])]
PC1 PC2
```

## 14) `dnakit.epigenome.plot_snapatac2_scatter3d_png()`

- **作用：** 画三维散点图，可按类别或数值着色。
- **处理方法：** 以 `show=False` 调用 `snapatac2.plotting.scatter3d(...)`，并在 Matplotlib 三维坐标轴上绘制这些点。三维绘制时每条轨迹使用一种颜色，因此数值颜色保留在返回的图数据中，但 PNG 中没有颜色条。
- **API：** `dnakit.epigenome.plot_snapatac2_scatter3d_png(x: 'Sequence[float]', y: 'Sequence[float]', z: 'Sequence[float]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, x_label: 'str' = 'X', y_label: 'str' = 'Y', z_label: 'str' = 'Z', color_label: 'str' = 'Color', marker_size: 'float' = 2.0, marker_opacity: 'float' = 0.5, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `x`、`y` 和 `z` 是等长的有限数值序列。`color`、各标签、`marker_size` 和 `marker_opacity` 的规则与二维散点图相同。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `points`、`marker_size` 和 `marker_opacity`；轨迹类型为 `scatter3d`。
- **示例代码：**

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

- **示例输出：**

```text
completed scatter3d {'points': 4, 'marker_size': 2.0, 'marker_opacity': 0.5}
True True
scatter3d [1.0, 2.0, 3.0, 4.0] [4.0, 3.0, 2.0, 1.0] [1.0, 1.0, 2.0, 2.0]
```

## 15) `dnakit.epigenome.plot_snapatac2_network_scores_png()`

- **作用：** 按到 TSS 的距离区间，绘制调控网络边得分的平均值。
- **处理方法：** 用节点类型和边构建 `rustworkx.PyDiGraph`，调用 `snapatac2.plotting.network_scores(graph, score_name=...)`，再用 Matplotlib 画出柱状图。
- **API：** `dnakit.epigenome.plot_snapatac2_network_scores_png(node_types: 'Sequence[str]', edges: 'Sequence[Mapping[str, Any]]', output_png: 'str', *, score_name: 'str', width: 'int' = 800, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。rustworkx 随 snapatac2 一起安装。
- **输入：** `node_types` 是非空字符串列表。`edges` 是非空的映射列表，含整数 `source` 和 `target` 节点索引，以及可选的数值 `distance`、`cor_score`、`regr_score`。每条边都需要 `distance`，且至少一条边有所选得分。`score_name` 为 `cor_score` 或 `regr_score`。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `score_name` 和 `edges`；轨迹为柱状图，x 是距离区间标签，y 是平均得分。
- **示例代码：**

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

- **示例输出：**

```text
completed network_scores {'score_name': 'cor_score', 'edges': 4}
True True
bar ['500 - 2K', '2K - 20K', '20K - 50K'] [0.5, 0.3, -0.05]
Distance to TSS (bp) | Average score
```

## 16) `dnakit.epigenome.plot_snapatac2_network_edge_stat_png()`

- **作用：** 为每一种源节点与目标节点类型组合绘制相关性得分的小提琴图。
- **处理方法：** 与网络得分图一样构建 `rustworkx.PyDiGraph`，调用 `snapatac2.plotting.network_edge_stat(graph)`，并用 Matplotlib 为每个 `source -> target` 类型组合画一个小提琴。只有一个值或没有离散度的组画成一个点。
- **API：** `dnakit.epigenome.plot_snapatac2_network_edge_stat_png(node_types: 'Sequence[str]', edges: 'Sequence[Mapping[str, Any]]', output_png: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。rustworkx 随 snapatac2 一起安装。
- **输入：** `node_types` 和 `edges` 的规则与网络得分图相同，但不要求 `distance`；至少一条边需要 `cor_score`。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `edges`；每个小提琴轨迹名为 `<源类型> -> <目标类型>`，其 `y` 为相关性得分。
- **示例代码：**

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

- **示例输出：**

```text
completed network_edge_stat {'edges': 4}
True True
violin peak -> gene [0.5, 0.3, -0.2, 0.1]
```

## 17) `dnakit.epigenome.plot_snapatac2_valid_cells_png()`

- **作用：** 在双对数坐标上绘制排序后的条形码计数，用于选择细胞判定阈值。
- **处理方法：** 调用 `snapatac2.plotting.valid_cells(counts, ...)`，它把计数从高到低排序，snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_valid_cells_png(counts: 'Sequence[float]', output_png: 'str', *, width: 'int' = 500, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `counts` 至少有两个为正的有限数值，每个条形码一个，顺序不限。`width` 和 `height` 为正整数。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `barcodes`；轨迹的 x 是排名，y 是排序后的计数。
- **示例代码：**

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

- **示例输出：**

```text
completed valid_cells {'barcodes': 6}
True True
[1.0, 2.0, 3.0, 4.0, 5.0, 6.0] [1000.0, 800.0, 500.0, 100.0, 50.0, 10.0]
log log
```

## 18) `dnakit.epigenome.plot_snapatac2_tsse_png()`

- **作用：** 绘制 TSS 富集分数与每个细胞唯一 fragment 数量的关系。
- **处理方法：** 把两个序列放入内存中 `AnnData` 的 `obs` 表，调用 `snapatac2.plotting.tsse(adata, min_fragment=...)`（密度等高线图），snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_tsse_png(n_fragment: 'Sequence[float]', tsse: 'Sequence[float]', output_png: 'str', *, min_fragment: 'int' = 500, width: 'int' = 500, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `n_fragment` 和 `tsse` 是等长、至少两个有限数值的序列，每个细胞一个。`min_fragment` 为非负整数。`width` 和 `height` 为正整数。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `cells` 和 `min_fragment`；`data["cells_plotted"]` 是 `n_fragment >= min_fragment` 的细胞数。
- **示例代码：**

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

- **示例输出：**

```text
completed tsse {'cells': 5, 'min_fragment': 500}
True True
4 ['contour']
Number of unique fragments | TSS enrichment score
```

## 19) `dnakit.epigenome.plot_snapatac2_frag_size_distr_png()`

- **作用：** 根据各 fragment 长度的计数绘制 fragment 大小分布。
- **处理方法：** 调用 `snapatac2.plotting.frag_size_distr(counts, max_recorded_size=...)`（计数对 fragment 大小的折线图），snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。在下面的示例中，索引 0 的值不被绘制，大小从 1 到 `max_recorded_size`。
- **API：** `dnakit.epigenome.plot_snapatac2_frag_size_distr_png(distribution: 'Sequence[float]', output_png: 'str', *, max_recorded_size: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `distribution` 至少有两个非负有限计数，索引即 fragment 长度。`max_recorded_size` 为正整数。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `max_recorded_size` 和 `sizes`（给出的计数个数）。
- **示例代码：**

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

- **示例输出：**

```text
completed frag_size_distr {'max_recorded_size': 9, 'sizes': 10}
True True
[1, 2, 3, 4, 5, 6, 7, 8, 9]
[0.0, 5.0, 10.0, 30.0, 20.0, 10.0, 5.0, 2.0, 1.0]
```

## 20) `dnakit.epigenome.plot_snapatac2_spectral_eigenvalues_png()`

- **作用：** 绘制谱嵌入的特征值以及 snapatac2 选出的拐点。
- **处理方法：** 把数值存入空 `AnnData` 的 `adata.uns["spectral_eigenvalue"]`，调用 `snapatac2.plotting.spectral_eigenvalues(adata)`，并读取 snapatac2 写入 `adata.uns["num_eigen"]` 的拐点。snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_spectral_eigenvalues_png(eigenvalues: 'Sequence[float]', output_png: 'str', *, width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `eigenvalues` 至少有两个有限数值。`width` 和 `height` 为正整数。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `eigenvalues`（个数）；`data["num_eigen"]` 是拐点位置，snapatac2 找不到拐点时为 `None`（例如平滑下降的序列）。
- **示例代码：**

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

- **示例输出：**

```text
completed spectral_eigenvalues {'eigenvalues': 12}
True True
elbow: 9
Component | Eigenvalue
```

## 21) `dnakit.epigenome.plot_snapatac2_regions_png()`

- **作用：** 绘制所选 peak 按组、经 RPKM 归一化的可及性。
- **处理方法：** 用细胞 × peak 矩阵构建内存中的 `AnnData`（peak 名为变量，分组标签在 `obs[groupby]`），调用 `snapatac2.plotting.regions(adata, groupby=..., peaks=...)`，它返回组与所选 peak 的热图。snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_regions_png(matrix: 'Sequence[Sequence[float]]', peak_names: 'Sequence[str]', cell_groups: 'Sequence[str]', peaks: 'Mapping[str, Sequence[str]]', output_png: 'str', *, groupby: 'str' = 'group', width: 'int' = 600, height: 'int' = 400, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `matrix` 是由有限数值组成的细胞 × peak 矩阵。`peak_names` 每列一个唯一名称，格式为 `chrom:start-end`。`cell_groups` 每行一个标签。`peaks` 把每个组名映射到取自 `peak_names` 的非空名称列表。`groupby` 是 `obs` 列名。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `cells`、`peaks` 和 `groupby`；`data["selected_peaks"]` 是 `peaks` 中列出的 peak 数。热图轨迹的 x 是组，y 是所选 peak，z 是归一化后的值。
- **示例代码：**

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

- **示例输出：**

```text
completed regions {'cells': 4, 'peaks': 3, 'groupby': 'group'}
True True
heatmap ['A', 'B'] ['chr1:100-200', 'chr1:300-400', 'chr2:10-90']
3 2 3
```

## 22) `dnakit.epigenome.plot_snapatac2_umap_png()`

- **作用：** 绘制二维或三维嵌入（如 UMAP），按分组或数值着色。
- **处理方法：** 以 `show=False` 调用 `snapatac2.plotting.umap(embedding, color=...)`。源码中，没有 `color` 的三列嵌入会让每个点获得同一个标签 `all cells`，因为 snapatac2 2.10 对未着色的三维嵌入会失败。snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_umap_png(embedding: 'Sequence[Sequence[float]]', output_png: 'str', *, color: 'Sequence[Any] | None' = None, marker_size: 'float | None' = None, marker_opacity: 'float' = 1.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `embedding` 是两列或三列的有限数值矩阵。`color` 每行一个值。`marker_size` 若给出必须为正；`marker_opacity` 在 0 到 1 之间。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `points` 和 `dimensions`。二维嵌入的坐标轴标题为 `UMAP-1` 和 `UMAP-2`。
- **示例代码：**

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

- **示例输出：**

```text
completed umap {'points': 4, 'dimensions': 2}
True True
[('a', [0.0, 1.0], [0.0, 1.0]), ('b', [2.0, 3.0], [0.0, 1.0])]
UMAP-1 UMAP-2
```

## 23) `dnakit.epigenome.plot_snapatac2_motif_enrichment_png()`

- **作用：** 把多个组的 motif 富集得分画成聚类热图。
- **处理方法：** 把每组的列转成 Polars `DataFrame`，调用 `snapatac2.plotting.motif_enrichment(tables, min_log_fc=..., max_fdr=...)`，它保留满足以下条件的 motif 并对其聚类：至少在一个组中 log2 倍数变化的绝对值达到 `min_log_fc`，且至少在一个组中校正 p 值不高于 `max_fdr`（不要求是同一个组）。snapatac2 先构建 Plotly 图；DNAKit 再用 Matplotlib 重新绘制该图（Plotly 自带的 PNG 导出需要 Chrome）并写出 PNG。
- **API：** `dnakit.epigenome.plot_snapatac2_motif_enrichment_png(enrichment: 'Mapping[str, Mapping[str, Sequence[Any]]]', output_png: 'str', *, min_log_fc: 'float' = 1.0, max_fdr: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。Polars 随 snapatac2 一起安装。
- **输入：** `enrichment` 把至少两个组名映射到列映射，列为等长的 `id`、`log2(fold change)`、`adjusted p-value` 和 `p-value`；各组必须以相同顺序列出相同的 motif id。`min_log_fc` 非负，`max_fdr` 在 (0, 1] 内。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`image_bytes`、`parameters` 和 `data`。`data` 包含绘制的 `traces`（类型、名称以及 x/y/z 数值；超过 20000 个值的数组会被 `{"omitted": true, "length": n}` 取代）、图标题、坐标轴标题与类型（`xaxis_title`、`xaxis_type`、`yaxis_title`、`yaxis_type`），图中有刻度文字时还有刻度文字。`parameters` 含 `groups`、`motifs`、`min_log_fc` 和 `max_fdr`。刻度文字为组名（x）和被保留的 motif（y）。
- **示例代码：**

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

- **示例输出：**

```text
completed motif_enrichment {'groups': 2, 'motifs': 3, 'min_log_fc': 1.0, 'max_fdr': 0.01}
True True
['A', 'B'] ['M1', 'M3']
```

## 24) `dnakit.epigenome.plot_snapatac2_coverage_png()`

- **作用：** 在 h5ad 文件的某个区域上，为分组细胞绘制 fragment 覆盖度轨迹。
- **处理方法：** 用 `snapatac2.read(..., backed="r")` 以只读方式打开文件，从 `obs[groupby]` 取标签，用 snapatac2 内部的 `_snapatac2.get_coverage(adata, region, groups)` 计算各组覆盖度，并用 Matplotlib 为每组画一条填充轨迹（y 轴标签为 `RPM`）。该函数不调用 snapatac2 的公开绘图函数，内部函数可能随 snapatac2 版本变化。
- **API：** `dnakit.epigenome.plot_snapatac2_coverage_png(input_h5ad: 'str', output_png: 'str', *, region: 'str', groupby: 'str', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 snapatac2（extra `dnakit[bio-snapatac2]`），以及 numpy（extra `dnakit[bio-analysis]`）和 Matplotlib（extra `dnakit[bio-plots]`）。不需要联网，也不需要 Chrome。
- **输入：** `input_h5ad` 是已存在、含 fragment 的 h5ad 文件（例如来自 `snapatac2.pp.import_fragments`）。`region` 为 `chrom:start-end`，且 start 小于 end。`groupby` 是 `obs` 列名。`output_png` 必须以 `.png` 结尾、尚不存在，且其所在目录必须存在，否则抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`（`coverage`）、`image_path`、`image_bytes`、`parameters`（`input_h5ad` 的绝对路径、`region`、`groupby`）和 `data`：`region`、`groupby`、`positions`（区域内位置数）和按组排序的 `tracks` 列表，每项含 `group`、`max` 和 `signal`（覆盖度数值，超过 20000 个值时为 `omitted` 标记）。
- **示例代码：**

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

- **示例输出：**

```text
completed coverage
True True
chr1:900-2000 cluster 1100
c1 2000.0 1100
c2 2000.0 1100
```
