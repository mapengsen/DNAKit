# gReLU：模型库、数据集、绘图与基序工具

这一组函数是 `dnakit.predictions` 中基于 [gReLU](https://github.com/Genentech/gReLU)（MIT 许可）的适配器。gReLU 是 Genentech 开发的序列到功能模型训练与解释库。每个函数都用调用者给出的参数调用一个 gReLU 函数，不重新实现其功能。这一组涵盖：基因组与区间工具（染色体末端检查、GC 匹配的背景区间、序列提取、变异序列、染色体大小、GTF 读取）；[Hugging Face](https://huggingface.co/Genentech) 上的 gReLU 模型库以及旧版 [Weights & Biases](https://wandb.ai/) 项目；JASPAR 基序（[pyjaspar](https://github.com/asntech/pyjaspar)）；tomtom 基序比较（[memelite](https://github.com/jmschrei/memesuite-lite)）；以及预测、归因和基序的绘图。基因组通过 [genomepy](https://github.com/vanheeringen-lab/genomepy) 获取。上述各包均为 MIT 或 Apache-2.0 许可。

用 `from dnakit.predictions import ...` 导入这些函数。依赖是可选的（额外依赖组 `dnakit[bio-grelu]`：`grelu`、`huggingface-hub`、`wandb`、`memelite`），只有调用函数时才会导入。缺少包时，函数抛出带安装提示的 `BackendUnavailableError`。需要访问远程服务的函数（基因组下载、Hugging Face 与 Weights & Biases 模型库）只在被调用时才联网。本页所有示例都在离线环境中运行，使用 gReLU 1.1.0。

以下三点适用于许多函数：

- 凡是 `genome` 参数为基因组名称（例如 `hg38`）时，gReLU 会在第一次使用时用 genomepy 安装它，这会下载基因组。如果 `genome` 是本地 FASTA 文件的路径，gReLU 直接使用该文件（需要 `.fai` 索引，pyfaidx 会在首次使用时生成，还需要染色体大小文件 `<文件>.sizes`）。示例使用这样的本地文件，因此可以离线运行。
- gReLU 自己的提示信息（例如 `Keeping 1 intervals`）会输出到标准输出。gReLU 的错误会以 `BackendExecutionError` 重新抛出，并保留上游的消息。
- 绘图函数把图像保存为 `.png`、`.pdf` 或 `.svg`（文件必须尚不存在，所在文件夹必须存在），并返回包含 `status`、`backend`、`version`、`function`、`image_path` 和 `data`（绘图输入的描述）的字典。gReLU 的若干绘图函数（在各自小节中标明）返回的是 plotnine 对象，必须先绘制再保存。DNAKit 直接保存当前 matplotlib 图形而不绘制，所以这些函数保存的图像是空白的（已用 gReLU 1.1.0 和 plotnine 0.15.8 核实）。

## 示例数据

示例使用一个含两条序列（40 bp 和 30 bp）的基因组 `tiny.fa`、它的染色体大小文件 `tiny.fa.sizes`（制表符分隔）、一个小的 GTF 文件 `genes.gtf`（制表符分隔）和一个含两个基序的 MEME 数据库 `db.meme`。其他输入在示例代码中构造。

**示例文件：** `tiny.fa`

```text
>chr1
ACGTTGCAAGGCTTAACCGGTATCGATCGGATCCAGTCAA
>chr2
TTGACCAGTCAGGCATCGATTCGAGGCCTA
```

**示例文件：** `tiny.fa.sizes`

```text
chr1	40
chr2	30
```

**示例文件：** `genes.gtf`

```text
chr1	ex	gene	101	400	.	+	.	gene_id "g1"; gene_name "GENEA";
chr1	ex	exon	101	200	.	+	.	gene_id "g1"; gene_name "GENEA";
chr1	ex	exon	301	400	.	+	.	gene_id "g1"; gene_name "GENEA";
chr2	ex	gene	51	150	.	-	.	gene_id "g2"; gene_name "GENEB";
```

**示例文件：** `db.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GATA_like
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.05 0.85
0.85 0.05 0.05 0.05

MOTIF CACGTG_like
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.05 0.85 0.05 0.05
0.85 0.05 0.05 0.05
0.05 0.85 0.05 0.05
0.05 0.05 0.85 0.05
0.05 0.05 0.05 0.85
0.05 0.05 0.85 0.05
```

## 1) `dnakit.predictions.grelu_check_chrom_ends()`

- **作用：** 检查基因组区间是否超出基因组的染色体末端。
- **处理方法：** 调用 gReLU 的 `grelu.data.preprocess.check_chrom_ends`。它标出起点为负的区间；给出基因组时，还标出终点大于 `read_sizes` 所得染色体大小的区间。不给基因组时只检查负起点，并由 gReLU 发出警告。所有区间都合法时返回空值，否则抛出错误并列出出错区间的行标签。
- **API：** `grelu_check_chrom_ends(data: 'Any', genome: 'str | None' = None, *, verbose: bool = True) -> 'None'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载（需要联网）；带 `.sizes` 文件的本地 FASTA 路径可离线使用。
- **输入：** `data` 是含 `chrom`、`start`、`end` 列的 pandas DataFrame（或在 `.var` 中保存区间的 AnnData 对象）。`genome` 是基因组名称或本地 FASTA 文件路径。
- **输出：** 所有区间都在染色体内时返回 `None`。否则抛出 `BackendExecutionError`，消息中列出出错区间的索引。
- **示例代码：**

```python
import pandas as pd

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_check_chrom_ends

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr2"], "start": [-5, 10, 20], "end": [10, 20, 35]}
)
print(grelu_check_chrom_ends(intervals.iloc[[1]], "tiny.fa", verbose=False))
try:
    grelu_check_chrom_ends(intervals, "tiny.fa", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, str(error).split("failed: ")[1])
```

- **示例输出：**

```text
None
BackendExecutionError Indices of intervals that extend beyond the chromosome ends: 0,2.
```

## 2) `dnakit.predictions.grelu_check_reference()`

- **作用：** 检查变异表中的参考等位基因是否与参考基因组一致。
- **处理方法：** 调用 gReLU 的 `grelu.variant.check_reference`。对每个变异，用 `intervals_to_strings` 取出 `pos`（1 起始）处的基因组碱基，并与 `ref` 比较。`ref` 等于 `null_string` 的行被跳过。不一致时以 Python `UserWarning` 报告并列出行标签，而不是抛出错误。
- **API：** `grelu_check_reference(variants: 'Any', *, genome: 'str' = 'hg38', null_string: 'str' = '-', verbose: bool = True) -> 'None'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；本地 FASTA 路径可离线使用。
- **输入：** `variants` 是含 `chrom`、`pos`、`ref`、`alt` 列的 DataFrame。`genome` 是基因组名称或 FASTA 路径（默认 `hg38`）。`null_string` 表示缺失的碱基（默认 `-`）。
- **输出：** `None`。请通过警告（例如 `warnings.catch_warnings(record=True)`）查看不一致的行。输入表保持不变。
- **示例代码：**

```python
import warnings

import pandas as pd

from dnakit.predictions import grelu_check_reference

variants = pd.DataFrame(
    {"chrom": ["chr1", "chr1"], "pos": [5, 6], "ref": ["T", "A"], "alt": ["G", "C"]}
)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    result = grelu_check_reference(variants, genome="tiny.fa", verbose=False)
messages = [str(w.message) for w in caught if "not centered" in str(w.message)]
print(result, len(messages))
print(messages[0])
```

- **示例输出：**

```text
None 1
Sequences are not centered on reference at 1 indices [1]
```

## 3) `dnakit.predictions.grelu_download_dataset()`

- **作用：** 从 Hugging Face 下载 gReLU 模型库的数据集文件，并返回本地路径。
- **处理方法：** 调用 gReLU 的 `grelu.resources.download_dataset`，它以 `repo_type="dataset"` 调用 `huggingface_hub.hf_hub_download`。文件保存在 Hugging Face 缓存中，返回缓存路径。
- **API：** `grelu_download_dataset(repo_id: 'str', filename: 'str' = 'data.h5ad', *, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接（公开仓库不需要密钥）。设置 `HF_HUB_OFFLINE=1` 后只使用本地缓存。
- **输入：** `repo_id` 是非空的数据集仓库 ID，例如 `Genentech/human-atac-catlas-data`；`filename` 是仓库中的文件（默认 `data.h5ad`）。
- **输出：** 本地文件路径字符串。下载失败时抛出 `BackendExecutionError`。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_download_dataset

try:
    grelu_download_dataset("example-org/missing-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "local cache" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_download_dataset

path = grelu_download_dataset("Genentech/human-atac-catlas-data", "data.h5ad")
```

## 4) `dnakit.predictions.grelu_download_model()`

- **作用：** 从 Hugging Face 下载 gReLU 模型库的模型检查点，并返回本地路径。
- **处理方法：** 调用 gReLU 的 `grelu.resources.download_model`，它调用 `huggingface_hub.hf_hub_download`。文件保存在 Hugging Face 缓存中。
- **API：** `grelu_download_model(repo_id: 'str', filename: 'str' = 'model.ckpt', *, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接（公开仓库不需要密钥）。设置 `HF_HUB_OFFLINE=1` 后只使用本地缓存。
- **输入：** `repo_id` 是非空的模型仓库 ID，例如 `Genentech/human-atac-catlas-model`；`filename` 是检查点文件（默认 `model.ckpt`）。
- **输出：** 本地文件路径字符串。失败时抛出 `BackendExecutionError`。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_download_model

try:
    grelu_download_model("example-org/missing-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "local cache" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_download_model

path = grelu_download_model("Genentech/human-atac-catlas-model")
```

## 5) `dnakit.predictions.grelu_filter_chrom_ends()`

- **作用：** 丢弃超出染色体末端的区间，也可以先从两端各去掉一段边距。
- **处理方法：** 调用 gReLU 的 `grelu.data.preprocess.filter_chrom_ends`。保留 `start >= pad` 的区间；给出基因组时，还要求 `end <= 染色体大小 - pad`。gReLU 会输出 `Keeping N intervals`。
- **API：** `grelu_filter_chrom_ends(data: 'Any', genome: 'str | None' = None, *, pad: 'int' = 0, inplace: 'bool' = False, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；带 `.sizes` 文件的本地 FASTA 路径可离线使用。
- **输入：** `data` 是含 `chrom`、`start`、`end` 的 DataFrame（或在 `.var` 中保存区间的 AnnData 对象）。`genome` 是基因组名称或 FASTA 路径。`pad` 是每个染色体末端要忽略的碱基数。`inplace=True` 时直接修改输入。
- **输出：** 相同格式的过滤后区间（除非 `inplace=True`，否则是新的 DataFrame）。
- **示例代码：**

```python
import pandas as pd

from dnakit.predictions import grelu_filter_chrom_ends

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr1"], "start": [3, 10, 33], "end": [10, 20, 38]}
)
print(len(grelu_filter_chrom_ends(intervals, "tiny.fa", verbose=False)))
kept = grelu_filter_chrom_ends(intervals, "tiny.fa", pad=5, verbose=False)
print(kept.to_dict("records"))
```

- **示例输出：**

```text
Keeping 3 intervals
3
Keeping 1 intervals
[{'chrom': 'chr1', 'start': 10, 'end': 20}]
```

## 6) `dnakit.predictions.grelu_get_base_models()`

- **作用：** 列出模型库中某个模型微调所基于的基础模型。
- **处理方法：** 调用 gReLU 的 `grelu.resources.get_base_models`，用 `HfApi().model_info` 读取模型卡的 `base_model` 字段。单个字符串会变成只有一个元素的列表；字段为空则返回空列表。
- **API：** `grelu_get_base_models(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** `repo_id` 是非空的模型仓库 ID。
- **输出：** 基础模型仓库 ID 的列表。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_base_models

try:
    grelu_get_base_models("Genentech/borzoi-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_get_base_models

bases = grelu_get_base_models("Genentech/human-atac-catlas-model")
```

## 7) `dnakit.predictions.grelu_get_dataset_info()`

- **作用：** 读取 Hugging Face 上某个数据集仓库的数据集卡元数据。
- **处理方法：** 调用 gReLU 的 `grelu.resources.get_dataset_info`，使用 `HfApi().dataset_info` 和 `list_repo_files(repo_type="dataset")`。
- **API：** `grelu_get_dataset_info(repo_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** `repo_id` 是非空的数据集仓库 ID。
- **输出：** 含 `id`、`tags`、`card_data`、`downloads`、`last_modified` 和 `files` 的字典。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_dataset_info

try:
    grelu_get_dataset_info("Genentech/human-atac-catlas-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_get_dataset_info

info = grelu_get_dataset_info("Genentech/human-atac-catlas-data")
print(info["files"])
```

## 8) `dnakit.predictions.grelu_get_datasets_by_model()`

- **作用：** 列出模型卡中与某个模型关联的数据集。
- **处理方法：** 调用 gReLU 的 `grelu.resources.get_datasets_by_model`，返回用 `HfApi().model_info` 读取的模型卡 `datasets` 字段（字段为空则返回空列表）。
- **API：** `grelu_get_datasets_by_model(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** `repo_id` 是非空的模型仓库 ID。
- **输出：** 数据集仓库 ID 的列表。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_datasets_by_model

try:
    grelu_get_datasets_by_model("Genentech/human-atac-catlas-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_get_datasets_by_model

datasets = grelu_get_datasets_by_model("Genentech/human-atac-catlas-model")
```

## 9) `dnakit.predictions.grelu_get_gc_matched_intervals()`

- **作用：** 采样 GC 含量与给定区间相匹配的背景区间。
- **处理方法：** 调用 gReLU 的 `grelu.data.preprocess.get_gc_matched_intervals`，它在基因组 FASTA 上以输入区间的长度调用 tangermeme 的 `extract_matching_loci`。之后 gReLU 去除黑名单区域：`blacklist=None` 时，它按基因组名称查找内置黑名单，找不到（例如基因组是 FASTA 路径）就跳过这一步并输出一条消息。
- **API：** `grelu_get_gc_matched_intervals(intervals: 'Any', genome: 'str', *, binwidth: 'float' = 0.1, chroms: 'str' = 'autosomes', blacklist: 'str | None' = None, seed: 'int | None' = None, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu（含 tangermeme）和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；本地 FASTA 路径可离线使用。
- **输入：** `intervals` 是含 `chrom`、`start`、`end` 且长度一致的 DataFrame。`genome` 是非空的基因组名称或 FASTA 路径。`binwidth` 是 GC 分辨率。`chroms` 是 `autosomes`（chr1 至 chr22）或其他 gReLU 简称；gReLU 也接受染色体名称列表，对没有 chr1 至 chr22 的基因组需要这样用。`blacklist` 是黑名单文件，`seed` 固定采样。
- **输出：** GC 匹配区间的 DataFrame（`chrom`、`start`、`end`）；候选足够时行数与输入区间数相同。
- **示例代码：**

```python
import contextlib
import io
import random

import pandas as pd

from dnakit.predictions import grelu_get_gc_matched_intervals

rng = random.Random(2)
with open("big.fa", "w") as handle:
    for chrom in ("chr1", "chr2"):
        sequence = "".join(rng.choice("ACGT") for _ in range(20000))
        handle.write(f">{chrom}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i:i + 60] + "\n")
with open("big.fa.sizes", "w") as handle:
    handle.write("chr1\t20000\nchr2\t20000\n")

positives = pd.DataFrame(
    {"chrom": ["chr1", "chr1", "chr2"], "start": [1000, 3000, 500], "end": [1100, 3100, 600]}
)
with contextlib.redirect_stdout(io.StringIO()):
    matched = grelu_get_gc_matched_intervals(
        positives, "big.fa", chroms=["chr1", "chr2"], seed=1, verbose=False
    )
print(len(matched), list(matched.columns), (matched.end - matched.start).unique().tolist())
```

- **示例输出：**

```text
3 ['chrom', 'start', 'end'] [100]
```

## 10) `dnakit.predictions.grelu_get_genome()`

- **作用：** 把基因组加载为 gReLU 基因组对象；给出名称时通过 genomepy 安装。
- **处理方法：** 调用 gReLU 的 `grelu.io.genome.get_genome`。如果参数是已存在的文件路径，返回 gReLU 的 `CustomGenome`（pyfaidx 的封装，提供 `get_seq`、`sizes_file` 和 `genome_file`）。否则，若该名称的基因组尚未安装，就用 `genomepy.install_genome` 安装（会下载），并返回 `genomepy.Genome`。
- **API：** `grelu_get_genome(genome: 'str', *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时需要联网；带 `.fai` 索引（pyfaidx 会生成）的本地 FASTA 路径可离线使用。
- **输入：** `genome` 是非空的基因组名称或本地 FASTA 文件路径。
- **输出：** 基因组对象。`get_seq(chrom, start, end)` 使用 1 起始的闭区间坐标。
- **示例代码：**

```python
from dnakit.predictions import grelu_get_genome

genome = grelu_get_genome("tiny.fa", verbose=False)
print(type(genome).__name__, genome.get_seq("chr2", 1, 8))
print(genome.sizes_file.endswith("tiny.fa.sizes"))
```

- **示例输出：**

```text
CustomGenome TTGACCAG
True
```

## 11) `dnakit.predictions.grelu_get_jaspar()`

- **作用：** 从 JASPAR 数据库获取基序，得到位置概率矩阵。
- **处理方法：** 调用 gReLU 的 `grelu.io.motifs.get_jaspar`，它打开 pyjaspar 的 `jaspardb(release)` 并调用 `fetch_motifs`。pyjaspar 以内置 SQLite 文件的形式附带各个 JASPAR 版本，所以测试中不需要网络。键为 `<matrix_id>_<name>`，值是形状为 (4, 长度) 的 NumPy 数组。`species` 接受分类 ID，也接受 `human` 和 `mouse`。注意：gReLU 1.1.0 总是获取 `CORE` 集合，会忽略 `collection` 参数。
- **API：** `grelu_get_jaspar(release: 'str' = 'JASPAR2024', collection: 'str' = 'CORE', *, tax_group: 'str | None' = None, species: 'str | None' = None, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 pyjaspar（grelu 的依赖；额外依赖组 `dnakit[bio-grelu]`）。使用内置版本时无需联网即可运行。
- **输入：** `release` 例如 `JASPAR2024`；`tax_group` 例如 `plants` 或 `vertebrates`；`species` 为分类 ID 字符串，例如 `3702` 表示拟南芥。
- **输出：** 从基序名称到形状为 (4, 长度) 的数组的字典。
- **示例代码：**

```python
from dnakit.predictions import grelu_get_jaspar

motifs = grelu_get_jaspar(tax_group="plants", species="3702", verbose=False)
print(type(motifs).__name__, len(motifs))
name = sorted(motifs)[0]
print(name, motifs[name].shape)
```

- **示例输出：**

```text
dict 591
MA0001.3_AGL3 (4, 10)
```

## 12) `dnakit.predictions.grelu_get_model_info()`

- **作用：** 读取 Hugging Face 上某个模型仓库的模型卡元数据。
- **处理方法：** 调用 gReLU 的 `grelu.resources.get_model_info`，使用 `HfApi().model_info` 和 `list_repo_files`。
- **API：** `grelu_get_model_info(repo_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** `repo_id` 是非空的模型仓库 ID。
- **输出：** 含 `id`、`tags`、`card_data`、`downloads`、`last_modified` 和 `files` 的字典。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_model_info

try:
    grelu_get_model_info("Genentech/human-atac-catlas-model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_get_model_info

info = grelu_get_model_info("Genentech/human-atac-catlas-model")
print(info["files"])
```

## 13) `dnakit.predictions.grelu_get_models_by_dataset()`

- **作用：** 列出 gReLU 集合中在某个数据集上训练的模型。
- **处理方法：** 调用 gReLU 的 `grelu.resources.get_models_by_dataset`。它列出模型库集合中的模型，并逐个读取其模型卡的 `datasets` 字段，所以每个模型发送一次请求。
- **API：** `grelu_get_models_by_dataset(repo_id: 'str', *, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** `repo_id` 是非空的数据集仓库 ID。
- **输出：** 模型仓库 ID 的列表。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_get_models_by_dataset

try:
    grelu_get_models_by_dataset("Genentech/human-atac-catlas-data", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_get_models_by_dataset

models = grelu_get_models_by_dataset("Genentech/human-atac-catlas-data")
```

## 14) `dnakit.predictions.grelu_intervals_to_strings()`

- **作用：** 从基因组中提取基因组区间的 DNA 序列。
- **处理方法：** 调用 gReLU 的 `grelu.sequence.format.intervals_to_strings`。区间为 0 起始的半开区间（`start`、`end`）；序列转为大写，若存在 `strand` 列，则 `strand` 为 `-` 的行取反向互补。
- **API：** `grelu_intervals_to_strings(intervals: 'Any', genome: 'str', *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；本地 FASTA 路径可离线使用。
- **输入：** `intervals` 是含 `chrom`、`start`、`end`（可含 `strand`）的 DataFrame，或表示单个区间的 Series 或字典。`genome` 是非空的基因组名称或 FASTA 路径。
- **输出：** 输入 DataFrame 时返回序列字符串列表；输入 Series 或字典时返回一个字符串。
- **示例代码：**

```python
import pandas as pd

from dnakit.predictions import grelu_intervals_to_strings

intervals = pd.DataFrame(
    {"chrom": ["chr1", "chr2"], "start": [10, 20], "end": [20, 30], "strand": ["+", "-"]}
)
print(grelu_intervals_to_strings(intervals, "tiny.fa", verbose=False))
print(grelu_intervals_to_strings(intervals.iloc[0], "tiny.fa", verbose=False))
```

- **示例输出：**

```text
['GCTTAACCGG', 'TAGGCCTCGA']
GCTTAACCGG
```

## 15) `dnakit.predictions.grelu_list_datasets()`

- **作用：** 列出 Hugging Face 上 gReLU 模型库的数据集仓库 ID。
- **处理方法：** 调用 gReLU 的 `grelu.resources.list_datasets`，用 `HfApi().get_collection` 读取模型库集合，保留 ID 以 `-data` 结尾的条目。
- **API：** `grelu_list_datasets(*, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** 无。
- **输出：** 数据集仓库 ID 的列表。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_list_datasets

try:
    grelu_list_datasets(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_list_datasets

datasets = grelu_list_datasets()
```

## 16) `dnakit.predictions.grelu_list_models()`

- **作用：** 列出 Hugging Face 上 gReLU 模型库的模型仓库 ID。
- **处理方法：** 调用 gReLU 的 `grelu.resources.list_models`，用 `HfApi().get_collection` 读取模型库集合，保留 ID 以 `-model` 结尾的条目。
- **API：** `grelu_list_models(*, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）以及网络连接。
- **输入：** 无。
- **输出：** 模型仓库 ID 的列表。
- **示例代码：**

```python
import os

os.environ["HF_HUB_OFFLINE"] = "1"

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_list_models

try:
    grelu_list_models(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "offline mode is enabled" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_list_models

models = grelu_list_models()
```

## 17) `dnakit.predictions.grelu_load_model()`

- **作用：** 从 Hugging Face 模型库加载 gReLU 模型。
- **处理方法：** 调用 gReLU 的 `grelu.resources.load_model`，它用 `download_model` 下载检查点，再用 `LightningModel.load_from_checkpoint(path, map_location=device)` 加载。`project` 和 `model_name` 参数仅用于给出清晰的错误：在 gReLU 1.1.0 中，传入其中任何一个都会抛出弃用错误，并指向 `grelu_wandb_load_model`。DNAKit 的文档字符串说本函数可以从 Weights & Biases 加载，但实际不能。
- **API：** `grelu_load_model(repo_id: 'str | None' = None, filename: 'str' = 'model.ckpt', *, device: 'str | int' = 'cpu', project: 'str | None' = None, model_name: 'str | None' = None, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 huggingface_hub（额外依赖组 `dnakit[bio-grelu]`）、首次下载时的网络连接，以及模型检查点（从 Hub 下载）。
- **输入：** `repo_id` 是模型仓库 ID（实际必需）；`filename` 是检查点文件（默认 `model.ckpt`）；`device` 是 `cpu` 或 GPU 序号。
- **输出：** gReLU 的 `LightningModel`。缺少 `repo_id` 或使用旧参数时抛出 `BackendExecutionError`。
- **示例代码：**

```python
from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_load_model

try:
    grelu_load_model(project="zoo", model_name="model", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "API has changed" in str(error))
try:
    grelu_load_model(verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "repo_id is required" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_load_model

model = grelu_load_model("Genentech/human-atac-catlas-model", device="cpu")
```

## 18) `dnakit.predictions.grelu_make_insertion_bigwig()`

- **作用：** 把 Tn5 插入位点的片段文件转换为分链的插入位点 bigWig。
- **处理方法：** 调用 gReLU 的 `grelu.data.preprocess.make_insertion_bigwig`。它运行一条 shell 管道（`cat` 或 `zcat`、`awk`、`sort`、用于 `chroms` 的可选 `grep`、`bedtools genomecov -bg -5`、`bedtools sort`），在 `tmp_dir` 中生成 bedGraph，再用基因组的染色体大小文件运行 `bedGraphToBigWig`，在 `out_dir` 中生成 bigWig，并删除 bedGraph。DNAKit 会先检查 `frag_file` 是否存在。
- **API：** `grelu_make_insertion_bigwig(frag_file: 'str', genome: 'str', *, out_prefix: 'str | None' = None, plus_shift: 'int' = 0, minus_shift: 'int' = 0, chroms: 'Sequence[str] | str | None' = None, tmp_dir: 'str' = './', out_dir: 'str' = './', verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`），以及 PATH 中的外部程序 `bedtools`、`awk`、`sort`、`cat`/`zcat` 和 UCSC 工具 `bedGraphToBigWig`。测试环境没有安装 `bedGraphToBigWig`，所以示例展示输入检查；基因组名称在首次使用时下载。
- **输入：** `frag_file` 是已存在的片段文件（类 BED，可为 `.gz`）；`genome` 是基因组名称或 FASTA 路径；`out_prefix` 是输出文件名（默认是片段文件名加 `.bw`；尽管名为前缀，实际用作完整文件名）；`plus_shift` 和 `minus_shift` 平移插入位点；`chroms` 限定染色体；`tmp_dir` 和 `out_dir` 是已存在的文件夹。
- **输出：** bigWig 文件路径字符串。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_make_insertion_bigwig

try:
    grelu_make_insertion_bigwig("missing_fragments.tsv", "tiny.fa", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] frag_file must be an existing file
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_make_insertion_bigwig

bigwig = grelu_make_insertion_bigwig(
    "fragments.tsv.gz", "hg38", out_prefix="insertions.bw", plus_shift=4, minus_shift=-5
)
```

## 19) `dnakit.predictions.grelu_plot_attention_matrix()`

- **作用：** 绘制 Transformer 注意力矩阵在序列窗口上的热图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_attention_matrix`，它用 `start_pos` 与 `end_pos` 之间的分箱坐标构建 seaborn 热图，并可为每个高亮区间画一个方框。图保存到 `output_path`。
- **API：** `grelu_plot_attention_matrix(attn: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int | None' = None, highlight_intervals: 'Any' = None, figsize: 'tuple[float, float]' = (5, 4), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`）和 matplotlib。
- **输入：** `attn` 是有限值的二维方阵；`output_path` 是已存在文件夹中的新 `.png`、`.pdf` 或 `.svg` 文件；`start_pos` 和 `end_pos` 给出基因组窗口（箱大小为 `(end_pos - start_pos) / n`）；`highlight_intervals` 是含 `start` 与 `end` 的 DataFrame。
- **输出：** 含 `status`、`backend`、`version`、`function`、`image_path` 和 `data` 的字典（`{"attn": {"shape": ..., "values": ...}}`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_attention_matrix

attn = np.eye(6) * 0.5 + 0.1
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_attention_matrix(attn, str(Path(folder) / "attention.png"), verbose=False)
    print(result["status"], result["function"], Path(result["image_path"]).is_file())
print(result["data"]["attn"]["shape"])
```

- **示例输出：**

```text
completed grelu.visualize.plot_attention_matrix True
[6, 6]
```

## 20) `dnakit.predictions.grelu_plot_attributions()`

- **作用：** 把单条序列的逐碱基归因分数绘制成序列标识图，并可加高亮。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_attributions`，它把 (4, L) 数组切成 `start_pos:end_pos`，并用 logomaker 画成标识图，每个字母的高度就是它的分数。用 `highlight_centers` 与 `highlight_width`，或 `highlight_positions` 添加高亮方框。默认的 `end_pos=-1` 会丢掉最后一个位置。
- **API：** `grelu_plot_attributions(attrs: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int' = -1, figsize: 'tuple[float, float]' = (20, 2), ticks: 'int' = 10, highlight_centers: 'Sequence[int] | None' = None, highlight_width: 'int | Sequence[int]' = 5, highlight_positions: 'Sequence[int] | None' = None, ylim: 'tuple[float, float] | None' = None, facecolor: 'str | None' = 'yellow', edgecolor: 'str | None' = None, alpha: 'float | None' = 0.15, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 logomaker 和 matplotlib。
- **输入：** `attrs` 是有限值的二维数组，形状 (4, L)，各行依次为 A、C、G、T。`output_path` 是新的图像文件。其余参数设置窗口、刻度、高亮位置、颜色和坐标轴范围。
- **输出：** 绘图字典（见上文）；`data["attrs"]` 保存形状和数值。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_attributions

attrs = np.zeros((4, 12))
attrs[2, 3:7] = [0.2, 0.9, 0.5, 0.7]
attrs[0, 7:10] = [0.4, 0.6, 0.3]
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_attributions(
        attrs, str(Path(folder) / "attributions.png"), end_pos=12,
        highlight_centers=[5], figsize=(6, 1.5), verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["attrs"]["shape"])
```

- **示例输出：**

```text
completed True
[4, 12]
```

## 21) `dnakit.predictions.grelu_plot_binary_preds()`

- **作用：** 按类别和任务，用箱线图绘制二分类标签的预测概率。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_binary_preds`，它把预测和标签展平成每个观测与位置一行，并为每个任务用 plotnine 画出按标签分组的预测箱线图。已知问题：gReLU 返回的 plotnine 对象没有被绘制，所以 DNAKit 保存的是空白图像（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_binary_preds(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `preds` 和 `labels` 是形状相同的有限值数组（观测数、任务数、长度）。`tasks` 选择任务索引或名称；`bins` 选择长度方向的位置；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 保存 `preds` 和 `labels` 的描述。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_binary_preds

rng = np.random.default_rng(0)
preds = rng.random((8, 2, 3))
labels = (rng.random((8, 2, 3)) > 0.5).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_binary_preds(preds, labels, str(Path(folder) / "binary.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"], result["data"]["labels"]["shape"])
```

- **示例输出：**

```text
completed True
[8, 2, 3] [8, 2, 3]
```

## 22) `dnakit.predictions.grelu_plot_calibration_curve()`

- **作用：** 绘制预测概率相对于观测二分类标签的校准曲线。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_calibration_curve`，它对概率分箱，并绘制每个箱中阳性的观测比例（`aggregate=True` 时在任务间聚合）。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_calibration_curve(probs: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, aggregate: 'bool' = True, figsize: 'tuple[float, float]' = (4, 3), show_legend: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `probs` 和 `labels` 是形状相同的有限值三维数组（观测数、任务数、长度）。DNAKit 也接受二维的 `probs`，但此时 gReLU 1.1.0 会以 `BackendExecutionError` 失败。`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 保存 `probs` 和 `labels` 的描述。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_calibration_curve

rng = np.random.default_rng(1)
probs = rng.random((20, 2, 3))
labels = (rng.random((20, 2, 3)) < probs).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_calibration_curve(
        probs, labels, str(Path(folder) / "calibration.png"), verbose=False
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["probs"]["shape"])
```

- **示例输出：**

```text
completed True
[20, 2, 3]
```

## 23) `dnakit.predictions.grelu_plot_distribution()`

- **作用：** 把一组数值的分布绘制成直方图或密度图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_distribution`，用 plotnine 构建直方图或密度图。DNAKit 函数的文档字符串还提到箱线图，但 gReLU 只支持 `histogram` 和 `density`；其他 `method` 会得到空图。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_distribution(values: 'Sequence[float]', output_path: 'str', *, title: 'str' = 'metric', method: 'str' = 'histogram', figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `values` 是一维数值序列；`title` 是坐标轴标签；`method` 是 `histogram` 或 `density`；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data["values"]` 是数值列表。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import grelu_plot_distribution

values = [0.1, 0.4, 0.35, 0.8, 0.6, 0.5, 0.45, 0.9]
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_distribution(
        values, str(Path(folder) / "distribution.png"), title="auroc", method="density",
        verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(len(result["data"]["values"]))
```

- **示例输出：**

```text
completed True
8
```

## 24) `dnakit.predictions.grelu_plot_evolution()`

- **作用：** 绘制定向进化各轮次中分数的变化。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_evolution`，它取 `iter`、`prediction_score`、`seq_score`、`total_score` 以及第九列之后的所有列，为每个分数按迭代轮次画箱线图（`seq_score` 全为零时，丢弃预测分数和序列分数）。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_evolution(df: 'Any', output_path: 'str', *, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `df` 是 `grelu.design.evolve` 产生格式的 DataFrame，含上述四个分数列，且总列数至少为九。`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data["df"]` 含 `columns`、`rows` 和 `values`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pandas as pd

from dnakit.predictions import grelu_plot_evolution

df = pd.DataFrame(
    {
        "iter": [0, 0, 1, 1],
        "prediction_score": [0.1, 0.2, 0.4, 0.5],
        "seq_score": [0.0, 0.0, 0.0, 0.0],
        "total_score": [0.1, 0.2, 0.4, 0.5],
        "seq": ["ACGT"] * 4,
        "chrom": ["x"] * 4,
        "start": [0] * 4,
        "end": [4] * 4,
        "strand": ["+"] * 4,
    }
)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_evolution(df, str(Path(folder) / "evolution.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["df"]["rows"], result["data"]["df"]["columns"][:2])
```

- **示例输出：**

```text
completed True
4 ['iter', 'prediction_score']
```

## 25) `dnakit.predictions.grelu_plot_gc_match()`

- **作用：** 绘制阳性区间与 GC 匹配的阴性区间的 GC 含量直方图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_gc_match`，它用 `gc_distribution` 在基因组上计算每组区间的 GC 分布，并用 plotnine 画出并列的条形图。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_gc_match(positives: 'Any', negatives: 'Any', output_path: 'str', *, binwidth: 'float' = 0.1, genome: 'str' = 'hg38', figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu、plotnine 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；本地 FASTA 路径可离线使用。
- **输入：** `positives` 和 `negatives` 是含 `chrom`、`start`、`end` 的 DataFrame；`binwidth` 是 GC 分辨率；`genome` 是基因组名称或 FASTA 路径（默认 `hg38`）；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 描述两个区间表（`columns`、`rows`、`values`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pandas as pd

from dnakit.predictions import grelu_plot_gc_match

positives = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 10], "end": [10, 20]})
negatives = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [20, 0], "end": [30, 10]})
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_gc_match(
        positives, negatives, str(Path(folder) / "gc.png"), binwidth=0.25,
        genome="tiny.fa", verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["positives"]["rows"], result["data"]["negatives"]["rows"])
```

- **示例输出：**

```text
completed True
2 2
```

## 26) `dnakit.predictions.grelu_plot_ism()`

- **作用：** 把体外饱和突变（ISM）表绘制成热图或序列标识图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_ISM`。`method="heatmap"` 时，对表的 `start_pos:end_pos` 区域画 seaborn 热图（行为 A、C、G、T，每个位置一列）。`method="logo"` 时，画每个位置负的平均突变效应的标识图，字母放在列名所给的参考碱基上。DNAKit 的文档字符串写的是“折线图”，但 gReLU 只支持 `heatmap` 和 `logo`；其他取值会以 `BackendExecutionError` 失败。
- **API：** `grelu_plot_ism(ism_preds: 'Any', output_path: 'str', *, start_pos: 'int | None' = None, end_pos: 'int | None' = None, figsize: 'tuple[float, float]' = (8, 1.5), method: 'str' = 'heatmap', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 seaborn、logomaker 和 matplotlib。
- **输入：** `ism_preds` 是以碱基 A、C、G、T 为行索引、每个序列位置一列的 DataFrame，列名为参考碱基（如 `grelu.model.interpret.ISM_predict` 所产生）。`start_pos` 和 `end_pos` 选择列；`method` 是 `heatmap` 或 `logo`；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data["ism_preds"]` 含 `columns`、`rows` 和 `values`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from dnakit.predictions import grelu_plot_ism

ism = pd.DataFrame(
    np.arange(24).reshape(4, 6) / 10, index=list("ACGT"), columns=list("ACGTAC")
)
with tempfile.TemporaryDirectory() as folder:
    heat = grelu_plot_ism(ism, str(Path(folder) / "ism.png"), verbose=False)
    logo = grelu_plot_ism(ism, str(Path(folder) / "ism_logo.svg"), method="logo", verbose=False)
    print(heat["status"], logo["status"], Path(logo["image_path"]).is_file())
print(heat["data"]["ism_preds"]["rows"], heat["data"]["ism_preds"]["columns"])
```

- **示例输出：**

```text
completed completed True
4 ['A', 'C', 'G', 'T', 'A', 'C']
```

## 27) `dnakit.predictions.grelu_plot_motif()`

- **作用：** 为四列位置矩阵表示的基序绘制序列标识图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_motif`，它把转置后的矩阵交给 tangermeme 的 `plot_pwm`。图保存到 `output_path`。
- **API：** `grelu_plot_motif(motif: 'Any', output_path: 'str', *, name: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 tangermeme 和 matplotlib。
- **输入：** `motif` 是有限值二维数组，形状 (长度, 4)，各列依次为 A、C、G、T（概率或分数）。`name` 是可选标题。`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data["motif"]` 保存形状和数值。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_motif

motif = np.array([[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1]])
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_motif(motif, str(Path(folder) / "motif.svg"), name="demo", verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["motif"]["shape"])
```

- **示例输出：**

```text
completed True
[3, 4]
```

## 28) `dnakit.predictions.grelu_plot_position_effect()`

- **作用：** 绘制模型输出如何随插入基序的位置而变化。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_position_effect`，它对每个位置计算均值、标准误和 95% 置信区间，并用 plotnine 绘制。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_position_effect(preds: 'Any', positions: 'Sequence[int]', output_path: 'str', *, title: 'str | None' = None, xlab: 'str | None' = None, figsize: 'tuple[float, float]' = (6, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine、scipy 和 matplotlib。
- **输入：** `preds` 是有限值二维数组（序列数、位置数），例如 `marginalize_pattern_spacing` 或 `shuffle_tiles` 的输出；`positions` 是每列一个值的列表；`title` 和 `xlab` 是可选标签；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 保存 `preds` 和 `positions`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_position_effect

preds = np.random.default_rng(0).random((5, 6))
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_position_effect(
        preds, [0, 10, 20, 30, 40, 50], str(Path(folder) / "position.png"),
        title="motif effect", xlab="position", verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"], result["data"]["positions"])
```

- **示例输出：**

```text
completed True
[5, 6] [0, 10, 20, 30, 40, 50]
```

## 29) `dnakit.predictions.grelu_plot_pred_distribution()`

- **作用：** 绘制每个任务中阳性和阴性标签对应的模型预测分布。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_pred_distribution`，它按任务收集预测和标签，并用 plotnine 画密度图。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_pred_distribution(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `preds` 和 `labels` 是形状相同的有限值三维数组（观测数、任务数、长度）；`tasks` 和 `bins` 选择任务和位置；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 保存 `preds` 和 `labels` 的描述。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_pred_distribution

rng = np.random.default_rng(2)
preds = rng.random((10, 2, 3))
labels = (rng.random((10, 2, 3)) > 0.5).astype(float)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_pred_distribution(
        preds, labels, str(Path(folder) / "pred_distribution.png"), tasks=[0], verbose=False
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["labels"]["shape"])
```

- **示例输出：**

```text
completed True
[10, 2, 3]
```

## 30) `dnakit.predictions.grelu_plot_pred_scatter()`

- **作用：** 为每个任务绘制预测值相对于实测标签的散点图。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_pred_scatter`，它按任务收集预测和标签，并用 plotnine 画散点图（`density=True` 时按密度着色）。已知问题：plotnine 对象在保存前没有被绘制，所以保存的图像是空白的（已用 gReLU 1.1.0 核实）。
- **API：** `grelu_plot_pred_scatter(preds: 'Any', labels: 'Any', output_path: 'str', *, tasks: 'Sequence[int] | Sequence[str] | None' = None, bins: 'Sequence[int] | None' = None, density: 'bool' = False, figsize: 'tuple[float, float]' = (4, 3), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`），以及 plotnine 和 matplotlib。
- **输入：** `preds` 和 `labels` 是形状相同的有限值三维数组（观测数、任务数、长度）；`tasks`、`bins` 和 `density` 控制绘制内容；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data` 保存 `preds` 和 `labels` 的描述。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_pred_scatter

rng = np.random.default_rng(3)
labels = rng.random((12, 2, 1))
preds = labels + rng.normal(scale=0.1, size=labels.shape)
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_pred_scatter(preds, labels, str(Path(folder) / "scatter.png"), verbose=False)
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["preds"]["shape"])
```

- **示例输出：**

```text
completed True
[12, 2, 1]
```

## 31) `dnakit.predictions.grelu_plot_tracks()`

- **作用：** 在序列窗口上绘制一条或多条预测或信号轨道。
- **处理方法：** 调用 gReLU 的 `grelu.visualize.plot_tracks`，用 matplotlib 在 `start_pos` 与 `end_pos` 之间为每条轨道画一个堆叠面板，可带面板标题和高亮区间。
- **API：** `grelu_plot_tracks(tracks: 'Any', output_path: 'str', *, start_pos: 'int' = 0, end_pos: 'int | None' = None, titles: 'Sequence[str] | None' = None, figsize: 'tuple[float, float]' = (20, 1.5), highlight_intervals: 'Any' = None, facecolor: 'str | None' = 'yellow', edgecolor: 'str | None' = None, alpha: 'float | None' = 0.15, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 grelu（额外依赖组 `dnakit[bio-grelu]`）和 matplotlib。
- **输入：** `tracks` 是有限值二维数组（轨道数、位置数）；`titles` 每条轨道一个标签；`highlight_intervals` 是含 `start` 和 `end` 的 DataFrame；`output_path` 是新的图像文件。
- **输出：** 绘图字典；`data["tracks"]` 保存形状和数值。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.predictions import grelu_plot_tracks

x = np.linspace(0, 6, 40)
tracks = np.vstack([np.sin(x) + 1, np.cos(x) + 1])
with tempfile.TemporaryDirectory() as folder:
    result = grelu_plot_tracks(
        tracks, str(Path(folder) / "tracks.png"), titles=["sin", "cos"], figsize=(6, 3),
        verbose=False,
    )
    print(result["status"], Path(result["image_path"]).is_file())
print(result["data"]["tracks"]["shape"])
```

- **示例输出：**

```text
completed True
[2, 40]
```

## 32) `dnakit.predictions.grelu_read_gtf()`

- **作用：** 把基因组注释中的基因特征读成表格。
- **处理方法：** 调用 gReLU 的 `grelu.io.genome.read_gtf`，它打开 `genomepy.Annotation(genome)` 并读取其 `named_gtf` 表。如果该基因组名称的注释尚未安装，gReLU 会用 genomepy 安装（需要下载，可能还需要 UCSC 工具）。genomepy 也接受 GTF 文件路径并直接读取，示例用的就是这种方式。表被重置为默认索引，`seqname` 改名为 `chrom`，再把 `chrom`、`start`、`end` 移到最前。`features` 过滤 `feature` 列。
- **API：** `grelu_read_gtf(genome: 'str', features: 'Sequence[str] | str | None' = None, *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时需要联网；本地 GTF 路径可离线使用。
- **输入：** `genome` 是非空的基因组名称或 GTF 文件路径；`features` 是特征名称或列表，例如 `exon`、`CDS` 或 `transcript`。
- **输出：** pandas DataFrame，列为 `chrom`、`start`、`end`、名称列、`source`、`feature`、`score`、`strand`、`frame` 和 `attribute`。genomepy 关于缺少同目录文件的提示写入日志，而不是标准输出。
- **示例代码：**

```python
from dnakit.predictions import grelu_read_gtf

table = grelu_read_gtf("genes.gtf", verbose=False)
print(table.shape, table.columns.tolist()[:3])
exons = grelu_read_gtf("genes.gtf", "exon", verbose=False)
print(exons[["chrom", "start", "end"]].to_dict("records"))
```

- **示例输出：**

```text
(4, 10) ['chrom', 'start', 'end']
[{'chrom': 'chr1', 'start': 101, 'end': 200}, {'chrom': 'chr1', 'start': 301, 'end': 400}]
```

## 33) `dnakit.predictions.grelu_read_sizes()`

- **作用：** 把基因组的染色体大小读成表格。
- **处理方法：** 调用 gReLU 的 `grelu.io.genome.read_sizes`，它用 `get_genome` 加载基因组（给出名称时用 genomepy 安装），并读取其 `sizes_file`。
- **API：** `grelu_read_sizes(genome: 'str' = 'hg38', *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时需要联网；带 `.sizes` 文件的本地 FASTA 路径可离线使用。
- **输入：** `genome` 是非空的基因组名称（默认 `hg38`）或 FASTA 路径。
- **输出：** 含 `chrom` 和 `size` 两列的 DataFrame。
- **示例代码：**

```python
from dnakit.predictions import grelu_read_sizes

sizes = grelu_read_sizes("tiny.fa", verbose=False)
print(sizes.to_dict("records"))
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'size': 40}, {'chrom': 'chr2', 'size': 30}]
```

## 34) `dnakit.predictions.grelu_run_tomtom()`

- **作用：** 用 tomtom 算法把查询基序与参考 MEME 数据库比较。
- **处理方法：** 调用 gReLU 的 `grelu.interpret.motifs.run_tomtom`，它使用 memelite 的纯 Python `tomtom`（不使用 MEME Suite 程序）。它读取 MEME 文件，对每个查询与目标的组合计算 p 值、偏移、重叠和方向，并加上 E 值（p 值乘以目标数）、q 值（Benjamini-Hochberg）和共有序列字符串。注意：DNAKit 检查每个基序是否有四列，其文档字符串也说行是位置，但 gReLU 1.1.0 把每个矩阵读作形状 (4, L)（行为 A、C、G、T）。矩阵被原样传递，所以长度不是 4 的 (长度, 4) 矩阵会被读错，而四行、其他长度的矩阵会被拒绝。只有四个位置时，请按 gReLU 的方向（行为 A、C、G、T）传入矩阵，如示例所示，并检查 `Query_consensus` 列。目标少于 25 个时，tomtom 还会发出警告，因为此时 p 值取决于目标集合。
- **API：** `grelu_run_tomtom(motifs: 'Mapping[str, Sequence[Sequence[float]]]', meme_file: 'str', *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu、memelite 和 statsmodels（额外依赖组 `dnakit[bio-grelu]`）。可离线运行。
- **输入：** `motifs` 是从名称到有限值四列概率矩阵的非空映射（见上面的注意事项）。`meme_file` 是已存在的 MEME 文件路径，或内置名称 `hocomoco_v12`、`hocomoco_v13` 和 `consensus` 之一。
- **输出：** 每个查询与目标组合一行的 DataFrame，列为 `Query_ID`、`Target_ID`、`Optimal_offset`、`p-value`、`E-value`、`q-value`、`Overlap`、`Query_consensus`、`Target_consensus` 和 `Orientation`。
- **示例代码：**

```python
import warnings

from dnakit.predictions import grelu_run_tomtom

warnings.simplefilter("ignore")
# rows A, C, G, T; four positions: G A T A
query = [
    [0.05, 0.85, 0.05, 0.85],
    [0.05, 0.05, 0.05, 0.05],
    [0.85, 0.05, 0.05, 0.05],
    [0.05, 0.05, 0.85, 0.05],
]
table = grelu_run_tomtom({"query": query}, "db.meme", verbose=False)
best = table.loc[table["p-value"].idxmin()]
print(table.shape, best["Target_ID"], best["Query_consensus"], best["Orientation"])
```

- **示例输出：**

```text
(2, 10) GATA_like GATA +
```

## 35) `dnakit.predictions.grelu_variant_to_seqs()`

- **作用：** 返回单核苷酸变异周围的参考序列和替代序列。
- **处理方法：** 调用 gReLU 的 `grelu.variant.variant_to_seqs`。它先运行 `check_reference`（不一致只发出警告），再生成以变异为中心、长度为 `seq_len` 的区间，提取参考序列，并在中心位置分别插入参考等位基因和替代等位基因。
- **API：** `grelu_variant_to_seqs(chrom: 'str', pos: 'int', ref: 'str', alt: 'str', genome: 'str', *, seq_len: 'int' = 1, verbose: bool = True) -> 'tuple[str, str]'`。
- **运行要求：** 需要可选包 grelu 和 genomepy（额外依赖组 `dnakit[bio-grelu]`）。基因组名称在首次使用时下载；本地 FASTA 路径可离线使用。
- **输入：** `chrom` 和 `genome` 是非空字符串；`pos` 是 1 起始的位置；`ref` 和 `alt` 是单个碱基；`seq_len` 是返回序列的长度（默认 1，此时只返回等位基因本身）。
- **输出：** 元组 `(参考序列, 替代序列)`。
- **示例代码：**

```python
from dnakit.predictions import grelu_variant_to_seqs

ref_seq, alt_seq = grelu_variant_to_seqs("chr1", 5, "T", "G", "tiny.fa", seq_len=5, verbose=False)
print(ref_seq, alt_seq)
```

- **示例输出：**

```text
GTTGC GTGGC
```

## 36) `dnakit.predictions.grelu_wandb_artifacts()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上某个项目的制品（artifact）。
- **处理方法：** 调用 `grelu.resources.wandb.artifacts`。它先登录 Weights & Biases 的 `host`（`wandb.login`，匿名访问），再列出项目 `grelu/<project>` 的制品类型（可只取等于 `type_is` 或包含 `type_contains` 的类型），并返回这些类型下各集合的名称。没有匹配的类型时会失败。
- **API：** `grelu_wandb_artifacts(project: 'str', *, host: 'str' = 'https://api.wandb.ai', type_is: 'str | None' = None, type_contains: 'str | None' = None, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接（会请求匿名访问；有 API 密钥时可事先登录）。
- **输入：** `project` 是 `grelu` 实体下的非空项目名；`host` 是 W&B API 地址；`type_is` 和 `type_contains` 过滤制品类型。
- **输出：** 制品集合名称的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_artifacts

try:
    grelu_wandb_artifacts("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError project must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_artifacts

names = grelu_wandb_artifacts("human-atac-catlas", type_contains="model")
```

## 37) `dnakit.predictions.grelu_wandb_datasets()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上某个项目的数据集。
- **处理方法：** 调用 `grelu.resources.wandb.datasets`，即限定为类型名包含 `dataset` 的 `artifacts`。
- **API：** `grelu_wandb_datasets(project: 'str', *, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `project` 是非空项目名；`host` 是 W&B API 地址。
- **输出：** 数据集制品名称的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_datasets

try:
    grelu_wandb_datasets("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError project must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_datasets

names = grelu_wandb_datasets("human-atac-catlas")
```

## 38) `dnakit.predictions.grelu_wandb_get_artifact()`

- **作用：** 按名称和别名获取旧版 gReLU 模型库中的 Weights & Biases 制品。
- **处理方法：** 调用 `grelu.resources.wandb.get_artifact`，它登录 Weights & Biases 并返回 `wandb.Api().artifact("grelu/<project>/<name>:<alias>")`。返回的是制品对象；调用 `artifact.download()` 之前不会下载文件。
- **API：** `grelu_wandb_get_artifact(name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `name` 和 `project` 是非空字符串；`alias` 是制品别名（默认 `latest`）；`host` 是 W&B API 地址。
- **输出：** `wandb` 制品对象。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_artifact

try:
    grelu_wandb_get_artifact("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError name must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_get_artifact

artifact = grelu_wandb_get_artifact("model", "human-atac-catlas")
```

## 39) `dnakit.predictions.grelu_wandb_get_dataset_by_model()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上与某个模型关联的数据集。
- **处理方法：** 调用 `grelu.resources.wandb.get_dataset_by_model`。它获取模型制品，找到记录该制品的运行，并返回该运行使用过的制品名称。
- **API：** `grelu_wandb_get_dataset_by_model(model_name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `model_name` 和 `project` 是非空字符串；`alias` 是模型别名；`host` 是 W&B API 地址。
- **输出：** 数据集制品名称的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_dataset_by_model

try:
    grelu_wandb_get_dataset_by_model("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError model_name must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_get_dataset_by_model

datasets = grelu_wandb_get_dataset_by_model("model", "human-atac-catlas")
```

## 40) `dnakit.predictions.grelu_wandb_get_model_by_dataset()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上用某个数据集训练的模型。
- **处理方法：** 调用 `grelu.resources.wandb.get_model_by_dataset`。它获取数据集制品，取第一个使用过它的运行，并返回该运行记录的制品名称。没有任何运行使用该数据集时会失败。
- **API：** `grelu_wandb_get_model_by_dataset(dataset_name: 'str', project: 'str', *, host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `dataset_name` 和 `project` 是非空字符串；`alias` 是数据集别名；`host` 是 W&B API 地址。
- **输出：** 模型制品名称的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_get_model_by_dataset

try:
    grelu_wandb_get_model_by_dataset("", "human-atac-catlas", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError dataset_name must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_get_model_by_dataset

models = grelu_wandb_get_model_by_dataset("dataset", "human-atac-catlas")
```

## 41) `dnakit.predictions.grelu_wandb_load_model()`

- **作用：** 从旧版 Weights & Biases 模型库的项目中加载 gReLU 模型。
- **处理方法：** 调用 `grelu.resources.wandb.load_model`。它获取模型制品，下载到临时文件夹，并用 `LightningModel.load_from_checkpoint` 加载 `checkpoint_file`。`weights_only=False` 是默认值，因为 gReLU 检查点包含的不只是张量；加载来源不可信的检查点可能执行代码。
- **API：** `grelu_wandb_load_model(project: 'str', model_name: 'str', *, device: 'str | int' = 'cpu', host: 'str' = 'https://api.wandb.ai', alias: 'str' = 'latest', checkpoint_file: 'str' = 'model.ckpt', weights_only: 'bool' = False, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）、访问 Weights & Biases 的网络连接，以及模型检查点（从制品下载）。
- **输入：** `project` 和 `model_name` 是非空字符串；`device` 是 `cpu` 或 GPU 序号；`alias` 是模型别名；`checkpoint_file` 是制品中的文件名；`weights_only` 传给检查点加载器。
- **输出：** gReLU 的 `LightningModel`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_load_model

try:
    grelu_wandb_load_model("human-atac-catlas", "", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError model_name must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_load_model

model = grelu_wandb_load_model("human-atac-catlas", "model", device="cpu")
```

## 42) `dnakit.predictions.grelu_wandb_models()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上某个项目的模型。
- **处理方法：** 调用 `grelu.resources.wandb.models`，即限定为类型名包含 `model` 的 `artifacts`。
- **API：** `grelu_wandb_models(project: 'str', *, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `project` 是非空项目名；`host` 是 W&B API 地址。
- **输出：** 模型制品名称的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_models

try:
    grelu_wandb_models("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError project must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_models

names = grelu_wandb_models("human-atac-catlas")
```

## 43) `dnakit.predictions.grelu_wandb_projects()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上的项目。
- **处理方法：** 调用 `grelu.resources.wandb.projects`。它登录 `host` 上的 Weights & Biases（`wandb.login`，先尝试匿名访问），并列出 `grelu` 实体下的项目。登录失败时 gReLU 抛出 `RuntimeError`，DNAKit 将其重新抛出为 `BackendExecutionError`。
- **API：** `grelu_wandb_projects(*, host: 'str' = 'https://api.wandb.ai', verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `host` 是 W&B API 地址（默认 `https://api.wandb.ai`）。
- **输出：** 项目名称的列表。
- **示例代码：**

```python
import os
import tempfile

folder = tempfile.mkdtemp()
os.environ["HOME"] = folder
for name in ("WANDB_CONFIG_DIR", "WANDB_CACHE_DIR", "WANDB_DIR"):
    os.environ[name] = folder
os.environ.pop("WANDB_API_KEY", None)

from dnakit.exceptions import BackendExecutionError
from dnakit.predictions import grelu_wandb_projects

try:
    grelu_wandb_projects(host="http://127.0.0.1:9", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, "is not configured" in str(error))
```

- **示例输出：**

```text
BackendExecutionError True
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_projects

projects = grelu_wandb_projects()
```

## 44) `dnakit.predictions.grelu_wandb_runs()`

- **作用：** 列出旧版 gReLU 模型库在 Weights & Biases 上某个项目的每个运行的一个属性。
- **处理方法：** 调用 `grelu.resources.wandb.runs`，它登录 Weights & Biases，对 `grelu/<project>` 的每个运行返回 `getattr(run, field)`，并可用传给 `wandb.Api().runs` 的 `filters` 映射过滤。
- **API：** `grelu_wandb_runs(project: 'str', *, host: 'str' = 'https://api.wandb.ai', field: 'str' = 'id', filters: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'list[Any]'`。
- **运行要求：** 需要可选包 grelu 和 wandb（额外依赖组 `dnakit[bio-grelu]`）以及访问 Weights & Biases 的网络连接。
- **输入：** `project` 是非空项目名；`field` 是运行属性，例如 `id`、`name` 或 `state`；`filters` 是 MongoDB 风格的过滤映射；`host` 是 W&B API 地址。
- **输出：** 每个运行一个值的列表。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import grelu_wandb_runs

try:
    grelu_wandb_runs("", field="name", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[1])
```

- **示例输出：**

```text
ConfigurationError project must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.predictions import grelu_wandb_runs

run_ids = grelu_wandb_runs("human-atac-catlas", field="id")
```
