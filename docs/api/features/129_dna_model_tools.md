# DNA 模型工具包：AlphaGenome、Enformer、Selene、tangermeme 与 genomic-benchmarks

这一组函数是五个可选包的适配器，它们都面向基于序列的深度学习模型：[AlphaGenome](https://github.com/google-deepmind/alphagenome)（Google DeepMind 的 SDK，用于托管的 AlphaGenome 模型及其绘图工具）、[enformer-pytorch](https://github.com/lucidrains/enformer-pytorch)（Enformer 模型）、[Selene](https://github.com/FunctionLab/selene)（selene-sdk，训练和应用序列模型的框架）、[tangermeme](https://github.com/jmschrei/tangermeme)（序列模型的归因与基序绘图）和 [genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks)（基因组序列的基准数据集）。它们可以连接 AlphaGenome 服务，选择和补齐预测轨道，绘制轨道、序列标志图、接触图、转录本和 ROC 曲线，提取 Enformer 嵌入，运行 Selene 流程和模型预测，计算归因，以及加载基准数据集。AlphaGenome、Enformer、Selene 和 genomic-benchmarks 的函数在 `dnakit.predictions` 中，tangermeme 的函数在 `dnakit.regulation` 中。示例用 `from dnakit.predictions import ...` 和 `from dnakit.regulation import ...` 导入。

依赖是可选的，只有在调用相应函数时才会导入。缺少某个包时，函数会抛出 `BackendUnavailableError` 并附带安装命令。包元数据给出这些 Python 代码的许可证：alphagenome 和 genomic-benchmarks 为 Apache-2.0，enformer-pytorch 和 tangermeme 为 MIT，selene-sdk 为 BSD；都不是 GPL。这些许可证只覆盖相应的包。托管的 AlphaGenome 服务有自己的使用条款并需要 API 密钥，预训练的 Enformer 权重和基准数据集则来自第三方服务器，各有各的条款。

`predict_selene_sequences()` 会反序列化完整的模型文件，Selene 流程函数会导入配置中指定的模型 Python 文件和类。请只使用可信的文件。

## 示例数据

Selene 流程示例（第 15、16、17 节）共用三个声明的文件：一个很小的模型定义、一个用 `analyze` 操作给 FASTA 文件中的序列打分的 YAML 配置，以及两条序列。其余示例都在代码中构造输入。

**示例文件：** `model.py`

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

**示例文件：** `seqs.fa`

```text
>s1
ACGTACGT
>s2
TTTTGGGG
```

**示例文件：** `config.yml`

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

- **作用：** 用 API 密钥连接 AlphaGenome Atlas 服务，并描述这个连接。
- **处理方法：** 用密钥、超时和地址调用 `alphagenome.atlas.atlas.create()`。alphagenome 会打开一条安全的 gRPC 通道并等待它就绪。函数返回客户端类名、地址和超时；客户端对象本身不返回（它保持通道打开），密钥也不会出现在结果中。
- **API：** `create_alphagenome_atlas_client(api_key: 'str', *, timeout: 'float | None' = None, address: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、网络连接和 AlphaGenome API 密钥。
- **输入：** `api_key` 为非空字符串。`timeout` 为正数秒，或 `None` 表示不限时等待。`address` 为可选的非空服务器地址，默认 `dns:///gdmscience.googleapis.com:443`。
- **输出：** `service`（`alphagenome.atlas`）、`source`（实际使用的地址）和 `result`（含 `client_type`、`address`、`timeout_seconds`）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import create_alphagenome_atlas_client

try:
    create_alphagenome_atlas_client("   ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **示例输出：**

```text
ConfigurationError api_key must be a non-empty string.
```

示例不会建立连接：密钥在使用网络之前就会被检查。安装依赖后的用法：

```python
from dnakit.predictions import create_alphagenome_atlas_client

info = create_alphagenome_atlas_client("YOUR_API_KEY", timeout=30.0)
print(info["result"]["client_type"])
```

## 2) `dnakit.predictions.create_alphagenome_dna_client()`

- **作用：** 用 API 密钥连接托管的 AlphaGenome DNA 模型服务，并描述这个客户端。
- **处理方法：** 把 `model_version` 映射为 alphagenome 的 `ModelVersion` 枚举，再调用 `alphagenome.models.dna_client.create()`。alphagenome 会打开一条安全的 gRPC 通道并等待它就绪。客户端对象不返回，密钥也不会出现在结果中。
- **API：** `create_alphagenome_dna_client(api_key: 'str', *, model_version: 'str | None' = None, timeout: 'float | None' = None, address: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、网络连接和 AlphaGenome API 密钥。
- **输入：** `api_key` 为非空字符串。`model_version` 为 `None`（使用 alphagenome 的默认值，报告为 `all_folds`）或 `all_folds`、`fold_0`、`fold_1`、`fold_2`、`fold_3` 之一。`timeout` 和 `address` 同 `create_alphagenome_atlas_client()`。
- **输出：** `service`（`alphagenome.models.dna_client`）、`source`（实际使用的地址）和 `result`（含 `client_type`、`address`、`model_version`、`model_version_source`（`requested` 或 `upstream default`）、`timeout_seconds`）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import create_alphagenome_dna_client

try:
    create_alphagenome_dna_client("key", model_version="fold_9", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **示例输出：**

```text
ConfigurationError model_version must be one of: all_folds, fold_0, fold_1, fold_2, fold_3.
```

示例不会建立连接：模型版本在使用网络之前就会被检查。安装依赖后的用法：

```python
from dnakit.predictions import create_alphagenome_dna_client

info = create_alphagenome_dna_client("YOUR_API_KEY", model_version="fold_0", timeout=30.0)
print(info["result"]["model_version"])
```

## 3) `dnakit.predictions.get_alphagenome_junctions_to_plot()`

- **作用：** 按 alphagenome 的规则，选出某条预测轨道在某条链上的剪接连接。
- **处理方法：** 用连接列表、数值矩阵和轨道名构建 alphagenome 的 `JunctionData`，再调用 `alphagenome.data.junction_data.get_junctions_to_plot()`。它保留位于指定链上、且该轨道的值不低于阈值的连接，并把 `k` 设为该值（保留两位小数）。本函数只选数据，绘图请用 `plot_alphagenome_sashimi()`。
- **API：** `get_alphagenome_junctions_to_plot(junctions: 'list[dict[str, Any]]', values: 'list[list[float]]', track_names: 'list[str]', name: 'str', strand: 'str', *, k_threshold: 'float | None' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）和 numpy（extra `dnakit[bio-analysis]`）。可离线运行。
- **输入：** `junctions` 为对象列表，含 `chromosome`、`start`、`end`（大于 `start`）和 `strand`（`+` 或 `-`）。`values` 每个连接一行、每条轨道一列。`track_names` 为不重复的名称，每列一个。`name` 为要选的轨道，`strand` 为 `+` 或 `-`。`k_threshold` 为数字，或 `None`（采用 alphagenome 的规则：最大值的 5 %）。
- **输出：** `track`、`strand`、`k_threshold`、`count` 和 `junctions`（含 `chromosome`、`start`、`end`、`strand`、`k` 的对象）。
- **示例代码：**

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

- **示例输出：**

```text
2
100 200 0.5
300 400 3.0
```

## 4) `dnakit.predictions.pad_alphagenome_track()`

- **作用：** 把一条“位置 × 通道”的轨道居中放进更长的窗口。
- **处理方法：** 调用 alphagenome 的 `pad_track()`：长度差为偶数时两侧补得一样多，为奇数时多出的一个位置放在末尾。函数把原始第一行的下标报告为 `offset`。
- **API：** `pad_alphagenome_track(track: 'list[list[float]]', new_len: 'int', *, value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）和 numpy（extra `dnakit[bio-analysis]`）。可离线运行。
- **输入：** `track` 为非空、矩形的有限数字行列表。`new_len` 为整数，不小于行数。`value` 为有限的填充数字。
- **输出：** `length`、`new_length`、`offset`、`value` 和 `padded`（补齐后的各行）。
- **示例代码：**

```python
from dnakit.predictions import pad_alphagenome_track

result = pad_alphagenome_track([[1.0, 2.0], [3.0, 4.0]], 5, value=-1.0, verbose=False)
print(result["length"], result["new_length"], result["offset"])
print(result["padded"])
```

- **示例输出：**

```text
2 5 1
[[-1.0, -1.0], [1.0, 2.0], [3.0, 4.0], [-1.0, -1.0], [-1.0, -1.0]]
```

## 5) `dnakit.predictions.plot_alphagenome_components()`

- **作用：** 把多个模型输出叠成多个组件，画在同一个基因组区间上并保存为 PNG。
- **处理方法：** 把每个组件对象转换为 alphagenome 组件（`Tracks`、`SeqLogo` 或 `EmptyComponent`），再对该区间调用 `alphagenome.visualization.plot_components.plot()`。图以 150 dpi 保存后关闭。
- **API：** `plot_alphagenome_components(components: 'list[dict[str, Any]]', output_path: 'str', *, chromosome: 'str', start: 'int', end: 'int', fig_width: 'float' = 20.0, fig_height_scale: 'float' = 1.0, title: 'str | None' = None, despine: 'bool' = True, despine_keep_bottom: 'bool' = False, hspace: 'float' = 0.3, xlabel: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `components` 为非空对象列表。`tracks` 对象含 `names`、`values`（每个碱基一行、每条轨道一列），可选 `strands`（`+`、`-` 或 `.`）、`track_height` 和 `filled`。`seq_logo` 对象含 `scores`（每个碱基一行、每个字母一列），可选 `alphabet`、`fig_height` 和 `ylabel`。`empty` 对象可带 `fig_height`。行数必须等于 `end - start`。这里不支持转录本和接触图组件。`output_path` 必须以 `.png` 结尾；缺少的文件夹会被创建，已有文件会被覆盖。
- **输出：** `path`（PNG 文件）和 `data`（含 `interval`、各组件的描述 `components` 与 `title`）。
- **示例代码：**

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

- **示例输出：**

```text
True
['tracks', 'seq_logo', 'empty'] {'chromosome': 'chr1', 'start': 100, 'end': 104}
```

## 6) `dnakit.predictions.plot_alphagenome_contact_map()`

- **作用：** 用对数色标画一张方形接触图并保存为 PNG。
- **处理方法：** 对矩阵调用 alphagenome 的 `plot_contact_map()`。色标范围为给定的 `vmin`、`vmax`，缺省时取矩阵中最小的正值和最大值；结果会报告实际使用的色标。
- **API：** `plot_alphagenome_contact_map(contact_map: 'list[list[float]]', output_path: 'str', *, vmin: 'float | None' = None, vmax: 'float | None' = None, square: 'bool' = True, cbar_shrink: 'float' = 0.4, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `contact_map` 为方阵，元素为有限数字，且至少含一个正值（或显式给出 `vmin`）。`vmin` 必须为正且小于 `vmax`。`square` 为布尔值，`cbar_shrink` 为正数。`output_path` 必须以 `.png` 结尾。
- **输出：** `path` 和 `data`（含矩阵、`vmin`、`vmax`、`square`、`cbar_shrink`）。
- **示例代码：**

```python
from dnakit.predictions import plot_alphagenome_contact_map

matrix = [[1.0, 0.5, 0.1], [0.5, 1.0, 0.4], [0.1, 0.4, 1.0]]
result = plot_alphagenome_contact_map(matrix, "contacts.png", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["vmin"], result["data"]["vmax"])
```

- **示例输出：**

```text
True
0.1 1.0
```

## 7) `dnakit.predictions.plot_alphagenome_sashimi()`

- **作用：** 把剪接连接画成 Sashimi 图并保存为 PNG。
- **处理方法：** 用读段数构建 alphagenome 的 `Junction` 对象，再用由 `seed` 播种的随机数生成器调用 `alphagenome.visualization.plot.sashimi_plot()`，因此弧线抖动和 PNG 可重现。结果还列出不低于 `filter_threshold` 的连接，即实际画出的连接。
- **API：** `plot_alphagenome_sashimi(junctions: 'list[dict[str, Any]]', output_path: 'str', *, interval: 'dict[str, Any] | None' = None, filter_threshold: 'float' = 0.01, annotate_counts: 'bool' = True, seed: 'int | None' = 0, color: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `junctions` 为对象列表，含 `chromosome`、`start`、`end`、`strand`（`+` 或 `-`）和 `k`（非负的读段数）；至少一个 `k` 要大于 0。`interval` 为可选对象（`chromosome`、`start`、`end`），用来限定计数标注的范围。`filter_threshold` 为非负数。`seed` 为整数或 `None`。`color` 为可选的颜色名。
- **输出：** `path` 和 `data`（含 `junctions`、`drawn`、`drawn_count`、`total_k`、`filter_threshold`、`annotate_counts`、`seed`、`ylim`）。
- **示例代码：**

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

- **示例输出：**

```text
True
2 7.001 0
```

## 8) `dnakit.predictions.plot_alphagenome_seqlogo()`

- **作用：** 根据每个位置的字母高度画序列标志图并保存为 PNG。
- **处理方法：** 在 matplotlib 坐标轴上调用 alphagenome 的 `seqlogo()`。正的高度叠在坐标轴上方，负的高度叠在下方。
- **API：** `plot_alphagenome_seqlogo(letter_heights: 'list[list[float]]', output_path: 'str', *, alphabet: 'str' = 'ACGT', one_based: 'bool' = True, start: 'int' = 0, letter_colors: 'str' = 'default', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `letter_heights` 每个位置一行，每个 `alphabet` 字母一列。`one_based` 和 `start` 决定 x 轴的编号。`letter_colors` 为 `default` 或 `factorbook`。`output_path` 必须以 `.png` 结尾。
- **输出：** `path` 和 `data`（含 `letter_heights`、`alphabet`、`one_based`、`start`、`letter_colors`、`positions`）。
- **示例代码：**

```python
from dnakit.predictions import plot_alphagenome_seqlogo

heights = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, -1]]
result = plot_alphagenome_seqlogo(heights, "seqlogo.png", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["positions"], result["data"]["alphabet"], result["data"]["letter_colors"])
```

- **示例输出：**

```text
True
3 ACGT default
```

## 9) `dnakit.predictions.plot_alphagenome_track()`

- **作用：** 画一条轨道并保存为 PNG，绘制方式由轨道的形状决定。
- **处理方法：** 调用 alphagenome 的 `plot_track()`。一维列表或单列的二维列表画成折线（`filled` 为真时填充）；四列的二维列表画成序列标志图；两列的二维列表画成两条折线（正链和负链）。所选的画法以 `rendering` 报告。
- **API：** `plot_alphagenome_track(track: 'list[Any]', output_path: 'str', *, x: 'list[float] | None' = None, legend: 'bool' = False, ylim: 'list[float] | None' = None, color: 'str | list[str] | None' = None, filled: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `track` 为一维列表，或有 1、2、4 列的二维列表。`x` 为可选列表，每行一个位置。`ylim` 为 `[low, high]`，`high` 要大于 `low`。`color` 为颜色或颜色列表。`legend` 和 `filled` 为布尔值。
- **输出：** `path` 和 `data`（含 `track`、`length`、`rendering`（`line`、`sequence_logo` 或 `two_strand_lines`）、`x`、`ylim`、`filled`）。
- **示例代码：**

```python
from dnakit.predictions import plot_alphagenome_track

result = plot_alphagenome_track([0.1, 0.5, 0.9, 0.4], "track.png", filled=True, verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["rendering"], result["data"]["length"], result["data"]["filled"])
```

- **示例输出：**

```text
True
line 4 True
```

## 10) `dnakit.predictions.plot_alphagenome_tracks()`

- **作用：** 把多条轨道画成上下堆叠的面板，放在同一个 PNG 中。
- **处理方法：** 对命名的轨道调用 alphagenome 的 `plot_tracks()`。每条轨道一个面板，y 轴可以各自独立（`ylim="auto"`）、共用同一个范围（`"same"`）或使用给定范围。
- **API：** `plot_alphagenome_tracks(tracks: 'dict[str, list[Any]]', output_path: 'str', *, x: 'list[float] | None' = None, title: 'str | None' = None, legend: 'bool' = False, fig_width: 'float' = 20.0, fig_track_height: 'float | dict[str, float]' = 1.5, ylim: 'str | list[float] | None' = 'auto', yticks_min_max_only: 'bool' = False, ylab: 'bool' = True, color: 'str | dict[str, str] | None' = None, horizontal_ylab: 'bool' = True, filled_tracks: 'list[str] | None' = None, despine: 'bool' = True, despine_keep_bottom: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）、numpy（extra `dnakit[bio-analysis]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `tracks` 把每个轨道名映射到一维或二维列表；所有轨道长度必须相同。`ylim` 为 `auto`、`same`、`[low, high]` 或 `None`。`fig_track_height` 为一个数字，或为每条轨道各给一个高度的字典。`color` 为一种颜色，或只包含已有轨道名的字典。`filled_tracks` 为轨道名列表。其余参数为布尔值、数字和标题。
- **输出：** `path` 和 `data`（含 `tracks`、`track_names`、`length`、`x`、`ylim`、`title`、`filled_tracks`）。
- **示例代码：**

```python
from dnakit.predictions import plot_alphagenome_tracks

tracks = {"liver": [0, 1, 2, 1], "brain": [2, 1, 0, 1]}
result = plot_alphagenome_tracks(
    tracks, "tracks.png", title="Two tracks", filled_tracks=["brain"], verbose=False
)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["track_names"], result["data"]["length"], result["data"]["filled_tracks"])
```

- **示例输出：**

```text
True
['liver', 'brain'] 4 ['brain']
```

## 11) `dnakit.predictions.plot_alphagenome_transcripts()`

- **作用：** 把转录本模型画在一个基因组区间上并保存为 PNG。
- **处理方法：** 用外显子和编码区间构建 alphagenome 的 `Transcript` 对象，再在 matplotlib 坐标轴上调用 `alphagenome.visualization.plot_transcripts.plot_transcripts()`。图的高度随转录本数量增加。
- **API：** `plot_alphagenome_transcripts(transcripts: 'list[dict[str, Any]]', output_path: 'str', *, chromosome: 'str', start: 'int', end: 'int', zero_origin: 'bool' = False, label_name: 'str | None' = None, plot_labels_once: 'bool' = False, fig_width: 'float' = 12.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 alphagenome（extra `dnakit[bio-alphagenome]`）和 matplotlib（extra `dnakit[bio-plots]`）。可离线运行。
- **输入：** `transcripts` 为非空对象列表，含 `transcript_id`、`strand`（`+` 或 `-`）、`exons`（`[start, end]` 对的列表），可选 `cds`（区间对；非编码转录本可省略或设为 null）和 `info`（由简单值构成的对象）。`chromosome`、`start`、`end` 给出区间。`label_name` 指定显示在每个转录本旁边的 `info` 键，所有转录本都必须有这个键。`zero_origin` 和 `plot_labels_once` 为布尔值。
- **输出：** `path` 和 `data`（含 `interval`、各转录本的描述 `transcripts`、`zero_origin`、`label_name`、`plot_labels_once`）。
- **示例代码：**

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

- **示例输出：**

```text
True
['T1', 'T2'] gene_name
```

## 12) `dnakit.predictions.download_pretrained_enformer()`

- **作用：** 从 Hugging Face hub 下载预训练的 Enformer 权重并保存到本地文件。
- **处理方法：** 调用 `enformer_pytorch.modeling_enformer.from_pretrained()`（仅当 `use_tf_gamma` 不是 `None` 时才传入它），用 `torch.save()` 保存 `model.state_dict()`，并在旁边把模型配置写成按键排序的 JSON。
- **API：** `download_pretrained_enformer(name: 'str', output_state_dict: 'str', *, use_tf_gamma: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 enformer-pytorch 和 torch（extra `dnakit[neural-enformer]`），首次调用需要联网（之后使用 Hugging Face 缓存）。函数在缺包时给出的安装提示写的是名为 `enformer` 的 extra，而 `pyproject.toml` 中实际的 extra 是 `neural-enformer`。
- **输入：** `name` 为非空的 Hugging Face 模型 ID，例如 `EleutherAI/enformer-official-rough`。`output_state_dict` 为已有文件夹中的新文件路径；`<output_state_dict>.config.json` 也不能已存在。`use_tf_gamma` 为 `None`、真或假。
- **输出：** `status`、`backend`、`model`、`state_dict` 和 `config`（两个文件路径）、`parameter_count`，以及 `output_heads`（排序后的输出头名称）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import download_pretrained_enformer

try:
    download_pretrained_enformer("", "enformer.pt", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **示例输出：**

```text
ConfigurationError name must be a Hugging Face model id such as EleutherAI/enformer-official-rough.
```

示例不会下载任何内容：模型名在使用网络之前就会被检查。安装依赖后的用法：

```python
from dnakit.predictions import download_pretrained_enformer

result = download_pretrained_enformer(
    "EleutherAI/enformer-official-rough", "enformer_official.pt"
)
print(result["parameter_count"], result["output_heads"])
```

## 13) `dnakit.predictions.enformer_embeddings()`

- **作用：** 用本地 state dict 计算 DNA 序列的 Enformer 嵌入，并写入 `.npy` 文件。
- **处理方法：** 构建 `Enformer(EnformerConfig(**config))`，在 CPU 上加载 state dict（`torch.load(..., weights_only=True)`），切换到 eval 模式，用 `enformer_pytorch.data.str_to_seq_indices()` 编码序列，再调用 `enformer_pytorch.finetune.get_enformer_embeddings()`，最后用 `numpy.save()` 保存数组。
- **API：** `enformer_embeddings(state_dict: 'str', sequences: 'Sequence[str]', output_npy: 'str', *, config: 'Mapping[str, Any]', freeze: 'bool' = False, train_layernorms_only: 'bool' = False, train_last_n_layers_only: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 enformer-pytorch 和 torch（extra `dnakit[neural-enformer]`）以及 numpy（extra `dnakit[bio-analysis]`）。权重已在磁盘上时可离线运行。函数在缺包时给出的安装提示写的是名为 `enformer` 的 extra，而 `pyproject.toml` 中实际的 extra 是 `neural-enformer`。
- **输入：** `state_dict` 为已有的权重文件，`config` 为非空映射，内容是训练这些权重时所用的 `EnformerConfig` 参数（可以读取 `download_pretrained_enformer()` 写出的 JSON 后传入）。`sequences` 为非空的字符串列表，只含 A、C、G、T、N，且长度一致。`output_npy` 为已有文件夹中的新文件路径。`freeze` 不能与 `train_layernorms_only` 或 `train_last_n_layers_only` 同时使用，后者不能为负。
- **输出：** `status`、`backend`、`output`（文件路径）、`shape`、`sequence_count`、`sequence_length`，以及三个训练选项。
- **示例代码：**

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

- **示例输出：**

```text
[2, 16, 96] 2 128
(2, 16, 96)
```

示例使用一个随机初始化的极小 Enformer，因此只打印形状。使用真实检查点时，请传入 `download_pretrained_enformer()` 保存的配置。

## 14) `dnakit.predictions.predict_selene_sequences()`

- **作用：** 在 CPU 上用完整的、训练好的 Selene 风格 torch 模型预测 DNA 序列。
- **处理方法：** 用 `torch.load(..., weights_only=False)` 加载文件，要求它是 `torch.nn.Module`，切换到 eval 模式，用 selene 的 `Genome.sequence_to_encoding()` 对序列做 one-hot 编码（N 变成均匀的 0.25 行），再对 (batch, length, 4) 数组调用 `selene_sdk.predict.model_predict.predict()`（不使用 CUDA）；selene 在调用模型前会把它转置为 (batch, 4, length)。
- **API：** `predict_selene_sequences(model_path: 'str', sequences: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`，它的依赖会安装 torch、numpy 和 matplotlib）。可离线运行。
- **输入：** `model_path` 为已有文件，由 `torch.save(model)` 写出完整模块；state dict 会被拒绝。文件会被反序列化，所以只能用可信文件。`sequences` 为一条 DNA 字符串，或最多 256 条、长度一致、不超过 10,000 个碱基的字符串列表。
- **输出：** `status`、`operation`、`model_path`、`sequence_count`、`length`、`shape`（batch, features）、`predictions`（嵌套列表）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
[3, 1] [[2.0], [0.0], [2.5]]
```

模型收到的是 (batch, 4, length) 张量，并统计每条序列中 A 的个数；第三条序列的 `N` 每个按 0.25 计，所以 `ANNA` 得 2.5。

## 15) `dnakit.predictions.run_selene_config_file()`

- **作用：** 运行 YAML 配置文件描述的 Selene 流程，并把输出写到输出文件夹下。
- **处理方法：** 用 selene 的 `load_path()` 加载 YAML（不实例化对象），把其中的 `output_dir` 键换成给定文件夹，再调用 `selene_sdk.utils.config_utils.parse_configs_and_run()`。这与 selene 的 `load_and_parse_configs_and_run()` 在读完文件后所做的相同，区别是运行文件夹中的配置副本是重新导出的 `configuration.yaml`，而不是原文的拷贝。selene 会把进度信息打印到标准输出。
- **API：** `run_selene_config_file(config_path: 'str', output_dir: 'str', *, create_subdirectory: 'bool' = True, lr: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`）。可离线运行。
- **输入：** `config_path` 为已有的 YAML 文件，顶层为映射，含最多 16 个名称的 `ops` 列表（取自 `train`、`evaluate`、`analyze`），以及指向已有模型文件的 `model: {path: ...}`。`output_dir` 为文件夹（不存在时创建）。`create_subdirectory` 为真时，运行结果放进新建的带时间戳的子文件夹。`lr` 为正数，仅当流程要训练时才需要。
- **输出：** `status`、`operation`、`source`（配置文件路径）、`operations`、`output_dir`、`run_directory`、`files`（运行文件夹中文件名的排序列表）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
['analyze'] ['configuration.yaml', 'model.py']
index	name	f1	f2
0	s1	5.00e-01	7.50e-01
1	s2	5.00e-01	7.50e-01
```

权重是固定的（权重全为 0，偏置为 0 和 ln 3），所以每条序列的两个特征分别是 sigmoid(0) = 0.5 和 sigmoid(ln 3) = 0.75。

## 16) `dnakit.predictions.run_selene_configs()`

- **作用：** 用解析后的配置字典运行 Selene 流程，并把输出写到输出文件夹下。
- **处理方法：** 复制字典，设置它的 `output_dir`，再调用 `selene_sdk.utils.config_utils.parse_configs_and_run()`。字典的结构与 selene 的 `load_path()` 返回的一致，因此用 `!obj:` 标签给出的对象在流程开始时仍未实例化。
- **API：** `run_selene_configs(configs: 'Mapping[str, Any]', output_dir: 'str', *, create_subdirectory: 'bool' = True, lr: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`）。可离线运行。
- **输入：** `configs` 为键是字符串的映射，含 `ops` 列表和 `model: {path: ...}`，同 `run_selene_config_file()`；调用者的对象不会被修改。`output_dir`、`create_subdirectory` 和 `lr` 同 `run_selene_config_file()`。当 `ops` 含 `train` 时，`configs` 要有 `lr`，或者传入 `lr` 参数。
- **输出：** 键与 `run_selene_config_file()` 相同，其中 `source` 为 `configs dict`。
- **示例代码：**

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

- **示例输出：**

```text
configs dict ['analyze'] ['configuration.yaml', 'model.py']
['seqs_predictions.tsv']
```

## 17) `dnakit.predictions.run_selene_operations()`

- **作用：** 按给定顺序运行 Selene 由配置定义的操作（`train`、`evaluate`、`analyze`）。
- **处理方法：** 用操作名、配置和输出文件夹调用 selene 的 `execute()`。与 `run_selene_configs()` 不同，它不会创建运行文件夹、复制模型或写配置文件；每个操作只读取它需要的配置段。
- **API：** `run_selene_operations(operations: 'Sequence[str]', configs: 'Mapping[str, Any]', output_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`）。可离线运行。
- **输入：** `operations` 为最多 16 个名称的列表，取自 `train`、`evaluate`、`analyze`。`configs` 为映射，含这些操作要读取的配置段（所有操作都需要 `model`，下面的示例还需要 `analyze_sequences` 和 `prediction`；列出 `train` 时还要有 `lr`）。`output_dir` 为文件夹（不存在时创建）。
- **输出：** `status`、`operation`、`operations`、`output_dir` 和 `provenance`。各操作写出的文件由操作自己报告，位于配置中指定的文件夹里。
- **示例代码：**

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

- **示例输出：**

```text
completed ['analyze']
['seqs_predictions.tsv']
```

## 18) `dnakit.predictions.selene_precision_recall_curves()`

- **作用：** 为每个模型特征写出一条精确率-召回率曲线，保存为 SVG 文件。
- **处理方法：** 按 selene 的规则选出特征，用 scikit-learn 计算每个特征的平均精确率，再调用 `selene_sdk.utils.performance_metrics.visualize_precision_recall_curves()`，它会在输出文件夹中写出 `precision_recall_curves.svg`。matplotlib 样式只在这次调用中生效，打开的图会被关闭。
- **API：** `selene_precision_recall_curves(prediction: 'Sequence[Sequence[float]]', target: 'Sequence[Sequence[int]]', output_dir: 'str', *, report_gt_feature_n_positives: 'int' = 50, style: 'str' = 'seaborn-v0_8-colorblind', fig_title: 'str' = 'Feature precision-recall curves', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`，它的依赖会安装 scikit-learn 和 matplotlib）。可离线运行。
- **输入：** `prediction` 和 `target` 为形状相同的 (样本, 特征) 表，至少两个样本；`target` 只含 0 和 1。只有两种标签都出现、且正例数大于 `report_gt_feature_n_positives` 的特征才会被画出。`output_dir` 为文件夹（不存在时创建），其中不能已有 `precision_recall_curves.svg`。`style` 为 matplotlib 样式名。
- **输出：** `status`、`operation`、`path`（SVG 文件）、`image_format`、`data`（`sample_count`、`feature_count`、`positives_per_feature`、`plotted_features`、`metric` = `average_precision`、`metric_per_plotted_feature`、`style`）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
precision_recall_curves.svg [0, 1]
[0.917, 1.0]
```

## 19) `dnakit.predictions.selene_roc_curves()`

- **作用：** 为每个模型特征写出一条 ROC 曲线，保存为 SVG 文件。
- **处理方法：** 按 selene 的规则选出特征，用 scikit-learn 计算每个特征的 ROC AUC，再调用 `selene_sdk.utils.performance_metrics.visualize_roc_curves()`，它会在输出文件夹中写出 `roc_curves.svg`。默认样式是 `seaborn-v0_8-colorblind`，因为 selene 自己的默认样式名在 matplotlib 中已不存在。
- **API：** `selene_roc_curves(prediction: 'Sequence[Sequence[float]]', target: 'Sequence[Sequence[int]]', output_dir: 'str', *, report_gt_feature_n_positives: 'int' = 50, style: 'str' = 'seaborn-v0_8-colorblind', fig_title: 'str' = 'Feature ROC curves', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`，它的依赖会安装 scikit-learn 和 matplotlib）。可离线运行。
- **输入：** 同 `selene_precision_recall_curves()`；文件为 `roc_curves.svg`。默认 `report_gt_feature_n_positives=50` 时，特征需要多于 50 个正例才会被画出，所以很小的表什么也画不出。
- **输出：** 键与 `selene_precision_recall_curves()` 相同，其中 `metric` = `roc_auc`。
- **示例代码：**

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

- **示例输出：**

```text
roc_curves.svg [3, 3]
[0, 1] [0.889, 1.0]
```

## 20) `dnakit.predictions.selene_score_heatmap()`

- **作用：** 把“位置 × 碱基”的计算机突变扫描（in silico mutagenesis）得分矩阵画成热图 PNG。
- **处理方法：** 在 matplotlib 坐标轴上调用 selene 的 `selene_sdk.interpret.vis.heatmap()`，序列类型为 selene 的 `Genome`。图中从上到下的行依次是碱基 T、G、C、A，与 selene 一致。
- **API：** `selene_score_heatmap(score_matrix: 'Sequence[Sequence[float]]', output_path: 'str', *, mask: 'Sequence[Sequence[bool]] | None' = None, cmap: 'str' = 'Blues', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`）。可离线运行。
- **输入：** `score_matrix` 最多 5,000 行（位置），每行四个数字，对应 A、C、G、T。`mask` 为可选的同形状网格，为 True 的格子会被隐藏。`cmap` 为 matplotlib 色图名。`output_path` 必须是新的 `.png` 路径。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`position_count`、`cmap`、`row_labels_top_to_bottom`、`score_matrix`、`mask`）和 `provenance`。
- **示例代码：**

```python
import os

from dnakit.predictions import selene_score_heatmap

scores = [[0.1, 0.2, 0.3, 0.4], [0.4, 0.1, 0.0, 0.2], [0.0, 0.5, 0.1, 0.1]]
mask = [[False, False, False, False], [False, True, False, False], [False, False, False, False]]
result = selene_score_heatmap(scores, "ism_heatmap.png", mask=mask, verbose=False)
print(os.path.basename(result["path"]), open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["position_count"], result["data"]["row_labels_top_to_bottom"])
```

- **示例输出：**

```text
ism_heatmap.png True
3 ['T', 'G', 'C', 'A']
```

## 21) `dnakit.predictions.plot_sequence_logo_scores()`

- **作用：** 根据“位置 × 碱基”的得分矩阵（例如归因得分）画 DNA 序列标志图并保存为 PNG。
- **处理方法：** 用普通的 matplotlib 图形块重现 selene 的 `sequence_logo()`，因为在当前安装的 matplotlib 下 selene 自己的调用会失败（它的字形路径效果用到了已被移除的方法）。每个位置的碱基按数值堆叠（`order="value"`，负值在零线以下），或按字母表顺序堆叠；每个条形都用 selene 自带的该碱基字形轮廓填充（`selene_sdk.interpret.vis._SVG_PATHS`）。
- **API：** `plot_sequence_logo_scores(score_matrix: 'Sequence[Sequence[float]]', image_path: 'str | Path', *, order: 'str' = 'value', width: 'float' = 1.0, color_scheme: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`）。可离线运行。
- **输入：** `score_matrix` 最多 5,000 行，每行四个有限数字，对应 A、C、G、T。`order` 为 `value` 或 `alpha`。`width` 为非负的条形宽度。`color_scheme` 为至少四种颜色的列表，依次对应 A、C、G、T；默认是 seaborn Set1 调色板的前四种颜色。`image_path` 必须是已有文件夹中的新 `.png` 路径。
- **输出：** `status`、`operation`、`image_path`、`image_format`、`order`、`bases` 和 `scores`（二者都是每个堆叠层一行、每个位置一列，按绘制顺序排列）以及 `provenance`。
- **示例代码：**

```python
from dnakit.predictions import plot_sequence_logo_scores

scores = [[0.2, 0.1, -0.1, 0.5], [0.0, 0.3, 0.1, 0.1], [-0.2, 0.4, 0.0, 0.2]]
result = plot_sequence_logo_scores(scores, "logo_scores.png", verbose=False)
print(open(result["image_path"], "rb").read(4) == b"\x89PNG", result["order"])
for bases, values in zip(result["bases"], result["scores"]):
    print("".join(bases), values)
```

- **示例输出：**

```text
True value
GAA [-0.1, 0.0, -0.2]
CGG [0.1, 0.1, 0.0]
ATT [0.2, 0.1, 0.2]
TCC [0.5, 0.3, 0.4]
```

打印出来的每一行是一个堆叠层，每一列是一个位置。在第一个位置上，堆叠从 G（-0.1，在零线以下）经 C（0.1）和 A（0.2）到 T（0.5）。

## 22) `dnakit.predictions.plot_variant_diffs_scatter()`

- **作用：** 把每个变异的最大绝对差值得分画成散点，并写入新的 HTML 文件。
- **处理方法：** 按染色体（selene 的 `sort_standard_chrs`）和位置给变异排序，在选定的特征列上取最大值，可选地只保留不低于某个百分位数的变异，再画出 plotly 散点，并用 `plotly.offline.plot()` 写出。selene 自己的函数在 plotly 6 下会失败，且在使用特征过滤时会把标签与得分错位配对；本函数改用 plotly 6 的 title 对象写布局，并按标签顺序取行和列。
- **API：** `plot_variant_diffs_scatter(scores: 'Sequence[Sequence[float]]', labels: 'Sequence[Sequence[Any]]', features: 'Sequence[str]', output_path: 'str | Path', *, feature_indices: 'Sequence[int] | None' = None, nth_percentile: 'float | None' = None, threshold_line: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 selene-sdk（extra `dnakit[bio-selene]`，它的依赖会安装 plotly 和 numpy）。可离线运行。
- **输入：** `scores` 为“变异 × 特征”的有限数字矩阵（最多 5,000 行），`features` 为各列的名称。`labels` 每行一条 `(chrom, pos, name, ref, alt)`。`feature_indices` 只保留这些列。`nth_percentile` 在 0 到 100 之间。`threshold_line` 为数字。`output_path` 必须是已有文件夹中的新 `.html` 路径。
- **输出：** `status`、`operation`、`output_path`、`output_format`、`features`（实际使用的特征名）、`variants`（按绘制顺序排列，含 `chrom`、`pos`、`ref`、`alt`、`max_abs_diff` 的对象）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
variant_scatter.html True
chr1 200 0.9
chr2 50 0.5
```

## 23) `dnakit.regulation.plot_attribution_logo()`

- **作用：** 根据“位置 × 碱基”的归因矩阵画序列标志图并保存为 PNG。
- **处理方法：** 对转置后的矩阵调用 tangermeme 的 `plot_logo()`（tangermeme 使用 (4, length) 布局），可带窗口、注释表和对称的 y 轴界限。界限在绘制后再设置到坐标轴上，因为没有注释时 tangermeme 会忽略它。
- **API：** `plot_attribution_logo(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, start: 'int | None' = None, end: 'int | None' = None, ylim: 'float | None' = None, annotations: 'Sequence[Mapping[str, Any]] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tangermeme 和 torch（extra `dnakit[bio-tangermeme]`）、matplotlib 和 pandas。可离线运行。
- **输入：** `attributions` 最多 10,000 行，每行四个有限数字，对应 A、C、G、T。`start` 和 `end` 选择窗口（`start` 小于 `end`，且都在序列范围内）。`ylim` 为正数界限，坐标轴范围为 `-ylim` 到 `ylim`。`annotations` 为最多 1,000 条记录的列表，含 `name`、`start`、`end`、`score`；tangermeme 只会画起点在窗口起点之后、终点在窗口终点之前的注释。`output_path` 必须是新的 `.png` 路径。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`position_count`、`start`、`end`、`ylim`（`[-ylim, ylim]` 或 `None`）、`alphabet`、`attributions`、`annotations`）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
True
5 [-1.0, 1.0] ['A', 'C', 'G', 'T']
```

## 24) `dnakit.regulation.plot_categorical_attribution_scatter()`

- **作用：** 绘制每个位置的归因总和，并按占优势的碱基着色，保存为 PNG。
- **处理方法：** 调用 tangermeme 的 `plot_categorical_scatter()`。y 值是四个碱基的总和，颜色取决于绝对值最大的碱基（A 绿、C 蓝、G 橙、T 红）。颜色是固定的，因为 tangermeme 的 `colors` 参数会失败。函数还会报告每个位置的优势碱基和总和。它适合长到标志图无法显示的序列。
- **API：** `plot_categorical_attribution_scatter(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, size: 'float' = 20.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tangermeme 和 torch（extra `dnakit[bio-tangermeme]`）和 matplotlib。可离线运行。
- **输入：** `attributions` 最多 10,000 行，每行四个有限数字，对应 A、C、G、T。`size` 为正的点大小。`output_path` 必须是新的 `.png` 路径。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`position_count`、`alphabet`、`dominant_base_index_by_position`（0 到 3 依次对应 A、C、G、T）、`summed_attribution_by_position`、`attributions`）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
True
[0, 1, 2, 3, 0]
[0.5, 0.8, -0.3, 0.9, 0.4]
```

## 25) `dnakit.regulation.plot_pwm_logo()`

- **作用：** 把位置权重矩阵画成信息量标志图并保存为 PNG。
- **处理方法：** 对转置后的矩阵调用 tangermeme 的 `plot_pwm()`，可带名称和对数保护项 `eps`。
- **API：** `plot_pwm_logo(pwm: 'Sequence[Sequence[float]]', output_path: 'str', *, name: 'str | None' = None, eps: 'float' = 1e-07, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tangermeme 和 torch（extra `dnakit[bio-tangermeme]`）、numpy 和 matplotlib。可离线运行。
- **输入：** `pwm` 最多 10,000 行，每行四个概率，对应 A、C、G、T；每行非负，且和为 1（误差不超过 1e-4）。`name` 为可选标题，最多 200 个字符。`eps` 为 0 到 1 之间的数。`output_path` 必须是新的 `.png` 路径。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`position_count`、`alphabet`、`name`、`pwm`）和 `provenance`。
- **示例代码：**

```python
from dnakit.regulation import plot_pwm_logo

pwm = [[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25], [0.0, 1.0, 0.0, 0.0]]
result = plot_pwm_logo(pwm, "pwm_logo.png", name="demo motif", verbose=False)
print(open(result["path"], "rb").read(4) == b"\x89PNG")
print(result["data"]["position_count"], result["data"]["name"])
```

- **示例输出：**

```text
True
3 demo motif
```

## 26) `dnakit.regulation.plot_sequence_attributions()`

- **作用：** 用 torch 模型计算 DNA 序列的 DeepLIFT-SHAP 归因，并为每条序列画一张标志图。
- **处理方法：** 对序列做 one-hot 编码，再调用 tangermeme 的 `plot_attributions()`，归因方法为 `deep_lift_shap`，参考序列为二核苷酸打乱序列（`dinucleotide_shuffle`，用 `random_state` 播种）。各标志图上下堆叠在一个 PNG 里。本函数不负责加载模型。
- **API：** `plot_sequence_attributions(model: 'Any', sequences: 'Sequence[str]', output_path: 'str', *, target: 'int' = 0, n_shuffles: 'int' = 20, random_state: 'int' = 0, hypothetical: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tangermeme 和 torch（extra `dnakit[bio-tangermeme]`）和 matplotlib。可离线运行。
- **输入：** `model` 为已处于 eval 模式的 `torch.nn.Module`。`sequences` 为 1 到 32 条 DNA 字符串的列表，长度一致（不超过 10,000 个碱基）。`target` 为要解释的模型输出的下标，非负。`n_shuffles` 在 1 到 200 之间，`random_state` 为非负整数，`hypothetical` 为布尔值。`output_path` 必须是新的 `.png` 路径。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`sequence_count`、`length`、`target`、`n_shuffles`、`random_state`、`hypothetical`、`sequences`，以及按（序列、位置、碱基）排列的 `attributions` 列表）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
True
2 5 2 5 4
```

## 27) `dnakit.regulation.render_interactive_logo_png()`

- **作用：** 把 tangermeme 的交互式归因标志图画成静态 PNG（悬停提示不会被保存）。
- **处理方法：** 在 matplotlib 坐标轴上调用 tangermeme 的 `interactive_logo()` 并保存图。tangermeme 的交互式标志图需要 mpld3；无法导入时，函数抛出 `BackendUnavailableError`，且不写任何文件。
- **API：** `render_interactive_logo_png(attributions: 'Sequence[Sequence[float]]', output_path: 'str', *, start: 'int | None' = None, end: 'int | None' = None, ylim: 'Sequence[float] | None' = None, annotations: 'Sequence[Mapping[str, Any]] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tangermeme、torch 和 mpld3，以及 matplotlib。`pyproject.toml` 中的 extra `dnakit[bio-tangermeme]` 列出了 mpld3，但函数的文档字符串和报错提示都说它不属于该 extra，出现这个报错时请用 `python -m pip install mpld3` 安装。检查本页时所用的环境没有安装 mpld3，所以示例展示的是这个报错。
- **输入：** 同 `plot_attribution_logo()`，区别是 `ylim` 为 `(low, high)` 对，且 `low` 小于 `high`。
- **输出：** `status`、`operation`、`path`、`image_format`、`data`（`position_count`、`start`、`end`、`ylim`、`alphabet`、`attributions`）和 `provenance`。
- **示例代码：**

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

- **示例输出：**

```text
BackendUnavailableError mpld3 could not be imported
```

上面的输出来自没有 mpld3 的环境。安装了 mpld3 时，示例会打印 `drawn 3`。安装依赖后的用法：

```python
from dnakit.regulation import render_interactive_logo_png

attributions = [[0.5, 0.0, 0.1, -0.1], [0.0, 0.7, 0.0, 0.1], [0.1, 0.0, -0.4, 0.0]]
result = render_interactive_logo_png(attributions, "interactive_logo.png", ylim=[-1.0, 1.0])
print(result["data"]["position_count"])
```

## 28) `dnakit.predictions.download_genomic_benchmark_cloud_cache()`

- **作用：** 从 Zenodo 云缓存下载预先构建好的完整序列基准数据集，并描述它。
- **处理方法：** 先检查（名称, 版本）是否在包的 `CLOUD_CACHE` 中，再调用 `genomic_benchmarks.loc2seq.cloud_caching.download_from_cloud_cache()`，它会下载一个 zip 压缩包并解压到 `dest_path` 的父文件夹。不在缓存中的名称会在建立任何连接之前被拒绝。使用 `force_download` 时，已有的 `dest_path` 会先被删除。
- **API：** `download_genomic_benchmark_cloud_cache(dataset_name: 'str', version: 'int', dest_path: 'str', *, force_download: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和网络连接。
- **输入：** `dataset_name` 与 `version` 必须是云缓存中存在的一对（例如 `demo_human_or_worm`，版本 0）。`dest_path` 为要创建的数据集文件夹，缺少的父文件夹会被创建。`force_download` 为布尔值。
- **输出：** `service`、`source`（Zenodo 记录的 URL）和 `result`（含 `dataset`、`local_path` 和 `summary`：`record_count`、`split_counts`、`label_counts`，以及含 `split`、`label`、`length`、`sequence` 的 `first_record`）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import download_genomic_benchmark_cloud_cache

try:
    download_genomic_benchmark_cloud_cache("not_a_dataset", 0, "full_data", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
```

- **示例输出：**

```text
ConfigurationError not_a_dataset version 0 is not in the genomic_benchmarks cloud cache.
```

示例不会下载任何内容：云缓存之外的名称会先被拒绝。安装依赖后的用法：

```python
from dnakit.predictions import download_genomic_benchmark_cloud_cache

result = download_genomic_benchmark_cloud_cache("demo_human_or_worm", 0, "data/demo_human_or_worm")
print(result["result"]["summary"]["record_count"])
```

## 29) `dnakit.predictions.download_genomic_benchmark_sequences()`

- **作用：** 构建区间列表型基准数据集的完整序列版本，并描述它。
- **处理方法：** 调用 `genomic_benchmarks.loc2seq.download_dataset()`：它读取数据集的 `metadata.yaml` 和区间表，从 `cache_path` 取参考基因组（缺失时下载），切出序列，并写出 `<dest_path>/<name>/<split>/<label>/<id>.txt`。启用云缓存时，如果该名称和版本有预先构建的副本就直接使用。随后函数按划分和标签统计文件数。
- **API：** `download_genomic_benchmark_sequences(interval_list_dataset: 'str', *, version: 'int | None' = None, dest_path: 'str | None' = None, cache_path: 'str | None' = None, force_download: 'bool' = False, use_cloud_cache: 'bool' = True, local_repo: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）。使用包内列出的基准名称，以及 `cache_path` 中还没有的参考基因组时，需要联网；下面的示例使用本地数据集文件夹和本地参考序列，所以可以离线运行。
- **输入：** `interval_list_dataset` 为基准名称，或在 `local_repo=True` 时为本地数据集文件夹的路径，其中含 `metadata.yaml` 和 `<split>/<class>.csv.gz` 文件。给出 `version` 时必须与元数据一致。`dest_path` 和 `cache_path` 默认位于 `~/.genomic_benchmarks` 下。`force_download` 强制重新下载参考序列。`use_cloud_cache` 和 `local_repo` 为布尔值。
- **输出：** `service`、`source`（数据集名称）和 `result`（含 `dataset`、`local_path` 和 `summary`：`record_count`、`split_counts`、按类别文件夹名统计的 `label_counts`、`first_record`）。
- **示例代码：**

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

- **示例输出：**

```text
mini_ds 4 {'train': 2, 'test': 2} {'pos': 4}
8 ACGTACGT
```

示例写入一条含两个序列的参考序列和一张含两行的区间表；第二个区间在负链上，所以它的序列是反向互补的。

## 30) `dnakit.predictions.get_genomic_benchmark_dataset()`

- **作用：** 加载一个命名的 genomic-benchmarks 分类数据集的某个划分，并描述它。
- **处理方法：** 调用 `genomic_benchmarks.dataset_getters.pytorch_datasets.get_dataset()`。如果 `~/.genomic_benchmarks/<dataset_name>` 不存在（或设置了 `force_download`），包会先下载数据集；`use_cloud_cache` 为真且名称在缓存中时，从 Zenodo 云缓存下载。函数统计记录数和标签数并返回第一条记录；不返回 PyTorch 数据集对象。标签是类别文件夹按文件系统列出顺序得到的下标，所以编号不保证按字母顺序。
- **API：** `get_genomic_benchmark_dataset(dataset_name: 'str', split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `dataset_name` 为非空的基准名称。`split` 为 `train` 或 `test`。`version` 为非负整数或 `None`。`force_download` 和 `use_cloud_cache` 为布尔值。
- **输出：** `service`、`source`（云缓存 URL；不在缓存中的名称则是 GitHub 数据集 URL）和 `result`（含 `dataset`、`split`、`local_path` 和 `summary`：`record_count`、`split_counts`、`label_counts`、`first_record`）。
- **示例代码：**

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

- **示例输出：**

```text
tiny_set train
4 {'train': 4} {'0': 2, '1': 2} 8
```

## 31) `dnakit.predictions.load_demo_coding_vs_intergenomic_seqs()`

- **作用：** 加载演示用的“编码序列与基因间序列”数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `DemoCodingVsIntergenomicSeqs` 获取函数（数据集 `demo_coding_vs_intergenomic_seqs`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_demo_coding_vs_intergenomic_seqs(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `demo_coding_vs_intergenomic_seqs`。
- **示例代码：**

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

- **示例输出：**

```text
demo_coding_vs_intergenomic_seqs https://zenodo.org/records/16605299
4 {'test': 4} {'0': 2, '1': 2}
```

## 32) `dnakit.predictions.load_demo_human_or_worm()`

- **作用：** 加载区分人与线虫序列的演示数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `DemoHumanOrWorm` 获取函数（数据集 `demo_human_or_worm`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_demo_human_or_worm(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `demo_human_or_worm`。
- **示例代码：**

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

- **示例输出：**

```text
demo_human_or_worm 6 {'0': 3, '1': 3}
8
```

## 33) `dnakit.predictions.load_demo_mouse_enhancers()`

- **作用：** 加载演示用的小鼠增强子数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `DemoMouseEnhancers` 获取函数。该数据集在上游以 `dummy_mouse_enhancers_ensembl` 存放，`~/.genomic_benchmarks` 下使用的就是这个文件夹名。汇总方式同 `get_genomic_benchmark_dataset()`。
- **API：** `load_demo_mouse_enhancers(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `dummy_mouse_enhancers_ensembl`。
- **示例代码：**

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

- **示例输出：**

```text
dummy_mouse_enhancers_ensembl 4 {'test': 4}
```

## 34) `dnakit.predictions.load_drosophila_enhancers_stark()`

- **作用：** 加载Stark 等人的果蝇增强子数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `DrosophilaEnhancersStark` 获取函数（数据集 `drosophila_enhancers_stark`）。这个上游函数没有 `use_cloud_cache` 参数，并且只要可能就总是使用云缓存，所以本函数也没有这个参数。汇总方式同 `get_genomic_benchmark_dataset()`。
- **API：** `load_drosophila_enhancers_stark(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download` 和 `version` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `drosophila_enhancers_stark`。
- **示例代码：**

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

- **示例输出：**

```text
drosophila_enhancers_stark 4 {'0': 2, '1': 2}
```

## 35) `dnakit.predictions.load_human_enhancers_cohn()`

- **作用：** 加载Cohn 等人的人类增强子数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `HumanEnhancersCohn` 获取函数（数据集 `human_enhancers_cohn`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_human_enhancers_cohn(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `human_enhancers_cohn`。
- **示例代码：**

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

- **示例输出：**

```text
human_enhancers_cohn 4 {'0': 2, '1': 2}
```

## 36) `dnakit.predictions.load_human_enhancers_ensembl()`

- **作用：** 加载由 Ensembl 调控区域构建的人类增强子数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `HumanEnhancersEnsembl` 获取函数（数据集 `human_enhancers_ensembl`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_human_enhancers_ensembl(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `human_enhancers_ensembl`。
- **示例代码：**

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

- **示例输出：**

```text
human_enhancers_ensembl 4 {'0': 2, '1': 2}
```

## 37) `dnakit.predictions.load_human_nontata_promoters()`

- **作用：** 加载人类无 TATA 框启动子数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `HumanNontataPromoters` 获取函数（数据集 `human_nontata_promoters`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_human_nontata_promoters(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `human_nontata_promoters`。
- **示例代码：**

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

- **示例输出：**

```text
human_nontata_promoters 4 {'0': 2, '1': 2}
```

## 38) `dnakit.predictions.load_human_ocr_ensembl()`

- **作用：** 加载由 Ensembl 构建的人类开放染色质区域数据集的某个划分，并描述它。
- **处理方法：** 调用包的 `HumanOcrEnsembl` 获取函数（数据集 `human_ocr_ensembl`），并像 `get_genomic_benchmark_dataset()` 那样汇总。数据集不在磁盘上时，会下载到 `~/.genomic_benchmarks`。
- **API：** `load_human_ocr_ensembl(split: 'str', *, force_download: 'bool' = False, version: 'int | None' = None, use_cloud_cache: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomic-benchmarks（extra `dnakit[bio-genomic-benchmarks]`）和 torch（该 extra 不会安装 torch）。除非数据集已在磁盘上，否则需要联网；下面的示例把一个很小的替身数据集放进一个充当主目录的文件夹，所以可以离线运行。
- **输入：** `split` 为 `train` 或 `test`。`force_download`、`version` 和 `use_cloud_cache` 的含义同 `get_genomic_benchmark_dataset()`。
- **输出：** 结构与 `get_genomic_benchmark_dataset()` 相同，其中 `dataset` 为 `human_ocr_ensembl`。
- **示例代码：**

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

- **示例输出：**

```text
human_ocr_ensembl 4 {'0': 2, '1': 2}
```
