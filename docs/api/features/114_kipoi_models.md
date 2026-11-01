# Kipoi：运行、验证和描述本地模型库中的模型

这一组 `dnakit.predictions` 中的函数是对 [kipoi](https://github.com/kipoi/kipoi)（基因组学模型库与 API）的适配器。它们加载 kipoi 的模型和数据加载器，描述并列出它们，把预测写成 TSV 文件，用模型自带的示例输入和期望预测来验证模型，并在当前 Python 环境或 Singularity 镜像中运行 `kipoi test`。几个小工具还暴露了 kipoi 的批处理、示例参数解析和预测比较。示例和大部分函数处理本地模型目录（`source="dir"`）；默认来源 `"kipoi"` 会克隆或更新 kipoi 基于 git 的模型库，因此需要联网。

kipoi 采用 MIT 许可（已在已安装包的元数据中核对）。这些函数用 `from dnakit.predictions import ...` 导入。依赖是可选的，只有调用函数时才会导入；缺少 kipoi 时函数会抛出 `BackendUnavailableError` 并给出安装命令 `python -m pip install kipoi`。Singularity 相关函数还会调用外部程序 `singularity`。

所有预测输出都写成 TSV，因为 kipoi 0.8.6 在 NumPy 2 下无法写 HDF5 文件。

## 示例数据

所有示例都使用同一个很小的本地模型，它同时也是示例的工作目录：一个返回 one-hot 编码序列 GC 比例的玩具模型，一个把每条 FASTA 记录前 20 个碱基做 one-hot 编码的数据加载器，以及一个两条记录的 FASTA 文件。需要模型库（装有多个模型文件夹的目录）的函数，其示例会把这些文件复制到临时目录的子文件夹中。

**示例文件：** `model.yaml`

```text
type: custom
args:
  file: model.py
  object: GcModel
default_dataloader: .
info:
  authors:
    - name: DNAKit docs
  doc: Toy model returning the GC fraction of one-hot encoded sequences.
  version: 0.1.0
  license: MIT
dependencies:
  conda:
    - python
  pip:
    - numpy
schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot DNA (A,C,G,T)
  targets:
    name: gc
    shape: (1,)
    doc: GC fraction
```

**示例文件：** `model.py`

```text
import numpy as np
from kipoi.model import BaseModel


class GcModel(BaseModel):
    def predict_on_batch(self, x):
        x = np.asarray(x)
        return x[:, :, 1:3].sum(axis=(1, 2)) / x.shape[1]
```

**示例文件：** `dataloader.yaml`

```text
type: Dataset
defined_as: dataloader.py::SeqDataset
args:
  fasta_file:
    doc: FASTA file
    example: example.fa
info:
  authors:
    - name: DNAKit docs
  doc: One-hot encodes each FASTA record
output_schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot
  metadata:
    name:
      type: str
      doc: record id
dependencies:
  pip:
    - numpy
```

**示例文件：** `dataloader.py`

```text
import numpy as np
from kipoi.data import Dataset


class SeqDataset(Dataset):
    def __init__(self, fasta_file):
        self.names, self.seqs = [], []
        for line in open(fasta_file):
            line = line.strip()
            if line.startswith(">"):
                self.names.append(line[1:].split()[0])
                self.seqs.append("")
            elif line:
                self.seqs[-1] += line.upper()

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        s = self.seqs[i][:20].ljust(20, "A")
        x = np.zeros((20, 4), dtype="float32")
        for j, c in enumerate(s):
            if c in "ACGT":
                x[j, "ACGT".index(c)] = 1
        return {"inputs": x, "metadata": {"name": self.names[i]}}
```

**示例文件：** `example.fa`

```text
>a
GGGGCCCCGGGGCCCCGGGG
>b
AAAAAAAAAATTTTTTTTTT
```

## 1) `dnakit.predictions.compare_kipoi_predictions()`

- **作用：** 用 kipoi 的比较方法，按给定的小数位数比较两个预测数组。
- **处理方法：** 把两个输入转成 NumPy 数组，调用 `kipoi.cli.main.compare_numpy_dict`（`kipoi test` 使用的比较）并设 `exact=False`。当处处满足 `abs(expected - observed) < 1.5 * 10**(-decimal)` 且形状一致时两个数组匹配。kipoi 抛出的 `AssertionError` 会转成 `matched=False` 和消息的第一行。
- **API：** `compare_kipoi_predictions(observed: 'Any', expected: 'Any', *, decimal: 'int' = 7, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** `observed` 和 `expected` 是数值数组或类数组；嵌套字典（object 数组）会以 `ConfigurationError` 拒绝。`decimal` 是不小于 0 的整数。
- **输出：** 字典，含 `matched`、`decimal`、`shape`、`expected_shape`、`mismatch`（kipoi 消息的第一行，或 `None`）和 `provenance`。
- **示例代码：**

```python
from dnakit.predictions import compare_kipoi_predictions

same = compare_kipoi_predictions([0.1, 0.2], [0.1, 0.2000001], decimal=5, verbose=False)
print(same["matched"], same["mismatch"])
different = compare_kipoi_predictions([0.1, 0.2], [0.1, 0.3], decimal=5, verbose=False)
print(different["matched"], different["mismatch"])
```

- **示例输出：**

```text
True None
False Arrays are not almost equal to 5 decimals
```

## 2) `dnakit.predictions.describe_kipoi_component()`

- **作用：** 读取本地目录中的 `model.yaml` 或 `dataloader.yaml` 并报告其描述。
- **处理方法：** 在目录上以 `which="model"` 或 `which="dataloader"` 调用 `kipoi.sources.load_component_descr`。目录按原样读取，不下载任何内容。
- **API：** `describe_kipoi_component(component_dir: 'str', *, which: 'str' = 'model', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** `component_dir` 必须是已存在的目录；`which` 为 `"model"` 或 `"dataloader"`，其他值会抛出 `ConfigurationError`。
- **输出：** `type`、`version`、`doc`、`license`，`inputs` 和 `targets` 模式（各含 `name`、`shape`、`doc`），数据加载器还有 `defined_as`，以及 `component_dir`、`which` 和 `provenance`。模型的模式来自 `schema`，数据加载器的来自 `output_schema`。
- **示例代码：**

```python
from dnakit.predictions import describe_kipoi_component

model = describe_kipoi_component(".", verbose=False)
print(model["type"], model["version"], model["license"])
print(model["inputs"]["name"], model["inputs"]["shape"], model["targets"]["name"])
loader = describe_kipoi_component(".", which="dataloader", verbose=False)
print(loader["type"], loader["defined_as"])
```

- **示例输出：**

```text
custom 0.1.0 MIT
seq [None, 4] gc
Dataset dataloader.py::SeqDataset
```

## 3) `dnakit.predictions.describe_kipoi_dataloader()`

- **作用：** 描述 kipoi 数据加载器：类型、定义它的类、输出模式和参数。
- **处理方法：** 调用 `kipoi.get_dataloader_descr`，并报告每个参数的文档、类型、是否可选和示例值。远程示例文件以 URL、md5 和名称报告，不会下载。
- **API：** `describe_kipoi_dataloader(dataloader: 'str', *, source: 'str' = 'kipoi', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；使用 `source="kipoi"` 时还需要联网来克隆或更新模型库，使用 `source="dir"` 则不需要。
- **输入：** `dataloader` 是来源中的名称，或在 `source="dir"` 时是相对于当前目录的目录。`source` 是 kipoi 来源名或 `"dir"`。
- **输出：** `dataloader`、`source`、`provenance`，以及 `description`，其中含 `type`、`defined_as`、`version`、`doc`、`license`、`inputs`、`targets` 和 `arguments` 映射。
- **示例代码：**

```python
from dnakit.predictions import describe_kipoi_dataloader

result = describe_kipoi_dataloader(".", source="dir", verbose=False)
description = result["description"]
print(description["type"], description["defined_as"])
print(description["inputs"]["shape"])
for name, argument in description["arguments"].items():
    print(name, argument["doc"], argument["optional"], argument["example"])
```

- **示例输出：**

```text
Dataset dataloader.py::SeqDataset
[None, 4]
fasta_file FASTA file False example.fa
```

## 4) `dnakit.predictions.get_kipoi_example_kwargs()`

- **作用：** 返回模型默认数据加载器的示例关键字参数，即 `kipoi get-example` 背后的信息。
- **处理方法：** 加载模型描述，找到其默认数据加载器的描述（以导入方式给出的会在模型目录内读取），并调用其 `download_example(output_dir=..., absolute_path=True, dry_run=...)`。远程示例文件会下载到 `output_dir`，除非 `dry_run=True`；像 `"example.fa"` 这样的普通示例值按 `dataloader.yaml` 中的写法返回。
- **API：** `get_kipoi_example_kwargs(model: 'str', *, source: 'str' = 'kipoi', output_dir: 'str' = 'example', dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；默认来源和远程示例文件需要联网，带普通示例值的本地 `source="dir"` 模型则不需要。
- **输入：** `model` 是模型名，或在 `source="dir"` 时是目录；`output_dir` 会被转成绝对路径；`dry_run` 必须是布尔值。
- **输出：** `model`、`source`、`output_dir`（绝对路径）、`dry_run`、`example_kwargs` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.predictions import get_kipoi_example_kwargs

result = get_kipoi_example_kwargs(".", source="dir", dry_run=True, verbose=False)
print(result["example_kwargs"], result["dry_run"], Path(result["output_dir"]).name)
```

- **示例输出：**

```text
{'fasta_file': 'example.fa'} True example
```

## 5) `dnakit.predictions.list_kipoi_directory_dataloaders()`

- **作用：** 列出本地 kipoi 模型库中的数据加载器。
- **处理方法：** 对由该目录构建的 `LocalSource` 调用 `kipoi.config.list_dataloaders`，因此不会访问任何已配置的网络来源。每行保留 `dataloader`、`version`、`type`、`doc`、`inputs`、`targets` 和 `license` 列。
- **API：** `list_kipoi_directory_dataloaders(directory: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** `directory` 必须是已存在的目录，其子文件夹中含有 `dataloader.yaml`。
- **输出：** `directory`、`total_dataloaders`、`dataloaders`（行字典列表）和 `provenance`。
- **示例代码：**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import list_kipoi_directory_dataloaders

zoo = Path(tempfile.mkdtemp())
(zoo / "gc").mkdir()
for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
    shutil.copy(name, zoo / "gc" / name)
result = list_kipoi_directory_dataloaders(str(zoo), verbose=False)
print(result["total_dataloaders"])
for row in result["dataloaders"]:
    print(row["dataloader"], row["type"], row["version"], row["inputs"])
```

- **示例输出：**

```text
1
gc Dataset 0.1 seq
```

## 6) `dnakit.predictions.list_kipoi_directory_models()`

- **作用：** 列出本地 kipoi 模型库中的模型。
- **处理方法：** 对由该目录构建的 `LocalSource` 调用 `kipoi.config.list_models`。每行保留 `model`、`version`、`type`、`doc`、`inputs`、`targets` 和 `license` 列。
- **API：** `list_kipoi_directory_models(directory: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** `directory` 必须是已存在的目录，其子文件夹中含有 `model.yaml`；空目录得到 0 个模型。
- **输出：** `directory`、`total_models`、`models`（行字典列表）和 `provenance`。
- **示例代码：**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import list_kipoi_directory_models

zoo = Path(tempfile.mkdtemp())
(zoo / "gc").mkdir()
for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
    shutil.copy(name, zoo / "gc" / name)
result = list_kipoi_directory_models(str(zoo), verbose=False)
print(result["total_models"])
for row in result["models"]:
    print(row["model"], row["type"], row["version"], row["license"])
```

- **示例输出：**

```text
1
gc custom 0.1.0 MIT
```

## 7) `dnakit.predictions.load_kipoi_dataloader()`

- **作用：** 返回 kipoi 数据加载器的类。
- **处理方法：** 调用 `kipoi.get_dataloader_factory`。用数据加载器参数实例化返回的类，并用 `batch_iter` 迭代。
- **API：** `load_kipoi_dataloader(dataloader: 'str', *, source: 'str' = 'kipoi', verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；默认来源需要联网，`source="dir"` 则不需要。
- **输入：** `dataloader` 是来源中的名称，或在 `source="dir"` 时是相对于当前目录的目录，例如 `"my_model/."`。未知的数据加载器会抛出 `BackendExecutionError`。
- **输出：** 数据加载器类（kipoi `Dataset` 的子类）。
- **示例代码：**

```python
from dnakit.predictions import load_kipoi_dataloader

loader_class = load_kipoi_dataloader(".", source="dir", verbose=False)
loader = loader_class(fasta_file="example.fa")
print(loader_class.__name__, len(loader))
for batch in loader.batch_iter(batch_size=2):
    print(sorted(batch), batch["inputs"].shape, batch["metadata"]["name"].tolist())
```

- **示例输出：**

```text
SeqDataset 2
['inputs', 'metadata'] (2, 20, 4) ['a', 'b']
```

## 8) `dnakit.predictions.load_kipoi_model()`

- **作用：** 加载 kipoi 模型及其默认数据加载器。
- **处理方法：** 用模型名、来源和 `with_dataloader` 调用 `kipoi.get_model`。返回的 kipoi 模型对象带有 `predict_on_batch`、`schema`、`info`、`default_dataloader` 以及 kipoi 附加的其他属性。
- **API：** `load_kipoi_model(model: 'str', *, source: 'str' = 'kipoi', with_dataloader: 'bool' = True, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；默认来源会克隆或更新 kipoi 的模型库，需要联网，`source="dir"` 不下载任何内容。
- **输入：** `model` 是模型名，或在 `source="dir"` 时是目录路径（相对于当前目录或绝对路径）。`with_dataloader` 必须是布尔值。
- **输出：** kipoi 模型对象。
- **示例代码：**

```python
import numpy as np

from dnakit.predictions import load_kipoi_model

model = load_kipoi_model(".", source="dir", verbose=False)
print(type(model).__name__, model.schema.inputs.shape)
x = np.zeros((2, 20, 4), dtype="float32")
x[0, :, 1] = 1
x[1, :, 0] = 1
print(model.predict_on_batch(x).tolist())
```

- **示例输出：**

```text
GcModel (None, 4)
[1.0, 0.0]
```

## 9) `dnakit.predictions.predict_kipoi_model()`

- **作用：** 用模型默认数据加载器的数据运行 kipoi 模型，并把预测写入 TSV 文件。
- **处理方法：** 执行 `kipoi predict` 的各步骤：加载模型，用校验过的 `dataloader_args` 实例化默认数据加载器，用 `batch_iter` 迭代，对每个批次调用 `predict_on_batch`，并在 `kipoi.cli.main.prepare_batch` 之后用 kipoi 的 TSV 写入器写出各批次。
- **API：** `predict_kipoi_model(model: 'str', output_file: 'str', *, source: 'str' = 'kipoi', dataloader_args: 'Mapping[str, Any] | None' = None, batch_size: 'int' = 32, num_workers: 'int' = 0, keep_inputs: 'bool' = False, keep_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；默认来源需要联网，`source="dir"` 则不需要。
- **输入：** `output_file` 必须以 `.tsv` 结尾、必须尚不存在，且其所在文件夹必须存在（拒绝 HDF5）。`dataloader_args` 把参数名映射到值；相对路径相对于当前目录。`batch_size` 至少为 1，`num_workers` 至少为 0；`keep_inputs` 和 `keep_metadata` 是布尔值。
- **输出：** `status`（文件写出时为 `completed`）、`model`、`source`、`dataloader`（类名）、`batches`、`samples`、`outputs`（写出的路径）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import predict_kipoi_model

output = Path(tempfile.mkdtemp()) / "predictions.tsv"
result = predict_kipoi_model(
    ".",
    str(output),
    source="dir",
    dataloader_args={"fasta_file": "example.fa"},
    verbose=False,
)
print(result["status"], result["dataloader"], result["batches"], result["samples"])
print(output.read_text().split())
```

- **示例输出：**

```text
completed SeqDataset 1 2
['preds', '1.0', '0.0']
```

## 10) `dnakit.predictions.validate_kipoi_model_example()`

- **作用：** 在本地模型目录的示例输入上运行模型，并可选地与期望预测比较，与 `kipoi test` 相同。
- **处理方法：** 以 `source="dir"` 加载模型，对数据加载器的示例参数调用其 `pipeline.predict_example`；当有期望时，用 `kipoi.readers.HDF5Reader` 读取 HDF5 文件，并按 `model.yaml` 中 `test.precision_decimal` 的精度，用 `compare_numpy_dict` 比较 `model.predict_on_batch(inputs)` 与存储的 `preds`。期望来自 `expect_file`，或 `model.yaml` 中本地的 `test.expect` 路径；远程的 `test.expect` 会在 `note` 中报告，不会下载。
- **API：** `validate_kipoi_model_example(model_dir: 'str', *, batch_size: 'int' = 32, output_file: 'str | None' = None, expect_file: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）（它会装上 h5py，读取 `expect_file` 时要用）。示例或期望是远程文件的模型需要联网来获取示例；这里的本地示例则不需要。
- **输入：** `model_dir` 必须是含 `model.yaml` 的已存在目录；`output_file`（可选）必须是新的 `.tsv` 路径；`expect_file`（可选）必须是已存在的 HDF5 文件，含 `inputs` 和 `preds` 数据集，即 `kipoi predict -o file.h5 --keep_inputs` 写出的格式，且只比较单个预测数组。
- **输出：** `status`（仅比较失败时为 `failed`）、`model_dir`、`samples`、`prediction_shape`、`outputs`、`expect_file`、`compared_samples`、`matched`（`True`、`False`，没有比较时为 `None`）、`precision_decimal`、`note` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import h5py
import numpy as np

from dnakit.predictions import validate_kipoi_model_example

plain = validate_kipoi_model_example(".", verbose=False)
print(plain["status"], plain["samples"], plain["prediction_shape"], plain["matched"])

inputs = np.zeros((2, 20, 4), dtype="float32")
inputs[0, :, 1] = 1
inputs[1, :, 0] = 1
expect = Path(tempfile.mkdtemp()) / "expected.h5"
with h5py.File(expect, "w") as handle:
    handle.create_dataset("inputs", data=inputs)
    handle.create_dataset("preds", data=np.array([1.0, 0.0]))
checked = validate_kipoi_model_example(".", expect_file=str(expect), verbose=False)
print(checked["status"], checked["matched"], checked["compared_samples"], checked["precision_decimal"])
with h5py.File(expect, "w") as handle:
    handle.create_dataset("inputs", data=inputs)
    handle.create_dataset("preds", data=np.array([0.0, 1.0]))
wrong = validate_kipoi_model_example(".", expect_file=str(expect), verbose=False)
print(wrong["status"], wrong["matched"])
```

- **示例输出：**

```text
completed 2 [2] None
completed True 2 7
failed False
```

## 11) `dnakit.predictions.predict_kipoi_model_in_singularity()`

- **作用：** 在显式指定的 Singularity 镜像中，对本地模型目录运行 `kipoi predict` 并写出 TSV 预测。
- **处理方法：** 构造命令 `singularity exec -B <dirs> <image> kipoi predict <model_dir> --source dir --batch_size <n> -o <output_file>`，给出参数时再加上以 JSON 表示的 `--dataloader_args`，并挂载模型目录、输出文件夹和本地文件参数所在的文件夹。命令通过 DNAKit 的外部命令执行器运行，超时 3600 秒；非零退出会抛出带日志尾部的 `BackendExecutionError`。`dry_run=True` 时只返回命令，不运行。
- **API：** `predict_kipoi_model_in_singularity(container_image: 'str', model_dir: 'str', output_file: 'str', *, dataloader_args: 'Mapping[str, Any] | None' = None, batch_size: 'int' = 32, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 真实运行需要 PATH 上的外部程序 `singularity`（或 Apptainer），以及含有 kipoi 的本地 `.sif` 镜像。缺少可执行文件会抛出 `BackendUnavailableError`。试运行既不需要该程序也不需要 kipoi；镜像路径必须是已存在的文件。
- **输入：** `container_image` 必须是已存在的文件，`model_dir` 是已存在的目录，`output_file` 是已存在文件夹中的新 `.tsv` 路径；`dataloader_args` 是映射或 `None`；`batch_size` 至少为 1。
- **输出：** `status`（`dry_run` 或 `completed`）、`operation`、`container_image`、`command`、`bind_directories`、`outputs`（真实运行时为写出的 TSV）、`return_code`、`log_tail` 和 `provenance`（真实运行时含 singularity 版本）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.predictions import predict_kipoi_model_in_singularity

folder = Path(tempfile.mkdtemp())
image = folder / "kipoi.sif"
image.write_text("placeholder image file")
result = predict_kipoi_model_in_singularity(
    str(image),
    ".",
    str(folder / "predictions.tsv"),
    dataloader_args={"fasta_file": "example.fa"},
    dry_run=True,
    verbose=False,
)
command = result["command"]
print(result["status"], command[:3], len(command))
print(command[5], command[6], command[-2:])
print(len(result["bind_directories"]), result["outputs"], result["return_code"])
```

- **示例输出：**

```text
dry_run ['singularity', 'exec', '-B'] 16
kipoi predict ['--dataloader_args', '{"fasta_file": "example.fa"}']
2 [] None
```

依赖安装后的用法（不执行）：

```python
from dnakit.predictions import predict_kipoi_model_in_singularity

result = predict_kipoi_model_in_singularity(
    "kipoi.sif",
    "my_model",
    "predictions.tsv",
    dataloader_args={"fasta_file": "reads.fa"},
)
print(result["status"], result["outputs"])
```

## 12) `dnakit.predictions.test_kipoi_model_in_singularity()`

- **作用：** 在显式指定的 Singularity 镜像中，对本地模型目录运行 `kipoi test`。
- **处理方法：** 构造并运行 `singularity exec -B <model dir> <image> kipoi test --batch_size <n> --source dir <model_dir>`，即 `kipoi test-source --singularity` 对每个模型运行的同一个 `kipoi test`。kipoi 自己的测试会把预测与本地的 `test.expect` 比较。运行超时 1800 秒。非零退出会抛出带日志尾部的 `BackendExecutionError`；`dry_run=True` 时改为返回命令。
- **API：** `test_kipoi_model_in_singularity(container_image: 'str', model_dir: 'str', *, batch_size: 'int' = 32, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 真实运行需要 PATH 上的外部程序 `singularity`（或 Apptainer），以及含有 kipoi 的本地 `.sif` 镜像。当前环境没有 singularity，所以下面的示例展示程序缺失时抛出的错误。
- **输入：** `container_image` 必须是已存在的文件，`model_dir` 是已存在的目录，`batch_size` 至少为 1，`dry_run` 是布尔值。
- **输出：** 与预测函数相同的字典：`status`、`operation`、`container_image`、`command`、`bind_directories`、`outputs`（空）、`return_code`、`log_tail` 和 `provenance`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.predictions import test_kipoi_model_in_singularity

image = Path(tempfile.mkdtemp()) / "kipoi.sif"
image.write_text("placeholder image file")
os.environ["PATH"] = ""
try:
    test_kipoi_model_in_singularity(str(image), ".", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "singularity executable was not found" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

依赖安装后的用法（不执行）：

```python
from dnakit.predictions import test_kipoi_model_in_singularity

result = test_kipoi_model_in_singularity("kipoi.sif", "my_model", batch_size=4)
print(result["status"], result["return_code"])
```

## 13) `dnakit.predictions.test_kipoi_model_in_current_env()`

- **作用：** 在当前 Python 环境中对本地模型目录运行 `kipoi test`，记录到警告时判为失败。
- **处理方法：** 以 `--source dir --batch_size <n>` 调用 kipoi 的 `cli_test`，并捕获标准输出和标准错误。预测 `dataloader.yaml` 中的示例输入，若 `model.yaml` 有本地的 `test.expect`，则按模型的精度比较。测试期间记录到的任何 `WARNING` 都会使调用抛出 `BackendExecutionError`，与 `kipoi.cli.source_test.test_model` 对其子进程的处理一致；kipoi 异常或非零退出同样会抛出。不会安装任何东西，也不会创建 conda 环境。
- **API：** `test_kipoi_model_in_current_env(model_dir: 'str', *, batch_size: 'int' = 32, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`），且模型自身的依赖必须已经可以导入。
- **输入：** `model_dir` 必须是已存在的目录；`batch_size` 至少为 1。
- **输出：** `status`、`operation`、`model_dir`、`batch_size`、`warnings`（成功时总是 0）、`log_tail`（捕获的输出，可能含进度条）和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.predictions import test_kipoi_model_in_current_env

result = test_kipoi_model_in_current_env(".", batch_size=2, verbose=False)
print(result["status"], result["batch_size"], result["warnings"], Path(result["model_dir"]).is_dir())
```

- **示例输出：**

```text
completed 2 0 True
```

## 14) `dnakit.predictions.prepare_kipoi_batch()`

- **作用：** 按 `kipoi predict` 的方式把模型预测附加到一个数据加载器批次上。
- **处理方法：** 把批次复制到新字典，并调用 `kipoi.cli.main.prepare_batch`。预测存放在 `"preds"` 下。除非 `keep_inputs` 为真，否则删除 `"inputs"` 和 `"targets"`；除非 `keep_metadata` 为真，否则删除 `"metadata"`。调用者的字典不会被修改。
- **API：** `prepare_kipoi_batch(dataloader_batch: 'dict[str, Any]', predictions: 'Any', *, keep_inputs: 'bool' = False, keep_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** `dataloader_batch` 是映射，例如 kipoi 数据加载器产出的一个批次（非映射会抛出 `ConfigurationError`）；`predictions` 是模型输出；`keep_inputs` 和 `keep_metadata` 是布尔值。
- **输出：** 新的字典，含 `preds`，并按标志位含 `inputs`、`targets` 和 `metadata`。
- **示例代码：**

```python
import numpy as np

from dnakit.predictions import prepare_kipoi_batch

batch = {
    "inputs": np.zeros((2, 3)),
    "targets": np.ones(2),
    "metadata": {"name": np.array(["a", "b"])},
}
predictions = np.array([0.5, 0.25])
print(sorted(prepare_kipoi_batch(batch, predictions, verbose=False)))
full = prepare_kipoi_batch(batch, predictions, keep_inputs=True, keep_metadata=True, verbose=False)
print(sorted(full), full["preds"].tolist())
print(sorted(batch))
```

- **示例输出：**

```text
['preds']
['inputs', 'metadata', 'preds', 'targets'] [0.5, 0.25]
['inputs', 'metadata', 'targets']
```

## 15) `dnakit.predictions.concatenate_kipoi_batches()`

- **作用：** 沿第一个轴按字段把数据加载器批次拼接起来。
- **处理方法：** 对批次列表调用 `kipoi.cli.main.numpy_collate_concat`：数组用 `numpy.concatenate` 拼接，嵌套字典或列表逐项拼接。结果的结构与一个批次相同，包含所有批次的全部样本。
- **API：** `concatenate_kipoi_batches(batches: 'Sequence[Any]', *, verbose: bool = True) -> 'Any'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）。
- **输入：** 非空的批次序列，且各批次结构相同。字符串、字节串或映射，或空序列，会抛出 `ConfigurationError`。
- **输出：** 一个结构相同的批次，按顺序包含所有输入批次的样本。
- **示例代码：**

```python
import numpy as np

from dnakit.predictions import concatenate_kipoi_batches

first = {"inputs": np.zeros((2, 3)), "metadata": {"name": np.array(["a", "b"])}}
second = {"inputs": np.ones((3, 3)), "metadata": {"name": np.array(["c", "d", "e"])}}
joined = concatenate_kipoi_batches([first, second], verbose=False)
print(joined["inputs"].shape, joined["inputs"][:, 0].tolist())
print(joined["metadata"]["name"].tolist())
```

- **示例输出：**

```text
(5, 3) [0.0, 0.0, 1.0, 1.0, 1.0]
['a', 'b', 'c', 'd', 'e']
```

## 16) `dnakit.predictions.kipoi_dataloader_example_kwargs()`

- **作用：** 解析本地 kipoi 数据加载器的示例参数。
- **处理方法：** 用 `kipoi.sources.load_component_descr` 加载 `dataloader.yaml`，并对其参数调用 `kipoi.specs.example_kwargs`。普通值按写法返回；远程示例文件变成 `cache_dir` 下的本地路径（缺失或 md5 不符时下载到那里，除非 `dry_run`），没有 `cache_dir` 时则是含 `url`、`md5` 和 `name` 的字典。示例中使用了 `dry_run` 或不给 `cache_dir`，所以不会下载任何内容。
- **API：** `kipoi_dataloader_example_kwargs(dataloader_dir: 'str', *, cache_dir: 'str | None' = None, absolute_path: 'bool' = True, dry_run: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；下载远程示例文件（给出 `cache_dir` 且不是 `dry_run`）需要联网。
- **输入：** `dataloader_dir` 必须是含 `dataloader.yaml` 的已存在目录；`cache_dir` 可选；`absolute_path` 和 `dry_run` 是布尔值。
- **输出：** `dataloader_dir`、`cache_dir`、`dry_run`、`example_kwargs` 和 `provenance`。给出 `cache_dir` 时，即使 `dry_run` 为真也会创建缓存文件夹。
- **示例代码：**

```python
import tempfile
import textwrap
from pathlib import Path

from dnakit.predictions import kipoi_dataloader_example_kwargs

print(kipoi_dataloader_example_kwargs(".", verbose=False)["example_kwargs"])

folder = Path(tempfile.mkdtemp())
(folder / "dataloader.py").write_text("from kipoi.data import Dataset\n\n\nclass Loader(Dataset):\n    pass\n")
(folder / "dataloader.yaml").write_text(textwrap.dedent("""\
    type: Dataset
    defined_as: dataloader.py::Loader
    args:
      fasta_file:
        doc: FASTA file
        example:
          url: https://example.org/data/example.fa
          md5: d41d8cd98f00b204e9800998ecf8427e
      window:
        doc: window size
        example: 10
    info:
      authors:
        - name: DNAKit docs
      doc: Remote example
    output_schema:
      inputs:
        name: seq
        shape: (None, 4)
        doc: one-hot
    """))
remote = kipoi_dataloader_example_kwargs(str(folder), verbose=False)["example_kwargs"]
print(remote["fasta_file"]["url"], remote["window"])
cache = folder / "cache"
planned = kipoi_dataloader_example_kwargs(str(folder), cache_dir=str(cache), dry_run=True, verbose=False)
print(Path(planned["example_kwargs"]["fasta_file"]).relative_to(cache), cache.exists())
```

- **示例输出：**

```text
{'fasta_file': 'example.fa'}
https://example.org/data/example.fa 10
fasta_file True
```

## 17) `dnakit.predictions.kipoi_dataloader_default_args()`

- **作用：** 解析本地 kipoi 数据加载器的默认参数。
- **处理方法：** 加载 `dataloader.yaml` 并调用 `kipoi.specs.download_default_args`，这是 kipoi 加载模型时执行的步骤。远程默认文件会下载到 `output_dir` 并替换为本地路径；其他非 `None` 的默认值按原样返回，已存在的文件路径变成绝对路径，布尔值保持为标志。数值、列表、元组和字典默认值绕过上游调用，因为 kipoi 0.8.6 会把每个非布尔默认值传给 `os.path.exists`，而它会把整数当作文件描述符。
- **API：** `kipoi_dataloader_default_args(dataloader_dir: 'str', output_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`）；远程默认文件需要联网，没有远程默认值的数据加载器不下载任何内容。
- **输入：** `dataloader_dir` 必须是含 `dataloader.yaml` 的已存在目录；`output_dir` 会被转成绝对路径。
- **输出：** `dataloader_dir`、`output_dir`、`overrides`（参数名到解析后默认值的映射）和 `provenance`。没有默认值的参数不会出现。注意数值、列表和字典默认值排在其他默认值之后。
- **示例代码：**

```python
import tempfile
import textwrap
from pathlib import Path

from dnakit.predictions import kipoi_dataloader_default_args

folder = Path(tempfile.mkdtemp())
(folder / "dataloader.py").write_text("from kipoi.data import Dataset\n\n\nclass Loader(Dataset):\n    pass\n")
(folder / "dataloader.yaml").write_text(textwrap.dedent("""\
    type: Dataset
    defined_as: dataloader.py::Loader
    args:
      fasta_file:
        doc: FASTA file
      window:
        doc: window size
        default: 20
      strict:
        doc: switch
        default: true
      label:
        doc: label
        default: x
    info:
      authors:
        - name: DNAKit docs
      doc: Defaults
    output_schema:
      inputs:
        name: seq
        shape: (None, 4)
        doc: one-hot
    """))
output = folder / "downloads"
result = kipoi_dataloader_default_args(str(folder), str(output), verbose=False)
print(dict(sorted(result["overrides"].items())), output.exists())
```

- **示例输出：**

```text
{'label': 'x', 'strict': True, 'window': 20} False
```

## 18) `dnakit.predictions.test_kipoi_source_models()`

- **作用：** 在当前 Python 环境中，对本地模型来源中的模型运行 `kipoi test`。
- **处理方法：** 列出 `source_dir` 的模型文件夹（或校验 `models` 中的名称），再对每个模型运行 `test_kipoi_model_in_current_env`，因此记录到 `WARNING` 就算失败。失败会收集在结果中，而不会中断运行，除非 `stop_on_failure` 为真。这相当于 `kipoi test-source` 中的预测部分；不使用 conda 环境和 git 范围。
- **API：** `test_kipoi_source_models(source_dir: 'str', *, models: 'Sequence[str] | None' = None, batch_size: 'int' = 4, stop_on_failure: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 kipoi（额外依赖 `dnakit[bio-kipoi]`），且模型自身的依赖必须已经可以导入。
- **输入：** `source_dir` 必须是已存在的目录；`models` 是可选的子文件夹名序列（绝对路径、`..` 和非文件夹会抛出 `ConfigurationError`，目录中没有模型时同样如此）；`batch_size` 至少为 1；`stop_on_failure` 是布尔值。
- **输出：** `status`（`passed` 或 `failed`）、`source_dir`、`batch_size`、`tested`、`passed`、`failed`（失败的名称）、`results`（每个模型一行，含 `model`、`status`、`warnings` 和 `error`）和 `provenance`。
- **示例代码：**

```python
import shutil
import tempfile
from pathlib import Path

from dnakit.predictions import test_kipoi_source_models

zoo = Path(tempfile.mkdtemp())
for folder in ["gc_a", "gc_b"]:
    (zoo / folder).mkdir()
    for name in ["model.yaml", "model.py", "dataloader.yaml", "dataloader.py", "example.fa"]:
        shutil.copy(name, zoo / folder / name)
result = test_kipoi_source_models(str(zoo), verbose=False)
print(result["status"], result["tested"], result["passed"], result["failed"])
for row in result["results"]:
    print(row["model"], row["status"], row["warnings"], row["error"])
```

- **示例输出：**

```text
passed 2 2 []
gc_a passed 0 None
gc_b passed 0 None
```
