# 工程化与扩展

统一管理 DNAKit 的计算后端、Python API、CLI、配置工作流、批量与并行计算、缓存、版本和复现信息。

## 1) Python API {#1-eng-003-python-api}

- **作用：** 通过稳定的 Python 对象、函数和结果类型调用 DNAKit 功能，便于在脚本、Notebook、测试和其他软件中组合分析流程。
- **处理方法：** 导入对应模块，通过类型化输入直接调用函数并读取返回对象。
- **API：**`dnakit` 顶层及 `dnakit.<domain>` 模块为命名空间入口，无直接调用参数；各函数参数见对应功能条目。
- **输入：**对应函数要求的 `DNASequence`、`DNARecord`、`DNASet` 或配置对象。
- **输出：** 各函数返回对应类型的结果对象；本例的 `length_features` 返回 `LengthResult`。
- **示例代码：**

```python
import dnakit
from dnakit.descriptors import length_features

sequence = dnakit.normalize(" acgt ").sequence
assert sequence is not None
print(length_features(sequence).symbol_length)  # 4
```

- **示例结果：**

```text
4
```

## 2) CLI {#2-eng-004-cli}

- **作用：** 通过命令行参数运行标准化、描述、指纹、搜索、比较等常用功能，并输出结构化 JSON，适合 shell 流程和自动化任务。
- **处理方法：** 解析子命令参数，调用相应功能并输出 JSON。
- **API：**`dnakit COMMAND[必须] ARGS[可选]`；每个子命令的完整参数以 `dnakit COMMAND --help` 为准。
- **输入：**子命令及其序列、文件或配置参数。
- **输出：** CLI 在标准输出打印 JSON；本例包含 `distance`、`method` 和 `mismatches`。
- **示例代码：**

```bash
dnakit compare ACGT ACGA --method hamming
```

- **示例结果：**

```text
{"costs": {"substitution": 1.0}, "distance": 1.0, "dp_cells": null, "edit_path": null, "exceeded_max_distance": false, "iupac_matching": "literal", "left_id": null, "left_length": 4, "max_cells": null, "max_distance": null, "method": "hamming", "mismatches": [{"left_symbol": "T", "position": 3, "right_symbol": "A"}], "name": "hamming_distance", "right_id": null, "right_length": 4}
```

## 3) 批量计算 {#3-eng-006}

- **作用：** 对大量 DNA 记录重复调用同一计算，按输入顺序收集成功结果和逐条错误，避免单个失败中断整个批次。
- **处理方法：** 依次调用 operation，并按输入顺序收集结果或错误。
- **API：**`dnakit.batch.run_batch(records[必须], operation[必须], name[必须], config[可选], progress[可选])`、`dnakit.batch.iter_batch(records[必须], operation[必须], config[可选], progress[可选])`；`config` 使用 `dnakit.batch.BatchConfig`。
- **输入：** `DNA([...])` 或 `DNARecord` 迭代器、批处理 callable 和可选错误、seed、resume 配置。`DNA` 会按记录顺序迭代，不会把整个集合当成一条序列。
- **输出：** 返回 `BatchResult`，`items` 按输入顺序保存每条结果或错误。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.batch import run_batch

records = DNA(["A", "CC"])
result = run_batch(
    records,
    lambda record, context: record.sequence.symbol_length,
    name="length",
)
print([item.value for item in result.items])  # [1, 2]
```

- **示例结果：**

```text
[1, 2]
```

## 4) 并行计算 {#4-eng-007}

- **作用：** 使用受控工作线程并行执行独立记录任务，同时保持结果与输入顺序一致，并统一传播取消、进度和错误信息。
- **处理方法：** 使用线程池执行独立任务，以 jobs 和 max_in_flight 限制并发。
- **API：**`dnakit.batch.run_batch(records[必须], operation[必须], name[必须], config[可选], progress[可选])`；`config` 使用 `dnakit.batch.BatchConfig`，本项设置 `execution_mode="thread"`。
- **输入：** `DNA([...])` 或记录迭代器、callable、worker 数和可选 `max_in_flight`。
- **输出：** 返回 `BatchResult`，包含成功数量 `success_count` 和逐条结果 `items`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.batch import BatchConfig, run_batch

records = DNA(["ACGT"] * 4)
result = run_batch(
    records,
    lambda record, context: (record.id, context.seed),
    name="threaded",
    config=BatchConfig(seed=7, jobs=2, execution_mode="thread", max_in_flight=2),
)
print(result.success_count)  # 4
```

- **示例结果：**

```text
4
```

## 5) 分块与流式处理 {#5-eng-008}

- **作用：** 以迭代器方式分块读取、转换和写出大文件，使内存只保留当前批次，并返回进度事件，适合超出内存容量的数据。
- **处理方法：** 逐条读取迭代器，每达到 chunk_size 条就返回一个分块。
- **API：**`dnakit.read(source[必须], format[可选], config[可选])`、`dnakit.RecordSource(iterator[必须], close_callback[可选], source_name[可选], format[可选])`、`dnakit.io.iter_chunks(values[必须], chunk_size[必须])`、`dnakit convert input_path[必须] output_path[必须] --input-format[可选] --output-format[可选] --overwrite[可选] --progress/--no-progress[可选]`。
- **输入：**文件路径或迭代器；可选格式、压缩和 chunk size。
- **输出：** `iter_chunks` 返回分块迭代器，每个分块是最多 `chunk_size` 个元素的元组。
- **示例代码：**

```python
from dnakit.io import iter_chunks

for chunk in iter_chunks(range(5), chunk_size=2):
    print(chunk)
# (0, 1), (2, 3), (4,)
```

- **示例结果：**

```text
(0, 1)
(2, 3)
(4,)
```

## 6) 缓存 {#6-eng-009}

- **作用：** 根据输入内容、功能名称和参数生成缓存键，复用已有计算结果并校验完整性，减少重复的高成本计算。
- **处理方法：** 对输入与参数生成缓存键，读写 JSON 并校验缓存内容。
- **API：**`dnakit.cache.CacheKey(namespace[必须], digest[必须], schema_version[可选])`、`dnakit.cache.CacheKey.from_components(namespace[必须], components[必须], schema_version[可选])`、`dnakit.cache.JSONCache(root[必须], max_entry_bytes[可选])`、`dnakit.cache.JSONCache.get(key[必须])`、`dnakit.cache.JSONCache.put(key[必须], value[必须])`、`dnakit.cache.JSONCache.invalidate(key[必须])`、`dnakit.cache.JSONCache.clear(namespace[可选])`。
- **输入：**专用缓存目录、namespace、输入/参数/算法和后端版本组件。
- **输出：** `get` 返回缓存中的值或 `None`；`put`、`invalidate` 和 `clear` 更新缓存。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.cache import CacheKey, JSONCache

with TemporaryDirectory() as directory:
    cache = JSONCache(Path(directory) / "cache")
    key = CacheKey.from_components("length", {"sequence": "ACGT", "algorithm": "v1"})
    cache.put(key, {"value": 4})
    print(cache.get(key))  # {"value": 4}
```

- **示例结果：**

```text
{'value': 4}
```

## 7) 错误与警告 {#7-eng-012}

- **作用：** 使用结构化异常和 Issue 统一返回错误码、上下文、严重程度及修复提示，使 CLI、API 和报告能够一致解释失败或数据问题。
- **处理方法：** 通过结构化异常字段返回错误码、上下文和修复提示。
- **API：**`dnakit.exceptions.DNAKitError(message[必须], code[可选], context[可选], hint[可选])`、`dnakit.core.Issue(code[必须], severity[必须], message[必须], location[可选], details[可选])`；CLI 非零退出码无调用参数。
- **输入：**非法输入、错误配置或缺少后端的请求。
- **输出：** 错误时抛出 `DNAKitError`，提供错误码、上下文和修复提示；校验问题以 `Issue` 表示。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.exceptions import DNAKitError

try:
    DNASequence("AX")
except DNAKitError as error:
    print(error.code)     # INVALID_ALPHABET
    print(error.context)
```

- **示例结果：**

```text
INVALID_ALPHABET
{'alphabet': 'strict', 'part_index': 0, 'part_offset': 1, 'symbol': 'X'}
```

## 8) Agent 与 MCP {#8-eng-013-agent-mcp}

- **作用：** 将稳定公开功能转换成 Agent 可搜索、可描述和可调用的 MCP 工具，同时复用现有 DNAKit 实现。
- **处理方法：** 从工具注册表读取 Schema，校验参数后分发到对应函数。
- **API：**`dnakit-mcp`、`dnakit.tools.default_tool_registry()`、`dnakit.tools.create_server()`。
- **输入：**工具名称、符合自动生成 Schema 的 JSON 参数；文件写入、模型下载和外部程序还需显式授权。
- **输出：** 返回工具目录或 MCP 调用结果，包含工具名称、参数或执行结果。
- **示例代码：**

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry()
result = registry.execute(
    "dnakit.thermodynamics.molecular_weight",
    {"sequence": "ACGT"},
)
print(result["value_dalton"])
```

- **示例结果：**

```text
1173.84
```

Agent 客户端配置见 [大模型工具 → MCP](../../llm_tools/mcp.md)。

<span id="custom-llm"></span>

通用大模型调用的自定义提示词、文件输入、服务商配置和安全边界见 [大模型工具 → Python API](../../llm_tools/python_api.md#custom-llm)；[Agent](../../llm_tools/agent.md#agent-tool-calling) 页面介绍一句话驱动工具调用，无需安装 MCP。

<span id="bundled-script-example"></span>[执行随包脚本](#dnakit.workflows.run)；[坐标转换组合示例](../../faq/tutorial_integrations.md#engineering)。

<span id="report-example"></span>[生成 PDF 报告](#dnakit.workflows.write_wes_report)；[完整文件示例](../../faq/tutorial_integrations.md#engineering)。

<span id="mcp-validation"></span>[验证范围说明](../../faq/tutorial_integrations.md#engineering)。

<!-- BEGIN GENERATED API: encode-credentials -->

<span id="encode-credentials"></span>

## 9) 管理 ENCODE API 凭据 {#dnakit.encode.manage_credentials}

- **作用：** 配置、查看状态或清除 ENCODE 凭据。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.manage_credentials(action, access_key=None, secret_key=None, *, verbose=True)`
- **输入：** 必须提供 `action`。
- **输出：** 返回 JSON 文本，报告所选凭据操作的结果。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.manage_credentials(
    action='check',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，报告所选凭据操作的结果。

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.encode.manage_credentials-parameters)

<!-- END GENERATED API: encode-credentials -->

<!-- BEGIN GENERATED API: execution -->

<span id="execution"></span>

## 10) 读取已执行 Notebook 的单元格与输出 {#dnakit.databases.read_executed_notebook}

- **作用：** 读取已执行 Notebook 的单元格与输出。
- **处理方法：** 读取本地已执行 notebook 中保存的单元格输出，并按长度上限截取。
- **API：** `dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 准备已经执行并保存输出的 `analysis_executed.ipynb`，或替换为自己的文件路径。
- **输出：** 返回字典，包含 `notebook_path`, `total_cells`, `cells`, `search` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.read_executed_notebook(
    notebook_path='analysis_executed.ipynb',
    max_output_chars=2000,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "notebook_path": "<string>",
  "total_cells": "<integer>",
  "cells": [
    {
      "idx": "<integer>",
      "type": "<string>",
      "source": "<string>",
      "output": "<string>"
    }
  ],
  "search": "<string>",
  "matching_cells": [
    "<value>"
  ],
  "n_matches": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.databases.read_executed_notebook-parameters)

## 11) 定位并校验随包资源文件 {#dnakit.workflows.resource_path}

- **作用：** 取得已校验的随包资源路径。
- **处理方法：** 定位随包资源并校验校验和，返回本地路径。
- **API：** `dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)`
- **输入：** 必须提供 `skill_id`, `resource`。
- **输出：** pathlib.Path；指向通过哈希校验的本地资源。
- **示例代码：**

```python
import dnakit

result = dnakit.workflows.resource_path(
    skill_id='genomic-coordinates',
    resource='scripts/convert_coords.py',
    verbose=False,
)
print(result)
```

- **示例结果：** pathlib.Path；指向通过哈希校验的本地资源。

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.workflows.resource_path-parameters)

## 12) 执行随包工作流并获取退出码与输出 {#dnakit.workflows.run}

- **作用：** 传递原脚本参数，执行随包工作流并返回真实退出码和输出。
- **处理方法：** 使用指定解释器执行随包脚本，收集退出码、日志和来源记录。
- **API：** `dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)`
- **输入：** 必须提供 `skill_id`, `script`。
- **输出：** 字典包含 skill、script、returncode、output、upstream_url、upstream_revision、license；returncode 非零表示脚本失败。
- **示例代码：**

```python
import dnakit

result = dnakit.workflows.run(
    skill_id='genomic-coordinates',
    script='scripts/convert_coords.py',
    arguments=['--from', 'bed', '--to', 'gff', 'chr1', '999', '1000', '--format', 'json'],
    timeout_seconds=60,
    verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

- **示例结果：** 字典包含 skill、script、returncode、output、upstream_url、upstream_revision、license；returncode 非零表示脚本失败。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "skill": "<str>",
  "script": "<str>",
  "returncode": "<int>",
  "output": "<str>",
  "upstream_url": "<str>"
}
```

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.workflows.run-parameters)

## 13) 执行组合流程的支持脚本 {#dnakit.workflows.run_support}

- **作用：** 执行组合流程依赖的随包支持脚本。
- **处理方法：** 执行指定辅助脚本，收集退出码与日志。
- **API：** `dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)`
- **输入：** 必须提供 `name`。
- **输出：** 字典包含 name、returncode、output；须检查退出码。
- **示例代码：**

```python
import dnakit

result = dnakit.workflows.run_support(
    name='gwas-lookup',
    arguments=['--help'],
    verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

- **示例结果：** 字典包含 name、returncode、output；须检查退出码。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "name": "<str>",
  "returncode": "<int>",
  "output": "<str>"
}
```

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.workflows.run_support-parameters)

## 14) 列出随包 Python、R 与 Shell 脚本 {#dnakit.workflows.scripts}

- **作用：** 列出技能包中可运行的 Python、R 和 Shell 脚本。
- **处理方法：** 读取随包资源目录，列出指定 Skill 的脚本路径。
- **API：** `dnakit.workflows.scripts(skill_id, *, verbose=True)`
- **输入：** 必须提供 `skill_id`。
- **输出：** 技能包中的脚本相对路径列表；可能为空。
- **示例代码：**

```python
import dnakit

result = dnakit.workflows.scripts(
    skill_id='genomic-coordinates',
    verbose=False,
)
print(result)
```

- **示例结果：** 技能包中的脚本相对路径列表；可能为空。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<relative script path>"
]
```

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.workflows.scripts-parameters)

## 15) 将已有 Markdown 报告渲染为 PDF {#dnakit.workflows.write_wes_report}

- **作用：** 把已有 Markdown 报告渲染为 PDF，并记录来源哈希。
- **处理方法：** 将准备好的 Markdown 排版为 PDF，并写入来源清单。
- **API：** `dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)`
- **输入：** 必须提供 `input_report`, `output_dir`。 先准备 `report.md` 并安装 `dnakit[reports]`；结果将写入 `reports/pdf`。
- **输出：** 返回报告路径和来源记录的字典；输出 wes_report.pdf 与 wes_report.manifest.json。需要 dnakit[reports]。
- **示例代码：**

```python
import dnakit

result = dnakit.workflows.write_wes_report(
    input_report='report.md',
    output_dir='reports/pdf',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回报告路径和来源记录的字典；输出 wes_report.pdf 与 wes_report.manifest.json。需要 dnakit[reports]。

[完整参数与使用说明](../../faq/api/14_engineering.md#dnakit.workflows.write_wes_report-parameters)

<!-- END GENERATED API: execution -->

## 16) 按依赖执行 Python 工作流 {#bio-59}

- **作用：** 串联 Python 函数并保留各步骤的状态和结果。
- **计算方法：** 检查依赖图并按拓扑顺序执行；失败步骤可阻断下游，独立步骤可继续。
- **API：** `dnakit.bioinformatics.workflows.Step(name, function, dependencies={}, parameters={})`; `dnakit.bioinformatics.workflows.run_workflow(steps, fail_fast=True)`
- **输入：** 命名步骤、实际 Python 函数和参数名→上游步骤名映射。
- **输出：** 返回按依赖顺序执行的工作流结果，包含各步骤状态和产物。
- **示例代码：**

```python
from dnakit.bioinformatics.workflows import Step, run_workflow

steps = [
    Step("sequence", lambda: "ACGT"),
    Step("length", lambda sequence: len(sequence), dependencies={"sequence": "sequence"}),
]
results = run_workflow(steps)
print([(step.name, step.status, step.value) for step in results])
```

- **示例结果：**

```text
[('sequence', 'completed', 'ACGT'), ('length', 'completed', 4)]
```

[依赖与实现说明](../../faq/bioinformatics.md#workflows)。
