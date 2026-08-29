# 后端、性能与可复现性

统一管理 DNAKit 的计算后端、Python API、CLI、配置工作流、批量与并行计算、缓存、版本和复现信息。

外部工具被注册不等于已经具备对应的科学计算执行器。

## 1) `ENG-001` 统一后端接口

- **作用：** 通过统一注册表发现、检查和调用 Primer3、NUPACK 等可选后端，返回版本、能力和一致的错误信息，避免各模块自行探测外部工具。
- **API**：`dnakit.backends.backend_registry.probe(backend_id[必须])`、`dnakit.backends.backend_registry.load(backend_id[必须], capability[可选])`、`dnakit.backends.BackendRegistry()`；领域级 Primer3 和 Dashing adapter 另见对应功能页。
- **输入**：后端 ID；可选 capability、显式可执行文件路径和资源限制。
- **示例代码**：

```python
from dnakit.backends import backend_registry

print(tuple(backend_registry))
primer3 = backend_registry.probe("primer3-cli")
print(primer3.available, primer3.version)
```

- **示例结果：**

```text
('blastn', 'dashing', 'mmseqs2', 'nupack', 'primer3-cli', 'sourmash')
False None
```

## 2) `ENG-002` 原生与外部标记

- **作用：** 为每个结果记录实现属于原生、外部适配、重新实现或新增方法，同时保存算法、后端、来源和许可信息，使计算过程可追溯。
- **API**：`dnakit.core.ImplementationInfo(label[可选], execution_mode[可选], origin_class[可选], license_expression[可选], citations[可选])`、`dnakit.core.Provenance(dnakit_version[可选], python_version[可选], platform[可选], dependency_versions[可选], implementation[可选], backend[可选], reference[可选])`、`dnakit.core.BackendInfo(name[必须], version[可选], executable_path[可选], package_location[可选], license_expression[可选], capabilities[可选], available[可选], metadata[可选])`。
- **输入**：实现标签；可选执行模式、来源、许可、引用和后端信息。
- **示例代码**：

```python
from dnakit.core import ImplementationInfo, Provenance

implementation = ImplementationInfo(
    label="adapter",
    execution_mode="hybrid",
    origin_class="integration",
    license_expression="GPL-2.0-or-later",
)
provenance = Provenance(implementation=implementation)
print(provenance.implementation.label)  # adapter
```

- **示例结果：**

```text
adapter
```

## 3) `ENG-003` Python API

- **作用：** 通过稳定的 Python 对象、函数和结果类型调用 DNAKit 功能，便于在脚本、Notebook、测试和其他软件中组合分析流程。
- **API**：`dnakit` 顶层及 `dnakit.<domain>` 模块为命名空间入口，无直接调用参数；各函数参数见对应功能条目。
- **输入**：对应函数要求的 `DNASequence`、`DNARecord`、`DNASet` 或配置对象。
- **示例代码**：

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

## 4) `ENG-004` CLI

- **作用：** 通过命令行参数运行标准化、描述、指纹、搜索、比较等常用功能，并输出结构化 JSON，适合 shell 流程和自动化任务。
- **API**：`dnakit COMMAND[必须] ARGS[可选]`；每个子命令的完整参数以 `dnakit COMMAND --help` 为准。
- **输入**：子命令及其序列、文件或配置参数。
- **示例代码**：

```bash
dnakit describe ACGTACGT
dnakit compare ACGT ACGA --method hamming
```

- **示例结果：**

```text
{"base_composition": {"ambiguity_policy": "ignore", "counts": {"A": 2, "C": 2, "G": 2, "T": 2}, "cross_gaps": false, "denominator": 8, "fractions": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "gap_count": 0, "ignored_ambiguity_count": 0, "method": "canonical_base_count", "name": "base_composition", "sequence_id": null, "unknown_gap_count": 0}, "complexity": {"ambiguity_policy": "ignore", "by_k": {"1": 1.0, "2": 0.5714285714285714, "3": 0.6666666666666666, "4": 0.8, "5": 1.0, "6": 1.0}, "cross_gaps": false, "formula": "product_k(unique_kmers/min(4**k,valid_kmer_positions))", "gap_count": 0, "max_observations": 10000000, "max_word_size": 6, "method": "vocabulary-observed-over-possible-product", "name": "linguistic_complexity", "observation_count": 33, "observed_by_k": {"1": 4, "2": 4, "3": 4, "4": 4, "5": 4, "6": 3}, "possible_by_k": {"1": 4, "2": 7, "3": 6, "4": 5, "5": 4, "6": 3}, "score": 0.3047619047619048, "sequence_id": null, "unknown_gap_count": 0}, "gc_at": {"ambiguity_policy": "ignore", "at_count": 4, "at_fraction": 0.5, "cross_gaps": false, "denominator": 8, "gap_count": 0, "gc_count": 4, "gc_fraction": 0.5, "ignored_ambiguity_count": 0, "method": "canonical_base_fraction", "name": "gc_at_content", "sequence_id": null, "unknown_gap_count": 0}, "repeat": {"ambiguity_policy": "ignore", "comparisons": 4, "cross_gaps": false, "denominator": 8, "gap_count": 0, "max_comparisons": 5000000, "max_unit_length": 20, "method": "maximal-exact-tandem-repeat-union", "min_repeats": 2, "min_unit_length": 1, "name": "exact_repeat_fraction", "repeat_count_by_unit": {"4": 1}, "repeat_fraction": 1.0, "repeated_base_count": 8, "runs": [{"repeat_count": 2, "symbol_end": 8, "symbol_start": 0, "unit": "ACGT", "unit_length": 4}], "sequence_id": null, "unknown_gap_count": 0}}
{"costs": {"substitution": 1.0}, "distance": 1.0, "dp_cells": null, "edit_path": null, "exceeded_max_distance": false, "iupac_matching": "literal", "left_id": null, "left_length": 4, "max_cells": null, "max_distance": null, "method": "hamming", "mismatches": [{"left_symbol": "T", "position": 3, "right_symbol": "A"}], "name": "hamming_distance", "right_id": null, "right_length": 4}
```

## 5) `ENG-005` 配置工作流

- **作用：** 按 YAML/JSON 中声明的输入、操作和依赖顺序执行多步骤流程，保存每步结果、错误和运行清单，便于复现批量分析。
- **API**：`dnakit workflow config_path[必须] --dry-run[可选] --resume[可选] --progress/--no-progress[可选]`、`dnakit.workflows.load_workflow(path[必须])`、`dnakit.workflows.run_workflow(path[必须], dry_run[可选], resume[可选], progress[可选])`。
- **输入**：`dnakit-workflow-v1` 配置和本地输入文件。
- **示例代码**：

```bash
DNAKIT_WORK_DIR="$(mktemp -d)"
cp examples/fixed_demo.fasta examples/advanced_workflow.yml "$DNAKIT_WORK_DIR/"
dnakit workflow "$DNAKIT_WORK_DIR/advanced_workflow.yml" --no-progress
```

- **示例结果：**

```text
{
  "status": "succeeded",
  "run_id": "fixed-demo",
  "artifacts": [
    "splits/train.fasta",
    "splits/test.fasta",
    "report.html",
    "run-manifest.json"
  ]
}
```

## 6) `ENG-006` 批量计算

- **作用：** 对大量 DNA 记录重复调用同一计算，按输入顺序收集成功结果和逐条错误，避免单个失败中断整个批次。
- **API**：`dnakit.batch.run_batch(records[必须], operation[必须], name[必须], config[可选], progress[可选])`、`dnakit.batch.iter_batch(records[必须], operation[必须], config[可选], progress[可选])`；`config` 使用 `dnakit.batch.BatchConfig`。
- **输入**：记录迭代器、批处理 callable 和可选错误、seed、resume 配置。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence
from dnakit.batch import run_batch

records = [
    DNARecord(DNASequence("A"), "a"),
    DNARecord(DNASequence("CC"), "b"),
]
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

## 7) `ENG-007` 并行计算

- **作用：** 使用受控工作线程并行执行独立记录任务，同时保持结果与输入顺序一致，并统一传播取消、进度和错误信息。
- **API**：`dnakit.batch.run_batch(records[必须], operation[必须], name[必须], config[可选], progress[可选])`；`config` 使用 `dnakit.batch.BatchConfig`，本项设置 `execution_mode="thread"`。
- **输入**：记录、callable、worker 数和可选 `max_in_flight`。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence
from dnakit.batch import BatchConfig, run_batch

records = [DNARecord(DNASequence("ACGT"), str(index)) for index in range(4)]
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

## 8) `ENG-008` 分块与流式处理

- **作用：** 以迭代器方式分块读取、转换和写出大文件，使内存只保留当前批次，并返回进度事件，适合超出内存容量的数据。
- **API**：`dnakit.read(source[必须], format[可选], config[可选])`、`dnakit.RecordSource(iterator[必须], close_callback[可选], source_name[可选], format[可选])`、`dnakit.io.iter_chunks(values[必须], chunk_size[必须])`、`dnakit convert input_path[必须] output_path[必须] --input-format[可选] --output-format[可选] --overwrite[可选] --progress/--no-progress[可选]`。
- **输入**：文件路径或迭代器；可选格式、压缩和 chunk size。
- **示例代码**：

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

## 9) `ENG-009` 缓存

- **作用：** 根据输入内容、功能名称和参数生成缓存键，复用已有计算结果并校验完整性，减少重复的高成本计算。
- **API**：`dnakit.cache.CacheKey(namespace[必须], digest[必须], schema_version[可选])`、`dnakit.cache.CacheKey.from_components(namespace[必须], components[必须], schema_version[可选])`、`dnakit.cache.JSONCache(root[必须], max_entry_bytes[可选])`、`dnakit.cache.JSONCache.get(key[必须])`、`dnakit.cache.JSONCache.put(key[必须], value[必须])`、`dnakit.cache.JSONCache.invalidate(key[必须])`、`dnakit.cache.JSONCache.clear(namespace[可选])`。
- **输入**：专用缓存目录、namespace、输入/参数/算法和后端版本组件。
- **示例代码**：

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

## 10) `ENG-010` 随机种子

- **作用：** 统一管理随机 seed 和稳定记录顺序，使随机划分、抽样、突变和聚类在相同配置下可重复，并把 seed 写入结果。
- **API**：`dnakit.datasets.split(records[必须], config[必须])`、`dnakit.batch.run_batch(records[必须], operation[必须], name[必须], config[可选], progress[可选])`；随机性由对应 `SplitConfig.seed[可选]`、`BatchConfig.seed[可选]` 等配置字段控制。
- **输入**：操作输入和整数 seed。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.datasets import SplitConfig, split

records = DNASet([DNARecord(DNASequence("A"), str(index)) for index in range(4)])
config = SplitConfig(ratios={"train": 0.5, "test": 0.5}, seed=7)
first = split(records, config=config)
second = split(records, config=config)
print(first.assignments == second.assignments)  # True
```

- **示例结果：**

```text
True
```

## 11) `ENG-011` 版本追踪

- **作用：** 生成运行清单，记录 DNAKit、Python、依赖、后端、参考库、输入及输出文件的版本与校验值，用于复现和审计一次分析。
- **API**：`dnakit.core.Provenance(dnakit_version[可选], python_version[可选], platform[可选], dependency_versions[可选], implementation[可选], backend[可选], reference[可选])`、`dnakit.core.BackendInfo(name[必须], version[可选], executable_path[可选], package_location[可选], license_expression[可选], capabilities[可选], available[可选], metadata[可选])`、`dnakit.evaluation.ReferenceLibrary(records[必须], name[必须], version[必须], source[必须], digest[必须], digest_scope[必须], date[必须], filters[必须], index_parameters[必须])`、`dnakit.core.RunManifest(run_id[必须], command[必须], resolved_config[必须], provenance[必须], started_at[必须], status[必须], seed[可选], seed_derivation[可选], inputs[可选], outputs[可选], finished_at[可选], issues[可选])`。
- **输入**：当前运行上下文；可选依赖版本和后端信息。
- **示例代码**：

```python
from dnakit.core import Provenance

provenance = Provenance(dependency_versions={"example-backend": "1.2.3"})
print(provenance.dnakit_version)
print(provenance.python_version)
print(provenance.platform)
```

- **示例结果：**

```text
0.1.1
3.10.16
Linux-6.18.33.2-microsoft-standard-WSL2-x86_64-with-glibc2.39
```

## 12) `ENG-012` 错误与警告

- **作用：** 使用结构化异常和 Issue 统一返回错误码、上下文、严重程度及修复提示，使 CLI、API 和报告能够一致解释失败或数据问题。
- **API**：`dnakit.exceptions.DNAKitError(message[必须], code[可选], context[可选], hint[可选])`、`dnakit.core.Issue(code[必须], severity[必须], message[必须], location[可选], details[可选])`；CLI 非零退出码无调用参数。
- **输入**：非法输入、错误配置或缺少后端的请求。
- **示例代码**：

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

## 13) `ENG-013` 单元测试

- **作用：** 使用自动化单元和集成测试验证 API 在正常、边界及错误输入下的返回值和异常，防止代码修改破坏既有行为。
- **API**：`python -m pytest paths[可选] options[可选]`；严格默认配置来自 `pyproject.toml`。
- **输入**：仓库测试集和当前开发环境。
- **示例代码**：

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src:. python -m pytest -q
```

- **示例结果：**

```text
945 passed, 1 skipped
```

## 14) `ENG-014` 一致性验证

- **作用：** 使用人工可核对的固定预期值、性质测试和可比外部工具结果验证算法正确性，并记录差异、容差和环境边界。
- **API**：`python -m validation.run_validation --output[可选] --skip-optional[可选] --force[可选]`。
- **输入**：固定验证夹具、对照依赖及新的报告路径。
- **示例代码**：

```bash
DNAKIT_VALIDATION_DIR="$(mktemp -d)"
PYTHONNOUSERSITE=1 python -m validation.run_validation \
  --output "$DNAKIT_VALIDATION_DIR/report.json"
```

- **示例结果：**

```text
Validation report written: .../report.json
```

## 15) `ENG-015` 性能 benchmark

- **作用：** 在固定输入、seed 和任务配置下测量运行时间及内存分配峰值，输出可审计报告，用于发现性能回退而非宣称普遍速度。
- **API**：`python -m benchmarks.benchmark_core --sizes[可选] --repeats[可选] --warmups[可选] --seed[可选] --tasks[可选] --implementations[可选] --output[可选] --force[可选]`。
- **输入**：规模、任务、实现、重复、预热、seed 和新报告路径。
- **示例代码**：

```bash
DNAKIT_BENCH_DIR="$(mktemp -d)"
PYTHONNOUSERSITE=1 python -m benchmarks.benchmark_core \
  --sizes 100,1000 \
  --repeats 1 \
  --warmups 1 \
  --seed 20260813 \
  --tasks construct,normalize,gc_content \
  --implementations dnakit \
  --output "$DNAKIT_BENCH_DIR/report.json"
```

- **示例结果：**

```text
Benchmark report written: .../report.json
```

## 16) `ENG-016` 文档与教程

- **作用：** 集中提供安装、教程、API、Notebook、workflow、FAQ、算法边界和验证证据，使用户能够找到正确入口并理解结果限制。
- **API**：`mkdocs serve options[可选]`、`mkdocs build options[可选]`；两个 Notebook 使用固定输入，无额外调用参数。
- **输入**：仓库文档源文件及文档开发依赖。
- **示例代码**：

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src:. mkdocs serve
```

- **示例结果：**

```text
INFO    -  [15:21:49] Serving on http://127.0.0.1:8000/
```

## 17) `ENG-017` 可选图形入口

- **作用：** 定义 Notebook、Web 或 Galaxy 等图形化入口可复用的 API 和结果边界；当前仅说明扩展位置，不表示已经提供完整 GUI。
- **API**：当前没有接受用户 DNA 输入的图形计算 API，因此没有可填写参数；只提供固定输入、预生成结果的静态页面。
- **输入**：当前无可执行输入接口。
- **示例代码**：只能预览静态文档，不能提交序列计算。

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src:. mkdocs serve
```

- **示例结果：**

```text
INFO    -  [15:21:49] Serving on http://127.0.0.1:8000/
```
