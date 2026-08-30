# 工程化与扩展

统一管理 DNAKit 的计算后端、Python API、CLI、配置工作流、批量与并行计算、缓存、版本和复现信息。

## 1) `ENG-003` Python API

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

## 2) `ENG-004` CLI

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

## 3) `ENG-006` 批量计算

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

## 4) `ENG-007` 并行计算

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

## 5) `ENG-008` 分块与流式处理

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

## 6) `ENG-009` 缓存

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

## 7) `ENG-012` 错误与警告

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
