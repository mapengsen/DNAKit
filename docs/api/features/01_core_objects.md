# 1. 核心数据对象

普通用户只需要记住一个核心对象：`dnakit.DNA(...)`。输入一条序列或多条序列、ID、拓扑、metadata 和 feature 都使用这个入口，返回类型始终是 `DNA`。

```python
import dnakit

# 一条序列
dna = dnakit.DNA("ACGT")

# 一条序列，并附加信息
annotated = dnakit.DNA(
    "ACGT",
    id="seq-1",
    topology="circular",
    metadata={"species": "synthetic"},
    features=[{"type": "motif", "start": 1, "end": 3}],
)

# 多条序列，仍然返回 DNA
dataset = dnakit.DNA(["ACGT", "TTAA"])
detailed = dnakit.DNA(
    [
        {"sequence": "ACGT", "id": "seq-1"},
        {"sequence": "TTAA", "id": "seq-2", "topology": "circular"},
    ]
)

print(type(dna).__name__, type(dataset).__name__)
print(dataset.ids, dataset[0].symbols)
```

```text
DNA DNA
('sequence_1', 'sequence_2') ACGT
```

使用规则只有三条：字符串表示一条序列；字符串列表表示多条序列；多条序列需要不同附加信息时，用包含 `sequence` 的字典列表。`data[0]` 和 `data[1:3]` 仍返回 `DNA`。下面的 `DNASequence`、`DNARecord` 和 `DNASet` 是内部明确分层及旧代码兼容对象，普通使用不需要分别构造。

## 1) CORE-001 DNA序列对象

- **作用：** 保存 DNA 序列值、字母表、线性或环状状态、单链或双链状态以及 Gap 信息，作为描述符、搜索、比对等计算的基础序列对象。
- **普通 API：** `dnakit.DNA(data[必须], alphabet[可选], topology[可选], strandedness[可选])`。
- **高级兼容对象：** `dnakit.DNASequence(parts[必须], alphabet[可选], topology[可选], strandedness[可选])`。
- **输入：** `DNA` 可直接接收原始字符串并自动标准化；高级 `DNASequence` 只接收已经标准化的序列。
- **示例代码：**

```python
import dnakit

seq = dnakit.DNA("ACGT", topology="linear", strandedness="double")
print(seq.symbols)
print(seq.symbol_length)
print(seq.record_count)
```

- **示例结果：**

```text
ACGT
4
1
```

- **限制：** `len(dna)` 表示记录数，不表示碱基数；碱基数使用 `symbol_length`。多记录对象读取 `symbols` 等单序列属性时会报错，应先用下标选一条。高级 `DNASequence` 仍只接受已规范化大写 DNA。

## 2) CORE-002 DNA记录对象

- **作用：** 在一条 DNA 序列上附加 ID、描述、功能区和样本信息，使计算结果能够追溯到具体记录及其来源。
- **普通 API：** `dnakit.DNA(data[必须], id[可选], description[可选], features[可选], metadata[可选], letter_annotations[可选])`。
- **高级兼容对象：** `dnakit.DNARecord(sequence[必须], id[必须], ...)`。
- **输入：** 普通入口只需序列；其余信息都是同一调用中的可选参数。未提供 ID 时自动生成 `sequence_1`。
- **示例代码：**

```python
import dnakit

record = dnakit.DNA(
    "ACGT",
    id="seq-1",
    description="示例序列",
    metadata={"species": "human"},
)
print(record.id, record.symbols, record.metadata["species"])
```

- **示例结果：**

```text
seq-1 ACGT human
```

- **限制：** 显式 ID 必须非空；feature 不得超出可解析坐标；metadata 必须可转为 JSON；每个 `letter_annotations` 数组必须是有限数值，且长度等于 `symbol_length`。

## 3) CORE-003 DNA数据集对象

- **作用：** 使用同一个 `DNA` 对象按固定顺序管理一条或多条记录，并支持下标和切片选择，作为普通用户统一的数据输入入口。
- **普通 API：** `dnakit.DNA(data[必须], name[可选], source[可选], version[可选], collection_metadata[可选])`；文件读取使用 `dnakit.read(..., mode="dna")`。
- **高级兼容对象：** `dnakit.DNASet(...)`、`DNASet.from_records(...)`、`DNASet.from_sequences(...)` 和 `read_set(...)`。
- **输入：** 简单多序列使用字符串列表；每条记录需要不同信息时使用字典列表。
- **示例代码：**

```python
import dnakit

dataset = dnakit.DNA(
    ["AC", "GT"],
    name="demo",
)
print(dataset.ids)
print(dataset[1].symbols)
```

- **示例结果：**

```text
('sequence_1', 'sequence_2')
GT
```

- **限制：** `DNA` 会在 `max_records` 上限内物化输入；字符串列表表示多条序列。单条 multipart 序列应使用 `{"parts": [...]}`，含显式 `Gap` 的片段列表也会被识别为一条 gapped 序列。筛选、去重、聚类和划分由相应领域模块提供。

## 4) CORE-004 特征对象

- **作用：** feature 是 `DNA` 中一条序列的附加注释，可标记 ORF、motif、限制酶位点和重复序列，并在编辑、导出和可视化时保留区域含义；它不是另一种 DNA 数据对象。
- **普通 API：** 在 `dnakit.DNA(..., features=[...])` 中直接传字典；简单区间使用 `type`、`start`、`end`，其余字段可选。
- **高级兼容对象：** `dnakit.DNAFeature(...)`、`dnakit.Interval(...)`；复合或未解析位置使用 `CompoundLocation`、`UnresolvedLocation`。
- **输入：** 普通 feature 字典必填 `type`、`start`、`end`；也可用 `location` 代替 `start/end`。
- **示例代码：**

```python
import dnakit

dna = dnakit.DNA(
    "ACGT",
    id="seq-1",
    features=[
        {
            "type": "motif",
            "start": 1,
            "end": 3,
            "id": "m1",
            "strand": "forward",
            "label": "示例位点",
        }
    ],
)
feature = dna.features[0]
print(feature.type, feature.location, feature.strand.value)
```

- **示例结果：**

```text
motif Interval(start=1, end=3) forward
```

- **限制：** 内部位置使用 0-based 半开区间；`phase` 仅允许 `0`、`1`、`2` 或 `None`。feature 属于对应记录，不能超出可解析的序列坐标。

## 5) CORE-005 Gap对象

- **作用：** 明确保存序列中已知或未知长度的缺口，防止缺失区域在坐标计算或序列拼接时被误当成连续碱基。
- **API：** `dnakit.DNA(parts[必须], ...)`、`dnakit.Gap(length[必须], kind[可选], crossable[可选], evidence[可选], metadata[可选])`、`dnakit.GapKind`；`DNASequence.from_fragments(...)` 为高级兼容入口。
- **输入：** 必填为正整数长度，或用 `None` 表示未知长度；可选 `kind`、`crossable`、`evidence`、`metadata`。
- **示例代码：**

```python
import dnakit

gap = dnakit.Gap(500, kind="scaffold", crossable=False, evidence=("paired-ends",))
seq = dnakit.DNA(["AC", gap, "GT"])
print(seq.symbol_length)
print(seq.coordinate_span)
```

- **示例结果：**

```text
4
504
```

- **限制：** 长度必须为正整数或 `None`；上下游片段由 `DNASequence.parts` 的位置表达。未知长度 Gap 的总坐标跨度为 `None`。

## 6) CORE-006 序列类型声明

- **作用：** 声明序列允许的字符、线性或环状形态，以及单链或双链类型，供输入校验和后续算法选择正确的处理规则。
- **API：** `dnakit.DNA(data[必须], alphabet[可选], topology[可选], strandedness[可选])`；`DNASequence`、`DNAAlphabet`、`Topology` 和 `Strandedness` 保留为高级类型。
- **输入：** 必填 DNA 内容；可选 `alphabet="strict"|"iupac"`、`topology="linear"|"circular"`、`strandedness="single"|"double"`。
- **示例代码：**

```python
import dnakit

seq = dnakit.DNA(
    "ACGN",
    alphabet="iupac",
    topology="circular",
    strandedness="double",
)
print(seq.alphabet.value, seq.topology.value, seq.strandedness.value)
```

- **示例结果：**

```text
iupac circular double
```

- **限制：** strict 仅允许 `A/C/G/T`；IUPAC 允许标准模糊符号。gapped 类型由 `parts` 中是否存在显式 `Gap` 派生，不是单独的字符串开关。构造时会核对内容与声明，且空序列不能声明为环状。

## 7) CORE-007 坐标系统

- **作用：** 把不同文件格式的坐标统一转换为 DNAKit 使用的 0-based 半开区间，避免区间截取和格式转换时出现一位偏差。
- **API：** `dnakit.core.ExternalInterval(start[必须], end[必须], system[必须], strand[可选])`、`dnakit.Interval(start[必须], end[必须])`、`dnakit.core.CompoundLocation(parts[必须])`、`dnakit.core.import_location(external[必须], sequence_length[可选])`、`dnakit.core.export_location(location[必须], target_system[必须], sequence_length[可选])`、`dnakit.core.reverse_strand_location(location[必须], sequence_length[必须])`。
- **输入：** 必填起点、终点、来源或目标坐标体系；可选 `strand`、`sequence_length`。跨环状原点时必须给出序列长度。
- **示例代码：**

```python
from dnakit.core import ExternalInterval, export_location, import_location

external = ExternalInterval(2, 8, system="1-based-closed", strand="forward")
internal = import_location(external)
(converted,) = export_location(internal, target_system="0-based-half-open")
print(internal)
print(converted.start, converted.end)
```

- **示例结果：**

```text
Interval(start=1, end=8)
1 8
```

- **限制：** 内部坐标固定为 0-based 半开；跨原点会返回 `CompoundLocation`。`import_location()` 不保留外部 `strand`，`export_location()` 固定返回 `Strand.UNKNOWN`，链信息应另存于 `DNAFeature.strand` 等字段。正负链的几何反转需显式调用 `reverse_strand_location(..., sequence_length=...)`；未解析位置不能导出，坐标转换不会替用户推断序列长度。
