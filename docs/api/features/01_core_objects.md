# 1. 核心数据对象

使用 `dnakit.DNA(...)` 创建一条或多条 DNA 序列。

## 1) DNA序列对象 {#1-core-001-dna}

- **作用：** 保存 DNA 序列值、字母表、线性或环状状态、单链或双链状态以及 Gap 信息，作为描述符、搜索、比对等计算的基础序列对象。
- **计算方法：** 标准化字符串后保存序列、拓扑和链类型；通过同一 DNA 对象访问记录与序列。
- **API：** `dnakit.DNA(data[必须], alphabet[可选], topology[可选], strandedness[可选])`。
- **输入：** `DNA` 可直接接收原始字符串并自动标准化；高级 `DNASequence` 只接收已经标准化的序列。
- **输出：** `DNA(...)` 返回 `DNASequence`，包含序列符号和长度。
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

## 2) DNA记录对象 {#2-core-002-dna}

- **作用：** 在一条 DNA 序列上附加 ID、描述、功能区和样本信息，使计算结果能够追溯到具体记录及其来源。
- **计算方法：** 将 ID、描述、区间注释和 metadata 绑定到同一条序列。
- **API：** `dnakit.DNA(data[必须], id[可选], description[可选], features[可选], metadata[可选], letter_annotations[可选])`。
- **输入：** 普通入口只需序列；其余信息都是同一调用中的可选参数。未提供 ID 时自动生成 `sequence_1`。
- **输出：** `DNA(...)` 返回 `DNARecord`，包含 ID、序列、注释和 metadata。
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

## 3) DNA数据集对象 {#3-core-003-dna}

- **作用：** 使用同一个 `DNA` 对象按固定顺序管理一条或多条记录，并支持下标和切片选择，作为普通用户统一的数据输入入口。
- **计算方法：** 按输入顺序构建记录集合；整数下标和切片均返回 DNA。
- **API：** `dnakit.DNA(data[必须], name[可选], source[可选], version[可选], collection_metadata[可选])`；文件读取使用 `dnakit.read(..., mode="dna")`。
- **输入：** 简单多序列使用字符串列表；每条记录需要不同信息时使用字典列表。
- **输出：** `DNA(...)` 返回 `DNASet`，包含有序记录和记录 ID。
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

## 4) 特征对象 {#4-core-004}

- **作用：** feature 是 `DNA` 中一条序列的附加注释，可标记 ORF、motif、限制酶位点和重复序列，并在编辑、导出和可视化时保留区域含义；它不是另一种 DNA 数据对象。
- **计算方法：** 把 feature 字典转换为带类型和 0-based 半开区间的注释。
- **API：** 在 `dnakit.DNA(..., features=[...])` 中直接传字典；简单区间使用 `type`、`start`、`end`，其余字段可选。
- **输入：** 普通 feature 字典必填 `type`、`start`、`end`；也可用 `location` 代替 `start/end`。
- **输出：** 返回包含 `DNAFeature` 的 `DNARecord`；特征提供类型、区间和链方向。
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

## 5) Gap对象 {#5-core-005-gap}

- **作用：** 明确保存序列中已知或未知长度的缺口，防止缺失区域在坐标计算或序列拼接时被误当成连续碱基。
- **计算方法：** 用独立 Gap 对象记录已知长度或未知长度，不将其自动替换为普通碱基。
- **API：** `dnakit.DNA(parts[必须], ...)`、`dnakit.Gap(length[必须], kind[可选], crossable[可选], evidence[可选], metadata[可选])`、`dnakit.GapKind`；`DNASequence.from_fragments(...)` 为高级兼容入口。
- **输入：** 必填为正整数长度，或用 `None` 表示未知长度；可选 `kind`、`crossable`、`evidence`、`metadata`。
- **输出：** `Gap` 表示缺口长度；包含它的 `DNASequence` 分别提供符号长度和坐标跨度。
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

## 6) 序列类型声明 {#6-core-006}

- **作用：** 声明序列允许的字符、线性或环状形态，以及单链或双链类型，供输入校验和后续算法选择正确的处理规则。
- **计算方法：** 用显式字母表、拓扑和链类型约束对象的有效状态。
- **API：** `dnakit.DNA(data[必须], alphabet[可选], topology[可选], strandedness[可选])`；`DNASequence`、`DNAAlphabet`、`Topology` 和 `Strandedness` 保留为高级类型。
- **输入：** 必填 DNA 内容；可选 `alphabet="strict"|"iupac"`、`topology="linear"|"circular"`、`strandedness="single"|"double"`。
- **输出：** 返回带 `alphabet`、`topology` 和 `strandedness` 声明的 `DNASequence`。
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

## 7) 坐标系统 {#7-core-007}

- **作用：** 把不同文件格式的坐标统一转换为 DNAKit 使用的 0-based 半开区间，避免区间截取和格式转换时出现一位偏差。
- **计算方法：** 按来源坐标体系换算端点；内部使用包含 start、不包含 end 的区间。
- **API：** `dnakit.core.ExternalInterval(start[必须], end[必须], system[必须], strand[可选])`、`dnakit.Interval(start[必须], end[必须])`、`dnakit.core.CompoundLocation(parts[必须])`、`dnakit.core.import_location(external[必须], sequence_length[可选])`、`dnakit.core.export_location(location[必须], target_system[必须], sequence_length[可选])`、`dnakit.core.reverse_strand_location(location[必须], sequence_length[必须])`。
- **输入：** 必填起点、终点、来源或目标坐标体系；可选 `strand`、`sequence_length`。跨环状原点时必须给出序列长度。
- **输出：** 返回坐标区间对象；转换后提供起点 `start` 和终点 `end`。
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

[FAQ：补充说明](../../faq/tutorial_basics.md#core-objects)

[FAQ：补充说明](../../faq/tutorial_basics.md#compatibility-01_core_objects)
