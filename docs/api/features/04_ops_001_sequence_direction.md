# 序列补全

对 DNA 序列进行逆序、互补或反向互补，从而得到方向转换后的新序列。

## 1) 逆序 {#1-ops-0011}

- **作用：** 只把 DNA 序列的字符顺序倒过来而不替换碱基，用于检查或构造相反的位置顺序。
- **处理方法：** 从末尾向前读取碱基，保留原字符。
- **API：** `dnakit.ops.reverse(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **输出：** 返回处理后的 `DNASequence`，序列内容可从 `symbols` 读取。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **示例结果：**

```text
CGCTAA
```

## 2) 互补 {#2-ops-0012}

- **作用：** 按碱基配对规则把每个 DNA 碱基替换为互补碱基，但保持原有位置顺序，用于构造对应互补链。
- **处理方法：** 按 DNA 碱基互补表逐位替换，位置顺序不变。
- **API：** `dnakit.ops.complement(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **输出：** 返回处理后的 `DNASequence`，序列内容可从 `symbols` 读取。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **示例结果：**

```text
TTAGCG
```

## 3) 反向互补 {#3-ops-0013}

- **作用：** 同时反转序列顺序并替换为互补碱基，得到另一条链的 5′→3′ 表示，可用于反向链搜索、方向统一和引物分析。
- **处理方法：** 反转碱基顺序后取互补，同步映射方向、feature 和逐碱基注释。
- **API：** `dnakit.ops.reverse_complement(value, output_path=None, *, feature_policy="preserve", read_config=None, write_config=None, sequence_column=None, show_progress=True)`。
- **输入：** DNA 文本、文件路径，或单记录 `DNA`、`DNASequence`、`DNARecord`。文本返回字符串，文件返回 `WriteResult`；输入和输出路径必须不同，默认不覆盖已有输出。[FAQ：详细说明](../../faq/tutorial_integrations.md#reverse-complement)。
- **输出：** 返回反向互补后的 `DNASequence` 或 `DNA`，文件输入可写出到指定路径。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **示例结果：**

```text
GCGATT
```

<!-- BEGIN GENERATED API: sequence-direction -->

<span id="sequence-direction"></span>

## 4) 计算 DNA 序列的反向互补 {#dnakit.databases.sequence_reverse_complement}

- **作用：** 计算 DNA 序列的反向互补。
- **处理方法：** 按 A↔T、C↔G、N↔N 互补后反转顺序。
- **API：** `dnakit.databases.sequence_reverse_complement(*, operation=None, sequence, verbose=True)`
- **输入：** 必须提供 `sequence`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sequence_reverse_complement(
    operation='reverse_complement',
    sequence='ATGCGATCG',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "original": "<string>",
    "reverse_complement": "<string>",
    "length": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/04_ops_001_sequence_direction.md#dnakit.databases.sequence_reverse_complement-parameters)

<!-- END GENERATED API: sequence-direction -->

## 5) 计算双链 DNA 的反向互补 {#bio-52-reverse}

- **作用：** 反转双链 DNA 的方向，并保留双链及末端结构。
- **计算方法：** 交换 Watson 链与 Crick 链，并调整两条链的相对位置。
- **API：** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.reverse_complement()`
- **输入：** 完全互补双链或 A/C/G/T 序列；不支持 nick、内部 gap 或错配。
- **输出：** 返回反向互补后的 `Duplex`，保留两条链及末端结构。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
print(molecule.reverse_complement().watson)
```

- **示例结果：**

```text
AACCGGTT
```

示例序列与它的反向互补相同。

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。
