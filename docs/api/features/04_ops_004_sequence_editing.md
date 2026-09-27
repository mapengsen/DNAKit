# 序列编辑

在指定位置插入、删除或替换 DNA 片段，切割双链 DNA、修剪序列两端或掩蔽指定区间，并按操作规则处理坐标与注释。

## 1) 插入、删除与替换序列

- **作用：** 在指定位置插入、删除或替换碱基，返回编辑后的序列、操作记录和坐标变化，并同步处理受影响的 feature。
- **处理方法：** 按位置或区间切分序列，插入、移除或替换片段，并更新受影响的坐标。
- **API：** `dnakit.ops.insert(dna[必须], position[必须], fragment[必须])`、`dnakit.ops.delete(dna[必须], start[必须], end[必须])`、`dnakit.ops.substitute(dna[必须], start[必须], end[必须], fragment[必须])`。
- **输入：** 普通用户传单记录 `DNA` 和位置/区间；插入或替换还需新 DNA 片段。三个函数都返回新的 `DNA`。
- **输出：** 三个函数均返回编辑后的新 `DNA`，可从 `symbols` 读取序列，原输入不变。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import delete, insert, substitute

dna = DNA("AACCGG", id="seq-1")
inserted = insert(dna, 2, "TT")
deleted = delete(dna, 2, 4)
replaced = substitute(dna, 2, 4, "TN")
print(inserted.symbols)
print(deleted.symbols)
print(replaced.symbols)
```

- **示例结果：**

```text
AATTCCGG
AAGG
AATNGG
```

<!-- BEGIN GENERATED API: digest -->

<span id="digest"></span>

## 2) 模拟限制性内切酶消化与片段生成 {#dnakit.databases.dna_virtual_digest}

- **作用：** 模拟限制性内切酶消化与片段生成。
- **处理方法：** 匹配限制酶识别位点，按切割偏移拆分线性或环状序列。
- **API：** `dnakit.databases.dna_virtual_digest(*, operation=None, sequence, enzymes=None, circular=None, verbose=True)`
- **输入：** 必须提供 `sequence`。
- **输出：** 返回字典，包含 `fragments`, `n_fragments`, `enzymes_used`, `cut_sites` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dna_virtual_digest(
    operation='virtual_digest',
    sequence='ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT',
    enzymes=['EcoRI', 'BamHI', 'HindIII'],
    circular=False,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "fragments": [
    {
      "sequence": "<string>",
      "length": "<integer>",
      "start": "<integer>",
      "end": "<integer>"
    }
  ],
  "n_fragments": "<integer>",
  "enzymes_used": [
    "<string>"
  ],
  "cut_sites": [
    {
      "enzyme": "<string>",
      "position": "<integer>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/04_ops_004_sequence_editing.md#dnakit.databases.dna_virtual_digest-parameters)

<!-- END GENERATED API: digest -->

## 3) 按两条链的位置切割双链 DNA {#bio-52-cut}

- **作用：** 切开双链 DNA，并检查生成片段的平端或黏端。
- **计算方法：** 显式保存两条 5′→3′ 链和 offset，按两条链的切割位置生成片段。
- **API：** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.cut(watson_position, crick_position)`
- **输入：** 线性互补双链或 A/C/G/T 序列；两条链的切点均使用 Watson 链从左到右的坐标，并位于链内部。不支持 nick、内部 gap 或错配。
- **输出：** 返回切开的两个 `Duplex` 片段及其末端信息。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
left, right = molecule.cut(3, 5)
print(left.ends)
```

- **示例结果：**

```text
{'left': {'polarity': 'blunt', 'sequence': ''}, 'right': {'polarity': '5prime', 'sequence': 'CG'}}
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 4) Trimming（修剪）：修剪序列两端 {#trimming}

Trimming 用于从 DNA 序列的左端和右端删除指定长度，返回缩短后的新序列。

- **作用：** 按指定长度删除序列两端，返回保留片段及同步后的注释。接头或低质量区域的长度需由调用方确定。
- **计算方法：** 按 left 和 right 的长度删除两端，保留中间序列，并按配置转换注释坐标。
- **API：** `dnakit.ops.trim(dna[必须], left[可选], right[可选], feature_policy[可选], letter_annotation_policy[可选])`。
- **输入：** `left` 是从左端删除的长度，`right` 是从右端删除的长度；两者都是非负整数。
- **输出：** 返回修剪后的新 `DNA`，不会原地修改输入。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import trim

dna = DNA("AACCGGTT", id="seq-1")
trimmed = trim(dna, left=1, right=2)
print(trimmed.symbols)
```

- **示例结果：**

```text
ACCGG
```

[FAQ：修剪、删除与截取](../../faq/tutorial_basics.md#trimming-comparison)

## 5) Masking（掩蔽）：掩蔽指定序列区间 {#masking}

Masking 用指定字符（默认为 `N`）替换 DNA 序列中一个或多个目标区间的碱基，保留原有序列长度和坐标。

- **作用：** 在不改变序列长度和坐标的情况下掩蔽低质量、低复杂度或指定区域，避免这些位置影响搜索、统计或建模。
- **计算方法：** 将给定半开区间内的碱基替换为 symbol，保留长度并按配置处理注释。
- **API：** `dnakit.ops.mask(dna[必须], intervals[必须], symbol[可选], feature_policy[可选], letter_annotation_policy[可选])`。
- **输入：** `intervals` 是一个或多个 0-based 半开区间 `(start, end)`；`symbol` 默认为 `N`。
- **输出：** 返回长度不变的掩蔽后新 `DNA`，不会原地修改输入。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import mask

dna = DNA("AACCGGTT", id="seq-1")
masked = mask(dna, [(2, 4), (6, 8)])
print(masked.symbols)
```

- **示例结果：**

```text
AANNGGNN
```
