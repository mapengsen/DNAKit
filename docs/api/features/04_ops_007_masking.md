# OPS-007 Masking（掩蔽）

Masking 用指定字符替换 DNA 序列中一个或多个目标区间的序列为任意字符N，保留原有序列长度和坐标。

- **作用：** 在不改变序列长度和坐标的情况下掩蔽低质量、低复杂度或指定区域，避免这些位置影响搜索、统计或建模。
- **API：** `dnakit.ops.mask(dna[必须], intervals[必须], symbol[可选], feature_policy[可选], letter_annotation_policy[可选])`。
- **输入：** `intervals` 是一个或多个 0-based 半开区间 `(start, end)`；`symbol` 默认为 `N`。
- **返回：** 长度不变的新 `DNA`，不会原地修改输入。

<span id="_1"></span>**示例**

```python
from dnakit import DNA
from dnakit.ops import mask

dna = DNA("AACCGGTT", id="seq-1")
masked = mask(dna, [(2, 4), (6, 8)])
print(masked.symbols)
```

```text
AANNGGNN
```
