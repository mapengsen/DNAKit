# 环状序列操作

旋转环状 DNA 的起点，或提取跨越环状原点的序列片段。

## 1) 旋转环状序列并提取跨原点片段

- **作用：** 在不改变环状 DNA 内容的情况下重新指定存储起点，或截取跨越原点的片段，并正确换算新旧坐标。
- **计算方法：** 按偏移量循环移动起点；canonical_origin 选择规范起点并同步变换注释坐标。
- **API：** `dnakit.ops.rotate(dna[必须], offset[必须])`、`dnakit.ops.canonical_origin(dna[必须])`，均返回新的 `DNA` 并同步注释。
- **输入：** 必填显式声明为 circular 的单记录 `DNA`；rotate 另需偏移。只提取序列值时可把 `dna.sequence` 传给 `circular_subsequence()`。
- **输出：** 返回旋转、规范化或跨原点提取后的 `DNA` 序列。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import canonical_origin, circular_subsequence, rotate

dna = DNA("GATTACA", topology="circular")
rotated = rotate(dna, 2)
canonical = canonical_origin(dna)
wrapped = circular_subsequence(dna.sequence, 5, 2)
print(rotated.symbols)
print(canonical.symbols)
print(wrapped.symbols)
```

- **示例结果：**

```text
TTACAGA
ACAGATT
CAGA
```

[FAQ：补充说明](../../faq/tutorial_basics.md#compatibility-04_ops_008_circular_sequence_operations)
