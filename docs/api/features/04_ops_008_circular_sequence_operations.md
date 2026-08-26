# OPS-008 环状序列操作

旋转环状 DNA 的起点，或提取跨越环状原点的序列片段。

- **作用：** 在不改变环状 DNA 内容的情况下重新指定存储起点，或截取跨越原点的片段，并正确换算新旧坐标。
- **普通 API：** `dnakit.ops.rotate(dna[必须], offset[必须])`、`dnakit.ops.canonical_origin(dna[必须])`，均返回新的 `DNA` 并同步注释。
- **高级序列值 API：** `dnakit.ops.circular_subsequence(sequence[必须], start[必须], end[必须], allow_gaps[可选])`。
- **输入：** 必填显式声明为 circular 的单记录 `DNA`；rotate 另需偏移。只提取序列值时可把 `dna.sequence` 传给 `circular_subsequence()`。
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

- **限制：** 线性序列会被拒绝；任何未知长度 Gap 都会令总坐标跨度不可解析，即使 origin 位于边界也不能旋转；已知 Gap 则不能把旋转起点放在其内部。`canonical_origin()` 采用正向字典序最小旋转，不考虑反向互补，且拒绝含 Gap 的序列。旧的 `rotate_record()`、`canonical_origin_record()` 只用于详细变更审计。
