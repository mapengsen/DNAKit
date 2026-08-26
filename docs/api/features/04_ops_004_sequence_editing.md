# OPS-004 序列编辑

在指定位置插入、删除或替换 DNA 片段，从而得到编辑后的序列并记录坐标变化。

- **作用：** 在指定位置插入、删除或替换碱基，返回编辑后的序列、操作记录和坐标变化，并同步处理受影响的 feature。
- **API：** `dnakit.ops.insert(dna[必须], position[必须], fragment[必须])`、`dnakit.ops.delete(dna[必须], start[必须], end[必须])`、`dnakit.ops.substitute(dna[必须], start[必须], end[必须], fragment[必须])`。
- **输入：** 普通用户传单记录 `DNA` 和位置/区间；插入或替换还需新 DNA 片段。三个函数都返回新的 `DNA`。
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

- **限制：** 不允许在 Gap 内部插入，也不允许删除或替换跨越 Gap；未知长度 Gap 无法解析编辑坐标。插入/替换片段必须是线性、无 Gap，并与源序列的 strandedness 一致。无后缀函数会同步 feature 和逐碱基注释；原有 `*_record()` 仅在需要详细变更审计时使用。
