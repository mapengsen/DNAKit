# OPS-003 子序列提取

从 DNA 序列中提取指定片段。

- **作用：** 按起止坐标截取 DNA 片段并保留来源范围，也支持环状序列跨原点截取，用于提取基因、功能区域或指定窗口。
- **API：** `dnakit.ops.subsequence(sequence[必须], start[必须], end[必须], allow_gaps[可选])`、`dnakit.ops.circular_subsequence(sequence[必须], start[必须], end[必须], allow_gaps[可选])`。
- **输入：** 必填 `DNASequence`、`start`、`end`；可选 `allow_gaps=True` 保留或截短已知 Gap。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import subsequence

seq = DNASequence("AACCGG")
selected = subsequence(seq, 1, 5)
print(selected.symbols)
```

- **示例结果：**

```text
ACCG
```
