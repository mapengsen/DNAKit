# OPS-003 子序列提取

按照线性或环状坐标从 DNA 序列中提取指定片段。

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

- **限制：** 默认拒绝跨 Gap 提取；已知 Gap 只有在 `allow_gaps=True` 时才保留，未知长度 Gap 和非法边界会报错。环状跨原点提取必须显式声明 `topology="circular"`。当前函数不直接接收 `DNAFeature` 或滑窗配置；feature 提取应传其已解析 `Interval.start/end`，滑窗应由调用方显式迭代坐标。
