# OPS-007 Trimming（修剪）

Trimming 用于从 DNA 序列的左端和右端删除指定长度，返回缩短后的新序列。

- **作用：** 从序列两端识别并删除接头、引物或低质量碱基，返回保留片段及修剪坐标，用于 reads 和扩增产物的预处理。
- **API：** `dnakit.ops.trim(dna[必须], left[可选], right[可选], feature_policy[可选], letter_annotation_policy[可选])`。
- **输入：** `left` 是从左端删除的长度，`right` 是从右端删除的长度；两者都是非负整数。
- **返回：** 新的 `DNA`，不会原地修改输入。

<span id="_1"></span>**与任意区间操作的区别**

| 需求 | 应使用的 API |
| --- | --- |
| 从序列左右两端删除碱基 | `trim(dna, left=..., right=...)` |
| 删除序列中任意一段区间 | `delete(dna, start, end)` |
| 只保留序列中任意一段区间 | `subsequence(dna, start, end)` |

`trim()` 专指两端修剪；任意内部位置的删除属于序列编辑，由 `delete()` 处理。

<span id="_2"></span>**示例**

```python
from dnakit import DNA
from dnakit.ops import trim

dna = DNA("AACCGGTT", id="seq-1")
trimmed = trim(dna, left=1, right=2)
print(trimmed.symbols)
```

```text
ACCGG
```

<span id="_3"></span>**限制**

- `left + right` 不能超过序列总长度。
- 可截短已知长度 Gap，但拒绝未知长度 Gap 和环状序列。
- 无后缀的 `trim()` 会同步 feature 和逐碱基注释；`trim_record()` 仅用于读取详细变更审计。
