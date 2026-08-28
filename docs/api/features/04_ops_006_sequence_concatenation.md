# OPS-006 序列拼接

将多个 DNA 序列进行拼接，或根据精确重叠区域进行去重叠拼接。

## 1. 普通序列拼接

- **作用：** 按输入顺序连接多条 DNA 片段，可在连接处插入 linker 或结构化 Gap，并返回各片段在新序列中的坐标，适合构建组合序列。
- **API：** `dnakit.ops.concat(sequences[必须], linker[可选], gap[可选])`。
- **输入：** 至少两条 `DNASequence` 或规范化字符串；`linker` 和 `gap` 不能同时使用。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import concat

joined = concat(
    [DNASequence("AA"), DNASequence("CC")],
    linker="T",
)
print(joined.symbols)
```

- **示例结果：**

```text
AATCC
```

## 2. 去重叠拼接

- **作用：** 根据两段 DNA 的末端重叠关系完成连接，使已匹配的重叠区域只保留一次，并报告接头位置和最终序列。
- **API：** `dnakit.ops.concat_overlap(sequences[必须], min_overlap[可选], max_overlap[可选])`。
- **输入：** 必须是两条线性、无显式 Gap 的 `DNASequence` 或规范化字符串。
- **规则：** 自动选择最长的“左片段后缀 = 右片段前缀”精确重叠；`min_overlap` 默认是 1，`max_overlap` 可限制搜索长度。IUPAC 字符按字面比较，不进行兼容性推断。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import concat_overlap

joined = concat_overlap(
    [DNASequence("AAACCC"), DNASequence("CCCGG")],
)
print(joined.symbols)
```

- **示例结果：**

```text
AAACCCGG
```
