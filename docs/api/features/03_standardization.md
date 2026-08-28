# 标准化

统一 DNA 序列的大小写、空白、`U/T` 和非法字符处理方式，从而得到标准化序列。

## 1) STD-001 字符标准化

- **作用：** 把原始文本整理成统一的 DNA 序列格式，按设置处理大小写、空白、`U`、模糊字符和非法字符，并记录具体修改，供后续计算使用一致输入。
- **API：** `dnakit.normalize(raw[必须], keep_ambiguous[可选], keep_u[可选], keep_other[可选], config[可选])`；`config` 使用 `dnakit.NormalizationConfig`。
- **输入：** 必填原始字符串、UTF-8 `bytes`、`DNASequence`，或由字符串和 `Gap` 组成的可迭代对象。
- **核心 args：**
  - `keep_ambiguous=True`：保留 `RYSWKMBDHVN` 模糊碱基；设为 `False` 时删除。
  - `keep_u=False`：默认删除 `U`；设为 `True` 时保留。
  - `keep_other=False`：默认删除除 `A/C/G/T`、IUPAC 模糊碱基和 `U` 之外的其他字符；设为 `True` 时保留。
- **默认行为：** 保留模糊碱基，删除 `U` 和其他非 DNA 字符。
- **示例代码：**

```python
from dnakit import normalize

result = normalize(
    " acn-uX\n",
    keep_ambiguous=True,
    keep_u=False,
    keep_other=False,
)
print(result.sequence.symbols)
print(result.sequence.alphabet.value)
print(
    [
        change.operation
        for change in result.changes
        if change.operation.startswith("delete_")
    ]
)
```

- **示例结果：**

```text
ACN
iupac
['delete_other', 'delete_u', 'delete_other']
```
