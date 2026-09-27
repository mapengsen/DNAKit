# 标准化

统一 DNA 序列的大小写、空白、`U/T` 和非法字符处理方式，从而得到标准化序列。

## 1) 字符标准化 {#1-std-001}

- **作用：** 把原始文本整理成统一的 DNA 序列格式，按设置处理大小写、空白、`U`、模糊字符和非法字符，并记录具体修改，供后续计算使用一致输入。
- **计算方法：** 按配置清理空白、统一大小写，并逐字符保留或删除 U、模糊碱基和其他符号；同时记录修改。
- **API：** `dnakit.normalize(raw[必须], keep_ambiguous[可选], keep_u[可选], keep_other[可选], config[可选])`；`config` 使用 `dnakit.NormalizationConfig`。
- **输入：** 必填原始字符串、UTF-8 `bytes`、`DNASequence`，或由字符串和 `Gap` 组成的可迭代对象。 默认 `keep_ambiguous=True`、`keep_u=False`、`keep_other=False`。
- **输出：** 返回 `NormalizationResult`，示例展示 `sequence`、`changes`。
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

[FAQ：补充说明](../../faq/tutorial_basics.md#normalization-options)
