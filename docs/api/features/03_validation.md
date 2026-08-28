# 合法性检查

检查 DNA 序列、记录或数据集是否符合字母表、长度、Gap、质量和元数据等规则，并返回结构化的问题列表。

只要发现一处错误，返回结果的 `is_valid` 就是 `False`；所有问题仍会保存在结构化 `issues` 中。

- **作用：** 一次检查字母表、非法符号、空序列、长度、Gap、模糊碱基比例、记录元数据、PHRED 质量、空集合和重复 ID，并返回逐项问题报告，用于在分析前筛除或修正无效数据。

- **API：** `dnakit.validate(value[必须], config[可选])`；单条输入返回 `ValidationReport`，集合输入返回 `DatasetValidationReport`。配置类型为 `dnakit.ValidationConfig` 或 `dnakit.DatasetValidationConfig`。
- **输入：** 普通用户直接传 `DNA`。其中只有一条记录时使用 `ValidationConfig` 并返回单条报告；包含多条记录时可使用 `DatasetValidationConfig`，也可直接传 `ValidationConfig` 作为每条记录的规则。旧核心对象仍兼容。
- **检查结果：** 单条结果的 `is_valid=False` 表示至少一条规则产生错误；集合结果的 `is_valid=False` 表示集合为空、任意记录不合法，或启用唯一 ID 检查时发现重复 ID。警告不会单独导致结果不合法。
- **示例代码：**

```python
from dnakit import DNA, ValidationConfig, validate

single = validate(
    DNA("ACGN", alphabet="iupac"),
    config=ValidationConfig(alphabet="strict"),
)
print(single.is_valid)
print([issue.code for issue in single.issues])

records = DNA(
    [
        {"sequence": "ACGT", "id": "valid"},
        {"sequence": "AC", "id": "short"},
    ]
)
collection = validate(records, config=ValidationConfig(min_length=3))
print(collection.is_valid)
print([issue.code for issue in collection.issues])
```

- **示例结果：**

```text
False
['STD_INVALID_SYMBOL']
False
['STD_INVALID_RECORDS']
```

- **命令行：** 原始文本命令仍会先标准化，再调用统一验证逻辑：

```bash
dnakit validate ACGT --sequence-length 4
```
