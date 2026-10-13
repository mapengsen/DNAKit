# G18c 复核：sbol3、screed 与 pyfaidx

本页记录 G18c 组的复核结果。台账中状态为 `other` 的行共 123 个：sbol3 118 个、screed 3 个、pyfaidx 2 个。复核后每一行都已有明确状态和一行原因，本组没有新增包装函数。三个模块 `sbol3_review_functions.py`、`screed_review_functions.py` 与 `pyfaidx_review_functions.py` 都不导出函数，供之后接入时使用。

- 本组没有可离线运行的 DNA 级新函数：sbol3 的 118 行是全局命名空间配置函数、类模块中的重新导出，以及安装目录路径辅助函数；screed 与 pyfaidx 的 5 行是 HAVA 示例格式、文件名与修改时间工具，以及一个与现有包装重复的 pair 写出函数。
- 已有的 DNA 级包装仍在原模块中：sbol3 的 `copy_sbol3_design`、`to_sbol3_display_id` 等在 `molecular_lab/sbol3_functions.py`；screed 与 pyfaidx 在 `sequence/screed_functions.py` 与 `sequence/pyfaidx_functions.py`。
- 测试文件 `others/tests/unit/test_review_G18c.py` 调用真实的上游接口或手工计算，验证每个分类依据。

## 复核结果

| 包 | 复核前 other | 新状态 | 新增 wrapped |
| --- | --- | --- | --- |
| sbol3 | 118 | internal 118 | 0 |
| screed | 3 | not-dna 2，duplicate-of 1（`format_fastx_records`） | 0 |
| pyfaidx | 2 | not-dna 2 | 0 |

- sbol3 的 118 行：每个类模块都重新导出 pySBOL3 的 `get_namespace`、`set_namespace` 与 `set_defaults`，共 38 个模块 × 3 = 114 行；加上顶层同名的 3 行；另有 `sbol3.document.data_path` 1 行。命名空间相关的 `get_namespace`、`set_namespace` 与 `set_defaults` 读写进程内的全局命名空间，`data_path` 只把相对路径拼到 sbol3 安装目录后面。按严格规则，全部标记为 `internal`，不包装。
- screed 的 `hava_iter` 与 `read_hava_sequences`：上游的 HAVA 是 API 参考用的示例格式，每条记录固定 6 行，字段名为占位词（`hava`、`quarzk` 等），没有 DNA 含义，标记为 `not-dna`。
- screed 的 `write_fastx_pair`：上游调用了未定义的 `write_record`，直接抛出 NameError。它要写出的内容与 `format_fastx_records` 相同（逐条写入两条记录），因此标记为 `duplicate-of:format_fastx_records`。
- pyfaidx 的 `getmtime_fsspec` 读取文件修改时间，`get_valid_filename` 清理文件名，两者都不处理序列数据，标记为 `not-dna`。

## 1) `format_fastx_records` 覆盖 `write_fastx_pair`

- **作用：** 把两条记录按顺序写成 FASTA 或 FASTQ 文本。它覆盖上游 `screed.screedRecord.write_fastx_pair` 的用途，而后者在上游不可用。
- **处理方法：** `format_fastx_records` 对每条记录调用 `screed.screedRecord.write_fastx`，输出与两次单条写入的结果相同。FASTQ 模式下每条记录都必须带与序列等长的质量串；上游只检查两条记录是否同时带质量。
- **API：** `format_fastx_records(records[必须], output_format[可选], max_records[可选], verbose[可选])`
- **输入：** 映射列表，每项含 `name`、`sequence`，FASTQ 还需要 `quality`；可选 `description`。
- **输出：** 字典，含 `record_count`、`text` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import format_fastx_records

pair = [
    {"name": "read1", "sequence": "ACGT"},
    {"name": "read2", "sequence": "GGCC"},
]
fasta = format_fastx_records(pair, output_format="fasta", verbose=False)
print(fasta["record_count"])
print(fasta["text"], end="")

fastq = format_fastx_records(
    [dict(pair[0], quality="IIII")],
    output_format="fastq",
    verbose=False,
)
print(fastq["text"], end="")
```

**示例输出：**

```text
2
>read1
ACGT
>read2
GGCC
@read1
ACGT
+
IIII
```

## 2) sbol3 命名空间辅助函数为什么不包装

- **作用：** 说明 `sbol3.component.set_namespace`、`sbol3.attachment.set_defaults` 等 114 个重新导出函数为什么不包装。它们与顶层 `sbol3.set_namespace`、`sbol3.set_defaults` 是同一对象。
- **处理方法：** 这些函数只读写 pySBOL3 进程内的全局命名空间。`set_namespace` 只接受 URL，否则抛出 ValueError；`set_defaults` 把命名空间恢复为 None。包装它们会让一次调用改变其他所有 SBOL3 操作，所以不包装。
- **API：** 无，本组不导出。
- **输入：** 无。
- **输出：** 无。

**示例代码：**

```python
import sbol3

sbol3.component.set_namespace("https://example.org/lab")
print(sbol3.get_namespace())
sbol3.attachment.set_defaults()
print(sbol3.get_namespace())
```

**示例输出：**

```text
https://example.org/lab
None
```

## 未包装的其余行

- sbol3：`document.data_path` 1 行，同样不包装；其余 117 行见第 2 节。
- screed：`hava_iter`、`read_hava_sequences` 两行为 HAVA 示例格式，`not-dna`。
- pyfaidx：`getmtime_fsspec`、`get_valid_filename` 两行为文件工具，`not-dna`。
