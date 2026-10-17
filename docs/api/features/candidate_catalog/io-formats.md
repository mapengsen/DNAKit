# 格式转换、压缩与文件管理

本类收录 **131** 项有文档的功能，均可通过正式的 DNAKit 命名空间调用。各格式小节中描述的 96 个函数都会计算真实结果，没有一个无条件返回 `status: "not_implemented"`：88 个 `<操作>_<格式>` 函数（8 种格式 PLINK `map` 与 `ped`、`paf`、`newick`、`nexus`、`phylip`、UCSC `twobit`、`zarr`，每种 11 个操作：`read`、`write`、`validate`、`stream`、`summarize`、`index`、`merge`、`split`、`filter`、`sort`、`convert`），7 个自动识别格式并调用对应函数的通用函数 `convert`、`filter`、`index`、`merge`、`sort`、`stream`、`summarize`，以及 `resume_download`。11 个 `*_zarr` 函数需要可选的 `zarr>=3` 和 `numpy` 包，只有在未安装它们时才返回 `status: "not_implemented"`（带 `required_dependency` 和 `message`）；`resume_download` 需要能访问所下载的服务器。各格式函数除非给出 `output`/`output_dir` 否则不写任何文件，并受 `max_bytes`/`max_records`/`max_cells`/`max_elements` 等上限约束。其余 35 个函数（压缩、校验和、清单、元数据表等）在各自的小节中说明。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) build_file_manifest

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `build_file_manifest` 功能。候选编号 CAP-0308，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `build_file_manifest` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.build_file_manifest(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `build_file_manifest` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import build_file_manifest

result = build_file_manifest(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0308",
  "operation": "build_file_manifest",
  "result_type": "dict",
  "result_preview": "{'path': None, 'size_bytes': 8, 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b', 'media_type': None}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) chunk_large_file

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `chunk_large_file` 功能。候选编号 CAP-0309，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `chunk_large_file` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.chunk_large_file(input[required], chunk_size[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`chunk_size`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `chunk_large_file` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import chunk_large_file

result = chunk_large_file(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0309",
  "operation": "chunk_large_file",
  "result_type": "list",
  "result_preview": "[b'>a\\nACGT\\n']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) clean_temporary_files

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `clean_temporary_files` 功能。候选编号 CAP-0310，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `clean_temporary_files` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.clean_temporary_files(input[required], suffixes[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`suffixes`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `clean_temporary_files` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import clean_temporary_files

result = clean_temporary_files(
    input=Path("nonexistent-temp-directory")
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0310",
  "operation": "clean_temporary_files",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) compress

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `compress` 功能。候选编号 CAP-0311，来自 25 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `compress` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（arrow、bam、bcf、bed、bedgraph、bigbed、bigwig、bim、cram、csv、embl、fam 等，另有 13 种）。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。已有同名函数证据：`compare`。

- **API：** `dnakit.io_formats.compress(input[required], level[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`level`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `compress` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import compress

result = compress(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0311",
  "operation": "compress",
  "result_type": "bytes",
  "result_preview": "<bytes length=28>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) compress_map

- **作用：** 用确定性的 gzip 压缩PLINK `.map` 标记表（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按PLINK `.map` 标记表解析或校验内容，二进制内容也逐字节保留；`decompress_map` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_map

envelope = compress_map('1\trs1\t0\t1000\n1\trs2\t0\t2000\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAAzPkLCo25DTgNDQwMOAyBHKMgBwjEAcAAo9oVBoAAAA=",
    "original_bytes": 26,
    "compressed_bytes": 41,
    "sha256_original": "0d63a7da7a1761ff406f984298d593a88a953100925de4daf0f0b345438301e9"
  }
}
```

## 6) compress_newick

- **作用：** 用确定性的 gzip 压缩Newick 树文件（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按Newick 树文件解析或校验内容，二进制内容也逐字节保留；`decompress_newick` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_newick

envelope = compress_newick('((A:1,B:2):1,C:3);', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA9PQcLQy1HGyMtIEUs5WxprWADsnB7wSAAAA",
    "original_bytes": 18,
    "compressed_bytes": 36,
    "sha256_original": "1687864f64fa679bf85d6ff2ba780918d7c19e3965bcf14f658c9ee9c568c968"
  }
}
```

## 7) compress_nexus

- **作用：** 用确定性的 gzip 压缩NEXUS 文件（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按NEXUS 文件解析或校验内容，二进制内容也逐字节保留；`decompress_nexus` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_nexus

envelope = compress_nexus('#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=2;\n  TAXLABELS A B;\nEND;\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA1P2c40IDeZycnX39FMIcYxwtOZSUHDx9HX1C/b09wtW8AOK2RqBBIEMH0cnV59gBUcFJ2suVz8Xay4AkYWK6z4AAAA=",
    "original_bytes": 62,
    "compressed_bytes": 77,
    "sha256_original": "80a4f6422311b951751d59815bb5a36465fe705e84fb840ecb8f30d8c4a36a8e"
  }
}
```

## 8) compress_paf

- **作用：** 用确定性的 gzip 压缩PAF 比对表（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按PAF 比对表解析或校验内容，二进制内容也逐字节保留；`decompress_paf` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_paf

envelope = compress_paf('q1\t100\t0\t100\t+\tt1\t500\t10\t110\t95\t100\t60\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAAys05DQ0MOA0AJPanCWGnKYGIA6nIRBbmoKFzQy4AKXfwIgnAAAA",
    "original_bytes": 39,
    "compressed_bytes": 48,
    "sha256_original": "b5f356a5606c2b16c7d2b21441ca4b368d794e6aa5cb9dc4cd87d84a373fcbc6"
  }
}
```

## 9) compress_ped

- **作用：** 用确定性的 gzip 压缩PLINK `.ped` 基因型表（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按PLINK `.ped` 基因型表解析或校验内容，二进制内容也逐字节保留；`decompress_ped` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_ped

envelope = compress_ped('F1 I1 0 0 1 -9 A A G G\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA3MzVPA0VDAAQkMFXUsFRyB0V3DnAgC6+u2xFwAAAA==",
    "original_bytes": 23,
    "compressed_bytes": 40,
    "sha256_original": "e51dd9783e6ea7af725321bed570ce6692918994b6b8d982c24cf918e3d6a858"
  }
}
```

## 10) compress_phylip

- **作用：** 用确定性的 gzip 压缩PHYLIP 比对（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按PHYLIP 比对解析或校验内容，二进制内容也逐字节保留；`decompress_phylip` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import compress_phylip

envelope = compress_phylip(' 2 4\nA ACGT\nB ACGA\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA1MwUjDhclRwdHYP4XICUY5cAGNU7OsTAAAA",
    "original_bytes": 19,
    "compressed_bytes": 36,
    "sha256_original": "64de1e5a4b3f8ab13596909fb54c08ddde79434041d40df2eb2b4ed435815d7d"
  }
}
```

## 11) compress_twobit

- **作用：** 用确定性的 gzip 压缩UCSC 2bit 序列文件（文本、字节或文件），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按UCSC 2bit 序列文件解析或校验内容，二进制内容也逐字节保留；`decompress_twobit` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。

- **API：** `dnakit.io_formats.compress_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节）（目录会抛出 `InputFormatError`，只有 zarr 函数接受目录）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import compress_twobit, write_twobit

blob = base64.b64decode(write_twobit({"chr1": "ACGTACGTNN"}, verbose=False)["result"]["data_base64"])
envelope = compress_twobit(blob, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA3NWd5RiAAJGBghgSc4oMpQEMrigghxAzMSAAHPmMAAAOxL/rzQAAAA=",
    "original_bytes": 52,
    "compressed_bytes": 50,
    "sha256_original": "9aa5d8bef92c1ae93359b6716557ed23609d3e7f4719ad399e765e40cb2475d1"
  }
}
```

## 12) compress_zarr

- **作用：** 用确定性的 gzip 压缩Zarr 存储（目录、文件、文本或字节），并以 base64 返回。

- **处理方法：** 先把输入转为字节（见输入），再按 `options['level']`（1-9，默认 6）做 gzip 压缩。gzip 头固定（修改时间为 0，不含文件名），所以相同的输入和级别总是得到相同的 `data`。结果做 base64 编码，并连同大小和未压缩字节的 SHA-256 一起返回。函数不会按Zarr 存储解析或校验内容，二进制内容也逐字节保留；`decompress_zarr` 可还原完全相同的字节。只有给出 `options['output']` 才会写文件：此时把原始 gzip 字节写入该文件（已存在的文件需要 `overwrite=True`）。目录会先被归档：其中的文件和子目录（路径相对于该目录）按排序顺序写入 tar 归档，时间和属主清零，因此归档可复现，被压缩的就是这个 tar 归档。目录内的符号链接和特殊文件会抛出 `InputFormatError`。

- **API：** `dnakit.io_formats.compress_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本（`str`，按 UTF-8 编码）、`bytes`/`bytearray`，或已存在文件的路径（`str` 或 `Path`，压缩其字节），或目录的路径（见处理方法）。单行字符串若是已存在的路径则按路径读取。不存在的 `Path`，或带有已知文件后缀（如 `.2bit`）但不存在的字符串，会抛出 `InputFormatError`（代码 `IO_FORMATS_FILE_NOT_FOUND`），而不是被当作文本压缩。不接受映射、列表和数字。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝大于该值的输入
    - `level`：1-9 的整数，默认 6：gzip 压缩级别
    - `output`：要写入的 .gz 文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`encoding`（`"gzip+base64"`）、`data`（base64 文本）、`original_bytes`（压缩前大小）、`compressed_bytes`（gzip 流大小）、`sha256_original`（未压缩字节的 SHA-256）、`archive`（`"tar"`，仅当输入是目录时），给出 `options['output']` 时还有 `output_path`。无法使用的输入抛出 `InputFormatError`；选项无效或含未知选项（如 `level` 不在 1-9 内）、输出文件已存在而未设 `overwrite=True`、输入超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json
import shutil
from pathlib import Path

from dnakit.io_formats import compress_zarr

store = Path("demo_store.zarr")
(store / "a").mkdir(parents=True, exist_ok=True)
(store / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')
(store / "a" / "zarr.json").write_text('{"zarr_format": 3, "node_type": "array"}')

envelope = compress_zarr(store, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
shutil.rmtree(store)
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA+3T0QqCMBiG4V3K2HGklnPQzcggC4KcTD2w6N6bHYUUQTATep+Tf2MH++Hjs4mILg2M1o8ZTOeLc5FttJA6/mpC9G1nvZRzfLVENrlY79en1tXR/hhDLfL8ff56mr8xJhMyjbbRkz/P/6rG/MuD82fbqZ3crqSq3b4qu6Gpwl2FVzuo26/3RBzx2/9V/7Wm/7P42P+jd31D/wEAAAAAAAAAAAAAWLI7rWWEdQAoAAA=",
    "original_bytes": 10240,
    "compressed_bytes": 185,
    "sha256_original": "ae3c6d9a98200511fcb14ad7fa368b2b6fee1b53c073c3b76128f516cb2a8421",
    "archive": "tar"
  }
}
```

## 13) convert

- **作用：** 把任意受支持格式（map、ped、paf、newick、nexus、phylip、twobit、zarr；FASTA 仅作来源）的文件转换为另一种格式。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `convert_ped`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 目标格式为 `options['to']`（必填），必须是所识别格式支持的目标之一。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.convert(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填），必须是所识别格式支持的目标之一
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = convert(PED, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 3,
    "source_format": "ped",
    "bytes": 101,
    "sha256": "288b2c84cc53c5a994f59ab2bcc18e128999dacfd1d9ffe7bf1629c930e9ea3f",
    "text": "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tV1\tV2\nF1\tI1\t0\t0\t1\t2\tA/A\tC/G\nF1\tI2\t0\t0\t2\t1\tA/G\t0/0\nF1\tI3\tI1\tI2\t1\t-9\tG/G\tC/C\n"
  }
}
```

## 14) convert_map

- **作用：** 把 PLINK `.map` 文件转换为其他格式（`map`、`bed`、`tsv` 或 `json`）。

- **处理方法：** `map` 写出规范化的制表符分隔文件；`bed` 为每个变异写出 `chrom`、`bp-1`、`bp`、`ID`（BED 为 0 起始、左闭右开，所以 1 起始的 bp 变成区间 [bp-1, bp)）；`tsv` 增加表头 `chrom id cm pos`；`json` 写出记录列表。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：map、bed、tsv、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = convert_map(MAP, options={"to": "bed"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "bed",
    "record_count": 4,
    "source_format": "map",
    "bytes": 46,
    "sha256": "c3b8e07bd896ef7348d8fd8270d0bf63a95fb78bc265fc8cd8488a90e1240e9f",
    "text": "1\t99\t100\trs1\n1\t49\t50\trs2\n2\t9\t10\trs3\nX\t4\t5\trs4\n"
  }
}
```

## 15) convert_newick

- **作用：** 把 Newick 文件转换为其他格式（`newick`、`nexus`、`tsv` 或 `json`）。

- **处理方法：** `newick` 写出规范化文件；`nexus` 写出含不同叶节点名的 TAXA 块和树名为 `tree1`、`tree2`……的 TREES 块；`tsv` 写出节点表 `tree, node_id, parent_id, name, length, is_tip`（树序号从 0 起，根节点的 parent 为空）；`json` 写出带扁平节点表的记录列表。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：newick、nexus、tsv、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = convert_newick(TREES, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 2,
    "source_format": "newick",
    "bytes": 188,
    "sha256": "45d02b5266f4e4b290b362ea416595353d61720d12bf5f014eba3e7b64903f3f",
    "text": "tree\tnode_id\tparent_id\tname\tlength\tis_tip\n0\t0\t\t\t\t0\n0\t1\t0\t\t0.5\t0\n0\t2\t1\tA\t1\t1\n0\t3\t1\tB\t2\t1\n0\t4\t0\t\t0.5\t0\n0\t5\t4\tC\t1\t1\n0\t6\t4\tD\t1\t1\n1\t0\t\t\t\t0\n1\t1\t0\tA\t\t1\n1\t2\t0\tB\t\t1\n1\t3\t0\t\t\t0\n1\t4\t3\tC\t\t1\n1\t5\t3\tD\t\t1\n"
  }
}
```

## 16) convert_nexus

- **作用：** 把 NEXUS 文件转换为其他格式（`nexus`、`phylip`、`fasta`、`twobit`、`newick` 或 `json`）。

- **处理方法：** `nexus` 写出规范化的 DATA 块和 TREES 块；`phylip` 和 `fasta` 写出矩阵（含空白的分类单元名无法写入 PHYLIP，会抛出 `InputFormatError`）；`twobit` 把核苷酸矩阵打包成 2bit 文件（非 ACGT 字母变成 N 块）；`newick` 每行写一棵树（不含树的文件抛出 `InputFormatError`）；`json` 写出含 `taxon_count`、`site_count`、`datatype`、`trees` 和 `records` 的对象。`options['to']` 必填。`phylip`、`fasta` 和 `twobit` 目标只保存矩阵：文件同时含 TREES 块时，照常转换矩阵，并返回一条说明树已被丢弃的警告。

- **API：** `dnakit.io_formats.convert_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：nexus、phylip、fasta、twobit、newick、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `layout`：'sequential'（默认）或 'interleaved'：PHYLIP/NEXUS 输出的排布
    - `block_size`：整数，默认 60：交错（interleaved）输出每块的位点数
    - `strict_names`：布尔，默认 False：写出严格的 10 字符 PHYLIP 名称（更长时报错）
    - `line_width`：整数，默认 60：FASTA 行宽（0 表示每条序列一行）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`（`twobit` 目标为 `data_base64`，并另有 `sequence_count`）；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。丢弃了树时（目标 `phylip`、`fasta`、`twobit`），`result` 还有 `warnings`（含一条消息的列表），同一条消息也在封装的 `warnings` 中。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""


envelope = convert_nexus(NEXUS, options={"to": "twobit"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [
    "The 1 tree(s) of the TREES block were dropped: the twobit format holds only the matrix."
  ],
  "result": {
    "format": "twobit",
    "record_count": 3,
    "sequence_count": 3,
    "warnings": [
      "The 1 tree(s) of the TREES block were dropped: the twobit format holds only the matrix."
    ],
    "source_format": "nexus",
    "bytes": 98,
    "sha256": "90833cd90314c8282040eb0e07e8cdb33abbfa3834000cbcbd4ba6e1ad8781e5",
    "data_base64": "QydBGgAAAAADAAAAAAAAAAFBJAAAAAFCNgAAAANDIGRIAAAACAAAAAAAAAAAAAAAAAAAAJycCAAAAAAAAAAAAAAAAAAAAJyeCAAAAAEAAAADAAAAAQAAAAAAAAAAAAAAnJw="
  }
}
```

## 17) convert_paf

- **作用：** 把 PAF 文件转换为其他格式（`paf`、`bed`、`tsv` 或 `json`）。

- **处理方法：** `paf` 写出规范化文件；`bed` 为每条比对写出目标区间（`target`、`target_start`、`target_end`、`query`、`mapq`、`strand`；无表头）；`tsv` 写出含 12 列以及每个标签名一列（按名称排序）的表头和标签取值（缺失时为空）；`json` 写出记录列表。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：paf、bed、tsv、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = convert_paf(PAF, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 2,
    "source_format": "paf",
    "bytes": 216,
    "sha256": "9fdb7fe2e77ddbbc6c4e27f7eec4b229553c7255da17988ad1e489b3afbf6752",
    "text": "query_name\tquery_length\tquery_start\tquery_end\tstrand\ttarget_name\ttarget_length\ttarget_start\ttarget_end\tmatches\tblock_length\tmapq\tcg\ttp\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\t90M\tP\nq2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\t\t\n"
  }
}
```

## 18) convert_ped

- **作用：** 把 PLINK `.ped` 文件转换为其他格式（`ped`、`tsv` 或 `json`）。

- **处理方法：** `ped` 写出规范化的制表符分隔文件；`tsv` 写出表头 `FID IID PAT MAT SEX PHENO V1..Vk`，基因型写成 `a1/a2`；`json` 写出含 `variant_count` 和 `records` 的对象。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：ped、tsv、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = convert_ped(PED, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 3,
    "source_format": "ped",
    "bytes": 101,
    "sha256": "288b2c84cc53c5a994f59ab2bcc18e128999dacfd1d9ffe7bf1629c930e9ea3f",
    "text": "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tV1\tV2\nF1\tI1\t0\t0\t1\t2\tA/A\tC/G\nF1\tI2\t0\t0\t2\t1\tA/G\t0/0\nF1\tI3\tI1\tI2\t1\t-9\tG/G\tC/C\n"
  }
}
```

## 19) convert_phylip

- **作用：** 把 PHYLIP 文件转换为其他格式（`phylip`、`nexus`、`fasta`、`twobit` 或 `json`）。

- **处理方法：** `phylip` 重新排布比对（`layout`、`block_size`、`strict_names`）；`nexus` 写出 DATA 块；`fasta` 写出按 `line_width` 折行的记录（0 表示每条序列一行）；`twobit` 把核苷酸序列打包成 2bit 文件（非 ACGT 字母变成 N 块）；`json` 写出记录。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：phylip、nexus、fasta、twobit、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）
    - `layout`：'sequential'（默认）或 'interleaved'：PHYLIP/NEXUS 输出的排布
    - `block_size`：整数，默认 60：交错（interleaved）输出每块的位点数
    - `strict_names`：布尔，默认 False：写出严格的 10 字符 PHYLIP 名称（更长时报错）
    - `line_width`：整数，默认 60：FASTA 行宽（0 表示每条序列一行）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`（`twobit` 目标为 `data_base64`，并另有 `sequence_count`）；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = convert_phylip(PHYLIP, options={"to": "fasta"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "fasta",
    "record_count": 3,
    "source_format": "phylip",
    "bytes": 36,
    "sha256": "5f65c0036a72aabb466c7dc332dea148b25953356b317de61c897fe465de5dea",
    "text": ">A\nACGTACGT\n>B\nACGTACGA\n>C\nACG-ACGT\n"
  }
}
```

## 20) convert_twobit

- **作用：** 把 2bit 文件转换为其他格式（`twobit`、`fasta` 或 `json`）。

- **处理方法：** `fasta` 按 `line_width` 折行写出序列（0 表示每条序列一行），除非 `soft_mask=False`，否则软屏蔽以小写保留；`json` 写出含 `name`、`length`、`n_blocks`、`mask_blocks` 和 `sequence` 的列表；`twobit` 重新编码文件（例如改变 `byte_order`）。`names` 把转换限制在部分序列；解码超过 `max_bases` 个碱基会被拒绝。FASTA 转 2bit 用 `convert(format='fasta', to='twobit')`。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `to`：目标格式名（必填）：twobit、fasta、json
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `line_width`：整数，默认 60：FASTA 行宽（0 表示一行）
    - `soft_mask`：布尔，默认 True：用小写表示软屏蔽碱基
    - `max_bases`：整数，默认 50000000：最多解码的碱基数
    - `names`：要转换的序列名列表（默认全部）
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`（目标格式）、`record_count`、`source_format`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`（`twobit` 目标为 `data_base64`，并另有 `sequence_count`）；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import convert_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = convert_twobit(TWOBIT, options={"to": "fasta", "line_width": 5}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "fasta",
    "record_count": 2,
    "source_format": "twobit",
    "bytes": 30,
    "sha256": "d7f761cab7931fa25acad2f48d3b4914cca20d16f8ba89ba150358615295ea28",
    "text": ">s1\nACGTN\nNacgt\n>s2\nGGCCA\nATT\n"
  }
}
```

## 21) convert_zarr

- **作用：** 把 Zarr 存储转换为其他格式（`zarr`、`json`、`tsv` 或 `npz`）。

- **处理方法：** `json` 写出元数据和数值；`tsv` 以制表符分隔、无表头地写出一个 1 维或 2 维数组（用 `array` 选择，存储含多个数组时必填）；`npz` 写出包含所有数组的 numpy 归档；`zarr` 把存储复制到 `output`（此目标必填）。转换超过 `max_elements` 个元素会被拒绝。`options['to']` 必填。

- **API：** `dnakit.io_formats.convert_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `to`：目标格式名（必填）：zarr、json、tsv、npz
    - `output`：输出路径（'zarr' 目标必填，其他目标可选）
    - `overwrite`：布尔，默认 False：替换 output 处已存在的 zarr 存储
    - `array`：'tsv' 目标使用的数组路径（1 维或 2 维）
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`json`、`tsv`、`npz` 目标为通常的 `format`（目标格式）、`record_count`、`bytes`、`sha256` 和 `text`（`npz` 为 `data_base64`），给出 `options['output']` 时改为返回 `output_path`（已存在的文件需要 `overwrite=True`）；`zarr` 目标为 `format`、`record_count` 和 `output_path`。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import convert_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = convert_zarr("demo.zarr", options={"to": "tsv", "array": "dosage"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 4,
    "bytes": 16,
    "sha256": "b593e6e0fd7a0d5575370d49623eabd2848d22f9592cd6b8f9e90695c6ec6d57",
    "text": "0\t1\n1\t1\n0\t0\n2\t1\n"
  }
}
```

## 22) decode_metadata_categorical

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `decode_metadata_categorical` 功能。候选编号 CAP-0329，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `decode_metadata_categorical` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.decode_metadata_categorical(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `decode_metadata_categorical` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import decode_metadata_categorical

result = decode_metadata_categorical(
    input={"data": {"group": [0, 1]}, "categories": {"group": ["A", "B"]}}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0329",
  "operation": "decode_metadata_categorical",
  "result_type": "dict",
  "result_preview": "{'group': ['A', 'B']}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) decompress

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `decompress` 功能。候选编号 CAP-0330，来自 25 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `decompress` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（arrow、bam、bcf、bed、bedgraph、bigbed、bigwig、bim、cram、csv、embl、fam 等，另有 13 种）。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.decompress(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `decompress` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import decompress

result = decompress(
    input=__import__("gzip").compress(b"ACGT")
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0330",
  "operation": "decompress",
  "result_type": "bytes",
  "result_preview": "<bytes length=4>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) decompress_map

- **作用：** 把 gzip 数据（`compress_map` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的PLINK `.map` 标记表的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做PLINK `.map` 标记表的解析或校验。

- **API：** `dnakit.io_formats.decompress_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_map` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_map

payload = 'H4sIAAAAAAAAAzPkLCo25DTgNDQwMOAyBHKMgBwjEAcAAo9oVBoAAAA='
envelope = decompress_map(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "original_bytes": 26,
    "compressed_bytes": 41,
    "sha256": "0d63a7da7a1761ff406f984298d593a88a953100925de4daf0f0b345438301e9",
    "text": "1\trs1\t0\t1000\n1\trs2\t0\t2000\n"
  }
}
```

## 25) decompress_newick

- **作用：** 把 gzip 数据（`compress_newick` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的Newick 树文件的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做Newick 树文件的解析或校验。

- **API：** `dnakit.io_formats.decompress_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_newick` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_newick

payload = 'H4sIAAAAAAAAA9PQcLQy1HGyMtIEUs5WxprWADsnB7wSAAAA'
envelope = decompress_newick(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "original_bytes": 18,
    "compressed_bytes": 36,
    "sha256": "1687864f64fa679bf85d6ff2ba780918d7c19e3965bcf14f658c9ee9c568c968",
    "text": "((A:1,B:2):1,C:3);"
  }
}
```

## 26) decompress_nexus

- **作用：** 把 gzip 数据（`compress_nexus` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的NEXUS 文件的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做NEXUS 文件的解析或校验。

- **API：** `dnakit.io_formats.decompress_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_nexus` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_nexus

payload = 'H4sIAAAAAAAAA1P2c40IDeZycnX39FMIcYxwtOZSUHDx9HX1C/b09wtW8AOK2RqBBIEMH0cnV59gBUcFJ2suVz8Xay4AkYWK6z4AAAA='
envelope = decompress_nexus(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "original_bytes": 62,
    "compressed_bytes": 77,
    "sha256": "80a4f6422311b951751d59815bb5a36465fe705e84fb840ecb8f30d8c4a36a8e",
    "text": "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=2;\n  TAXLABELS A B;\nEND;\n"
  }
}
```

## 27) decompress_paf

- **作用：** 把 gzip 数据（`compress_paf` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的PAF 比对表的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做PAF 比对表的解析或校验。

- **API：** `dnakit.io_formats.decompress_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_paf` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_paf

payload = 'H4sIAAAAAAAAAys05DQ0MOA0AJPanCWGnKYGIA6nIRBbmoKFzQy4AKXfwIgnAAAA'
envelope = decompress_paf(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "original_bytes": 39,
    "compressed_bytes": 48,
    "sha256": "b5f356a5606c2b16c7d2b21441ca4b368d794e6aa5cb9dc4cd87d84a373fcbc6",
    "text": "q1\t100\t0\t100\t+\tt1\t500\t10\t110\t95\t100\t60\n"
  }
}
```

## 28) decompress_ped

- **作用：** 把 gzip 数据（`compress_ped` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的PLINK `.ped` 基因型表的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做PLINK `.ped` 基因型表的解析或校验。

- **API：** `dnakit.io_formats.decompress_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_ped` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_ped

payload = 'H4sIAAAAAAAAA3MzVPA0VDAAQkMFXUsFRyB0V3DnAgC6+u2xFwAAAA=='
envelope = decompress_ped(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "original_bytes": 23,
    "compressed_bytes": 40,
    "sha256": "e51dd9783e6ea7af725321bed570ce6692918994b6b8d982c24cf918e3d6a858",
    "text": "F1 I1 0 0 1 -9 A A G G\n"
  }
}
```

## 29) decompress_phylip

- **作用：** 把 gzip 数据（`compress_phylip` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的PHYLIP 比对的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做PHYLIP 比对的解析或校验。

- **API：** `dnakit.io_formats.decompress_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_phylip` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_phylip

payload = 'H4sIAAAAAAAAA1MwUjDhclRwdHYP4XICUY5cAGNU7OsTAAAA'
envelope = decompress_phylip(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "original_bytes": 19,
    "compressed_bytes": 36,
    "sha256": "64de1e5a4b3f8ab13596909fb54c08ddde79434041d40df2eb2b4ed435815d7d",
    "text": " 2 4\nA ACGT\nB ACGA\n"
  }
}
```

## 30) decompress_twobit

- **作用：** 把 gzip 数据（`compress_twobit` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的UCSC 2bit 序列文件的完全相同的字节。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做UCSC 2bit 序列文件的解析或校验。

- **API：** `dnakit.io_formats.decompress_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_twobit` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：写入解压字节的文件路径（可选；不给则内联返回）
    - `overwrite`：布尔值，默认 False：允许替换已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import decompress_twobit

payload = 'H4sIAAAAAAAAA3NWd5RiAAJGBghgSc4oMpQEMrigghxAzMSAAHPmMAAAOxL/rzQAAAA='
envelope = decompress_twobit(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "original_bytes": 52,
    "compressed_bytes": 50,
    "sha256": "9aa5d8bef92c1ae93359b6716557ed23609d3e7f4719ad399e765e40cb2475d1",
    "data_base64": "QydBGgAAAAABAAAAAAAAAARjaHIxGQAAAAoAAAABAAAACAAAAAIAAAAAAAAAAAAAAJycAA=="
  }
}
```

## 31) decompress_zarr

- **作用：** 把 gzip 数据（`compress_zarr` 的 `data`、gzip 字节或 `.gz` 文件）解码还原为原来的Zarr 存储的完全相同的字节（目录归档会解压到 `options['output']`）。

- **处理方法：** 先把输入转为 gzip 流：以 gzip 魔数（`1f 8b`）开头的字节原样使用，其他输入必须是 gzip 流的合法 base64 文本（忽略空白）。解压受 `max_bytes` 限制（接受串联的 gzip 成员；数据被截断或损坏、末尾有多余数据都会报错）。字节是合法 UTF-8 时以 `text` 返回，否则以 `data_base64` 返回，因此往返逐字节一致；给出 `options['output']` 时改为把字节写入该文件（已存在的文件需要 `overwrite=True`）。不做Zarr 存储的解析或校验。若解压出的字节是 tar 归档（`compress_zarr` 对目录生成的就是这种归档），`options['output']` 指定解压到的目录（除非 `overwrite=True`，否则必须是新目录或空目录），并安全地解包：先检查所有成员再写入任何内容，绝对路径、`..` 分量、链接、设备和其他特殊成员都会被拒绝并抛出 `InputFormatError`（代码 `IO_FORMATS_UNSAFE_ARCHIVE`）。不给 `output` 时不写任何文件，tar 归档以 `data_base64` 返回。

- **API：** `dnakit.io_formats.decompress_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：base64 文本，即 `compress_zarr` 的 `result["data"]`；gzip `bytes`；或 `.gz` 文件的路径（`str` 或 `Path`，原始 gzip，或保存 base64 文本的文件）。不接受整个结果封装或其 `result` 映射。不存在的 `Path`，或带有已知后缀（如 `.gz`）但不存在的字符串，会抛出 `InputFormatError`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝输入或解压后数据大于该值
    - `output`：解包归档的目录（除非 overwrite=True，否则必须是新目录或空目录），对其他数据则是文件路径（可选；不给则不写任何文件）
    - `overwrite`：布尔值，默认 False：允许非空目录（已有文件会被替换）或已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`original_bytes`（解压后大小）、`compressed_bytes`、`sha256`（解压后字节的 SHA-256），以及内容 `text`（合法 UTF-8）或 `data_base64`，给出 `options['output']` 时为 `output_path`。对于 tar 归档，`result` 含 `archive: "tar"`，以及 `output_path` 加 `file_count`（解出的文件数）或 `data_base64`（归档本身）。输入不是合法 base64、不是 gzip、被截断或已损坏时抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_BASE64`、`IO_FORMATS_BAD_GZIP`）并附提示；不再返回带失败结果的 `status: "ok"`。选项无效或含未知选项、输出已存在而未设 `overwrite=True`、数据超过 `max_bytes` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json
import shutil

from dnakit.io_formats import decompress_zarr

payload = 'H4sIAAAAAAAAA+3T0QqCMBiG4V3K2HGklnPQzcggC4KcTD2w6N6bHYUUQTATep+Tf2MH++Hjs4mILg2M1o8ZTOeLc5FttJA6/mpC9G1nvZRzfLVENrlY79en1tXR/hhDLfL8ff56mr8xJhMyjbbRkz/P/6rG/MuD82fbqZ3crqSq3b4qu6Gpwl2FVzuo26/3RBzx2/9V/7Wm/7P42P+jd31D/wEAAAAAAAAAAAAAWLI7rWWEdQAoAAA='
envelope = decompress_zarr(payload, options={"output": "restored_store"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
shutil.rmtree("restored_store")
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "original_bytes": 10240,
    "compressed_bytes": 185,
    "sha256": "ae3c6d9a98200511fcb14ad7fa368b2b6fee1b53c073c3b76128f516cb2a8421",
    "archive": "tar",
    "file_count": 2,
    "output_path": "restored_store"
  }
}
```

## 32) detect_compression

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `detect_compression` 功能。候选编号 CAP-0339，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `detect_compression` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.detect_compression(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `detect_compression` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import detect_compression

result = detect_compression(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0339",
  "operation": "detect_compression",
  "result_type": "str",
  "result_preview": "'none'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) detect_file_format

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `detect_file_format` 功能。候选编号 CAP-0340，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `detect_file_format` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.detect_file_format(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `detect_file_format` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import detect_file_format

result = detect_file_format(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0340",
  "operation": "detect_file_format",
  "result_type": "str",
  "result_preview": "'fasta'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) encode_metadata_categorical

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `encode_metadata_categorical` 功能。候选编号 CAP-0342，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `encode_metadata_categorical` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.encode_metadata_categorical(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `encode_metadata_categorical` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import encode_metadata_categorical

result = encode_metadata_categorical(
    input={"group": ["A", "B", "A"]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0342",
  "operation": "encode_metadata_categorical",
  "result_type": "dict",
  "result_preview": "{'data': {'group': [0, 1, 0]}, 'categories': {'group': ['A', 'B']}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) estimate_storage_size

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `estimate_storage_size` 功能。候选编号 CAP-0343，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `estimate_storage_size` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.estimate_storage_size(input[required], compression[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`compression`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `estimate_storage_size` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import estimate_storage_size

result = estimate_storage_size(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0343",
  "operation": "estimate_storage_size",
  "result_type": "int",
  "result_preview": "8",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) filter

- **作用：** 用对应格式专用过滤函数的条件过滤任意受支持格式的文件。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `filter_paf`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.filter(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = filter(PAF, options={"min_mapq": 30}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 55,
    "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
    "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 37) filter_map

- **作用：** 按染色体、变异 ID、位置和遗传距离过滤 PLINK `.map` 文件，并返回过滤后的 `.map` 内容。

- **处理方法：** 变异必须通过所有给定条件（逻辑与）才被保留：染色体在/不在给定列表中、ID 在/不在给定列表中、`region`（`chrom:start-end`，1 起始、含两端，染色体相同且 `start <= bp <= end`）、bp 的 `min_position`/`max_position`，以及遗传距离的 `min_cm`/`max_cm`。保持文件顺序；不给任何条件时全部通过。

- **API：** `dnakit.io_formats.filter_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `chromosomes`：要保留的染色体名列表
    - `exclude_chromosomes`：要丢弃的染色体名列表
    - `ids`：要保留的变异 ID 列表
    - `exclude_ids`：要丢弃的变异 ID 列表
    - `region`：'chrom:start-end'：保留 start <= bp <= end 的变异（1 起始，含两端）
    - `min_position`：整数：保留 bp >= 该值的变异
    - `max_position`：整数：保留 bp <= 该值的变异
    - `min_cm`：浮点数：保留遗传距离 >= 该值的变异
    - `max_cm`：浮点数：保留遗传距离 <= 该值的变异

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = filter_map(MAP, options={"chromosomes": ["1", "X"], "max_position": 90}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 2,
    "input_count": 4,
    "output_count": 2,
    "removed_count": 2,
    "bytes": 23,
    "sha256": "7175112a87148a8e643383f0e65de0b0832f3778411172821c9bf7c5e8da8881",
    "text": "1\trs2\t0.5\t50\nX\trs4\t0\t5\n"
  }
}
```

## 38) filter_newick

- **作用：** 按叶节点名、短枝和叶节点数过滤 Newick 文件，并返回过滤后的 Newick 内容。

- **处理方法：** 对每棵树依次执行：`keep_tips`/`drop_tips` 把树剪枝到需要的叶节点（没有剩余叶节点的树被丢弃；只剩一个子节点的内部节点被合并，枝长相加）；`collapse_below` 把短于阈值的内部枝收缩为多分支节点；然后 `min_tips`/`max_tips` 丢弃剩余叶节点数不在范围内的树。

- **API：** `dnakit.io_formats.filter_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `keep_tips`：要保留的叶节点名列表（其余叶节点被剪除；单子节点被合并，枝长相加）
    - `drop_tips`：要剪除的叶节点名列表
    - `collapse_below`：浮点数：把短于该值的内部枝收缩为多分支节点
    - `min_tips`：整数：丢弃叶节点数更少的树（剪枝之后判断）
    - `max_tips`：整数：丢弃叶节点数更多的树（剪枝之后判断）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = filter_newick(TREES, options={"drop_tips": ["D"]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "removed_count": 0,
    "bytes": 32,
    "sha256": "f2260e09bd860329f784c631272292e6da494b022ae3ee88587218dff8d2e942",
    "text": "((A:1,B:2):0.5,C:1.5);\n(A,B,C);\n"
  }
}
```

## 39) filter_nexus

- **作用：** 按分类单元名以及每个分类单元/每个位点的缺失数据过滤 NEXUS 文件，并把树剪枝到保留的分类单元。

- **处理方法：** 分类单元必须通过所有给定条件才保留：在 `keep` 中、不在 `exclude` 中、匹配 `name_pattern`（`re.search`，最长 256 个字符），且缺口/缺失/歧义比例不超过 `max_missing_fraction`。随后删去在保留的分类单元上缺失比例超过 `max_site_missing_fraction` 的比对列。分类单元过滤作用于矩阵、TAXA 标签以及（通过剪去叶节点）每棵树（没有剩余叶节点的树被丢弃）。

- **API：** `dnakit.io_formats.filter_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `keep`：要保留的分类单元名列表
    - `exclude`：要丢弃的分类单元名列表
    - `name_pattern`：分类单元名必须匹配的正则表达式（re.search，最长 256 个字符）
    - `max_missing_fraction`：[0,1] 内的浮点数：丢弃缺口/缺失/歧义比例超过它的分类单元
    - `max_site_missing_fraction`：[0,1] 内的浮点数：丢弃缺失比例（在保留的分类单元上计算）超过它的比对列

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = filter_nexus(NEXUS, options={"exclude": ["B"]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 3,
    "input_count": 4,
    "output_count": 3,
    "removed_count": 1,
    "bytes": 207,
    "sha256": "a4c11f842c4c2c54c81268d8380d449f4dfe72e33af3929d45d73aa6c735684e",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=2 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGTACGT\n    'C d'  ACG-ACGT\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] (A:2,'C d':2);\nEND;\n"
  }
}
```

## 40) filter_paf

- **作用：** 按比对质量、一致性、长度、链方向、名称、区间和标签条件过滤 PAF 文件，并返回过滤后的 PAF 内容。

- **处理方法：** 比对必须通过所有给定条件（逻辑与）才被保留：`min_mapq`、`min_matches`、`min_block_length`、`min_identity`（matches / block_length）、`min_query_coverage`（(query_end - query_start) / query_length）、`strand`、要保留或丢弃的查询/目标序列名列表、`target_region`/`query_region`（`name:start-end`，比对与该序列上 0 起始、左闭右开的区间重叠即保留）以及 `tags`（所列标签都必须存在且取值相同，按文本比较）。保持输入顺序。

- **API：** `dnakit.io_formats.filter_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `min_mapq`：整数：保留比对质量 >= 该值的记录
    - `min_identity`：浮点数：保留 matches/block_length >= 该值的记录
    - `min_block_length`：整数：保留比对块长度 >= 该值的记录
    - `min_matches`：整数：保留匹配碱基数 >= 该值的记录
    - `min_query_coverage`：浮点数：保留 (query_end-query_start)/query_length >= 该值的记录
    - `strand`：'+' 或 '-'
    - `query_names`：要保留的查询序列名列表
    - `exclude_query_names`：要丢弃的查询序列名列表
    - `target_names`：要保留的目标序列名列表
    - `exclude_target_names`：要丢弃的目标序列名列表
    - `target_region`：'name:start-end'：保留与该 0 起始、左闭右开的目标区间重叠的比对
    - `query_region`：'name:start-end'：保留与该 0 起始、左闭右开的查询区间重叠的比对
    - `tags`：映射 {tag: value}：保留带有所有这些标签且取值相同的比对（如 {'tp': 'P'}）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = filter_paf(PAF, options={"min_mapq": 30, "tags": {"tp": "P"}}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 55,
    "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
    "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 41) filter_ped

- **作用：** 按家系、个体、性别、表型和基因型缺失条件过滤 PLINK `.ped` 文件，并返回过滤后的 `.ped` 内容。

- **处理方法：** 个体过滤条件以逻辑与组合：保留或丢弃的家系/个体、`sex`（1/`male`、2/`female`、0/`unknown`）、`phenotype` 类别（`case` = 2、`control` = 1、`missing` = 0 或 -9、`observed` = 非缺失、`quantitative` = 其他任何值）、`founders_only`（父母均为 `0`）以及 `max_missing_rate`（每个个体缺失基因型的比例）。随后变异过滤从所有保留的行中删去基因型列：`variant_indices` 保留列出的 0 起始列，`max_variant_missing_rate` 丢弃缺失比例（在通过其他过滤条件的个体上计算）超过阈值的列。对应的 `.map` 文件需要由调用者做同样的列筛选。

- **API：** `dnakit.io_formats.filter_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `families`：要保留的家系 ID 列表
    - `exclude_families`：要丢弃的家系 ID 列表
    - `individuals`：要保留的个体 ID 列表
    - `exclude_individuals`：要丢弃的个体 ID 列表
    - `sex`：只保留一种性别：1/'male'、2/'female' 或 0/'unknown'
    - `phenotype`：只保留一类：'case'（2）、'control'（1）、'missing'（0 或 -9）、'observed'（非缺失）或 'quantitative'（其他任何值）
    - `founders_only`：布尔，默认 False：只保留父、母 ID 均为 0 的个体（奠基者）
    - `max_missing_rate`：[0,1] 内的浮点数：丢弃缺失基因型比例超过它的个体
    - `variant_indices`：要保留的 0 起始变异列下标列表（对 .map 文件应做同样的筛选）
    - `max_variant_missing_rate`：[0,1] 内的浮点数：丢弃缺失基因型比例超过它的变异列（在通过其他过滤条件的个体上计算）
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、`variant_count`、`kept_variant_indices`（保留下来的 0 起始变异列）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = filter_ped(PED, options={"founders_only": True, "max_variant_missing_rate": 0.4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 2,
    "variant_count": 1,
    "input_count": 3,
    "output_count": 2,
    "removed_count": 1,
    "kept_variant_indices": [
      0
    ],
    "bytes": 36,
    "sha256": "8c55765c8c50404f7b360d6c86a60f6191faa94deb8837c14f9e6b013f6d0f56",
    "text": "F1\tI1\t0\t0\t1\t2\tA\tA\nF1\tI2\t0\t0\t2\t1\tA\tG\n"
  }
}
```

## 42) filter_phylip

- **作用：** 按分类单元名以及每个分类单元/每个位点的缺失数据过滤 PHYLIP 文件。

- **处理方法：** 分类单元必须通过所有给定条件才保留：在 `keep` 中、不在 `exclude` 中、匹配 `name_pattern`（`re.search`，最长 256 个字符），且缺口/缺失/歧义比例不超过 `max_missing_fraction`。随后删去在保留的分类单元上缺失比例超过 `max_site_missing_fraction` 的比对列。

- **API：** `dnakit.io_formats.filter_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）
    - `keep`：要保留的分类单元名列表
    - `exclude`：要丢弃的分类单元名列表
    - `name_pattern`：分类单元名必须匹配的正则表达式（re.search，最长 256 个字符）
    - `max_missing_fraction`：[0,1] 内的浮点数：丢弃缺口/缺失/歧义比例超过它的分类单元
    - `max_site_missing_fraction`：[0,1] 内的浮点数：丢弃缺失比例（在保留的分类单元上计算）超过它的比对列

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = filter_phylip(PHYLIP, options={"max_site_missing_fraction": 0.2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 3,
    "output_count": 3,
    "removed_count": 0,
    "bytes": 62,
    "sha256": "89eeb0f0715875534c807a2240cb826134cee6d673b22d4197eb4d981a042300",
    "text": " 3 7\nA          ACGACGT\nB          ACGACGA\nC          ACGACGT\n"
  }
}
```

## 43) filter_twobit

- **作用：** 按序列名、长度和 N 比例过滤 2bit 文件，并返回过滤后的 2bit 内容。

- **处理方法：** 序列必须通过所有给定条件才被整条保留（记录按位原样复制，不重新编码）：在 `names` 中、不在 `exclude` 中、匹配 `name_pattern`（`re.search`，最长 256 个字符）、长度介于 `min_length` 和 `max_length` 之间，且 N 比例（N 碱基数 / 长度）不超过 `max_n_fraction`。

- **API：** `dnakit.io_formats.filter_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `names`：要保留的序列名列表
    - `exclude`：要丢弃的序列名列表
    - `name_pattern`：序列名必须匹配的正则表达式（re.search，最长 256 个字符）
    - `min_length`：整数：保留至少有这么多碱基的序列
    - `max_length`：整数：序列至多这么多碱基
    - `max_n_fraction`：[0,1] 内的浮点数：丢弃 N 比例超过它的序列
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`removed_count`、`sequence_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `data_base64`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `data_base64`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import filter_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = filter_twobit(TWOBIT, options={"min_length": 9}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 1,
    "sequence_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 58,
    "sha256": "2c8722ce953b5a3db475f68d7f8f16d9f09381af10b7e1472ac129e3eed9f6cf",
    "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMRcAAAAKAAAAAQAAAAQAAAACAAAAAQAAAAYAAAAEAAAAAAAAAJwJwA=="
  }
}
```

## 44) filter_zarr

- **作用：** 按数组路径和第一轴上的行过滤 Zarr 存储，写出新存储。

- **处理方法：** 把 `arrays`（路径或分组前缀；默认全部）选中且不在 `exclude` 中的数组复制到 `output`（必填）处的新存储；给出 `rows` 时只复制 0 轴上选中的行（`'start:stop[:step]'`、`{start, stop, step}` 或下标列表；所选的所有至少 1 维的数组在 0 轴上长度必须相同）。复制超过 `max_elements` 个元素会被拒绝。

- **API：** `dnakit.io_formats.filter_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `output`：新 zarr 存储的目录路径（必填）
    - `overwrite`：布尔，默认 False：替换 output 处已存在的 zarr 存储
    - `arrays`：要保留的数组路径（或分组前缀）
    - `exclude`：要丢弃的数组路径（或分组前缀）
    - `rows`：沿 0 轴保留的行：'start:stop[:step]'、{start, stop, step} 或下标列表（所选的所有至少 1 维的数组在 0 轴上长度必须相同）
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`（写入的数组个数）、`kept_rows`（给出 `rows` 时）和 `output_path`（不返回内容）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import filter_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = filter_zarr("demo.zarr", options={"output": "filtered.zarr", "overwrite": True, "arrays": ["dosage"], "rows": "1:3"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 1,
    "kept_rows": 2,
    "output_path": "filtered.zarr"
  }
}
```

## 45) get_file_checksum

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `get_file_checksum` 功能。候选编号 CAP-0353，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `get_file_checksum` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.get_file_checksum(input[required], algorithm[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`algorithm`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `get_file_checksum` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import get_file_checksum

result = get_file_checksum(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0353",
  "operation": "get_file_checksum",
  "result_type": "str",
  "result_preview": "'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) index

- **作用：** 为任意受支持格式的文件建立查找索引。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `index_nexus`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.index(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = index(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C d": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    },
    "trees": {
      "t1": 0
    }
  }
}
```

## 47) index_map

- **作用：** 为 PLINK `.map` 文件建立查找索引：变异 ID 对应的行号和字节偏移。

- **处理方法：** 把每个变异 ID 映射到它在（解压后）文本中的从 1 开始的行号、字节偏移、字节长度及其染色体和位置，并列出每条染色体的变异数、bp 范围和首行行号，从而无需重新解析即可定位变异。ID 出现在多行时，`entries` 只保留它的最后一行（`record_count` 仍统计每一行）；这类 ID 列在 `duplicate_ids` 中，并返回一条警告。

- **API：** `dnakit.io_formats.index_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"id"`）、`record_count`、`entries`（每个 ID 的 `line`、`offset`、`length`、`chrom`、`pos`）、`chromosomes`（每条染色体的 `variant_count`、`min_pos`、`max_pos`、`first_line`）和 `duplicate_ids`（出现在多行的 ID，按首次出现的顺序；没有则为空）。ID 重复时 `result` 还有 `warnings`（含一条消息的列表），同一条消息也在封装的 `warnings` 中。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
    "X\trs1\t0\t7\n"
)

envelope = index_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [
    "1 variant ID(s) occur on several lines (rs1): each entry keeps the last line of its ID, while record_count (5) counts every line."
  ],
  "result": {
    "format": "map",
    "key": "id",
    "record_count": 5,
    "entries": {
      "rs1": {
        "line": 5,
        "offset": 46,
        "length": 10,
        "chrom": "X",
        "pos": 7
      },
      "rs2": {
        "line": 2,
        "offset": 12,
        "length": 13,
        "chrom": "1",
        "pos": 50
      },
      "rs3": {
        "line": 3,
        "offset": 25,
        "length": 11,
        "chrom": "2",
        "pos": 10
      },
      "rs4": {
        "line": 4,
        "offset": 36,
        "length": 10,
        "chrom": "X",
        "pos": 5
      }
    },
    "chromosomes": {
      "1": {
        "variant_count": 2,
        "min_pos": 50,
        "max_pos": 100,
        "first_line": 1
      },
      "2": {
        "variant_count": 1,
        "min_pos": 10,
        "max_pos": 10,
        "first_line": 3
      },
      "X": {
        "variant_count": 2,
        "min_pos": 5,
        "max_pos": 7,
        "first_line": 4
      }
    },
    "duplicate_ids": [
      "rs1"
    ],
    "warnings": [
      "1 variant ID(s) occur on several lines (rs1): each entry keeps the last line of its ID, while record_count (5) counts every line."
    ]
  }
}
```

## 48) index_newick

- **作用：** 为 Newick 文件建立查找索引：树的字节偏移，以及包含每个叶节点的树。

- **处理方法：** 给出每棵树的从 0 开始的序号、在（解压后）文本中的字节偏移与字节长度、叶节点数，并对每个叶节点名列出包含它的树的序号。

- **API：** `dnakit.io_formats.index_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"tree index"`）、`record_count`、`trees`（每棵树的 `index`、`offset`、`length`、`tip_count`）和 `tips`（叶节点名到树序号列表）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = index_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "key": "tree index",
    "record_count": 2,
    "trees": [
      {
        "index": 0,
        "offset": 0,
        "length": 30,
        "tip_count": 4
      },
      {
        "index": 1,
        "offset": 31,
        "length": 12,
        "tip_count": 4
      }
    ],
    "tips": {
      "A": [
        0,
        1
      ],
      "B": [
        0,
        1
      ],
      "C": [
        0,
        1
      ],
      "D": [
        0,
        1
      ]
    }
  }
}
```

## 49) index_nexus

- **作用：** 为 NEXUS 文件建立查找索引：分类单元和树。

- **处理方法：** 把每个分类单元名映射到它的行序号、序列长度和非缺失残基数，把每个树名映射到它在树中的位置。NEXUS 矩阵可能是交错排布，所以不给字节偏移。

- **API：** `dnakit.io_formats.index_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"taxon name"`）、`record_count`、`entries`（每个分类单元的 `ordinal`、`length`、`residues`）和 `trees`（树名到位置）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = index_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C d": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    },
    "trees": {
      "t1": 0
    }
  }
}
```

## 50) index_paf

- **作用：** 为 PAF 文件建立查找索引：查询序列（或目标序列）名对应的行号、字节偏移和区间。

- **处理方法：** `by='query'`（默认）或 `by='target'` 时，把每个序列名映射到它的比对列表，每项含该行在（解压后）文本中的从 1 开始的行号、字节偏移、字节长度，以及该序列上的比对区间（`start`、`end`，同 PAF 为 0 起始、左闭右开）；列表按起点、终点和行号排序。

- **API：** `dnakit.io_formats.index_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `by`：'query'（默认）或 'target'：建立索引所用的序列名

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"query_name"` 或 `"target_name"`）、`record_count` 和 `entries`（名称到 `{line, offset, length, start, end}` 的列表）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = index_paf(PAF, options={"by": "target"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "key": "target_name",
    "record_count": 2,
    "entries": {
      "t1": [
        {
          "line": 1,
          "offset": 0,
          "length": 55,
          "start": 100,
          "end": 190
        }
      ],
      "t2": [
        {
          "line": 2,
          "offset": 55,
          "length": 34,
          "start": 0,
          "end": 50
        }
      ]
    }
  }
}
```

## 51) index_ped

- **作用：** 为 PLINK `.ped` 文件建立查找索引：个体对应的行号和字节偏移。

- **处理方法：** 把 `family:individual` 映射到它在（解压后）文本中的从 1 开始的行号、字节偏移和字节长度，并列出每个家系的个体。

- **API：** `dnakit.io_formats.index_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"family:individual"`）、`record_count`、`entries`（每个个体的 `line`、`offset`、`length`）和 `families`（家系 ID 到个体 ID 列表）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = index_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "key": "family:individual",
    "record_count": 3,
    "entries": {
      "F1:I1": {
        "line": 1,
        "offset": 0,
        "length": 22
      },
      "F1:I2": {
        "line": 2,
        "offset": 22,
        "length": 22
      },
      "F1:I3": {
        "line": 3,
        "offset": 44,
        "length": 25
      }
    },
    "families": {
      "F1": [
        "I1",
        "I2",
        "I3"
      ]
    }
  }
}
```

## 52) index_phylip

- **作用：** 为 PHYLIP 文件建立查找索引：分类单元的行序号、长度和残基数。

- **处理方法：** 把每个分类单元名映射到它的行序号、序列长度和非缺失残基数。交错排布的 PHYLIP 并不是每行一个分类单元，所以不给字节偏移。

- **API：** `dnakit.io_formats.index_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"taxon name"`）、`record_count` 和 `entries`（每个分类单元的 `ordinal`、`length`、`residues`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = index_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    }
  }
}
```

## 53) index_twobit

- **作用：** 为 2bit 文件建立查找索引：序列名及其文件偏移。

- **处理方法：** 只读取文件头和每条序列的记录头：每条序列的记录文件偏移、长度、N 块数与屏蔽块数，以及其打包数据的偏移。

- **API：** `dnakit.io_formats.index_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：


- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"sequence name"`）、`record_count` 和 `entries`（每条序列的 `offset`、`length`、`n_block_count`、`mask_block_count`、`packed_offset`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import index_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = index_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "key": "sequence name",
    "record_count": 2,
    "entries": {
      "s1": {
        "offset": 30,
        "length": 10,
        "n_block_count": 1,
        "mask_block_count": 1,
        "packed_offset": 62
      },
      "s2": {
        "offset": 65,
        "length": 8,
        "n_block_count": 0,
        "mask_block_count": 0,
        "packed_offset": 81
      }
    }
  }
}
```

## 54) index_zarr

- **作用：** 为 Zarr 存储建立查找索引：数组路径及其形状和分块网格。

- **处理方法：** 列出每个数组路径及其类型、形状、dtype、分块形状、分块网格和分块数（已初始化的和总数）。不读取分块数据。

- **API：** `dnakit.io_formats.index_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：


- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`key`（`"array path"`）、`record_count`、`entries`（路径到 `{shape, dtype, chunks, chunk_grid, nchunks, nchunks_initialized}`）、`groups` 和 `path`。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import index_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = index_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "key": "array path",
    "record_count": 2,
    "entries": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1
      },
      "pos": {
        "shape": [
          4
        ],
        "dtype": "int64",
        "chunks": [
          4
        ],
        "chunk_grid": [
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1
      }
    },
    "groups": [],
    "path": "demo.zarr"
  }
}
```

## 55) list_file_records

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `list_file_records` 功能。候选编号 CAP-0363，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `list_file_records` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.list_file_records(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `list_file_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import list_file_records

result = list_file_records(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0363",
  "operation": "list_file_records",
  "result_type": "list",
  "result_preview": "['>a', 'ACGT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 56) merge

- **作用：** 把同一种受支持格式的多个文件合并为一个。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `merge_phylip`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 输入为非空列表；由第一个元素决定格式，所有元素都必须是该格式。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.merge(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：同一种受支持格式的文件组成的非空列表（或 `{'inputs': [...]}`），每个文件的形式同对应的格式专用函数。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = merge([PHYLIP, PHYLIP], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 2,
    "bytes": 90,
    "sha256": "6ef4966d5f176941e22c4ba195bd108ca859a24e6fcce8f468255ce2bbc51031",
    "text": " 3 16\nA          ACGTACGTACGTACGT\nB          ACGTACGAACGTACGA\nC          ACG-ACGTACG-ACGT\n"
  }
}
```

## 57) merge_map

- **作用：** 把多个 `.map` 文件合并为一个 PLINK `.map` 文件。

- **处理方法：** 变异按输入顺序拼接。变异 ID 重复时，`duplicates='error'` 抛出 `InputFormatError`（`IO_FORMATS_DUPLICATE_ID`），`duplicates='first'` 则保留第一次出现的记录。`sort=True` 时合并结果按染色体（自然顺序）和 bp、再按 ID 排序。

- **API：** `dnakit.io_formats.merge_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `duplicates`：'error'（默认）或 'first'：变异 ID 重复出现时如何处理
    - `sort`：布尔，默认 False：把合并后的变异按染色体和 bp 排序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)
MAP_B = "2\trs3\t0\t10\n3\trs5\t0\t77\n"

envelope = merge_map([MAP, MAP_B], options={"duplicates": "first", "sort": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 5,
    "input_count": 2,
    "bytes": 57,
    "sha256": "fbdbc7656e124d28d44c918b51a2a367f6266809d601f238fcb2153f6689a483",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\n3\trs5\t0\t77\nX\trs4\t0\t5\n"
  }
}
```

## 58) merge_metadata_tables

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `merge_metadata_tables` 功能。候选编号 CAP-0366，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge_metadata_tables` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`merge_metadata`。

- **API：** `dnakit.io_formats.merge_metadata_tables(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `merge_metadata_tables` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import merge_metadata_tables

result = merge_metadata_tables(
    input=[[{"sample_id": "sample1"}], [{"sample_id": "sample2"}]]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0366",
  "operation": "merge_metadata_tables",
  "result_type": "list",
  "result_preview": "[{'sample_id': 'sample1'}, {'sample_id': 'sample2'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 59) merge_newick

- **作用：** 把多个 Newick 文件合并为一个多树 Newick 文件。

- **处理方法：** 所有输入的树按输入顺序拼接，每棵树一行规范化 Newick。`deduplicate=True` 时丢弃规范化 Newick 与前面某棵树相同的树。

- **API：** `dnakit.io_formats.merge_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `deduplicate`：布尔，默认 False：丢弃规范化 Newick 与前面某棵树相同的树

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = merge_newick([TREES, TREES], options={"deduplicate": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "bytes": 44,
    "sha256": "060f56a376c85eddacd58d52e718831c6b84a379f89394849619c36ed07806e7",
    "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"
  }
}
```

## 60) merge_nexus

- **作用：** 把多个 NEXUS 文件合并为一个 NEXUS 文件，按分类单元连接位点或追加分类单元。

- **处理方法：** `mode='concatenate'`（默认）构建超矩阵：分类单元取各输入按首次出现顺序的并集，同名分类单元的位点相连，某输入中缺少的分类单元用缺失符按该输入的位点数填充。`mode='append'` 追加所有输入的分类单元（名称重复抛出 `InputFormatError`；所有输入的位点数必须相同）。数据类型必须一致；各输入的树原样带入（重名时改名），不做剪枝，也不扩展到新增的分类单元。

- **API：** `dnakit.io_formats.merge_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `mode`：'concatenate'（默认：把同名分类单元的位点连接成超矩阵，某输入中缺少的分类单元用缺失字符填充）或 'append'（追加分类单元；位点数必须相等）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""
NEXUS_B = "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=3;\n  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n  MATRIX\n    A TTT\n    D GGG\n  ;\nEND;\n"

envelope = merge_nexus([NEXUS, NEXUS_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 5,
    "input_count": 2,
    "bytes": 268,
    "sha256": "18d97702dbe400e22c8bcaf861e10a7269a6320d60863079fe9414e00736a258",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=4 NCHAR=11;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGTACGTTTT\n    B      ACGTACGA???\n    'C d'  ACG-ACGT???\n    D      ????????GGG\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] ((A:1,B:1):1,'C d':2);\nEND;\n"
  }
}
```

## 61) merge_paf

- **作用：** 把多个 PAF 文件合并为一个 PAF 文件。

- **处理方法：** 比对按输入顺序拼接。`deduplicate=True` 时丢弃与前面完全相同的行（所有字段和标签都相同）。`sort_by`（`sort_paf` 的某个键）用于在合并之后排序。

- **API：** `dnakit.io_formats.merge_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `deduplicate`：布尔，默认 False：丢弃与前面完全相同的比对行
    - `sort_by`：可选排序键（同 sort_paf），应用于合并后的记录

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = merge_paf([PAF, PAF], options={"deduplicate": True, "sort_by": "mapq"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "input_count": 2,
    "bytes": 89,
    "sha256": "5ced75fd261d6ae80944159151cb702eb1c8258e15679a6fba89f91a2282bf99",
    "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 62) merge_ped

- **作用：** 把变异相同的多个 `.ped` 文件合并为一个 PLINK `.ped` 文件。

- **处理方法：** 各行按输入顺序拼接。所有文件的变异数必须相同，否则抛出 `InputFormatError`。（家系, 个体）重复时，`duplicates='error'` 抛出 `InputFormatError`（`IO_FORMATS_DUPLICATE_ID`），`duplicates='first'` 则保留第一行。函数不检查各文件的变异是否是同一批变异且顺序一致，这由调用者负责。

- **API：** `dnakit.io_formats.merge_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `duplicates`：'error'（默认）或 'first'：（家系, 个体）对重复出现时如何处理
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、`variant_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)
PED_B = "F2 J1 0 0 2 1 C C T T\n"

envelope = merge_ped([PED, PED_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 4,
    "variant_count": 2,
    "input_count": 2,
    "bytes": 91,
    "sha256": "32ddebc1aec0e9da60a191a2b5a8775a2b4c7151cf4d3be89f05c343daca7a90",
    "text": "F1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\nF2\tJ1\t0\t0\t2\t1\tC\tC\tT\tT\n"
  }
}
```

## 63) merge_phylip

- **作用：** 把多个 PHYLIP 比对合并为一个 PHYLIP 文件，按分类单元连接位点或追加分类单元。

- **处理方法：** `mode='concatenate'`（默认）构建超矩阵：分类单元取各输入按首次出现顺序的并集，同名分类单元的位点相连，某输入中缺少的分类单元用 `?` 按该输入的位点数填充。`mode='append'` 追加所有输入的分类单元（名称重复抛出 `InputFormatError`；所有输入的位点数必须相同）。数据类型必须一致。

- **API：** `dnakit.io_formats.merge_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）
    - `mode`：'concatenate'（默认：把同名分类单元的位点连接成超矩阵，某输入中缺少的分类单元用缺失字符填充）或 'append'（追加分类单元；位点数必须相等）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"
PHYLIP_B = " 2 3\nA  TTT\nD  GGG\n"

envelope = merge_phylip([PHYLIP, PHYLIP_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 4,
    "input_count": 2,
    "bytes": 98,
    "sha256": "613e33924aea618ca1521ba246371fbe922448952251c676c6d39c4a77ead034",
    "text": " 4 11\nA          ACGTACGTTTT\nB          ACGTACGA???\nC          ACG-ACGT???\nD          ????????GGG\n"
  }
}
```

## 64) merge_twobit

- **作用：** 把多个 2bit 文件合并为一个 2bit 文件。

- **处理方法：** 按输入顺序把所有输入的序列连接成一个版本 0 的文件；打包记录原样复制、不重新编码。名称重复时，`duplicates='error'` 抛出 `InputFormatError`（`IO_FORMATS_DUPLICATE_ID`），`duplicates='first'` 则保留第一条序列。`byte_order` 指定结果的字节序。

- **API：** `dnakit.io_formats.merge_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `duplicates`：'error'（默认）或 'first'：名称重复出现时如何处理
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`（合并的输入个数）、`sequence_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `data_base64`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `data_base64`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import merge_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = merge_twobit([TWOBIT, TWOBIT], options={"duplicates": "first"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "input_count": 2,
    "bytes": 83,
    "sha256": "12ae7a3592d786f9ae950efd4bc62172c1006c9fedd5815f798de219e356f5a1",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMR4AAAACczJBAAAACgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcAIAAAAAAAAAAAAAAAAAAAA9aA="
  }
}
```

## 65) merge_zarr

- **作用：** 把多个 Zarr 存储按相同路径的数组拼接，合并为一个存储。

- **处理方法：** 所有输入存储必须含相同的数组路径；相同路径的数组沿 `axis`（默认 0）按输入顺序拼接，除该轴外形状必须一致，dtype 用 `numpy.result_type` 提升。结果是 `output`（必填；`overwrite=True` 时替换已有存储）处的新存储。合并超过 `max_elements` 个元素会被拒绝。

- **API：** `dnakit.io_formats.merge_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：非空列表（或 `{'inputs': [...]}`），每一项都是Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `output`：新 zarr 存储的目录路径（必填）
    - `overwrite`：布尔，默认 False：替换 output 处已存在的 zarr 存储
    - `axis`：整数，默认 0：拼接所沿的轴
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`（写入的数组个数）、`input_count`、`axis` 和 `output_path`（不返回内容）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import merge_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = merge_zarr(["demo.zarr", "demo.zarr"], options={"output": "merged.zarr", "overwrite": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 2,
    "input_count": 2,
    "axis": 0,
    "output_path": "merged.zarr"
  }
}
```

## 66) normalize_sample_metadata

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `normalize_sample_metadata` 功能。候选编号 CAP-0374，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `normalize_sample_metadata` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.normalize_sample_metadata(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `normalize_sample_metadata` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import normalize_sample_metadata

result = normalize_sample_metadata(
    input=[{"sample_id": "sample1", "group": "A"}]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0374",
  "operation": "normalize_sample_metadata",
  "result_type": "list",
  "result_preview": "[{'sample_id': 'sample1', 'group': 'A'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 67) read_map

- **作用：** 把 PLINK `.map` 文件解析为结构化的变异记录。

- **处理方法：** 每个非空行必须有四个以空白分隔的列：染色体、变异 ID、遗传距离 cM（有限数值）和碱基位置 bp（>= 0 的整数）。空行被跳过，格式错误的行抛出 `InputFormatError`，记录保持文件中的顺序。最多返回 `max_records` 个变异；更大的文件请用 `stream_map`。

- **API：** `dnakit.io_formats.read_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count` 和 `records`；每条记录为 `{chrom, id, cm, pos}`（`cm` 为浮点数，`pos` 为整数）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = read_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "records": [
      {
        "chrom": "1",
        "id": "rs1",
        "cm": 0.0,
        "pos": 100
      },
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      },
      {
        "chrom": "2",
        "id": "rs3",
        "cm": 0.0,
        "pos": 10
      },
      {
        "chrom": "X",
        "id": "rs4",
        "cm": 0.0,
        "pos": 5
      }
    ]
  }
}
```

## 68) read_newick

- **作用：** 把 Newick 树文件解析为结构化的树。

- **处理方法：** 解析每棵以 `;` 结尾的树（标签、枝长、`[...]` 注释）得到节点树；记录含重新序列化的规范化 Newick 字符串、叶节点名和先序排列的扁平节点表（`id`、`parent`、`name`、`length`，根节点的 `parent` 为 `null`）。语法错误抛出 `InputFormatError`；每棵树最多 `max_nodes` 个节点，文件最多 `max_records` 棵树。

- **API：** `dnakit.io_formats.read_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）
    - `max_nodes`：整数，默认 2000000：每棵树的最大节点数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count` 和 `records`；每条记录为 `{index, newick, tips, nodes}`，含规范化的 Newick 字符串、叶节点名和先序排列的扁平节点表 `{id, parent, name, length}`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = read_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "records": [
      {
        "index": 0,
        "newick": "((A:1,B:2):0.5,(C:1,D:1):0.5);",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": null,
            "length": 0.5
          },
          {
            "id": 2,
            "parent": 1,
            "name": "A",
            "length": 1.0
          },
          {
            "id": 3,
            "parent": 1,
            "name": "B",
            "length": 2.0
          },
          {
            "id": 4,
            "parent": 0,
            "name": null,
            "length": 0.5
          },
          {
            "id": 5,
            "parent": 4,
            "name": "C",
            "length": 1.0
          },
          {
            "id": 6,
            "parent": 4,
            "name": "D",
            "length": 1.0
          }
        ]
      },
      {
        "index": 1,
        "newick": "(A,B,(C,D));",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": "A",
            "length": null
          },
          {
            "id": 2,
            "parent": 0,
            "name": "B",
            "length": null
          },
          {
            "id": 3,
            "parent": 0,
            "name": null,
            "length": null
          },
          {
            "id": 4,
            "parent": 3,
            "name": "C",
            "length": null
          },
          {
            "id": 5,
            "parent": 3,
            "name": "D",
            "length": null
          }
        ]
      }
    ]
  }
}
```

## 69) read_nexus

- **作用：** 把 NEXUS 文件解析为结构化的序列和树。

- **处理方法：** 解析 TAXA、DATA/CHARACTERS 和 TREES 块（忽略其他块）：分类单元标签、矩阵（顺序或交错排布，展开 `matchchar`）、声明的数据类型、缺失符和缺口符，以及树（应用 TRANSLATE 表，`[&R]`/`[&U]` 记为 `rooted` 标志）。`[]` 中的注释被丢弃。带引号的分类单元名（`'C d'`）会去掉引号。

- **API：** `dnakit.io_formats.read_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`taxon_count`、`site_count`、`datatype`、`taxa`、`tree_count`、`missing`、`gap`、`layout` 和 `records`（序列记录 `{kind: "sequence", name, sequence}`，其后是树记录 `{kind: "tree", name, rooted, newick}`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = read_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 4,
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "taxa": [
      "A",
      "B",
      "C d"
    ],
    "tree_count": 1,
    "missing": "?",
    "gap": "-",
    "layout": "sequential",
    "records": [
      {
        "kind": "sequence",
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "kind": "sequence",
        "name": "B",
        "sequence": "ACGTACGA"
      },
      {
        "kind": "sequence",
        "name": "C d",
        "sequence": "ACG-ACGT"
      },
      {
        "kind": "tree",
        "name": "t1",
        "rooted": true,
        "newick": "((A:1,B:1):1,'C d':2);"
      }
    ]
  }
}
```

## 70) read_paf

- **作用：** 把 PAF 比对文件解析为结构化的比对记录。

- **处理方法：** 每个非空行按制表符切成 12 个必需列（查询名/长度/起点/终点、链方向、目标名/长度/起点/终点、匹配碱基数、比对块长度、比对质量）；其余列是可选的类 SAM 标签 `TAG:TYPE:VALUE`，以 `tags` 返回，并保留其声明的类型。格式错误的行抛出 `InputFormatError`。最多返回 `max_records` 条比对。

- **API：** `dnakit.io_formats.read_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count` 和 `records`；每条记录含 12 个 PAF 字段以及 `tags`（标签到 `{type, value}`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = read_paf(PAF, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "records": [
      {
        "query_name": "q1",
        "query_length": 100,
        "query_start": 0,
        "query_end": 90,
        "strand": "+",
        "target_name": "t1",
        "target_length": 1000,
        "target_start": 100,
        "target_end": 190,
        "matches": 80,
        "block_length": 90,
        "mapq": 60,
        "tags": {
          "tp": {
            "type": "A",
            "value": "P"
          },
          "cg": {
            "type": "Z",
            "value": "90M"
          }
        }
      },
      {
        "query_name": "q2",
        "query_length": 50,
        "query_start": 0,
        "query_end": 50,
        "strand": "-",
        "target_name": "t2",
        "target_length": 500,
        "target_start": 0,
        "target_end": 50,
        "matches": 40,
        "block_length": 50,
        "mapq": 10,
        "tags": {}
      }
    ]
  }
}
```

## 71) read_ped

- **作用：** 把 PLINK `.ped` 文件解析为结构化的个体记录。

- **处理方法：** 每行有 `6 + 2k` 个以空白分隔的列（家系、个体、父亲、母亲、性别、表型，随后每个变异两个等位基因），且每行的 `k` 相同。`sex` 读作整数；`pheno` 保持原字符串。基因型矩阵大小受 `max_cells`（个体 x 变异）限制，行数受 `max_records` 限制。

- **API：** `dnakit.io_formats.read_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）
    - `max_cells`：整数，默认 20000000：拒绝基因型单元格（个体 x 变异）更多的文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`variant_count` 和 `records`；每条记录为 `{fid, iid, pat, mat, sex, pheno, genotypes}`，`genotypes` 是 `[等位基因1, 等位基因2]` 对的列表（`sex` 为整数，`pheno` 为字符串；等位基因 `0` 表示缺失）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = read_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 3,
    "variant_count": 2,
    "records": [
      {
        "fid": "F1",
        "iid": "I1",
        "pat": "0",
        "mat": "0",
        "sex": 1,
        "pheno": "2",
        "genotypes": [
          [
            "A",
            "A"
          ],
          [
            "C",
            "G"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I2",
        "pat": "0",
        "mat": "0",
        "sex": 2,
        "pheno": "1",
        "genotypes": [
          [
            "A",
            "G"
          ],
          [
            "0",
            "0"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I3",
        "pat": "I1",
        "mat": "I2",
        "sex": 1,
        "pheno": "-9",
        "genotypes": [
          [
            "G",
            "G"
          ],
          [
            "C",
            "C"
          ]
        ]
      }
    ]
  }
}
```

## 72) read_phylip

- **作用：** 把 PHYLIP 比对文件解析为结构化的序列。

- **处理方法：** 读取宽松 PHYLIP（名称是第一个以空白分隔的词；残基块中可含空格）或严格 PHYLIP（`input_strict_names=True`：名称占前 10 个字符）；除非给出 `input_layout`，否则自动识别顺序和交错排布。首行 `ntax nchar` 必须与数据一致。`datatype` 指定用于校验的字符集。

- **API：** `dnakit.io_formats.read_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `max_records`：整数，默认 200000：记录数超过该值时抛出异常而不是返回（请改用 stream_* 分页）
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`taxon_count`、`site_count`、`datatype` 和 `records`（每条记录为 `{name, sequence}`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = read_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "records": [
      {
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "name": "B",
        "sequence": "ACGTACGA"
      },
      {
        "name": "C",
        "sequence": "ACG-ACGT"
      }
    ]
  }
}
```

## 73) read_twobit

- **作用：** 把 UCSC 2bit 文件解码为序列，或解码其中一条序列的一个区间。

- **处理方法：** 把打包的 2 位碱基解码为文本：N 块变为 `N`，屏蔽块变为小写（`soft_mask=False` 时输出大写）。`names` 选择序列（默认全部）；`start`/`end` 给出恰好一条指定序列的 0 起始、左闭右开区间 `[start, end)`，否则抛出 `ConfigurationError`。解码超过 `max_bases` 个碱基会被拒绝。

- **API：** `dnakit.io_formats.read_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `names`：要读取的序列名列表（默认全部）
    - `start`：整数：0 起始的开始位置（仅限单个序列名）
    - `end`：整数：0 起始、不含终点的结束位置（仅限单个序列名）
    - `soft_mask`：布尔，默认 True：把被屏蔽的碱基输出为小写
    - `max_bases`：整数，默认 50000000：超过该碱基数则拒绝解码

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count` 和 `records`；每条记录为 `{name, start, end, length, sequence}`，其中 `length` 为整条序列的长度，`[start, end)` 为解码的区间。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import read_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = read_twobit(TWOBIT, options={"names": ["s1"], "start": 2, "end": 8}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "records": [
      {
        "name": "s1",
        "start": 2,
        "end": 8,
        "length": 10,
        "sequence": "GTNNac"
      }
    ],
    "record_count": 1
  }
}
```

## 74) read_zarr

- **作用：** 把 Zarr 目录存储读取为结构化的数组。

- **处理方法：** 以只读方式打开存储，对每个数组返回形状、dtype、分块形状、分块网格、分块数、压缩器、zarr 格式、未压缩大小和属性，并且（除非 `include_data=False`）以嵌套列表返回数值。`array` 选择一个数组（默认全部）；`start`/`stop` 在单个数组的 0 轴上选取行。读取超过 `max_elements` 个元素会被拒绝。

- **API：** `dnakit.io_formats.read_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `array`：要读取的数组路径（默认全部数组）
    - `include_data`：布尔，默认 True：返回元素值
    - `start`：整数：沿 0 轴的第一行（仅限单个数组）
    - `stop`：整数：结束行（不含）
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`path`、`arrays`（路径到 `{shape, dtype, chunks, chunk_grid, nchunks, nchunks_initialized, compressors, zarr_format, nbytes, attrs, data}`；`include_data=False` 时没有 `data`）、`groups`、`attrs`、`record_count` 和 `records`（相同的元数据，不含 `data`，带 `path`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import read_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = read_zarr("demo.zarr", options={"array": "dosage", "start": 1, "stop": 3}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "path": "demo.zarr",
    "arrays": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {},
        "data": [
          [
            1,
            1
          ],
          [
            0,
            0
          ]
        ]
      }
    },
    "groups": [],
    "attrs": {},
    "record_count": 1,
    "records": [
      {
        "path": "dosage",
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {}
      }
    ]
  }
}
```

## 75) repack_compressed_file

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `repack_compressed_file` 功能。候选编号 CAP-0384，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `repack_compressed_file` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.repack_compressed_file(input[required], level[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`level`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `repack_compressed_file` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import repack_compressed_file

result = repack_compressed_file(
    input=__import__("gzip").compress(b"ACGT")
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0384",
  "operation": "repack_compressed_file",
  "result_type": "bytes",
  "result_preview": "<bytes length=24>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 76) resume_download

- **作用：** 通过 HTTP(S) 下载文件，用 Range 请求续传中断的传输。

- **处理方法：** 先把 URL 下载到 `<output>.part`。若该分片文件已有字节，则保留它们，并用 `Range: bytes=<n>-` 请求其余部分，同时用保存在 `<output>.part.meta` 中的 ETag（或 Last-Modified）通过 `If-Range` 校验。服务器返回 `200`（忽略了 Range 或资源已变化）时从字节 0 重新下载（`restarted: true`）；`416` 且总大小一致表示分片文件已经完整。`206` 响应必须从请求的偏移开始。收到的大小等于声明的大小，并且（若给出 `expected_sha256`）摘要一致（不一致时删除分片文件并抛出 `DownloadError`）后，分片文件被重命名为 `output`。传输中途断开时保留 `.part` 文件并抛出 `DownloadError`（`DOWNLOAD_LENGTH_MISMATCH` 或 `DOWNLOAD_NETWORK_ERROR`），再次调用即可续传。只接受 `http` 和 `https` URL（重定向也一样）；超过 `max_bytes` 抛出 `DownloadError`。这是本页中唯一使用网络的函数。

- **API：** `dnakit.io_formats.resume_download(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：URL（`http://` 或 `https://`）字符串，或映射 `{'url': ..., 'output': ...}`。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `url`：http(s) URL（也可作为 input 传入）
    - `output`：目标文件路径（input 不是带 'output' 的映射时必填）
    - `overwrite`：布尔，默认 False：替换已存在的目标文件
    - `chunk_size`：整数，默认 65536：每步读取的字节数
    - `timeout`：秒数，默认 30：套接字超时
    - `max_bytes`：整数，默认 1073741824：拒绝下载超过该字节数的资源
    - `expected_sha256`：完成的文件必须具有的十六进制摘要（不一致时删除未完成的文件）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`output_path`、`bytes`（完成文件的大小）、`sha256`、`restarted`（下载是否不得不从字节 0 重新开始）、`complete`（成功时恒为 `true`）和 `resumed_from`（传输从哪个字节偏移续传，全新下载为 0）。HTTP 或网络故障抛出 `DownloadError`；选项无效抛出 `ConfigurationError`；非 http(s) 的 URL 抛出 `InputFormatError`。

- **示例代码：**

```python
import hashlib
import http.server
import json
import threading
from pathlib import Path

from dnakit.io_formats import resume_download

PAYLOAD = bytes(range(256)) * 40  # the remote file: 10240 bytes


class RangeHandler(http.server.BaseHTTPRequestHandler):
    """A tiny local server that honours 'Range: bytes=<n>-' requests."""

    def log_message(self, *args):
        pass

    def do_GET(self):
        start = 0
        header = self.headers.get("Range")
        if header:
            start = int(header.split("=")[1].rstrip("-"))
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(PAYLOAD) - 1}/{len(PAYLOAD)}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(PAYLOAD) - start))
        self.send_header("ETag", '"v1"')
        self.end_headers()
        self.wfile.write(PAYLOAD[start:])


server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f"http://127.0.0.1:{server.server_address[1]}/file.bin"

# an interrupted earlier transfer: 3000 bytes are already stored in out.bin.part
Path("out.bin.part").write_bytes(PAYLOAD[:3000])
Path("out.bin.part.meta").write_text(json.dumps({"etag": '"v1"'}))

envelope = resume_download(
    url,
    options={"output": "out.bin", "expected_sha256": hashlib.sha256(PAYLOAD).hexdigest()},
    verbose=False,
)
server.shutdown()
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
assert Path("out.bin").read_bytes() == PAYLOAD  # the finished file equals the remote file
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "output_path": "out.bin",
    "bytes": 10240,
    "sha256": "e96760a87768717bcebcfd25ddc7d46b4dbc95a4b0014def080c08539f7d90d0",
    "restarted": false,
    "complete": true,
    "resumed_from": 3000
  }
}
```

## 77) sort

- **作用：** 用对应格式专用排序函数的顺序对任意受支持格式的文件排序。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `sort_map`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.sort(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = sort(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 46,
    "sha256": "4299c54d1a304472628eaa9018e807d127e591bc315c46a7272a56910661277b",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
  }
}
```

## 78) sort_map

- **作用：** 按染色体和位置（或按 ID、遗传距离）对 PLINK `.map` 文件中的变异排序。

- **处理方法：** `key='position'`（默认）依次按染色体自然顺序（1..22、X、Y、XY、M/MT，其他名称按字母序）、bp、ID 排序；`key='id'` 先按 ID，再按染色体和 bp；`key='cm'` 先按遗传距离，再按染色体和 bp。`reverse=True` 使顺序倒过来。

- **API：** `dnakit.io_formats.sort_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `key`：'position'（默认：染色体自然顺序，然后 bp，然后 ID）、'id' 或 'cm'
    - `reverse`：布尔，默认 False：倒序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = sort_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 46,
    "sha256": "4299c54d1a304472628eaa9018e807d127e591bc315c46a7272a56910661277b",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
  }
}
```

## 79) sort_newick

- **作用：** 对 Newick 文件中每棵树的分支排序（梯形化或字母序）。

- **处理方法：** 对每个节点的子节点重新排序，树的拓扑不变：`order='ascending'`（默认）把叶节点少的分支排在前（梯形化），`'descending'` 相反，`'alphabetical'` 按各分支最小的叶节点名排序；并列时按该名称决定。文件中树的顺序保持不变。

- **API：** `dnakit.io_formats.sort_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `order`：'ascending'（默认：叶节点少的分支在前）、'descending' 或 'alphabetical'（按每个分支最小的叶节点名）；并列时按该名称决定

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = sort_newick(TREES, options={"order": "descending"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 44,
    "sha256": "654ac0beaf1d726fa4cc69552def8d54c2df37576b814570a46d7560004d7032",
    "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n((C,D),A,B);\n"
  }
}
```

## 80) sort_nexus

- **作用：** 按名称或缺失数据比例对 NEXUS 文件的分类单元（矩阵行）排序。

- **处理方法：** `key='name'`（默认）按分类单元名排序；`key='missing'` 按缺口/缺失/歧义字符的比例、再按名称排序。`reverse=True` 使顺序倒过来。位点和树保持不变。

- **API：** `dnakit.io_formats.sort_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `key`：'name'（默认）或 'missing'
    - `reverse`：布尔，默认 False：倒序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = sort_nexus(NEXUS, options={"key": "missing", "reverse": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 235,
    "sha256": "d105b30c92592d6e0bc698e1fd80bfebaadaf5d21acb7a30b952b2b45dd18f83",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    'C d'  ACG-ACGT\n    B      ACGTACGA\n    A      ACGTACGT\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] ((A:1,B:1):1,'C d':2);\nEND;\n"
  }
}
```

## 81) sort_paf

- **作用：** 按目标序列、查询序列、比对质量、一致性、块长度或匹配数对 PAF 文件中的比对排序。

- **处理方法：** `key='target'`（默认）按目标名、目标起点和终点排序；`query` 按查询名、查询起点和终点；`mapq`、`identity`（matches / block_length）、`block_length` 和 `matches` 按对应数值。`reverse=True` 使顺序倒过来。

- **API：** `dnakit.io_formats.sort_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `key`：'target'（默认：目标名、起点、终点）、'query'（查询名、起点、终点）、'mapq'、'identity'、'block_length' 或 'matches'
    - `reverse`：布尔，默认 False：倒序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = sort_paf(PAF, options={"key": "mapq"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 89,
    "sha256": "5ced75fd261d6ae80944159151cb702eb1c8258e15679a6fba89f91a2282bf99",
    "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 82) sort_ped

- **作用：** 按家系、个体 ID、性别或表型对 PLINK `.ped` 文件中的个体排序。

- **处理方法：** `key='family'`（默认）先按家系 ID 再按个体 ID；`individual` 先按个体 ID 再按家系 ID；`sex` 按性别码，再按家系和个体；`phenotype` 按表型的数值，再按家系和个体（解析器要求每行的表型都是数值）。ID 按字符串比较。`reverse=True` 使顺序倒过来。

- **API：** `dnakit.io_formats.sort_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `key`：'family'（默认：先家系 ID 后个体 ID）、'individual'、'sex' 或 'phenotype'
    - `reverse`：布尔，默认 False：倒序
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`variant_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = sort_ped(PED, options={"key": "phenotype"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 3,
    "variant_count": 2,
    "input_count": 3,
    "output_count": 3,
    "bytes": 69,
    "sha256": "e96809700ea4ad197aae18b2174944b94d7ec502c4febe6d1549d56d1d78f358",
    "text": "F1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\n"
  }
}
```

## 83) sort_phylip

- **作用：** 按名称或缺失数据比例对 PHYLIP 文件的分类单元（行）排序。

- **处理方法：** `key='name'`（默认）按分类单元名排序；`key='missing'` 按缺口/缺失/歧义字符的比例、再按名称排序。`reverse=True` 使顺序倒过来。位点保持不变。

- **API：** `dnakit.io_formats.sort_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）
    - `key`：'name'（默认）或 'missing'（先按缺失比例，再按名称）
    - `reverse`：布尔，默认 False：倒序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = sort_phylip(PHYLIP, options={"key": "missing", "reverse": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 3,
    "output_count": 3,
    "bytes": 65,
    "sha256": "05bc6f3276e92b5bef45beaa26c3ffeb9d2af1c36f3e0ec87386a1152c9f7677",
    "text": " 3 8\nC          ACG-ACGT\nB          ACGTACGA\nA          ACGTACGT\n"
  }
}
```

## 84) sort_twobit

- **作用：** 按名称、染色体自然顺序或长度对 2bit 文件中的序列排序。

- **处理方法：** 输出索引中的序列被重新排序，打包记录按位原样复制。`key='name'`（默认）是普通名称顺序；`'natural'` 是染色体顺序（`chr1`..`chr22`、`X`、`Y`、`M`，其他名称按字母序）；`'length'` 按序列长度。`reverse=True` 使顺序倒过来。

- **API：** `dnakit.io_formats.sort_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output`：输出文件路径（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `key`：'name'（默认）、'natural'（chr1..22、X、Y、M 顺序）或 'length'
    - `reverse`：布尔，默认 False：倒序
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`input_count`、`output_count`、`sequence_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `data_base64`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `data_base64`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import sort_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = sort_twobit(TWOBIT, options={"key": "length"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 83,
    "sha256": "651601b14986cb13fa6b5342f96e6fbb9ed4474860233892def1156920d74518",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMh4AAAACczEwAAAACAAAAAAAAAAAAAAAAAAAAPWgCgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcA="
  }
}
```

## 85) sort_zarr

- **作用：** 按键数组的值对 Zarr 存储中所有数组的行（0 轴）排序，写出新存储。

- **处理方法：** `by`（必填）指定的一维数组的值给出稳定排序（`reverse=True` 时反转升序结果）；第一维长度等于键数组长度的每个数组沿 0 轴按该顺序置换，其他数组原样复制。结果是 `output`（必填）处的新存储。排序超过 `max_elements` 个元素会被拒绝。

- **API：** `dnakit.io_formats.sort_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `output`：新 zarr 存储的目录路径（必填）
    - `overwrite`：布尔，默认 False：替换 output 处已存在的 zarr 存储
    - `by`：1 维键数组的路径（必填）；各行按其值稳定排序
    - `reverse`：布尔，默认 False：反转升序结果
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`（写入的数组个数）、`key`、`permuted_arrays`、`copied_unchanged` 和 `output_path`（不返回内容）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import sort_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = sort_zarr("demo.zarr", options={"output": "sorted.zarr", "overwrite": True, "by": "pos"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 2,
    "key": "pos",
    "permuted_arrays": [
      "dosage",
      "pos"
    ],
    "copied_unchanged": [],
    "output_path": "sorted.zarr"
  }
}
```

## 86) split_map

- **作用：** 把 PLINK `.map` 文件按染色体拆分，或拆成每块 N 个变异。

- **处理方法：** `by='chromosome'`（默认）时每条染色体一个分片，按自然顺序（1..22、X、Y、XY、M/MT、其他）排列并以染色体名命名。`by='size'` 时按每块 `size` 个变异连续切分，命名为 `0001`、`0002`……，此时 `size` 必填。每个分片都是完整的 `.map` 文件。

- **API：** `dnakit.io_formats.split_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `by`：'chromosome'（默认，每条染色体一个分片）或 'size'
    - `size`：整数：by='size' 时每个分片的变异数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = split_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "part_count": 3,
    "parts": [
      {
        "name": "1",
        "record_count": 2,
        "bytes": 25,
        "sha256": "40be7a615faec1b5f39695fe44a7f975a678b81bb09a0ac0aa2608389cbefd6e",
        "text": "1\trs1\t0\t100\n1\trs2\t0.5\t50\n"
      },
      {
        "name": "2",
        "record_count": 1,
        "bytes": 11,
        "sha256": "fc5af41f3a73f2c542af0c2b9e7be799905f0ae6459c62f6162daa2088473cbd",
        "text": "2\trs3\t0\t10\n"
      },
      {
        "name": "X",
        "record_count": 1,
        "bytes": 10,
        "sha256": "ac4cc4dc84e656a5e8cf48bab597529c20e7a0af3d275866c00dccc7470660bb",
        "text": "X\trs4\t0\t5\n"
      }
    ]
  }
}
```

## 87) split_newick

- **作用：** 把 Newick 文件拆成每个含 N 棵树的文件。

- **处理方法：** 每连续 `trees_per_part` 棵树（默认 1）构成一个分片，命名为 `0001`、`0002`……；每个分片都是完整的 Newick 文件。

- **API：** `dnakit.io_formats.split_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `trees_per_part`：整数，默认 1：每个输出分片的树数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = split_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "part_count": 2,
    "parts": [
      {
        "name": "0001",
        "record_count": 1,
        "bytes": 31,
        "sha256": "2443c552d427d6bcde12f7e9de1eea8a374ca83fae68f6b5be7f2cdb9de0b008",
        "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n"
      },
      {
        "name": "0002",
        "record_count": 1,
        "bytes": 13,
        "sha256": "2cead55aa9365339f55531d03a9973092b4c1d2a66503e615296d53c826612b6",
        "text": "(A,B,(C,D));\n"
      }
    ]
  }
}
```

## 88) split_nexus

- **作用：** 按位点块、分类单元组、命名的位点分区或树组拆分 NEXUS 文件。

- **处理方法：** 必须恰好给出 `site_size`（连续的位点块，分片名 `sites_<起>-<止>`）、`taxa_size`（连续的分类单元组，`taxa_0001`……）、`partitions`（标签到 `start-end[,start-end...]`，1 起始、含两端的位点区间，每个标签一个分片）或 `trees_per_part`（树组，`trees_0001`……；文件必须含树）之一。比对分片只含矩阵，树分片含 TAXA 块和树。每个分片都是完整的 NEXUS 文件。

- **API：** `dnakit.io_formats.split_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `site_size`：整数：按这么多个位点切成连续的块
    - `taxa_size`：整数：按这么多个分类单元切成连续的组
    - `partitions`：映射 标签 -> 'start-end[,start-end...]'（1 起始、含两端的位点区间）
    - `trees_per_part`：整数：每个分片的树数（无矩阵的文档）；site_size、taxa_size、partitions、trees_per_part 四者必须恰好给一个

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = split_nexus(NEXUS, options={"site_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "part_count": 2,
    "parts": [
      {
        "name": "sites_1-4",
        "record_count": 3,
        "bytes": 162,
        "sha256": "f730e6e646487331208af5613fdcf90ac1bcd885e831dc963b4c860a4132b190",
        "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=4;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGT\n    B      ACGT\n    'C d'  ACG-\n    ;\nEND;\n"
      },
      {
        "name": "sites_5-8",
        "record_count": 3,
        "bytes": 162,
        "sha256": "4f5f5fcb316930bf01371cd7d93e8f613a41f263e25d7d12b1a2522983797e59",
        "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=4;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGT\n    B      ACGA\n    'C d'  ACGT\n    ;\nEND;\n"
      }
    ]
  }
}
```

## 89) split_paf

- **作用：** 把 PAF 文件按查询序列、目标序列或链方向拆分，或拆成每块 N 条比对。

- **处理方法：** `by='query'`（默认）、`'target'` 或 `'strand'` 时每个不同取值一个分片，按首次出现的顺序排列并以该取值命名；`by='size'` 时按每块 `size` 条比对连续切分，命名为 `0001`、`0002`……（此时 `size` 必填）。每个分片都是完整的 PAF 文件。

- **API：** `dnakit.io_formats.split_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `by`：'query'（默认）、'target'、'strand' 或 'size'
    - `size`：整数：by='size' 时每个分片的比对数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = split_paf(PAF, options={"by": "strand"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "part_count": 2,
    "parts": [
      {
        "name": "+",
        "record_count": 1,
        "bytes": 55,
        "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
        "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
      },
      {
        "name": "-",
        "record_count": 1,
        "bytes": 34,
        "sha256": "5bec66460b177fb4f757cf63b7dc30daee9f94c45e100137e47512b6dc9b5bd3",
        "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
      }
    ]
  }
}
```

## 90) split_ped

- **作用：** 把 PLINK `.ped` 文件按家系拆分，或拆成每块 N 个个体。

- **处理方法：** `by='family'`（默认）时每个家系一个分片，按首次出现的顺序排列并以家系 ID 命名。`by='size'` 时按每块 `size` 个个体连续切分，命名为 `0001`、`0002`……，此时 `size` 必填。每个分片保留全部变异列，都是完整的 `.ped` 文件。

- **API：** `dnakit.io_formats.split_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `by`：'family'（默认，每个家系一个分片）或 'size'
    - `size`：整数：by='size' 时每个分片的个体数
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`variant_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = split_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "part_count": 1,
    "parts": [
      {
        "name": "F1",
        "record_count": 3,
        "variant_count": 2,
        "bytes": 69,
        "sha256": "5083df15f1e129b12dfd6258db03b28858011c33cce5a983955ed2f3ace7a2b7",
        "text": "F1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\n"
      }
    ]
  }
}
```

## 91) split_phylip

- **作用：** 按位点块、分类单元组或命名的位点分区拆分 PHYLIP 文件。

- **处理方法：** 必须恰好给出 `site_size`（连续的位点块，分片名 `sites_<起>-<止>`）、`taxa_size`（连续的分类单元组，`taxa_0001`……）或 `partitions`（标签到 `start-end[,start-end...]`，1 起始、含两端的位点区间，每个标签一个分片）之一。PHYLIP 文件不含树，所以 NEXUS 和 Newick 拆分函数的 `trees_per_part` 选项在这里会被拒绝，并抛出说明 PHYLIP 没有树的 `ConfigurationError`。每个分片都是完整的 PHYLIP 文件。

- **API：** `dnakit.io_formats.split_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）
    - `site_size`：整数：按这么多个位点切成连续的块
    - `taxa_size`：整数：按这么多个分类单元切成连续的组
    - `partitions`：映射 标签 -> 'start-end[,start-end...]'（1 起始、含两端的位点区间）。site_size、taxa_size、partitions 三者必须恰好给一个

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = split_phylip(PHYLIP, options={"site_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "part_count": 2,
    "parts": [
      {
        "name": "sites_1-4",
        "record_count": 3,
        "bytes": 53,
        "sha256": "878fecc7ee77357a027fd67f53ca9a126c676fd3fb4de77519f8d7a15695fcd2",
        "text": " 3 4\nA          ACGT\nB          ACGT\nC          ACG-\n"
      },
      {
        "name": "sites_5-8",
        "record_count": 3,
        "bytes": 53,
        "sha256": "05cbfd292180c3c0d25dc5286ea79a552121d97f8d753d01a20fb56bd2c1995d",
        "text": " 3 4\nA          ACGT\nB          ACGA\nC          ACGT\n"
      }
    ]
  }
}
```

## 92) split_twobit

- **作用：** 把 2bit 文件拆成每个含 N 条序列的 2bit 文件。

- **处理方法：** 每连续 `sequences_per_part` 条序列（默认 1）写入一个单独的版本 0 文件；每个分片都是完整的 2bit 文件。`sequences_per_part=1` 时分片以序列名命名，否则命名为 `0001`、`0002`……。

- **API：** `dnakit.io_formats.split_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `output_dir`：把分片写成 part_0001_<label><扩展名> 的目录（可选）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `sequences_per_part`：整数，默认 1：每个输出分片的序列数
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`；每个分片含 `name`、`record_count`、`sequence_count`、`bytes`、`sha256` 及其内容（`text`，二进制格式为 `data_base64`）；给出 `options['output_dir']` 时每个分片写入 `part_0001_<name><扩展名>`、`part_0002_...`，并返回 `output_path` 而不再返回内容。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import split_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = split_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "part_count": 2,
    "parts": [
      {
        "name": "s1",
        "record_count": 1,
        "sequence_count": 1,
        "bytes": 58,
        "sha256": "2c8722ce953b5a3db475f68d7f8f16d9f09381af10b7e1472ac129e3eed9f6cf",
        "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMRcAAAAKAAAAAQAAAAQAAAACAAAAAQAAAAYAAAAEAAAAAAAAAJwJwA=="
      },
      {
        "name": "s2",
        "record_count": 1,
        "sequence_count": 1,
        "bytes": 41,
        "sha256": "99bac5b233c34dd7e0e7a506a9b1b176814d4680c658ebf52792b0681b440e44",
        "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMhcAAAAIAAAAAAAAAAAAAAAAAAAA9aA="
      }
    ]
  }
}
```

## 93) split_zarr

- **作用：** 沿 0 轴把 Zarr 存储拆成每个含 N 行的存储。

- **处理方法：** 所有数组的 0 轴被切成每块 `rows_per_part` 行（必填）的连续块，每块一个存储，写入 `output_dir`（必填），文件名为 `part_0001_rows_<起>-<止>.zarr` 等；所有至少 1 维的数组在 0 轴上长度必须相同（否则抛出 `ConfigurationError`；可先用 `filter_zarr` 选出兼容的数组），没有维度的数组复制到每个分片。拆分超过 `max_elements` 个元素会被拒绝。

- **API：** `dnakit.io_formats.split_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `output_dir`：接收 part_0001_<label>.zarr 存储的目录（必填）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `rows_per_part`：整数（必填）：每个分片沿 0 轴的行数
    - `max_elements`：整数，默认 5000000：拒绝载入多于该数量的数组元素

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`part_count` 和 `parts`（每个分片的 `name`（如 `rows_0-3`）、`record_count`（数组个数）、`row_count` 和 `output_path`；分片总是写入 `output_dir`）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import split_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = split_zarr("demo.zarr", options={"rows_per_part": 3, "output_dir": "parts"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "part_count": 2,
    "parts": [
      {
        "name": "rows_0-3",
        "record_count": 2,
        "row_count": 3,
        "output_path": "parts/part_0001_rows_0-3.zarr"
      },
      {
        "name": "rows_3-4",
        "record_count": 2,
        "row_count": 1,
        "output_path": "parts/part_0002_rows_3-4.zarr"
      }
    ]
  }
}
```

## 94) stream

- **作用：** 返回任意受支持格式文件中的一页记录，并给出下一页的游标。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `stream_map`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.stream(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = stream(MAP, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "chrom": "1",
        "id": "rs1",
        "cm": 0.0,
        "pos": 100
      },
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      }
    ]
  }
}
```

## 95) stream_map

- **作用：** 返回 PLINK `.map` 文件中的一页变异记录，并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 限制）并返回下标区间 `[start, start + limit)` 的记录；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（每条记录为 `{chrom, id, cm, pos}`）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = stream_map(MAP, options={"start": 1, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "start": 1,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 3,
    "done": false,
    "records": [
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      },
      {
        "chrom": "2",
        "id": "rs3",
        "cm": 0.0,
        "pos": 10
      }
    ]
  }
}
```

## 96) stream_newick

- **作用：** 返回 Newick 文件中的一页树，并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 和 `max_nodes` 限制）并返回下标区间 `[start, start + limit)` 的树；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（每条记录为 `{index, newick, tips, nodes}`，含规范化的 Newick 字符串、叶节点名和先序排列的扁平节点表 `{id, parent, name, length}`）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = stream_newick(TREES, options={"start": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "start": 1,
    "limit": 1000,
    "returned": 1,
    "total": 2,
    "next_start": null,
    "done": true,
    "records": [
      {
        "index": 1,
        "newick": "(A,B,(C,D));",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": "A",
            "length": null
          },
          {
            "id": 2,
            "parent": 0,
            "name": "B",
            "length": null
          },
          {
            "id": 3,
            "parent": 0,
            "name": null,
            "length": null
          },
          {
            "id": 4,
            "parent": 3,
            "name": "C",
            "length": null
          },
          {
            "id": 5,
            "parent": 3,
            "name": "D",
            "length": null
          }
        ]
      }
    ]
  }
}
```

## 97) stream_nexus

- **作用：** 返回 NEXUS 文件中的一页序列（之后是树），并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 限制）并返回合并列表（先矩阵各行，后树）中下标区间 `[start, start + limit)` 的记录；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（序列记录 `{kind: "sequence", name, sequence}`，其后是树记录 `{kind: "tree", name, rooted, newick}`）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = stream_nexus(NEXUS, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "kind": "sequence",
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "kind": "sequence",
        "name": "B",
        "sequence": "ACGTACGA"
      }
    ]
  }
}
```

## 98) stream_paf

- **作用：** 返回 PAF 文件中的一页比对记录，并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 限制）并返回下标区间 `[start, start + limit)` 的记录；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（每条记录含 12 个 PAF 字段以及 `tags`（标签到 `{type, value}`））；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = stream_paf(PAF, options={"limit": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "start": 0,
    "limit": 1,
    "returned": 1,
    "total": 2,
    "next_start": 1,
    "done": false,
    "records": [
      {
        "query_name": "q1",
        "query_length": 100,
        "query_start": 0,
        "query_end": 90,
        "strand": "+",
        "target_name": "t1",
        "target_length": 1000,
        "target_start": 100,
        "target_end": 190,
        "matches": 80,
        "block_length": 90,
        "mapq": 60,
        "tags": {
          "tp": {
            "type": "A",
            "value": "P"
          },
          "cg": {
            "type": "Z",
            "value": "90M"
          }
        }
      }
    ]
  }
}
```

## 99) stream_ped

- **作用：** 返回 PLINK `.ped` 文件中的一页个体记录，并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 和 `max_cells` 限制）并返回下标区间 `[start, start + limit)` 的记录；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（每条记录为 `{fid, iid, pat, mat, sex, pheno, genotypes}`，`genotypes` 是 `[等位基因1, 等位基因2]` 对的列表）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = stream_ped(PED, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 3,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "fid": "F1",
        "iid": "I1",
        "pat": "0",
        "mat": "0",
        "sex": 1,
        "pheno": "2",
        "genotypes": [
          [
            "A",
            "A"
          ],
          [
            "C",
            "G"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I2",
        "pat": "0",
        "mat": "0",
        "sex": 2,
        "pheno": "1",
        "genotypes": [
          [
            "A",
            "G"
          ],
          [
            "0",
            "0"
          ]
        ]
      }
    ]
  }
}
```

## 100) stream_phylip

- **作用：** 返回 PHYLIP 文件中的一页序列，并给出下一页的游标。

- **处理方法：** 解析文件（受 `max_bytes` 限制）并返回下标区间 `[start, start + limit)` 的序列；以 `start=next_start` 反复调用，直到 `done` 为 true。

- **API：** `dnakit.io_formats.stream_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（每条记录为 `{name, sequence}`）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = stream_phylip(PHYLIP, options={"start": 1, "limit": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "start": 1,
    "limit": 1,
    "returned": 1,
    "total": 3,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "name": "B",
        "sequence": "ACGTACGA"
      }
    ]
  }
}
```

## 101) stream_records

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `stream_records` 功能。候选编号 CAP-0411，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `stream_records` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`trim_record`。

- **API：** `dnakit.io_formats.stream_records(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `stream_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import stream_records

result = stream_records(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0411",
  "operation": "stream_records",
  "result_type": "list",
  "result_preview": "['>a', 'ACGT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 102) stream_twobit

- **作用：** 返回 2bit 文件中的一页序列窗口，并给出下一页的游标。

- **处理方法：** 每条序列被切成 `window` 个碱基（默认 10000）的窗口，所有序列的窗口按文件顺序构成记录列表；`start` 和 `limit` 以窗口计，只解码本页的窗口（每页最多 `max_bases` 个碱基）。

- **API：** `dnakit.io_formats.stream_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）
    - `window`：整数，默认 10000：每个窗口的碱基数（分页的单位）
    - `soft_mask`：布尔，默认 True：用小写表示软屏蔽碱基
    - `max_bases`：整数，默认 50000000：每页最多碱基数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`start`、`limit`、`returned`、`total`（窗口总数）、`window`、`next_start`、`done` 和 `records`；每条记录是一个窗口的 `{name, start, end, sequence}`（最后一页的 `next_start` 为 `null`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import stream_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = stream_twobit(TWOBIT, options={"window": 6, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "window": 6,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "name": "s1",
        "start": 0,
        "end": 6,
        "sequence": "ACGTNN"
      },
      {
        "name": "s1",
        "start": 6,
        "end": 10,
        "sequence": "acgt"
      }
    ]
  }
}
```

## 103) stream_zarr

- **作用：** 返回 Zarr 存储中某个数组的一页行，并给出下一页的游标。

- **处理方法：** 沿 0 轴逐页读取某个数组的行：每条记录是一行（嵌套列表），`start`/`limit` 以行计。存储含多个数组时必须给出 `array`。每页最多读取 `max_elements` 个元素。

- **API：** `dnakit.io_formats.stream_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `start`：整数，默认 0：本页第一条记录的下标
    - `limit`：整数，默认 1000：每页记录数（最大 100000）
    - `array`：要分页读取的数组路径（存储中有多个数组时必填）
    - `max_elements`：整数，默认 5000000：每页最多元素数

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`array`、`start`、`limit`、`returned`、`total`、`next_start`、`done` 和 `records`（数组的行，为嵌套列表）；`next_start` 是取下一页时传给 `options['start']` 的值，最后一页为 `null`，此时 `done` 为 `true`。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import stream_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = stream_zarr("demo.zarr", options={"array": "dosage", "start": 2, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "array": "dosage",
    "start": 2,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": null,
    "done": true,
    "records": [
      [
        0,
        0
      ],
      [
        2,
        1
      ]
    ]
  }
}
```

## 104) summarize

- **作用：** 汇总任意受支持格式的文件。

- **处理方法：** 格式取自 `options['format']`（`map`、`newick`、`nexus`、`paf`、`ped`、`phylip`、`twobit` 或 `zarr`；`fasta` 仅可作为 `convert` 的来源），否则自动识别：先看路径（后缀 `.map`、`.ped`、`.paf`、`.nwk`/`.newick`/`.tre`/`.tree`、`.nex`/`.nexus`/`.nxs`、`.phy`/`.phylip`、`.2bit`、`.fa`/`.fasta`/`.fna`，末尾的 `.gz` 被忽略；目录视为 zarr），否则看内容的前 64 KiB（2bit 签名、`#NEXUS`、以 `(` 开头为 Newick、以 `>` 开头为 FASTA、首行至少 12 个制表符分隔列且第 2 到 4 列为整数为 PAF、`ntax nchar` 头为 PHYLIP、四列且末列为整数为 MAP、至少六列且列数为偶数为 PED）。随后由与对应格式专用函数（例如 `summarize_twobit`）相同的引擎执行，因此选项和 `result` 都与该函数一致。无法识别格式时抛出 `InputFormatError`（`IO_FORMATS_UNKNOWN_FORMAT`）。 对 FASTA 输入执行 `convert` 以外的操作会抛出 `ConfigurationError`。

- **API：** `dnakit.io_formats.summarize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意受支持格式的文件，形式同对应的格式专用函数（路径、文本或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。除下列键外，还接受所识别格式对应的格式专用函数的选项。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `format`：输入格式名（map、newick、nexus、paf、ped、phylip、twobit、zarr）；省略时按后缀/内容自动识别

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：与所用格式专用函数的 `result` 完全相同，其 `format` 字段说明用了哪一个（FASTA 输入的 `convert` 中 `source_format` 为 `fasta`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import summarize, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = summarize(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "sequence_count": 2,
    "total_length": 18,
    "n_bases": 2,
    "soft_masked_bases": 4,
    "gc_fraction": 0.5,
    "version": 0,
    "byte_order": "little",
    "sequences": [
      {
        "name": "s1",
        "length": 10,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 2,
        "soft_masked": 4,
        "gc_fraction": 0.5
      },
      {
        "name": "s2",
        "length": 8,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 0,
        "soft_masked": 0,
        "gc_fraction": 0.5
      }
    ]
  }
}
```

## 105) summarize_map

- **作用：** 汇总 PLINK `.map` 文件：变异数与染色体数、bp 范围和重复 ID。

- **处理方法：** 按染色体统计变异数及其最小/最大 bp（染色体按自然顺序 1..22、X、Y、XY、M/MT、其他排列），统计出现多次的 ID 个数、cM 的取值范围、cM 是否全为 0，以及文件是否已按染色体和位置排序。

- **API：** `dnakit.io_formats.summarize_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`variant_count`、`chromosome_count`、`chromosomes`（每条染色体的 `variant_count`、`min_pos`、`max_pos`）、`duplicate_id_count`、`cm_min`、`cm_max`、`all_cm_zero` 和 `sorted`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = summarize_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "variant_count": 4,
    "chromosome_count": 3,
    "chromosomes": {
      "1": {
        "variant_count": 2,
        "min_pos": 50,
        "max_pos": 100
      },
      "2": {
        "variant_count": 1,
        "min_pos": 10,
        "max_pos": 10
      },
      "X": {
        "variant_count": 1,
        "min_pos": 5,
        "max_pos": 5
      }
    },
    "duplicate_id_count": 0,
    "cm_min": 0.0,
    "cm_max": 0.5,
    "all_cm_zero": false,
    "sorted": false
  }
}
```

## 106) summarize_newick

- **作用：** 汇总 Newick 文件：每棵树的节点数、深度、高度、总枝长和分辨率。

- **处理方法：** 每棵树统计：叶节点数、内部节点数和节点总数；根的度数；`rooted`（根恰有两个子节点）；`fully_resolved` 标志；以边数计的深度；高度 = 根到叶的最长路径长度（缺失的枝长按 0 计）；总枝长（没有任何枝长时为 `null`）；负枝长个数、无名叶节点数和重复叶节点名个数。

- **API：** `dnakit.io_formats.summarize_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`tree_count`、`trees`（每棵树的 `tip_count`、`internal_node_count`、`node_count`、`root_degree`、`rooted`、`fully_resolved`、`max_depth_edges`、`height`、`total_branch_length`、`has_branch_lengths`、`negative_branch_count`、`unnamed_tip_count`、`duplicate_tip_count`）、`tip_count_min` 和 `tip_count_max`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = summarize_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "tree_count": 2,
    "trees": [
      {
        "tip_count": 4,
        "internal_node_count": 3,
        "node_count": 7,
        "root_degree": 2,
        "rooted": true,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 2.5,
        "total_branch_length": 6.0,
        "has_branch_lengths": true,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      },
      {
        "tip_count": 4,
        "internal_node_count": 2,
        "node_count": 6,
        "root_degree": 3,
        "rooted": false,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 0.0,
        "total_branch_length": null,
        "has_branch_lengths": false,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      }
    ],
    "tip_count_min": 4,
    "tip_count_max": 4
  }
}
```

## 107) summarize_nexus

- **作用：** 汇总 NEXUS 文件：比对统计（恒定位点、可变位点、简约信息位点、组成）和树的统计。

- **处理方法：** 比对统计同 `summarize_phylip`（一列的无歧义字符只有一种状态为恒定位点，多种状态为可变位点，至少有两种状态各在至少两个分类单元中出现为简约信息位点；缺口、缺失符和歧义码被忽略，核苷酸中 `U` 按 `T` 计），并附上每棵树的 `summarize_newick` 统计。

- **API：** `dnakit.io_formats.summarize_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`taxon_count`、`site_count`、`datatype`、`equal_length`、`constant_sites`、`variable_sites`、`parsimony_informative_sites`、`all_missing_sites`、`total_characters`、`missing_fraction`、`taxon_missing_fraction`（`min`、`mean`、`max`）、`duplicate_sequence_count`、`composition`、`gc_content`（核苷酸）；有树时还有 `tree_count` 和 `trees`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = summarize_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "equal_length": true,
    "constant_sites": 7,
    "variable_sites": 1,
    "parsimony_informative_sites": 0,
    "all_missing_sites": 0,
    "total_characters": 24,
    "missing_fraction": 0.041666666666666664,
    "taxon_missing_fraction": {
      "min": 0.0,
      "mean": 0.041666666666666664,
      "max": 0.125
    },
    "duplicate_sequence_count": 0,
    "composition": {
      "-": 1,
      "A": 7,
      "C": 6,
      "G": 6,
      "T": 4
    },
    "gc_content": 0.5217391304347826,
    "tree_count": 1,
    "trees": [
      {
        "tip_count": 3,
        "internal_node_count": 2,
        "node_count": 5,
        "root_degree": 2,
        "rooted": true,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 2.0,
        "total_branch_length": 5.0,
        "has_branch_lengths": true,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      }
    ]
  }
}
```

## 108) summarize_paf

- **作用：** 汇总 PAF 文件：比对数、一致性、比对碱基数、链方向和比对质量直方图。

- **处理方法：** 统计比对数以及不同的查询序列/目标序列个数；比对碱基数 = 各比对块长度之和；单条比对的一致性 = matches / block_length，给出均值、最小值和最大值，`overall_identity` = 总匹配数 / 总块长度（按碱基加权）；`mapq` 是每个比对质量取值对应的比对条数。空输入时各比值为 `null`。

- **API：** `dnakit.io_formats.summarize_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`alignment_count`、`query_count`、`target_count`、`strand`（`+`、`-`）、`aligned_bases`、`matching_bases`、`overall_identity`、`mean_identity`、`min_identity`、`max_identity`、`mapq`（取值到条数）和 `mean_mapq`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = summarize_paf(PAF, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "alignment_count": 2,
    "query_count": 2,
    "target_count": 2,
    "strand": {
      "+": 1,
      "-": 1
    },
    "aligned_bases": 140,
    "matching_bases": 120,
    "overall_identity": 0.8571428571428571,
    "mean_identity": 0.8444444444444444,
    "min_identity": 0.8,
    "max_identity": 0.8888888888888888,
    "mapq": {
      "10": 1,
      "60": 1
    },
    "mean_mapq": 35.0
  }
}
```

## 109) summarize_ped

- **作用：** 汇总 PLINK `.ped` 文件：个体数、家系数、变异数、性别、表型类别和基因型缺失率。

- **处理方法：** 统计个体数、家系数、变异数和奠基者数（父母均为 `0`）；性别分为男（1）、女（2）、未知（其他）；表型分为病例（2）、对照（1）、缺失（0 或 -9）、数量性状（其他任何值）。任一等位基因为 `0` 即该基因型为缺失；`genotype_missing_rate` 为缺失基因型数 / 全部基因型数，`individual_missing_rate` 给出各个体缺失率的最小值、均值和最大值。

- **API：** `dnakit.io_formats.summarize_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`individual_count`、`family_count`、`variant_count`、`founder_count`、`sex`（`male`、`female`、`unknown`）、`phenotype`（`case`、`control`、`missing`、`quantitative`）、`genotype_missing_rate` 和 `individual_missing_rate`（`min`、`mean`、`max`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = summarize_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "individual_count": 3,
    "family_count": 1,
    "variant_count": 2,
    "founder_count": 2,
    "sex": {
      "male": 2,
      "female": 1,
      "unknown": 0
    },
    "phenotype": {
      "case": 1,
      "control": 1,
      "missing": 1,
      "quantitative": 0
    },
    "genotype_missing_rate": 0.16666666666666666,
    "individual_missing_rate": {
      "min": 0.0,
      "mean": 0.16666666666666666,
      "max": 0.5
    }
  }
}
```

## 110) summarize_phylip

- **作用：** 汇总 PHYLIP 文件：比对统计（恒定位点、可变位点、简约信息位点、组成、缺失数据）。

- **处理方法：** 一列的无歧义字符只有一种状态为恒定位点，有多种状态为可变位点，至少有两种不同的无歧义状态各在至少两个分类单元中出现为简约信息位点（缺口、缺失符和歧义码被忽略；DNA/RNA 中 `U` 按 `T` 计；没有任何无歧义字符的列计入 `all_missing_sites`）。核苷酸数据的 `gc_content` = (G + C) / (A + C + G + T + U)；`missing_fraction` 和各分类单元的比例统计 `datatype` 对应的缺口、缺失和歧义字符；同时统计重复序列。

- **API：** `dnakit.io_formats.summarize_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：输入（解压后）字节数上限
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`taxon_count`、`site_count`、`datatype`、`equal_length`、`constant_sites`、`variable_sites`、`parsimony_informative_sites`、`all_missing_sites`、`total_characters`、`missing_fraction`、`taxon_missing_fraction`（`min`、`mean`、`max`）、`duplicate_sequence_count`、`composition`（字符到计数）和 `gc_content`（核苷酸）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = summarize_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "equal_length": true,
    "constant_sites": 7,
    "variable_sites": 1,
    "parsimony_informative_sites": 0,
    "all_missing_sites": 0,
    "total_characters": 24,
    "missing_fraction": 0.041666666666666664,
    "taxon_missing_fraction": {
      "min": 0.0,
      "mean": 0.041666666666666664,
      "max": 0.125
    },
    "duplicate_sequence_count": 0,
    "composition": {
      "-": 1,
      "A": 7,
      "C": 6,
      "G": 6,
      "T": 4
    },
    "gc_content": 0.5217391304347826
  }
}
```

## 111) summarize_twobit

- **作用：** 汇总 2bit 文件：每条序列的长度、碱基组成、N 与软屏蔽碱基数、GC 比例。

- **处理方法：** 每条序列：长度、A、C、G、T（不区分大小写）和 N 的个数、软屏蔽（小写）碱基数，以及 `gc_fraction = (G + C) / (A + C + G + T)`；文件合计把各序列的长度、N 碱基数和屏蔽碱基数相加，并在所有序列上计算 `gc_fraction`。

- **API：** `dnakit.io_formats.summarize_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：


- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`sequence_count`、`total_length`、`n_bases`、`soft_masked_bases`、`gc_fraction`、`version`、`byte_order` 和 `sequences`（每条序列的 `name`、`length`、`A`、`C`、`G`、`T`、`N`、`soft_masked`、`gc_fraction`）。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import summarize_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = summarize_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "sequence_count": 2,
    "total_length": 18,
    "n_bases": 2,
    "soft_masked_bases": 4,
    "gc_fraction": 0.5,
    "version": 0,
    "byte_order": "little",
    "sequences": [
      {
        "name": "s1",
        "length": 10,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 2,
        "soft_masked": 4,
        "gc_fraction": 0.5
      },
      {
        "name": "s2",
        "length": 8,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 0,
        "soft_masked": 0,
        "gc_fraction": 0.5
      }
    ]
  }
}
```

## 112) summarize_zarr

- **作用：** 汇总 Zarr 存储：数组形状、dtype、分块和总量。

- **处理方法：** 统计分组数和数组数，并汇总元素数、未压缩字节数和已初始化的分块数；每个数组一项元数据。不读取分块数据。

- **API：** `dnakit.io_formats.summarize_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：


- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`path`、`kind`（`"group"` 或 `"array"`）、`array_count`、`group_count`、`total_elements`、`total_nbytes`、`initialized_chunks`、`total_chunks`、`attrs` 和 `arrays`（每个数组的形状、dtype、分块、分块网格、分块数、压缩器、zarr 格式、`nbytes`、`attrs`）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import summarize_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = summarize_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "path": "demo.zarr",
    "kind": "group",
    "array_count": 2,
    "group_count": 0,
    "total_elements": 12,
    "total_nbytes": 96,
    "initialized_chunks": 2,
    "total_chunks": 2,
    "attrs": {},
    "arrays": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {}
      },
      "pos": {
        "shape": [
          4
        ],
        "dtype": "int64",
        "chunks": [
          4
        ],
        "chunk_grid": [
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 32,
        "attrs": {}
      }
    }
  }
}
```

## 113) validate

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `validate` 功能。候选编号 CAP-0423，来自 25 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（arrow、bam、bcf、bed、bedgraph、bigbed、bigwig、bim、cram、csv、embl、fam 等，另有 13 种）。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`generic_operation_name`。已有同名函数证据：`dnakit.standardize.validate|validate|validate_set`。

- **API：** `dnakit.io_formats.validate(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `validate` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import validate

result = validate(
    input='>a\nACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0423",
  "operation": "validate",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'size_bytes': 8, 'line_count': 2, 'error': None}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 114) validate_map

- **作用：** 校验 PLINK `.map` 文件并报告错误与警告。

- **处理方法：** 每个非空行需要四个以空白分隔的列：染色体、变异 ID、遗传距离 cM（有限数值）和碱基位置 bp（>= 0 的整数）。一行的每个问题各占一条记录，形如 `line N: field 'NAME' (column C) has invalid value 'VALUE': WHY`（例如 cM 为 `x` 时是 `field 'cM' (column 3) has invalid value 'x': not a finite number`，cM 和 bp 都无效的行给出两条记录）；列数不对的行给出一条 `line N: expected 4 columns ...`。变异 ID 重复属于错误；（染色体, bp）重复属于警告。

- **API：** `dnakit.io_formats.validate_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.map` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`（`errors` 为空时为 true）、`errors`、`warnings`（`line N: ...` 形式的消息列表）和 `variant_count`（成功解析的变异数）。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_map

BAD = "1\trs1\t0\t100\n1\trs1\t0\t200\n2\trs3\tx\tabc\n2\trs9\t0\t10\n"

envelope = validate_map(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "valid": false,
    "errors": [
      "line 3: field 'cM' (column 3) has invalid value 'x': not a finite number",
      "line 3: field 'bp' (column 4) has invalid value 'abc': not an integer"
    ],
    "warnings": []
  }
}
```

## 115) validate_metadata_schema

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `validate_metadata_schema` 功能。候选编号 CAP-0425，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate_metadata_schema` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`validate_metadata`。

- **API：** `dnakit.io_formats.validate_metadata_schema(input[required], required[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`required`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `validate_metadata_schema` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import validate_metadata_schema

result = validate_metadata_schema(
    input={"sample_id": "sample1"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0425",
  "operation": "validate_metadata_schema",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'missing': [], 'keys': ['sample_id']}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 116) validate_newick

- **作用：** 校验 Newick 树文件并报告错误与警告。

- **处理方法：** 每棵以 `;` 结尾的树必须括号平衡、标签格式正确、枝长为有限数值（语法错误以 `tree N: ...` 报告）。叶节点名重复、叶节点无名属于错误；负枝长和单子（一元）内部节点属于警告。

- **API：** `dnakit.io_formats.validate_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 树文件（一棵或多棵以 `;` 结尾的树），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`（`tree N: ...` 形式的消息列表）、`tree_count` 和 `tip_count`（所有树的叶节点总数）。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_newick

BAD = "((A:1,A:2):0.5,(C:-1,D:1));\n"

envelope = validate_newick(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "valid": false,
    "errors": [
      "tree 1: 1 duplicated tip name(s)"
    ],
    "warnings": [
      "tree 1: 1 negative branch length(s)"
    ],
    "tree_count": 1,
    "tip_count": 4
  }
}
```

## 117) validate_nexus

- **作用：** 校验 NEXUS 文件并报告错误与警告。

- **处理方法：** 检查 `#NEXUS` 头、块结构、声明的 NTAX/NCHAR 与矩阵是否一致、各行等长、分类单元唯一，以及核苷酸/蛋白数据的字符集；树中不在矩阵里的叶节点以警告报告。

- **API：** `dnakit.io_formats.validate_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：NEXUS 文件（`#NEXUS`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`、`taxon_count`、`site_count` 和 `tree_count`。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_nexus

BAD = "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=5;\n  FORMAT DATATYPE=DNA;\n  MATRIX\n    A ACGTA\n    B ACG\n  ;\nEND;\n"

envelope = validate_nexus(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "valid": false,
    "errors": [
      "sequences have different lengths [3, 5]; an alignment needs equal lengths",
      "header declares 5 characters but sequences have length [3, 5]"
    ],
    "warnings": [],
    "taxon_count": 2,
    "site_count": 5,
    "tree_count": 0
  }
}
```

## 118) validate_paf

- **作用：** 校验 PAF 文件并报告错误与警告。

- **处理方法：** 每个非空行需要 12 个以制表符分隔的必需列，长度/坐标/计数为整数，查询和目标上均满足 `start <= end <= length`，`matches <= block_length`，链方向为 `+` 或 `-`，`mapq` 在 [0, 255] 内；可选列必须形如 `TAG:TYPE:VALUE`。一行的每个问题各占一条记录（单个字段为 `line N: field 'NAME' (column C) has invalid value 'VALUE': WHY`，字段之间的关系为 `line N: query coordinates must satisfy ...`），所以一行可以产生多条记录；少于 12 列的行给出一条记录。

- **API：** `dnakit.io_formats.validate_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PAF（Pairwise mApping Format）文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`（`line N: ...` 形式的消息列表），文件有效时还有 `alignment_count`。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_paf

BAD = "q1\t100\t0\t90\t*\tt1\t1000\t100\t190\t95\t90\t300\n"

envelope = validate_paf(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "valid": false,
    "errors": [
      "line 1: field 'strand' (column 5) has invalid value '*': must be '+' or '-'",
      "line 1: matches must satisfy 0 <= matches <= block_length",
      "line 1: field 'mapq' (column 12) has invalid value '300': must be in [0, 255]"
    ],
    "warnings": []
  }
}
```

## 119) validate_ped

- **作用：** 校验 PLINK `.ped` 文件并报告错误与警告。

- **处理方法：** 每行需要 `6 + 2k` 个以空白分隔的列，且每行的 `k` 相同（否则报 `line N: has ... variants but earlier rows have ...`），性别为整数、表型为数值。一行的每个问题各占一条记录，形如 `line N: field 'NAME' (column C) has invalid value 'VALUE': WHY`（因此性别和表型都无效的行给出两条记录，变异数不符也不会掩盖它们）。（家系, 个体）重复属于错误；性别码不在 {0, 1, 2} 内、以及父母不是同一家系中的个体属于警告；等位基因 `0` 表示缺失。给出 `variant_count` 时，文件的变异数必须等于该值，否则报错。

- **API：** `dnakit.io_formats.validate_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PLINK `.ped` 文件，文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `variant_count`：整数：期望的变异数（例如对应 .map 文件的行数）
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`（`line N: ...` 形式的消息列表）、`individual_count` 和 `variant_count`。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 two 1x A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = validate_ped(PED, options={"variant_count": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "valid": false,
    "errors": [
      "line 2: field 'sex' (column 5) has invalid value 'two': not an integer",
      "line 2: field 'phenotype' (column 6) has invalid value '1x': not a finite number"
    ],
    "warnings": []
  }
}
```

## 120) validate_phylip

- **作用：** 校验 PHYLIP 比对文件并报告错误与警告。

- **处理方法：** 检查首行计数与各行是否一致（若各行无法在顺序或交错排布下构成 `ntax` 个恰有 `nchar` 个字符的分类单元，则报一条错误）、序列是否等长、分类单元名是否唯一，以及（按所选 `datatype`）每个字符是否属于该字符集。

- **API：** `dnakit.io_formats.validate_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：PHYLIP 比对文件（首行 `ntax nchar`），文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的输入
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于校验字符集和统计缺失字符（核苷酸为缺口、?、N；蛋白为 X）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：PHYLIP 输入的排布
    - `input_strict_names`：布尔，默认 False：输入名称占每行前 10 个字符（否则取第一个以空白分隔的词）

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`，文件可解析时还有 `taxon_count`、`site_count` 和 `tree_count`（PHYLIP 恒为 0）。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_phylip

BAD = " 3 4\nA  ACGT\nA  ACGA\nC  ACJT\n"

envelope = validate_phylip(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "valid": false,
    "errors": [
      "duplicate taxon names: ['A']",
      "C: characters ['J'] are not valid for datatype 'dna' (set the datatype option if this is a different alphabet)"
    ],
    "warnings": [],
    "taxon_count": 3,
    "site_count": 4,
    "tree_count": 0
  }
}
```

## 121) validate_twobit

- **作用：** 校验 UCSC 2bit 文件并报告错误与警告。

- **处理方法：** 检查签名（两种字节序）、版本 0、保留字段、索引偏移、名称唯一且非空、N 块和屏蔽块有序、非空、不重叠且落在序列内，以及打包数据是否完整地落在文件内。完全不是 2bit 的文件会作为错误报告（`sequence_count` 为 0）。

- **API：** `dnakit.io_formats.validate_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：UCSC 2bit 基因组文件，可以是路径（`str` 或 `Path`）或原始 `bytes`（不是已存在路径的 `str` 会被拒绝，因为文本无法容纳二进制文件）；两种字节序都接受。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：


- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`、`sequence_count`，头部可读时还有 `byte_order` 和 `file_size`。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。

- **示例代码：**

```python
import base64
import json

from dnakit.io_formats import validate_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = validate_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "valid": true,
    "errors": [],
    "warnings": [],
    "sequence_count": 2,
    "byte_order": "little",
    "file_size": 83
  }
}
```

## 122) validate_zarr

- **作用：** 校验 Zarr 目录存储并报告错误与警告。

- **处理方法：** 打开存储，检查每个数组的形状、分块形状和分块数，并解码每个已初始化数组的第一个分块，以发现损坏的编解码器或被截断的分块文件（大于 `max_chunk_bytes` 的首个分块不解码）。

- **API：** `dnakit.io_formats.validate_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Zarr 目录存储的路径（`str` 或 `Path`；zarr v2 或 v3；需要可选的 `zarr>=3` 和 `numpy` 包）；不接受文本内容、bytes 和 gzip。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_chunk_bytes`：整数，默认 67108864：大于该字节数的首个分块不做解码

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`valid`、`errors`、`warnings`、`array_count`、`group_count` 和 `path`。内容问题记录在 `errors`/`warnings` 中而不是抛出；完全无法读取的输入（文件不存在、gzip 或 UTF-8 无效）仍抛出 `InputFormatError`；选项无效时抛出 `ConfigurationError`。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。

- **示例代码：**

```python
import json

from dnakit.io_formats import validate_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = validate_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "valid": true,
    "errors": [],
    "warnings": [],
    "array_count": 2,
    "group_count": 0,
    "path": "demo.zarr"
  }
}
```

## 123) verify_checksum

- **作用：** 为“格式转换、压缩与文件管理”类分析提供可复用的 `verify_checksum` 功能。候选编号 CAP-0433，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `verify_checksum` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/io_formats/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.io_formats.verify_checksum(input[required], expected[required], algorithm[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`、`expected`；可选参数：`algorithm`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `verify_checksum` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.io_formats import verify_checksum

result = verify_checksum(
    input='>a\nACGT\n',
    expected=__import__("hashlib").sha256(b">a\nACGT\n").hexdigest()
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0433",
  "operation": "verify_checksum",
  "result_type": "bool",
  "result_preview": "True",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 124) write_map

- **作用：** 把变异记录序列化为 PLINK `.map` 文件，或规范化已有的 `.map` 内容。

- **处理方法：** 每条记录被整理成四列，以制表符分隔、每个变异一行、保持输入顺序（cM 最多 15 位有效数字）。已有的 `.map` 文本或路径会先被解析再重新序列化，因此以空白分隔的输入会变为制表符分隔。除非给出 `output`，否则不写磁盘。

- **API：** `dnakit.io_formats.write_map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录列表，每条为字典 `{'chrom', 'id', 'cm', 'pos'}` 或列表 `[chrom, id, cm, bp]` / `[chrom, id, bp]`（此时 `cm` 取 0），也可以是待规范化的现有 `.map` 文本/路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_map


envelope = write_map([["1", "rs1", 100], {"chrom": "2", "id": "rs2", "cm": 0.5, "pos": 10}], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 2,
    "bytes": 25,
    "sha256": "a0fc688e87258f9344dc06823654e51300a21c589a3421303fcada7d66454b49",
    "text": "1\trs1\t0\t100\n2\trs2\t0.5\t10\n"
  }
}
```

## 125) write_newick

- **作用：** 把树序列化为 Newick 文件，或规范化已有的 Newick 内容。

- **处理方法：** 按给定描述（解析 Newick 文本、扁平节点表或嵌套字典）构造每棵树，并按输入顺序写成以 `;` 结尾的规范化 Newick 行。扁平节点表中每个节点需要 `id` 和 `parent`（根节点的 `parent` 为 `null`），`name` 和 `length` 可选。

- **API：** `dnakit.io_formats.write_newick(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：Newick 文本/路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）、`{'newick': ...}` 字典、扁平的 `{'nodes': [{'id', 'parent', 'name', 'length'}, ...]}` 字典、嵌套的 `{'name', 'length', 'children'}` 字典，或这些对象的列表。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_newick


envelope = write_newick({"name": "root", "children": [{"name": "A", "length": 1}, {"name": "B", "length": 2}]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 1,
    "bytes": 15,
    "sha256": "d189445c127e5c47ea6623d600519888cb3f22fa38218bca1f679a1f5e099c38",
    "text": "(A:1,B:2)root;\n"
  }
}
```

## 126) write_nexus

- **作用：** 把序列序列化为 NEXUS 文件，或规范化已有的 NEXUS 内容。

- **处理方法：** 写出带一个 DATA 块（`DIMENSIONS NTAX/NCHAR`、`FORMAT DATATYPE=... MISSING=? GAP=-`、矩阵；需要时给名称加引号）的 `#NEXUS` 文件，输入含树时再加 TREES 块。`layout='interleaved'` 时矩阵按每块 `block_size` 个位点交错写出。FASTA、PHYLIP 和 NEXUS 输入会先被解析；普通映射和记录列表使用 `datatype`（默认 `dna`）。长度不等的序列抛出 `InputFormatError`。

- **API：** `dnakit.io_formats.write_nexus(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA/PHYLIP/NEXUS 文本或路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）、`{name: sequence}` 映射，或 `{'name', 'sequence'}` 记录列表。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `layout`：'sequential'（默认）或 'interleaved'：PHYLIP/NEXUS 输出的排布
    - `block_size`：整数，默认 60：交错（interleaved）输出每块的位点数
    - `strict_names`：布尔，默认 False：写出严格的 10 字符 PHYLIP 名称（更长时报错）
    - `line_width`：整数，默认 60：FASTA 行宽（0 表示每条序列一行）
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：用于未声明类型的记录

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_nexus


envelope = write_nexus({"A": "ACGTACGT", "B": "ACGTACGA"}, options={"layout": "interleaved", "block_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 2,
    "bytes": 178,
    "sha256": "90285f9c9dce3e7d66a73c016e18a2e474be23f9498a6210624311c0f37385f5",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=2 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;\n    MATRIX\n    A  ACGT\n    B  ACGT\n\n    A  ACGT\n    B  ACGA\n    ;\nEND;\n"
  }
}
```

## 127) write_paf

- **作用：** 把比对记录序列化为 PAF 文件，或规范化已有的 PAF 内容。

- **处理方法：** 每条记录写成一个制表符分隔的行：12 个必需字段，后接其标签 `TAG:TYPE:VALUE`。输入映射中以普通值给出的标签（`{'tp': 'P'}`）以类型 `Z` 写出。已有的 PAF 文本或路径会被解析后重新序列化。会检查每条记录的坐标（查询和目标上均需 `0 <= start <= end <= length`），不合格的记录抛出 `InputFormatError`。

- **API：** `dnakit.io_formats.write_paf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：含 12 个必需 PAF 字段（`query_name`、`query_length`、`query_start`、`query_end`、`strand`、`target_name`、`target_length`、`target_start`、`target_end`、`matches`、`block_length`、`mapq`）及可选 `tags` 的字典列表，也可以是待规范化的现有 PAF 文本/路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_paf

RECORD = {
    "query_name": "q", "query_length": 10, "query_start": 0, "query_end": 10, "strand": "+",
    "target_name": "t", "target_length": 20, "target_start": 0, "target_end": 10,
    "matches": 9, "block_length": 10, "mapq": 5, "tags": {"tp": "P"},
}

envelope = write_paf([RECORD], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "bytes": 36,
    "sha256": "65ccfbff32154adca8949b2e2d4e75fb2fd8c528ca3f9274613f4f05b6df05a9",
    "text": "q\t10\t0\t10\t+\tt\t20\t0\t10\t9\t10\t5\ttp:Z:P\n"
  }
}
```

## 128) write_ped

- **作用：** 把个体记录序列化为 PLINK `.ped` 文件，或规范化已有的 `.ped` 内容。

- **处理方法：** 每个个体写成 `fid iid pat mat sex pheno`，其后是每个变异的两个等位基因，以制表符分隔、每个个体一行、保持输入顺序。已有的 `.ped` 文本或路径先被解析再重新序列化（以空白分隔的输入变为制表符分隔）。基因型矩阵大小受 `max_cells` 限制。

- **API：** `dnakit.io_formats.write_ped(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：个体列表，每个为字典 `{'fid', 'iid', 'pat', 'mat', 'sex', 'pheno', 'genotypes': [[a1, a2], ...]}`（`fid` 和 `iid` 必填；`pat`、`mat`、`sex` 默认为 0，`pheno` 默认为 -9，`genotypes` 默认为空），也可以是待规范化的现有 `.ped` 文本/路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `max_cells`：整数，默认 20000000：基因型单元格（个体 x 变异）数上限

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`variant_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_ped


envelope = write_ped([{"fid": "F1", "iid": "1", "sex": 1, "pheno": 2, "genotypes": [["A", "C"], ["G", "G"]]}], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 1,
    "variant_count": 2,
    "bytes": 21,
    "sha256": "df689e3f321aa8c782f823e51956fa08c4a41c841436038a2bde33b30e8f4944",
    "text": "F1\t1\t0\t0\t1\t2\tA\tC\tG\tG\n"
  }
}
```

## 129) write_phylip

- **作用：** 把序列序列化为 PHYLIP 比对文件，或规范化已有的 PHYLIP 内容。

- **处理方法：** 写出首行 `ntax nchar`，每个分类单元一行：默认是宽松 PHYLIP（名称补足到 10 个字符，再加一个空格），`strict_names=True` 时是严格 PHYLIP（名称字段恰好 10 个字符；名称更长则抛出 `InputFormatError`）。含空白的名称无法写出。`layout='interleaved'` 时按每块 `block_size` 个位点写出。序列必须等长（否则抛出 `InputFormatError`）；序列会按 `datatype` 的字符集检查（`dna`（默认）和 `rna`：IUPAC 核苷酸编码、`-`、`?`、`.`、`*` 以及文档自己的缺口/缺失字符，不区分大小写；`protein`：氨基酸字母加 B、Z、X、J、U、O、`*`、`-`、`?`、`.`；`any`：不检查），字符集之外的字符会抛出 `InputFormatError`（代码 `IO_FORMATS_BAD_ALPHABET`），并列出这些字符及其所在的序列。

- **API：** `dnakit.io_formats.write_phylip(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA/PHYLIP/NEXUS 文本或路径（文件路径（`str` 或 `Path`），或内容本身（`str`/`bytes`）；字符串若是已存在的文件则按路径读取；gzip 数据按魔数识别并解压，受 `max_bytes` 限制）、`{name: sequence}` 映射，或 `{'name', 'sequence'}` 记录列表。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `datatype`：'dna'（默认）、'rna'、'protein' 或 'any'：序列所要检查的字符集（核苷酸为 IUPAC 编码、缺口、?、N；蛋白为 X；'any' 不检查）
    - `layout`：'sequential'（默认）或 'interleaved'：PHYLIP/NEXUS 输出的排布
    - `block_size`：整数，默认 60：交错（interleaved）输出每块的位点数
    - `strict_names`：布尔，默认 False：写出严格的 10 字符 PHYLIP 名称（更长时报错）
    - `line_width`：整数，默认 60：FASTA 行宽（0 表示每条序列一行）
    - `input_layout`：'auto'（默认）、'sequential' 或 'interleaved'：输入为 PHYLIP 文本时的排布
    - `input_strict_names`：布尔，默认 False：输入为 PHYLIP 文本时名称是否占每行前 10 个字符

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `text`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `text`。内容格式错误或字符超出字符集时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_phylip


envelope = write_phylip(">a\nACGT\n>b\nACGA\n", options={"strict_names": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 2,
    "bytes": 35,
    "sha256": "162355e387d5b8a7ccc32c599fa049f4b8bdd217ad19b419ba040139e4875603",
    "text": " 2 4\na         ACGT\nb         ACGA\n"
  }
}
```

## 130) write_twobit

- **作用：** 把序列序列化为 UCSC 2bit 文件，或重新编码已有的 2bit 文件。

- **处理方法：** 把序列打包成版本 0 的 2bit 文件：A、C、G、T 各占 2 位，其他字母变成 N 块，小写字母变成屏蔽块（所以 R、Y 等 IUPAC 码按 N 存储）。`byte_order` 选择小端（默认）或大端。序列名必须唯一且非空。

- **API：** `dnakit.io_formats.write_twobit(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本/路径、`{name: sequence}` 映射、`{'name', 'sequence'}` 记录列表，或要重新编码的现有 2bit 文件（路径或 bytes）。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `max_bytes`：整数，默认 268435456：拒绝超过该字节数的文本/文件输入
    - `output`：输出文件路径（可选；不给时内容随结果直接返回）
    - `overwrite`：布尔，默认 False：允许覆盖已存在的输出文件
    - `byte_order`：'little'（默认）或 'big'：所写 2bit 文件的字节序

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`、`sequence_count`、序列化输出的 `bytes` 与 `sha256`，以及内容 `data_base64`；给出 `options['output']` 时写入文件（已存在的文件需要 `overwrite=True`），并返回 `output_path` 而不再返回 `data_base64`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_twobit


envelope = write_twobit(">s1\nACGTNNacgt\n>s2\nGGCC\n", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "bytes": 82,
    "sha256": "a75737f8defa5dad84385216e2f7b9091991fce11a1d9b3ebb45b4d15849423f",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMR4AAAACczJBAAAACgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcAEAAAAAAAAAAAAAAAAAAAA9Q=="
  }
}
```

## 131) write_zarr

- **作用：** 用具名数组创建 Zarr 目录存储。

- **处理方法：** 在 `output` 处创建 zarr v3（默认）或 v2（`zarr_format=2`）目录存储，每个名称一个数组；`chunks` 设置分块形状（默认由 zarr 选择），`attrs` 设置根分组的属性。路径已存在时需要 `overwrite=True`；`output` 必填。其他 dtype（例如字符串）的数组会抛出 `InputFormatError`。

- **API：** `dnakit.io_formats.write_zarr(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{数组名: 类数组}` 映射、`{'name', 'data'}` 记录列表，或单个类数组对象（存为 `options['array']`，默认 `data`）；支持的 dtype 为 bool、int、uint 和 float。`options` 是可选映射，未知键抛出 `ConfigurationError`。`verbose`（默认 `True`）会打印进度消息；示例传入 `verbose=False`。选项键（默认值如所述）：

    - `output`：新 zarr 存储的目录路径（必填）
    - `overwrite`：布尔，默认 False：替换 output 处已存在的 zarr 存储
    - `array`：输入为单个类数组对象时的数组名（默认 'data'）
    - `chunks`：整数或整数列表：分块形状（默认由 zarr 选择）
    - `attrs`：存为根分组属性的映射
    - `zarr_format`：2 或 3（默认 3）：zarr 格式版本

- **输出：** 调用返回 DNAKit 结果封装，包含 `status`（`"ok"`）、`module`（`"io_formats"`）、`operation`、`provenance`（`implementation: "native"`，`backend: "local"`）、`input_summary`、`warnings` 和 `result`（示例只打印 `status`、`warnings` 和 `result`）。`result`：`format`、`record_count`（写入的数组个数）、`zarr_format` 和 `output_path`（存储写入 `options['output']`，不返回内容）。如果无法导入 `zarr`（>= 3）或 `numpy`，函数返回 `status: "not_implemented"`，带有 `required_dependency: "pip install 'zarr>=3' numpy"` 和 `message`（`zarr and numpy are required: <导入错误>`，或 `zarr <版本> is installed but this adapter needs zarr >= 3.`），而没有 `result`。内容格式错误时抛出 `InputFormatError`；选项无效、含未知选项或输入超过大小限制时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.io_formats import write_zarr


envelope = write_zarr({"pos": [100, 200, 300]}, options={"output": "written.zarr", "overwrite": True, "chunks": 2, "attrs": {"source": "demo"}}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **示例输出：**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 1,
    "zarr_format": 3,
    "output_path": "written.zarr"
  }
}
```

