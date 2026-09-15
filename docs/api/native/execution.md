# 接口发现、工作流与报告

本页收录 **12 个新增公开函数**。 [完整索引](index.md) · [调用约定与完整示例](usage.md)

每个函数均列出当前 Python 签名、参数、返回值和已登记的执行要求。提供方说明及字段描述保留源码原文。`None` 默认值在生成的数据库包装中通常表示不向服务传该参数；实际必填项与组合约束以输入 Schema 为准。坐标系和基因组版本按各接口说明使用。

scripts、resource_path、describe、search、list_operations 用于本地发现或资源读取。run、run_support 会执行脚本，可能联网或写文件；解释器、科学软件和参考数据需要按具体任务准备。原有严格 YAML 工作流 run_workflow 与这里的脚本执行入口是不同函数。

## 函数目录

- [`dnakit.databases.call`](#dnakit.databases.call)
- [`dnakit.databases.describe`](#dnakit.databases.describe)
- [`dnakit.databases.execute`](#dnakit.databases.execute)
- [`dnakit.databases.execute_tool`](#dnakit.databases.execute_tool)
- [`dnakit.databases.list_operations`](#dnakit.databases.list_operations)
- [`dnakit.databases.read_executed_notebook`](#dnakit.databases.read_executed_notebook)
- [`dnakit.databases.search`](#dnakit.databases.search)
- [`dnakit.workflows.resource_path`](#dnakit.workflows.resource_path)
- [`dnakit.workflows.run`](#dnakit.workflows.run)
- [`dnakit.workflows.run_support`](#dnakit.workflows.run_support)
- [`dnakit.workflows.scripts`](#dnakit.workflows.scripts)
- [`dnakit.workflows.write_wes_report`](#dnakit.workflows.write_wes_report)

## `dnakit.databases.call` {#dnakit.databases.call}

根据操作名和参数字典调用数据库接口。

```text
dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | str | 必填 | 本地操作名或登记的上游别名。 |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `max_result_bytes` | int | 20000000 | 序列化结果的字节上限，超限抛出异常。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "补充说明（源码原文）"

    Call a bundled DNAKit database adapter without installing ToolUniverse.

    An operation may use its native snake_case name or recorded upstream name.
    Input JSON is checked against the pinned service schema; returned provider
    envelopes preserve the service's fields, pagination and evidence details.
    Failure raises an exception, including provider error envelopes. Models,
    reference data, API credentials and specialist libraries remain external.

## `dnakit.databases.describe` {#dnakit.databases.describe}

查看单个接口的实际输入 Schema、上游别名和实现信息。

```text
dnakit.databases.describe(operation, *, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | str | 必填 | 本地操作名或登记的上游别名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

字典；input_schema 为实际输入约束，config 保存提供方配置，upstream_name 为原接口名。

??? info "补充说明（源码原文）"

    Return the function's actual upstream API schema and implementation provenance.

## `dnakit.databases.execute` {#dnakit.databases.execute}

使用操作名或结构化请求分发一次调用。

```text
dnakit.databases.execute(request, arguments=None, *, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `request` | str &#124; Mapping&#91;str, object&#93; | 必填 | 操作名，或包含 name、arguments 的字典。 |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

返回类型为 object，具体结构由所调用操作决定；不保证统一的 data 层。提供方错误会抛出异常。

??? info "补充说明（源码原文）"

    Execute an operation name plus arguments, or a JSON name/arguments request.

## `dnakit.databases.execute_tool` {#dnakit.databases.execute_tool}

根据工具名与参数调用本地 DNAKit 数据库操作；禁止递归调用自身。

```text
dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `tool_name` | string | 必填（输入 Schema） | Name of the tool to execute (e.g., &#x27;example_tool_name&#x27;) |
| `arguments` | object &#124; None | None | Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {&quot;param1&quot;: &quot;value1&quot;, &quot;param2&quot;: 5}. JSON string example: &quot;{\&quot;param1\&quot;: \&quot;value1\&quot;, \&quot;param2\&quot;: 5}&quot;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;object&quot;, &quot;additionalProperties&quot;: true}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "tool_name": {
          "type": "string",
          "description": "Name of the tool to execute (e.g., 'example_tool_name')"
        },
        "arguments": {
          "oneOf": [
            {
              "type": "object",
              "additionalProperties": true
            },
            {
              "type": "string"
            }
          ],
          "description": "Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {\"param1\": \"value1\", \"param2\": 5}. JSON string example: \"{\\\"param1\\\": \\\"value1\\\", \\\"param2\\\": 5}\"."
        }
      },
      "required": [
        "tool_name"
      ]
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "type": "object",
      "description": "Tool execution result as a structured object (dict). The actual content depends on the executed tool. May contain error information if execution failed. Returns structured data instead of JSON string for programmatic use.",
      "additionalProperties": true
    }
    ```

**执行与来源：** `network`；上游接口名 `execute_tool`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("execute_tool", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.list_operations` {#dnakit.databases.list_operations}

离线列出全部数据库操作名与上游别名。

```text
dnakit.databases.list_operations(*, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

列表，每项包含 name、upstream_name。

??? info "补充说明（源码原文）"

    List native operation names and recorded upstream aliases without loading providers.

## `dnakit.databases.read_executed_notebook` {#dnakit.databases.read_executed_notebook}

Read a pre-executed Jupyter notebook (.ipynb, often named '*_executed.ipynb') and return its cells with printed outputs as structured JSON. Use this BEFORE reimplementing an analysis when the data folder ships with an executed notebook — the notebook already ran the authoritative pipeline with the exact package versions, filters, and thresholds that match the reference answers. Reimplementing with pydeseq2/scanpy/etc. produces different numbers. Accepts either a folder (auto-discovers the notebook) or a direct path. Optionally filters cells by a comma-separated search term list or a regex so the agent sees only the cells relevant to its question.

```text
dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `data_folder` | string &#124; null | None | Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path. |
| `notebook_path` | string &#124; null | None | Direct path to a specific .ipynb file. Provide either this OR data_folder. |
| `search` | string | None | Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in &#x27;matching_cells&#x27;. Case-insensitive. Example: &#x27;upregulated,log2FoldChange&#x27; or &#x27;padj.*0\.05&#x27; |
| `regex` | boolean | None | Treat &#x27;search&#x27; as a Python regex instead of a comma-separated term list (default: false) |
| `max_output_chars` | integer | None | Truncate each cell source/output to this many characters (default 400) |
| `include_source` | boolean | None | Include cell source code in the response (default true). Set false to halve response size when you only need outputs. |
| `verbose` | bool | True | 运行提示开关。 |

**返回值**

返回提供方结果，经 DNAKit 参数处理、错误检查和大小限制；下方为提供方登记的结果 Schema，属于结构参考，不是本次实时请求输出。不同接口的 data 层与分页字段可能不同。

??? info "完整输入 Schema（含别名与条件约束）"

    ```json
    {
      "type": "object",
      "properties": {
        "data_folder": {
          "type": [
            "string",
            "null"
          ],
          "description": "Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path."
        },
        "notebook_path": {
          "type": [
            "string",
            "null"
          ],
          "description": "Direct path to a specific .ipynb file. Provide either this OR data_folder."
        },
        "search": {
          "type": "string",
          "description": "Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in 'matching_cells'. Case-insensitive. Example: 'upregulated,log2FoldChange' or 'padj.*0\\.05'"
        },
        "regex": {
          "type": "boolean",
          "description": "Treat 'search' as a Python regex instead of a comma-separated term list (default: false)"
        },
        "max_output_chars": {
          "type": "integer",
          "description": "Truncate each cell source/output to this many characters (default 400)"
        },
        "include_source": {
          "type": "boolean",
          "description": "Include cell source code in the response (default true). Set false to halve response size when you only need outputs."
        }
      },
      "required": []
    }
    ```

??? info "提供方登记的结果 Schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "notebook_path": {
              "type": "string"
            },
            "total_cells": {
              "type": "integer"
            },
            "cells": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "idx": {
                    "type": "integer"
                  },
                  "type": {
                    "type": "string"
                  },
                  "source": {
                    "type": "string"
                  },
                  "output": {
                    "type": "string"
                  }
                }
              }
            },
            "search": {
              "type": "string"
            },
            "matching_cells": {
              "type": "array"
            },
            "n_matches": {
              "type": "integer"
            }
          },
          "required": [
            "notebook_path",
            "total_cells",
            "cells"
          ]
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**执行与来源：** `network`；上游接口名 `read_executed_notebook`。

**查看参数示例（离线，不调用提供方）：**

```python
import dnakit
info = dnakit.databases.describe("read_executed_notebook", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.search` {#dnakit.databases.search}

按名称、任务关键词或描述离线搜索数据库接口。

```text
dnakit.databases.search(query, *, limit=10, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | str | 必填 | 本地目录搜索词，匹配名称与描述；使用英文功能词或函数名。 |
| `limit` | int | 10 | 最多返回的匹配条数；数据库目录搜索允许 1–100。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

列表，每项包含 name、description、upstream_name。

??? info "补充说明（源码原文）"

    Find bundled scientific API operations locally by task, name or description.

## `dnakit.workflows.resource_path` {#dnakit.workflows.resource_path}

取得已校验的随包资源路径。

```text
dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `resource` | str | 必填 | 技能包内相对路径，读取前校验哈希。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

pathlib.Path；指向通过哈希校验的本地资源。

??? info "补充说明（源码原文）"

    Locate a checksum-verified bundled file, including resources of dependent skills.

## `dnakit.workflows.run` {#dnakit.workflows.run}

传递原脚本参数，执行随包工作流并返回真实退出码和输出。

```text
dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `script` | str | 必填 | 技能包内 .py、.R 或 .sh 脚本相对路径。 |
| `arguments` | list&#91;str&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `workdir` | str &#124; Path &#124; None | None | 子进程工作目录；None 使用当前目录。 |
| `timeout_seconds` | float | 300.0 | 子进程执行超时秒数。 |
| `max_output_bytes` | int | 1000000 | 最多捕获的子进程输出字节数。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

字典包含 skill、script、returncode、output、upstream_url、upstream_revision、license；returncode 非零表示脚本失败。

??? info "补充说明（源码原文）"

    Run a bundled workflow with its original CLI options and true scientific backends.

    The caller chooses an explicit script and working directory. Python uses the
    current environment; R, shell tools, trained models and data retain their own
    requirements. This executes code and may write analysis outputs or access APIs.
    Nonzero exit status is returned, not interpreted as scientific success.

## `dnakit.workflows.run_support` {#dnakit.workflows.run_support}

执行组合流程依赖的随包支持脚本。

```text
dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `name` | str | 必填 | 支持流程名称，允许值见本条说明。 |
| `arguments` | list&#91;str&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `workdir` | str &#124; Path &#124; None | None | 子进程工作目录；None 使用当前目录。 |
| `timeout_seconds` | float | 300 | 子进程执行超时秒数。 |
| `max_output_bytes` | int | 1000000 | 最多捕获的子进程输出字节数。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

字典包含 name、returncode、output；须检查退出码。

??? info "补充说明（源码原文）"

    Run a licensed support stage of a composite skill with its documented CLI options.

    Names: clinical-variant-reporter, clinical-trial-finder, omics-target-evidence-mapper,
    gwas-lookup, fine-mapping.
    Source implementations are bundled in DNAKit; scientific dependencies remain external.

## `dnakit.workflows.scripts` {#dnakit.workflows.scripts}

列出技能包中可运行的 Python、R 和 Shell 脚本。

```text
dnakit.workflows.scripts(skill_id, *, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

技能包中的脚本相对路径列表；可能为空。

??? info "补充说明（源码原文）"

    List bundled Python, R and shell resources; consult their usage before execution.

## `dnakit.workflows.write_wes_report` {#dnakit.workflows.write_wes_report}

把已有 Markdown 报告渲染为 PDF，并记录来源哈希。

```text
dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)
```

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `input_report` | str &#124; Path | 必填 | 已有 UTF-8 Markdown 报告路径，最多 5 MB。 |
| `output_dir` | str &#124; Path | 必填 | 报告输出目录，生成 PDF 和来源清单。 |
| `language` | str | &#x27;en&#x27; | 报告标签语言 en 或 es；不翻译输入正文。 |
| `title` | str &#124; None | None | 自定义报告标题；None 使用语言对应默认标题。 |
| `logo_left` | str &#124; Path &#124; None | None | 左侧 Logo 图片路径。 |
| `logo_right` | str &#124; Path &#124; None | None | 右侧 Logo 图片路径。 |
| `font_path` | str &#124; Path &#124; None | None | 可选 TTF 字体路径。 |
| `overwrite` | bool | False | 是否覆盖已存在的报告文件。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |

**返回值**

返回报告路径和来源记录的字典；输出 wes_report.pdf 与 wes_report.manifest.json。需要 dnakit[reports]。

??? info "补充说明（源码原文）"

    Render supplied Markdown sections and tables as a PDF, with a source hash.

    language selects English or Spanish report labels; input content is preserved
    without translation or new clinical interpretations. Optional logos and a TTF
    font are supported. Requires the reportlab rendering backend (dnakit[reports]).
    This implementation was written for DNAKit; no proprietary report code is used.

