# 外部基因组分析与 MCP：参数与补充说明

作用、计算方法与示例见[功能教程](../../api/features/39_external_genomics.md)。论文、作者和许可证证据见[来源 FAQ](../genomics_sources.md)。

## 目录状态与检索约定 {#discovery}

`search_projects` 与 `list_analyses` 接受可选的 `query` 和 `domain`，多个关键词需要同时匹配。`describe_analysis(analysis)` 接受已登记的分析 ID，返回输入键名、参数类型与范围、安装说明、输出文件和来源。

`reference_only` 表示发现了来源，但尚无该项目的本地接口。`available` 还需区分原生函数、外部适配器、MCP 或工作流文档；收录不代表上游全部能力均已实现或完成科学验证。初始缺口字段描述扩展前的情况，`remaining_gaps` 描述当前限制。

## 外部分析的参数与运行约定 {#external-analysis}

使用前按具体分析单独准备软件、参考数据库、模型权重和凭据，统一参考版本、样本标识和等位基因方向。输入前缀必须对应完整文件组，例如 PLINK 的 BED/BIM/FAM。

| 参数 | 使用约定 |
| --- | --- |
| `analysis` | 必填；`list_analyses` 中的分析 ID。 |
| `inputs` | 必填；字典键必须与 `describe_analysis` 的输入定义完全一致，值为文件路径或前缀。 |
| `output_directory` | 必填；尚不存在的输出目录。 |
| `parameters` | 可选；仅接受登记的参数，默认值和范围见 `describe_analysis`；不接受任意命令行参数。 |
| `threads` | 默认 `1`，整数范围 `1–1024`。 |
| `executable_paths` | 可选；按后端程序名指定路径，如 `{"plink2": "/path/to/plink2"}`。 |
| `timeout_seconds` | 仅执行接口；默认 `3600.0`，范围 `(0, 86400]` 秒。 |
| `max_artifact_bytes` | 仅执行接口；默认 `1000000000`，限制单个产物文件的字节数。 |
| `verbose` | 默认 `True`；控制运行提示和进度显示。 |

`plan_analysis` 校验输入并生成命令，不写文件；计划会标记尚未找到的可执行程序。`run_analysis` 要求后端可用，在新输出目录保存命令、输入、参数、来源、日志与 `dnakit-run.json`。执行失败仍写失败记录；`completed` 表示程序退出与文件契约通过，科学结论仍需 QC。

PLINK 四项已用真实 PLINK 2 和固定随机种子的人工数据运行。各补充适配器的 `validation` 区分接口审查、人工数据实际运行和未运行，见[335 条审核](../user_catalog.md)。FastANI 无匹配、Barrnap 无基因、Cas-OFFinder 无命中等有效空结果会明确允许；其他必须产物为空会判失败。`modkit-pileup` 的参考 FASTA 需预先建立 `.fai`。

## MCP 安装与连接参数 {#mcp-connection}

在源码环境安装 MCP 客户端支持：

```bash
python -m pip install -e ".[agent]"
```

各第三方服务需按其上游说明单独安装、配置。`list_mcp_servers` 返回启动方式和安装信息，不读取凭据或启动服务；功能教程中的服务目录保留逐项配置。客户端不会自动安装服务，目录收录也不代表本机服务可用或完成科学验证。

| 参数 | 使用约定 |
| --- | --- |
| `server_id` | 必填；服务目录中的 ID。 |
| `tool_name` | 单次工具调用必填；必须取自实时工具清单。 |
| `arguments` | 单次工具调用必填；符合工具 `inputSchema` 的 JSON 对象。 |
| `transport` | 默认 `"stdio"`；`"http"` 仅用于目录登记的作者公开端点。 |
| `executable_path` | 可选；单独安装的服务可执行程序或解释器路径。 |
| `workdir` | 可选；运行目录，源码型服务通常需要指向项目目录。 |
| `environment_variables` | 可选；将服务变量名映射到当前环境中已有的变量名，值不是凭据原文。 |
| `timeout_seconds` | 默认 `60.0` 秒；控制请求或整个工作流会话的超时。 |
| `max_response_bytes` | 默认 `4194304`；SDK 解码后的响应大小限制。 |
| `max_tools` | 默认 `5000`；允许的工具清单大小上限。 |
| `verbose` | 默认 `True`；控制运行提示和进度显示。 |

`list_mcp_tools` 启动 stdio 服务或连接登记的 HTTP 服务，获取真实 schema。`call_mcp_tool` 校验参数后才执行工具；远端 `isError` 会保留为 `remote_error`，已知转发凭据会从结果中脱敏。每次独立调用建立并关闭一个会话，服务内存状态不会跨调用保留。响应大小在 SDK 解码后检查，不是网络流量限额。

本地真实 stdio 协议、参数验证、会话状态和错误中止已验证；补充审核还对 6 个作者 HTTP 服务完成真实工具发现。服务科学分析、账户和模型环境的验证边界见[补充清单 FAQ](../user_catalog.md)。

## 同一会话的调用列表 {#mcp-workflow}

`call_mcp_workflow` 接受必填的 `server_id` 与 `calls`，其余连接参数与上节一致。`calls` 包含 `1–100` 项，每项只能包含 `tool_name` 和 `arguments`。工具参数显式传入，不执行脚本，也不自动将上一步输出代入下一步。

每步重新读取 schema 并校验参数；任一步远端报错即停止，返回已执行步骤（包括失败步骤）。配置或连接错误会抛异常。此前已执行的上传或计算不会回滚。`timeout_seconds` 限制整个会话，响应字节限制作用于累计解码结果；按实际工具说明判断联网、文件与计算资源使用。

## DNAKit MCP 入口 {#dnakit-mcp}

安装客户端支持后启动统一入口：

```bash
dnakit-mcp
```

先用 `search_dnakit_tools` 搜索 `genomics` 或 `integrations`，再用 `describe_dnakit_tool` 查看完整签名与 schema。外部分析和远端 MCP 调用标记为 `external_process`；统一执行入口使用已有的 `allow_side_effects=true` 调用机制。

## 完整 API 签名

::: dnakit.integrations
    options:
      show_source: false
