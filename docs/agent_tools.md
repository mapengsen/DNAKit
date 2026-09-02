# Agent 与 MCP 工具

DNAKit 可以通过 MCP 被 Codex、Claude、Cursor、VS Code 等支持 MCP 的 Agent 调用。Agent
不直接加载全部 Python 函数，而是先搜索工具目录、查看参数 Schema，再执行选中的公开 API。

## 1) 安装

安装 DNAKit 和 Agent 可选依赖：

```bash
python -m pip install "dnakit[agent]"
```

在源码开发环境中安装：

```bash
python -m pip install -e ".[agent]"
```

深度学习预测、数据库查询或 ToolUniverse 后端仍需安装各自的可选依赖。例如：

```bash
python -m pip install "dnakit[agent,external-tools,neural,neural-enformer]"
```

## 2) 启动服务

```bash
dnakit-mcp
```

该命令默认使用本地 `stdio` MCP 传输。它会等待 Agent 客户端连接，因此直接运行时不会输出普通终端内容。

## 3) Agent 配置

在 Agent 的 MCP 配置中加入：

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/.venv/bin/dnakit-mcp"
    }
  }
}
```

Windows 虚拟环境通常使用：

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "D:\\project\\.venv\\Scripts\\dnakit-mcp.exe"
    }
  }
}
```

必须使用 Agent 实际能够访问的绝对路径。不同客户端的配置文件位置不同，但启动命令相同。

## 4) Agent 可见入口

MCP 服务只暴露六个紧凑入口，避免一次发送数百个工具 Schema：

1. `dnakit_catalog_stats`：查看工具总数；
2. `list_dnakit_categories`：查看功能分类；
3. `list_dnakit_tools`：分页列出工具；
4. `search_dnakit_tools`：按功能搜索；
5. `describe_dnakit_tool`：查看完整输入 Schema 和安全属性；
6. `call_dnakit_tool`：执行具体 DNAKit 功能。

当前源码会自动发现 326 个公开函数，其中 323 个可以直接通过 JSON 参数调用。两个必须接收 Python
回调的底层批处理函数，以及任意命令执行底层函数保持 Python-only。普通批量任务可由 Agent 多次调用具体
DNAKit 功能完成。

## 5) 调用示例

用户可以直接对 Agent 说：

```text
使用 DNAKit 计算 ACGTACGT 的单链分子量。
```

Agent 的典型调用过程是：

```text
search_dnakit_tools("molecular weight")
describe_dnakit_tool("dnakit.thermodynamics.molecular_weight")
call_dnakit_tool(
    tool_name="dnakit.thermodynamics.molecular_weight",
    arguments={"sequence": "ACGTACGT", "strand": "single"}
)
```

返回值是有界、可序列化的 JSON，并保留算法版本、适用边界和 provenance 等 DNAKit 审计信息。

## 6) 输入与安全规则

- `DNASequence` 可以直接使用序列字符串；`DNARecord` 使用包含 `sequence` 和 `id` 的对象。
- 配置类使用普通 JSON 对象，Agent 可通过工具 Schema 查看字段和默认值。
- 27 个本地任务 checkpoint 可通过 `dnakit.predictions.predict_enformer_benchmark` 调用；它属于模型工具，Agent 执行前仍需显式确认并提供 checkpoint 路径。
- Python callback、已打开文件、活动 backend 实例等对象不能通过 JSON 传入。
- 文件写入、模型下载和外部程序功能默认拒绝执行；确认目标及许可后，需要显式设置
  `allow_side_effects=true`。
- 数据库查询和模型功能仍受网络、API、模型权重、GPU、许可证及可选依赖限制。
- 单次参数最大 2 MB，单次返回默认最大 5 MB；大型结果应使用 DNAKit 文件输出功能。
