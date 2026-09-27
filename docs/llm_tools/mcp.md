# MCP

## 推荐：一句话安装

先安装并登录一个支持 MCP、能读取文档和执行本地命令的 Agent，例如 Codex、Claude Code 或
Gemini CLI。将下面这句话复制到 **Agent 对话框**：

```text
请阅读 https://mapengsen.github.io/DNAKit/setup.md，帮我安装并配置 DNAKit MCP，然后验证工具调用。
```

Agent 会完成环境检查、安装、客户端配置和工具验证。需要登录或重启客户端时，按提示操作。
这个入门示例使用本地计算和默认关键词搜索，无需额外 API key。

## 工作方式

DNAKit 可以通过 MCP 被 Codex、Claude、Cursor、VS Code 等支持 MCP 的 Agent 调用。Agent
不直接加载全部 Python 函数，而是先搜索工具目录、查看参数 Schema，再执行选中的公开 API。

本页使用本地 `stdio` 模式，MCP Server 是运行在用户电脑上的程序。用户安装 DNAKit、在 Agent 中配置启动命令后，由 Agent
启动服务并调用本地函数；无需由 DNAKit 作者提供计算服务器。模型负责决策，实际计算、文件读写和
模型权重加载由服务所在环境完成，工具也可以请求外部 API。

使用 DeepSeek 等模型 API 时，需要选择同时支持该 API 和 MCP 的 Agent，或自行实现 MCP 客户端。
模型服务能否切换取决于 Agent。DNAKit 的内置 `run_agent()` 则直接调用本地工具目录，不经过 MCP，
用法见 [Agent](agent.md)。

工具生成、说明优化或兼容修复、本地结果缓存和知识库检索（RAG）由内置 `run_agent()` 的[对应配置](agent.md#agent-actions)开启，不是 MCP 新增动作，也不会自动为其他 MCP 客户端开启这些能力。已验证并登记到 DNAKit 的生成工具由现有搜索、描述和调用入口使用；独立运行的 MCP 服务需重启后加载更新。

下面的完整流程是：**安装 → 配置客户端 → 检查连接 → 输入自然语言任务**。

## 1) 安装

准备 Python 3.10 及以上版本，并安装、登录所选 Agent。Codex 用户可以先参考
[OpenAI 官方 CLI 安装说明](https://learn.chatgpt.com/docs/codex/cli)。

在 Linux 或 macOS 的专用环境中安装 DNAKit 和 MCP 依赖：

```bash
python3 -m venv "$HOME/.venvs/dnakit-mcp"
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install --upgrade pip
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install "dnakit[agent]"
```

`[agent]` 安装 MCP 接入依赖。训练、绘图、深度模型等任务还需要在**同一环境**安装对应的可选依赖、
外部程序或模型权重，见[安装说明](../installation.md)。

测试本地源码时，在仓库根目录将最后一条安装命令的 `"dnakit[agent]"` 替换为 `-e ".[agent]"`。
安装环境必须能被客户端访问。

<span id="codex-setup"></span>

## 2) Codex 接入与验证

### 2.1 注册本地服务

在已经安装 Codex 的同一个环境中执行：

```bash
codex mcp add dnakit -- "$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
codex mcp list
```

`codex mcp list` 中出现 `dnakit` 表示配置已登记；实际连接和执行还需按下文检查。

### 2.2 设置等待时间

打开 Codex 用户配置 `~/.codex/config.toml`。在刚生成的 `[mcp_servers.dnakit]` 区块中，
紧接 `command` 行加入下面两行；若已存在同名设置，修改其值即可：

```toml
startup_timeout_sec = 60
tool_timeout_sec = 3600
```

启动等待设为 60 秒，适用于导入工具目录较慢的环境。单次工具执行上限示例为一小时，按任务调整；
它不保证任务一定在时限内完成。配置格式和参数见
[OpenAI 官方 MCP 文档](https://learn.chatgpt.com/docs/extend/mcp?surface=cli)。

### 2.3 启动并完成第一次调用

在自己的分析工作目录运行：

```bash
codex
```

首次使用按提示登录。已经打开的会话在修改配置后重新启动。在 Codex 交互界面输入：

```text
/mcp
```

确认 `dnakit` 已连接、能看到本页列出的八个入口，然后发送：

> 请使用 DNAKit MCP 计算 AACG 的反向互补序列。请实际调用工具，并报告使用的工具名称和返回结果。

预期过程为 `search_dnakit_tools` → `describe_dnakit_tool` → `call_dnakit_tool`，
选中的功能为 `dnakit.ops.reverse_complement`，结果中的序列应为 **`CGTT`**。
结果可能包含结构化记录，不能要求整个返回值只是一段字符串。该示例无需深度模型权重，
默认关键词搜索也无需 embedding 或额外的 DeepSeek API key；对话模型由所用 Agent 提供。

### 2.4 Codex IDE 扩展

同一 Codex 主机环境中的 CLI 和 IDE 扩展共享 MCP 配置。按上面配置后重启扩展，
在齿轮菜单的 **MCP servers** 中检查连接。也可以通过 **Add server** 选择 **STDIO**，
填写实际的 `dnakit-mcp` 绝对路径后保存并重启扩展。

## 3) 其他 Agent 客户端

对于使用 `mcpServers` JSON 配置格式的客户端，可以在其指定的配置文件中加入：

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/dnakit-mcp"
    }
  }
}
```

将 `command` 替换为安装环境中可执行文件的绝对路径。配置文件位置、启用和重启方式由客户端决定；
上面的 JSON 不是 Codex `config.toml` 的格式。连接后同样用 `AACG` 示例验证实际调用。

## 4) 服务如何启动、在哪里执行

完成客户端配置后，**Agent 自动启动 `dnakit-mcp` 子进程，通过标准输入／输出交换 MCP 消息**。
通常无需另开终端启动服务。函数在服务所在机器的 Python 环境中运行，使用该环境的文件和计算资源。

手动执行下面的命令仅适合检查启动错误，需使用同一安装环境的可执行文件：

```bash
"$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
```

它可能输出启动日志并等待协议输入，不会出现普通聊天界面。另一终端里单独启动的进程不会自动
被 Agent 接管。本入口默认不监听 HTTP 端口；`http://127.0.0.1:8080/DNAKit/...` 这类文档预览
地址是文档站点地址，不能填作 DNAKit 的 MCP 服务 URL。

## 5) Agent MCP 可见入口

MCP 服务只暴露八个紧凑入口，避免一次发送数百个工具 Schema：

1. `dnakit_catalog_stats`：查看工具总数；
2. `list_dnakit_categories`：查看功能分类；
3. `list_dnakit_tools`：分页列出工具；
4. `search_dnakit_tools`：按功能搜索，默认关键词匹配，可选云端语义／混合检索；
5. `describe_dnakit_tool`：查看完整输入 Schema 和安全属性；
6. `call_dnakit_tool`：执行具体 DNAKit 功能。
7. `search_dnakit_skills`：在本地检索已收录的工作流技能摘要。
8. `read_dnakit_skill`：按需、分页读取指定技能或辅助文本。

安装、依赖与引用见 [Skills 教程](skills.md)。发现技能不会自动安装其底层工具；只有选中的工作流内容需要进入对话。

可发现和可调用的函数数量随安装版本变化，以 `dnakit_catalog_stats` 返回的统计为准。
需要 Python 回调等活动对象的接口可能保持 Python-only；普通批量任务可由 Agent 多次调用具体
DNAKit 功能完成。目录中的函数数与论文中的描述符、算法等能力条目数采用不同计数口径。

工具目录从 `PUBLIC_TOOL_MODULES` 中各模块的 `__all__` 动态发现函数，执行时调用相同的 Python
函数对象，不维护独立的 Agent 版本。修改已有函数实现、签名或文档后，重启 Python／MCP 服务即可
在同一源码环境中生效；已经运行的进程有导入与目录缓存。若 Agent 使用另一个虚拟环境或普通安装包，
需要更新那个环境中的 DNAKit。新增函数须公开导出，新增模块还须加入工具模块列表；需要回调等活动
Python 对象的接口仍可能无法通过 JSON 调用。

`input_schema.description` 保留函数说明（最多 8,000 字符），让仅传入 Schema 的 LLM 也能看到
参数的语义约束。搜索词与调用参数不同：搜索可用 `download human hg38 reference genome fasta`，
但 `download_genome` 的 `query` 应为 `hg38`；不要把搜索词直接复制为组装标识。

<span id="tool-search-modes"></span>

## 6) 只获取候选工具，不执行任务

客户端只调用 `search_dnakit_tools` 即可获取候选，例如下面的 MCP `tools/call` 参数：

```json
{
  "name": "search_dnakit_tools",
  "arguments": {"query": "download human hg38 reference genome", "limit": 3}
}
```

响应包含 `ok`、`count` 和 `tools`，不会下载 hg38 或执行候选工具。需要参数说明时再调用 `describe_dnakit_tool`；是否调用 `call_dnakit_tool` 由客户端决定。默认 `mode="keyword"`，不需要 embedding API key，也不调用聊天 LLM。

## 7) 可选的云端 embedding 与混合检索

在 MCP 服务进程中配置独立的 embedding 服务即可启用 `semantic` 和 `hybrid` 模式。
服务端配置或缓存路径不由 LLM 参数指定。以下地址、模型、密钥须替换为服务商提供的实际值，
也可以通过客户端支持的环境变量／密钥管理功能传入。

Codex 在现有服务配置中添加环境变量子表：

```toml
[mcp_servers.dnakit.env]
DNAKIT_EMBEDDING_BASE_URL = "https://your-embedding-service.example/v1"
DNAKIT_EMBEDDING_MODEL = "your-text-embedding-model"
DNAKIT_EMBEDDING_API_KEY = "your-embedding-api-key"
```

使用 `mcpServers` JSON 格式的其他客户端：

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/.venv/bin/dnakit-mcp",
      "env": {
        "DNAKIT_EMBEDDING_BASE_URL": "https://your-embedding-service.example/v1",
        "DNAKIT_EMBEDDING_MODEL": "your-text-embedding-model",
        "DNAKIT_EMBEDDING_API_KEY": "your-embedding-api-key"
      }
    }
  }
}
```

启动时不生成向量；即使配置了这些变量，**默认搜索仍然是关键词检索**。客户端明确选择混合检索：

```json
{
  "name": "search_dnakit_tools",
  "arguments": {"query": "下载人类 hg38 参考基因组", "mode": "hybrid", "limit": 3}
}
```

将 `mode` 改为 `semantic` 即为纯向量检索。两者只返回候选，但会将查询和工具描述发送到 embedding API，并在运行 MCP 服务的机器上缓存工具向量。可选环境变量 `DNAKIT_EMBEDDING_CACHE_DIR` 指定该机器上的缓存目录。模型或 API 配置只在服务进程启动时读取，修改后需要重启。缺少配置／密钥或 API 调用失败时，返回 `ok=false` 和具体错误，不会静默切换检索方式。

需要自定义批次、维度、超时或关闭磁盘缓存时，可使用 Python 创建服务：

```python
import os

from dnakit.tools import AgentToolRegistry, EmbeddingConfig, create_server

registry = AgentToolRegistry(embedding_config=EmbeddingConfig(
    model=os.environ["DNAKIT_EMBEDDING_MODEL"],
    base_url=os.environ["DNAKIT_EMBEDDING_BASE_URL"],
    cache_dir=None,
))
create_server(registry).run(show_banner=False)
```

Python、内置 Agent 和 MCP 使用同一个 `registry.search()`。模式说明、缓存更新及其他参数见 [Python API：关键词、云端 embedding 与混合检索](python_api.md#tool-search-modes)。

## 8) 常见问题

| 现象 | 检查方式 |
| --- | --- |
| 找不到 `dnakit-mcp` | 确认在配置所指的 Python 环境安装了 `dnakit[agent]`；使用可执行文件的绝对路径。 |
| 配置已登记但连接失败 | 查看 `/mcp` 或客户端启动日志；核对运行环境、路径和 `startup_timeout_sec`。 |
| 手动运行后一直等待 | 本地 `stdio` 服务需要客户端交换协议消息；通过客户端配置启动后再发送自然语言任务。 |
| 提示缺少依赖、外部程序或权重 | 在 MCP 服务使用的环境中安装该功能的依赖，按错误信息指定模型或外部程序路径。 |
| 文件写入或下载被拒绝 | 在任务中明确需要执行的操作和输出路径，由 Agent 按任务授权在 `call_dnakit_tool` 中设置 `allow_side_effects=true`。 |
| 修改代码后仍是旧行为 | 确认服务环境安装的是目标版本或可编辑源码，并重启客户端中的 MCP 服务。 |
