# 大模型工具

根据使用方式选择对应页面：

- [Python API](llm_tools/python_api.md)：只获取候选工具、云端 embedding／混合检索与通用大模型调用，无需安装 MCP。默认检索为本地关键词匹配。
- [MCP](llm_tools/mcp.md)：安装依赖、启动服务、配置客户端与使用八个工具入口。
- [Skills](llm_tools/skills.md)：搜索 DNA 工作流、安装统一入口或单个技能、检查依赖和引用来源。
- [Agent](llm_tools/agent.md)：使用内置 `run_agent()` 自动调用工具与补充 Python，按需生成工具、优化说明或兼容修复源码，并配置本地结果缓存、断点续跑和知识库检索（RAG）。

<span id="custom-llm"></span>
<span id="llm-file-input"></span>

[通用大模型调用](llm_tools/python_api.md#custom-llm)与[本地文件输入](llm_tools/python_api.md#llm-file-input)的完整说明见 Python API 页面。
