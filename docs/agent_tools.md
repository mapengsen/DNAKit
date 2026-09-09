# Agent 与 MCP 工具

DNAKit 可以通过 MCP 被 Codex、Claude、Cursor、VS Code 等支持 MCP 的 Agent 调用。Agent
不直接加载全部 Python 函数，而是先搜索工具目录、查看参数 Schema，再执行选中的公开 API。

本页也收录可直接嵌入 Python 流程的[通用大模型调用](#custom-llm)，包括[指定本地文件路径](#llm-file-input)。该功能不需要启动 MCP 服务，也不需要安装 Agent 可选依赖。

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

<span id="custom-llm"></span>

## 7) 通用大模型调用（Python API）

`dnakit.llm` 可插入任意 Python 处理流程：自定义提示词 + 当前步骤的数据 → 大模型输出 → 继续本地代码。它不局限于 DNA 或训练文件结构分析，也不需要使用 DNAKit CLI。

公开入口为 `LLMConfig`、`LLMResult`、同步函数 `call_llm()` 和异步函数 `acall_llm()`。仅依赖 Python 标准库和 DNAKit 核心，不需要安装 OpenAI SDK、PyTorch 或训练扩展。首版使用 OpenAI-compatible Chat Completions 协议；不自动适配其他厂商的私有协议。

### 支持的服务商

不是只能使用 DeepSeek。`provider` 预设支持下面六家服务商，并保留 `custom` 自定义兼容接口；同步和异步调用使用同一套配置。模型名称必须由用户填写，不会自动选用收费模型。

| 服务商                                                                                                             | `provider`   | 默认 API base URL                                           | 默认密钥环境变量       |
| ------------------------------------------------------------------------------------------------------------------ | -------------- | ----------------------------------------------------------- | ---------------------- |
| DeepSeek                                                                                                           | `deepseek`   | `https://api.deepseek.com`                                | `DEEPSEEK_API_KEY`   |
| [Google Gemini](https://ai.google.dev/gemini-api/docs/openai)                                                       | `gemini`     | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY`     |
| [千问 / 阿里云百炼](https://help.aliyun.com/zh/model-studio/qwen-api-via-openai-chat-completions)                   | `qwen`       | `https://dashscope.aliyuncs.com/compatible-mode/v1`       | `DASHSCOPE_API_KEY`  |
| [Kimi / Moonshot](https://platform.kimi.com/docs/get-api-key)                                                       | `kimi`       | `https://api.moonshot.cn/v1`                              | `MOONSHOT_API_KEY`   |
| [OpenRouter](https://openrouter.ai/docs/quickstart)                                                                 | `openrouter` | `https://openrouter.ai/api/v1`                            | `OPENROUTER_API_KEY` |
| [OpenAI](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create) | `openai`     | `https://api.openai.com/v1`                               | `OPENAI_API_KEY`     |
| 自定义兼容接口                                                                                                     | `custom`     | 必须显式填写`base_url`                                    | `DNAKIT_LLM_API_KEY` |

```python
from dnakit.llm import LLMConfig, available_llm_providers, call_llm

print([provider.name for provider in available_llm_providers()])  # 仅本地查询
config = LLMConfig(
    provider="gemini",  # 可换为 qwen、kimi、openrouter、openai、deepseek
    model="填写你账户支持的实际模型名",
)
result = call_llm("用中文总结这段内容。", input_data="Your input text", config=config)
print(result.text)
```

`available_llm_providers()` 返回 `LLMProvider` 预设对象；`get_llm_provider(name)` 解析单个预设。支持 `google` / `google-gemini`、`dashscope` / `alibaba`、`moonshot` 和 `openai-compatible` 别名。错误的 provider 名称会在请求前报错，不会回退至 DeepSeek。

兼容与地域注意事项：

- Gemini 使用 Google 官方的 **OpenAI 兼容入口**和 Gemini API key，不是原生 `generateContent` URL，也不包含 Vertex AI 的服务账号认证。
- Qwen 预设为北京旧域名，官方仍保留支持。可用 `base_url` 指定账户所属地域/业务空间地址，例如 `https://<实际WorkspaceId>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` 或新加坡的 `https://<实际WorkspaceId>.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1`。密钥与地域必须匹配。
- Kimi 国际账户可覆盖为 `base_url="https://api.moonshot.ai/v1"`，并使用对应平台密钥。
- OpenRouter 的 `model` 使用其目录中的完整模型 ID；选择服务商不等于所有模型都支持 JSON、思考模式或相同参数。
- 显式填写的 `base_url`、`api_key_env` 和 `api_key` 优先于预设。为兼容已有代码，不填写 `provider` 时仍使用 DeepSeek 地址及 `DNAKIT_LLM_API_KEY`，不会悄悄改读其他环境变量。

### 模型扩展参数

`extra_body` 将有限、可序列化的 JSON 字段合并到 HTTP 请求体顶层，可传服务商的思考参数等。它不能覆盖 `model`、`messages`、`stream`、token 限额、输出格式或认证字段，也不能开启尚未支持的工具/多模态请求。参数是否有效仍由具体模型决定；默认不注入模型专属参数。

```python
config = LLMConfig(
    provider="qwen",
    model="填写支持关闭思考的实际模型名",
    extra_body={"enable_thinking": False},
)
```

对于要求关闭思考才能使用非流式或 JSON 输出的 Qwen 模型，可按模型文档显式配置上例。OpenAI/Kimi 的部分模型可传 `extra_body={"reasoning_effort": "low"}`；Gemini 的对应 Google 扩展字段需要保留 HTTP 层嵌套，例如 `extra_body={"extra_body": {"google": {"thinking_config": {"thinking_level": "low"}}}}`。这些只是特定参数形式的例子，不保证适用于同一服务商的所有模型。

### 文本输入与输出

```python
import os

from dnakit.llm import LLMConfig, call_llm

config = LLMConfig(
    provider="deepseek",  # 可切换为 gemini、qwen、kimi、openrouter 或 openai
    model=os.environ["DNAKIT_LLM_MODEL"],  # 使用你账户支持的实际模型名
    api_key_env="DNAKIT_LLM_API_KEY",  # 也支持 api_key=os.environ["...照你的变量名..."]
    max_output_tokens=1024,
    timeout_seconds=60,
)

result = call_llm(
    "请将下面这段文字翻译成中文，只返回译文。",
    input_data="The sample passed the quality check.",
    config=config,
)
translated_text = result.text
print(translated_text)
print(dict(result.usage))  # 服务商返回的 token 用量，不是费用估算
```

在 PyCharm 的运行配置中设置两个环境变量即可；无需把密钥写进代码。`base_url` 是 API 地址，不是 `platform.deepseek.com/usage` 之类的账户页面。模型名称按服务商当前可用列表配置，DNAKit 不替用户选择或切换模型。

<span id="llm-file-input"></span>

### 直接指定本地文件路径

`call_llm()` 和 `acall_llm()` 均支持 `file_path`，可直接读取并处理一个本地文本文件，无需手动 `read_text()`：

```python
from dnakit.llm import LLMConfig, call_llm

result = call_llm(
    "用中文概括文件内容，只根据原文回答。",
    file_path="/mnt/d/Pycharm_workspace/DNA/DNAKit/LICENSE",
    config=LLMConfig(
        provider="gemini",  # 也可使用其他已支持的服务商
        model="填写实际模型名",
        api_key_env="GEMINI_API_KEY",
    ),
)
print(result.text)
```

文件路径支持字符串或 `Path`、绝对路径或相对于当前工作目录的路径，包含空格也可使用；Python 字符串中的下划线不需要加反斜杠转义。Windows 原生 Python 应使用 Windows 路径，例如 `r"D:\Pycharm_workspace\DNA\DNAKit\LICENSE"`；`/mnt/d/...` 是 WSL/Linux 路径。

**指定 `file_path` 表示允许发送该文件的完整文本内容。** 请求只附带文件名，不自动附带绝对目录；不会读取邻近文件、扫描目录、递归读取文件中的链接、执行代码或自动分批请求。`input_data` 与 `file_path` 不能同时使用。原先 `input_data="/path/to/file"` 的行为保持不变：它仍只是一个字符串，不会触发读取。

可以先完全在本地预览：

```python
from dnakit.llm import read_llm_file

text = read_llm_file("LICENSE")  # 仅读取，不联网
print(text[:200])
```

默认 `file_encoding="utf-8-sig"` 支持 UTF-8 及 UTF-8 BOM；其他编码须显式指定，例如 `file_encoding="gb18030"` 或 `"utf-16"`。`read_llm_file()` 对应参数名为 `encoding`。

默认 `max_file_bytes=64000` 限制原文件字节数，`read_llm_file()` 对应参数名为 `max_bytes`。完整请求还受 `LLMConfig.max_request_bytes` 限制，后者包括提示词、文件名、JSON 转义和其他请求字段；提高文件上限不等于提高请求上限。超限时在请求前报错，不会偷偷截断、抽样或拆成多次付费调用。

支持无后缀 LICENSE、TXT、Markdown、CSV/TSV、JSON/JSONL、YAML、XML、Python 等可解码文本；表格与代码按原文发送，不执行或自动清洗。PDF、Word/Excel 二进制文档、图片、压缩包等不能直接处理，需先自行提取文本再传入 `input_data`。空文件、目录、设备文件、错误编码和超限文件会明确报错。完整示例为 `others/examples/custom_llm_file.py`。

### JSON 输入与结构化输出

```python
rows = [
    {"id": "a", "note": "quality check passed"},
    {"id": "b", "note": "measurement needs repeating"},
]
result = call_llm(
    '将 note 翻译为中文。返回 JSON 对象 {"items": [{"id": "原 ID", "note_zh": "译文"}]}。',
    input_data=rows,
    output_format="json",
    system_prompt="你是实验记录翻译助手。保留原始 ID，不补充未提供的实验结果。",
    config=config,
)
output = result.data  # 普通 dict；文本模式下为 None
if output is None or not isinstance(output.get("items"), list):
    raise ValueError("输出缺少预期的 items 列表")
# 在这里继续验证记录 ID、字段类型等，再交给下一步 Python 代码。
```

JSON 模式请求 `response_format={"type": "json_object"}` 并校验返回语法、非有限数值、重复键和顶层对象类型；不保证业务字段符合预期，也不保证模型判断正确。需要数组时，将其放在对象字段内，例如 `{"items": [...]}`。所选模型需要支持 JSON mode；服务商拒绝请求或返回非法 JSON 时会明确报错，DNAKit 不会静默切换输出模式或模型。

`LLMResult` 提供 `text`、`data`、`model`、`response_id`、`finish_reason`、`usage` 和 `to_dict()`。`data` 每次返回独立的普通字典，修改它不会改变原结果。结果不保留请求提示词或原始输入，但模型输出本身仍可能包含敏感内容。

### 任意 Python 步骤与异步程序

你可以把上一步产生的字符串、有限数值、布尔值、列表、元组或字符串键字典传入 `input_data`。`None` 表示不额外附加数据，只使用提示词。DataFrame、NumPy 数组、自定义对象需要由调用方先明确转换；例如只共享表格样本时使用 `df.head(5).to_dict(orient="records")`，并检查缺失值和非有限数值。

```python
from dnakit.llm import acall_llm


async def summarize_step(current_result, config):
    result = await acall_llm(
        "用两句话概括结果，不推断未提供的信息。",
        input_data=current_result,
        config=config,
    )
    return result.text
```

`acall_llm()` 在线程中进行 HTTP 调用，不阻塞 asyncio 事件循环。取消 await 不代表已提交的请求被取消，服务商仍可能完成并计费。首次连接、每次调用均无自动重试、无隐式会话历史；循环或并发调用的次数由你的 Python 代码控制。
