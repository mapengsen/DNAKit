# Python API

本页介绍只获取候选工具、可选的云端 embedding／混合检索，以及通用大模型调用。均可直接在 Python 中使用，无需安装或启动 MCP，也无需安装 Agent 可选依赖。工具检索默认使用本地关键词匹配，不需要 API key。

安装 DNAKit 核心包即可：

```bash
python -m pip install dnakit
```

如需接入支持 MCP 的客户端，请参阅 [MCP](mcp.md)。

<span id="agent-tool-calling"></span>

自然语言驱动的 `run_agent()`、Python 补充执行、工具生成与优化、本地结果缓存及知识库检索（RAG）见 [Agent](agent.md#agent-actions) 页面，其中提供开关、参数、示例与执行限制。

<span id="tool-discovery"></span>

## 1) 只获取候选工具，不执行任务

其他程序只想知道“可以使用哪些工具”时，直接调用 `registry.search()`，无需运行 `run_agent()`：

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry()
tools = registry.search("download human hg38 reference genome", limit=3)
print([tool.to_summary() for tool in tools])
```

返回按相关性排序的 `AgentToolSpec` 元组；`to_summary()` 提供工具名称、功能描述、分类等信息。这里仅查询本地工具目录，不调用聊天模型、不下载 hg38，也不执行候选工具。关键词模式建议使用英文功能词或函数名称，`limit` 默认为 20，可设置为 1–100。

需要查看完整参数 Schema 时，单独调用：

```python
info = registry.describe("dnakit.references.download_genome")
print(info["input_schema"])
```

工具执行是另一个步骤，由调用方决定是否使用 `registry.execute()`。单独运行 `search()` 或 `describe()` 不会触发它。

<span id="tool-search-modes"></span>

## 2) 关键词、云端 embedding 与混合检索

六种模式共用 `registry.search()`，返回格式相同，**默认始终是 `keyword`**。即使配置了 embedding 服务，省略 `mode` 也不会联网。

| `mode` | 检索方式 | 是否调用 embedding API |
| --- | --- | --- |
| `keyword`（默认） | 匹配工具名称、功能短说明、中英文分类 | 否 |
| `semantic` | 将查询和上述工具描述转成向量，按余弦相似度排序 | 是 |
| `hybrid` | 合并关键词排名和语义排名 | 是 |
| `bm25` | 本地 Okapi BM25 词项检索 | 否 |
| `bm25_hybrid` | BM25 与语义排名的等权 RRF 融合 | 是 |
| `rerank` | BM25＋Dense 召回固定候选池，再通过独立重排 API 排序 | 是，另调用重排 API |

混合检索使用等权倒数排名融合（RRF，`k=60`），完整函数名或短函数名的精确匹配优先。分类 `category` 和 `include_python_only` 过滤在检索前生效。空白查询用于列出工具，或分类过滤后没有候选时，不请求 embedding API。向量模式返回最相近的候选，不保证候选一定能完成需求；执行前仍需查看 Schema。

### 配置独立的云端 embedding 服务

支持 OpenAI-compatible `POST /embeddings` 文本向量接口；不需要安装 OpenAI SDK、PyTorch、向量数据库或本地 embedding 模型。可接入提供该兼容接口的服务，例如[阿里云百炼文本向量 API](https://help.aliyun.com/zh/model-studio/embedding)。私有协议不自动适配。

先在终端或 PyCharm 运行配置中填写以下环境变量，将占位值替换为自己账户的实际配置：

```bash
export DNAKIT_EMBEDDING_BASE_URL="https://your-embedding-service.example/v1"
export DNAKIT_EMBEDDING_MODEL="your-text-embedding-model"
export DNAKIT_EMBEDDING_API_KEY="your-embedding-api-key"
```

`BASE_URL` 填 API 前缀，DNAKit 自动追加一次 `/embeddings`，也接受完整的该接口地址。选择支持中英文语义检索的文本 embedding 模型。这里的配置与 `LLMConfig`、`DEEPSEEK_API_KEY` 相互独立，不会自动使用聊天模型或聊天密钥。

```python
import os

from dnakit.tools import AgentToolRegistry, EmbeddingConfig

registry = AgentToolRegistry(embedding_config=EmbeddingConfig(
    model=os.environ["DNAKIT_EMBEDDING_MODEL"],
    base_url=os.environ["DNAKIT_EMBEDDING_BASE_URL"],
))
tools = registry.search("下载人类 hg38 参考基因组", mode="semantic", limit=3)
print([tool.to_summary() for tool in tools])
```

仅获取混合检索候选：

```python
tools = registry.search("下载人类 hg38 参考基因组", mode="hybrid", limit=3)
print([tool.to_summary() for tool in tools])
```

回到默认关键词检索只需省略 `mode`，仍然不执行任务：

```python
tools = registry.search("download hg38 genome", limit=3)
```

### 缓存与请求参数

首次语义／混合检索先请求查询向量，再分批生成尚未缓存的工具向量。后续相同工具描述可复用本地缓存，每次搜索仍请求查询向量。缓存根据接口地址、模型、指定维度和工具描述区分；修改工具名称、短说明或分类后，重启进程会为变化的描述重新生成向量。服务商在同一个模型 ID 下更新权重时，应手动清理对应缓存，避免混用新旧向量。

| `EmbeddingConfig` 参数 | 默认值与用途 |
| --- | --- |
| `model`、`base_url` | 必填；文本 embedding 模型和 API 地址 |
| `api_key` | 不填写时读取 `api_key_env`，不会读取聊天密钥 |
| `api_key_env` | `DNAKIT_EMBEDDING_API_KEY`；可改为服务商的环境变量名 |
| `cache_dir` | `$XDG_CACHE_HOME/dnakit/tool_embeddings`，未设置时使用 `~/.cache/dnakit/tool_embeddings`；`None` 表示仅内存缓存 |
| `batch_size` | `10`；按服务商限制调整每批工具描述条数 |
| `dimensions` | `None`，使用模型默认维度；仅在服务商支持时指定 |
| `timeout_seconds` | `60`；单次 HTTP 请求超时 |
| `max_request_bytes` | `1_000_000`；单批请求字节上限 |
| `max_response_bytes` | `16_000_000`；单批响应字节上限 |

云端只接收搜索文本和工具的名称、短说明、分类；这些接口不会读取输入路径指向的 DNA 文件。向量与相似度排序保存在／运行于调用 registry 的机器上。磁盘缓存不保存 API key 或查询原文。构造配置、创建 registry 和默认关键词检索均不创建缓存、不请求云端。明确选择语义／混合模式后需要相应配置；缺少密钥、网络错误、超限或非法向量会明确报错，不会静默切回关键词模式，也不会自动重试付费请求。

将这个 registry 用于 `run_agent()` 的示例见 [Agent：配置工具检索](agent.md#tool-search-modes)。

MCP 也使用相同检索实现，参阅 [MCP 的只获取候选与云端检索配置](mcp.md#tool-search-modes)。

### BM25、BM25＋Dense 与检索后重排

`bm25` 不需要模型、API key 或新增依赖，可以直接运行：

```python
from dnakit.tools import AgentToolRegistry, BM25Config

registry = AgentToolRegistry(bm25_config=BM25Config(k1=1.2, b=0.75))
tools = registry.search("download human reference genome", mode="bm25", limit=10)
print([tool.name for tool in tools])
```

BM25 使用正值 IDF `log(1 + (N - df + 0.5) / (df + 0.5))`；默认 `k1=1.2`、`b=0.75`。统计量在过滤后的候选目录内计算，仅返回得分大于零的工具，同分按完整工具名排序。文本来自工具名、检索说明（默认短说明）和中英文类别，与语义检索使用相同字段。英文转小写，下划线／连字符标识符同时保留整体和分词；中文文档使用单字及相邻双字，多个连续汉字的查询使用相邻双字。没有自动翻译、词干化或停用词删除。中文分词不会把中文描述翻译成工具的英文说明。

`bm25_hybrid` 复用 `EmbeddingConfig`，以等权 RRF（`k=60`）融合完整 BM25 正分排名和完整语义排名，最后截取 `limit`。它不使用原 `hybrid` 的工具名称精确匹配优先规则。

`rerank` 先运行 `bm25_hybrid`，再将前 `candidate_limit` 个候选发送到重排 API。使用 [Cohere Rerank 风格协议](https://docs.cohere.com/v2/reference/rerank)：`POST /rerank`，请求包含 `model`、`query`、字符串数组 `documents`、`top_n`；响应包含 `results`，每项具有 `index` 和数值 `relevance_score`。其他协议需要另行适配，不能直接把聊天或 embedding 接口地址当作重排地址。

```bash
export DNAKIT_RERANK_BASE_URL="https://your-reranker-service.example/v2"
export DNAKIT_RERANK_MODEL="your-reranker-model"
export DNAKIT_RERANK_API_KEY="your-reranker-api-key"
```

配置好本节前面的 embedding 环境变量后：

```python
import os

from dnakit.tools import AgentToolRegistry, BM25Config, EmbeddingConfig, RerankerConfig

registry = AgentToolRegistry(
    bm25_config=BM25Config(k1=1.2, b=0.75),
    embedding_config=EmbeddingConfig(
        model=os.environ["DNAKIT_EMBEDDING_MODEL"],
        base_url=os.environ["DNAKIT_EMBEDDING_BASE_URL"],
    ),
    reranker_config=RerankerConfig(
        model=os.environ["DNAKIT_RERANK_MODEL"],
        base_url=os.environ["DNAKIT_RERANK_BASE_URL"],
        candidate_limit=50,
    ),
)
for mode in ("bm25", "bm25_hybrid", "rerank"):
    tools = registry.search("下载人类 hg38 参考基因组", mode=mode, limit=10)
    print(mode, [tool.name for tool in tools])
```

`RerankerConfig` 的 `model`、`base_url` 必填；地址可为 API 前缀或完整 `/rerank` 地址。默认独立读取 `DNAKIT_RERANK_API_KEY`，可用 `api_key_env` 或隐藏于 repr 的 `api_key` 覆盖，不借用聊天或 embedding 密钥。`candidate_limit` 默认 50，允许 1–1000；`limit` 仍为 1–100，且不能超过候选池上限。可设置 `timeout_seconds=60`、`max_request_bytes=1_000_000`、`max_response_bytes=2_000_000`。

重排每次对整个固定候选池评分，再截取输出；相同分数保留召回顺序，不随 Top-1／5／10 改变候选池。所有索引、结果数量及有限数值均被核对，缺失或重复结果报错。缺少配置、密钥、网络失败或超限均明确报错，不自动降级或重试。重排不执行候选工具，也不持久化查询、重排分数或密钥；其 embedding 缓存仍按前述规则保存。`run_agent(..., registry=registry, search_mode="rerank", search_limit=10)` 使用同一实现；工具目录刷新会保留这些配置。知识库 RAG 的模式仍单独配置。

<span id="custom-llm"></span>

## 3) 通用大模型调用

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
