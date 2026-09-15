# Python API

This page covers natural-language tool execution, candidate-only discovery, optional cloud embedding/hybrid search, and general-purpose LLM calls. All run directly in Python without an MCP server or the Agent optional dependencies. Tool search defaults to local keyword matching and needs no API key.

Install DNAKit core:

```bash
python -m pip install dnakit
```

For MCP-capable clients, see [MCP](mcp.md).

<span id="agent-tool-calling"></span>

## 1) Natural-language tool execution (no MCP installation required)

`run_agent()` sends a natural-language request to the LLM and loops through tool search, schema
inspection, execution, and summarization. The user script does not need to hard-code specific tools
or a workflow:

```python
import argparse
from pathlib import Path

from dnakit.llm import LLMConfig
from dnakit.tools import run_agent

parser = argparse.ArgumentParser(description="Use natural language to call DNAKit tools with DeepSeek")
parser.add_argument("--model", default="deepseek-flash", help="DeepSeek model name")
parser.add_argument("--api-key", help="API key; reads DEEPSEEK_API_KEY when omitted")
parser.add_argument("--output", type=Path, help="Optional output directory")
parser.add_argument("--request", help="Natural-language task; defaults to downloading hg38")
parser.add_argument("--max-steps", type=int, default=50, help="Maximum number of LLM requests")
args = parser.parse_args()

config = LLMConfig(
    provider="deepseek",
    model=args.model,
    api_key=args.api_key,
    extra_body={"thinking": {"type": "disabled"}},
)
request = args.request or (
    f"Download the human hg38 reference genome FASTA to "
    f"{args.output or Path(__file__).resolve().parents[2] / 'data/fig3'} without overwriting existing files."
)
request += f"\nOutput directory: {args.output.expanduser().resolve()}" if args.output else ""
result = run_agent(request, config=config, allow_side_effects=True, max_steps=args.max_steps)
print(result.text)
```

Run the example with:

```bash
python write_paper/figure/fig3/01_download_hg38.py \
  --api-key "key.........." \
  --model deepseek-v4-pro \
  --request "Download the human hg38 reference genome FASTA to write_paper/data/fig3"
```

### Resolve missing URLs and pass files between tools

A task may specify only a dataset name, subset and output directory. `run_agent()` already
uses `call_llm()` at every step. To obtain missing external information, the model can discover:

| Tool | Purpose |
| --- | --- |
| `dnakit.search.web_search` | Search the internet for candidate pages, returning titles and URLs |
| `dnakit.search.read_webpage` | Read an HTTPS page or text/JSON API, returning text and observed links |
| `dnakit.download.download_url` | Save a confirmed file URL to a local directory |

The Agent's `search` action still searches the installed tool catalog. Web search uses Tavily
when `TAVILY_API_KEY` is set, otherwise keyless DuckDuckGo HTML; `provider` can select explicitly.
No additional Python dependency is needed. The search key stays in the process environment.
A DuckDuckGo challenge is reported as an error, never as an empty successful search. Configure
Tavily or read a known candidate page to verify it. Web content is untrusted data, not Agent
instructions. No page scripts run and no unobserved download links are invented by the reader.

```python
from dnakit.search import read_webpage, web_search

candidates = web_search("NT revised H2AFZ dataset")
page = read_webpage(candidates.records[0]["url"], match="H2AFZ")
print(page.to_dict())
```

The reader defaults to 8,000 text characters and 30 links. `match` filters by substring;
`offset` and `link_offset` paginate, with explicit truncation flags and continuation offsets.
Binary resources belong in the downloader. `call_llm(file_path=...)` remains a text-file input.

The existing ML trainer accepts `CustomPropertyTrainingConfig(cv_folds=10, ...)`. It validates
all folds of a training-only source, saves validation-fold membership, out-of-fold predictions
and metrics, then refits the final model on all source rows. This is supported only with
`finetune_method="machine_learning"`; split ratios are ignored and mapped split columns are
rejected. An independent test file is used only for subsequent prediction. Each fold fits its
own preprocessing and model using training rows alone.
`predict_custom_properties(..., output_path="predictions.csv")` saves predicted labels locally.
Its `to_dict()`/Agent result then contains a file summary; Python prediction/probability attributes
remain available. Feed the CSV to `plot_csv_distributions` and save the plot with `save_image`.

Agent `call` arguments can reference earlier successful results, for example
`{"artifact": {"$ref": {"step": 12, "path": ["label"]}}}`. Step numbers start at 1;
the path selects result keys/list indices, with an empty path selecting the whole result.
References resolve locally, so the model need not regenerate SVG or copy large payloads.
Direct `registry.execute()` calls do not interpret these Agent references.

Automatic field detection checks candidate ID uniqueness for built-in readers. Duplicate names
use stable source-row IDs without dropping samples, with an inspection warning. Explicitly
mapped `DataMapping.id_field` values still require unique IDs.

### Python fallback and failure recovery

When tools cannot directly exchange data or a processing step is missing, the same
Agent can write and execute Python, then resume registered tool calls. Enable it explicitly:

```python
from dnakit.tools import PythonExecutionConfig, run_agent

result = run_agent(
    request,
    config=config,
    allow_side_effects=True,
    code_execution=PythonExecutionConfig(
        workspace="write_paper/data/my_run",
        timeout_seconds=600,
        max_attempts=5,
    ),
)
print(result.text)
print(result.execution_dir)
```

The model prefers existing tools. After finding a capability gap or receiving an error,
it can choose the `execute_python` runtime action with `reason`, complete `code`,
`inputs`, and `expected_outputs`. Input `$ref` values resolve locally into a JSON
file passed as `sys.argv[1]`; the model need not copy sequences, predictions or SVG
into source code. Subsequent tool calls can reference execution output paths.

Scripts run in fresh subprocesses with the current Python interpreter. Nonzero exit,
timeout, or missing/empty declared output files produce failed observations that the
model can use to correct its code. Identical code, inputs and output lists reuse the
previous observation. Inspection scripts count toward `max_attempts`; all actions
also count toward `max_steps`. Defaults are 120 seconds per execution and 5 attempts.
Each stdout/stderr tail sent to the model is limited to 16,000 bytes; full logs stay
on disk. POSIX cleanup terminates the process group on timeout, interruption or script
exit. On Windows, only the directly launched Python process is terminated.

Audit files live under `<workspace>/.dnakit_agent/run_*/`: `history.jsonl` records
actions and observations, and `summary.json` records status and actual execution count.
Each `attempt_*/` contains `script.py`, `inputs.json`, `request.json` (reason and output
paths), `stdout.log`, `stderr.log`, and `result.json` (exit code, time, file checks).
Automatic code attempts are distinct from human interventions.

This is **not a sandbox**: scripts have the current user's filesystem, network and
environment access. The workspace is a cwd, not an access boundary. Execution is
disabled by default and also requires caller-owned `allow_side_effects=True`.
File existence is not scientific validation: scripts should assert task-specific
row counts, columns and sequence lengths. Missing credentials, dependencies or
external resources can still block completion.

The Fig. 3 `02_h2afz_svm.py` example and Fig. 2 `agent_comparison/run_task.py --mode dnakit`
support `--allow-code-execution`, `--code-timeout-seconds`, and `--max-code-attempts`.
The native Chat comparison keeps its original behavior.

<span id="tool-discovery"></span>

## 2) Get candidate tools without executing a task

Programs that only need to know which tools are available can call `registry.search()` directly, without running `run_agent()`:

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry()
tools = registry.search("download human hg38 reference genome", limit=3)
print([tool.to_summary() for tool in tools])
```

The result is a relevance-ranked tuple of `AgentToolSpec` objects; `to_summary()` includes the name, description, category, and other metadata. This searches the local catalog without calling a chat model, downloading hg38, or executing candidates. For keyword mode, use English capability terms or function names. `limit` defaults to 20 and accepts 1–100.

Inspect the complete parameter schema separately:

```python
info = registry.describe("dnakit.references.download_genome")
print(info["input_schema"])
```

Execution is a separate step: the caller decides whether to use `registry.execute()`. Neither `search()` nor `describe()` invokes it.

<span id="tool-search-modes"></span>

## 3) Keyword, cloud embedding, and hybrid search

All three modes use `registry.search()` and return the same format. **The default is always `keyword`**: omitting `mode` makes no network requests even when an embedding service is configured.

| `mode` | Retrieval method | Calls the embedding API? |
| --- | --- | --- |
| `keyword` (default) | Match tool names, short descriptions, and bilingual categories | No |
| `semantic` | Embed the query and those tool descriptions, then rank by cosine similarity | Yes |
| `hybrid` | Combine keyword and semantic rankings | Yes |

Hybrid search uses equally weighted reciprocal rank fusion (RRF, `k=60`), prioritizing exact full or short function names. `category` and `include_python_only` filters apply before retrieval. Blank queries list tools without embedding requests; an empty filtered catalog also makes no request. Vector search returns the nearest candidates, which are not guaranteed to satisfy the task; inspect their schemas before execution.

### Configure an independent cloud embedding service

The client supports OpenAI-compatible text-vector `POST /embeddings` endpoints without an OpenAI SDK, PyTorch, vector database, or local embedding model. Services offering this interface include [Alibaba Cloud Model Studio](https://help.aliyun.com/zh/model-studio/embedding). Proprietary protocols are not adapted automatically.

Set these environment variables in your terminal or PyCharm run configuration, replacing the placeholders with your account's actual settings:

```bash
export DNAKIT_EMBEDDING_BASE_URL="https://your-embedding-service.example/v1"
export DNAKIT_EMBEDDING_MODEL="your-text-embedding-model"
export DNAKIT_EMBEDDING_API_KEY="your-embedding-api-key"
```

`BASE_URL` is the API prefix; DNAKit appends `/embeddings` once and also accepts the complete endpoint URL. Choose a text embedding model that supports Chinese/English semantic retrieval. These settings are independent of `LLMConfig` and `DEEPSEEK_API_KEY`; the chat model and its credential are never selected implicitly.

```python
import os

from dnakit.tools import AgentToolRegistry, EmbeddingConfig

registry = AgentToolRegistry(embedding_config=EmbeddingConfig(
    model=os.environ["DNAKIT_EMBEDDING_MODEL"],
    base_url=os.environ["DNAKIT_EMBEDDING_BASE_URL"],
))
tools = registry.search("Download the human hg38 reference genome", mode="semantic", limit=3)
print([tool.to_summary() for tool in tools])
```

Get hybrid-search candidates without execution:

```python
tools = registry.search("Download the human hg38 reference genome", mode="hybrid", limit=3)
print([tool.to_summary() for tool in tools])
```

Omit `mode` to return to default keyword search, still without executing a task:

```python
tools = registry.search("download hg38 genome", limit=3)
```

### Caching and request options

The first semantic/hybrid search requests the query vector, then generates uncached tool vectors in batches. Unchanged descriptions reuse the local cache; each search still requests its query vector. Cache keys distinguish endpoint, model, requested dimensions, and tool text. After tool names, short descriptions, or categories change, restart the process to regenerate only changed descriptions. If a provider changes weights under the same model ID, clear the corresponding cache to avoid mixing old and new vectors.

| `EmbeddingConfig` option | Default and purpose |
| --- | --- |
| `model`, `base_url` | Required text embedding model and API address |
| `api_key` | Omitted: read `api_key_env`; never borrow the chat credential |
| `api_key_env` | `DNAKIT_EMBEDDING_API_KEY`; can name a provider-specific variable |
| `cache_dir` | `$XDG_CACHE_HOME/dnakit/tool_embeddings`, or `~/.cache/dnakit/tool_embeddings` when unset; `None` uses memory only |
| `batch_size` | `10`; adjust descriptions per batch to provider limits |
| `dimensions` | `None`, the model's default; specify only when supported by the provider |
| `timeout_seconds` | `60` per HTTP request |
| `max_request_bytes` | `1_000_000` per request batch |
| `max_response_bytes` | `16_000_000` per response batch |

The cloud receives only query text and tool names, short descriptions, and categories. These interfaces do not open DNA files referenced by input paths. Vectors are cached and similarity ranking runs on the machine hosting the registry. Disk caches contain neither API keys nor raw queries. Configuration/registry construction and default keyword search create no cache files and make no cloud requests. Explicit semantic/hybrid searches require the corresponding configuration. Missing keys, network errors, limits, and invalid vectors produce errors without silently switching to keywords or automatically retrying paid requests.

### Use the same retrieval in the built-in Agent

Pass the configured `registry` above to `run_agent()`. The chat model still selects tools, generates arguments, and summarizes results; the embedding model only produces vectors:

```python
from dnakit.llm import LLMConfig
from dnakit.tools import run_agent

result = run_agent(
    "Calculate the molecular weight of ACGT.",
    config=LLMConfig(provider="deepseek", model="deepseek-flash"),
    registry=registry,
    search_mode="hybrid",
    search_limit=10,
)
print(result.text)
```

`run_agent()` defaults to `search_mode="keyword"` and `search_limit=5`; candidate limits accept 1–100. The model cannot change the caller's mode or limit. `max_steps=50` bounds chat LLM turns, excluding independent embedding API requests. This example executes the selected molecular-weight tool; the earlier `registry.search()` examples only return candidates.

MCP shares the same implementation; see [candidate-only discovery and cloud search configuration](mcp.md#tool-search-modes).

<span id="custom-llm"></span>

## 4) General-purpose LLM calls

`dnakit.llm` can be inserted into any Python workflow: custom prompt + data from the current step
→ LLM output → further local processing. It is not limited to DNA or training-file structure
analysis and does not require the DNAKit CLI.

The public entry points are `LLMConfig`, `LLMResult`, the synchronous `call_llm()`, and the
asynchronous `acall_llm()`. They depend only on the Python standard library and DNAKit core, with
no OpenAI SDK, PyTorch, or training extensions required. The initial implementation uses the
OpenAI-compatible Chat Completions protocol and does not automatically adapt proprietary provider
protocols.

### Supported providers

DeepSeek is not the only option. The `provider` presets support the six providers below, plus a
`custom` compatible endpoint. Synchronous and asynchronous calls use the same configuration. Users
must supply the model name; no paid model is selected automatically.

| Provider | `provider` | Default API base URL | Default key environment variable |
| --- | --- | --- | --- |
| DeepSeek | `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| [Google Gemini](https://ai.google.dev/gemini-api/docs/openai) | `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| [Qwen / Alibaba Cloud Model Studio](https://help.aliyun.com/zh/model-studio/qwen-api-via-openai-chat-completions) | `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| [Kimi / Moonshot](https://platform.kimi.com/docs/get-api-key) | `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| [OpenRouter](https://openrouter.ai/docs/quickstart) | `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| [OpenAI](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create) | `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| Custom compatible endpoint | `custom` | Explicit `base_url` required | `DNAKIT_LLM_API_KEY` |

```python
from dnakit.llm import LLMConfig, available_llm_providers, call_llm

print([provider.name for provider in available_llm_providers()])  # Local lookup only.
config = LLMConfig(
    provider="gemini",  # Or qwen, kimi, openrouter, openai, deepseek.
    model="your-actual-model-id",
)
result = call_llm("Summarize this text in Chinese.", input_data="Your input text", config=config)
print(result.text)
```

`available_llm_providers()` returns `LLMProvider` preset objects; `get_llm_provider(name)` resolves
one preset. Supported aliases include `google` / `google-gemini`, `dashscope` / `alibaba`,
`moonshot`, and `openai-compatible`. An invalid provider name fails before the request instead of
falling back to DeepSeek.

Compatibility and regional considerations:

- Gemini uses Google's official **OpenAI-compatible endpoint** and a Gemini API key, not the native `generateContent` URL or Vertex AI service-account authentication.
- The Qwen preset uses the legacy Beijing URL, which remains officially supported. Set `base_url` to the address for your account's region/workspace, such as `https://<yourWorkspaceId>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` or `https://<yourWorkspaceId>.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1` for Singapore. The key must match the region.
- International Kimi accounts can override `base_url="https://api.moonshot.ai/v1"` and use the corresponding platform key.
- OpenRouter's `model` must be the complete model ID from its catalog. Selecting a provider does not mean every model supports JSON, thinking mode, or the same parameters.
- Explicit `base_url`, `api_key_env`, and `api_key` values override presets. For backward compatibility, omitting `provider` still uses the DeepSeek URL and `DNAKIT_LLM_API_KEY`, without silently switching environment variables.

### Model-specific request options

`extra_body` merges bounded, JSON-serializable fields into the top level of the HTTP request body,
including provider-specific thinking options. It cannot override `model`, `messages`, `stream`,
token limits, output format, or authentication fields, and cannot enable unsupported tool or
multimodal requests. The model determines whether a parameter is valid; no model-specific options
are injected by default.

```python
config = LLMConfig(
    provider="qwen",
    model="a-model-that-supports-disabling-thinking",
    extra_body={"enable_thinking": False},
)
```

For Qwen models that require thinking to be disabled for non-streaming or JSON output, configure
the example above explicitly according to the model documentation. Some OpenAI/Kimi models accept
`extra_body={"reasoning_effort": "low"}`. Gemini's Google-specific extension must retain its HTTP
nesting, for example `extra_body={"extra_body": {"google": {"thinking_config": {"thinking_level": "low"}}}}`.
These illustrate parameter formats for particular models, not settings guaranteed to work across
all models from a provider.

### Text input and output

```python
import os

from dnakit.llm import LLMConfig, call_llm

config = LLMConfig(
    provider="deepseek",  # Or gemini, qwen, kimi, openrouter, openai.
    model=os.environ["DNAKIT_LLM_MODEL"],  # An actual model supported by your account.
    api_key_env="DNAKIT_LLM_API_KEY",  # Or api_key=os.environ["YOUR_KEY_VARIABLE"].
    max_output_tokens=1024,
    timeout_seconds=60,
)

result = call_llm(
    "Translate the following text into Chinese. Return only the translation.",
    input_data="The sample passed the quality check.",
    config=config,
)
translated_text = result.text
print(translated_text)
print(dict(result.usage))  # Provider-reported token usage, not an estimated cost.
```

Set the two environment variables in your PyCharm run configuration; there is no need to write
keys into the code. `base_url` is the API address, not an account page such as
`platform.deepseek.com/usage`. Configure an actual model currently available from the provider;
DNAKit does not select or switch models for you.

<span id="llm-file-input"></span>

### Explicit local-file paths

Both `call_llm()` and `acall_llm()` accept `file_path` to read and process a local text file directly,
without a separate `read_text()` call:

```python
from dnakit.llm import LLMConfig, call_llm

result = call_llm(
    "Summarize the file in Chinese, using only its contents.",
    file_path="/mnt/d/Pycharm_workspace/DNA/DNAKit/LICENSE",
    config=LLMConfig(
        provider="gemini",  # Other supported providers can also be used.
        model="your-actual-model-id",
        api_key_env="GEMINI_API_KEY",
    ),
)
print(result.text)
```

Paths can be strings or `Path` objects, absolute or relative to the current working directory,
and may contain spaces. Underscores do not need backslash escapes inside Python strings. Use
native Windows paths such as `r"D:\Pycharm_workspace\DNA\DNAKit\LICENSE"` with Windows Python;
`/mnt/d/...` is a WSL/Linux path.

**Specifying `file_path` authorizes sending the complete text of that file.** The request includes
the filename, but does not automatically include the absolute directory. It does not read nearby
files, scan directories, recursively follow links in the file, execute code, or automatically
split the input into multiple requests. `input_data` and `file_path` cannot be used together.
The original `input_data="/path/to/file"` behavior remains unchanged: it is a literal string and
does not trigger file reading.

You can preview the text entirely locally first:

```python
from dnakit.llm import read_llm_file

text = read_llm_file("LICENSE")  # Local read only, no network request.
print(text[:200])
```

The default `file_encoding="utf-8-sig"` supports UTF-8 with or without a BOM. Specify other encodings
explicitly, such as `file_encoding="gb18030"` or `"utf-16"`. The corresponding `read_llm_file()`
parameter is named `encoding`.

The default `max_file_bytes=64000` limits raw file size; the corresponding `read_llm_file()`
parameter is `max_bytes`. The complete request is also subject to `LLMConfig.max_request_bytes`,
which includes the prompt, filename, JSON escaping, and other request fields. Raising the file
limit does not raise the request limit. Oversized inputs fail before transmission rather than
being silently truncated, sampled, or split into multiple paid calls.

Supported inputs include extensionless LICENSE files, TXT, Markdown, CSV/TSV, JSON/JSONL, YAML,
XML, Python, and other decodable text. Tables and code are sent as text without execution or
automatic cleaning. PDF, binary Word/Excel documents, images, and archives cannot be processed
directly; extract their text separately and pass it through `input_data`. Empty files, directories,
device files, invalid encodings, and oversized files produce explicit errors. The complete example
is `others/examples/custom_llm_file.py`.

### JSON input and structured output

```python
rows = [
    {"id": "a", "note": "quality check passed"},
    {"id": "b", "note": "measurement needs repeating"},
]
result = call_llm(
    'Translate note into Chinese. Return JSON {"items": [{"id": "original ID", "note_zh": "translation"}]}.',
    input_data=rows,
    output_format="json",
    system_prompt="You translate experimental notes. Preserve original IDs and do not invent findings.",
    config=config,
)
output = result.data  # An ordinary dict; None in text mode.
if output is None or not isinstance(output.get("items"), list):
    raise ValueError("Expected an items array in the output")
# Continue validating record IDs and field types before passing data to the next Python step.
```

JSON mode requests `response_format={"type": "json_object"}` and validates JSON syntax, non-finite
values, duplicate keys, and the top-level object type. It does not guarantee that business fields
match expectations or that the model's judgment is correct. Put arrays inside an object, such as
`{"items": [...]}`. The selected model must support JSON mode. Provider rejections or invalid JSON
produce explicit errors; DNAKit does not silently switch output modes or models.

`LLMResult` provides `text`, `data`, `model`, `response_id`, `finish_reason`, `usage`, and `to_dict()`.
Each access to `data` returns an independent ordinary dictionary; modifying it does not change the
original result. The result does not retain the request prompt or original input, although model
output itself may still contain sensitive information.

### Arbitrary Python steps and asynchronous programs

Pass strings, finite numbers, booleans, lists, tuples, or string-keyed dictionaries from a previous
step through `input_data`. `None` means no additional data is attached; only the prompt is used.
Callers must explicitly convert DataFrames, NumPy arrays, and custom objects. For example, to share
only a table sample, use `df.head(5).to_dict(orient="records")` and check for missing and non-finite
values.

```python
from dnakit.llm import acall_llm


async def summarize_step(current_result, config):
    result = await acall_llm(
        "Summarize the result in two sentences without inferring information not provided.",
        input_data=current_result,
        config=config,
    )
    return result.text
```

`acall_llm()` runs the HTTP request in a worker thread without blocking the asyncio event loop.
Cancelling the await does not cancel an already submitted request; the provider may still complete
it and charge for it. There are no automatic retries or implicit conversation history on either
the initial connection or subsequent calls. Your Python code controls the number of looped or
concurrent calls.
