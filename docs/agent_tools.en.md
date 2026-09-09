# Agent and MCP tools

DNAKit can be called from MCP-capable Agents such as Codex, Claude, Cursor, and VS Code. The
Agent searches a compact catalog, inspects the selected input schema, and then calls the existing
public Python API.

This page also covers [general-purpose LLM calls](#custom-llm) within Python workflows, including [local-file input](#llm-file-input). These calls do not require an MCP server or the Agent extra.

## 1) Installation

Install DNAKit with the Agent extra:

```bash
python -m pip install "dnakit[agent]"
```

For a source checkout:

```bash
python -m pip install -e ".[agent]"
```

Deep-learning prediction, database access, and ToolUniverse-backed capabilities retain their own
optional dependencies. For example:

```bash
python -m pip install "dnakit[agent,external-tools,neural,neural-enformer]"
```

## 2) Start the server

```bash
dnakit-mcp
```

The command uses the local `stdio` MCP transport. It waits for an Agent client, so it does not
print ordinary terminal output when run directly.

## 3) Agent configuration

Add the following server to the Agent's MCP configuration:

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/.venv/bin/dnakit-mcp"
    }
  }
}
```

A Windows virtual environment normally uses:

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "D:\\project\\.venv\\Scripts\\dnakit-mcp.exe"
    }
  }
}
```

Use an absolute executable path visible to the Agent. Configuration file locations vary by client,
but the launch command is the same.

## 4) Compact MCP surface

The server exposes only six discovery and execution tools, rather than sending hundreds of schemas
to the model at once:

1. `dnakit_catalog_stats` reports catalog counts;
2. `list_dnakit_categories` lists capability groups;
3. `list_dnakit_tools` returns a bounded page;
4. `search_dnakit_tools` searches by capability;
5. `describe_dnakit_tool` returns the complete schema and safety metadata;
6. `call_dnakit_tool` executes one selected DNAKit function.

The current source discovers 326 public functions, of which 323 accept JSON-adaptable inputs. Two
low-level batch functions that require Python callbacks and the arbitrary command executor remain
Python-only. An Agent can perform ordinary batching by making repeated calls to a selected DNAKit
function.

## 5) Example

Ask the Agent:

```text
Use DNAKit to calculate the single-stranded molecular weight of ACGTACGT.
```

The expected tool sequence is:

```text
search_dnakit_tools("molecular weight")
describe_dnakit_tool("dnakit.thermodynamics.molecular_weight")
call_dnakit_tool(
    tool_name="dnakit.thermodynamics.molecular_weight",
    arguments={"sequence": "ACGTACGT", "strand": "single"}
)
```

The response is bounded JSON and retains DNAKit algorithm versions, applicability statements, and
provenance.

## 6) Inputs and safety

- A `DNASequence` accepts sequence text; a `DNARecord` accepts an object with `sequence` and `id`.
- Configuration classes use ordinary JSON objects whose fields are documented by the tool schema.
- The 27 local task checkpoints are callable through `dnakit.predictions.predict_enformer_benchmark`; it is a model tool, so an Agent still needs explicit confirmation and a checkpoint path.
- Python callbacks, open file handles, and live backend instances cannot be passed through JSON.
- File writes, model downloads, and external programs are denied by default and require an explicit
  `allow_side_effects=true` after review.
- Database and model tools remain subject to network access, APIs, weights, GPU availability,
  licenses, and optional dependencies.
- Arguments are limited to 2 MB and results to 5 MB by default; use file-output tools for large
  results.

## 7) General-purpose LLM calls (Python API) {#custom-llm}

`dnakit.llm` inserts a custom prompt/input/output step into any Python workflow, independently of DNA training or file-schema analysis. Use `LLMConfig`, `LLMResult`, `call_llm()` and `acall_llm()`; no OpenAI SDK, neural dependencies, or CLI arguments are required. Version one targets OpenAI-compatible Chat Completions endpoints, not provider-specific protocols.

### Provider presets

Select a provider and your actual model ID; the URL and key environment variable are filled automatically. Models are never selected or switched implicitly.

| Provider | `provider` | Default base URL | Default key environment variable |
| --- | --- | --- | --- |
| DeepSeek | `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| [Google Gemini](https://ai.google.dev/gemini-api/docs/openai) | `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| [Qwen / Model Studio](https://help.aliyun.com/zh/model-studio/qwen-api-via-openai-chat-completions) | `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| [Kimi / Moonshot](https://platform.kimi.com/docs/get-api-key) | `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| [OpenRouter](https://openrouter.ai/docs/quickstart) | `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| [OpenAI](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create) | `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| Custom compatible endpoint | `custom` | Explicit `base_url` required | `DNAKIT_LLM_API_KEY` |

```python
from dnakit.llm import LLMConfig, available_llm_providers, call_llm

print([provider.name for provider in available_llm_providers()])  # No network request.
config = LLMConfig(provider="gemini", model="your-actual-model-id")
result = call_llm("Summarize this text.", input_data="Your input text", config=config)
print(result.text)
```

`available_llm_providers()` returns immutable `LLMProvider` presets; `get_llm_provider(name)` resolves a preset. Aliases include `google`/`google-gemini`, `dashscope`/`alibaba`, `moonshot`, and `openai-compatible`. Unknown names fail before transmission rather than falling back to DeepSeek.

Explicit `base_url`, `api_key_env`, and `api_key` override presets. Omitting `provider` preserves the original DeepSeek URL and `DNAKIT_LLM_API_KEY` behavior. Both sync and async APIs use the same configuration.

Gemini uses Google's OpenAI-compatible endpoint and Gemini API key, not native `generateContent` or Vertex AI service-account authentication. Qwen's default is the still-supported Beijing legacy URL; override it with your region/workspace URL, such as `https://<WorkspaceId>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` or `https://<WorkspaceId>.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1`, with a matching regional key. International Kimi accounts can override `base_url="https://api.moonshot.ai/v1"`. OpenRouter expects a complete model ID from its catalog. A provider preset does not promise that every model supports every feature.

### Model-specific request options

`extra_body` merges bounded JSON fields into the HTTP request body without replacing core fields, authentication, or enabling unsupported tool/multimodal requests. For Qwen models that require non-thinking mode for synchronous/JSON output, configure `extra_body={"enable_thinking": False}` where supported. Some OpenAI/Kimi models accept `extra_body={"reasoning_effort": "low"}`. Gemini's Google-specific extension retains its HTTP nesting, for example `extra_body={"extra_body": {"google": {"thinking_config": {"thinking_level": "low"}}}}`. These are parameter-shape examples, not universal model settings; no model-specific options are injected by default.

### Text and JSON examples

```python
import os

from dnakit.llm import LLMConfig, call_llm

config = LLMConfig(
    provider="deepseek",  # Or gemini, qwen, kimi, openrouter, openai.
    model=os.environ["DNAKIT_LLM_MODEL"],
    api_key_env="DNAKIT_LLM_API_KEY",
    max_output_tokens=1024,
    timeout_seconds=60,
)
result = call_llm(
    "Translate the input into Chinese. Return only the translation.",
    input_data="The sample passed the quality check.",
    config=config,
)
print(result.text)
print(dict(result.usage))
```

Configure the two environment variables in your IDE, using an actual model identifier supported by your account. `base_url` must be the API address, not an account/usage webpage. Explicit `api_key=...` is also supported, but do not commit credentials to source control.

### Explicit local-file input {#llm-file-input}

Both `call_llm()` and `acall_llm()` accept `file_path`, so a local text file can be processed directly:

```python
from dnakit.llm import LLMConfig, call_llm, read_llm_file

print(read_llm_file("LICENSE")[:200])  # Local preview only, no API request.
result = call_llm(
    "Summarize the file using only its contents.",
    file_path="/mnt/d/Pycharm_workspace/DNA/DNAKit/LICENSE",
    config=LLMConfig(provider="gemini", model="your-actual-model-id"),
)
print(result.text)
```

Paths can be strings or `Path` objects, absolute or relative to the working directory. Underscores do not need escaping inside Python strings. Use native Windows paths such as `r"D:\Pycharm_workspace\DNA\DNAKit\LICENSE"` with Windows Python; `/mnt/d/...` is a WSL/Linux path.

**Using `file_path` authorizes sending that complete file's text**, together with its basename, but not its absolute directory. No neighboring files, folders, or referenced links are read. `file_path` and non-None `input_data` are mutually exclusive. Path-looking strings passed as `input_data` remain literal text, preserving the original privacy boundary.

The default `file_encoding="utf-8-sig"` handles UTF-8 with or without BOM; select other encodings explicitly, for example `"gb18030"` or `"utf-16"`. The default `max_file_bytes=64000` caps raw file bytes. `read_llm_file()` uses the corresponding names `encoding` and `max_bytes` and never calls an API. `LLMConfig.max_request_bytes` separately limits the complete serialized request, including JSON escaping, prompt, and filename. Oversized inputs fail before transmission rather than being silently truncated, sampled, or split into paid requests.

Text is accepted regardless of extension, including LICENSE, TXT, Markdown, CSV/TSV, JSON/JSONL, YAML, XML, and source code. Files are sent as text without automatic cleaning or execution. PDF, binary Office documents, images, and archives require separate text extraction before passing `input_data`. Empty files, directories, special devices, invalid encodings, and oversized files are explicit errors. The complete editable example is `others/examples/custom_llm_file.py`.

Use JSON mode for downstream Python data processing:

```python
result = call_llm(
    'Translate each note into Chinese. Return {"items": [{"id": "original ID", "note_zh": "text"}]}.',
    input_data=[{"id": "a", "note": "quality check passed"}],
    output_format="json",
    system_prompt="Translate faithfully. Preserve IDs and do not invent experimental findings.",
    config=config,
)
output = result.data
if output is None or not isinstance(output.get("items"), list):
    raise ValueError("Expected an items array")
# Validate IDs and field types before passing the output to your next Python step.
```

`LLMResult` exposes `text`, `data`, `model`, `response_id`, `finish_reason`, `usage`, and `to_dict()`. `data` returns a fresh ordinary dict in JSON mode and `None` in text mode. JSON mode requests `response_format={"type": "json_object"}` and rejects malformed JSON, duplicate keys, non-finite values, and non-object roots. It does not validate your business schema or factual correctness. Put arrays inside an object, such as `{"items": [...]}`. Select a model supporting JSON mode; provider rejections and invalid JSON cause explicit errors, with no silent output-mode fallback or model switching by DNAKit.

`input_data` may be text, finite numbers, booleans, lists/tuples, or string-keyed mappings. `None` means no additional inline input. Explicitly convert and sample tables/arrays first, for example `df.head(5).to_dict(orient="records")`, checking missing/non-finite values. Only the explicit `file_path` parameter reads a file. Paths passed as `input_data` are never implicitly opened; file handles, DataFrames, and arbitrary objects still require conversion.

For asyncio programs, use `await acall_llm(prompt, input_data=current_result, config=config)`. It runs the HTTP call in a worker thread without blocking the event loop. Cancelling the await does not cancel an already submitted provider request, which may still complete and be billed.

Important defaults and boundaries:

- `temperature=None` omits the parameter; supported models can use values from `0` to `2`.
- `max_output_tokens=2048` uses `token_limit_parameter="auto"`: the OpenAI preset or a direct `api.openai.com` URL selects `max_completion_tokens`, while other presets select `max_tokens`. Either field can be explicitly requested. Set `instruction_role="developer"` where required. The new/legacy token-field distinction follows the official OpenAI reference; actual parameter support remains model-specific.
- The default socket timeout is 60 seconds, the full serialized request limit is 100,000 bytes, and the response limit is 2,000,000 bytes. Oversized requests fail before transmission, without silently truncating the input.
- This generic API sends **all explicitly supplied prompt/input data** within the configured limit, including a file's text when `file_path` is used. Unlike training schema analysis, it does not sample automatically. There is no workspace scanning, provider Files API attachment upload, conversation history, prompt cache, or automatic persistence; file text is sent in the same Chat Completions request.
- Credentials are put only in the Authorization header, not config repr, result metadata, or error text. Remote endpoints require HTTPS; HTTP is allowed only for localhost. All redirects are refused. Outputs can still contain sensitive information supplied by the caller; do not log them casually or export credentials using `dataclasses.asdict(config)`.
- Errors such as `LLM_OUTPUT_TRUNCATED`, `LLM_REFUSED`, `INVALID_LLM_JSON_OUTPUT`, and `LLM_API_ERROR` are distinguishable through `DNAKitError.code`. Raw request/provider error text is not included. There are no automatic retries that might duplicate requests or charges.
- Version one supports text/JSON only, not image/audio attachments, streamed tokens, or automatic tool calls. Model-generated Python, SQL, and shell are never executed; prompts are not a security sandbox.

The editable IDE example is `others/examples/custom_llm.py`.
