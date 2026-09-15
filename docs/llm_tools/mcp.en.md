# MCP

## Recommended: install with one prompt

First install and sign in to an MCP-capable Agent that can read documentation and run local commands,
such as Codex, Claude Code, or Gemini CLI. Paste this into the **Agent's chat**:

```text
Read https://mapengsen.github.io/DNAKit/setup.md, install and configure DNAKit MCP for me, and verify a tool call.
```

The Agent checks the environment, installs DNAKit, configures the client, and verifies a tool call.
Follow its instructions if sign-in or a client restart is needed. This example uses local computation
and default keyword search, so no additional API key is required.

## How it works

DNAKit can be called by MCP-capable Agents such as Codex, Claude, Cursor, and VS Code. Rather than
loading every Python function at once, the Agent searches the tool catalog, inspects parameter
schemas, and executes the selected public API.

This page uses local `stdio` transport: the MCP server is a program on the user's machine. After installing DNAKit and configuring
the launch command, the Agent starts that program and calls local functions. The DNAKit author does
not need to host the user's computation. The model makes decisions; computation, file access, and
model-weight loading take place in the server's environment. Tools can also call external APIs.

To use a model API such as DeepSeek, choose an Agent that supports both that API and MCP, or implement
an MCP client. Model-provider selection depends on the Agent. DNAKit's built-in `run_agent()` calls
the local registry directly without MCP; see [Python API](python_api.md).

The complete workflow is **install → configure the client → check the connection → submit a task**.

## 1) Installation

Use Python 3.10 or later, and install and sign in to your chosen Agent. Codex users can start with
the [official OpenAI CLI installation guide](https://learn.chatgpt.com/docs/codex/cli).

On Linux or macOS, install DNAKit and its MCP dependencies in a dedicated environment:

```bash
python3 -m venv "$HOME/.venvs/dnakit-mcp"
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install --upgrade pip
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install "dnakit[agent]"
```

`[agent]` installs MCP integration dependencies. Training, plotting, and deep-model tasks require
their optional dependencies, external programs, or model weights in the **same environment**;
see [Installation](../installation.md).

To test a local checkout, run the final installation command from the repository root with
`-e ".[agent]"` instead of `"dnakit[agent]"`. The client must be able to execute the configured path.

<span id="codex-setup"></span>

## 2) Connect and verify Codex

### 2.1 Register the local server

Run these commands in the environment where Codex is installed:

```bash
codex mcp add dnakit -- "$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
codex mcp list
```

Seeing `dnakit` in `codex mcp list` confirms registration. Check the actual connection and execution
as described below.

### 2.2 Configure timeouts

Open the Codex user configuration, `~/.codex/config.toml`. In the newly created `[mcp_servers.dnakit]`
table, add these settings immediately after `command`, or update them if already present:

```toml
startup_timeout_sec = 60
tool_timeout_sec = 3600
```

A 60-second startup timeout accommodates slower catalog imports. The example tool timeout allows
up to one hour per call; adjust it for your task. It does not guarantee completion within that time.
See the [official OpenAI MCP documentation](https://learn.chatgpt.com/docs/extend/mcp?surface=cli)
for configuration formats and options.

### 2.3 Start Codex and make the first call

Run this from your analysis working directory:

```bash
codex
```

Follow the sign-in prompt on first use. Restart an existing session after changing configuration.
In the Codex interactive interface, enter:

```text
/mcp
```

Check that `dnakit` is connected and exposes the eight entry points listed on this page, then ask:

> Use DNAKit MCP to compute the reverse complement of AACG. Actually call the tools and report the selected tool name and its result.

The expected flow is `search_dnakit_tools` → `describe_dnakit_tool` → `call_dnakit_tool`, selecting
`dnakit.ops.reverse_complement`. The sequence in the result should be **`CGTT`**. The result can
contain structured records, so do not require the entire response to be a plain string.
This example needs no deep-model weights. Default keyword search also needs no embedding or
additional DeepSeek API key; the Agent supplies the conversation model.

### 2.4 Codex IDE extension

The CLI and IDE extension share MCP configuration on the same Codex host. After configuring the
server above, restart the extension and check **MCP servers** in its gear menu. Alternatively,
choose **Add server**, select **STDIO**, enter the actual absolute `dnakit-mcp` path, save, and
restart the extension.

## 3) Other Agent clients

For clients that use the `mcpServers` JSON format, add this entry to the configuration file
specified by that client:

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/dnakit-mcp"
    }
  }
}
```

Replace `command` with the executable's absolute path in the installation environment. File location,
enablement, and restart procedures depend on the client. This JSON is not the format of Codex's
`config.toml`. After connecting, use the same `AACG` example to verify an actual call.

## 4) Server startup and execution environment

Once configured, **the Agent starts `dnakit-mcp` as a child process and exchanges MCP messages over
standard input/output**. Normally, no separate terminal is needed to start the server. Functions run
in the Python environment on the server's machine, using its files and computing resources.

Running the command manually is useful only for checking startup errors;
use the executable from the same installation environment:

```bash
"$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
```

It may print startup logs and wait for protocol input; it does not open a chat interface. A process
started in another terminal is not automatically adopted by the Agent. This entry point does not
listen on an HTTP port by default. A documentation-preview address such as
`http://127.0.0.1:8080/DNAKit/...` is a documentation URL, not a DNAKit MCP server URL.

## 5) MCP entry points available to the Agent

The MCP server exposes eight compact entry points, avoiding sending hundreds of tool schemas
at once:

1. `dnakit_catalog_stats`: view tool counts;
2. `list_dnakit_categories`: view capability categories;
3. `list_dnakit_tools`: list tools with pagination;
4. `search_dnakit_tools`: search by capability; keyword matching by default, with optional cloud semantic/hybrid retrieval;
5. `describe_dnakit_tool`: inspect the complete input schema and safety properties;
6. `call_dnakit_tool`: execute a specific DNAKit function.
7. `search_dnakit_skills`: search curated workflow metadata locally.
8. `read_dnakit_skill`: read a selected workflow or resource, with pagination.

See the [Skills tutorial](skills.md) for installation, dependencies, and citations. Skill discovery does not install upstream tools; only selected workflow text enters the conversation.

Discoverable and callable function counts depend on the installed version; use
`dnakit_catalog_stats` for the current counts. Interfaces requiring live Python objects, such as
callbacks, may remain Python-only. Agents can perform ordinary batch tasks by repeatedly calling
specific functions. Catalog function counts and the paper's descriptor/algorithm capability counts
use different counting units.

The catalog dynamically discovers functions exported through `__all__` in `PUBLIC_TOOL_MODULES`
and executes those same Python function objects; no separate Agent implementation is maintained.
After changing an existing implementation, signature, or docstring, restart Python/the MCP server
to load the changes in the same source environment. Running processes retain import and catalog
caches. If the Agent uses another virtual environment or a regular package installation, update
DNAKit in that environment. New functions must be publicly exported, and new modules must also be
added to the tool module list. Interfaces requiring callbacks or other live Python objects may
still be unavailable through JSON.

`input_schema.description` retains function documentation (up to 8,000 characters), so an LLM
receiving only the schema can still see parameter semantics. Search terms and execution arguments
are different: search for `download human hg38 reference genome fasta`, but pass `hg38` as
`download_genome`'s `query`. Do not copy search keywords directly into the assembly identifier.

<span id="tool-search-modes"></span>

## 6) Get candidates without executing a task

Clients can call only `search_dnakit_tools` to obtain candidates. For example, use these MCP `tools/call` parameters:

```json
{
  "name": "search_dnakit_tools",
  "arguments": {"query": "download human hg38 reference genome", "limit": 3}
}
```

The response includes `ok`, `count`, and `tools`; it does not download hg38 or execute candidates. Call `describe_dnakit_tool` if parameter documentation is needed. The client decides whether to invoke `call_dnakit_tool`. The default `mode="keyword"` needs no embedding API key and makes no chat LLM calls.

## 7) Optional cloud embedding and hybrid search

Configure an independent embedding service in the MCP server process to enable `semantic` and
`hybrid` modes. The LLM cannot specify server configuration or cache paths in tool arguments.
Replace the address, model, and credential below with actual provider settings, or supply them
through your client's environment/credential management.

For Codex, add an environment subtable to the existing server configuration:

```toml
[mcp_servers.dnakit.env]
DNAKIT_EMBEDDING_BASE_URL = "https://your-embedding-service.example/v1"
DNAKIT_EMBEDDING_MODEL = "your-text-embedding-model"
DNAKIT_EMBEDDING_API_KEY = "your-embedding-api-key"
```

For other clients using the `mcpServers` JSON format:

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

Startup does not generate embeddings. **Keyword search remains the default**, even with these variables configured. A client explicitly requests hybrid search as follows:

```json
{
  "name": "search_dnakit_tools",
  "arguments": {"query": "Download the human hg38 reference genome", "mode": "hybrid", "limit": 3}
}
```

Set `mode` to `semantic` for pure vector retrieval. Both modes return candidates only, but send query/tool descriptions to the embedding API and cache tool vectors on the machine running the MCP server. The optional `DNAKIT_EMBEDDING_CACHE_DIR` variable selects a cache directory on that machine. Model/API settings are read at server startup; restart after changing them. Missing configuration/credentials or API failures return `ok=false` with a specific error, without silently changing retrieval modes.

To customize batch size, dimensions, or timeouts, or disable disk caching, create the server in Python:

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

Python, the built-in Agent, and MCP share `registry.search()`. See [Python API: keyword, cloud embedding, and hybrid search](python_api.md#tool-search-modes) for ranking modes, cache updates, and other settings.

## 8) Troubleshooting

| Symptom | What to check |
| --- | --- |
| `dnakit-mcp` cannot be found | Install `dnakit[agent]` in the Python environment named by the configuration and use the executable's absolute path. |
| Registered server fails to connect | Check `/mcp` or the client startup log; verify the execution environment, paths, and `startup_timeout_sec`. |
| Manual startup keeps waiting | A local `stdio` server needs a client to exchange protocol messages. Configure the client to start it, then submit your task there. |
| Missing dependency, external program, or weights | Install the capability's dependencies in the MCP environment and provide the model or executable path requested by the error. |
| File writes or downloads are denied | State the requested operation and output path so the Agent can pass `allow_side_effects=true` to `call_dnakit_tool` within the task's authorization. |
| Code edits have no effect | Verify that the server environment uses the intended release or editable checkout, then restart the client's MCP server. |
