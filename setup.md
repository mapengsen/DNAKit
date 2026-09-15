# Set up DNAKit MCP

This file is an installation guide for an Agent acting on a user's request to
install DNAKit MCP. Respond in the user's language. Carry out the setup using
your available local tools; do not ask the user to paste the installation prompt
again. If a required step needs user login or a client restart, explain exactly
what remains to be done.

DNAKit provides DNA sequence analysis tools through a local stdio MCP server.
The client starts the server as a child process, and DNAKit functions run in that
process's Python environment. The initial example uses keyword tool search and
local computation, with no additional API key, model weights, or external program.

## 1. Inspect the environment

- Identify the current MCP client and the environment where it runs. Ask which
  client to configure only if the session does not establish that.
- Check for Python 3.10 or later, an existing compatible DNAKit virtual
  environment, and an existing MCP server named `dnakit`.
- Reuse an environment that already meets the requirements. Preserve other MCP
  servers and user settings when merging or updating the `dnakit` configuration.
- Use the user's requested installation source. For a requested local checkout,
  verify that its `pyproject.toml` identifies the project as `dnakit`.

## 2. Install DNAKit

For a regular installation on Linux or macOS, use a dedicated environment:

```bash
python3 -m venv "$HOME/.venvs/dnakit-mcp"
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install --upgrade pip
"$HOME/.venvs/dnakit-mcp/bin/python" -m pip install "dnakit[agent]"
```

For a requested local checkout, run the installation from its repository root,
using the chosen environment's Python and `-m pip install -e ".[agent]"`.
Resolve the actual absolute path of that environment's `dnakit-mcp` executable.
The commands below assume the dedicated environment above; substitute the actual
path if reusing another environment.

Confirm the installed versions of `dnakit` and `fastmcp`, and that the executable
exists. If the installed release lacks the MCP entry point or required tools,
report the version mismatch rather than claiming installation succeeded.

## 3. Configure the current client

Choose only the section for the user's client. Register the server if missing,
or update its existing entry. Do not launch a separate background server: stdio
connections are managed by the client. The documentation website is not an MCP
HTTP endpoint.

### Codex

```bash
codex mcp add dnakit -- "$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
codex mcp list
```

In `~/.codex/config.toml`, add or update the following setting within the existing
`[mcp_servers.dnakit]` table, alongside `command`:

```toml
startup_timeout_sec = 60
```

Restart an existing Codex CLI session or IDE extension to load the configuration
if needed. In the CLI, `/mcp` shows the active servers. The IDE extension exposes
MCP servers in its settings. CLI and IDE configuration is shared on the same
Codex host.

Reference: https://learn.chatgpt.com/docs/extend/mcp?surface=cli

### Claude Code

```bash
claude mcp add --scope user --transport stdio dnakit -- "$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
```

When starting a new session from a terminal, allow 60 seconds for MCP startup:

```bash
MCP_TIMEOUT=60000 claude
```

Use `/mcp` inside the session to check the connection. If already working inside
Claude Code, have the user restart with this environment setting when required;
do not start a nested interactive session as a substitute for loading tools.

Reference: https://code.claude.com/docs/en/mcp

### Gemini CLI

```bash
gemini mcp add --scope user --transport stdio --timeout 60000 dnakit "$HOME/.venvs/dnakit-mcp/bin/dnakit-mcp"
gemini mcp list
```

Start or restart Gemini CLI with `gemini` and check `/mcp` in the session. Complete
any required project trust interaction in the client before testing local tools.

Reference: https://geminicli.com/docs/tools/mcp-server/

### Other MCP clients

Consult the client's official documentation for its local stdio configuration.
For clients using `mcpServers` JSON, merge an entry like this with the actual
absolute executable path:

```json
{
  "mcpServers": {
    "dnakit": {
      "command": "/absolute/path/to/dnakit-mcp"
    }
  }
}
```

This JSON is not the format of Codex's `config.toml`.

## 4. Verify an actual MCP tool call

Check that the connected server exposes these six entry points:

- `dnakit_catalog_stats`
- `list_dnakit_categories`
- `list_dnakit_tools`
- `search_dnakit_tools`
- `describe_dnakit_tool`
- `call_dnakit_tool`

Then perform the following through MCP:

1. Use `search_dnakit_tools` to search for `reverse complement`.
2. Use `describe_dnakit_tool` to inspect the actual input schema for
   `dnakit.ops.reverse_complement`.
3. Follow that schema to process the sequence `AACG` with `call_dnakit_tool`.
4. Confirm the call succeeded and the sequence in the returned records is `CGTT`.

Inspect the current schema instead of guessing arguments. A direct Python call
alone does not verify MCP connectivity, and the full result can be a structured
record rather than a plain string.

If the current session needs restarting to load the new tools, report
"configuration complete; tool verification pending restart" and give the user
this prompt in their language for the new session:

> Use DNAKit MCP to compute the reverse complement of AACG. Search for a suitable
> tool, inspect its parameters, actually call it, and report the tool name and result.

## 5. Report the outcome

Report the installation environment, executable path, configuration file changed,
and the actual verification result. Distinguish a registered configuration from
a connected server and a successful tool call. If something fails, inspect the
client's startup log, executable path, installed dependencies, and startup timeout,
then describe any remaining problem accurately.

Keep default keyword search for this initial setup. Configure cloud embeddings,
service credentials, optional dependencies, or model weights only when a user's
later task needs them.
