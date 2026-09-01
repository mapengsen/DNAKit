# Agent and MCP tools

DNAKit can be called from MCP-capable Agents such as Codex, Claude, Cursor, and VS Code. The
Agent searches a compact catalog, inspects the selected input schema, and then calls the existing
public Python API.

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
python -m pip install "dnakit[agent,external-tools,neural]"
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

The current source discovers 322 public functions, of which 319 accept JSON-adaptable inputs. Two
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
- Python callbacks, open file handles, and live backend instances cannot be passed through JSON.
- File writes, model downloads, and external programs are denied by default and require an explicit
  `allow_side_effects=true` after review.
- Database and model tools remain subject to network access, APIs, weights, GPU availability,
  licenses, and optional dependencies.
- Arguments are limited to 2 MB and results to 5 MB by default; use file-output tools for large
  results.
