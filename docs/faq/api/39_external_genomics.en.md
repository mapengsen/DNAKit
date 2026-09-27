# External genomics analyses and MCP: parameters and notes

See the [feature guide](../../api/features/39_external_genomics.md) for purposes, methods and examples, and the [source FAQ](../genomics_sources.md) for papers, authors and license evidence.

## Discovery and catalog status {#discovery}

`search_projects` and `list_analyses` accept optional `query` and `domain` filters; all query terms must match. `describe_analysis(analysis)` accepts a registered analysis ID and returns input keys, parameter types and ranges, installation instructions, outputs and sources.

`reference_only` means a discovered project without a registered project-specific interface. For `available`, check whether the record represents a native function, executable adapter, MCP bridge or bundled skill. Inclusion does not mean all upstream methods were implemented or scientifically validated. Initial gap fields describe the repository before this expansion; `remaining_gaps` describes current limits.

## External analysis parameters and execution {#external-analysis}

Prepare the software, references, model weights and credentials required by the selected analysis separately. Harmonize reference builds, sample identifiers and allele orientation. Input prefixes must refer to complete file sets, such as PLINK BED/BIM/FAM.

| Parameter | Contract |
| --- | --- |
| `analysis` | Required; an ID returned by `list_analyses`. |
| `inputs` | Required; mapping with exactly the keys declared by `describe_analysis`, and file paths or prefixes as values. |
| `output_directory` | Required; a directory that does not yet exist. |
| `parameters` | Optional; registered parameters only. Defaults and ranges come from `describe_analysis`; arbitrary command-line arguments are not accepted. |
| `threads` | Default `1`; integer in `1–1024`. |
| `executable_paths` | Optional program-to-path overrides, such as `{"plink2": "/path/to/plink2"}`. |
| `timeout_seconds` | Execution only; default `3600.0`, range `(0, 86400]` seconds. |
| `max_artifact_bytes` | Execution only; default `1000000000`, limiting each artifact file in bytes. |
| `verbose` | Default `True`; controls runtime messages and progress display. |

`plan_analysis` validates inputs and builds commands without writing files; missing executables remain explicit in the plan. `run_analysis` requires the backend and a new output directory. It saves commands, inputs, parameters, sources, logs and `dnakit-run.json`, including failed runs. A `completed` process/file contract still requires scientific QC.

The four PLINK adapters were exercised with real PLINK 2 and deterministic synthetic data. Each additional adapter records its interface and synthetic-backend execution status in `validation`; see the [335-project review](../user_catalog.md). Valid empty FastANI/Barrnap/Cas-OFFinder results are explicitly allowed; other required empty outputs fail. `modkit-pileup` requires a pre-existing FASTA `.fai`.

## MCP installation and connection parameters {#mcp-connection}

Install MCP client support in a source checkout:

```bash
python -m pip install -e ".[agent]"
```

Install and configure each third-party service separately from its upstream instructions. `list_mcp_servers` returns installation and launch recipes without reading credentials or starting services; the feature guide retains the per-service catalog. The client does not install servers, and catalog inclusion does not establish local availability or scientific validation.

| Parameter | Contract |
| --- | --- |
| `server_id` | Required; a service catalog ID. |
| `tool_name` | Required for a single tool call; obtained from the live tool list. |
| `arguments` | Required for a single tool call; a JSON object matching the tool's `inputSchema`. |
| `transport` | Default `"stdio"`; `"http"` uses only cataloged author-documented endpoints. |
| `executable_path` | Optional path to the separately installed server executable or interpreter. |
| `workdir` | Optional working directory, commonly required for source checkouts. |
| `environment_variables` | Optional mapping from server variable names to existing local environment variable names, never credential values. |
| `timeout_seconds` | Default `60.0` seconds; request or complete workflow session timeout. |
| `max_response_bytes` | Default `4194304`; response limit after SDK decoding. |
| `max_tools` | Default `5000`; maximum tool list size. |
| `verbose` | Default `True`; controls runtime messages and progress display. |

`list_mcp_tools` starts a stdio service or connects to a registered HTTP service to obtain live schemas. `call_mcp_tool` validates arguments before invocation; remote `isError` is retained as `remote_error`, and known forwarded credential values are redacted. Each standalone call opens and closes one session, so in-memory state does not persist between calls. Response size is checked after SDK decoding, not as a wire-level traffic bound.

Tests exercise real local stdio discovery, argument validation, session state and failure handling; six author-hosted HTTP services also passed live tool discovery in the supplemental review. See the [review FAQ](../user_catalog.md) for the limits of scientific execution, account and model-environment validation.

## Calls within one session {#mcp-workflow}

`call_mcp_workflow` requires `server_id` and `calls`, with the same connection options as above. `calls` contains `1–100` entries, each with exactly `tool_name` and `arguments`. Arguments are explicit: no scripts run and no previous output is substituted automatically.

Schemas are refreshed and validated before every step. A remote error stops execution and returns all executed steps, including the failed step. Configuration and connection errors raise exceptions. Earlier uploads or computation are not rolled back. `timeout_seconds` bounds the entire session, and the response size limit applies to accumulated decoded results. Consult the selected tools for network, filesystem and compute effects.

## DNAKit MCP entry point {#dnakit-mcp}

After installing client support, start the unified server:

```bash
dnakit-mcp
```

Search `genomics` or `integrations` with `search_dnakit_tools`, then inspect full signatures and schemas with `describe_dnakit_tool`. External analysis and remote MCP calls have the `external_process` effect and use the existing `allow_side_effects=true` invocation mechanism.

## Full API signatures

::: dnakit.integrations
    options:
      show_source: false
