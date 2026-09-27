---
name: dnakit-external-genomics
description: Discover reviewed DNA projects, run registered external genomics analyses and connect to registered scientific MCP servers through DNAKit. Use when a workflow needs an optional executable or external service.
license: MIT
---

# External genomic analyses and MCP

## Select an integration

Use `dnakit.integrations.search_projects(query)` for the dated literature/project survey, `list_analyses(query)` for executable analysis recipes, and `list_mcp_servers()` for supported MCP launch configurations. A project search hit may be research-only; inspect its integration status. Do not infer a callable API from a paper or invent a server endpoint.

## External software

1. Read `describe_analysis(analysis_id)` for exact input keys, required indexes, parameter schema, installation and original sources. Provision that software and reference data separately.
2. Call `plan_analysis(analysis_id, inputs=..., output_directory=..., parameters=...)`. Review the resolved command and expected artifacts. A plan creates no files and performs no backend analysis.
3. Use `run_analysis` with the same arguments and a new output directory. Threads, wall time and output bounds are explicit options. Failure manifests preserve status and logs; do not treat file existence as success. Inspect expected outputs and scientific QC after execution.

## Third-party MCP

Install `pip install "dnakit[agent]"` for the client, then separately install the selected server from its recorded upstream. Credentials are supplied by documented environment variable names. Some similarly named PyPI packages are unrelated: follow the catalogue's actual installation source.

Use `list_mcp_tools(server_id, ...)` to inspect live schemas before `call_mcp_tool(server_id, tool_name, arguments, ...)`. Use `call_mcp_workflow(server_id, calls, ...)` for stateful sequences in one connection, such as loading an oracle and then predicting. Each standalone call opens and closes its own session. A remote error stops the workflow; report the failed step. Call only tools needed for the user's authorized analysis and account for uploads or compute jobs described by their schema.

Use the user's supplied organism, reference build and input files. Report inputs, parameters, backend versions, outputs and validation actually performed. Catalogue inclusion is not evidence that a dependency is installed or that an analysis has run.

Cite the projects returned by `describe_analysis` / `list_mcp_servers` and the methods actually used. BOLE is a cited paper/service until its public source or documented API is available; this integration does not fabricate a BOLE server.
