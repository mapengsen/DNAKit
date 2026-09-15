# API discovery, workflows and reports

This page documents **12 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

scripts, resource_path, describe, search and list_operations support local discovery or resource reads. run and run_support execute scripts which may access networks or write files; prepare their interpreters, scientific software and reference data. The existing strict YAML run_workflow is separate from these script runners.

## Function index

- [`dnakit.databases.call`](#dnakit.databases.call)
- [`dnakit.databases.describe`](#dnakit.databases.describe)
- [`dnakit.databases.execute`](#dnakit.databases.execute)
- [`dnakit.databases.execute_tool`](#dnakit.databases.execute_tool)
- [`dnakit.databases.list_operations`](#dnakit.databases.list_operations)
- [`dnakit.databases.read_executed_notebook`](#dnakit.databases.read_executed_notebook)
- [`dnakit.databases.search`](#dnakit.databases.search)
- [`dnakit.workflows.resource_path`](#dnakit.workflows.resource_path)
- [`dnakit.workflows.run`](#dnakit.workflows.run)
- [`dnakit.workflows.run_support`](#dnakit.workflows.run_support)
- [`dnakit.workflows.scripts`](#dnakit.workflows.scripts)
- [`dnakit.workflows.write_wes_report`](#dnakit.workflows.write_wes_report)

## `dnakit.databases.call` {#dnakit.databases.call}

Call a bundled DNAKit database adapter without installing ToolUniverse.

An operation may use its native snake_case name or recorded upstream name.
Input JSON is checked against the pinned service schema; returned provider
envelopes preserve the service's fields, pagination and evidence details.
Failure raises an exception, including provider error envelopes. Models,
reference data, API credentials and specialist libraries remain external.

```text
dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | str | Required | Native operation name or registered upstream alias. |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `max_result_bytes` | int | 20000000 | Serialized result byte limit; exceeding it raises an exception. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Returns object with the called operation's result structure; a common data envelope is not guaranteed. Provider errors raise exceptions.

## `dnakit.databases.describe` {#dnakit.databases.describe}

Return the function's actual upstream API schema and implementation provenance.

```text
dnakit.databases.describe(operation, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | str | Required | Native operation name or registered upstream alias. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Mapping with actual input_schema, provider config and upstream_name.

## `dnakit.databases.execute` {#dnakit.databases.execute}

Execute an operation name plus arguments, or a JSON name/arguments request.

```text
dnakit.databases.execute(request, arguments=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `request` | str &#124; Mapping&#91;str, object&#93; | Required | Operation name or a mapping with name and arguments. |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Returns object with the called operation's result structure; a common data envelope is not guaranteed. Provider errors raise exceptions.

## `dnakit.databases.execute_tool` {#dnakit.databases.execute_tool}

Dispatch a tool name and arguments to a local DNAKit database operation; recursive self-calls are rejected.

```text
dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `tool_name` | string | Required by input schema | Name of the tool to execute (e.g., &#x27;example_tool_name&#x27;) |
| `arguments` | object &#124; None | None | Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {&quot;param1&quot;: &quot;value1&quot;, &quot;param2&quot;: 5}. JSON string example: &quot;{\&quot;param1\&quot;: \&quot;value1\&quot;, \&quot;param2\&quot;: 5}&quot;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;object&quot;, &quot;additionalProperties&quot;: true}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "tool_name": {
          "type": "string",
          "description": "Name of the tool to execute (e.g., 'example_tool_name')"
        },
        "arguments": {
          "oneOf": [
            {
              "type": "object",
              "additionalProperties": true
            },
            {
              "type": "string"
            }
          ],
          "description": "Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {\"param1\": \"value1\", \"param2\": 5}. JSON string example: \"{\\\"param1\\\": \\\"value1\\\", \\\"param2\\\": 5}\"."
        }
      },
      "required": [
        "tool_name"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Tool execution result as a structured object (dict). The actual content depends on the executed tool. May contain error information if execution failed. Returns structured data instead of JSON string for programmatic use.",
      "additionalProperties": true
    }
    ```

**Execution and provenance:** `network`；upstream name `execute_tool`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("execute_tool", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.list_operations` {#dnakit.databases.list_operations}

List native operation names and recorded upstream aliases without loading providers.

```text
dnakit.databases.list_operations(*, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

List of name and upstream_name mappings.

## `dnakit.databases.read_executed_notebook` {#dnakit.databases.read_executed_notebook}

Read a pre-executed Jupyter notebook (.ipynb, often named '*_executed.ipynb') and return its cells with printed outputs as structured JSON. Use this BEFORE reimplementing an analysis when the data folder ships with an executed notebook — the notebook already ran the authoritative pipeline with the exact package versions, filters, and thresholds that match the reference answers. Reimplementing with pydeseq2/scanpy/etc. produces different numbers. Accepts either a folder (auto-discovers the notebook) or a direct path. Optionally filters cells by a comma-separated search term list or a regex so the agent sees only the cells relevant to its question.

```text
dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `data_folder` | string &#124; null | None | Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path. |
| `notebook_path` | string &#124; null | None | Direct path to a specific .ipynb file. Provide either this OR data_folder. |
| `search` | string | None | Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in &#x27;matching_cells&#x27;. Case-insensitive. Example: &#x27;upregulated,log2FoldChange&#x27; or &#x27;padj.*0\.05&#x27; |
| `regex` | boolean | None | Treat &#x27;search&#x27; as a Python regex instead of a comma-separated term list (default: false) |
| `max_output_chars` | integer | None | Truncate each cell source/output to this many characters (default 400) |
| `include_source` | boolean | None | Include cell source code in the response (default true). Set false to halve response size when you only need outputs. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "data_folder": {
          "type": [
            "string",
            "null"
          ],
          "description": "Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path."
        },
        "notebook_path": {
          "type": [
            "string",
            "null"
          ],
          "description": "Direct path to a specific .ipynb file. Provide either this OR data_folder."
        },
        "search": {
          "type": "string",
          "description": "Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in 'matching_cells'. Case-insensitive. Example: 'upregulated,log2FoldChange' or 'padj.*0\\.05'"
        },
        "regex": {
          "type": "boolean",
          "description": "Treat 'search' as a Python regex instead of a comma-separated term list (default: false)"
        },
        "max_output_chars": {
          "type": "integer",
          "description": "Truncate each cell source/output to this many characters (default 400)"
        },
        "include_source": {
          "type": "boolean",
          "description": "Include cell source code in the response (default true). Set false to halve response size when you only need outputs."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "notebook_path": {
              "type": "string"
            },
            "total_cells": {
              "type": "integer"
            },
            "cells": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "idx": {
                    "type": "integer"
                  },
                  "type": {
                    "type": "string"
                  },
                  "source": {
                    "type": "string"
                  },
                  "output": {
                    "type": "string"
                  }
                }
              }
            },
            "search": {
              "type": "string"
            },
            "matching_cells": {
              "type": "array"
            },
            "n_matches": {
              "type": "integer"
            }
          },
          "required": [
            "notebook_path",
            "total_cells",
            "cells"
          ]
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `read_executed_notebook`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("read_executed_notebook", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.search` {#dnakit.databases.search}

Find bundled scientific API operations locally by task, name or description.

```text
dnakit.databases.search(query, *, limit=10, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | str | Required | Local name/description search; use English task terms or a function name. |
| `limit` | int | 10 | Maximum matches; database catalog search allows 1–100. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

List of name, description and upstream_name mappings.

## `dnakit.workflows.resource_path` {#dnakit.workflows.resource_path}

Locate a checksum-verified bundled file, including resources of dependent skills.

```text
dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `resource` | str | Required | Relative resource path within the skill; verified before use. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

pathlib.Path to a checksum-verified local resource.

## `dnakit.workflows.run` {#dnakit.workflows.run}

Run a bundled workflow with its original CLI options and true scientific backends.

The caller chooses an explicit script and working directory. Python uses the
current environment; R, shell tools, trained models and data retain their own
requirements. This executes code and may write analysis outputs or access APIs.
Nonzero exit status is returned, not interpreted as scientific success.

```text
dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `script` | str | Required | Relative .py, .R or .sh script path inside the bundle. |
| `arguments` | list&#91;str&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `workdir` | str &#124; Path &#124; None | None | Subprocess working directory; None uses the current directory. |
| `timeout_seconds` | float | 300.0 | Subprocess timeout in seconds. |
| `max_output_bytes` | int | 1000000 | Maximum captured subprocess output bytes. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Mapping with skill, script, returncode, output, upstream_url, upstream_revision and license; nonzero returncode indicates failure.

## `dnakit.workflows.run_support` {#dnakit.workflows.run_support}

Run a licensed support stage of a composite skill with its documented CLI options.

Names: clinical-variant-reporter, clinical-trial-finder, omics-target-evidence-mapper,
gwas-lookup, fine-mapping.
Source implementations are bundled in DNAKit; scientific dependencies remain external.

```text
dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `name` | str | Required | Support workflow name; allowed values are described below. |
| `arguments` | list&#91;str&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `workdir` | str &#124; Path &#124; None | None | Subprocess working directory; None uses the current directory. |
| `timeout_seconds` | float | 300 | Subprocess timeout in seconds. |
| `max_output_bytes` | int | 1000000 | Maximum captured subprocess output bytes. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Mapping with name, returncode and output; check the exit code.

## `dnakit.workflows.scripts` {#dnakit.workflows.scripts}

List bundled Python, R and shell resources; consult their usage before execution.

```text
dnakit.workflows.scripts(skill_id, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

List of relative script paths; may be empty.

## `dnakit.workflows.write_wes_report` {#dnakit.workflows.write_wes_report}

Render supplied Markdown sections and tables as a PDF, with a source hash.

language selects English or Spanish report labels; input content is preserved
without translation or new clinical interpretations. Optional logos and a TTF
font are supported. Requires the reportlab rendering backend (dnakit[reports]).
This implementation was written for DNAKit; no proprietary report code is used.

```text
dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `input_report` | str &#124; Path | Required | Existing UTF-8 Markdown report path, at most 5 MB. |
| `output_dir` | str &#124; Path | Required | Report output directory for the PDF and provenance manifest. |
| `language` | str | &#x27;en&#x27; | Report label language en or es; input text is not translated. |
| `title` | str &#124; None | None | Custom report title; None uses the language default. |
| `logo_left` | str &#124; Path &#124; None | None | Left logo image path. |
| `logo_right` | str &#124; Path &#124; None | None | Right logo image path. |
| `font_path` | str &#124; Path &#124; None | None | Optional TTF font path. |
| `overwrite` | bool | False | Allow overwriting existing report files. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Returns a mapping of report paths and provenance; writes wes_report.pdf and wes_report.manifest.json. Requires dnakit[reports].

