# Discover and call database APIs: parameters and notes

See the [tutorial](../../api/features/22_database_query_discovery.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

scripts, resource_path, describe, search and list_operations support local discovery or resource reads. run and run_support execute scripts which may access networks or write files; prepare their interpreters, scientific software and reference data. The existing strict YAML run_workflow is separate from these script runners.

## Call a database operation by name and parameters {#dnakit.databases.call-parameters}

`dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | str | Required | Native operation name or registered upstream alias. |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `max_result_bytes` | int | 20000000 | Serialized result byte limit; exceeding it raises an exception. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Inspect database input schemas and provenance {#dnakit.databases.describe-parameters}

`dnakit.databases.describe(operation, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | str | Required | Native operation name or registered upstream alias. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Execute a structured database request {#dnakit.databases.execute-parameters}

`dnakit.databases.execute(request, arguments=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `request` | str &#124; Mapping&#91;str, object&#93; | Required | Operation name or a mapping with name and arguments. |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Dispatch a database operation using a compatible tool name {#dnakit.databases.execute_tool-parameters}

`dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `tool_name` | string | Required | Name of the tool to execute (e.g., &#x27;example_tool_name&#x27;) |
| `arguments` | object &#124; None | None | Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {&quot;param1&quot;: &quot;value1&quot;, &quot;param2&quot;: 5}. JSON string example: &quot;{\&quot;param1\&quot;: \&quot;value1\&quot;, \&quot;param2\&quot;: 5}&quot;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;object&quot;, &quot;additionalProperties&quot;: true}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | Progress flag. |


## List database operations and aliases offline {#dnakit.databases.list_operations-parameters}

`dnakit.databases.list_operations(*, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Discover database operations by task keywords offline {#dnakit.databases.search-parameters}

`dnakit.databases.search(query, *, limit=10, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | str | Required | Local name/description search; use English task terms or a function name. |
| `limit` | int | 10 | Maximum matches; database catalog search allows 1–100. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |
