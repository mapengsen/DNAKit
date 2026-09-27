# Discover and call database APIs

Find and call database APIs by name or keyword.

<!-- BEGIN GENERATED API: discovery -->

<span id="discovery"></span>

## 1) Call a database operation by name and parameters {#dnakit.databases.call}

- **Purpose:** Call a database operation by name and parameters.
- **Method:** Validate the operation and arguments, then dispatch to its adapter and preserve the returned fields.
- **API：** `dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)`
- **Input:** Provide `operation`.
- **Output:** Returns the selected operation's response; this example retrieves the Ensembl BRCA1 record.
- **Example code:**

```python
import dnakit

result = dnakit.databases.call(
    operation='ensembl_lookup_gene',
    arguments={'id': 'BRCA1', 'species': 'homo_sapiens'},
    verbose=False,
)
print(result)
```

- **Example output:** Returns the selected operation's response; this example retrieves the Ensembl BRCA1 record.

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.call-parameters)

## 2) Inspect database input schemas and provenance {#dnakit.databases.describe}

- **Purpose:** Inspect database input schemas and provenance.
- **Method:** Read parameter constraints and provenance from the local operation catalog, without network access.
- **API：** `dnakit.databases.describe(operation, *, verbose=True)`
- **Input:** Provide `operation`.
- **Output:** Mapping with actual input_schema, provider config and upstream_name.
- **Example code:**

```python
import dnakit

result = dnakit.databases.describe(
    operation='gnomad_get_sv_by_region',
    verbose=False,
)
print(result)
```

- **Example output:** Mapping with actual input_schema, provider config and upstream_name.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "name": "<str>",
  "upstream_name": "<str>",
  "input_schema": "<dict>",
  "config": "<dict>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.describe-parameters)

## 3) Execute a structured database request {#dnakit.databases.execute}

- **Purpose:** Execute a structured database request.
- **Method:** Parse an operation name or name/arguments request and dispatch to the database adapter.
- **API：** `dnakit.databases.execute(request, arguments=None, *, verbose=True)`
- **Input:** Provide `request`.
- **Output:** Returns the requested operation's response; this example retrieves the Ensembl BRCA1 record.
- **Example code:**

```python
import dnakit

example_request = {'name': 'ensembl_lookup_gene',
 'arguments': {'id': 'BRCA1', 'species': 'homo_sapiens'}}

result = dnakit.databases.execute(
    request=example_request,
    verbose=False,
)
print(result)
```

- **Example output:** Returns the requested operation's response; this example retrieves the Ensembl BRCA1 record.

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.execute-parameters)

## 4) Dispatch a database operation using a compatible tool name {#dnakit.databases.execute_tool}

- **Purpose:** Dispatch a tool name and arguments to a local DNAKit database operation; recursive self-calls are rejected.
- **Method:** Use `Execute` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)`
- **Input:** Provide `tool_name`.
- **Output:** Returns the result of the operation selected by tool_name, using that operation's fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.execute_tool(
    tool_name='list_tools',
    arguments={'mode': 'names', 'limit': '3'},
    verbose=False,
)
print(result)
```

- **Example output:** Returns the result of the operation selected by tool_name, using that operation's fields.

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.execute_tool-parameters)

## 5) List database operations and aliases offline {#dnakit.databases.list_operations}

- **Purpose:** List database operations and aliases offline.
- **Method:** List function names and upstream aliases from the local registry.
- **API：** `dnakit.databases.list_operations(*, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** List of name and upstream_name mappings.
- **Example code:**

```python
import dnakit

result = dnakit.databases.list_operations(
    verbose=False,
)
print(result)
```

- **Example output:** List of name and upstream_name mappings.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "name": "<str>",
    "upstream_name": "<str>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.list_operations-parameters)

## 6) Discover database operations by task keywords offline {#dnakit.databases.search}

- **Purpose:** Discover database operations by task keywords offline.
- **Method:** Match keywords against local operation names and descriptions and rank by match count.
- **API：** `dnakit.databases.search(query, *, limit=10, verbose=True)`
- **Input:** Provide `query`.
- **Output:** List of name, description and upstream_name mappings.
- **Example code:**

```python
import dnakit

result = dnakit.databases.search(
    query='gnomad',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** List of name, description and upstream_name mappings.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "name": "<str>",
    "description": "<str>",
    "upstream_name": "<str>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_discovery.md#dnakit.databases.search-parameters)

<!-- END GENERATED API: discovery -->
