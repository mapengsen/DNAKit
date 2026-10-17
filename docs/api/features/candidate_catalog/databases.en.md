# Databases

This category documents **8** functions of the `databases` module. 5 of them (`cache`, `download`, `export`, `fetch`, `list_versions`) always return `status: "not_implemented"`: each would need a named database, a release catalogue or a cache layout that DNAKit does not define, and each message points to the working alternative. `query` calls one bundled database operation by name and `map_ids` maps gene identifiers through MyGene.info; both contact remote services, and `map_ids` returns `not_implemented` until the optional package `mygene` is installed. `search` lists matching bundled operations. Example outputs are real and were produced offline; `input_summary` is abbreviated to `{...}`.

## 1) cache

- **Purpose:** Not implemented: no database-result cache is defined for the database adapters.

- **Method:** No computation is performed. The database adapters keep no cache, so there is nothing for this function to populate or inspect.

- **API:** `dnakit.databases.cache(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a database-result cache with a defined layout and expiry policy") and a `message` (The database adapters keep no cache; there is nothing for this function to populate or inspect.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.databases import cache

result = cache("uniprot")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "cache",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a database-result cache with a defined layout and expiry policy",
  "message": "The database adapters keep no cache; there is nothing for this function to populate or inspect."
}
```

## 2) download

- **Purpose:** Not implemented: a generic download needs a named database release; use the specific download or fetch functions.

- **Method:** No computation is performed and nothing is downloaded. There is no generic download across the bundled databases; `dnakit.databases.search` finds the specific operation, which `query` can then call.

- **API:** `dnakit.databases.download(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a named database release and network access") and a `message` (No generic database download exists; use dnakit.databases.search to find the specific operation.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.databases import download

result = download({"database": "clinvar", "release": "2024-01"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "download",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a named database release and network access",
  "message": "No generic database download exists; use dnakit.databases.search to find the specific operation."
}
```

## 3) export

- **Purpose:** Not implemented: there is no generic export of database content.

- **Method:** No computation is performed and no file is written. To keep the result of a database call, pass the result of `query` to the functions of the `io_formats` module.

- **API:** `dnakit.databases.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a database snapshot and a target format") and a `message` (No generic database export exists; export the result of a query with the io functions.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.databases import export

result = export({"database": "clinvar"}, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "export",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a database snapshot and a target format",
  "message": "No generic database export exists; export the result of a query with the io functions."
}
```

## 4) fetch

- **Purpose:** Not implemented: `fetch` does not name a database; use `query` or the specific `fetch_*` functions.

- **Method:** No computation is performed and no request is sent. Call `dnakit.databases.query` with a bundled operation name instead, or one of the specific functions such as `fetch_ensembl_gene_info`.

- **API:** `dnakit.databases.fetch(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a named database and network access") and a `message` (fetch does not name a database; call dnakit.databases.query with an operation name instead.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.databases import fetch

result = fetch("ENSG00000012048")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "fetch",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a named database and network access",
  "message": "fetch does not name a database; call dnakit.databases.query with an operation name instead."
}
```

## 5) list_versions

- **Purpose:** Not implemented: DNAKit keeps no catalogue of database release versions.

- **Method:** No computation is performed. No DNAKit function lists database releases; the provider operations return whatever release the remote service currently serves.

- **API:** `dnakit.databases.list_versions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a release catalogue of the database") and a `message` (No DNAKit function lists database releases.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.databases import list_versions

result = list_versions("clinvar")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "list_versions",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a release catalogue of the database",
  "message": "No DNAKit function lists database releases."
}
```

## 6) map_ids

- **Purpose:** Map gene identifiers to MyGene.info records by delegating to `dnakit.databases.map_gene_identifiers`.

- **Method:** The identifiers are sent to the MyGene.info service, so the call needs network access and the optional package `mygene`. Without `mygene` the function returns `status: "not_implemented"` with `required_dependency: "pip install mygene"` and the import error as `message` (this is what the example shows). The identifiers are not mapped locally.

- **API:** `dnakit.databases.map_ids(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of identifiers or a comma-separated string (an empty list or another type raises `InputFormatError`). `options`: `scopes` (identifier type, default `"entrezgene"`), `species` (default `"human"`), `fields` (default `"symbol,name,taxid"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result` (`provenance.implementation` is `dnakit.databases.map_gene_identifiers`) when `mygene` is installed and the service answers; `result` is the `map_gene_identifiers` document with one row per identifier under `results`. Without `mygene` the envelope has `status: "not_implemented"`, `required_dependency` and `message` instead, as in the example (the example was run in an environment without `mygene`, and no request was sent).

- **Example code:**

```python
import json
from dnakit.databases import map_ids

result = map_ids("1017, 7157", options={"scopes": "entrezgene", "species": "human"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "databases",
  "operation": "map_ids",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install mygene",
  "message": "[BACKEND_UNAVAILABLE] mygene could not be imported: No module named 'mygene' Hint: Install the optional dependency with: python -m pip install \"dnakit[bio-mygene]\""
}
```

## 7) query

- **Purpose:** Call one bundled DNAKit database operation by name through `dnakit.databases.call`.

- **Method:** The operation name must be one of the bundled operations (see `dnakit.databases.list_operations` and `search`); the arguments are validated against that operation's input schema and the call is dispatched to the provider implementation, which contacts the remote database, so a successful call needs network access. An unknown operation raises `ConfigurationError` (`UNKNOWN_DATABASE_OPERATION`); a provider failure raises `BackendExecutionError` (`DATABASE_PROVIDER_ERROR`). `query` does not itself interpret or cache the returned document. No successful remote call is shown below because the example must run offline; it shows the error for an unknown operation.

- **API:** `dnakit.databases.query(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the operation name (arguments in `options["arguments"]`), or a mapping `{"name": operation, "arguments": {...}}` (other keys raise `InputFormatError`); `arguments` must be a mapping (`ConfigurationError` otherwise). `options`: `arguments`, `max_result_bytes` (default 20,000,000, bounds the result).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result` (`provenance.implementation` is `dnakit.databases.call`); `result` is the provider's JSON document, whose fields depend on the operation. The example prints the error of an unknown operation (`code` and `message` of the `ConfigurationError`) instead of an envelope.

- **Example code:**

```python
import json
from dnakit.databases import query
from dnakit.exceptions import ConfigurationError

try:
    query("definitely_unknown_operation")
except ConfigurationError as exc:
    print(json.dumps({"error": type(exc).__name__, "code": exc.code, "message": exc.message}, indent=2))
    # result = query({"name": "pdb_get_structure", "arguments": {"pdb_id": "1ABC"}})  # needs network
```

- **Example output:**

```json
{
  "error": "ConfigurationError",
  "code": "UNKNOWN_DATABASE_OPERATION",
  "message": "Unknown DNAKit database operation: definitely_unknown_operation"
}
```

## 8) search

- **Purpose:** List the bundled database and web-API operations whose name or description match a query. It searches the local catalogue of operations (currently 450), not the databases themselves, and makes no network request.

- **Method:** `dnakit.databases.search` is `dnakit.databases._runtime.search`, a plain keyword matcher and not part of the shared candidate runtime (the generic `search(input, options)` placeholder in `src/dnakit/databases/_candidates.py` is shadowed by it and is not the public function). The query is lower-cased and split into distinct word tokens (`\w+`). For each operation the text `name + upstream_name + description` is lower-cased, and the score is the number of query tokens that occur in it as substrings (so a short token such as `a` matches nearly everything). Operations with a score above 0 are sorted by score (descending), then by name, and the first `limit` are returned. There is no stemming, synonym handling or relevance model. Use `dnakit.databases.describe(name)` for the parameters of one operation and `query` or `call` to run it (these may contact the remote service).

- **API:** `dnakit.databases.search(query[required], limit[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `query`: non-empty string. `limit`: integer from 1 to 100, default 10 (booleans are rejected). `verbose`: print progress messages, default `True`. A wrong `query` or `limit` raises `ConfigurationError` (code `CONFIGURATION_ERROR`, message `Provide a nonempty query and limit between 1 and 100.`).

- **Output:** Unlike the candidate-runtime functions this does not return an envelope (no `status`, `provenance` or `input_summary`): it returns a list of mappings, best match first, each with `name` (DNAKit operation name), `description` and `upstream_name` (name of the original API operation). The list is empty when nothing matches.

- **Example code:**

```python
import json
from dnakit.databases import search

result = search("clinvar variant", limit=2)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
[
  {
    "name": "clin_var_get_clinical_significance",
    "description": "Get clinical significance information for a variant from ClinVar. Returns pathogenicity classification and clinical interpretations.",
    "upstream_name": "ClinVar_get_clinical_significance"
  },
  {
    "name": "clin_var_get_variant_details",
    "description": "Get variant summary information from ClinVar by variant ID. Returns accession, title, genes, clinical significance, review status, chromosome location, and variation name. For the same variant, this provides a different field view than ClinVar_get_clinical_significance.",
    "upstream_name": "ClinVar_get_variant_details"
  }
]
```
