# Engineering and extensions

Unified management of DNAKit's computing backend, Python API, CLI, configuration workflow, batch and parallel computing, cache, version and recurrence information.

## 1) Python API {#1-eng-003-python-api}

- **Purpose:** Call DNAKit functions through stable Python objects, functions, and result types, making it easy to combine analysis workflows in scripts, notebooks, tests, and other software.
- **Method:** Import the relevant module, call its function and inspect the returned object.
- **API:** The `dnakit` top-level and `dnakit.<domain>` modules are namespace entries and have no direct calling parameters; see the corresponding function entries for each function parameter.
- **Input:** `DNASequence`, `DNARecord`, `DNASet` or configuration object required by the corresponding function.
- **Output:** Each function returns its documented result type; here `length_features` returns `LengthResult`.
- **Example code:**

```python
import dnakit
from dnakit.descriptors import length_features

sequence = dnakit.normalize(" acgt ").sequence
assert sequence is not None
print(length_features(sequence).symbol_length)  # 4
```

- **Example output:**

```text
4
```

## 2) CLI {#2-eng-004-cli}

- **Purpose:** Run common functions such as standardization, description, fingerprinting, search, and comparison through command line parameters, and output structured JSON, suitable for shell processes and automation tasks.
- **Method:** Parse subcommand arguments, execute the corresponding function and print JSON.
- **API:** `dnakit COMMAND[required] ARGS[optional]`; the complete parameters of each subcommand are subject to `dnakit COMMAND --help`.
- **Input:** subcommand and its sequence, file or configuration parameters.
- **Output:** The CLI prints JSON to standard output; this example includes `distance`, `method`, and `mismatches`.
- **Example code:**

```bash
dnakit compare ACGT ACGA --method hamming
```

- **Example output:**

```text
{"costs": {"substitution": 1.0}, "distance": 1.0, "dp_cells": null, "edit_path": null, "exceeded_max_distance": false, "iupac_matching": "literal", "left_id": null, "left_length": 4, "max_cells": null, "max_distance": null, "method": "hamming", "mismatches": [{"left_symbol": "T", "position": 3, "right_symbol": "A"}], "name": "hamming_distance", "right_id": null, "right_length": 4}
```

## 3) Batch calculation {#3-eng-006-batch-calculation}

- **Purpose:** Repeatedly call the same calculation on a large number of DNA records, collecting successful results and errors one by one in the order of input, avoiding a single failure interrupting the entire batch.
- **Method:** Apply the operation to each record and collect results or errors in input order.
- **API:** `dnakit.batch.run_batch(records[required], operation[required], name[required], config[optional], progress[optional])`, `dnakit.batch.iter_batch(records[required], operation[required], config[optional], progress[optional])`; `config` uses `dnakit.batch.BatchConfig`.
- **Input:** `DNA([...])` or a `DNARecord` iterator, a batch callable, and optional error, seed, and resume configuration. `DNA` iterates over its records in order.
- **Output:** Returns `BatchResult`; `items` holds each value or error in input order.
- **Example code:**

```python
from dnakit import DNA
from dnakit.batch import run_batch

records = DNA(["A", "CC"])
result = run_batch(
    records,
    lambda record, context: record.sequence.symbol_length,
    name="length",
)
print([item.value for item in result.items])  # [1, 2]
```

- **Example output:**

```text
[1, 2]
```

## 4) Parallel calculation {#4-eng-007-parallel-calculation}

- **Purpose:** Use controlled worker threads to execute independent recording tasks in parallel, while keeping the results consistent with the input order, and uniformly propagating cancellation, progress and error information.
- **Method:** Run independent tasks in a thread pool bounded by jobs and max_in_flight.
- **API:** `dnakit.batch.run_batch(records[required], operation[required], name[required], config[optional], progress[optional])`; `config` uses `dnakit.batch.BatchConfig`, this setting is `execution_mode="thread"`.
- **Input:** `DNA([...])` or a record iterator, callable, worker count, and optional `max_in_flight`.
- **Output:** Returns `BatchResult` with `success_count` and per-record `items`.
- **Example code:**

```python
from dnakit import DNA
from dnakit.batch import BatchConfig, run_batch

records = DNA(["ACGT"] * 4)
result = run_batch(
    records,
    lambda record, context: (record.id, context.seed),
    name="threaded",
    config=BatchConfig(seed=7, jobs=2, execution_mode="thread", max_in_flight=2),
)
print(result.success_count)  # 4
```

- **Example output:**

```text
4
```

## 5) Chunking and streaming {#5-eng-008-chunking-and-streaming}

- **Purpose:** Read, convert and write large files in blocks using iterators, so that only the current batch is retained in the memory, and progress events are returned, which is suitable for data that exceeds the memory capacity.
- **Method:** Read an iterator and yield a chunk whenever chunk_size items are collected.
- **API:** `dnakit.read(source[required], format[optional], config[optional])`, `dnakit.RecordSource(iterator[required], close_callback[optional], source_name[optional], format[optional])`, `dnakit.io.iter_chunks(values[required], chunk_size[required])`, `dnakit convert input_path[required] output_path[required] --input-format[optional] --output-format[optional] --overwrite[optional] --progress/--no-progress[optional]`.
- **Input:** file path or iterator; optional format, compression and chunk size.
- **Output:** `iter_chunks` returns an iterator of tuples with at most `chunk_size` items each.
- **Example code:**

```python
from dnakit.io import iter_chunks

for chunk in iter_chunks(range(5), chunk_size=2):
    print(chunk)
# (0, 1), (2, 3), (4,)
```

- **Example output:**

```text
(0, 1)
(2, 3)
(4,)
```

## 6) cache {#6-eng-009-cache}

- **Purpose:** Generate cache keys based on input content, function names and parameters, reuse existing calculation results and verify integrity, reducing repeated high-cost calculations.
- **Method:** Derive keys from input and parameters, then read or write integrity-checked JSON entries.
- **API:** `dnakit.cache.CacheKey(namespace[required], digest[required], schema_version[optional])`, `dnakit.cache.CacheKey.from_components(namespace[required], components[required], schema_version[optional])`, `dnakit.cache.JSONCache(root[required], max_entry_bytes[optional])`, `dnakit.cache.JSONCache.get(key[required])`, `dnakit.cache.JSONCache.put(key[required], value[required])`, `dnakit.cache.JSONCache.invalidate(key[required])`, `dnakit.cache.JSONCache.clear(namespace[optional])`.
- **Input:** Private cache directory, namespace, input/parameters/algorithm and backend version components.
- **Output:** `get` returns a cached value or `None`; `put`, `invalidate`, and `clear` update the cache.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.cache import CacheKey, JSONCache

with TemporaryDirectory() as directory:
    cache = JSONCache(Path(directory) / "cache")
    key = CacheKey.from_components("length", {"sequence": "ACGT", "algorithm": "v1"})
    cache.put(key, {"value": 4})
    print(cache.get(key))  # {"value": 4}
```

- **Example output:**

```text
{'value': 4}
```

## 7) Errors and Warnings {#7-eng-012-errors-and-warnings}

- **Purpose:** Use structured exceptions and Issues to uniformly return error codes, context, severity and repair tips, so that CLI, API and reports can consistently explain failures or data problems.
- **Method:** Return error codes, context and repair hints in structured exceptions.
- **API:** `dnakit.exceptions.DNAKitError(message[required], code[optional], context[optional], hint[optional])`, `dnakit.core.Issue(code[required], severity[required], message[required], location[optional], details[optional])`; CLI non-zero exit code without calling parameters.
- **Input:** Illegal input, misconfigured or missing backend request.
- **Output:** Errors raise `DNAKitError` with a code, context, and hint; validation problems use `Issue`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.exceptions import DNAKitError

try:
    DNASequence("AX")
except DNAKitError as error:
    print(error.code)     # INVALID_ALPHABET
    print(error.context)
```

- **Example output:**

```text
INVALID_ALPHABET
{'alphabet': 'strict', 'part_index': 0, 'part_offset': 1, 'symbol': 'X'}
```

## 8) Agents and MCP {#8-eng-013-agents-and-mcp}

- **Purpose:** Convert stable public functions into searchable, inspectable, and callable MCP tools while reusing the existing DNAKit implementations.
- **Method:** Read tool schemas from the registry, validate arguments and dispatch to the function.
- **API:** `dnakit-mcp`, `dnakit.tools.default_tool_registry()`, and `dnakit.tools.create_server()`.
- **Input:** A tool name and JSON arguments matching the generated schema; file writes, model downloads, and external programs also require explicit authorization.
- **Output:** Returns tool catalog or MCP call results with tool names, parameters, or execution results.
- **Example code:**

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry()
result = registry.execute(
    "dnakit.thermodynamics.molecular_weight",
    {"sequence": "ACGT"},
)
print(result["value_dalton"])
```

- **Example output:**

```text
1173.84
```

See [LLM tools → MCP](../../llm_tools/mcp.md) for client configuration.

<span id="custom-llm"></span>

See [LLM tools → Python API](../../llm_tools/python_api.md#custom-llm) for general-purpose LLM calls, including custom prompts, file input, provider configuration, and safety boundaries. The [Agent](../../llm_tools/agent.md#agent-tool-calling) page covers natural-language tool execution without installing MCP.

<span id="bundled-script-example"></span>[Run a bundled script](#dnakit.workflows.run); [coordinate conversion example](../../faq/tutorial_integrations.md#engineering).

<span id="report-example"></span>[Render a PDF report](#dnakit.workflows.write_wes_report); [complete file example](../../faq/tutorial_integrations.md#engineering).

<span id="mcp-validation"></span>[Validation scope](../../faq/tutorial_integrations.md#engineering).

<!-- BEGIN GENERATED API: encode-credentials -->

<span id="encode-credentials"></span>

## 9) Manage ENCODE API credentials {#dnakit.encode.manage_credentials}

- **Purpose:** Manage ENCODE API credentials.
- **Method:** Read or update local ENCODE tracking records and format the selected operation's result.
- **API：** `dnakit.encode.manage_credentials(action, access_key=None, secret_key=None, *, verbose=True)`
- **Input:** Provide `action`.
- **Output:** Python return type: str.  JSON with action result.
- **Example code:**

```python
import dnakit

result = dnakit.encode.manage_credentials(
    action='check',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with action result.

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.encode.manage_credentials-parameters)

<!-- END GENERATED API: encode-credentials -->

<!-- BEGIN GENERATED API: execution -->

<span id="execution"></span>

## 10) Read cells and outputs from an executed notebook {#dnakit.databases.read_executed_notebook}

- **Purpose:** Read cells and outputs from an executed notebook.
- **Method:** Read saved cell outputs from a local executed notebook and apply the output-length limit.
- **API：** `dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Prepare an executed notebook with saved outputs at `analysis_executed.ipynb`, or supply your own path.
- **Output:** Returns a mapping with `notebook_path`, `total_cells`, `cells`, `search` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.read_executed_notebook(
    notebook_path='analysis_executed.ipynb',
    max_output_chars=2000,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "notebook_path": "<string>",
  "total_cells": "<integer>",
  "cells": [
    {
      "idx": "<integer>",
      "type": "<string>",
      "source": "<string>",
      "output": "<string>"
    }
  ],
  "search": "<string>",
  "matching_cells": [
    "<value>"
  ],
  "n_matches": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.databases.read_executed_notebook-parameters)

## 11) Locate and verify a bundled resource file {#dnakit.workflows.resource_path}

- **Purpose:** Locate and verify a bundled resource file.
- **Method:** Locate a bundled resource, verify its checksum and return its local path.
- **API：** `dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)`
- **Input:** Provide `skill_id`, `resource`.
- **Output:** pathlib.Path to a checksum-verified local resource.
- **Example code:**

```python
import dnakit

result = dnakit.workflows.resource_path(
    skill_id='genomic-coordinates',
    resource='scripts/convert_coords.py',
    verbose=False,
)
print(result)
```

- **Example output:** pathlib.Path to a checksum-verified local resource.

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.workflows.resource_path-parameters)

## 12) Run a bundled workflow and retrieve its exit code and output {#dnakit.workflows.run}

- **Purpose:** Run a bundled workflow and retrieve its exit code and output.
- **Method:** Execute the bundled script with the selected interpreter and collect its exit code, output and provenance.
- **API：** `dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)`
- **Input:** Provide `skill_id`, `script`.
- **Output:** Mapping with skill, script, returncode, output, upstream_url, upstream_revision and license; nonzero returncode indicates failure.
- **Example code:**

```python
import dnakit

result = dnakit.workflows.run(
    skill_id='genomic-coordinates',
    script='scripts/convert_coords.py',
    arguments=['--from', 'bed', '--to', 'gff', 'chr1', '999', '1000', '--format', 'json'],
    timeout_seconds=60,
    verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

- **Example output:** Mapping with skill, script, returncode, output, upstream_url, upstream_revision and license; nonzero returncode indicates failure.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "skill": "<str>",
  "script": "<str>",
  "returncode": "<int>",
  "output": "<str>",
  "upstream_url": "<str>"
}
```

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.workflows.run-parameters)

## 13) Run a support script for a composite workflow {#dnakit.workflows.run_support}

- **Purpose:** Run a support script for a composite workflow.
- **Method:** Execute the selected support script and collect its exit code and output.
- **API：** `dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)`
- **Input:** Provide `name`.
- **Output:** Mapping with name, returncode and output; check the exit code.
- **Example code:**

```python
import dnakit

result = dnakit.workflows.run_support(
    name='gwas-lookup',
    arguments=['--help'],
    verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

- **Example output:** Mapping with name, returncode and output; check the exit code.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "name": "<str>",
  "returncode": "<int>",
  "output": "<str>"
}
```

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.workflows.run_support-parameters)

## 14) List bundled Python, R and shell scripts {#dnakit.workflows.scripts}

- **Purpose:** List bundled Python, R and shell scripts.
- **Method:** Read bundled resources and list script paths for the selected skill.
- **API：** `dnakit.workflows.scripts(skill_id, *, verbose=True)`
- **Input:** Provide `skill_id`.
- **Output:** List of relative script paths; may be empty.
- **Example code:**

```python
import dnakit

result = dnakit.workflows.scripts(
    skill_id='genomic-coordinates',
    verbose=False,
)
print(result)
```

- **Example output:** List of relative script paths; may be empty.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<relative script path>"
]
```

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.workflows.scripts-parameters)

## 15) Render an existing Markdown report as a PDF {#dnakit.workflows.write_wes_report}

- **Purpose:** Render an existing Markdown report as a PDF.
- **Method:** Render the prepared Markdown as a PDF and write a provenance manifest.
- **API：** `dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)`
- **Input:** Provide `input_report`, `output_dir`. Prepare `report.md` and install `dnakit[reports]`; output is written to `reports/pdf`.
- **Output:** Returns a mapping of report paths and provenance; writes wes_report.pdf and wes_report.manifest.json. Requires dnakit[reports].
- **Example code:**

```python
import dnakit

result = dnakit.workflows.write_wes_report(
    input_report='report.md',
    output_dir='reports/pdf',
    verbose=False,
)
print(result)
```

- **Example output:** Returns a mapping of report paths and provenance; writes wes_report.pdf and wes_report.manifest.json. Requires dnakit[reports].

[Full parameters and usage notes](../../faq/api/14_engineering.md#dnakit.workflows.write_wes_report-parameters)

<!-- END GENERATED API: execution -->

## 16) Execute Python workflows by dependency {#bio-59}

- **Purpose:** Chain Python functions while retaining each step's status and result.
- **Method:** Validate dependencies and run in topological order; failures can block downstream steps while independent steps continue.
- **API:** `dnakit.bioinformatics.workflows.Step(name, function, dependencies={}, parameters={})`; `dnakit.bioinformatics.workflows.run_workflow(steps, fail_fast=True)`
- **Input:** Named steps, Python callables and parameter-name→upstream-step-name mappings.
- **Output:** Returns workflow results with each dependency-ordered step status and artifact.
- **Example code:**

```python
from dnakit.bioinformatics.workflows import Step, run_workflow

steps = [
    Step("sequence", lambda: "ACGT"),
    Step("length", lambda sequence: len(sequence), dependencies={"sequence": "sequence"}),
]
results = run_workflow(steps)
print([(step.name, step.status, step.value) for step in results])
```

- **Example output:**

```text
[('sequence', 'completed', 'ACGT'), ('length', 'completed', 4)]
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#workflows).
