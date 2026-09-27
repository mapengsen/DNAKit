# Engineering and extensions: parameters and notes

See the [tutorial](../../api/features/14_engineering.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials].

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

scripts, resource_path, describe, search and list_operations support local discovery or resource reads. run and run_support execute scripts which may access networks or write files; prepare their interpreters, scientific software and reference data. The existing strict YAML run_workflow is separate from these script runners.

## Read cells and outputs from an executed notebook {#dnakit.databases.read_executed_notebook-parameters}

`dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `data_folder` | string &#124; null | None | Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path. |
| `notebook_path` | string &#124; null | None | Direct path to a specific .ipynb file. Provide either this OR data_folder. |
| `search` | string | None | Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in &#x27;matching_cells&#x27;. Case-insensitive. Example: &#x27;upregulated,log2FoldChange&#x27; or &#x27;padj.*0\.05&#x27; |
| `regex` | boolean | None | Treat &#x27;search&#x27; as a Python regex instead of a comma-separated term list (default: false) |
| `max_output_chars` | integer | None | Truncate each cell source/output to this many characters (default 400) |
| `include_source` | boolean | None | Include cell source code in the response (default true). Set false to halve response size when you only need outputs. |
| `verbose` | bool | True | Progress flag. |


## Manage ENCODE API credentials {#dnakit.encode.manage_credentials-parameters}

`dnakit.encode.manage_credentials(action, access_key=None, secret_key=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `action` | Literal&#91;&#x27;store&#x27;, &#x27;check&#x27;, &#x27;clear&#x27;&#93; | Required | What to do: - &quot;store&quot;: Save new credentials (requires access_key and secret_key) - &quot;check&quot;: Check if credentials are configured - &quot;clear&quot;: Remove stored credentials |
| `access_key` | str &#124; None | None | Your ENCODE access key (only for action=&quot;store&quot;) |
| `secret_key` | str &#124; None | None | Your ENCODE secret key (only for action=&quot;store&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Locate and verify a bundled resource file {#dnakit.workflows.resource_path-parameters}

`dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `resource` | str | Required | Relative resource path within the skill; verified before use. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Run a bundled workflow and retrieve its exit code and output {#dnakit.workflows.run-parameters}

`dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `script` | str | Required | Relative .py, .R or .sh script path inside the bundle. |
| `arguments` | list&#91;str&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `workdir` | str &#124; Path &#124; None | None | Subprocess working directory; None uses the current directory. |
| `timeout_seconds` | float | 300.0 | Subprocess timeout in seconds. |
| `max_output_bytes` | int | 1000000 | Maximum captured subprocess output bytes. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Run a support script for a composite workflow {#dnakit.workflows.run_support-parameters}

`dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `name` | str | Required | Support workflow name; allowed values are described below. |
| `arguments` | list&#91;str&#93; &#124; None | None | Parameter mapping for database calls; CLI string list for workflows. |
| `workdir` | str &#124; Path &#124; None | None | Subprocess working directory; None uses the current directory. |
| `timeout_seconds` | float | 300 | Subprocess timeout in seconds. |
| `max_output_bytes` | int | 1000000 | Maximum captured subprocess output bytes. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## List bundled Python, R and shell scripts {#dnakit.workflows.scripts-parameters}

`dnakit.workflows.scripts(skill_id, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `skill_id` | str | Required | Bundled skill ID or registered legacy alias. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Render an existing Markdown report as a PDF {#dnakit.workflows.write_wes_report-parameters}

`dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)`

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
