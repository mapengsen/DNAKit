# Pipeline-definition tools: WDL, CWL, Snakemake, nf-core and Galaxy

This group of functions is an optional-dependency adapter set in `dnakit.workflows`, using [miniwdl](https://github.com/chanzuckerberg/miniwdl), [cwltool](https://github.com/common-workflow-language/cwltool), [Snakemake](https://github.com/snakemake/snakemake), [nf-core tools](https://github.com/nf-core/tools) and [bioblend](https://github.com/galaxyproject/bioblend). They are pipeline-definition tools, not DNA algorithms: they parse and validate WDL, CWL, Snakefile and nf-core checkouts, list tasks, steps, rules or processes, plan a Snakemake dry run, and read tool, history and workflow lists from a Galaxy server through its API. They do not read or compute any sequence result, and they do not execute pipelines.

All dependencies are optional and are imported only when the matching function is called. A missing package raises an error with the install command; the extras are `bio-wdl`, `bio-cwl`, `bio-snakemake`, `bio-nf-core` and `bio-bioblend`.

**Example file:** `pipe.wdl`
```text
version 1.0
task Count {
  input { File bam }
  command <<< wc -l ~{bam} >>>
  output { Int lines = 1 }
}
workflow Pipeline {
  input {
    File reads
    String sample
  }
  call Count { input: bam = reads }
  output { Int lines = Count.lines }
}
```

## 1) `dnakit.workflows.parse_wdl_document()`

- **Purpose:** Parse a WDL document and list its tasks, inputs, outputs and workflow calls.
- **Method:** Read the file with miniwdl's `WDL.load()` and type-check it; no task is run. Inputs without a default are marked required. Tasks and calls are sorted by name.
- **API:** `dnakit.workflows.parse_wdl_document(source, *, max_bytes=1000000)`.
- **Input:** A WDL file path; `max_bytes` is the file size limit, and a larger file raises an error. Imported files are read by miniwdl.
- **Output:** A dictionary with `wdl_version`, `imports`, `tasks`, `workflow` and `provenance`.
- **Example code:**

```python
from dnakit.workflows import parse_wdl_document

result = parse_wdl_document("pipe.wdl")
print([task["name"] for task in result["tasks"]])
print([(item["name"], item["required"]) for item in result["workflow"]["inputs"]])
print(result["workflow"]["calls"])
```

- **Example output:**

```text
['Count']
[('reads', True), ('sample', True)]
[{'name': 'Count', 'task': 'Count'}]
```

## 2) `dnakit.workflows.check_wdl_inputs()`

- **Purpose:** Check a set of JSON inputs against the input declarations of a WDL workflow.
- **Method:** Compare the input names with the required inputs, then use miniwdl's `Value.from_json()` to check that each value can be built as the declared type. `File` values are checked as strings only; the file is not opened.
- **API:** `dnakit.workflows.check_wdl_inputs(source, inputs, *, max_bytes=1000000)`.
- **Input:** A WDL file path; `inputs` maps input names to values. Names may be `call.input` for call-level inputs.
- **Output:** `valid`, `missing` (required inputs that were not given), `unknown` (names the workflow does not declare), `invalid` (values of the wrong type) and `provenance`.
- **Example code:**

```python
from dnakit.workflows import check_wdl_inputs

result = check_wdl_inputs("pipe.wdl", {"reads": 3, "sample": "S1", "extra": 1})
print(result["valid"], result["missing"], result["unknown"])
print([item["name"] for item in result["invalid"]])
```

- **Example output:**

```text
False [] ['extra']
['reads']
```

## 3) `dnakit.workflows.read_cwl_workflow()`

- **Purpose:** Load a CWL Workflow with cwltool and list its inputs, steps and outputs.
- **Method:** Call cwltool's `load_tool()`, which validates the document with schema-salad; no step is executed. Steps are sorted by id because cwltool's internal order depends on `PYTHONHASHSEED`. Each step lists its input sources and output names.
- **API:** `dnakit.workflows.read_cwl_workflow(source, *, max_bytes=1000000)`.
- **Input:** A CWL workflow file path whose `class` is `Workflow`. Tool files referenced by `run:` are resolved relative to the workflow file.
- **Output:** A dictionary with `inputs`, `steps`, `outputs` and `provenance`. Each item in `steps` has `name`, `run`, `run_class`, `inputs` and `outputs`.
- **Example file:** `wf.cwl`
```text
cwlVersion: v1.2
class: Workflow
inputs:
  bam: File
outputs:
  lines:
    type: File
    outputSource: count/lines
steps:
  count:
    run: count.cwl
    in:
      bam: bam
    out: [lines]
```

- **Example file:** `count.cwl`
```text
cwlVersion: v1.2
class: CommandLineTool
baseCommand: [wc, -l]
inputs:
  bam:
    type: File
    inputBinding: {position: 1}
outputs:
  lines:
    type: stdout
stdout: count.txt
```

- **Example code:**

```python
from dnakit.workflows import read_cwl_workflow

result = read_cwl_workflow("wf.cwl")
step = result["steps"][0]
print(step["name"], step["run"], step["run_class"])
print(step["inputs"], result["outputs"][0]["sources"])
```

- **Example output:**

```text
count count.cwl CommandLineTool
[{'name': 'bam', 'sources': ['bam']}] ['count/lines']
```

## 4) `dnakit.workflows.read_snakefile()`

- **Purpose:** Parse a Snakefile and list its rules with their inputs and outputs.
- **Method:** Parse the file and its `include` directives through the Snakemake API; no rule is run. Parsing happens in a scratch directory, so no `.snakemake/` folder is written next to the source. Rules are sorted by name.
- **API:** `dnakit.workflows.read_snakefile(source, *, max_bytes=1000000)`.
- **Input:** A Snakefile path; `max_bytes` is the file size limit.
- **Output:** A dictionary with `default_target`, `rules` (each item has `name`, `inputs` and `outputs`) and `provenance`.
- **Example file:** `Snakefile`
```text
rule all:
    input:
        expand("counts/{sample}.txt", sample=["A", "B"])

rule align:
    input:
        "reads/{sample}.fastq"
    output:
        "bam/{sample}.bam"
    shell:
        "cp {input} {output}"

rule count:
    input:
        "bam/{sample}.bam"
    output:
        "counts/{sample}.txt"
    shell:
        "wc -l {input} > {output}"
```

- **Example code:**

```python
from dnakit.workflows import read_snakefile

result = read_snakefile("Snakefile")
print(result["default_target"])
for rule in result["rules"]:
    print(rule["name"], rule["inputs"], rule["outputs"])
```

- **Example output:**

```text
all
align ['reads/{sample}.fastq'] ['bam/{sample}.bam']
all ['counts/A.txt', 'counts/B.txt'] []
count ['bam/{sample}.bam'] ['counts/{sample}.txt']
```

## 5) `dnakit.workflows.plan_snakemake_jobs()`

- **Purpose:** Dry-run a Snakefile and list the jobs needed to reach the targets.
- **Method:** Build the job DAG with the Snakemake API and keep only the jobs that still have to run, so existing outputs are not listed. No rule is run and no lock is taken. `workdir` defaults to the directory of the Snakefile, where inputs are looked up; Snakemake writes a `.snakemake/` folder with a log there. Jobs are sorted by rule, outputs and inputs.
- **API:** `dnakit.workflows.plan_snakemake_jobs(source, *, targets=(), workdir=None, max_bytes=1000000, max_jobs=10000)`.
- **Input:** A Snakefile path; `targets` is a list of target files or rule names, and an empty list uses the default target; `max_jobs` limits the number of planned jobs.
- **Output:** A dictionary with `job_count`, `jobs_by_rule`, `jobs` (each item has `rule`, `inputs`, `outputs` and `wildcards`) and `provenance`.
- **Example code:**

```python
from pathlib import Path
from dnakit.workflows import plan_snakemake_jobs

Path("reads").mkdir(exist_ok=True)
Path("bam").mkdir(exist_ok=True)
for sample in ("A", "B"):
    Path(f"reads/{sample}.fastq").write_text("ACGT\n")
Path("bam/A.bam").write_text("bam\n")

result = plan_snakemake_jobs("Snakefile")
print(result["job_count"], result["jobs_by_rule"])
for job in result["jobs"]:
    print(job["rule"], job["outputs"])
```

- **Example output:**

```text
4 {'align': 1, 'all': 1, 'count': 2}
align ['bam/B.bam']
all []
count ['counts/A.txt']
count ['counts/B.txt']
```

## 6) `dnakit.workflows.list_nfcore_processes()`

- **Purpose:** List the Nextflow processes in an nf-core pipeline checkout.
- **Method:** Use nf_core's `NFCoreComponent` to parse the `main.nf` files under `modules/nf-core` and `modules/local`, and the flat `.nf` files directly in either folder. Everything is offline: GitHub is not contacted and no repository is cloned. Only process names and paths are returned; containers and channel shapes need the `nextflow` command, so they are not read.
- **API:** `dnakit.workflows.list_nfcore_processes(source, *, max_components=500, max_bytes=1000000)`.
- **Input:** A pipeline directory that contains a `modules/` subdirectory; `max_components` limits the number of components, and `max_bytes` limits the size of each `main.nf`.
- **Output:** A dictionary with `component_count`, `process_count`, `components` (each item has `name`, `kind`, `process` and `path`) and `provenance`. `kind` is `nf-core` or `local`.
- **Example code:**

```python
from pathlib import Path
from dnakit.workflows import list_nfcore_processes

Path("pipe/modules/nf-core/fastqc").mkdir(parents=True, exist_ok=True)
Path("pipe/modules/nf-core/fastqc/main.nf").write_text(
    'process FASTQC {\n    script:\n    """\n    fastqc\n    """\n}\n'
)
Path("pipe/modules/local").mkdir(parents=True, exist_ok=True)
Path("pipe/modules/local/count_lines.nf").write_text(
    'process COUNT_LINES {\n    script:\n    """\n    wc -l\n    """\n}\n'
)

result = list_nfcore_processes("pipe")
print(result["process_count"])
for row in result["components"]:
    print(row["kind"], row["name"], row["process"], row["path"])
```

- **Example output:**

```text
2
local count_lines COUNT_LINES modules/local/count_lines.nf
nf-core fastqc FASTQC modules/nf-core/fastqc/main.nf
```

## 7) `dnakit.workflows.list_galaxy_tools()`

- **Purpose:** List the tools installed on a Galaxy server.
- **Method:** Call bioblend's `GalaxyInstance.tools.get_tools()`, which requests `/api/tools`, and sort the result by `id`. The API key is sent in the `x-api-key` request header and never appears in the output. The `url` must not contain credentials, a query string or a fragment.
- **API:** `dnakit.workflows.list_galaxy_tools(url, *, api_key=None, name=None, max_items=1000)`.
- **Input:** A Galaxy server address (`http://` or `https://`); `api_key` is optional; `name` filters by tool name; `max_items` is the largest number of items accepted.
- **Output:** A dictionary with `server`, `count`, `items` (each item has `id`, `name`, `version` and `section`) and `provenance`.
- **Example code (not run locally):**

```python
import os
from dnakit.workflows import list_galaxy_tools

result = list_galaxy_tools(
    "https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"], name="Concatenate"
)
print(result["count"], [item["id"] for item in result["items"]])
```

- **Example output (not run locally):**

```text
(Not run locally: needs a reachable Galaxy server and an API key. Prints the number of matching tools and their IDs.)
```

## 8) `dnakit.workflows.list_galaxy_histories()`

- **Purpose:** List the non-deleted histories of the current Galaxy user.
- **Method:** Call bioblend's `GalaxyInstance.histories.get_histories()`, which requests `/api/histories`, and sort the result by `id`. Authentication works as in `list_galaxy_tools()`.
- **API:** `dnakit.workflows.list_galaxy_histories(url, *, api_key=None, name=None, max_items=1000)`.
- **Input:** As in `list_galaxy_tools()`; `name` filters by history name.
- **Output:** A dictionary with `server`, `count`, `items` (each item has `id`, `name`, `count`, `deleted` and `published`) and `provenance`.
- **Example code (not run locally):**

```python
import os
from dnakit.workflows import list_galaxy_histories

result = list_galaxy_histories("https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"])
print(result["count"], [item["name"] for item in result["items"]])
```

- **Example output (not run locally):**

```text
(Not run locally: needs a reachable Galaxy server and an API key. Prints the number of histories and their names.)
```

## 9) `dnakit.workflows.list_galaxy_workflows()`

- **Purpose:** List the workflows stored on a Galaxy server, which are Galaxy's pipeline definitions.
- **Method:** Call bioblend's `GalaxyInstance.workflows.get_workflows()`, which requests `/api/workflows`, and sort the result by `id`. The function only reads the list; it does not invoke a workflow.
- **API:** `dnakit.workflows.list_galaxy_workflows(url, *, api_key=None, name=None, max_items=1000)`.
- **Input:** As in `list_galaxy_tools()`; `name` filters by workflow name.
- **Output:** A dictionary with `server`, `count`, `items` (each item has `id`, `name` and `published`) and `provenance`.
- **Example code (not run locally):**

```python
import os
from dnakit.workflows import list_galaxy_workflows

result = list_galaxy_workflows("https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"])
print(result["count"], [item["name"] for item in result["items"]])
```

- **Example output (not run locally):**

```text
(Not run locally: needs a reachable Galaxy server and an API key. Prints the number of workflows and their names.)
```

## How these functions are verified

- The test file is `others/tests/unit/test_workflow_backends.py`, with 25 tests. In an environment with miniwdl, cwltool, snakemake, nf-core and bioblend installed, 24 pass; the live Galaxy test runs only when `DNAKIT_LIVE_TESTS=1` is set and is skipped by default. Without those packages, tests that need an upstream package skip, and the missing-package tests still run.
- The WDL and CWL tests write hand-written definitions into `tmp_path` and compare tasks, inputs, steps and sources with the expected values. The CWL workflow listing and the Snakemake plan are identical under different `PYTHONHASHSEED` values.
- The Snakemake dry-run test compares the planned job count (4) with a hand calculation, checks that the dry run creates no output files, and checks that missing inputs and unknown targets raise errors.
- The nf-core test builds a pipeline checkout in a temporary directory and runs offline.
- The Galaxy tests replace `requests.get` with a fake server, so bioblend's client code runs unchanged. They check that the API key appears only in the request header, that `max_items` and HTTP errors raise errors, and that no network is used.
- The results are checked against hand-written definitions; they are not compared line by line with the upstream command-line output.

## Current scope

- Wrapped: WDL parsing and input checks (miniwdl); CWL workflow loading and step listing (cwltool); Snakefile rule listing and dry-run job planning (Snakemake API); nf-core checkout process listing (nf_core, offline); Galaxy tool, history and workflow listing (bioblend).
- Not wrapped: pipeline execution (`miniwdl run`, cwltool execution, Snakemake execution, Nextflow runs). Galaxy workflow invocation, data upload and tool runs are not wrapped because they change server state. nf-core commands such as `pipelines create`, `lint` and `sync` are not wrapped either.
- Limits: only local files are supported, with no remote WDL imports and no URL-form CWL inputs. nf-core output contains process names and paths only. Snakemake dry runs use the upstream private attribute `WorkflowApi._workflow`, so `bio-snakemake` is limited to a version range. A dry run writes `.snakemake/` into the working directory but runs no rule.
