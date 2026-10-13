# Three-dimensional genome

This category contains **9** documented candidate capabilities; 9 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 9 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) aggregate

- **Purpose:** Document a reusable aggregate capability for Three-dimensional genome workflows. Candidate CAP-0252 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the aggregate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.genome3d.aggregate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `aggregate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import aggregate

result = aggregate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0252",
  "operation": "aggregate",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'aggregate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) annotate

- **Purpose:** Document a reusable annotate capability for Three-dimensional genome workflows. Candidate CAP-0253 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.genome3d.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import annotate

result = annotate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0253",
  "operation": "annotate",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'annotate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438948...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) calculate

- **Purpose:** Document a reusable calculate capability for Three-dimensional genome workflows. Candidate CAP-0254 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.genome3d.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import calculate

result = calculate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0254",
  "operation": "calculate",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'calculate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) compare

- **Purpose:** Document a reusable compare capability for Three-dimensional genome workflows. Candidate CAP-0255 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `compare|dnakit.similarity.compare`.

- **API:** `dnakit.genome3d.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import compare

result = compare(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0255",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'compare', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389481...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) export

- **Purpose:** Document a reusable export capability for Three-dimensional genome workflows. Candidate CAP-0256 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the export operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.genome3d.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `export` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import export

result = export(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0256",
  "operation": "export",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'export', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) normalize

- **Purpose:** Document a reusable normalize capability for Three-dimensional genome workflows. Candidate CAP-0258 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.standardize.normalize|normalize|vcf_normalize`.

- **API:** `dnakit.genome3d.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `normalize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import normalize

result = normalize(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0258",
  "operation": "normalize",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'normalize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) plot

- **Purpose:** Document a reusable plot capability for Three-dimensional genome workflows. Candidate CAP-0259 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.genome3d.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import plot

result = plot(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0259",
  "operation": "plot",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'plot', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) read

- **Purpose:** Document a reusable read capability for Three-dimensional genome workflows. Candidate CAP-0260 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.io.read|read|reads`.

- **API:** `dnakit.genome3d.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import read

result = read(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0260",
  "operation": "read",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'read', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) write

- **Purpose:** Document a reusable write capability for Three-dimensional genome workflows. Candidate CAP-0261 comes from 12 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (chiarseq, compartment, contact_matrix, distance_decay, enhancer_promoter_contact, hic, hichip, insulation, loop, microc, placseq, tad). The implementation is available at `src/dnakit/genome3d/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.io.write|write`.

- **API:** `dnakit.genome3d.write(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.genome3d import write

result = write(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0261",
  "operation": "write",
  "result_type": "dict",
  "result_preview": "{'module': 'genome3d', 'operation': 'write', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
